// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/track_parametrization.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/edm/track_state_collection.hpp"
#include "traccc/edm/track_state_helpers.hpp"
#include "traccc/finding/measurement_selector.hpp"
#include "traccc/fitting/details/cholesky_solve.hpp"
#include "traccc/fitting/kalman_filter/two_filters_smoother.hpp"
#include "traccc/fitting/status_codes.hpp"

// Detray include(s).
#include <detray/definitions/geometry.hpp>
#include <detray/geometry/identifier.hpp>

// VecMem include(s).
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <random>
#include <utility>

namespace {

using algebra_t = traccc::default_algebra;
using scalar_t = traccc::scalar;
using bound_vector_t = traccc::bound_vector<algebra_t>;
using bound_matrix_t = traccc::bound_matrix<algebra_t>;
using bound_params_t = traccc::bound_track_parameters<algebra_t>;

constexpr std::size_t N{traccc::e_bound_size};

using dmat = std::array<std::array<double, N>, N>;
using dvec = std::array<double, N>;

/// Typical standard deviations of the bound parameters (loc0, loc1, phi,
/// theta, q/p, time): The covariances span many orders of magnitude
constexpr dvec sigmas{0.1, 0.2, 1e-3, 2e-3, 1e-3, 1.};

/// Invert a matrix in double precision (Gauss-Jordan, full pivoting is not
/// needed for the well-conditioned test matrices)
dmat invert(dmat a) {
  dmat inv{};
  for (std::size_t i = 0; i < N; ++i) {
    inv[i][i] = 1.;
  }
  for (std::size_t c = 0; c < N; ++c) {
    // Partial pivoting
    std::size_t p = c;
    for (std::size_t r = c + 1; r < N; ++r) {
      if (std::abs(a[r][c]) > std::abs(a[p][c])) {
        p = r;
      }
    }
    std::swap(a[c], a[p]);
    std::swap(inv[c], inv[p]);
    const double d = a[c][c];
    for (std::size_t k = 0; k < N; ++k) {
      a[c][k] /= d;
      inv[c][k] /= d;
    }
    for (std::size_t r = 0; r < N; ++r) {
      if (r != c) {
        const double f = a[r][c];
        for (std::size_t k = 0; k < N; ++k) {
          a[r][k] -= f * a[c][k];
          inv[r][k] -= f * inv[c][k];
        }
      }
    }
  }
  return inv;
}

dmat add(const dmat& a, const dmat& b) {
  dmat r{};
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < N; ++j) {
      r[i][j] = a[i][j] + b[i][j];
    }
  }
  return r;
}

dvec mult(const dmat& a, const dvec& v) {
  dvec r{};
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < N; ++j) {
      r[i] += a[i][j] * v[j];
    }
  }
  return r;
}

/// Random symmetric positive definite covariance with realistic scales and
/// correlations
dmat random_covariance(std::mt19937& gen, const double corr = 0.5) {
  std::uniform_real_distribution<double> dist(-1., 1.);
  dmat b{};
  for (auto& row : b) {
    for (auto& e : row) {
      e = dist(gen);
    }
  }
  dmat c{};
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < N; ++j) {
      double s = 0.;
      for (std::size_t k = 0; k < N; ++k) {
        s += b[i][k] * b[j][k];
      }
      c[i][j] = corr * s / static_cast<double>(N) + (i == j ? 1. : 0.);
      c[i][j] *= sigmas[i] * sigmas[j];
    }
  }
  return c;
}

bound_matrix_t to_float(const dmat& m) {
  bound_matrix_t r = traccc::matrix::zero<bound_matrix_t>();
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < N; ++j) {
      traccc::getter::element(r, i, j) = static_cast<scalar_t>(m[i][j]);
    }
  }
  return r;
}

bound_vector_t to_float(const dvec& v) {
  bound_vector_t r = traccc::matrix::zero<bound_vector_t>();
  for (std::size_t i = 0; i < N; ++i) {
    traccc::getter::element(r, i, 0) = static_cast<scalar_t>(v[i]);
  }
  return r;
}

detray::geometry::identifier test_surface() {
  detray::geometry::identifier sf{};
  sf.set_volume(1u)
      .set_id(detray::surface_id::e_sensitive)
      .set_index(2u)
      .set_transform(2u);
  return sf;
}

/// Bound track parameters from double precision values
bound_params_t make_params(const detray::geometry::identifier sf, const dvec& v,
                           const dmat& c) {
  bound_params_t p{};
  p.set_surface_link(sf);
  p.set_vector(to_float(v));
  p.set_covariance(to_float(c));
  return p;
}

/// Result of a smoother call
struct smoother_result {
  traccc::kalman_fitter_status status;
  bound_vector_t vec;
  bound_matrix_t cov;
};

/// Run the two-filters smoother on a single 2D measurement
smoother_result run_smoother(const dvec& fw_vec, const dmat& fw_cov,
                             const dvec& pred_vec, const dmat& pred_cov) {
  vecmem::host_memory_resource mr;
  const detray::geometry::identifier sf = test_surface();

  // A 2D measurement close to the forward filtered position
  traccc::edm::measurement_collection::host measurements{mr};
  measurements.push_back({{static_cast<float>(fw_vec[0] + 0.01),
                           static_cast<float>(fw_vec[1] - 0.02)},
                          {0.01f, 0.04f},
                          2u,
                          0.f,
                          0.f,
                          0u,
                          sf,
                          {0u, 1u},
                          0u});
  const auto meas_view = vecmem::get_data(measurements);
  const traccc::edm::measurement_collection::const_device meas_device{
      meas_view};

  traccc::edm::track_state_collection<algebra_t>::host states{mr};
  states.push_back(traccc::edm::make_track_state<algebra_t>(meas_device, 0u));
  auto state = states.at(0u);
  state.filtered_params() = make_params(sf, fw_vec, fw_cov);

  bound_params_t bound_params = make_params(sf, pred_vec, pred_cov);

  const traccc::kalman_fitter_status status =
      traccc::two_filters_smoother<algebra_t>{}(
          state, measurements.at(0u), bound_params,
          traccc::measurement_selector::config{}, false);

  return {status, state.smoothed_params().vector(),
          state.smoothed_params().covariance()};
}

/// Reference: Eq (3.38) in information form, in double precision
std::pair<dvec, dmat> reference(const dvec& fw_vec, const dmat& fw_cov,
                                const dvec& pred_vec, const dmat& pred_cov) {
  const dmat fw_inv = invert(fw_cov);
  const dmat pred_inv = invert(pred_cov);
  const dmat smt_cov = invert(add(fw_inv, pred_inv));
  const dvec a = mult(fw_inv, fw_vec);
  const dvec b = mult(pred_inv, pred_vec);
  dvec sum{};
  for (std::size_t i = 0; i < N; ++i) {
    sum[i] = a[i] + b[i];
  }
  return {mult(smt_cov, sum), smt_cov};
}

/// Compare a smoother result with the reference, relative to the standard
/// deviations of the reference
void expect_near_reference(const smoother_result& res, const dvec& ref_vec,
                           const dmat& ref_cov, const double tol) {
  ASSERT_EQ(res.status, traccc::kalman_fitter_status::SUCCESS);
  for (std::size_t i = 0; i < N; ++i) {
    const double sig_i = std::sqrt(ref_cov[i][i]);
    double diff = traccc::getter::element(res.vec, i, 0) - ref_vec[i];
    if (i == traccc::e_bound_phi) {
      diff = std::remainder(diff, 2. * std::numbers::pi);
    }
    EXPECT_LT(std::abs(diff) / sig_i, tol) << "vector element " << i;
    for (std::size_t j = 0; j < N; ++j) {
      const double sig_j = std::sqrt(ref_cov[j][j]);
      const double cdiff =
          traccc::getter::element(res.cov, i, j) - ref_cov[i][j];
      EXPECT_LT(std::abs(cdiff) / (sig_i * sig_j), tol)
          << "covariance element " << i << ", " << j;
    }
  }
}

}  // namespace

/// The solver reproduces the right-hand side for SPD matrices
TEST(cholesky_solve, spd_system) {
  std::mt19937 gen(42u);
  std::uniform_real_distribution<double> dist(-1., 1.);

  for (unsigned int n = 0u; n < 20u; ++n) {
    const dmat a = random_covariance(gen, 2.);
    dmat x{};
    for (auto& row : x) {
      for (auto& e : row) {
        e = dist(gen);
      }
    }

    // B = A X in double precision
    dmat b{};
    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = 0; j < N; ++j) {
        for (std::size_t k = 0; k < N; ++k) {
          b[i][j] += a[i][k] * x[k][j];
        }
      }
    }

    bound_matrix_t sol = to_float(b);
    ASSERT_TRUE(
        (traccc::details::cholesky_solve<algebra_t, N, N>(to_float(a), sol)));

    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = 0; j < N; ++j) {
        EXPECT_NEAR(traccc::getter::element(sol, i, j), x[i][j],
                    1e-3 / sigmas[i])
            << "element " << i << ", " << j;
      }
    }
  }
}

/// The solver fails for matrices that are not positive definite
TEST(cholesky_solve, not_positive_definite) {
  std::mt19937 gen(7u);
  dmat a = random_covariance(gen);
  // Negative variance
  a[3][3] = -a[3][3];

  bound_matrix_t x = traccc::matrix::identity<bound_matrix_t>();
  EXPECT_FALSE(
      (traccc::details::cholesky_solve<algebra_t, N, N>(to_float(a), x)));

  // Singular matrix
  bound_matrix_t zero = traccc::matrix::zero<bound_matrix_t>();
  EXPECT_FALSE((traccc::details::cholesky_solve<algebra_t, N, N>(zero, x)));
}

/// The smoothed parameters match the information form of Eq (3.38)
TEST(two_filters_smoother, matches_information_form) {
  std::mt19937 gen(1234u);
  std::normal_distribution<double> normal(0., 1.);

  for (unsigned int n = 0u; n < 50u; ++n) {
    const dmat fw_cov = random_covariance(gen);
    const dmat pred_cov = random_covariance(gen);

    const dvec fw_vec{1.5, -20., 0.3, 1.2, -0.5, 2.};
    dvec pred_vec{};
    for (std::size_t i = 0; i < N; ++i) {
      pred_vec[i] = fw_vec[i] + sigmas[i] * normal(gen);
    }

    const auto [ref_vec, ref_cov] =
        reference(fw_vec, fw_cov, pred_vec, pred_cov);
    const smoother_result res =
        run_smoother(fw_vec, fw_cov, pred_vec, pred_cov);

    expect_near_reference(res, ref_vec, ref_cov, 1e-3);
  }
}

/// The phi difference is wrapped: tracks close to phi = +-pi are smoothed
/// correctly
TEST(two_filters_smoother, phi_wrap_around) {
  std::mt19937 gen(99u);
  const dmat fw_cov = random_covariance(gen);
  const dmat pred_cov = random_covariance(gen);

  const double pi = std::numbers::pi;
  const dvec fw_vec{1.5, -20., pi - 0.001, 1.2, -0.5, 2.};
  // Same direction as above, but expressed across the phi boundary
  const dvec pred_vec{1.5, -20., -pi + 0.001, 1.2, -0.5, 2.};

  // Reference with the unwrapped predicted angle (pi + 0.001)
  dvec pred_unwrapped = pred_vec;
  pred_unwrapped[traccc::e_bound_phi] = pi + 0.001;
  const auto [ref_vec, ref_cov] =
      reference(fw_vec, fw_cov, pred_unwrapped, pred_cov);

  const smoother_result res = run_smoother(fw_vec, fw_cov, pred_vec, pred_cov);

  expect_near_reference(res, ref_vec, ref_cov, 1e-3);

  // The smoothed phi stays at the boundary and in the valid range
  const scalar_t phi = traccc::getter::element(res.vec, traccc::e_bound_phi, 0);
  EXPECT_GT(std::abs(phi), pi - 0.01);
  EXPECT_LE(std::abs(phi), pi);
}
