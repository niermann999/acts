// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/definitions/hints.hpp"
#include "traccc/definitions/math.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"

// System include(s).
#include <cstddef>

namespace traccc::details {

/// Solve the linear system A X = B in place, with A symmetric positive
/// definite (Cholesky decomposition A = L L^T)
///
/// In contrast to a general inversion with pivoting, all element indices are
/// known at compile time. Therefore, the matrices can stay in registers on
/// GPUs. Only the lower triangle of @c A is read.
///
/// @tparam algebra_t The algebra type
/// @tparam N         The number of rows and columns of A
/// @tparam M         The number of right-hand side columns
///
/// @param[in]     A The symmetric positive definite matrix
/// @param[in,out] X The right-hand side B on input, the solution on output
///
/// @returns false if A is not (numerically) positive definite. In this case,
///          the content of @c X is undefined.
///
/// @note The sizes cannot be deduced through the matrix alias: Call as
///       @c cholesky_solve<algebra_t, N, M>(A, X)
template <typename algebra_t, std::size_t N, std::size_t M>
TRACCC_HOST_DEVICE inline bool cholesky_solve(
    const detray::dmatrix<algebra_t, N, N>& A,
    detray::dmatrix<algebra_t, N, M>& X) {
  using scalar_t = detray::dscalar<algebra_t>;

  // Lower triangle of L and the inverse of its diagonal
  detray::dmatrix<algebra_t, N, N> L{};
  scalar_t inv_diag[N];

  // Cholesky-Banachiewicz algorithm
  TRACCC_PRAGMA_UNROLL
  for (std::size_t j = 0; j < N; ++j) {
    scalar_t d = getter::element(A, j, j);
    TRACCC_PRAGMA_UNROLL
    for (std::size_t k = 0; k < j; ++k) {
      d -= getter::element(L, j, k) * getter::element(L, j, k);
    }
    // Also catches NaN
    if (!(d > 0.f)) {
      return false;
    }
    const scalar_t l_jj = math::sqrt(d);
    getter::element(L, j, j) = l_jj;
    inv_diag[j] = 1.f / l_jj;

    TRACCC_PRAGMA_UNROLL
    for (std::size_t i = j + 1; i < N; ++i) {
      scalar_t s = getter::element(A, i, j);
      TRACCC_PRAGMA_UNROLL
      for (std::size_t k = 0; k < j; ++k) {
        s -= getter::element(L, i, k) * getter::element(L, j, k);
      }
      getter::element(L, i, j) = s * inv_diag[j];
    }
  }

  TRACCC_PRAGMA_UNROLL
  for (std::size_t c = 0; c < M; ++c) {
    // Forward substitution: L Y = B
    TRACCC_PRAGMA_UNROLL
    for (std::size_t i = 0; i < N; ++i) {
      scalar_t s = getter::element(X, i, c);
      TRACCC_PRAGMA_UNROLL
      for (std::size_t k = 0; k < i; ++k) {
        s -= getter::element(L, i, k) * getter::element(X, k, c);
      }
      getter::element(X, i, c) = s * inv_diag[i];
    }
    // Back substitution: L^T X = Y
    TRACCC_PRAGMA_UNROLL
    for (std::size_t ii = 0; ii < N; ++ii) {
      const std::size_t i = N - 1 - ii;
      scalar_t s = getter::element(X, i, c);
      TRACCC_PRAGMA_UNROLL
      for (std::size_t k = i + 1; k < N; ++k) {
        s -= getter::element(L, k, i) * getter::element(X, k, c);
      }
      getter::element(X, i, c) = s * inv_diag[i];
    }
  }

  return true;
}

}  // namespace traccc::details
