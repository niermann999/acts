// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include <type_traits>

#include "traccc/definitions/math.hpp"
#include "traccc/edm/measurement_helpers.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/edm/track_state_helpers.hpp"
#include "traccc/finding/measurement_selector.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/matrix_helpers.hpp"
#include "traccc/utils/prob.hpp"
#include "traccc/utils/subspace.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/jagged_device_vector.hpp>

namespace traccc::device {

namespace details {

/// Assemble the surface sequence of a track from the surface descriptors
/// that were recorded per departure link during the propagation
///
/// The surfaces between the link @c L and its previous link are stored at
/// the previous link: @c [material..., sensitive surface of L]. For the
/// links of the first step, the seed surface is the first entry. The ranges
/// are collected from the last measurement back to the seed and the sequence
/// is reversed at the end.
///
/// The backward fit with the direct navigator expects a track state on every
/// sensitive surface in the sequence. Therefore, the sensitive surfaces of
/// holes are not added. The same applies to sensitive surfaces that the CKF
/// aborter passed without a measurement search (closer than the minimum step
/// length to the departure surface).
///
/// @tparam surface_t The surface descriptor type of the detector
///
/// @param tip       The tip link of the track
/// @param links     The candidate links
/// @param n_meas    The number of measurements
/// @param payload   The build tracks payload with the link surface data
/// @param sequence  The output surface sequence of the track
///
/// @returns false if the sequence overflowed or is empty
template <typename surface_t>
TRACCC_HOST_DEVICE inline bool fill_surface_sequence(
    const unsigned int tip,
    const vecmem::device_vector<const candidate_link>& links,
    const unsigned int n_meas, const build_tracks_payload& payload,
    vecmem::device_vector<surface_t> sequence) {
  const auto* link_surfaces =
      static_cast<const surface_t*>(payload.link_surfaces_ptr);
  const vecmem::device_vector<const unsigned int> link_n_surfaces(
      payload.link_n_surfaces_view);
  const unsigned int max_n_sf = payload.max_surfaces_per_step;

  sequence.clear();

  // Skip the holes after the last measurement
  candidate_link L = links.at(tip);
  while (L.meas_idx >= n_meas && L.step != 0u) {
    L = links.at(L.previous_candidate_idx);
  }

  // Collect the surfaces of every step in reverse order
  while (L.step != 0u) {
    const unsigned int prev_idx = L.previous_candidate_idx;
    const candidate_link P = links.at(prev_idx);
    const unsigned int n_sf = math::min(link_n_surfaces.at(prev_idx), max_n_sf);

    // Leave out the sensitive surface of L, if it is a hole, and the seed
    // surface, if the first link is a hole
    const unsigned int first =
        (P.step == 0u && P.meas_idx >= n_meas) ? math::min(1u, n_sf) : 0u;
    const unsigned int last =
        (L.meas_idx >= n_meas && n_sf > first) ? n_sf - 1u : n_sf;

    for (unsigned int i = last; i > first; --i) {
      const unsigned int sf_idx = i - 1u;
      const surface_t& sf_desc = link_surfaces[prev_idx * max_n_sf + sf_idx];

      // Only the sensitive surface of L and the seed surface have a track
      // state: Leave out all other sensitive surfaces
      const bool has_track_state =
          (sf_idx == n_sf - 1u) || (P.step == 0u && sf_idx == 0u);
      if (sf_desc.is_sensitive() && !has_track_state) {
        continue;
      }

      if (sequence.size() == sequence.capacity()) {
        TRACCC_WARNING_DEVICE("Surface sequence overflow for track with tip %d",
                              tip);
        sequence.clear();
        return false;
      }
      sequence.push_back(sf_desc);
    }

    L = P;
  }

  // Forward order for the direct navigator
  const unsigned int n_total = sequence.size();
  for (unsigned int i = 0u; i < n_total / 2u; ++i) {
    const surface_t tmp = sequence[i];
    sequence[i] = sequence[n_total - 1u - i];
    sequence[n_total - 1u - i] = tmp;
  }

  return n_total > 0u;
}

}  // namespace details

template <typename surface_t>
TRACCC_HOST_DEVICE inline void build_tracks(
    const global_index_t globalIndex, bool run_mbf,
    const measurement_selector::config calib_cfg,
    const build_tracks_payload& payload,
    vecmem::data::jagged_vector_view<surface_t> track_surfaces_view) {
  const edm::measurement_collection::const_device measurements(
      payload.tracks_view.measurements);

  const bound_track_parameters_collection_types::const_device seeds(
      payload.seeds_view);

  const vecmem::device_vector<const candidate_link> links(payload.links_view);

  const vecmem::device_vector<const unsigned int> tips(payload.tips_view);

  edm::track_collection<default_algebra>::device track_candidates(
      payload.tracks_view.tracks);
  edm::track_state_collection<default_algebra>::device track_states(
      payload.tracks_view.states);
  bound_track_parameters_collection_types::device link_predicted_params(
      payload.link_predicted_parameter_view);
  bound_track_parameters_collection_types::device link_filtered_params(
      payload.link_filtered_parameter_view);

  if (globalIndex >= tips.size()) {
    return;
  }

  const vecmem::device_vector<const unsigned int> tip_to_output_map(
      payload.tip_to_output_map);
  const unsigned int output_idx = tip_to_output_map.size()
                                      ? tip_to_output_map.at(globalIndex)
                                      : globalIndex;

  if (output_idx == std::numeric_limits<unsigned int>::max()) {
    return;
  }

  const unsigned int tip = tips.at(globalIndex);
  edm::track track = track_candidates.at(output_idx);

  // Get the link corresponding to tip
  unsigned int link_idx = tip;
  candidate_link L = links.at(link_idx);

  const unsigned int n_meas = measurements.size();

  // Surface sequence, fit order and liveness for the Kalman smoother
  if (payload.link_surfaces_ptr != nullptr) {
    vecmem::jagged_device_vector<surface_t> track_surfaces(track_surfaces_view);
    const bool is_valid = details::fill_surface_sequence<surface_t>(
        tip, links, n_meas, payload, track_surfaces.at(output_idx));

    vecmem::device_vector<device::sort_key> track_sort_keys(
        payload.track_sort_keys_view);
    vecmem::device_vector<unsigned int> track_indices(
        payload.track_indices_view);
    vecmem::device_vector<unsigned int> track_liveness(
        payload.track_liveness_view);

    // Key = The number of measurements (same as in the Kalman fitter)
    track_sort_keys.at(output_idx) =
        static_cast<device::sort_key>(track.constituent_links().size());
    track_indices.at(output_idx) = output_idx;
    track_liveness.at(output_idx) = is_valid ? 1u : 0u;
  }

  // Both smoothers need the filtered track states
  const bool make_states = (link_filtered_params.capacity() > 0u);

  // Track summary variables
  scalar ndf_sum = 0.f;
  scalar chi2_sum = 0.f;

  bound_matrix<default_algebra> big_lambda_tilde =
      matrix::zero<bound_matrix<default_algebra>>();
  bound_vector<default_algebra> small_lambda_tilde =
      matrix::zero<bound_vector<default_algebra>>();
  bound_matrix<default_algebra> accumulated_jacobian =
      matrix::identity<bound_matrix<default_algebra>>();

  bool mbf_first_processed = false;

  // Reversely iterate to fill the track candidates
  for (auto it = track.constituent_links().rbegin();
       it != track.constituent_links().rend(); it++) {
    while (L.meas_idx >= n_meas && L.step != 0u) {
      if (run_mbf) {
        accumulated_jacobian =
            accumulated_jacobian * payload.jacobian_ptr[link_idx];
      }

      link_idx = L.previous_candidate_idx;
      L = links.at(link_idx);
    }

    assert(L.meas_idx < n_meas);

    if (make_states) {
      const unsigned int track_state_index = track_states.push_back(
          edm::make_track_state<default_algebra>(measurements, L.meas_idx));
      auto track_state = track_states.at(track_state_index);

      track_state.set_hole(false);
      track_state.filtered_chi2() = L.chi2;
      track_state.filtered_params() = link_filtered_params.at(link_idx);

      *it = {edm::track_constituent_link::track_state, track_state_index};

      // MBF smoothing (the Kalman smoother runs after the track finding)
      if (run_mbf) {
        track_state.set_smoothed(true);

        // TODO: The fact that we store the chi2 three times is nonsense.
        track_state.smoothed_chi2() = L.chi2;
        track_state.backward_chi2() = L.chi2;

        const bound_track_parameters<>& predicted_params =
            link_predicted_params.at(link_idx);

        const bound_track_parameters<>::vector_type& predicted_vec =
            predicted_params.vector();
        const bound_track_parameters<>::covariance_type& predicted_covariance =
            predicted_params.covariance();

        bound_vector<default_algebra> small_lambda_hat;
        bound_matrix<default_algebra> big_lambda_hat;

        if (!mbf_first_processed) {
          small_lambda_hat = matrix::zero<bound_vector<default_algebra>>();
          big_lambda_hat = matrix::zero<bound_matrix<default_algebra>>();
          mbf_first_processed = true;
        } else {
          small_lambda_hat =
              matrix::transpose(accumulated_jacobian) * small_lambda_tilde;
          big_lambda_hat = matrix::transpose(accumulated_jacobian) *
                           big_lambda_tilde * accumulated_jacobian;
        }

        const edm::measurement meas = measurements.at(L.meas_idx);

        // Measurement data on surface
        const detray::dmatrix<default_algebra, 2, 1> meas_local =
            measurement_selector::calibrated_measurement_position<
                default_algebra, 2>(meas, calib_cfg);

        // Spatial resolution (Measurement covariance)
        const detray::dmatrix<default_algebra, 2, 2> V =
            measurement_selector::calibrated_measurement_covariance<
                default_algebra, 2>(meas, calib_cfg);

        // TODO: Does not work for line surfaces!
        const detray::dmatrix<default_algebra, 2, e_bound_size> H =
            measurement_selector::observation_model<default_algebra, 2>(
                meas, predicted_params, false);

        const auto S_inv = masked_inverse<default_algebra>(
            H * predicted_covariance * matrix::transpose(H) + V,
            meas.dimensions());

        const auto K = predicted_covariance * matrix::transpose(H) * S_inv;
        const auto C_hat =
            matrix::identity<detray::dmatrix<default_algebra, e_bound_size,
                                             e_bound_size>>() -
            K * H;
        const auto y = meas_local - H * predicted_vec;

        small_lambda_tilde = -1.f * matrix::transpose(H) * S_inv * y +
                             matrix::transpose(C_hat) * small_lambda_hat;
        big_lambda_tilde = matrix::transpose(H) * S_inv * H +
                           matrix::transpose(C_hat) * big_lambda_hat * C_hat;

        track_state.smoothed_params().vector() =
            predicted_vec - predicted_covariance * small_lambda_tilde;
        track_state.smoothed_params().covariance() =
            predicted_covariance -
            (predicted_covariance * big_lambda_tilde * predicted_covariance);

        accumulated_jacobian = payload.jacobian_ptr[link_idx];
      }
    } else {
      *it = {edm::track_constituent_link::measurement, L.meas_idx};
    }

    // Sanity check on chi2
    assert(L.chi2 < std::numeric_limits<traccc::scalar>::max());
    assert(L.chi2 >= 0.f);

    ndf_sum += static_cast<scalar>(measurements.at(L.meas_idx).dimensions());
    chi2_sum += L.chi2;

    // Break the loop if the iterator is at the first candidate and fill the
    // seed and track quality
    if (it == track.constituent_links().rend() - 1) {
      if (run_mbf) {
        track.fit_outcome() = track_fit_outcome::SUCCESS;
        track.params() = track_states.at(it->index).smoothed_params();
      } else {
        track.fit_outcome() = track_fit_outcome::UNKNOWN;
        track.params() = seeds.at(L.seed_idx);
      }
      track.ndf() = ndf_sum - 5.f;
      track.chi2() = chi2_sum;
      track.pval() = prob(track.chi2(), track.ndf());
      track.nholes() = L.n_skipped;
    } else {
      link_idx = L.previous_candidate_idx;
      L = links.at(link_idx);
    }
  }

#ifndef NDEBUG
  // Assert that we did not make any duplicate track states.
  for (const auto& i : track.constituent_links()) {
    for (const auto& j : track.constituent_links()) {
      if (i.index != j.index) {
        // TODO: Re-enable me!
        // assert(measurements.at(i->index).identifier() !=
        //       measurement.at(j->index).identifier());
      }
    }
  }
#endif
}

}  // namespace traccc::device
