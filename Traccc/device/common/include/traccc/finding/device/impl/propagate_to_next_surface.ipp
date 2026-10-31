// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// @TODO: Remove this once alpaka compiles without warning.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic warning "-Wmaybe-uninitialized"
#endif

// Project include(s).
#include "traccc/finding/details/combinatorial_kalman_filter_types.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/particle.hpp"
#include "traccc/utils/propagation.hpp"

// Detray include(s).
#include <detray/utils/tuple_helpers.hpp>

namespace traccc::device {

namespace details {

/// Implementation of the propagation to the next surface for both smoother
/// variants
///
/// @tparam propagator_t    The propagator type
/// @tparam bfield_t        The magnetic field type
/// @tparam kalman_smoother Record the surfaces for the Kalman smoother
///                         (otherwise, accumulate the Jacobians if the MBF
///                         smoother is configured)
template <typename propagator_t, typename bfield_t, bool kalman_smoother>
TRACCC_HOST_DEVICE inline void propagate_to_next_surface_impl(
    const global_index_t globalIndex, const finding_config& cfg,
    const typename propagator_t::detector_type* const det_data_ptr,
    const bfield_t& field_data,
    const propagate_to_next_surface_payload& payload) {
  using algebra_t = typename propagator_t::detector_type::algebra_type;
  using scalar_t = detray::dscalar<algebra_t>;

  if (globalIndex >= payload.n_in_params) {
    return;
  }

  // Theta id
  vecmem::device_vector<const unsigned int> param_ids(payload.param_ids_view);

  const unsigned int param_id = param_ids.at(globalIndex);

  // Links
  vecmem::device_vector<const candidate_link> links(payload.links_view);

  const unsigned int link_idx = payload.prev_links_idx + param_id;
  const candidate_link link = links.at(link_idx);
  assert(link.step == payload.step);
  const unsigned int n_cands = link.step + 1 - link.n_skipped;

  // Parameter liveness
  vecmem::device_vector<unsigned int> params_liveness(
      payload.params_liveness_view);

  // tips
  vecmem::device_vector<unsigned int> tips(payload.tips_view);
  vecmem::device_vector<unsigned int> tip_lengths(payload.tip_lengths_view);

  // Parameters
  bound_track_parameters_collection_types::device params(payload.params_view);

  if (params_liveness.at(param_id) == 0u) {
    return;
  }

  // Input bound track parameter
  const bound_track_parameters<> in_par = params.at(param_id);

  // Create propagator
  propagator_t propagator(cfg.propagation);

  // Create propagator state
  typename propagator_t::state propagation(in_par, field_data, *det_data_ptr);
  propagation.set_particle(
      traccc::detail::correct_particle_hypothesis(cfg.ptc_hypothesis, in_par));
  propagation.stepping()
      .template set_constraint<detray::step::constraint::e_accuracy>(
          cfg.propagation.stepping.step_constraint);

  // Actor states
  // Pathlimit aborter
  typename detray::actor::pathlimit_aborter<scalar_t>::state aborter_state{};
  // Parameter updater
  typename detray::actor::parameter_updater_state<algebra_t> updater_state{
      cfg.propagation, in_par};
  // CKF-interactor
  traccc::details::ckf_interactor_t::state interactor_state;
  // Momentum aborter
  typename detray::actor::momentum_aborter<scalar_t>::state
      momentum_aborter_state{};
  // CKF aborter
  typename ckf_aborter::state ckf_aborter_state{};

  /*
   * If we are running the MBF smoother, we need to accumulate the Jacobians
   * between the two sensitives multiplicatively. To this end, we ask the
   * parameter transporter to multiply the Jacobians into this matrix, which
   * is set to the multiplicative identity.
   */
  if constexpr (!kalman_smoother) {
    if (cfg.run_smoother == smoother_type::e_mbf) {
      assert(payload.tmp_jacobian_view.ptr() != nullptr);

      vecmem::device_vector<bound_matrix<default_algebra>> tmp_jacobian(
          payload.tmp_jacobian_view);
      tmp_jacobian.at(param_id) =
          matrix::identity<bound_matrix<default_algebra>>();
      updater_state.set_full_jacobian(&(tmp_jacobian.at(param_id)));
    }
  }

  /*
   * If we are running the Kalman smoother, the parameter transporter records
   * the descriptors of all sensitive and material surfaces between the
   * departure surface (excluded) and the next sensitive surface (included).
   * They are stored per departure link and shared by all of its branches. At
   * the first step, also the seed surface has to be added, if it has a track
   * state (i.e. the link is not a hole).
   */
  using surface_t = typename propagator_t::detector_type::surface_type;
  // Recomputed where needed, so that they are not live during the propagation
  [[maybe_unused]] auto get_link_surfaces = [&payload, link_idx]() {
    return static_cast<surface_t*>(payload.link_surfaces_ptr) +
           link_idx * payload.max_surfaces_per_step;
  };
  [[maybe_unused]] auto is_seed_surface_recorded = [&payload, &link]() {
    return (link.step == 0u) && (link.meas_idx < payload.n_measurements);
  };

  if constexpr (kalman_smoother) {
    assert(payload.link_surfaces_ptr != nullptr);
    assert(payload.max_surfaces_per_step > 0u);

    updater_state.set_surface_sequence(get_link_surfaces(),
                                       payload.max_surfaces_per_step);

    if (is_seed_surface_recorded()) {
      updater_state.template record_surface<surface_t>(
          det_data_ptr->surface(in_par.surface_link()));
    }
  }

  // Notify the KF and material interaction only at the first propagation
  // initialization. For all subsequent initializations, they will have run
  // already on the previous step
  updater_state.notify_on_initial(link.step == 0);
  momentum_aborter_state.min_pT(static_cast<scalar_t>(cfg.min_pT));
  momentum_aborter_state.min_p(static_cast<scalar_t>(cfg.min_p));
  ckf_aborter_state.min_step_length = cfg.min_step_length_for_next_surface;
  ckf_aborter_state.max_count = cfg.max_step_counts_for_next_surface;

  // Propagate to the next surface
  propagator.propagate(
      propagation, detray::tie(aborter_state, updater_state, interactor_state,
                               momentum_aborter_state, ckf_aborter_state));

  // If a surface found, add the parameter for the next step
  if (ckf_aborter_state.success) {
    assert(propagation.navigation().is_on_sensitive());
    assert(!updater_state.bound_params().is_invalid());

    params[param_id] = updater_state.bound_params();
    params_liveness[param_id] = 1u;

    const scalar theta = params[param_id].theta();
    if (theta <= 0.f || theta >= 2.f * constant<traccc::scalar>::pi) {
      TRACCC_ERROR_DEVICE("Theta is zero after propagation");
      params_liveness[param_id] = 0u;
    }

    if (!std::isfinite(params[param_id].phi())) {
      TRACCC_ERROR_DEVICE(
          "Phi is infinite after propagation (Matrix inversion)");
      params_liveness[param_id] = 0u;
    }

    if (math::fabs(params[param_id].qop()) == 0.f) {
      TRACCC_ERROR_DEVICE("q/p is zero after propagation");
      params_liveness[param_id] = 0u;
    }
  } else {
    params_liveness[param_id] = 0u;
  }

  if constexpr (kalman_smoother) {
    // Do not continue tracks with an incomplete surface sequence
    if (updater_state.surface_sequence_overflow()) {
      TRACCC_WARNING_DEVICE(
          "Surface sequence overflow: %d surfaces found (capacity %d)",
          updater_state.n_recorded_surfaces(), payload.max_surfaces_per_step);
      params_liveness[param_id] = 0u;
    }

    link_surface_range range{};

    // The range is only used if the track continues: Remove the sensitive
    // surfaces without a track state once per link (i.e. the ones that the
    // CKF aborter passed without a measurement search). Only the seed surface
    // and the destination surface (last entry) keep their sensitive surface.
    if (params_liveness[param_id] != 0u) {
      surface_t* const link_surfaces = get_link_surfaces();
      const unsigned int n_recorded = updater_state.n_recorded_surfaces();
      const unsigned int first = is_seed_surface_recorded() ? 1u : 0u;
      assert(n_recorded > first);

      unsigned int n_kept = first;
      for (unsigned int i = first; i < n_recorded - 1u; ++i) {
        if (!link_surfaces[i].is_sensitive()) {
          link_surfaces[n_kept++] = link_surfaces[i];
        }
      }
      link_surfaces[n_kept++] = link_surfaces[n_recorded - 1u];
      range.n_surfaces = n_kept;
    }

    // The surfaces of this link follow the surfaces of the previous link.
    // If this link is a hole, its sensitive surface (last entry of the
    // previous range) is left out
    vecmem::device_vector<link_surface_range> link_surface_ranges(
        payload.link_surface_ranges_view);
    if (link.step != 0u) {
      const link_surface_range prev_range =
          link_surface_ranges.at(link.previous_candidate_idx);
      const unsigned int is_hole =
          (link.meas_idx >= payload.n_measurements) ? 1u : 0u;
      assert(prev_range.n_surfaces >= is_hole);
      range.offset = prev_range.offset + prev_range.n_surfaces - is_hole;
    }
    link_surface_ranges.at(link_idx) = range;
  }

  if (params_liveness[param_id] == 0 &&
      n_cands >= cfg.min_track_candidates_per_track) {
    TRACCC_VERBOSE_DEVICE("Create tip: No next sensitive found");
    auto tip_pos = tips.push_back(link_idx);
    tip_lengths.at(tip_pos) = n_cands;
  }
}

}  // namespace details

template <typename propagator_t, typename bfield_t>
TRACCC_HOST_DEVICE inline void propagate_to_next_surface(
    const global_index_t globalIndex, const finding_config& cfg,
    const typename propagator_t::detector_type* const det_data_ptr,
    const bfield_t& field_data,
    const propagate_to_next_surface_payload& payload) {
  details::propagate_to_next_surface_impl<propagator_t, bfield_t, false>(
      globalIndex, cfg, det_data_ptr, field_data, payload);
}

template <typename propagator_t, typename bfield_t>
TRACCC_HOST_DEVICE inline void propagate_to_next_surface_kalman(
    const global_index_t globalIndex, const finding_config& cfg,
    const typename propagator_t::detector_type* const det_data_ptr,
    const bfield_t& field_data,
    const propagate_to_next_surface_payload& payload) {
  details::propagate_to_next_surface_impl<propagator_t, bfield_t, true>(
      globalIndex, cfg, det_data_ptr, field_data, payload);
}

}  // namespace traccc::device
