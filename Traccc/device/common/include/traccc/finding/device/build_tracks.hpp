// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/device/sort_key.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/finding/candidate_link.hpp"
#include "traccc/finding/device/link_surface_range.hpp"
#include "traccc/finding/finding_config.hpp"
#include "traccc/finding/measurement_selector.hpp"

// VecMem include(s).
#include <vecmem/containers/data/jagged_vector_view.hpp>
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// (Event Data) Payload for the @c traccc::device::build_tracks function
struct build_tracks_payload {
  /**
   * @brief View object to the vector of measurements
   */
  bound_track_parameters_collection_types::const_view seeds_view;

  /**
   * @brief View object to the vector of candidate links
   */
  vecmem::data::vector_view<const candidate_link> links_view;

  /**
   * @brief View object to the vector of tips
   */
  vecmem::data::vector_view<const unsigned int> tips_view;

  /**
   * @brief View object to the vector of track candidates
   */
  edm::track_container<default_algebra>::view tracks_view;

  /**
   * @brief Optional mapping from tip index to output index
   */
  vecmem::data::vector_view<const unsigned int> tip_to_output_map;
  bound_matrix<default_algebra>* jacobian_ptr = nullptr;
  bound_track_parameters_collection_types::view link_predicted_parameter_view;
  bound_track_parameters_collection_types::view link_filtered_parameter_view;

  /**
   * @brief Type-erased pointer to the surface descriptors per departure link
   * (detector surface type, @c max_surfaces_per_step entries per link). The
   * surface sequences per track are only filled if this pointer is set.
   */
  const void* link_surfaces_ptr = nullptr;

  /**
   * @brief View object to the recorded surface range per departure link
   */
  vecmem::data::vector_view<const link_surface_range>
      link_surface_ranges_view{};

  /**
   * @brief Maximum number of surfaces that can be recorded per link
   */
  unsigned int max_surfaces_per_step = 0u;

  /**
   * @brief View object to the sort keys of the output tracks for the Kalman
   * smoother (number of measurements). Only filled with the surface sequences
   */
  vecmem::data::vector_view<device::sort_key> track_sort_keys_view{};

  /**
   * @brief View object to the output track indices for the Kalman smoother
   */
  vecmem::data::vector_view<unsigned int> track_indices_view{};

  /**
   * @brief View object to the liveness of the output tracks for the Kalman
   * smoother (0 if the surface sequence is not valid)
   */
  vecmem::data::vector_view<unsigned int> track_liveness_view{};
};

/// Function for building full tracks from the link container:
/// The full tracks are built using the link container and tip link container.
/// Since every link holds an information of the link from the previous step,
/// we can build a full track by iterating from a tip link backwardly.
///
/// If the surface descriptors were recorded during the propagation, the
/// surface sequence of every track is assembled for the Kalman smoother.
///
/// @tparam surface_t The surface descriptor type of the detector
///
/// @param[in] globalIndex         The index of the current thread
/// @param[in] run_mbf             Whether to run the MBF smoother
/// @param[in] calib_cfg           The measurement calibration configuration
/// @param[in,out] payload         The function call payload
/// @param[out] track_surfaces_view The surface sequences per output track
///
template <typename surface_t>
TRACCC_HOST_DEVICE inline void build_tracks(
    global_index_t globalIndex, bool run_mbf,
    const measurement_selector::config calib_cfg,
    const build_tracks_payload& payload,
    vecmem::data::jagged_vector_view<surface_t> track_surfaces_view);

}  // namespace traccc::device

// Include the implementation.
#include "./impl/build_tracks.ipp"
