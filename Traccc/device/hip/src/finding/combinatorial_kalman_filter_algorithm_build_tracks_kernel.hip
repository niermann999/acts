// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "../utils/global_index.hpp"
#include "../utils/hip_error_handling.hpp"
#include "../utils/utils.hpp"
#include "traccc/hip/finding/combinatorial_kalman_filter_algorithm.hpp"

// Project include(s).
#include "traccc/finding/device/build_tracks.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/detector_buffer.hpp"
#include "traccc/geometry/detector_type_list.hpp"

// HIP include(s).
#include <hip/hip_runtime.h>

namespace traccc::hip {
namespace kernels {

/// HIP kernel for running @c traccc::device::build_tracks
template <typename surface_t>
__global__ void build_tracks(
    const bool run_mbf, const measurement_selector::config calib_cfg,
    const device::build_tracks_payload payload,
    const vecmem::data::jagged_vector_view<surface_t> track_surfaces_view) {
  device::build_tracks<surface_t>(details::global_index1(), run_mbf, calib_cfg,
                                  payload, track_surfaces_view);
}

}  // namespace kernels

void combinatorial_kalman_filter_algorithm::build_tracks_kernel(
    unsigned int n_threads, bool run_mbf_smoother,
    const measurement_selector::config& calib_cfg, const detector_buffer& det,
    const move_only_any& track_surfaces,
    const device::build_tracks_payload& payload) const {
  // Establish the kernel launch parameters.
  const unsigned int deviceThreads = warp_size() * 2;
  const unsigned int deviceBlocks =
      (n_threads + deviceThreads - 1) / deviceThreads;

  // Launch the kernel for the appropriate detector type.
  detector_buffer_visitor<detector_type_list>(
      det, [&]<detray::concepts::detector detector_t>(
               const detray::detector_view_t<detector_t>&) {
        using surface_t =
            typename detray::detector_device_t<detector_t>::surface_type;

        hipLaunchKernelGGL(
            kernels::build_tracks<surface_t>, deviceBlocks, deviceThreads, 0,
            details::get_stream(stream()), run_mbf_smoother, calib_cfg, payload,
            get_track_surfaces_view<surface_t>(track_surfaces, payload));
        TRACCC_HIP_ERROR_CHECK(hipGetLastError());
      });
}

}  // namespace traccc::hip
