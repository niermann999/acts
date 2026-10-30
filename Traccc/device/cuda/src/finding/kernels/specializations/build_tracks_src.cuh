// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "../../../utils/global_index.hpp"
#include "../build_tracks.cuh"

// Project include(s).
#include "traccc/finding/device/build_tracks.hpp"

namespace traccc::cuda {
namespace kernels {

template <typename surface_t>
__global__ void build_tracks(
    const __grid_constant__ bool run_mbf,
    const __grid_constant__ measurement_selector::config calib_cfg,
    const __grid_constant__ device::build_tracks_payload payload,
    const __grid_constant__ vecmem::data::jagged_vector_view<surface_t>
        track_surfaces_view) {
  device::build_tracks<surface_t>(details::global_index1(), run_mbf, calib_cfg,
                                  payload, track_surfaces_view);
}

}  // namespace kernels

template <detray::concepts::detector detector_t>
void build_tracks(
    const dim3& grid_size, const dim3& block_size, std::size_t shared_mem_size,
    const cudaStream_t& stream, bool run_mbf,
    const measurement_selector::config& calib_cfg,
    const device::build_tracks_payload& payload,
    vecmem::data::jagged_vector_view<typename detector_t::surface_type>
        track_surfaces_view) {
  kernels::build_tracks<typename detector_t::surface_type>
      <<<grid_size, block_size, shared_mem_size, stream>>>(
          run_mbf, calib_cfg, payload, track_surfaces_view);
}

}  // namespace traccc::cuda
