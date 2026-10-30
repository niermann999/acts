// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/finding/device/build_tracks.hpp"
#include "traccc/finding/finding_config.hpp"

// Detray include(s).
#include <detray/core/concepts.hpp>

// VecMem include(s).
#include <vecmem/containers/data/jagged_vector_view.hpp>

// CUDA include(s).
#include <cuda_runtime.h>

namespace traccc::cuda {

template <detray::concepts::detector detector_t>
void build_tracks(
    const dim3& grid_size, const dim3& block_size, std::size_t shared_mem_size,
    const cudaStream_t& stream, bool run_mbf,
    const measurement_selector::config& calib_cfg,
    const device::build_tracks_payload& payload,
    vecmem::data::jagged_vector_view<typename detector_t::surface_type>
        track_surfaces_view);

}  // namespace traccc::cuda
