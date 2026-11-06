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
#include "../fit_backward.hpp"

// Project include(s).
#include "traccc/fitting/device/fit_backward.hpp"

// System include(s).
#include <cassert>

namespace traccc::cuda {
namespace kernels {

template <typename fitter_t>
__global__ __launch_bounds__(128) void fit_backward(
    const __grid_constant__ fitting_config cfg,
    const __grid_constant__ device::fit_payload payload,
    const device::fit_tpayload<
        typename fitter_t::detector_type::const_view_type,
        typename fitter_t::bfield_type, typename fitter_t::surface_type>*
        tpayload) {
  // The device detector object is built once in global memory
  assert(tpayload->det_ptr != nullptr);
  const auto* __restrict__ const det_ptr =
      static_cast<const typename fitter_t::detector_type*>(tpayload->det_ptr);

  device::fit_backward<fitter_t>(details::global_index1(), cfg, payload,
                                 *tpayload, *det_ptr);
}

}  // namespace kernels

template <typename fitter_t>
void fit_backward(const dim3& grid_size, const dim3& block_size,
                  std::size_t shared_mem_size, const cudaStream_t& stream,
                  const fitting_config& cfg, const device::fit_payload& payload,
                  const device::fit_tpayload<
                      typename fitter_t::detector_type::const_view_type,
                      typename fitter_t::bfield_type,
                      typename fitter_t::surface_type>* tpayload) {
  kernels::fit_backward<fitter_t>
      <<<grid_size, block_size, shared_mem_size, stream>>>(cfg, payload,
                                                           tpayload);
}

}  // namespace traccc::cuda
