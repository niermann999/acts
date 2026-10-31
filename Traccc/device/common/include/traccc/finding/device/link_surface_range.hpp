// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace traccc::device {

/// Range of surface descriptors that the propagation recorded for a
/// departure link of the CKF (used for the Kalman smoother)
///
/// The surfaces between a link @c P and its next sensitive surface are stored
/// at @c P in the order of the propagation: @c [material..., sensitive].
/// At the first step, the seed surface is the first entry, if the link is not
/// a hole. Sensitive surfaces without a track state are removed after the
/// propagation, so the range can be copied into the surface sequence of the
/// track without further checks.
struct link_surface_range {
  /// Position of the first surface of this range in the surface sequence of
  /// every track that contains the link
  unsigned int offset = 0u;
  /// Number of surface descriptors in the range
  unsigned int n_surfaces = 0u;
};

}  // namespace traccc::device
