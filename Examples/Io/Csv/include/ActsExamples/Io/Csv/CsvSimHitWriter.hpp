// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <limits>
#include <string>

namespace ActsExamples {

/// Write out a simhit collection before detector digitization in comma-
/// separated-value format.
///
/// This writes one file per event containing information about the
/// global space points, momenta (before and after interaction )and hit index
/// into the configured output directory. By default it writes to the
/// current working directory. Files are named using the following schema
///
///     event000000001-<stem>.csv
///     event000000002-<stem>.csv
///     ...
///
/// and each line in the file corresponds to one simhit.
class CsvSimHitWriter final : public WriterT<SimHitContainer>> {
 public:
  struct Config {
    /// Which simulated (truth) hits collection to use.
    std::string inputSimulatedHits;
    /// Where to place output files
    std::string outputDir;
    /// Output filename stem.
    std::string outputStem;
    /// Number of decimal digits for floating point precision in output.
    size_t outputPrecision = std::numeric_limits<float>::max_digits10;
  };

  /// Construct the cluster writer.
  ///
  /// @params cfg is the configuration object
  /// @params lvl is the logging level
  CsvSimHitWriter(const Config& cfg, Acts::Logging::Level lvl);

 protected:
  /// Type-specific write implementation.
  ///
  /// @param[in] ctx is the algorithm context
  /// @param[in] simHits are the simhits to be written
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const SimHitContainer& simHits) final override;

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
