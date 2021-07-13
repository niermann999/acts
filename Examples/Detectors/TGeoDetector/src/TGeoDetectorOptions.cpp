// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TGeoDetector/TGeoDetectorOptions.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <cstdlib>
#include <fstream>
#include <vector>

#include <boost/program_options.hpp>

void ActsExamples::Options::addTGeoGeometryOptions(Description& desc) {
  using boost::program_options::value;

  // due to the way program options handles options that can occur multiple
  // times, all options of a logical block must always be present.
  //
  // each detector volume configuration is one logical block which can
  // be repeated as many times as there are usable detector volumes.
  //
  // required per-volume options:
  //
  //   --geo-tgeo-volume InnerPixels
  //   --geo-tgeo-sfbin-r-tolerance 5:5
  //   --geo-tgeo-sfbin-phi-tolerance 0.025:0.025
  //   --geo-tgeo-sfbin-z-tolerance 5:5
  //   --geo-tgeo-nlayers 0  # boolean switch whether there are negative layers
  //   --geo-tgeo-clayers 1  # boolean switch whether there are central layers
  //   --geo-tgeo-players 0  # boolean switch whether there are positive layers
  //
  // within each volume there can be negative/central/positive layers depending
  // on the which `--geo-tgeo-{n,c,p}layers` flags are set to true. if any of
  // them are set, they must be followed by the corresponding layer option. if
  // the `*layers` option is false, **no** further options **must** be set.
  //
  // examples: negative and central layers, but not positive layers
  //
  //   --geo-tgeo-nlayers 1
  //   --geo-tgeo-nvolume-name Pixel::Pixel
  //   --geo-tgeo-nmodule-name Pixel::siLog
  //   --geo-tgeo-nmodule-axes YZX
  //   --geo-tgeo-nlayer-r-range 0:135
  //   --geo-tgeo-nlayer-z-range -3000:-250
  //   --geo-tgeo-nlayer-r-split -1.
  //   --geo-tgeo-nlayer-z-split 10.
  //
  //   --geo-tgeo-clayers 1
  //   --geo-tgeo-cvolume-name Pixel::Pixel
  //   --geo-tgeo-cmodule-name Pixel::siLog
  //   --geo-tgeo-cmodule-axes YZX
  //   --geo-tgeo-clayer-r-range 0:135
  //   --geo-tgeo-clayer-z-range -250:250
  //   --geo-tgeo-clayer-r-split 5.
  //   --geo-tgeo-clayer-z-split -1.
  //
  //   --geo-tgeo-players 0
  //   # no --geo-tgeo-{cvolume,cmodule,clayer}* options
  //
  auto opt = desc.add_options();
  // required global options
  opt("geo-tgeo-filename", value<std::string>()->default_value(""),
      "Root file name.");
  opt("geo-tgeo-jsonconfig", value<std::string>()->default_value(""),
      "Json config file name.");
}

std::vector<double> ActsExamples::Options::readBeampipeBuilderParam(const std::string& path) {

  nlohmann::json djson;
  if (path.empty()) {
    return {};
  }
  std::ifstream infile(path, std::ifstream::in | std::ifstream::binary);
  // rely on exception for error handling
  infile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  infile >> djson;

  return djson["geo-tgeo-beampipe-parameters"].get<std::vector<double>>();
}

void ActsExamples::Options::from_json(const nlohmann::json& j,
                             Acts::TGeoLayerBuilder::SplitterConfig& msc) {

  msc.moduleType = j["geo-tgeo-splitter-type"];
  msc.splitParamMap = j["geo-tgeo-splitter-map"].get<std::map<std::string, std::vector<unsigned int>>>();
}

void ActsExamples::Options::from_json(const nlohmann::json& j,
                             Acts::TGeoLayerBuilder::LayerConfig& psc) {

  psc.volumeName  = j["geo-tgeo-volume-name"];
  psc.sensorNames = j["geo-tgeo-module-name"].get<std::vector<std::string>>();
  psc.localAxes   = j["geo-tgeo-module-axes"];
  auto r_range = j["geo-tgeo-layer-r-range"].get<std::pair<double, double>>();
  auto z_range = j["geo-tgeo-layer-z-range"].get<std::pair<double, double>>();
  psc.parseRanges = {{Acts::binR, r_range}, {Acts::binZ, z_range}};
  double r_split = j["geo-tgeo-layer-r-split"];
  double z_split = j["geo-tgeo-layer-z-split"];
  if (0 < r_split) {
      psc.splitConfigs.emplace_back(Acts::binR, r_split);
  }
  if (0 < z_split) {
      psc.splitConfigs.emplace_back(Acts::binZ, z_split);
  }
}


std::vector<Acts::TGeoLayerBuilder::Config>
ActsExamples::Options::readTGeoLayerBuilderConfigs(const std::string& path) {
  std::vector<Acts::TGeoLayerBuilder::Config> detLayerConfigs = {};

  nlohmann::json djson;
  if (path.empty()) {
    return detLayerConfigs;
  }
  std::ifstream infile(path, std::ifstream::in | std::ifstream::binary);
  // rely on exception for error handling
  infile.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  infile >> djson;

  double unitScalor = djson["geo-tgeo-unit-scalor"];
  std::string worldVolume = djson["geo-tgeo-worldvolume"];

  for (const auto& volume : djson["LayerConfigs"]["Volumes"]) {
    Acts::TGeoLayerBuilder::Config layerBuilderConfig;
    // subdetector selection
    std::string subDetector = volume["geo-tgeo-volume"];
    layerBuilderConfig.configurationName = subDetector;
    layerBuilderConfig.unit = unitScalor;

    // configure surface autobinning
    auto binToleranceR =
      volume["geo-tgeo-sfbin-r-tolerance"].get<std::vector<double>>();
    auto binToleranceZ =
      volume["geo-tgeo-sfbin-z-tolerance"].get<std::vector<double>>();
    auto binTolerancePhi =
      volume["geo-tgeo-sfbin-phi-tolerance"].get<std::vector<double>>();

    std::vector<std::pair<double, double>> binTolerances(
        static_cast<size_t>(Acts::binValues), {0., 0.});
    binTolerances[Acts::binR] = {binToleranceR[0],
                                 binToleranceR[1]};
    binTolerances[Acts::binZ] = {binToleranceZ[0],
                                 binToleranceZ[1]};
    binTolerances[Acts::binPhi] = {binTolerancePhi[0],
                                   binTolerancePhi[1]};
    layerBuilderConfig.autoSurfaceBinning = true;
    layerBuilderConfig.surfaceBinMatcher =
        Acts::SurfaceBinningMatcher(binTolerances);

    auto isLayers = std::vector<bool>{volume["geo-tgeo-nlayers"], 
                                      volume["geo-tgeo-clayers"],
                                      volume["geo-tgeo-players"]};
    size_t ncp = 0;
    for (const auto& layer : volume["Layers"]) {
      if (not isLayers[ncp]) {
        ncp++;
        continue;
      }
      Acts::TGeoLayerBuilder::LayerConfig lConfig;
      from_json(layer, lConfig);

      for(const auto& splitter : layer["Splitters"]) {
        Acts::TGeoLayerBuilder::SplitterConfig sConfig;
        from_json(splitter, sConfig);
        lConfig.splitterConfigs.push_back(sConfig);
      }
      layerBuilderConfig.layerConfigurations[ncp++].push_back(lConfig);
    }

    detLayerConfigs.push_back(layerBuilderConfig);
  }

  return detLayerConfigs;
}
