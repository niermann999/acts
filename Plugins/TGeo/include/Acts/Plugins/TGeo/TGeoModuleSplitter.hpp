// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/TGeo/TGeoLayerBuilder.hpp"
#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"

#include <map>

namespace Acts {

class TGeoDetectorElement;

class TGeoModuleSplitter : public TGeoLayerBuilder {

  public:

  /// Constructor
  TGeoModuleSplitter() = default;

  TGeoModuleSplitter(std::shared_ptr<TGeoLayerBuilder> layerBuilder) :
    m_layerBuilder(layerBuilder) {}

  /// If applicable, returns a split detector element 
  inline std::vector<std::shared_ptr<const TGeoDetectorElement>> buildElements(std::shared_ptr<const TGeoDetectorElement> detElement) {
    return m_layerBuilder->buildElements(detElement);
  }

  private:

  /// Instance of wrapped layer builder
  std::shared_ptr<TGeoLayerBuilder> m_layerBuilder;
};
}  // namespace Acts