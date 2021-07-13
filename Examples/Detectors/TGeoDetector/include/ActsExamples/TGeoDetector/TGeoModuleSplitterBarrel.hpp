// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/TGeo/TGeoModuleSplitter.hpp"

namespace Acts {

class AddTGeoModuleSplitterBarrel : public TGeoModuleSplitter {

using splitParamMap = std::map<std::string, std::vector<unsigned int>>;

public:

  /// Constructor from another splitter to add splitting of a different shape
  AddTGeoModuleSplitterBarrel(std::shared_ptr<TGeoLayerBuilder> moduleSplitter, splitParamMap map) 
            : m_moduleSplitter(moduleSplitter),  m_paramMap(map) {}

  /// If applicable, returns a split detector element 
  inline std::vector<std::shared_ptr<const TGeoDetectorElement>> buildElements(std::shared_ptr<const TGeoDetectorElement> detElement) {
    return m_moduleSplitter->buildElements(detElement);
  }

private:
    std::shared_ptr<TGeoLayerBuilder> m_moduleSplitter;
    splitParamMap m_paramMap;
};
} // namespace Acts