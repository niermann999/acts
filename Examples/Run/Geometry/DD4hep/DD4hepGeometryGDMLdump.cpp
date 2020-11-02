// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/DD4hepDetectorOptions.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepGeometryService.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"

#include <boost/program_options.hpp>

using namespace ActsExamples;

int main(int argc, char* argv[]) {
  // Setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addOutputOptions(desc);
  Options::addDD4hepOptions(desc);
  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  // Setup the DD4hep detector
  auto dd4hepCfg = Options::readDD4hepConfig<po::variables_map>(vm);
  auto geometrySvc = std::make_shared<DD4hep::DD4hepGeometryService>(dd4hepCfg);
  auto dd4hepDet = *geometrySvc->lcdd();

  // Call dd4hep GDML plugin to convert the geomtery and save it to disk
  // see
  try {
    dd4hepDet.apply("DD4hepGeometry2GDML", 0, "outfile");
    return EXIT_SUCCESS;
  }
  catch(const std::exception& e)  {
      dd4hep::except("RunPlugin","++ Exception while executing plugin %s:\n\t\t%s",
                     name ? name : "<unknown>", e.what());
  }
  catch(...)  {
      dd4hep::except("RunPlugin","++ UNKNOWN Exception while executing plugin %s.",name ? name : "<unknown>");
  }
  ::exit(EINVAL);
  return EINVAL;
}
