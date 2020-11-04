// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/DD4hepDetectorOptions.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepGeometryService.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <boost/program_options.hpp>

#include <DD4hep/Detector.h>

using namespace ActsExamples;

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
  // Setup and parse options
  auto desc = Options::makeDefaultOptions();
  desc.add_options()
      ("dd4hep-input",
      po::value<std::vector<std::string>>()->multitoken()->default_value(
          {"file:Detectors/DD4hepDetector/compact/OpenDataDetector/"
           "OpenDataDetector.xml"}),
      "The locations of the input DD4hep files, use 'file:foo.xml'. In case "
      "you want to read in multiple files, just seperate the strings by "
      "space.")
    ("output-dir", po::value<std::string>()->default_value(""),
                   "Output directory location.");
  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  // Setup the DD4hep detector
  dd4hep::Detector& dd4hepDet = dd4hep::Detector::getInstance();

  // DD4hep detector xml inputfiles
  std::vector<std::string> geoFiles = vm["dd4hep-input"].as<std::vector<std::string>>();
  const char* inputFiles[] = {geoFiles[0].c_str()};

  // Output directory
  std::string outputDir = vm["output-dir"].template as<std::string>();
  std::string outputFile = joinPaths(outputDir,"dd4hep_detector.gdml");
  // Due to dd4hep argument parsing, an output file name cannot be at 0th position
  const char* outputFiles[] = {0,outputFile.c_str()};
  
  // Logging
  auto logLevel = ActsExamples::Options::readLogLevel(vm);
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("DD4hep2Gdml", logLevel));

  // Call dd4hep GDML plugin to convert the geomtery and save it to disk
  // see dd4hep/UtilityApps/src/run_plugin.h
  try {
    dd4hepDet.apply("DD4hep_CompactLoader",1,(char**)inputFiles);
    dd4hepDet.apply("DD4hepGeometry2GDML", 1, (char**)&outputFiles[1]);
    ACTS_INFO("Completed conversion");
    return EXIT_SUCCESS;
  }
  catch(const std::exception& e)  {
    ACTS_FATAL("DD4hep plugin error: " + std::string(e.what()));
  }
  catch(...)  {
    ACTS_FATAL("DD4hep plugin error: <unknown>");
  }
  return EXIT_FAILURE;
}
