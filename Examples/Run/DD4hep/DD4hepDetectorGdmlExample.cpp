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
#include <fstream>
#include <memory>

#include <DD4hep/Detector.h>

using namespace ActsExamples;

int main(int argc, char* argv[]) {
  // Setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addOutputOptions(desc);
  Options::addDD4hepOptions(desc);
  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  // Setup the DD4hep detector
  dd4hep::Detector& dd4hepDet = dd4hep::Detector::getInstance();

  // DD4hep detector xml inputfile
  //std::string inputFile = vm["dd4hep-input"].as<read_strings>();
  std::string inputFile = "/afs/cern.ch/user/j/jonierma/Project_Files/src/ACTS/Examples/Detectors/DD4hepDetector/compact/OpenDataDetector/OpenDataDetector.xml";
  const char* inputFiles[] = {inputFile.c_str()};

  // Output directory
  //std::string outputDir = vm["output-dir"].template as<std::string>();
  //std::string outputFile = "/eos/user/j/jonierma/data/geometry/OpenDataDetector/DD4hep_detector.gdml";
  std::string outputFile = "dd4hep_detector.gdml";
  const char* outputFiles[] = {0,outputFile.c_str()};
  
  // Logfile
  std::ofstream logfile("log.txt");
  std::unique_ptr<const Acts::Logger> myLogger
      = Acts::getDefaultLogger("MyLogger", Acts::Logging::INFO, &logfile);
  //ACTS_LOCAL_LOGGER(myLogger);
  // Call dd4hep GDML plugin to convert the geomtery and save it to disk
  // see dd4hep/UtilityApps/src/run_plugin.h
  try {
    dd4hepDet.apply("DD4hep_CompactLoader",1,(char**)inputFiles);
    dd4hepDet.apply("DD4hepGeometry2GDML", 1, (char**)&outputFiles[1]);
    //ACTS_INFO("Completed conversion");
    logfile.close();
    return EXIT_SUCCESS;
  }
  catch(const std::exception& e)  {
    //dd4hep::except("DD4hep: RunPlugin","++ Exception while executing plugin "
    //               "%s:\n\t\t%s", name ? name : "<unknown>", e.what());
    //ACTS_INFO("DD4hep plugin error: " + std::string(e.what()));
    std::cout << "DD4hep exception: " << std::string(e.what()) << std::endl;
  }
  catch(...)  {
    //dd4hep::except("DD4hep: RunPlugin","++ UNKNOWN Exception while executing "
    //               "plugin %s.",name ? name : "<unknown>");
    //ACTS_INFO("DD4hep plugin error: <unknown>");
    std::cout << "DD4hep exception: <unknown>" << std::endl;
  }
  logfile.close();
  return EXIT_FAILURE;
}
