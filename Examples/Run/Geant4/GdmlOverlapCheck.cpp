// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/program_options.hpp>

#include "ActsExamples/DD4hepDetector/DD4hepDetectorOptions.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepGeometryService.hpp"
#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"
#include "ActsExamples/Geant4DD4hep/DD4hepDetectorConstruction.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"

#include <G4GDMLParser.hh>
#include "G4UImanager.hh"

using namespace ActsExamples;
using namespace ActsExamples;

int main(int argc, char* argv[]) {
  // Setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addDD4hepOptions(desc);
  desc.add_options()(
      "gdml-file",
      boost::program_options::value<std::string>()->default_value(""),
      "GDML detector file.");

  auto vm = Options::parse(desc, argc, argv);


  if (vm.empty()) {
    return EXIT_FAILURE;
  }
  auto gdmlFile = vm["gdml-file"].as<std::string>();

  std::unique_ptr<G4VUserDetectorConstruction> g4detector;
  if (gdmlFile.compare("") != 0) {
    g4detector = std::make_unique<GdmlDetectorConstruction>(gdmlFile, true);
  }
  else {
    auto dd4hepCfg = Options::readDD4hepConfig<po::variables_map>(vm);
    auto geometrySvc = std::make_shared<DD4hep::DD4hepGeometryService>(dd4hepCfg);
    g4detector =
      std::make_unique<DD4hepDetectorConstruction>(*geometrySvc->lcdd());
  }
  // Setup the GDML detector which automatically runs the overlap check
  G4VPhysicalVolume* world = g4detector->Construct();
  
  //G4GDMLParser parser;
  //parser.Write("dd4hep_detector_alt.gdml", world);
  
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4String command = "/geometry/test/run ";
  UImanager->ApplyCommand(command);
  return EXIT_SUCCESS;
}
