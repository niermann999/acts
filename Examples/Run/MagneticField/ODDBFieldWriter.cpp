// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/BFieldMapUtils.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SolenoidBField.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Plugins/BField/BFieldOptions.hpp"

#include <string>

#include <boost/program_options.hpp>

#include "BFieldWritingBase.hpp"

namespace po = boost::program_options;

using namespace Acts::UnitLiterals;

template <class T>
struct always_false : std::false_type {};

/// The main executable
///
/// Creates an InterpolatedBFieldMap from a txt or csv file and writes out the
/// grid points and values of the map into root format. The Field can then be
/// displayed using the root script printBField.cpp

int main(int argc, char* argv[]) {
  using boost::program_options::value;

  // setup and parse options
  auto desc = ActsExamples::Options::makeDefaultOptions();
  ActsExamples::Options::addBFieldOptions(desc);
  desc.add_options()("bf-file-out",
                     value<std::string>()->default_value("BFieldOut.root"),
                     "Set this name for an output root file.")(
      "bf-map-out", value<std::string>()->default_value("bField"),
      "Set this name for the tree in the out file.")(
      "bf-out-rz", value<bool>()->default_value(false),
      "Please set this flag to true, if you want to print out the field map in "
      "cylinder coordinates (r,z). The default are cartesian coordinates "
      "(x,y,z). ")(
      "bf-rRange", value<read_range>()->multitoken(),
      "[optional] range which the bfield map should be written out in either r "
      "(cylinder "
      "coordinates) or x/y (cartesian coordinates)  in [mm]. In case no value "
      "is handed over the whole map will be written out. Please "
      "hand over by simply seperating the values by space")(
      "bf-zRange", value<read_range>()->multitoken(),
      "[optional] range which the bfield map should be written out in z in "
      "[mm].In case no value is handed over for 'bf-rRange' and 'bf-zRange the "
      "whole map will be written out. "
      "Please hand over by simply seperating the values by space")(
      "bf-rBins", value<size_t>()->default_value(200),
      "[optional] The number of bins in r. This parameter only needs to be "
      "specified if 'bf-rRange' and 'bf-zRange' are given.")(
      "bf-ZBins", value<size_t>()->default_value(300),
      "[optional] The number of bins in z. This parameter only needs to be "
      "specified if 'bf-rRange' and 'bf-zRange' are given.")(
      "bf-PhiBins", value<size_t>()->default_value(100),
      "[optional] The number of bins in phi. This parameter only needs to be "
      "specified if 'bf-rRange' and 'bf-zRange' are given and 'bf-out-rz' is "
      "turned on.");
  auto vm = ActsExamples::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  //auto bFieldVar = ActsExamples::Options::readBField(vm);
  const double L = 2.85_m;
  const double R = (1.200 + 1.140) * 0.5_m;
  const size_t nCoils = 1200;
  const double bMagCenter = 2.5_T;
  const size_t nBinsR = 150;
  const size_t nBinsZ = 200;

  double rMin = -0.1;
  double rMax = R * 2.;
  double zMin = 2 * (-L / 2.);
  double zMax = 2 * (L / 2.);
  
  Acts::SolenoidBField bSolenoidField({R, L, nCoils, bMagCenter});
  std::cout << "Building interpolated field map" << std::endl;
  auto mapper = Acts::solenoidFieldMapper({rMin, rMax}, {zMin, zMax},
                                          {nBinsR, nBinsZ}, bSolenoidField);
  using BField_t = Acts::InterpolatedBFieldMap<decltype(mapper)>;

  BField_t::Config cfg(std::move(mapper));
  auto bFieldMap = BField_t(std::move(cfg));
  /*using Cache = typename BField_t::Cache;
  Acts::MagneticFieldContext mctx{};*/

  ActsExamples::BField::writeField<BField_t>(vm, std::make_shared<const BField_t>(bFieldMap));
  
  /*return std::visit(
      [&](auto& bField) -> int {
        using field_type =
            typename std::decay_t<decltype(bField)>::element_type;
        if constexpr (!std::is_same_v<field_type, InterpolatedBFieldMap2D> &&
                      !std::is_same_v<field_type, InterpolatedBFieldMap3D>) {
          std::cout << "Bfield map could not be read. Exiting." << std::endl;
          return EXIT_FAILURE;
        } else {
          ActsExamples::BField::writeField<field_type>(vm, bField);
          return EXIT_SUCCESS;
        }
      },
      bFieldVar);*/
}
