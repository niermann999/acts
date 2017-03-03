// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialMapping.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Plugins/MaterialPlugins/MaterialMapping.hpp"
#include "ACTS/EventData/NeutralParameters.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Plugins/MaterialPlugins/SurfaceMaterialRecord.hpp"
#include "ACTS/Surfaces/CylinderBounds.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/Helpers.hpp"

Acts::MaterialMapping::MaterialMapping(const Config&           cnf,
                                       std::unique_ptr<Logger> log)
  : m_cnf(cnf)
  , m_logger(std::move(log))
  , m_surfaceMaterialRecords()

{
  // check if extrapolation engine is given
  if (!m_cnf.extrapolationEngine) {
    ACTS_ERROR("[!]No extrapolation engine given!");
  } else
    ACTS_INFO("Extrapolation engine successfully retrieved!");
}

Acts::MaterialMapping::~MaterialMapping()
{
}

void
Acts::MaterialMapping::setLogger(std::unique_ptr<Logger> newLogger)
{
  m_logger = std::move(newLogger);
}

void
Acts::MaterialMapping::mapMaterialTrackRecord(const MaterialTrackRecord& matTrackRec)
{
  // access the parameters
  double   theta  = matTrackRec.theta();
  double   phi    = matTrackRec.phi();
  auto     spos   = matTrackRec.position();
  Vector3D vertex(spos.x,spos.y,spos.z);
  // get the steps
  std::vector<MaterialStep> materialSteps = matTrackRec.materialSteps();
  
  // let's extrapolate through the ACTS detector and record all surfaces
  // that have a material proxy
  if (materialSteps.size()) {
    // get the number of materialsteps
    ACTS_DEBUG("Successfuly retrieved " << materialSteps.size()
                                        << "materialSteps");
    // propagate through the detector and collect the layers hit in the given
    // direction eta phi
    // calculate the direction in cartesian coordinates
    Vector3D direction(cos(phi)*sin(theta), sin(phi)*sin(theta),cos(theta));
    // create the beginning neutral parameters to extrapolate through the
    // geometry
    // std::unique_ptr<ActsSymMatrixD<NGlobalPars>> cov;
    Acts::NeutralCurvilinearParameters startParameters(
        nullptr, vertex, direction);
    // create a neutral extrapolation cell and configure it:
    // - to collect surfaces with a SurfaceMaterialProxy
    // - to step at the detector boundary 
    // - to run in a FATRAS type approach 
    ExtrapolationCell<NeutralParameters> ecc(startParameters);
    ecc.addConfigurationMode(ExtrapolationMode::StopAtBoundary);
    ecc.addConfigurationMode(ExtrapolationMode::FATRAS);
    ecc.addConfigurationMode(ExtrapolationMode::CollectSensitive);
    ecc.addConfigurationMode(ExtrapolationMode::CollectMaterial);
    // call the extrapolation engine
    // screen output
    ACTS_DEBUG("===> forward extrapolation - collecting material layers <<===");
    // call the extrapolation engine
    Acts::ExtrapolationCode eCode = m_cnf.extrapolationEngine->extrapolate(ecc);
    // end parameter, if there
    if (eCode.isSuccess()) {
      // number of layers hit
      size_t nSurfacesHit = ecc.extrapolationSteps.size();
      ACTS_VERBOSE("[+] Extrapolation to layers did succeed and found "
                   << nLayersHit
                   << " layers.");
      // find all the intersected material - remember the last parameters
      std::unique_ptr<const Acts::NeutralParameters> parameters = nullptr;
      // loop over the collected information
      for (auto& es : ecc.extrapolationSteps) {
        if (es.configuration.checkMode(ExtrapolationMode::CollectMaterial)){
          // get the parameters
          auto parameters = es.parameters;          

        }
        // continue if we have parameters
      }  // loop over extrapolationsteps
    }    // extrapolation success
  }      // stepCollection
  return bool(layersAndHits.size());
}

void
Acts::MaterialMapping::associateLayerMaterial(
    const MaterialTrackRecord& matTrackRec,
    std::vector<std::pair<const Acts::Layer*, Acts::Vector3D>>& layersAndHits)
{
  // now go through the material step collection and find best fitting layer
  // layers are ordered, hence you can move the starting point along
  std::vector<std::pair<const Acts::Layer*, Acts::Vector3D>>::iterator
      currentLayer
      = layersAndHits.begin();
  // access the material steps of this track record
  std::vector<MaterialStep> materialSteps = matTrackRec.materialSteps();

  // create object which connects layer with the original material step and its
  // assigned position on the layer
  std::map<const Acts::Layer*,
           std::pair<const Vector3D, std::vector<MaterialStep>>>
      layersPosAndSteps;

  // loop through hits and find the closest layer, the start point moves
  // outwards as we go
  for (auto& step : materialSteps) {
    ACTS_VERBOSE("[L] starting from layer "
                 << std::distance(layersAndHits.begin(), currentLayer)
                 << " from layer collection for this step.");
    // step length and position
    Acts::Vector3D pos(step.position().x, step.position().y, step.position().z);
    // now find the closest layer
    // if the currentlayer is the last layer and the hit is still inside ->
    // assign & check if the layers before have been assigned the right way -
    // reassign in case another layer fits better
    if (currentLayer != std::prev(layersAndHits.end())) {
      // search through the layers - this is the reference distance for
      // projection
      double currentDistance = (pos - currentLayer->second).mag();
      ACTS_VERBOSE(
          "  - current distance is " << currentDistance << " from "
                                     << Acts::toString(pos)
                                     << " and "
                                     << Acts::toString(currentLayer->second));
      // check if other layer is more suitable
      for (std::vector<std::pair<const Acts::Layer*, Acts::Vector3D>>::iterator
               testLayer
           = std::next(currentLayer);
           testLayer != layersAndHits.end();
           ++testLayer) {
        // calculate the distance to the testlayer
        double testDistance = (pos - testLayer->second).mag();
        ACTS_VERBOSE("[L] Testing layer "
                     << std::distance(layersAndHits.begin(), testLayer)
                     << " from layer collection for this step.");
        ACTS_VERBOSE(
            " - test distance is " << testDistance << " from "
                                   << Acts::toString(pos)
                                   << " and "
                                   << Acts::toString(testLayer->second));
        if (testDistance < currentDistance) {
          ACTS_VERBOSE("[L] Skipping over to current layer "
                       << std::distance(layersAndHits.begin(), testLayer)
                       << " because "
                       << testDistance
                       << " < "
                       << currentDistance);
          // the test distance did shrink - update currentlayer
          currentLayer    = testLayer;
          currentDistance = testDistance;
        }  // testdistance < currentDistance
        else
          break;  // stick to the layer you have
      }           // check for better fitting layers
    }             // if last layer
    // the current layer *should* be correct now
    const Acts::Layer* assignedLayer = currentLayer->first;
    // correct material thickness with pathcorrection
    double theta = matTrackRec.theta();
    double phi   = matTrackRec.phi();
    // calculate the direction in cartesian coordinates
    Acts::Vector3D direction(
        cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
    // access the path correction of the associated material surface
    double pathCorrection
        = assignedLayer->materialSurface()->pathCorrection(pos, direction);

    // create material Properties
    const Acts::MaterialProperties* layerMaterialProperties
        = new MaterialProperties(step.material().material(),
                                 step.material().thickness() / pathCorrection);
    // correct also the thickness of the material step
    Acts::MaterialStep updatedStep(*layerMaterialProperties, step.position());
    // fill the current material step and its assigned position
    // first check if layer is already there
    auto layerPosSteps = layersPosAndSteps.find(assignedLayer);
    // just fill in material step if layer is already there
    // otherwise create new entry
    if (layerPosSteps != layersPosAndSteps.end())
      layerPosSteps->second.second.push_back(updatedStep);
    else
      layersPosAndSteps.emplace(
          assignedLayer,
          std::make_pair(currentLayer->second,
                         std::vector<MaterialStep>{updatedStep}));

    // associate the hit
    ACTS_VERBOSE("[L] Now associate hit " << Acts::toString(pos) << " at "
                                          << currentLayer->second);
  }

  // associate the steps
  for (auto& layer : layersPosAndSteps) {
    associateHit(layer.first, layer.second.first, layer.second.second);
  }
}

void
Acts::MaterialMapping::associateHit(
    const Layer*                           layer,
    const Acts::Vector3D&                  position,
    const std::vector<Acts::MaterialStep>& layerMaterialSteps)
{
  auto layerRecord = m_surfaceMaterialRecords.find(layer);
  // if layer was not present already create new Material Record
  if (layerRecord == m_surfaceMaterialRecords.end()) {
    // get the bin utility
    const Acts::BinUtility* binUtility = layer->material()->binUtility();
    // create the material record
    ACTS_VERBOSE("[L] Creating new Layer record, for layer  "
                 << layer->geoID()
                 << " at position "
                 << Acts::toString(position));
    m_surfaceMaterialRecords[layer] = SurfaceMaterialRecord(binUtility);
  }
  ACTS_VERBOSE("[L] Add new layer material properties  at position "
               << Acts::toString(position));
  // add material to record, if record exists already
  m_surfaceMaterialRecords[layer].addLayerMaterialProperties(position,
                                                   layerMaterialSteps);
}

void
Acts::MaterialMapping::averageLayerMaterial()
{
  ACTS_VERBOSE(m_surfaceMaterialRecords.size() << " SurfaceMaterialRecords to be averaged");
  // average over the layer material
  for (auto& layRecord : m_surfaceMaterialRecords) {
    layRecord.second.averageMaterial();
  }
}

void
Acts::MaterialMapping::finalizeLayerMaterial()
{
  ACTS_VERBOSE(m_surfaceMaterialRecords.size()
               << " SurfaceMaterialRecords to be finalized");
  // finally set the material of the layers
  for (auto& layRecord : m_surfaceMaterialRecords) {
    // @todo check with Julia how to fix this
    //layRecord.first->materialSurface()->setAssociatedMaterial(
    //    layRecord.second.layerMaterial());
  }
}
