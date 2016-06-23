// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ConeLayer.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_LAYERS_CONELAYER_H
#define ACTS_LAYERS_CONELAYER_H

#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Surfaces/ConeSurface.hpp"
#include <algorithm>

namespace Acts {

class ConeBounds;
class OverlapDescriptor;
class ApproachDescriptor;

/// 
/// @class ConeLayer
/// 
/// Class to describe a conical detector layer for tracking, it inhertis from
/// both, Layer base class and ConeSurface class
/// 
class ConeLayer : virtual public ConeSurface, public Layer
{
public:
  /// Factory for shared layer
  /// 
  /// @param transform is the 3D transform that poisitions the layer in 3D frame
  /// @param cbounds is the conical bound description
  /// @param thickness is the layer thickness along the normal axis
  /// @param od is the overlap description for the sensitive layers
  /// @param ad is the approach descriptor for navigation towards the layer
  /// @param laytyp is the layer type
  ///
  /// @TODO chage od and ad to unique_ptr
  ///
  /// @return is a shared pointer to a layer
  static LayerPtr
  create(std::shared_ptr<Transform3D>      transform,
         std::shared_ptr<const ConeBounds> cbounds,
         std::unique_ptr<SurfaceArray>     surfaceArray,
         double                            thickness = 0.,
         OverlapDescriptor*                od        = nullptr,
         ApproachDescriptor*               ad        = nullptr,
         LayerType                         laytyp    = Acts::active)
  {
    return LayerPtr(new ConeLayer(
        transform, cbounds, std::move(surfaceArray), thickness, od, ad, laytyp));
  }

  /// Factory for shared layer with shift 
  ///
  /// @param cla is the source clone layer
  /// @param shift is the additional shift applied after copying
  ///
  /// @return is a shared pointer to a layer
  static LayerPtr
  create(const ConeLayer& cla, const Transform3D& shift)
  {
    return LayerPtr(new ConeLayer(cla, shift));
  }

  /// Factory for shared layer with shift - clone
  ///
  /// @param shift is the additional shift applied after copying
  ///
  /// @return is a shared pointer to a layer
  LayerPtr
  cloneWithShift(const Transform3D& shift) const override
  {
    return ConeLayer::create(*this, shift);
  }

  /// Default Constructor - delete
  ConeLayer() = delete;
    
  /// Copy constructor of ConeLayer - delete
  ConeLayer(const ConeLayer& cla) = delete;

  /// Assignment operator for ConeLayers - delete
  ConeLayer&
  operator=(const ConeLayer&)
      = delete;

  /// Destructor
  virtual ~ConeLayer() {}
  
  /// Transforms the layer into a Surface representation for extrapolation 
  const ConeSurface&
  surfaceRepresentation() const override;

protected:
  /// Private constructor with arguments
  /// 
  /// @param transform is the 3D transform that poisitions the layer in 3D frame
  /// @param cbounds is the conical bound description
  /// @param thickness is the layer thickness along the normal axis
  /// @param od is the overlap description for the sensitive layers
  /// @param ad is the approach descriptor for navigation towards the layer
  /// @param laytyp is the layer type
  ///
  /// @TODO chage od and ad to unique_ptr
  ConeLayer(std::shared_ptr<Transform3D>      transform,
            std::shared_ptr<const ConeBounds> cbounds,
            std::unique_ptr<SurfaceArray>     surfaceArray,
            double                            thickness = 0.,
            OverlapDescriptor*                od        = nullptr,
            ApproachDescriptor*               ad        = nullptr,
            LayerType                         laytyp    = Acts::active);

  /// Private copy constructor with shift, called by create(args*)
  ///
  /// @param cla is the source clone layer
  /// @param shift is the additional shift applied after copying
  ConeLayer(const ConeLayer& cla, const Transform3D& shift);
};

}  // end of namespace

#endif  // TRKGEOMETY_CONELAYER_H