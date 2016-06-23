// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// Layer.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETECTOR_LAYER_H
#define ACTS_DETECTOR_LAYER_H 1

// Core module
#include "ACTS/EventData/NeutralParameters.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Utilities/ApproachDescriptor.hpp"
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/GeometryObject.hpp"
#include "ACTS/Utilities/Intersection.hpp"
#include "ACTS/Utilities/OverlapDescriptor.hpp"
#include "ACTS/Volumes/AbstractVolume.hpp"

namespace Acts {

class Surface;
class SurfaceMaterial;
class MaterialProperties;
class BinUtility;
class Volume;
class VolumeBounds;
class TrackingVolume;
class DetachedTrackingVolume;
class ApproachDescriptor;
class ICompatibilityEstimator;

typedef ObjectIntersection<Surface> SurfaceIntersection;

// master typedef
class Layer;
typedef std::shared_ptr<const Layer> LayerPtr;
typedef std::pair<const Layer*, const Layer*> NextLayers;

/// 
/// @enum LayerType
/// 
/// For readability
/// 
enum LayerType { passive = 0, active = 1 };
 
/// @class Layer
/// 
/// Base Class for a Detector Layer in the Tracking realm.
/// An actual implemented Detector Layer inheriting from this base
/// class has to inherit from a specific type of Surface as well.
/// In addition, a Layer can carry:
/// 
/// - a SurfaceArray of Surfaces holding the actual detector elements or
/// subSurfaces.
/// - SurfaceMaterial for Surface-based materialUpdates
/// - an OverlapDescriptor (mainly used for blind extrapolation)
/// - a pointer to the TrackingVolume (can only be set by such)
/// - an active/passive code :
/// 0      - activ
/// 1      - passive
/// [....] - other
/// 
/// The search type for compatible surfaces on a layer is [ the higher the
/// number, the faster ]:
/// --------- Layer internal ------------------------------------------------
/// -1     - untested: provide all layer surfaces to the extrapolation engine
///               - does not work with endSurface, will be increased to 0 if
/// endSurface is given
///               - debug mode only !
///  0     - test all on intersection and provide to the extrapolation engine
/// --------- Overlap descriptor --------------------------------------------
///  1     - provide bin surface and registered neighbours and bin mates
///               - does not work with endSurface, will be increased to 2 if
/// endSurface is given
///  2     - as 1 but with intersection test @TODO compatibility test
///  3     - provide bin surface and next bin surfaces (if differ)
///               - does not work with endSurface, will be increased to 4 if
/// endSurface is given
///  4     - as 3 but with intersection test @TODO compatibility test
///  5     - whatever the overlap descriptor returns with this
/// 
/// 
class Layer : public virtual GeometryObject
{
  /// Declare the TrackingVolume as a friend, to be able to register previous,
  /// next and set the enclosing TrackingVolume
  friend class TrackingVolume;

  /// Declare the DetachedTrackingVolume as a friend to be able to register it
  friend class DetachedTrackingVolume;

public:
  /// Clone at a with a shift - this is the only clone allowed 
  /// @param shift is the additional transform applied after cloning
  virtual LayerPtr
  cloneWithShift(const Transform3D& shift) const = 0;

  /// Destructor
  virtual ~Layer();

  /// Assignment operator - forbidden, layer assignment must not be ambiguous 
  Layer&
  operator=(const Layer& lay)
      = delete;

  /// Return the entire SurfaceArray, returns a nullptr if no SurfaceArray 
  const SurfaceArray*
  surfaceArray() const;

  /// Transforms the layer into a Surface representation for extrapolation 
  /// @note the layer can be hosting many surfaces, but this is the global
  /// one to which one can extrapolate
  virtual const Surface&
  surfaceRepresentation() const = 0;

  /// Return the Thickness of the Layer
  /// this is by definition along the normal vector of the surfaceRepresentation
  double
  thickness() const;

  /// templated onLayer() method
  /// @paramt parameters are the templated (charged/neutral) on layer check
  /// @param bchk is the boundary check directive
  template <class T>
  bool
  onLayer(const T&             parameters,
          const BoundaryCheck& bchk = true) const;

  /// geometrical isOnLayer() method
  /// @note using isOnSurface() with Layer specific tolerance 
  /// @param gpos is the gobal position to be checked
  /// @param bchk is the boundary check directive       
  virtual bool
  isOnLayer(const Vector3D&      gp,
            const BoundaryCheck& bchk = true) const;

  /// Return method for the overlap descriptor, can be nullptr        
  const OverlapDescriptor*
  overlapDescriptor() const;

  /// Return method for the approach descriptor, can be nullptr
  const ApproachDescriptor*
  approachDescriptor() const;

  ///  Surface seen on approach 
  /// for surcfaces without sub structure, this is the surfaceRepresentation
  /// @param gpos is the global position to start the approach from 
  /// @param dir is the direction from where to attempt the approach
  /// @param pdir is the direction prescription
  /// @param bhck is the boundary check directive
  /// @param resolveSubSurfaces is a boolean directive whether to resolve structure or not
  /// @note reasons for resolving are: collect & find material, collect & find sensitivex
  /// @ice is a (future) compatibility estimator that could be used to modify the straight line approach
  virtual const SurfaceIntersection
  surfaceOnApproach(const Vector3D&                gpos,
                    const Vector3D&                dir,
                    PropDirection                  pdir,
                    const BoundaryCheck&           bchk,
                    bool                           resolveSubSurfaces = false,
                    const ICompatibilityEstimator* ice = nullptr) const;

  ///  get compatible surfaces starting from charged parameters
  ///  returns back the compatible surfaces either with straight line estimation,
  ///  or (@TODO later) with a compatiblityEstimator.
  ///  - if start/end surface are given, surfaces are provided in between (start
  /// & end excluded)
  ///  - the boolean indicates if the surfaces are direction ordered
  ///
  /// @param cSurfaces are the retrun surface intersections
  /// @param pars are the (charged) track parameters for the search
  /// @param pdir is the propagation direction prescription
  /// @param bchk is the boundary check directive                  
  /// @param collectSensitive is the prescription to find the sensitive surfaces
  /// @param collectPassive is the prescription to find all passive surfaces
  /// @param searchType is the level of depth for the search
  /// @param startSurface is an (optional) start surface for the search: excluded in return
  /// @param endSurface is an (optional) end surface for the search: excluded in return
  /// @ice is a (future) compatibility estimator that could be used to modify the straight line approach                               
  virtual bool
  compatibleSurfaces(std::vector<SurfaceIntersection>& cSurfaces,
                     const TrackParameters&            pars,
                     PropDirection                     pdir,
                     const BoundaryCheck&              bchk,
                     bool                              collectSensitive,
                     bool                              collectPassive,
                     int                               searchType,
                     const Surface*                    startSurface = nullptr,
                     const Surface*                    endSurface   = nullptr,
                     const ICompatibilityEstimator*    ice = nullptr) const;

  ///  get compatible surfaces starting from neutral parameters
  ///  returns back the compatible surfaces either with straight line estimation,
  ///  or (@TODO later) with a compatiblityEstimator.
  ///  - if start/end surface are given, surfaces are provided in between (start
  /// & end excluded)
  ///  - the boolean indicates if the surfaces are direction ordered
  ///
  /// @param cSurfaces are the retrun surface intersections
  /// @param pars are the (charged) track parameters for the search
  /// @param pdir is the propagation direction prescription
  /// @param bchk is the boundary check directive                  
  /// @param collectSensitive is the prescription to find the sensitive surfaces
  /// @param collectPassive is the prescription to find all passive surfaces
  /// @param searchType is the level of depth for the search
  /// @param startSurface is an (optional) start surface for the search: excluded in return
  /// @param endSurface is an (optional) end surface for the search: excluded in return
  /// @ice is a (future) compatibility estimator that could be used to modify the straight line approach  
  virtual bool
  compatibleSurfaces(std::vector<SurfaceIntersection>& cSurfaces,
                     const NeutralParameters&          pars,
                     PropDirection                     pdir,
                     const BoundaryCheck&              bchk,
                     bool                              collectSensitive,
                     bool                              collectPassive,
                     int                               searchType,
                     const Surface*                    startSurface = nullptr,
                     const Surface*                    endSurface   = nullptr,
                     const ICompatibilityEstimator*    ice = nullptr) const;

  /// Has sub-structure method:
  /// @note sub-structure depending on :
  ///   (a) only when required to resolve sub surfaces for sensitive hits
  ///   (b) also material is ordered with sub structure 
  virtual bool
  hasSubStructure(bool resolveSensitive = false) const;

  //  Boolean check method if layer has material:
  // - checks if any of the layer surfaces has material:
  // - can be approach surfaces or layer surface
  virtual bool
  hasMaterial() const;

  /// Boolean check method if layer has sensitive surfaces 
  /// @note checks if a surfaceArray is present
  virtual bool
  hasSensitive() const;

  /// Fast navigation to next layer
  /// @param gpos is the start position for the search
  /// @param dir is the direction for the search
  const Layer*
  nextLayer(const Vector3D& pos, const Vector3D& dir) const;

  /// get the confining TrackingVolume 
  const TrackingVolume*
  enclosingTrackingVolume() const;

  /// get the confining DetachedTrackingVolume 
  const DetachedTrackingVolume*
  enclosingDetachedTrackingVolume() const;

  /// register Volume associated to the layer - if you want to do that by hand
  void
  registerRepresentingVolume(const AbstractVolume* theVol) const;

  ///  return the abstract volume that represents the layer
  const AbstractVolume*
  representingVolume() const;

protected:
  /// Default Constructor
  Layer();

  /// Copy Constructor 
  Layer(const Layer& lay);

  /// Constructor with pointer to SurfaceArray (passing ownership)
  /// @param surfaceArray is the array of sensitive surfaces
  /// @param thickness is the normal thickness of the Layer
  /// @param od overlap descriptor (@TODO change to unique_ptr)  
  /// @param ad approach descriptor (@TODO change to unique_ptr)  
  /// 
  Layer(std::unique_ptr<SurfaceArray> surfaceArray,
        double                        thickness = 0.,
        OverlapDescriptor*            od        = nullptr,
        ApproachDescriptor*           ad        = nullptr,
        LayerType                     ltype     = Acts::passive);

  /// get compatible surfaces starting from charged parameters - forward call
  /// from explicit methods
  template <class T>
  bool
  getCompatibleSurfaces(std::vector<SurfaceIntersection>& cSurfaces,
                        const T&                          pars,
                        PropDirection                     pdir,
                        const BoundaryCheck&              bchk,
                        bool                              collectSensitive,
                        bool                              collectPassive,
                        int                               searchType,
                        const Surface*                    startSurface = nullptr,
                        const Surface*                    endSurface   = nullptr,
                        const ICompatibilityEstimator*    ice = nullptr) const;

  /// test compatible surface - checking directly for intersection & collection
  /// geometrical test compatible surface method
  void
  testCompatibleSurface(std::vector<SurfaceIntersection>& cSurfaces,
                        const Surface&                    surface,
                        const Vector3D&                   gpos,
                        const Vector3D&                   dir,
                        PropDirection                     pdir,
                        const BoundaryCheck&              bchk,
                        double                            maxPathLength,
                        bool                              collectSensitive,
                        bool                              collectPassive,
                        bool                              intersectionTest,
                        const Surface*                    startSurface = nullptr,
                        const Surface*                    endSurface   = nullptr,
                        const ICompatibilityEstimator*    ice = nullptr) const;

  ///  private method to set enclosing TrackingVolume, called by friend class
  /// only
  ///  optionally, the layer can be resized to the dimensions of the
  /// TrackingVolume
  ///  - Bounds of the Surface are resized
  ///  - MaterialProperties dimensions are resized
  ///  - SubSurface array boundaries are NOT resized
  /// 
  /// @param tvol is the tracking volume the layer is confined
  void
  encloseTrackingVolume(const TrackingVolume& tvol) const;

  ///  private method to set the enclosed detached TV,
  /// called by friend class only 
  /// @param dtvol is the detached tracking volume the layer is confined
  void
  encloseDetachedTrackingVolume(const DetachedTrackingVolume& tvol) const;

  //// the previous Layer according to BinGenUtils
  mutable NextLayers                    m_nextLayers; 
  //// A binutility to find the next layer
  mutable const BinUtility*             m_nextLayerUtility; 
  //// SurfaceArray on this layer Surface
  std::unique_ptr<SurfaceArray>         m_surfaceArray;    
  //// thickness of the Layer
  double                                m_layerThickness;  
  //// descriptor for overlap/next surface
  OverlapDescriptor*                    m_overlapDescriptor;
  //// descriptor for surface on approach
  //// @note this is a mutable member, since resize may trigger to 
  mutable ApproachDescriptor*           m_approachDescriptor;
  //// the enclosing TrackingVolume
  //// @note this is a mutable member since it's set after the layer creation
  mutable const TrackingVolume*         m_enclosingTrackingVolume;  
  //// the eventual enclosing detached Tracking volume
  //// @note this is a mutable member since it is set after layer creation
  mutable const DetachedTrackingVolume* m_enclosingDetachedTrackingVolume;  
  //// Representing Volume
  //// can be used as appraoch suraces 
  mutable const AbstractVolume*         m_representingVolume; 
  //// make a passive/active divisio
  LayerType                             m_layerType;
};

inline const SurfaceArray*
Layer::surfaceArray() const
{
  return m_surfaceArray.get();
}

inline double
Layer::thickness() const
{
  return m_layerThickness;
}

inline const OverlapDescriptor*
Layer::overlapDescriptor() const
{
  return m_overlapDescriptor;
}

inline const TrackingVolume*
Layer::enclosingTrackingVolume() const
{
  return m_enclosingTrackingVolume;
}

inline void
Layer::encloseTrackingVolume(const TrackingVolume& tvol) const
{
  m_enclosingTrackingVolume = &(tvol);
}

inline const DetachedTrackingVolume*
Layer::enclosingDetachedTrackingVolume() const
{
  return m_enclosingDetachedTrackingVolume;
}

inline void
Layer::encloseDetachedTrackingVolume(const DetachedTrackingVolume& tvol) const
{
  m_enclosingDetachedTrackingVolume = &(tvol);
}

inline const AbstractVolume*
Layer::representingVolume() const
{
  return m_representingVolume;
}

inline const Layer*
Layer::nextLayer(const Vector3D& gp, const Vector3D& mom) const
{
  // no binutility -> no chance to find out the direction
  if (!m_nextLayerUtility) return nullptr;
  return (m_nextLayerUtility->nextDirection(gp, mom) < 0) ? m_nextLayers.first
                                                          : m_nextLayers.second;
}

inline void
Layer::registerRepresentingVolume(const AbstractVolume* theVol) const
{
  delete m_representingVolume;
  m_representingVolume = theVol;
}

#include "ACTS/Layers/detail/Layer.ipp"

/// Layers are constructedd with shared_ptr factories, hence the layer array is
/// describes as:
typedef BinnedArray<LayerPtr> LayerArray;

}  // end of namespace

#endif  // ACTS_DETECTOR_LAYER_H