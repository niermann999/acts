// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// VolumeBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_VOLUMES_VOLUMEBOUNDS_H
#define ACTS_VOLUMES_VOLUMEBOUNDS_H 1

#include <iomanip>
#include <iostream>
#include <memory>
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/BinningType.hpp"

namespace Acts {

class Surface;
class Volume;

class VolumeBounds;
typedef std::shared_ptr<const VolumeBounds> VolumeBoundsPtr;

/// @class VolumeBounds
///
/// Pure Absract Base Class for Volume bounds.
///
/// Acts::VolumeBounds are a set of up to six confining Surfaces that are stored in
/// a std::vector.
/// Each type of Acts::VolumeBounds has to implement a decomposeToSurfaces() and a
/// inside() method.
///
/// The orientation of the Surfaces are in a way that the normal vector points to
/// the oustide world.
///
/// The Volume, retrieving a set of Surfaces from the VolumeBounds, can turn the
/// Surfaces into BoundarySurfaces.

class VolumeBounds
{
public:
  
  /// Default Constructor*/
  VolumeBounds() {}
  
  /// Destructor
  virtual ~VolumeBounds() {}
  
  ///  clone() method to make deep copy in Volume copy constructor and for
  /// assigment operator  of the Surface class.
  virtual VolumeBounds*
  clone() const = 0;

  /// Checking if position given in volume frame is inside
  /// @param gpos is the global position to be checked
  /// @param tol is the tolerance applied for the inside check
  virtual bool
  inside(const Vector3D& gpos, double tol = 0.) const = 0;

  /// Method to decompose the Bounds into Surfaces
  /// the Volume can turn them into BoundarySurfaces
  /// @param transform is the 3D transform to be applied to the boundary surfaces
  ///        to position them in 3D space
  /// @note this is factory method
  virtual const std::vector<const Surface*>
  decomposeToSurfaces(std::shared_ptr<Transform3D> transform) const = 0;

  /// Binning offset - overloaded for some R-binning types 
  virtual Vector3D
  binningOffset(BinningValue bValue) const;

  /// Binning borders in double 
  virtual double
  binningBorder(BinningValue bValue) const;
  
  /// Output Method for std::ostream, to be overloaded by child classes 
  virtual std::ostream&
  dump(std::ostream& sl) const = 0;
};

/// Binning offset - overloaded for some R-binning types 
inline Vector3D VolumeBounds::binningOffset(BinningValue) const
{  // standard offset is 0.,0.,0.
  return Vector3D(0., 0., 0.);
}

inline double VolumeBounds::binningBorder(BinningValue) const
{
  return 0.;
}

/// Overload of << operator for std::ostream for debug output
std::ostream&
operator<<(std::ostream& sl, const VolumeBounds& vb);

}  // end of namespace

#endif  // ACTS_VOLUMES_VOLUMEBOUNDS_H