// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DetectorElementBase.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETELEMENTBASE_DETELEMENTBASE_H
#define ACTS_DETELEMENTBASE_DETELEMENTBASE_H 1

#ifdef ACTS_GEOMETRY_DETELEMENT_PLUGIN
#include ACTS_GEOMETRY_DETELEMENT_PLUGIN
#else

#include <memory>
#include <vector>
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Identifier.hpp"

namespace Acts {

class Surface;
class SurfaceBounds;

/// @class DetectorElementBase
///
/// This is the base class for all tracking detector elements
/// with read-out relevant information.
///
/// If a DetectorElement has a second element (or even a triple setup)
/// that would naturally fall into the same bin, one can register that as a
/// binmember.
///
/// DetectorElements close by can be registered as neighbours as this will help
/// the navigation.
///
class DetectorElementBase
{
public:
  /// Constructor 
  DetectorElementBase() {}
  
  /// virtual Destructor 
  virtual ~DetectorElementBase() {}
  
  /// Identifier 
  virtual Identifier
  identify() const = 0;

  ///  Return local to global transform 
  /// (optionally associated with an identifier)
  virtual const Transform3D&
  transform(const Identifier& identifier = Identifier()) const = 0;

  ///  Return surface association
  /// (optionally associated with an identifier)
  virtual const Surface&
  surface(const Identifier& identifier = Identifier()) const = 0;

  ///  Return full list of associated surfaces
  virtual const std::vector< std::shared_ptr<const Surface> >&
  surfaces() const = 0;

  /// Returns the thickness of the module
  virtual double
  thickness() const = 0;

  /// Fast access to bin members
  const std::vector<const DetectorElementBase*>&
  binmembers() const;

  /// Reigster the bin members
  void
  registerBinmembers(std::vector<const DetectorElementBase*>& binmembers) const;

  /// Fast access to neighbours
  const std::vector<const DetectorElementBase*>&
  neighbours() const;

  /// Reigster the neighbours
  void
  registerNeighbours(std::vector<const DetectorElementBase*>& neighbours) const;

private:
  mutable std::vector<const DetectorElementBase*> m_binmembers;
  mutable std::vector<const DetectorElementBase*> m_neighbours;
};

inline const std::vector<const DetectorElementBase*>&
DetectorElementBase::binmembers() const
{
  return m_binmembers;
}

inline const std::vector<const DetectorElementBase*>&
DetectorElementBase::neighbours() const
{
  return m_neighbours;
}

inline void
DetectorElementBase::registerBinmembers(
    std::vector<const DetectorElementBase*>& binmembers) const
{
  for (auto& bmember : binmembers) {
    // only fill if it's not yet registered
    if (find(m_binmembers.begin(), m_binmembers.end(), bmember)
        == m_binmembers.end())
      m_binmembers.push_back(bmember);
  }
}

inline void
DetectorElementBase::registerNeighbours(
    std::vector<const DetectorElementBase*>& neighbours) const
{
  for (auto& neighbour : neighbours) {
    // only fill if it's not yet registered
    if (find(m_neighbours.begin(), m_neighbours.end(), neighbour)
        == m_neighbours.end())
      m_neighbours.push_back(neighbour);
  }
}

}  // end of ns
#endif

#endif  // ACTS_GEOMETRY_DETELEMENT_PLUGIN