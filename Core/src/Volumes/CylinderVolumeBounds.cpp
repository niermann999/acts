// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderVolumeBounds.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Volumes/CylinderVolumeBounds.hpp"
#include "ACTS/Surfaces/CylinderBounds.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Surfaces/DiscSurface.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include <iostream>
#include <math.h>

double Acts::CylinderVolumeBounds::s_numericalStable = 10e-2;

Acts::CylinderVolumeBounds::CylinderVolumeBounds()
  : VolumeBounds(), m_valueStore(4, 0.)
{
}

Acts::CylinderVolumeBounds::CylinderVolumeBounds(double radius, double halez)
  : VolumeBounds(), m_valueStore(4, 0.)
{
  m_valueStore.at(bv_innerRadius)   = 0.;
  m_valueStore.at(bv_outerRadius)   = fabs(radius);
  m_valueStore.at(bv_halfPhiSector) = M_PI;
  m_valueStore.at(bv_halfZ)         = fabs(halez);
}

Acts::CylinderVolumeBounds::CylinderVolumeBounds(double rinner,
                                                 double router,
                                                 double halez)
  : VolumeBounds(), m_valueStore(4, 0.)
{
  m_valueStore.at(bv_innerRadius)   = fabs(rinner);
  m_valueStore.at(bv_outerRadius)   = fabs(router);
  m_valueStore.at(bv_halfPhiSector) = M_PI;
  m_valueStore.at(bv_halfZ)         = fabs(halez);
}

Acts::CylinderVolumeBounds::CylinderVolumeBounds(double rinner,
                                                 double router,
                                                 double haphi,
                                                 double halez)
  : VolumeBounds(), m_valueStore(4, 0.)
{
  m_valueStore.at(bv_innerRadius)   = fabs(rinner);
  m_valueStore.at(bv_outerRadius)   = fabs(router);
  m_valueStore.at(bv_halfPhiSector) = fabs(haphi);
  m_valueStore.at(bv_halfZ)         = fabs(halez);
}

Acts::CylinderVolumeBounds::CylinderVolumeBounds(
    const CylinderVolumeBounds& cylbo)
  : VolumeBounds(), m_valueStore(cylbo.m_valueStore)
{
}

Acts::CylinderVolumeBounds::~CylinderVolumeBounds()
{
}

Acts::CylinderVolumeBounds&
Acts::CylinderVolumeBounds::operator=(const CylinderVolumeBounds& cylbo)
{
  if (this != &cylbo) m_valueStore = cylbo.m_valueStore;
  return *this;
}

const std::vector<const Acts::Surface*>
Acts::CylinderVolumeBounds::decomposeToSurfaces(
    std::shared_ptr<Transform3D> transformPtr) const
{
  std::vector<const Surface*> rSurfaces;
  rSurfaces.reserve(6);

  // set the transform
  Transform3D transform = (transformPtr == nullptr)
      ? Transform3D::Identity()
      : (*transformPtr.get());
  Transform3D*     tTransform = nullptr;
  RotationMatrix3D discRot(transform.rotation());
  Vector3D         cylCenter(transform.translation());

  // bottom Disc (negative z)
  RotationMatrix3D bottomDiscRot;
  bottomDiscRot.col(0) = discRot.col(1);
  bottomDiscRot.col(1) = discRot.col(0);
  bottomDiscRot.col(2) = -discRot.col(2);

  std::shared_ptr<const DiscBounds> dBounds = discBounds();
  tTransform = new Transform3D(
      transform * AngleAxis3D(M_PI, Vector3D(1., 0., 0.))
      * Translation3D(Vector3D(0., 0., halflengthZ())));
  rSurfaces.push_back(new DiscSurface(
      std::shared_ptr<Transform3D>(tTransform), dBounds));
  // top Disc (positive z)
  tTransform = new Transform3D(
      discRot
      * Translation3D(cylCenter + halflengthZ() * discRot.col(2)));
  rSurfaces.push_back(new DiscSurface(
      std::shared_ptr<Transform3D>(tTransform), dBounds));

  // outer Cylinder - shares the transform
  rSurfaces.push_back(
      new CylinderSurface(transformPtr, outerCylinderBounds()));

  // innermost Cylinder
  if (innerRadius() > s_numericalStable)
    rSurfaces.push_back(
        new CylinderSurface(transformPtr, innerCylinderBounds()));

  // the cylinder is sectoral
  if (fabs(halfPhiSector() - M_PI) > s_numericalStable) {
    std::shared_ptr<const PlanarBounds> sp12Bounds = sectorPlaneBounds();
    // sectorPlane 1 (negative phi)
    Transform3D* sp1Transform = new Transform3D(
        transform
        * AngleAxis3D(-halfPhiSector(), Vector3D(0., 0., 1.))
        * Translation3D(Vector3D(mediumRadius(), 0., 0.))
        * AngleAxis3D(M_PI / 2, Vector3D(1., 0., 0.)));
    rSurfaces.push_back(new PlaneSurface(
        std::shared_ptr<Transform3D>(sp1Transform), sp12Bounds));
    // sectorPlane 2 (positive phi)
    Transform3D* sp2Transform = new Transform3D(
        transform
        * AngleAxis3D(halfPhiSector(), Vector3D(0., 0., 1.))
        * Translation3D(Vector3D(mediumRadius(), 0., 0.))
        * AngleAxis3D(-M_PI / 2, Vector3D(1., 0., 0.)));
    rSurfaces.push_back(new PlaneSurface(
        std::shared_ptr<Transform3D>(sp2Transform), sp12Bounds));
  }
  return rSurfaces;
}

std::shared_ptr<const Acts::CylinderBounds>
Acts::CylinderVolumeBounds::innerCylinderBounds() const
{
  return std::make_shared<const CylinderBounds>(m_valueStore.at(bv_innerRadius),
                                                m_valueStore.at(bv_halfPhiSector),
                                                m_valueStore.at(bv_halfZ));
}

std::shared_ptr<const Acts::CylinderBounds>
Acts::CylinderVolumeBounds::outerCylinderBounds() const
{
  return std::make_shared<const CylinderBounds>(m_valueStore.at(bv_outerRadius),
                                                m_valueStore.at(bv_halfPhiSector),
                                                m_valueStore.at(bv_halfZ));
}

std::shared_ptr<const Acts::DiscBounds>
Acts::CylinderVolumeBounds::discBounds() const
{
  return std::shared_ptr<const DiscBounds>(new RadialBounds(m_valueStore.at(bv_innerRadius),
                                                            m_valueStore.at(bv_outerRadius),
                                                            m_valueStore.at(bv_halfPhiSector)));
}

std::shared_ptr<const Acts::PlanarBounds>
Acts::CylinderVolumeBounds::sectorPlaneBounds() const
{
  return std::shared_ptr<const PlanarBounds>(new RectangleBounds(0.5 * (m_valueStore.at(bv_outerRadius)
                                                - m_valueStore.at(bv_innerRadius)),
                                                m_valueStore.at(bv_halfZ)));
}

std::ostream&
Acts::CylinderVolumeBounds::dump(std::ostream& sl) const
{
  return dumpT<std::ostream>(sl);
}