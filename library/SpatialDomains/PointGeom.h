////////////////////////////////////////////////////////////////////////////////
//
//  File:  PointGeom.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Point geometry information
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_POINTGEOM_H
#define NEKTAR_SPATIALDOMAINS_POINTGEOM_H

#include <StdRegions/StdRegions.hpp>
#include <SpatialDomains/Geometry0D.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
#include <memory>

namespace Nektar
{

namespace SpatialDomains
{

class PointGeom;
class SegGeom;

// shorthand for boost pointer
typedef std::shared_ptr<PointGeom> PointGeomSharedPtr;
typedef std::map<int, PointGeomSharedPtr> PointGeomMap;

class PointGeom : public Geometry0D,
                  public NekPoint<NekDouble>,
                  public std::enable_shared_from_this<PointGeom>
{
public:
    SPATIAL_DOMAINS_EXPORT PointGeom();
    SPATIAL_DOMAINS_EXPORT PointGeom(const int coordim,
                                     const int vid,
                                     NekDouble x,
                                     NekDouble y,
                                     NekDouble z);
    SPATIAL_DOMAINS_EXPORT PointGeom(const PointGeom &T);

    SPATIAL_DOMAINS_EXPORT ~PointGeom();

    SPATIAL_DOMAINS_EXPORT int GetVid()
    {
        return m_globalID;
    }

    SPATIAL_DOMAINS_EXPORT void GetCoords(NekDouble &x,
                                          NekDouble &y,
                                          NekDouble &z);
    SPATIAL_DOMAINS_EXPORT void GetCoords(Array<OneD, NekDouble> &coords);
    SPATIAL_DOMAINS_EXPORT void UpdatePosition(NekDouble x,
                                               NekDouble y,
                                               NekDouble z);

    SPATIAL_DOMAINS_EXPORT void Mult(PointGeom &a, PointGeom &b);
    SPATIAL_DOMAINS_EXPORT void Add(PointGeom &a, PointGeom &b);
    SPATIAL_DOMAINS_EXPORT void Sub(PointGeom &a, PointGeom &b);
    SPATIAL_DOMAINS_EXPORT void Rotate (PointGeom& a, int dir, NekDouble angle);
    SPATIAL_DOMAINS_EXPORT NekDouble dist(PointGeom &a);
    SPATIAL_DOMAINS_EXPORT NekDouble dot(PointGeom &a);

    SPATIAL_DOMAINS_EXPORT friend bool operator==(const PointGeom &x,
                                                  const PointGeom &y);
    SPATIAL_DOMAINS_EXPORT friend bool operator==(const PointGeom &x,
                                                  const PointGeom *y);
    SPATIAL_DOMAINS_EXPORT friend bool operator==(const PointGeom *x,
                                                  const PointGeom &y);
    SPATIAL_DOMAINS_EXPORT friend bool operator!=(const PointGeom &x,
                                                  const PointGeom &y);
    SPATIAL_DOMAINS_EXPORT friend bool operator!=(const PointGeom &x,
                                                  const PointGeom *y);
    SPATIAL_DOMAINS_EXPORT friend bool operator!=(const PointGeom *x,
                                                  const PointGeom &y);

protected:
    virtual void v_GenGeomFactors();
    virtual PointGeomSharedPtr v_GetVertex(int i) const;
};

}
}

#endif // NEKTAR_SPATIALDOMAINS_POINTGEOM_H
