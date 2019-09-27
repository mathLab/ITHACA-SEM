////////////////////////////////////////////////////////////////////////////////
//
//  File: Geometry2D.h
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
//  Description: 2D geometry information
//
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_GEOMETRY2D_H
#define NEKTAR_SPATIALDOMAINS_GEOMETRY2D_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdExpansion2D.h>

#include <SpatialDomains/Geometry.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>

namespace Nektar
{

namespace SpatialDomains
{

class Geometry0D;
class Geometry1D;
class Geometry2D;
class PointGeom;
class SegGeom;

// shorthand for boost pointer
typedef std::shared_ptr<PointGeom> PointGeomSharedPtr;
typedef std::shared_ptr<Geometry0D> Geometry0DSharedPtr;
typedef std::shared_ptr<Geometry1D> Geometry1DSharedPtr;
typedef std::shared_ptr<Geometry2D> Geometry2DSharedPtr;
typedef std::shared_ptr<SegGeom> SegGeomSharedPtr;
typedef std::vector<Geometry2DSharedPtr> Geometry2DVector;
typedef std::vector<PointGeomSharedPtr> PointGeomVector;
typedef std::vector<SegGeomSharedPtr> SegGeomVector;

/// 2D geometry information
class Geometry2D : public Geometry
{
public:
    SPATIAL_DOMAINS_EXPORT Geometry2D();
    SPATIAL_DOMAINS_EXPORT Geometry2D(const int coordim, CurveSharedPtr curve);
    SPATIAL_DOMAINS_EXPORT virtual ~Geometry2D();

    SPATIAL_DOMAINS_EXPORT static const int kDim = 2;

protected:
    PointGeomVector                      m_verts;
    SegGeomVector                        m_edges;
    std::vector<StdRegions::Orientation> m_eorient;
    CurveSharedPtr                       m_curve;

    void NewtonIterationForLocCoord(const Array<OneD, const NekDouble> &coords,
                                    const Array<OneD, const NekDouble> &ptsx,
                                    const Array<OneD, const NekDouble> &ptsy,
                                    Array<OneD, NekDouble> &Lcoords,
                                    NekDouble &resid);

private:
    //---------------------------------------
    // Helper functions
    //---------------------------------------
    virtual int v_GetShapeDim() const;
    virtual PointGeomSharedPtr v_GetVertex(int i) const;
    virtual Geometry1DSharedPtr v_GetEdge(int i) const;
    virtual int v_GetNumVerts() const;
    virtual int v_GetNumEdges() const;
    virtual StdRegions::Orientation v_GetEorient(const int i) const;
};

}
}

#endif // NEKTAR_SPATIALDOMAINS_GEOMETRY2D_H
