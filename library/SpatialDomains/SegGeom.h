////////////////////////////////////////////////////////////////////////////////
//
//  File:  SegGeom.h
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
//  Description: Segment geometry information
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_SEGGEOM_H
#define NEKTAR_SPATIALDOMAINS_SEGGEOM_H

#include <StdRegions/StdRegions.hpp>
#include <LibUtilities/Foundations/Basis.h>
#include <SpatialDomains/Geometry1D.h>
#include <SpatialDomains/PointGeom.h>
#include <SpatialDomains/Curve.hpp>
#include <SpatialDomains/SpatialDomainsDeclspec.h>

namespace Nektar
{

namespace SpatialDomains
{
class SegGeom;
typedef std::shared_ptr<SegGeom> SegGeomSharedPtr;
typedef std::map<int, SegGeomSharedPtr> SegGeomMap;

class SegGeom : public Geometry1D
{
public:
    SPATIAL_DOMAINS_EXPORT SegGeom();
    SPATIAL_DOMAINS_EXPORT SegGeom(
        int id,
        const int coordim,
        const PointGeomSharedPtr vertex[],
        const CurveSharedPtr curve = CurveSharedPtr());

    SPATIAL_DOMAINS_EXPORT SegGeom(const SegGeom &in);

    SPATIAL_DOMAINS_EXPORT SegGeomSharedPtr GenerateOneSpaceDimGeom(void);

    SPATIAL_DOMAINS_EXPORT ~SegGeom();

    SPATIAL_DOMAINS_EXPORT static StdRegions::Orientation GetEdgeOrientation(
        const SegGeom &edge1, const SegGeom &edge2);

    inline SPATIAL_DOMAINS_EXPORT CurveSharedPtr GetCurve()
    {
        return m_curve;
    }

    SPATIAL_DOMAINS_EXPORT static const int kNverts = 2;
    SPATIAL_DOMAINS_EXPORT static const int kNedges = 1;

protected:
    SpatialDomains::PointGeomSharedPtr m_verts[kNverts];
    StdRegions::Orientation m_porient[kNverts];

    virtual PointGeomSharedPtr v_GetVertex(const int i) const;
    virtual LibUtilities::ShapeType v_GetShapeType() const;
    virtual NekDouble v_GetLocCoords(
        const Array<OneD, const NekDouble> &coords,
        Array<OneD, NekDouble> &Lcoords);
    virtual void v_GenGeomFactors();
    virtual void v_FillGeom();
    virtual void v_Reset(CurveMap &curvedEdges, CurveMap &curvedFaces);
    virtual void v_Setup();

    virtual NekDouble v_GetCoord(
        const int i, const Array<OneD, const NekDouble> &Lcoord);
    virtual int v_GetNumVerts() const;
    virtual int v_GetNumEdges() const;
    virtual bool v_ContainsPoint(
        const Array<OneD, const NekDouble> &gloCoord,
        Array<OneD, NekDouble> &locCoord,
        NekDouble tol,
        NekDouble &resid);

private:
    /// Boolean indicating whether object owns the data
    CurveSharedPtr m_curve;

    void SetUpXmap();
};
} // end of namespace
} // end of namespace

#endif // NEKTAR_SPATIALDOMAINS_SEGGEOM_H
