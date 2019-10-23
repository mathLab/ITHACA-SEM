////////////////////////////////////////////////////////////////////////////////
//
//  File: QuadGeom.h
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
//  Description:
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_QUADGEOM_H
#define NEKTAR_SPATIALDOMAINS_QUADGEOM_H

#include <StdRegions/StdRegions.hpp>

#include <SpatialDomains/Geometry2D.h>
#include <SpatialDomains/PointGeom.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
#include <SpatialDomains/SegGeom.h>

namespace Nektar
{
namespace SpatialDomains
{

class QuadGeom;
struct Curve;
typedef std::shared_ptr<Curve> CurveSharedPtr;
typedef std::shared_ptr<QuadGeom> QuadGeomSharedPtr;
typedef std::map<int, QuadGeomSharedPtr> QuadGeomMap;

class QuadGeom : public Geometry2D
{
public:
    SPATIAL_DOMAINS_EXPORT QuadGeom();
    SPATIAL_DOMAINS_EXPORT QuadGeom(const QuadGeom &in);
    SPATIAL_DOMAINS_EXPORT QuadGeom(
        const int id,
        const SegGeomSharedPtr edges[],
        const CurveSharedPtr curve = CurveSharedPtr());
    SPATIAL_DOMAINS_EXPORT ~QuadGeom();

    /// Get the orientation of face1.
    SPATIAL_DOMAINS_EXPORT static StdRegions::Orientation GetFaceOrientation(
                    const QuadGeom &face1, const QuadGeom &face2,
                    bool doRot = false, int dir = 0, NekDouble angle = 0.0,
                    NekDouble tol = 1e-8);
    SPATIAL_DOMAINS_EXPORT static StdRegions::Orientation GetFaceOrientation(
                    const PointGeomVector &face1, const PointGeomVector &face2,
                    bool doRot = false, int dir = 0, NekDouble angle = 0.0,
                    NekDouble tol = 1e-8);

    SPATIAL_DOMAINS_EXPORT static const int kNverts = 4;
    SPATIAL_DOMAINS_EXPORT static const int kNedges = 4;
    SPATIAL_DOMAINS_EXPORT static const std::string XMLElementType;

protected:
    SPATIAL_DOMAINS_EXPORT virtual NekDouble v_GetCoord(
        const int i, const Array<OneD, const NekDouble> &Lcoord);
    SPATIAL_DOMAINS_EXPORT void v_GenGeomFactors();
    SPATIAL_DOMAINS_EXPORT virtual void v_FillGeom();
    SPATIAL_DOMAINS_EXPORT virtual NekDouble v_GetLocCoords(
        const Array<OneD, const NekDouble> &coords,
        Array<OneD, NekDouble> &Lcoords);
    SPATIAL_DOMAINS_EXPORT virtual bool v_ContainsPoint(
        const Array<OneD, const NekDouble> &gloCoord,
        Array<OneD, NekDouble> &locCoord,
        NekDouble tol,
        NekDouble &resid);
    SPATIAL_DOMAINS_EXPORT virtual void v_Reset(CurveMap &curvedEdges,
                                                CurveMap &curvedFaces);
    SPATIAL_DOMAINS_EXPORT virtual void v_Setup();

private:
    void SetUpXmap();
};

}
}

#endif // NEKTAR_SPATIALDOMAINS_QUADGEOM_H
