////////////////////////////////////////////////////////////////////////////////
//
//  File: PrismGeom.h
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
//  Description: Prismatic geometry definition.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_PRISMGEOM_H
#define NEKTAR_SPATIALDOMAINS_PRISMGEOM_H

#include <SpatialDomains/Geometry3D.h>
#include <SpatialDomains/QuadGeom.h>
#include <SpatialDomains/TriGeom.h>

namespace Nektar
{
namespace SpatialDomains
{

class PrismGeom : public Geometry3D
{
public:
    SPATIAL_DOMAINS_EXPORT PrismGeom();
    SPATIAL_DOMAINS_EXPORT PrismGeom(int id, const Geometry2DSharedPtr faces[]);
    SPATIAL_DOMAINS_EXPORT ~PrismGeom();

    SPATIAL_DOMAINS_EXPORT static const int kNverts = 6;
    SPATIAL_DOMAINS_EXPORT static const int kNedges = 9;
    SPATIAL_DOMAINS_EXPORT static const int kNqfaces = 3;
    SPATIAL_DOMAINS_EXPORT static const int kNtfaces = 2;
    SPATIAL_DOMAINS_EXPORT static const int kNfaces = kNqfaces + kNtfaces;
    SPATIAL_DOMAINS_EXPORT static const std::string XMLElementType;

protected:
    virtual void v_GenGeomFactors();
    virtual NekDouble v_GetLocCoords(const Array<OneD, const NekDouble> &coords,
                                     Array<OneD, NekDouble> &Lcoords);
    virtual bool v_ContainsPoint(const Array<OneD, const NekDouble> &gloCoord,
                                 Array<OneD, NekDouble> &locCoord,
                                 NekDouble tol,
                                 NekDouble &resid);
    virtual int v_GetVertexEdgeMap(const int i, const int j) const;
    virtual int v_GetVertexFaceMap(const int i, const int j) const;
    virtual int v_GetEdgeFaceMap(const int i, const int j) const;
    virtual int v_GetDir(const int faceidx, const int facedir) const;
    virtual void v_Reset(CurveMap &curvedEdges, CurveMap &curvedFaces);
    virtual void v_Setup();

private:
    void SetUpLocalEdges();
    void SetUpLocalVertices();
    void SetUpEdgeOrientation();
    void SetUpFaceOrientation();
    void SetUpXmap();

    static const unsigned int VertexEdgeConnectivity[6][3];
    static const unsigned int VertexFaceConnectivity[6][3];
    static const unsigned int EdgeFaceConnectivity[9][2];
};

typedef std::shared_ptr<PrismGeom> PrismGeomSharedPtr;
typedef std::map<int, PrismGeomSharedPtr> PrismGeomMap;
}
}

#endif // NEKTAR_SPATIALDOMAINS_PRISMGEOM_H
