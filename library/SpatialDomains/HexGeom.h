////////////////////////////////////////////////////////////////////////////////
//
//  File: HexGeom.h
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
//  Description: Hexahedral geometry information.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_HEXGEOM_H
#define NEKTAR_SPATIALDOMAINS_HEXGEOM_H

#include <SpatialDomains/Geometry3D.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>

namespace Nektar
{
namespace SpatialDomains
{

class QuadGeom;
typedef std::shared_ptr<QuadGeom> QuadGeomSharedPtr;

class HexGeom : public Geometry3D
{
public:
    SPATIAL_DOMAINS_EXPORT HexGeom();
    SPATIAL_DOMAINS_EXPORT HexGeom(int id, const QuadGeomSharedPtr faces[]);
    SPATIAL_DOMAINS_EXPORT ~HexGeom();

    SPATIAL_DOMAINS_EXPORT static const int kNverts = 8;
    SPATIAL_DOMAINS_EXPORT static const int kNedges = 12;
    SPATIAL_DOMAINS_EXPORT static const int kNqfaces = 6;
    SPATIAL_DOMAINS_EXPORT static const int kNtfaces = 0;
    SPATIAL_DOMAINS_EXPORT static const int kNfaces = kNqfaces + kNtfaces;
    SPATIAL_DOMAINS_EXPORT static const std::string XMLElementType;

protected:
    virtual void v_GenGeomFactors();
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

    static const unsigned int VertexEdgeConnectivity[8][3];
    static const unsigned int VertexFaceConnectivity[8][3];
    static const unsigned int EdgeFaceConnectivity[12][2];
};

typedef std::shared_ptr<HexGeom> HexGeomSharedPtr;
typedef std::map<int, HexGeomSharedPtr> HexGeomMap;
}
}

#endif
