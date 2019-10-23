////////////////////////////////////////////////////////////////////////////////
//
//  File: TetGeom.h
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
//  Description: Tetrahedral geometry information.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_TETGEOM
#define NEKTAR_SPATIALDOMAINS_TETGEOM

#include <SpatialDomains/SpatialDomainsDeclspec.h>
#include <SpatialDomains/Geometry3D.h>
#include <SpatialDomains/TriGeom.h>

namespace Nektar
{
namespace SpatialDomains
{

class TetGeom : public Geometry3D
{
public:
    SPATIAL_DOMAINS_EXPORT TetGeom();
    SPATIAL_DOMAINS_EXPORT TetGeom(int id, const TriGeomSharedPtr faces[]);
    SPATIAL_DOMAINS_EXPORT ~TetGeom();

    SPATIAL_DOMAINS_EXPORT static const int kNverts = 4;
    SPATIAL_DOMAINS_EXPORT static const int kNedges = 6;
    SPATIAL_DOMAINS_EXPORT static const int kNqfaces = 0;
    SPATIAL_DOMAINS_EXPORT static const int kNtfaces = 4;
    SPATIAL_DOMAINS_EXPORT static const int kNfaces = kNqfaces + kNtfaces;
    SPATIAL_DOMAINS_EXPORT static const std::string XMLElementType;

protected:
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
    virtual void v_GenGeomFactors();

private:
    void SetUpLocalEdges();
    void SetUpLocalVertices();
    void SetUpEdgeOrientation();
    void SetUpFaceOrientation();
    void SetUpXmap();

    static const unsigned int VertexEdgeConnectivity[4][3];
    static const unsigned int VertexFaceConnectivity[4][3];
    static const unsigned int EdgeFaceConnectivity[6][2];
};

typedef std::shared_ptr<TetGeom> TetGeomSharedPtr;
typedef std::vector<TetGeomSharedPtr> TetGeomVector;
typedef std::map<int, TetGeomSharedPtr> TetGeomMap;


}
}

#endif // NEKTAR_SPATIALDOMAINS_TETGEOM
