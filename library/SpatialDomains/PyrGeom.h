////////////////////////////////////////////////////////////////////////////////
//
//  File: PyrGeom.h
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
//  Description: Pyramidic geometry information.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_PYRGEOM_H
#define NEKTAR_SPATIALDOMAINS_PYRGEOM_H

#include <SpatialDomains/Geometry3D.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>

namespace Nektar
{
namespace SpatialDomains
{

class PyrGeom : public Geometry3D
{
public:
    SPATIAL_DOMAINS_EXPORT PyrGeom();
    SPATIAL_DOMAINS_EXPORT PyrGeom(int id, const Geometry2DSharedPtr faces[]);
    SPATIAL_DOMAINS_EXPORT ~PyrGeom();

    SPATIAL_DOMAINS_EXPORT static const int kNverts = 5;
    SPATIAL_DOMAINS_EXPORT static const int kNedges = 8;
    SPATIAL_DOMAINS_EXPORT static const int kNqfaces = 1;
    SPATIAL_DOMAINS_EXPORT static const int kNtfaces = 4;
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
    virtual int v_GetDir(const int faceidx, const int facedir) const;
    virtual void v_Reset(CurveMap &curvedEdges, CurveMap &curvedFaces);
    virtual void v_Setup();

private:
    void SetUpLocalEdges();
    void SetUpLocalVertices();
    void SetUpEdgeOrientation();
    void SetUpFaceOrientation();
    void SetUpXmap();
};

typedef std::shared_ptr<PyrGeom> PyrGeomSharedPtr;
typedef std::map<int, PyrGeomSharedPtr> PyrGeomMap;
}
}


#endif // NEKTAR_SPATIALDOMAINS_PYRGEOM_H
