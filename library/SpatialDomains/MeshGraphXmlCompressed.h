////////////////////////////////////////////////////////////////////////////////
//
//  File:  Geometry.h
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
//  Description:  This file contains the base class specification for the
//                Geometry class.
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_MGXMLCOM_H
#define NEKTAR_SPATIALDOMAINS_MGXMLCOM_H

#include "MeshGraphXml.h"

namespace Nektar
{
namespace SpatialDomains
{

class MeshGraphXmlCompressed : public MeshGraphXml
{
public:
    MeshGraphXmlCompressed()
    {
    }

    virtual ~MeshGraphXmlCompressed()
    {
    }

    static MeshGraphSharedPtr create()
    {
        return MemoryManager<MeshGraphXmlCompressed>::AllocateSharedPtr();
    }

    static std::string className;

private:
    void ReadVertices();
    void ReadCurves();

    void ReadEdges();
    void ReadFaces();

    void ReadElements1D();
    void ReadElements2D();
    void ReadElements3D();

    void WriteVertices(TiXmlElement *geomTag, PointGeomMap &verts);
    void WriteEdges(TiXmlElement *geomTag, SegGeomMap &edges);
    void WriteTris(TiXmlElement *faceTag, TriGeomMap &tris);
    void WriteQuads(TiXmlElement *faceTag, QuadGeomMap &quads);
    void WriteHexs(TiXmlElement *elmtTag, HexGeomMap &hexs);
    void WritePrisms(TiXmlElement *elmtTag, PrismGeomMap & pris);
    void WritePyrs(TiXmlElement *elmtTag, PyrGeomMap & pyrs);
    void WriteTets(TiXmlElement *elmtTag, TetGeomMap & tets);
    void WriteCurves(TiXmlElement *geomTag, CurveMap &edges, CurveMap &faces);
};

} // end of namespace
} // end of namespace

#endif
