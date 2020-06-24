////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputNekpp.h
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
//  Description: GMSH converter.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_NEKMESH_OUTPUTNEKPP
#define UTILITIES_NEKMESH_OUTPUTNEKPP

#include <NekMesh/Module/Module.h>

namespace Nektar
{
namespace NekMesh
{

/// Converter for Gmsh files.
class OutputNekpp : public NekMesh::OutputModule
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(NekMesh::MeshSharedPtr m)
    {
        return MemoryManager<OutputNekpp>::AllocateSharedPtr(m);
    }
    static NekMesh::ModuleKey className1, className2;

    OutputNekpp(NekMesh::MeshSharedPtr m);
    virtual ~OutputNekpp();

    /// Write mesh to output file.
    virtual void Process();

private:
    LibUtilities::Interpreter m_strEval;

    void TransferVertices(SpatialDomains::MeshGraphSharedPtr graph);
    void TransferEdges(
        SpatialDomains::MeshGraphSharedPtr graph,
        std::unordered_map<int, SpatialDomains::SegGeomSharedPtr> &edgeMap);
    void TransferFaces(
        SpatialDomains::MeshGraphSharedPtr graph,
        std::unordered_map<int, SpatialDomains::SegGeomSharedPtr> &edgeMap);
    void TransferElements(SpatialDomains::MeshGraphSharedPtr graph);
    void TransferCurves(SpatialDomains::MeshGraphSharedPtr graph);
    void TransferComposites(SpatialDomains::MeshGraphSharedPtr graph);
    void TransferDomain(SpatialDomains::MeshGraphSharedPtr graph);
};
}
}

#endif
