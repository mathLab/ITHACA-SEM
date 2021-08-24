////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessCurvedEdges.h
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
//  Description: Abstract base class for creating curved edges on boundaries.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_NEKMESH_PROCESSCURVEDEDGES
#define UTILITIES_NEKMESH_PROCESSCURVEDEDGES

#include <NekMesh/Module/Module.h>

namespace Nektar
{
namespace NekMesh
{

class ProcessCurvedEdges : public NekMesh::ProcessModule
{
public:
    ProcessCurvedEdges(NekMesh::MeshSharedPtr m);
    virtual ~ProcessCurvedEdges();

    /// Write mesh to output file.
    virtual void Process();
    void GenerateEdgeNodes(NekMesh::EdgeSharedPtr edge)
    {
        v_GenerateEdgeNodes( edge);
    }

    virtual std::string GetModuleName()
    {
        return "ProcessCurvedEdges";
    }
protected:
    virtual void v_GenerateEdgeNodes(NekMesh::EdgeSharedPtr edge) = 0;
};
}
}

#endif
