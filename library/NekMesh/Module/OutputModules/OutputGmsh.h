////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputGmsh.h
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

#ifndef UTILITIES_NEKMESH_OUTPUTGMSH
#define UTILITIES_NEKMESH_OUTPUTGMSH

#include <NekMesh/Module/Module.h>
#include <LibUtilities/BasicUtils/HashUtils.hpp>

namespace Nektar
{
namespace NekMesh
{

struct ElmtConfigHash : std::unary_function<NekMesh::ElmtConfig, std::size_t>
{
    std::size_t operator()(NekMesh::ElmtConfig const& el) const
    {
        return hash_combine(
            (int)el.m_e, el.m_faceNodes, el.m_volumeNodes, el.m_order);
    }
};

/// Converter for Gmsh files.
class OutputGmsh : public NekMesh::OutputModule
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(NekMesh::MeshSharedPtr m) {
        return MemoryManager<OutputGmsh>::AllocateSharedPtr(m);
    }
    static NekMesh::ModuleKey className;

    OutputGmsh(NekMesh::MeshSharedPtr m);
    virtual ~OutputGmsh();

    /// Write mesh to output file.
    virtual void Process();
    
    virtual std::string GetModuleName()
    {
        return "OutputGmsh";
    }

private:
    std::unordered_map<NekMesh::ElmtConfig, unsigned int, ElmtConfigHash> elmMap;
};

}
}

#endif
