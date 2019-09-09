////////////////////////////////////////////////////////////////////////////////
//
//  File: CFIMesh.h
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
//  Description: class which extracts an exisiting mesh from CFI source.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESHUTILS_CFIMESH_CFIMESH
#define NEKMESHUTILS_CFIMESH_CFIMESH

#include <NekMeshUtils/Module/Module.h>

#include <NekMeshUtils/CADSystem/CFI/CADSystemCFI.h>

namespace Nektar
{
namespace NekMeshUtils
{

/**
 * @brief class containing all surface meshing routines methods and classes
 */
class CFIMesh : public ProcessModule
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(MeshSharedPtr m)
    {
        return MemoryManager<CFIMesh>::AllocateSharedPtr(m);
    }
    static ModuleKey className;

    CFIMesh(MeshSharedPtr m);
    virtual ~CFIMesh();

    virtual void Process();

private:

    CADSystemCFISharedPtr m_cad;
    cfi::Model *m_model;
    std::map<std::string, int> m_nameToCurveId;
    std::map<std::string, int> m_nameToFaceId;
    std::map<std::string, int> m_nameToVertId;
};
}
}

#endif
