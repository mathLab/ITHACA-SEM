////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessLinear.h
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
//  Description: Removes high-order information from the mesh.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_NEKMESH_PROCESSLINEAR
#define UTILITIES_NEKMESH_PROCESSLINEAR

#include <NekMeshUtils/Module/Module.h>

namespace Nektar
{
namespace Utilities
{
/**
 * @brief This processing module removes all the high-order information
 * from the mesh leaving just the linear elements
 */
class ProcessLinear : public NekMeshUtils::ProcessModule
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(NekMeshUtils::MeshSharedPtr m)
    {
        return MemoryManager<ProcessLinear>::AllocateSharedPtr(m);
    }
    static NekMeshUtils::ModuleKey className;

    ProcessLinear(NekMeshUtils::MeshSharedPtr m);
    virtual ~ProcessLinear();

    /// Write mesh to output file.
    virtual void Process();
private:
    bool Invalid(NekMeshUtils::ElementSharedPtr el, NekDouble thr);
};
}
}

#endif
