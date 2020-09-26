////////////////////////////////////////////////////////////////////////////////
//
//  File: SurfaceMeshing.h
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
//  Description: class containing all surfacemeshing routines and classes.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESHUTILS_SURFACEMESHING_SURFACEMESH
#define NEKMESHUTILS_SURFACEMESHING_SURFACEMESH

#include <NekMeshUtils/Module/Module.h>
#include <NekMeshUtils/SurfaceMeshing/FaceMesh.h>
#include <NekMeshUtils/SurfaceMeshing/CurveMesh.h>

namespace Nektar
{
namespace NekMeshUtils
{

/**
 * @brief class containing all surface meshing routines methods and classes
 */
class SurfaceMesh : public ProcessModule
{
public:

    /// Creates an instance of this class
    static std::shared_ptr<Module> create(MeshSharedPtr m)
    {
        return MemoryManager<SurfaceMesh>::AllocateSharedPtr(m);
    }
    static ModuleKey className;

    SurfaceMesh(MeshSharedPtr m);
    virtual ~SurfaceMesh();

    virtual void Process();

private:
    void Report();
    /// map of individual surface meshes from parametric surfaces
    std::map<int, FaceMeshSharedPtr> m_facemeshes;
    /// map of individual curve meshes of the curves in the domain
    std::map<int, CurveMeshSharedPtr> m_curvemeshes;
};

}
}

#endif
