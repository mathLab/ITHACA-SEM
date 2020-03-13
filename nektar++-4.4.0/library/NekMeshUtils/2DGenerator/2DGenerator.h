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
//  License for the specific language governing rights and limitations under
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

#ifndef NEKMESHUTILS_2D_2D
#define NEKMESHUTILS_2D_2D

#include <NekMeshUtils/Module/Module.h>
#include <NekMeshUtils/SurfaceMeshing/CurveMesh.h>
#include <NekMeshUtils/SurfaceMeshing/FaceMesh.h>

#include <LibUtilities/Interpreter/AnalyticExpressionEvaluator.hpp>

namespace Nektar
{
namespace NekMeshUtils
{

/**
 * @brief class containing all surface meshing routines methods and classes
 */
class Generator2D : public ProcessModule
{
public:
    /// Creates an instance of this class
    static boost::shared_ptr<Module> create(MeshSharedPtr m)
    {
        return MemoryManager<Generator2D>::AllocateSharedPtr(m);
    }
    static ModuleKey className;

    Generator2D(MeshSharedPtr m);
    virtual ~Generator2D();

    virtual void Process();

private:

    void MakeBLPrep();

    void MakeBL(int faceid);

    void Report();
    /// map of individual surface meshes from parametric surfaces
    std::map<int, FaceMeshSharedPtr> m_facemeshes;
    /// map of individual curve meshes of the curves in the domain
    std::map<int, CurveMeshSharedPtr> m_curvemeshes;

    std::vector<unsigned int> m_blCurves;
    LibUtilities::AnalyticExpressionEvaluator m_thickness;
    int m_thickness_ID;
    std::map<NodeSharedPtr, std::vector<EdgeSharedPtr> > m_nodesToEdge;
};
}
}

#endif
