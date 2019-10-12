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

#ifndef NEKMESHUTILS_2D_2D
#define NEKMESHUTILS_2D_2D

#include <NekMeshUtils/Module/Module.h>
#include <NekMeshUtils/SurfaceMeshing/CurveMesh.h>
#include <NekMeshUtils/SurfaceMeshing/FaceMesh.h>

#include <LibUtilities/Interpreter/Interpreter.h>

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
    static std::shared_ptr<Module> create(MeshSharedPtr m)
    {
        return MemoryManager<Generator2D>::AllocateSharedPtr(m);
    }
    static ModuleKey className;

    Generator2D(MeshSharedPtr m);

    virtual ~Generator2D();

    virtual void Process();

private:
    void FindBLEnds();

    void MakeBLPrep();

    void PeriodicPrep();

    void MakePeriodic();

    void MakeBL(int faceid);

    void Report();
    /// map of individual surface meshes from parametric surfaces
    std::map<int, FaceMeshSharedPtr> m_facemeshes;
    /// map of individual curve meshes of the curves in the domain
    std::map<int, CurveMeshSharedPtr> m_curvemeshes;
    /// map of periodic curve pairs
    std::map<unsigned, unsigned> m_periodicPairs;
    /// list of curves needing a BL
    std::vector<unsigned int> m_blCurves;
    /// map of curves and Bl ends: 0, 1 or 2 (for both)
    std::map<unsigned, unsigned> m_blends;
    /// list of BL edges
    std::vector<EdgeSharedPtr> m_blEdges;
    /// BL thickness expression
    LibUtilities::Interpreter m_thickness;
    /// BL thickness expression ID
    int m_thickness_ID;
    /// map of BL curve nodes to adjacent edges
    std::map<NodeSharedPtr, std::vector<EdgeSharedPtr>> m_nodesToEdge;
};
}
}

#endif
