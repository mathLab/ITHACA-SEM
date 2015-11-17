////////////////////////////////////////////////////////////////////////////////
//
//  File: SurfaceMesh.h
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
//  Description: class for indivdual surface meshes
//
////////////////////////////////////////////////////////////////////////////////


#ifndef MESHUTILS_SURFACEMESHING_FACEMESH
#define MESHUTILS_SURFACEMESHING_FACEMESH

#include <MeshUtils/MeshElements/MeshElements.h>
#include <MeshUtils/CADSystem/CADSurf.h>
#include <MeshUtils/Octree/Octree.h>
#include <MeshUtils/SurfaceMeshing/CurveMesh.h>

namespace Nektar
{
namespace MeshUtils
{

class FaceMesh
{
public:
    friend class MemoryManager<FaceMesh>;

    /**
     * @brief Default constructor
     */
    FaceMesh(   const int id,
                MeshSharedPtr m,
                CADSurfSharedPtr cad,
                OctreeSharedPtr oct,
                const std::map<int, CurveMeshSharedPtr> &cmeshes)
                    : m_mesh(m), m_cadsurf(cad), m_octree(oct),
                      m_curvemeshes(cmeshes),m_id(id)

    {
        m_edgeloops = m_cadsurf->GetEdges();
    };

    /**
     * @brief mesh exectuation command
     */
    void Mesh();

    void Smoothing();

    void DiagonalSwap();


private:

    /**
     * @brief Calculate the paramter plane streching factor
     */
    void Stretching();

    void BuildLocalMesh();

    void OptimiseLocalMesh();

    /**
     * @brief Validate the surface mesh base on the octree and real
     * dimensions of the edges
     */
    bool Validate();

    /**
     * @brief Get the boundries of the surface and extracts the nodes from
     * the curve meshes in the correct order
     */
    void OrientateCurves();

    /**
     * @brief addes a new stiener point to the triangulation for meshing
     */
    void AddNewPoint(Array<OneD, NekDouble> uv);

    ///mesh pointer
    MeshSharedPtr m_mesh;
    /// CAD surface
    CADSurfSharedPtr m_cadsurf;
    /// Octree object
    OctreeSharedPtr m_octree;
    /// Map of the curve meshes which bound the surfaces
    std::map<int, CurveMeshSharedPtr> m_curvemeshes;
    /// data structure containing the edges, their order and oreientation for the surface
    std::vector<EdgeLoop> m_edgeloops;
    /// id of the surface mesh
    int m_id;
    /// list of boundary nodes in their order loops
    std::vector<std::vector<NodeSharedPtr> > orderedLoops;
    /// list of stiener points in the triangulation
    std::vector<NodeSharedPtr> m_stienerpoints;
    /// paramter plane and real space aspect ratios
    NekDouble pasr,asr;
    /// triangle connectiviities
    std::vector<std::vector<NodeSharedPtr> > m_connec;

    NodeSet m_localNodes;
    EdgeSet m_localEdges;
    std::vector<ElementSharedPtr> m_localElements;
};

typedef boost::shared_ptr<FaceMesh> FaceMeshSharedPtr;

}
}

#endif
