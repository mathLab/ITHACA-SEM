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


#ifndef NEKTAR_MESHUTILS_SURFACEMESHING_SURFACEMESH_H
#define NEKTAR_MESHUTILS_SURFACEMESHING_SURFACEMESH_H

#include <boost/shared_ptr.hpp>

#include <MeshUtils/MeshElem.hpp>
#include <LibUtilities/CADSystem/CADSurf.h>
#include <MeshUtils/Octree/Octree.h>
#include <MeshUtils/SurfaceMeshing/CurveMesh.h>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>


namespace Nektar {
namespace MeshUtils {

class SurfaceMesh
{
    public:
        friend class MemoryManager<SurfaceMesh>;

        /**
         * @brief Default constructor
         */
        SurfaceMesh(const int id,
                    const bool verb,
                    const LibUtilities::CADSurfSharedPtr &cad,
                    const OctreeSharedPtr &oct,
                    const std::map<int, CurveMeshSharedPtr> &cmeshes)
                        : m_verbose(verb), m_cadsurf(cad), m_octree(oct),
                          m_curvemeshes(cmeshes),m_id(id)

        {
            m_edges = m_cadsurf->GetEdges();
        };

        /**
         * @brief mesh exectuation command
         */
        void Mesh(std::map<int, MeshNodeSharedPtr> &Nodes,
                  std::map<int, MeshEdgeSharedPtr> &Edges,
                  std::map<int, MeshTriSharedPtr> &Tris);

        /**
         * @brief Print report of the surface mesh to screen
         */
        void Report();


    private:

        /**
         * @brief Calculate the paramter plane streching factor
         */
        void Stretching();

        /**
         * @brief Validate the surface mesh base on the octree and real
         * dimensions of the edges
         */
        bool Validate(std::map<int, MeshNodeSharedPtr> &Nodes);

        /**
         * @brief Get the boundries of the surface and extracts the nodes from
         * the curve meshes in the correct order
         */
        void OrientateCurves(std::map<int, MeshNodeSharedPtr> &Nodes);

        /**
         * @brief addes a new stiener point to the triangulation for meshing
         */
        void AddNewPoint(Array<OneD, NekDouble> uv,
                         std::map<int, MeshNodeSharedPtr> &Nodes);

        ///verbosity
        bool m_verbose;
        /// CAD surface
        LibUtilities::CADSurfSharedPtr m_cadsurf;
        /// Octree object
        OctreeSharedPtr m_octree;
        /// Map of the curve meshes which bound the surfaces
        std::map<int, CurveMeshSharedPtr> m_curvemeshes;
        /// data structure containing the edges, their order and oreientation for the surface
        std::vector<std::vector<std::pair<int,int> > > m_edges;
        /// id of the surface mesh
        int m_id;
        /// list of boundary nodes in their order loops
        std::vector<std::vector<int> > orderedLoops;
        /// center coords of the loops
        std::vector<std::vector<NekDouble> > m_centers;
        /// list of stiener points in the triangulation
        std::vector<int> m_stienerpoints;
        /// number of triangles and points
        int numtri, nump;
        /// list of node connectivities forming triangles
        Array<OneD, Array<OneD, int> > Connec;
        /// paramter plane and real space aspect ratios
        NekDouble pasr,asr;
};

typedef boost::shared_ptr<SurfaceMesh> SurfaceMeshSharedPtr;

}
}

#endif
