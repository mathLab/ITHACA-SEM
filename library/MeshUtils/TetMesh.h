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
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_MESHUTILS_TETMESH_TETMESH_H
#define NEKTAR_LIB_UTILITIES_MESHUTILS_TETMESH_TETMESH_H

#include <boost/shared_ptr.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <MeshUtils/Octree.h>
#include <MeshUtils/SurfaceMeshing/SurfaceMeshing.h>
#include <MeshUtils/MeshElem.hpp>

namespace Nektar{
namespace MeshUtils{

    class TetMesh
    {
    public:
        friend class MemoryManager<TetMesh>;

        TetMesh(const int id,
                const bool verb,
                const OctreeSharedPtr &oct,
                const SurfaceMeshingSharedPtr &sm)
                    : m_id(id), m_verbose(verb), m_octree(oct),
                      m_surfacemesh(sm)
        {
        };

        void Mesh();

        void Get(std::map<int, MeshNodeSharedPtr> &n,
                 std::map<int, MeshEdgeSharedPtr> &e,
                 std::map<int, MeshTriSharedPtr> &ti,
                 std::map<int, MeshTetSharedPtr> &te)
        {
            n = Nodes;
            e = Edges;
            ti = Tris;
            te = Tets;
        }

    private:

        bool Validate(std::map<int, MeshNodeSharedPtr> &Nodes);

        void AddNewPoint(Array<OneD, NekDouble> loc,
                                      std::map<int, MeshNodeSharedPtr> &Nodes,
                                      std::vector<int> tetnodes);

        int m_id;
        bool m_verbose;
        OctreeSharedPtr m_octree;
        SurfaceMeshingSharedPtr m_surfacemesh;

        std::vector<int> nodesintris;
        std::vector<int> m_stienerpoints;

        std::vector<int> nodesaddedinthisvalidate;

        int numtet;
        Array<OneD, Array<OneD, int> > tetconnect;

        std::map<int, MeshNodeSharedPtr> Nodes;
        std::map<int, MeshEdgeSharedPtr> Edges;
        std::map<int, MeshTriSharedPtr> Tris;
        std::map<int, MeshTetSharedPtr> Tets;
    };

    typedef boost::shared_ptr<TetMesh> TetMeshSharedPtr;

}
}

#endif
