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


#ifndef NEKTAR_LIB_UTILITIES_MESHUTILS_SURFACEMESH_SURFACEMESH_H
#define NEKTAR_LIB_UTILITIES_MESHUTILS_SURFACEMESH_SURFACEMESH_H

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

        SurfaceMesh(const int id,
                    const bool verb,
                    const LibUtilities::CADSurfSharedPtr &cad,
                    const OctreeSharedPtr &oct,
                    const std::map<int, CurveMeshSharedPtr> &cmeshes,
                    const int &order)
                        : m_verbose(verb), m_cadsurf(cad), m_octree(oct),
                          m_curvemeshes(cmeshes),m_order(order),m_id(id)

        {
            m_edges = m_cadsurf->GetEdges();
        };

        void Mesh(std::map<int, MeshNodeSharedPtr> &Nodes,
                  std::map<int, MeshEdgeSharedPtr> &Edges,
                  std::map<int, MeshTriSharedPtr> &Tris);

        void Report();


    private:

        void Stretching();

        bool Validate(std::map<int, MeshNodeSharedPtr> &Nodes);

        void OrientateCurves(std::map<int, MeshNodeSharedPtr> &Nodes);

        void AddNewPoint(Array<OneD, NekDouble> uv,
                         std::map<int, MeshNodeSharedPtr> &Nodes);

        bool m_verbose;
        LibUtilities::CADSurfSharedPtr m_cadsurf;
        OctreeSharedPtr m_octree;
        std::map<int, CurveMeshSharedPtr> m_curvemeshes;

        std::vector<std::vector<std::pair<int,int> > > m_edges;

        int m_order;

        int m_id;

        std::vector<std::vector<int> > orderedLoops;
        std::vector<std::vector<NekDouble> > m_centers;

        std::vector<int> m_stienerpoints;

        int numtri, nump;
        Array<OneD, Array<OneD, int> > Connec;

        NekDouble pasr,asr;

    };

    typedef boost::shared_ptr<SurfaceMesh> SurfaceMeshSharedPtr;
}
}

#endif
