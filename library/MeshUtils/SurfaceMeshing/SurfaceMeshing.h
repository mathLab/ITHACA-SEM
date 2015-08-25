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


#ifndef NEKTAR_LIB_UTILITIES_MESHUTILS_SURFACEMESHING_SURFACEMESHING_H
#define NEKTAR_LIB_UTILITIES_MESHUTILS_SURFACEMESHING_SURFACEMESHING_H

#include <boost/shared_ptr.hpp>

#include <LibUtilities/CADSystem/CADSystem.h>
#include <MeshUtils/Octree/Octree.h>
#include <MeshUtils/SurfaceMeshing/SurfaceMesh.h>
#include <MeshUtils/SurfaceMeshing/CurveMesh.h>
#include <MeshUtils/MeshElem.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

namespace Nektar {
namespace MeshUtils {

    class SurfaceMeshing
    {
    public:
        friend class MemoryManager<SurfaceMeshing>;

        SurfaceMeshing(bool verbose,
                       const LibUtilities::CADSystemSharedPtr &cad,
                       const OctreeSharedPtr &octree,
                       const int &order) :
                           m_cad(cad),m_octree(octree),
                           m_order(order),m_verbose(verbose)
        {
        };

        void Mesh();

        void HOSurf();

        void Get(std::map<int, MeshNodeSharedPtr> &n,
                 std::map<int, MeshEdgeSharedPtr> &e,
                 std::map<int, MeshTriSharedPtr> &t)
        {
            t = Tris;
            e = Edges;
            n = Nodes;
        }

        int GetNLinear(){return nodeinlinearmesh;}

    private:

        Array<OneD, NekDouble> EdgeGrad(Array<OneD, NekDouble> uv1,
                                        Array<OneD, NekDouble> uv2,
                                        Array<OneD, NekDouble> uvx, int surf, bool &valid);

        NekDouble EdgeF(NekDouble ux, NekDouble vx,
                        std::vector<Array<OneD,NekDouble> > bcs, int surf);

        NekDouble FaceF(NekDouble ux, NekDouble vx,
                        std::vector<Array<OneD,NekDouble> > bcs, int surf);

        Array<OneD, NekDouble> FaceGrad(Array<OneD, NekDouble> uvx,
                                        std::vector<Array<OneD, NekDouble> > bcs,
                                        int surf, bool &valid);

        void Find1DBounds(NekDouble &a, NekDouble &b,
                          Array<OneD, NekDouble> uvi,
                          Array<OneD, NekDouble> df,
                          Array<OneD, NekDouble> bounds);

        NekDouble BrentOpti(NekDouble ax, NekDouble bx,
                            NekDouble cx, NekDouble &fx,
                            NekDouble tol, int surf,
                            Array<OneD, NekDouble> uvi,
                            Array<OneD, NekDouble> df,
                            Array<OneD, NekDouble> bounds,
                            std::vector<Array<OneD,NekDouble> > bcs,
                            NekDouble (SurfaceMeshing::*GetF)(
                            NekDouble, NekDouble,
                            std::vector<Array<OneD,NekDouble> >, int));

        void Validate();
        void Optimise();
        //void EnergyEval(double *p, double *x, int m, int n, void *data);

        LibUtilities::CADSystemSharedPtr m_cad;
        OctreeSharedPtr m_octree;

        std::map<int, SurfaceMeshSharedPtr> m_surfacemeshes;
        std::map<int, CurveMeshSharedPtr> m_curvemeshes;

        int m_order;

        Array<OneD, NekDouble> W; //a way to pass the node weights to face evaluation

        bool m_verbose;

        int nodeinlinearmesh;

        std::map<int, MeshNodeSharedPtr> Nodes;
        std::map<int, MeshEdgeSharedPtr> Edges;
        std::map<int, MeshTriSharedPtr> Tris;


    };

    typedef boost::shared_ptr<SurfaceMeshing> SurfaceMeshingSharedPtr;
}
}

#endif
