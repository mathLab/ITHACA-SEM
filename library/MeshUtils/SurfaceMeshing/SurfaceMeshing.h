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


#ifndef NEKTAR_MESHUTILS_SURFACEMESHING_SURFACEMESHING_H
#define NEKTAR_MESHUTILS_SURFACEMESHING_SURFACEMESHING_H

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

/**
 * @brief class containing all surface meshing routines methods and classes
 */
class SurfaceMeshing
{
    public:
        friend class MemoryManager<SurfaceMeshing>;

        /**
         * @brief Default constructor, requires the cad and octree objects to
         * begin
         */
        SurfaceMeshing(bool verbose,
                       const LibUtilities::CADSystemSharedPtr &cad,
                       const OctreeSharedPtr &octree,
                       const int &order) :
                           m_cad(cad),m_octree(octree),
                           m_order(order),m_verbose(verbose)
        {
        };

        /**
         * @brief Run all linear meshing routines
         */
        void Mesh();

        /**
         * @brief run all high-order surface meshing routines
         */
        void HOSurf();

        /**
         * @brief extract mesh objects
         */
        void Get(std::map<int, MeshNodeSharedPtr> &n,
                 std::map<int, MeshEdgeSharedPtr> &e,
                 std::map<int, MeshTriSharedPtr> &t)
        {
            t = Tris;
            e = Edges;
            n = Nodes;
        }

    private:

        /**
         * @brief get the gadient of the edge spring energy optimsation function
         */
        Array<OneD, NekDouble> EdgeGrad(NekDouble ux, NekDouble vx,
                                        std::vector<Array<OneD,NekDouble> > bcs,
                                        std::vector<NekDouble> weights,
                                        int surf, bool &valid);

        /**
         * @brief get the value of the edge spring energy
         */
        NekDouble EdgeF(NekDouble ux, NekDouble vx,
                        std::vector<Array<OneD,NekDouble> > bcs,
                        std::vector<NekDouble> weights,
                        int surf, bool &valid);

        /**
         * @brief get the value of face spring energy
         */
        NekDouble FaceF(NekDouble ux, NekDouble vx,
                        std::vector<Array<OneD,NekDouble> > bcs,
                        std::vector<NekDouble> weights,
                        int surf, bool &valid);

        /**
         * @brief get the gadient of the face spring energy optimsation function
         */
        Array<OneD, NekDouble> FaceGrad(NekDouble ux, NekDouble vx,
                                        std::vector<Array<OneD,NekDouble> > bcs,
                                        std::vector<NekDouble> weights,
                                        int surf, bool &valid);

        /**
         * @brief get the bracketing for brent optimisation
         */
        void Find1DBounds(NekDouble &a, NekDouble &b,
                          Array<OneD, NekDouble> uvi,
                          Array<OneD, NekDouble> df,
                          Array<OneD, NekDouble> bounds);

        /**
         * @brief perform brent optimsation
         */
        NekDouble BrentOpti(NekDouble ax, NekDouble bx,
                            NekDouble cx, NekDouble &fx,
                            NekDouble tol, int surf,
                            Array<OneD, NekDouble> uvi,
                            Array<OneD, NekDouble> df,
                            Array<OneD, NekDouble> bounds,
                            std::vector<Array<OneD,NekDouble> > bcs,
                            std::vector<NekDouble> weights,
                            NekDouble (SurfaceMeshing::*GetF)(
                            NekDouble, NekDouble,
                            std::vector<Array<OneD,NekDouble> >,
                            std::vector<NekDouble>, int, bool &));

        /**
         * @brief Validate the linear surface mesh
         */
        void Validate();

        /**
         * @brief Optimise the linear surface mesh using spring relaxation
         * and edge swapping
         */
        void Optimise();

        /// CAD object
        LibUtilities::CADSystemSharedPtr m_cad;
        /// Octree object
        OctreeSharedPtr m_octree;
        /// map of individual surface meshes from parametric surfaces
        std::map<int, SurfaceMeshSharedPtr> m_surfacemeshes;
        /// map of individual curve meshes of the curves in the domain
        std::map<int, CurveMeshSharedPtr> m_curvemeshes;
        /// order of the high-order mesh to be created
        int m_order;
        /// verbosity of the routines
        bool m_verbose;
        /// map of mesh nodes
        std::map<int, MeshNodeSharedPtr> Nodes;
        /// map of mesh edges
        std::map<int, MeshEdgeSharedPtr> Edges;
        /// map of mesh triangles
        std::map<int, MeshTriSharedPtr> Tris;
};

typedef boost::shared_ptr<SurfaceMeshing> SurfaceMeshingSharedPtr;

}
}

#endif
