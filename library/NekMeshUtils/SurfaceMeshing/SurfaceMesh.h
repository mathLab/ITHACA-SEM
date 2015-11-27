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


#ifndef NekMeshUtils_SURFACEMESHING_SURFACEMESH
#define NekMeshUtils_SURFACEMESHING_SURFACEMESH

#include <NekMeshUtils/MeshElements/MeshElements.h>
#include <NekMeshUtils/CADSystem/CADSystem.h>
#include <NekMeshUtils/Octree/Octree.h>
#include <NekMeshUtils/SurfaceMeshing/FaceMesh.h>
#include <NekMeshUtils/SurfaceMeshing/CurveMesh.h>


namespace Nektar
{
namespace NekMeshUtils
{

/**
 * @brief class containing all surface meshing routines methods and classes
 */
class SurfaceMesh
{
    public:
        friend class MemoryManager<SurfaceMesh>;

        /**
         * @brief Default constructor, requires the cad and octree objects to
         * begin
         */
        SurfaceMesh(   MeshSharedPtr m,
                       CADSystemSharedPtr cad,
                       OctreeSharedPtr octree,
                       std::vector<unsigned int> sy, NekDouble b) :
                          m_mesh(m), m_cad(cad), m_octree(octree), m_symsurfs(sy), m_bl(b)
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
         * @brief Validate the linear surface mesh
         */
        void Validate();

        /**
         * @brief Print brief information to screen
         */
        void Report();


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
                            NekDouble (SurfaceMesh::*GetF)(
                            NekDouble, NekDouble,
                            std::vector<Array<OneD,NekDouble> >,
                            std::vector<NekDouble>, int, bool &));

        /// mesh object
        MeshSharedPtr m_mesh;
        /// CAD object
        CADSystemSharedPtr m_cad;
        /// Octree object
        OctreeSharedPtr m_octree;
        /// map of individual surface meshes from parametric surfaces
        std::map<int, FaceMeshSharedPtr> m_facemeshes;
        /// map of individual curve meshes of the curves in the domain
        std::map<int, CurveMeshSharedPtr> m_curvemeshes;

        std::vector<unsigned int> m_symsurfs;
        NekDouble m_bl;
};

typedef boost::shared_ptr<SurfaceMesh> SurfaceMeshSharedPtr;

}
}

#endif
