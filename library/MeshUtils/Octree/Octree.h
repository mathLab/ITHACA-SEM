////////////////////////////////////////////////////////////////////////////////
//
//  File: Octree.h
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
//  Description: octree object header
//
////////////////////////////////////////////////////////////////////////////////


#ifndef NEKTAR_MESHUTILS_OCTREE_OCTREE_H
#define NEKTAR_MESHUTILS_OCTREE_OCTREE_H

#include <boost/shared_ptr.hpp>

#include <LibUtilities/CADSystem/CADSystem.h>
#include <MeshUtils/Octree/CurvaturePoint.hpp>
#include <MeshUtils/Octree/Octant.h>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

namespace Nektar {
namespace MeshUtils {

/**
 * @brief class for octree
 *
 * This class contains the routines to generate and query a automatically
 * generated set of mesh spacing parameters based on the CAD
 */
class Octree
{
    public:
        friend class MemoryManager<Octree>;

        /**
         * @brief Defualt constructor
         *
         * @param cad CAD object
         * @param ver bool verbose
         */
        Octree(const LibUtilities::CADSystemSharedPtr &cad,
               const bool ver) : m_cad(cad), m_verbose(ver)
        {
        };

        /**
         * @brief executes octree building routines
         *
         * @param min minimum delta to be found in the octree
         * @param max maximum delta to be found in the Octree
         * @param eps curvature sensivity parameter
         */
        void Build(const NekDouble min, const NekDouble max,
                   const NekDouble eps);

        /**
         * @brief once constructed queryies the octree based on x,y,z location
         * to get a mesh spacing
         *
         * @param loc array of x,y,z
         * @return mesh spacing parameter
         * @todo improve this algorithm for both robustness and completness,
         * functions just fine regardless
         */
        NekDouble Query(Array<OneD, NekDouble> loc);

        /**
         * @brief returns the miminum spacing in the octree (for meshing purposes)
         *
         * @return miminum delta in octree
         */
        NekDouble GetMinDelta(){return m_minDelta;}

    private:

        /**
         * @brief Smooths specification over all octants to a
         * gradation criteria
         */
        void SmoothAllOctants();
        
        /**
         * @brief gets an optimum number of curvature sampling points and
         * calculates the curavture at these points
         */
        void CompileCuravturePointList();

        /**
         * @brief Recursive alorithm which divides and creates new octants based
         * on the geometry
         */
        void InitialSubDivide(int parent);

        /**
         * @brief Recursive alorithm which subdivides octants so that the
         * neighbours differ by no more that one level
         */
        void SubDivideByLevel();

        /**
         * @brief Subdivision step for smoothoctants()
         */
        void SubDivideLevel(int parent);

        /**
         * @brief Smooths specification over the surface octants to a
         * gradation criteria
         */
        void SmoothSurfaceOctants();

        /**
         * @brief takes the mesh specification from surface octants and
         * progates that through the domain so all octants have a specification
         * using gradiation crieteria
         */
        void PropagateDomain();

        /**
         * @brief estimates the number of elements to be creted in the mesh
         */
        int CountElemt();

        /**
         * @brief Calculates the difference in delta divided by the difference
         * in location between two octants i and j
         */
        NekDouble ddx(int i, int j);

        /// minimum delta in the octree
        NekDouble m_minDelta;
        /// maximum delta in the octree
        NekDouble m_maxDelta;
        /// curavture sensivity paramter
        NekDouble m_eps;

        /// cad object
        LibUtilities::CADSystemSharedPtr m_cad;
        /// verbose output?
        bool m_verbose;
        /// max and min dimensions of the domain 6 varibles
        Array<OneD, NekDouble> BoundingBox;
        /// list of curvature sample points
        std::vector<CurvaturePointSharedPtr> m_cpList;
        /// list of octants
        std::vector<OctantSharedPtr> OctantList;
        /// number which do not need subdivding, when this is 0 octree complete
        int m_totNotDividing;

};

typedef boost::shared_ptr<Octree> OctreeSharedPtr;

}
}

#endif
