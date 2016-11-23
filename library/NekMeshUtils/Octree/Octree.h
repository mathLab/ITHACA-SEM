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

#ifndef NEKTAR_MESHUTILS_OCTREE_OCTREE
#define NEKTAR_MESHUTILS_OCTREE_OCTREE

#include "SourcePoint.hpp"
#include "Octant.h"
#include <NekMeshUtils/MeshElements/Mesh.h>

namespace Nektar
{
namespace NekMeshUtils
{

/**
 * @brief class for octree
 *
 * This class contains the routines to generate and query a automatically
 * generated set of mesh spacing parameters based on the CAD
 */
class Octree
{
public:

    Octree(MeshSharedPtr m) : m_mesh(m)
    {
        m_udsfileset = false;
    }

    Octree()
    {
    }

    /**
     * @brief builds the octree based on curvature sampling and user defined
     * spacing
     */
    void Process();

    /**
     * @brief once constructed queryies the octree based on x,y,z location
     * to get a mesh spacing
     *
     * @param loc array of x,y,z
     * @return mesh spacing parameter
     */
    NekDouble Query(Array<OneD, NekDouble> loc);

    /**
     * @brief returns the miminum spacing in the octree (for meshing purposes)
     *
     * @return miminum delta in octree
     */
    NekDouble GetMinDelta()
    {
        return m_minDelta;
    }

    /**
     * @brief sets the parameters used for curvature sampling
     *
     * @param min minimum spacing to be found in the mesh
     * @param max maximum spacing to be found in the mesh
     * @param eps curvature sensivity relating radius of curvature to spacing
     */
    void SetParameters(NekDouble &min, NekDouble &max, NekDouble &ep)
    {
        m_minDelta = min;
        m_maxDelta = max;
        m_eps = ep;
    }

    /**
     * @brief populates the mesh m with a invalid hexahedral mesh based on the
     *        octree, used for visualisation
     * @param nm name of the mesh file to be made
     */
    void WriteOctree(std::string nm);

    /**
     * @brief informs the octree there is a user defined spacing file
     *
     * @param nm name of the user defined spacing file
     */
    void UDS(std::string nm)
    {
        m_udsfile = nm;
        m_udsfileset = true;
    }

private:

    /**
     * @brief Smooths specification over all octants to a gradation criteria
     */
    void SmoothAllOctants();

    /**
     * @brief gets an optimum number of curvature sampling points and
     * calculates the curavture at these points
     */
    void CompileSourcePointList();

    /**
     * @brief Function which initiates and controls the subdivision process
     */
    void SubDivide();

    /**
     * @brief Smooths specification over the surface encompasing octants to a
     *        gradation criteria
     */
    void SmoothSurfaceOctants();

    /**
     * @brief takes the mesh specification from surface octants and
     *        progates that through the domain so all octants have a
     * specification
     *        using gradiation crieteria
     */
    void PropagateDomain();

    /**
     * @brief estimates the number of elements to be created in the mesh
     */
    int CountElemt();

    /**
     * @brief Calculates the difference in delta divided by the difference
     *        in location between two octants i and j
     */
    NekDouble ddx(OctantSharedPtr i, OctantSharedPtr j);

    /**
     * @brief Looks over all leaf octants and checks that their neigbour
     *        assigments are valid
     */
    bool VerifyNeigbours();

    /// minimum delta in the octree
    NekDouble m_minDelta;
    /// maximum delta in the octree
    NekDouble m_maxDelta;
    /// curavture sensivity paramter
    NekDouble m_eps;
    /// x,y,z location of the center of the octree
    Array<OneD, NekDouble> m_centroid;
    /// physical size of the octree
    NekDouble m_dim;
    /// list of source points
    std::vector<SPBaseSharedPtr> m_SPList;
    /// list of leaf octants
    std::vector<OctantSharedPtr> m_octants;
    /// master octant for searching
    OctantSharedPtr m_masteroct;
    /// number of octants made, used for id index
    int m_numoct;
    /// Mesh object
    MeshSharedPtr m_mesh;
    /// user defined spacing has been set
    bool m_udsfileset;
    /// name of the user defined spacing file
    std::string m_udsfile;
};
typedef boost::shared_ptr<Octree> OctreeSharedPtr;

}
}

#endif
