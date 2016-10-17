////////////////////////////////////////////////////////////////////////////////
//
//  File: CADSurf.h
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
//  Description: CAD object surface.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NekMeshUtils_CADSYSTEM_CADSURF
#define NekMeshUtils_CADSYSTEM_CADSURF

#include <NekMeshUtils/CADSystem/CADObj.h>
#include <NekMeshUtils/CADSystem/CADVert.h>
#include <NekMeshUtils/CADSystem/CADSystem.h>

namespace Nektar
{
namespace NekMeshUtils
{

class CADCurve;
typedef boost::shared_ptr<CADCurve> CADCurveSharedPtr;

/**
 * @brief base class for a cad surface
 */

class CADSurf : public CADObj
{
public:
    friend class MemoryManager<CADSurf>;

    /**
     * @brief Default constructor.
     */
    CADSurf()
    {
    }

    ~CADSurf()
    {
    }

    /**
     * @brief Get the loop structures which bound the cad surface
     */
    std::vector<EdgeLoop> GetEdges()
    {
        return m_edges;
    }

    /**
     * @brief Get the limits of the parametric space for the surface.
     *
     * @return Array of 4 entries with parametric umin,umax,vmin,vmax.
     */
    virtual Array<OneD, NekDouble> GetBounds() = 0;

    /**
     * @brief Get the normal vector at parametric point u,v.
     *
     * @param uv Array of u and v parametric coords.
     * @return Array of xyz components of normal vector.
     */
    virtual Array<OneD, NekDouble> N    (Array<OneD, NekDouble> uv) = 0;

    /**
     * @brief Get the set of first derivatives at parametric point u,v
     *
     * @param uv Array of u and v parametric coords.
     * @return Array of xyz copmonents of first derivatives.
     */
    virtual Array<OneD, NekDouble> D1   (Array<OneD, NekDouble> uv) = 0;

    /**
     * @brief Get the set of second derivatives at parametric point u,v
     *
     * @param uv array of u and v parametric coords
     * @return array of xyz copmonents of second derivatives
     */
    virtual Array<OneD, NekDouble> D2   (Array<OneD, NekDouble> uv) = 0;

    /**
     * @brief Get the x,y,z at parametric point u,v.
     *
     * @param uv Array of u and v parametric coords.
     * @return Array of xyz location.
     */
    virtual Array<OneD, NekDouble> P    (Array<OneD, NekDouble> uv) = 0;

    /**
     * @brief Performs a reverse look up to find u,v and x,y,z.
     *
     * @param p Array of xyz location
     * @return The parametric location of xyz on this surface
     */
    virtual Array<OneD, NekDouble> locuv(Array<OneD, NekDouble> p) = 0;

    /**
     * @brief does unconstrained locuv to project point from anywhere
     * and calculate the distance between the orthonormal projection to the surface
     * and the point
     */
    virtual NekDouble DistanceTo(Array<OneD, NekDouble> p) = 0;

    /**
     * @brief takes a point from anywhere find the nearest surface point and its
     * uv
     */
    virtual void ProjectTo(Array<OneD, NekDouble> &tp, Array<OneD, NekDouble> &uv) = 0;

    /**
     * @brief returns curvature at point uv
     */
    virtual NekDouble Curvature(Array<OneD, NekDouble> uv) = 0;

    /**
     * @brief sets the flag to reverse the normal for this suface,
     * this is determined in cadsystem and ensures all surface normals,
     * point internaly
     */
    void SetReverseNomral()
    {
        m_correctNormal = false;
    }

    /**
     * @brief query reversed normal
     */
    bool IsReversedNormal()
    {
        return !m_correctNormal;
    }

protected:
    /// normal
    bool m_correctNormal;
    /// List of bounding edges in loops with orientation.
    std::vector<EdgeLoop> m_edges;

    /// Function which tests the the value of uv used is within the surface
    virtual void Test(Array<OneD, NekDouble> uv) = 0;
};

typedef boost::shared_ptr<CADSurf> CADSurfSharedPtr;

typedef LibUtilities::NekFactory<std::string, CADSurf> CADSurfFactory;

CADSurfFactory& GetCADSurfFactory();

}
}

#endif
