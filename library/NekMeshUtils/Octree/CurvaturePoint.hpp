////////////////////////////////////////////////////////////////////////////////
//
//  File: Curavturepoint.hpp
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
//  Description: class and methods of curvature sampling point
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NekMeshUtils_OCTREE_CURAVTUREPOINT
#define NekMeshUtils_OCTREE_CURAVTUREPOINT

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

namespace Nektar
{
namespace NekMeshUtils
{

/**
 * @brief class for a curvature samlping Point
 */
class CurvaturePoint
{
public:
    friend class MemoryManager<CurvaturePoint>;

    /**
     * @brief constructor for a valid point (has radius of curvature)
     */
    CurvaturePoint(int i,
                   Array<OneD, NekDouble> uv,
                   Array<OneD, NekDouble> l,
                   NekDouble d,
                   bool bnd = true)
        : sid(i), m_uv(uv), m_loc(l), m_delta(d), m_boundary(bnd)
    {
        m_valid = true;
    }

    /**
     * @brief constructor for a invalid point
     */
    CurvaturePoint(int i, Array<OneD, NekDouble> uv, Array<OneD, NekDouble> l)
        : sid(i), m_uv(uv), m_loc(l)
    {
        m_delta    = -1;
        m_valid    = false;
        m_boundary = true;
    }

    /**
     * @brief return bool on whether point is valid
     */
    bool IsValid()
    {
        return m_valid;
    }

    bool Isboundary()
    {
        return m_boundary;
    }

    /**
     * @brief get mesh spacing paramter
     */
    NekDouble GetDelta()
    {
        if (m_valid)
        {
            return m_delta;
        }
        else
        {
            return -1;
        }
    }

    /**
     * @brief get location of point
     */
    Array<OneD, NekDouble> GetLoc()
    {
        return m_loc;
    }

    /**
     * @brief gets the corresponding cad information for the point
     */
    void GetCAD(int &surf, Array<OneD, NekDouble> &uv)
    {
        surf = sid;
        uv   = m_uv;
    }

    void SetDelta(NekDouble i)
    {
        m_delta = i;
    }

private:
    /// surf id
    int sid;
    /// uv coord on surf
    Array<OneD, NekDouble> m_uv;
    /// x,y,z location
    Array<OneD, NekDouble> m_loc;
    /// normal vector of surface at point
    NekDouble m_delta;
    /// valid point or not
    bool m_valid;

    bool m_boundary;
};

typedef boost::shared_ptr<CurvaturePoint> CurvaturePointSharedPtr;
}
}

#endif
