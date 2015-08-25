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

#ifndef NEKTAR_MESHUTILS_OCTREE_CURAVTUREPOINT_H
#define NEKTAR_MESHUTILS_OCTREE_CURAVTUREPOINT_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

namespace Nektar {
namespace MeshUtils {

/**
 * @brief class for a curvature samlping Point
 */
class CurvaturePoint {

    public:
        friend class MemoryManager<CurvaturePoint>;

        /**
         * @brief constructor for a valid point (has radius of curvature)
         */
        CurvaturePoint(Array<OneD, NekDouble> l, NekDouble R,
                       Array<OneD, NekDouble> n) : m_loc(l), m_norm(n),
                                                   m_radius(R)
        {
            m_valid = true;
        }

        /**
         * @brief constructor for a invalid point
         */
        CurvaturePoint(Array<OneD, NekDouble> l,
                       Array<OneD, NekDouble> n) : m_loc(l), m_norm(n)
        {
            m_delta = -1;
            m_valid = false;
        }

        /**
         * @brief calculate mesh spacing delta from radius of curvature
         */
        void Process(NekDouble min, NekDouble max, NekDouble eps)
        {
            if(m_valid)
            {
                m_delta = 2.0*m_radius*sqrt(eps*(2.0-eps));

                if(m_delta>max)
                {
                    m_delta = max;
                }
                if(m_delta<min)
                {
                    m_delta = min;
                }
            }
        }

        /**
         * @brief return bool on whether point is valid
         */
        bool IsValid(){return m_valid;}

        /**
         * @brief get mesh spacing paramter
         */
        NekDouble GetDelta()
        {
            if(m_valid)
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
        Array<OneD, NekDouble> GetLoc(){return m_loc;}

        /**
         * @brief get normal vector
         */
        Array<OneD, NekDouble> GetNormal(){return m_norm;}

    private:

        /// x,y,z location
        Array<OneD, NekDouble> m_loc;
        ///normal vector of surface at point
        Array<OneD, NekDouble> m_norm;
        /// radius of curvature at point
        NekDouble m_radius;
        /// mesh spacing parameter at point
        NekDouble m_delta;
        /// valid point or not
        bool m_valid;
};

typedef boost::shared_ptr<CurvaturePoint> CurvaturePointSharedPtr;

}
}

#endif
