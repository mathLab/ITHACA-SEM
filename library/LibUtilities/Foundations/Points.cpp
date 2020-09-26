///////////////////////////////////////////////////////////////////////////////
//
// File Points.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
// 
// Description: 1D Points definitions 
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Points.h>

namespace Nektar
{
    namespace LibUtilities 
    {
        /**
        * @class PointsKey
        * Specification for a set of points. This includes the total number of
        * points, as well as their distribution.
        */

        bool operator==(const PointsKey &lhs, const PointsKey &rhs)
        {
            return (lhs.m_numpoints == rhs.m_numpoints &&
                    lhs.m_pointstype == rhs.m_pointstype);
        }

        bool operator<(const PointsKey &lhs, const PointsKey &rhs)
        {
            if (lhs.m_pointstype < rhs.m_pointstype)
            {
                return true;
            }

            if (lhs.m_pointstype > rhs.m_pointstype)
            {
                return false;
            }

            if (lhs.m_numpoints < rhs.m_numpoints)
            {
                return true;
            }
            
            if (lhs.m_numpoints > rhs.m_numpoints)
            {
                return false;
            }
            
            if(lhs.m_factor < rhs.m_factor)
            {
                return true;
            }
            
            if(lhs.m_factor > rhs.m_factor)
            {
                return true;
            }

            return false;
        }

        bool PointsKey::opLess::operator()(const PointsKey &lhs, const PointsKey &rhs) const
        {
            return (lhs.m_pointstype < rhs.m_pointstype);
        }

        std::ostream& operator<<(std::ostream& os, const PointsKey& rhs)
        {
            os << "NumPoints: " << rhs.GetNumPoints() << " PointsType: " << kPointsTypeStr[rhs.GetPointsType()] << std::endl;

            return os;
        }
        
        
        /**
         * @class Nektar::LibUtilities::Points
         * This encapsulates a set of points, specified by a PointKey. The 
         * class stores not only the point coordinates, but also the 
         * integration weights and derivative matrix coefficients. Memory is 
         * allocated from the memory pool if in use.
         */

    }
}

