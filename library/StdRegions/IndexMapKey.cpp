///////////////////////////////////////////////////////////////////////////////
//
// File IndexMapKey.cpp
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
// Description: Definition of IndexMapKey
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/IndexMapKey.h>
#include <StdRegions/StdRegions.hpp>

namespace Nektar
{
    namespace StdRegions
    {
        IndexMapKey::IndexMapKey(
            const StdRegions::IndexMapType  indexmapType,
            const LibUtilities::ShapeType   shapeType,
            const unsigned short            p, 
            const unsigned short            q,
            const unsigned short            r,
            const unsigned short            entityID,
            const StdRegions::Orientation   orientation)
            : m_indexMapType(indexmapType),
              m_shapeType(shapeType),
              m_p(p),
              m_q(q),
              m_r(r),
              m_entityID(entityID),
              m_orientation(orientation)
        {
        }

        IndexMapKey::IndexMapKey(const IndexMapKey& rhs,
                                 const StdRegions::IndexMapType indexmapType):
            m_indexMapType (indexmapType),
            m_shapeType(rhs.m_shapeType),
            m_p            (rhs.m_p),
            m_q            (rhs.m_q),
            m_r            (rhs.m_r),
            m_entityID     (rhs.m_entityID),
            m_orientation  (rhs.m_orientation)
        {
        }

        IndexMapKey::IndexMapKey(const IndexMapKey& rhs):
            m_indexMapType (rhs.m_indexMapType),
            m_shapeType(rhs.m_shapeType),
            m_p            (rhs.m_p),
            m_q            (rhs.m_q),
            m_r            (rhs.m_r),
            m_entityID     (rhs.m_entityID),
            m_orientation  (rhs.m_orientation)
        {
        }

        bool IndexMapKey::opLess::operator()(const IndexMapKey &lhs, 
                                             const IndexMapKey &rhs) const
        {        
            return (lhs.m_indexMapType < rhs.m_indexMapType);
        }

        bool operator<(const IndexMapKey &lhs, const IndexMapKey &rhs)
        {   
            if(lhs.m_indexMapType < rhs.m_indexMapType)
            {
                return true;
            }

            if(lhs.m_indexMapType > rhs.m_indexMapType)
            {
                return false;
            }
            
            if(lhs.m_shapeType < rhs.m_shapeType)
            {
                return true;
            }
            if(lhs.m_shapeType > rhs.m_shapeType)
            {
                return false;
            }

            if(lhs.m_p < rhs.m_p)
            {
                return true;
            }
            if(lhs.m_p > rhs.m_p)
            {
                return false;
            }

            if(lhs.m_q < rhs.m_q)
            {
                return true;
            }
            if(lhs.m_q > rhs.m_q)
            {
                return false;
            }

            if(lhs.m_r < rhs.m_r)
            {
                return true;
            }
            if(lhs.m_r > rhs.m_r)
            {
                return false;
            }

            if(lhs.m_entityID < rhs.m_entityID)
            {
                return true;
            }
            if(lhs.m_entityID > rhs.m_entityID)
            {
                return false;
            }

            if(lhs.m_orientation < rhs.m_orientation)
            {
                return true;
            }
            if(lhs.m_orientation > rhs.m_orientation)
            {
                return false;
            }

            return false;
        }

        bool operator==(const IndexMapKey &lhs, const IndexMapKey &rhs)
        {
            if(lhs.m_indexMapType != rhs.m_indexMapType)
            {
                return false;
            }

            if(lhs.m_shapeType != rhs.m_shapeType)
            {
                return false;
            }

            if(lhs.m_p != rhs.m_p)
            {
                return false;
            }

            if(lhs.m_q != rhs.m_q)
            {
                return false;
            }

            if(lhs.m_r != rhs.m_r)
            {
                return false;
            }

            if(lhs.m_entityID != rhs.m_entityID)
            {
                return false;
            }

            if(lhs.m_orientation != rhs.m_orientation)
            {
                return false;
            }

            return true;
        }

        std::ostream& operator<<(std::ostream& os, const IndexMapKey& rhs)
        {
            os << "IndexMapType: " << IndexMapTypeMap[rhs.GetIndexMapType()] 
               << std::endl;
            return os;
        }
    } // end StdRegion namespace
} // end Nektar namespace
