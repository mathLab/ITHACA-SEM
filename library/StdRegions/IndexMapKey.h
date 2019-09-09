///////////////////////////////////////////////////////////////////////////////
//
// File IndexMapKey.h
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
// Description: Headers for Index Maps
//
///////////////////////////////////////////////////////////////////////////////

#ifndef INDEXMAPKEY_H
#define INDEXMAPKEY_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdRegionsDeclspec.h>
#include <ostream>
#include <memory>

namespace Nektar
{
    namespace StdRegions
    {
        struct IndexValue
        {
            unsigned short index;
            short sign;
        };
		
        typedef Array<OneD, IndexValue> IndexMapValues;
		
        class IndexMapKey
        {
        public:
            STD_REGIONS_EXPORT IndexMapKey(
                const StdRegions::IndexMapType  indexmapType,
                const LibUtilities::ShapeType   shapeType,
                const unsigned short            p, 
                const unsigned short            q,
                const unsigned short            r,
                const unsigned short            entityID    = 0,
                const StdRegions::Orientation   orientation = eNoOrientation);
            
            STD_REGIONS_EXPORT IndexMapKey(const IndexMapKey& rhs,const StdRegions::IndexMapType indexmapType);

            STD_REGIONS_EXPORT IndexMapKey(const IndexMapKey& rhs);

            virtual ~IndexMapKey()
            {
            }

            // Used to lookup the create function in NekManager.
            struct opLess
            {
                STD_REGIONS_EXPORT bool operator()(const IndexMapKey &lhs, const IndexMapKey &rhs) const;
            };

            // Used for finding value given the key in NekManager.
            STD_REGIONS_EXPORT friend bool operator<(const IndexMapKey &lhs, const IndexMapKey &rhs);
            STD_REGIONS_EXPORT friend bool operator==(const IndexMapKey &lhs, const IndexMapKey &rhs);
            STD_REGIONS_EXPORT friend bool opLess::operator()(const IndexMapKey &lhs, const IndexMapKey &rhs) const;

            IndexMapType GetIndexMapType() const
            {
                return m_indexMapType;
            }
			
            Orientation GetIndexOrientation() const
            {
                return m_orientation;
            }
			
            int GetIndexEntity() const
            {
                return m_entityID;
            }

        protected:
			
            IndexMapType  m_indexMapType;
            
            LibUtilities::ShapeType m_shapeType;
			
            unsigned short m_p;
            unsigned short m_q;
            unsigned short m_r;
			
            unsigned short m_entityID;
			
            Orientation m_orientation;
			
        private:
			
            IndexMapKey();
        };
        //=======================================================================================

        STD_REGIONS_EXPORT std::ostream& operator<<(std::ostream& os, const IndexMapKey& rhs);

        typedef  std::shared_ptr<IndexMapKey> IndexMapKeySharedPtr;
        typedef  std::shared_ptr<IndexMapValues> IndexMapValuesSharedPtr;
    } // end of namespace
} // end of namespace

#endif
