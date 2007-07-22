///////////////////////////////////////////////////////////////////////////////
//
// File: ConsistentObjectAccess.hpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: A wrapper around objects or pointers to objects which give
// the same interface for both types (i.e., -> will work for both).
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_CONSISTENT_ACCESS_OBJECT_HPP
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_CONSISTENT_ACCESS_OBJECT_HPP

#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/call_traits.hpp>

#include <iostream>

namespace Nektar
{
    template<typename DataType>
    class ConsistentObjectAccess
    {            
        public:
            ConsistentObjectAccess() :
                m_obj()
            {
            }
            
                template<typename ParamType1, typename ParamType2, typename ParamType3>
                ConsistentObjectAccess(const ParamType1& p1,
                                    const ParamType2& p2,
                                    const ParamType3& p3) :
                    m_obj(p1, p2, p3)
                {
                }
            
            explicit ConsistentObjectAccess(typename boost::call_traits<DataType>::const_reference rhs) :
                m_obj(rhs)
            {
            }
            
            ConsistentObjectAccess(const ConsistentObjectAccess<DataType>& rhs) :
                m_obj(rhs.m_obj)
            {
            }
            
            ConsistentObjectAccess<DataType>& operator=(const ConsistentObjectAccess<DataType>& rhs)
            {
                m_obj = rhs.m_obj;
                return *this;
            }
            
            ConsistentObjectAccess<DataType>& operator=(typename boost::call_traits<DataType>::const_reference rhs) 
            {
                m_obj = rhs;
                return *this;
            }
    
            ~ConsistentObjectAccess() 
            {
            }
    
            const DataType* operator->() const 
            {
                return &m_obj;
            }
    
            DataType* operator->() 
            {
                return &m_obj;
            }
            
            DataType& operator*()
            {
                return m_obj;
            }
    
            const DataType& operator*() const
            {
                return m_obj;
            }
    
            ConsistentObjectAccess<DataType> operator+=(const ConsistentObjectAccess<DataType>& rhs)
            {
                m_obj += rhs.m_obj;
                return *this;
            }
    
            ConsistentObjectAccess<DataType> operator*=(const ConsistentObjectAccess<DataType>& rhs)
            {
                m_obj *= rhs.m_obj;
                return *this;
            }
            
            
        private:
            DataType m_obj;
    };
    
    template<typename DataType>
    class ConsistentObjectAccess<boost::shared_ptr<DataType> >
    {            
        public:
            ConsistentObjectAccess() :
                m_obj(new DataType())
            {
            }
            
            template<typename ParamType1, typename ParamType2, typename ParamType3>
            ConsistentObjectAccess(typename boost::call_traits<ParamType1>::const_reference p1,
                                    typename boost::call_traits<ParamType2>::const_reference p2,
                                    typename boost::call_traits<ParamType3>::const_reference p3) :
                m_obj(new DataType(p1, p2, p3))
            {
            }
            
            explicit ConsistentObjectAccess(boost::shared_ptr<DataType> rhs) :
                m_obj(rhs)
            {
            }
                
            ConsistentObjectAccess(const ConsistentObjectAccess<DataType>& rhs) :
                m_obj(new DataType(*rhs.m_obj.get()))
            {
            }
            
            ConsistentObjectAccess<DataType>& operator=(const ConsistentObjectAccess<DataType>& rhs)
            {
                m_obj = boost::shared_ptr<DataType>(new DataType(*rhs.m_obj.get()));
                return *this;
            }
    
            ConsistentObjectAccess<DataType>& operator=(boost::shared_ptr<DataType> rhs)
            {
                m_obj = rhs;
                return *this;
            }
    
            ~ConsistentObjectAccess() 
            {
            }

            boost::shared_ptr<DataType> operator->() 
            {
                return m_obj;
            }
    
            boost::shared_ptr<const DataType> operator->() const 
            {
                return m_obj;
            }
    
            DataType& operator*()
            {
                return *m_obj;
            }
    
            const DataType& operator*() const
            {
                return *m_obj;
            }
            
            ConsistentObjectAccess<DataType> operator+=(const ConsistentObjectAccess<DataType>& rhs)
            {
                *(m_obj) += *(rhs.m_obj);
                return *this;
            }
    
            ConsistentObjectAccess<DataType> operator*=(const ConsistentObjectAccess<DataType>& rhs)
            {
                *(m_obj) *= *(rhs.m_obj);
                return *this;
            }
            
        private:
            boost::shared_ptr<DataType> m_obj;
    };
    
    template<typename DataType>
    bool operator==(const ConsistentObjectAccess<DataType>& lhs, const ConsistentObjectAccess<DataType>& rhs)
    {
        return *lhs == *rhs;
    }
            
    template<typename DataType>
    bool operator==(const ConsistentObjectAccess<DataType>& lhs, typename boost::call_traits<DataType>::const_reference rhs)
    {
        return *lhs == rhs;
    }
    
    template<typename DataType>
    bool operator==(typename boost::call_traits<DataType>::const_reference lhs, const ConsistentObjectAccess<DataType>& rhs)
    {
        return lhs == *rhs;
    }
    
    template<typename DataType>
    bool operator!=(const ConsistentObjectAccess<DataType>& lhs, const ConsistentObjectAccess<DataType>& rhs)
    {
        return !(lhs == rhs);
    }
            
    template<typename DataType>
    bool operator!=(const ConsistentObjectAccess<DataType>& lhs, typename boost::call_traits<DataType>::const_reference rhs)
    {
        return !(lhs == rhs);
    }
    
    template<typename DataType>
    bool operator!=(typename boost::call_traits<DataType>::const_reference lhs, const ConsistentObjectAccess<DataType>& rhs)
    {
        return !(lhs == rhs);
    }
    
    template<typename DataType>
    std::ostream& operator<<(std::ostream& os, const ConsistentObjectAccess<DataType>& rhs)
    {
        os << *rhs;
        return os;
    }

};
    
#endif //NEKTAR_LIB_UTILITIES_BASIC_UTILS_CONSISTENT_ACCESS_OBJECT_HPP

/**
    $Log: ConsistentObjectAccess.hpp,v $
    Revision 1.3  2007/07/20 00:39:36  bnelson
    Replaced boost::shared_ptr with Nektar::ptr

    Revision 1.2  2007/01/29 01:35:17  bnelson
    Removed memory manager requirements.

    Revision 1.1  2006/11/06 17:06:20  bnelson
    *** empty log message ***

 **/
 
