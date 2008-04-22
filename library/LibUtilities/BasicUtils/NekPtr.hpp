///////////////////////////////////////////////////////////////////////////////
//
// File NekPtr.hpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Scientific Computing and Imaging Institute,
// University of Utah (USA) and Department of Aeronautics, Imperial
// College London (UK).
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
// Description: A lightweight shared pointer.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_NEK_PTR_HPP
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_NEK_PTR_HPP

#include <LibUtilities/BasicUtils/ArrayPolicies.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

namespace Nektar
{
    template<typename T>
    class NekPtr
    {
        public:
            NekPtr() :
                m_data(0),
                m_count(0),
                m_size(0)
            {
            }
            
            explicit NekPtr(T* d) :
                m_data(d),
                m_count(MemoryManager<unsigned int>::RawAllocate(1)),
                m_size(1)
            {
                *m_count = 1;
            }
            
            NekPtr(T* d, unsigned int s) :
                m_data(d),
                m_count(MemoryManager<unsigned int>::RawAllocate(1)),
                m_size(s)
            {
                *m_count = 1;
            }
             
            NekPtr(const NekPtr<T>& rhs) :
                m_data(rhs.m_data),
                m_count(rhs.m_count),
                m_size(rhs.m_size)
            {
                if( m_count != 0 )
                {
                    *m_count += 1;
                }
            }
            
            NekPtr<T>& operator=(const NekPtr<T>& rhs)
            {
                NekPtr temp(rhs);
                Swap(temp);
                return *this;
            }
             
            ~NekPtr()
            {
                if( m_count == 0 )
                {
                    return;
                }
                
                *m_count -= 1;
                if( *m_count == 0 )
                {
                    ArrayDestructionPolicy<T>::Destroy(m_data, m_size);
                    MemoryManager<T>::RawDeallocate(m_data, m_size);
                    MemoryManager<unsigned int>::RawDeallocate(m_count, 1);
                }
            }
            
            T* get() { return m_data; }
            const T* get() const { return m_data; }
            
        private:
            void Swap(NekPtr<T>& rhs)
            {
                std::swap(m_data, rhs.m_data);
                std::swap(m_count, rhs.m_count);
                std::swap(m_size, rhs.m_size);
            }
            
            T* m_data;
            unsigned int* m_count;
            int m_size;
    };

}


#endif //NEKTAR_LIB_UTILITIES_BASIC_UTILS_NEK_PTR_HPP



