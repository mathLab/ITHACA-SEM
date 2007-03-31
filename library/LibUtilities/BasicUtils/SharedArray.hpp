///////////////////////////////////////////////////////////////////////////////
//
// File SharedArray.hpp
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
// Description: 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_SHARED_ARRAY_HPP
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_SHARED_ARRAY_HPP

// This is a copy of boost::shared_array with minor modifications.

#include <boost/config.hpp>   // for broken compiler workarounds

#include <boost/assert.hpp>
#include <boost/checked_delete.hpp>

#include <boost/detail/shared_count.hpp>
#include <boost/detail/workaround.hpp>
#include <boost/bind.hpp>

#include <cstddef>            // for std::ptrdiff_t
#include <algorithm>          // for std::swap
#include <functional>         // for std::less

namespace Nektar
{

    /// \brief Wraps a reference counted T* and deletes it when it is no longer referenced.
    ///
    /// SharedArray is a replacement for boost::shared_array.  It addresses the following 
    /// problem with boost::shared_array
    ///
    /// A "const boost::shared_array<T>" still allows you to modify the internal data
    /// held by the array.  In other words, it is possible to do the following:
    ///
    /// const boost::shared_array<double>& foo();
    ///
    /// const boost::shared_array<double>& a = foo();
    /// a[0] = 1.0;
    ///
    /// With SharedArrays, a "const SharedArray" will not allow you to modify the 
    /// underlying data.  This has some unusual implications.  Consider the following:
    ///
    /// const SharedArray<double>& foo();
    ///
    /// SharedArray<double> obj = foo();
    ///
    /// The variable obj now contains a <em>copy</em> of the SharedArray returned by foo.
    /// Why?  If it were not a copy, then, since obj is not a constant variable, we could
    /// do the following:
    ///
    /// obj[0] = 1.0;
    ///
    /// This is an error - we shouldn't be able to modify the underlying array since we returned 
    /// it as a constant.
    ///
    /// The solution is to write your call as follows:
    ///
    /// const SharedArray<double>& obj = foo();
    ///
    /// Note that this is a reference and therefore the internal reference count is not updated
    /// and the reference will be invalid if the SharedArray returned by foo is ever destroyed.  In 
    /// other words, it behaves just like any other C++ reference object.
    ///
    /// If you are returning a non-const shared array, such as
    ///
    /// SharedArray<double> bar();
    ///
    /// Then the following:
    ///
    /// SharedArray<double> obj = bar();
    ///
    /// Creates a new shared array which refers to the original array and updates the internal 
    /// reference count just as the boost::shared_array does.
    ///
    /// Additionally, it is no longer possible to assign a const SharedArray to a non-const SharedArray.
    ///
    /// const SharedArray<double> rhs;
    /// SharedArray<double> lhs;
    ///
    /// lhs = rhs; // Compile error.
    ///
    /// To do this, call the Assign method:
    ///
    /// lhs.Assign(rhs);
    ///
    /// This creates a copy of the underlying array and gives it to lhs.
    ///
    /// An additional features of a SharedArray is the ability to create a new shared array
    /// with an offset:
    ///
    /// SharedArray<double> coeffs = MemoryManager::AllocateSharedArray<double>(10);
    /// SharedArray<double> offset = coefffs+5; // offset[5] = coeffs[0]
    /// 
    /// coeffs and offset point to the same underlying array, so if coeffs goes out of scope the array
    /// is not deleted until offset goes out of scope.  
    template<class T> 
    class SharedArray
    {
        private:

            // Borland 5.5.1 specific workarounds
            typedef boost::checked_array_deleter<T> deleter;
            typedef SharedArray<T> this_type;

        public:
            typedef T element_type;

            SharedArray() : 
                px(0), 
                pn((T*)0, deleter()),
                m_offset(0),
                m_size(0)
            {
            }

            explicit SharedArray(unsigned int s) :
                px(new T[s]),
                pn(px, deleter()),
                m_offset(0),
                m_size(0)
            {
            }

            SharedArray(T* p, unsigned int s) : 
                px(p), 
                pn(p, deleter()),
                m_offset(0),
                m_size(s)
            {
            }
            
            SharedArray(SharedArray<T>& rhs, unsigned int offset) :
                px(rhs.px),
                pn(rhs.pn),
                m_offset(offset),
                m_size(rhs.m_size)
            {
            }

            //
            // Requirements: D's copy constructor must not throw
            //
            // shared_array will release p by calling d(p)
            //

            template<class D> 
            SharedArray(T* p, unsigned int s, D d) : 
                px(p), 
                pn(p, d),
                m_offset(0),
                m_size(s)
            {
            }

            //  generated destructor is fine
            SharedArray(SharedArray<T>& rhs) :
                px(rhs.px),
                pn(rhs.pn),
                m_offset(rhs.m_offset),
                m_size(rhs.m_size)
            {
            }

            SharedArray<T>& operator=(SharedArray<T>& rhs)
            {
                if( this != &rhs )
                {
                    px = rhs.px;
                    pn = rhs.pn;
                    m_offset = rhs.m_offset;
                    m_size = rhs.m_size;
                }

                return *this;
            }

            void reset()
            {
                this_type().swap(*this);
            }

            void reset(T* p, unsigned int size)
            {
                BOOST_ASSERT(p == 0 || p != px);
                this_type(p, size).swap(*this);
            }

            template <class D> void reset(T * p, unsigned int size, D d)
            {
                this_type(p, size, d).swap(*this);
            }

            T& operator[] (std::ptrdiff_t i)  // never throws
            {
                BOOST_ASSERT(px != 0);
                BOOST_ASSERT(i >= 0);
                return px[i + m_offset];
            }

            const T& operator[] (std::ptrdiff_t i) const  // never throws
            {
                BOOST_ASSERT(px != 0);
                BOOST_ASSERT(i >= 0);
                return px[i + m_offset];
            }
            
            T* get() // never throws
            {
                return px + m_offset;
            }

            const T* get() const 
            {
                return px + m_offset;
            }
            
            unsigned int GetSize() const { return m_size; }

            // implicit conversion to "bool"

        #if defined(__SUNPRO_CC) && BOOST_WORKAROUND(__SUNPRO_CC, <= 0x530)

            operator bool () const
            {
                return px != 0;
            }

        #elif defined(__MWERKS__) && BOOST_WORKAROUND(__MWERKS__, BOOST_TESTED_AT(0x3003))
            typedef T * (this_type::*unspecified_bool_type)() const;
            
            operator unspecified_bool_type() const // never throws
            {
                return px == 0? 0: &this_type::get;
            }

        #else 

            typedef T * this_type::*unspecified_bool_type;

            operator unspecified_bool_type() const // never throws
            {
                return px == 0? 0: &this_type::px;
            }

        #endif

            bool operator! () const // never throws
            {
                return px == 0;
            }

            bool unique() const // never throws
            {
                return pn.unique();
            }

            long use_count() const // never throws
            {
                return pn.use_count();
            }

            void swap(SharedArray<T>& other) // never throws
            {
                std::swap(px, other.px);
                pn.swap(other.pn);
                std::swap(m_offset, other.m_offset);
                std::swap(m_size, other.m_size);
            }

            typedef T* iterator;
            typedef const T* const_iterator;

            iterator begin() { return px + m_offset; }
            iterator end() { return px + m_size; }

            const_iterator begin() const { return px + m_offset; }
            const_iterator end() const { return px + m_size; }

        private:
            //SharedArray<T>& operator=(const SharedArray<T>& rhs);

            T* px;                     // contained pointer
            boost::detail::shared_count pn;    // reference counter
            unsigned int m_offset;
            unsigned int m_size;

    };  

    template<class T> inline bool operator==(SharedArray<T> const & a, SharedArray<T> const & b) // never throws
    {
        return a.get() == b.get();
    }

    template<class T> inline bool operator!=(SharedArray<T> const & a, SharedArray<T> const & b) // never throws
    {
        return a.get() != b.get();
    }

    template<class T> inline bool operator<(SharedArray<T> const & a, SharedArray<T> const & b) // never throws
    {
        return std::less<T*>()(a.get(), b.get());
    }

    template<class T> void swap(SharedArray<T> & a, SharedArray<T> & b) // never throws
    {
        a.swap(b);
    }

    template<typename T>
    SharedArray<T> operator+(SharedArray<T>& lhs, unsigned int offset)
    {
        return SharedArray<T>(lhs, offset);
    }

    template<typename T>
    SharedArray<T> operator+(unsigned int offset, SharedArray<T>& rhs)
    {
        return SharedArray<T>(rhs, offset);
    }
}

#endif //NEKTAR_LIB_UTILITIES_BASIC_UTILS_SHARED_ARRAY_HPP
