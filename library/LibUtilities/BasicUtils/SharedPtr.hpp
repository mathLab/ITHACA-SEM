///////////////////////////////////////////////////////////////////////////////
//
// File SharedPtr.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_SHARED_PTR_HPP
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_SHARED_PTR_HPP    

#include <boost/config.hpp>   // for broken compiler workarounds

#include <boost/assert.hpp>
#include <boost/checked_delete.hpp>
#include <boost/throw_exception.hpp>
#include <boost/detail/shared_count.hpp>
#include <boost/detail/workaround.hpp>

#include <memory>               // for std::auto_ptr
#include <algorithm>            // for std::swap
#include <functional>           // for std::less
#include <typeinfo>             // for std::bad_cast
#include <iosfwd>               // for std::basic_ostream

#ifdef BOOST_MSVC  // moved here to work around VC++ compiler crash
    #pragma warning(push)
    #pragma warning(disable:4284) // odd return type for operator->
#endif

namespace Nektar
{
    struct static_cast_tag {};
    struct const_cast_tag {};
    struct dynamic_cast_tag {};
    struct polymorphic_cast_tag {};

    template<class T>
    class const_ptr
    {
        private:
            typedef const_ptr<T> this_type;

        public:
            typedef T element_type;
            typedef T value_type;
            typedef T* pointer;
            typedef T& reference;
            typedef const T& const_reference;

            const_ptr() : 
                px(0), 
                pn()
            {
            }

            template<class Y>
            explicit const_ptr(Y* p ) : 
                px(p), 
                pn(p)
            {
            }

                
            template<class Y, class D> 
            const_ptr(Y* p, const D& d) : 
                px(p), 
                pn(p, d)
            {
            }

            template<class Y>
            const_ptr(const const_ptr<Y>& r) : 
                px(r.px), 
                pn(r.pn)
            {
            }

            const_ptr(const const_ptr<T>& rhs) :
                px(rhs.px),
                pn(rhs.pn)
            {
            }

            template<class Y>
            const_ptr(const const_ptr<Y>& r, static_cast_tag) : 
                px(static_cast<element_type*>(r.px)), 
                pn(r.pn)
            {
            }

            template<class Y>
            const_ptr(const const_ptr<Y>& r, const_cast_tag) : 
                px(const_cast<element_type*>(r.px)), 
                pn(r.pn)
            {
            }

            template<class Y>
            const_ptr(const const_ptr<Y>& r, dynamic_cast_tag) : 
                px(dynamic_cast<element_type *>(r.px)), 
                pn(r.pn)
            {
                if(px == 0) // need to allocate new counter -- the cast failed
                {
                    pn = boost::detail::shared_count();
                }
            }

            template<class Y>
            const_ptr(const const_ptr<Y>& r, polymorphic_cast_tag) : 
                px(dynamic_cast<element_type *>(r.px)), 
                pn(r.pn)
            {
                if(px == 0)
                {
                    boost::throw_exception(std::bad_cast());
                }
            }

            template<class Y>
            const const_ptr<T>& operator=(const const_ptr<Y>& r)
            {
                px = r.px;
                pn = r.pn; 
                return *this;
            }

            const const_ptr<T>& operator=(const const_ptr<T>& r)
            {
                px = r.px;
                pn = r.pn; 
                return *this;
            }

            void reset()
            {
                this_type().swap(*this);
            }

            template<class Y> 
            void reset(Y* p)
            {
                BOOST_ASSERT(p == 0 || p != px); // catch self-reset errors
                this_type(p).swap(*this);
            }

            template<class Y, class D> 
            void reset(Y* p, const D& d)
            {
                this_type(p, d).swap(*this);
            }

            const_reference operator*() const 
            {
                BOOST_ASSERT(px != 0);
                return *px;
            }

            const T* operator->() const
            {
                BOOST_ASSERT(px != 0);
                return px;
            }
            
            const T* get() const
            {
                return px;
            }

            operator bool() const
            {
                return px != 0;
            }

            // operator! is redundant, but some compilers need it
            bool operator!() const
            {
                return px == 0;
            }

            bool unique() const
            {
                return pn.unique();
            }

            long use_count() const
            {
                return pn.use_count();
            }

            void swap(const_ptr<T>& other)
            {
                std::swap(px, other.px);
                pn.swap(other.pn);
            }

            template<class Y> bool 
            internal_less(const const_ptr<Y>& rhs) const
            {
                return pn < rhs.pn;
            }

        private:
            T* px;                     
            boost::detail::shared_count pn;

    };  

    template<class T, class U> 
    inline bool operator==(const const_ptr<T>& a, const const_ptr<U>& b)
    {
        return a.get() == b.get();
    }

    template<class T, class U> 
    inline bool operator!=(const const_ptr<T>& a, const const_ptr<U>& b)
    {
        return a.get() != b.get();
    }


    template<class T, class U> 
    inline bool operator<(const const_ptr<T>& a, const const_ptr<U>& b)
    {
        return a.internal_less(b);
    }

    template<class T> 
    inline void swap(const const_ptr<T>& a, const const_ptr<T>& b)
    {
        a.swap(b);
    }

    template<class T, class U> 
    const_ptr<T> static_pointer_cast(const const_ptr<U>& r)
    {
        return const_ptr<T>(r, static_cast_tag());
    }

    template<class T, class U> 
    const_ptr<T> const_pointer_cast(const const_ptr<U>& r)
    {
        return const_ptr<T>(r, const_cast_tag());
    }

    template<class T, class U> \
    const_ptr<T> dynamic_pointer_cast(const const_ptr<U>& r)
    {
        return const_ptr<T>(r, dynamic_cast_tag());
    }


    template<class T> 
    inline const T* get_pointer(const const_ptr<T>& p)
    {
        return p.get();
    }
    

    template<class T>
    class ptr : public const_ptr<T>
    {
        public:
            typedef const_ptr<T> BaseType;
            typedef ptr<T> ThisType;
            using BaseType::element_type;
            using BaseType::value_type;
            using BaseType::pointer;
            using BaseType::reference;

            ptr() : 
                BaseType()
            {
            }

            template<class Y>
            explicit ptr(Y* p ) : 
                BaseType(p)
            {
            }

                
            template<class Y, class D> 
            ptr(Y* p, const D& d) : 
                BaseType(p, d)
            {
            }

            template<class Y>
            ptr(const ptr<Y>& r) : 
                BaseType(r)
            {
            }

            ptr(const ptr<T>& rhs) :
                BaseType(rhs)
            {
            }

            template<class Y>
            ptr(const ptr<Y>& r, static_cast_tag t) : 
                BaseType(r, t)
            {
            }

            template<class Y>
            ptr(const ptr<Y>& r, const_cast_tag t) : 
                BaseType(r, t)
            {
            }

            template<class Y>
            ptr(const const_ptr<Y>& r, dynamic_cast_tag t) : 
                BaseType(r, t)
            {
            }

            template<class Y>
            ptr(const const_ptr<Y>& r, polymorphic_cast_tag t) : 
                BaseType(r, t)
            {
            }

            template<class Y>
            const ptr<T>& operator=(const ptr<Y>& r)
            {
                BaseType::operator=(r);
                return *this;
            }

            const ptr<T>& operator=(const ptr<T>& r)
            {
                BaseType::operator=(r);
                return *this;
            }

            reference operator*() 
            {
                BOOST_ASSERT(px != 0);
                return *px;
            }

            T* operator->()
            {
                BOOST_ASSERT(px != 0);
                return px;
            }
            
            T* get()
            {
                return px;
            }

            using BaseType::operator*;
            using BaseType::operator->;
            using BaseType::get;
    };

}

#endif //NEKTAR_LIB_UTILITIES_BASIC_UTILS_SHARED_PTR_HPP
