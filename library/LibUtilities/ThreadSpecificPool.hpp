////////////////////////////////////////////////////////////////////////////////
//
//  File: ThreadSpecificPool.hpp
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
//  Description:
//
//
////////////////////////////////////////////////////////////////////////////////


#ifndef NEKATAR_LIB_UTILITES_THREAD_SPECIFIC_POOL_HPP
#define NEKATAR_LIB_UTILITES_THREAD_SPECIFIC_POOL_HPP

#include <boost/thread/tss.hpp>
#include <boost/pool/pool.hpp>

namespace Nektar
{
    /// \brief A memory pool which exists on a thread by thread basis.
    /// \param ByteSize The number of bytes in each chunk allocated by the pool.
    ///
    /// Provides a simple, thread specific memory pool that is based on byte size.
    /// The pool allocates and deallocates raw memory - the user is responsible for
    /// calling appropriate constructors/destructors when allocating objects.
    ///
    /// Example:
    ///
    /// \code
    /// ThreadSpecificPool<sizeof(TestClass)> pool;
    /// void* memory = pool.allocate();
    ///
    /// // Construct the object in the memory returned by the pool.
    /// TestClass* t = new (memory) TestClass;
    ///
    /// // Do stuff with t.
    ///
    /// // Destruct t and return it.
    /// t->~TestClass();
    /// pool.deallocate(t);
    /// \endcode
    template<unsigned int ByteSize>
    class ThreadSpecificPool
    {
        public:
            /// \brief Allocate a block of memory of size ByteSize.
            /// \throw std::bad_alloc if memory is exhausted.
            static void* Allocate()
            {
                // Check to see if we've created the pool for this thread.
                // If not, create one.  Note that the get method here will
                // always return NULL the first time allocate is called on
                // each thread.
                if( !m_pool.get() )
                {
                    m_pool.reset(new boost::pool<>(ByteSize));
                }

                return m_pool->malloc();
            }

            /// \brief Deallocate memory claimed by an earlier call to allocate.
            ///
            /// \attention It is an error to deallocate memory not allocated
            /// from this pool.  Doing this will result in undefined behavior.
            static void Deallocate(void* p)
            {
                if( m_pool.get() )
                {
                    m_pool->free(p);
                }
            }


        private:
            static boost::thread_specific_ptr<boost::pool<> > m_pool;
    };

    template<unsigned int ByteSize>
    boost::thread_specific_ptr<boost::pool<> > ThreadSpecificPool<ByteSize>::m_pool;
}



#endif //NEKATAR_LIB_UTILITES_THREAD_SPECIFIC_POOL_HPP

/**
    $Log: ThreadSpecificPool.hpp,v $
    Revision 1.1  2006/02/23 07:53:23  bnelson
    *** empty log message ***

**/

