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
#include <boost/thread/mutex.hpp>

#include <loki/Singleton.h>
#include <map>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

#include <cstring>

namespace Nektar
{
    namespace detail
    {
        /// \internal
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
        class ThreadSpecificPool
        {
            public:
                ThreadSpecificPool(size_t ByteSize) :
                    m_pool(),
                    m_blockSize(ByteSize),
                    m_mutex()
                {
                    // We can do the new in the constructor list because the thread specific 
                    // pointer doesn't have a supporting constructor.
                    m_pool = new boost::pool<>(ByteSize);
                }

                ~ThreadSpecificPool()
                {
                    // The documentation isn't particularly clear if delete needs to be called manually
                    // or if the thread specific pointer will call delete for me.  Looking through the 
                    // boost code doesn't make it any clearer. 
                }

                /// \brief Allocate a block of memory of size ByteSize.
                /// \throw std::bad_alloc if memory is exhausted.
                void* Allocate()
                {
                    boost::mutex::scoped_lock l(m_mutex);
                    void* result = m_pool->malloc();

#if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)
                    memset(result, 0, m_blockSize);
#endif //defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)

                    return result;
                }

                /// \brief Deallocate memory claimed by an earlier call to allocate.
                ///
                /// \attention It is an error to deallocate memory not allocated
                /// from this pool.  Doing this will result in undefined behavior.
                void Deallocate(const void* p)
                {
                    boost::mutex::scoped_lock l(m_mutex);
#if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)
                    // The idea here is to fill the returned memory with some known
                    // pattern, then detect that pattern on the allocate.  If the 
                    // pattern is no longer there then some memory corruption has 
                    // occurred.  However, I'm not sure how to distinguish between first
                    // time allocations and repeat allocations.

                    //memset(p, '+', m_pool->get_requested_size());
#endif //defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)

                    m_pool->free(const_cast<void*>(p));
                }


            private:
                //boost::thread_specific_ptr<boost::pool<> > m_pool;
                boost::pool<>* m_pool;
                size_t m_blockSize;
                boost::mutex m_mutex;
        };
    }

    class MemPool
    {
        public:
            typedef std::map<size_t, boost::shared_ptr<detail::ThreadSpecificPool> > PoolMapType;
            
        public:
            MemPool() :
                m_fourBytePool(4),
                m_pools(),
                m_upperBound(1024)
            {
                // The m_pools data member stores a collection of thread specific pools of varying size.  All memory requests
                // up to and including the largest pool size will be allocated from a pool (note that this means you may receive
                // more memory than you asked for).  For example, if there is a pool for 8 bytes and the next largest pool is 32 
                // bytes, then a request for 10 bytes will return a 32 byte chunk of memory from the 32 byte pool.
                
                typedef PoolMapType::value_type PairType;
                m_pools.insert(PairType(8, boost::shared_ptr<detail::ThreadSpecificPool>(new detail::ThreadSpecificPool(8))));
                m_pools.insert(PairType(16, boost::shared_ptr<detail::ThreadSpecificPool>(new detail::ThreadSpecificPool(16))));
                m_pools.insert(PairType(32, boost::shared_ptr<detail::ThreadSpecificPool>(new detail::ThreadSpecificPool(32))));
                m_pools.insert(PairType(64, boost::shared_ptr<detail::ThreadSpecificPool>(new detail::ThreadSpecificPool(64))));
                m_pools.insert(PairType(128, boost::shared_ptr<detail::ThreadSpecificPool>(new detail::ThreadSpecificPool(128))));
                m_pools.insert(PairType(256, boost::shared_ptr<detail::ThreadSpecificPool>(new detail::ThreadSpecificPool(256))));
                m_pools.insert(PairType(512, boost::shared_ptr<detail::ThreadSpecificPool>(new detail::ThreadSpecificPool(512))));
                m_pools.insert(PairType(1024, boost::shared_ptr<detail::ThreadSpecificPool>(new detail::ThreadSpecificPool(1024))));
            }
            
            ~MemPool()
            {
            }
            
            /// \brief Allocate a block of memory of size ByteSize.
            /// \throw std::bad_alloc if memory is exhausted.
            /// \param bytes The number of bytes to allocate.
            ///
            /// If the bytes parameter specifies a size that is handled by memory pools then the memory 
            /// is allocated from the pool.  Otherwise the memory is allocated with a call to new.
            ///
            /// Important: All memory allocated from this method must be returned to the pool
            /// via the Deallocate method.  Deleting pointers allocated from the memory pool with the 
            /// delete operator will result in undefined behavior.
            void* Allocate(size_t bytes)
            {
                if( bytes <= 4 )
                {
                    return m_fourBytePool.Allocate();
                }
                else if( bytes > m_upperBound )
                {
                    return ::operator new(bytes);
                }
                else
                {
                    PoolMapType::iterator iter = m_pools.lower_bound(bytes);
                    ASSERTL1(iter != m_pools.end(), "The memory manager is mishandling a memory request for " +
                        boost::lexical_cast<std::string>(bytes) + " bytes of memory.");
                    
                    return (*iter).second->Allocate();
                }
            }

            /// \brief Deallocate memory claimed by an earlier call to allocate.
            ///
            /// \attention It is an error to deallocate memory not allocated
            /// from this pool.  Doing this will result in undefined behavior.
            void Deallocate(void* p, size_t bytes)
            {
                if( bytes <= 4 )
                {
                    m_fourBytePool.Deallocate(p);
                }
                else if( bytes > m_upperBound )
                {
                    ::operator delete(p);
                }
                else
                {
                    PoolMapType::iterator iter = m_pools.lower_bound(bytes);
                    ASSERTL1(iter != m_pools.end(), "The memory manager is mishandling a memory request for " +
                        boost::lexical_cast<std::string>(bytes) + " bytes of memory.");
                    
                    (*iter).second->Deallocate(p);
                }
            }
            
        private:
            detail::ThreadSpecificPool m_fourBytePool;
            std::map<size_t, boost::shared_ptr<detail::ThreadSpecificPool> > m_pools;
            size_t m_upperBound;
    };

    LIB_UTILITIES_EXPORT MemPool& GetMemoryPool();

}



#endif //NEKATAR_LIB_UTILITES_THREAD_SPECIFIC_POOL_HPP

