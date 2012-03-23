////////////////////////////////////////////////////////////////////////////////
//
//  File: NekMemoryManager.hpp
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
//  A memory manager that allocates memory from thread specific pools.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_NEK_MEMORY_MANAGER_H
#define NEKTAR_LIB_UTILITIES_NEK_MEMORY_MANAGER_H

#include <LibUtilities/Memory/ThreadSpecificPool.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

#include <boost/mpl/contains.hpp>
#include <boost/mpl/list_c.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/bind.hpp>
#include <boost/type_traits.hpp>

#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/preprocessor/repetition/enum_trailing.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>
#include <boost/preprocessor/repetition/enum_trailing_params.hpp>
#include <boost/preprocessor/control/if.hpp>
#include <boost/preprocessor/facilities/empty.hpp>
#include <boost/preprocessor/comparison/equal.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>
#include <boost/preprocessor/repetition/enum_trailing_binary_params.hpp>

#include <vector>

#ifdef max
#undef max
#endif

namespace Nektar
{

#ifndef NEKTAR_MAX_MEMORY_MANAGER_CONSTRUCTOR_ARGS
#define NEKTAR_MAX_MEMORY_MANAGER_CONSTRUCTOR_ARGS 15
#endif //NEKTAR_MAX_MEMORY_MANAGER_CONSTRUCTOR_ARGS

    /// @brief General purpose memory allocation routines with the ability
    ///        to allocate from thread specific memory pools.
    ///
    /// If compiled with NEKTAR_MEMORY_POOL_ENABLED, the MemoryManager
    /// allocates from thread specific memory pools for small objects.
    /// Large objects are managed with the system supplied new/delete.
    /// These memory pools provide faster allocation and deallocation
    /// of small objects (particularly useful for shared pointers which
    /// allocate many 4 byte objects).
    ///
    /// @warning All memory allocated from the memory manager must be returned
    /// to the memory manager.  Calling delete on memory allocated from the
    /// manager will likely cause undefined behavior.  A particularly subtle
    /// violation of this rule occurs when giving memory allocated from the
    /// manager to a shared pointer.
    /// @code
    /// boost::shared_ptr<Obj> f(MemoryManager<Obj>::Allocate());
    /// @endcode
    /// Shared pointers call delete when they go out of scope, so this line of
    /// code will cause problems.  Instead, you should call the
    /// AllocateSharedPtr method:
    /// @code
    /// boost::shared_ptr<Obj> f = MemoryManager<Obj>::AllocateSharedPtr();
    /// @endcode
    template<typename DataType>
    class MemoryManager
    {
        private:
            template<typename ObjectType, typename CustomDeallocator>
            class DeallocateSharedPtr
            {
                public:
                    explicit DeallocateSharedPtr(const CustomDeallocator& d) :
                        m_dealloc(d)
                    {
                    }

                    void operator()(ObjectType*& m) const
                    {
                        m_dealloc();
                        MemoryManager<ObjectType>::Deallocate(m);
                    }

                private:
                    CustomDeallocator m_dealloc;

            };

            class DefaultCustomDeallocator
            {
                public:
                    void operator()() const
                    {
                    }
            };

        public:
            /// @brief Deallocate a pointer allocated by
            /// MemoryManager::Allocate.
            /// @note Results are undefined if called with a pointer to
            /// something that was not allocated with the memory manager.
            ///
            /// Use this method to deallocate a pointer you have allocated from
            /// the MemoryManager using the Allocate method.
            ///
            /// Example:
            /// @code
            /// CustObj* c = MemoryManager::Allocate<CustObj>();
            /// MemoryManager::Deallocate(c);
            /// @endcode
            static void Deallocate(DataType*& data)
            {
                #ifdef NEKTAR_MEMORY_POOL_ENABLED
                    data->~DataType();
                    GetMemoryPool().Deallocate(data, sizeof(DataType));
                #else //NEKTAR_MEMORY_POOL_ENABLED
                    delete data;
                #endif //NEKTAR_MEMORY_POOL_ENABLED

                data = NULL;
            }

            #ifdef NEKTAR_MEMORY_POOL_ENABLED
                /// @brief Allocates a single object from the memory pool.
                /// @throws unknown If the object throws an exception during
                /// construction, this method will catch it, release the memory
                /// back to the pool, then rethrow it.
                ///
                /// The allocated object must be returned to the memory pool
                /// via Deallocate.
                static DataType* Allocate()
                {
                    DataType* result = static_cast<DataType*>(GetMemoryPool().Allocate(sizeof(DataType)));

                    if( result )
                    {
                        try
                        {
                            new (result) DataType();
                        }
                        catch(...)
                        {
                            GetMemoryPool().Deallocate(result,
                                                              sizeof(DataType));
                            throw;
                        }
                    }

                    return result;
                }
                #define ALLOCATE_METHOD_GENERATOR(z, i, methodName) \
                template<BOOST_PP_ENUM_PARAMS(i, typename Arg)> \
                static DataType* methodName(BOOST_PP_ENUM_BINARY_PARAMS(i, const Arg, & arg)) \
                { \
                    DataType* result = static_cast<DataType*>(GetMemoryPool().Allocate(sizeof(DataType))); \
                    \
                    if( result ) \
                    { \
                        try \
                        { \
                            new (result) DataType(BOOST_PP_ENUM_PARAMS(i, arg)); \
                        } \
                        catch(...) \
                        { \
                            GetMemoryPool().Deallocate(result, sizeof(DataType)); \
                            throw; \
                        } \
                    } \
                    \
                    return result; \
                }
            #else //NEKTAR_MEMORY_POOL_ENABLED
                /// @brief Allocates a single object from the memory pool.
                /// @throws unknown Any exception thrown by DataType's default
                /// constructor will propogate through this method.
                ///
                /// The allocated object must be returned to the memory pool
                /// via Deallocate.
                static DataType* Allocate()
                {
                    return new DataType();
                }

                #define ALLOCATE_METHOD_GENERATOR(z, i, methodName) \
                template<BOOST_PP_ENUM_PARAMS(i, typename Arg)> \
                static DataType* methodName(BOOST_PP_ENUM_BINARY_PARAMS(i, const Arg, & arg)) \
                { \
                    return new DataType(BOOST_PP_ENUM_PARAMS(i, arg)); \
                }
            #endif //NEKTAR_MEMORY_POOL_ENABLED

            BOOST_PP_REPEAT_FROM_TO(1, NEKTAR_MAX_MEMORY_MANAGER_CONSTRUCTOR_ARGS, ALLOCATE_METHOD_GENERATOR, Allocate);

            /// @brief Allocate a shared pointer from the memory pool.
            ///
            /// The shared pointer does not need to be returned to the memory
            /// pool. When the reference count to this object reaches 0, the
            /// shared pointer will automatically return the memory.
            static boost::shared_ptr<DataType> AllocateSharedPtr()
            {
                return AllocateSharedPtrD(DefaultCustomDeallocator());
            }

            /// @def ALLOCATE_SHARED_PTR_METHOD_GENERATOR
            /// @brief Generator for allocating shared pointers to objects with
            /// constructors having a variable number of parameters.
            ///
            /// Uses the Boost preprocessor macros
            /// - BOOST_PP_ENUM_PARAMS(idx, text) which generates a list of
            /// parameters, in this case: "typename Arg0, typename Arg1, ..."
            /// used for the template definition.
            /// - BOOST_PP_ENUM_BINARY_PARAMS(idx, type, name) which generates
            /// a list of parameters and variables, in this case:
            /// "Arg0& arg0, Arg1& arg1, ..." used for the function prototype.
            ///
            /// @note All the parameter lists are references so whenever the
            /// MemoryManager::AllocateSharedPtr(...) is used, all parameters
            /// must be passed <b>by reference</b>. Consequently, fundamental
            /// datatype parameters must be defined locally first. So
            /// @code
            ///   MemoryManager::AllocateSharedPtr(Var, true);
            /// @endcode
            /// must be replaced by
            /// @code
            ///   bool flag = true;
            ///   MemoryManager::AllocateSharedPtr(Var, flag);
            /// @endcode
            #define ALLOCATE_SHARED_PTR_METHOD_GENERATOR(z, i, methodName) \
            template<BOOST_PP_ENUM_PARAMS(i, typename Arg)> \
            static boost::shared_ptr<DataType> methodName(BOOST_PP_ENUM_BINARY_PARAMS(i, const Arg, & arg)) \
            { \
                return AllocateSharedPtrD(DefaultCustomDeallocator(), BOOST_PP_ENUM_PARAMS(i, arg)); \
            }

            #define ALLOCATE_SHARED_PTR_METHOD_WITH_DEALLOCATOR_GENERATOR(z, i, methodName) \
            template<typename DeallocatorType BOOST_PP_ENUM_TRAILING_PARAMS(i, typename Arg)> \
            static boost::shared_ptr<DataType> methodName(const DeallocatorType& d BOOST_PP_ENUM_TRAILING_BINARY_PARAMS(i, Arg, & arg)) \
            { \
                DataType* data = Allocate(BOOST_PP_ENUM_PARAMS(i, arg)); \
                return boost::shared_ptr<DataType>(data, DeallocateSharedPtr<DataType, DeallocatorType>(d)); \
            }

            BOOST_PP_REPEAT_FROM_TO(0, NEKTAR_MAX_MEMORY_MANAGER_CONSTRUCTOR_ARGS, ALLOCATE_SHARED_PTR_METHOD_WITH_DEALLOCATOR_GENERATOR, AllocateSharedPtrD);
            BOOST_PP_REPEAT_FROM_TO(1, NEKTAR_MAX_MEMORY_MANAGER_CONSTRUCTOR_ARGS, ALLOCATE_SHARED_PTR_METHOD_GENERATOR, AllocateSharedPtr);


            /// \brief Allocates a chunk of raw, uninitialized memory, capable of holding NumberOfElements objects.
            /// \param NumberOfElements The number of elements the array should be capable of holding.
            ///
            /// This method is not meant to be called by client code.  Use Array instead.
            /// Any memory allocated from this method must be returned to the memory pool
            /// via RawDeallocate.  Failure to do so will result in memory leaks and undefined
            /// behavior.
            static DataType* RawAllocate(unsigned int NumberOfElements)
            {
                #ifdef NEKTAR_MEMORY_POOL_ENABLED
                    return static_cast<DataType*>(GetMemoryPool().Allocate(sizeof(DataType)*NumberOfElements));
                #else //NEKTAR_MEMORY_POOL_ENABLED
                    return static_cast<DataType*>(::operator new(NumberOfElements * sizeof(DataType)));
                #endif //NEKTAR_MEMORY_POOL_ENABLED
            }


            /// \brief Deallocates memory allocated from RawAllocate.
            /// \param array A pointer to the memory returned from RawAllocate.
            /// \param NumberOfElements The number of object held in the array.
            ///
            /// This method is not meant to be called by client code.  Use Array instead.
            /// Only memory allocated via RawAllocate should be returned to the pool here.
            static void RawDeallocate(DataType* array, unsigned int NumberOfElements)
            {
                #ifdef NEKTAR_MEMORY_POOL_ENABLED
                    GetMemoryPool().Deallocate(array, sizeof(DataType)*NumberOfElements);
                #else //NEKTAR_MEMORY_POOL_ENABLED
                    ::operator delete(array);
                #endif //NEKTAR_MEMORY_POOL_ENABLED

            }

            /////////////////////////////////////////////////////////////////
            ///\name Allocator Interface
            /// The allocator interface allows a MemoryManager object to be used
            /// in any object that allows an allocator parameter, such as STL
            /// containers.
            /////////////////////////////////////////////////////////////////
            typedef DataType value_type;
            typedef size_t size_type;
            typedef ptrdiff_t difference_type;
            typedef DataType* pointer;
            typedef const DataType* const_pointer;
            typedef DataType& reference;
            typedef const DataType& const_reference;

            MemoryManager() {}
            template<typename T>
            MemoryManager(const MemoryManager<T>& rhs) {}
            ~MemoryManager() {}

            pointer address(reference r) const { return &r; }
            const_pointer address(const_reference r) const { return &r; }

            pointer allocate(size_type n, std::allocator<void>::const_pointer hint = 0)//typename MemoryManager<void>::pointer hint = 0)
            {
                return RawAllocate(n);
            }

            void deallocate(pointer p, size_type n)
            {
                return RawDeallocate(p, n);
            }

            void construct(pointer p, const_reference val)
            {
                new(p) DataType(val);
            }

            void destroy(pointer p)
            {
                p->~DataType();
            }

            size_type max_size()
            {
                return std::numeric_limits<size_type>::max()/sizeof(DataType);
            }

            template<typename U>
            struct rebind
            {
                typedef MemoryManager<U> other;
            };

            /////////////////////////////////////////////////////////////////
            ///@}
            /////////////////////////////////////////////////////////////////

        private:

    };

    template<typename DataType>
    bool operator==(const MemoryManager<DataType>& lhs, const MemoryManager<DataType>& rhs)
    {
        return true;
    }

    template<typename DataType>
    bool operator!=(const MemoryManager<DataType>& lhs, const MemoryManager<DataType>& rhs)
    {
        return !(lhs == rhs);
    }

}

#endif //NEKTAR_LIB_UTILITIES_NEK_MEMORY_MANAGER_H


/**
    $Log: NekMemoryManager.hpp,v $
    Revision 1.20  2009/09/24 10:46:45  cbiotto
    Increasing NEKTAR_MAX_MEMORY_MANAGER_CONSTRUCTOR_ARGS to 11

    Revision 1.19  2008/06/10 06:00:37  bnelson
    Updated documentation.

    Revision 1.18  2008/05/16 05:43:19  bnelson
    Updated the memory manager so it is faster choosing the allocator to use, doesn't use the pools for anything larger than 1024 bytes, and doesn't issue a warning for large allocations.

    Revision 1.17  2008/04/06 22:29:13  bnelson
    Reimplmented shared arrays for speed improvements.

    Revision 1.16  2008/01/02 18:17:22  bnelson
    Removed commented code.

    Revision 1.15  2007/08/25 04:21:49  bnelson
    *** empty log message ***

    Revision 1.14  2007/07/27 00:22:37  bnelson
    Memory manager now accepts non-const parameters to the allocate methods.

    Revision 1.13  2007/07/22 23:03:28  bnelson
    Backed out Nektar::ptr.

    Revision 1.12  2007/07/20 00:19:42  bnelson
    Replaced boost::shared_ptr with Nektar::ptr

    Revision 1.11  2007/05/14 23:50:16  bnelson
    Updated MemoryManager to implement the allocator interface.

    Revision 1.10  2007/04/29 00:30:05  jfrazier
    Converted tmp space methods to return 1D multi_arrays.

    Revision 1.9  2007/04/06 04:36:21  bnelson
    Updated for const-correctness.

    Revision 1.8  2007/03/29 18:42:58  bnelson
    Replaced boost::shared_array with Nektar::SharedArray and fixed several problems where the compile time array size was being used instead of the run time array size.

    Revision 1.7  2007/03/21 16:13:19  sherwin
    Fixed double to NekDouble casting in Array<OneD, NekDouble>

    Revision 1.6  2007/03/20 16:58:41  sherwin
    Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D

    Revision 1.5  2007/01/29 01:27:49  bnelson
    Added additional Allocate methods which take more parameters for the constructor.

    Removed the requirement for the user to specify, for each class, whether it can use the memory manager or not.  The memory manager is now a cmake level variable.

    Revision 1.4  2006/10/30 05:09:36  bnelson
    Fixed an error deleting object in an array that were not being constructed.

    Revision 1.3  2006/08/25 01:30:36  bnelson
    Added allocation of raw arrays.

    Revision 1.2  2006/06/01 13:44:29  kirby
    *** empty log message ***

    Revision 1.1  2006/06/01 09:17:24  kirby
    *** empty log message ***

    Revision 1.6  2006/05/30 14:00:03  sherwin
    Updates to make MultiRegions and its Demos work

    Revision 1.5  2006/05/18 04:23:57  bnelson
    Added allocation functions that pass arguments to the constructors of the objects being created.

    Revision 1.4  2006/05/15 04:13:36  bnelson
    no message

    Revision 1.3  2006/05/14 21:31:49  bnelson
    Modified the upper bound on static shared array allocation.

    Revision 1.2  2006/05/06 20:36:16  sherwin
    Modifications to get LocalRegions/Project1D working

    Revision 1.1  2006/05/04 18:57:43  kirby
    *** empty log message ***

    Revision 1.5  2006/04/25 20:17:39  jfrazier
    Fixed a .Net issue with the deallocator function passed to shared_array.

    Revision 1.4  2006/04/01 22:00:11  sherwin
    Changed definition of ASSERT

    Revision 1.3  2006/03/21 09:21:31  sherwin
    Introduced NekMemoryManager

    Revision 1.2  2006/02/26 21:11:40  bnelson
    Added a variable sized array allocator.

    Revision 1.1  2006/02/23 07:53:23  bnelson
    *** empty log message ***

**/
