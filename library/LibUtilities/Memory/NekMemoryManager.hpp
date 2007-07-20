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
//#include <LibUtilities/BasicUtils/SharedArray.hpp>
//#include <boost/shared_array.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

#include <boost/mpl/contains.hpp>
#include <boost/mpl/list_c.hpp>
#include <LibUtilities/BasicUtils/SharedPtr.hpp>
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
#define NEKTAR_MAX_MEMORY_MANAGER_CONSTRUCTOR_ARGS 9
#endif //NEKTAR_MAX_MEMORY_MANAGER_CONSTRUCTOR_ARGS

    /// \brief General purpose memory allocation routines with the ability
    ///        to allocate from thread specific memory pools.
    ///
    /// Allows the allocation of pointers, shared pointers, arrays, and
    /// shared arrays.  Each of these options can also be allocated with
    /// the default C++ new/delete or with thread specific memory pools.
    /// These memory pools alow multiple threads to use
    /// the same manager object without blocking while waiting for memory
    /// allocation from the system.
    ///
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

           
            /// \brief Deallocate a pointer allocated MemoryManager::Allocate.
            /// \note Results are undefined if called with a pointer to something that 
            ///       was not allocated with the memory manager.
            ///
            /// Use this method to deallocate a pointer you have allocated from the 
            /// MemoryManager.
            ///
            /// Example:
            /// \code
            /// CustObj* c = MemoryManager::Allocate<CustObj>();
            /// MemoryManager::Deallocate(c);
            /// \endcode
            static void Deallocate(DataType*& data)
            {
                #ifdef NEKTAR_MEMORY_POOL_ENABLED
                    data->~DataType();
                    //ThreadSpecificPool<sizeof(DataType)>::Deallocate(data);
                    MemPool<sizeof(DataType)>::Type::Instance().Deallocate(data);
                #else //NEKTAR_MEMORY_POOL_ENABLED
                    delete data;
                #endif //NEKTAR_MEMORY_POOL_ENABLED
                
                data = NULL;
            }

           
//            template<typename Arg1>
//             static DataType* Allocate(const Arg1& arg1) 
//             { 
//                 DataType* result = static_cast<DataType*>(ThreadSpecificPool<sizeof(DataType)>::Allocate()); 
//                 
//                 if( result ) 
//                 { 
//                     try 
//                     { 
//                         new (result) DataType(arg1);BOOST_PP_ENUM_PARAMS(i, arg)); \
//                     } \
//                     catch(...) \
//                     { \
//                         ThreadSpecificPool<sizeof(DataType)>::Deallocate(result); \
//                         throw; \
//                     } \
//                 } \
//                 \
//                 return result; \
//             }
                
            #ifdef NEKTAR_MEMORY_POOL_ENABLED
                static DataType* Allocate() 
                { 
                    //DataType* result = static_cast<DataType*>(ThreadSpecificPool<sizeof(DataType)>::Allocate()); 
                    DataType* result = static_cast<DataType*>(MemPool<sizeof(DataType)>::Type::Instance().Allocate());
                    
                    if( result ) 
                    { 
                        try 
                        { 
                            new (result) DataType();
                        } 
                        catch(...) 
                        { 
                            //ThreadSpecificPool<sizeof(DataType)>::Deallocate(result); 
                            MemPool<sizeof(DataType)>::Type::Instance().Deallocate(result);
                            throw; 
                        } 
                    } 
                    
                    return result; 
                }
                #define ALLOCATE_METHOD_GENERATOR(z, i, methodName) \
                /* \brief test1 */ \
                template<BOOST_PP_ENUM_PARAMS(i, typename Arg)> \
                static DataType* methodName(BOOST_PP_ENUM_BINARY_PARAMS(i, const Arg, & arg)) \
                { \
                    DataType* result = static_cast<DataType*>(MemPool<sizeof(DataType)>::Type::Instance().Allocate()); \
                    \
                    if( result ) \
                    { \
                        try \
                        { \
                            new (result) DataType(BOOST_PP_ENUM_PARAMS(i, arg)); \
                        } \
                        catch(...) \
                        { \
                            MemPool<sizeof(DataType)>::Type::Instance().Deallocate(result); \
                            throw; \
                        } \
                    } \
                    \
                    return result; \
                }
            #else //NEKTAR_MEMORY_POOL_ENABLED
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

            static ptr<DataType> AllocateSharedPtr()
            {
                return AllocateSharedPtrD(DefaultCustomDeallocator());
                //DataType* data = Allocate(); 
                //return ptr<DataType>(data, &MemoryManager::DeallocateSharedPtr); 
            }
            
            #define ALLOCATE_SHARED_PTR_METHOD_GENERATOR(z, i, methodName) \
            template<BOOST_PP_ENUM_PARAMS(i, typename Arg)> \
            static ptr<DataType> methodName(BOOST_PP_ENUM_BINARY_PARAMS(i, const Arg, & arg)) \
            { \
                return AllocateSharedPtrD(DefaultCustomDeallocator(), BOOST_PP_ENUM_PARAMS(i, arg)); \
            }

            #define ALLOCATE_SHARED_PTR_METHOD_WITH_DEALLOCATOR_GENERATOR(z, i, methodName) \
            template<typename DeallocatorType BOOST_PP_ENUM_TRAILING_PARAMS(i, typename Arg)> \
            static ptr<DataType> methodName(const DeallocatorType& d BOOST_PP_ENUM_TRAILING_BINARY_PARAMS(i, const Arg, & arg)) \
            { \
                DataType* data = Allocate(BOOST_PP_ENUM_PARAMS(i, arg)); \
                return ptr<DataType>(data, DeallocateSharedPtr<DataType, DeallocatorType>(d)); \
            }
            
            BOOST_PP_REPEAT_FROM_TO(0, NEKTAR_MAX_MEMORY_MANAGER_CONSTRUCTOR_ARGS, ALLOCATE_SHARED_PTR_METHOD_WITH_DEALLOCATOR_GENERATOR, AllocateSharedPtrD);
            BOOST_PP_REPEAT_FROM_TO(1, NEKTAR_MAX_MEMORY_MANAGER_CONSTRUCTOR_ARGS, ALLOCATE_SHARED_PTR_METHOD_GENERATOR, AllocateSharedPtr);




            /// \brief Allocates a chunk of raw, uninitialized memory, capable of holding NumberOfElements objects.
            template<unsigned int NumberOfElements>
            static DataType* RawAllocate()
            {
                BOOST_STATIC_ASSERT(NumberOfElements > 0);
                
                #ifdef NEKTAR_MEMORY_POOL_ENABLED
                    return static_cast<DataType*>(MemPool<sizeof(DataType)*NumberOfElements>::Type::Instance().Allocate());
                #else //NEKTAR_MEMORY_POOL_ENABLED
                    return static_cast<DataType*>(::operator new(NumberOfElements * sizeof(DataType)));
                #endif //NEKTAR_MEMORY_POOL_ENABLED
            }
            
            /// \brief Allocates a chunk of raw, uninitialized memory, capable of holding NumberOfElements objects.
            static DataType* RawAllocate(unsigned int NumberOfElements)
            {
                // Even if 0 elements are requested, we should allocate some memory.
                if( NumberOfElements <= 1 )
                {
                    return RawAllocate<1>();
                }
                else if( NumberOfElements < 10 )
                {
                    return RawAllocate<10>();
                }
                else if( NumberOfElements < 20 )
                {
                    return RawAllocate<20>();
                }
                else if( NumberOfElements < 30 )
                {
                    return RawAllocate<30>();
                }
                else if( NumberOfElements < 40 )
                {
                    return RawAllocate<40>();
                }
                else if( NumberOfElements < 50 )
                {
                    return RawAllocate<50>();
                }
                else if( NumberOfElements < 60 )
                {
                    return RawAllocate<60>();
                }
                else
                {
                    return static_cast<DataType*>(::operator new(NumberOfElements * sizeof(DataType)));
                }
            }
            
            /// \brief Deallocates arrays of fundamental data types.
            template<unsigned int NumberOfElements>
            static void RawDeallocate(DataType* data)
            {
                BOOST_STATIC_ASSERT(NumberOfElements > 0);
                
                #ifdef NEKTAR_MEMORY_POOL_ENABLED
                    MemPool<sizeof(DataType)*NumberOfElements>::Type::Instance().Deallocate(data);
                #else //NEKTAR_MEMORY_POOL_ENABLED
                    ::operator delete(data);
                #endif //NEKTAR_MEMORY_POOL_ENABLED
            }
                    

            static void RawDeallocate(DataType* array, unsigned int NumberOfElements)
            {
                if( NumberOfElements <= 1 )
                {
                    RawDeallocate<1>(array);
                }
                else if( NumberOfElements < 10 )
                {
                    RawDeallocate<10>(array);
                }
                else if( NumberOfElements < 20 )
                {
                    RawDeallocate<20>(array);
                }
                else if( NumberOfElements < 30 )
                {
                    RawDeallocate<30>(array);
                }
                else if( NumberOfElements < 40 )
                {
                    RawDeallocate<40>(array);
                }
                else if( NumberOfElements < 50 )
                {
                    RawDeallocate<50>(array);
                }
                else if( NumberOfElements < 60 )
                {
                    RawDeallocate<60>(array);
                }
                else
                {
                    ::operator delete(array);
                }
            }
            
            /////////////////////////////////////////////////////////////////
            /// @[ Allocator Interface
            /// \todo I don't think this is the correct syntax for doxygen sections.
            /////////////////////////////////////////////////////////////////
            typedef DataType value_type;
            typedef size_t size_type;
            typedef ptrdiff_t difference_type;
            typedef DataType* pointer;
            typedef const DataType* const_pointer;
            typedef DataType& reference;
            typedef const DataType& const_reference;
            
            MemoryManager() {}
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
