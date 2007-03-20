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
#include <boost/shared_array.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/bind.hpp>
#include <boost/type_traits.hpp>

#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>  
#include <boost/preprocessor/repetition/enum_trailing.hpp>  
#include <boost/preprocessor/repetition/enum_binary_params.hpp>
#include <boost/preprocessor/repetition/enum_trailing_params.hpp>

#include <vector>
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
    class MemoryManager
    {

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
            template<typename DataType>
            static void Deallocate(DataType*& data)
            {
#ifdef NEKTAR_MEMORY_POOL_ENABLED
                data->~DataType();
                ThreadSpecificPool<sizeof(DataType)>::Deallocate(data);
#else //NEKTAR_MEMORY_POOL_ENABLED
                delete data;
#endif //NEKTAR_MEMORY_POOL_ENABLED
                data = NULL;
            }

           
#ifdef NEKTAR_MEMORY_POOL_ENABLED
            #define ALLOCATE_METHOD_GENERATOR(z, i, methodName) \
            /* \brief test1 */ \
            template<typename DataType BOOST_PP_ENUM_TRAILING_PARAMS(i, typename Arg)> \
            static DataType* methodName(BOOST_PP_ENUM_BINARY_PARAMS(i, const Arg, & arg)) \
            { \
                DataType* result = static_cast<DataType*>(ThreadSpecificPool<sizeof(DataType)>::Allocate()); \
                \
                if( result ) \
                { \
                    try \
                    { \
                        new (result) DataType(BOOST_PP_ENUM_PARAMS(i, arg)); \
                    } \
                    catch(...) \
                    { \
                        ThreadSpecificPool<sizeof(DataType)>::Deallocate(result); \
                        throw; \
                    } \
                } \
                \
                return result; \
            }
#else //NEKTAR_MEMORY_POOL_ENABLED
            #define ALLOCATE_METHOD_GENERATOR(z, i, methodName) \
            template<typename DataType BOOST_PP_ENUM_TRAILING_PARAMS(i, typename Arg)> \
            static DataType* methodName(BOOST_PP_ENUM_BINARY_PARAMS(i, const Arg, & arg)) \
            { \
                return new DataType(BOOST_PP_ENUM_PARAMS(i, arg)); \
            }
#endif //NEKTAR_MEMORY_POOL_ENABLED

            BOOST_PP_REPEAT_FROM_TO(0, NEKTAR_MAX_MEMORY_MANAGER_CONSTRUCTOR_ARGS, ALLOCATE_METHOD_GENERATOR, Allocate)

            #define ALLOCATE_SHARED_PTR_METHOD_GENERATOR(z, i, methodName) \
            template<typename DataType BOOST_PP_ENUM_TRAILING_PARAMS(i, typename Arg)> \
            static boost::shared_ptr<DataType> methodName(BOOST_PP_ENUM_BINARY_PARAMS(i, const Arg, & arg)) \
            { \
                DataType* data = Allocate<DataType>(BOOST_PP_ENUM_PARAMS(i, arg)); \
                return boost::shared_ptr<DataType>(data, &MemoryManager::DeallocateSharedPtr<DataType>); \
            }

            BOOST_PP_REPEAT_FROM_TO(0, NEKTAR_MAX_MEMORY_MANAGER_CONSTRUCTOR_ARGS, ALLOCATE_SHARED_PTR_METHOD_GENERATOR, AllocateSharedPtr)


            /// \brief Deallocate an array allocated via MemoryManager::AllocateArray.
            template<unsigned int ArraySize, typename DataType>
            static void DeallocateArray(DataType*& data)
            {
                DeallocateArray<ArraySize, DataType>(data, ArraySize);
            }

            /// \brief Allocate an array.
            ///
            /// Every pointer obtained through this method must be deallocated with 
            /// MemoryManager::DeallocateArray.
            template<unsigned int ArraySize, typename DataType>
            static DataType* AllocateArray(unsigned int itemsToCreate = ArraySize,
                                           typename boost::disable_if<boost::is_fundamental<DataType> >::type* p = NULL)
            {
                BOOST_STATIC_ASSERT(ArraySize > 0);
#ifdef NEKTAR_MEMORY_POOL_ENABLED
                DataType* result = static_cast<DataType*>(ThreadSpecificPool<sizeof(DataType)*ArraySize>::Allocate());

                if( result )
                {
                    unsigned int nextObjectToCreate = 0;
                    try
                    {
                        for(unsigned int i = 0; i < itemsToCreate; ++i)
                        {
                            DataType* memLocation = &result[i];
                            new (memLocation) DataType;
                            ++nextObjectToCreate;
                        }
                    }
                    catch(...)
                    {
                        for(unsigned int i = 0; i < nextObjectToCreate; ++i)
                        {
                            DataType* memLocation = &result[i];
                            memLocation->~DataType();
                        }

                        ThreadSpecificPool<sizeof(DataType)>::Deallocate(result);
                        throw;
                    }
                }

                return result;
#else //NEKTAR_MEMORY_POOL_ENABLED
                return new DataType[itemsToCreate];
#endif //NEKTAR_MEMORY_POOL_ENABLED
            }
           
            /// \brief Creates an array of fundamental type using the memory pool.
            template<unsigned int ArraySize, typename DataType>
            static DataType* AllocateArray(unsigned int itemsToCreate = ArraySize,
                                           typename boost::enable_if<boost::is_fundamental<DataType> >::type* p = NULL)
            {
                BOOST_STATIC_ASSERT(ArraySize > 0);
#ifdef NEKTAR_MEMORY_POOL_ENABLED
                return static_cast<DataType*>(ThreadSpecificPool<sizeof(DataType)*ArraySize>::Allocate());
#else //NEKTAR_MEMORY_POOL_ENABLED
                return new DataType[itemsToCreate];
#endif //NEKTAR_MEMORY_POOL_ENABLED
            }

            /// \brief Creates an array.
            template<typename DataType>
            static DataType* AllocateArray(unsigned int arraySize)
            {
                if( arraySize < 10 )
                {
                    return MemoryManager::AllocateArray<10, DataType>(arraySize);
                }
                else if( arraySize < 20 )
                {
                    return MemoryManager::AllocateArray<20, DataType>(arraySize);
                }
                else if( arraySize < 30 )
                {
                    return MemoryManager::AllocateArray<30, DataType>(arraySize);
                }
                else if( arraySize < 40 )
                {
                    return MemoryManager::AllocateArray<40, DataType>(arraySize);
                }
                else if( arraySize < 50 )
                {
                    return MemoryManager::AllocateArray<50, DataType>(arraySize);
                }
                else if( arraySize < 60 )
                {
                    return MemoryManager::AllocateArray<60, DataType>(arraySize);
                }
                else
                {
                    return new DataType[arraySize];
                }
            }

            /// \brief Deallocates an array.
            template<typename DataType>
            static void DeallocateArray(DataType*& array, unsigned int arraySize)
            {
                if( arraySize < 10 )
                {
                    MemoryManager::DeallocateArray<10, DataType>(array, arraySize);
                }
                else if( arraySize < 20 )
                {
                    MemoryManager::DeallocateArray<20, DataType>(array, arraySize);
                }
                else if( arraySize < 30 )
                {
                    MemoryManager::DeallocateArray<30, DataType>(array, arraySize);
                }
                else if( arraySize < 40 )
                {
                    MemoryManager::DeallocateArray<40, DataType>(array, arraySize);
                }
                else if( arraySize < 50 )
                {
                    MemoryManager::DeallocateArray<50, DataType>(array, arraySize);
                }
                else if( arraySize < 60 )
                {
                    MemoryManager::DeallocateArray<60, DataType>(array, arraySize);
                }
                else
                {
                    delete[] array;
                }

                array = 0;
            }


            /// \brief Allocate a shared array of given size and type.
            template<unsigned int ArraySize, typename DataType>
            static boost::shared_array<DataType> AllocateSharedArray()
            {
                return AllocateSharedArray<ArraySize, DataType>(ArraySize);
            }

            template<typename DataType>
            static boost::shared_array<DataType> AllocateSharedArray(unsigned int arraySize)
            {
                if( arraySize < 10 )
                {
                    return MemoryManager::AllocateSharedArray<10, DataType>(arraySize);
                }
                else if( arraySize < 20 )
                {
                    return MemoryManager::AllocateSharedArray<20, DataType>(arraySize);
                }
                else if( arraySize < 30 )
                {
                    return MemoryManager::AllocateSharedArray<30, DataType>(arraySize);
                }
                else if( arraySize < 40 )
                {
                    return MemoryManager::AllocateSharedArray<40, DataType>(arraySize);
                }
                else if( arraySize < 50 )
                {
                    return MemoryManager::AllocateSharedArray<50, DataType>(arraySize);
                }
                else if( arraySize < 60 )
                {
                    return MemoryManager::AllocateSharedArray<60, DataType>(arraySize);
                }
                else
                {
                    return boost::shared_array<DataType>(new DataType[arraySize]);
                }
            }

        private:
            /// \brief Custom deletion policy for shared pointers.
            template<typename DataType>
            static void DeallocateSharedPtr(DataType* data)
            {
                Deallocate<DataType>(data);
            }

            /// \brief Custom deletion policy for shared arrays.
            template<unsigned int ArraySize, typename DataType>
            static void DeallocateSharedArray(DataType* data, unsigned int arraySize)
            {
                BOOST_STATIC_ASSERT(ArraySize > 0);
                DeallocateArray<ArraySize, DataType>(data, arraySize);
            }

            template<unsigned int ArraySize, typename DataType>
            static boost::shared_array<DataType> AllocateSharedArray(unsigned int arraySize)
            {
                DataType* data = AllocateArray<ArraySize, DataType>(arraySize);
                //void (*deallocator)(DataType *) = &MemoryManager::DeallocateSharedArray<ArraySize, DataType>;
                //return boost::shared_array<DataType>(data, deallocator);
                return boost::shared_array<DataType>(data, boost::bind(MemoryManager::DeallocateSharedArray<ArraySize, DataType>, _1, arraySize));
            }

            /// \brief Deallocates arrays of fundamental data types.
            template<unsigned int ArraySize, typename DataType>
            static void DeallocateArray(DataType*& data, unsigned int numToDelete,
                    typename boost::enable_if<boost::is_fundamental<DataType> >::type* p = NULL)
            {
                BOOST_STATIC_ASSERT(ArraySize > 0);
#ifdef NEKTAR_MEMORY_POOL_ENABLED
                ThreadSpecificPool<sizeof(DataType)*ArraySize>::Deallocate(data);
#else //NEKTAR_MEMORY_POOL_ENABLED
                delete [] data;
#endif //NEKTAR_MEMORY_POOL_ENABLED
                data = NULL;
            }

            /// \brief Deallocate a pointer allocated by the Memory Manager with a memory pool.
            /// \param numToDelete The number of objects in this memory chunk on which the destructor must be called.
            template<unsigned int ArraySize, typename DataType>
            static void DeallocateArray(DataType*& data, unsigned int numToDelete,
                typename boost::disable_if<boost::is_fundamental<DataType> >::type* p = NULL)
            {
#ifdef NEKTAR_MEMORY_POOL_ENABLED
                BOOST_STATIC_ASSERT(ArraySize > 0);
                for(unsigned int i = 0; i < numToDelete; ++i)
                {
                    DataType* memLocation = &data[i];
                    memLocation->~DataType();
                }

                ThreadSpecificPool<sizeof(DataType)*sizeof(ArraySize)>::Deallocate(data);
#else //NEKTAR_MEMORY_POOL_ENABLED
                delete [] data;
#endif //NEKTAR_MEMORY_POOL_ENABLED
                data = NULL;
            }

    };

    inline NekDoubleSharedArray  GetDoubleTmpSpace(const int size)
    {
	return MemoryManager::AllocateSharedArray<double>(size);
    }
    
    inline NekIntSharedArray  GetIntTmpSpace(const int size)
    {
	return MemoryManager::AllocateSharedArray<int>(size);
    }
    
}

#endif //NEKTAR_LIB_UTILITIES_NEK_MEMORY_MANAGER_H


/**
    $Log: NekMemoryManager.hpp,v $
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
