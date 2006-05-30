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

#include <LibUtilities/ThreadSpecificPool.hpp>
#include <LibUtilities/ErrorUtil.hpp>

#include <boost/mpl/contains.hpp>
#include <boost/mpl/list_c.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/bind.hpp>
#include <boost/type_traits.hpp>

namespace Nektar
{
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
    /// If you wish to use the memory pool with user defined data types, then
    /// you must create a static const MemoryPoolEnabler data member in your
    /// class and set it to eEnabled.  Once this is done all uses of the
    /// MemoryManager with your class will use the memory pool.  If you do not
    /// wish to use the memory manager then set this data member to eDisabled.
    ///
    /// For example:
    /// \code
    ///
    /// // This class doesn't use the pool.
    /// class Disabled
    /// {
    ///     public:
    ///         Disabled() {}
    ///
    ///         static const MemoryManager::MemoryPoolEnabler MemoryPoolEnabled = MemoryManager::eDisabled;
    /// };
    ///
    /// // This class uses the pool.
    /// class Enabled
    /// {
    ///     public:
    ///         Enabled() {}
    ///
    ///         static const MemoryManager::MemoryPoolEnabler MemoryPoolEnabled = MemoryManager::eEnabled;
    /// };
    /// \endcode
    ///
    ///
    class MemoryManager
    {

        public:

            /// \brief Used by clients to mark classes that should use memory pools.
            enum MemoryPoolEnabler
            {
                eDisabled,
                eEnabled
            };

            ////////////////////////////////////////////////////////////////////
            /// \name Object Allocation
            /// Routines to allocate and deallocate individual objects.
            /// \{
            ////////////////////////////////////////////////////////////////////

            /// \brief Deallocate a pointer allocated by the
            /// MemoryManager without a memory pool.
            template<typename DataType>
            static void Deallocate(DataType*& data,
                    typename boost::disable_if_c<DataType::MemoryPoolEnabled>::type* p = NULL)
            {
                delete data;
                data = NULL;
            }

            /// \brief Deallocate a pointer allocated by the Memory
            /// Manager with a memory pool.
            template<typename DataType>
            static void Deallocate(DataType*& data,
                    typename boost::enable_if_c<DataType::MemoryPoolEnabled>::type* p = NULL)
            {
                data->~DataType();
                ThreadSpecificPool<sizeof(DataType)>::Deallocate(data);
                data = NULL;
            }

            /// \brief Allocate a pointer to an object of type
            /// DataType*.
            ///
            /// This method based heavily on the boost object pool code.
            template<typename DataType>
            static DataType* Allocate(
                typename boost::enable_if_c<DataType::MemoryPoolEnabled>::type* p = NULL)
            {
                DataType* result = static_cast<DataType*>(ThreadSpecificPool<sizeof(DataType)>::Allocate());

                if( result )
                {
                    try
                    {
                        // This is placement new, constructing the
                        // object inside the memory returned by the
                        // pool.
                        new (result) DataType;
                    }
                    catch(...)
                    {
                        // Clean up the memory since the object didn't
                        // get created successfully.  Note that we
                        // don't call the destructor here because the
                        // object wasn't fully created.
                        ThreadSpecificPool<sizeof(DataType)>::Deallocate(result);
                        throw;
                    }
                }

                return result;
            }

            template<typename DataType, typename Arg1Type>
            static DataType* Allocate(const Arg1Type& arg1,
                typename boost::enable_if_c<DataType::MemoryPoolEnabled>::type* p = NULL)
            {
                DataType* result = static_cast<DataType*>(ThreadSpecificPool<sizeof(DataType)>::Allocate());

                if( result )
                {
                    try
                    {
                        // This is placement new, constructing the
                        // object inside the memory returned by the
                        // pool.
                        new (result) DataType(arg1);
                    }
                    catch(...)
                    {
                        // Clean up the memory since the object didn't
                        // get created successfully.  Note that we
                        // don't call the destructor here because the
                        // object wasn't fully created.
                        ThreadSpecificPool<sizeof(DataType)>::Deallocate(result);
                        throw;
                    }
                }

                return result;
            }

            template<typename DataType, typename Arg1Type, typename Arg2Type>
            static DataType* Allocate(const Arg1Type& arg1, const Arg2Type& arg2,
                typename boost::enable_if_c<DataType::MemoryPoolEnabled>::type* p = NULL)
            {
                DataType* result = static_cast<DataType*>(ThreadSpecificPool<sizeof(DataType)>::Allocate());

                if( result )
                {
                    try
                    {
                        // This is placement new, constructing the
                        // object inside the memory returned by the
                        // pool.
                        new (result) DataType(arg1, arg2);
                    }
                    catch(...)
                    {
                        // Clean up the memory since the object didn't
                        // get created successfully.  Note that we
                        // don't call the destructor here because the
                        // object wasn't fully created.
                        ThreadSpecificPool<sizeof(DataType)>::Deallocate(result);
                        throw;
                    }
                }

                return result;
            }

            /// \brief Allocate a pointer without using the memory
            /// pool.
            template<typename DataType>
            static DataType* Allocate(
                typename boost::disable_if_c<DataType::MemoryPoolEnabled>::type* p = NULL)
            {
                return new DataType();
            }


            ////////////////////////////////////////////////////////////////////
            /// \}
            ////////////////////////////////////////////////////////////////////

            ////////////////////////////////////////////////////////////////////
            /// \name Shared Pointer Allocation
            /// Routines to allocate a single object in a shared pointer.  In
            /// all cases the shared pointer will correctly release memory,
            /// regardless of how the memory was allocated.
            /// \{
            ////////////////////////////////////////////////////////////////////

            /// \brief Custom deletion policy for shared pointers.
            template<typename DataType>
            static void DeallocateSharedPtr(DataType* data)
            {
                Deallocate<DataType>(data);
            }

            /// \brief Allocate an object in a shared pointer.
            template<typename DataType>
            static boost::shared_ptr<DataType> AllocateSharedPtr()
            {
                DataType* data = Allocate<DataType>();
                return boost::shared_ptr<DataType>(data, &MemoryManager::DeallocateSharedPtr<DataType>);
            }

            template<typename DataType, typename Arg1Type>
                    static boost::shared_ptr<DataType> AllocateSharedPtr(const Arg1Type& arg1)
            {
                DataType* data = Allocate<DataType>(arg1);
                return boost::shared_ptr<DataType>(data, &MemoryManager::DeallocateSharedPtr<DataType>);
            }

            template<typename DataType, typename Arg1Type, typename Arg2Type>
                    static boost::shared_ptr<DataType> AllocateSharedPtr(const Arg1Type& arg1,
                                                                        const Arg2Type& arg2)
            {
                DataType* data = Allocate<DataType>(arg1, arg2);
                return boost::shared_ptr<DataType>(data, &MemoryManager::DeallocateSharedPtr<DataType>);
            }


            ////////////////////////////////////////////////////////////////////
            /// \}
            ////////////////////////////////////////////////////////////////////


            ////////////////////////////////////////////////////////////////////
            /// \name Array Allocation
            /// Allocates an array of user specified data.
            ///
            /// Care must be taken when deleting arrays.  The ArraySize template
            /// parameter must match the ArraySize parameter used when allocating
            /// the array.  A mismatch will cause undefined behavior.  For this
            /// reason, the use of the shared_array versions is preferred.
            /// \{
            ////////////////////////////////////////////////////////////////////

            /// \brief Deallocate a pointer allocated by the MemoryManager without a memory pool.
            /// \param data The data to be deleted.  This parameter will be set to
            ///             NULL when the method is finished.
            template<unsigned int ArraySize, typename DataType>
            static void DeallocateArray(DataType*& data,
                    typename boost::disable_if_c<DataType::MemoryPoolEnabled>::type* p = NULL)
            {
                BOOST_STATIC_ASSERT(ArraySize > 0);
                delete [] data;
                data = NULL;
            }

            /// \brief Deallocates arrays of fundamental data types.
            template<unsigned int ArraySize, typename DataType>
            static void DeallocateArray(DataType*& data,
                    typename boost::enable_if<boost::is_fundamental<DataType> >::type* p = NULL)
            {
                BOOST_STATIC_ASSERT(ArraySize > 0);
                ThreadSpecificPool<sizeof(DataType)*ArraySize>::Deallocate(data);
                data = NULL;
            }

            /// \brief Deallocate a pointer allocated by the Memory Manager with a memory pool.
            template<unsigned int ArraySize, typename DataType>
            static void DeallocateArray(DataType*& data,
                    typename boost::enable_if_c<DataType::MemoryPoolEnabled>::type* p = NULL)
            {
                BOOST_STATIC_ASSERT(ArraySize > 0);
                for(unsigned int i = 0; i < ArraySize; ++i)
                {
                    DataType* memLocation = &data[i];
                    memLocation->~DataType();
                }

                // Dangerous if I pass in the wrong array size.
                ThreadSpecificPool<sizeof(DataType)*sizeof(ArraySize)>::Deallocate(data);
                data = NULL;
            }

            /// \brief Allocate ArraySize values in an array of
            /// DataType using the memory pool.
            template<unsigned int ArraySize, typename DataType>
            static DataType* AllocateArray(
                typename boost::enable_if_c<DataType::MemoryPoolEnabled>::type* p = NULL)
            {
                BOOST_STATIC_ASSERT(ArraySize > 0);
                DataType* result = static_cast<DataType*>(ThreadSpecificPool<sizeof(DataType)*ArraySize>::Allocate());

                if( result )
                {
                    unsigned int nextObjectToCreate = 0;
                    try
                    {
                        for(unsigned int i = 0; i < ArraySize; ++i)
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
            }

            /// \brief Create an array without using the memory pool.
            template<unsigned int ArraySize, typename DataType>
            static DataType* AllocateArray(
                typename boost::disable_if_c<DataType::MemoryPoolEnabled>::type* p = NULL)
            {
                BOOST_STATIC_ASSERT(ArraySize > 0);
                return new DataType[ArraySize];
            }

            /// \brief Creates an array of fundamental type using the memory pool.
            template<unsigned int ArraySize, typename DataType>
            static DataType* AllocateArray(
                typename boost::enable_if<boost::is_fundamental<DataType> >::type* p = NULL)
            {
                BOOST_STATIC_ASSERT(ArraySize > 0);
                return static_cast<DataType*>(ThreadSpecificPool<sizeof(DataType)*ArraySize>::Allocate());
            }

            ////////////////////////////////////////////////////////////////////
            /// \}
            ////////////////////////////////////////////////////////////////////

            ////////////////////////////////////////////////////////////////////
            /// \name Shared Array Allocation
            /// \{
            ////////////////////////////////////////////////////////////////////

            /// \brief Custom deletion policy for shared arrays.
            template<unsigned int ArraySize, typename DataType>
            static void DeallocateSharedArray(DataType* data)
            {
                BOOST_STATIC_ASSERT(ArraySize > 0);
                DeallocateArray<ArraySize, DataType>(data);
            }

            /// \brief Allocate a shared array of given size and type.
            template<unsigned int ArraySize, typename DataType>
            static boost::shared_array<DataType> AllocateSharedArray()
            {
                DataType* data = AllocateArray<ArraySize, DataType>();
                void (*deallocator)(DataType *) = &MemoryManager::DeallocateSharedArray<ArraySize, DataType>;
                return boost::shared_array<DataType>(data, deallocator);
            }

            template<typename DataType>
            static boost::shared_array<DataType> AllocateSharedArray(unsigned int arraySize)
            {
                if( arraySize < 10 )
                {
                    return MemoryManager::AllocateSharedArray<10, DataType>();
                }
                else if( arraySize < 20 )
                {
                    return MemoryManager::AllocateSharedArray<20, DataType>();
                }
                else if( arraySize < 30 )
                {
                    return MemoryManager::AllocateSharedArray<30, DataType>();
                }
                else if( arraySize < 40 )
                {
		  return MemoryManager::AllocateSharedArray<40, DataType>();
                }
                else if( arraySize < 50 )
                {
                    return MemoryManager::AllocateSharedArray<50, DataType>();
                }
                else if( arraySize < 60 )
                {
                    return MemoryManager::AllocateSharedArray<60, DataType>();
                }
                else
                {
                    return boost::shared_array<DataType>(new DataType[arraySize]);
                }
            }

      ////////////////////////////////////////////////////////////////////
      /// \}
      ////////////////////////////////////////////////////////////////////
    };

  // typedef boost access for use in temporary memory allocation
  typedef boost::shared_array<double> BstShrDArray;
  typedef boost::shared_array<int>    BstShrIArray;

  inline BstShrDArray  GetDoubleTmpSpace(const int size)
  {
    return MemoryManager::AllocateSharedArray<double>(size);
  }

  inline BstShrIArray  GetIntTmpSpace(const int size)
  {
    return MemoryManager::AllocateSharedArray<int>(size);
  }
  
}

#endif //NEKTAR_LIB_UTILITIES_NEK_MEMORY_MANAGER_H


/**
    $Log: NekMemoryManager.hpp,v $
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
