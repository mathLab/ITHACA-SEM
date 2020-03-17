///////////////////////////////////////////////////////////////////////////////
//
// File ArrayPolicies.hpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Scientific Computing and Imaging Institute,
// University of Utah (USA) and Department of Aeronautics, Imperial
// College London (UK).
//
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
// Description: Implementations of the array creation and destruction policies.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_ARRAY_POLICIES_HPP
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_ARRAY_POLICIES_HPP

#include <type_traits>
#include <memory>

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/multi_array.hpp>

namespace Nektar
{
    template<typename ObjectType, typename enabled = void>
    class ArrayInitializationPolicy;
    
    /// \internal
    /// \brief Does nothing.
    template<typename ObjectType>
    class ArrayInitializationPolicy<
        ObjectType,
        typename std::enable_if<std::is_fundamental<ObjectType>::value>::type >
    {
        public:
            static void Initialize(ObjectType* data, size_t itemsToCreate)
            {
                boost::ignore_unused(data, itemsToCreate);
            }
            
            static void Initialize(ObjectType* data, size_t itemsToCreate, const ObjectType& initValue)
            {
                std::fill_n(data, itemsToCreate, initValue);
            }

            static void Initialize(ObjectType* data, size_t itemsToCreate, const ObjectType* initValue)
            {
                std::copy(initValue, initValue + itemsToCreate, data);
            }
    };
    
    /// \internal
    /// \brief Default initializes all elements.
    ///
    /// \internal 
    /// The user calls Initialize, which immediately calls DoInitialization.  There reason
    /// for this separation is because the code for creating an array of default initialized 
    /// elements and an array of copy constructed elements is identical except for the actual 
    /// call to placement new.
    ///
    /// 
    template<typename ObjectType>
    class ArrayInitializationPolicy<
        ObjectType,
        typename std::enable_if<!std::is_fundamental<ObjectType>::value>::type>
    {
        public:
            /// \brief Initalize each element in the array with ObjectType's default constructor.
            /// \param data The array of values to populate.
            /// \param itemsToCreate The size of data.
            static void Initialize(ObjectType* data, size_t itemsToCreate)
            {
                DoInitialization(
                    data, itemsToCreate,
                    [](ObjectType *element) { new (element) ObjectType; });
            }
            
            /// \brief Initalize each element in the array with ObjectType's copy constructor.
            /// \param data The array of values to populate.
            /// \param itemsToCreate The size of data.
            /// \param initValue The inital value each element in data will have.
            static void Initialize(ObjectType* data, size_t itemsToCreate, const ObjectType& initValue)
            {
                DoInitialization(
                    data, itemsToCreate,
                    [&](ObjectType *element) { new (element) ObjectType(initValue); });
            }

            static void Initialize(ObjectType* data, size_t itemsToCreate, const ObjectType* initValue)
            {
                DoInitialization(
                    data, itemsToCreate,
                    [&](ObjectType *element) { new (element) ObjectType(*initValue); initValue++; });
            }
            
            private:
                template<typename CreateType>
                static void DoInitialization(ObjectType* data, size_t itemsToCreate, const CreateType &f)
                {
                    size_t nextObjectToCreate = 0;
                    try
                    {
                        for(size_t i = 0; i < itemsToCreate; ++i)
                        {
                            ObjectType* memLocation = &data[i];
                            f(memLocation);
                            ++nextObjectToCreate;
                        }
                    }
                    catch(...)
                    {
                        for(size_t i = 0; i < nextObjectToCreate; ++i)
                        {
                            ObjectType* memLocation = &data[nextObjectToCreate - i - 1];
                            memLocation->~ObjectType();
                        }
                        throw;
                    }
                }
        };
    
    
    template<typename ObjectType, typename enabled = void>
    class ArrayDestructionPolicy;
    
    template<typename ObjectType>
    class ArrayDestructionPolicy<ObjectType,
                                 typename std::enable_if<std::is_fundamental<ObjectType>::value>::type>
    {
        public:
            static void Destroy(ObjectType* data, size_t itemsToDestroy)
            {
                boost::ignore_unused(data, itemsToDestroy);
            }
    };
    
    template<typename ObjectType>
    class ArrayDestructionPolicy<ObjectType,
                                 typename std::enable_if<!std::is_fundamental<ObjectType>::value>::type>
    {
        public:
            static void Destroy(ObjectType* data, size_t itemsToDestroy)
            {
                for(size_t i = 0; i < itemsToDestroy; ++i)
                {
                    ObjectType* memLocation = &data[itemsToDestroy - i - 1];
                    memLocation->~ObjectType();
                }
            }
    };

    template<typename Dim, typename DataType, typename ExtentListType>
    std::shared_ptr<boost::multi_array_ref<DataType, Dim::Value> > 
    CreateStorage(const ExtentListType& extent)
    {
        typedef boost::multi_array_ref<DataType, Dim::Value> ArrayType;
        size_t size = std::accumulate(extent.begin(), extent.end(), 1, 
            std::multiplies<size_t>());
        DataType* storage = MemoryManager<DataType>::RawAllocate(size);
        return MemoryManager<ArrayType>::AllocateSharedPtrD(
            [=](boost::multi_array_ref<DataType, Dim::Value> *ptr) {
                boost::ignore_unused(ptr);
                ArrayDestructionPolicy<DataType>::Destroy(storage, size);
                MemoryManager<DataType>::RawDeallocate(storage, size);
            },
            storage, extent);
    }
    
    template<typename DataType>
    std::shared_ptr<boost::multi_array_ref<DataType, 1> >
    CreateStorage(size_t d1)
    {
        std::vector<size_t> extents = { d1 };
        return CreateStorage<OneD, DataType>(extents);
    } 
    
    template<typename DataType>
    std::shared_ptr<boost::multi_array_ref<DataType, 2> >
    CreateStorage(size_t d1, size_t d2)
    {
        std::vector<size_t> extents = { d1, d2 };
        return CreateStorage<TwoD, DataType>(extents);
    }
    
    template<typename DataType>
    std::shared_ptr<boost::multi_array_ref<DataType, 3> >
    CreateStorage(size_t d1, size_t d2, size_t d3)
    {
        std::vector<size_t> extents = { d1, d2, d3 };
        return CreateStorage<ThreeD, DataType>(extents);
    }
}

#endif //NEKTAR_LIB_UTILITIES_BASIC_UTILS_ARRAY_POLICIES_HPP
