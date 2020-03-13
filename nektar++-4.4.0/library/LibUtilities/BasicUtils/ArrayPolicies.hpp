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
// Description: Implementations of the array creation and destruction policies.
//
///////////////////////////////////////////////////////////////////////////////


#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_ARRAY_POLICIES_HPP
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_ARRAY_POLICIES_HPP

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>
#include <boost/bind.hpp>

namespace Nektar
{
    template<typename ObjectType, typename enabled = void>
    class ArrayInitializationPolicy;
    
    /// \internal
    /// \brief Does nothing.
    template<typename ObjectType>
    class ArrayInitializationPolicy<ObjectType, 
                                typename boost::enable_if<boost::is_fundamental<ObjectType> >::type >
    {
        public:
            static void Initialize(ObjectType* data, unsigned int itemsToCreate)
            {
            }
            
            static void Initialize(ObjectType* data, unsigned int itemsToCreate, const ObjectType& initValue)
            {
                std::fill_n(data, itemsToCreate, initValue);
            }

            static void Initialize(ObjectType* data, unsigned int itemsToCreate, const ObjectType* initValue)
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
    class ArrayInitializationPolicy<ObjectType, 
                                typename boost::disable_if<boost::is_fundamental<ObjectType> >::type >
    {
        public:
            /// \brief Initalize each element in the array with ObjectType's default constructor.
            /// \param data The array of values to populate.
            /// \param itemsToCreate The size of data.
            static void Initialize(ObjectType* data, unsigned int itemsToCreate)
            {
                DoInitialization(data, itemsToCreate, 
                                    boost::bind(&ArrayInitializationPolicy<ObjectType>::DefaultConstructionWithPlacementNew, _1));
            }
            
            /// \brief Initalize each element in the array with ObjectType's copy constructor.
            /// \param data The array of values to populate.
            /// \param itemsToCreate The size of data.
            /// \param initValue The inital value each element in data will have.
            static void Initialize(ObjectType* data, unsigned int itemsToCreate, const ObjectType& initValue)
            {
                DoInitialization(data, itemsToCreate, 
                                    boost::bind(&ArrayInitializationPolicy<ObjectType>::CopyConstructionWithPlacementNew, _1, boost::ref(initValue)));
            }

            static void Initialize(ObjectType* data, unsigned int itemsToCreate, const ObjectType* initValue)
            {
                DoInitialization(data, itemsToCreate, 
                        boost::bind(&ArrayInitializationPolicy<ObjectType>::CopyConstructionFromArray, _1, boost::ref(initValue)));
            }
            
            private:
                template<typename CreateType>
                static void DoInitialization(ObjectType* data, unsigned int itemsToCreate, const CreateType& f)
                {
                    unsigned int nextObjectToCreate = 0;
                    try
                    {
                        for(unsigned int i = 0; i < itemsToCreate; ++i)
                        {
                            ObjectType* memLocation = &data[i];
                            f(memLocation);
                            ++nextObjectToCreate;
                        }
                    }
                    catch(...)
                    {
                        for(unsigned int i = 0; i < nextObjectToCreate; ++i)
                        {
                            ObjectType* memLocation = &data[nextObjectToCreate - i - 1];
                            memLocation->~ObjectType();
                        }
                        throw;
                    }
                }
                
                static void DefaultConstructionWithPlacementNew(ObjectType* element)
                {
                    new (element) ObjectType;
                }
                
                static void CopyConstructionWithPlacementNew(ObjectType* element, const ObjectType& initValue)
                {
                    new (element) ObjectType(initValue);
                }

                static void CopyConstructionFromArray(ObjectType* element, const ObjectType*& rhs)
                {
                    new (element) ObjectType(*rhs);
                    rhs += 1;
                }
        };
    
    
    template<typename ObjectType, typename enabled = void>
    class ArrayDestructionPolicy;
    
    template<typename ObjectType>
    class ArrayDestructionPolicy<ObjectType, 
                            typename boost::enable_if<boost::is_fundamental<ObjectType> >::type >
    {
        public:
            static void Destroy(ObjectType* data, unsigned int itemsToDestroy)
            {
            }
    };
    
    template<typename ObjectType>
    class ArrayDestructionPolicy<ObjectType, 
                            typename boost::disable_if<boost::is_fundamental<ObjectType> >::type >
    {
        public:
            static void Destroy(ObjectType* data, unsigned int itemsToDestroy)
            {
                for(unsigned int i = 0; i < itemsToDestroy; ++i)
                {
                    ObjectType* memLocation = &data[itemsToDestroy - i - 1];
                    memLocation->~ObjectType();
                }
            }
    };
    
    template<typename DataType>
    void DeleteStorage(DataType* data, unsigned int num)
    {
        ArrayDestructionPolicy<DataType>::Destroy(data, num);
        MemoryManager<DataType>::RawDeallocate(data, num);
    }
    
    template<typename Dim, typename DataType, typename ExtentListType>
    boost::shared_ptr<boost::multi_array_ref<DataType, Dim::Value> > 
    CreateStorage(const ExtentListType& extent)
    {
        typedef boost::multi_array_ref<DataType, Dim::Value> ArrayType;
        unsigned int size = std::accumulate(extent.begin(), extent.end(), 1, 
            std::multiplies<unsigned int>());
        DataType* storage = MemoryManager<DataType>::RawAllocate(size);
        return MemoryManager<ArrayType>::AllocateSharedPtrD(
                boost::bind(DeleteStorage<DataType>, storage, size),
                storage, extent);
    }
    
    template<typename DataType>
    boost::shared_ptr<boost::multi_array_ref<DataType, 1> >
    CreateStorage(unsigned int d1)
    {
        std::vector<unsigned int> extents(1, d1);
        return CreateStorage<OneD, DataType>(extents);
    } 
    
    template<typename DataType>
    boost::shared_ptr<boost::multi_array_ref<DataType, 2> >
    CreateStorage(unsigned int d1, unsigned int d2)
    {
        unsigned int vals[]  = {d1, d2};
        std::vector<unsigned int> extents(vals, vals+2);
        return CreateStorage<TwoD, DataType>(extents);
    }
    
    template<typename DataType>
    boost::shared_ptr<boost::multi_array_ref<DataType, 3> >
    CreateStorage(unsigned int d1, unsigned int d2, unsigned int d3)
    {
        unsigned int vals[]  = {d1, d2, d3};
        std::vector<unsigned int> extents(vals, vals+3);
        return CreateStorage<ThreeD, DataType>(extents);
    }
}

#endif //NEKTAR_LIB_UTILITIES_BASIC_UTILS_ARRAY_POLICIES_HPP
