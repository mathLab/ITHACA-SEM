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

#include <LibUtilities/BasicUtils/mojo.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/multi_array.hpp>
#include <boost/shared_ptr.hpp>

namespace Nektar
{
    // Forward declaration for a ConstArray constructor.
    template<Dimension Dim, typename DataType>
    class Array;

    /// \brief An array of unchangeable elements.
    /// \param Dim The array's dimensionality.
    /// \param DataType The data type held at each position in the array.
    ///
    /// ConstArray is a thin wrapper around boost::multi_array.  It supports
    /// much of the multi_array interface (with a notable exception of views)
    /// including the arbitrary data type and dimension.
    ///
    /// ConstArray sets itself apart from multi_array by not allowing any 
    /// modification of the underlying data.  For modifiable arrays, \see Array.
    ///
    /// ConstArray is most often used as the return value of a method when 
    /// the user needs access to an array of data but should not be allowed 
    /// to modify it.
    ///
    ///\code
    /// class Sample
    /// {
    ///     public:
    ///         ConstArray<OneD, double>& getData() const { return m_data; }
    ///        
    ///     private:
    ///         Array<OneD, double> m_data;
    /// };
    ///\endcode
    ///
    /// In this example, each instance of Sample contains an array.  The getData method 
    /// gives the user access to the array values, but does not allow modification of those
    /// values.
    ///\section Usage Test
    /// 
    /// For those with a strong C background, the easiest way to look at an array is to 
    /// transform "Array<OneD, DataType>" to "DataType*".  Therfore, in any code where you 
    /// would have specified "DataType*" you will want to use "Array<OneD, DataType>", and 
    /// where you would have used 
    ///
    /// <TABLE>
    /// <TR>
    ///     <TH>Native C Parameter Type</TH>
    ///     <TH>Array Parameter Type</TH>
    ///     <TH>Description</TH>
    /// </TR>
    /// <TR>
    ///     <TD>DataType*</TD>
    ///     <TD>Array&lt;OneD, DataType&gt; </TD>
    ///     <TD></TD>
    /// </TR>
    /// <TR>
    ///     <TD>DataType*</TD>
    ///     <TD>Array&lt;OneD, DataType&gt; </TD>
    ///     <TD></TD>
    /// </TR>
    /// </TABLE>
    /// const Array<>& is kind of weird.  You can't change it directly, but you can create a copy
    /// of it and change that.  So the const can be ignored as needed.  Mostly useful for cases 
    /// where you want to be able to change the array, and you need to accept temporaries.
    ///
    ///\section Efficiency Efficiency Considerations
    ///
    ///\code
    /// // Instead of this:
    /// ConstArray<OneD, NekDouble> lessEfficient(points->GetZ());
    ///
    /// // Do this
    /// ConstArray<OneD, NekDouble>& moreEfficient = points->GetZ();
    ///\endcode
    ///
    /// If you know the size, don't create an empty array and then populate it on another line.
    template<Dimension Dim, typename DataType>
    class ConstArray
    {
        public:
            typedef boost::multi_array_ref<DataType, Dim> ArrayType;

            typedef typename ArrayType::const_reference const_reference;
            typedef typename ArrayType::index index;
            typedef typename ArrayType::const_iterator const_iterator;
            typedef typename ArrayType::element element;
            typedef typename ArrayType::size_type size_type;
            
        public:
            template<typename ObjectType, typename enabled = void>
            class InitializationPolicy;
            
            /// \brief Does nothing.
            template<typename ObjectType>
            class InitializationPolicy<ObjectType, 
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
            class InitializationPolicy<ObjectType, 
                                       typename boost::disable_if<boost::is_fundamental<ObjectType> >::type >
            {
                public:
                    /// \brief Initalize each element in the array with ObjectType's default constructor.
                    /// \param data The array of values to populate.
                    /// \param itemsToCreate The size of data.
                    static void Initialize(ObjectType* data, unsigned int itemsToCreate)
                    {
                        DoInitialization(data, itemsToCreate, 
                                         boost::bind(&InitializationPolicy<ObjectType>::DefaultConstructionWithPlacementNew, _1));
                    }
                    
                    /// \brief Initalize each element in the array with ObjectType's copy constructor.
                    /// \param data The array of values to populate.
                    /// \param itemsToCreate The size of data.
                    /// \param initValue The inital value each element in data will have.
                    static void Initialize(ObjectType* data, unsigned int itemsToCreate, const ObjectType& initValue)
                    {
                        DoInitialization(data, itemsToCreate, 
                                         boost::bind(&InitializationPolicy<ObjectType>::CopyConstructionWithPlacementNew, _1, boost::ref(initValue)));
                    }

                    static void Initialize(ObjectType* data, unsigned int itemsToCreate, const ObjectType* initValue)
                    {
                        DoInitialization(data, itemsToCreate, 
                                boost::bind(&InitializationPolicy<ObjectType>::CopyConstructionFromArray, _1, boost::ref(initValue)));
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
                                    memLocation->~DataType();
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
            class DestructionPolicy;
            
            template<typename ObjectType>
            class DestructionPolicy<ObjectType, 
                                    typename boost::enable_if<boost::is_fundamental<ObjectType> >::type >
            {
                public:
                    static void Destroy(ObjectType* data, unsigned int itemsToDestroy)
                    {
                    }
            };
            
            template<typename ObjectType>
            class DestructionPolicy<ObjectType, 
                                    typename boost::disable_if<boost::is_fundamental<ObjectType> >::type >
            {
                public:
                    static void Destroy(ObjectType* data, unsigned int itemsToDestroy)
                    {
                        for(unsigned int i = 0; i < itemsToDestroy; ++i)
                        {
                            DataType* memLocation = &data[itemsToDestroy - i - 1];
                            memLocation->~DataType();
                        }
                    }
            };
            
        public:
            /// \brief Constructs an empty array.
            ///
            /// An empty array will still have a small amount of memory allocated for it, typically around 
            /// 4 bytes.  
            /// \post size() == 0
            /// \post begin() == end()
            ConstArray() :
                m_data(),
                m_offset(0)
            {
                std::vector<unsigned int> extents(Dim, 0);
                CreateStorage(extents);
            }
            
            // We could have used the boost pre-processor library to automate the generation 
            // of the constructors that take varying parameters based on the dimension (similar 
            // to how it is done in the MemoryManager.  However, since we don't anticipate that 
            // anyone will create arrays with a dimension greater than 3, the code is clearer by
            // spelling it out explicitly.
            
            /// \brief Constructs a 1 dimensional array.  The elements of the array are not initialized.
            explicit ConstArray(unsigned int arraySize) :
                m_data(),
                m_offset(0)
            {
                BOOST_STATIC_ASSERT(Dim==1);
                std::vector<unsigned int> extents(Dim, arraySize);
                CreateStorage(extents);
                InitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements());
            }
            
            /// \brief Constructs a 2 dimensional array.  The elements of the array are not initialized.
            ConstArray(unsigned int dim1Size, unsigned int dim2Size) :
                m_data(),
                m_offset(0)
            {
                BOOST_STATIC_ASSERT(Dim==2);
                unsigned int vals[] = {dim1Size, dim2Size};
                std::vector<unsigned int> extents(vals, vals+2);
                CreateStorage(extents);
                InitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements());
            }
            
            /// \brief Constructs a 3 dimensional array.  The elements of the array are not initialized.
            ConstArray(unsigned int dim1Size, unsigned int dim2Size, unsigned int dim3Size) :
                m_data(),
                m_offset(0)
            {
                BOOST_STATIC_ASSERT(Dim==3);
                unsigned int vals[] = {dim1Size, dim2Size, dim3Size};
                std::vector<unsigned int> extents(vals, vals+3);
                CreateStorage(extents);
                InitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements());
            }
            
            /// \brief Constructs a 1D array.  Every element is given the value initValue.
            ConstArray(unsigned int arraySize, const DataType& initValue) :
                m_data(),
                m_offset(0)
            {
                BOOST_STATIC_ASSERT(Dim == 1);
                std::vector<unsigned int> extents(Dim, arraySize);
                CreateStorage(extents);
                InitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements(), initValue);
            }
            
            ConstArray(unsigned int dim1Size, unsigned int dim2Size, const DataType& initValue) :
                m_data(),
                m_offset(0)
            {
                BOOST_STATIC_ASSERT(Dim==2);
                unsigned int vals[] = {dim1Size, dim2Size};
                std::vector<unsigned int> extents(vals, vals+2);
                CreateStorage(extents);
                InitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements(), initValue);
            }
            
            ConstArray(unsigned int dim1Size, unsigned int dim2Size, unsigned int dim3Size, const DataType& initValue) :
                m_data(),
                m_offset(0)
            {
                BOOST_STATIC_ASSERT(Dim==3);
                unsigned int vals[] = {dim1Size, dim2Size, dim3Size};
                std::vector<unsigned int> extents(vals, vals+3);
                CreateStorage(extents);
                InitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements(), initValue);
            }
            
            ConstArray(unsigned int arraySize, const DataType* d) :
                m_data(),
                m_offset(0)
            {
                BOOST_STATIC_ASSERT(Dim==1);
                std::vector<unsigned int> extents(Dim, arraySize);
                CreateStorage(extents);
                InitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements(), d);
            }
             
            ConstArray(const ConstArray<Dim, DataType>& rhs) :
                m_data(rhs.m_data),
                m_offset(rhs.m_offset)
            {
            }
             
            ConstArray(const ConstArray<Dim, DataType>& rhs, unsigned int offset) :
                m_data(rhs.m_data),
                m_offset(offset)
            {
                BOOST_STATIC_ASSERT(Dim == 1);
            }
            
            ConstArray(const Array<Dim, DataType>& rhs) :
                m_data(rhs.m_data),
                m_offset(rhs.m_offset)
            {
            }
            
            ConstArray<Dim, DataType>& operator=(const ConstArray<Dim, DataType>& rhs)
            {
                ConstArray<Dim, DataType> temp(rhs);
                Swap(temp);
                return *this;
            }

            const_reference operator[](index i) const 
            { 
                ASSERTL1(i < num_elements(), "Array Bounds Error.");
                return (*m_data)[i+m_offset]; 
            }
        
            const_iterator begin() const 
            {
                const_iterator result = m_data->begin();
                result += m_offset;
                return result;
            }
            
            const_iterator end() const 
            { 
                return m_data->end(); 
            }
            
            const element* get() const 
            { 
                return m_data->data() + m_offset; 
            }
            
            const element* data() const 
            { 
                return m_data->data() + m_offset; 
            }
            
//             size_type size() const 
//             { 
//                 return m_data->size() - m_offset; 
//             }
            
            size_type num_dimensions() const 
            { 
                return m_data->num_dimensions(); 
            }
            
            bool operator==(const ConstArray<Dim, DataType>& rhs) const
            {
                return m_data == rhs.m_data;
            }
            
            bool operator!=(const ConstArray<Dim, DataType>& rhs) const
            {
                return !(*this == rhs);
            }
            
            const size_type* shape() const { return m_data->shape(); }
            
            size_type num_elements() const { return m_data->num_elements() - m_offset; }

        protected:
            boost::shared_ptr<ArrayType>& GetData() { return m_data; }
            const boost::shared_ptr<ArrayType>& GetData() const { return m_data; }
            unsigned int GetOffset() const { return m_offset; }
            
        private:
            static void DeleteStorage(DataType* data, unsigned int num)
            {
                DestructionPolicy<DataType>::Destroy(data, num);
                MemoryManager<DataType>::RawDeallocate(data, num);
            }
            
            template<typename ExtentListType>
            void CreateStorage(const ExtentListType& extent)
            {
                unsigned int size = std::accumulate(extent.begin(), extent.end(), 1, 
                    std::multiplies<unsigned int>());
                DataType* storage = MemoryManager<DataType>::RawAllocate(size);
                m_data = MemoryManager<ArrayType>::AllocateSharedPtrD(
                        boost::bind(&ConstArray<Dim, DataType>::DeleteStorage, storage, size),
                        storage, extent);
            }
            
            void Swap(ConstArray<Dim, DataType>& rhs)
            {
                std::swap(m_data, rhs.m_data);
                std::swap(m_offset, rhs.m_offset);
            }
            
            boost::shared_ptr<ArrayType> m_data;
            unsigned int m_offset;
    };
    
    template<Dimension Dim, typename DataType>
    class Array : public ConstArray<Dim, DataType>
    {
        public:
            typedef ConstArray<Dim, DataType> BaseType;
            typedef typename BaseType::ArrayType ArrayType;
            typedef typename BaseType::index index;
            
            typedef typename ArrayType::iterator iterator;
            typedef typename ArrayType::reference reference;
            typedef typename ArrayType::const_reference const_reference;
            typedef typename ArrayType::element element;
            
        public:
            Array() :
                BaseType()
            {
            }
            
            explicit Array(unsigned int size) :
                BaseType(size)
            {
            }
            
            Array(unsigned int dim1Size, unsigned int dim2Size) :
                BaseType(dim1Size, dim2Size)
            {
            }
            
            Array(unsigned int dim1Size, unsigned int dim2Size, unsigned int dim3Size) :
                BaseType(dim1Size, dim2Size, dim3Size)
            {
            }
            
            Array(unsigned int arraySize, const DataType& initValue) :
                BaseType(arraySize, initValue)
            {
            }
            
            Array(unsigned int dim1Size, unsigned int dim2Size, const DataType& initValue) :
                BaseType(dim1Size, dim2Size, initValue)
            {
            }
            
            Array(unsigned int dim1Size, unsigned int dim2Size, unsigned int dim3Size, const DataType& initValue) :
                BaseType(dim1Size, dim2Size, dim3Size, initValue)
            {
            }

            
            Array(unsigned int arraySize, const DataType* d) :
                BaseType(arraySize, d)
            {
            }
            
            Array(const Array<Dim, DataType>& rhs) :
                BaseType(rhs)
            {
            }
            
            Array(const Array<Dim, DataType>& rhs, unsigned int offset) :
                BaseType(rhs, offset)
            {
            }
            
            explicit Array(const ConstArray<Dim, DataType>& rhs) :
                BaseType(rhs.num_elements())
            {
                std::copy(rhs.begin(), rhs.end(), begin());
            }
            
            Array<Dim, DataType>& operator=(const Array<Dim, DataType>& rhs)
            {
                ConstArray<Dim, DataType>::operator=(rhs);
                return *this;
            }
            
            element* get() { return this->GetData()->data(); }
            element* data() { return this->GetData()->data(); }
            
            const element* get() const { return this->GetData()->data(); }
            const element* data() const { return this->GetData()->data(); }
            
            reference operator[](index i) { return (*(this->GetData()))[i+this->GetOffset()]; }
            const_reference operator[](index i) const { return (*(this->GetData()))[i+this->GetOffset()]; }
            
            iterator begin()
            {
                iterator result = this->GetData()->begin();
                result += this->GetOffset();
                return result;
            }
            
            iterator end() { return this->GetData()->end(); }
            
            using BaseType::begin;
            using BaseType::end;
            
        private:
    };
        
        
    template<Dimension Dim, typename DataType>
    ConstArray<Dim, DataType> operator+(const ConstArray<Dim, DataType>& lhs, unsigned int offset)
    {
        return ConstArray<Dim, DataType>(lhs, offset);
    }
    
    template<Dimension Dim, typename DataType>
    Array<Dim, DataType> operator+(const Array<Dim, DataType>& lhs, unsigned int offset)
    {
        return Array<Dim, DataType>(lhs, offset);
    }
    
    template<typename DataType>
    void CopyArray(const ConstArray<OneD, DataType>& source, Array<OneD, DataType>& dest)
    {
        if( dest.num_elements() != source.num_elements() )
        {
            dest = Array<OneD, DataType>(source.num_elements());
        }
        
        std::copy(source.data(), source.data() + source.num_elements(), dest.data());
    }
    
    static Array<OneD, NekDouble> NullNekDouble1DArray;
}

#endif //NEKTAR_LIB_UTILITIES_BASIC_UTILS_SHARED_ARRAY_HPP
