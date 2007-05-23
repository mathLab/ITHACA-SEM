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
#include <LibUtilities/BasicUtils/ArrayPolicies.hpp>

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
    class ArrayImpl
    {
        public:
            typedef boost::multi_array_ref<DataType, Dim> ArrayType;

            typedef typename ArrayType::const_reference const_reference;
            typedef typename ArrayType::index index;
            typedef typename ArrayType::const_iterator const_iterator;
            typedef typename ArrayType::element element;
            typedef typename ArrayType::size_type size_type;
            typedef typename ArrayType::reference reference;
            typedef typename ArrayType::iterator iterator;
            
        public:
            ArrayImpl() :
                m_data()
            {
            }
            
            ArrayImpl(const ArrayImpl<Dim, DataType>& rhs) :
                m_data(rhs.m_data)
            {
            }
            
            ArrayImpl<Dim, DataType>& operator=(const ArrayImpl<Dim, DataType>& rhs)
            {
                if( m_data != rhs.m_data )
                {
                    m_data = rhs.m_data;
                }
                return *this;
            }
            
            const_reference operator[](index i) const 
            { 
                ASSERTL1(i < num_elements(), "Array Bounds Error.");
                return (*m_data)[i]; 
            }
            
            reference operator[](index i) 
            { 
                ASSERTL1(i < num_elements(), "Array Bounds Error.");
                return (*m_data)[i]; 
            }
        
            const_iterator begin() const 
            {
                const_iterator result = m_data->begin();
                return result;
            }
            
            iterator begin() 
            {
                iterator result = m_data->begin();
                return result;
            }
            
            const_iterator end() const 
            { 
                return m_data->end(); 
            }
            
            iterator end() 
            { 
                return m_data->end(); 
            }
            
            const element* get() const 
            { 
                return m_data->data(); 
            }
            
            element* get() 
            { 
                return m_data->data(); 
            }
            
            const element* data() const 
            { 
                return m_data->data(); 
            }
            
            element* data() 
            { 
                return m_data->data(); 
            }
            
//             size_type size() const 
//             { 
//                 return m_data->size() - m_offset; 
//             }
            
            size_type num_dimensions() const 
            { 
                return m_data->num_dimensions(); 
            }
            
            
            
            const size_type* shape() const { return m_data->shape(); }
            
            size_type num_elements() const { return m_data->num_elements(); }
            
        protected:
            static void DeleteStorage(DataType* data, unsigned int num)
            {
                ArrayDestructionPolicy<DataType>::Destroy(data, num);
                MemoryManager<DataType>::RawDeallocate(data, num);
            }
            
            template<typename ExtentListType>
            void CreateStorage(const ExtentListType& extent)
            {
                unsigned int size = std::accumulate(extent.begin(), extent.end(), 1, 
                    std::multiplies<unsigned int>());
                DataType* storage = MemoryManager<DataType>::RawAllocate(size);
                m_data = MemoryManager<ArrayType>::AllocateSharedPtrD(
                        boost::bind(&ArrayImpl<Dim, DataType>::DeleteStorage, storage, size),
                        storage, extent);
            }

            boost::shared_ptr<ArrayType> m_data;
    };
    
    
    template<Dimension Dim, typename DataType>
    class ConstArray;
    
    template<typename DataType>
    class ConstArray<OneD, DataType> : protected ArrayImpl<OneD, DataType>
    {
        public:
            typedef ArrayImpl<OneD, DataType> BaseType;
            typedef typename BaseType::const_reference const_reference;
            typedef typename BaseType::reference reference;
            typedef typename BaseType::iterator iterator;
            typedef typename BaseType::const_iterator const_iterator;
            typedef typename BaseType::index index;
            typedef typename BaseType::element element;
            typedef typename BaseType::size_type size_type;
            
        public:
            /// \brief Constructs an empty array.
            ///
            /// An empty array will still have a small amount of memory allocated for it, typically around 
            /// 4 bytes.  
            /// \post size() == 0
            /// \post begin() == end()
            ConstArray() :
                BaseType(),
                m_offset(0)
            {
                std::vector<unsigned int> extents(1, 0);
                this->CreateStorage(extents);
            }
            
            /// \brief Constructs a 1D array.  Every element is given the value initValue.
            ConstArray(unsigned int arraySize) :
                BaseType(),
                m_offset(0)
            {
                std::vector<unsigned int> extents(1, arraySize);
                this->CreateStorage(extents);
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements());
            }
            
            /// \brief Constructs a 1D array.  Every element is given the value initValue.
            ConstArray(unsigned int arraySize, const DataType& initValue) :
                BaseType(),
                m_offset(0)
            {
                std::vector<unsigned int> extents(1, arraySize);
                this->CreateStorage(extents);
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements(), initValue);
            }
            
            ConstArray(unsigned int arraySize, const DataType* d) :
                BaseType(),
                m_offset(0)
            {
                std::vector<unsigned int> extents(1, arraySize);
                this->CreateStorage(extents);
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements(), d);
            }
            
            ConstArray(const ConstArray<OneD, DataType>& rhs) :
                BaseType(rhs),
                m_offset(rhs.m_offset)
            {
            }
            
            ConstArray(const ConstArray<OneD, DataType>& rhs, unsigned int off) :
                BaseType(rhs),
                m_offset(off)
            {
            }
            
            ConstArray<OneD, DataType>& operator=(const ConstArray<OneD, DataType>& rhs)
            {
                BaseType::operator=(rhs);
                m_offset = rhs.m_offset;
                return *this;
            }
            
            const_reference operator[](index i) const 
            {
                return BaseType::operator[](i+m_offset);
            }
        
            const_iterator begin() const 
            {
                const_iterator result = BaseType::begin();
                result += m_offset;
                return result;
            }
            
            const_iterator end() const 
            { 
                return BaseType::end();
            }
            
            const element* get() const 
            { 
                return BaseType::data() + m_offset;
            }
            
            const element* data() const 
            { 
                return BaseType::data() + m_offset;
            }
            
            
            size_type num_elements() const 
            { 
                return BaseType::num_elements() - m_offset;
            }

            using BaseType::shape;
            using BaseType::num_dimensions;
            
            template<typename T>
            friend bool operator==(const ConstArray<OneD, T>&, const ConstArray<OneD, T>&);
            
            unsigned int GetOffset() const { return m_offset; }
        private:
            using BaseType::m_data;
            unsigned int m_offset;
    };
    
    template<typename DataType>
    class ConstArray<TwoD, DataType> : protected ArrayImpl<TwoD, DataType>
    {
        public:
            typedef ArrayImpl<TwoD, DataType> BaseType;
            
        public:
            ConstArray() :
                BaseType()
            {
                std::vector<unsigned int> extents(1, 0);
                this->CreateStorage(extents);
            }
            
            /// \brief Constructs a 2 dimensional array.  The elements of the array are not initialized.
            ConstArray(unsigned int dim1Size, unsigned int dim2Size) :
                BaseType()
            {
                unsigned int vals[] = {dim1Size, dim2Size};
                std::vector<unsigned int> extents(vals, vals+2);
                this->CreateStorage(extents);
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements());
            }
            
            ConstArray(unsigned int dim1Size, unsigned int dim2Size, const DataType& initValue) :
                BaseType()
            {
                unsigned int vals[] = {dim1Size, dim2Size};
                std::vector<unsigned int> extents(vals, vals+2);
                this->CreateStorage(extents);
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements(), initValue);
            }
            
            ConstArray(const ConstArray<TwoD, DataType>& rhs) :
                BaseType(rhs)
            {
            }
            
            ConstArray<TwoD, DataType>& operator=(const ConstArray<TwoD, DataType>& rhs)
            {
                BaseType::operator=(rhs);
                return *this;
            }
            
            using BaseType::begin;
            using BaseType::end;
            using BaseType::operator[];
            using BaseType::data;
            using BaseType::get;
            using BaseType::num_elements;
            using BaseType::shape;
            using BaseType::num_dimensions;
            
        private:
            using BaseType::m_data;
            
            
    };
    
    template<typename DataType>
    class ConstArray<ThreeD, DataType> : protected ArrayImpl<ThreeD, DataType>
    {
        public:
            typedef ArrayImpl<ThreeD, DataType> BaseType;
            
        public:
            ConstArray() :
                BaseType()
            {
                std::vector<unsigned int> extents(1, 0);
                this->CreateStorage(extents);
            }
            
            /// \brief Constructs a 3 dimensional array.  The elements of the array are not initialized.
            ConstArray(unsigned int dim1Size, unsigned int dim2Size, unsigned int dim3Size) :
                BaseType()
            {
                unsigned int vals[] = {dim1Size, dim2Size, dim3Size};
                std::vector<unsigned int> extents(vals, vals+3);
                this->CreateStorage(extents);
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements());
            }
            
            ConstArray(unsigned int dim1Size, unsigned int dim2Size, unsigned int dim3Size, const DataType& initValue) :
                BaseType()
            {
                unsigned int vals[] = {dim1Size, dim2Size, dim3Size};
                std::vector<unsigned int> extents(vals, vals+3);
                this->CreateStorage(extents);
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements(), initValue);
            }
            
            ConstArray(const ConstArray<ThreeD, DataType>& rhs) :
                BaseType(rhs)
            {
            }
            
            ConstArray<ThreeD, DataType>& operator=(const ConstArray<ThreeD, DataType>& rhs)
            {
                BaseType::operator=(rhs);
                return *this;
            }
            
            using BaseType::begin;
            using BaseType::end;
            using BaseType::operator[];
            using BaseType::data;
            using BaseType::get;
            using BaseType::num_elements;
            using BaseType::shape;
            using BaseType::num_dimensions;
            
        private:
            using BaseType::m_data;
    };

    template<Dimension Dim, typename DataType>
    class Array;
    
    
    /// Misc notes.
    ///
    /// Throught the 1D Array class you will see things like "using BaseType::begin" and 
    /// "using BaseType::end".  This is necessary to bring the methods from the ConstArray
    /// into scope in Array class.  Typically this is not necessary, but since we have 
    /// method names which match those in the base class, the base class names are hidden.
    /// Therefore, we have to explicitly bring them into scope to use them.
    template<typename DataType>
    class Array<OneD, DataType> : public ConstArray<OneD, DataType>
    {
        public:
            typedef ConstArray<OneD, DataType> BaseType;
            typedef typename ArrayImpl<OneD, DataType>::ArrayType ArrayType;
            typedef typename ArrayType::index index;
            
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
            
            Array(unsigned int arraySize, const DataType& initValue) :
                BaseType(arraySize, initValue)
            {
            }
            

            Array(unsigned int arraySize, const DataType* d) :
                BaseType(arraySize, d)
            {
            }
            
            Array(const Array<OneD, DataType>& rhs) :
                BaseType(rhs)
            {
            }
            
            Array(const Array<OneD, DataType>& rhs, unsigned int offset) :
                BaseType(rhs, offset)
            {
            }
            
            Array<OneD, DataType>& operator=(const Array<OneD, DataType>& rhs)
            {
                ConstArray<OneD, DataType>::operator=(rhs);
                return *this;
            }
            
            using BaseType::get;
            element* get() 
            {
                return ArrayImpl<OneD, DataType>::get() + this->GetOffset(); 
            }
            
            using BaseType::data;
            element* data() 
            { 
                return ArrayImpl<OneD, DataType>::data() + this->GetOffset(); 
            }
            
            using BaseType::operator[];
            reference operator[](index i) 
            { 
                return ArrayImpl<OneD, DataType>::operator[](i + this->GetOffset());
            }
            
            using BaseType::begin;
            iterator begin()
            {
                iterator result = ArrayImpl<OneD, DataType>::begin();
                result += this->GetOffset();
                return result;
            }
            
            using BaseType::end;
            iterator end() 
            { 
                return ArrayImpl<OneD, DataType>::end();
            }

        private:
    };
        
        
    template<typename DataType>
    class Array<TwoD, DataType> : public ConstArray<TwoD, DataType>
    {
        public:
            typedef ConstArray<TwoD, DataType> BaseType;
            
        public:
            Array() :
                BaseType()
            {
            }
            
            Array(unsigned int dim1Size, unsigned int dim2Size) :
                BaseType(dim1Size, dim2Size)
            {
            }
            
            
            Array(unsigned int dim1Size, unsigned int dim2Size, const DataType& initValue) :
                BaseType(dim1Size, dim2Size, initValue)
            {
            }
            
            Array(const Array<TwoD, DataType>& rhs) :
                BaseType(rhs)
            {
            }
            
            Array<TwoD, DataType>& operator=(const Array<TwoD, DataType>& rhs)
            {
                BaseType::operator=(rhs);
                return *this;
            }
    };
            
    template<typename DataType>
    class Array<ThreeD, DataType> : public ConstArray<ThreeD, DataType>
    {
        public:
            typedef ConstArray<ThreeD, DataType> BaseType;
            
        public:
            Array() :
                BaseType()
            {
            }
            
            Array(unsigned int dim1Size, unsigned int dim2Size, unsigned int dim3Size) :
                BaseType(dim1Size, dim2Size, dim3Size)
            {
            }
            
            Array(unsigned int dim1Size, unsigned int dim2Size, unsigned int dim3Size, const DataType& initValue) :
                BaseType(dim1Size, dim2Size, dim3Size, initValue)
            {
            }
            
            Array(const Array<ThreeD, DataType>& rhs) :
                BaseType(rhs)
            {
            }
            
            Array<ThreeD, DataType>& operator=(const Array<ThreeD, DataType>& rhs)
            {
                BaseType::operator=(rhs);
                return *this;
            }
    };
    
    template<typename DataType>
    bool operator==(const ConstArray<OneD, DataType>& lhs, const ConstArray<OneD, DataType>& rhs) 
    {
        return lhs.m_data == rhs.m_data;
    }
    
    template<typename DataType>
    bool operator!=(const ConstArray<OneD, DataType>& lhs, const ConstArray<OneD, DataType>& rhs) 
    {
        return !(lhs == rhs);
    }

    template<typename DataType>
    ConstArray<OneD, DataType> operator+(const ConstArray<OneD, DataType>& lhs, unsigned int offset)
    {
        return ConstArray<OneD, DataType>(lhs, offset);
    }
    
    template<typename DataType>
    Array<OneD, DataType> operator+(const Array<OneD, DataType>& lhs, unsigned int offset)
    {
        return Array<OneD, DataType>(lhs, offset);
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
