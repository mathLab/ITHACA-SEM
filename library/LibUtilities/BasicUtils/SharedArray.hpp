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

#include <LibUtilities/BasicUtils/ArrayPolicies.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <LibUtilities/BasicUtils/NekPtr.hpp>
#include <LibUtilities/BasicConst/NektarUnivConsts.hpp>

#include <boost/multi_array.hpp>
#include <boost/shared_ptr.hpp>

namespace Nektar
{
    class LinearSystem;

    // Forward declaration for a ConstArray constructor.
    template<typename Dim, typename DataType>
    class Array;

    /// \brief 1D Array of constant elements with garbage collection and bounds checking.
    template<typename DataType>
    class Array<OneD, const DataType>
    {
        public:
            typedef DataType* ArrayType;
            typedef const DataType& const_reference;
            typedef DataType& reference;

            typedef const DataType* const_iterator;
            typedef DataType* iterator;

            typedef DataType element;
            typedef unsigned int size_type;


        public:
            /// \brief Creates an empty array.
            Array() :
                m_size(0),
                m_capacity(0),
                m_data(0),
                m_count(0),
                m_offset(0)
            {
                CreateStorage(m_capacity);
            }

            /// \brief Creates an array of size dim1Size.
            ///
            /// If DataType is a fundamental type (double, int, etc.), then the allocated array is
            /// uninitialized.  If it is any other type, each element is initialized with DataType's default
            /// constructor.
            explicit Array(unsigned int dim1Size) :
                m_size(dim1Size),
                m_capacity(dim1Size),
                m_data(0),
                m_count(0),
                m_offset(0)
            {
                CreateStorage(m_capacity);
                ArrayInitializationPolicy<DataType>::Initialize(m_data + 1, m_capacity);
            }

             /// \brief Creates a 1D array with each element
             /// initialized to an initial value.
            /// \param dim1Size The array's size.
            /// \param initValue Each element's initial value.
            ///
            /// If DataType is a fundamental type (double, int, etc.),
            /// then the initial value is copied directly into each
            /// element.  Otherwise, the DataType's copy constructor
            /// is used to initialize each element.
            Array(unsigned int dim1Size, const DataType& initValue) :
                m_size(dim1Size),
                m_capacity(dim1Size),
                m_data(0),
                m_count(0),
                m_offset(0)
            {
                CreateStorage(m_capacity);
                ArrayInitializationPolicy<DataType>::Initialize(m_data + 1, m_capacity, initValue);
            }

            /// \brief Creates a 1D array a copies data into it.
            /// \param dim1Size the array's size.
            /// \param data The data to copy.
            ///
            /// If DataType is a fundamental type (double, int, etc.), then data is copied
            /// directly into the underlying storage.  Otherwise, the DataType's copy constructor
            /// is used to copy each element.
            Array(unsigned int dim1Size, const DataType* data) :
                m_size(dim1Size),
                m_capacity(dim1Size),
                m_data(0),
                m_count(0),
                m_offset(0)
            {
                CreateStorage(m_capacity);
                ArrayInitializationPolicy<DataType>::Initialize(m_data + 1, m_capacity, data);
            }

            /// \brief Creates a 1D array that references rhs.
            /// \param dim1Size The size of the array.  This is useful
            ///                 when you want this array to reference
            ///                 a subset of the elements in rhs.
            /// \param rhs      Array to reference.
            /// This constructor creates an array that references rhs.
            /// Any changes to rhs will be reflected in this array.
            /// The memory for the array will only be deallocated when
            /// both rhs and this array have gone out of scope.
            Array(unsigned int dim1Size, const Array<OneD, const DataType>& rhs) :
                m_size(dim1Size),
                m_capacity(rhs.m_capacity),
                m_data(rhs.m_data),
                m_count(rhs.m_count),
                m_offset(rhs.m_offset)
            {
                *m_count += 1;
                ASSERTL0(m_size <= rhs.num_elements(), "Requested size is larger than input array size.");
            }

            /// \brief Creates a reference to rhs.
            Array(const Array<OneD, const DataType>& rhs) :
                m_size(rhs.m_size),
                m_capacity(rhs.m_capacity),
                m_data(rhs.m_data),
                m_count(rhs.m_count),
                m_offset(rhs.m_offset)
            {
                *m_count += 1;
            }

            ~Array()
            {
                if( m_count == 0 )
                {
                    return;
                }

                *m_count -= 1;
                if( *m_count == 0 )
                {
                    ArrayDestructionPolicy<DataType>::Destroy(m_data+1, m_capacity);
                    MemoryManager<DataType>::RawDeallocate(m_data, m_capacity+1);
                }
            }

            /// \brief Creates a reference to rhs.
            Array<OneD, const DataType>& operator=(const Array<OneD, const DataType>& rhs)
            {
                *m_count -= 1;
                if( *m_count == 0 )
                {
                    ArrayDestructionPolicy<DataType>::Destroy(m_data+1, m_capacity);
                    MemoryManager<DataType>::RawDeallocate(m_data, m_capacity+1);
                }

                m_data = rhs.m_data;
                m_capacity = rhs.m_capacity;
                m_count = rhs.m_count;
                *m_count += 1;
                m_offset = rhs.m_offset;
                m_size = rhs.m_size;
                return *this;
            }

            const_iterator begin() const { return m_data + m_offset + 1; }
            const_iterator end() const { return m_data + m_offset + m_size + 1; }

            const_reference operator[](unsigned int i) const
            {
                ASSERTL1(static_cast<size_type>(i) < m_size, (std::string("Element ") +
                    boost::lexical_cast<std::string>(i) + std::string(" requested in an array of size ") +
                    boost::lexical_cast<std::string>(m_size)));
                return *(m_data + i + m_offset + 1);
            }

            /// \brief Returns a c-style pointer to the underlying array.
            const element* get() const { return m_data+m_offset+1; }

            /// \brief Returns a c-style pointer to the underlying array.
            const element* data() const { return m_data+m_offset+1; }

            /// \brief Returns 1.
            size_type num_dimensions() const { return 1; }

            /// \brief Returns the array's size.
            size_type num_elements() const { return m_size; }

            size_type capacity() const { return m_capacity; }

            /// \brief Returns the array's offset.
            unsigned int GetOffset() const { return m_offset; }

            /// \brief Returns true is this array and rhs overlap.
            bool Overlaps(const Array<OneD, const DataType>& rhs) const
            {
                const element* start = get();
                const element* end = start + m_size;

                const element* rhs_start = rhs.get();
                const element* rhs_end = rhs_start + rhs.num_elements();

                return (rhs_start >= start && rhs_start <= end) ||
                       (rhs_end >= start && rhs_end <= end);
            }

            template<typename T1, typename T2>
            friend bool operator==(const Array<OneD, T1>&, const Array<OneD, T2>&);
            LIB_UTILITIES_EXPORT friend bool operator==(const Array<OneD, NekDouble>&, const Array<OneD, NekDouble>&);
            LIB_UTILITIES_EXPORT friend bool IsEqual(const Array<OneD, const NekDouble>&,const Array<OneD, const NekDouble>&,NekDouble);

            /// \brief Creates an array with a specified offset.
            ///
            /// The return value will reference the same array as lhs,
            /// but with an offset.
            ///
            /// For example, in the following:
            /// \code
            /// Array<OneD, const double> result = anArray + 10;
            /// \endcode
            /// result[0] == anArray[10];
            template<typename T>
            friend Array<OneD, T> operator+(const Array<OneD, T>& lhs, unsigned int offset);

            template<typename T>
            friend Array<OneD, T> operator+(unsigned int offset, const Array<OneD, T>& rhs);

        protected:
            unsigned int m_size;
            unsigned int m_capacity;
            //boost::shared_ptr<DataType> m_data;
            //NekPtr<DataType> m_data;
            DataType* m_data;
            unsigned int* m_count;
            unsigned int m_offset;


        private:
        //            struct DestroyArray
        //            {
        //                DestroyArray(unsigned int elements) :
        //                    m_elements(elements) {}
        //
        //                void operator()(DataType* p)
        //                {
        //                    ArrayDestructionPolicy<DataType>::Destroy(p, m_elements);
        //                    MemoryManager<DataType>::RawDeallocate(p, m_elements);
        //                }
        //                unsigned int m_elements;
        //            };
        //
        // boost::shared_ptr<DataType>
        // NekPtr<DataType>
        void
            CreateStorage(unsigned int size)
            {
                DataType* storage = MemoryManager<DataType>::RawAllocate(size+1);
                m_data = storage;
                m_count = (unsigned int*)storage;
                *m_count = 1;
                //return NekPtr<DataType>(storage, size);
                //return boost::shared_ptr<DataType>(storage,
                        //boost::bind(DeleteStorage<DataType>, storage, size) );
                //        DestroyArray(size), MemoryManager<DataType>() );
            }

            template<typename T>
            static Array<OneD, T> CreateWithOffset(const Array<OneD, T>& rhs, unsigned int offset)
            {
                Array<OneD, T> result(rhs);
                result.m_offset += offset;
                result.m_size = rhs.m_size - offset;
                return result;
            }

    };


    /// \brief 2D array with garbage collection and bounds checking.
    template<typename DataType>
    class Array<TwoD, const DataType>
    {
        public:
            typedef boost::multi_array_ref<DataType, 2> ArrayType;
            typedef typename ArrayType::const_reference const_reference;
            typedef typename ArrayType::reference reference;

            typedef typename ArrayType::index index;
            typedef typename ArrayType::const_iterator const_iterator;
            typedef typename ArrayType::iterator iterator;

            typedef typename ArrayType::element element;
            typedef typename ArrayType::size_type size_type;



        public:
            Array() :
                m_data(CreateStorage<DataType>(0, 0))
            {
            }

            /// \brief Constructs a 3 dimensional array.  The elements of the array are not initialized.
            Array(unsigned int dim1Size, unsigned int dim2Size) :
                m_data(CreateStorage<DataType>(dim1Size, dim2Size))
            {
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements());
            }

            Array(unsigned int dim1Size, unsigned int dim2Size, const DataType& initValue) :
                m_data(CreateStorage<DataType>(dim1Size, dim2Size))
            {
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements(), initValue);
            }

            Array(unsigned int dim1Size, unsigned int dim2Size, const DataType* data) :
                m_data(CreateStorage<DataType>(dim1Size, dim2Size))
            {
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements(), data);
            }

            Array(const Array<TwoD, const DataType>& rhs) :
                m_data(rhs.m_data)
            {
            }

            Array<TwoD, const DataType>& operator=(const Array<TwoD, const DataType>& rhs)
            {
                m_data = rhs.m_data;
                return *this;
            }

            template<typename T>
            friend bool operator==(const Array<TwoD, T>&, const Array<TwoD, T>&);
            LIB_UTILITIES_EXPORT friend bool operator==(const Array<TwoD, NekDouble>&, const Array<TwoD, NekDouble>&);
            LIB_UTILITIES_EXPORT friend bool IsEqual(const Array<TwoD, const NekDouble>&,const Array<TwoD, const NekDouble>&,NekDouble);

            const_iterator begin() const { return m_data->begin(); }
            const_iterator end() const { return m_data->end(); }
            const_reference operator[](index i) const { return (*m_data)[i]; }
            const element* get() const { return m_data->data(); }
            const element* data() const { return m_data->data(); }
            size_type num_dimensions() const { return m_data->num_dimensions(); }
            const size_type* shape() const { return m_data->shape(); }
            size_type num_elements() const { return m_data->num_elements(); }

            size_type GetRows() const { return m_data->shape()[0]; }
            size_type GetColumns() const { return m_data->shape()[1]; }

        protected:
            boost::shared_ptr<ArrayType> m_data;

        private:

    };

    /// \brief 3D array with garbage collection and bounds checking.
    template<typename DataType>
    class Array<ThreeD, const DataType>
    {
        public:
            typedef boost::multi_array_ref<DataType, 3> ArrayType;
            typedef typename ArrayType::const_reference const_reference;
            typedef typename ArrayType::reference reference;

            typedef typename ArrayType::index index;
            typedef typename ArrayType::const_iterator const_iterator;
            typedef typename ArrayType::iterator iterator;

            typedef typename ArrayType::element element;
            typedef typename ArrayType::size_type size_type;



        public:
            Array() :
                m_data(CreateStorage<DataType>(0, 0, 0))
            {
            }

            /// \brief Constructs a 3 dimensional array.  The elements of the array are not initialized.
            Array(unsigned int dim1Size, unsigned int dim2Size, unsigned int dim3Size) :
                m_data(CreateStorage<DataType>(dim1Size, dim2Size, dim3Size))
            {
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements());
            }

            Array(unsigned int dim1Size, unsigned int dim2Size, unsigned int dim3Size, const DataType& initValue) :
                m_data(CreateStorage<DataType>(dim1Size, dim2Size, dim3Size))
            {
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements(), initValue);
            }

            Array(const Array<ThreeD, const DataType>& rhs) :
                m_data(rhs.m_data)
            {
            }

            Array<ThreeD, const DataType>& operator=(const Array<ThreeD, const DataType>& rhs)
            {
                m_data = rhs.m_data;
                return *this;
            }

            const_iterator begin() const { return m_data->begin(); }
            const_iterator end() const { return m_data->end(); }
            const_reference operator[](index i) const { return (*m_data)[i]; }
            const element* get() const { return m_data->data(); }
            const element* data() const { return m_data->data(); }
            size_type num_dimensions() const { return m_data->num_dimensions(); }
            const size_type* shape() const { return m_data->shape(); }
            size_type num_elements() const { return m_data->num_elements(); }

        protected:
            boost::shared_ptr<ArrayType> m_data;

        private:

    };


    enum AllowWrappingOfConstArrays
    {
        eVECTOR_WRAPPER
    };

    /// \brief 1D Array
    ///
    /// \ref pageNekArrays
    ///
    /// Misc notes.
    ///
    /// Throught the 1D Array class you will see things like "using BaseType::begin" and
    /// "using BaseType::end".  This is necessary to bring the methods from the ConstArray
    /// into scope in Array class.  Typically this is not necessary, but since we have
    /// method names which match those in the base class, the base class names are hidden.
    /// Therefore, we have to explicitly bring them into scope to use them.
    template<typename DataType>
    class Array<OneD, DataType> : public Array<OneD, const DataType>
    {
        public:
            typedef Array<OneD, const DataType> BaseType;
            typedef typename BaseType::iterator iterator;
            typedef typename BaseType::reference reference;
            typedef typename BaseType::size_type size_type;
            typedef typename BaseType::element element;

        public:
            Array() :
                BaseType()
            {
            }

            explicit Array(unsigned int dim1Size) :
                BaseType(dim1Size)
            {
            }

            Array(unsigned int dim1Size, const DataType& initValue) :
                BaseType(dim1Size, initValue)
            {
            }

            Array(unsigned int dim1Size, const DataType* data) :
                BaseType(dim1Size, data)
            {
            }

            Array(unsigned int dim1Size, const Array<OneD, DataType>& rhs) :
                BaseType(dim1Size, rhs)
            {
            }

            Array(unsigned int dim1Size, const Array<OneD, const DataType>& rhs) :
                BaseType(dim1Size, rhs.data())
            {
            }

            Array(const Array<OneD, DataType>& rhs) :
                BaseType(rhs)
            {
            }

            Array(const Array<OneD, const DataType>& rhs) :
                BaseType(rhs.num_elements(), rhs.data())
            {
            }

            Array<OneD, DataType>& operator=(const Array<OneD, DataType>& rhs)
            {
                BaseType::operator=(rhs);
                return *this;
            }

            static Array<OneD, DataType> CreateWithOffset(const Array<OneD, DataType>& rhs, unsigned int offset)
            {
                Array<OneD, DataType> result(rhs);
                result.m_offset += offset;
                result.m_size = rhs.m_size - offset;
                return result;
            }

            using BaseType::begin;
            iterator begin() { return this->m_data + 1 +this->m_offset; }

            using BaseType::end;
            iterator end() { return this->m_data + 1 + this->m_offset + this->m_size; }

            using BaseType::operator[];
            reference operator[](unsigned int i)
            {
                ASSERTL1(static_cast<size_type>(i) < this->num_elements(), (std::string("Element ") +
                    boost::lexical_cast<std::string>(i) + std::string(" requested in an array of size ") +
                    boost::lexical_cast<std::string>(this->num_elements())));
                return (get())[i];
            }


            using BaseType::get;
            element* get() { return this->m_data + 1 +this->m_offset; }

            using BaseType::data;
            element* data() { return this->m_data + 1 +this->m_offset; }

            template<typename T1>
            friend class NekVector;

            template<typename T1, typename T3>
            friend class NekMatrix;

            friend class LinearSystem;

        protected:
            Array(const Array<OneD, const DataType>& rhs, AllowWrappingOfConstArrays a) :
                BaseType(rhs)
            {
            }

            void ChangeSize(unsigned int newSize)
            {
                ASSERTL1(newSize <= this->m_capacity, "Can't change an array size to something larger than its capacity.");
                this->m_size = newSize;
            }

        private:

    };

    /// \brief A 2D array.
    template<typename DataType>
    class Array<TwoD, DataType> : public Array<TwoD, const DataType>
    {
        public:
            typedef Array<TwoD, const DataType> BaseType;
            typedef typename BaseType::iterator iterator;
            typedef typename BaseType::reference reference;
            typedef typename BaseType::index index;
            typedef typename BaseType::size_type size_type;
            typedef typename BaseType::element element;

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

            Array(unsigned int dim1Size, unsigned int dim2Size, const DataType* data) :
                BaseType(dim1Size, dim2Size, data)
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

            using BaseType::begin;
            iterator begin() { return this->m_data->begin(); }

            using BaseType::end;
            iterator end() { return this->m_data->end(); }

            using BaseType::operator[];
            reference operator[](index i) { return (*this->m_data)[i]; }

            using BaseType::get;
            element* get() { return this->m_data->data(); }

            using BaseType::data;
            element* data() { return this->m_data->data(); }

        private:

    };

    /// \brief A 3D array.
    template<typename DataType>
    class Array<ThreeD, DataType> : public Array<ThreeD, const DataType>
    {
        public:
            typedef Array<ThreeD, const DataType> BaseType;
            typedef typename BaseType::iterator iterator;
            typedef typename BaseType::reference reference;
            typedef typename BaseType::index index;
            typedef typename BaseType::size_type size_type;
            typedef typename BaseType::element element;

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

            using BaseType::begin;
            iterator begin() { return this->m_data->begin(); }

            using BaseType::end;
            iterator end() { return this->m_data->end(); }

            using BaseType::operator[];
            reference operator[](index i) { return (*this->m_data)[i]; }

            using BaseType::get;
            element* get() { return this->m_data->data(); }

            using BaseType::data;
            element* data() { return this->m_data->data(); }

        private:

    };

    LIB_UTILITIES_EXPORT bool IsEqual(const Array<OneD, const NekDouble>& lhs,
                 const Array<OneD, const NekDouble>& rhs,
                 NekDouble tol = NekConstants::kNekZeroTol);
    LIB_UTILITIES_EXPORT bool operator==(const Array<OneD, NekDouble>& lhs, const Array<OneD, NekDouble>& rhs);

    template<typename T1, typename T2>
    bool operator==(const Array<OneD, T1>& lhs, const Array<OneD, T2>& rhs)
    {
        if( lhs.num_elements() != rhs.num_elements() )
        {
            return false;
        }

        if( lhs.data() == rhs.data() )
        {
            return true;
        }

        for(unsigned int i = 0; i < lhs.num_elements(); ++i)
        {
            if( lhs[i] != rhs[i] )
            {
                return false;
            }
        }

        return true;
    }

    template<typename T1, typename T2>
    bool operator!=(const Array<OneD, T1>& lhs, const Array<OneD, T2>& rhs)
    {
        return !(lhs == rhs);
    }

    template<typename DataType>
    Array<OneD, DataType> operator+(const Array<OneD, DataType>& lhs, unsigned int offset)
    {
        return Array<OneD, const DataType>::CreateWithOffset(lhs, offset);
    }

    template<typename DataType>
    Array<OneD, DataType> operator+(unsigned int offset, const Array<OneD, DataType>& rhs)
    {
        return Array<OneD, const DataType>::CreateWithOffset(rhs, offset);
    }

//    template<typename DataType>
//    Array<OneD, DataType> operator+(const Array<OneD, DataType>& lhs, unsigned int offset)
//    {
//        return Array<OneD, DataType>::CreateWithOffset(lhs, offset);
//    }

    template<typename ConstDataType, typename DataType>
    void CopyArray(const Array<OneD, ConstDataType>& source, Array<OneD, DataType>& dest)
    {
        if( dest.num_elements() != source.num_elements() )
        {
            dest = Array<OneD, DataType>(source.num_elements());
        }

        std::copy(source.data(), source.data() + source.num_elements(), dest.data());
    }

    template<typename ConstDataType, typename DataType>
    void CopyArrayN(const Array<OneD, ConstDataType>& source, Array<OneD, DataType>& dest, unsigned int n)
    {
        if( dest.num_elements() != n )
        {
            dest = Array<OneD, DataType>(n);
        }

        std::copy(source.data(), source.data() + n, dest.data());
    }

    static Array<OneD, NekDouble> NullNekDouble1DArray;
    static Array<OneD, Array<OneD, NekDouble> > NullNekDoubleArrayofArray;

    LIB_UTILITIES_EXPORT bool IsEqual(const Array<TwoD, const NekDouble>& lhs,
                 const Array<TwoD, const NekDouble>& rhs,
                 NekDouble tol = NekConstants::kNekZeroTol);
    LIB_UTILITIES_EXPORT bool operator==(const Array<TwoD, NekDouble>& lhs, const Array<TwoD, NekDouble>& rhs) ;

    template<typename DataType>
    bool operator==(const Array<TwoD, DataType>& lhs, const Array<TwoD, DataType>& rhs)
    {
        return *lhs.m_data == *rhs.m_data;
    }

    template<typename DataType>
    bool operator!=(const Array<TwoD, DataType>& lhs, const Array<TwoD, DataType>& rhs)
    {
        return !(lhs == rhs);
    }


}

#endif //NEKTAR_LIB_UTILITIES_BASIC_UTILS_SHARED_ARRAY_HPP
