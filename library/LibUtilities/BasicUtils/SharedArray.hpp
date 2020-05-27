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
#include <LibUtilities/BasicUtils/Deprecated.hpp>
#include <LibUtilities/BasicUtils/RealComparison.hpp>

#include <boost/core/ignore_unused.hpp>
#include <boost/multi_array.hpp>

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
#ifdef WITH_PYTHON
        struct PythonInfo {
            void *m_pyObject; // Underlying PyObject pointer
            void (*m_callback)(void *); // Callback
        };
#endif
        public:
            typedef DataType* ArrayType;
            typedef const DataType& const_reference;
            typedef DataType& reference;

            typedef const DataType* const_iterator;
            typedef DataType* iterator;

            typedef DataType element;
            typedef size_t size_type;


        public:
            /// \brief Creates an empty array.
            Array() :
#ifdef WITH_PYTHON
                m_pythonInfo(nullptr),
#endif
                m_size( 0 ),
                m_capacity( 0 ),
                m_data( nullptr ),
                m_count( nullptr ),
                m_offset( 0 )
            {
                CreateStorage(m_capacity);
            }

            /// \brief Creates an array of size dim1Size.
            ///
            /// If DataType is a fundamental type (double, int, etc.), then the allocated array is
            /// uninitialized.  If it is any other type, each element is initialized with DataType's default
            /// constructor.
            explicit Array(size_type dim1Size) :
#ifdef WITH_PYTHON
                m_pythonInfo(nullptr),
#endif
                m_size( dim1Size ),
                m_capacity( dim1Size ),
                m_data( nullptr ),
                m_count( nullptr ),
                m_offset( 0 )
            {
                CreateStorage(m_capacity);
                ArrayInitializationPolicy<DataType>::Initialize( m_data, m_capacity );
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
            Array(size_type dim1Size, const DataType& initValue) :
#ifdef WITH_PYTHON
                m_pythonInfo(nullptr),
#endif
                m_size( dim1Size ),
                m_capacity( dim1Size ),
                m_data( nullptr ),
                m_count( nullptr ),
                m_offset( 0 )
            {
                CreateStorage(m_capacity);
                ArrayInitializationPolicy<DataType>::Initialize( m_data, m_capacity, initValue );
            }

            /// \brief Creates a 1D array a copies data into it.
            /// \param dim1Size the array's size.
            /// \param data The data to copy.
            ///
            /// If DataType is a fundamental type (double, int, etc.), then data is copied
            /// directly into the underlying storage.  Otherwise, the DataType's copy constructor
            /// is used to copy each element.
            Array(size_type dim1Size, const DataType* data) :
#ifdef WITH_PYTHON
                m_pythonInfo(nullptr),
#endif
                m_size( dim1Size ),
                m_capacity( dim1Size ),
                m_data( nullptr ),
                m_count( nullptr ),
                m_offset( 0 )
            {
                CreateStorage(m_capacity);
                ArrayInitializationPolicy<DataType>::Initialize( m_data, m_capacity, data );
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
            Array(size_type dim1Size, const Array<OneD, const DataType>& rhs) :
#ifdef WITH_PYTHON
                m_pythonInfo(rhs.m_pythonInfo),
#endif
                m_size(dim1Size),
                m_capacity(rhs.m_capacity),
                m_data(rhs.m_data),
                m_count(rhs.m_count),
                m_offset(rhs.m_offset)
            {
                *m_count += 1;
                ASSERTL0(m_size <= rhs.size(), "Requested size is \
                    larger than input array size.");
            }

#ifdef WITH_PYTHON
            /// \brief Creates a 1D array a copies data into it.
            /// \param dim1Size the array's size.
            /// \param data The data to reference.
            /// \param memory_pointer Pointer to the memory address of the array
            /// \param python_decrement Pointer to decrementer
            Array(size_type dim1Size, DataType* data, void* memory_pointer,
                    void (*python_decrement)(void *)) :
                m_size( dim1Size ),
                m_capacity( dim1Size ),
                m_data( data ),
                m_count( nullptr ),
                m_offset( 0 )
            {
                m_count = new size_type();
                *m_count = 1;

                m_pythonInfo = new PythonInfo *();
                *m_pythonInfo = new PythonInfo();
                (*m_pythonInfo)->m_callback = python_decrement;
                (*m_pythonInfo)->m_pyObject = memory_pointer;
            }
#endif

            /// \brief Creates a reference to rhs.
            Array(const Array<OneD, const DataType>& rhs) :
#ifdef WITH_PYTHON
                m_pythonInfo(rhs.m_pythonInfo),
#endif
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
                if( m_count == nullptr )
                {
                    return;
                }

                *m_count -= 1;
                if( *m_count == 0 )
                {
#ifdef WITH_PYTHON
                    if (*m_pythonInfo == nullptr)
                    {
                        ArrayDestructionPolicy<DataType>::Destroy( m_data, m_capacity );
                        MemoryManager<DataType>::RawDeallocate( m_data, m_capacity );
                    }
                    else
                    {
                        (*m_pythonInfo)->m_callback((*m_pythonInfo)->m_pyObject);
                        delete *m_pythonInfo;
                    }

                    delete m_pythonInfo;

#else

                    ArrayDestructionPolicy<DataType>::Destroy( m_data, m_capacity );
                    MemoryManager<DataType>::RawDeallocate( m_data, m_capacity );

#endif

                    delete m_count; // Clean up the memory used for the reference count.
                }
            }

            /// \brief Creates a reference to rhs.
            Array<OneD, const DataType>& operator=(const Array<OneD, const DataType>& rhs)
            {
                *m_count -= 1;
                if( *m_count == 0 )
                {
#ifdef WITH_PYTHON
                    if (*m_pythonInfo == nullptr)
                    {
                        ArrayDestructionPolicy<DataType>::Destroy( m_data, m_capacity );
                        MemoryManager<DataType>::RawDeallocate( m_data, m_capacity );
                    }
                    else if ((*rhs.m_pythonInfo) != nullptr && (*m_pythonInfo)->m_pyObject != (*rhs.m_pythonInfo)->m_pyObject)
                    {
                        (*m_pythonInfo)->m_callback((*m_pythonInfo)->m_pyObject);
                        delete *m_pythonInfo;
                    }

                    delete m_pythonInfo;
#else

                    ArrayDestructionPolicy<DataType>::Destroy( m_data, m_capacity );
                    MemoryManager<DataType>::RawDeallocate( m_data, m_capacity );
#endif
                    delete m_count; // Clean up the memory used for the reference count.
                }

                m_data = rhs.m_data;
                m_capacity = rhs.m_capacity;
                m_count = rhs.m_count;
                *m_count += 1;
                m_offset = rhs.m_offset;
                m_size = rhs.m_size;
#ifdef WITH_PYTHON
                m_pythonInfo = rhs.m_pythonInfo;
#endif
                return *this;
            }

            const_iterator begin() const { return m_data + m_offset; }
            const_iterator end() const { return m_data + m_offset + m_size; }

            const_reference operator[](size_type i) const
            {
                ASSERTL1(i < m_size,
                         std::string("Element ") + std::to_string(i) +
                         std::string(" requested in an array of size ") +
                         std::to_string(m_size));
                return *( m_data + i + m_offset );
            }

            /// \brief Returns a c-style pointer to the underlying array.
            const element* get() const { return m_data + m_offset; }

            /// \brief Returns a c-style pointer to the underlying array.
            const element* data() const { return m_data + m_offset; }

            /// \brief Returns 1.
            size_type num_dimensions() const { return 1; }

            /// \brief Returns the array's size.
            size_type size() const { return m_size; }

            /// \brief Returns the array's size.
            /// Deprecated
            DEPRECATED(5.1.0, size) size_type num_elements() const
            {
                WARNINGL1(false,
                          "member function num_elements() is deprecated, "
                          "use size() instead.");
                return m_size;
            }

            /// \brief Returns the array's capacity.
            size_type capacity() const { return m_capacity; }

            /// \brief Returns the array's offset.
            size_type GetOffset() const { return m_offset; }

            /// \brief Returns the array's reference counter.
            size_type GetCount() const { return *m_count; }

            /// \brief Returns true is this array and rhs overlap.
            bool Overlaps(const Array<OneD, const DataType>& rhs) const
            {
                const element* start = get();
                const element* end = start + m_size;

                const element* rhs_start = rhs.get();
                const element* rhs_end = rhs_start + rhs.size();

                return (rhs_start >= start && rhs_start <= end) ||
                       (rhs_end >= start && rhs_end <= end);
            }

#ifdef WITH_PYTHON
            bool IsPythonArray()
            {
                return *m_pythonInfo != nullptr;
            }

            void ToPythonArray(void* memory_pointer, void (*python_decrement)(void *))
            {
                *m_pythonInfo = new PythonInfo();
                (*m_pythonInfo)->m_callback = python_decrement;
                (*m_pythonInfo)->m_pyObject = memory_pointer;
            }
#endif

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
            friend Array<OneD, T> operator+(const Array<OneD, T>& lhs,
                typename Array<OneD, T>::size_type offset);

            template<typename T>
            friend Array<OneD, T> operator+(typename Array<OneD, T>::size_type offset,
                const Array<OneD, T>& rhs);

        protected:

#ifdef WITH_PYTHON
            PythonInfo **m_pythonInfo;
#endif

            size_type m_size;
            size_type m_capacity;
            DataType* m_data;

            // m_count points to an integer used as a reference count to this array's data (m_data).
            // Previously, the reference count was stored in the first 4 bytes of the m_data array.
            size_type* m_count;

            size_type m_offset;

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
        void
            CreateStorage( size_type size )
            {
                DataType* storage = MemoryManager<DataType>::RawAllocate( size );
                m_data = storage;

                // Allocate an integer to hold the reference count.  Note 1, all arrays that share this array's
                // data (ie, point to m_data) will also share the m_count data.  Note 2, previously m_count
                // pointed to "(unsigned int*)storage".
                m_count = new size_type();
                *m_count = 1;
#ifdef WITH_PYTHON
                m_pythonInfo = new PythonInfo*();
                *m_pythonInfo = nullptr;
#endif
            }

            template<typename T>
            static Array<OneD, T> CreateWithOffset(const Array<OneD, T>& rhs, size_type offset)
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

            /// \brief Constructs a 2 dimensional array.  The elements of the array are not initialized.
            Array(size_type dim1Size, size_type dim2Size) :
                m_data(CreateStorage<DataType>(dim1Size, dim2Size))
            {
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(),
                    m_data->num_elements());
            }

            Array(size_type dim1Size, size_type dim2Size, const DataType& initValue) :
                m_data(CreateStorage<DataType>(dim1Size, dim2Size))
            {
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(),
                    m_data->num_elements(), initValue);
            }

            Array(size_type dim1Size, size_type dim2Size, const DataType* data) :
                m_data(CreateStorage<DataType>(dim1Size, dim2Size))
            {
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(),
                    m_data->num_elements(), data);
            }

            Array(const Array<TwoD, const DataType>& rhs) :
                m_data(rhs.m_data)
            {
            }

            Array<TwoD, const DataType>& operator=(const Array<TwoD,
                const DataType>& rhs)
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
            // m_data is a shared_ptr to a boost::multi_array_ref
            size_type size() const { return m_data->num_elements(); }
            // deprecated interface
            DEPRECATED(5.1.0, size) size_type num_elements() const
            {
                WARNINGL1(false,
                          "member function num_elements() is deprecated, "
                          "use size() instead.");
                return m_data->num_elements();
            }

            size_type GetRows() const { return m_data->shape()[0]; }
            size_type GetColumns() const { return m_data->shape()[1]; }

        protected:
            std::shared_ptr<ArrayType> m_data;

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
    /// Through out the 1D Array class you will see things like "using BaseType::begin" and
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

            explicit Array(size_type dim1Size) :
                BaseType(dim1Size)
            {
            }

            Array(size_type dim1Size, const DataType& initValue) :
                BaseType(dim1Size, initValue)
            {
            }

            Array(size_type dim1Size, const DataType* data) :
                BaseType(dim1Size, data)
            {
            }

            Array(size_type dim1Size, const Array<OneD, DataType>& rhs) :
                BaseType(dim1Size, rhs)
            {
            }

            Array(size_type dim1Size, const Array<OneD, const DataType>& rhs) :
                BaseType(dim1Size, rhs.data())
            {
            }

            Array(const Array<OneD, DataType>& rhs) :
                BaseType(rhs)
            {
            }

            Array(const Array<OneD, const DataType>& rhs) :
                BaseType(rhs.size(), rhs.data())
            {
            }

#ifdef WITH_PYTHON
            Array(size_type dim1Size, DataType* data, void* memory_pointer, void (*python_decrement)(void *)) :
                BaseType(dim1Size, data, memory_pointer, python_decrement)
            {
            }

            void ToPythonArray(void* memory_pointer, void (*python_decrement)(void *))
            {
                BaseType::ToPythonArray(memory_pointer, python_decrement);
            }
#endif

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

            typename Array<OneD, const DataType>::const_iterator
                     begin() const { return Array<OneD, const DataType>::begin(); }
            iterator begin() { return this->m_data + this->m_offset; }

            using BaseType::end;
            iterator end() { return this->m_data + this->m_offset + this->m_size; }

            using BaseType::operator[];
            reference operator[](size_type i)
            {
                ASSERTL1(static_cast<size_type>(i) < this->size(),
                         std::string("Element ") + std::to_string(i) +
                         std::string(" requested in an array of size ") +
                         std::to_string(this->size()));
                return (get())[i];
            }


            using BaseType::get;
            element* get() { return this->m_data + this->m_offset; }

            using BaseType::data;
            element* data() { return this->m_data + this->m_offset; }

            template<typename T1>
            friend class NekVector;

            template<typename T1, typename T3>
            friend class NekMatrix;

            friend class LinearSystem;

        protected:
            Array(const Array<OneD, const DataType>& rhs, AllowWrappingOfConstArrays a) :
                BaseType(rhs)
            {
                boost::ignore_unused(a);
            }

            void ChangeSize(size_type newSize)
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

            Array(size_type dim1Size, size_type dim2Size) :
                BaseType(dim1Size, dim2Size)
            {
            }

            Array(size_type dim1Size, size_type dim2Size, const DataType& initValue) :
                BaseType(dim1Size, dim2Size, initValue)
            {
            }

            Array(size_type dim1Size, size_type dim2Size, const DataType* data) :
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

    // compare whatever
    template<typename T>
    inline bool IsEqualImpl(const T& lhs, const T& rhs, std::false_type)
    {
        return lhs == rhs;
    }

    // compare floating point value
    template<typename T>
    inline bool IsEqualImpl(const T& lhs, const T& rhs, std::true_type)
    {
        return LibUtilities::IsRealEqual(lhs, rhs);
    }

    template<typename T>
    inline bool IsEqual(const T& lhs, const T& rhs)
    {
        return IsEqualImpl(lhs, rhs, std::is_floating_point<T>());
    }

    template<typename T1, typename T2>
    bool operator==(const Array<OneD, T1>& lhs,
                    const Array<OneD, T2>& rhs)
    {
        if (lhs.size() != rhs.size())
        {
            return false;
        }

        if (lhs.data() == rhs.data())
        {
            return true;
        }

        typename Array<OneD, T1>::size_type size_value(lhs.size());
        for (typename Array<OneD, T1>::size_type i = 0; i < size_value; ++i)
        {
            if (!IsEqual(lhs[i], rhs[i]))
            {
                return false;
            }
        }

        return true;
    }

    template<typename T1, typename T2>
    bool operator!=(const Array<OneD, T1>& lhs,
                    const Array<OneD, T2>& rhs)
    {
        return !(lhs == rhs);
    }

    template<typename DataType>
    Array<OneD, DataType> operator+(const Array<OneD, DataType>& lhs,
        typename Array<OneD, DataType>::size_type offset)
    {
        return Array<OneD, const DataType>::CreateWithOffset(lhs, offset);
    }

    template<typename DataType>
    Array<OneD, DataType> operator+(typename Array<OneD, DataType>::size_type offset,
        const Array<OneD, DataType>& rhs)
    {
        return Array<OneD, const DataType>::CreateWithOffset(rhs, offset);
    }

    template<typename ConstDataType, typename DataType>
    void CopyArray(const Array<OneD, ConstDataType>& source, Array<OneD, DataType>& dest)
    {
        if( dest.size() != source.size() )
        {
            dest = Array<OneD, DataType>(source.size());
        }

        std::copy(source.data(), source.data() + source.size(), dest.data());
    }

    template<typename ConstDataType, typename DataType>
    void CopyArrayN(const Array<OneD, ConstDataType>& source, Array<OneD,
        DataType>& dest, typename Array<OneD, DataType>::size_type n)
    {
        if( dest.size() != n )
        {
            dest = Array<OneD, DataType>(n);
        }

        std::copy(source.data(), source.data() + n, dest.data());
    }

    static Array<OneD, int> NullInt1DArray;
    static Array<OneD, NekDouble> NullNekDouble1DArray;
    static Array<OneD, Array<OneD, NekDouble> > NullNekDoubleArrayofArray;
    static Array<OneD, Array<OneD, Array<OneD, NekDouble> > > 
            NullNekDoubleArrayofArrayofArray;

    template<class T>
    using TensorOfArray1D = Array<OneD, T>;
    template<class T>
    using TensorOfArray2D = Array<OneD, Array<OneD, T>>;
    template<class T>
    using TensorOfArray3D = Array<OneD, Array<OneD, Array<OneD, T>>>;

    template<typename T1, typename T2>
    bool operator==(const Array<TwoD, T1>& lhs,
                    const Array<TwoD, T2>& rhs)
    {
        if ( (lhs.GetRows() != rhs.GetRows()) ||
            (lhs.GetColumns() != rhs.GetColumns()) )
        {
            return false;
        }

        if ( lhs.data() == rhs.data() )
        {
            return true;
        }

        for (typename Array<OneD, T1>::size_type i = 0; i < lhs.GetRows(); ++i)
        {
            for (typename Array<OneD, T1>::size_type j = 0;
                j < lhs.GetColumns(); ++j)
            {
                if (!IsEqual(lhs[i][j], rhs[i][j]))
                {
                    return false;
                }
            }
        }

        return true;
    }

    template<typename T1, typename T2>
    bool operator!=(const Array<TwoD, T1>& lhs,
                    const Array<TwoD, T2>& rhs)
    {
        return !(lhs == rhs);
    }

    namespace LibUtilities
    {
        static std::vector<NekDouble> NullNekDoubleVector;
        static std::vector<unsigned int> NullUnsignedIntVector;
        static std::vector<std::vector<NekDouble> > NullVectorNekDoubleVector
            = { NullNekDoubleVector };
    }
}

#endif //NEKTAR_LIB_UTILITIES_BASIC_UTILS_SHARED_ARRAY_HPP
