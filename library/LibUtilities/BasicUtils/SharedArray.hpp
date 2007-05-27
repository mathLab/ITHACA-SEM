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

#include <boost/multi_array.hpp>
#include <boost/shared_ptr.hpp>

namespace Nektar
{
    
    /// \page Arrays
    ///
    /// An Array object is a thin wrapper around native arrays.  Arrays provide 
    /// all the functionality of native arrays, with the additional benefits of 
    /// automatic use of the Nektar memory pool, automatic memory allocation and 
    /// deallocation, and easy creation and use of multi-dimensional arrays.
    ///
    /// The two Array classes are ConstArray and Array.  Each takes two template parameters.
    /// The first is a Dimension enumeration and indicates the array's dimensionality.
    /// The second is the type of data held in the array.
    ///
    /// A ConstArray provides constant access to an array.  A ConstArray can, through 
    /// its assignment operator and copy constructor, access many different arrays, but 
    /// it is impossible to modify the data.
    ///
    /// NumDimensions
    ///
    /// DataType - Requirements of DataType.
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
    /// const ConstArray<OneD, NekDouble>& moreEfficient = points->GetZ();
    ///\endcode
    ///
    /// If you know the size, don't create an empty array and then populate it on another line.
    ///
    ///\section ArrayCreation Array Creation
    ///
    /// Each array comes with three types of constructors.
    ///
    ///\li Empty Array Constructor.  This constructor takes one parameter for each dimension, then creates 
    /// an empty array of the appropriate size.  If the array data type is a fundamental C++ type
    /// (such as double, int, char, etc.) then the array will be uninitialized.  If the array is any 
    /// other type each element will be initialized with the data type's default constructor.  
    ///
    ///\code
    /// class Double 
    /// {
    ///     public:
    ///         Double() : m_value(0.0) {}
    ///     private:
    ///         double m_value;
    /// };
    ///
    /// class Integer
    /// {
    ///     public: 
    ///         Integer(int i) : m_value(i) {}
    ///     private:
    ///         int m_value;
    /// };
    ///
    /// void ExampleFunction() 
    /// {
    ///     // Creates a 10 element uninitialized array.
    ///     Array<OneD, NekDouble> a(10);
    ///
    ///     // Creates a 10 element array initialized with the Double default constructor.
    ///     // For each element i in the array, b[i].m_value == 0.0
    ///     // The default constructor is called 10 times.
    ///     Array<OneD, Double> b(10);
    ///
    ///     // Creates a 10x20 array with a total of 200 elements using Double's default constructor.
    ///     Array<TwoD, Double> c(10, 20);
    ///
    ///     // This is an error.  Class Integer doesn't have a default constructor.
    ///     Array<OneD, Integer> c(15);
    /// }
    ///\endcode
    ///
    ///\li Initial Value Constructor.  This constructor takes the apprpriate number of size parameters 
    /// for the array's dimension, then an additional parameter which specifies the initial value for 
    /// each element.  If the array's data type is a fundamental type (such as double, int, etc.), then
    /// the initial value is copied directly into the array.  Otherwise, each element is created using the 
    /// data type's copy constructor.
    ///
    ///\code
    /// class Double 
    /// {
    ///     public:
    ///         Double(double v) : m_value(v) {}
    ///         Double(const Double& rhs) : m_value(rhs.m_value) {}
    ///     private:
    ///         double m_value;
    /// };
    ///
    /// class Integer
    /// {
    ///     public: 
    ///         Integer(int i) : m_value(i) {}
    ///     private:
    ///         Integer(const Integer& rhs);
    ///         int m_value;
    /// };
    ///
    /// void ExampleInitialValueConstructors() 
    /// {
    ///     // Copies the value 8.9 into each a's 10 elements.
    ///     Array<OneD, NekDouble> a(10, 8.9);
    ///
    ///     // Copies d into each of b's 10 elements.
    ///     Double d(1.7);
    ///     Array<OneD, Double> b(10, d);
    ///
    ///     // Does not compile since Integer does not have a public copy constructor.
    ///     Integer i(2);
    ///     Array<OneD, Integer> c(10, i);
    /// }
    ///\endcode
    ///
    
    
    
    // Forward declaration for a ConstArray constructor.
    template<Dimension Dim, typename DataType>
    class Array;    
    
    template<Dimension Dim, typename DataType>
    class ConstArray;
    
    /// \brief A 1D immutable array.  See \ref Arrays for details.
    template<typename DataType>
    class ConstArray<OneD, DataType>
    {
        public:
            typedef boost::multi_array_ref<DataType, 1> ArrayType;
            typedef typename ArrayType::const_reference const_reference;
            typedef typename ArrayType::reference reference;
            
            typedef typename ArrayType::index index;
            typedef typename ArrayType::const_iterator const_iterator;
            typedef typename ArrayType::iterator iterator;
            
            typedef typename ArrayType::element element;
            typedef typename ArrayType::size_type size_type;
            
            
        public:
            /// \brief Creates an empty, size 0 array.
            ConstArray() :
                m_data(CreateStorage<DataType>(0)),
                m_offset(0)
            {
            }
            
            /// \brief Constructs an empty 1 dimensional array.
            /// \param dim1Size The array's size.
            ///
            /// If DataType is a fundamental type (double, int, etc.), then the allocated array is 
            /// uninitialized.  If it is any other type, each element is initialized with DataType's default
            /// constructor.
            explicit ConstArray(unsigned int dim1Size) :
                m_data(CreateStorage<DataType>(dim1Size)),
                m_offset(0)
            {
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements());
            }
            
            /// \brief Creates a 1D array with each element initialized to an initial value.
            /// \param dim1Size The array's size.
            /// \param initValue Each element's initial value.
            ///
            /// If DataType is a fundamental type (double, int, etc.), then the initial value 
            /// is copied directly into each element.  Otherwise, the DataType's copy constructor
            /// is used to initialize each element.
            ConstArray(unsigned int dim1Size, const DataType& initValue) :
                m_data(CreateStorage<DataType>(dim1Size)),
                m_offset(0)
            {
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements(), initValue);
            }
            
            /// \brief Creates a 1D array with a copy the given raw array.
            /// \param dim1Size the array's size.
            /// \param data The data to copy.  
            ///
            /// If DataType is a fundamental type (double, int, etc.), then data is copied
            /// directly into the underlying storage.  Otherwise, the DataType's copy constructor
            /// is used to copy each element.
            ConstArray(unsigned int dim1Size, const DataType* data) :
                m_data(CreateStorage<DataType>(dim1Size)),
                m_offset(0)
            {
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements(), data);
            }
            
            ConstArray(const ConstArray<OneD, DataType>& rhs) :
                m_data(rhs.m_data),
                m_offset(rhs.m_offset)
            {
            }
            
            ConstArray<OneD, DataType>& operator=(const ConstArray<OneD, DataType>& rhs)
            {
                m_data = rhs.m_data;
                m_offset = rhs.m_offset;
                return *this;
            }
            
            const_iterator begin() const { return m_data->begin() + m_offset; }
            const_iterator end() const { return m_data->end(); }
            
            const_reference operator[](index i) const 
            {
                ASSERTL1(i < num_elements(), (std::string("Element ") +
                    boost::lexical_cast<std::string>(i) + std::string(" requested in an array of size ") +
                    boost::lexical_cast<std::string>(num_elements())));
                return (*m_data)[i+m_offset]; 
            }
            
            const element* get() const { return m_data->data()+m_offset; }
            const element* data() const { return m_data->data()+m_offset; }
            size_type num_dimensions() const { return m_data->num_dimensions(); }
            //const size_type* shape() const { return m_data->shape(); }
            size_type num_elements() const { return m_data->num_elements()-m_offset; }
            unsigned int GetOffset() const { return m_offset; }
            
            template<typename T>
            friend bool operator==(const ConstArray<OneD, T>&, const ConstArray<OneD, T>&);
            
            template<typename T>
            friend ConstArray<OneD, T> operator+(const ConstArray<OneD, T>& lhs, unsigned int offset);
   
            template<typename T>
            friend ConstArray<OneD, T> operator+(unsigned int offset, const ConstArray<OneD, T>& rhs);
            
        protected:
            boost::shared_ptr<ArrayType> m_data;
            unsigned int m_offset;
            
        private:
            static ConstArray<OneD, DataType> CreateWithOffset(const ConstArray<OneD, DataType>& rhs, unsigned int offset)
            {
                ConstArray<OneD, DataType> result(rhs);
                result.m_offset = offset;
                return result;
            }
            
    };
    
    template<typename DataType>
    class ConstArray<TwoD, DataType>
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
            ConstArray() :
                m_data(CreateStorage<DataType>(0, 0))
            {
            }
            
            /// \brief Constructs a 3 dimensional array.  The elements of the array are not initialized.
            ConstArray(unsigned int dim1Size, unsigned int dim2Size) :
                m_data(CreateStorage<DataType>(dim1Size, dim2Size))
            {
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements());
            }
            
            ConstArray(unsigned int dim1Size, unsigned int dim2Size, const DataType& initValue) :
                m_data(CreateStorage<DataType>(dim1Size, dim2Size))
            {
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements(), initValue);
            }
            
            ConstArray(unsigned int dim1Size, unsigned int dim2Size, const DataType* data) :
                m_data(CreateStorage<DataType>(dim1Size, dim2Size))
            {
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements(), data);
            }
            
            ConstArray(const ConstArray<TwoD, DataType>& rhs) :
                m_data(rhs.m_data)
            {
            }
            
            ConstArray<TwoD, DataType>& operator=(const ConstArray<TwoD, DataType>& rhs)
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

            size_type GetRows() const { return m_data->shape()[0]; }
            size_type GetColumns() const { return m_data->shape()[1]; }
            
        protected:
            boost::shared_ptr<ArrayType> m_data;
            
        private:
            
    };
    
    template<typename DataType>
    class ConstArray<ThreeD, DataType>
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
            ConstArray() :
                m_data(CreateStorage<DataType>(0, 0, 0))
            {
            }
            
            /// \brief Constructs a 3 dimensional array.  The elements of the array are not initialized.
            ConstArray(unsigned int dim1Size, unsigned int dim2Size, unsigned int dim3Size) :
                m_data(CreateStorage<DataType>(dim1Size, dim2Size, dim3Size))
            {
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements());
            }
            
            ConstArray(unsigned int dim1Size, unsigned int dim2Size, unsigned int dim3Size, const DataType& initValue) :
                m_data(CreateStorage<DataType>(dim1Size, dim2Size, dim3Size))
            {
                ArrayInitializationPolicy<DataType>::Initialize(m_data->data(), m_data->num_elements(), initValue);
            }
            
            ConstArray(const ConstArray<ThreeD, DataType>& rhs) :
                m_data(rhs.m_data)
            {
            }
            
            ConstArray<ThreeD, DataType>& operator=(const ConstArray<ThreeD, DataType>& rhs)
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

    template<Dimension Dim, typename DataType>
    class Array;
    
    
    /// \brief 1D Array
    ///
    /// \ref Arrayp
    ///
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
            
            Array(unsigned int dim1Size) :
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
            
            
            Array(const Array<TwoD, DataType>& rhs) :
                BaseType(rhs)
            {
            }
            
            Array(const ConstArray<TwoD, DataType>& rhs) :
                BaseType(rhs)
            {
            }
            
            Array<TwoD, DataType>& operator=(const Array<TwoD, DataType>& rhs)
            {
                BaseType::operator=(rhs);
                return *this;
            }
            
            static Array<OneD, DataType> CreateWithOffset(const Array<OneD, DataType>& rhs, unsigned int offset)
            {
                Array<OneD, DataType> result(rhs);
                result.m_offset = offset;
                return result;
            }
            
            using BaseType::begin;
            iterator begin() { return this->m_data->begin()+this->m_offset; }
            
            using BaseType::end;
            iterator end() { return this->m_data->end(); }
            
            using BaseType::operator[];
            reference operator[](index i)
            {
                ASSERTL1(i < this->num_elements(), (std::string("Element ") +
                    boost::lexical_cast<std::string>(i) + std::string(" requested in an array of size ") +
                    boost::lexical_cast<std::string>(this->num_elements())));
                return (*this->m_data)[i+this->m_offset]; 
            }
                
            
            using BaseType::get;
            element* get() { return this->m_data->data()+this->m_offset; }
            
            using BaseType::data;
            element* data() { return this->m_data->data()+this->m_offset; }
            
        private:
            
    };

    /// \brief A 2D array.
    template<typename DataType>
    class Array<TwoD, DataType> : public ConstArray<TwoD, DataType>
    {
        public:
            typedef ConstArray<TwoD, DataType> BaseType;
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
    class Array<ThreeD, DataType> : public ConstArray<ThreeD, DataType>
    {
        public:
            typedef ConstArray<ThreeD, DataType> BaseType;
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
        return ConstArray<OneD, DataType>::CreateWithOffset(lhs, offset);
    }

    template<typename DataType>
    ConstArray<OneD, DataType> operator+(unsigned int offset, const ConstArray<OneD, DataType>& rhs)
    {
        return ConstArray<OneD, DataType>::CreateWithOffset(rhs, offset);
    }
    
    template<typename DataType>
    Array<OneD, DataType> operator+(const Array<OneD, DataType>& lhs, unsigned int offset)
    {
        return Array<OneD, DataType>::CreateWithOffset(lhs, offset);
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
