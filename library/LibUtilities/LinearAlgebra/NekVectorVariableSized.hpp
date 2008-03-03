///////////////////////////////////////////////////////////////////////////////
//
// File: NekVectorVariableSized.hpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_VECTOR_VARIABLE_SIZED_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_VECTOR_VARIABLE_SIZED_HPP

//#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
//#include <LibUtilities/LinearAlgebra/NekPoint.hpp>
//#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>
//#include <LibUtilities/LinearAlgebra/NekVectorMetadata.hpp>

#include <LibUtilities/LinearAlgebra/PointerWrapper.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/Memory/DeleteNothing.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>

//#include <functional>
//#include <algorithm>
//#include <math.h>
//
//#include <boost/call_traits.hpp>
//#include <boost/type_traits.hpp>
//#include <boost/shared_array.hpp>


namespace Nektar
{
//     template<typename MatrixDataType, typename DataType, typename StorageType, typename Type, unsigned int space>
//     class BinaryExpressionTraits<NekMatrix<MatrixDataType, StorageType, Type>, NekVector<DataType, 0, space>, MultiplyOp>
//     {
//         public:
//             typedef NekVector<DataType, 0, space> ResultType;
//     };
        
    
    template<typename DataType, typename space>
    class NekVector<const DataType, VariableSizedVector, space>
    {
        public:
            /// \brief Creates an empty vector.
            NekVector() :
                m_size(0),
                m_data(),
                m_wrapperType(eCopy)
            {
            }
            
            /// \brief Creates a vector of given size.  The elements are not initialized.
            explicit NekVector(unsigned int size) :
                m_size(size),
                m_data(size),
                m_wrapperType(eCopy)
            {
            }

            /// \brief Creates a vector with given size and initial value.
            NekVector(unsigned int size, typename boost::call_traits<DataType>::const_reference a) :
                m_size(size),
                m_data(size),
                m_wrapperType(eCopy)
            {
                std::fill_n(m_data.get(), m_size, a);
            }


            explicit NekVector(const std::string& vectorValues) :
                m_size(0),
                m_data(),
                m_wrapperType(eCopy)
            {
                try
                {
                    std::vector<DataType> values = FromString<DataType>(vectorValues);
                    m_size = values.size();
                    m_data = Array<OneD, DataType>(m_size);
                    std::copy(values.begin(), values.end(), m_data);

                    ASSERTL0(m_size > 0, "Error converting string values to vector");
                }
                catch(std::runtime_error& e)
                {
                    NEKERROR(ErrorUtil::efatal, e.what());
                }
            }

            NekVector(typename boost::call_traits<DataType>::const_reference x,
                      typename boost::call_traits<DataType>::const_reference y,
                      typename boost::call_traits<DataType>::const_reference z) :
                m_size(3),
                m_data(m_size),
                m_wrapperType(eCopy)
            {
                m_data[0] = x;
                m_data[1] = y;
                m_data[2] = z;
            }
            
            NekVector(const NekVector<const DataType, VariableSizedVector, space>& rhs) :
                m_size(rhs.GetDimension()),
                m_data(rhs.m_data),
                m_wrapperType(rhs.m_wrapperType)
            {
                if( m_wrapperType = eCopy )
                {
                    m_data = Array<OneD, DataType>(m_size);
                    std::copy(rhs.begin(), rhs.end(), m_data.get());
                }
            }
            
//            explicit NekVector(const ConstArray<OneD, DataType>& ptr) :
//                m_size(ptr.num_elements()),
//                m_data(m_size),
//                m_wrapperType(eCopy)
//            {
//                CopyArray(ptr, m_data);
//            }

//            NekVector(unsigned int size, const ConstArray<OneD, DataType>& ptr) :
//                m_size(size),
//                m_data(size),
//                m_wrapperType(eCopy)
//            {
//                ASSERTL0(size <= ptr.num_elements(), "Attempting to populate a vector of size " +
//                    boost::lexical_cast<std::string>(size) + " but the incoming array only has " +
//                    boost::lexical_cast<std::string>(ptr.num_elements()) + " elements.");
//
//                std::copy(ptr.begin(), ptr.begin()+size, m_data.begin());
//            }
            
            NekVector(unsigned int size, const DataType* const ptr) :
                m_size(size),
                m_data(size, ptr),
                m_wrapperType(eCopy)
            {
            }
            
            explicit NekVector(const Array<OneD, DataType>& ptr, PointerWrapper h = eCopy) :
                m_size(ptr.num_elements()),
                m_data(ptr),
                m_wrapperType(h)
            {
                if( h == eCopy )
                {
                    m_data = Array<OneD, DataType>(m_size);
                    CopyArray(ptr, m_data);
                }
            }

            NekVector(unsigned int size, Array<OneD, DataType>& ptr, PointerWrapper h = eCopy) :
                m_size(size),
                m_data(ptr),
                m_wrapperType(h)
            {
                if( h == eCopy )
                {
                    ASSERTL0(size <= ptr.num_elements(), "Attempting to populate a vector of size " +
                        boost::lexical_cast<std::string>(size) + " but the incoming array only has " +
                        boost::lexical_cast<std::string>(ptr.num_elements()) + " elements.");

                    m_data = Array<OneD, DataType>(size);
                    std::copy(ptr.begin(), ptr.begin()+size, m_data.begin());
                }
            }
            
            explicit NekVector(const ConstArray<OneD, DataType>& ptr, PointerWrapper h = eCopy) :
                m_size(ptr.num_elements()),
                m_data(ptr, eVECTOR_WRAPPER),
                m_wrapperType(h)
            {
                if( h == eCopy )
                {
                    m_data = Array<OneD, DataType>(m_size);
                    CopyArray(ptr, m_data);
                }
            }

            NekVector(unsigned int size, const ConstArray<OneD, DataType>& ptr, PointerWrapper h = eCopy) :
                m_size(size),
                m_data(ptr, eVECTOR_WRAPPER),
                m_wrapperType(h)
            {
                if( h == eCopy )
                {
                    ASSERTL0(size <= ptr.num_elements(), "Attempting to populate a vector of size " +
                        boost::lexical_cast<std::string>(size) + " but the incoming array only has " +
                        boost::lexical_cast<std::string>(ptr.num_elements()) + " elements.");

                    m_data = Array<OneD, DataType>(size);
                    std::copy(ptr.begin(), ptr.begin()+size, m_data.begin());
                }
            }

            
            #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
            template<typename ExpressionPolicyType>
            NekVector(const Expression<ExpressionPolicyType>& rhs) :
                m_size(rhs.GetMetadata().Rows),
                m_data(rhs.GetMetadata().Rows),
                m_wrapperType(eCopy)
            {
                /// TODO Make sure this works correctly with eWrapper
                BOOST_MPL_ASSERT(( boost::is_same<typename Expression<ExpressionPolicyType>::ResultType, NekVector<const DataType, VariableSizedVector, space> > ));
                rhs.Evaluate(*this);
            }
            #endif //NEKTAR_USE_EXPRESSION_TEMPLATES

            virtual ~NekVector() {}

            NekVector<const DataType, VariableSizedVector, space>& operator=(const NekVector<const DataType, VariableSizedVector, space>& rhs)
            {
                if( m_wrapperType == eCopy  )
                {
                    // If the current vector is a copy, then regardless of the rhs type 
                    // we just copy over the values, resizing if needed.
                    if( GetDimension() != rhs.GetDimension() )
                    {
                        m_size = rhs.GetDimension();
                        m_data = Array<OneD, DataType>(m_size);
                    }
                }
                else if( m_wrapperType == eWrapper )
                {
                    // If the current vector is wrapped, then just copy over the top,
                    // but the sizes of the two vectors must be the same.
                    ASSERTL0(GetDimension() == rhs.GetDimension(), "Wrapped NekVectors must have the same dimension in operator=");
                }

                std::copy(rhs.begin(), rhs.end(), m_data.get());
                return *this;
            }
            
            /// \brief Returns the number of dimensions for the point.
            inline unsigned int GetDimension() const
            {
                return m_size;
            }

            inline unsigned int GetRows() const
            {
                return m_size;
            }

            
            const DataType* GetRawPtr() const
            {
                return m_data.get();
            }
     
            const ConstArray<OneD, DataType>& GetPtr() const { return m_data; }
     
            typedef const DataType* const_iterator;
            const_iterator begin() const { return GetRawPtr(); }
            const_iterator end() const { return GetRawPtr() + GetDimension(); }
             
            typename boost::call_traits<DataType>::const_reference operator()(unsigned int i) const
            {
                ASSERTL1(( i >= 0) && (i < GetDimension()), "Invalid access to m_data via parenthesis operator");
                return m_data[i];
            }

            typename boost::call_traits<DataType>::const_reference operator[](unsigned int i) const
            {
                return m_data[i];
            }

            typename boost::call_traits<DataType>::const_reference x() const
            {
                ASSERTL1( GetDimension() >= 1, "Invalid use of NekVector::x");
                return (*this)(0);
            }

            typename boost::call_traits<DataType>::const_reference y() const
            {
                ASSERTL1( GetDimension() >= 2, "Invalid use of NekVector::y");
                return (*this)(1);
            }

            typename boost::call_traits<DataType>::const_reference z() const
            {
                ASSERTL1( GetDimension() >= 3, "Invalid use of NekVector::z");
                return (*this)(2);
            }

            NekVector<const DataType, VariableSizedVector, space> operator-() const { return Negate(*this); }

            DataType Magnitude() const { return Nektar::Magnitude(*this); }

            DataType Dot(const NekVector<DataType, VariableSizedVector, space>& rhs) const { return Nektar::Dot(*this, rhs); }

            NekVector<DataType, VariableSizedVector, space> Cross(const NekVector<DataType, VariableSizedVector, space>& rhs) const
            {
                return Nektar::Cross(*this, rhs);
            }
            
            std::string AsString() const { return Nektar::AsString(*this); }

            // Norms
            DataType L1Norm() const { return Nektar::L1Norm(*this); }
            DataType L2Norm() const { return Nektar::L2Norm(*this); }
            DataType InfinityNorm() const { return Nektar::InfinityNorm(*this); }
            
                         
        protected:
            NekVector(const NekVectorMetadata& m) :
                m_size(m.Rows),
                m_data(m.Rows),
                m_wrapperType(eCopy)
            {
            }
            
            inline Array<OneD, DataType>& GetData() { return m_data; }
            inline void SetSize(unsigned int s) { m_size = s; }
            inline void SetWrapperType(PointerWrapper p) { m_wrapperType = p; }
            inline void SetData(const Array<OneD, DataType>& newData) { m_data = newData; }
            
        private:
            unsigned int m_size;
            Array<OneD, DataType> m_data;
            PointerWrapper m_wrapperType;
    };
    
    // \param DataType The type of data held by each element of the vector.
    // \param dim The number of elements in the vector.  If set to 0, the vector
    //            will have a variable number of elements.
    // \param space The space of the vector.
    template<typename DataType, typename space>
    class NekVector<DataType, VariableSizedVector, space> : public NekVector<const DataType, VariableSizedVector, space>
    {
        public:
            typedef NekVector<const DataType, VariableSizedVector, space> BaseType;
            
        public:
            /// \brief Creates an empty vector.
            NekVector() : BaseType() {}

            /// \brief Creates a vector of given size.  The elements are not initialized.
            explicit NekVector(unsigned int size) : BaseType(size) {}
                
            /// \brief Creates a vector with given size and initial value.
            NekVector(unsigned int size, typename boost::call_traits<DataType>::const_reference a) :
                BaseType(size, a) {}
                

            explicit NekVector(const std::string& vectorValues) :
                BaseType(vectorValues) {}
                
            NekVector(typename boost::call_traits<DataType>::const_reference x,
                      typename boost::call_traits<DataType>::const_reference y,
                      typename boost::call_traits<DataType>::const_reference z) :
                      BaseType(x, y, z) {}

            NekVector(const NekVector<DataType, VariableSizedVector, space>& rhs) :
                BaseType(rhs) {}
                
            explicit NekVector(const ConstArray<OneD, DataType>& ptr) :
                BaseType(ptr) {}
                

            NekVector(unsigned int size, const ConstArray<OneD, DataType>& ptr) :
                BaseType(size, ptr) {}
                            
            NekVector(unsigned int size, const DataType* const ptr) :
                BaseType(size, ptr) {}
                
            explicit NekVector(const Array<OneD, DataType>& ptr, PointerWrapper h = eCopy) :
                BaseType(ptr, h) {}
            
            NekVector(unsigned int size, const Array<OneD, DataType>& ptr, PointerWrapper h = eCopy) :
                BaseType(size, ptr, h) {}
                
            #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
            template<typename ExpressionPolicyType>
            NekVector(const Expression<ExpressionPolicyType>& rhs) :
                BaseType(rhs.GetMetadata()) 
            {
                /// TODO Make sure this works correctly with eWrapper
                BOOST_MPL_ASSERT(( boost::is_same<typename Expression<ExpressionPolicyType>::ResultType, NekVector<DataType, VariableSizedVector, space> > ));
                rhs.Evaluate(*this);
            }
            #endif //NEKTAR_USE_EXPRESSION_TEMPLATES
                    
            virtual ~NekVector() {}

            #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
            template<typename ExpressionPolicyType>
            NekVector<DataType, VariableSizedVector, space>& operator=(const Expression<ExpressionPolicyType>& rhs)
            {
                /// TODO Make sure this works correctly with eWrapper.
                BOOST_MPL_ASSERT(( boost::is_same<typename Expression<ExpressionPolicyType>::ResultType, NekVector<DataType, VariableSizedVector, space> > ));

                this->SetSize(rhs.GetMetadata().Rows);
                if( this->GetData().num_elements() != this->GetDimension() )
                {
                    this->SetData(Array<OneD, DataType>(this->GetDimension()));
                }
                this->SetWrapperType(eCopy);

                rhs.Evaluate(*this);
                return *this;            
            }
            #endif //NEKTAR_USE_EXPRESSION_TEMPLATES

            NekVector<DataType, VariableSizedVector, space>& operator=(const NekVector<DataType, VariableSizedVector, space>& rhs)
            {
                BaseType::operator=(rhs);
                return *this;
            }

            using BaseType::GetRawPtr;
            DataType* GetRawPtr()
            {
                return this->GetData().get();
            }
            
            using BaseType::GetPtr;
            Array<OneD, DataType>& GetPtr() { return this->GetData(); }
     
            typedef DataType* iterator;
            
            using BaseType::begin;
            using BaseType::end;
            
            iterator begin() { return GetRawPtr(); }
            iterator end() { return GetRawPtr() + this->GetDimension(); }
            
            using BaseType::operator();
            
            /// \brief Returns i^{th} element.
            /// \param i The element to return.
            /// \pre i < dim
            /// \return A reference to the i^{th} element.
            ///
            /// Retrieves the i^{th} element.  Since it returns a reference you may
            /// assign a new value (i.e., p(2) = 3.2;)
            ///
            /// This operator performs range checking.
            typename boost::call_traits<DataType>::reference operator()(unsigned int i)
            {
                ASSERTL1((i >= 0) && (i < GetDimension()), "Invalid access to m_data via parenthesis operator");
                return this->GetData()[i];
            }

            using BaseType::operator[];
            typename boost::call_traits<DataType>::reference operator[](unsigned int i)
            {
                return this->GetData()[i];
            }

            typename boost::call_traits<DataType>::reference x()
            {
                ASSERTL1(GetDimension() >= 1, "Invalid use of NekVector::x");
                return (*this)(0);
            }

            typename boost::call_traits<DataType>::reference y()
            {
                ASSERTL1(GetDimension() >= 2, "Invalid use of NekVector::y");
                return (*this)(1);
            }

            typename boost::call_traits<DataType>::reference z()
            {
                ASSERTL1(GetDimension() >= 3, "Invalid use of NekVector::z");
                return (*this)(2);
            }

            void SetX(typename boost::call_traits<DataType>::const_reference val)
            {
                ASSERTL1(GetDimension() >= 1, "Invalid use of NekVector::SetX");
                this->GetData()[0] = val;
            }

            void SetY(typename boost::call_traits<DataType>::const_reference val)
            {
                ASSERTL1(GetDimension() >= 2, "Invalid use of NekVector::SetX");
                this->GetData()[1] = val;
            }

            void SetZ(typename boost::call_traits<DataType>::const_reference val)
            {
                ASSERTL1(GetDimension() >= 3, "Invalid use of NekVector::SetX");
                this->GetData()[2] = val;
            }

            /// Arithmetic Routines

            // Unitary operators

            NekVector<DataType, VariableSizedVector, space>& operator+=(const NekVector<DataType, VariableSizedVector, space>& rhs)
            {
                PlusEqual(*this, rhs);
                return *this;
            }

            NekVector<DataType, VariableSizedVector, space>& operator-=(const NekVector<DataType, VariableSizedVector, space>& rhs)
            {
                MinusEqual(*this, rhs);
                return *this;
            }

            NekVector<DataType, VariableSizedVector, space>& operator*=(typename boost::call_traits<DataType>::const_reference rhs)
            {
                TimesEqual(*this, rhs);
                return *this;
            }
            
            NekVector<DataType, VariableSizedVector, space>& operator/=(typename boost::call_traits<DataType>::const_reference rhs)
            {
                DivideEqual(*this, rhs);
                return *this;
            }

            
            void Normalize() { return Nektar::Normalize(*this); }
            
            
        private:
            // Prevents accidental use of wrapped mode around ConstArrays.
            NekVector(const ConstArray<OneD, DataType>& ptr, PointerWrapper h);

    };    
    
}

#endif // NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_VECTOR_VARIABLE_SIZED_HPP
