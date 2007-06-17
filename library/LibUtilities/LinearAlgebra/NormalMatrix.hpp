///////////////////////////////////////////////////////////////////////////////
//
// File: NormalMatrix.hpp
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
// Description: Interface classes for matrices
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NORMAL_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NORMAL_MATRIX_HPP

#include <LibUtilities/LinearAlgebra/MatrixBase.hpp>
#include <LibUtilities/LinearAlgebra/FullMatrixStoragePolicy.hpp>
#include <LibUtilities/LinearAlgebra/DiagonalMatrixStoragePolicy.hpp>
#include <LibUtilities/BasicUtils/OperatorGenerators.hpp>
#include <LibUtilities/BasicUtils/BinaryExpressionTraits.hpp>

namespace Nektar
{
    // Define the expression traits for a matrix before defining the matrix or the 
    // OperatorGenerator won't work.
    
    /// \brief Gives all expressions involving the same storage types the same result type.
    template<typename DataType, typename LhsStorageType, typename RhsStorageType, typename OpType>
    class BinaryExpressionTraits<NekMatrix<DataType, LhsStorageType, StandardMatrixTag>, NekMatrix<DataType, RhsStorageType, StandardMatrixTag>, OpType>
    {
        public:
            typedef NekMatrix<DataType, FullMatrixTag, StandardMatrixTag> ResultType;
    };
    
    template<typename DataType, typename StorageType, typename OpType>
    class BinaryExpressionTraits<NekMatrix<DataType, StorageType, StandardMatrixTag>, NekMatrix<DataType, StorageType, StandardMatrixTag>, OpType>
    {
        public:
            typedef NekMatrix<DataType, StorageType, StandardMatrixTag> ResultType;
    };
    
    
    template<typename DataType, typename StorageType>
    class NekMatrix<DataType, StorageType, StandardMatrixTag> : public Matrix<DataType>,
                                                                public OperatorGeneratorL3R3<NekMatrix, NekMatrix, AddOp>,
                                                                public OperatorGeneratorL3R3<NekMatrix, NekMatrix, MultiplyOp>
    {
        public:
            typedef Matrix<DataType> BaseType;
            typedef NekMatrix<DataType, StorageType, StandardMatrixTag> ThisType;
            typedef MatrixStoragePolicy<DataType, StorageType> StoragePolicy;
            
            typedef typename StoragePolicy::GetValueReturnType GetValueType;
            typedef typename boost::call_traits<DataType>::const_reference ConstGetValueType;
            
        public:
            NekMatrix() :
                BaseType(0, 0),
                m_data(StoragePolicy::Initialize()),
                m_wrapperType(eCopy)
            {
            }
            
            NekMatrix(unsigned int rows, unsigned int columns) :
                BaseType(rows, columns),
                m_data(StoragePolicy::Initialize(rows, columns)),
                m_wrapperType(eCopy)
            {
            }
            
            NekMatrix(unsigned int rows, unsigned int columns, typename boost::call_traits<DataType>::const_reference initValue) :
                BaseType(rows, columns),
                m_data(StoragePolicy::Initialize(rows, columns, initValue)),
                m_wrapperType(eCopy)
            {
            }
            
            NekMatrix(unsigned int rows, unsigned int columns, const DataType* data) :
                BaseType(rows, columns),
                m_data(StoragePolicy::Initialize(rows, columns, data)),
                m_wrapperType(eCopy)
            {
            }
            
            NekMatrix(unsigned int rows, unsigned int columns, const ConstArray<OneD, DataType>& d) :
                BaseType(rows, columns),
                m_data(StoragePolicy::Initialize(rows, columns, d)),
                m_wrapperType(eCopy)
            {
            }
            
            NekMatrix(unsigned int rows, unsigned int columns, const Array<OneD, DataType>& d, PointerWrapper wrapperType = eCopy) :
                BaseType(rows, columns),
                m_data(wrapperType == eCopy ? StoragePolicy::Initialize(rows, columns, d) : d),
                m_wrapperType(wrapperType)
            {
            }
            
            NekMatrix(const ThisType& rhs) :
                BaseType(rhs),
                m_data(),
                m_wrapperType(rhs.m_wrapperType)
            {
                if( m_wrapperType == eCopy )
                {
                    m_data = StoragePolicy::Initialize(this->GetRows(), this->GetColumns());
                    CopyArray(rhs.m_data, m_data);
                }
                else
                {
                    m_data = rhs.m_data;
                }
            }
            
            MatrixStorage GetStorageType() const 
            {
                return static_cast<MatrixStorage>(ConvertToMatrixStorageEnum<StorageType>::Value);
            }
            
            // TODO - Copy constructors from other types of matrices.
            
            ThisType& operator=(const ThisType& rhs)
            {
                if( this != &rhs )
                {
                    BaseType::operator=(rhs);
                    if( m_wrapperType == eCopy  )
                    {
                        // If the current vector is a matrix, then regardless of the rhs type 
                        // we just copy over the values, resizing if needed.
                        if( m_data.num_elements() != rhs.m_data.num_elements() )
                        {
                            m_data = Array<OneD, DataType>(rhs.m_data.num_elements());
                        }
                    }
                    else if( m_wrapperType == eWrapper )
                    {
                        // If the current matrix is wrapped, then just copy over the top,
                        // but the sizes of the two matrices must be the same.
                        ASSERTL0(m_data.num_elements() == rhs.m_data.num_elements(), "Wrapped NekMatrices must have the same dimension in operator=");
                    }
    
                    CopyArray(rhs.m_data, m_data);
                }
                
                return *this;
            }
            
            
            ConstGetValueType operator()(unsigned int row, unsigned int column) const
            {
                ASSERTL2(row < this->GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetRows()) +
                    std::string(" rows"));
                ASSERTL2(column < this->GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetColumns()) +
                    std::string(" columns"));
                    
                return StoragePolicy::GetValue(this->GetRows(), this->GetColumns(), row, column, m_data);    
            }
            
            GetValueType operator()(unsigned int row, unsigned int column)
            {
                ASSERTL2(row < this->GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetRows()) +
                    std::string(" rows"));
                ASSERTL2(column < this->GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetColumns()) +
                    std::string(" columns"));
                    
                return StoragePolicy::GetValue(this->GetRows(), this->GetColumns(), row, column, m_data);    
            }

            void SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d)
            {
                ASSERTL2(row < this->GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetRows()) +
                    std::string(" rows"));
                ASSERTL2(column < this->GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetColumns()) +
                    std::string(" columns"));
                StoragePolicy::SetValue(this->GetRows(), this->GetColumns(), row, column, m_data, d);
            }
            
            typename boost::call_traits<DataType>::const_reference GetValue(unsigned int row, unsigned int column) const
            {
                ASSERTL2(row < this->GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetRows()) +
                    std::string(" rows"));
                ASSERTL2(column < this->GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetColumns()) +
                    std::string(" columns"));
                return StoragePolicy::GetValue(this->GetRows(), this->GetColumns(), row, column, m_data);
            }
            
            Array<OneD, DataType>& GetPtr()
            {
                return m_data;
            }
            
            const ConstArray<OneD, DataType>& GetPtr() const
            {
                return m_data;
            }
            
            DataType* GetRawPtr()
            {
                return m_data.data();
            }
            
            const DataType* GetRawPtr() const
            {
                return m_data.data();
            }
            
            typedef DataType* iterator;
            typedef const DataType* const_iterator;
            
            iterator begin() { return m_data.data(); }
            iterator end() { return m_data.data() + m_data.num_elements(); }
            
            const_iterator begin() const { return m_data.data(); }
            const_iterator end() const { return m_data.data() + m_data.num_elements(); }
            
            unsigned int GetStorageSize() const
            {
                return m_data.num_elements();
            }
            
            bool operator==(const NekMatrix<DataType, StorageType, StandardMatrixTag>& rhs) const
            {
                if( GetStorageSize() != rhs.GetStorageSize() )
                {
                    return false;
                }
                
                return std::equal(begin(), end(), rhs.begin());
            }
            
            
       protected:
            
            
        private:
            virtual typename boost::call_traits<DataType>::value_type v_GetValue(unsigned int row, unsigned int column) const 
            {
                return ThisType::operator()(row, column);
            }
            
            virtual unsigned int v_GetStorageSize() const 
            {
                return ThisType::GetStorageSize();
            }
            
            virtual MatrixStorage v_GetStorageType() const
            {
                return ThisType::GetStorageType();
            }
            
            virtual void v_SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d)
            {
                return ThisType::SetValue(row, column, d);
            }
            
            Array<OneD, DataType> m_data;
            PointerWrapper m_wrapperType;
    };
    

    //GENERATE_ADDITION_OPERATOR_L3R3(NekMatrix, NekMatrix);
    //GENERATE_MULTIPLICATION_OPERATOR_L3R3(NekMatrix, NekMatrix);
}


#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NORMAL_MATRIX_HPP
