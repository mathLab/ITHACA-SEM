///////////////////////////////////////////////////////////////////////////////
//
// File: StandardMatrix.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_STANDARD_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_STANDARD_MATRIX_HPP

#include <LibUtilities/LinearAlgebra/MatrixBase.hpp>
#include <LibUtilities/LinearAlgebra/PointerWrapper.h>
#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrixMetadata.hpp>
#include <LibUtilities/LinearAlgebra/MatrixFuncs.h>

namespace Nektar
{
    /// \brief Standard Matrix
    /// \param DataType The type stored in each element.
    ///
    /// Matrices are stored in column major order to make it easier to interoperate with 
    /// Blas and Lapack.
    template<typename DataType>
    class NekMatrix<const DataType, StandardMatrixTag> : public Matrix<DataType>
    {
        public:
            typedef Matrix<DataType> BaseType;
            typedef NekMatrix<const DataType, StandardMatrixTag> ThisType;
            typedef const DataType NumberType;
          
            public:
                            
                template<typename T, typename MatrixType>
                class iterator_impl
                {
                    public:
                        typedef T value_type;
                        typedef std::input_iterator_tag iterator_category;
                        typedef unsigned int difference_type;
                        typedef typename boost::call_traits<value_type>::reference reference;
                        typedef typename boost::call_traits<value_type>::const_reference const_reference;
                        typedef typename boost::remove_reference<value_type>::type* pointer;

                    public:
                        iterator_impl(pointer d, pointer e, bool isEnd = false) :
                            m_data(d),
                            m_end(e),
                            m_curRow(std::numeric_limits<unsigned int>::max()),
                            m_curColumn(std::numeric_limits<unsigned int>::max()),
                            m_matrix(NULL),
                            m_curIndex(std::numeric_limits<unsigned int>::max()),
                            m_transpose('N')
                        {
                            if( isEnd )
                            {
                                m_data = m_end;
                            }
                        }

                        iterator_impl(MatrixType* m, char transpose, bool isEnd = false) :
                            m_data(NULL),
                            m_end(NULL),
                            m_curRow(0),
                            m_curColumn(0),
                            m_matrix(m),
                            m_curIndex(0),
                            m_transpose(transpose)
                        {
                            if( isEnd )
                            {
                                m_curRow = std::numeric_limits<unsigned int>::max();
                                m_curColumn = std::numeric_limits<unsigned int>::max();
                                m_curIndex = std::numeric_limits<unsigned int>::max();
                            }
                        }
                        
                        iterator_impl(const iterator_impl<T, MatrixType>& rhs) :
                            m_data(rhs.m_data),
                            m_end(rhs.m_end),
                            m_curRow(rhs.m_curRow),
                            m_curColumn(rhs.m_curColumn),
                            m_matrix(rhs.m_matrix),
                            m_curIndex(rhs.m_curIndex),
                            m_transpose(rhs.m_transpose)
                        {
                        }

                        iterator_impl<T, MatrixType>& operator=(const iterator_impl<T, MatrixType>& rhs)
                        {
                            m_data = rhs.m_data;
                            m_end = rhs.m_end;
                            m_curRow = rhs.m_curRow;
                            m_curColumn = rhs.m_curColumn;
                            m_matrix = rhs.m_matrix;
                            m_curIndex = rhs.m_curIndex;
                            m_transpose = rhs.m_transpose;
                            return *this;
                        }

                        reference operator*()
                        {
                            if( m_data )
                            {
                                ASSERTL1(m_data < m_end, "Attempt to dereference matrix iterator after its end.");
                                return *m_data;
                            }
                            else
                            {
                                return m_matrix->GetPtr()[m_curIndex];
                            }
                        }

                        const_reference operator*() const
                        {
                            if( m_data )
                            {
                                ASSERTL1(m_data < m_end, "Attempt to dereference matrix iterator after its end.");
                                return *m_data;
                            }
                            else
                            {
                                return m_matrix->GetPtr(m_curIndex);
                            }
                        }

                        /// \brief Prefix increment operator.
                        iterator_impl<T, MatrixType>& operator++()
                        {
                            if( m_data )
                            {
                                ++m_data;
                            }
                            else
                            {
                                boost::tie(m_curRow, m_curColumn) = 
                                    m_matrix->Advance(m_curRow, m_curColumn, m_transpose);
                                if( m_curRow == std::numeric_limits<unsigned int>::max() )
                                {
                                    m_curIndex = m_curRow;
                                }
                                else
                                {
                                    m_curIndex = m_matrix->CalculateIndex(m_curRow, m_curColumn, m_transpose);
                                }
                            }
                            return *this;
                        }

                        /// \postfix increment operator.
                        iterator_impl<T, MatrixType> operator++(int)
                        {
                            iterator_impl<T, MatrixType> result = *this;
                            ++(*this);
                            return result;
                        }

                        bool operator==(const iterator_impl<T, MatrixType>& rhs)
                        {
                            return m_data == rhs.m_data &&
                                   m_end == rhs.m_end &&
                                   m_curRow == rhs.m_curRow &&
                                   m_curColumn == rhs.m_curColumn &&
                                   m_matrix == rhs.m_matrix &&
                                   m_curIndex == rhs.m_curIndex &&
                                   m_transpose == rhs.m_transpose;
                        }

                        bool operator!=(const iterator_impl<T, MatrixType>& rhs)
                        {
                            return !(*this == rhs);
                        }

                    private:
                        // Used when the matrix is not transposed
                        T* m_data;
                        T* m_end;

                        // Used when the matrix is transposed.
                        unsigned int m_curRow;
                        unsigned int m_curColumn;
                        MatrixType* m_matrix;
                        unsigned int m_curIndex;
                        char m_transpose;
            };

        public:
            /// \brief Creates an empty matrix.
            NekMatrix() :
                BaseType(0, 0),
                m_data(),
                m_wrapperType(eCopy),
                m_storagePolicy(eFULL),
                m_numberOfSuperDiagonals(std::numeric_limits<unsigned int>::max()),
                m_numberOfSubDiagonals(std::numeric_limits<unsigned int>::max())
            {
                m_data = Array<OneD, DataType>(GetRequiredStorageSize());
            }
            
            /// \brief Creates a rows by columns matrix.
            /// \brief rows The number of rows in the matrix.
            /// \brief columns The number of columns in the matrix.
            NekMatrix(unsigned int rows, unsigned int columns, MatrixStorage policy = eFULL,
                      unsigned int subDiagonals = std::numeric_limits<unsigned int>::max(),
                      unsigned int superDiagonals = std::numeric_limits<unsigned int>::max()) :
                BaseType(rows, columns),
                m_data(),
                m_wrapperType(eCopy),
                m_storagePolicy(policy),
                m_numberOfSuperDiagonals(superDiagonals),
                m_numberOfSubDiagonals(subDiagonals)
            {
                m_data = Array<OneD, DataType>(GetRequiredStorageSize());
            }

            /// \brief Creates a rows by columns matrix.
            /// \brief rows The number of rows in the matrix.
            /// \brief columns The number of columns in the matrix.
            /// \brief initValue The value used to initialize each element.
            /// \brief policySpecificData Data the is specific to the storage type assigned to this matrix.
            ///                           Check MatrixStoragePolicy<StorageType> to see if the StoragePolicy
            ///                           for this matrix has policy specific data.
            NekMatrix(unsigned int rows, unsigned int columns, typename boost::call_traits<DataType>::const_reference initValue,
                      MatrixStorage policy = eFULL,
                      unsigned int subDiagonals = std::numeric_limits<unsigned int>::max(),
                      unsigned int superDiagonals = std::numeric_limits<unsigned int>::max()) :
                BaseType(rows, columns),
                m_data(),
                m_wrapperType(eCopy),
                m_storagePolicy(policy),
                m_numberOfSuperDiagonals(superDiagonals),
                m_numberOfSubDiagonals(subDiagonals)
            {
                m_data = Array<OneD, DataType>(GetRequiredStorageSize());
                std::fill(m_data.begin(), m_data.end(), initValue);
            }
            
            /// \brief Creates a rows by columns matrix.
            /// \brief rows The number of rows in the matrix.
            /// \brief columns The number of columns in the matrix.
            /// \brief data An array of data use to initialize the matrix.
            /// \brief policySpecificData Data the is specific to the storage type assigned to this matrix.
            ///                           Check MatrixStoragePolicy<StorageType> to see if the StoragePolicy
            ///                           for this matrix has policy specific data.
            NekMatrix(unsigned int rows, unsigned int columns, const DataType* data,
                      MatrixStorage policy = eFULL,
                      unsigned int subDiagonals = std::numeric_limits<unsigned int>::max(),
                      unsigned int superDiagonals = std::numeric_limits<unsigned int>::max()) :
                BaseType(rows, columns),
                m_data(),
                m_wrapperType(eCopy),
                m_storagePolicy(policy),
                m_numberOfSuperDiagonals(superDiagonals),
                m_numberOfSubDiagonals(subDiagonals)
            {
                unsigned int size = GetRequiredStorageSize();
                m_data = Array<OneD, DataType>(size);
                std::copy(data, data + size, m_data.begin());
            }
            
            /// \brief Creates a rows by columns matrix.
            /// \brief rows The number of rows in the matrix.
            /// \brief columns The number of columns in the matrix.
            /// \brief d An array of data used to initialize the matrix.  Values from d are copied into the matrix.
            /// \brief policySpecificData Data the is specific to the storage type assigned to this matrix.
            ///                           Check MatrixStoragePolicy<StorageType> to see if the StoragePolicy
            ///                           for this matrix has policy specific data.
            NekMatrix(unsigned int rows, unsigned int columns, const Array<OneD, const DataType>& d,
                      MatrixStorage policy = eFULL,
                      unsigned int subDiagonals = std::numeric_limits<unsigned int>::max(),
                      unsigned int superDiagonals = std::numeric_limits<unsigned int>::max()) :
                BaseType(rows, columns),
                m_data(),
                m_wrapperType(eCopy),
                m_storagePolicy(policy),
                m_numberOfSuperDiagonals(superDiagonals),
                m_numberOfSubDiagonals(subDiagonals)
            {
                m_data = Array<OneD, DataType>(GetRequiredStorageSize());
                CopyArrayN(d, m_data, m_data.size());
            }
            
            /// \brief Creates a rows by columns matrix.
            /// \brief rows The number of rows in the matrix.
            /// \brief columns The number of columns in the matrix.
            /// \brief d An array of data used to initialize the matrix.  If wrapperType is eCopy, then
            ///          each element is copied from d to the matrix.  If wrapperType is eWrapper, then 
            ///          the matrix uses d directly as its matrix data and no copies are made.
            /// \brief policySpecificData Data the is specific to the storage type assigned to this matrix.
            ///                           Check MatrixStoragePolicy<StorageType> to see if the StoragePolicy
            ///                           for this matrix has policy specific data.
            NekMatrix(unsigned int rows, unsigned int columns, const Array<OneD, const DataType>& d, PointerWrapper wrapperType = eCopy,
                      MatrixStorage policy = eFULL,
                      unsigned int subDiagonals = std::numeric_limits<unsigned int>::max(),
                      unsigned int superDiagonals = std::numeric_limits<unsigned int>::max()) :
                BaseType(rows, columns),
                m_data(),
                m_wrapperType(wrapperType),
                m_storagePolicy(policy),
                m_numberOfSuperDiagonals(superDiagonals),
                m_numberOfSubDiagonals(subDiagonals)
            {
                if( wrapperType == eWrapper )
                {
                     m_data = Array<OneD, DataType>(d, eVECTOR_WRAPPER);
                }
                else
                {
                    m_data = Array<OneD, DataType>(GetRequiredStorageSize());
                    CopyArrayN(d, m_data, m_data.num_elements());
                }
            }
            
            NekMatrix(const ThisType& rhs) :
                BaseType(rhs),
                m_data(),
                m_wrapperType(rhs.m_wrapperType),
                m_storagePolicy(rhs.m_storagePolicy),
                m_numberOfSuperDiagonals(rhs.m_numberOfSuperDiagonals),
                m_numberOfSubDiagonals(rhs.m_numberOfSubDiagonals)
            {
                if( m_wrapperType == eWrapper )
                {
                    m_data = rhs.m_data;
                }
                else
                {
                    m_data = Array<OneD, DataType>(GetRequiredStorageSize());
                    CopyArrayN(rhs.m_data, m_data, m_data.num_elements());
                }
            }
            
            #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
                explicit NekMatrix(const NekMatrixMetadata& d) :
                    BaseType(d.Rows, d.Columns),
                    m_data(),
                    m_wrapperType(eCopy),
                    m_storagePolicy(eFULL),
                    m_numberOfSuperDiagonals(std::numeric_limits<unsigned int>::max()),
                    m_numberOfSubDiagonals(std::numeric_limits<unsigned int>::max())
                {
                    m_data = Array<OneD, DataType>(d.RequiredStorageSize);
                }

                template<typename ExpressionPolicyType>
                NekMatrix(const Expression<ExpressionPolicyType>& rhs) :
                    BaseType(rhs.GetMetadata().Rows, rhs.GetMetadata().Columns),
                    m_wrapperType(eCopy),
                    m_storagePolicy(eFULL),
                    m_numberOfSuperDiagonals(std::numeric_limits<unsigned int>::max()),
                    m_numberOfSubDiagonals(std::numeric_limits<unsigned int>::max())
                {
                    BOOST_MPL_ASSERT(( boost::is_same<typename Expression<ExpressionPolicyType>::ResultType, NekMatrix<const DataType, StandardMatrixTag> > ));
                    m_data = Array<OneD, DataType>(rhs.GetMetadata().RequiredStorageSize);
                    rhs.Evaluate(*this);
                }
            #endif //NEKTAR_USE_EXPRESSION_TEMPLATES

            
            MatrixStorage GetType() const { return m_storagePolicy; }
           
            ThisType& operator=(const ThisType& rhs)
            {
                if( this == &rhs )
                {
                    return *this;
                }

                BaseType::operator=(rhs);
                m_storagePolicy = rhs.m_storagePolicy;
                m_numberOfSubDiagonals = rhs.m_numberOfSubDiagonals;
                m_numberOfSuperDiagonals = rhs.m_numberOfSuperDiagonals;
                
                ResizeDataArrayIfNeeded();
                
                unsigned int requiredStorageSize = GetRequiredStorageSize();
                std::copy(rhs.m_data.data(), rhs.m_data.data() + requiredStorageSize, m_data.data());
                
                return *this;
            }
            
            #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
                template<typename ExpressionPolicyType>
                ThisType& operator=(const Expression<ExpressionPolicyType>& rhs)
                {
                    BOOST_MPL_ASSERT(( boost::is_same<typename Expression<ExpressionPolicyType>::ResultType, NekMatrix<const DataType, StandardMatrixTag> > ));
                    m_storagePolicy = eFULL;
                    m_numberOfSubDiagonals = std::numeric_limits<unsigned int>::max();
                    m_numberOfSuperDiagonals = std::numeric_limits<unsigned int>::max();
                    
                    if( this->GetRows() != rhs.GetMetadata().Rows ||
                        this->GetColumns() != rhs.GetMetadata().Columns )
                    {
                        Resize(rhs.GetMetadata().Rows, rhs.GetMetadata().Columns);
                        ResizeDataArrayIfNeeded(rhs.GetMetadata());
                    }

                    this->SetTransposeFlag('N');
                    rhs.Evaluate(*this);
                    return *this;
                }
            #endif //NEKTAR_USE_EXPRESSION_TEMPLATES
            
            /// \brief Returns the element value at the given row and column.
            typename boost::call_traits<DataType>::const_reference operator()(unsigned int row, unsigned int column) const
            {
                ASSERTL2(row < this->GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetRows()) +
                    std::string(" rows"));
                ASSERTL2(column < this->GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetColumns()) +
                    std::string(" columns"));
                
                return this->GetValue(row, column, this->GetTransposeFlag());
            }

            /// \brief Returns the element value at the given row and column.  
            /// \brief row The element's row.
            /// \brief column The element's column.
            /// \brief transpose If transpose = 'N', then the return value is element [row, column].
            ///                  If transpose = 'T', then the return value is element [column, row].
            typename boost::call_traits<DataType>::const_reference operator()(unsigned int row, unsigned int column, char transpose) const
            {    
                return this->GetValue(row, column, transpose);
            }
            
            unsigned int CalculateIndex(unsigned int row, unsigned int col, const char transpose) const
            {
                unsigned int numRows = this->GetSize()[0];
                unsigned int numColumns = this->GetSize()[1];
                if(transpose == 'T' )
                {
                    std::swap(row, col);
                }
                switch(m_storagePolicy)
                {
                    case eFULL:
                        return FullMatrixFuncs::CalculateIndex(numRows, numColumns, row, col);
                        break;
                    case eDIAGONAL:
                        return DiagonalMatrixFuncs::CalculateIndex(row, col);
                        break;
                    case eUPPER_TRIANGULAR:
                        return UpperTriangularMatrixFuncs::CalculateIndex(row, col);                        
                        break;
                    case eLOWER_TRIANGULAR:
                        return LowerTriangularMatrixFuncs::CalculateIndex(numColumns, row, col);
                        break;
                    case eSYMMETRIC:
                        return SymmetricMatrixFuncs::CalculateIndex(row, col);
                        break;
                    case eBANDED:
                        return BandedMatrixFuncs::CalculateIndex(numRows, numColumns,
                            row, col, m_numberOfSubDiagonals, m_numberOfSuperDiagonals);
                        break;
                    case eSYMMETRIC_BANDED:
                        NEKERROR(ErrorUtil::efatal, "Not yet implemented.");
                        break;
                    case eUPPER_TRIANGULAR_BANDED:
                        NEKERROR(ErrorUtil::efatal, "Not yet implemented.");
                        break;
                    case eLOWER_TRIANGULAR_BANDED:
                        NEKERROR(ErrorUtil::efatal, "Not yet implemented.");
                        break;
                        
                    default:
                        NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
                }
                
                return std::numeric_limits<unsigned int>::max();
            }
            
            /// \brief Returns the element value at the given row and column.
            typename boost::call_traits<DataType>::const_reference GetValue(unsigned int row, unsigned int column) const
            {
                ASSERTL2(row < this->GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetRows()) +
                    std::string(" rows"));
                ASSERTL2(column < this->GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetColumns()) +
                    std::string(" columns"));
                    
                return GetValue(row, column, this->GetTransposeFlag());
            }

            /// \brief Returns the element value at the given row and column.  
            /// \brief row The element's row.
            /// \brief column The element's column.
            /// \brief transpose If transpose = 'N', then the return value is element [row, column].
            ///                  If transpose = 'T', then the return value is element [column, row].
            typename boost::call_traits<DataType>::const_reference GetValue(unsigned int row, unsigned int column, char transpose) const
            {
                static DataType defaultReturnValue;
                unsigned int index = CalculateIndex(row, column, transpose);
                if( index != std::numeric_limits<unsigned int>::max() )
                {
                    return m_data[index];
                }
                else
                {
                    return defaultReturnValue;
                }
            }
            
            const Array<OneD, const DataType>& GetPtr() const
            {
                return m_data;
            }
            
            /// \brief Returns the scaling used by this matrix.
            ///
            /// Since this matrix is not a scaled matrix, this method always returns 1.
            DataType Scale() const
            {
                return DataType(1);
            }

            const DataType* GetRawPtr() const
            {
                return m_data.data();
            }
            
            typedef iterator_impl<const DataType, const ThisType> const_iterator;
            
            
            const_iterator begin() const
            {
                return begin(this->GetTransposeFlag());
            }

            const_iterator begin(char transpose) const
            { 
                if( transpose == 'N' )
                {
                    return const_iterator(m_data.data(), m_data.data() + m_data.num_elements());
                }
                else
                {
                    return const_iterator(this, transpose);
                }
            }

            const_iterator end() const
            {
                return end(this->GetTransposeFlag());
            }

            const_iterator end(char transpose) const
            {
                if( transpose == 'N' )
                {
                    return const_iterator(m_data.data(), m_data.data() + m_data.num_elements(), true);
                }
                else
                {
                    return const_iterator(this, transpose, true);
                }
            }
            
            
            unsigned int GetStorageSize() const
            {
                return m_data.num_elements();
            }
            
             // Banded matrices only.
            // Get the specified number of sub diagonals, or calculate the 
            // number if this object has not been initialized.
            unsigned int GetNumberOfSubDiagonals() const
            {
                if( m_numberOfSubDiagonals != std::numeric_limits<unsigned int>::max() )
                {
                    return m_numberOfSubDiagonals;
                }
                else if( this->GetRows() > 0 )
                {
                    return this->GetRows()-1;
                }
                else
                {
                    return 0;
                }
            }

            unsigned int GetNumberOfSuperDiagonals() const
            {
                if( m_numberOfSuperDiagonals != std::numeric_limits<unsigned int>::max() )
                {
                    return m_numberOfSuperDiagonals;
                }
                else if( this->GetRows() > 0 )
                {
                    return this->GetRows()-1;
                }
                else
                {
                    return 0;
                }
            }
            
            unsigned int CalculateNumberOfRows() const
            {
                return GetNumberOfSubDiagonals() + GetNumberOfSuperDiagonals() + 1;
            }
            /// \brief Returns true if the this matrix and rhs are equivalent.
            ///
            /// Two matrices are equivalent if they have the same size and each element
            /// is the same.
            bool operator==(const NekMatrix<DataType, StandardMatrixTag>& rhs) const
            {                
                if( this->GetRows() != rhs.GetRows() ||
                    this->GetColumns() != rhs.GetColumns() )
                {
                    return false;
                }

                if( this->GetTransposeFlag() == rhs.GetTransposeFlag() )
                {
                    return std::equal(begin(), end(), rhs.begin());
                }
                else
                {
                    for(unsigned int i = 0; i < this->GetRows(); ++i)
                    {
                        for(unsigned int j = 0; j < this->GetColumns(); ++j)
                        {
                            if( (*this)(i,j) != rhs(i,j) )
                            {
                                return false;
                            }
                        }
                    }
                }

                return true;
            }
                   
            PointerWrapper GetWrapperType() const { return m_wrapperType; }

            // The following methods purposefully hide the method with the same name in the base class.
            // We don't need to call the virtual function if we have an actual pointer to the 
            // derived class.
//            unsigned int GetRows() const
//            {
//                return BaseType::GetRowsForTranspose(this->GetRawTransposeFlag());
//            }
//            
//            unsigned int GetColumns() const
//            {
//                return BaseType::GetColumnsForTranspose(this->GetRawTransposeFlag());
//            }
            
            char GetTransposeFlag() const 
            {
                return this->GetRawTransposeFlag();
            }  
            
            
            boost::tuples::tuple<unsigned int, unsigned int> 
            Advance(unsigned int curRow, unsigned int curColumn) const
            {
                return Advance(curRow, curColumn, this->GetTransposeFlag());
            }
            
            boost::tuples::tuple<unsigned int, unsigned int> 
            Advance(unsigned int curRow, unsigned int curColumn, char transpose) const
            {
                unsigned int numRows = this->GetTransposedRows(transpose);
                unsigned int numColumns = this->GetTransposedColumns(transpose);
                
                switch(m_storagePolicy)
                {
                    case eFULL:
                        return FullMatrixFuncs::Advance(
                            numRows, numColumns, curRow, curColumn);
                        break;
                    case eDIAGONAL:
                        return DiagonalMatrixFuncs::Advance(
                            numRows, numColumns, curRow, curColumn);
                        break;
                    case eUPPER_TRIANGULAR:
                        return UpperTriangularMatrixFuncs::Advance(
                            numRows, numColumns, curRow, curColumn);
                        break;
                        
                    case eLOWER_TRIANGULAR:
                        return LowerTriangularMatrixFuncs::Advance(
                            numRows, numColumns, curRow, curColumn);
                        break;
                        
                    case eSYMMETRIC:
                        return SymmetricMatrixFuncs::Advance(
                            numRows, numColumns, curRow, curColumn);
                        break;
                    case eBANDED:
                        return BandedMatrixFuncs::Advance(
                            numRows, numColumns, curRow, curColumn);
                        break;
                    case eSYMMETRIC_BANDED:
                        NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
                        break;
                    case eUPPER_TRIANGULAR_BANDED:
                        NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
                        break;
                    case eLOWER_TRIANGULAR_BANDED:
                        NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
                        break;
                        
                    default:
                        NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
                }
                return boost::tuples::tuple<unsigned int, unsigned int>(curRow, curColumn);
            }
         
        protected:
            Array<OneD, DataType>& GetData() { return m_data; }
            void ResizeDataArrayIfNeeded(unsigned int requiredStorageSize)
            {
                if( m_wrapperType == eCopy  )
                {
                    // If the current vector is a matrix, then regardless of the rhs type 
                    // we just copy over the values, resizing if needed.
                    if( m_data.num_elements() < requiredStorageSize )
                    {
                        m_data = Array<OneD, DataType>(requiredStorageSize);
                    }
                }
                else if( m_wrapperType == eWrapper )
                {
                    // If the current matrix is wrapped, then just copy over the top,
                    // but the sizes of the two matrices must be the same.
                    ASSERTL0(m_data.num_elements() >= requiredStorageSize, "Wrapped NekMatrices must have the same dimension in operator=");
                }
            }
            
            void ResizeDataArrayIfNeeded()
            {
                unsigned int requiredStorageSize = GetRequiredStorageSize();
                ResizeDataArrayIfNeeded(requiredStorageSize);
            }
            
            void ResizeDataArrayIfNeeded(const NekMatrixMetadata& data)
            {
                ResizeDataArrayIfNeeded(data.RequiredStorageSize);
            }

        private:
            unsigned int GetRequiredStorageSize()
            {
                switch(m_storagePolicy)
                {
                    case eFULL:
                        return FullMatrixFuncs::GetRequiredStorageSize(this->GetRows(), this->GetColumns());
                        break;
                    case eDIAGONAL:
                        return DiagonalMatrixFuncs::GetRequiredStorageSize(this->GetRows(), this->GetColumns());
                        break;
                    case eUPPER_TRIANGULAR:
                        return UpperTriangularMatrixFuncs::GetRequiredStorageSize(this->GetRows(), this->GetColumns());
                        break;
                    case eLOWER_TRIANGULAR:
                        return LowerTriangularMatrixFuncs::GetRequiredStorageSize(this->GetRows(), this->GetColumns());
                        break;
                    case eSYMMETRIC:
                        return SymmetricMatrixFuncs::GetRequiredStorageSize(this->GetRows(), this->GetColumns());
                        break;
                    case eBANDED:
                        return BandedMatrixFuncs::GetRequiredStorageSize(this->GetRows(), this->GetColumns(),
                            m_numberOfSubDiagonals, m_numberOfSuperDiagonals);
                        break;
                    case eSYMMETRIC_BANDED:
                        NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
                        break;
                    case eUPPER_TRIANGULAR_BANDED:
                        NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
                        break;
                    case eLOWER_TRIANGULAR_BANDED:
                        NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
                        break;
                        
                    default:
                        NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
                }

                return 0;
            }
            
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
                return ThisType::GetType();
            }
            
            // We need to rethink class structure a little.  This shouldn't be necessary.
            virtual void v_SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d)
            {
            }
            
            Array<OneD, DataType> m_data;
            PointerWrapper m_wrapperType;
            MatrixStorage m_storagePolicy;
            
            // Only used by banded matrices.
            unsigned int m_numberOfSuperDiagonals;
            unsigned int m_numberOfSubDiagonals;
    };

     
    template<typename DataType>
    class NekMatrix<DataType, StandardMatrixTag> : public NekMatrix<const DataType, StandardMatrixTag>
    {
        public:
            typedef NekMatrix<const DataType, StandardMatrixTag> BaseType;
            typedef NekMatrix<DataType, StandardMatrixTag> ThisType;
            typedef DataType NumberType;

        public:
            NekMatrix() :
                BaseType()
            {
            }
            
            NekMatrix(unsigned int rows, unsigned int columns,
                      MatrixStorage policy = eFULL,
                      unsigned int subDiagonals = std::numeric_limits<unsigned int>::max(),
                      unsigned int superDiagonals = std::numeric_limits<unsigned int>::max()) :
                BaseType(rows, columns, policy, subDiagonals, superDiagonals)
            {
            }

            NekMatrix(unsigned int rows, unsigned int columns, typename boost::call_traits<DataType>::const_reference initValue,
                      MatrixStorage policy = eFULL,
                      unsigned int subDiagonals = std::numeric_limits<unsigned int>::max(),
                      unsigned int superDiagonals = std::numeric_limits<unsigned int>::max()) :
                BaseType(rows, columns, initValue, policy, subDiagonals, superDiagonals)
            {
            }
            
            NekMatrix(unsigned int rows, unsigned int columns, const DataType* data,
                      MatrixStorage policy = eFULL,
                      unsigned int subDiagonals = std::numeric_limits<unsigned int>::max(),
                      unsigned int superDiagonals = std::numeric_limits<unsigned int>::max()) :
                BaseType(rows, columns, data, policy, subDiagonals, superDiagonals)
            {
            }
            
            NekMatrix(unsigned int rows, unsigned int columns, const Array<OneD, const DataType>& d,
                      MatrixStorage policy = eFULL,
                      unsigned int subDiagonals = std::numeric_limits<unsigned int>::max(),
                      unsigned int superDiagonals = std::numeric_limits<unsigned int>::max()) :
                BaseType(rows, columns, d, policy, subDiagonals, superDiagonals)
            {
            }
            
            NekMatrix(unsigned int rows, unsigned int columns, const Array<OneD, DataType>& d, PointerWrapper wrapperType = eCopy,
                      MatrixStorage policy = eFULL,
                      unsigned int subDiagonals = std::numeric_limits<unsigned int>::max(),
                      unsigned int superDiagonals = std::numeric_limits<unsigned int>::max()) :
                BaseType(rows, columns, d, wrapperType, policy, subDiagonals, superDiagonals)
            {

            }
            
            NekMatrix(const ThisType& rhs) :
                BaseType(rhs)
            {
            }
            
            #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
                explicit NekMatrix(const NekMatrixMetadata& d) :
                    BaseType(d)
                {
                }

                template<typename ExpressionPolicyType>
                NekMatrix(const Expression<ExpressionPolicyType>& rhs) :
                    BaseType(rhs.GetMetadata())
                {
                    BOOST_MPL_ASSERT(( boost::is_same<typename Expression<ExpressionPolicyType>::ResultType, NekMatrix<DataType, StandardMatrixTag> > ));
                    rhs.Evaluate(*this);
                }
            #endif //NEKTAR_USE_EXPRESSION_TEMPLATES

            
            // TODO - Copy constructors from other types of matrices.
            
            ThisType& operator=(const ThisType& rhs)
            {
                if( this == &rhs )
                {
                    return *this;
                }

                BaseType::operator=(rhs);
                return *this;
            }
            
            #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
                template<typename ExpressionPolicyType>
                ThisType& operator=(const Expression<ExpressionPolicyType>& rhs)
                {
                    BOOST_MPL_ASSERT(( 
                        boost::mpl::or_
                        <
                            boost::is_same<typename Expression<ExpressionPolicyType>::ResultType, NekMatrix<DataType, StandardMatrixTag> >,
                            boost::is_same<typename Expression<ExpressionPolicyType>::ResultType, NekMatrix<const DataType, StandardMatrixTag> >
                        > ));
                    
                    if( this->GetRows() != rhs.GetMetadata().Rows ||
                        this->GetColumns() != rhs.GetMetadata().Columns )
                    {
                        Resize(rhs.GetMetadata().Rows, rhs.GetMetadata().Columns);
                        this->ResizeDataArrayIfNeeded(rhs.GetMetadata());
                    }

                    this->SetTransposeFlag('N');
                    rhs.Evaluate(*this);
                    return *this;
                }
            #endif //NEKTAR_USE_EXPRESSION_TEMPLATES

            class Proxy
            {
                public:
                    Proxy() : m_value(defaultReturnValue) {}                
                    explicit Proxy(DataType& value) : m_value(value) {}
                    Proxy(const Proxy& rhs) : m_value(rhs.m_value) {}
                    Proxy& operator=(const Proxy& rhs)
                    {
                        m_value = rhs.m_value;
                        return *this;
                    }
                    
                    DataType& operator*() { return m_value; }
                    const DataType& operator*() const { return m_value; }
                    operator DataType&() { return m_value; }
                    operator const DataType&() const { return m_value; }
                    void operator=(const DataType& newValue)
                    {
                        if( &m_value != &defaultReturnValue )
                        {
                            m_value = newValue;
                        }
                    }
                    
                private:
                    DataType& m_value;
                    static DataType defaultReturnValue;
            };

            
            using BaseType::operator();
            Proxy operator()(unsigned int row, unsigned int column)
            {
                ASSERTL2(row < this->GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetRows()) +
                    std::string(" rows"));
                ASSERTL2(column < this->GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetColumns()) +
                    std::string(" columns"));
                    
                return (*this)(row, column, this->GetTransposeFlag());
            }
            
                        
            Proxy operator()(unsigned int row, unsigned int column, char transpose)
            {   
                unsigned int index = this->CalculateIndex(row, column, transpose);
                if( index != std::numeric_limits<unsigned int>::max() )
                {
                    return Proxy(this->GetData()[index]);
                }
                else
                {
                    return Proxy();
                }    
//                static DataType defaultReturnValue;
//                unsigned int index = CalculateIndex(row, column, transpose);
//                if( index != std::numeric_limits<unsigned int>::max() )
//                {
//                    return GetData()[index];
//                }
//                else
//                {
//                    // Reset the default in case someone overwrites it.
//                    defaultReturnValue = DataType();
//                    return defaultReturnValue;
//                }
            }
            
            void SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d)
            {
                ASSERTL2(row < this->GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetRows()) +
                    std::string(" rows"));
                ASSERTL2(column < this->GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetColumns()) +
                    std::string(" columns"));
                SetValue(row, column, d, this->GetTransposeFlag());
            }

            void SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d, char transpose)
            {
                unsigned int index = this->CalculateIndex(row, column, transpose);
                if( index != std::numeric_limits<unsigned int>::max() )
                {
                    this->GetData()[index] = d;
                }
                else
                {
                    NEKERROR(ErrorUtil::efatal, "Can't assign values into zeroed elements of a special array.");
                }
            }

            using BaseType::GetPtr;
            Array<OneD, DataType>& GetPtr()
            {
                return this->GetData();
            }

            using BaseType::GetRawPtr;
            DataType* GetRawPtr()
            {
                return this->GetData().data();
            }

            typedef typename BaseType::template iterator_impl<DataType, ThisType> iterator;
            
            using BaseType::begin;
            using BaseType::end;
            iterator begin()
            {
                return begin(this->GetTransposeFlag());
            }

            iterator begin(char transpose) 
            { 
                if( transpose == 'N' )
                {
                    return iterator(this->GetData().data(), this->GetData().data() + this->GetData().num_elements());
                }
                else
                {
                    return iterator(this, transpose);
                }
            }

            iterator end()
            {
                return end(this->GetTransposeFlag());
            }

            iterator end(char transpose)
            {
                if( transpose == 'N' )
                {
                    return iterator(this->GetData().data(), this->GetData().data() + this->GetData().num_elements(), true);
                }
                else
                {
                    return iterator(this, transpose, true);
                }
            }

            void Invert()
            {
                ASSERTL0(this->GetRows()==this->GetColumns(), "Only square matrices can be inverted.");
                ASSERTL0(this->GetTransposeFlag()=='N', "Only untransposed matrices may be inverted.");

                switch(this->GetType())
                {
                    case eFULL:
                        FullMatrixFuncs::Invert(this->GetRows(), this->GetColumns(),
                            this->GetData(), this->GetTransposeFlag());
                        break;
                    case eDIAGONAL:
                        DiagonalMatrixFuncs::Invert(this->GetRows(), this->GetColumns(),
                            this->GetData());
                        break;
                    case eUPPER_TRIANGULAR:
                    case eLOWER_TRIANGULAR:
                    case eSYMMETRIC:
                    case eBANDED:
                        NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
                        break;
                    case eSYMMETRIC_BANDED:
                        NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
                        break;
                    case eUPPER_TRIANGULAR_BANDED:
                        NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
                        break;
                    case eLOWER_TRIANGULAR_BANDED:
                        NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
                        break;
                        
                    default:
                        NEKERROR(ErrorUtil::efatal, "Unhandled matrix type");
                } 
            }
                
            
        protected:
            
            
        private:
            virtual void v_SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d)
            {
                return ThisType::SetValue(row, column, d);
            }
    };
    
    template<typename DataType>
    DataType NekMatrix<DataType, StandardMatrixTag>::Proxy::defaultReturnValue;
            
    template<typename DataType>
    NekMatrix<DataType, StandardMatrixTag>
    Transpose(NekMatrix<DataType, StandardMatrixTag>& rhs)
    {
        NekMatrix<DataType, StandardMatrixTag> result(rhs.GetRows(), rhs.GetColumns(), 
            rhs.GetPtr(), eWrapper, rhs.GetType(), rhs.GetNumberOfSubDiagonals(), 
            rhs.GetNumberOfSuperDiagonals());
        result.Transpose();
        return result;
    }
    

}


#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_STANDARD_MATRIX_HPP
