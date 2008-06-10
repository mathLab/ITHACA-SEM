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
#include <LibUtilities/LinearAlgebra/FullMatrixStoragePolicy.hpp>
#include <LibUtilities/LinearAlgebra/DiagonalMatrixStoragePolicy.hpp>
#include <LibUtilities/LinearAlgebra/TriangularMatrixStoragePolicy.hpp>
#include <LibUtilities/LinearAlgebra/PointerWrapper.h>
#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrixMetadata.hpp>

namespace Nektar
{
    /// \brief Standard Matrix
    /// \param DataType The type stored in each element.
    /// \param StorageType The matrix storage type.  Valid values are DiagonalMatrixTag,
    ///                    FullMatrixTag, BandedMatrixTag, UpperTriangularMatrixTag,
    ///                    LowerTriangularMatrixTag, and SymmetricMatrixTag.
    ///
    /// Matrices are stored in column major order to make it easier to interoperate with 
    /// Blas and Lapack.
    template<typename DataType, typename StorageType>
    class NekMatrix<const DataType, StorageType, StandardMatrixTag> : public Matrix<DataType>
    {
        public:
            typedef Matrix<DataType> BaseType;
            typedef NekMatrix<const DataType, StorageType, StandardMatrixTag> ThisType;
            typedef MatrixStoragePolicy<DataType, StorageType> StoragePolicy;
            typedef const DataType NumberType;
            typedef typename StoragePolicy::PolicySpecificDataHolderType PolicySpecificDataHolderType;

            typedef typename StoragePolicy::GetValueReturnType GetValueType;
            typedef typename boost::call_traits<DataType>::const_reference ConstGetValueType;
            
            public:
                template<typename T, typename MatrixType>
                class iterator_impl
                {
                    public:
                        typedef typename StoragePolicy::template reference<T>::type value_type;
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
                            m_transpose(transpose)
                        {
                            if( isEnd )
                            {
                                m_curRow = std::numeric_limits<unsigned int>::max();
                                m_curColumn = std::numeric_limits<unsigned int>::max();
                            }
                        }

                        iterator_impl(const iterator_impl<T, MatrixType>& rhs) :
                            m_data(rhs.m_data),
                            m_end(rhs.m_end),
                            m_curRow(rhs.m_curRow),
                            m_curColumn(rhs.m_curColumn),
                            m_matrix(rhs.m_matrix),
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
                                return (*m_matrix)(m_curRow, m_curColumn, m_transpose);
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
                                return (*m_matrix)(m_curRow, m_curColumn, m_transpose);
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
                                    StoragePolicy::Advance(m_matrix->GetRowsForTranspose(m_transpose), m_matrix->GetColumnsForTranspose(m_transpose),
                                        m_curRow, m_curColumn, 
                                        m_matrix->GetPolicySpecificDataHolderType());
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
                        char m_transpose;
            };

        public:
            /// \brief Creates an empty matrix.
            NekMatrix() :
                BaseType(0, 0),
                m_data(StoragePolicy::Initialize()),
                m_wrapperType(eCopy),
                m_policySpecificData()
            {
            }
            
            /// \brief Creates a rows by columns matrix.
            /// \brief rows The number of rows in the matrix.
            /// \brief columns The number of columns in the matrix.
            /// \brief policySpecificData Data the is specific to the storage type assigned to this matrix.
            ///                           Check MatrixStoragePolicy<StorageType> to see if the StoragePolicy
            ///                           for this matrix has policy specific data.
            NekMatrix(unsigned int rows, unsigned int columns, const PolicySpecificDataHolderType& policySpecificData = PolicySpecificDataHolderType()) :
                BaseType(rows, columns),
                m_data(StoragePolicy::Initialize(rows, columns, policySpecificData)),
                m_wrapperType(eCopy),
                m_policySpecificData(policySpecificData)
            {
            }

            /// \brief Creates a rows by columns matrix.
            /// \brief rows The number of rows in the matrix.
            /// \brief columns The number of columns in the matrix.
            /// \brief initValue The value used to initialize each element.
            /// \brief policySpecificData Data the is specific to the storage type assigned to this matrix.
            ///                           Check MatrixStoragePolicy<StorageType> to see if the StoragePolicy
            ///                           for this matrix has policy specific data.
            NekMatrix(unsigned int rows, unsigned int columns, typename boost::call_traits<DataType>::const_reference initValue,
                      const PolicySpecificDataHolderType& policySpecificData = PolicySpecificDataHolderType()) :
                BaseType(rows, columns),
                m_data(StoragePolicy::Initialize(rows, columns, initValue, policySpecificData)),
                m_wrapperType(eCopy),
                m_policySpecificData(policySpecificData)
            {
            }
            
            /// \brief Creates a rows by columns matrix.
            /// \brief rows The number of rows in the matrix.
            /// \brief columns The number of columns in the matrix.
            /// \brief data An array of data use to initialize the matrix.
            /// \brief policySpecificData Data the is specific to the storage type assigned to this matrix.
            ///                           Check MatrixStoragePolicy<StorageType> to see if the StoragePolicy
            ///                           for this matrix has policy specific data.
            NekMatrix(unsigned int rows, unsigned int columns, const DataType* data,
                      const PolicySpecificDataHolderType& policySpecificData = PolicySpecificDataHolderType()) :
                BaseType(rows, columns),
                m_data(StoragePolicy::Initialize(rows, columns, data, policySpecificData)),
                m_wrapperType(eCopy),
                m_policySpecificData(policySpecificData)
            {
            }
            
            /// \brief Creates a rows by columns matrix.
            /// \brief rows The number of rows in the matrix.
            /// \brief columns The number of columns in the matrix.
            /// \brief d An array of data used to initialize the matrix.  Values from d are copied into the matrix.
            /// \brief policySpecificData Data the is specific to the storage type assigned to this matrix.
            ///                           Check MatrixStoragePolicy<StorageType> to see if the StoragePolicy
            ///                           for this matrix has policy specific data.
            NekMatrix(unsigned int rows, unsigned int columns, const Array<OneD, const DataType>& d,
                      const PolicySpecificDataHolderType& policySpecificData = PolicySpecificDataHolderType()) :
                BaseType(rows, columns),
                m_data(StoragePolicy::Initialize(rows, columns, d, policySpecificData)),
                m_wrapperType(eCopy),
                m_policySpecificData(policySpecificData)
            {
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
                      const PolicySpecificDataHolderType& policySpecificData = PolicySpecificDataHolderType()) :
                BaseType(rows, columns),
                m_data(),
                m_wrapperType(wrapperType),
                m_policySpecificData(policySpecificData)
            {
                if( wrapperType == eWrapper )
                {
                     m_data = Array<OneD, DataType>(d, eVECTOR_WRAPPER);
                }
                else
                {
                    m_data = StoragePolicy::Initialize(rows, columns, d, policySpecificData);
                }
            }
            
            NekMatrix(const ThisType& rhs) :
                BaseType(rhs),
                m_data(),
                m_wrapperType(rhs.m_wrapperType),
                m_policySpecificData(rhs.m_policySpecificData)
            {
                if( m_wrapperType == eWrapper )
                {
                    m_data = rhs.m_data;
                }
                else
                {
                    m_data = StoragePolicy::Initialize(this->GetRows(), this->GetColumns(), m_policySpecificData);
                    CopyArray(rhs.m_data, m_data);
                }
            }
            
            #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
                explicit NekMatrix(const NekMatrixMetadata& d) :
                    BaseType(d.Rows, d.Columns),
                    m_data(),
                    m_wrapperType(eCopy),
                    m_policySpecificData()
                {
                    m_data = StoragePolicy::Initialize(this->GetRows(), this->GetColumns(), m_policySpecificData);
                }

                template<typename ExpressionPolicyType>
                NekMatrix(const Expression<ExpressionPolicyType>& rhs) :
                    BaseType(rhs.GetMetadata().Rows, rhs.GetMetadata().Columns),
                    m_data(),
                    m_wrapperType(eCopy),
                    m_policySpecificData()
                {
                    BOOST_MPL_ASSERT(( boost::is_same<typename Expression<ExpressionPolicyType>::ResultType, NekMatrix<const DataType, StorageType, StandardMatrixTag> > ));
                    m_data = StoragePolicy::Initialize(this->GetRows(), this->GetColumns(), m_policySpecificData);
                    rhs.Evaluate(*this);
                }
            #endif //NEKTAR_USE_EXPRESSION_TEMPLATES

            MatrixStorage GetStorageType() const 
            {
                return static_cast<MatrixStorage>(ConvertToMatrixStorageEnum<StorageType>::Value);
            }
            
            // TODO - Copy constructors from other types of matrices.
            
            ThisType& operator=(const ThisType& rhs)
            {
                if( this == &rhs )
                {
                    return *this;
                }

                BaseType::operator=(rhs);
                m_policySpecificData = rhs.m_policySpecificData;
                
                ResizeDataArrayIfNeeded(this->GetRows(), this->GetColumns(), m_policySpecificData);
                
                unsigned int requiredStorageSize = StoragePolicy::GetRequiredStorageSize(this->GetRows(), this->GetColumns(), m_policySpecificData);
                std::copy(rhs.m_data.data(), rhs.m_data.data() + requiredStorageSize, m_data.data());
                
                return *this;
            }
            
            #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
                template<typename ExpressionPolicyType>
                ThisType& operator=(const Expression<ExpressionPolicyType>& rhs)
                {
                    BOOST_MPL_ASSERT(( boost::is_same<typename Expression<ExpressionPolicyType>::ResultType, NekMatrix<const DataType, StorageType, StandardMatrixTag> > ));
                    m_policySpecificData = PolicySpecificDataHolderType();
                    if( this->GetRows() != rhs.GetMetadata().Rows ||
                        this->GetColumns() != rhs.GetMetadata().Columns )
                    {
                        Resize(rhs.GetMetadata().Rows, rhs.GetMetadata().Columns);
                        ResizeDataArrayIfNeeded(this->GetRows(), this->GetColumns(), m_policySpecificData);
                    }

                    this->SetTransposeFlag('N');
                    rhs.Evaluate(*this);
                    return *this;
                }
            #endif //NEKTAR_USE_EXPRESSION_TEMPLATES
            
            /// \brief Returns the element value at the given row and column.
            ConstGetValueType operator()(unsigned int row, unsigned int column) const
            {
                ASSERTL2(row < this->GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetRows()) +
                    std::string(" rows"));
                ASSERTL2(column < this->GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetColumns()) +
                    std::string(" columns"));
                    
                return (*this)(row, column, this->GetTransposeFlag());
            }

            /// \brief Returns the element value at the given row and column.  
            /// \brief row The element's row.
            /// \brief column The element's column.
            /// \brief transpose If transpose = 'N', then the return value is element [row, column].
            ///                  If transpose = 'T', then the return value is element [column, row].
            ConstGetValueType operator()(unsigned int row, unsigned int column, char transpose) const
            {       
                return StoragePolicy::GetValue(this->GetRowsForTranspose(transpose), this->GetColumnsForTranspose(transpose), row, column, m_data, transpose, m_policySpecificData);    
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
                return StoragePolicy::GetValue(this->GetRowsForTranspose(transpose), this->GetColumnsForTranspose(transpose), row, column, m_data, transpose, m_policySpecificData);
            }

            
            const PolicySpecificDataHolderType& GetPolicySpecificDataHolderType() const 
            {
                return m_policySpecificData;
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
            
            //typedef DataType* iterator;
            //typedef const DataType* const_iterator;
            //iterator begin() { return m_data.data(); }
            //iterator end() { return m_data.data() + m_data.num_elements(); }
            
            //const_iterator begin() const { return m_data.data(); }
            //const_iterator end() const { return m_data.data() + m_data.num_elements(); }

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
            
            /// \brief Returns true if the this matrix and rhs are equivalent.
            ///
            /// Two matrices are equivalent if they have the same size and each element
            /// is the same.
            bool operator==(const NekMatrix<DataType, StorageType, StandardMatrixTag>& rhs) const
            {
                if( GetStorageSize() != rhs.GetStorageSize() )
                {
                    return false;
                }
                
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
            unsigned int GetRows() const
            {
                return BaseType::GetRowsForTranspose(this->GetRawTransposeFlag());
            }
            
            unsigned int GetColumns() const
            {
                return BaseType::GetColumnsForTranspose(this->GetRawTransposeFlag());
            }
            
            char GetTransposeFlag() const 
            {
                return this->GetRawTransposeFlag();
            }  
            
            
        protected:
            Array<OneD, DataType>& GetData() { return m_data; }
            
            PolicySpecificDataHolderType& GetPolicySpecificDataHolderTypeByReference()  
            {
                return m_policySpecificData;
            }
            
            void ResizeDataArrayIfNeeded(unsigned int newRows, unsigned int newColumns, const PolicySpecificDataHolderType& policySpecificData)
            {
                unsigned int requiredStorageSize = StoragePolicy::GetRequiredStorageSize(newRows, newColumns, policySpecificData);

                if( m_wrapperType == eCopy  )
                {
                    // If the current vector is a matrix, then regardless of the rhs type 
                    // we just copy over the values, resizing if needed.
                    if( m_data.num_elements() < requiredStorageSize )
                    {
                        m_data = StoragePolicy::Initialize(newRows, newColumns, policySpecificData);
                    }
                }
                else if( m_wrapperType == eWrapper )
                {
                    // If the current matrix is wrapped, then just copy over the top,
                    // but the sizes of the two matrices must be the same.
                    ASSERTL0(m_data.num_elements() >= requiredStorageSize, "Wrapped NekMatrices must have the same dimension in operator=");
                }
            }

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
            
            // We need to rethink class structure a little.  This shouldn't be necessary.
            virtual void v_SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d)
            {
            }
            
            Array<OneD, DataType> m_data;
            PointerWrapper m_wrapperType;
            PolicySpecificDataHolderType m_policySpecificData;
    };

     
    template<typename DataType, typename StorageType>
    class NekMatrix<DataType, StorageType, StandardMatrixTag> : public NekMatrix<const DataType, StorageType, StandardMatrixTag>
    {
        public:
            typedef NekMatrix<const DataType, StorageType, StandardMatrixTag> BaseType;
            typedef NekMatrix<DataType, StorageType, StandardMatrixTag> ThisType;
            typedef typename BaseType::StoragePolicy StoragePolicy;
            typedef DataType NumberType;
            typedef typename StoragePolicy::PolicySpecificDataHolderType PolicySpecificDataHolderType;

            typedef typename StoragePolicy::GetValueReturnType GetValueType;
            typedef typename boost::call_traits<DataType>::const_reference ConstGetValueType;
            
        public:
            NekMatrix() :
                BaseType()
            {
            }
            
            NekMatrix(unsigned int rows, unsigned int columns, const PolicySpecificDataHolderType& policySpecificData = PolicySpecificDataHolderType()) :
                BaseType(rows, columns, policySpecificData)
            {
            }

            NekMatrix(unsigned int rows, unsigned int columns, typename boost::call_traits<DataType>::const_reference initValue,
                      const PolicySpecificDataHolderType& policySpecificData = PolicySpecificDataHolderType()) :
                BaseType(rows, columns, initValue, policySpecificData)
            {
            }
            
            NekMatrix(unsigned int rows, unsigned int columns, const DataType* data,
                      const PolicySpecificDataHolderType& policySpecificData = PolicySpecificDataHolderType()) :
                BaseType(rows, columns, data, policySpecificData)
            {
            }
            
            NekMatrix(unsigned int rows, unsigned int columns, const Array<OneD, const DataType>& d,
                      const PolicySpecificDataHolderType& policySpecificData = PolicySpecificDataHolderType()) :
                BaseType(rows, columns, d, policySpecificData)
            {
            }
            
            NekMatrix(unsigned int rows, unsigned int columns, const Array<OneD, DataType>& d, PointerWrapper wrapperType = eCopy,
                      const PolicySpecificDataHolderType& policySpecificData = PolicySpecificDataHolderType()) :
                BaseType(rows, columns, d, wrapperType, policySpecificData)
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
                    BOOST_MPL_ASSERT(( boost::is_same<typename Expression<ExpressionPolicyType>::ResultType, NekMatrix<DataType, StorageType, StandardMatrixTag> > ));
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
                            boost::is_same<typename Expression<ExpressionPolicyType>::ResultType, NekMatrix<DataType, StorageType, StandardMatrixTag> >,
                            boost::is_same<typename Expression<ExpressionPolicyType>::ResultType, NekMatrix<const DataType, StorageType, StandardMatrixTag> >
                        > ));
                    
                    this->GetPolicySpecificDataHolderTypeByReference() = PolicySpecificDataHolderType();
                    if( this->GetRows() != rhs.GetMetadata().Rows ||
                        this->GetColumns() != rhs.GetMetadata().Columns )
                    {
                        Resize(rhs.GetMetadata().Rows, rhs.GetMetadata().Columns);
                        this->ResizeDataArrayIfNeeded(this->GetRows(), this->GetColumns(), this->GetPolicySpecificDataHolderType());
                    }

                    this->SetTransposeFlag('N');
                    rhs.Evaluate(*this);
                    return *this;
                }
            #endif //NEKTAR_USE_EXPRESSION_TEMPLATES
            
            
            using BaseType::operator();
            GetValueType operator()(unsigned int row, unsigned int column)
            {
                ASSERTL2(row < this->GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetRows()) +
                    std::string(" rows"));
                ASSERTL2(column < this->GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetColumns()) +
                    std::string(" columns"));
                    
                return (*this)(row, column, this->GetTransposeFlag());
            }
            
            GetValueType operator()(unsigned int row, unsigned int column, char transpose)
            {       
                return StoragePolicy::GetValue(this->GetRowsForTranspose(transpose), this->GetColumnsForTranspose(transpose), row, column, this->GetData(), transpose, this->GetPolicySpecificDataHolderType());    
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
                StoragePolicy::SetValue(this->GetRowsForTranspose(transpose), this->GetColumnsForTranspose(transpose), row, column, this->GetData(), d, transpose, this->GetPolicySpecificDataHolderType());
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

            
            //typedef DataType* iterator;
            //typedef const DataType* const_iterator;
            //iterator begin() { return m_data.data(); }
            //iterator end() { return m_data.data() + m_data.num_elements(); }
            
            //const_iterator begin() const { return m_data.data(); }
            //const_iterator end() const { return m_data.data() + m_data.num_elements(); }


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
                StoragePolicy::Invert(this->GetRows(), this->GetColumns(), this->GetData(), this->GetTransposeFlag(), this->GetPolicySpecificDataHolderType());
            }
                
            
        protected:
            
            
        private:
            virtual void v_SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d)
            {
                return ThisType::SetValue(row, column, d);
            }
    };

    template<typename DataType, typename StorageType>
    NekMatrix<DataType, StorageType, StandardMatrixTag>
    Transpose(NekMatrix<DataType, StorageType, StandardMatrixTag>& rhs)
    {
        NekMatrix<DataType, StorageType, StandardMatrixTag> result(rhs.GetRows(), rhs.GetColumns(), 
            rhs.GetPtr(), eWrapper, rhs.GetPolicySpecificDataHolderType());
        result.Transpose();
        return result;
    }
    

}


#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_STANDARD_MATRIX_HPP
