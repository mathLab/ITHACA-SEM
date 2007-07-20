///////////////////////////////////////////////////////////////////////////////
//
// File: BlockMatrix.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BLOCK_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BLOCK_MATRIX_HPP

#include <LibUtilities/LinearAlgebra/MatrixBase.hpp>
#include <LibUtilities/BasicUtils/BinaryExpressionTraits.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LinearAlgebra/MatrixTraits.hpp>
#include <LibUtilities/BasicUtils/OperatorGenerators.hpp>
#include <LibUtilities/LinearAlgebra/MatrixTraits.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>

#include <LibUtilities/BasicUtils/SharedPtr.hpp>

namespace Nektar
{
    template<typename DataType, typename StorageType, typename InnerMatrixType>
    class NekMatrix<NekMatrix<DataType, StorageType, InnerMatrixType>, StorageType, BlockMatrixTag> : public Matrix<DataType>
    {
        public:
            typedef Matrix<DataType> BaseType;
            typedef NekMatrix<DataType, StorageType, InnerMatrixType> InnerType;
            typedef NekMatrix<InnerType, StorageType, BlockMatrixTag> ThisType;
            typedef typename NekMatrix<DataType, StorageType, InnerMatrixType>::NumberType NumberType;

            // Each inner matrix type can possible return references or value types from GetValue.
            // Query the type here to find out.
            typedef typename InnerType::GetValueType GetValueType;
            typedef typename InnerType::ConstGetValueType ConstGetValueType;
                
        public:
            template<typename MatrixType>
            class iterator_base
            {
                public:
                    typedef typename MatrixType::InnerType IteratorInnerType;
                    
                public:                   
                    iterator_base(MatrixType& m, unsigned int curRow, unsigned int curCol) :
                        m_matrix(m),
                        m_curRow(curRow),
                        m_curColumn(curCol)
                    {
                    }
                    
                    iterator_base(MatrixType& m) :
                        m_matrix(m),
                        m_curRow(std::numeric_limits<unsigned int>::max()),
                        m_curColumn(std::numeric_limits<unsigned int>::max())
                    {
                    }
                    
                    iterator_base(const iterator_base<MatrixType>& rhs) :
                        m_matrix(rhs.m_matrix),
                        m_curRow(rhs.m_curRow),
                        m_curColumn(rhs.m_curColumn)
                    {
                    }
                    
                    void operator++()
                    {
                        if( m_curColumn != std::numeric_limits<unsigned int>::max() )
                        {
                            ++m_curColumn;
                            
                            if( m_curColumn >= m_matrix.GetColumns() )
                            {
                                m_curColumn = 0;
                                ++m_curRow;
                                
                                if( m_curRow >= m_matrix.GetRows() )
                                {
                                    m_curColumn = std::numeric_limits<unsigned int>::max();
                                    m_curRow = std::numeric_limits<unsigned int>::max();
                                }
                            }
                            
//                             if( m_curColumn != std::numeric_limits<unsigned int>::max() )
//                             {
//                                 m_curBlock = m_matrix.GetBlock(m_curRow, m_curColumn);
//                             }
                        }
                    }
                    
                    NumberType operator*()
                    {
                        return m_matrix(m_curRow, m_curColumn);
                    }
                                        
                    bool operator==(const iterator_base<MatrixType>& rhs)
                    {
                        return m_curRow == rhs.m_curRow && m_curColumn == rhs.m_curColumn;
                    }
                    
                    bool operator!=(const iterator_base<MatrixType>& rhs)
                    {
                        return !(*this == rhs);
                    }
                    
                private:
                    iterator_base<MatrixType>& operator=(const iterator_base<MatrixType>& rhs);
                                
                    MatrixType& m_matrix;
                    //ptr<IteratorInnerType> m_curBlock;
                    unsigned int m_curRow;
                    unsigned int m_curColumn;
            };
            
            typedef iterator_base<ThisType> iterator;
            typedef iterator_base<const ThisType> const_iterator;
                                                            
        public:
            NekMatrix() :
                BaseType(0,0),
                m_data(),
                m_rowSizes(),
                m_columnSizes(),
                m_storageSize(),
                m_numberOfBlockRows(0),
                m_numberOfBlockColumns(0)
            {
            }
            
            NekMatrix(unsigned int numberOfBlockRows, unsigned int numberOfBlockColumns,
                      unsigned int rowsPerBlock, unsigned int columnsPerBlock) :
                BaseType(numberOfBlockRows*rowsPerBlock, numberOfBlockColumns*columnsPerBlock),
                m_data(numberOfBlockRows, numberOfBlockColumns, ptr<InnerType>()),
                m_rowSizes(numberOfBlockRows),
                m_columnSizes(numberOfBlockColumns),
                m_storageSize(this->GetRows()*this->GetColumns()),
                m_numberOfBlockRows(numberOfBlockRows),
                m_numberOfBlockColumns(numberOfBlockColumns)
            {
                for(unsigned int i = 1; i <= numberOfBlockRows; ++i)
                {
                    m_rowSizes[i-1] = i*rowsPerBlock-1;
                }
                
                for(unsigned int i = 1; i <= numberOfBlockColumns; ++i)
                {
                    m_columnSizes[i-1] = i*columnsPerBlock-1;
                }
            }
            
            NekMatrix(unsigned int numberOfBlockRows, unsigned int numberOfBlockColumns,
                      unsigned int* rowsPerBlock, unsigned int* columnsPerBlock) :
                BaseType(std::accumulate(rowsPerBlock, rowsPerBlock + numberOfBlockRows, 0),
                         std::accumulate(columnsPerBlock, columnsPerBlock + numberOfBlockColumns, 0)),
                m_data(numberOfBlockRows, numberOfBlockColumns, ptr<InnerType>()),
                m_rowSizes(numberOfBlockRows),
                m_columnSizes(numberOfBlockColumns),
                m_storageSize(this->GetRows()*this->GetColumns()),
                m_numberOfBlockRows(numberOfBlockRows),
                m_numberOfBlockColumns(numberOfBlockColumns)
            {
                m_rowSizes[0] = rowsPerBlock[0] - 1;
                for(unsigned int i = 1; i < numberOfBlockRows; ++i)
                {
                    m_rowSizes[i] = rowsPerBlock[i] + m_rowSizes[i-1];
                }
                
                m_columnSizes[0] = columnsPerBlock[0] - 1;
                for(unsigned int i = 1; i < numberOfBlockColumns; ++i)
                {
                    m_columnSizes[i] = columnsPerBlock[i] + m_columnSizes[i-1];
                }
            }
                
            ptr<const InnerType> GetBlock(unsigned int row, unsigned int column) const        
            {
                ASSERTL2(row < m_numberOfBlockRows, std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockRows) +
                    std::string(" rows"));
                ASSERTL2(column < m_numberOfBlockColumns, std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockColumns) +
                    std::string(" columns"));
                return m_data[row][column];
            }
            
            ptr<InnerType> GetBlock(unsigned int row, unsigned int column)       
            {
                ASSERTL2(row < m_numberOfBlockRows, std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockRows) +
                    std::string(" rows"));
                ASSERTL2(column < m_numberOfBlockColumns, std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockColumns) +
                    std::string(" columns"));
                return m_data[row][column];
            }
            
            void SetBlock(unsigned int row, unsigned int column, ptr<InnerType>& m)
            {
                ASSERTL2(row < m_numberOfBlockRows, std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockRows) +
                    std::string(" rows"));
                ASSERTL2(column < m_numberOfBlockColumns, std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockColumns) +
                    std::string(" columns"));
                m_data[row][column] = m;
            }
            
            
            
            
            ConstGetValueType operator()(unsigned int row, unsigned int col) const
            {
                ASSERTL2(row < this->GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetRows()) +
                    std::string(" rows"));
                ASSERTL2(column < this->GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetColumns()) +
                    std::string(" columns"));
                    
                unsigned int blockRow = std::lower_bound(m_rowSizes.begin(), m_rowSizes.end(), row) - m_rowSizes.begin();
                unsigned int blockColumn = std::lower_bound(m_columnSizes.begin(), m_columnSizes.end(), col) - m_columnSizes.begin();
                unsigned int actualRow = row-(m_rowSizes[blockRow]-1);
                unsigned int actualCol = col-(m_columnSizes[blockColumn]-1);
                
                ASSERTL2(GetBlock(blockRow, blockColumn), std::string("Attempting to access block (") +
                    boost::lexical_cast<std::string>(blockRow) + std::string(", ") + 
                    boost::lexical_cast<std::string>(blockColumn) + std::string(") of a block matrix but it is null."));
                    
                return GetBlock(blockRow, blockColumn)->operator()(actualRow, actualCol);
            }
            
            GetValueType operator()(unsigned int row, unsigned int col)
            {
                ASSERTL2(row < this->GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetRows()) +
                    std::string(" rows"));
                ASSERTL2(column < this->GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetColumns()) +
                    std::string(" columns"));
                    
                unsigned int blockRow = std::lower_bound(m_rowSizes.begin(), m_rowSizes.end(), row) - m_rowSizes.begin();
                unsigned int blockColumn = std::lower_bound(m_columnSizes.begin(), m_columnSizes.end(), col) - m_columnSizes.begin();
                unsigned int actualRow = row-(m_rowSizes[blockRow]-1);
                unsigned int actualCol = col-(m_columnSizes[blockColumn]-1);
                
                ASSERTL2(GetBlock(blockRow, blockColumn), std::string("Attempting to access block (") +
                    boost::lexical_cast<std::string>(blockRow) + std::string(", ") + 
                    boost::lexical_cast<std::string>(blockColumn) + std::string(") of a block matrix but it is null."));
                    
                return GetBlock(blockRow, blockColumn)->operator()(actualRow, actualCol);
            }
            
            unsigned int GetStorageSize() const 
            {
                return m_storageSize;
            }
             
            MatrixStorage GetStorageType() const
            {
                return static_cast<MatrixStorage>(ConvertToMatrixStorageEnum<StorageType>::Value);
            }            
            
            unsigned int GetNumberOfBlockRows() const { return m_numberOfBlockRows; }
            unsigned int GetNumberOfBlockColumns() const { return m_numberOfBlockColumns; }
            
            iterator begin() { return iterator(*this, 0, 0); }
            iterator end() { return iterator(*this); }
            const_iterator begin() const { return const_iterator(*this, 0, 0); }
            const_iterator end() const { return const_iterator(*this); }
                    
        public:
        
        private:
            virtual typename boost::call_traits<DataType>::value_type v_GetValue(unsigned int row, unsigned int column) const 
            {
                return (*this)(row, column);
            }
            
            virtual unsigned int v_GetStorageSize() const 
            {
                return this->GetStorageSize();
            }
            
            virtual MatrixStorage v_GetStorageType() const
            {
                return this->GetStorageType();
            }
            
            virtual void v_SetValue(unsigned int row, unsigned int column, typename boost::call_traits<DataType>::const_reference d)
            {
                (*this)(row, column) = d;
            }
            
            Array<TwoD, ptr<InnerType> > m_data;
            Array<OneD, unsigned int> m_rowSizes;
            Array<OneD, unsigned int> m_columnSizes;
            unsigned int m_storageSize;
            unsigned int m_numberOfBlockRows;
            unsigned int m_numberOfBlockColumns;            
    };
    
}


#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BLOCK_MATRIX_HPP
