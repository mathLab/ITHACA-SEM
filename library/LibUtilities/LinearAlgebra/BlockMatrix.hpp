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
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/LinearAlgebra/FullMatrixStoragePolicy.hpp>
#include <LibUtilities/LinearAlgebra/DiagonalMatrixStoragePolicy.hpp>
#include <LibUtilities/LinearAlgebra/TriangularMatrixStoragePolicy.hpp>
#include <LibUtilities/LinearAlgebra/BandedMatrixStoragePolicy.hpp>

#include <boost/shared_ptr.hpp>

namespace Nektar
{
    template<typename DataType, typename InnerStorageType, typename InnerMatrixType, typename StorageType>
    class NekMatrix<NekMatrix<DataType, InnerStorageType, InnerMatrixType>, StorageType, BlockMatrixTag> : public ConstMatrix<typename NekMatrix<DataType, InnerStorageType, InnerMatrixType>::NumberType>
    {
        public:
            typedef NekMatrix<DataType, InnerStorageType, InnerMatrixType> InnerType;
            typedef NekMatrix<InnerType, StorageType, BlockMatrixTag> ThisType;
            typedef typename InnerType::NumberType NumberType;
            typedef ConstMatrix<NumberType> BaseType;

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

                    // TODO
                    // This won't work if we want to specify the data as banded, 
                    // because the data type for banded requires consturctor paramters.
                    // It should work for everything else.  We may need to slightly
                    // rethink the template parameters.
                    typedef MatrixStoragePolicy<NumberType, StorageType> StoragePolicy;

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
                        if( m_curRow != std::numeric_limits<unsigned int>::max() )
                        {
                            boost::tie(m_curRow, m_curColumn) = StoragePolicy::Advance(
                                m_matrix.GetRows(), m_matrix.GetColumns(), m_curRow, m_curColumn, m_data);
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
                    //boost::shared_ptr<IteratorInnerType> m_curBlock;
                    unsigned int m_curRow;
                    unsigned int m_curColumn;
                    typename StoragePolicy::PolicySpecificDataHolderType m_data;
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
                m_data(numberOfBlockRows, numberOfBlockColumns, boost::shared_ptr<InnerType>()),
                m_rowSizes(numberOfBlockRows),
                m_columnSizes(numberOfBlockColumns),
                m_storageSize(0),
                m_numberOfBlockRows(numberOfBlockRows),
                m_numberOfBlockColumns(numberOfBlockColumns)
            {
                m_storageSize = this->GetRows()*this->GetColumns();
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
                m_data(numberOfBlockRows, numberOfBlockColumns, boost::shared_ptr<InnerType>()),
                m_rowSizes(numberOfBlockRows),
                m_columnSizes(numberOfBlockColumns),
                m_storageSize(0),
                m_numberOfBlockRows(numberOfBlockRows),
                m_numberOfBlockColumns(numberOfBlockColumns)
            {
                Initialize(rowsPerBlock, columnsPerBlock);
            }
            
            NekMatrix(unsigned int numberOfBlockRows, unsigned int numberOfBlockColumns,
                      const Array<OneD, const unsigned int>& rowsPerBlock, const Array<OneD, const unsigned int>& columnsPerBlock) :
                BaseType(std::accumulate(rowsPerBlock.data(), rowsPerBlock.data() + numberOfBlockRows, 0),
                         std::accumulate(columnsPerBlock.data(), columnsPerBlock.data() + numberOfBlockColumns, 0)),
                m_data(numberOfBlockRows, numberOfBlockColumns, boost::shared_ptr<InnerType>()),
                m_rowSizes(numberOfBlockRows),
                m_columnSizes(numberOfBlockColumns),
                m_storageSize(0),
                m_numberOfBlockRows(numberOfBlockRows),
                m_numberOfBlockColumns(numberOfBlockColumns)
            {
                Initialize(rowsPerBlock.data(), columnsPerBlock.data());
            }

            NekMatrix(const Array<OneD, const unsigned int>& rowsPerBlock,
                      const Array<OneD, const unsigned int>& columnsPerBlock) :
                BaseType(std::accumulate(rowsPerBlock.begin(), rowsPerBlock.end(), 0),
                         std::accumulate(columnsPerBlock.begin(), columnsPerBlock.end(), 0)),
                m_data(rowsPerBlock.num_elements(), columnsPerBlock.num_elements(), boost::shared_ptr<InnerType>()),
                m_rowSizes(rowsPerBlock.num_elements()),
                m_columnSizes(columnsPerBlock.num_elements()),
                m_storageSize(0),
                m_numberOfBlockRows(rowsPerBlock.num_elements()),
                m_numberOfBlockColumns(columnsPerBlock.num_elements())
            {
                Initialize(rowsPerBlock.data(), columnsPerBlock.data());
            }
                
            boost::shared_ptr<const InnerType> GetBlock(unsigned int row, unsigned int column) const        
            {
                ASSERTL2(row < m_numberOfBlockRows, std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockRows) +
                    std::string(" rows"));
                ASSERTL2(column < m_numberOfBlockColumns, std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockColumns) +
                    std::string(" columns"));
                return m_data[row][column];
            }
            
            boost::shared_ptr<InnerType> GetBlock(unsigned int row, unsigned int column)       
            {
                ASSERTL2(row < m_numberOfBlockRows, std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockRows) +
                    std::string(" rows"));
                ASSERTL2(column < m_numberOfBlockColumns, std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockColumns) +
                    std::string(" columns"));
                return m_data[row][column];
            }
            
            void SetBlock(unsigned int row, unsigned int column, boost::shared_ptr<InnerType>& m)
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
                ASSERTL2(col < this->GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(col) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetColumns()) +
                    std::string(" columns"));
                    
                unsigned int blockRow = std::lower_bound(m_rowSizes.begin(), m_rowSizes.end(), row) - m_rowSizes.begin();
                unsigned int blockColumn = std::lower_bound(m_columnSizes.begin(), m_columnSizes.end(), col) - m_columnSizes.begin();
                unsigned int actualRow = row;
                if( blockRow > 0 )
                {
                    actualRow = row-(m_rowSizes[blockRow-1])-1;
                }

                unsigned int actualCol = col;
                if( blockColumn > 0 )
                {
                    actualCol = col-(m_columnSizes[blockColumn-1])-1;
                }
                    
                const boost::shared_ptr<const InnerType> block = GetBlock(blockRow, blockColumn);
                if( block )
                {
                    return (*block)(actualRow, actualCol);
                }
                else
                {
                    return m_zeroElement;
                }
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
            
            unsigned int GetNumberOfRowsInBlockRow(unsigned int blockRow) const
            {
                ASSERTL2(blockRow < m_numberOfBlockRows, std::string("Block Row ") + boost::lexical_cast<std::string>(blockRow) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockRows) +
                    std::string(" block rows"));
                if( blockRow == 0 )
                {
                    return m_rowSizes[blockRow]+1;
                }
                else
                {
                    return m_rowSizes[blockRow] - m_rowSizes[blockRow-1];
                }
                
            }

            unsigned int GetNumberOfColumnsInBlockColumn(unsigned int blockCol) const
            {
                ASSERTL2(blockCol < m_numberOfBlockColumns, std::string("Block column ") + boost::lexical_cast<std::string>(blockCol) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockColumns) +
                    std::string(" block columns"));
                if( blockCol == 0 )
                {
                    return m_columnSizes[blockCol]+1;
                }
                else
                {
                    return m_columnSizes[blockCol] - m_columnSizes[blockCol-1];
                }
            }

            iterator begin() { return iterator(*this, 0, 0); }
            iterator end() { return iterator(*this); }
            const_iterator begin() const { return const_iterator(*this, 0, 0); }
            const_iterator end() const { return const_iterator(*this); }
                    
        public:
        
        private:
            void Initialize(const unsigned int* rowsPerBlock, const unsigned int* columnsPerBlock) 
            {
                m_storageSize = this->GetRows()*this->GetColumns();
                m_rowSizes[0] = rowsPerBlock[0] - 1;
                for(unsigned int i = 1; i < m_numberOfBlockRows; ++i)
                {
                    m_rowSizes[i] = rowsPerBlock[i] + m_rowSizes[i-1];
                }
                
                m_columnSizes[0] = columnsPerBlock[0] - 1;
                for(unsigned int i = 1; i < m_numberOfBlockColumns; ++i)
                {
                    m_columnSizes[i] = columnsPerBlock[i] + m_columnSizes[i-1];
                }
            }

            virtual typename boost::call_traits<NumberType>::value_type v_GetValue(unsigned int row, unsigned int column) const 
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
            
            Array<TwoD, boost::shared_ptr<InnerType> > m_data;
            Array<OneD, unsigned int> m_rowSizes;
            Array<OneD, unsigned int> m_columnSizes;
            unsigned int m_storageSize;
            unsigned int m_numberOfBlockRows;
            unsigned int m_numberOfBlockColumns; 
            static NumberType m_zeroElement;
    };

    template<typename DataType, typename InnerStorageType, typename InnerMatrixType, typename StorageType>
    typename NekMatrix<NekMatrix<DataType, InnerStorageType, InnerMatrixType>, StorageType, BlockMatrixTag>::NumberType
    NekMatrix<NekMatrix<DataType, InnerStorageType, InnerMatrixType>, StorageType, BlockMatrixTag>::m_zeroElement(0);
}


#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BLOCK_MATRIX_HPP
