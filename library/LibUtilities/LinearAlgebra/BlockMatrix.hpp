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
#include <LibUtilities/LinearAlgebra/MatrixFuncs.h>

#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>

namespace Nektar
{
    template<typename DataType, typename InnerMatrixType>
    class NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag> : public ConstMatrix<typename NekMatrix<DataType, InnerMatrixType>::NumberType>
    {
        public:
            typedef NekMatrix<DataType, InnerMatrixType> InnerType;
            typedef NekMatrix<InnerType, BlockMatrixTag> ThisType;
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
                    typedef FullMatrixFuncs StoragePolicy;

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
                                m_matrix.GetRows(), m_matrix.GetColumns(), m_curRow, m_curColumn);
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
            };
            
            typedef iterator_base<ThisType> iterator;
            typedef iterator_base<const ThisType> const_iterator;
                                                            
        public:
            explicit NekMatrix(MatrixStorage type = eFULL) :
                BaseType(0,0),
                m_data(),
                m_rowSizes(),
                m_columnSizes(),
                m_storageSize(),
                m_numberOfBlockRows(0),
                m_numberOfBlockColumns(0),
                m_storageType(type)
            {
            }
            
            NekMatrix(unsigned int numberOfBlockRows, unsigned int numberOfBlockColumns,
                      unsigned int rowsPerBlock, unsigned int columnsPerBlock,
                      MatrixStorage type = eFULL) :
                BaseType(numberOfBlockRows*rowsPerBlock, numberOfBlockColumns*columnsPerBlock),
                m_data(),
                m_rowSizes(numberOfBlockRows),
                m_columnSizes(numberOfBlockColumns),
                m_storageSize(0),
                m_numberOfBlockRows(numberOfBlockRows),
                m_numberOfBlockColumns(numberOfBlockColumns),
                m_storageType(type)
            {
                m_storageSize = GetRequiredStorageSize();
                m_data = Array<OneD, boost::shared_ptr<InnerType> >(m_storageSize, boost::shared_ptr<InnerType>());
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
                      unsigned int* rowsPerBlock, unsigned int* columnsPerBlock,
                      MatrixStorage type = eFULL) :
                BaseType(std::accumulate(rowsPerBlock, rowsPerBlock + numberOfBlockRows, 0),
                         std::accumulate(columnsPerBlock, columnsPerBlock + numberOfBlockColumns, 0)),
                m_data(),
                m_rowSizes(numberOfBlockRows),
                m_columnSizes(numberOfBlockColumns),
                m_storageSize(0),
                m_numberOfBlockRows(numberOfBlockRows),
                m_numberOfBlockColumns(numberOfBlockColumns),
                m_storageType(type)
            {
                m_storageSize = GetRequiredStorageSize();
                m_data = Array<OneD, boost::shared_ptr<InnerType> >(m_storageSize, boost::shared_ptr<InnerType>());
                Initialize(rowsPerBlock, columnsPerBlock);
            }
            
            NekMatrix(unsigned int numberOfBlockRows, unsigned int numberOfBlockColumns,
                      const Array<OneD, const unsigned int>& rowsPerBlock, const Array<OneD, const unsigned int>& columnsPerBlock,
                      MatrixStorage type = eFULL) :
                BaseType(std::accumulate(rowsPerBlock.data(), rowsPerBlock.data() + numberOfBlockRows, 0),
                         std::accumulate(columnsPerBlock.data(), columnsPerBlock.data() + numberOfBlockColumns, 0)),
                m_data(),
                m_rowSizes(numberOfBlockRows),
                m_columnSizes(numberOfBlockColumns),
                m_storageSize(0),
                m_numberOfBlockRows(numberOfBlockRows),
                m_numberOfBlockColumns(numberOfBlockColumns),
                m_storageType(type)
            {
                m_storageSize = GetRequiredStorageSize();
                m_data = Array<OneD, boost::shared_ptr<InnerType> >(m_storageSize, boost::shared_ptr<InnerType>());
                Initialize(rowsPerBlock.data(), columnsPerBlock.data());
            }

            NekMatrix(const Array<OneD, const unsigned int>& rowsPerBlock,
                      const Array<OneD, const unsigned int>& columnsPerBlock,
                      MatrixStorage type = eFULL) :
                BaseType(std::accumulate(rowsPerBlock.begin(), rowsPerBlock.end(), 0),
                         std::accumulate(columnsPerBlock.begin(), columnsPerBlock.end(), 0)),
                m_data(),
                m_rowSizes(rowsPerBlock.num_elements()),
                m_columnSizes(columnsPerBlock.num_elements()),
                m_storageSize(0),
                m_numberOfBlockRows(rowsPerBlock.num_elements()),
                m_numberOfBlockColumns(columnsPerBlock.num_elements()),
                m_storageType(type)
            {
                m_storageSize = GetRequiredStorageSize();
                m_data = Array<OneD, boost::shared_ptr<InnerType> >(m_storageSize, boost::shared_ptr<InnerType>());
                Initialize(rowsPerBlock.data(), columnsPerBlock.data());
            }
            
            NekMatrix(const ThisType& rhs) :
                BaseType(rhs),
                m_data(rhs.m_data.num_elements()),
                m_rowSizes(rhs.m_rowSizes),
                m_columnSizes(rhs.m_columnSizes),
                m_storageSize(rhs.m_storageSize),
                m_numberOfBlockRows(rhs.m_numberOfBlockRows),
                m_numberOfBlockColumns(rhs.m_numberOfBlockColumns),
                m_storageType(rhs.m_storageType)
            {
                for(unsigned int i = 0; i < rhs.m_data.num_elements(); ++i)
                {
                    m_data[i] = InnerType::CreateWrapper(rhs.m_data[i]);
                }
            }
            
            
            unsigned int GetRequiredStorageSize() const
            {
                return BaseType::GetRequiredStorageSize(this->GetStorageType(), this->GetNumberOfBlockRows(),
                    this->GetNumberOfBlockColumns());
            }
            
            unsigned int CalculateBlockIndex(unsigned int row, unsigned int column) const
            {
                unsigned int blockRows = GetNumberOfBlockRows();
                unsigned int blockCols = GetNumberOfBlockColumns();
                
                if( this->GetTransposeFlag() == 'T' )
                {
                    std::swap(blockRows, blockCols);
                }
                
                return BaseType::CalculateIndex(this->GetStorageType(), 
                    row, column, blockRows, blockCols, this->GetTransposeFlag());
            }
            
            const InnerType* GetBlockPtr(unsigned int row, unsigned int column) const
            {
                ASSERTL2(row < m_numberOfBlockRows, std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockRows) +
                    std::string(" rows"));
                ASSERTL2(column < m_numberOfBlockColumns, std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockColumns) +
                    std::string(" columns"));
                return m_data[CalculateBlockIndex(row, column)].get();
            }
            
            boost::shared_ptr<const InnerType> GetBlock(unsigned int row, unsigned int column) const        
            {
                ASSERTL2(row < m_numberOfBlockRows, std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockRows) +
                    std::string(" rows"));
                ASSERTL2(column < m_numberOfBlockColumns, std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockColumns) +
                    std::string(" columns"));
                return m_data[CalculateBlockIndex(row, column)];
            }
            
            boost::shared_ptr<InnerType>& GetBlock(unsigned int row, unsigned int column)       
            {
                ASSERTL2(row < m_numberOfBlockRows, std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockRows) +
                    std::string(" rows"));
                ASSERTL2(column < m_numberOfBlockColumns, std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockColumns) +
                    std::string(" columns"));
                return m_data[CalculateBlockIndex(row, column)];
            }
            
            void SetBlock(unsigned int row, unsigned int column, boost::shared_ptr<InnerType>& m)
            {
                ASSERTL2(row < m_numberOfBlockRows, std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockRows) +
                    std::string(" rows"));
                ASSERTL2(column < m_numberOfBlockColumns, std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                    std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockColumns) +
                    std::string(" columns"));
                m_data[CalculateBlockIndex(row, column)] = InnerType::CreateWrapper(m);
            }
            
            
            
            
            ConstGetValueType operator()(unsigned int row, unsigned int col) const
            {
                ASSERTL2(row < this->GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetRows()) +
                    std::string(" rows"));
                ASSERTL2(col < this->GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(col) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(this->GetColumns()) +
                    std::string(" columns"));
                
                
                const Array<OneD, unsigned int>* rowSizes = &m_rowSizes;
                const Array<OneD, unsigned int>* columnSizes = &m_columnSizes;
                
                if( this->GetTransposeFlag() == 'T' )
                {
                    std::swap(rowSizes, columnSizes);
                }
                
                unsigned int blockRow = std::lower_bound(rowSizes->begin(), rowSizes->end(), row) - rowSizes->begin();
                unsigned int blockColumn = std::lower_bound(columnSizes->begin(), columnSizes->end(), col) - columnSizes->begin();
                const boost::shared_ptr<const InnerType> block = GetBlock(blockRow, blockColumn);
                
                unsigned int actualRow = row;
                if( blockRow > 0 )
                {
                    actualRow = row-((*rowSizes)[blockRow-1])-1;
                }

                unsigned int actualCol = col;
                if( blockColumn > 0 )
                {
                    actualCol = col-((*columnSizes)[blockColumn-1])-1;
                }
                    
                
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
            
            MatrixStorage GetType() const
            {
                return m_storageType;
            }
            
            unsigned int GetNumberOfBlockRows() const 
            {
                if( this->GetTransposeFlag() == 'N' )
                {
                    return m_numberOfBlockRows; 
                }
                else
                {
                    return m_numberOfBlockColumns;
                }
            }
            
            unsigned int GetNumberOfBlockColumns() const 
            {
                if( this->GetTransposeFlag() == 'N' )
                {
                    return m_numberOfBlockColumns; 
                }
                else
                {
                    return m_numberOfBlockRows;
                }
            }
            
            unsigned int GetNumberOfRowsInBlockRow(unsigned int blockRow) const
            {
                if( this->GetTransposeFlag() == 'N' )
                {
                    return GetNumberOfElementsInBlock(blockRow, m_numberOfBlockRows, m_rowSizes);
                }
                else
                {
                    return GetNumberOfElementsInBlock(blockRow, m_numberOfBlockColumns, m_columnSizes);
                }
            }

            unsigned int GetNumberOfColumnsInBlockColumn(unsigned int blockCol) const
            {
                if( this->GetTransposeFlag() == 'T' )
                {
                    return GetNumberOfElementsInBlock(blockCol, m_numberOfBlockRows, m_rowSizes);
                }
                else
                {
                    return GetNumberOfElementsInBlock(blockCol, m_numberOfBlockColumns, m_columnSizes);
                }
            }

            iterator begin() { return iterator(*this, 0, 0); }
            iterator end() { return iterator(*this); }
            const_iterator begin() const { return const_iterator(*this, 0, 0); }
            const_iterator end() const { return const_iterator(*this); }
                    
        public:
            static ThisType CreateWrapper(const ThisType& rhs)
            {
                return ThisType(rhs);
            }
            
            static boost::shared_ptr<ThisType> CreateWrapper(const boost::shared_ptr<ThisType>& rhs)
            {
                return boost::shared_ptr<ThisType>(new ThisType(*rhs));
            }

        private:
            
            static unsigned int GetNumberOfElementsInBlock(unsigned int block, unsigned int totalBlocks, const Array<OneD, unsigned int>& sizes)
            {
                ASSERTL2(block < totalBlocks, std::string("Block Element ") + boost::lexical_cast<std::string>(block) + 
                    std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(totalBlocks) +
                    std::string(" blocks."));
                if( block == 0 )
                {
                    return sizes[block]+1;
                }
                else
                {
                    return sizes[block] - sizes[block-1];
                }
            }
            
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
                return this->GetType();
            }
            
            virtual void v_Transpose()
            {
                BOOST_FOREACH(boost::shared_ptr<InnerType> ptr, m_data)
                {
                    if( ptr.get() != 0 )
                    {
                        ptr->Transpose();
                    }
                }
            }
            
            //Array<TwoD, boost::shared_ptr<InnerType> > m_data;
            Array<OneD, boost::shared_ptr<InnerType> > m_data;
            Array<OneD, unsigned int> m_rowSizes;
            Array<OneD, unsigned int> m_columnSizes;
            unsigned int m_storageSize;
            unsigned int m_numberOfBlockRows;
            unsigned int m_numberOfBlockColumns;
            MatrixStorage m_storageType; 
            static NumberType m_zeroElement;
    };

    template<typename DataType, typename InnerMatrixType>
    typename NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::NumberType
    NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::m_zeroElement(0);
    
    template<typename InnerMatrixType>
    NekMatrix<InnerMatrixType, BlockMatrixTag>
    Transpose(NekMatrix<InnerMatrixType, BlockMatrixTag>& rhs)
    {
        NekMatrix<InnerMatrixType, BlockMatrixTag> result(rhs);
        result.Transpose();
        return result;
    }
}


#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BLOCK_MATRIX_HPP
