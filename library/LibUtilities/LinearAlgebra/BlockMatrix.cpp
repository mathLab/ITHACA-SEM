///////////////////////////////////////////////////////////////////////////////
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

#include <LibUtilities/LinearAlgebra/StandardMatrix.hpp>
#include <LibUtilities/LinearAlgebra/ScaledMatrix.hpp>
#include <LibUtilities/LinearAlgebra/BlockMatrix.hpp>


namespace Nektar
{
    template<typename DataType, typename InnerMatrixType>
    NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::NekMatrix(MatrixStorage type) :
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

    template<typename DataType, typename InnerMatrixType>
    NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::NekMatrix(unsigned int numberOfBlockRows, unsigned int numberOfBlockColumns,
              unsigned int rowsPerBlock, unsigned int columnsPerBlock,
              MatrixStorage type) :
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

    template<typename DataType, typename InnerMatrixType>
    NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::NekMatrix(unsigned int numberOfBlockRows, unsigned int numberOfBlockColumns,
              const unsigned int* rowsPerBlock, const unsigned int* columnsPerBlock,
              MatrixStorage type) :
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

    template<typename DataType, typename InnerMatrixType>
    NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::NekMatrix(unsigned int numberOfBlockRows, unsigned int numberOfBlockColumns,
              const Array<OneD, const unsigned int>& rowsPerBlock, const Array<OneD, const unsigned int>& columnsPerBlock,
              MatrixStorage type) :
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

    template<typename DataType, typename InnerMatrixType>
    NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::NekMatrix(const Array<OneD, const unsigned int>& rowsPerBlock,
              const Array<OneD, const unsigned int>& columnsPerBlock,
              MatrixStorage type) :
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

    template<typename DataType, typename InnerMatrixType>
    NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::NekMatrix(const NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>& rhs) :
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

    template<typename DataType, typename InnerMatrixType>
    unsigned int NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::GetRequiredStorageSize() const
    {
        return BaseType::GetRequiredStorageSize(this->GetStorageType(), this->GetNumberOfBlockRows(),
            this->GetNumberOfBlockColumns());
    }

    template<typename DataType, typename InnerMatrixType>
    unsigned int NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::CalculateBlockIndex(unsigned int row, unsigned int column) const
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

    template<typename DataType, typename InnerMatrixType>
    const typename NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::InnerType*
    NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::GetBlockPtr(unsigned int row, unsigned int column) const
    {
        ASSERTL2(row < m_numberOfBlockRows, std::string("Row ") + boost::lexical_cast<std::string>(row) +
            std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockRows) +
            std::string(" rows"));
        ASSERTL2(column < m_numberOfBlockColumns, std::string("Column ") + boost::lexical_cast<std::string>(column) +
            std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockColumns) +
            std::string(" columns"));
        int x = CalculateBlockIndex(row,column);
        if (x == -1)
        {
            return 0;
        }
        else
        {
            return m_data[x].get();
        }
    }

    template<typename DataType, typename InnerMatrixType>
    boost::shared_ptr<const typename NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::InnerType>
    NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::GetBlock(unsigned int row, unsigned int column) const
    {
        ASSERTL2(row < m_numberOfBlockRows, std::string("Row ") + boost::lexical_cast<std::string>(row) +
            std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockRows) +
            std::string(" rows"));
        ASSERTL2(column < m_numberOfBlockColumns, std::string("Column ") + boost::lexical_cast<std::string>(column) +
            std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockColumns) +
            std::string(" columns"));
        int x = CalculateBlockIndex(row,column);
        if (x < 0)
        {
            return boost::shared_ptr<const InnerType>();
        }
        else
        {
            return m_data[x];
        }
    }

    template<typename DataType, typename InnerMatrixType>
    boost::shared_ptr<typename NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::InnerType>&
    NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::GetBlock(unsigned int row, unsigned int column)
    {
        ASSERTL2(row < m_numberOfBlockRows, std::string("Row ") + boost::lexical_cast<std::string>(row) +
            std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockRows) +
            std::string(" rows"));
        ASSERTL2(column < m_numberOfBlockColumns, std::string("Column ") + boost::lexical_cast<std::string>(column) +
            std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockColumns) +
            std::string(" columns"));
        int x = CalculateBlockIndex(row,column);
        if (x == -1)
        {
            return m_nullBlockPtr;
        }
        else
        {
            return m_data[x];
        }
    }

    template<typename DataType, typename InnerMatrixType>
    void NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::SetBlock(unsigned int row, unsigned int column, boost::shared_ptr<InnerType>& m)
    {
        ASSERTL2(row < m_numberOfBlockRows, std::string("Row ") + boost::lexical_cast<std::string>(row) +
            std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockRows) +
            std::string(" rows"));
        ASSERTL2(column < m_numberOfBlockColumns, std::string("Column ") + boost::lexical_cast<std::string>(column) +
            std::string(" requested in a block matrix with a maximum of ") + boost::lexical_cast<std::string>(m_numberOfBlockColumns) +
            std::string(" columns"));
        m_data[CalculateBlockIndex(row, column)] = InnerType::CreateWrapper(m);
    }



    template<typename DataType, typename InnerMatrixType>
    typename NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::ConstGetValueType
    NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::operator()(unsigned int row, unsigned int col) const
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

    template<typename DataType, typename InnerMatrixType>
    unsigned int NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::GetStorageSize() const
    {
        return m_storageSize;
    }

    template<typename DataType, typename InnerMatrixType>
    MatrixStorage NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::GetType() const
    {
        return m_storageType;
    }

    template<typename DataType, typename InnerMatrixType>
    unsigned int NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::GetNumberOfBlockRows() const
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

    template<typename DataType, typename InnerMatrixType>
    unsigned int NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::GetNumberOfBlockColumns() const
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

    template<typename DataType, typename InnerMatrixType>
    unsigned int NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::GetNumberOfRowsInBlockRow(unsigned int blockRow) const
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

    template<typename DataType, typename InnerMatrixType>
    unsigned int NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::GetNumberOfColumnsInBlockColumn(unsigned int blockCol) const
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

    template<typename DataType, typename InnerMatrixType>
    typename NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::iterator NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::begin() { return iterator(*this, 0, 0); }

    template<typename DataType, typename InnerMatrixType>
    typename NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::iterator NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::end() { return iterator(*this); }

    template<typename DataType, typename InnerMatrixType>
    typename NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::const_iterator NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::begin() const { return const_iterator(*this, 0, 0); }

    template<typename DataType, typename InnerMatrixType>
    typename NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::const_iterator NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::end() const { return const_iterator(*this); }


    template<typename DataType, typename InnerMatrixType>
    NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag> NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::CreateWrapper(const NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>& rhs)
    {
        return ThisType(rhs);
    }

    template<typename DataType, typename InnerMatrixType>
    boost::shared_ptr<NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag> > NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::CreateWrapper(const boost::shared_ptr<NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag> >& rhs)
    {
        return boost::shared_ptr<ThisType>(new ThisType(*rhs));
    }


    template<typename DataType, typename InnerMatrixType>
    unsigned int NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::GetNumberOfElementsInBlock(unsigned int block, unsigned int totalBlocks, const Array<OneD, unsigned int>& sizes)
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

    template<typename DataType, typename InnerMatrixType>
    void NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::Initialize(const unsigned int* rowsPerBlock, const unsigned int* columnsPerBlock)
    {
        m_storageSize = this->GetRows()*this->GetColumns();
        if (this->GetRows() > 0)
        {
            m_rowSizes[0] = rowsPerBlock[0] - 1;
            for(unsigned int i = 1; i < m_numberOfBlockRows; ++i)
            {
                m_rowSizes[i] = rowsPerBlock[i] + m_rowSizes[i-1];
            }
        }
        if (this->GetColumns() > 0)
        {
            m_columnSizes[0] = columnsPerBlock[0] - 1;
            for(unsigned int i = 1; i < m_numberOfBlockColumns; ++i)
            {
                m_columnSizes[i] = columnsPerBlock[i] + m_columnSizes[i-1];
            }
        }
    }

    template<typename DataType, typename InnerMatrixType>
    typename boost::call_traits<typename NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::NumberType>::value_type
    NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::v_GetValue(unsigned int row, unsigned int column) const
    {
        return (*this)(row, column);
    }

    template<typename DataType, typename InnerMatrixType>
    unsigned int NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::v_GetStorageSize() const
    {
        return this->GetStorageSize();
    }

    template<typename DataType, typename InnerMatrixType>
    MatrixStorage NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::v_GetStorageType() const
    {
        return this->GetType();
    }

    template<typename DataType, typename InnerMatrixType>
    void NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::v_Transpose()
    {
        BOOST_FOREACH(boost::shared_ptr<InnerType> ptr, m_data)
        {
            if( ptr.get() != 0 )
            {
                ptr->Transpose();
            }
        }
    }


//    template<typename DataType, typename InnerMatrixType>
//    template<typename MatrixType>
//    NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::iterator_base<MatrixType>::iterator_base(MatrixType& m, unsigned int curRow, unsigned int curCol) :
//        m_matrix(m),
//        m_curRow(curRow),
//        m_curColumn(curCol)
//    {
//    }

//    template<typename DataType, typename InnerMatrixType>
//    template<typename MatrixType>
//    NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::iterator_base<MatrixType>::iterator_base(MatrixType& m) :
//        m_matrix(m),
//        m_curRow(std::numeric_limits<unsigned int>::max()),
//        m_curColumn(std::numeric_limits<unsigned int>::max())
//    {
//    }

//    template<typename DataType, typename InnerMatrixType>
//    template<typename MatrixType>
//    NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::iterator_base<MatrixType>::iterator_base(const NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::iterator_base<MatrixType>& rhs) :
//        m_matrix(rhs.m_matrix),
//        m_curRow(rhs.m_curRow),
//        m_curColumn(rhs.m_curColumn)
//    {
//    }

//    template<typename DataType, typename InnerMatrixType>
//    template<typename MatrixType>
//    void NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::iterator_base<MatrixType>::operator++()
//    {
//        if( m_curRow != std::numeric_limits<unsigned int>::max() )
//        {
//            boost::tie(m_curRow, m_curColumn) = StoragePolicy::Advance(
//                m_matrix.GetRows(), m_matrix.GetColumns(), m_curRow, m_curColumn);
//        }
//    }

//    template<typename DataType, typename InnerMatrixType>
//    template<typename MatrixType>
//    typename NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::NumberType
//    NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::iterator_base<MatrixType>::operator*()
//    {
//        return m_matrix(m_curRow, m_curColumn);
//    }

//    template<typename DataType, typename InnerMatrixType>
//    template<typename MatrixType>
//    bool NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::iterator_base<MatrixType>::operator==(const NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::iterator_base<MatrixType>& rhs)
//    {
//        return m_curRow == rhs.m_curRow && m_curColumn == rhs.m_curColumn;
//    }

//    template<typename DataType, typename InnerMatrixType>
//    template<typename MatrixType>
//    bool NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::iterator_base<MatrixType>::operator!=(const NekMatrix<NekMatrix<DataType, InnerMatrixType>, BlockMatrixTag>::iterator_base<MatrixType>& rhs)
//    {
//        return !(*this == rhs);
//    }


    template LIB_UTILITIES_EXPORT class NekMatrix< NekMatrix<NekDouble, StandardMatrixTag>, BlockMatrixTag>;
    template LIB_UTILITIES_EXPORT class NekMatrix<NekMatrix< NekMatrix<NekDouble, StandardMatrixTag>, BlockMatrixTag>, BlockMatrixTag>;
    template LIB_UTILITIES_EXPORT class NekMatrix<NekMatrix< NekMatrix<NekDouble, StandardMatrixTag>, ScaledMatrixTag>, BlockMatrixTag>;
    template LIB_UTILITIES_EXPORT class NekMatrix<NekMatrix<NekMatrix< NekMatrix<NekDouble, StandardMatrixTag>, ScaledMatrixTag>, BlockMatrixTag>, BlockMatrixTag>;
}

