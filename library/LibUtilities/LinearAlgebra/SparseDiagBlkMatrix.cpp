///////////////////////////////////////////////////////////////////////////////
//
// File: SparseDiagBlkMatrix.hpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
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
// Description: Diagonal block sparse matrix class templated by underlying sparse
//              storage format
//
///////////////////////////////////////////////////////////////////////////////

#include <map>
#include <vector>
#include <utility>
#include <algorithm>
#include <fstream>

#include <LibUtilities/LinearAlgebra/MatrixStorageType.h>
#include <LibUtilities/LinearAlgebra/NistSparseDescriptors.hpp>
#include <LibUtilities/LinearAlgebra/SparseMatrixFwd.hpp>
#include <LibUtilities/LinearAlgebra/SparseDiagBlkMatrix.hpp>
#include <LibUtilities/LinearAlgebra/StorageSmvBsr.hpp>

namespace Nektar
{

    template<typename SparseStorageType>
    NekSparseDiagBlkMatrix<SparseStorageType>::
                NekSparseDiagBlkMatrix(const SparseStorageSharedPtrVector&
                                             sparseStoragePtrVector):
        m_rows(0),
        m_cols(0),
        m_rowoffset(sparseStoragePtrVector.size()+1, 0.0),
        m_mulCallsCounter(0),
        m_submatrix(sparseStoragePtrVector.size(),sparseStoragePtrVector)
    {
        for (int i = 0; i < sparseStoragePtrVector.size(); i++)
        {
            const IndexType rows    = sparseStoragePtrVector[i]->GetRows();
            m_rows += rows;
            m_cols += sparseStoragePtrVector[i]->GetColumns();
            m_rowoffset[i+1] = m_rowoffset[i] + rows;
        }
    }

    template<typename SparseStorageType>
    NekSparseDiagBlkMatrix<SparseStorageType>::NekSparseDiagBlkMatrix(const NekSparseDiagBlkMatrix& src):
        m_rows(src.m_rows),
        m_cols(src.m_cols),
        m_rowoffset(src.m_rowoffset),
        m_mulCallsCounter(src.m_mulCallsCounter),
        m_submatrix(src.m_submatrix)
    {
    }

    template<typename SparseStorageType>
    NekSparseDiagBlkMatrix<SparseStorageType>::~NekSparseDiagBlkMatrix()
    {
    }

    template<typename SparseStorageType>
    IndexType NekSparseDiagBlkMatrix<SparseStorageType>::GetRows() const
    {
        return m_rows;
    }

    template<typename SparseStorageType>
    IndexType NekSparseDiagBlkMatrix<SparseStorageType>::GetColumns() const
    {
        return m_cols;
    }

    /// number of rows at i-th submatrix
    template<typename SparseStorageType>
    IndexType NekSparseDiagBlkMatrix<SparseStorageType>::GetRows(int i) const
    {
        return m_submatrix[i]->GetRows();
    }

    /// number of columns at i-th submatrix
    template<typename SparseStorageType>
    IndexType NekSparseDiagBlkMatrix<SparseStorageType>::GetColumns(int i) const
    {
        return m_submatrix[i]->GetColumns();
    }

    template<typename SparseStorageType>
    IndexType NekSparseDiagBlkMatrix<SparseStorageType>::GetNumberOfMatrixBlocks() const
    {
        return m_submatrix.size();
    }

    template<typename SparseStorageType>
    IndexType NekSparseDiagBlkMatrix<SparseStorageType>::GetNumNonZeroEntries()
    {
        IndexType nnz = 0;
        for (int i = 0; i < m_submatrix.size(); i++)
        {
            nnz += m_submatrix[i]->GetNumNonZeroEntries();
        }
        return nnz;
    }

    // nnz of i-th CSR matrix
    template<typename SparseStorageType>
    IndexType NekSparseDiagBlkMatrix<SparseStorageType>::GetNumNonZeroEntries(int i) const
    {
        return m_submatrix[i]->GetNumNonZeroEntries();
    }

    template<typename SparseStorageType>
    typename NekSparseDiagBlkMatrix<SparseStorageType>::DataType NekSparseDiagBlkMatrix<SparseStorageType>::GetFillInRatio(int i) const
    {
        return m_submatrix[i]->GetFillInRatio();
    }

    template<typename SparseStorageType>
    typename NekSparseDiagBlkMatrix<SparseStorageType>::DataType NekSparseDiagBlkMatrix<SparseStorageType>::GetFillInRatio() const
    {
        IndexType stored = 0;
        IndexType nnz = 0;
        for (int i = 0; i < m_submatrix.size(); i++)
        {
            stored += m_submatrix[i]->GetNumStoredDoubles();
            nnz    += m_submatrix[i]->GetNumNonZeroEntries();
        }
        return (DataType)stored/(DataType)nnz;
    }


    template<typename SparseStorageType>
    typename boost::call_traits<typename SparseStorageType::DataType>::const_reference NekSparseDiagBlkMatrix<SparseStorageType>::operator()(IndexType block, IndexType loc_row, IndexType loc_column) const
    {
        return m_submatrix[block]->GetValue(loc_row, loc_column);
    }

    template<typename SparseStorageType>
    typename boost::call_traits<typename SparseStorageType::DataType>::const_reference
            NekSparseDiagBlkMatrix<SparseStorageType>::operator()
                         (IndexType glob_row, IndexType glob_column) const
    {
        IndexType i = 0;
        static DataType defaultReturnValue = 0;

        signed int local_row = glob_row;
        signed int local_col = glob_column;

        while ((local_row >= m_submatrix[i]->GetRows()) &&
               (local_col >= 0))
        {
            local_row -= m_submatrix[i]->GetRows();
            local_col -= m_submatrix[i]->GetColumns();
            i++;
        }
        /// \todo double check, might be a bug when local_col > local column
        if (local_col < 0) return defaultReturnValue;

        return m_submatrix[i]->GetValue(local_row, local_col);
    }


    template<typename SparseStorageType>
    void NekSparseDiagBlkMatrix<SparseStorageType>::Multiply(
                  const DataVectorType  &in,
                        DataVectorType &out)
    {
        for (int i = 0; i < m_submatrix.size(); ++i)
        {
            m_submatrix[i]->Multiply(&in[m_rowoffset[i]], &out[m_rowoffset[i]]);
        }
        m_mulCallsCounter++;
    }

    template<typename SparseStorageType>
    void NekSparseDiagBlkMatrix<SparseStorageType>::Multiply(
                                  const DataType* in, DataType* out)
    {
        for (int i = 0; i < m_submatrix.size(); ++i)
        {
            m_submatrix[i]->Multiply(&in[m_rowoffset[i]], &out[m_rowoffset[i]]);
        }
        m_mulCallsCounter++;
    }

    template<typename SparseStorageType>
    void NekSparseDiagBlkMatrix<SparseStorageType>::
        MultiplySubMatrix(const IndexType blockNum,
                                DataType* in,
                                DataType* out)
    {
        m_submatrix[blockNum]->Multiply(in + m_rowoffset[blockNum], out + m_rowoffset[blockNum]);
        // m_mulCallsCounter++;
    }

    template<typename SparseStorageType>
    size_t NekSparseDiagBlkMatrix<SparseStorageType>::GetMemoryFootprint()
    {
        size_t bytes =
                       sizeof(IndexType)*2 +      // sizes
                       sizeof(unsigned long)  +   // mulCallsCounter
                       sizeof(IndexVector) +
                       sizeof(IndexType)*m_rowoffset.capacity() +
                       sizeof(SparseStorageSharedPtrVector) +
                       sizeof(SparseStorageSharedPtr)*m_submatrix.capacity();
        for (int i = 0; i < m_submatrix.size(); i++)
        {
            bytes += m_submatrix[i]->GetMemoryUsage();
        }
        return bytes;
    }

    template<typename SparseStorageType>
    size_t NekSparseDiagBlkMatrix<SparseStorageType>::GetMemoryFootprint(IndexType i) const
    {
        return m_submatrix[i]->GetMemoryUsage();
    }

    template<typename SparseStorageType>
    unsigned long NekSparseDiagBlkMatrix<SparseStorageType>::GetMulCallsCounter() const
    {
        return m_mulCallsCounter;
    }

    template<typename SparseStorageType>
    typename SparseStorageType::DataType  NekSparseDiagBlkMatrix<SparseStorageType>::GetAvgRowDensity()
    {
        DataType avgRowDensity = 0.0;
        for (int i = 0; i < m_submatrix.size(); i++)
        {
            avgRowDensity += (DataType) m_submatrix[i]->GetNumNonZeroEntries() /
                             (DataType) m_submatrix[i]->GetRows();
        }
        return avgRowDensity / m_submatrix.size();
    }

    template<typename SparseStorageType>
    typename SparseStorageType::DataType  NekSparseDiagBlkMatrix<SparseStorageType>::GetAvgRowDensity(IndexType  i) const
    {
        return (DataType) m_submatrix[i]->GetNumNonZeroEntries() /
               (DataType) m_submatrix[i]->GetRows();
    }

    template<typename SparseStorageType>
    IndexType NekSparseDiagBlkMatrix<SparseStorageType>::GetBandwidth()
    {
        IndexType bandwidth = 0;
        for (int i = 0; i < m_submatrix.size(); i++)
        {
            typename SparseStorageType::const_iterator entry = m_submatrix[i]->begin();
            for (; entry != m_submatrix[i]->end(); ++entry)
            {
                bandwidth = (std::max)(static_cast<int>(bandwidth),
                    2*abs( static_cast<int>(entry->first.first - entry->first.second)+1));
            }
        }
        return bandwidth;
    }

    template<typename SparseStorageType>
    IndexType NekSparseDiagBlkMatrix<SparseStorageType>::GetBandwidth(IndexType  i)
    {
        IndexType bandwidth = 0;
        typename SparseStorageType::const_iterator entry = m_submatrix[i]->begin();
        for (; entry != m_submatrix[i]->end(); ++entry)
        {
            bandwidth = (std::max)(static_cast<int>(bandwidth),
                2*abs( static_cast<int>(entry->first.first - entry->first.second)+1));
        }
        return bandwidth;
    }

    template<typename SparseStorageType>
    COOMatTypeSharedPtr NekSparseDiagBlkMatrix<SparseStorageType>::GetCooStorage(IndexType i)
    {
        COOMatTypeSharedPtr coo (new COOMatType());
        typename SparseStorageType::const_iterator entry = m_submatrix[i]->begin();
        for (; entry != m_submatrix[i]->end(); entry++)
        {
            IndexType loc_row = entry->first.first;
            IndexType loc_col = entry->first.second;
            (*coo)[std::make_pair(loc_row, loc_col) ] = entry->second;
        }
        return coo;
    }

    template<typename SparseStorageType>
    COOMatTypeSharedPtr NekSparseDiagBlkMatrix<SparseStorageType>::GetCooStorage()
    {
        COOMatTypeSharedPtr coo (new COOMatType());
        IndexType row_offset = 0;
        IndexType col_offset = 0;
        for (IndexType i = 0; i < m_submatrix.size(); i++)
        {
            typename SparseStorageType::const_iterator entry = m_submatrix[i]->begin();
            for (; entry != m_submatrix[i]->end(); entry++)
            {
                IndexType loc_row = entry->first.first;
                IndexType loc_col = entry->first.second;
                (*coo)[std::make_pair(loc_row + row_offset, loc_col + col_offset) ] = entry->second;
            }
            row_offset += m_submatrix[i]->GetRows();
            col_offset += m_submatrix[i]->GetColumns();
        }
        return coo;
    }

    // Generates a matrix each element of which is a number of
    // nonzero entries to block-matrix associated with *this object.
    // E.g. for unit 6x6 matrix and block size = 2 it will generate
    // 3x3 matrix with elements 2 along the diagonal.
    //
    // The output can be visualised with Python's matplotlib's spy().
    //
    template<typename SparseStorageType>
    void NekSparseDiagBlkMatrix<SparseStorageType>::
            writeSubmatrixSparsityPatternTo(std::ostream& out,
                                   const IndexType subMatrixIdx,
                                   const IndexType blockSize)
    {
        const int matRows = m_submatrix[subMatrixIdx]->GetRows();
        const int gridRows = matRows / blockSize + (matRows % blockSize > 0);
        const int gridCols = gridRows;

        std::vector< std::vector<int> > grid (gridRows);
        for (int row = 0; row < gridRows; row++)
        {
            grid[row].resize(gridCols,0.0);
        }

        typename SparseStorageType::const_iterator entry =
                                    m_submatrix[subMatrixIdx]->begin();
        typename SparseStorageType::const_iterator stop  =
                                    m_submatrix[subMatrixIdx]->end();
        for (; entry != stop; ++entry)
        {
            const IndexType row = entry->first.first;
            const IndexType col = entry->first.second;
            const int gridRow = row / blockSize;
            const int gridCol = col / blockSize;
            grid[gridRow][gridCol]++;
        }

        for (int row = 0; row < gridRows; row++)
        {
            for (int col = 0; col < gridCols; col++)
            {
                out << grid[row][col] << " ";
            }
            out << std::endl;
        }
    }

    template<typename SparseStorageType>
    void NekSparseDiagBlkMatrix<SparseStorageType>::
            writeSparsityPatternTo(std::ostream& out,
                                   const IndexType blockSize)
    {
        const int matRows = GetRows();
        const int gridRows = matRows / blockSize + (matRows % blockSize > 0);
        const int gridCols = gridRows;

        std::vector< std::vector<int> > grid (gridRows);
        for (int row = 0; row < gridRows; row++)
        {
            grid[row].resize(gridCols,0.0);
        }

        IndexType row_offset = 0;
        IndexType col_offset = 0;
        for (int i = 0; i < m_submatrix.size(); i++)
        {
            typename SparseStorageType::const_iterator entry = m_submatrix[i]->begin();
            typename SparseStorageType::const_iterator stop  = m_submatrix[i]->end();
            for (; entry != stop; ++entry)
            {
                const IndexType row = entry->first.first  + row_offset;
                const IndexType col = entry->first.second + col_offset;
                const int gridRow = row / blockSize;
                const int gridCol = col / blockSize;
                grid[gridRow][gridCol]++;
            }
            row_offset += m_submatrix[i]->GetRows();
            col_offset += m_submatrix[i]->GetColumns();
        }

        for (int row = 0; row < gridRows; row++)
        {
            for (int col = 0; col < gridCols; col++)
            {
                out << grid[row][col] << " ";
            }
            out << std::endl;
        }
    }

    /// Complementary routine to the previous. It generates
    /// exact non-zero pattern of a given block matrix entry
    /// to *this object. E.g., for a 6x6 matrix defined as
    /// A(i,i):=i, A(i,j):=0, block size = 2 and
    /// blk_row=blk_col=2 it produces 2x2 matrix with 5 and 6
    /// along the diagonal.
    template<typename SparseStorageType>
    void NekSparseDiagBlkMatrix<SparseStorageType>::
            writeSubmatrixBlockSparsityPatternTo(std::ostream& out,
                                        const IndexType subMatrixIdx,
                                        const IndexType blk_row,
                                        const IndexType blk_col,
                                              IndexType blockSize)
    {
        blockSize = (std::min)(blockSize, m_submatrix[subMatrixIdx]->GetRows());
        std::vector< std::vector<int> > grid (blockSize);
        for (int row = 0; row < blockSize; row++)
        {
            grid[row].resize(blockSize,0.0);
        }

        typename SparseStorageType::const_iterator entry =
                                    m_submatrix[subMatrixIdx]->begin();
        typename SparseStorageType::const_iterator stop  =
                                    m_submatrix[subMatrixIdx]->end();
        for (; entry != stop; ++entry)
        {
            const IndexType row = entry->first.first;
            const IndexType col = entry->first.second;

            if (blk_row != row / blockSize ) continue;
            if (blk_col != col / blockSize ) continue;
            grid[row % blockSize][col % blockSize]++;
        }

        for (int row = 0; row < blockSize; row++)
        {
            for (int col = 0; col < blockSize; col++)
            {
                out << grid[row][col] << " ";
            }
            out << std::endl;
        }
    }

    template<typename SparseStorageType>
    void NekSparseDiagBlkMatrix<SparseStorageType>::writeBlockSparsityPatternTo(
                std::ostream& out,
                const IndexType blk_row,
                const IndexType blk_col,
                IndexType blockSize)
    {
        blockSize = (std::min)(blockSize, GetRows());
        std::vector< std::vector<int> > grid (blockSize);
        for (int row = 0; row < blockSize; row++)
        {
            grid[row].resize(blockSize,0.0);
        }

        IndexType row_offset = 0;
        IndexType col_offset = 0;
        for (int i = 0; i < m_submatrix.size(); i++)
        {
            typename SparseStorageType::const_iterator entry =
                                        m_submatrix[i]->begin();
            typename SparseStorageType::const_iterator stop  =
                                        m_submatrix[i]->end();
            for (; entry != stop; ++entry)
            {
                const IndexType row = entry->first.first  + row_offset;
                const IndexType col = entry->first.second + col_offset;

                if (blk_row != row / blockSize ) continue;
                if (blk_col != col / blockSize ) continue;
                grid[row % blockSize][col % blockSize]++;
            }
            row_offset += m_submatrix[i]->GetRows();
            col_offset += m_submatrix[i]->GetColumns();
        }

        for (int row = 0; row < blockSize; row++)
        {
            for (int col = 0; col < blockSize; col++)
            {
                out << grid[row][col] << " ";
            }
            out << std::endl;
        }
    }



    // explicit instantiation
    template class NekSparseDiagBlkMatrix<StorageSmvBsr<NekDouble> >;


} // namespace

