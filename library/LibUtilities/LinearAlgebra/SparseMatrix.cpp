///////////////////////////////////////////////////////////////////////////////
//
// File: SparseMatrix.cpp
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
// Description: Generic sparse matrix class templated by underlying sparse
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
#include <LibUtilities/LinearAlgebra/SparseMatrix.hpp>
#include <LibUtilities/LinearAlgebra/StorageSmvBsr.hpp>

using std::min;
using std::max;

namespace Nektar
{

    template<typename SparseStorageType>
    NekSparseMatrix<SparseStorageType>::NekSparseMatrix(const SparseStorageSharedPtr& sparseStoragePtr):
        m_mulCallsCounter(0),
        m_sparseStorage(sparseStoragePtr)
    {
    }

    template<typename SparseStorageType>
    NekSparseMatrix<SparseStorageType>::NekSparseMatrix(const NekSparseMatrix& src):
        m_mulCallsCounter(src.m_mulCallsCounter),
        m_sparseStorage(src.m_sparseStorage)
    {
    }

    template<typename SparseStorageType>
    NekSparseMatrix<SparseStorageType>::~NekSparseMatrix()
    {
    }


    template<typename SparseStorageType>
    IndexType NekSparseMatrix<SparseStorageType>::GetRows() const
    {
        return m_sparseStorage->GetRows();
    }

    template<typename SparseStorageType>
    IndexType NekSparseMatrix<SparseStorageType>::GetColumns() const
    {
        return m_sparseStorage->GetColumns();
    }

    template<typename SparseStorageType>
    IndexType NekSparseMatrix<SparseStorageType>::GetNumNonZeroEntries() const
    {
        return m_sparseStorage->GetNumNonZeroEntries();
    }


    template<typename SparseStorageType>
    const typename NekSparseMatrix<SparseStorageType>::DataType NekSparseMatrix<SparseStorageType>::GetFillInRatio() const
    {
        return m_sparseStorage->GetFillInRatio();
    }



    template<typename SparseStorageType>
    typename boost::call_traits<typename SparseStorageType::DataType>::const_reference NekSparseMatrix<SparseStorageType>::operator()(unsigned int row, unsigned int column) const
    {
        return m_sparseStorage->GetValue(row, column);
    }

    template<typename SparseStorageType>
    typename SparseStorageType::const_iterator NekSparseMatrix<SparseStorageType>::begin() const
    {
        return m_sparseStorage->begin();
    }

    template<typename SparseStorageType>
    typename SparseStorageType::const_iterator NekSparseMatrix<SparseStorageType>::end() const
    {
        return m_sparseStorage->end();
    }


    template<typename SparseStorageType>
    void NekSparseMatrix<SparseStorageType>::Multiply(const DataVectorType &in,
                        DataVectorType &out)
    {
        m_sparseStorage->Multiply(in,out);
        m_mulCallsCounter++;
    }

    template<typename SparseStorageType>
    void NekSparseMatrix<SparseStorageType>::Multiply(const DataType* in,
                        DataType* out)
    {
        m_sparseStorage->Multiply(in,out);
        m_mulCallsCounter++;
    }

    template<typename SparseStorageType>
    size_t NekSparseMatrix<SparseStorageType>::GetMemoryFootprint() const
    {
        return m_sparseStorage->GetMemoryUsage() +
               sizeof(SparseStorageSharedPtr) +
               sizeof(unsigned long); // mulCallsCounter
    }

    template<typename SparseStorageType>
    unsigned long NekSparseMatrix<SparseStorageType>::GetMulCallsCounter() const
    {
        return m_mulCallsCounter;
    }

    template<typename SparseStorageType>
    const typename SparseStorageType::DataType NekSparseMatrix<SparseStorageType>::GetAvgRowDensity() const
    {
        return (DataType) m_sparseStorage->GetNumNonZeroEntries() /
               (DataType) m_sparseStorage->GetRows();
    }

    template<typename SparseStorageType>
    IndexType NekSparseMatrix<SparseStorageType>::GetBandwidth()
    {
        int bandwidth = 0;

        typename SparseStorageType::const_iterator entry = m_sparseStorage->begin();

        for (; entry != m_sparseStorage->end(); ++entry)
        {
            bandwidth = (std::max)(bandwidth, 2*std::abs((int)(entry->first.first - entry->first.second+1)));
        }
        return bandwidth;
    }

    template<typename SparseStorageType>
    COOMatTypeSharedPtr NekSparseMatrix<SparseStorageType>::GetCooStorage()
    {
        COOMatTypeSharedPtr coo (new COOMatType());

        typename SparseStorageType::const_iterator entry = m_sparseStorage->begin();
        for (; entry != m_sparseStorage->end(); entry++)
        {
            (*coo)[std::make_pair(entry->first.first, entry->first.second) ] = entry->second;
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
    void NekSparseMatrix<SparseStorageType>::writeSparsityPatternTo(std::ostream& out, IndexType blockSize)
    {
        const int matRows = m_sparseStorage->GetRows();
        const int gridRows = matRows / blockSize + (matRows % blockSize > 0);
        const int gridCols = gridRows;

        std::vector< std::vector<int> > grid (gridRows);
        for (int row = 0; row < gridRows; row++)
        {
            grid[row].resize(gridCols,0);
        }

        typename SparseStorageType::const_iterator entry = m_sparseStorage->begin();
        for (; entry != m_sparseStorage->end(); entry++)
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

    /// Complementary routine to the previous. It generates
    /// exact non-zero pattern of a given block matrix entry
    /// to *this object. E.g., for a 6x6 matrix defined as
    /// A(i,i):=i, A(i,j):=0, block size = 2 and 
    /// blk_row=blk_col=2 it produces 2x2 matrix with 5 and 6
    /// along the diagonal.
    template<typename SparseStorageType>
    void NekSparseMatrix<SparseStorageType>::writeBlockSparsityPatternTo(std::ostream& out,
                    const IndexType blk_row, const IndexType blk_col, IndexType blockSize)
    {
        blockSize = (std::min)(blockSize, m_sparseStorage->GetRows());
        std::vector< std::vector<int> > grid (blockSize);
        for (int row = 0; row < blockSize; row++)
        {
            grid[row].resize(blockSize,0);
        }

        typename SparseStorageType::const_iterator entry = m_sparseStorage->begin();
        for (; entry != m_sparseStorage->end(); entry++)
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



    // explicit instantiation
    template class NekSparseMatrix<StorageSmvBsr<NekDouble> >;


} // namespace

