///////////////////////////////////////////////////////////////////////////////
//
// File: StorageNistBsr.cpp
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
// Description: Interface to NIST BSR sparse matrix storage.
//
///////////////////////////////////////////////////////////////////////////////

#include <map>
#include <vector>
#include <utility>
#include <fstream>
#include <cmath>

#include <LibUtilities/BasicConst/NektarUnivConsts.hpp>
#include <LibUtilities/LinearAlgebra/NistSparseBlas.hpp>
#include <LibUtilities/LinearAlgebra/StorageNistBsr.hpp>
#include <LibUtilities/LinearAlgebra/NistSparseDescriptors.hpp>

#include <boost/lexical_cast.hpp>

namespace Nektar
{

    template<typename DataType>
    StorageNistBsr<DataType>::const_iterator::const_iterator(
                    MatrixStorage       matType,
                    IndexType           begin,
                    IndexType           end,
                    IndexType           blkDim,
                    const DataVectorType&     val,
                    const IndexVectorType&    indx,
                    const IndexVectorType&    pntr) :
        m_matType(matType),
        m_iter(),
        m_begin(begin),
        m_end(end),
        m_blkDim(blkDim),
        m_val(val),
        m_indx(indx),
        m_pntr(pntr)
    {
        m_iter.nnzindex = begin;

        // Here we rewind the iterator to the 'begin'-th
        // nonzero double value. Consecutive nonzeros are
        // to lie within dense block and when it's done,
        // jump to the next dense block.
        // NB: m_val stores some explicit zeros that need to be skept.
        if (begin < end)
        {
            // determine offset of 'begin'-th nonzero value

            m_iter.storageindex = 0;
            m_iter.nnzindex     = 0;
            while(m_iter.nnzindex < begin)
            {
                // explicit zero?
                if (m_val[m_iter.storageindex] > NekConstants::kNekSparseNonZeroTol)
                {
                    m_iter.nnzindex++;
                }
                m_iter.storageindex++;
                if (m_iter.storageindex >= m_val.num_elements())
                {
                    std::cout << "const_iterator: 'begin' out stored values bounds" << std::endl;
                    throw 1;
                }
            }
            m_iter.second   = m_val[m_iter.storageindex];

            CoordType c = storageIndexToFullCoord(m_iter.storageindex);
            m_iter.first.first  = c.first;
            m_iter.first.second = c.second;
        }
    }

    template<typename DataType>
    CoordType StorageNistBsr<DataType>::const_iterator::storageIndexToFullCoord(IndexType storageIndex)
    {
        // find local row and column indices within this block

        const IndexType elms       = m_blkDim*m_blkDim;
        const IndexType block      = storageIndex / elms;
        const IndexType loc_offset = storageIndex % elms;
        const IndexType loc_row    = loc_offset % m_blkDim;
        const IndexType loc_col    = loc_offset / m_blkDim;

        // find block row and block column

        IndexType block_col = m_indx[block+1]-1;
        IndexType block_row = 0;
        while(block >= m_pntr[block_row+1]-1)
        {
            block_row++;
        }

        // full matrix coordinates of this nonzero entry

        CoordType coord;
        coord.first  = block_row*m_blkDim + loc_row;
        coord.second = block_col*m_blkDim + loc_col;
        return coord;
    }


    template<typename DataType>
    StorageNistBsr<DataType>::const_iterator::const_iterator(const const_iterator& src):
        m_matType(src.m_matType),
        m_iter(src.m_iter),
        m_begin(src.m_begin),
        m_end(src.m_end),
        m_val(src.m_val),
        m_indx(src.m_indx),
        m_pntr(src.m_pntr)
    {
    }

    template<typename DataType>
    StorageNistBsr<DataType>::const_iterator::~const_iterator()
    {
    }

    template<typename DataType>
    typename StorageNistBsr<DataType>::const_iterator StorageNistBsr<DataType>::const_iterator::operator++(int)
    {
        const_iterator out = *this;
        forward();
        return out;
    }

    template<typename DataType>
    typename StorageNistBsr<DataType>::const_iterator& StorageNistBsr<DataType>::const_iterator::operator++()
    {
        forward();
        return *this;
    }

    template<typename DataType>
    const typename StorageNistBsr<DataType>::const_iterator::IterType& StorageNistBsr<DataType>::const_iterator::operator*()
    {
        return m_iter;
    }

    template<typename DataType>
    const typename StorageNistBsr<DataType>::const_iterator::IterType* StorageNistBsr<DataType>::const_iterator::operator->()
    {
        return &m_iter;
    }

    template<typename DataType>
    const bool StorageNistBsr<DataType>::const_iterator::operator==(const const_iterator& rhs)
    {
        return m_iter.nnzindex == rhs.m_iter.nnzindex;
    }

    template<typename DataType>
    const bool StorageNistBsr<DataType>::const_iterator::operator!=(const const_iterator& rhs)
    {
        return !(m_iter.nnzindex == rhs.m_iter.nnzindex);
    }

    template<typename DataType>
    void StorageNistBsr<DataType>::const_iterator::forward()
    {
        while((m_iter.storageindex+1 < m_val.num_elements()) &&
              (m_val[++m_iter.storageindex] <= NekConstants::kNekSparseNonZeroTol));

        m_iter.nnzindex++;

        if (m_iter.storageindex >= m_val.num_elements())
        {
            m_iter.nnzindex = m_end;
            return;
        }

        m_iter.second   = m_val[m_iter.storageindex];

        CoordType c = storageIndexToFullCoord(m_iter.storageindex);
        m_iter.first.first  = c.first;
        m_iter.first.second = c.second;
    }





    template<typename DataType>
    StorageNistBsr<DataType>::StorageNistBsr(
                    const IndexType  blkRows,
                    const IndexType  blkCols,
                    const IndexType  blkDim,
                    const BCOMatType&   bcoMat,
                    const MatrixStorage matType):
        m_matType (matType),
        m_blkRows (blkRows),
        m_blkCols (blkCols),
        m_blkDim  (blkDim),
        m_bnnz    (bcoMat.size()),
        m_nnz     (0),
        m_val     (m_bnnz * blkDim*blkDim),
        m_indx    (m_bnnz+1),
        m_pntr    (blkRows+1)
    {
        if (matType != Nektar::eFULL)
        {
            /// \todo: - iterators over strictly lower-triangular part
            ///        - number off-diagonal elements calculation
            ///        - clear distinction between stored and factual nonzeros
            ///          (row density)
            std::cout << "matrix type not implemented" << std::endl;
            throw 1;
        }

        processBcoInput(blkRows,blkCols,blkDim,bcoMat);
    }


    template<typename DataType>
    StorageNistBsr<DataType>::StorageNistBsr(const StorageNistBsr& src):
        m_matType(src.m_matType),
        m_blkRows (src.m_blkRows),
        m_blkCols (src.m_blkCols),
        m_blkDim(src.m_blkDim),
        m_bnnz(src.m_bnnz),
        m_nnz(src.m_nnz),
        m_val(src.m_val),
        m_indx(src.m_indx),
        m_pntr(src.m_pntr)
    {
    }

    template<typename DataType>
    StorageNistBsr<DataType>::~StorageNistBsr()
    {
    }


    template<typename DataType>
    const IndexType StorageNistBsr<DataType>::GetRows() const
    {
        return m_blkRows*m_blkDim;
    }

    template<typename DataType>
    const IndexType StorageNistBsr<DataType>::GetColumns() const
    {
        return m_blkCols*m_blkDim;
    }

    template<typename DataType>
    const IndexType StorageNistBsr<DataType>::GetNumNonZeroEntries() const
    {
        return m_nnz;
    }

    template<typename DataType>
    const IndexType StorageNistBsr<DataType>::GetBlkSize() const
    {
        return m_blkDim;
    }

    template<typename DataType>
    const IndexType StorageNistBsr<DataType>::GetNumStoredDoubles() const
    {
        return m_bnnz*m_blkDim*m_blkDim;
    }


    template<typename DataType>
    const DataType StorageNistBsr<DataType>::GetFillInRatio() const
    {
        return (DataType)(m_bnnz*m_blkDim*m_blkDim)/(DataType)m_nnz;
    }


    template<typename DataType>
    const size_t StorageNistBsr<DataType>::GetMemoryUsage(IndexType nnz, IndexType nRows) const
    {
        return sizeof(DataType) *m_val.capacity()   +
               sizeof(IndexType)*m_indx.capacity() +
               sizeof(IndexType)*m_pntr.capacity()  +
               sizeof(IndexType)*5   + //< blkRows + blkCols + blkDim + nnz + bnnz
               sizeof(MatrixStorage);
    }


    template<typename DataType>
    const typename boost::call_traits<DataType>::const_reference StorageNistBsr<DataType>::GetValue(IndexType grow, IndexType gcolumn) const
    {
        IndexType  brow = grow    / m_blkDim;
        IndexType  bcol = gcolumn / m_blkDim;
        IndexType  lrow = grow    % m_blkDim;
        IndexType  lcol = gcolumn % m_blkDim;

        // rewind to the first entry of the first
        // block in the current block row
        IndexType  offset = (m_pntr[brow]-1)*m_blkDim*m_blkDim;

        IndexType i;
        static DataType defaultReturnValue;
        for( i = m_pntr[brow]; i < m_pntr[brow+1]; i++)
        {
            if(bcol == m_indx[i]-1)
            {
                return m_val[offset+lrow + lcol*m_blkDim];
            }
            offset += m_blkDim*m_blkDim;
        }

        return defaultReturnValue;
    }


    template<typename DataType>
    typename StorageNistBsr<DataType>::const_iterator StorageNistBsr<DataType>::begin() const
    {
        return const_iterator(m_matType, 0, m_nnz, m_blkDim, m_val, m_indx, m_pntr);
    }

    template<typename DataType>
    typename StorageNistBsr<DataType>::const_iterator StorageNistBsr<DataType>::end() const
    {
        return const_iterator(m_matType, m_nnz, m_nnz, m_blkDim, m_val, m_indx, m_pntr);
    }



    // General non-symmetric zero-based BSR multiply
    // C = A*B where A is matrix and B is vector.
    // No scaling. Previous content of C is discarded.
    template<typename DataType>
    void StorageNistBsr<DataType>::Multiply(const DataVectorType &in,
                                                  DataVectorType &out)
    {
        const IndexType m = m_blkRows*m_blkDim;
        SparseBlas::Dbsrmm(0,           // no transpose
                           m_blkRows,
                           1,           // multiplying to a vector
                           m_blkCols,
                           1.0,         // alpha
                           NistSpBlasDescra[m_matType],
                           &m_val[0], // nonzero entries
                           (int*)&m_indx[0], // column indices
                           (int*)&m_pntr[0],
                           (int*)&m_pntr[0]+1,
                           m_blkDim,
                           &in[0],
                           m,
                           0.0,         // beta
                           &out[0],
                           m,
                           NULL,
                           0);
    }

    template<typename DataType>
    void StorageNistBsr<DataType>::Multiply(const DataType* in,
                                                  DataType* out)
    {
        const IndexType m = m_blkRows*m_blkDim;
        SparseBlas::Dbsrmm(0,           // no transpose
                           m_blkRows,
                           1,           // multiplying to a vector
                           m_blkCols,
                           1.0,         // alpha
                           NistSpBlasDescra[m_matType],
                           &m_val[0], // nonzero entries
                           (int*)&m_indx[0], // column indices
                           (int*)&m_pntr[0],
                           (int*)&m_pntr[0]+1,
                           m_blkDim,
                           in,
                           m,
                           0.0,         // beta
                           out,
                           m,
                           NULL,
                           0);
    }


    /// General non-symmetric BSR multiply using light NIST Sparse BLAS API
    template<typename DataType>
    void StorageNistBsr<DataType>::MultiplyLight(
            const DataVectorType &in,
                  DataVectorType &out)
    {
        ASSERTL0(false, "StorageNistBsr::MultiplyLight() routine hasn't been tested to work");

        //using light API of NIST Sparse Blas
        SparseBlas::BSR_VecMult_CAB_double(
                            m_blkRows,
                            m_blkCols,
                            &m_val[0], // nonzero entries
                            (int*)&m_indx[0], // column indices
                            (int*)&m_pntr[0],
                            (int*)&m_pntr[0]+1,
                            m_blkDim,
                            &in[0],
                            &out[0],
                            0);
    }


    /// Symmetric BSR multiply using light NIST Sparse BLAS API
    template<typename DataType>
    void StorageNistBsr<DataType>::MultiplyLightSymm(
                const DataVectorType &in,
                      DataVectorType &out)
    {
        ASSERTL0(false, "StorageNistBsr::MultiplyLightSymm() routine hasn't been tested to work");

        SparseBlas::BSRsymm_VecMult_CAB_double(
                            m_blkRows,
                            m_blkCols,
                            &m_val[0], // nonzero entries
                            (int*)&m_indx[0], // column indices
                            (int*)&m_pntr[0],
                            (int*)&m_pntr[0]+1,
                            m_blkDim,
                            &in[0],
                            &out[0],
                            0);
    }

    // converts input vector of BCO blocks
    // to the internal Nist BSR representation
    template<typename DataType>
    void StorageNistBsr<DataType>::processBcoInput(
                        const IndexType  blkRows,
                        const IndexType  blkColumns,
                        const IndexType  blkDim,
                        const BCOMatType&   bcoMat)
    {
        IndexType i;
        BCOMatTypeConstIt  entry;
        IndexType rowcoord;
        IndexType colcoord;
        BCOEntryType       value;

        std::vector<IndexType> tmp(blkRows+1,0);

        // calculate the number of entries on each row
        // and store the result in tmp
        for(entry = bcoMat.begin(); entry != bcoMat.end(); entry++)
        {
            rowcoord = (entry->first).first;
            tmp[rowcoord]++;

            colcoord = (entry->first).second;
            value    =  entry->second;
        }
        // Based upon this information, fill the array m_pntr
        // which basically contains the offset of each row's 
        // first entry in the other arrays m_val and m_indx
        m_pntr[0] = 1;
        for(i = 0; i < blkRows; i++)
        {
            m_pntr[i+1] = m_pntr[i] + tmp[i]; 
        }

        // Copy the values of m_pntr into tmp as this will be needed 
        // in the following step
        std::copy(&m_pntr[0],&m_pntr[0]+blkRows+1,&tmp[0]);

        // Now, fill in index and value entries.
        for(entry = bcoMat.begin(); entry != bcoMat.end(); entry++)
        {
            rowcoord = (entry->first).first;
            colcoord = (entry->first).second;
            value    =  entry->second;
            int blkSize = blkDim*blkDim;

            for (i = 0; i < blkSize; i++)
            {
                m_val [ blkSize*(tmp[rowcoord]-1) + i ] = value[i];
                if (std::abs(value[i]) > NekConstants::kNekSparseNonZeroTol) m_nnz++;
            }

            m_indx[ tmp[rowcoord] ] = colcoord+1;
            tmp[rowcoord]++;
        }
    }


    // explicit instantiation
    template class StorageNistBsr<NekDouble>;


} // namespace

