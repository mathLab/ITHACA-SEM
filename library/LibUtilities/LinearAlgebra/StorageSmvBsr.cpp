///////////////////////////////////////////////////////////////////////////////
//
// File: StorageSmvBsr.cpp
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
// Description: 0-based sparse BSR storage class with own unrolled multiply
//              kernels.
//
///////////////////////////////////////////////////////////////////////////////

#include <map>
#include <vector>
#include <utility>
#include <fstream>
#include <cmath>

#include <LibUtilities/BasicConst/NektarUnivConsts.hpp>
#include <LibUtilities/LinearAlgebra/Blas.hpp>
#include <LibUtilities/LinearAlgebra/StorageSmvBsr.hpp>
#include <LibUtilities/LinearAlgebra/NistSparseDescriptors.hpp>

namespace Nektar
{

    template<typename DataType>
    StorageSmvBsr<DataType>::const_iterator::const_iterator(
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
                if (m_iter.storageindex >= m_val.size())
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
    CoordType StorageSmvBsr<DataType>::const_iterator::storageIndexToFullCoord(IndexType storageIndex)
    {
        // find local row and column indices within this block

        const IndexType elms       = m_blkDim*m_blkDim;
        const IndexType block      = storageIndex / elms;
        const IndexType loc_offset = storageIndex % elms;
        const IndexType loc_row    = loc_offset % m_blkDim;
        const IndexType loc_col    = loc_offset / m_blkDim;

        // find block row and block column

        IndexType block_col = m_indx[block];
        IndexType block_row = 0;
        while(block >= m_pntr[block_row+1])
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
    StorageSmvBsr<DataType>::const_iterator::const_iterator(const const_iterator& src):
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
    StorageSmvBsr<DataType>::const_iterator::~const_iterator()
    {
    }

    template<typename DataType>
    typename StorageSmvBsr<DataType>::const_iterator StorageSmvBsr<DataType>::const_iterator::operator++(int)
    {
        const_iterator out = *this;
        forward();
        return out;
    }

    template<typename DataType>
    typename StorageSmvBsr<DataType>::const_iterator& StorageSmvBsr<DataType>::const_iterator::operator++()
    {
        forward();
        return *this;
    }

    template<typename DataType>
    const typename StorageSmvBsr<DataType>::const_iterator::IterType& StorageSmvBsr<DataType>::const_iterator::operator*()
    {
        return m_iter;
    }

    template<typename DataType>
    const typename StorageSmvBsr<DataType>::const_iterator::IterType* StorageSmvBsr<DataType>::const_iterator::operator->()
    {
        return &m_iter;
    }

    template<typename DataType>
    bool StorageSmvBsr<DataType>::const_iterator::operator==(const const_iterator& rhs)
    {
        return m_iter.nnzindex == rhs.m_iter.nnzindex;
    }

    template<typename DataType>
    bool StorageSmvBsr<DataType>::const_iterator::operator!=(const const_iterator& rhs)
    {
        return !(m_iter.nnzindex == rhs.m_iter.nnzindex);
    }

    template<typename DataType>
    void StorageSmvBsr<DataType>::const_iterator::forward()
    {
        while((m_iter.storageindex+1 < m_val.size()) &&
              (m_val[++m_iter.storageindex] <= NekConstants::kNekSparseNonZeroTol));

        m_iter.nnzindex++;

        if (m_iter.storageindex >= m_val.size())
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
    StorageSmvBsr<DataType>::StorageSmvBsr(
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

        processBcoInput(blkRows,blkDim,bcoMat);
    }


    template<typename DataType>
    StorageSmvBsr<DataType>::StorageSmvBsr(const StorageSmvBsr& src):
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
    StorageSmvBsr<DataType>::~StorageSmvBsr()
    {
    }


    template<typename DataType>
    IndexType StorageSmvBsr<DataType>::GetRows() const
    {
        return m_blkRows*m_blkDim;
    }

    template<typename DataType>
    IndexType StorageSmvBsr<DataType>::GetColumns() const
    {
        return m_blkCols*m_blkDim;
    }

    template<typename DataType>
    IndexType StorageSmvBsr<DataType>::GetNumNonZeroEntries() const
    {
        return m_nnz;
    }

    template<typename DataType>
    IndexType StorageSmvBsr<DataType>::GetBlkSize() const
    {
        return m_blkDim;
    }


    template<typename DataType>
    IndexType StorageSmvBsr<DataType>::GetNumStoredDoubles() const
    {
        return m_bnnz*m_blkDim*m_blkDim;
    }

    template<typename DataType>
    DataType StorageSmvBsr<DataType>::GetFillInRatio() const
    {
        return (DataType)(m_bnnz*m_blkDim*m_blkDim)/(DataType)m_nnz;
    }


    template<typename DataType>
    size_t StorageSmvBsr<DataType>::GetMemoryUsage() const
    {
        return sizeof(DataType) *m_val.capacity()   +
               sizeof(IndexType)*m_indx.capacity() +
               sizeof(IndexType)*m_pntr.capacity()  +
               sizeof(IndexType)*5   + //< blkRows + blkCols + blkDim + nnz + bnnz
               sizeof(MatrixStorage);
    }


    template<typename DataType>
    const DataType &StorageSmvBsr<DataType>::GetValue(IndexType grow, IndexType gcolumn) const
    {
        IndexType  brow = grow    / m_blkDim;
        IndexType  bcol = gcolumn / m_blkDim;
        IndexType  lrow = grow    % m_blkDim;
        IndexType  lcol = gcolumn % m_blkDim;

        // rewind to the first entry of the first
        // block in the current block row
        IndexType  offset = m_pntr[brow]*m_blkDim*m_blkDim;

        IndexType i;
        static DataType defaultReturnValue;
        for( i = m_pntr[brow]; i < m_pntr[brow+1]; i++)
        {
            if(bcol == m_indx[i])
            {
                return m_val[offset+lrow + lcol*m_blkDim];
            }
            offset += m_blkDim*m_blkDim;
        }

        return defaultReturnValue;
    }


    template<typename DataType>
    typename StorageSmvBsr<DataType>::const_iterator StorageSmvBsr<DataType>::begin() const
    {
        return const_iterator(m_matType, 0, m_nnz, m_blkDim, m_val, m_indx, m_pntr);
    }

    template<typename DataType>
    typename StorageSmvBsr<DataType>::const_iterator StorageSmvBsr<DataType>::end() const
    {
        return const_iterator(m_matType, m_nnz, m_nnz, m_blkDim, m_val, m_indx, m_pntr);
    }



    // General non-symmetric zero-based BSR multiply
    // C = A*B where A is matrix and B is vector.
    // No scaling. Previous content of C is discarded.
    template<typename DataType>
    void StorageSmvBsr<DataType>::Multiply(
            const DataVectorType &in,
                  DataVectorType &out)
    {
        const double* b = &in[0];
              double* c = &out[0];
        const double* val = &m_val[0];
        const int* bindx  = (int*)&m_indx[0];
        const int* bpntrb = (int*)&m_pntr[0];
        const int* bpntre = (int*)&m_pntr[0]+1;
        const int  mb = m_blkRows;

        switch(m_blkDim)
        {
        case 1:  Multiply_1x1(mb,val,bindx,bpntrb,bpntre,b,c); return;
        case 2:  Multiply_2x2(mb,val,bindx,bpntrb,bpntre,b,c); return;
        default:
            {
                Multiply_generic(mb,val,bindx,bpntrb,bpntre,b,c); return;
            }
        }
    }


    template<typename DataType>
    void StorageSmvBsr<DataType>::Multiply(
            const DataType*  in,
                  DataType*  out)
    {
        const double* b = &in[0];
              double* c = &out[0];
        const double* val = &m_val[0];
        const int* bindx  = (int*)&m_indx[0];
        const int* bpntrb = (int*)&m_pntr[0];
        const int* bpntre = (int*)&m_pntr[0]+1;
        const int  mb = m_blkRows;

        switch(m_blkDim)
        {
        case 1:  Multiply_1x1(mb,val,bindx,bpntrb,bpntre,b,c); return;
        case 2:  Multiply_2x2(mb,val,bindx,bpntrb,bpntre,b,c); return;
        default:
            {
                Multiply_generic(mb,val,bindx,bpntrb,bpntre,b,c); return;
            }
        }
    }


    // This function is defined for backward compatibility with
    // NIST storage classes
    template<typename DataType>
    void StorageSmvBsr<DataType>::MultiplyLight(
            const DataVectorType &in,
                  DataVectorType &out)
    {
        const double* b = &in[0];
              double* c = &out[0];
        const double* val = &m_val[0];
        const int* bindx  = (int*)&m_indx[0];
        const int* bpntrb = (int*)&m_pntr[0];
        const int* bpntre = (int*)&m_pntr[0]+1;
        const int  mb = m_blkRows;

        switch(m_blkDim)
        {
        case 1:  Multiply_1x1(mb,val,bindx,bpntrb,bpntre,b,c); return;
        case 2:  Multiply_2x2(mb,val,bindx,bpntrb,bpntre,b,c); return;
        default:
            {
                Multiply_generic(mb,val,bindx,bpntrb,bpntre,b,c); return;
            }
        }
    }

    /// Zero-based CSR multiply.
    /// Essentially this is slightly modified copy-paste from
    /// NIST Sparse Blas 0.9b routine CSR_VecMult_CAB_double()
    template<typename DataType>
    void StorageSmvBsr<DataType>::Multiply_1x1(
            const int mb,
            const double* val,
            const int* bindx,
            const int* bpntrb,
            const int* bpntre,
            const double* b,
                  double* c)
    {
        for (int i=0;i!=mb;i++)
        {
            double t = 0;
            int jb = bpntrb[i];
            int je = bpntre[i];
            for (int j=jb;j!=je;j++)
            {
                t += b[bindx[j]] * (*val++);
            }
            c[i] = t;
        }
    }

    /// Zero-based BSR multiply unrolled for 2x2 blocks.
    /// Essentially this is slightly optimised copy-paste from
    /// NIST Sparse Blas 0.9b routine BSR_VecMult_BAB_double()
    template<typename DataType>
    void StorageSmvBsr<DataType>::Multiply_2x2(
            const int mb,
            const double* val,
            const int* bindx,
            const int* bpntrb,
            const int* bpntre,
            const double* b,
                  double* c)
    {
        const int lb = 2;

        const double *pval = val;
        double *pc=c;

        for (int i=0;i!=mb;i++)
        {
            int jb = bpntrb[i];
            int je = bpntre[i];
            pc[0] = 0;
            pc[1] = 0;
            for (int j=jb;j!=je;j++)
            {
                int bs=bindx[j]*lb;
                const double *pb = &b[bs];

                pc[0] += pb[0] * pval[0] + pb[1] * pval[2];
                pc[1] += pb[0] * pval[1] + pb[1] * pval[3];
                pval += 4;
            }
            pc += 2;
        }
    }

    /// Zero-based BSR multiply unrolled for 3x3 blocks
    template<typename DataType>
    void StorageSmvBsr<DataType>::Multiply_3x3(
            const int mb,
            const double* val,
            const int* bindx,
            const int* bpntrb,
            const int* bpntre,
            const double* b,
                  double* c)
    {
        const int lb = 3;

        const double *pval = val;
        double *pc=c;

        for (int i=0;i!=mb;i++)
        {
            int jb = bpntrb[i];
            int je = bpntre[i];
            pc[0] = 0;
            pc[1] = 0;
            pc[2] = 0;
            for (int j=jb;j!=je;j++)
            {
                int bs=bindx[j]*lb;
                const double *pb = &b[bs];

                pc[0] += pb[0] * pval[0] + pb[1] * pval[3] + pb[2] * pval[6];
                pc[1] += pb[0] * pval[1] + pb[1] * pval[4] + pb[2] * pval[7];
                pc[2] += pb[0] * pval[2] + pb[1] * pval[5] + pb[2] * pval[8];
                pval += 9;
            }
            pc += 3;
        }
    }

    /// Zero-based BSR multiply unrolled for 4x4 blocks
    template<typename DataType>
    void StorageSmvBsr<DataType>::Multiply_4x4(
            const int mb,
            const double* val,
            const int* bindx,
            const int* bpntrb,
            const int* bpntre,
            const double* b,
                  double* c)
    {
        const int lb = 4;

        const double *pval = val;
        double *pc=c;

        for (int i=0;i!=mb;i++)
        {
            int jb = bpntrb[i];
            int je = bpntre[i];
            pc[0] = 0;
            pc[1] = 0;
            pc[2] = 0;
            pc[3] = 0;
            for (int j=jb;j!=je;j++)
            {
                int bs=bindx[j]*lb;
                const double *pb = &b[bs];

                pc[0] += pb[0] * pval[0] + pb[1] * pval[4] + pb[2] * pval[ 8] + pb[3] * pval[12];
                pc[1] += pb[0] * pval[1] + pb[1] * pval[5] + pb[2] * pval[ 9] + pb[3] * pval[13];
                pc[2] += pb[0] * pval[2] + pb[1] * pval[6] + pb[2] * pval[10] + pb[3] * pval[14];
                pc[3] += pb[0] * pval[3] + pb[1] * pval[7] + pb[2] * pval[11] + pb[3] * pval[15];
                pval += 16;
            }
            pc += 4;
        }
    }

    /// Generic zero-based BSR multiply for higher matrix ranks
    template<typename DataType>
    void StorageSmvBsr<DataType>::Multiply_generic(
            const int mb,
            const double* val,
            const int* bindx,
            const int* bpntrb,
            const int* bpntre,
            const double* b,
                  double* c)
    {
        const int lb = m_blkDim;
        const double *pval = val;
        const int mm=lb*lb;
        double *pc=c;
        for (int i=0;i!=mb*lb;i++) *pc++ = 0;

        pc=c;
        for (int i=0;i!=mb;i++)
        {
            int jb = bpntrb[i];
            int je = bpntre[i];
            for (int j=jb;j!=je;j++)
            {
                Blas::Dgemv('N',lb,lb,1.0,pval,lb,&b[bindx[j]*lb],1,1.0,pc,1);
                pval+=mm;
            }
            pc += lb;
        }
    }



    // converts input vector of BCO blocks
    // to the internal BSR representation
    template<typename DataType>
    void StorageSmvBsr<DataType>::processBcoInput(
                        const IndexType  blkRows,
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
        }
        // Based upon this information, fill the array m_pntr
        // which basically contains the offset of each row's
        // first entry in the other arrays m_val and m_indx
        m_pntr[0] = 0;
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
                m_val [ blkSize*(tmp[rowcoord]) + i ] = value[i];
                if (std::abs(value[i]) > NekConstants::kNekSparseNonZeroTol) m_nnz++;
            }

            m_indx[ tmp[rowcoord] ] = colcoord;
            tmp[rowcoord]++;
        }
    }


    // explicit instantiation
    template class StorageSmvBsr<NekDouble>;


} // namespace

