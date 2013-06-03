///////////////////////////////////////////////////////////////////////////////
//
// File: StorageNistCsr.cpp
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
// Description: Interface to NIST CSR sparse matrix storage.
//
///////////////////////////////////////////////////////////////////////////////

#include <map>
#include <vector>
#include <utility>
#include <fstream>

#include <LibUtilities/LinearAlgebra/NistSparseBlas.hpp>
#include <LibUtilities/LinearAlgebra/StorageNistCsr.hpp>
#include <LibUtilities/LinearAlgebra/NistSparseDescriptors.hpp>

#include <boost/lexical_cast.hpp>

namespace Nektar
{

    template<typename DataType>
    StorageNistCsr<DataType>::const_iterator::const_iterator(
                    MatrixStorage       matType,
                    IndexType           begin,
                    IndexType           end,
                    const DataVectorType&     val,
                    const IndexVectorType&    indx,
                    const IndexVectorType&    pntr) :
        m_matType(matType),
        m_iter(),
        m_begin(begin),
        m_end(end),
        m_val(val),
        m_indx(indx),
        m_pntr(pntr)
    {
        m_iter.nnzindex = begin;
        if (begin < end)
        {
            IndexType& row  = m_iter.first.first;
            IndexType& col  = m_iter.first.second;
            m_iter.second   = m_val[begin];
            col             = m_indx[begin];
            row             = 0;
            while(begin >= m_pntr[row+1])
            {
                row++;
            }
        }
    }

    template<typename DataType>
    StorageNistCsr<DataType>::const_iterator::const_iterator(const const_iterator& src):
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
    StorageNistCsr<DataType>::const_iterator::~const_iterator()
    {
    }

    template<typename DataType>
    typename StorageNistCsr<DataType>::const_iterator StorageNistCsr<DataType>::const_iterator::operator++(int)
    {
        const_iterator out = *this;
        forward();
        return out;
    }

    template<typename DataType>
    typename StorageNistCsr<DataType>::const_iterator& StorageNistCsr<DataType>::const_iterator::operator++()
    {
        forward();
        return *this;
    }

    template<typename DataType>
    const typename StorageNistCsr<DataType>::const_iterator::IterType& StorageNistCsr<DataType>::const_iterator::operator*()
    {
        return m_iter;
    }

    template<typename DataType>
    const typename StorageNistCsr<DataType>::const_iterator::IterType* StorageNistCsr<DataType>::const_iterator::operator->()
    {
        return &m_iter;
    }

    template<typename DataType>
    const bool StorageNistCsr<DataType>::const_iterator::operator==(const const_iterator& rhs)
    {
        return m_iter.nnzindex == rhs.m_iter.nnzindex;
    }

    template<typename DataType>
    const bool StorageNistCsr<DataType>::const_iterator::operator!=(const const_iterator& rhs)
    {
        return !(m_iter.nnzindex == rhs.m_iter.nnzindex);
    }

    template<typename DataType>
    void StorageNistCsr<DataType>::const_iterator::forward()
    {
        IndexType& row = m_iter.first.first;
        IndexType& col = m_iter.first.second;

        m_iter.nnzindex++;

        while((m_pntr.num_elements() > row+1) && (m_pntr[row+1] == m_iter.nnzindex))
        {
            row++;
        }
        if (m_indx.num_elements() > m_iter.nnzindex)
        {
            col = m_indx[m_iter.nnzindex];
            m_iter.second = m_val[m_iter.nnzindex];
        }
    }


    template<typename DataType>
    StorageNistCsr<DataType>::StorageNistCsr(const StorageNistCsr& src):
        m_matType(src.m_matType),
        m_rows(src.m_rows),
        m_columns(src.m_columns),
        m_nnz(src.m_nnz),
        m_val(src.m_val),
        m_indx(src.m_indx),
        m_pntr(src.m_pntr)
    {
    }

    template<typename DataType>
    StorageNistCsr<DataType>::~StorageNistCsr()
    {
    }


    template<typename DataType>
    StorageNistCsr<DataType>::StorageNistCsr(
                    const IndexType rows,
                    const IndexType columns,
                    const COOMatType&  cooMat,
                    const MatrixStorage matType):
        m_matType(matType),
        m_rows(rows),
        m_columns(columns),
        m_nnz(cooMat.size()),
        m_val(m_nnz),
        m_indx(m_nnz),
        m_pntr(rows+1)
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


        processCooInput(rows,columns,cooMat);

    }


    template<typename DataType>
    const IndexType StorageNistCsr<DataType>::GetRows() const
    {
        return m_rows;
    }

    template<typename DataType>
    const IndexType StorageNistCsr<DataType>::GetColumns() const
    {
        return m_columns;
    }

    template<typename DataType>
    const IndexType StorageNistCsr<DataType>::GetNumNonZeroEntries() const
    {
        return m_nnz;
    }

    template<typename DataType>
    const IndexType StorageNistCsr<DataType>::GetNumStoredDoubles() const
    {
        return m_nnz;
    }

    template<typename DataType>
    const IndexType StorageNistCsr<DataType>::GetBlkSize() const
    {
        return 1;
    }


    template<typename DataType>
    const DataType StorageNistCsr<DataType>::GetFillInRatio() const
    {
        return 1.0;
    }


    template<typename DataType>
    const size_t StorageNistCsr<DataType>::GetMemoryUsage(IndexType nnz, IndexType nRows) const
    {
        return sizeof(DataType) *m_val.capacity()   +
               sizeof(IndexType)*m_indx.capacity() +
               sizeof(IndexType)*m_pntr.capacity()  +
               sizeof(IndexType)*3   + //<  sizes + nnz
               sizeof(MatrixStorage);
    }


    template<typename DataType>
    const typename boost::call_traits<DataType>::const_reference StorageNistCsr<DataType>::GetValue(IndexType row, IndexType column) const
    {
/*
        ASSERTL1(row < GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                 std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(GetRows()) +
                 std::string(" rows"));
        ASSERTL1(column < GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                 std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(GetColumns()) +
                 std::string(" columns"));
*/
        IndexType i;
        static DataType defaultReturnValue;
        for( i = m_pntr[row]; i < m_pntr[row+1]; i++)
        {
            if(column == m_indx[i])
            {
                return m_val[i];
            }
        }

        return defaultReturnValue;
    }


    template<typename DataType>
    typename StorageNistCsr<DataType>::const_iterator StorageNistCsr<DataType>::begin() const
    {
        return const_iterator(m_matType, 0, m_nnz, m_val, m_indx, m_pntr);
    }

    template<typename DataType>
    typename StorageNistCsr<DataType>::const_iterator StorageNistCsr<DataType>::end() const
    {
        return const_iterator(m_matType, m_nnz, m_nnz, m_val, m_indx, m_pntr);
    }


    template<typename DataType>
    void StorageNistCsr<DataType>::Multiply(const DataVectorType &in,
                                                  DataVectorType &out)
    {
        const IndexType m = m_rows;
        SparseBlas::Dcsrmm(0,           // no transpose
                           m,
                           1,           // multiplying to a vector
                           m_columns,
                           1.0,         // alpha
                           NistSpBlasDescra[m_matType],
                           &m_val[0], // nonzero entries
                           (int*)&m_indx[0], // column indices
                           (int*)&m_pntr[0],
                           (int*)&m_pntr[0]+1,
                           &in[0],
                           m,
                           0.0,         // beta
                           &out[0],
                           m,
                           NULL,
                           0);
    }

    template<typename DataType>
    void StorageNistCsr<DataType>::Multiply(const DataType* in,
                                                  DataType* out)
    {
        const IndexType m = m_rows;
        SparseBlas::Dcsrmm(0,           // no transpose
                           m,
                           1,           // multiplying to a vector
                           m_columns,
                           1.0,         // alpha
                           NistSpBlasDescra[m_matType],
                           &m_val[0], // nonzero entries
                           (int*)&m_indx[0], // column indices
                           (int*)&m_pntr[0],
                           (int*)&m_pntr[0]+1,
                           in,
                           m,
                           0.0,         // beta
                           out,
                           m,
                           NULL,
                           0);
    }


    // General unsymmetric CSR multiply
    template<typename DataType>
    void StorageNistCsr<DataType>::MultiplyLight(
            const DataVectorType &in,
                  DataVectorType &out)
    {
        //using light API of NIST Sparse Blas
        SparseBlas::CSR_VecMult_CAB_double(
                            m_rows,
                            m_columns,
                            &m_val[0], // nonzero entries
                            (int*)&m_indx[0], // column indices
                            (int*)&m_pntr[0],
                            (int*)&m_pntr[0]+1,
                            &in[0],
                            &out[0],
                            0);
    }

    // Symmetric CSR multiply --- works correctly only
    // for lower-triangular stored matrix. For upper-triangular
    // need to use CSC* call.
    template<typename DataType>
    void StorageNistCsr<DataType>::MultiplyLightSymm(const DataVectorType &in,
                                      DataVectorType &out)
    {
        //using light API of NIST Sparse Blas
        SparseBlas::CSRsymm_VecMult_CAB_double(
                            m_rows,
                            m_columns,
                            &m_val[0], // nonzero entries
                            (int*)&m_indx[0], // column indices
                            (int*)&m_pntr[0],
                            (int*)&m_pntr[0]+1,
                            &in[0],
                            &out[0],
                            0);
    }


    // converts input vector of COO blocks
    // to the internal Nist CSR representation
    template<typename DataType>
    void StorageNistCsr<DataType>::processCooInput(
                    const IndexType  nRows,
                    const IndexType  nColumns,
                    const COOMatType&   cooMat)
    {
        IndexType i;
        COOMatTypeConstIt entry;
        int rowcoord;
        int colcoord;
        DataType     value;

        std::vector<IndexType> tmp(nRows+1,0);

        m_rows = nRows;
        m_columns = nColumns;

        // calculate the number of entries on each row
        // and store the result in tmp
        for(entry = cooMat.begin(); entry != cooMat.end(); entry++)
        {
            rowcoord = (entry->first).first;
            tmp[rowcoord]++;
        }
        // Based upon this information, fill the array m_pntr
        // which basically contains the offset of each row's 
        // first entry in the other arrays m_val and m_indx
        m_pntr[0] = 0;
        for(i = 0; i < nRows; i++)
        {
            m_pntr[i+1] = m_pntr[i] + tmp[i]; 
        }
        // Copy the values of m_pntr into tmp as this will be needed 
        // in the following step
        std::copy(&m_pntr[0],&m_pntr[0]+nRows+1,&tmp[0]);
        // Now, fill in index and value entries.
        for(entry = cooMat.begin(); entry != cooMat.end(); entry++)
        {
            rowcoord = (entry->first).first;
            colcoord = (entry->first).second;
            value    =  entry->second;

            m_val [ tmp[rowcoord] ] = value;
            m_indx[ tmp[rowcoord] ] = colcoord;
            tmp[rowcoord]++;
        }
    }



    // explicit instantiation
    template class StorageNistCsr<NekDouble>;


} // namespace


