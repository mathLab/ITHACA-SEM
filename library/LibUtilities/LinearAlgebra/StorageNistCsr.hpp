///////////////////////////////////////////////////////////////////////////////
//
// File: StorageNistCsr.hpp
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

///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_STORAGE_NIST_CSR_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_STORAGE_NIST_CSR_HPP

#include <map>
#include <vector>
#include <utility>
#include <fstream>

#include <LibUtilities/LinearAlgebra/MatrixStorageType.h>
#include <LibUtilities/LinearAlgebra/SparseMatrixFwd.hpp>

#include <boost/call_traits.hpp>

namespace Nektar
{

    /*
     *  This class is an interface to NIST CSR sparse matrix storage.
     *  The CSR (Compressed Sparse Row) sparse format assumes sparse
     *  matrix is a collection of individual non-zero entries with
     *  additional indexing: one index array stores column indices of
     *  non-zero entries row by row; another one stores offsets of each
     *  row in the first indexing array. Indexing is 0-based.
     *
     *  The constructor takes input matrix in coordinate storage
     *  (COO) sparse format.
     *
     *  Multiply kernels are provided by NIST Sparse BLAS v 0.9
     *
     */

    template<typename T>
    class StorageNistCsr
    {

    public:
        typedef T                             DataType;
        typedef Array<OneD, DataType>         DataVectorType;
        typedef Array<OneD, const DataType>   ConstDataVectorType;
        typedef Array<OneD, IndexType>        IndexVectorType;

        /// \internal
        /// \brief Forward iterator through nonzero elements of the matrix
        ///        that mimics forward iteration of COOMatType.
        ///        It's a dirty hack, not a real iterator in the C++ sense.
        class const_iterator
        {
            struct IterType
            {
                CoordType first;      //< (row, column)
                DataType  second;     //< value
                IndexType nnzindex;   //< index of this nnz entry
            };

            public:
                const_iterator( MatrixStorage       matType,
                          IndexType           begin,
                          IndexType           end,
                          const DataVectorType&     val,
                          const IndexVectorType&    indx,
                          const IndexVectorType&    pntr);
                const_iterator(const const_iterator& src);
                ~const_iterator();

                const_iterator  operator++(int);
                const_iterator& operator++();
                const IterType& operator*();
                const IterType* operator->();
                const bool operator==(const const_iterator& rhs);
                const bool operator!=(const const_iterator& rhs);

                // one way conversion: iterator -> const_iterator
                // operator const_iterator<T const, Tag>() const;

            private:
                void forward();

                MatrixStorage               m_matType;
                IterType                    m_iter;
                IndexType                   m_begin;
                IndexType                   m_end;
                const DataVectorType&       m_val;
                const IndexVectorType&      m_indx;
                const IndexVectorType&      m_pntr;
        };


    public:
        // Construct a CSR sparse matrix based on input COO storage
        LIB_UTILITIES_EXPORT StorageNistCsr( const IndexType rows,
                        const IndexType columns,
                        const COOMatType&  cooMat,
                        const MatrixStorage matType = eFULL);

        // Copy constructor
        LIB_UTILITIES_EXPORT StorageNistCsr(const StorageNistCsr& src);

        LIB_UTILITIES_EXPORT ~StorageNistCsr();


        LIB_UTILITIES_EXPORT const IndexType GetRows() const;
        LIB_UTILITIES_EXPORT const IndexType GetColumns() const;
        // NumNonZeroEntries == NumStoredDoubles for CSR, but
        // not for all sparse storages.
        LIB_UTILITIES_EXPORT const IndexType GetNumNonZeroEntries() const;
        LIB_UTILITIES_EXPORT const IndexType GetNumStoredDoubles() const;
        LIB_UTILITIES_EXPORT const IndexType GetBlkSize() const;
        LIB_UTILITIES_EXPORT const DataType  GetFillInRatio() const;
        LIB_UTILITIES_EXPORT const size_t GetMemoryUsage(IndexType nnz, IndexType nRows) const;

        LIB_UTILITIES_EXPORT const_iterator begin() const;
        LIB_UTILITIES_EXPORT const_iterator end() const;

        LIB_UTILITIES_EXPORT const typename boost::call_traits<DataType>::const_reference
                GetValue(IndexType row, IndexType column) const;

        LIB_UTILITIES_EXPORT void Multiply(const DataType* in,
                            DataType* out);
        LIB_UTILITIES_EXPORT void Multiply(const DataVectorType &in,
                            DataVectorType &out);
        LIB_UTILITIES_EXPORT void MultiplyLight(const DataVectorType &in,
                                 DataVectorType &out);
        LIB_UTILITIES_EXPORT void MultiplyLightSymm(const DataVectorType &in,
                                     DataVectorType &out);


    protected:

        // converts input vector of COO blocks
        // to the internal Nist CSR representation
        void processCooInput(
                        const IndexType  nRows,
                        const IndexType  nColumns,
                        const COOMatType&  cooMat);


        MatrixStorage   m_matType;
        unsigned long   m_mulCallsCounter;

        IndexType       m_rows;
        IndexType       m_columns;
        IndexType       m_nnz; //< number of nonzero entries in the sparse matrix

        DataVectorType  m_val;  //< values of non-zero entries
        IndexVectorType m_indx; //< column indices of non-zero entries
        IndexVectorType m_pntr; //< m_pntr(i) contains index in m_val of first non-zero element in row i

    private:

    };

} // namespace

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_STORAGE_NIST_CSR_HPP
