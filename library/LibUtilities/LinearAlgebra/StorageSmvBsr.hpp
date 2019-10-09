///////////////////////////////////////////////////////////////////////////////
//
// File: StorageSmvBsr.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_STORAGE_SMV_BSR_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_STORAGE_SMV_BSR_HPP

#include <map>
#include <vector>
#include <utility>
#include <fstream>

#include <LibUtilities/LinearAlgebra/MatrixStorageType.h>
#include <LibUtilities/LinearAlgebra/SparseMatrixFwd.hpp>

namespace Nektar
{
    /*
     *  Zero-based BSR (Block Sparse Row) storage class with its sparse
     *  multiply kernels built upon its own dense unrolled multiply kernels
     *  up to 4x4 matrices. When matrix is larger than or 4x4, the
     *  multiply kernel calls dense dgemv from BLAS.
     *
     *  The BSR sparse format assumes sparse matrix is a CSR collection of
     *  dense square blocks of same size. In contrast with Nist BSR class
     *  this one uses zero-based storage. The constructor takes input matrix in
     *  block coordinate storage (BCO) sparse format.
     *
     */

    template<typename T>
    class StorageSmvBsr
    {

    public:
        typedef T                             DataType;
        typedef Array<OneD, DataType>         DataVectorType;
        typedef Array<OneD, const DataType>   ConstDataVectorType;
        typedef Array<OneD, IndexType>        IndexVectorType;

        typedef void (*MultiplyKernel)(const double*, const double*, double*);

        /// \internal
        /// \brief Forward iterator through nonzero (double) elements of the matrix
        ///        that mimics forward iteration of BCOMatType.
        ///        It's a dirty hack, not a real iterator in the C++ sense.
        class const_iterator
        {
            struct IterType
            {
                CoordType  first;       //< (row, column)
                DataType   second;      //< value
                IndexType  nnzindex;    //< index of this nnz entry
                IndexType  storageindex;//< offset of this nnz entry in the storage
            };

            public:
                const_iterator( MatrixStorage       matType,
                          IndexType           begin,
                          IndexType           end,
                          IndexType           blkDim,
                          const DataVectorType&     val,
                          const IndexVectorType&    indx,
                          const IndexVectorType&    pntr);
                const_iterator(const const_iterator& src);
                ~const_iterator();

                const_iterator operator++(int);
                const_iterator& operator++();
                const IterType& operator*();
                const IterType* operator->();
                bool operator==(const const_iterator& rhs);
                bool operator!=(const const_iterator& rhs);

                // one way conversion: iterator -> const_iterator
                // operator const_iterator<T const, Tag>() const;

            private:
                void forward();
                CoordType storageIndexToFullCoord(IndexType storageIndex);

                MatrixStorage               m_matType;
                IterType                    m_iter;
                IndexType                   m_begin;
                IndexType                   m_end;
                IndexType                   m_blkDim;
                const DataVectorType&       m_val;
                const IndexVectorType&      m_indx;
                const IndexVectorType&      m_pntr;
        };


    public:
        // Constructs zero-based BSR sparse matrix based on input BCO storage
        LIB_UTILITIES_EXPORT StorageSmvBsr( const IndexType  blkRows,
                             const IndexType  blkCols,
                             const IndexType  blkDim,
                             const BCOMatType&   bcoMat,
                             const MatrixStorage matType = eFULL);

        // Copy constructor
        LIB_UTILITIES_EXPORT StorageSmvBsr(const StorageSmvBsr& src);

        LIB_UTILITIES_EXPORT ~StorageSmvBsr();

        LIB_UTILITIES_EXPORT IndexType GetRows() const;
        LIB_UTILITIES_EXPORT IndexType GetColumns() const;
        LIB_UTILITIES_EXPORT IndexType GetNumNonZeroEntries() const;
        LIB_UTILITIES_EXPORT IndexType GetNumStoredDoubles() const;
        LIB_UTILITIES_EXPORT IndexType GetBlkSize() const;
        LIB_UTILITIES_EXPORT DataType  GetFillInRatio() const;
        LIB_UTILITIES_EXPORT size_t GetMemoryUsage() const;

        LIB_UTILITIES_EXPORT const_iterator begin() const;
        LIB_UTILITIES_EXPORT const_iterator end() const;

        LIB_UTILITIES_EXPORT const DataType &GetValue(
            IndexType row, IndexType column) const;

        LIB_UTILITIES_EXPORT void Multiply(const DataType* in,
                            DataType* out);
        LIB_UTILITIES_EXPORT void Multiply(const DataVectorType &in,
                            DataVectorType &out);
        LIB_UTILITIES_EXPORT void MultiplyLight(const DataVectorType &in,
                                 DataVectorType &out);


    protected:

        // converts input vector of BCO blocks
        // to the internal zero-based BSR representation
        void processBcoInput(
                        const IndexType  blkRows,
                        const IndexType  blkDim,
                        const BCOMatType&   bcoMat);


        void Multiply_1x1(const int mb, const double* val,
                    const int* bindx, const int* bpntrb, const int* bpntre,
                    const double* b, double* c);

        void Multiply_2x2(const int mb, const double* val,
                    const int* bindx, const int* bpntrb, const int* bpntre,
                    const double* b, double* c);

        void Multiply_3x3(const int mb, const double* val,
                    const int* bindx, const int* bpntrb, const int* bpntre,
                    const double* b, double* c);

        void Multiply_4x4(const int mb, const double* val,
                    const int* bindx, const int* bpntrb, const int* bpntre,
                    const double* b, double* c);

        void Multiply_generic(const int mb, const double* val,
                    const int* bindx, const int* bpntrb, const int* bpntre,
                    const double* b, double* c);

        // interface to lowest level LibSMV multiply kernels
        MultiplyKernel   m_mvKernel;

        MatrixStorage    m_matType;

        IndexType        m_blkRows;  // number of block rows
        IndexType        m_blkCols;  // number of block columns
        IndexType        m_blkDim;   // rank of each block

        IndexType        m_bnnz;  //< number of nonzero blocks
        IndexType        m_nnz;   //< number of factual nonzero entries in the sparse matrix

        DataVectorType   m_val;  // values of non-zero entries
        IndexVectorType  m_indx; // column indices of non-zero entries
        IndexVectorType  m_pntr; // m_pntr(i) contains index in m_val of first non-zero element in row i

    private:

    };



} // namespace

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_STORAGE_SMV_BSR_HPP
