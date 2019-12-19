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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SPARSE_DIAG_BLK_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SPARSE_DIAG_BLK_MATRIX_HPP

#include <map>
#include <vector>
#include <utility>
#include <algorithm>
#include <fstream>

#include <LibUtilities/LinearAlgebra/SparseMatrixFwd.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <boost/call_traits.hpp>


namespace Nektar
{

    /*
     * This is a class-container to diagonal block matrix
     * with elements being sparse matrices. The type of
     * sparse entries is defined with template parameter.
     *
     */
    template<typename SparseStorageType>
    class NekSparseDiagBlkMatrix
    {
    public:

        typedef SparseStorageType                    StorageType;
        typedef typename SparseStorageType::DataType DataType;
        typedef std::shared_ptr<SparseStorageType>   SparseStorageSharedPtr;
        typedef Array<OneD, DataType>                DataVectorType;
        typedef Array<OneD, const DataType>          ConstDataVectorType;
        typedef Array<OneD, SparseStorageSharedPtr>  SparseStorageSharedPtrVector;


        LIB_UTILITIES_EXPORT NekSparseDiagBlkMatrix(const SparseStorageSharedPtrVector& sparseStoragePtrVector);
        LIB_UTILITIES_EXPORT NekSparseDiagBlkMatrix(const NekSparseDiagBlkMatrix& src);
        LIB_UTILITIES_EXPORT ~NekSparseDiagBlkMatrix();

        LIB_UTILITIES_EXPORT IndexType GetRows() const;
        LIB_UTILITIES_EXPORT IndexType GetColumns() const;
        LIB_UTILITIES_EXPORT IndexType GetNumNonZeroEntries();
        LIB_UTILITIES_EXPORT DataType  GetFillInRatio() const;

        LIB_UTILITIES_EXPORT IndexType GetRows(int i) const;
        LIB_UTILITIES_EXPORT IndexType GetColumns(int i) const;
        LIB_UTILITIES_EXPORT IndexType GetNumberOfMatrixBlocks() const;
        LIB_UTILITIES_EXPORT IndexType GetNumNonZeroEntries(int i) const;
        LIB_UTILITIES_EXPORT DataType  GetFillInRatio(int i) const;


        LIB_UTILITIES_EXPORT typename boost::call_traits<DataType>::const_reference
                 operator()(const IndexType row, const IndexType column) const;
        LIB_UTILITIES_EXPORT typename boost::call_traits<DataType>::const_reference
                 operator()(const IndexType block, const IndexType row, const IndexType column) const;

        //typename SparseStorageType::const_iterator begin() const;
        //typename SparseStorageType::const_iterator end() const;

        LIB_UTILITIES_EXPORT void Multiply(const DataVectorType &in,
                            DataVectorType &out);
        LIB_UTILITIES_EXPORT void Multiply(const DataType* in,
                            DataType* out);
        LIB_UTILITIES_EXPORT void MultiplySubMatrix( const IndexType blockNum,
                                DataType* in,
                                DataType* out);

        LIB_UTILITIES_EXPORT size_t GetMemoryFootprint();
        LIB_UTILITIES_EXPORT size_t GetMemoryFootprint(IndexType i) const;

        LIB_UTILITIES_EXPORT unsigned long GetMulCallsCounter() const;
        LIB_UTILITIES_EXPORT DataType   GetAvgRowDensity();
        LIB_UTILITIES_EXPORT DataType   GetAvgRowDensity(IndexType i) const;
        LIB_UTILITIES_EXPORT IndexType  GetBandwidth();
        LIB_UTILITIES_EXPORT IndexType  GetBandwidth(IndexType i);
        LIB_UTILITIES_EXPORT COOMatTypeSharedPtr GetCooStorage();
        LIB_UTILITIES_EXPORT COOMatTypeSharedPtr GetCooStorage(IndexType i);

        LIB_UTILITIES_EXPORT void writeSparsityPatternTo(std::ostream& out, IndexType blockSize = 64);
        LIB_UTILITIES_EXPORT void writeSubmatrixSparsityPatternTo(std::ostream& out, const IndexType subMatrixIdx, IndexType blockSize = 64);
        LIB_UTILITIES_EXPORT void writeBlockSparsityPatternTo(std::ostream& out,
                        const IndexType blk_row = 0, const IndexType blk_col = 0, IndexType blockSize = 64);
        LIB_UTILITIES_EXPORT void writeSubmatrixBlockSparsityPatternTo(std::ostream& out, const IndexType subMatrixIdx,
                        const IndexType blk_row = 0, const IndexType blk_col = 0, IndexType blockSize = 64);

    protected:

        IndexType      m_rows;
        IndexType      m_cols;
        IndexVector    m_rowoffset;
        unsigned long  m_mulCallsCounter;
        SparseStorageSharedPtrVector m_submatrix;

    private:

    };

}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SPARSE_DIAG_BLK_MATRIX_HPP
