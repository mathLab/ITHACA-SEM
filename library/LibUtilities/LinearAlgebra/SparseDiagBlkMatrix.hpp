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
// Description: generic sparse matrix class templated by underlying sparse
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
#include <LibUtilities/BasicUtils/SharedArray.hpp> ///< \todo: forward declare

#include <boost/call_traits.hpp>


namespace Nektar
{

    // Construct a CSR sparse matrix based on input matrix in
    // coordinate storage (COO) sparse format.
    // This COO sparse matrix is given as a map< pair<int,int>, NekDouble>
    // where the pair refers to the coordinate of the non-zero entry
    // and the NekDouble contains its value.
    // The constructor now converts from COO storage to CSR storage.
    //
    // Generic, symmetric, diagonal and lower-/upper-triangular 
    // input matrix type identifiers respected.
    // For symmetric matrix type input COO storage should define its
    // upper-triangular part. For all other matrix properties the
    // conversion from COO to CSR is straightforward.
    //
    template<typename SparseStorageType>
    class NekSparseDiagBlkMatrix
    {
    public:

        typedef SparseStorageType                    StorageType;
        typedef typename SparseStorageType::DataType DataType;
        typedef boost::shared_ptr<SparseStorageType> SparseStorageSharedPtr;
        typedef Array<OneD, DataType>                DataVectorType;
        typedef Array<OneD, const DataType>          ConstDataVectorType;
        typedef Array<OneD, SparseStorageSharedPtr>  SparseStorageSharedPtrVector;


        NekSparseDiagBlkMatrix(const SparseStorageSharedPtrVector& sparseStoragePtrVector);
        NekSparseDiagBlkMatrix(const NekSparseDiagBlkMatrix& src);
        ~NekSparseDiagBlkMatrix();

        const IndexType GetRows() const;
        const IndexType GetColumns() const;
        const IndexType GetNumNonZeroEntries();
        const DataType  GetFillInRatio() const;

        const IndexType GetRows(int i) const;
        const IndexType GetColumns(int i) const;
        const IndexType GetNumberOfMatrixBlocks() const;
        const IndexType GetNumNonZeroEntries(int i) const;
        const DataType  GetFillInRatio(int i) const;


        typename boost::call_traits<DataType>::const_reference
                 operator()(const IndexType row, const IndexType column) const;
        typename boost::call_traits<DataType>::const_reference
                 operator()(const IndexType block, const IndexType row, const IndexType column) const;

        //typename SparseStorageType::const_iterator begin() const;
        //typename SparseStorageType::const_iterator end() const;

        void Multiply(const DataVectorType &in,
                            DataVectorType &out);
        void Multiply(const DataType* in,
                            DataType* out);
        void MultiplySubMatrix( const IndexType blockNum,
                                DataType* in,
                                DataType* out);

        const size_t GetMemoryFootprint();
        const size_t GetMemoryFootprint(IndexType i) const;

        const unsigned long GetMulCallsCounter() const;
        const DataType   GetAvgRowDensity();
        const DataType   GetAvgRowDensity(IndexType i) const;
        const IndexType  GetBandwidth();
        const IndexType  GetBandwidth(IndexType i);
        COOMatTypeSharedPtr GetCooStorage();
        COOMatTypeSharedPtr GetCooStorage(IndexType i);

        void writeSparsityPatternTo(std::ostream& out, IndexType blockSize = 64);
        void writeSubmatrixSparsityPatternTo(std::ostream& out, const IndexType subMatrixIdx, IndexType blockSize = 64);
        void writeBlockSparsityPatternTo(std::ostream& out,
                        const IndexType blk_row = 0, const IndexType blk_col = 0, IndexType blockSize = 64);
        void writeSubmatrixBlockSparsityPatternTo(std::ostream& out, const IndexType subMatrixIdx,
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
