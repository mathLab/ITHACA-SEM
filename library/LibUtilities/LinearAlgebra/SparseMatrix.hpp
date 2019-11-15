///////////////////////////////////////////////////////////////////////////////
//
// File: SparseMatrix.hpp
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
// Description: generic sparse matrix class templated by underlying sparse
//              storage format
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SPARSE_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SPARSE_MATRIX_HPP

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
     * This is a class-container to a single sparse matrices.
     * The type of sparse entry is defined with template parameter.
     *
     */
    template<typename SparseStorageType>
    class NekSparseMatrix
    {
    public:

        typedef SparseStorageType                    StorageType;
        typedef typename SparseStorageType::DataType DataType;
        typedef std::shared_ptr<SparseStorageType>   SparseStorageSharedPtr;
        typedef Array<OneD, DataType>                DataVectorType;
        typedef Array<OneD, const DataType>          ConstDataVectorType;


        LIB_UTILITIES_EXPORT NekSparseMatrix(const SparseStorageSharedPtr& sparseStoragePtr);
        LIB_UTILITIES_EXPORT NekSparseMatrix(const NekSparseMatrix& src);
        LIB_UTILITIES_EXPORT ~NekSparseMatrix();

        LIB_UTILITIES_EXPORT IndexType GetRows() const;
        LIB_UTILITIES_EXPORT IndexType GetColumns() const;
        LIB_UTILITIES_EXPORT IndexType GetNumNonZeroEntries() const;

        LIB_UTILITIES_EXPORT const DataType  GetFillInRatio() const;
        LIB_UTILITIES_EXPORT size_t          GetMemoryFootprint() const;
        LIB_UTILITIES_EXPORT unsigned long   GetMulCallsCounter() const;
        LIB_UTILITIES_EXPORT const DataType  GetAvgRowDensity() const;
        LIB_UTILITIES_EXPORT IndexType       GetBandwidth();
        LIB_UTILITIES_EXPORT COOMatTypeSharedPtr GetCooStorage();



        LIB_UTILITIES_EXPORT typename boost::call_traits<DataType>::const_reference
                 operator()(const IndexType row, const IndexType column) const;
        LIB_UTILITIES_EXPORT typename SparseStorageType::const_iterator begin() const;
        LIB_UTILITIES_EXPORT typename SparseStorageType::const_iterator end() const;

        LIB_UTILITIES_EXPORT void Multiply(const DataVectorType &in,
                            DataVectorType &out);
        LIB_UTILITIES_EXPORT void Multiply(const DataType* in,
                            DataType* out);


        LIB_UTILITIES_EXPORT void writeSparsityPatternTo(std::ostream& out, IndexType blockSize = 64);
        LIB_UTILITIES_EXPORT void writeBlockSparsityPatternTo(std::ostream& out,
                        const IndexType blk_row = 0, const IndexType blk_col = 0, IndexType blockSize = 64);

    protected:

        unsigned long           m_mulCallsCounter;
        SparseStorageSharedPtr  m_sparseStorage;

    private:

    };

}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SPARSE_MATRIX_HPP
