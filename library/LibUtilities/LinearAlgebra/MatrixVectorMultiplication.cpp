///////////////////////////////////////////////////////////////////////////////
//
// File: MatrixVectorMultiplication.hpp
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#include <type_traits>

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/LinearAlgebra/ExplicitInstantiation.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

#include <LibUtilities/LinearAlgebra/Blas.hpp>
#include <LibUtilities/LinearAlgebra/CanGetRawPtr.hpp>

#include <LibUtilities/LinearAlgebra/StandardMatrix.hpp>
#include <LibUtilities/LinearAlgebra/BlockMatrix.hpp>
#include <LibUtilities/LinearAlgebra/ScaledMatrix.hpp>
#include <LibUtilities/LinearAlgebra/Blas.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>
#include <LibUtilities/LinearAlgebra/MatrixOperations.hpp>

namespace Nektar
{

    template<typename DataType, typename LhsDataType, typename MatrixType>
    void NekMultiplyUnspecializedMatrixType(DataType* result,
                                            const NekMatrix<LhsDataType, MatrixType>& lhs,
                                            const DataType* rhs)
    {
        for(unsigned int i = 0; i < lhs.GetRows(); ++i)
        {
            DataType accum = DataType(0);
            for(unsigned int j = 0; j < lhs.GetColumns(); ++j)
            {
                accum += lhs(i,j)*rhs[j];
            }
            result[i] = accum;
        }
    }

    template<typename DataType, typename LhsDataType, typename MatrixType>
    void NekMultiplyDiagonalMatrix(DataType* result,
                                   const NekMatrix<LhsDataType, MatrixType>& lhs,
                                   const DataType* rhs)
    {
        int n = lhs.GetRows();
        for(unsigned int i = 0; i < n; ++i)
        {
            result[i] = lhs(i,i)*rhs[i];
        }
    }

    template<typename DataType, typename LhsDataType>
    void NekMultiplyDiagonalMatrix(DataType* result,
                                   const NekMatrix<LhsDataType, StandardMatrixTag>& lhs,
                                   const DataType* rhs)
    {
        int n = lhs.GetRows();
        const DataType* mat_ptr = lhs.GetRawPtr();
        Vmath::Vmul(n, mat_ptr, 1, rhs, 1, result, 1);
    }

    template<typename DataType, typename LhsDataType, typename MatrixType>
    void NekMultiplyBandedMatrix(DataType* result,
                    const NekMatrix<LhsDataType, MatrixType>& lhs,
                    const DataType* rhs,
                    typename std::enable_if<
                                 CanGetRawPtr<NekMatrix<LhsDataType, MatrixType> >::value >::type* p=0)
    {
        boost::ignore_unused(p);

        int m  = lhs.GetRows();
        int n  = lhs.GetColumns();
        int kl = lhs.GetNumberOfSubDiagonals();
        int ku = lhs.GetNumberOfSuperDiagonals();
        DataType alpha = lhs.Scale();
        const DataType* a = lhs.GetRawPtr();
        int lda = kl + ku + 1;
        const DataType* x = rhs;
        int incx = 1;
        DataType beta = 0.0;
        DataType* y = result;
        int incy = 1;
        Blas::Gbmv('N', m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);

    }

    template<typename DataType, typename LhsDataType>
    void NekMultiplyBandedMatrix(DataType* result,
                    const NekMatrix<LhsDataType, BlockMatrixTag>& lhs,
                    const DataType* rhs,
                    typename std::enable_if<
                                 !CanGetRawPtr<NekMatrix<LhsDataType, BlockMatrixTag> >::value>::type* p=0)
    {
        boost::ignore_unused(result, lhs, rhs, p);

        NEKERROR(ErrorUtil::efatal, "Banded block matrix multiplication not yet implemented");
    }

    template<typename DataType, typename LhsInnerMatrixType>
    void FullBlockMatrixMultiply(NekVector<DataType>& result,
                     const NekMatrix<LhsInnerMatrixType, BlockMatrixTag>& lhs,
                     const NekVector<DataType>& rhs)
    {
        unsigned int numberOfBlockRows = lhs.GetNumberOfBlockRows();
        unsigned int numberOfBlockColumns = lhs.GetNumberOfBlockColumns();
        DataType* result_ptr = result.GetRawPtr();
        const DataType* rhs_ptr = rhs.GetRawPtr();
        
        for(unsigned int i = 0; i < result.GetDimension(); ++i)
        {
            result_ptr[i] = DataType(0);
        }
        Array<OneD, DataType> temp(result.GetDimension());
        DataType* temp_ptr = temp.get();
        
        unsigned int curResultRow = 0;
        for(unsigned int blockRow = 0; blockRow < numberOfBlockRows; ++blockRow)
        {
            unsigned int rowsInBlock = lhs.GetNumberOfRowsInBlockRow(blockRow);

            if( blockRow != 0 )
            {
                curResultRow += lhs.GetNumberOfRowsInBlockRow(blockRow-1);
            }

            if( rowsInBlock == 0 )
            {
                continue;
            }

            DataType* resultWrapper = result_ptr + curResultRow;
            
            unsigned int curWrapperRow = 0;
            for(unsigned int blockColumn = 0; blockColumn < numberOfBlockColumns; ++blockColumn)
            {
                if( blockColumn != 0 )
                {
                    curWrapperRow += lhs.GetNumberOfColumnsInBlockColumn(blockColumn-1);
                }

                //const std::shared_ptr<const LhsInnerMatrixType>& block = lhs.GetBlock(blockRow, blockColumn);
                const LhsInnerMatrixType* block = lhs.GetBlockPtr(blockRow, blockColumn);
                if( !block )
                {
                    continue;
                }

                unsigned int columnsInBlock = lhs.GetNumberOfColumnsInBlockColumn(blockColumn);
                if( columnsInBlock == 0 )
                {
                    continue;
                }

                const DataType* rhsWrapper = rhs_ptr + curWrapperRow;
                Multiply(temp_ptr, *block, rhsWrapper);
                for(unsigned int i = 0; i < rowsInBlock; ++i)
                {
                    resultWrapper[i] += temp_ptr[i];
                }
            }
        }
    }

    template<typename DataType, typename LhsInnerMatrixType>
    void DiagonalBlockMatrixMultiply(NekVector<DataType>& result,
                     const NekMatrix<LhsInnerMatrixType, BlockMatrixTag>& lhs,
                     const NekVector<DataType>& rhs)
    {
        unsigned int numberOfBlockRows = lhs.GetNumberOfBlockRows();
        DataType* result_ptr = result.GetRawPtr();
        const DataType* rhs_ptr = rhs.GetRawPtr();
        
        Array<OneD, unsigned int> rowSizes;
        Array<OneD, unsigned int> colSizes;
        lhs.GetBlockSizes(rowSizes, colSizes);

        unsigned int curResultRow = 0;
        unsigned int curWrapperRow = 0;
        unsigned int rowsInBlock = 0;
        unsigned int columnsInBlock = 0;
        for(unsigned int blockRow = 0; blockRow < numberOfBlockRows; ++blockRow)
        {
            curResultRow  += rowsInBlock;
            curWrapperRow += columnsInBlock;
            if ( blockRow == 0)
            {
                rowsInBlock    = rowSizes[blockRow] + 1;
                columnsInBlock = colSizes[blockRow] + 1;
            }
            else
            {
                rowsInBlock    = rowSizes[blockRow] - rowSizes[blockRow-1];
                columnsInBlock = colSizes[blockRow] - colSizes[blockRow-1];
            }

            if( rowsInBlock == 0)
            {
                continue;
            }
            if( columnsInBlock == 0)
            {
                std::fill(result.begin()+curResultRow,
                         result.begin()+curResultRow + rowsInBlock, 0.0);
                continue;
            }

            const LhsInnerMatrixType* block = lhs.GetBlockPtr(blockRow, blockRow);
            if( !block )
            {
                continue;
            }

            DataType* resultWrapper = result_ptr + curResultRow;            
            const DataType* rhsWrapper = rhs_ptr + curWrapperRow;
            Multiply(resultWrapper, *block, rhsWrapper);
            //resultWrapper = (*block)*rhsWrapper;
        }
        curResultRow  += rowsInBlock;
        if (curResultRow < result.GetRows())
        {
            std::fill(result.begin()+curResultRow, result.end(), 0.0);
        }
    }

    void DiagonalBlockFullScalMatrixMultiply(NekVector<double>& result,
                     const NekMatrix<NekMatrix<NekMatrix<NekDouble, StandardMatrixTag>, ScaledMatrixTag>, BlockMatrixTag>& lhs,
                     const NekVector<double>& rhs)
    {
        unsigned int numberOfBlockRows = lhs.GetNumberOfBlockRows();
        double* result_ptr = result.GetRawPtr();
        const double* rhs_ptr = rhs.GetRawPtr();
        
        Array<OneD, unsigned int> rowSizes;
        Array<OneD, unsigned int> colSizes;
        lhs.GetBlockSizes(rowSizes, colSizes);

        unsigned int curResultRow = 0;
        unsigned int curWrapperRow = 0;
        unsigned int rowsInBlock = 0;
        unsigned int columnsInBlock = 0;
        for(unsigned int blockRow = 0; blockRow < numberOfBlockRows; ++blockRow)
        {
            curResultRow  += rowsInBlock;
            curWrapperRow += columnsInBlock;
            if ( blockRow == 0)
            {
                rowsInBlock    = rowSizes[blockRow] + 1;
                columnsInBlock = colSizes[blockRow] + 1;
            }
            else
            {
                rowsInBlock    = rowSizes[blockRow] - rowSizes[blockRow-1];
                columnsInBlock = colSizes[blockRow] - colSizes[blockRow-1];
            }

            if( rowsInBlock == 0)
            {
                continue;
            }
            if( columnsInBlock == 0)
            {
                std::fill(result.begin()+curResultRow,
                         result.begin()+curResultRow + rowsInBlock, 0.0);
                continue;
            }

            const DNekScalMat* block = lhs.GetBlockPtr(blockRow, blockRow);
            if( !block )
            {
                continue;
            }

            double* resultWrapper = result_ptr + curResultRow;
            const double* rhsWrapper = rhs_ptr + curWrapperRow;

            // Multiply
            const unsigned int* size = block->GetSize();
            Blas::Gemv('N', size[0], size[1], block->Scale(),
                        block->GetRawPtr(), size[0], rhsWrapper, 1,
                        0.0, resultWrapper, 1);
        }
        curResultRow  += rowsInBlock;
        if (curResultRow < result.GetRows())
        {
            std::fill(result.begin()+curResultRow, result.end(), 0.0);
        }
    }

    void DiagonalBlockFullScalMatrixMultiply(NekVector<NekSingle>& result,
                     const NekMatrix<NekMatrix<NekMatrix<NekSingle, StandardMatrixTag>, ScaledMatrixTag>, BlockMatrixTag>& lhs,
                     const NekVector<NekSingle>& rhs)
    {
        unsigned int numberOfBlockRows = lhs.GetNumberOfBlockRows();
        NekSingle* result_ptr = result.GetRawPtr();
        const NekSingle* rhs_ptr = rhs.GetRawPtr();
        
        Array<OneD, unsigned int> rowSizes;
        Array<OneD, unsigned int> colSizes;
        lhs.GetBlockSizes(rowSizes, colSizes);

        unsigned int curResultRow = 0;
        unsigned int curWrapperRow = 0;
        unsigned int rowsInBlock = 0;
        unsigned int columnsInBlock = 0;
        for(unsigned int blockRow = 0; blockRow < numberOfBlockRows; ++blockRow)
        {
            curResultRow  += rowsInBlock;
            curWrapperRow += columnsInBlock;
            if ( blockRow == 0)
            {
                rowsInBlock    = rowSizes[blockRow] + 1;
                columnsInBlock = colSizes[blockRow] + 1;
            }
            else
            {
                rowsInBlock    = rowSizes[blockRow] - rowSizes[blockRow-1];
                columnsInBlock = colSizes[blockRow] - colSizes[blockRow-1];
            }

            if( rowsInBlock == 0)
            {
                continue;
            }
            if( columnsInBlock == 0)
            {
                std::fill(result.begin()+curResultRow,
                         result.begin()+curResultRow + rowsInBlock, 0.0);
                continue;
            }

            const SNekScalMat* block = lhs.GetBlockPtr(blockRow, blockRow);
            if( !block )
            {
                continue;
            }

            NekSingle* resultWrapper = result_ptr + curResultRow;
            const NekSingle* rhsWrapper = rhs_ptr + curWrapperRow;

            // Multiply
            const unsigned int* size = block->GetSize();
            Blas::Gemv('N', size[0], size[1], block->Scale(),
                        block->GetRawPtr(), size[0], rhsWrapper, 1,
                        0.0, resultWrapper, 1);
        }
        curResultRow  += rowsInBlock;
        if (curResultRow < result.GetRows())
        {
            std::fill(result.begin()+curResultRow, result.end(), 0.0);
        }
    }

    template<typename DataType>
    void NekMultiplyLowerTriangularMatrix(DataType* result,
                     const NekMatrix<DataType, StandardMatrixTag>& lhs,
                     const DataType* rhs)
    {
        int vectorSize = lhs.GetColumns();
        std::copy(rhs, rhs+vectorSize, result);
        int n = lhs.GetRows();
        const DataType* a = lhs.GetRawPtr();
        DataType* x = result;
        int incx = 1;
        
        Blas::Tpmv('L', lhs.GetTransposeFlag(), 'N', n, a, x, incx);
    }

    template<typename DataType>
    void NekMultiplyLowerTriangularMatrix(DataType* result,
                     const NekMatrix<NekMatrix<DataType, StandardMatrixTag>, ScaledMatrixTag>& lhs,
                     const DataType* rhs)
    {
        NekMultiplyLowerTriangularMatrix(result, *lhs.GetOwnedMatrix(), rhs);
        
        for(unsigned int i = 0; i < lhs.GetColumns(); ++i)
        {
            result[i] *= lhs.Scale();
        }
    }

    template<typename DataType, typename LhsDataType, typename MatrixType>
    void NekMultiplyLowerTriangularMatrix(DataType* result,
                    const NekMatrix<LhsDataType, MatrixType>& lhs,
                    const DataType* rhs)
    {
       for(unsigned int i = 0; i < lhs.GetRows(); ++i)
       {
           DataType accum = DataType(0);
           for(unsigned int j = 0; j <= i; ++j)
           {
               accum += lhs(i,j)*rhs[j];
           }
           result[i] = accum;
       }
    }

    template<typename DataType>
    void NekMultiplyUpperTriangularMatrix(DataType* result,
                     const NekMatrix<DataType, StandardMatrixTag>& lhs,
                     const DataType* rhs)
    {
        int vectorSize = lhs.GetColumns();
        std::copy(rhs, rhs+vectorSize, result);
        int n = lhs.GetRows();
        const DataType* a = lhs.GetRawPtr();
        DataType* x = result;
        int incx = 1;
        
        Blas::Tpmv('U', lhs.GetTransposeFlag(), 'N', n, a, x, incx);
    }

    template<typename DataType>
    void NekMultiplyUpperTriangularMatrix(DataType* result,
                     const NekMatrix<NekMatrix<DataType, StandardMatrixTag>, ScaledMatrixTag>& lhs,
                     const DataType* rhs)
    {
        NekMultiplyUpperTriangularMatrix(result, *lhs.GetOwnedMatrix(), rhs);
        
        for(unsigned int i = 0; i < lhs.GetColumns(); ++i)
        {
            result[i] *= lhs.Scale();
        }
    }
    
    template<typename DataType, typename LhsDataType, typename MatrixType>
    void NekMultiplyUpperTriangularMatrix(DataType* result,
                    const NekMatrix<LhsDataType, MatrixType>& lhs,
                    const DataType* rhs)
    {
       for(unsigned int i = 0; i < lhs.GetRows(); ++i)
       {
           DataType accum = DataType(0);
           for(unsigned int j = i; j < lhs.GetColumns(); ++j)
           {
               accum += lhs(i,j)*rhs[j];
           }
           result[i] = accum;
       }
    }

    template<typename DataType, typename InnerMatrixType, typename MatrixTag>
    void NekMultiplySymmetricMatrix(DataType* result, const NekMatrix<InnerMatrixType, MatrixTag>& lhs, const DataType* rhs,
                                    typename std::enable_if<CanGetRawPtr<NekMatrix<InnerMatrixType, MatrixTag> >::value >::type* p=0)
    {
        boost::ignore_unused(p);

        const unsigned int* size = lhs.GetSize();

        DataType alpha = lhs.Scale();
        const DataType* a = lhs.GetRawPtr();
        const DataType* x = rhs;
        int incx = 1;
        DataType beta = 0.0;
        DataType* y = result;
        int incy = 1;

        Blas::Spmv('U', size[0], alpha, a, x, incx, beta, y, incy);
    }

    template<typename DataType, typename InnerMatrixType, typename MatrixTag>
    void NekMultiplySymmetricMatrix(DataType* result, const NekMatrix<InnerMatrixType, MatrixTag>& lhs, const DataType* rhs,
                                    typename std::enable_if<!CanGetRawPtr<NekMatrix<InnerMatrixType, MatrixTag> >::value >::type* p = 0)
    {
        boost::ignore_unused(p);

        NekMultiplyUnspecializedMatrixType(result, lhs, rhs);
    }

    template<typename DataType, typename InnerMatrixType, typename MatrixTag>
    void NekMultiplyFullMatrix(DataType* result, const NekMatrix<InnerMatrixType, MatrixTag>& lhs, const DataType* rhs,
                               typename std::enable_if<CanGetRawPtr<NekMatrix<InnerMatrixType, MatrixTag> >::value >::type* p=0)
    {
        boost::ignore_unused(p);

        const unsigned int* size = lhs.GetSize();

        char t = lhs.GetTransposeFlag();

        DataType alpha = lhs.Scale();
        const DataType* a = lhs.GetRawPtr();
        int lda = size[0];
        const DataType* x = rhs;
        int incx = 1;
        DataType beta = 0.0;
        DataType* y = result;
        int incy = 1;
        
        Blas::Gemv(t, size[0], size[1], alpha, a, lda, x, incx, beta, y, incy);
    }

    template<typename DataType, typename InnerMatrixType, typename MatrixTag>
    void NekMultiplyFullMatrix(DataType* result, const NekMatrix<InnerMatrixType, MatrixTag>& lhs, const DataType* rhs,
        typename std::enable_if<!CanGetRawPtr<NekMatrix<InnerMatrixType, MatrixTag> >::value>::type* p = 0)
    {
        boost::ignore_unused(p);

        NekMultiplyUnspecializedMatrixType(result, lhs, rhs);
    }

    template<typename DataType, typename LhsDataType, typename MatrixType>
    void Multiply(DataType* result,
                  const NekMatrix<LhsDataType, MatrixType>& lhs,
                  const DataType* rhs)
    {
        switch(lhs.GetType())
        {
            case eFULL:
                NekMultiplyFullMatrix(result, lhs, rhs);
                break;
            case eDIAGONAL:
                NekMultiplyDiagonalMatrix(result, lhs, rhs);
                break;
            case eUPPER_TRIANGULAR:
                NekMultiplyUpperTriangularMatrix(result, lhs, rhs);
                break;
            case eLOWER_TRIANGULAR:
                NekMultiplyLowerTriangularMatrix(result, lhs, rhs);
                break;
            case eSYMMETRIC:
                NekMultiplySymmetricMatrix(result, lhs, rhs);
                break;
            case eBANDED:
                NekMultiplyBandedMatrix(result, lhs, rhs);
                break;
            case eSYMMETRIC_BANDED:
            case eUPPER_TRIANGULAR_BANDED:
            case eLOWER_TRIANGULAR_BANDED:
            default:
                NekMultiplyUnspecializedMatrixType(result, lhs, rhs);
        } 
    }

    template<typename DataType, typename LhsDataType, typename MatrixType>
    void Multiply(NekVector<DataType>& result,
                  const NekMatrix<LhsDataType, MatrixType>& lhs,
                  const NekVector<DataType>& rhs)
    {
                       
        ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
            std::to_string(lhs.GetColumns()) + 
            std::string(" and a right side vector with row count ") + 
            std::to_string(rhs.GetRows()) + std::string(" can't be multiplied."));
        Multiply(result.GetRawPtr(), lhs, rhs.GetRawPtr());
    }

    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_SINGLE_MATRIX(Multiply, NEKTAR_STANDARD_AND_SCALED_MATRICES, (1, (void)), (1, (NekVector<NekDouble>&)), (1,(const NekVector<NekDouble>&)))
    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_SINGLE_MATRIX(Multiply, NEKTAR_STANDARD_AND_SCALED_MATRICES_SINGLE, (1, (void)), (1, (NekVector<NekSingle>&)), (1,(const NekVector<NekSingle>&)))

    template<typename DataType, typename LhsInnerMatrixType>
    void Multiply(NekVector<DataType>& result,
                  const NekMatrix<LhsInnerMatrixType, BlockMatrixTag>& lhs,
                  const NekVector<DataType>& rhs)
    {
        if( lhs.GetStorageType() == eDIAGONAL )
        {
            DiagonalBlockMatrixMultiply(result, lhs, rhs);
        }
        else
        {
            FullBlockMatrixMultiply(result, lhs, rhs);
        }
    }

    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_SINGLE_MATRIX(Multiply, NEKTAR_BLOCK_MATRIX_TYPES, (1, (void)), (1, (NekVector<NekDouble>&)), (1,(const NekVector<NekDouble>&)))
    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_SINGLE_MATRIX(Multiply, NEKTAR_BLOCK_MATRIX_TYPES_SINGLE, (1, (void)), (1, (NekVector<NekSingle>&)), (1,(const NekVector<NekSingle>&)))

    template<typename DataType, typename LhsDataType, typename MatrixType>
    NekVector<DataType> 
    Multiply(const NekMatrix<LhsDataType, MatrixType>& lhs,
             const NekVector<DataType>& rhs)
    {
       NekVector<DataType> result(lhs.GetRows(), DataType(0));
       Multiply(result, lhs, rhs);
       return result;
    }

    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_SINGLE_MATRIX(Multiply, NEKTAR_ALL_MATRIX_TYPES, (1, (NekVector<NekDouble>)), (0, ()), (1,(const NekVector<NekDouble>&)))
    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_SINGLE_MATRIX(Multiply, NEKTAR_ALL_MATRIX_TYPES_SINGLE, (1, (NekVector<NekSingle>)), (0, ()), (1,(const NekVector<NekSingle>&)))

}
