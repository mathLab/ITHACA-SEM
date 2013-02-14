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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/ExplicitInstantiation.h>
#include <LibUtilities/LinearAlgebra/MatrixVectorMultiplication.hpp>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/RawType.hpp>

#include <LibUtilities/LinearAlgebra/Blas.hpp>
#include <LibUtilities/LinearAlgebra/CanGetRawPtr.hpp>

#include <ExpressionTemplates/CommutativeTraits.hpp>
#include <boost/type_traits.hpp>

#include <LibUtilities/LinearAlgebra/StandardMatrix.hpp>
#include <LibUtilities/LinearAlgebra/BlockMatrix.hpp>
#include <LibUtilities/LinearAlgebra/ScaledMatrix.hpp>
#include <LibUtilities/LinearAlgebra/Blas.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>
#include <LibUtilities/BasicUtils/OperatorGenerators.hpp>
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

    template<typename LhsDataType, typename MatrixType>
    void NekMultiplyBandedMatrix(NekDouble* result,
                    const NekMatrix<LhsDataType, MatrixType>& lhs,
                    const NekDouble* rhs,
                    typename boost::enable_if
                    <
                        CanGetRawPtr<NekMatrix<LhsDataType, MatrixType> >
                    >::type* p=0)
    {
        int m = lhs.GetRows();
        int n = lhs.GetColumns();
        int kl = lhs.GetNumberOfSubDiagonals();
        int ku = lhs.GetNumberOfSuperDiagonals();
        double alpha = lhs.Scale();
        const double* a = lhs.GetRawPtr();
        int lda = kl + ku + 1;
        const double* x = rhs;
        int incx = 1;
        double beta = 0.0;
        double* y = result;
        int incy = 1;
        Blas::Dgbmv('N', m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);

    }

    template<typename DataType, typename LhsDataType>
    void NekMultiplyBandedMatrix(DataType* result,
                    const NekMatrix<LhsDataType, BlockMatrixTag>& lhs,
                    const DataType* rhs,
                    typename boost::enable_if
                    <
                        boost::mpl::not_<CanGetRawPtr<NekMatrix<LhsDataType, BlockMatrixTag> > >
                    >::type* p=0)
    {
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

                //const boost::shared_ptr<const LhsInnerMatrixType>& block = lhs.GetBlock(blockRow, blockColumn);
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

    template<typename LhsInnerMatrixType>
    void DiagonalBlockMatrixMultiply(NekVector<double>& result,
                     const NekMatrix<LhsInnerMatrixType, BlockMatrixTag>& lhs,
                     const NekVector<double>& rhs)
    {
        unsigned int numberOfBlockRows = lhs.GetNumberOfBlockRows();
        double* result_ptr = result.GetRawPtr();
        const double* rhs_ptr = rhs.GetRawPtr();
        
        std::fill(result.begin(), result.end(), 0.0);
        
        unsigned int curResultRow = 0;
        unsigned int curWrapperRow = 0;
        for(unsigned int blockRow = 0; blockRow < numberOfBlockRows; ++blockRow)
        {

            if( blockRow != 0 )
            {
                curResultRow += lhs.GetNumberOfRowsInBlockRow(blockRow-1);
            }

            unsigned int blockColumn = blockRow;
            if( blockColumn != 0 )
            {
                curWrapperRow += lhs.GetNumberOfColumnsInBlockColumn(blockColumn-1);
            }

            unsigned int rowsInBlock = lhs.GetNumberOfRowsInBlockRow(blockRow);
            if( rowsInBlock == 0 )
            {
                continue;
            }

            unsigned int columnsInBlock = lhs.GetNumberOfColumnsInBlockColumn(blockColumn);
            if( columnsInBlock == 0 )
            {
                continue;
            }

            const LhsInnerMatrixType* block = lhs.GetBlockPtr(blockRow, blockColumn);
            if( !block )
            {
                continue;
            }

            double* resultWrapper = result_ptr + curResultRow;            
            const double* rhsWrapper = rhs_ptr + curWrapperRow;
            Multiply(resultWrapper, *block, rhsWrapper);
            //resultWrapper = (*block)*rhsWrapper;
        }
    }

    void NekMultiplyLowerTriangularMatrix(NekDouble* result,
                     const NekMatrix<NekDouble, StandardMatrixTag>& lhs,
                     const NekDouble* rhs)
    {
        int vectorSize = lhs.GetColumns();
        std::copy(rhs, rhs+vectorSize, result);
        int n = lhs.GetRows();
        const double* a = lhs.GetRawPtr();
        double* x = result;
        int incx = 1;
        
        Blas::Dtpmv('L', lhs.GetTransposeFlag(), 'N', n, a, x, incx);
    }

    void NekMultiplyLowerTriangularMatrix(NekDouble* result,
                     const NekMatrix<NekMatrix<NekDouble, StandardMatrixTag>, ScaledMatrixTag>& lhs,
                     const NekDouble* rhs)
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

    void NekMultiplyUpperTriangularMatrix(NekDouble* result,
                     const NekMatrix<NekDouble, StandardMatrixTag>& lhs,
                     const NekDouble* rhs)
    {
        int vectorSize = lhs.GetColumns();
        std::copy(rhs, rhs+vectorSize, result);
        int n = lhs.GetRows();
        const double* a = lhs.GetRawPtr();
        double* x = result;
        int incx = 1;
        
        Blas::Dtpmv('U', lhs.GetTransposeFlag(), 'N', n, a, x, incx);
    }

    void NekMultiplyUpperTriangularMatrix(NekDouble* result,
                     const NekMatrix<NekMatrix<NekDouble, StandardMatrixTag>, ScaledMatrixTag>& lhs,
                     const NekDouble* rhs)
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

    template<typename InnerMatrixType, typename MatrixTag>
    void NekMultiplyFullMatrix(NekDouble* result, const NekMatrix<InnerMatrixType, MatrixTag>& lhs, const NekDouble* rhs,
        typename boost::enable_if
        <
            CanGetRawPtr<NekMatrix<InnerMatrixType, MatrixTag> >
        >::type* p=0)
    {
        int m = lhs.GetRows();
        int n = lhs.GetColumns();

        char t = lhs.GetTransposeFlag();
        if( t == 'T' )
        {
            std::swap(m, n);
        }

        double alpha = lhs.Scale();
        const double* a = lhs.GetRawPtr();
        int lda = m;
        const double* x = rhs;
        int incx = 1;
        double beta = 0.0;
        double* y = result;
        int incy = 1;
        
        Blas::Dgemv(t, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }

    template<typename InnerMatrixType, typename MatrixTag>
    void NekMultiplyFullMatrix(NekDouble* result, const NekMatrix<InnerMatrixType, MatrixTag>& lhs, const NekDouble* rhs,
        typename boost::enable_if
        <
            boost::mpl::not_<CanGetRawPtr<NekMatrix<InnerMatrixType, MatrixTag> > >
        >::type* p = 0)
    {
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
                NekMultiplyUnspecializedMatrixType(result, lhs, rhs);
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
            boost::lexical_cast<std::string>(lhs.GetColumns()) + 
            std::string(" and a right side vector with row count ") + 
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));
        Multiply(result.GetRawPtr(), lhs, rhs.GetRawPtr());
    }

    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_SINGLE_MATRIX(Multiply, NEKTAR_STANDARD_AND_SCALED_MATRICES, (1, (void)), (1, (NekVector<NekDouble>&)), (1,(const NekVector<NekDouble>&)));

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

    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_SINGLE_MATRIX(Multiply, NEKTAR_BLOCK_MATRIX_TYPES, (1, (void)), (1, (NekVector<NekDouble>&)), (1,(const NekVector<NekDouble>&)));

    template<typename DataType, typename LhsDataType, typename MatrixType>
    NekVector<DataType> 
    Multiply(const NekMatrix<LhsDataType, MatrixType>& lhs,
             const NekVector<DataType>& rhs)
    {
       NekVector<DataType> result(lhs.GetRows(), DataType(0));
       Multiply(result, lhs, rhs);
       return result;
    }

    NEKTAR_GENERATE_EXPLICIT_FUNCTION_INSTANTIATION_SINGLE_MATRIX(Multiply, NEKTAR_ALL_MATRIX_TYPES, (1, (NekVector<NekDouble>)), (0, ()), (1,(const NekVector<NekDouble>&)));


}
