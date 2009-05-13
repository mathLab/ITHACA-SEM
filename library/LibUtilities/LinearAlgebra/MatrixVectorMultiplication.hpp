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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_VECTOR_MULTIPLICATION_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_VECTOR_MULTIPLICATION_HPP

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/LinearAlgebra/Blas.hpp>
#include <LibUtilities/LinearAlgebra/CanGetRawPtr.hpp>

namespace Nektar
{

    /// \brief The default matrix vector multiplication when a specialized algorithm is 
    ///        not available.
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
        for(unsigned int i = 0; i < lhs.GetRows(); ++i)
        {
            result[i] = lhs(i,i)*rhs[i];
        }
    }
     
    /// \brief Matrix vector multiplication using blas.
    template<typename MatrixType, typename LhsDataType>
    void NekMultiplyFullMatrix(double* result,
                               const NekMatrix<LhsDataType, MatrixType>& lhs,
                               const double* rhs,
                               typename boost::enable_if
                               <
                                    boost::mpl::and_
                                    <
                                        boost::is_same
                                        <
                                            typename RawType<typename NekMatrix<LhsDataType, MatrixType>::NumberType>::type,
                                            double
                                        >,
                                        CanGetRawPtr<NekMatrix<LhsDataType, MatrixType> >
                                    >
                               >::type* p = 0)
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
        
    template<typename DataType, typename LhsDataType, typename MatrixType>
    void NekMultiplyFullMatrix(DataType* result,
                               const NekMatrix<LhsDataType, MatrixType>& lhs,
                               const DataType* rhs)
    {
        NekMultiplyUnspecializedMatrixType(result, lhs, rhs);
    }
    
    template<typename LhsDataType, typename MatrixType>
    void NekMultiplyBandedMatrix(double* result,
                    const NekMatrix<LhsDataType, MatrixType>& lhs,
                    const double* rhs,
                    typename boost::enable_if
                    <
                        boost::mpl::and_
                        <
                            boost::is_same
                            <
                                typename RawType
                                <
                                    typename NekMatrix<LhsDataType, MatrixType>::NumberType
                                >::type, 
                                double
                            >,
                            CanGetRawPtr<NekMatrix<LhsDataType, MatrixType> >
                        >
                    >::type* p = 0)
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


    template<typename DataType, typename LhsDataType, typename MatrixType>
    void NekMultiplyBandedMatrix(DataType* result,
                    const NekMatrix<LhsDataType, MatrixType>& lhs,
                    const DataType* rhs,
                    typename boost::enable_if
                    <
                        boost::mpl::not_<CanGetRawPtr<NekMatrix<LhsDataType, MatrixType> > >
                    >::type* p = 0)
    {
        NEKERROR(ErrorUtil::efatal, "Banded block matrix multiplication not yet implemented");
    }
    
    template<typename DataType, typename LhsDataType, typename MatrixType>
    void NekMultiplyBandedMatrix(DataType* result,
                    const NekMatrix<LhsDataType, MatrixType>& lhs,
                    const DataType* rhs,
                    typename boost::enable_if
                    <
                        CanGetRawPtr<NekMatrix<LhsDataType, MatrixType> >
                    >::type* p = 0)
    {
        unsigned int subDiagonals = lhs.GetNumberOfSubDiagonals();
        unsigned int superDiagonals = lhs.GetNumberOfSuperDiagonals();
        unsigned int packedRows = subDiagonals+superDiagonals+1;

        const DataType* rawData = lhs.GetRawPtr();
        DataType scale = lhs.Scale();
        for(unsigned int i = 0; i < lhs.GetRows(); ++i)
        {
            result[i] = DataType(0);
        }
        
        for(unsigned int i = 0; i < lhs.GetRows(); ++i)
        {          
            unsigned int start = i <= superDiagonals ? superDiagonals-i : 0;
            unsigned int end = packedRows;
            if( lhs.GetColumns() - i <= subDiagonals )
            {
                end -= subDiagonals - (lhs.GetColumns() - i) + 1;
            }
            
            unsigned int resultOffset = (i <= superDiagonals) ? 0 : i - superDiagonals;

            for(unsigned int j = start; j < end; ++j, ++resultOffset)
            {
               result[resultOffset] += rawData[i*packedRows+j]*scale*rhs[i];
            }
        }
    }
    
    template<typename LhsDataType>
    void NekMultiplyUpperTriangularMatrix(double* result,
                     const NekMatrix<LhsDataType, StandardMatrixTag>& lhs,
                     const double* rhs,
                     typename boost::enable_if
                     <
                        boost::mpl::and_
                        <
                            boost::is_same<double, typename RawType<LhsDataType>::type>,
                            CanGetRawPtr<NekMatrix<LhsDataType, StandardMatrixTag> >
                        >
                     >::type* p = 0)
    {
        int vectorSize = lhs.GetColumns();
        std::copy(rhs, rhs+vectorSize, result);
        int n = lhs.GetRows();
        const double* a = lhs.GetRawPtr();
        double* x = result;
        int incx = 1;
        
        Blas::Dtpmv('U', lhs.GetTransposeFlag(), 'N', n, a, x, incx);
    }

    template<typename dim, typename space, typename LhsDataType>
    void NekMultiplyUpperTriangularMatrix(double* result,
                     const NekMatrix<LhsDataType, ScaledMatrixTag>& lhs,
                     const double* rhs)
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

    template<typename LhsDataType>
    void NekMultiplyLowerTriangularMatrix(double* result,
                     const NekMatrix<LhsDataType, StandardMatrixTag>& lhs,
                     const double* rhs,
                     typename boost::enable_if
                     <
                        boost::mpl::and_
                        <
                            boost::is_same<double, typename RawType<LhsDataType>::type>,
                            CanGetRawPtr<NekMatrix<LhsDataType, StandardMatrixTag> >
                        >
                     >::type* p = 0)
    {
        int vectorSize = lhs.GetColumns();
        std::copy(rhs, rhs+vectorSize, result);
        int n = lhs.GetRows();
        const double* a = lhs.GetRawPtr();
        double* x = result;
        int incx = 1;
        
        Blas::Dtpmv('L', lhs.GetTransposeFlag(), 'N', n, a, x, incx);
    }

    template<typename LhsStorageType>
    void NekMultiplyLowerTriangularMatrix(double* result,
                     const NekMatrix<LhsStorageType, ScaledMatrixTag>& lhs,
                     const double* rhs)
    {
        NekMultiplyLowerTriangularMatrix(result, *lhs.GetOwnedMatrix(), rhs);
        
        for(unsigned int i = 0; i < lhs.GetColumns(); ++i)
        {
            result[i] *= lhs.Scale();
        }
    }

    template<typename DataType, typename LhsDataType, typename MatrixType>
    void NekMultiply(DataType* result,
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
    
    /// \brief Matrix vector multiplication for normal and scaled matrices.
    template<typename DataType, typename LhsDataType, typename MatrixType, typename dim, typename space>
    void NekMultiply(NekVector<DataType, VariableSizedVector, space>& result,
                    const NekMatrix<LhsDataType, MatrixType>& lhs,
                    const NekVector<const DataType, dim, space>& rhs)
    {
                       
        ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + 
            std::string(" and a right side vector with row count ") + 
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));
        NekMultiply(result.GetRawPtr(), lhs, rhs.GetRawPtr());
    }

    template<typename DataType, typename LhsInnerMatrixType, typename dim, typename space>
    void FullBlockMatrixMultiply(NekVector<DataType, VariableSizedVector, space>& result,
                     const NekMatrix<LhsInnerMatrixType, BlockMatrixTag>& lhs,
                     const NekVector<const DataType, dim, space>& rhs)
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
                NekMultiply(temp_ptr, *block, rhsWrapper);
                for(unsigned int i = 0; i < rowsInBlock; ++i)
                {
                    resultWrapper[i] += temp_ptr[i];
                }
            }
        }
    }
    
    template<typename LhsInnerMatrixType, typename dim, typename space>
    void DiagonalBlockMatrixMultiply(NekVector<double, VariableSizedVector, space>& result,
                     const NekMatrix<LhsInnerMatrixType, BlockMatrixTag>& lhs,
                     const NekVector<const double, dim, space>& rhs)
    {
        unsigned int numberOfBlockRows = lhs.GetNumberOfBlockRows();
        double* result_ptr = result.GetRawPtr();
        const double* rhs_ptr = rhs.GetRawPtr();
        
        unsigned int curResultRow = 0;
        unsigned int curWrapperRow = 0;
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

            double* resultWrapper = result_ptr + curResultRow;
            
            unsigned int blockColumn = blockRow;

            if( blockColumn != 0 )
            {
                curWrapperRow += lhs.GetNumberOfColumnsInBlockColumn(blockColumn-1);
            }

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

            const double* rhsWrapper = rhs_ptr + curWrapperRow;
            NekMultiply(resultWrapper, *block, rhsWrapper);
            //resultWrapper = (*block)*rhsWrapper;
        }
    }
    
    template<typename DataType, typename LhsInnerMatrixType, typename dim, typename space>
    void DiagonalBlockMatrixMultiply(NekVector<DataType, VariableSizedVector, space>& result,
                     const NekMatrix<LhsInnerMatrixType, BlockMatrixTag>& lhs,
                     const NekVector<const DataType, dim, space>& rhs)
    {
        unsigned int numberOfBlockRows = lhs.GetNumberOfBlockRows();
                
        unsigned int curResultRow = 0;
        unsigned int curWrapperRow = 0;
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

            NekVector<DataType, VariableSizedVector, space> resultWrapper(rowsInBlock, result.GetPtr() + curResultRow, eWrapper);
            
            unsigned int blockColumn = blockRow;

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

            NekVector<const DataType, VariableSizedVector, space> rhsWrapper(columnsInBlock, rhs.GetPtr() + curWrapperRow, eWrapper);
            resultWrapper = (*block)*rhsWrapper;
        }
    }
    
    /// \brief Matrix/Vector multiplicatin for block matrices.
    template<typename DataType, typename LhsInnerMatrixType, typename dim, typename space>
    void NekMultiply(NekVector<DataType, VariableSizedVector, space>& result,
                     const NekMatrix<LhsInnerMatrixType, BlockMatrixTag>& lhs,
                     const NekVector<const DataType, dim, space>& rhs)
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



   
    /// We need two version of NekMultiply, one for constant sized vectors and one for 
    /// variable sized vectors because the constructors to initialize the vector to 0 are 
    /// different in each case.
    ///
    /// Note that I can't make a compile time decision about the size of the result vector
    /// because the matrix dimensions are runtime only.
    template<typename DataType, typename LhsDataType, typename MatrixType, typename dim, typename space>
    NekVector<typename boost::remove_const<DataType>::type, VariableSizedVector, space> 
    NekMultiply(const NekMatrix<LhsDataType, MatrixType>& lhs,
                const NekVector<DataType, dim, space>& rhs)
    {
       NekVector<typename boost::remove_const<DataType>::type, VariableSizedVector, space> result(lhs.GetRows(), DataType(0));
       NekMultiply(result, lhs, rhs);
       return result;
    }

    template<typename DataType, typename LhsDataType, typename MatrixType, typename space>
    NekVector<typename boost::remove_const<DataType>::type, VariableSizedVector, space> 
    NekMultiply(const NekMatrix<LhsDataType, MatrixType>& lhs,
                const NekVector<DataType, VariableSizedVector, space>& rhs)
    {
       NekVector<typename boost::remove_const<DataType>::type, VariableSizedVector, space> result(lhs.GetRows(), DataType(0));
       NekMultiply(result, lhs, rhs);
       return result;
    }
}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_VECTOR_MULTIPLICATION_HPP
