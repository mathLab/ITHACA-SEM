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

namespace Nektar
{

    /// \brief The default matrix vector multiplication when a specialized algorithm is 
    ///        not available.
    template<typename DataType, typename LhsDataType, typename MatrixType, typename dim, typename space>
    void NekMultiplyUnspecializedMatrixType(NekVector<DataType, VariableSizedVector, space>& result,
                                            const NekMatrix<LhsDataType, MatrixType>& lhs,
                                            const NekVector<const DataType, dim, space>& rhs)
    {
        for(unsigned int i = 0; i < lhs.GetRows(); ++i)
        {
            DataType accum = DataType(0);
            for(unsigned int j = 0; j < lhs.GetColumns(); ++j)
            {
                accum += lhs(i,j)*rhs(j);
            }
            result[i] = accum;
        }
    }
    
    template<typename DataType, typename LhsDataType, typename MatrixType, typename dim, typename space>
    void NekMultiplyDiagonalMatrix(NekVector<DataType, VariableSizedVector, space>& result,
                                   const NekMatrix<LhsDataType, MatrixType>& lhs,
                                   const NekVector<const DataType, dim, space>& rhs)
    {
        for(unsigned int i = 0; i < lhs.GetRows(); ++i)
        {
            result[i] = lhs(i,i)*rhs(i);
        }
    }
     
    /// \brief Matrix vector multiplication using blas.
    template<typename dim, typename space, typename MatrixType, typename LhsDataType>
    void NekMultiplyFullMatrix(NekVector<double, VariableSizedVector, space>& result,
                               const NekMatrix<LhsDataType, MatrixType>& lhs,
                               const NekVector<const double, dim, space>& rhs)
    {
        ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
           boost::lexical_cast<std::string>(lhs.GetColumns()) + 
           std::string(" and a right side vector with row count ") + 
           boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

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
        const double* x = rhs.GetRawPtr();
        int incx = 1;
        double beta = 0.0;
        double* y = result.GetRawPtr();
        int incy = 1;
        
        Blas::Dgemv(t, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }
        
    template<typename DataType, typename LhsDataType, typename MatrixType, typename dim, typename space>
    void NekMultiplyFullMatrix(NekVector<DataType, VariableSizedVector, space>& result,
                               const NekMatrix<LhsDataType, MatrixType>& lhs,
                               const NekVector<const DataType, dim, space>& rhs)
    {
        NekMultiplyUnspecializedMatrixType(result, lhs, rhs);
    }
    
    template<typename LhsDataType, typename MatrixType, typename dim, typename space>
    void NekMultiplyBandedMatrix(NekVector<double, VariableSizedVector, space>& result,
                    const NekMatrix<LhsDataType, MatrixType>& lhs,
                    const NekVector<const double, dim, space>& rhs)
    {
        ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
           boost::lexical_cast<std::string>(lhs.GetColumns()) + 
           std::string(" and a right side vector with row count ") + 
           boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));
   
        int m = lhs.GetRows();
        int n = lhs.GetColumns();
        int kl = lhs.GetNumberOfSubDiagonals();
        int ku = lhs.GetNumberOfSuperDiagonals();
        double alpha = lhs.Scale();
        const double* a = lhs.GetRawPtr();
        int lda = kl + ku + 1;
        const double* x = rhs.GetRawPtr();
        int incx = 1;
        double beta = 0.0;
        double* y = result.GetRawPtr();
        int incy = 1;
        Blas::Dgbmv('N', m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);

    }

    template<typename DataType, typename LhsDataType, typename MatrixType, typename dim, typename space>
    void NekMultiplyBandedMatrix(NekVector<DataType, VariableSizedVector, space>& result,
                    const NekMatrix<LhsDataType, MatrixType>& lhs,
                    const NekVector<const DataType, dim, space>& rhs)
    {
        ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
           boost::lexical_cast<std::string>(lhs.GetColumns()) + 
           std::string(" and a right side vector with row count ") + 
           boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

        unsigned int subDiagonals = lhs.GetNumberOfSubDiagonals();
        unsigned int superDiagonals = lhs.GetNumberOfSuperDiagonals();
        unsigned int packedRows = subDiagonals+superDiagonals+1;

        const DataType* rawData = lhs.GetRawPtr();
        DataType scale = lhs.Scale();
        for(unsigned int i = 0; i < result.GetDimension(); ++i)
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
    
    template<typename LhsDataType, typename dim, typename space>
    void NekMultiplyUpperTriangularMatrix(NekVector<double, VariableSizedVector, space>& result,
                     const NekMatrix<LhsDataType, StandardMatrixTag>& lhs,
                     const NekVector<const double, dim, space>& rhs,
                     typename boost::enable_if<boost::is_same<double, typename RawType<LhsDataType>::type> >::type* p = 0)
    {
        ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
           boost::lexical_cast<std::string>(lhs.GetColumns()) + 
           std::string(" and a right side vector with row count ") + 
           boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

        result = rhs;
        int n = lhs.GetRows();
        const double* a = lhs.GetRawPtr();
        double* x = result.GetRawPtr();
        int incx = 1;
        
        Blas::Dtpmv('U', lhs.GetTransposeFlag(), 'N', n, a, x, incx);
    }

    template<typename dim, typename space, typename LhsDataType>
    void NekMultiplyUpperTriangularMatrix(NekVector<double, VariableSizedVector, space>& result,
                     const NekMatrix<LhsDataType, ScaledMatrixTag>& lhs,
                     const NekVector<const double, dim, space>& rhs)
    {
        ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
           boost::lexical_cast<std::string>(lhs.GetColumns()) + 
           std::string(" and a right side vector with row count ") + 
           boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

        NekMultiplyUpperTriangularMatrix(result, *lhs.GetOwnedMatrix(), rhs);
        result *= lhs.Scale();
    }
    
    template<typename DataType, typename LhsDataType, typename MatrixType, typename dim, typename space>
    void NekMultiplyUpperTriangularMatrix(NekVector<DataType, VariableSizedVector, space>& result,
                    const NekMatrix<LhsDataType, MatrixType>& lhs,
                    const NekVector<const DataType, dim, space>& rhs)
    {
       ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
           boost::lexical_cast<std::string>(lhs.GetColumns()) + 
           std::string(" and a right side vector with row count ") + 
           boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

       for(unsigned int i = 0; i < lhs.GetRows(); ++i)
       {
           DataType accum = DataType(0);
           for(unsigned int j = i; j < lhs.GetColumns(); ++j)
           {
               accum += lhs(i,j)*rhs(j);
           }
           result[i] = accum;
       }
    }


    template<typename DataType, typename LhsDataType, typename MatrixType, typename dim, typename space>
    void NekMultiplyLowerTriangularMatrix(NekVector<DataType, VariableSizedVector, space>& result,
                    const NekMatrix<LhsDataType, MatrixType>& lhs,
                    const NekVector<const DataType, dim, space>& rhs)
    {
       ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
           boost::lexical_cast<std::string>(lhs.GetColumns()) + 
           std::string(" and a right side vector with row count ") + 
           boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

       for(unsigned int i = 0; i < lhs.GetRows(); ++i)
       {
           DataType accum = DataType(0);
           for(unsigned int j = 0; j <= i; ++j)
           {
               accum += lhs(i,j)*rhs(j);
           }
           result[i] = accum;
       }
    }

    template<typename LhsDataType, typename dim, typename space>
    void NekMultiplyLowerTriangularMatrix(NekVector<double, VariableSizedVector, space>& result,
                     const NekMatrix<LhsDataType, StandardMatrixTag>& lhs,
                     const NekVector<const double, dim, space>& rhs,
                     typename boost::enable_if
                     <
                        boost::is_same<double, typename RawType<LhsDataType>::type>
                     >::type* p = 0)
    {
        ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
           boost::lexical_cast<std::string>(lhs.GetColumns()) + 
           std::string(" and a right side vector with row count ") + 
           boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

        result = rhs;
        int n = lhs.GetRows();
        const double* a = lhs.GetRawPtr();
        double* x = result.GetRawPtr();
        int incx = 1;
        
        Blas::Dtpmv('L', lhs.GetTransposeFlag(), 'N', n, a, x, incx);
    }

    template<typename dim, typename space, typename LhsStorageType>
    void NekMultiplyLowerTriangularMatrix(NekVector<double, VariableSizedVector, space>& result,
                     const NekMatrix<LhsStorageType, ScaledMatrixTag>& lhs,
                     const NekVector<const double, dim, space>& rhs)
    {
        ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
           boost::lexical_cast<std::string>(lhs.GetColumns()) + 
           std::string(" and a right side vector with row count ") + 
           boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

        NekMultiplyLowerTriangularMatrix(result, *lhs.GetOwnedMatrix(), rhs);
        result *= lhs.Scale();
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

    /// \brief Matrix/Vector multiplicatin for block matrices.
    template<typename DataType, typename LhsInnerMatrixType, typename dim, typename space>
    void NekMultiply(NekVector<DataType, VariableSizedVector, space>& result,
                     const NekMatrix<LhsInnerMatrixType, BlockMatrixTag>& lhs,
                     const NekVector<const DataType, dim, space>& rhs)
    {
        unsigned int numberOfBlockRows = lhs.GetNumberOfBlockRows();
        unsigned int numberOfBlockColumns = lhs.GetNumberOfBlockColumns();

        for(unsigned int i = 0; i < result.GetDimension(); ++i)
        {
            result[i] = DataType(0);
        }
        
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

            NekVector<DataType, VariableSizedVector, space> resultWrapper(rowsInBlock, result.GetPtr() + curResultRow, eWrapper);

            unsigned int curWrapperRow = 0;
            for(unsigned int blockColumn = 0; blockColumn < numberOfBlockColumns; ++blockColumn)
            {
                if( blockColumn != 0 )
                {
                    curWrapperRow += lhs.GetNumberOfColumnsInBlockColumn(blockColumn-1);
                }

                const boost::shared_ptr<const LhsInnerMatrixType>& block = lhs.GetBlock(blockRow, blockColumn);
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
                
                resultWrapper += (*block)*rhsWrapper;
            }
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

#endif NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_VECTOR_MULTIPLICATION_HPP
