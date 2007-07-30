///////////////////////////////////////////////////////////////////////////////
//
// File: MatrixOperations.hpp
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
// Description: Defines the global functions needed for matrix operations.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_OPERATIONS_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_OPERATIONS_HPP

// Since this file defines all of the operations for all combination of matrix types, 
// we have to include all matrix specializations first.
#include <LibUtilities/LinearAlgebra/NormalMatrix.hpp>
#include <LibUtilities/LinearAlgebra/BlockMatrix.hpp>
#include <LibUtilities/LinearAlgebra/ScaledMatrix.hpp>
#include <LibUtilities/LinearAlgebra/Blas.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>

namespace Nektar
{
    ////////////////////////////////////////////////////////////////////////////////////
	// Multiplication
	////////////////////////////////////////////////////////////////////////////////////
    template<typename DataType, typename LhsDataType, typename LhsStorageType, typename LhsMatrixType, typename RhsType>
    void NekMultiply(NekMatrix<DataType, LhsStorageType, StandardMatrixTag>& result,
                     const NekMatrix<LhsDataType, LhsStorageType, LhsMatrixType>& lhs,
                     const RhsType& rhs)
    {
        // TODO - optimize for the different matrix types.
        result = NekMatrix<DataType, LhsStorageType, StandardMatrixTag>(lhs.GetRows(), lhs.GetColumns());
        for(unsigned int i = 0; i < lhs.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < lhs.GetColumns(); ++j)
            {
                result(i,j) = lhs(i,j)*rhs;
            }
        }
    }

    template<typename DataType, typename LhsDataType, typename LhsStorageType, typename LhsMatrixType, typename RhsType>
    void NekMultiply(boost::shared_ptr<NekMatrix<DataType, LhsStorageType, StandardMatrixTag> >& result,
                     const boost::shared_ptr<const NekMatrix<LhsDataType, LhsStorageType, LhsMatrixType> >& lhs,
                     const RhsType& rhs)
    {
        NekMultiply(*result, *lhs, rhs);
    }

    template<typename DataType, typename LhsType, typename RhsDataType, typename RhsStorageType, typename RhsMatrixType>
    void NekMultiply(NekMatrix<DataType, RhsStorageType, StandardMatrixTag>& result,
                     const LhsType& rhs,
                     const NekMatrix<RhsDataType, RhsStorageType, RhsMatrixType>& lhs)
                     
    {
        return NekMultiply(result, lhs, rhs);
    }
    
   
    template<typename DataType>
    void NekMultiply(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                            const NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& lhs,
                            const NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& rhs)
    {
        ASSERTL0(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + 
            std::string(" and a right side matrix with row count ") + 
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));
           
        result = NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>(lhs.GetRows(), rhs.GetColumns());

        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < result.GetColumns(); ++j)
            {
                DataType t = DataType(0);

                // Set the result(i,j) element.
                for(unsigned int k = 0; k < lhs.GetColumns(); ++k)
                {
                    t += lhs(i,k)*rhs(k,j);
                }
                result(i,j) = t;
            }
        }
    }
    
    #ifdef NEKTAR_USING_BLAS    
    /// \brief Floating point specialization when blas is in use.
    ///
    /// The signature of this function should match the unspecialized NekMultiply.
    void NekMultiply(NekMatrix<double, FullMatrixTag, StandardMatrixTag>& result,
                     const NekMatrix<double, FullMatrixTag, StandardMatrixTag>& lhs,
                     const NekMatrix<double, FullMatrixTag, StandardMatrixTag>& rhs);
    #endif //NEKTAR_USING_BLAS

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Matrix-Vector multiplication
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    template<typename DataType, typename LhsDataType, typename MatrixType, unsigned int dim, unsigned int space>
    void NekMultiply(NekVector<DataType, dim, space>& result,
                    const NekMatrix<LhsDataType, FullMatrixTag, MatrixType>& lhs,
                    const NekVector<DataType, dim, space>& rhs)
    {
       ASSERTL0(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
           boost::lexical_cast<std::string>(lhs.GetColumns()) + 
           std::string(" and a right side vector with row count ") + 
           boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

       result = NekVector<DataType, dim, space>(lhs.GetRows());
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

    #ifdef NEKTAR_USING_BLAS    
    /// \brief Floating point specialization when blas is in use.
    template<unsigned int dim, unsigned int space>
    void NekMultiply(NekVector<double, dim, space>& result,
                    const NekMatrix<double, FullMatrixTag, StandardMatrixTag>& lhs,
                    const NekVector<double, dim, space>& rhs)
    {
        ASSERTL0(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
           boost::lexical_cast<std::string>(lhs.GetColumns()) + 
           std::string(" and a right side vector with row count ") + 
           boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

        result = NekVector<double, dim, space>(lhs.GetRows());

        int m = lhs.GetRows();
        int n = lhs.GetColumns();
        double alpha = 1.0;
        const double* a = lhs.GetRawPtr();
        int lda = m;
        const double* x = rhs.GetPtr();
        int incx = 1;
        double beta = 0.0;
        double* y = result.GetPtr();
        int incy = 1;
        
        Blas::Dgemv('T', n, m, alpha, a, n, x, incx, beta, y, incy);
    }
    #endif //NEKTAR_USING_BLAS

    template<typename DataType, typename LhsDataType, typename MatrixType, unsigned int dim, unsigned int space>
    void NekMultiply(NekVector<DataType, dim, space>& result,
                    const NekMatrix<LhsDataType, DiagonalMatrixTag, MatrixType>& lhs,
                    const NekVector<DataType, dim, space>& rhs)
    {
       ASSERTL0(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
           boost::lexical_cast<std::string>(lhs.GetColumns()) + 
           std::string(" and a right side vector with row count ") + 
           boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

       result = NekVector<DataType, dim, space>(lhs.GetRows());
       for(unsigned int i = 0; i < lhs.GetRows(); ++i)
       {
           result[i] = lhs(i,i)*rhs(i);
       }
    }


	///////////////////////////////////////////////////////////////////
	// Addition
	///////////////////////////////////////////////////////////////////
    template<typename ResultType, typename LhsType, typename RhsType>
    void AddMatricesWithSameStorageType(ResultType& result, const LhsType& lhs, const RhsType& rhs)
    {
        typename ResultType::iterator result_iter = result.begin();
        typename RhsType::const_iterator rhs_iter = rhs.begin();
        typename LhsType::const_iterator lhs_iter = lhs.begin();
        while( result_iter != result.end() )
        {
            (*result_iter) = (*lhs_iter) + (*rhs_iter);
            ++result_iter;
            ++rhs_iter;
            ++lhs_iter;
        }
    }

	template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void NekAdd(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL0(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL0(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));
            
        result = NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>(lhs.GetRows(), lhs.GetColumns());
        return AddMatricesWithSameStorageType(result, lhs, rhs);
    }
    
    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void NekAdd(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                        const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
                        const NekMatrix<RhsDataType, DiagonalMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL0(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL0(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));
        
        result = NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>(lhs);
        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            result(i,i) += rhs(i,i);
        }
    }
    
    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void NekAdd(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                        const NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>& lhs,
                        const NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL0(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL0(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));
        
        result = NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>(rhs);
        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            result(i,i) += lhs(i,i);
        }
    }
    
    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void NekAdd(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, UpperTriangularMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL0(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL0(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));

        result = NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>(rhs);
        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            for(unsigned int j = i; j < result.GetColumns(); ++j)
            {
                result(i,i) += lhs(i,i);
            }
        }
    }

    // UT = UT + UT
    // Line 22
    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void NekAdd(NekMatrix<DataType, UpperTriangularMatrixTag, StandardMatrixTag>& result,
                const NekMatrix<LhsDataType, UpperTriangularMatrixTag, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, UpperTriangularMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL0(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL0(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));
            
        result = NekMatrix<DataType, UpperTriangularMatrixTag, StandardMatrixTag>(lhs.GetRows(), lhs.GetColumns());
        return AddMatricesWithSameStorageType(result, lhs, rhs);
    }

    // LT = LT + LT
    // Line 32
    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void NekAdd(NekMatrix<DataType, LowerTriangularMatrixTag, StandardMatrixTag>& result,
                const NekMatrix<LhsDataType, LowerTriangularMatrixTag, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, LowerTriangularMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL0(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL0(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));
            
        result = NekMatrix<DataType, LowerTriangularMatrixTag, StandardMatrixTag>(lhs.GetRows(), lhs.GetColumns());
        return AddMatricesWithSameStorageType(result, lhs, rhs);
    }

	////////////////////////////////////////////////////////////////////////////////////
	// Subtraction
	////////////////////////////////////////////////////////////////////////////////////
	template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void NekSubtract(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL0(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be subtracted."));
        ASSERTL0(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be subtracted."));
            
        result = NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>(lhs.GetRows(), lhs.GetColumns());
        typename NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>::iterator result_iter = result.begin();
        typename NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>::const_iterator rhs_iter = rhs.begin();
        typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::const_iterator lhs_iter = lhs.begin();
        while( result_iter != result.end() )
        {
            (*result_iter) = (*lhs_iter) - (*rhs_iter);
            ++result_iter;
            ++rhs_iter;
            ++lhs_iter;
        }
    }
    
    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void NekSubtract(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                        const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
                        const NekMatrix<RhsDataType, DiagonalMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL0(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be subtracted."));
        ASSERTL0(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be subtracted."));
        
        result = NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>(lhs);
        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            result(i,i) -= rhs(i,i);
        }
    }
    
    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void NekSubtract(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                        const NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>& lhs,
                        const NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL0(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be subtracted."));
        ASSERTL0(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be subtracted."));
        
        result = NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>(rhs);
        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            result(i,i) -= lhs(i,i);
        }
    }
    
    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void NekSubtract(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
                        const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
                        const NekMatrix<RhsDataType, UpperTriangularMatrixTag, RhsMatrixType>& rhs)
    {
        ASSERTL0(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be subtracted."));
        ASSERTL0(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be subtracted."));

        result = NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>(rhs);
        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            for(unsigned int j = i; j < result.GetColumns(); ++j)
            {
                result(i,i) -= lhs(i,i);
            }
        }
    }

}
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_OPERATIONS_HPP

