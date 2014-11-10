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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_OPERATIONS_DECLARATIONS_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_OPERATIONS_DECLARATIONS_HPP

// Since this file defines all of the operations for all combination of matrix types, 
// we have to include all matrix specializations first.
#include <LibUtilities/LinearAlgebra/StandardMatrix.hpp>
#include <LibUtilities/LinearAlgebra/BlockMatrix.hpp>
#include <LibUtilities/LinearAlgebra/ScaledMatrix.hpp>
#include <LibUtilities/LinearAlgebra/Blas.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>
#include <LibUtilities/BasicUtils/OperatorGenerators.hpp>
#include <LibUtilities/BasicUtils/RawType.hpp>
#include <LibUtilities/LinearAlgebra/MatrixVectorMultiplication.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrixMetadata.hpp>
#include <LibUtilities/BasicUtils/RawType.hpp>
#include <LibUtilities/LinearAlgebra/CanGetRawPtr.hpp>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>

#include <string>

namespace Nektar
{
    ////////////////////////////////////////////////////////////////////////////////////
    // Matrix-Vector Multiplication
    ////////////////////////////////////////////////////////////////////////////////////
    template<typename DataType, typename LhsDataType, typename MatrixType>
    NekVector<DataType> 
    Multiply(const NekMatrix<LhsDataType, MatrixType>& lhs,
             const NekVector<DataType>& rhs);
    
    template<typename DataType, typename LhsDataType, typename MatrixType>
    void Multiply(NekVector<DataType>& result,
                  const NekMatrix<LhsDataType, MatrixType>& lhs,
                  const NekVector<DataType>& rhs);
    
    template<typename DataType, typename LhsInnerMatrixType>
    void Multiply(NekVector<DataType>& result,
                  const NekMatrix<LhsInnerMatrixType, BlockMatrixTag>& lhs,
                  const NekVector<DataType>& rhs);

    ////////////////////////////////////////////////////////////////////////////////////
    // Matrix-Constant Multiplication
    ////////////////////////////////////////////////////////////////////////////////////
    template<typename ResultDataType, typename LhsDataType, typename LhsMatrixType>
    void Multiply(NekMatrix<ResultDataType, StandardMatrixTag>& result,
                     const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                     const ResultDataType& rhs);

    template<typename DataType, typename LhsDataType, typename LhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType, StandardMatrixTag>
    Multiply(const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                const DataType& rhs);

    template<typename RhsDataType, typename RhsMatrixType, typename ResultDataType>
    void Multiply(NekMatrix<ResultDataType, StandardMatrixTag>& result,
                     const ResultDataType& lhs,
                     const NekMatrix<RhsDataType, RhsMatrixType>& rhs);

    template<typename DataType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<RhsDataType, RhsMatrixType>::NumberType, StandardMatrixTag>
    Multiply(const DataType& lhs,
                const NekMatrix<RhsDataType, RhsMatrixType>& rhs);
    
    template<typename LhsDataType>
    void MultiplyEqual(NekMatrix<LhsDataType, StandardMatrixTag>& lhs,
                       typename boost::call_traits<LhsDataType>::const_reference rhs);
    
    
    ///////////////////////////////////////////////////////////////////
    // Matrix-Matrix Multipliation
    //////////////////////////////////////////////////////////////////
    
    template<typename LhsDataType, typename RhsDataType,
             typename LhsMatrixType, typename RhsMatrixType>
    void NekMultiplyFullMatrixFullMatrix(NekMatrix<NekDouble, StandardMatrixTag>& result,
                                         const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                                         const NekMatrix<RhsDataType, RhsMatrixType>& rhs,
                                         typename boost::enable_if
                                         <
                                            boost::mpl::and_
                                            <
                                                CanGetRawPtr<NekMatrix<LhsDataType, LhsMatrixType> >,
                                                CanGetRawPtr<NekMatrix<RhsDataType, RhsMatrixType> >
                                            >
                                         >::type* p = 0)
    {
        ASSERTL1(lhs.GetType() == eFULL && rhs.GetType() == eFULL, "Only full matrices are supported.");

        unsigned int M = lhs.GetRows();
        unsigned int N = rhs.GetColumns();
        unsigned int K = lhs.GetColumns();

        unsigned int LDA = M;
        if( lhs.GetTransposeFlag() == 'T' )
        {
            LDA = K;
        }

        unsigned int LDB = K;
        if( rhs.GetTransposeFlag() == 'T' )
        {
            LDB = N;
        }

        Blas::Dgemm(lhs.GetTransposeFlag(), rhs.GetTransposeFlag(), M, N, K,
                    lhs.Scale()*rhs.Scale(), lhs.GetRawPtr(), LDA,
                    rhs.GetRawPtr(), LDB, 0.0,
                    result.GetRawPtr(), result.GetRows());
    }
    
    
    
    template<typename LhsDataType, typename RhsDataType, typename DataType,
             typename LhsMatrixType, typename RhsMatrixType>
    void Multiply(NekMatrix<DataType, StandardMatrixTag>& result,
                     const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                     const NekMatrix<RhsDataType, RhsMatrixType>& rhs);
        
    template<typename RhsInnerType, typename RhsMatrixType>
    void MultiplyEqual(NekMatrix<double, StandardMatrixTag>& result,
                          const NekMatrix<RhsInnerType, RhsMatrixType>& rhs,
                          typename boost::enable_if
                          <
                            boost::mpl::and_
                            <
                                boost::is_same<typename RawType<typename NekMatrix<RhsInnerType, RhsMatrixType>::NumberType>::type, double>,
                                CanGetRawPtr<NekMatrix<RhsInnerType, RhsMatrixType> >
                            >
                          >::type* t = 0)
    {
        ASSERTL0(result.GetType() == eFULL && rhs.GetType() == eFULL, "Only full matrices supported.");
        unsigned int M = result.GetRows();
        unsigned int N = rhs.GetColumns();
        unsigned int K = result.GetColumns();

        unsigned int LDA = M;
        if( result.GetTransposeFlag() == 'T' )
        {
            LDA = K;
        }

        unsigned int LDB = K;
        if( rhs.GetTransposeFlag() == 'T' )
        {
            LDB = N;
        }
        double scale = rhs.Scale();
        Array<OneD, double>& buf = result.GetTempSpace();
        Blas::Dgemm(result.GetTransposeFlag(), rhs.GetTransposeFlag(), M, N, K,
            scale, result.GetRawPtr(), LDA, rhs.GetRawPtr(), LDB, 0.0,
            buf.data(), result.GetRows());
        result.SetSize(result.GetRows(), rhs.GetColumns());
        result.SwapTempAndDataBuffers();
    }
    
    template<typename DataType, typename RhsInnerType, typename RhsMatrixType>
    void MultiplyEqual(NekMatrix<DataType, StandardMatrixTag>& result,
                          const NekMatrix<RhsInnerType, RhsMatrixType>& rhs,
                          typename boost::enable_if
                          <
                            boost::mpl::or_
                            <
                                boost::mpl::not_<boost::is_same<typename RawType<typename NekMatrix<RhsInnerType, RhsMatrixType>::NumberType>::type, double> >,
                                boost::mpl::not_<CanGetRawPtr<NekMatrix<RhsInnerType, RhsMatrixType> > >
                            >
                          >::type* t = 0)
    {
        ASSERTL1(result.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
            boost::lexical_cast<std::string>(result.GetColumns()) + 
            std::string(" and a right side matrix with row count ") + 
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));
        NekMatrix<DataType, StandardMatrixTag> temp(result.GetRows(), result.GetColumns());
        
        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < result.GetColumns(); ++j)
            {
                DataType t = DataType(0);

                // Set the result(i,j) element.
                for(unsigned int k = 0; k < result.GetColumns(); ++k)
                {
                    t += result(i,k)*rhs(k,j);
                }
                temp(i,j) = t;
            }
        }
        
        result = temp;
    }

	template<typename LhsDataType, typename RhsDataType,
             typename LhsMatrixType, typename RhsMatrixType>
    NekMatrix<typename boost::remove_const<typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType>::type, StandardMatrixTag> 
    Multiply(const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        typedef typename boost::remove_const<typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType>::type NumberType;
        NekMatrix<NumberType, StandardMatrixTag> result(lhs.GetRows(), rhs.GetColumns());
        Multiply(result, lhs, rhs);
        return result;
    }



    ///////////////////////////////////////////////////////////////////
    // Addition
    ///////////////////////////////////////////////////////////////////

    template<typename DataType, typename RhsDataType, typename RhsMatrixType>
    void AddEqual(NekMatrix<DataType, StandardMatrixTag>& result,
                     const NekMatrix<RhsDataType, RhsMatrixType>& rhs);


    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void Add(NekMatrix<DataType, StandardMatrixTag>& result,
                const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, RhsMatrixType>& rhs);
    
            
    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType, StandardMatrixTag> 
    Add(const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
        const NekMatrix<RhsDataType, RhsMatrixType>& rhs);

    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void AddNegatedLhs(NekMatrix<DataType, StandardMatrixTag>& result,
                       const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                       const NekMatrix<RhsDataType, RhsMatrixType>& rhs);

    template<typename DataType, typename RhsDataType, typename RhsMatrixType>
    void AddEqualNegatedLhs(NekMatrix<DataType, StandardMatrixTag>& result,
                     const NekMatrix<RhsDataType, RhsMatrixType>& rhs);
    

    ////////////////////////////////////////////////////////////////////////////////////
    // Subtraction
    ////////////////////////////////////////////////////////////////////////////////////
    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void Subtract(NekMatrix<DataType, StandardMatrixTag>& result,
                const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, RhsMatrixType>& rhs);

    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void SubtractNegatedLhs(NekMatrix<DataType, StandardMatrixTag>& result,
                const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, RhsMatrixType>& rhs);

    template<typename DataType, typename RhsDataType, typename RhsMatrixType>
    void SubtractEqual(NekMatrix<DataType, StandardMatrixTag>& result,
                          const NekMatrix<RhsDataType, RhsMatrixType>& rhs);

    template<typename DataType, typename RhsDataType, typename RhsMatrixType>
    void SubtractEqualNegatedLhs(NekMatrix<DataType, StandardMatrixTag>& result,
                          const NekMatrix<RhsDataType, RhsMatrixType>& rhs);
    
    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType, StandardMatrixTag> 
    Subtract(const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, RhsMatrixType>& rhs);
}


#endif 

