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
#include <ExpressionTemplates/Operators.hpp>
#include <LibUtilities/LinearAlgebra/IsDgemmTraits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>

#include <string>

namespace Nektar
{
    
        
    ////////////////////////////////////////////////////////////////////////////////////
    // Matrix-Constant Multiplication
    ////////////////////////////////////////////////////////////////////////////////////
    template<typename ResultDataType, typename LhsDataType, typename LhsMatrixType>
    void Multiply(NekMatrix<ResultDataType, StandardMatrixTag>& result,
                     const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                     const ResultDataType& rhs)
    {
        // TODO - optimize for the different matrix types.
        for(unsigned int i = 0; i < lhs.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < lhs.GetColumns(); ++j)
            {
                result(i,j) = lhs(i,j)*rhs;
            }
        }
    }
    
    template<typename DataType, typename LhsDataType, typename LhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType, StandardMatrixTag> 
    Multiply(const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                const DataType& rhs)
    {
        typedef typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType ResultDataType;
        NekMatrix<ResultDataType, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
        Multiply(result, lhs, rhs);
        return result;
    }

    template<typename RhsDataType, typename RhsMatrixType, typename ResultDataType>
    void Multiply(NekMatrix<ResultDataType, StandardMatrixTag>& result,
                     const ResultDataType& lhs,
                     const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
                     
    {
        Multiply(result, rhs, lhs);
    }
    
    template<typename DataType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<RhsDataType, RhsMatrixType>::NumberType, StandardMatrixTag> 
    Multiply(const DataType& lhs,
                const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
                     
    {
        return Multiply(rhs, lhs);
    }
    
    template<typename LhsDataType>
    void MultiplyEqual(NekMatrix<LhsDataType, StandardMatrixTag>& lhs,
                          typename boost::call_traits<LhsDataType>::const_reference rhs)
    {
        typedef typename NekMatrix<LhsDataType, StandardMatrixTag>::iterator iterator;
        for( iterator iter = lhs.begin(); iter != lhs.end(); ++iter)
        {
            *iter *= rhs;
        }
    }
    
    
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
            result.GetRawPtr(), lhs.GetRows());
    }
    
    template<typename ResultType, typename LhsDataType, typename RhsDataType,
             typename LhsMatrixType, typename RhsMatrixType>
    void NekMultiplyDefaultImpl(NekMatrix<ResultType, StandardMatrixTag>& result,
                                         const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                                         const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + 
            std::string(" and a right side matrix with row count ") + 
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

        for(unsigned int i = 0; i < result.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < result.GetColumns(); ++j)
            {
                ResultType t = ResultType(0);

                // Set the result(i,j) element.
                for(unsigned int k = 0; k < lhs.GetColumns(); ++k)
                {
                    t += lhs(i,k)*rhs(k,j);
                }
                result(i,j) = t;
            }
        }
    }
     
    template<typename ResultType, typename LhsDataType, typename RhsDataType,
             typename LhsMatrixType, typename RhsMatrixType>
    void NekMultiplyFullMatrixFullMatrix(NekMatrix<ResultType, StandardMatrixTag>& result,
                                         const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                                         const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        NekMultiplyDefaultImpl(result, lhs, rhs);
    }
    
    template<typename LhsDataType, typename RhsDataType, typename DataType,
             typename LhsMatrixType, typename RhsMatrixType>
    void Multiply(NekMatrix<DataType, StandardMatrixTag>& result,
                     const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                     const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + 
            std::string(" and a right side matrix with row count ") + 
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

        result.SetSize(lhs.GetRows(), rhs.GetColumns());
        if( lhs.GetType() == eFULL && rhs.GetType() == eFULL)
        {
            NekMultiplyFullMatrixFullMatrix(result, lhs, rhs);
        }
        else
        {
            NekMultiplyDefaultImpl(result, lhs, rhs);
        }
    }
        
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
    template<typename ResultType, typename RhsType>
    void AddEqualMatricesWithSameStorageType(ResultType& result, const RhsType& rhs)
    {
        //typename ResultType::NumberType* r = result.GetRawPtr();
        //const typename RhsType::NumberType* rhs_buf = rhs.GetRawPtr();
        //unsigned int numElements = result.GetStorageSize();
        //for(unsigned int i = 0; i < numElements; ++i)
        //{
        //    r[i] += rhs_buf[i];
        //}
        typename ResultType::iterator result_iter = result.begin();
        typename RhsType::const_iterator rhs_iter = rhs.begin();
        while( result_iter != result.end() )
        {
            (*result_iter) += (*rhs_iter);
            ++result_iter;
            ++rhs_iter;
        }
    }

    
    template<typename ResultType, typename LhsType, typename RhsType>
    void AddMatricesWithSameStorageType(ResultType& result, const LhsType& lhs, const RhsType& rhs)
    {
        //if( lhs.GetTransposeFlag() == rhs.GetTransposeFlag() )
        //{
        //    // This is much more efficient if we don't need to worry about stepping through 
        //    // the data in different ways for each operand.
        //    typename ResultType::NumberType* r = result.GetRawPtr();
        //    const typename LhsType::NumberType* lhs_buf = lhs.GetRawPtr();
        //    const typename RhsType::NumberType* rhs_buf = rhs.GetRawPtr();
        //    unsigned int end = result.GetStorageSize();
        //    for(unsigned int i = 0; i < end; ++i)
        //    {
        //        r[i] = lhs_buf[i] + rhs_buf[i];
        //    }
        //}
        //else
        //{
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
        //}
    }


    template<typename DataType, typename RhsDataType, typename RhsMatrixType>
    void AddEqual(NekMatrix<DataType, StandardMatrixTag>& result,
                     const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        ASSERTL1(result.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(result.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL1(result.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(result.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));

        //return AddEqualMatricesWithSameStorageType(result, rhs);
        for(unsigned int i = 0; i < rhs.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < rhs.GetColumns(); ++j)
            {
                result(i,j) += rhs(i,j);
            }
        }
    }



    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void Add(NekMatrix<DataType, StandardMatrixTag>& result,
                const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));

        //AddMatricesWithSameStorageType(result, lhs, rhs);
        for(unsigned int i = 0; i < lhs.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < lhs.GetColumns(); ++j)
            {
                result(i,j) = lhs(i,j) + rhs(i,j);
            }
        }
    }

    
            
    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType, StandardMatrixTag> 
    Add(const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
           const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        typedef typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType NumberType;
        NekMatrix<NumberType, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
        Add(result, lhs, rhs);
        return result;
    }
    
//    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    void Add(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
//                        const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
//                        const NekMatrix<RhsDataType, DiagonalMatrixTag, RhsMatrixType>& rhs)
//    {
//        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
//            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
//            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
//        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
//            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
//            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));
//        
//        result = lhs;
//        for(unsigned int i = 0; i < result.GetRows(); ++i)
//        {
//            result(i,i) += rhs(i,i);
//        }
//    }
    
//    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    NekMatrix<typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType, FullMatrixTag, StandardMatrixTag> 
//    Add(const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
//           const NekMatrix<RhsDataType, DiagonalMatrixTag, RhsMatrixType>& rhs)
//    {
//        typedef typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType NumberType;
//        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
//        Add(result, lhs, rhs);
//        return result;
//    }
//    
//    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    void Add(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
//                        const NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>& lhs,
//                        const NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>& rhs)
//    {
//        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
//            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
//            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
//        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
//            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
//            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));
//
//        result = rhs;
//        for(unsigned int i = 0; i < result.GetRows(); ++i)
//        {
//            result(i,i) += lhs(i,i);
//        }
//    }
//    
//    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    NekMatrix<typename NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>::NumberType, FullMatrixTag, StandardMatrixTag> 
//    Add(const NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>& lhs,
//           const NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>& rhs)
//    {
//        typedef typename NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>::NumberType NumberType;
//        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
//        Add(result, lhs, rhs);
//        return result;
//    }
//    
//    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    void Add(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
//                const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
//                const NekMatrix<RhsDataType, UpperTriangularMatrixTag, RhsMatrixType>& rhs)
//    {
//        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
//            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
//            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
//        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
//            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
//            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));
//
//        result = lhs;
//        for(unsigned int i = 0; i < result.GetRows(); ++i)
//        {
//            for(unsigned int j = i; j < result.GetColumns(); ++j)
//            {
//                result(i,i) += lhs(i,i);
//            }
//        }
//    }
//    
//    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    NekMatrix<typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType, FullMatrixTag, StandardMatrixTag> 
//    Add(const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
//           const NekMatrix<RhsDataType, UpperTriangularMatrixTag, RhsMatrixType>& rhs)
//    {
//        typedef typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType NumberType;
//        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
//        Add(result, lhs, rhs);
//        return result;
//    }
//
//    // UT = UT + UT
//    // Line 22
//    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    void Add(NekMatrix<DataType, UpperTriangularMatrixTag, StandardMatrixTag>& result,
//                const NekMatrix<LhsDataType, UpperTriangularMatrixTag, LhsMatrixType>& lhs,
//                const NekMatrix<RhsDataType, UpperTriangularMatrixTag, RhsMatrixType>& rhs)
//    {
//        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
//            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
//            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
//        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
//            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
//            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));
//
//        return AddMatricesWithSameStorageType(result, lhs, rhs);
//    }
//
//    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    NekMatrix<typename NekMatrix<LhsDataType, UpperTriangularMatrixTag, LhsMatrixType>::NumberType, UpperTriangularMatrixTag, StandardMatrixTag> 
//    Add(const NekMatrix<LhsDataType, UpperTriangularMatrixTag, LhsMatrixType>& lhs,
//           const NekMatrix<RhsDataType, UpperTriangularMatrixTag, RhsMatrixType>& rhs)
//    {
//        typedef typename NekMatrix<LhsDataType, UpperTriangularMatrixTag, LhsMatrixType>::NumberType NumberType;
//        NekMatrix<NumberType, UpperTriangularMatrixTag, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
//        Add(result, lhs, rhs);
//        return result;
//    }
//
//    // LT = LT + LT
//    // Line 32
//    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    void Add(NekMatrix<DataType, LowerTriangularMatrixTag, StandardMatrixTag>& result,
//                const NekMatrix<LhsDataType, LowerTriangularMatrixTag, LhsMatrixType>& lhs,
//                const NekMatrix<RhsDataType, LowerTriangularMatrixTag, RhsMatrixType>& rhs)
//    {
//        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
//            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
//            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be added."));
//        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
//            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
//            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be added."));
//            
//        return AddMatricesWithSameStorageType(result, lhs, rhs);
//    }
//    
//    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    NekMatrix<typename NekMatrix<LhsDataType, LowerTriangularMatrixTag, LhsMatrixType>::NumberType, LowerTriangularMatrixTag, StandardMatrixTag> 
//    Add(const NekMatrix<LhsDataType, LowerTriangularMatrixTag, LhsMatrixType>& lhs,
//           const NekMatrix<RhsDataType, LowerTriangularMatrixTag, RhsMatrixType>& rhs)
//    {
//        typedef typename NekMatrix<LhsDataType, LowerTriangularMatrixTag, LhsMatrixType>::NumberType NumberType;
//        NekMatrix<NumberType, LowerTriangularMatrixTag, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
//        Add(result, lhs, rhs);
//        return result;
//    }

    ////////////////////////////////////////////////////////////////////////////////////
    // Subtraction
    ////////////////////////////////////////////////////////////////////////////////////
    template<typename ResultType, typename RhsType>
    void SubtractEqualMatricesWithSameStorageType(ResultType& result, const RhsType& rhs)
    {
        typename ResultType::iterator result_iter = result.begin();
        typename RhsType::const_iterator rhs_iter = rhs.begin();
        while( result_iter != result.end() )
        {
            (*result_iter) -= (*rhs_iter);
            ++result_iter;
            ++rhs_iter;
        }
    }


    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void Subtract(NekMatrix<DataType, StandardMatrixTag>& result,
                const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be subtracted."));
        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be subtracted."));

//        typename NekMatrix<DataType, StandardMatrixTag>::iterator result_iter = result.begin();
//        typename NekMatrix<RhsDataType, RhsMatrixType>::const_iterator rhs_iter = rhs.begin();
//        typename NekMatrix<LhsDataType, LhsMatrixType>::const_iterator lhs_iter = lhs.begin();
//        while( result_iter != result.end() )
//        {
//            (*result_iter) = (*lhs_iter) - (*rhs_iter);
//            ++result_iter;
//            ++rhs_iter;
//            ++lhs_iter;
//        }
        for(unsigned int i = 0; i < lhs.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < lhs.GetColumns(); ++j)
            {
                result(i,j) = lhs(i,j) - rhs(i,j);
            }
        }
    }

    template<typename DataType, typename RhsDataType, typename RhsMatrixType>
    void SubtractEqual(NekMatrix<DataType, StandardMatrixTag>& result,
                          const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        ASSERTL1(result.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
            boost::lexical_cast<std::string>(result.GetRows()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be subtracted."));
        ASSERTL1(result.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
            boost::lexical_cast<std::string>(result.GetColumns()) + std::string(" and ") +
            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be subtracted."));

        //SubtractEqualMatricesWithSameStorageType(result, rhs);
        for(unsigned int i = 0; i < rhs.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < rhs.GetColumns(); ++j)
            {
                result(i,j) -= rhs(i,j);
            }
        }
    }
    
    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType, StandardMatrixTag> 
    Subtract(const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        typedef typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType NumberType;
        NekMatrix<NumberType, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
        Subtract(result, lhs, rhs);
        return result;
    }
    
//    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    void Subtract(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
//                        const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
//                        const NekMatrix<RhsDataType, DiagonalMatrixTag, RhsMatrixType>& rhs)
//    {
//        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
//            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
//            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be subtracted."));
//        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
//            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
//            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be subtracted."));
//        
//        result = lhs;
//        for(unsigned int i = 0; i < result.GetRows(); ++i)
//        {
//            result(i,i) -= rhs(i,i);
//        }
//    }
//    
//    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    NekMatrix<typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType, FullMatrixTag, StandardMatrixTag> 
//    Subtract(const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
//                const NekMatrix<RhsDataType, DiagonalMatrixTag, RhsMatrixType>& rhs)
//    {
//        typedef typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType NumberType;
//        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
//        Subtract(result, lhs, rhs);
//        return result;
//    }
//
//    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    void Subtract(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
//                        const NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>& lhs,
//                        const NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>& rhs)
//    {
//        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
//            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
//            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be subtracted."));
//        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
//            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
//            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be subtracted."));
//        
//        result = rhs;
//        for(unsigned int i = 0; i < result.GetRows(); ++i)
//        {
//            result(i,i) -= lhs(i,i);
//        }
//    }
//    
//    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    NekMatrix<typename NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>::NumberType, FullMatrixTag, StandardMatrixTag> 
//    Subtract(const NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>& lhs,
//                const NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>& rhs)
//    {
//        typedef typename NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>::NumberType NumberType;
//        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(rhs.GetRows(), rhs.GetColuomns());
//        Subtract(result, lhs, rhs);
//        return result;
//    }
//    
//    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    void Subtract(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
//                        const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
//                        const NekMatrix<RhsDataType, UpperTriangularMatrixTag, RhsMatrixType>& rhs)
//    {
//        ASSERTL1(lhs.GetRows() == rhs.GetRows(), std::string("Matrices with different row counts  ") + 
//            boost::lexical_cast<std::string>(lhs.GetRows()) + std::string(" and ") +
//            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be subtracted."));
//        ASSERTL1(lhs.GetColumns() == rhs.GetColumns(), std::string("Matrices with different column counts  ") + 
//            boost::lexical_cast<std::string>(lhs.GetColumns()) + std::string(" and ") +
//            boost::lexical_cast<std::string>(rhs.GetColumns()) + std::string(" can't be subtracted."));
//
//        result = rhs;
//        for(unsigned int i = 0; i < result.GetRows(); ++i)
//        {
//            for(unsigned int j = i; j < result.GetColumns(); ++j)
//            {
//                result(i,i) -= lhs(i,i);
//            }
//        }
//    }
//    
//    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    NekMatrix<typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType, FullMatrixTag, StandardMatrixTag> 
//    Subtract(const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
//                const NekMatrix<RhsDataType, UpperTriangularMatrixTag, RhsMatrixType>& rhs)
//    {
//        typedef typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType NumberType;
//        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(rhs.GetRows(), rhs.GetColumns());
//        Subtract(result, lhs, rhs);
//        return result;
//    }



    template<typename ADataType, typename BDataType, typename CDataType, 
             typename AMatrixType, typename BMatrixType, typename CMatrixType>
    void Dgemm(NekMatrix<double, StandardMatrixTag>& result, 
               double alpha, const NekMatrix<ADataType, AMatrixType>& A, const NekMatrix<BDataType, BMatrixType>& B,
               double beta, const NekMatrix<CDataType, CMatrixType>& C)
    {
        result = C;
        unsigned int M = A.GetRows();
        unsigned int N = B.GetColumns();
        unsigned int K = A.GetColumns();

        unsigned int LDA = M;
        if( A.GetTransposeFlag() == 'T' )
        {
            LDA = K;
        }

        unsigned int LDB = K;
        if( B.GetTransposeFlag() == 'T' )
        {
            LDB = N;
        }

        Blas::Dgemm(A.GetTransposeFlag(), B.GetTransposeFlag(), M, N, K,
            A.Scale()*B.Scale()*alpha, A.GetRawPtr(), LDA,
            B.GetRawPtr(), LDB, beta*result.Scale(),
            result.GetRawPtr(), M);
    }

    template<typename ADataType, typename BDataType, typename CDataType, 
             typename AMatrixType, typename BMatrixType, typename CMatrixType>
    void Dgemm(NekMatrix<double, StandardMatrixTag>& result, 
               const NekMatrix<ADataType, AMatrixType>& A, double alpha, const NekMatrix<BDataType, BMatrixType>& B,
               double beta, const NekMatrix<CDataType, CMatrixType>& C)
    {
        Dgemm(result, alpha, A, B, beta, C);
    }

    template<typename ADataType, typename BDataType, typename CDataType, 
             typename AMatrixType, typename BMatrixType, typename CMatrixType>
    void Dgemm(NekMatrix<double, StandardMatrixTag>& result, 
               const NekMatrix<ADataType, AMatrixType>& A, const NekMatrix<BDataType, BMatrixType>& B, double alpha, 
               double beta, const NekMatrix<CDataType, CMatrixType>& C)
    {
        Dgemm(result, alpha, A, B, beta, C);
    }

    template<typename L, typename R, typename IndicesType, unsigned int index>
    struct BinaryBinaryEvaluateNodeOverride<L, AddOp, R, IndicesType, index,
        typename boost::enable_if
        <
            boost::mpl::and_
            <
                IsDgemmLeftSide<L>,
                IsDgemmRightSide<R>
            >
        >::type > : public boost::true_type 
    {
        template<typename ResultType, typename ArgumentVectorType>
        static void Evaluate(ResultType& accumulator, const ArgumentVectorType& args)
        {
            Dgemm(accumulator, 
                IsDgemmLeftSide<L>::GetValues<IndicesType, index>::GetAlpha(args),
                IsDgemmLeftSide<L>::GetValues<IndicesType, index>::GetA(args),
                IsDgemmLeftSide<L>::GetValues<IndicesType, index>::GetB(args),
                IsDgemmRightSide<R>::GetValues<IndicesType, index + L::TotalCount>::GetBeta(args),
                IsDgemmRightSide<R>::GetValues<IndicesType, index + L::TotalCount>::GetC(args) );
        }
    };


    // TODO - This case also catches chained multipliation, which we don't want.
    // But we do want it to catch a*A*B.
    //template<typename L, typename R, typename IndicesType, unsigned int index>
    //struct BinaryBinaryEvaluateNodeOverride<L, MultiplyOp, R, IndicesType, index,
    //    typename boost::enable_if
    //    <
    //        IsDgemmLeftSide<Node<L, MultiplyOp, R> >
    //    >::type > : public boost::true_type 
    //{
    //    typedef IsDgemmLeftSide<Node<L, MultiplyOp, R> > Accessor;

    //    template<typename ResultType, typename ArgumentVectorType>
    //    static void Evaluate(ResultType& accumulator, const ArgumentVectorType& args)
    //    {
    //        Dgemm(accumulator, 
    //            Accessor::GetValues<IndicesType, index>::GetAlpha(args),
    //            Accessor::GetValues<IndicesType, index>::GetA(args),
    //            Accessor::GetValues<IndicesType, index>::GetB(args),
    //            0.0, accumulator);
    //    }
    //};

    //    
    //    namespace Impl
    //    {
    //        template<typename LhsExpressionPolicyType, typename RhsExpressionPolicyType, typename ResultType, 
    //                 template <typename, typename> class OpType, 
    //                 template <typename, typename> class ParentOpType = BinaryNullOp,
    //                 typename enabled = void>
    //        struct DgemmBinaryExpressionEvaluator;
    //                            
    //        // AB + C
    //        template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //        struct DgemmBinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                                            MultiplyOp, 
    //                                            ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                        >,
    //                                        ConstantExpressionPolicy<NekMatrix<T3, M3> >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        AddOp, BinaryNullOp,
    //                                        typename boost::enable_if
    //                                        <
    //                                            boost::mpl::and_
    //                                            <
    //                                                CanGetRawPtr<NekMatrix<T1, M1> >,
    //                                                CanGetRawPtr<NekMatrix<T2, M2> >,
    //                                                CanGetRawPtr<NekMatrix<T3, M3> >
    //                                            >
    //                                        >::type>
    //        {
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                        MultiplyOp, 
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    > LhsType;
    //            typedef ConstantExpressionPolicy<NekMatrix<T3, M3> > RhsType;
    //            
    //            static void Eval(const Expression<LhsType>& lhs, 
    //                             const Expression<RhsType>& rhs,
    //                             Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //            {
    //                const NekMatrix<T1, M1>& a = *LhsType::Left(*lhs);
    //                const NekMatrix<T2, M2>& b = *LhsType::Right(*lhs);
    //                const NekMatrix<T3, M3>& c = *rhs;
    //                
    //                if( a.GetType() == eFULL && b.GetType() == eFULL && c.GetType() == eFULL )
    //                {
    //                    Dgemm(*result, 1.0, a, b, 1.0, c);
    //                }
    //                else
    //                {
    //                    typedef typename LhsType::ResultType LhsResultType;
    //                    typedef typename RhsType::ResultType RhsResultType;
    //                    lhs.Evaluate(result);
    //                    AddOp<LhsResultType, RhsResultType>::ApplyEqual(result, *rhs);
    //                }
    //            }
    //        };
    //        
    //        // AB + C - in case there is a block matrix in the expression.
    //        template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //        struct DgemmBinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                                            MultiplyOp, 
    //                                            ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                        >,
    //                                        ConstantExpressionPolicy<NekMatrix<T3, M3> >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        AddOp, BinaryNullOp,
    //                                        typename boost::enable_if
    //                                        <
    //                                            boost::mpl::or_
    //                                            <
    //                                                boost::mpl::not_<CanGetRawPtr<NekMatrix<T1, M1> > >,
    //                                                boost::mpl::not_<CanGetRawPtr<NekMatrix<T2, M2> > >,
    //                                                boost::mpl::not_<CanGetRawPtr<NekMatrix<T3, M3> > >
    //                                            >
    //                                        >::type>
    //        {
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                        MultiplyOp, 
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    > LhsType;
    //            typedef ConstantExpressionPolicy<NekMatrix<T3, M3> > RhsType;
    //            
    //            static void Eval(const Expression<LhsType>& lhs, 
    //                             const Expression<RhsType>& rhs,
    //                             Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //            {
    //                typedef typename LhsType::ResultType LhsResultType;
    //                typedef typename RhsType::ResultType RhsResultType;
    //                lhs.Evaluate(result);
    //                AddOp<LhsResultType, RhsResultType>::ApplyEqual(result, *rhs);
    //            }
    //        };


    //        // aAB + C
    //        template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //        struct DgemmBinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                        <
    //                                            BinaryExpressionPolicy
    //                                            <
    //                                                ConstantExpressionPolicy<double>,
    //                                                MultiplyOp,
    //                                                ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                                            >,
    //                                            MultiplyOp, 
    //                                            ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                        >,
    //                                        ConstantExpressionPolicy<NekMatrix<T3, M3> >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        AddOp, BinaryNullOp,
    //                                        typename boost::enable_if
    //                                        <
    //                                            boost::mpl::and_
    //                                            <
    //                                                CanGetRawPtr<NekMatrix<T1, M1> >,
    //                                                CanGetRawPtr<NekMatrix<T2, M2> >,
    //                                                CanGetRawPtr<NekMatrix<T3, M3> >
    //                                            >
    //                                        >::type>
    //        {
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                    > LhsLhsType;
    //                    
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        LhsLhsType,
    //                        MultiplyOp, 
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    > LhsType;
    //            typedef ConstantExpressionPolicy<NekMatrix<T3, M3> > RhsType;
    //            
    //            static void Eval(const Expression<LhsType>& lhs, 
    //                             const Expression<RhsType>& rhs,
    //                             Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //            {
    //                const Expression<LhsLhsType>& lhsExpression = LhsType::Left(*lhs);
    //                double alpha = *LhsLhsType::Left(*lhsExpression);
    //                const NekMatrix<T1, M1>& a = *LhsLhsType::Right(*lhsExpression);
    //                const NekMatrix<T2, M2>& b = *LhsType::Right(*lhs);
    //                const NekMatrix<T3, M3>& c = *rhs;
    //                
    //                if( a.GetType() == eFULL && b.GetType() == eFULL && c.GetType() == eFULL )
    //                {
    //                    Dgemm(*result, alpha, a, b, 1.0, c);
    //                }
    //                else
    //                {
    //                    typedef typename LhsType::ResultType LhsResultType;
    //                    typedef typename RhsType::ResultType RhsResultType;
    //                    lhs.Evaluate(result);
    //                    AddOp<LhsResultType, RhsResultType>::ApplyEqual(result, *rhs);
    //                }
    //            }
    //        };
    //        
    //        // aAB + C - in case there is a block matrix in the expression.
    //        template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //        struct DgemmBinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                        <
    //                                            BinaryExpressionPolicy
    //                                            <
    //                                                ConstantExpressionPolicy<double>,
    //                                                MultiplyOp,
    //                                                ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                                            >,
    //                                            MultiplyOp, 
    //                                            ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                        >,
    //                                        ConstantExpressionPolicy<NekMatrix<T3, M3> >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        AddOp, BinaryNullOp,
    //                                        typename boost::enable_if
    //                                        <
    //                                            boost::mpl::or_
    //                                            <
    //                                                boost::mpl::not_<CanGetRawPtr<NekMatrix<T1, M1> > >,
    //                                                boost::mpl::not_<CanGetRawPtr<NekMatrix<T2, M2> > >,
    //                                                boost::mpl::not_<CanGetRawPtr<NekMatrix<T3, M3> > >
    //                                            >
    //                                        >::type>
    //        {
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                    > LhsLhsType;
    //                    
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        LhsLhsType,
    //                        MultiplyOp, 
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    > LhsType;
    //            typedef ConstantExpressionPolicy<NekMatrix<T3, M3> > RhsType;
    //            
    //            static void Eval(const Expression<LhsType>& lhs, 
    //                             const Expression<RhsType>& rhs,
    //                             Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //            {
    //                typedef typename LhsType::ResultType LhsResultType;
    //                typedef typename RhsType::ResultType RhsResultType;
    //                lhs.Evaluate(result);
    //                AddOp<LhsResultType, RhsResultType>::ApplyEqual(result, *rhs);
    //            }
    //        };


    //        // AB + bC
    //        template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //        struct DgemmBinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                                            MultiplyOp, 
    //                                            ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                        >,
    //                                        BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<double>,
    //                                            MultiplyOp,
    //                                            ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                                        >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        AddOp, BinaryNullOp,
    //                                        typename boost::enable_if
    //                                        <
    //                                            boost::mpl::and_
    //                                            <
    //                                                CanGetRawPtr<NekMatrix<T1, M1> >,
    //                                                CanGetRawPtr<NekMatrix<T2, M2> >,
    //                                                CanGetRawPtr<NekMatrix<T3, M3> >
    //                                            >
    //                                        >::type>
    //        {
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                        MultiplyOp, 
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    > LhsType;
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                    > RhsType;
    //            
    //            static void Eval(const Expression<LhsType>& lhs, 
    //                             const Expression<RhsType>& rhs,
    //                             Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //            {
    //                const NekMatrix<T1, M1>& a = *LhsType::Left(*lhs);
    //                const NekMatrix<T2, M2>& b = *LhsType::Right(*lhs);
    //                double beta = *RhsType::Left(*rhs);
    //                const NekMatrix<T3, M3>& c = *RhsType::Right(*rhs);
    //                
    //                if( a.GetType() == eFULL && b.GetType() == eFULL && c.GetType() == eFULL )
    //                {
    //                    Dgemm(*result, 1.0, a, b, beta, c);
    //                }
    //                else
    //                {
    //                    lhs.Evaluate(result);
    //                    NekMatrix<double> rhsTemp = rhs.Evaluate();
    //                    AddOp<NekMatrix<double>, NekMatrix<double> >::ApplyEqual(result, rhsTemp);
    //                }
    //            }
    //        };
    //        
    //        // AB + bC - in case there is a block matrix in the expression.
    //        template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //        struct DgemmBinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                                            MultiplyOp, 
    //                                            ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                        >,
    //                                        BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<double>,
    //                                            MultiplyOp,
    //                                            ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                                        >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        AddOp, BinaryNullOp,
    //                                        typename boost::enable_if
    //                                        <
    //                                            boost::mpl::or_
    //                                            <
    //                                                boost::mpl::not_<CanGetRawPtr<NekMatrix<T1, M1> > >,
    //                                                boost::mpl::not_<CanGetRawPtr<NekMatrix<T2, M2> > >,
    //                                                boost::mpl::not_<CanGetRawPtr<NekMatrix<T3, M3> > >
    //                                            >
    //                                        >::type>
    //        {
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                        MultiplyOp, 
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    > LhsType;
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                    > RhsType;
    //            
    //            static void Eval(const Expression<LhsType>& lhs, 
    //                             const Expression<RhsType>& rhs,
    //                             Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //            {
    //                lhs.Evaluate(result);
    //                NekMatrix<double> rhsTemp = rhs.Evaluate();
    //                AddOp<NekMatrix<double>, NekMatrix<double> >::ApplyEqual(result, rhsTemp);
    //            }
    //        };
    //        
    //        // aAB + bC
    //        template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //        struct DgemmBinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                        <
    //                                            BinaryExpressionPolicy
    //                                            <
    //                                                ConstantExpressionPolicy<double>,
    //                                                MultiplyOp,
    //                                                ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                                            >,
    //                                            MultiplyOp,
    //                                            ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                        >,
    //                                        BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<double>,
    //                                            MultiplyOp,
    //                                            ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                                        >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        AddOp, BinaryNullOp,
    //                                        typename boost::enable_if
    //                                        <
    //                                            boost::mpl::and_
    //                                            <
    //                                                CanGetRawPtr<NekMatrix<T1, M1> >,
    //                                                CanGetRawPtr<NekMatrix<T2, M2> >,
    //                                                CanGetRawPtr<NekMatrix<T3, M3> >
    //                                            >
    //                                        >::type>
    //        {
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                    > LhsLhsType;
    //                    
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        LhsLhsType,
    //                        MultiplyOp, 
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    > LhsType;
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                    > RhsType;
    //            
    //            static void Eval(const Expression<LhsType>& lhs, 
    //                             const Expression<RhsType>& rhs,
    //                             Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //            {
    //                const Expression<LhsLhsType>& lhsExpression = LhsType::Left(*lhs);
    //                double alpha = *LhsLhsType::Left(*lhsExpression);
    //                const NekMatrix<T1, M1>& a = *LhsLhsType::Right(*lhsExpression);
    //                const NekMatrix<T2, M2>& b = *LhsType::Right(*lhs);
    //                const double beta = *RhsType::Left(*rhs);
    //                const NekMatrix<T3, M3>& c = *RhsType::Right(*rhs);
    //                
    //                if( a.GetType() == eFULL && b.GetType() == eFULL && c.GetType() == eFULL )
    //                {
    //                    Dgemm(*result, alpha, a, b, beta, c);
    //                }
    //                else
    //                {
    //                    typedef typename LhsType::ResultType LhsResultType;
    //                    typedef typename RhsType::ResultType RhsResultType;
    //                    lhs.Evaluate(result);
    //                    NekMatrix<double> rhsTemp = rhs.Evaluate();
    //                    AddOp<LhsResultType, RhsResultType>::ApplyEqual(result, rhsTemp);
    //                }
    //            }
    //        };
    //        
    //        // aAB + bC - in case there is a block matrix in the expression.
    //        template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //        struct DgemmBinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                        <
    //                                            BinaryExpressionPolicy
    //                                            <
    //                                                ConstantExpressionPolicy<double>,
    //                                                MultiplyOp,
    //                                                ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                                            >,
    //                                            MultiplyOp,
    //                                            ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                        >,
    //                                        BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<double>,
    //                                            MultiplyOp,
    //                                            ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                                        >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        AddOp, BinaryNullOp,
    //                                        typename boost::enable_if
    //                                        <
    //                                            boost::mpl::or_
    //                                            <
    //                                                boost::mpl::not_<CanGetRawPtr<NekMatrix<T1, M1> > >,
    //                                                boost::mpl::not_<CanGetRawPtr<NekMatrix<T2, M2> > >,
    //                                                boost::mpl::not_<CanGetRawPtr<NekMatrix<T3, M3> > >
    //                                            >
    //                                        >::type>
    //        {
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                    > LhsLhsType;
    //                    
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        LhsLhsType,
    //                        MultiplyOp, 
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    > LhsType;
    //             typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                    > RhsType;
    //            
    //            static void Eval(const Expression<LhsType>& lhs, 
    //                             const Expression<RhsType>& rhs,
    //                             Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //            {
    //                typedef typename LhsType::ResultType LhsResultType;
    //                typedef typename RhsType::ResultType RhsResultType;
    //                lhs.Evaluate(result);
    //                NekMatrix<double> rhsTemp = rhs.Evaluate();
    //                AddOp<LhsResultType, RhsResultType>::ApplyEqual(result, rhsTemp);
    //            }
    //        };
    //        
    //        
    //        // a(AB) + bC
    //        template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //        struct DgemmBinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<double>,
    //                                            MultiplyOp,
    //                                            BinaryExpressionPolicy
    //                                            <
    //                                                ConstantExpressionPolicy<NekMatrix<T1, M1> >,
    //                                                MultiplyOp,
    //                                                ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                            >
    //                                        >,
    //                                        BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<double>,
    //                                            MultiplyOp,
    //                                            ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                                        >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        AddOp, BinaryNullOp,
    //                                        typename boost::enable_if
    //                                        <
    //                                            boost::mpl::and_
    //                                            <
    //                                                CanGetRawPtr<NekMatrix<T1, M1> >,
    //                                                CanGetRawPtr<NekMatrix<T2, M2> >,
    //                                                CanGetRawPtr<NekMatrix<T3, M3> >
    //                                            >
    //                                        >::type>
    //        {
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >,
    //                        MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    > LhsRhsType;
    //                    
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        MultiplyOp, 
    //                        LhsRhsType
    //                    > LhsType;
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                    > RhsType;
    //            
    //            static void Eval(const Expression<LhsType>& lhs, 
    //                             const Expression<RhsType>& rhs,
    //                             Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //            {
    //                double alpha = *LhsType::Left(*lhs);
    //                const Expression<LhsRhsType>& rhsExpression = LhsType::Right(*lhs);
    //                const NekMatrix<T1, M1>& a = *LhsRhsType::Left(*rhsExpression);
    //                const NekMatrix<T2, M2>& b = *LhsRhsType::Right(*rhsExpression);
    //                const double beta = *RhsType::Left(*rhs);
    //                const NekMatrix<T3, M3>& c = *RhsType::Right(*rhs);
    //                
    //                if( a.GetType() == eFULL && b.GetType() == eFULL && c.GetType() == eFULL )
    //                {
    //                    Dgemm(*result, alpha, a, b, beta, c);
    //                }
    //                else
    //                {
    //                    typedef typename LhsType::ResultType LhsResultType;
    //                    typedef typename RhsType::ResultType RhsResultType;
    //                    lhs.Evaluate(result);
    //                    NekMatrix<double> rhsTemp = rhs.Evaluate();
    //                    AddOp<LhsResultType, RhsResultType>::ApplyEqual(result, rhsTemp);
    //                }
    //            }
    //        };
    //        
    //        // a(AB) + bC - in case there is a block matrix in the expression.
    //        template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //        struct DgemmBinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<double>,
    //                                            MultiplyOp,
    //                                            BinaryExpressionPolicy
    //                                            <
    //                                                ConstantExpressionPolicy<NekMatrix<T1, M1> >,
    //                                                MultiplyOp,
    //                                                ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                            >
    //                                        >,
    //                                        BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<double>,
    //                                            MultiplyOp,
    //                                            ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                                        >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        AddOp, BinaryNullOp,
    //                                        typename boost::enable_if
    //                                        <
    //                                            boost::mpl::or_
    //                                            <
    //                                                boost::mpl::not_<CanGetRawPtr<NekMatrix<T1, M1> > >,
    //                                                boost::mpl::not_<CanGetRawPtr<NekMatrix<T2, M2> > >,
    //                                                boost::mpl::not_<CanGetRawPtr<NekMatrix<T3, M3> > >
    //                                            >
    //                                        >::type>
    //        {
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >,
    //                        MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    > LhsRhsType;
    //                    
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        MultiplyOp, 
    //                        LhsRhsType
    //                    > LhsType;
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                    > RhsType;
    //            
    //            static void Eval(const Expression<LhsType>& lhs, 
    //                             const Expression<RhsType>& rhs,
    //                             Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //            {
    //                typedef typename LhsType::ResultType LhsResultType;
    //                typedef typename RhsType::ResultType RhsResultType;
    //                lhs.Evaluate(result);
    //                NekMatrix<double> rhsTemp = rhs.Evaluate();
    //                AddOp<LhsResultType, RhsResultType>::ApplyEqual(result, rhsTemp);
    //            }
    //        };
    //        
    //        // a(AB) + C
    //        template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //        struct DgemmBinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<double>,
    //                                            MultiplyOp,
    //                                            BinaryExpressionPolicy
    //                                            <
    //                                                ConstantExpressionPolicy<NekMatrix<T1, M1> >,
    //                                                MultiplyOp,
    //                                                ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                            >
    //                                        >,
    //                                        ConstantExpressionPolicy<NekMatrix<T3, M3> >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        AddOp, BinaryNullOp,
    //                                        typename boost::enable_if
    //                                        <
    //                                            boost::mpl::and_
    //                                            <
    //                                                CanGetRawPtr<NekMatrix<T1, M1> >,
    //                                                CanGetRawPtr<NekMatrix<T2, M2> >,
    //                                                CanGetRawPtr<NekMatrix<T3, M3> >
    //                                            >
    //                                        >::type>
    //        {
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >,
    //                        MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    > LhsRhsType;
    //                    
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        MultiplyOp, 
    //                        LhsRhsType
    //                    > LhsType;
    //            typedef ConstantExpressionPolicy<NekMatrix<T3, M3> > RhsType;
    //            
    //            static void Eval(const Expression<LhsType>& lhs, 
    //                             const Expression<RhsType>& rhs,
    //                             Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //            {
    //                double alpha = *LhsType::Left(*lhs);
    //                const Expression<LhsRhsType>& rhsExpression = LhsType::Right(*lhs);
    //                const NekMatrix<T1, M1>& a = *LhsRhsType::Left(*rhsExpression);
    //                const NekMatrix<T2, M2>& b = *LhsRhsType::Right(*rhsExpression);
    //                const NekMatrix<T3, M3>& c = *rhs;
    //                
    //                if( a.GetType() == eFULL && b.GetType() == eFULL && c.GetType() == eFULL )
    //                {
    //                    Dgemm(*result, alpha, a, b, 1.0, c);
    //                }
    //                else
    //                {
    //                    typedef typename LhsType::ResultType LhsResultType;
    //                    typedef typename RhsType::ResultType RhsResultType;
    //                    lhs.Evaluate(result);
    //                    const NekMatrix<T3, M3>& rhsTemp = *rhs;
    //                    AddOp<LhsResultType, RhsResultType>::ApplyEqual(result, rhsTemp);
    //                }
    //            }
    //        };
    //        
    //        // a(AB) + C - in case there is a block matrix in the expression.
    //        template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //        struct DgemmBinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<double>,
    //                                            MultiplyOp,
    //                                            BinaryExpressionPolicy
    //                                            <
    //                                                ConstantExpressionPolicy<NekMatrix<T1, M1> >,
    //                                                MultiplyOp,
    //                                                ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                            >
    //                                        >,
    //                                        ConstantExpressionPolicy<NekMatrix<T3, M3> >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        AddOp, BinaryNullOp,
    //                                        typename boost::enable_if
    //                                        <
    //                                            boost::mpl::or_
    //                                            <
    //                                                boost::mpl::not_<CanGetRawPtr<NekMatrix<T1, M1> > >,
    //                                                boost::mpl::not_<CanGetRawPtr<NekMatrix<T2, M2> > >,
    //                                                boost::mpl::not_<CanGetRawPtr<NekMatrix<T3, M3> > >
    //                                            >
    //                                        >::type>
    //        {
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >,
    //                        MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    > LhsRhsType;
    //                    
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        MultiplyOp, 
    //                        LhsRhsType
    //                    > LhsType;
    //            typedef ConstantExpressionPolicy<NekMatrix<T3, M3> > RhsType;
    //            
    //            static void Eval(const Expression<LhsType>& lhs, 
    //                             const Expression<RhsType>& rhs,
    //                             Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //            {
    //                typedef typename LhsType::ResultType LhsResultType;
    //                typedef typename RhsType::ResultType RhsResultType;
    //                lhs.Evaluate(result);
    //                const NekMatrix<T3, M3>& rhsTemp = *rhs;
    //                AddOp<LhsResultType, RhsResultType>::ApplyEqual(result, rhsTemp);
    //            }
    //        };


    //    }
    //    
    //    
    //    // AB + C
    //    template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //    struct BinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                    <
    //                                        ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                                        MultiplyOp, 
    //                                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                    >,
    //                                    ConstantExpressionPolicy<NekMatrix<T3, M3> >,
    //                                    NekMatrix<double, StandardMatrixTag>,
    //                                    AddOp, BinaryNullOp 
    //                                    >
    //    {
    //        typedef BinaryExpressionPolicy
    //                <
    //                    ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                    MultiplyOp, 
    //                    ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                > LhsType;
    //        typedef ConstantExpressionPolicy<NekMatrix<T3, M3> > RhsType;
    //        
    //        static void Eval(const Expression<LhsType>& lhs, 
    //                         const Expression<RhsType>& rhs,
    //                         Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //        {
    //            return Impl::DgemmBinaryExpressionEvaluator<LhsType, RhsType,
    //                                    NekMatrix<double, StandardMatrixTag>,
    //                                    AddOp, BinaryNullOp >::Eval(lhs, rhs, result);
    //        }
    //    };

    //    // aAB + C
    //    template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //    struct BinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                    <
    //                                        ConstantExpressionPolicy<double>,
    //                                        MultiplyOp,
    //                                        BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                                            MultiplyOp, 
    //                                            ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                        >
    //                                    >,
    //                                    ConstantExpressionPolicy<NekMatrix<T3, M3> >,
    //                                    NekMatrix<double, StandardMatrixTag>,
    //                                    AddOp, BinaryNullOp 
    //                                    >
    //    {
    //        typedef BinaryExpressionPolicy
    //                <
    //                    ConstantExpressionPolicy<double>,
    //                    MultiplyOp,
    //                    BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                        MultiplyOp, 
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    >
    //                > LhsType;
    //        typedef ConstantExpressionPolicy<NekMatrix<T3, M3> > RhsType;
    //        
    //        static void Eval(const Expression<LhsType>& lhs, 
    //                         const Expression<RhsType>& rhs,
    //                         Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //        {
    //            return Impl::DgemmBinaryExpressionEvaluator<LhsType, RhsType,
    //                                    NekMatrix<double, StandardMatrixTag>,
    //                                    AddOp, BinaryNullOp >::Eval(lhs, rhs, result);
    //        }
    //    };


    //    // aAB + C
    //    template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //    struct BinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                    <
    //                                        BinaryExpressionPolicy
    //                                        <                                            
    //                                            ConstantExpressionPolicy<double>,
    //                                            MultiplyOp,
    //                                            ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                                        >,
    //                                        MultiplyOp, 
    //                                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                    >,
    //                                    ConstantExpressionPolicy<NekMatrix<T3, M3> >,
    //                                    NekMatrix<double, StandardMatrixTag>,
    //                                    AddOp, BinaryNullOp 
    //                                    >
    //    {
    //        typedef BinaryExpressionPolicy
    //                <
    //                    BinaryExpressionPolicy
    //                    <                                            
    //                        ConstantExpressionPolicy<double>,
    //                        MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                    >,
    //                    MultiplyOp, 
    //                    ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                > LhsType;
    //        typedef ConstantExpressionPolicy<NekMatrix<T3, M3> > RhsType;
    //        
    //        static void Eval(const Expression<LhsType>& lhs, 
    //                         const Expression<RhsType>& rhs,
    //                         Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //        {
    //            return Impl::DgemmBinaryExpressionEvaluator<LhsType, RhsType,
    //                                    NekMatrix<double, StandardMatrixTag>,
    //                                    AddOp, BinaryNullOp >::Eval(lhs, rhs, result);
    //        }
    //    };
    //    
    //    // AB + bC
    //    template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //    struct BinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                    <
    //                                        ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                                        MultiplyOp, 
    //                                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                    >,
    //                                    BinaryExpressionPolicy
    //                                    <
    //                                        ConstantExpressionPolicy<double>,
    //                                        MultiplyOp,
    //                                        ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                                    >,
    //                                    NekMatrix<double, StandardMatrixTag>,
    //                                    AddOp, BinaryNullOp>
    //    {
    //        typedef BinaryExpressionPolicy
    //                <
    //                    ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                    MultiplyOp, 
    //                    ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                > LhsType;
    //        typedef BinaryExpressionPolicy
    //                <
    //                    ConstantExpressionPolicy<double>,
    //                    MultiplyOp,
    //                    ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                > RhsType;
    //        
    //        static void Eval(const Expression<LhsType>& lhs, 
    //                         const Expression<RhsType>& rhs,
    //                         Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //        {
    //            return Impl::DgemmBinaryExpressionEvaluator<LhsType, RhsType,
    //                                NekMatrix<double, StandardMatrixTag>,
    //                    AddOp, BinaryNullOp >::Eval(lhs, rhs, result);
    //        }
    //    };

    //    // aAB + bC
    //    template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //    struct BinaryExpressionEvaluator<BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                                    MultiplyOp,
    //                                        BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                                            MultiplyOp, 
    //                                            ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                        >
    //                                    >,
    //                                    BinaryExpressionPolicy
    //                                    <
    //                                        ConstantExpressionPolicy<double>,
    //                                        MultiplyOp,
    //                                        ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                                    >,
    //                                    NekMatrix<double, StandardMatrixTag>,
    //                                    AddOp, BinaryNullOp 
    //                                    >
    //    {
    //        typedef BinaryExpressionPolicy
    //                <
    //                    ConstantExpressionPolicy<double>,
    //                    MultiplyOp,
    //                    BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                        MultiplyOp, 
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    >
    //                > LhsType;
    //        typedef BinaryExpressionPolicy
    //                <
    //                    ConstantExpressionPolicy<double>,
    //                    MultiplyOp,
    //                    ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                > RhsType;
    //        
    //        static void Eval(const Expression<LhsType>& lhs, 
    //                         const Expression<RhsType>& rhs,
    //                         Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //        {
    //            return Impl::DgemmBinaryExpressionEvaluator<LhsType, RhsType,
    //                                    NekMatrix<double, StandardMatrixTag>,
    //                                    AddOp, BinaryNullOp >::Eval(lhs, rhs, result);
    //        }
    //    };

    //    // (aA)B + bC
    //    template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //    struct BinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                    <
    //                                        BinaryExpressionPolicy
    //                                        <                                            
    //                                            ConstantExpressionPolicy<double>,
    //                                            MultiplyOp,
    //                                            ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                                        >,
    //                                        MultiplyOp, 
    //                                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                    >,
    //                                    BinaryExpressionPolicy
    //                                    <
    //                                        ConstantExpressionPolicy<double>,
    //                                        MultiplyOp,
    //                                        ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                                    >,
    //                                    NekMatrix<double, StandardMatrixTag>,
    //                                    AddOp, BinaryNullOp 
    //                                    >
    //    {
    //        typedef BinaryExpressionPolicy
    //                <
    //                    BinaryExpressionPolicy
    //                    <                                            
    //                        ConstantExpressionPolicy<double>,
    //                        MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                    >,
    //                    MultiplyOp, 
    //                    ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                > LhsType;
    //        typedef BinaryExpressionPolicy
    //                <
    //                    ConstantExpressionPolicy<double>,
    //                    MultiplyOp,
    //                    ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                > RhsType;
    //        
    //        static void Eval(const Expression<LhsType>& lhs, 
    //                         const Expression<RhsType>& rhs,
    //                         Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //        {
    //            return Impl::DgemmBinaryExpressionEvaluator<LhsType, RhsType,
    //                                    NekMatrix<double, StandardMatrixTag>,
    //                                    AddOp, BinaryNullOp >::Eval(lhs, rhs, result);
    //        }
    //    };


    GENERATE_MULTIPLICATION_OPERATOR(NekMatrix, 2, NekMatrix, 2);

    GENERATE_MULTIPLICATION_OPERATOR(NekMatrix, 2, NekDouble, 0);
    GENERATE_MULTIPLICATION_OPERATOR(NekDouble, 0, NekMatrix, 2);
    GENERATE_MULTIPLICATION_OPERATOR(NekMatrix, 2, NekVector, 3);
    
    GENERATE_DIVISION_OPERATOR(NekMatrix, 2, NekMatrix, 2);
    GENERATE_ADDITION_OPERATOR(NekMatrix, 2, NekMatrix, 2);
    GENERATE_SUBTRACTION_OPERATOR(NekMatrix, 2, NekMatrix, 2);
            
}
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_OPERATIONS_HPP

