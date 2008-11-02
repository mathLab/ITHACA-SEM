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
#include <boost/utility/enable_if.hpp>

#include <string>

namespace Nektar
{
    ////////////////////////////////////////////////////////////////////////////////////
    // Multiplication
    ////////////////////////////////////////////////////////////////////////////////////
    template<typename ResultDataType, typename LhsDataType, typename LhsMatrixType>
    void NekMultiply(NekMatrix<ResultDataType, StandardMatrixTag>& result,
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
    NekMultiply(const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                const DataType& rhs)
    {
        typedef typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType ResultDataType;
        NekMatrix<ResultDataType, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
        NekMultiply(result, lhs, rhs);
        return result;
    }

    template<typename RhsDataType, typename RhsMatrixType, typename ResultDataType>
    void NekMultiply(NekMatrix<ResultDataType, StandardMatrixTag>& result,
                     const ResultDataType& lhs,
                     const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
                     
    {
        NekMultiply(result, rhs, lhs);
    }
    
    template<typename DataType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<RhsDataType, RhsMatrixType>::NumberType, StandardMatrixTag> 
    NekMultiply(const DataType& lhs,
                const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
                     
    {
        return NekMultiply(rhs, lhs);
    }
    
   
    template<typename LhsDataType, typename RhsDataType, typename DataType,
             typename LhsMatrixType, typename RhsMatrixType>
    void NekMultiply(NekMatrix<DataType, StandardMatrixTag>& result,
                     const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                     const NekMatrix<RhsDataType, RhsMatrixType>& rhs
//                     #ifdef NEKTAR_USING_BLAS
//                     // This is necessary to force double matrices to use the blas specific versions below.
//                     // Without this check, NekMatrix<double> would not use the versions below because
//                     // this is a better match than NekMatrix<const double>
//                     , typename boost::disable_if
//                     <
//                        boost::mpl::and_
//                        <
//                            boost::is_same<double, typename RawType<typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType>::type>,
//                            boost::is_same<double, typename RawType<typename NekMatrix<RhsDataType, RhsMatrixType>::NumberType>::type>
//                        >
//                     >::type* p = 0
//                     #endif //NEKTAR_USING_BLAS
                     )
    {
        ASSERTL1(lhs.GetColumns() == rhs.GetRows(), std::string("A left side matrix with column count ") + 
            boost::lexical_cast<std::string>(lhs.GetColumns()) + 
            std::string(" and a right side matrix with row count ") + 
            boost::lexical_cast<std::string>(rhs.GetRows()) + std::string(" can't be multiplied."));

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
        
        
    template<typename LhsDataType>
    void NekMultiplyEqual(NekMatrix<LhsDataType, StandardMatrixTag>& lhs,
                          typename boost::call_traits<LhsDataType>::const_reference rhs)
    {
        typedef typename NekMatrix<LhsDataType, StandardMatrixTag>::iterator iterator;
        for( iterator iter = lhs.begin(); iter != lhs.end(); ++iter)
        {
            *iter *= rhs;
        }
    }

#ifdef NEKTAR_USING_BLAS    
//    /// \brief Floating point specialization when blas is in use.
//    ///
//    /// The signature of this function should match the unspecialized NekMultiply.
//    void NekMultiply(NekMatrix<double, FullMatrixTag, StandardMatrixTag>& result,
//                     const NekMatrix<const double, FullMatrixTag, StandardMatrixTag>& lhs,
//                     const NekMatrix<const double, FullMatrixTag, StandardMatrixTag>& rhs);
//
//    void NekMultiply(NekMatrix<double, FullMatrixTag, StandardMatrixTag>& result,
//                     const NekMatrix<NekMatrix<const double, FullMatrixTag, StandardMatrixTag>, FullMatrixTag, ScaledMatrixTag>& lhs,
//                     const NekMatrix<NekMatrix<const double, FullMatrixTag, StandardMatrixTag>, FullMatrixTag, ScaledMatrixTag>& rhs);
//
//    void NekMultiply(NekMatrix<double, FullMatrixTag, StandardMatrixTag>& result,
//                     const NekMatrix<const double, FullMatrixTag, StandardMatrixTag>& lhs,
//                     const NekMatrix<NekMatrix<const double, FullMatrixTag, StandardMatrixTag>, FullMatrixTag, ScaledMatrixTag>& rhs);
//
//    void NekMultiply(NekMatrix<double, FullMatrixTag, StandardMatrixTag>& result,
//                     const NekMatrix<NekMatrix<const double, FullMatrixTag, StandardMatrixTag>, FullMatrixTag, ScaledMatrixTag>& lhs,
//                     const NekMatrix<const double, FullMatrixTag, StandardMatrixTag>& rhs);
//                     
    void NekMultiplyEqual(NekMatrix<double, StandardMatrixTag>& result,
                          const NekMatrix<const double, StandardMatrixTag>& rhs);
                          
                          
                     
#endif //NEKTAR_USING_BLAS

	template<typename LhsDataType, typename RhsDataType,
             typename LhsMatrixType, typename RhsMatrixType>
    NekMatrix<typename boost::remove_const<typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType>::type, StandardMatrixTag> 
    NekMultiply(const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        typedef typename boost::remove_const<typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType>::type NumberType;
        NekMatrix<NumberType, StandardMatrixTag> result(lhs.GetRows(), rhs.GetColumns());
        NekMultiply(result, lhs, rhs);
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
    void NekAddEqual(NekMatrix<DataType, StandardMatrixTag>& result,
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
    void NekAdd(NekMatrix<DataType, StandardMatrixTag>& result,
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
    NekAdd(const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
           const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        typedef typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType NumberType;
        NekMatrix<NumberType, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
        NekAdd(result, lhs, rhs);
        return result;
    }
    
//    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    void NekAdd(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
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
//    NekAdd(const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
//           const NekMatrix<RhsDataType, DiagonalMatrixTag, RhsMatrixType>& rhs)
//    {
//        typedef typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType NumberType;
//        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
//        NekAdd(result, lhs, rhs);
//        return result;
//    }
//    
//    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    void NekAdd(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
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
//    NekAdd(const NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>& lhs,
//           const NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>& rhs)
//    {
//        typedef typename NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>::NumberType NumberType;
//        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
//        NekAdd(result, lhs, rhs);
//        return result;
//    }
//    
//    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    void NekAdd(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
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
//    NekAdd(const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
//           const NekMatrix<RhsDataType, UpperTriangularMatrixTag, RhsMatrixType>& rhs)
//    {
//        typedef typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType NumberType;
//        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
//        NekAdd(result, lhs, rhs);
//        return result;
//    }
//
//    // UT = UT + UT
//    // Line 22
//    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    void NekAdd(NekMatrix<DataType, UpperTriangularMatrixTag, StandardMatrixTag>& result,
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
//    NekAdd(const NekMatrix<LhsDataType, UpperTriangularMatrixTag, LhsMatrixType>& lhs,
//           const NekMatrix<RhsDataType, UpperTriangularMatrixTag, RhsMatrixType>& rhs)
//    {
//        typedef typename NekMatrix<LhsDataType, UpperTriangularMatrixTag, LhsMatrixType>::NumberType NumberType;
//        NekMatrix<NumberType, UpperTriangularMatrixTag, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
//        NekAdd(result, lhs, rhs);
//        return result;
//    }
//
//    // LT = LT + LT
//    // Line 32
//    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    void NekAdd(NekMatrix<DataType, LowerTriangularMatrixTag, StandardMatrixTag>& result,
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
//    NekAdd(const NekMatrix<LhsDataType, LowerTriangularMatrixTag, LhsMatrixType>& lhs,
//           const NekMatrix<RhsDataType, LowerTriangularMatrixTag, RhsMatrixType>& rhs)
//    {
//        typedef typename NekMatrix<LhsDataType, LowerTriangularMatrixTag, LhsMatrixType>::NumberType NumberType;
//        NekMatrix<NumberType, LowerTriangularMatrixTag, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
//        NekAdd(result, lhs, rhs);
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
    void NekSubtract(NekMatrix<DataType, StandardMatrixTag>& result,
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
    void NekSubtractEqual(NekMatrix<DataType, StandardMatrixTag>& result,
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
    NekSubtract(const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, RhsMatrixType>& rhs)
    {
        typedef typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType NumberType;
        NekMatrix<NumberType, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
        NekSubtract(result, lhs, rhs);
        return result;
    }
    
//    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    void NekSubtract(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
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
//    NekSubtract(const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
//                const NekMatrix<RhsDataType, DiagonalMatrixTag, RhsMatrixType>& rhs)
//    {
//        typedef typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType NumberType;
//        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(lhs.GetRows(), lhs.GetColumns());
//        NekSubtract(result, lhs, rhs);
//        return result;
//    }
//
//    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    void NekSubtract(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
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
//    NekSubtract(const NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>& lhs,
//                const NekMatrix<RhsDataType, FullMatrixTag, RhsMatrixType>& rhs)
//    {
//        typedef typename NekMatrix<LhsDataType, DiagonalMatrixTag, LhsMatrixType>::NumberType NumberType;
//        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(rhs.GetRows(), rhs.GetColuomns());
//        NekSubtract(result, lhs, rhs);
//        return result;
//    }
//    
//    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
//    void NekSubtract(NekMatrix<DataType, FullMatrixTag, StandardMatrixTag>& result,
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
//    NekSubtract(const NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>& lhs,
//                const NekMatrix<RhsDataType, UpperTriangularMatrixTag, RhsMatrixType>& rhs)
//    {
//        typedef typename NekMatrix<LhsDataType, FullMatrixTag, LhsMatrixType>::NumberType NumberType;
//        NekMatrix<NumberType, FullMatrixTag, StandardMatrixTag> result(rhs.GetRows(), rhs.GetColumns());
//        NekSubtract(result, lhs, rhs);
//        return result;
//    }


    
    GENERATE_MULTIPLICATION_OPERATOR(NekMatrix, 2, NekMatrix, 2);
//    template<typename L1, typename L2, typename L3, typename R1, typename R2, typename R3>
//    typename Nektar::MultiplicationTraits<NekMatrix<L1, L2, L3>, NekMatrix<R1, R2, R3> >::ResultType
//    operator*(const NekMatrix<L1, L2, L3>& lhs, const NekMatrix<R1, R2, R3>& rhs)
//    {
//        return NekMultiply(lhs, rhs);
//    }
                    
//    // TODO - Any update possible for this case?  The constant value must be the same as the matrix 
//    // number type.  Seems pretty custom to me.
//    template<typename DataType, typename LhsDataType, typename StorageType, typename MatrixType>
//    typename MultiplicationTraits<NekMatrix<LhsDataType, StorageType, MatrixType>, DataType>::ResultType
//    operator*(const NekMatrix<LhsDataType, StorageType, MatrixType>& lhs, const DataType& rhs)
//    {
//        return NekMultiply(lhs, rhs);
//    }
//                        
//    template<typename DataType, typename RhsDataType, typename StorageType, typename MatrixType>
//    typename MultiplicationTraits<DataType, NekMatrix<RhsDataType, StorageType, MatrixType> >::ResultType
//    operator*(const DataType& lhs, const NekMatrix<RhsDataType, StorageType, MatrixType>& rhs)
//    {
//        return NekMultiply(lhs, rhs);
//    }
    GENERATE_MULTIPLICATION_OPERATOR(NekMatrix, 2, NekDouble, 0);
    GENERATE_MULTIPLICATION_OPERATOR(NekDouble, 0, NekMatrix, 2);
    GENERATE_MULTIPLICATION_OPERATOR(NekMatrix, 2, NekVector, 3);
    
    GENERATE_DIVISION_OPERATOR(NekMatrix, 2, NekMatrix, 2);
    GENERATE_ADDITION_OPERATOR(NekMatrix, 2, NekMatrix, 2);
    GENERATE_SUBTRACTION_OPERATOR(NekMatrix, 2, NekMatrix, 2);
    
    // TODO - Either update the GENERATE macros to allow non-type template parameters,
    // or put an if/else here for expression templates.
//    template<typename MatrixDataType, typename StorageType, typename MatrixType, 
//             typename VectorDataType, typename dim, typename space>
//    typename MultiplicationTraits<NekMatrix<MatrixDataType, StorageType, MatrixType>, NekVector<VectorDataType, dim, space> >::ResultType
//    operator*(const NekMatrix<MatrixDataType, StorageType, MatrixType>& lhs,
//              const NekVector<VectorDataType, dim, space>& rhs)
//    {
//        return NekMultiply(lhs, rhs);
//    }          
}
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_OPERATIONS_HPP

