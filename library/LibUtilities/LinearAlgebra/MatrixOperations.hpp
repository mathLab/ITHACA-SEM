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
    

    ////////////////////////////////////////////////////////////////////////////////////
    // Subtraction
    ////////////////////////////////////////////////////////////////////////////////////
    template<typename DataType, typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    void Subtract(NekMatrix<DataType, StandardMatrixTag>& result,
                const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, RhsMatrixType>& rhs);

    template<typename DataType, typename RhsDataType, typename RhsMatrixType>
    void SubtractEqual(NekMatrix<DataType, StandardMatrixTag>& result,
                          const NekMatrix<RhsDataType, RhsMatrixType>& rhs);
    
    template<typename LhsDataType, typename LhsMatrixType, typename RhsDataType, typename RhsMatrixType>
    NekMatrix<typename NekMatrix<LhsDataType, LhsMatrixType>::NumberType, StandardMatrixTag> 
    Subtract(const NekMatrix<LhsDataType, LhsMatrixType>& lhs,
                const NekMatrix<RhsDataType, RhsMatrixType>& rhs);
    

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
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
                    result.GetRawPtr(), result.GetRows());
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
#endif

}


#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
namespace expt
{
    template<typename L, typename R, typename IndicesType, unsigned int index>
    struct BinaryBinaryEvaluateNodeOverride<L, expt::AddOp, R, IndicesType, index,
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
                IsDgemmLeftSide<L>::template GetValues<IndicesType, index>::GetAlpha(args),
                IsDgemmLeftSide<L>::template GetValues<IndicesType, index>::GetA(args),
                IsDgemmLeftSide<L>::template GetValues<IndicesType, index>::GetB(args),
                IsDgemmRightSide<R>::template GetValues<IndicesType, index + L::TotalCount>::GetBeta(args),
                IsDgemmRightSide<R>::template GetValues<IndicesType, index + L::TotalCount>::GetC(args) );
        }
    };  
}
#endif

namespace Nektar
{
    // TODO - This case also catches chained multipliation, which we don't want.
    // But we do want it to catch a*A*B.
    //template<typename L, typename R, typename IndicesType, unsigned int index>
    //struct BinaryBinaryEvaluateNodeOverride<L, expt::MultiplyOp, R, IndicesType, index,
    //    typename boost::enable_if
    //    <
    //        IsDgemmLeftSide<Node<L, expt::MultiplyOp, R> >
    //    >::type > : public boost::true_type 
    //{
    //    typedef IsDgemmLeftSide<Node<L, expt::MultiplyOp, R> > Accessor;

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
    //                                            expt::MultiplyOp, 
    //                                            ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                        >,
    //                                        ConstantExpressionPolicy<NekMatrix<T3, M3> >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        expt::AddOp, BinaryNullOp,
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
    //                        expt::MultiplyOp, 
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
    //                    expt::AddOp<LhsResultType, RhsResultType>::ApplyEqual(result, *rhs);
    //                }
    //            }
    //        };
    //        
    //        // AB + C - in case there is a block matrix in the expression.
    //        template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //        struct DgemmBinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                                            expt::MultiplyOp, 
    //                                            ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                        >,
    //                                        ConstantExpressionPolicy<NekMatrix<T3, M3> >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        expt::AddOp, BinaryNullOp,
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
    //                        expt::MultiplyOp, 
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
    //                expt::AddOp<LhsResultType, RhsResultType>::ApplyEqual(result, *rhs);
    //            }
    //        };


    //        // aAB + C
    //        template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //        struct DgemmBinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                        <
    //                                            BinaryExpressionPolicy
    //                                            <
    //                                                ConstantExpressionPolicy<double>,
    //                                                expt::MultiplyOp,
    //                                                ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                                            >,
    //                                            expt::MultiplyOp, 
    //                                            ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                        >,
    //                                        ConstantExpressionPolicy<NekMatrix<T3, M3> >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        expt::AddOp, BinaryNullOp,
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
    //                        expt::MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                    > LhsLhsType;
    //                    
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        LhsLhsType,
    //                        expt::MultiplyOp, 
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
    //                    expt::AddOp<LhsResultType, RhsResultType>::ApplyEqual(result, *rhs);
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
    //                                                expt::MultiplyOp,
    //                                                ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                                            >,
    //                                            expt::MultiplyOp, 
    //                                            ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                        >,
    //                                        ConstantExpressionPolicy<NekMatrix<T3, M3> >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        expt::AddOp, BinaryNullOp,
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
    //                        expt::MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                    > LhsLhsType;
    //                    
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        LhsLhsType,
    //                        expt::MultiplyOp, 
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
    //                expt::AddOp<LhsResultType, RhsResultType>::ApplyEqual(result, *rhs);
    //            }
    //        };


    //        // AB + bC
    //        template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //        struct DgemmBinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                                            expt::MultiplyOp, 
    //                                            ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                        >,
    //                                        BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<double>,
    //                                            expt::MultiplyOp,
    //                                            ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                                        >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        expt::AddOp, BinaryNullOp,
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
    //                        expt::MultiplyOp, 
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    > LhsType;
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        expt::MultiplyOp,
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
    //                    expt::AddOp<NekMatrix<double>, NekMatrix<double> >::ApplyEqual(result, rhsTemp);
    //                }
    //            }
    //        };
    //        
    //        // AB + bC - in case there is a block matrix in the expression.
    //        template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //        struct DgemmBinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                                            expt::MultiplyOp, 
    //                                            ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                        >,
    //                                        BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<double>,
    //                                            expt::MultiplyOp,
    //                                            ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                                        >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        expt::AddOp, BinaryNullOp,
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
    //                        expt::MultiplyOp, 
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    > LhsType;
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        expt::MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                    > RhsType;
    //            
    //            static void Eval(const Expression<LhsType>& lhs, 
    //                             const Expression<RhsType>& rhs,
    //                             Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //            {
    //                lhs.Evaluate(result);
    //                NekMatrix<double> rhsTemp = rhs.Evaluate();
    //                expt::AddOp<NekMatrix<double>, NekMatrix<double> >::ApplyEqual(result, rhsTemp);
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
    //                                                expt::MultiplyOp,
    //                                                ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                                            >,
    //                                            expt::MultiplyOp,
    //                                            ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                        >,
    //                                        BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<double>,
    //                                            expt::MultiplyOp,
    //                                            ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                                        >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        expt::AddOp, BinaryNullOp,
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
    //                        expt::MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                    > LhsLhsType;
    //                    
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        LhsLhsType,
    //                        expt::MultiplyOp, 
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    > LhsType;
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        expt::MultiplyOp,
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
    //                    expt::AddOp<LhsResultType, RhsResultType>::ApplyEqual(result, rhsTemp);
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
    //                                                expt::MultiplyOp,
    //                                                ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                                            >,
    //                                            expt::MultiplyOp,
    //                                            ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                        >,
    //                                        BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<double>,
    //                                            expt::MultiplyOp,
    //                                            ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                                        >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        expt::AddOp, BinaryNullOp,
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
    //                        expt::MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                    > LhsLhsType;
    //                    
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        LhsLhsType,
    //                        expt::MultiplyOp, 
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    > LhsType;
    //             typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        expt::MultiplyOp,
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
    //                expt::AddOp<LhsResultType, RhsResultType>::ApplyEqual(result, rhsTemp);
    //            }
    //        };
    //        
    //        
    //        // a(AB) + bC
    //        template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //        struct DgemmBinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<double>,
    //                                            expt::MultiplyOp,
    //                                            BinaryExpressionPolicy
    //                                            <
    //                                                ConstantExpressionPolicy<NekMatrix<T1, M1> >,
    //                                                expt::MultiplyOp,
    //                                                ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                            >
    //                                        >,
    //                                        BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<double>,
    //                                            expt::MultiplyOp,
    //                                            ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                                        >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        expt::AddOp, BinaryNullOp,
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
    //                        expt::MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    > LhsRhsType;
    //                    
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        expt::MultiplyOp, 
    //                        LhsRhsType
    //                    > LhsType;
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        expt::MultiplyOp,
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
    //                    expt::AddOp<LhsResultType, RhsResultType>::ApplyEqual(result, rhsTemp);
    //                }
    //            }
    //        };
    //        
    //        // a(AB) + bC - in case there is a block matrix in the expression.
    //        template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //        struct DgemmBinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<double>,
    //                                            expt::MultiplyOp,
    //                                            BinaryExpressionPolicy
    //                                            <
    //                                                ConstantExpressionPolicy<NekMatrix<T1, M1> >,
    //                                                expt::MultiplyOp,
    //                                                ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                            >
    //                                        >,
    //                                        BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<double>,
    //                                            expt::MultiplyOp,
    //                                            ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                                        >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        expt::AddOp, BinaryNullOp,
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
    //                        expt::MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    > LhsRhsType;
    //                    
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        expt::MultiplyOp, 
    //                        LhsRhsType
    //                    > LhsType;
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        expt::MultiplyOp,
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
    //                expt::AddOp<LhsResultType, RhsResultType>::ApplyEqual(result, rhsTemp);
    //            }
    //        };
    //        
    //        // a(AB) + C
    //        template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //        struct DgemmBinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<double>,
    //                                            expt::MultiplyOp,
    //                                            BinaryExpressionPolicy
    //                                            <
    //                                                ConstantExpressionPolicy<NekMatrix<T1, M1> >,
    //                                                expt::MultiplyOp,
    //                                                ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                            >
    //                                        >,
    //                                        ConstantExpressionPolicy<NekMatrix<T3, M3> >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        expt::AddOp, BinaryNullOp,
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
    //                        expt::MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    > LhsRhsType;
    //                    
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        expt::MultiplyOp, 
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
    //                    expt::AddOp<LhsResultType, RhsResultType>::ApplyEqual(result, rhsTemp);
    //                }
    //            }
    //        };
    //        
    //        // a(AB) + C - in case there is a block matrix in the expression.
    //        template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //        struct DgemmBinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<double>,
    //                                            expt::MultiplyOp,
    //                                            BinaryExpressionPolicy
    //                                            <
    //                                                ConstantExpressionPolicy<NekMatrix<T1, M1> >,
    //                                                expt::MultiplyOp,
    //                                                ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                            >
    //                                        >,
    //                                        ConstantExpressionPolicy<NekMatrix<T3, M3> >,
    //                                        NekMatrix<double, StandardMatrixTag>,
    //                                        expt::AddOp, BinaryNullOp,
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
    //                        expt::MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    > LhsRhsType;
    //                    
    //            typedef BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                        expt::MultiplyOp, 
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
    //                expt::AddOp<LhsResultType, RhsResultType>::ApplyEqual(result, rhsTemp);
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
    //                                        expt::MultiplyOp, 
    //                                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                    >,
    //                                    ConstantExpressionPolicy<NekMatrix<T3, M3> >,
    //                                    NekMatrix<double, StandardMatrixTag>,
    //                                    expt::AddOp, BinaryNullOp 
    //                                    >
    //    {
    //        typedef BinaryExpressionPolicy
    //                <
    //                    ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                    expt::MultiplyOp, 
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
    //                                    expt::AddOp, BinaryNullOp >::Eval(lhs, rhs, result);
    //        }
    //    };

    //    // aAB + C
    //    template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //    struct BinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                    <
    //                                        ConstantExpressionPolicy<double>,
    //                                        expt::MultiplyOp,
    //                                        BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                                            expt::MultiplyOp, 
    //                                            ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                        >
    //                                    >,
    //                                    ConstantExpressionPolicy<NekMatrix<T3, M3> >,
    //                                    NekMatrix<double, StandardMatrixTag>,
    //                                    expt::AddOp, BinaryNullOp 
    //                                    >
    //    {
    //        typedef BinaryExpressionPolicy
    //                <
    //                    ConstantExpressionPolicy<double>,
    //                    expt::MultiplyOp,
    //                    BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                        expt::MultiplyOp, 
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
    //                                    expt::AddOp, BinaryNullOp >::Eval(lhs, rhs, result);
    //        }
    //    };


    //    // aAB + C
    //    template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //    struct BinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                    <
    //                                        BinaryExpressionPolicy
    //                                        <                                            
    //                                            ConstantExpressionPolicy<double>,
    //                                            expt::MultiplyOp,
    //                                            ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                                        >,
    //                                        expt::MultiplyOp, 
    //                                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                    >,
    //                                    ConstantExpressionPolicy<NekMatrix<T3, M3> >,
    //                                    NekMatrix<double, StandardMatrixTag>,
    //                                    expt::AddOp, BinaryNullOp 
    //                                    >
    //    {
    //        typedef BinaryExpressionPolicy
    //                <
    //                    BinaryExpressionPolicy
    //                    <                                            
    //                        ConstantExpressionPolicy<double>,
    //                        expt::MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                    >,
    //                    expt::MultiplyOp, 
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
    //                                    expt::AddOp, BinaryNullOp >::Eval(lhs, rhs, result);
    //        }
    //    };
    //    
    //    // AB + bC
    //    template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //    struct BinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                    <
    //                                        ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                                        expt::MultiplyOp, 
    //                                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                    >,
    //                                    BinaryExpressionPolicy
    //                                    <
    //                                        ConstantExpressionPolicy<double>,
    //                                        expt::MultiplyOp,
    //                                        ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                                    >,
    //                                    NekMatrix<double, StandardMatrixTag>,
    //                                    expt::AddOp, BinaryNullOp>
    //    {
    //        typedef BinaryExpressionPolicy
    //                <
    //                    ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                    expt::MultiplyOp, 
    //                    ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                > LhsType;
    //        typedef BinaryExpressionPolicy
    //                <
    //                    ConstantExpressionPolicy<double>,
    //                    expt::MultiplyOp,
    //                    ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                > RhsType;
    //        
    //        static void Eval(const Expression<LhsType>& lhs, 
    //                         const Expression<RhsType>& rhs,
    //                         Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //        {
    //            return Impl::DgemmBinaryExpressionEvaluator<LhsType, RhsType,
    //                                NekMatrix<double, StandardMatrixTag>,
    //                    expt::AddOp, BinaryNullOp >::Eval(lhs, rhs, result);
    //        }
    //    };

    //    // aAB + bC
    //    template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //    struct BinaryExpressionEvaluator<BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<double>,
    //                                    expt::MultiplyOp,
    //                                        BinaryExpressionPolicy
    //                                        <
    //                                            ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                                            expt::MultiplyOp, 
    //                                            ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                        >
    //                                    >,
    //                                    BinaryExpressionPolicy
    //                                    <
    //                                        ConstantExpressionPolicy<double>,
    //                                        expt::MultiplyOp,
    //                                        ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                                    >,
    //                                    NekMatrix<double, StandardMatrixTag>,
    //                                    expt::AddOp, BinaryNullOp 
    //                                    >
    //    {
    //        typedef BinaryExpressionPolicy
    //                <
    //                    ConstantExpressionPolicy<double>,
    //                    expt::MultiplyOp,
    //                    BinaryExpressionPolicy
    //                    <
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >, 
    //                        expt::MultiplyOp, 
    //                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                    >
    //                > LhsType;
    //        typedef BinaryExpressionPolicy
    //                <
    //                    ConstantExpressionPolicy<double>,
    //                    expt::MultiplyOp,
    //                    ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                > RhsType;
    //        
    //        static void Eval(const Expression<LhsType>& lhs, 
    //                         const Expression<RhsType>& rhs,
    //                         Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //        {
    //            return Impl::DgemmBinaryExpressionEvaluator<LhsType, RhsType,
    //                                    NekMatrix<double, StandardMatrixTag>,
    //                                    expt::AddOp, BinaryNullOp >::Eval(lhs, rhs, result);
    //        }
    //    };

    //    // (aA)B + bC
    //    template<typename T1, typename M1, typename T2, typename M2, typename T3, typename M3>
    //    struct BinaryExpressionEvaluator<BinaryExpressionPolicy
    //                                    <
    //                                        BinaryExpressionPolicy
    //                                        <                                            
    //                                            ConstantExpressionPolicy<double>,
    //                                            expt::MultiplyOp,
    //                                            ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                                        >,
    //                                        expt::MultiplyOp, 
    //                                        ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                                    >,
    //                                    BinaryExpressionPolicy
    //                                    <
    //                                        ConstantExpressionPolicy<double>,
    //                                        expt::MultiplyOp,
    //                                        ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                                    >,
    //                                    NekMatrix<double, StandardMatrixTag>,
    //                                    expt::AddOp, BinaryNullOp 
    //                                    >
    //    {
    //        typedef BinaryExpressionPolicy
    //                <
    //                    BinaryExpressionPolicy
    //                    <                                            
    //                        ConstantExpressionPolicy<double>,
    //                        expt::MultiplyOp,
    //                        ConstantExpressionPolicy<NekMatrix<T1, M1> >
    //                    >,
    //                    expt::MultiplyOp, 
    //                    ConstantExpressionPolicy<NekMatrix<T2, M2> > 
    //                > LhsType;
    //        typedef BinaryExpressionPolicy
    //                <
    //                    ConstantExpressionPolicy<double>,
    //                    expt::MultiplyOp,
    //                    ConstantExpressionPolicy<NekMatrix<T3, M3> >
    //                > RhsType;
    //        
    //        static void Eval(const Expression<LhsType>& lhs, 
    //                         const Expression<RhsType>& rhs,
    //                         Accumulator<NekMatrix<double, StandardMatrixTag> >& result)
    //        {
    //            return Impl::DgemmBinaryExpressionEvaluator<LhsType, RhsType,
    //                                    NekMatrix<double, StandardMatrixTag>,
    //                                    expt::AddOp, BinaryNullOp >::Eval(lhs, rhs, result);
    //        }
    //    };


    GENERATE_MULTIPLICATION_OPERATOR(NekMatrix, 2, NekMatrix, 2);

    GENERATE_MULTIPLICATION_OPERATOR(NekMatrix, 2, NekDouble, 0);
    GENERATE_MULTIPLICATION_OPERATOR(NekDouble, 0, NekMatrix, 2);
    GENERATE_MULTIPLICATION_OPERATOR(NekMatrix, 2, NekVector, 1);
    
    GENERATE_DIVISION_OPERATOR(NekMatrix, 2, NekMatrix, 2);
    GENERATE_ADDITION_OPERATOR(NekMatrix, 2, NekMatrix, 2);
    GENERATE_SUBTRACTION_OPERATOR(NekMatrix, 2, NekMatrix, 2);
            
}
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_MATRIX_OPERATIONS_HPP

