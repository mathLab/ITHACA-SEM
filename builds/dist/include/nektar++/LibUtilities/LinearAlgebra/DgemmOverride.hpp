///////////////////////////////////////////////////////////////////////////////
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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_DGEMM_OVERRIDE_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_DGEMM_OVERRIDE_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <ExpressionTemplates/ExpressionTemplates.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <LibUtilities/LinearAlgebra/CanGetRawPtr.hpp>
#include <LibUtilities/LinearAlgebra/Blas.hpp>
#include <LibUtilities/LinearAlgebra/BlasOverrideUtil.hpp>

namespace Nektar
{
    template<typename ADataType, typename BDataType, 
             typename AMatrixType, typename BMatrixType>
    void Dgemm(NekMatrix<double, StandardMatrixTag>& result, 
               double alpha, const NekMatrix<ADataType, AMatrixType>& A, const NekMatrix<BDataType, BMatrixType>& B)
    {
        if( A.GetType() != eFULL || B.GetType() != eFULL )
        {
            Multiply(result, A, B);
            MultiplyEqual(result, alpha);
            return;
        }

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
            B.GetRawPtr(), LDB, 0.0,
                    result.GetRawPtr(), result.GetRows());
    }

    template<typename ADataType, typename BDataType, 
             typename AMatrixType, typename BMatrixType>
    void Dgemm(NekMatrix<double, StandardMatrixTag>& result, 
               const NekMatrix<ADataType, AMatrixType>& A, double alpha, const NekMatrix<BDataType, BMatrixType>& B)
    {
        Dgemm(result, alpha, A, B);
    }

    template<typename ADataType, typename BDataType, 
             typename AMatrixType, typename BMatrixType>
    void Dgemm(NekMatrix<double, StandardMatrixTag>& result, 
               const NekMatrix<ADataType, AMatrixType>& A, const NekMatrix<BDataType, BMatrixType>& B, double alpha)
    {
        Dgemm(result, alpha, A, B);
    }

    template<typename ADataType, typename BDataType, typename CDataType, 
             typename AMatrixType, typename BMatrixType, typename CMatrixType>
    void Dgemm(NekMatrix<double, StandardMatrixTag>& result, 
               double alpha, const NekMatrix<ADataType, AMatrixType>& A, const NekMatrix<BDataType, BMatrixType>& B,
               double beta, const NekMatrix<CDataType, CMatrixType>& C)
    {
        if( A.GetType() != eFULL || B.GetType() != eFULL || C.GetType() != eFULL )
        {
            Multiply(result, A, B);
            MultiplyEqual(result, alpha);
            NekMatrix<double> temp = beta*C;
            AddEqual(result, temp);
            return;
        }

        if( C.GetTransposeFlag() == 'T' )
        {
            // In this case, the data array in C is not ordered correctly, 
            // so we need to copy it into result to get the right ordering.
            result.SetSize(C.GetRows(), C.GetColumns());

            for(unsigned int i = 0; i < C.GetRows(); ++i)
            {
                for(unsigned int j = 0; j < C.GetColumns(); ++j)
                {
                    result(i,j) = C(i,j);
                }
            }
        }
        else
        {
            result = C;
        }

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
}

namespace expt
{
    namespace impl
    {
        // The DgemmNodeEvaluator creates an evaluator for the requested node, or the 
        // child of a negate node.  The negation operation will occur as part of the dgemm
        // call, so we do not need to do it up front.
        template<typename NodeType, typename IndicesType, unsigned int index>
        struct DgemmNodeEvaluator
        {
            typedef EvaluateNodeWithTemporaryIfNeeded<NodeType, IndicesType, index> Evaluator;
            typedef typename Evaluator::ResultType Type;
            static double GetScale() { return 1.0; }
        };

        template<typename LhsType, typename IndicesType, unsigned int index>
        struct DgemmNodeEvaluator<Node<LhsType, NegateOp>, IndicesType, index>
        {
            typedef EvaluateNodeWithTemporaryIfNeeded<LhsType, IndicesType, index> Evaluator;
            typedef typename Evaluator::ResultType Type;
            static double GetScale() { return -1.0; }
        };
    }

    namespace impl
    {
        // To handle alpha*A*B - beta*C, this class provides the scale depending on the operator.
        template<typename OpType>
        struct BetaScale
        {
            static double GetScale() { return 1.0; }
        };

        template<>
        struct BetaScale<SubtractOp>
        {
            static double GetScale() { return -1.0; }
        };

        template<typename T>
        struct IsDouble : public boost::false_type{};

        template<>
        struct IsDouble<double> : public boost::true_type {};
        
    }

    namespace impl
    {
        template<typename NodeType, typename IndicesType, unsigned int index, typename enabled=void>
        struct AlphaABParameterAccessImpl : public boost::false_type {};

        template<typename A1, typename A2, typename A3, typename IndicesType, unsigned int index>
        struct AlphaABParameterAccessImpl< Node<Node<A1, MultiplyOp, A2>, MultiplyOp, A3>, IndicesType, index,
            typename boost::enable_if
            <
                Test3ArgumentAssociativeNode<Node<Node<A1, MultiplyOp, A2>, MultiplyOp, A3>, IsDouble, MultiplyOp,
                                             Nektar::CanGetRawPtr, MultiplyOp, Nektar::CanGetRawPtr>
                //boost::mpl::and_    
                //<
                //    boost::is_same<typename A1::ResultType, double>,
                //    CanGetRawPtr<typename A2::ResultType>,
                //    CanGetRawPtr<typename A3::ResultType>
                //>
            >::type> : public boost::true_type
        {
            typedef DgemmNodeEvaluator<A1, IndicesType, index> AlphaWrappedEvaluator;
            typedef typename AlphaWrappedEvaluator::Evaluator AlphaEvaluator;

            static const unsigned int A2Index = index + A1::TotalCount;
            typedef DgemmNodeEvaluator<A2, IndicesType, A2Index> AWrappedEvaluator;
            typedef typename AWrappedEvaluator::Evaluator AEvaluator;

            static const unsigned int A3Index = A2Index + A2::TotalCount;
            typedef DgemmNodeEvaluator<A3, IndicesType, A3Index> BWrappedEvaluator;
            typedef typename BWrappedEvaluator::Evaluator BEvaluator;

            template<typename ArgumentVectorType>
            static double GetAlpha(const ArgumentVectorType& args)
            {
                return AlphaEvaluator::Evaluate(args) * AWrappedEvaluator::GetScale() * BWrappedEvaluator::GetScale();
            }
        };

        template<typename A1, typename A2, typename A3, typename IndicesType, unsigned int index>
        struct AlphaABParameterAccessImpl< Node<Node<A1, MultiplyOp, A2>, MultiplyOp, A3>, IndicesType, index,
            typename boost::enable_if
            <
                Test3ArgumentAssociativeNode<Node<Node<A1, MultiplyOp, A2>, MultiplyOp, A3>, Nektar::CanGetRawPtr, MultiplyOp,
                                             IsDouble, MultiplyOp, Nektar::CanGetRawPtr>
                //boost::mpl::and_    
                //<
                //    CanGetRawPtr<typename A1::ResultType>,
                //    boost::is_same<typename A2::ResultType, double>,
                //    CanGetRawPtr<typename A3::ResultType>
                //>
            >::type> : public boost::true_type
        {
            typedef DgemmNodeEvaluator<A1, IndicesType, index> AWrappedEvaluator;
            typedef typename AWrappedEvaluator::Evaluator AEvaluator;

            static const unsigned int A2Index = index + A1::TotalCount;
            typedef DgemmNodeEvaluator<A2, IndicesType, A2Index> AlphaWrappedEvaluator;
            typedef typename AlphaWrappedEvaluator::Evaluator AlphaEvaluator;

            static const unsigned int A3Index = A2Index + A2::TotalCount;
            typedef DgemmNodeEvaluator<A3, IndicesType, A3Index> BWrappedEvaluator;
            typedef typename BWrappedEvaluator::Evaluator BEvaluator;

            template<typename ArgumentVectorType>
            static double GetAlpha(const ArgumentVectorType& args)
            {
                return AlphaEvaluator::Evaluate(args) * AWrappedEvaluator::GetScale() * BWrappedEvaluator::GetScale();
            }
        };

        template<typename A1, typename A2, typename A3, typename IndicesType, unsigned int index>
        struct AlphaABParameterAccessImpl< Node<Node<A1, MultiplyOp, A2>, MultiplyOp, A3>, IndicesType, index,
            typename boost::enable_if
            <
                Test3ArgumentAssociativeNode<Node<Node<A1, MultiplyOp, A2>, MultiplyOp, A3>, Nektar::CanGetRawPtr, MultiplyOp,
                                             Nektar::CanGetRawPtr, MultiplyOp, IsDouble>
                //boost::mpl::and_    
                //<
                //    CanGetRawPtr<typename A1::ResultType>,
                //    CanGetRawPtr<typename A2::ResultType>,
                //    boost::is_same<typename A3::ResultType, double>
                //>
            >::type> : public boost::true_type
        {
            typedef DgemmNodeEvaluator<A1, IndicesType, index> AWrappedEvaluator;
            typedef typename AWrappedEvaluator::Evaluator AEvaluator;

            static const unsigned int A2Index = index + A1::TotalCount;
            typedef DgemmNodeEvaluator<A2, IndicesType, A2Index> BWrappedEvaluator;
            typedef typename BWrappedEvaluator::Evaluator BEvaluator;

            static const unsigned int A3Index = A2Index + A2::TotalCount;
            typedef DgemmNodeEvaluator<A3, IndicesType, A3Index> AlphaWrappedEvaluator;
            typedef typename AlphaWrappedEvaluator::Evaluator AlphaEvaluator;

            template<typename ArgumentVectorType>
            static double GetAlpha(const ArgumentVectorType& args)
            {
                return AlphaEvaluator::Evaluate(args) * AWrappedEvaluator::GetScale() * BWrappedEvaluator::GetScale();
            }
        };

        // Plain M*M
        template<typename A1, typename A2, typename IndicesType, unsigned int index>
        struct AlphaABParameterAccessImpl< Node<A1, MultiplyOp, A2>, IndicesType, index,
            typename boost::enable_if
            <
                boost::mpl::and_
                <
                    TestBinaryNode<Node<A1, MultiplyOp, A2>, Nektar::CanGetRawPtr, MultiplyOp, Nektar::CanGetRawPtr>,
                    boost::mpl::not_<Test3ArgumentAssociativeNode<Node<A1, MultiplyOp, A2>, Nektar::CanGetRawPtr, MultiplyOp, Nektar::CanGetRawPtr, MultiplyOp, IsDouble> >,
                    boost::mpl::not_<Test3ArgumentAssociativeNode<Node<A1, MultiplyOp, A2>, Nektar::CanGetRawPtr, MultiplyOp, IsDouble, MultiplyOp, Nektar::CanGetRawPtr> >,
                    boost::mpl::not_<Test3ArgumentAssociativeNode<Node<A1, MultiplyOp, A2>, IsDouble, MultiplyOp, Nektar::CanGetRawPtr, MultiplyOp, Nektar::CanGetRawPtr> >
                >
                //boost::mpl::and_
                //<
                //    CanGetRawPtr<typename A1::ResultType>,
                //    CanGetRawPtr<typename A2::ResultType>,
                //    boost::is_same<typename A3::ResultType, double>
                //>
            >::type> : public boost::true_type
        {
            typedef DgemmNodeEvaluator<A1, IndicesType, index> AWrappedEvaluator;
            typedef typename AWrappedEvaluator::Evaluator AEvaluator;

            static const unsigned int A2Index = index + A1::TotalCount;
            typedef DgemmNodeEvaluator<A2, IndicesType, A2Index> BWrappedEvaluator;
            typedef typename BWrappedEvaluator::Evaluator BEvaluator;

            template<typename ArgumentVectorType>
            static double GetAlpha(const ArgumentVectorType& args)
            {
                return 1.0 * AWrappedEvaluator::GetScale() * BWrappedEvaluator::GetScale();
            }
        };

        // The separation here allows us to separate the negatin from the AlphaABNode.
        template<typename NodeType, typename IndicesType, unsigned int index>
        struct AlphaABParameterAccess : public AlphaABParameterAccessImpl<NodeType, IndicesType, index> {};
        
        template<typename T, typename IndicesType, unsigned int index>
        struct AlphaABParameterAccess<Node<T, NegateOp, void>, IndicesType, index> : public AlphaABParameterAccessImpl<T, IndicesType, index> 
        {
            typedef Node<T, NegateOp, void> NodeType;

            template<typename ArgumentVectorType>
            static double GetAlpha(const ArgumentVectorType& args)
            {
                return -AlphaABParameterAccessImpl<T, IndicesType, index>::GetAlpha(args);
            }
        };
        

        template<typename NodeType, typename IndicesType, unsigned int index, typename enabled=void>
        struct BetaCParameterAccessImpl : public boost::false_type {};

        template<typename T, typename IndicesType, unsigned int index>
        struct BetaCParameterAccessImpl< Node<T, void, void>, IndicesType, index> : public boost::true_type
        {
            typedef Node<T, void, void> NodeType;
            typedef DgemmNodeEvaluator<NodeType, IndicesType, index> CWrappedEvaluator;
            typedef typename CWrappedEvaluator::Evaluator CEvaluator;

            template<typename ArgumentVectorType>
            static double GetBeta(const ArgumentVectorType& args)
            {
                return 1.0;
            }
        };

        template<typename L, typename R, typename IndicesType, unsigned int index>
        struct BetaCParameterAccessImpl< Node<L, expt::MultiplyOp, R>, IndicesType, index, 
            typename boost::enable_if
            <
                boost::mpl::and_    
                <
                Nektar::CanGetRawPtr<typename L::ResultType>,
                    boost::is_same<typename R::ResultType, double>
                >
            >::type> : public boost::true_type
        {
            typedef DgemmNodeEvaluator<L, IndicesType, index> CWrappedEvaluator;
            typedef typename CWrappedEvaluator::Evaluator CEvaluator;

            static const unsigned int nextIndex = index + L::TotalCount;
            typedef DgemmNodeEvaluator<R, IndicesType, nextIndex> BetaWrappedEvaluator;
            typedef typename BetaWrappedEvaluator::Evaluator BetaEvaluator;
            
            template<typename ArgumentVectorType>
            static double GetBeta(const ArgumentVectorType& args)
            {
                return BetaEvaluator::Evaluate(args) * CWrappedEvaluator::GetScale();
            }
        };

        template<typename L, typename R, typename IndicesType, unsigned int index>
        struct BetaCParameterAccessImpl< Node<L, expt::MultiplyOp, R>, IndicesType, index, 
            typename boost::enable_if
            <
                boost::mpl::and_    
                <
                    Nektar::CanGetRawPtr<typename R::ResultType>,
                    boost::is_same<typename L::ResultType, double>
                >
            >::type> : public boost::true_type
        {
            typedef DgemmNodeEvaluator<L, IndicesType, index> BetaWrappedEvaluator;
            typedef typename BetaWrappedEvaluator::Evaluator BetaEvaluator;

            static const unsigned int nextIndex = index + L::TotalCount;
            typedef DgemmNodeEvaluator<R, IndicesType, nextIndex> CWrappedEvaluator;
            typedef typename CWrappedEvaluator::Evaluator CEvaluator;

            template<typename ArgumentVectorType>
            static double GetBeta(const ArgumentVectorType& args)
            {
                return BetaEvaluator::Evaluate(args) * CWrappedEvaluator::GetScale();
            }
        };

        template<typename NodeType, typename IndicesType, unsigned int index>
        struct BetaCParameterAccess : public BetaCParameterAccessImpl<NodeType, IndicesType, index> {};

        template<typename T, typename IndicesType, unsigned int index>
        struct BetaCParameterAccess<Node<T, NegateOp, void>, IndicesType, index> : public BetaCParameterAccessImpl<T, IndicesType, index> 
        {
            template<typename ArgumentVectorType>
            static double GetBeta(const ArgumentVectorType& args)
            {
                return -BetaCParameterAccessImpl<T, IndicesType, index>::GetBeta(args);
            }
        };

    }

    // Three cases to handle.
    // 1. alpha*A*B +/- beta*C
    // 2. beta*C +/- alpha*A*B - Implement through commutative transform.
    // 3. alpha*A*B +/- beta*C*D - This has to be considered separately because it matches both 1 and 2.
    // 4. alpha*A*B

    ///////////////////////////////////////////////////////////
    // Case 1: alpha*A*B +/- beta*C
    // Case 3 as well.
    ///////////////////////////////////////////////////////////
    template<typename L, typename OpType, typename R, typename IndicesType, unsigned int index>
    struct BinaryBinaryEvaluateNodeOverride<L, OpType, R, IndicesType, index,
        typename boost::enable_if
        <
            boost::mpl::and_
            <
                IsAdditiveOperator<OpType>,
                impl::AlphaABParameterAccess<L, IndicesType, index>,
                impl::BetaCParameterAccess<R, IndicesType, index + L::TotalCount>
                //IsAlphaABNode<L>,
                //IsBetaCNode<R>
            >
        >::type> : public boost::true_type
    {
        typedef Node<L, OpType, R> NodeType;
        typedef impl::AlphaABParameterAccess<L, IndicesType, index> LhsAccess;
        static const unsigned int rhsIndex = index + L::TotalCount;
        typedef impl::BetaCParameterAccess<R, IndicesType, rhsIndex> RhsAccess;
        typedef typename LhsAccess::AEvaluator AEvaluator;
        typedef typename LhsAccess::BEvaluator BEvaluator;
        typedef typename RhsAccess::CEvaluator CEvaluator;

        template<typename ResultType, typename ArgumentVectorType>
        static void Evaluate(ResultType& accumulator, const ArgumentVectorType& args)
        {
            typename AEvaluator::ResultType a = AEvaluator::Evaluate(args);
            typename BEvaluator::ResultType b = BEvaluator::Evaluate(args);
            typename CEvaluator::ResultType c = CEvaluator::Evaluate(args);

            double alpha = LhsAccess::GetAlpha(args);
            double beta = RhsAccess::GetBeta(args);

            beta *= impl::BetaScale<OpType>::GetScale();
            Nektar::Dgemm(accumulator, alpha, a, b, beta, c);
        }
    };

    ////////////////////////////////////////////////////////////
    // Case 2 - beta*C +/- alpha*A*B
    ////////////////////////////////////////////////////////////
    template<typename L, typename OpType, typename R, typename IndicesType, unsigned int index>
    struct BinaryBinaryEvaluateNodeOverride<L, OpType, R, IndicesType, index,
        typename boost::enable_if
        <
            boost::mpl::and_
            <
                IsAdditiveOperator<OpType>,
                impl::AlphaABParameterAccess<R, IndicesType, index + L::TotalCount>,
                impl::BetaCParameterAccess<L, IndicesType, index>,
                boost::mpl::not_<impl::AlphaABParameterAccess<L, IndicesType, index> >
                //IsAlphaABNode<R>,
                //IsBetaCNode<L>,
                //boost::mpl::not_<IsAlphaABNode<L> >
            >
        >::type> : public boost::true_type
    {
        typedef Node<L, OpType, R> NodeType;

        typedef impl::BetaCParameterAccess<L, IndicesType, index> LhsAccess;
        static const unsigned int rhsIndex = index + L::TotalCount;
        typedef impl::AlphaABParameterAccess<R, IndicesType, rhsIndex> RhsAccess;
        
        typedef typename RhsAccess::AEvaluator AEvaluator;
        typedef typename RhsAccess::BEvaluator BEvaluator;
        typedef typename LhsAccess::CEvaluator CEvaluator;

        template<typename ResultType, typename ArgumentVectorType>
        static void Evaluate(ResultType& accumulator, const ArgumentVectorType& args)
        {
            typename AEvaluator::ResultType a = AEvaluator::Evaluate(args);
            typename BEvaluator::ResultType b = BEvaluator::Evaluate(args);
            typename CEvaluator::ResultType c = CEvaluator::Evaluate(args);

            double alpha = RhsAccess::GetAlpha(args);
            double beta = LhsAccess::GetBeta(args);

            alpha *= impl::BetaScale<OpType>::GetScale();
            Nektar::Dgemm(accumulator, alpha, a, b, beta, c);
        }
    };

    ////////////////////////////////////////////////////////////
    // Case 4 - alpha*A*B
    ////////////////////////////////////////////////////////////
    template<typename L, typename OpType, typename R, typename IndicesType, unsigned int index>
    struct BinaryBinaryEvaluateNodeOverride<L, OpType, R, IndicesType, index,
        typename boost::enable_if
        <
            //IsAlphaABNode<Node<L, OpType, R> >
            impl::AlphaABParameterAccess<Node<L, OpType, R>, IndicesType, index>
        >::type> : public boost::true_type
    {
        typedef Node<L, OpType, R> NodeType;
        typedef impl::AlphaABParameterAccess<NodeType, IndicesType, index> LhsAccess;
        typedef typename LhsAccess::AEvaluator AEvaluator;
        typedef typename LhsAccess::BEvaluator BEvaluator;

        template<typename ResultType, typename ArgumentVectorType>
        static void Evaluate(ResultType& accumulator, const ArgumentVectorType& args)
        {
            typename AEvaluator::ResultType a = AEvaluator::Evaluate(args);
            typename BEvaluator::ResultType b = BEvaluator::Evaluate(args);

            double alpha = LhsAccess::GetAlpha(args);
            Nektar::Dgemm(accumulator, alpha, a, b);
        }
    };


}

#endif
#endif
