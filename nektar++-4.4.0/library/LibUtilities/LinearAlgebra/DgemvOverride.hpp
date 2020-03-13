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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_DGEMV_OVERRIDE_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_DGEMV_OVERRIDE_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <ExpressionTemplates/ExpressionTemplates.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <LibUtilities/LinearAlgebra/CanGetRawPtr.hpp>
#include <LibUtilities/LinearAlgebra/Blas.hpp>
#include <LibUtilities/LinearAlgebra/MatrixVectorMultiplication.hpp>

namespace Nektar
{
    template<typename ADataType, typename AMatrixType>
    void Dgemv(NekVector<double>& result, 
               double alpha, const NekMatrix<ADataType, AMatrixType>& A, const NekVector<double>& x)
    {
        if( A.GetType() != eFULL )
        {
            Multiply(result, A, x);
            MultiplyEqual(result, alpha);
            return;
        }

        unsigned int M = A.GetRows();
        unsigned int N = A.GetColumns();

        char t = A.GetTransposeFlag();
        if( t == 'T' )
        {
            std::swap(M,N);
        }

        int lda = M;

        Blas::Dgemv(t, M, N, alpha*A.Scale(), A.GetRawPtr(), lda, x.GetRawPtr(), 1, 0.0, result.GetRawPtr(), 1);
    }

    template<typename ADataType, typename AMatrixType>
    void Dgemv(NekVector<double>& result, 
               double alpha, const NekMatrix<ADataType, AMatrixType>& A, const NekVector<double>& x,
               double beta, const NekVector<double>& y)
    {
        if( A.GetType() != eFULL)
        {
            Multiply(result, A, x);
            MultiplyEqual(result, alpha);
            NekVector<double> temp = beta*y;
            AddEqual(result, temp);
            return;
        }

        result = y;
        unsigned int M = A.GetRows();
        unsigned int N = A.GetColumns();

        char t = A.GetTransposeFlag();
        if( t == 'T' )
        {
            std::swap(M,N);
        }

        int lda = M;
        Blas::Dgemv(t, M, N, alpha*A.Scale(), A.GetRawPtr(), lda, x.GetRawPtr(), 1, beta, result.GetRawPtr(), 1);
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
        struct DgemvNodeEvaluator
        {
            typedef EvaluateNodeWithTemporaryIfNeeded<NodeType, IndicesType, index> Evaluator;
            typedef typename Evaluator::ResultType Type;
            static double GetScale() { return 1.0; }
        };

        template<typename LhsType, typename IndicesType, unsigned int index>
        struct DgemvNodeEvaluator<Node<LhsType, NegateOp>, IndicesType, index>
        {
            typedef EvaluateNodeWithTemporaryIfNeeded<LhsType, IndicesType, index> Evaluator;
            typedef typename Evaluator::ResultType Type;
            static double GetScale() { return -1.0; }
        };
    }

    namespace dgemv_impl
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
    }

    namespace impl
    {
        template<typename NodeType, typename IndicesType, unsigned int index, typename enabled=void>
        struct AlphaAXParameterAccessImpl : public boost::false_type {};

        template<typename A1, typename A2, typename A3, typename IndicesType, unsigned int index>
        struct AlphaAXParameterAccessImpl< Node<Node<A1, MultiplyOp, A2>, MultiplyOp, A3>, IndicesType, index,
            typename boost::enable_if
            <
                Test3ArgumentAssociativeNode<Node<Node<A1, MultiplyOp, A2>, MultiplyOp, A3>, IsDouble, MultiplyOp,
                                            Nektar::CanGetRawPtr, MultiplyOp, Nektar::IsVector>
//                boost::mpl::and_
//                <
//                    boost::is_same<typename A1::ResultType, double>,
//                    CanGetRawPtr<typename A2::ResultType>,
//                    boost::is_same<typename A3::ResultType, NekVector<double> >
//                >
            >::type> : public boost::true_type
        {
            typedef DgemvNodeEvaluator<A1, IndicesType, index> AlphaWrappedEvaluator;
            typedef typename AlphaWrappedEvaluator::Evaluator AlphaEvaluator;

            static const unsigned int A2Index = index + A1::TotalCount;
            typedef DgemvNodeEvaluator<A2, IndicesType, A2Index> AWrappedEvaluator;
            typedef typename AWrappedEvaluator::Evaluator AEvaluator;

            static const unsigned int A3Index = A2Index + A2::TotalCount;
            typedef DgemvNodeEvaluator<A3, IndicesType, A3Index> XWrappedEvaluator;
            typedef typename XWrappedEvaluator::Evaluator BEvaluator;

            template<typename ArgumentVectorType>
            static double GetAlpha(const ArgumentVectorType& args)
            {
                return AlphaEvaluator::Evaluate(args) * AWrappedEvaluator::GetScale() * XWrappedEvaluator::GetScale();
            }
        };

        template<typename A1, typename A2, typename A3, typename IndicesType, unsigned int index>
        struct AlphaAXParameterAccessImpl< Node<Node<A1, MultiplyOp, A2>, MultiplyOp, A3>, IndicesType, index,
            typename boost::enable_if
            <
                Test3ArgumentAssociativeNode<Node<Node<A1, MultiplyOp, A2>, MultiplyOp, A3>, Nektar::CanGetRawPtr, MultiplyOp,
                                             IsDouble, MultiplyOp, Nektar::IsVector>
//                boost::mpl::and_
//                <
//                    CanGetRawPtr<typename A1::ResultType>,
//                    boost::is_same<typename A2::ResultType, double>,
//                    boost::is_same<typename A3::ResultType, NekVector<double> >
//                >
            >::type> : public boost::true_type 
        {
            typedef DgemvNodeEvaluator<A1, IndicesType, index> AWrappedEvaluator;
            typedef typename AWrappedEvaluator::Evaluator AEvaluator;

            static const unsigned int A2Index = index + A1::TotalCount;
            typedef DgemvNodeEvaluator<A2, IndicesType, A2Index> AlphaWrappedEvaluator;
            typedef typename AlphaWrappedEvaluator::Evaluator AlphaEvaluator;

            static const unsigned int A3Index = A2Index + A2::TotalCount;
            typedef DgemvNodeEvaluator<A3, IndicesType, A3Index> XWrappedEvaluator;
            typedef typename XWrappedEvaluator::Evaluator BEvaluator;

            template<typename ArgumentVectorType>
            static double GetAlpha(const ArgumentVectorType& args)
            {
                return AlphaEvaluator::Evaluate(args) * AWrappedEvaluator::GetScale() * XWrappedEvaluator::GetScale();
            }
        };

        template<typename A1, typename A2, typename A3, typename IndicesType, unsigned int index>
        struct AlphaAXParameterAccessImpl< Node<Node<A1, MultiplyOp, A2>, MultiplyOp, A3>, IndicesType, index,
            typename boost::enable_if
            <
                Test3ArgumentAssociativeNode<Node<Node<A1, MultiplyOp, A2>, MultiplyOp, A3>, Nektar::CanGetRawPtr, MultiplyOp,
                                             Nektar::IsVector, MultiplyOp, IsDouble>
//                boost::mpl::and_
//                <
//                    CanGetRawPtr<typename A1::ResultType>,
//                    boost::is_same<typename A2::ResultType, NekVector<double> >,
//                    boost::is_same<typename A3::ResultType, double>
//                >
            >::type> : public boost::true_type 
        {
            typedef DgemvNodeEvaluator<A1, IndicesType, index> AWrappedEvaluator;
            typedef typename AWrappedEvaluator::Evaluator AEvaluator;

            static const unsigned int A2Index = index + A1::TotalCount;
            typedef DgemvNodeEvaluator<A2, IndicesType, A2Index> XWrappedEvaluator;
            typedef typename XWrappedEvaluator::Evaluator XEvaluator;

            static const unsigned int A3Index = A2Index + A2::TotalCount;
            typedef DgemvNodeEvaluator<A3, IndicesType, A3Index> AlphaWrappedEvaluator;
            typedef typename AlphaWrappedEvaluator::Evaluator AlphaEvaluator;

            template<typename ArgumentVectorType>
            static double GetAlpha(const ArgumentVectorType& args)
            {
                return AlphaEvaluator::Evaluate(args) * AWrappedEvaluator::GetScale() * XWrappedEvaluator::GetScale();
            }
        };

        // A*B alone.
        template<typename A1, typename A2, typename IndicesType, unsigned int index>
        struct AlphaAXParameterAccessImpl< Node<A1, MultiplyOp, A2>, IndicesType, index,
            typename boost::enable_if
            <
                boost::mpl::and_    
                <
                    TestBinaryNode<Node<A1, MultiplyOp, A2>, Nektar::CanGetRawPtr, MultiplyOp, Nektar::IsVector>,
                    boost::mpl::not_<Test3ArgumentAssociativeNode<Node<A1, MultiplyOp, A2>, IsDouble, MultiplyOp, Nektar::CanGetRawPtr, MultiplyOp, Nektar::IsVector> >,
                    boost::mpl::not_<Test3ArgumentAssociativeNode<Node<A1, MultiplyOp, A2>, Nektar::CanGetRawPtr, MultiplyOp, IsDouble, MultiplyOp, Nektar::IsVector> >,
                    boost::mpl::not_<Test3ArgumentAssociativeNode<Node<A1, MultiplyOp, A2>, Nektar::CanGetRawPtr, MultiplyOp, Nektar::IsVector, MultiplyOp, IsDouble> >
//                    CanGetRawPtr<typename A1::ResultType>,
//                    boost::is_same<typename A2::ResultType, NekVector<double> >
                >
            >::type> : public boost::true_type 
        {
            typedef DgemvNodeEvaluator<A1, IndicesType, index> AWrappedEvaluator;
            typedef typename AWrappedEvaluator::Evaluator AEvaluator;

            static const unsigned int A2Index = index + A1::TotalCount;
            typedef DgemvNodeEvaluator<A2, IndicesType, A2Index> XWrappedEvaluator;
            typedef typename XWrappedEvaluator::Evaluator XEvaluator;

            template<typename ArgumentVectorType>
            static double GetAlpha(const ArgumentVectorType& args)
            {
                return 1.0 * AWrappedEvaluator::GetScale() * XWrappedEvaluator::GetScale();
            }
        };

        template<typename NodeType, typename IndicesType, unsigned int index>
        struct AlphaAXParameterAccess : public AlphaAXParameterAccessImpl<NodeType, IndicesType, index> {};
        
        template<typename T, typename IndicesType, unsigned int index>
        struct AlphaAXParameterAccess<Node<T, NegateOp, void>, IndicesType, index> : public AlphaAXParameterAccessImpl<T, IndicesType, index> 
        {
            typedef Node<T, NegateOp, void> NodeType;

            template<typename ArgumentVectorType>
            static double GetAlpha(const ArgumentVectorType& args)
            {
                return -AlphaAXParameterAccessImpl<T, IndicesType, index>::GetAlpha(args);
            }
        };
        

        template<typename NodeType, typename IndicesType, unsigned int index, typename enabled=void>
        struct BetaYParameterAccessImpl : public boost::false_type {};

        //template<typename T, typename IndicesType, unsigned int index>
        //struct BetaYParameterAccessImpl< Node<T, void, void>, IndicesType, index> : public boost::true_type 
        //{
        //    typedef Node<T, void, void> NodeType;
        //    typedef DgemvNodeEvaluator<NodeType, IndicesType, index> YWrappedEvaluator;
        //    typedef typename YWrappedEvaluator::Evaluator YEvaluator;

        //    template<typename ArgumentVectorType>
        //    static double GetBeta(const ArgumentVectorType& args)
        //    {
        //        return 1.0;
        //    }
        //};

        template<typename L, typename R, typename IndicesType, unsigned int index>
        struct BetaYParameterAccessImpl< Node<L, expt::MultiplyOp, R>, IndicesType, index, 
            typename boost::enable_if
            <
                boost::mpl::and_    
                <
                    boost::is_same<typename L::ResultType, Nektar::NekVector<double> >,
                    boost::is_same<typename R::ResultType, double>
                >
            >::type> : public boost::true_type 
        {
            typedef DgemvNodeEvaluator<L, IndicesType, index> YWrappedEvaluator;
            typedef typename YWrappedEvaluator::Evaluator YEvaluator;

            static const unsigned int nextIndex = index + L::TotalCount;
            typedef DgemvNodeEvaluator<R, IndicesType, nextIndex> BetaWrappedEvaluator;
            typedef typename BetaWrappedEvaluator::Evaluator BetaEvaluator;
            
            template<typename ArgumentVectorType>
            static double GetBeta(const ArgumentVectorType& args)
            {
                return BetaEvaluator::Evaluate(args) * YWrappedEvaluator::GetScale();
            }
        };

        template<typename L, typename R, typename IndicesType, unsigned int index>
        struct BetaYParameterAccessImpl< Node<L, expt::MultiplyOp, R>, IndicesType, index, 
            typename boost::enable_if
            <
                boost::mpl::and_    
                <
                    boost::is_same<typename R::ResultType, Nektar::NekVector<double> >,
                    boost::is_same<typename L::ResultType, double>
                >
            >::type> : public boost::true_type 
        {
            typedef DgemvNodeEvaluator<L, IndicesType, index> BetaWrappedEvaluator;
            typedef typename BetaWrappedEvaluator::Evaluator BetaEvaluator;

            static const unsigned int nextIndex = index + L::TotalCount;
            typedef DgemvNodeEvaluator<R, IndicesType, nextIndex> YWrappedEvaluator;
            typedef typename YWrappedEvaluator::Evaluator YEvaluator;

            template<typename ArgumentVectorType>
            static double GetBeta(const ArgumentVectorType& args)
            {
                return BetaEvaluator::Evaluate(args) * YWrappedEvaluator::GetScale();
            }
        };

        template<typename NodeType, typename IndicesType, unsigned int index>
        struct BetaYParameterAccessImpl< NodeType, IndicesType, index, 
            typename boost::enable_if
            <
                boost::is_same<typename NodeType::ResultType, Nektar::NekVector<double> >
            >::type> : public boost::true_type 
        {
            typedef DgemvNodeEvaluator<NodeType, IndicesType, index> YWrappedEvaluator;
            typedef typename YWrappedEvaluator::Evaluator YEvaluator;

            template<typename ArgumentVectorType>
            static double GetBeta(const ArgumentVectorType& args)
            {
                return 1.0 * YWrappedEvaluator::GetScale();
            }
        };

        template<typename NodeType, typename IndicesType, unsigned int index>
        struct BetaYParameterAccess : public BetaYParameterAccessImpl<NodeType, IndicesType, index> {};

        template<typename T, typename IndicesType, unsigned int index>
        struct BetaYParameterAccess<Node<T, NegateOp, void>, IndicesType, index> : public BetaYParameterAccessImpl<T, IndicesType, index> 
        {
            template<typename ArgumentVectorType>
            static double GetBeta(const ArgumentVectorType& args)
            {
                return -BetaYParameterAccessImpl<T, IndicesType, index>::GetBeta(args);
            }
        };

    }

    // Cases to handle
    // 1. alpha*A*x +- beta*y
    // 2. beta*y +- alpha*A*x
    // 3. alpha*A*x + beta*B*y
    // 4. alpha*A*x

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
                impl::AlphaAXParameterAccess<L, IndicesType, index>,
                impl::BetaYParameterAccess<R, IndicesType, index + L::TotalCount>
                //IsAlphaABNode<L>,
                //IsBetaCNode<R>
            >
        >::type> : public boost::true_type
    {
        typedef Node<L, OpType, R> NodeType;
        typedef impl::AlphaAXParameterAccess<L, IndicesType, index> LhsAccess;
        static const unsigned int rhsIndex = index + L::TotalCount;
        typedef impl::BetaYParameterAccess<R, IndicesType, rhsIndex> RhsAccess;
        typedef typename LhsAccess::AEvaluator AEvaluator;
        typedef typename LhsAccess::XEvaluator XEvaluator;
        typedef typename RhsAccess::YEvaluator YEvaluator;

        template<typename ResultType, typename ArgumentVectorType>
        static void Evaluate(ResultType& accumulator, const ArgumentVectorType& args)
        {
            typename AEvaluator::ResultType a = AEvaluator::Evaluate(args);
            typename XEvaluator::ResultType x = XEvaluator::Evaluate(args);
            typename YEvaluator::ResultType y = YEvaluator::Evaluate(args);

            double alpha = LhsAccess::GetAlpha(args);
            double beta = RhsAccess::GetBeta(args);

            beta *= impl::BetaScale<OpType>::GetScale();
            Nektar::Dgemv(accumulator, alpha, a, x, beta, y);
        }
    };

    ////////////////////////////////////////////////////////////
    // Case 2 - beta*B +/- alpha*A*X
    ////////////////////////////////////////////////////////////
    template<typename L, typename OpType, typename R, typename IndicesType, unsigned int index>
    struct BinaryBinaryEvaluateNodeOverride<L, OpType, R, IndicesType, index,
        typename boost::enable_if
        <
            boost::mpl::and_
            <
                impl::AlphaAXParameterAccess<R, IndicesType, index + L::TotalCount>,
                impl::BetaYParameterAccess<L, IndicesType, index>,
                boost::mpl::not_<impl::AlphaAXParameterAccess<L, IndicesType, index> >
                //IsAlphaABNode<R>,
                //IsBetaCNode<L>,
                //boost::mpl::not_<IsAlphaABNode<L> >
            >
        >::type> : public boost::true_type
    {
        typedef Node<L, OpType, R> NodeType;

        typedef impl::BetaYParameterAccess<L, IndicesType, index> LhsAccess;
        static const unsigned int rhsIndex = index + L::TotalCount;
        typedef impl::AlphaAXParameterAccess<R, IndicesType, rhsIndex> RhsAccess;
        
        typedef typename RhsAccess::AEvaluator AEvaluator;
        typedef typename RhsAccess::XEvaluator XEvaluator;
        typedef typename LhsAccess::YEvaluator YEvaluator;

        template<typename ResultType, typename ArgumentVectorType>
        static void Evaluate(ResultType& accumulator, const ArgumentVectorType& args)
        {
            typename AEvaluator::ResultType a = AEvaluator::Evaluate(args);
            typename XEvaluator::ResultType x = XEvaluator::Evaluate(args);
            typename YEvaluator::ResultType y = YEvaluator::Evaluate(args);

            double alpha = RhsAccess::GetAlpha(args);
            double beta = LhsAccess::GetBeta(args);

            alpha *= impl::BetaScale<OpType>::GetScale();
            Nektar::Dgemv(accumulator, alpha, a, x, beta, y);
        }
    };

    ////////////////////////////////////////////////////////////
    // Case 4 - alpha*A*B
    ////////////////////////////////////////////////////////////
    template<typename L, typename OpType, typename R, typename IndicesType, unsigned int index>
    struct BinaryBinaryEvaluateNodeOverride<L, OpType, R, IndicesType, index,
        typename boost::enable_if
        <
            impl::AlphaAXParameterAccess<Node<L, OpType, R>, IndicesType, index>
            //IsAlphaABNode<Node<L, OpType, R> >
        >::type> : public boost::true_type
    {
        typedef Node<L, OpType, R> NodeType;
        typedef impl::AlphaAXParameterAccess<NodeType, IndicesType, index> LhsAccess;
        typedef typename LhsAccess::AEvaluator AEvaluator;
        typedef typename LhsAccess::XEvaluator XEvaluator;

        template<typename ResultType, typename ArgumentVectorType>
        static void Evaluate(ResultType& accumulator, const ArgumentVectorType& args)
        {
            typename AEvaluator::ResultType a = AEvaluator::Evaluate(args);
            typename XEvaluator::ResultType x = XEvaluator::Evaluate(args);

            double alpha = LhsAccess::GetAlpha(args);
            Nektar::Dgemm(accumulator, alpha, a, x);
        }
    };


}

#endif
#endif
