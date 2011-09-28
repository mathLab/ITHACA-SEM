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
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_IS_DGEMM_TRAITS_HPP
#define NEKTAR_LIB_UTILITIES_IS_DGEMM_TRAITS_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

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
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>

#include <string>

namespace Nektar
{
    template<typename NodeType, typename enabled = void>
    struct IsDgemmRightSide : public boost::false_type {};

    template<typename L, typename R>
    struct IsDgemmRightSide<expt::Node<L, expt::MultiplyOp, R>,
        typename boost::enable_if
        <
            boost::mpl::and_
            <
                boost::is_same<typename L::ResultType, double>,
                CanGetRawPtr<typename R::ResultType >
            >
        >::type > : public boost::true_type
    {
        template<typename IndicesType, unsigned int index>
        struct GetValues
        {
            typedef expt::EvaluateNodeWithTemporaryIfNeeded<L, IndicesType, index> LEvaluator;

            static const unsigned int RIndex = index + L::TotalCount;
            typedef expt::EvaluateNodeWithTemporaryIfNeeded<R, IndicesType, RIndex> REvaluator;

            template<typename ArgumentVectorType>
            static double GetBeta(const ArgumentVectorType& args)
            {
                return LEvaluator::Evaluate(args);
            }

            template<typename ArgumentVectorType>
            static typename REvaluator::ResultType GetC(const ArgumentVectorType& args)
            {
                return REvaluator::Evaluate(args);
            }
        };
        
    };

    template<typename L, typename R>
    struct IsDgemmRightSide<expt::Node<L, expt::MultiplyOp, R>,
        typename boost::enable_if
        <
            boost::mpl::and_
            <
                boost::is_same<typename R::ResultType, double>,
                CanGetRawPtr<typename L::ResultType >
            >
        >::type > : public boost::true_type
    {
        template<typename IndicesType, unsigned int index>
        struct GetValues
        {
            typedef expt::EvaluateNodeWithTemporaryIfNeeded<L, IndicesType, index> LEvaluator;

            static const unsigned int RIndex = index + L::TotalCount;
            typedef expt::EvaluateNodeWithTemporaryIfNeeded<R, IndicesType, RIndex> REvaluator;

            template<typename ArgumentVectorType>
            static double GetBeta(const ArgumentVectorType& args)
            {
                return REvaluator::Evaluate(args);
            }

            template<typename ArgumentVectorType>
            static typename REvaluator::ResultType GetC(const ArgumentVectorType& args)
            {
                return LEvaluator::Evaluate(args);
            }
        };
        
    };

    template<typename T, typename Op>
    struct IsDgemmRightSide<expt::Node<T, Op, void> ,
        typename boost::enable_if
        < CanGetRawPtr<typename expt::Node<T, Op, void>::ResultType > >::type > : public boost::true_type
    {
        template<typename IndicesType, unsigned int index>
        struct GetValues
        {
            typedef expt::EvaluateNodeWithTemporaryIfNeeded< expt::Node<T, Op, void>, IndicesType, index> Evaluator;

            template<typename ArgumentVectorType>
            static double GetBeta(const ArgumentVectorType& args)
            {
                return 1.0;
            }

            template<typename ArgumentVectorType>
            static typename Evaluator::ResultType GetC(const ArgumentVectorType& args)
            {
                return Evaluator::Evaluate(args);
            }
        };
        
    };

    template<typename NodeType, typename enabled = void>
    struct IsDgemmLeftSide : public boost::false_type {};

    template<typename A1, typename A2, typename A3>
    struct IsDgemmLeftSide< expt::Node< expt::Node<A1, expt::MultiplyOp, A2>, expt::MultiplyOp, A3> ,
        typename boost::enable_if
        <
            boost::mpl::and_
            <
                boost::is_same<double, typename A1::ResultType>,
                CanGetRawPtr<typename A2::ResultType>,
                CanGetRawPtr<typename A3::ResultType>
            >
        >::type > : public boost::true_type
    {
        template<typename IndicesType, unsigned int index>
        struct GetValues
        {
            typedef expt::EvaluateNodeWithTemporaryIfNeeded<A1, IndicesType, index> A1Evaluator;

            static const unsigned int A2Index = index + A1::TotalCount;
            typedef expt::EvaluateNodeWithTemporaryIfNeeded<A2, IndicesType, A2Index> A2Evaluator;

            static const unsigned int A3Index = A2Index + A2::TotalCount;
            typedef expt::EvaluateNodeWithTemporaryIfNeeded<A3, IndicesType, A3Index> A3Evaluator;

            template<typename ArgumentVectorType>
            static typename A1Evaluator::ResultType GetAlpha(const ArgumentVectorType& args)
            {
                return A1Evaluator::Evaluate(args);
            }

            template<typename ArgumentVectorType>
            static typename A2Evaluator::ResultType GetA(const ArgumentVectorType& args)
            {
                return A2Evaluator::Evaluate(args);
            }

            template<typename ArgumentVectorType>
            static typename A3Evaluator::ResultType GetB(const ArgumentVectorType& args)
            {
                return A3Evaluator::Evaluate(args);
            }
        };
    };

    template<typename A1, typename A2, typename A3>
    struct IsDgemmLeftSide< expt::Node< expt::Node<A1, expt::MultiplyOp, A2>, expt::MultiplyOp, A3> ,
        typename boost::enable_if
        <
            boost::mpl::and_
            <
                boost::is_same<double, typename A2::ResultType>,
                CanGetRawPtr<typename A1::ResultType>,
                CanGetRawPtr<typename A3::ResultType>
            >
        >::type > : public boost::true_type
    {
        template<typename IndicesType, unsigned int index>
        struct GetValues
        {
            typedef expt::EvaluateNodeWithTemporaryIfNeeded<A1, IndicesType, index> A1Evaluator;

            static const unsigned int A2Index = index + A1::TotalCount;
            typedef expt::EvaluateNodeWithTemporaryIfNeeded<A2, IndicesType, A2Index> A2Evaluator;

            static const unsigned int A3Index = A2Index + A2::TotalCount;
            typedef expt::EvaluateNodeWithTemporaryIfNeeded<A3, IndicesType, A3Index> A3Evaluator;

            template<typename ArgumentVectorType>
            static typename A2Evaluator::ResultType GetAlpha(const ArgumentVectorType& args)
            {
                return A2Evaluator::Evaluate(args);
            }

            template<typename ArgumentVectorType>
            static typename A1Evaluator::ResultType GetA(const ArgumentVectorType& args)
            {
                return A1Evaluator::Evaluate(args);
            }

            template<typename ArgumentVectorType>
            static typename A3Evaluator::ResultType GetB(const ArgumentVectorType& args)
            {
                return A3Evaluator::Evaluate(args);
            }
        };
    };

    template<typename A1, typename A2, typename A3>
    struct IsDgemmLeftSide< expt::Node< expt::Node<A1, expt::MultiplyOp, A2>, expt::MultiplyOp, A3> ,
        typename boost::enable_if
        <
            boost::mpl::and_
            <
                boost::is_same<double, typename A3::ResultType>,
                CanGetRawPtr<typename A1::ResultType>,
                CanGetRawPtr<typename A2::ResultType>
            >
        >::type > : public boost::true_type
    {
        template<typename IndicesType, unsigned int index>
        struct GetValues
        {
            typedef expt::EvaluateNodeWithTemporaryIfNeeded<A1, IndicesType, index> A1Evaluator;

            static const unsigned int A2Index = index + A1::TotalCount;
            typedef expt::EvaluateNodeWithTemporaryIfNeeded<A2, IndicesType, A2Index> A2Evaluator;

            static const unsigned int A3Index = A2Index + A2::TotalCount;
            typedef expt::EvaluateNodeWithTemporaryIfNeeded<A3, IndicesType, A3Index> A3Evaluator;

            template<typename ArgumentVectorType>
            static typename A3Evaluator::ResultType GetAlpha(const ArgumentVectorType& args)
            {
                return A3Evaluator::Evaluate(args);
            }

            template<typename ArgumentVectorType>
            static typename A1Evaluator::ResultType GetA(const ArgumentVectorType& args)
            {
                return A1Evaluator::Evaluate(args);
            }

            template<typename ArgumentVectorType>
            static typename A2Evaluator::ResultType GetB(const ArgumentVectorType& args)
            {
                return A2Evaluator::Evaluate(args);
            }
        };
    };


    // Three matrices multiplied together constitutes a left hand side of 
    // a dgemm call (one subexpression will need to be evaluated first).
    template<typename A1, typename A2, typename A3>
    struct IsDgemmLeftSide< expt::Node< expt::Node<A1, expt::MultiplyOp, A2>, expt::MultiplyOp, A3> ,
        typename boost::enable_if
        <
            boost::mpl::and_
            <
                CanGetRawPtr<typename expt::Node<A1, expt::MultiplyOp, A2>::ResultType>,
                CanGetRawPtr<typename A3::ResultType>,
                boost::mpl::not_<boost::is_same<typename A1::ResultType, double> >,
                boost::mpl::not_<boost::is_same<typename A2::ResultType, double> >
            >
        >::type > : public boost::true_type
    {
        template<typename IndicesType, unsigned int index>
        struct GetValues
        {
            typedef expt::Node<A1, expt::MultiplyOp, A2> LhsNode;
            typedef expt::EvaluateNodeWithTemporaryIfNeeded<LhsNode, IndicesType, index> A1Evaluator;

            static const unsigned int A3Index = index + LhsNode::TotalCount;
            typedef expt::EvaluateNodeWithTemporaryIfNeeded<A3, IndicesType, A3Index> A3Evaluator;

            template<typename ArgumentVectorType>
            static double GetAlpha(const ArgumentVectorType& args)
            {
                return 1.0;
            }

            template<typename ArgumentVectorType>
            static typename A1Evaluator::ResultType GetA(const ArgumentVectorType& args)
            {
                return A1Evaluator::Evaluate(args);
            }

            template<typename ArgumentVectorType>
            static typename A3Evaluator::ResultType GetB(const ArgumentVectorType& args)
            {
                return A3Evaluator::Evaluate(args);
            }
        };
    };

    template<typename T, typename R>
    struct IsDgemmLeftSide< expt::Node< expt::Node<T, void, void>, expt::MultiplyOp, expt::Node<R, void, void> > ,
        typename boost::enable_if
        <
            boost::mpl::and_
            <
                CanGetRawPtr<typename expt::Node<T, void, void>::ResultType>,
                CanGetRawPtr<typename expt::Node<R, void, void>::ResultType>
            >
        >::type > : public boost::true_type
    {
        template<typename IndicesType, unsigned int index>
        struct GetValues
        {
            typedef expt::Node<T, void, void> A1;
            typedef expt::Node<R, void, void> A2;

            typedef expt::EvaluateNodeWithTemporaryIfNeeded<A1, IndicesType, index> A1Evaluator;

            static const unsigned int A2Index = index + A1::TotalCount;
            typedef expt::EvaluateNodeWithTemporaryIfNeeded<A2, IndicesType, A2Index> A2Evaluator;

            template<typename ArgumentVectorType>
            static double GetAlpha(const ArgumentVectorType& args)
            {
                return 1.0;
            }

            template<typename ArgumentVectorType>
            static typename A1Evaluator::ResultType GetA(const ArgumentVectorType& args)
            {
                return A1Evaluator::Evaluate(args);
            }

            template<typename ArgumentVectorType>
            static typename A2Evaluator::ResultType GetB(const ArgumentVectorType& args)
            {
                return A2Evaluator::Evaluate(args);
            }
        };
    };

}

#endif
#endif