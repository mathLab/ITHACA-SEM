///////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_EXPRESSION_TEMPLATES_EXPRESSION_EVALUATOR_HPP
#define NEKTAR_EXPRESSION_TEMPLATES_EXPRESSION_EVALUATOR_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include "Node.hpp"
#include "RemoveAllUnecessaryTemporaries.hpp"
#include <boost/version.hpp>
#include <boost/fusion/algorithm/iteration/accumulate.hpp>
#include <boost/fusion/include/accumulate.hpp>


namespace Nektar
{
    template<typename DataType, typename NodeType, typename Indices, unsigned int StartIndex>
    struct CreateFromTree
    {
        template<typename ArgVectorType>
        static DataType Apply(const ArgVectorType& tree)
        {
            return DataType();
        }
    };

    // Evaluate node - all temporaries that can be avoided have been.
    // Here are the rules.
    //
    // Constant Node.
    // The only time we evaluate a constant node on its own is when going down the 
    // left side of the evaluation, so we just assign it.  
    //
    // Unary Node
    template<typename NodeType, typename IndicesType, unsigned int index>
    struct EvaluateNode;
    
    
    //////////////////////////////////////////////
    // Constant 
    //////////////////////////////////////////////
    template<typename Type, typename IndicesType, unsigned int index>
    struct EvaluateNode<Node<Type, void, void>, IndicesType, index>
    {
        static const unsigned int MappedIndex = boost::mpl::at_c<IndicesType, index>::type::value;
        template<typename ResultType, typename ArgumentVectorType>
        static void Evaluate(ResultType& accumulator, const ArgumentVectorType& args)
        {
            // Will there ever be a case where accumulator and args[index] aren't the 
            // same type and it isn't caught before here?
            accumulator = boost::fusion::at_c<MappedIndex>(args);
        }
    };
    


    //////////////////////////////////////////////////////////////////
    // Binary Nodes
    //////////////////////////////////////////////////////////////////
    // - Right constant - no temporary.
    // - Right unary or binary - temporary for the right.
    // - Left result type is different than accumulator - temporary for left.
    // 
    // Four combinations.
    // Left No Temp, Right no temp.  A = L, OpEqual(A, R).  Left is same result type, right is constant
    // Left temp, right no temp.  T = L, A = Op(T, R).  Left is not same rsult type, right is constant
    // Left no temp, right temp.  T = R, A = L, OpEqual(A, R).  Left is same rsult type, right is not constant
    // Left temp, right temp. T1 = L, t2 = r, A = Op(T1, t2).  Left is not same result type, right is not constant.
    

    template<typename NodeType, typename IndicesType, unsigned int index>
    struct EvaluateNodeWithTemporaryIfNeeded;
    
    template<typename Type, typename IndicesType, unsigned int index>
    struct EvaluateNodeWithTemporaryIfNeeded<Node<Type, void, void>, IndicesType, index>
    {
        typedef const Type& ResultType;
        static const unsigned int MappedIndex = boost::mpl::at_c<IndicesType, index>::type::value;

        template<typename ArgumentVectorType>
        static ResultType Evaluate(const ArgumentVectorType& args)
        {
            return boost::fusion::at_c<MappedIndex>(args);
        }
    };
    
    template<typename Arg1, typename Op, typename Arg2, typename IndicesType, unsigned int index>
    struct EvaluateNodeWithTemporaryIfNeeded<Node<Arg1, Op, Arg2>, IndicesType, index>
    {
        typedef Node<Arg1, Op, Arg2> NodeType;
        typedef typename NodeType::ResultType ResultType;
        
        template<typename ArgumentVectorType>
        static ResultType Evaluate(const ArgumentVectorType& args)
        {
            ResultType temp = CreateFromTree<ResultType, NodeType, IndicesType, index>::Apply(args);
            EvaluateNode<NodeType, IndicesType, index>::Evaluate(temp, args);
            return temp;
        }
    };

    
    // Binary specialization when both children are leaf nodes.  Using the general purpose
    // algorithm results in unecessary temporaries.
    template<typename L, typename Op, typename R, typename IndicesType, unsigned int index>
    struct EvaluateNode<Node<Node<L, void, void>, Op, Node<R, void, void> >, IndicesType, index>
    {
        static const unsigned int LhsMappedIndex = boost::mpl::at_c<IndicesType, index>::type::value;
        static const unsigned int RhsMappedIndex = boost::mpl::at_c<IndicesType, index+1>::type::value;

        template<typename ResultType, typename ArgumentVectorType>
        static void Evaluate(ResultType& accumulator, const ArgumentVectorType& args)
        {
            Op::Op(accumulator, boost::fusion::at_c<LhsMappedIndex>(args),
                boost::fusion::at_c<RhsMappedIndex>(args));

            // Doesn't compile for m*v
            //accumulator = boost::fusion::at_c<LhsMappedIndex>(args);
            //Op::OpEqual(accumulator, boost::fusion::at_c<RhsMappedIndex>(args));

            // Causes a temp.
            //accumulator = Op::Op(boost::fusion::at_c<LhsMappedIndex>(args),
            //    boost::fusion::at_c<RhsMappedIndex>(args));
        }
    };

    // To override, provide a specialized version of this class.
    template<typename LhsType, typename Op, typename RhsType, typename IndicesType, unsigned int index, typename enabled = void>
    struct BinaryBinaryEvaluateNodeOverride : public boost::false_type {};

    template<typename LhsType, typename Op, typename RhsType, typename IndicesType, unsigned int index>
    struct EvaluateNode<Node<LhsType, Op, RhsType >, IndicesType, index>
    {
        static const int rhsNodeIndex = index + LhsType::TotalCount;
        typedef typename LhsType::ResultType LhsResultType;
        typedef EvaluateNodeWithTemporaryIfNeeded<RhsType, IndicesType, rhsNodeIndex> EvaluateNodeWithTemporaryIfNeededType;
        typedef typename EvaluateNodeWithTemporaryIfNeededType::ResultType RhsTempType;    
    
        typedef Node<LhsType, Op, RhsType> Expression;
        
        template<typename ResultType, typename ArgumentVectorType>
        static void Evaluate(ResultType& accumulator, const ArgumentVectorType& args, 
            typename boost::enable_if
            <
                boost::mpl::and_
                <
                    boost::mpl::not_<BinaryBinaryEvaluateNodeOverride<LhsType, Op, RhsType, IndicesType, index> >,
                    boost::is_same<ResultType, LhsResultType>
                >
            >::type* dummy = 0)
        {
            EvaluateNode<LhsType, IndicesType, index>::Evaluate(accumulator, args);
            RhsTempType rhs = EvaluateNodeWithTemporaryIfNeededType::Evaluate(args);
            Op::OpEqual(accumulator, rhs);
        }
        
        template<typename ResultType, typename ArgumentVectorType>
        static void Evaluate(ResultType& accumulator, const ArgumentVectorType& args, 
            typename boost::enable_if
            <
                boost::mpl::and_
                <
                    boost::mpl::not_<BinaryBinaryEvaluateNodeOverride<LhsType, Op, RhsType, IndicesType, index> >,
                    boost::mpl::not_<boost::is_same<ResultType, LhsResultType> >
                >
            >::type* dummy = 0)
        {
            LhsResultType temp = CreateFromTree<LhsResultType, LhsType, IndicesType, index>::Apply(args);
            EvaluateNode<LhsType, IndicesType, index>::Evaluate(temp, args);
            RhsTempType rhs = EvaluateNodeWithTemporaryIfNeededType::Evaluate(args);
            accumulator = Op::Op(temp, rhs);
        }

        template<typename ResultType, typename ArgumentVectorType>
        static void Evaluate(ResultType& accumulator, const ArgumentVectorType& args, 
            typename boost::enable_if
            <
                boost::mpl::and_
                <
                    BinaryBinaryEvaluateNodeOverride<LhsType, Op, RhsType, IndicesType, index>,

                    // Note that this condition should always be true, so we don't include a 
                    // version where ResultType and Expression::ResultType aren't the same, but we
                    // need it so that some parameter of enable_if relies on the function template 
                    // parameters.
                    boost::is_same<ResultType, typename Expression::ResultType>
                >
            >::type* dummy = 0)
        {
            return BinaryBinaryEvaluateNodeOverride<LhsType, Op, RhsType, IndicesType, index>::Evaluate(accumulator, args);
        }
    };

    template<typename T, typename R, typename enabled = void>
    struct IsAlias
    {
        static bool Apply(const T& lhs, const R& rhs)
        {
            return false;
        }
    };

    template<typename T>
    struct IsAlias<T, T>
    {
        static bool Apply(const T& lhs, const T& rhs)
        {
            return &lhs == &rhs;
        }
    };

    struct ExpressionEvaluator
    {
        template<typename T>
        struct ContainsAliasAccumulator
        {
            ContainsAliasAccumulator(const T& v) :
                value(v)
            {
            }

#if BOOST_VERSION < 104200
            template<typename ElementType>
            unsigned int operator()(const ElementType& rhs, const unsigned int& accum) const
#else
            template<typename ElementType>
            unsigned int operator()(const unsigned int& accum, const ElementType& rhs) const
#endif
            {
                return accum + IsAlias<T, ElementType>::Apply(value, rhs);
            }

            typedef int result_type;
            typedef int result;

            const T& value;
        };

        template<typename Expression>
        static bool ContainsAlias(const Expression& expression, typename Expression::ResultType& accum)
        {
            ContainsAliasAccumulator<typename Expression::ResultType> containsAliasAccumulator(accum);
            int numAliases = boost::fusion::accumulate(expression.GetData(), 0, containsAliasAccumulator);
            return numAliases > 0;
        }

        template<typename Expression>
        static void Evaluate(const Expression& expression, typename Expression::ResultType& accum)
        {
            typedef typename Expression::Indices Indices;

            // Perform the optimizations on the parse three.
            typedef typename RemoveUnecessaryTemporaries<Expression, Indices>::TransformedNodeType OptimizedParseTree;
            typedef typename RemoveUnecessaryTemporaries<Expression, Indices>::TransformedIndices TransformedIndices;

            if( ContainsAlias(expression, accum) )
            {
                typename Expression::ResultType temp = CreateFromTree<typename Expression::ResultType, OptimizedParseTree, TransformedIndices, 0>::Apply(expression.GetData());
                EvaluateNode<OptimizedParseTree, TransformedIndices, 0>::Evaluate(temp, expression.GetData());
                accum = temp;
            }
            else
            {
                EvaluateNode<OptimizedParseTree, TransformedIndices, 0>::Evaluate(accum, expression.GetData());
            }
        } 
    };
}

#endif
#endif //NEKTAR_EXPRESSION_TEMPLATES_EXPRESSION_EVALUATOR_HPP

