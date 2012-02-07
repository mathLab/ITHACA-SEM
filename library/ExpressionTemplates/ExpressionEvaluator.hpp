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

#ifndef EXPRESSION_TEMPLATES_EXPRESSION_EVALUATOR_HPP
#define EXPRESSION_TEMPLATES_EXPRESSION_EVALUATOR_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <ExpressionTemplates/Node.hpp>
#include <ExpressionTemplates/RemoveAllUnecessaryTemporaries.hpp>
#include <ExpressionTemplates/CreateFromTree.hpp>

#include <boost/fusion/algorithm/iteration/accumulate.hpp>
#include <boost/fusion/include/accumulate.hpp>
#include <boost/version.hpp>

namespace expt
{
    /// \brief EvaluateNode is responsible for evaluating the tree.
    ///
    /// The behavior of this class is different for each type of operator type,
    /// so all of the work is performed in the specializations below.
    template<typename NodeType, typename IndicesType, unsigned int index=0>
    struct EvaluateNode;
    
    
    //////////////////////////////////////////////
    // Constant 
    //////////////////////////////////////////////

    /// \brief Evaluates a constant node.
    template<typename Type, typename IndicesType, unsigned int index>
    struct EvaluateNode<Node<Type, void, void>, IndicesType, index>
    {
        static const unsigned int MappedIndex = boost::mpl::at_c<IndicesType, index>::type::value;

        template<typename ResultType, typename ArgumentVectorType>
        static void Evaluate(ResultType& accumulator, const ArgumentVectorType& args)
        {
            accumulator = boost::fusion::at_c<MappedIndex>(args);
        }
    };
    
    //////////////////////////////////////////////
    // Unary 
    //////////////////////////////////////////////
    template<typename ChildType, typename OpType, typename IndicesType, unsigned int index>
    struct EvaluateNode<Node<ChildType, OpType, void>, IndicesType, index>
    {
        template<typename ResultType, typename ArgumentVectorType>
        static void Evaluate(ResultType& accumulator, const ArgumentVectorType& args)
        {
            // Evaluate the child, then apply the unary operation.
            EvaluateNode<ChildType, IndicesType, index>::Evaluate(accumulator, args);
            OpType::Op(accumulator);
        }
    };

    //////////////////////////////////////////////////////////////////
    // Binary Nodes
    //////////////////////////////////////////////////////////////////
    namespace impl
    {
        // During binary node evaluation, there are time when the only difference between 
        // two evaluation strategies is whether or not the right child is constant and 
        // does not need a temporary, or is a unary/binary node and does need a temporary.
        //
        // This class takes care of these two scenarios.  
        template<typename NodeType, typename IndicesType, unsigned int index>
        struct EvaluateNodeWithTemporaryIfNeeded;
        
        // Constant nodes don't need temporaries, so we can just return a reference.
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
        
        // Anything else requires the construction of a temporary.
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
    }
    

    // To override the default behavior during binary node evaluation, create a specialization of 
    // this class with an Evaluate method defined as below.
    template<typename LhsType, typename Op, typename RhsType, typename IndicesType, unsigned int index, typename enabled = void>
    struct BinaryBinaryEvaluateNodeOverride : public boost::false_type {};

    namespace impl
    {
        template<typename LhsType, typename Op, typename RhsType, typename IndicesType, unsigned int index, typename enabled=void>
        struct EvaluateBinaryNode
        {
            static const int rhsNodeIndex = index + LhsType::TotalCount;
            typedef typename LhsType::ResultType LhsResultType;
            typedef EvaluateNodeWithTemporaryIfNeeded<RhsType, IndicesType, rhsNodeIndex> EvaluateNodeWithTemporaryIfNeededType;
            typedef typename EvaluateNodeWithTemporaryIfNeededType::ResultType RhsTempType;    
        
            typedef Node<LhsType, Op, RhsType> Expression;
            
            // No override and the lhs and result type are the same, so the accumulator 
            // can be passed to the lhs.
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
            
            // No override and the lhs is not the same type as the result, so we need 
            // to create a temporary to evaluate the lhs.
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
                Op::Op(accumulator, temp, rhs);
            }

            // There is a user supplied override that will be called.
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

        // Binary, two constant children.  A+B
        template<typename L, typename Op, typename R, typename IndicesType, unsigned int index>
        struct EvaluateBinaryNode<Node<L, void, void>, Op, Node<R, void, void>, IndicesType, index>
        {
            static const unsigned int LhsMappedIndex = boost::mpl::at_c<IndicesType, index>::type::value;
            static const unsigned int RhsMappedIndex = boost::mpl::at_c<IndicesType, index+1>::type::value;

            template<typename ResultType, typename ArgumentVectorType>
            static void Evaluate(ResultType& accumulator, const ArgumentVectorType& args)
            {
                Op::Op(accumulator, boost::fusion::at_c<LhsMappedIndex>(args),
                    boost::fusion::at_c<RhsMappedIndex>(args));
            }
        };
    }

    template<typename LhsType, typename OpType, typename RhsType, typename IndicesType, unsigned int index>
    struct EvaluateNode<Node<LhsType, OpType, RhsType>, IndicesType, index>
    {
        template<typename ResultType, typename ArgumentVectorType>
        static void Evaluate(ResultType& accumulator, const ArgumentVectorType& args)
        {
            impl::EvaluateBinaryNode<LhsType, OpType, RhsType, IndicesType, index>::Evaluate(accumulator, args);
        }
    };

    // During evaluation, we need to detect scenarios such as 
    // A = A*B + A
    //
    // If we evaluate A*B first, then the accumulator value will be overwritten and the 
    // expression will be incorrect.  When an alias is detected, the expression evaluator
    // creates a temporary for the accumulator.
    //
    // Users can specialize this if the default behaviors are not adequate.
    template<typename T, typename R>
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
        /// \brief A class that is used with boost::fusion::accumulate to iterate an expression's
        /// argument list to determine if there are any aliases.
        template<typename T>
        struct ContainsAliasAccumulator
        {
            ContainsAliasAccumulator(const T& v) :
                value(v)
            {
            }

            // The ordering of the parameters passed to this class in boost::fusion::accumulate
            // changed in boost 1.42
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

        /// \brief Determines if the accumulator is aliased anywhere in the expression.
        template<typename Expression>
        static bool ContainsAlias(const Expression& expression, typename Expression::ResultType& accum)
        {
            ContainsAliasAccumulator<typename Expression::ResultType> containsAliasAccumulator(accum);
            int numAliases = boost::fusion::accumulate(expression.GetData(), 0, containsAliasAccumulator);
            return numAliases > 0;
        }

        /// \brief Evaluates an expression, creating a new object as the accumulator.
        ///
        /// It is not always possible to modify the class definitions of the objects for which 
        /// expression templates will be used.  In those cases, user code can call this evaluate 
        /// method to evaluate the expression.
        template<typename Expression>
        static typename Expression::ResultType Evaluate(const Expression& expression)
        {
            typedef typename Expression::Indices Indices;

            // Perform the optimizations on the parse three.
            typedef typename RemoveUnecessaryTemporaries<Expression>::TransformedNodeType OptimizedParseTree;
            typedef typename RemoveUnecessaryTemporaries<Expression>::TransformedIndicesType TransformedIndicesType;

            typename Expression::ResultType result = CreateFromTree<typename Expression::ResultType, OptimizedParseTree, TransformedIndicesType, 0>::Apply(expression.GetData());
            EvaluateNode<OptimizedParseTree, TransformedIndicesType>::Evaluate(result, expression.GetData());
            return result;
        }

        /// \brief This method evaluates the expression, storing the results in accum.
        ///
        /// It is expected that this method is called from within a class constructor or assignment operator, 
        /// with accum = *this.
        template<typename Expression>
        static void Evaluate(const Expression& expression, typename Expression::ResultType& accum)
        {
            typedef typename Expression::Indices Indices;

            // Perform the optimizations on the parse three.
            typedef typename RemoveUnecessaryTemporaries<Expression>::TransformedNodeType OptimizedParseTree;
            typedef typename RemoveUnecessaryTemporaries<Expression>::TransformedIndicesType TransformedIndicesType;

            if( ContainsAlias(expression, accum) )
            {
                typename Expression::ResultType temp = CreateFromTree<typename Expression::ResultType, OptimizedParseTree, TransformedIndicesType, 0>::Apply(expression.GetData());
                EvaluateNode<OptimizedParseTree, TransformedIndicesType>::Evaluate(temp, expression.GetData());
                accum = temp;
            }
            else
            {
                EvaluateNode<OptimizedParseTree, TransformedIndicesType>::Evaluate(accum, expression.GetData());
            }
        } 

        template<typename Expression>
        static void EvaluateWithoutAliasingCheck(const Expression& expression, typename Expression::ResultType& accum)
        {
            typedef typename Expression::Indices Indices;

            // Perform the optimizations on the parse three.
            typedef typename RemoveUnecessaryTemporaries<Expression>::TransformedNodeType OptimizedParseTree;
            typedef typename RemoveUnecessaryTemporaries<Expression>::TransformedIndicesType TransformedIndicesType;

            EvaluateNode<OptimizedParseTree, TransformedIndicesType>::Evaluate(accum, expression.GetData());
        } 
    };
}

#endif
#endif //EXPRESSION_TEMPLATES_EXPRESSION_EVALUATOR_HPP

