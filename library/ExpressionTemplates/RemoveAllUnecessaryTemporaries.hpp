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

#ifndef EXPRESSION_TEMPLATES_REMOVE_ALL_UNECESSARY_TEMPORARIES_H
#define EXPRESSION_TEMPLATES_REMOVE_ALL_UNECESSARY_TEMPORARIES_H
#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <ExpressionTemplates/AssociativeTraits.hpp>
#include <ExpressionTemplates/CommutativeTraits.hpp>
#include <ExpressionTemplates/CommutativeTransform.hpp>
#include <ExpressionTemplates/AssociativeTransform.hpp>
#include <ExpressionTemplates/PushDownUnaryNodes.hpp>
#include <ExpressionTemplates/ForwardInverseTransform.hpp>
#include <ExpressionTemplates/BackwardInverseTransform.hpp>
#include <ExpressionTemplates/SortAssociativeCommutativeClusters.hpp>
#include <ExpressionTemplates/PerformCommutativeTransformIfNeeded.hpp>

namespace expt
{
    namespace impl
    {
        // Updates the parse tree so there are no unecessary temporaries.
        template<typename NodeType, typename IndicesType, unsigned int StartIndex = 0, typename enabled = void>
        struct RemoveUnecessaryTemporariesInternal;

        // Constant node (Base Case 0) - don't need to do anything.
        template<typename T, typename IndicesType, unsigned int StartIndex>
        struct RemoveUnecessaryTemporariesInternal<expt::Node<T, void, void>, IndicesType, StartIndex>
        {
            typedef expt::Node<T, void, void> TransformedNodeType;
            typedef IndicesType TransformedIndicesType;
        };

        // Unary node 
        template<typename T, typename UnaryOp, typename IndicesType, unsigned int StartIndex>
        struct RemoveUnecessaryTemporariesInternal<expt::Node<T, UnaryOp, void>, IndicesType, StartIndex>
        {          
            typedef typename RemoveUnecessaryTemporariesInternal<T, IndicesType, StartIndex>::TransformedNodeType ChildNodeType;
            typedef typename RemoveUnecessaryTemporariesInternal<T, IndicesType, StartIndex>::TransformedIndicesType ChildIndicesType;

            typedef expt::Node<ChildNodeType, UnaryOp, void> TransformedNodeType;
            typedef ChildIndicesType TransformedIndicesType;
        };
    

        // Binary Node
        template<typename L, typename OpType, typename R, typename IndicesType, unsigned int StartIndex>
        struct RemoveUnecessaryTemporariesInternal<expt::Node<L, OpType, R>, IndicesType, StartIndex>
        {
            static const int RightStart = StartIndex + L::TotalCount;

            // Optimize the right child.
            typedef typename RemoveUnecessaryTemporariesInternal<R, IndicesType, RightStart>::TransformedNodeType RightNode0Type;
            typedef typename RemoveUnecessaryTemporariesInternal<R, IndicesType, RightStart>::TransformedIndicesType Indices0Type;
            typedef typename expt::Node<L, OpType, RightNode0Type> Tree0Type;

            // If Op and ROp are associative, apply the associative transformation.
            // Note that the specializations for the AssociativeTransform object only apply when the right 
            // child is a binary node with the appropriate operator, so this transformation will not mistakenly
            // turn a constant right child into a binary right child.
            typedef typename AssociativeTransform<Tree0Type>::TransformedNodeType Tree1Type;
            typedef typename Tree1Type::Left LeftNode1Type;
            typedef Indices0Type Indices1Type;

            // Optimize the left child.
            typedef typename RemoveUnecessaryTemporariesInternal<LeftNode1Type, Indices1Type, StartIndex>::TransformedNodeType LeftNode2Type;
            typedef typename RemoveUnecessaryTemporariesInternal<LeftNode1Type, Indices1Type, StartIndex>::TransformedIndicesType Indices2Type;
            typedef expt::Node<LeftNode2Type, OpType, typename Tree1Type::Right> Tree2Type;

            // Apply commutative if Op is commutative and the number of temps for the rhs is greater 
            // than the number for the left.
            typedef typename PerformCommutativeTransformIfNeeded<Tree2Type, Indices2Type, StartIndex>::TransformedNodeType Tree3Type;
            typedef typename PerformCommutativeTransformIfNeeded<Tree2Type, Indices2Type, StartIndex>::TransformedIndicesType Indices3Type;

            typedef Tree3Type TransformedNodeType;
            typedef Indices3Type TransformedIndicesType;
        };
    }

    template<typename NodeType>
    struct RemoveUnecessaryTemporaries
    {
        typedef typename NodeType::Indices IndicesType;

        // Step 1 - Iterate the tree and evaluate all unary nodes that can be potentially simplified.  It may not 
        // reduce the number of nodes in the tree, but it does push the unary nodes down towards the leaves of the 
        // parse tree, making it easier to group operators.
        typedef typename PushDownUnaryNodes<NodeType>::Type T0;
        typedef IndicesType I0;

        // Step 2 - Now that we have pushed unary operators towards the leaves, 
        // do a backward associative transform to catch cases such as 
        // A + (-B) -> A-B.
        typedef typename BackwardInverseTransform<T0, I0, 0>::TransformedNodeType T1;
        typedef typename BackwardInverseTransform<T0, I0, 0>::TransformedIndicesType I1;

        // Now forward transform.
        typedef typename ForwardInverseTransform<T1>::Type T2;
        typedef I1 I2;

        // Now we can run the optimization procedures.
        typedef typename impl::RemoveUnecessaryTemporariesInternal<T2, I2, 0>::TransformedNodeType Tree0Type;
        typedef typename impl::RemoveUnecessaryTemporariesInternal<T2, I2, 0>::TransformedIndicesType Indices0Type;

        typedef typename SwapACNodesIfNeeded<Tree0Type, Indices0Type, 0>::TransformedNodeType Tree1Type;
        typedef typename SwapACNodesIfNeeded<Tree0Type, Indices0Type, 0>::TransformedIndicesType Indices1Type;

        typedef typename BackwardInverseTransform<Tree1Type, Indices1Type, 0>::TransformedNodeType Tree2Type;
        typedef typename BackwardInverseTransform<Tree1Type, Indices1Type, 0>::TransformedIndicesType Indices2Type;

        typedef Tree2Type TransformedNodeType;
        typedef Indices2Type TransformedIndicesType;
    };
}

#endif
#endif

