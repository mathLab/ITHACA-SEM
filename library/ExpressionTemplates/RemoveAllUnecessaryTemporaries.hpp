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

#ifndef NEKTAR_EXPRESSION_TEMPLATES_REMOVE_ALL_UNECESSARY_TEMPORARIES_H
#define NEKTAR_EXPRESSION_TEMPLATES_REMOVE_ALL_UNECESSARY_TEMPORARIES_H

#include <ExpressionTemplates/AssociativeTraits.hpp>
#include <ExpressionTemplates/CommutativeTraits.hpp>
#include <ExpressionTemplates/CommutativeTransform.hpp>
#include <ExpressionTemplates/AssociativeTransform.hpp>

namespace Nektar
{
    template<typename NodeType, typename IndicesType, unsigned int StartIndex = 0>
    class RemoveUnecessaryTemporaries;

    template<typename NodeType, typename IndicesType, unsigned int StartIndex>
    class Default
    {
        private:
            template<typename InnerNodeType, typename InnerIndicesType, unsigned int InnerStartIndex>
            struct Apply;

            template<typename T, typename InnerIndicesType, unsigned int InnerStartIndex>
            struct Apply<Node<T, void, void>, InnerIndicesType, InnerStartIndex>
            {
                typedef Node<T, void, void> TransformedNodeType;
                typedef InnerIndicesType TransformedIndices;
            };

            template<typename T, typename Op, typename InnerIndicesType, unsigned int InnerStartIndex>
            struct Apply<Node<T, Op, void>, InnerIndicesType, InnerStartIndex>
            {
                // TODO - still need to finish the unary cases.
            };

            template<typename L, typename Op, typename R, typename InnerIndicesType, unsigned int InnerStartIndex>
            struct Apply<Node<L, Op, R>, InnerIndicesType, InnerStartIndex>
            {
                typedef typename RemoveUnecessaryTemporaries<L, InnerIndicesType, InnerStartIndex>::TransformedNodeType OptimizedLeft;
                typedef typename RemoveUnecessaryTemporaries<L, InnerIndicesType, InnerStartIndex>::TransformedIndices FirstOptimizedIndices;

                static const int RightStart = InnerStartIndex + L::TotalCount;
                typedef typename RemoveUnecessaryTemporaries<R, FirstOptimizedIndices, RightStart>::TransformedNodeType OptimizedRight;
                typedef typename RemoveUnecessaryTemporaries<R, FirstOptimizedIndices, RightStart>::TransformedIndices SecondOptimizedIndices;

                typedef Node<OptimizedLeft, Op, OptimizedRight> TransformedNodeType;
                typedef SecondOptimizedIndices TransformedIndices;
            };

        public:

            typedef typename Apply<NodeType, IndicesType, StartIndex>::TransformedNodeType TransformedNodeType;
            typedef typename Apply<NodeType, IndicesType, StartIndex>::TransformedIndices TransformedIndices;
    };

    // RemoveUnecessaryTemporaries is a parse tree, compile time optimizer.
    template<typename NodeType, typename IndicesType, unsigned int StartIndex>
    class RemoveUnecessaryTemporaries : public Default<NodeType, IndicesType, StartIndex>
    {
    };

    // Constant-Constant specialization.  
    template<typename LeftType, typename Op, typename RightType, typename IndicesType, unsigned int StartIndex>
    class RemoveUnecessaryTemporaries<Node<Node<LeftType, void, void>, Op, Node<RightType, void, void> >, IndicesType, StartIndex > :
        public Default<Node<Node<LeftType, void, void>, Op, Node<RightType, void, void> >, IndicesType, StartIndex > {};

    // Left anything, right constant.  Nothing do do.
    template<typename L1, typename LOp, typename L2, typename Op, typename R, typename Indices, unsigned int StartIndex>
    class RemoveUnecessaryTemporaries<Node<Node<L1, LOp, L2>, Op, Node<R, void, void> >, Indices, StartIndex > :
        public Default<Node<Node<L1, LOp, L2>, Op, Node<R, void, void> >, Indices, StartIndex > {};

    // Case 1 - Left constant - Right anything.
    template<typename LeftType, typename Op, typename RightArg1, typename RightOp, typename RightArg2, typename IndicesType, unsigned int StartIndex>
    class RemoveUnecessaryTemporaries<Node<Node<LeftType, void, void>, Op, Node<RightArg1, RightOp, RightArg2> >, IndicesType, StartIndex >
    {
        private:
            template<typename NodeType, typename enabled = void>
            struct Apply : public Default<NodeType, IndicesType, StartIndex> {};

            template<typename LeftNode, typename InnerOp, typename RightNode>
            struct Apply<Node<LeftNode, InnerOp, RightNode>,
                    typename boost::enable_if
                        <
                            CommutativeTraits<typename LeftNode::ResultType, InnerOp, typename RightNode::ResultType>
                        >::type>
            {
                typedef Node<LeftNode, InnerOp, RightNode> InnerNodeType;
                typedef CommutativeTransform<InnerNodeType, IndicesType, StartIndex> Transform;

                typedef typename Transform::TransformedNodeType CommutativeNodeType;
                typedef typename Transform::TransformedIndicesType CommutativeIndices;

                typedef typename RemoveUnecessaryTemporaries<CommutativeNodeType, CommutativeIndices, StartIndex>::TransformedNodeType TransformedNodeType;
                typedef typename RemoveUnecessaryTemporaries<CommutativeNodeType, CommutativeIndices, StartIndex>::TransformedIndices TransformedIndices;
            };
        
        public:
            typedef Node<Node<LeftType, void, void>, Op, Node<RightArg1, RightOp, RightArg2> > NodeType;
            typedef typename Apply<NodeType>::TransformedNodeType TransformedNodeType;
            typedef typename Apply<NodeType>::TransformedIndices TransformedIndices;
    };

    // Case 2 - Left and Right nodes are arbitrary.
    template<typename L1, typename LOp, typename L2, typename Op, typename R1, typename ROp, typename R2, typename Indices, unsigned int StartIndex>
    struct RemoveUnecessaryTemporaries<Node<Node<L1, LOp, L2>, Op, Node<R1, ROp, R2> >, Indices, StartIndex >
    {
        template<typename NodeType, typename InnerIndicesType, typename enabled = void>
        struct Apply : public Default<NodeType, InnerIndicesType, StartIndex> {};

        // case 2.1 - If the node is associative, then apply the associative transformation to get a new 
        // tree with a new root and re-apply to the new root.
        template<typename InnerLhsNode, typename InnerOpType, typename InnerR1, typename InnerROpType, typename InnerR2, typename InnerIndicesType>
        struct Apply<Node<InnerLhsNode, InnerOpType, Node<InnerR1, InnerROpType, InnerR2> >, InnerIndicesType,
                typename boost::enable_if
                    <
                        TreeIsAssociative<Node<InnerLhsNode, InnerOpType, Node<InnerR1, InnerROpType, InnerR2> > >
                    >::type>
        {
            typedef Node<InnerLhsNode, InnerOpType, Node<InnerR1, InnerROpType, InnerR2> > NodeType;
            typedef typename AssociativeTransform<NodeType>::TransformedNodeType NewTree;

            typedef typename RemoveUnecessaryTemporaries<NewTree, InnerIndicesType, StartIndex>::TransformedNodeType TransformedNodeType;
            typedef typename RemoveUnecessaryTemporaries<NewTree, InnerIndicesType, StartIndex>::TransformedIndices TransformedIndices;
        };

        // case 2.2 - If not associative, but commutative, and the swapped tree is associative, 
        // do a commutative, followed by an associative.
        template<typename InnerL1, typename InnerLOp, typename InnerL2, typename InnerOp, typename InnerR1, typename InnerROp, typename InnerR2, typename InnerIndicesType>
        struct Apply<Node<Node<InnerL1, InnerLOp, InnerL2>, InnerOp, Node<InnerR1, InnerROp, InnerR2> >, InnerIndicesType,
            typename boost::enable_if
            <
                boost::mpl::and_
                <
                    boost::mpl::not_<TreeIsAssociative<Node<Node<InnerL1, InnerLOp, InnerL2>, InnerOp, Node<InnerR1, InnerROp, InnerR2> > > >,
                    CommutativeTraits<typename Node<InnerL1, InnerLOp, InnerL2>::ResultType, InnerOp, typename Node<InnerR1, InnerROp, InnerR2>::ResultType>,
                    TreeIsAssociative<Node<Node<InnerL1, InnerLOp, InnerL2>, InnerOp, Node<InnerR1, InnerROp, InnerR2> > >
                >
            >::type >
        {
            typedef Node<Node<InnerL1, InnerLOp, InnerL2>, InnerOp, Node<InnerR1, InnerROp, InnerR2> > InnerNodeType;
            typedef typename CommutativeTransform<InnerNodeType, Indices, StartIndex>::TransformedNodeType CommutativeTransformedNodeType;
            typedef typename CommutativeTransform<InnerNodeType, Indices, StartIndex>::TransformedIndicesType CommutativeIndicesType;

            typedef typename AssociativeTransform<CommutativeTransformedNodeType>::TransformedNodeType AssociativeTransformNodeType;

            typedef typename RemoveUnecessaryTemporaries<AssociativeTransformNodeType, CommutativeIndicesType, StartIndex>::TransformedNodeType TransformedNodeType;
            typedef typename RemoveUnecessaryTemporaries<AssociativeTransformNodeType, CommutativeIndicesType, StartIndex>::TransformedIndices TransformedIndices;
        };

        // Step 1 - Optimize both side independently.
        typedef Node<L1, LOp, L2> LhsNodeType;
        typedef Node<R1, ROp, R2> RhsNodeType;
        
        typedef typename RemoveUnecessaryTemporaries<LhsNodeType, Indices, StartIndex>::TransformedNodeType TransformedLhsNodeType;
        typedef typename RemoveUnecessaryTemporaries<LhsNodeType, Indices, StartIndex>::TransformedIndices TransformedIndicesFirstPass;

        static const int RightStart = StartIndex + LhsNodeType::TotalCount;
        typedef typename RemoveUnecessaryTemporaries<RhsNodeType, TransformedIndicesFirstPass, RightStart>::TransformedNodeType TransformedRhsNodeType;
        typedef typename RemoveUnecessaryTemporaries<RhsNodeType, TransformedIndicesFirstPass, RightStart>::TransformedIndices TransformedIndicesSecondPass;

        // Step 2 - Check if we can avoid a temporary through a full associative tree transform, and do 
        // it if we can.
        typedef Node<TransformedLhsNodeType, Op, TransformedRhsNodeType > NodeType;
        typedef typename Apply<NodeType, TransformedIndicesSecondPass>::TransformedNodeType TransformedNodeType;
        typedef typename Apply<NodeType, TransformedIndicesSecondPass>::TransformedIndices TransformedIndices;
    };

}

#endif //NEKTAR_EXPRESSION_TEMPLATES_REMOVE_ALL_UNECESSARY_TEMPORARIES_H
