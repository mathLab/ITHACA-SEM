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

#ifndef EXPRESSION_TEMPLATES_SORT_ASSOCIATIVE_COMMUTATIVE_CLUSTERS_HPP
#define EXPRESSION_TEMPLATES_SORT_ASSOCIATIVE_COMMUTATIVE_CLUSTERS_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <ExpressionTemplates/Node.hpp>
#include <ExpressionTemplates/AssociativeTransform.hpp>
#include <ExpressionTemplates/CommutativeTransform.hpp>
#include <boost/mpl/comparison.hpp>

namespace expt
{

    /// \brief Finds the furthest left node of the tree rooted at this node.
    ///
    /// When optimizing AC clusters, we need to swap the furthest left node 
    /// and the node requiring the most temporaries to minimize the number of 
    /// temporaries required by the cluster.  This metafunction finds the 
    /// furthest left node.
    template<typename NodeType, typename ACOpType, typename enabled = void>
    struct FindFurthestLeftNode
    {
        typedef NodeType Type;
    };

    // Continue searching if it is a binary node on the left with the same operator type
    // and it is commutative.  
    template<typename L1, typename LOp, typename L2, typename ACOpType, typename RhsType>
    struct FindFurthestLeftNode<Node<Node<L1, LOp, L2>, ACOpType, RhsType>, ACOpType,
            typename boost::enable_if
            <
                boost::mpl::and_
                <
                    boost::is_same<LOp, ACOpType>,
                    CommutativeTraits<typename Node<L1, LOp, L2>::ResultType, ACOpType, typename RhsType::ResultType>
                >
            >::type>
    {
        typedef typename  FindFurthestLeftNode< Node<L1, ACOpType, L2>, ACOpType>::Type Type;
    };

    template<typename L1, typename LOp, typename L2, typename ACOpType, typename RhsType>
    struct FindFurthestLeftNode<Node<Node<L1, LOp, L2>, ACOpType, RhsType>, ACOpType,
            typename boost::disable_if
            <
                boost::mpl::and_
                <
                    boost::is_same<LOp, ACOpType>,
                    CommutativeTraits<typename Node<L1, LOp, L2>::ResultType, ACOpType, typename RhsType::ResultType>
                >
            >::type>
    {
        typedef Node<L1, LOp, L2> Type;
    };


    template<typename NodeType, typename ACOpType, unsigned int index = 0, typename enabled = void>
    struct FindNodeWithMostTemporaries;

    template<typename T, typename ACOpType, unsigned int StartingIndex>
    struct FindNodeWithMostTemporaries<Node<T, void, void>, ACOpType, StartingIndex>
    {
        typedef Node<T, void, void> Type;
        static const unsigned int Index = StartingIndex;
    };

    template<typename T, typename TOp, typename ACOpType, unsigned int StartingIndex>
    struct FindNodeWithMostTemporaries<Node<T, TOp, void>, ACOpType, StartingIndex>
    {
        typedef Node<T, TOp, void> Type;
        static const unsigned int Index = StartingIndex;
    };

    template<typename LhsType, typename RhsType, typename ACOpType, unsigned int StartingIndex>
    struct FindNodeWithMostTemporaries<Node<LhsType, ACOpType, RhsType>, ACOpType, StartingIndex,
            typename boost::enable_if
            <
                CommutativeTraits<typename LhsType::ResultType, ACOpType, typename RhsType::ResultType>
            >::type>
    {
        typedef typename FindNodeWithMostTemporaries<LhsType, ACOpType, StartingIndex>::Type LhsLargest;
        static const unsigned int LhsIndex = FindNodeWithMostTemporaries<LhsType, ACOpType, StartingIndex>::Index;

        static const unsigned int rhsStartingIndex = LhsType::TotalCount + StartingIndex;
        typedef typename FindNodeWithMostTemporaries<RhsType, ACOpType, rhsStartingIndex>::Type RhsLargest;
        static const unsigned int RhsIndex = FindNodeWithMostTemporaries<RhsType, ACOpType, rhsStartingIndex>::Index;

        static const unsigned int lhsTempCount = (LhsIndex == StartingIndex ? 0 :  TemporaryCount<LhsLargest>::AsRhs + 1);
        static const unsigned int rhsTempCount = TemporaryCount<RhsLargest>::AsRhs;

        static const unsigned int test = lhsTempCount >= rhsTempCount;
        typedef typename boost::mpl::if_c< test, LhsLargest, RhsLargest>::type Type;
        static const unsigned int Index = (test ? LhsIndex : RhsIndex);
    };

    template<typename LhsType, typename RhsType, typename ACOpType, unsigned int StartingIndex>
    struct FindNodeWithMostTemporaries<Node<LhsType, ACOpType, RhsType>, ACOpType, StartingIndex,
            typename boost::disable_if
            <
                CommutativeTraits<typename LhsType::ResultType, ACOpType, typename RhsType::ResultType>
            >::type>
    {
        typedef Node<LhsType, ACOpType, RhsType> Type;
        static const unsigned int Index = StartingIndex;
    };

    template<typename LhsType, typename OpType, typename RhsType, typename ACOpType, unsigned int StartingIndex>
    struct FindNodeWithMostTemporaries<Node<LhsType, OpType, RhsType>, ACOpType, StartingIndex,
            typename boost::enable_if
            <
                boost::mpl::and_
                <
                    IsBinaryNode<Node<LhsType, OpType, RhsType> >,
                    boost::mpl::not_<boost::is_same<OpType, ACOpType> >
                >
            >::type>
    {
        typedef Node<LhsType, OpType, RhsType> Type;
        static const unsigned int Index = StartingIndex;
    };

    namespace impl
    {
        template<typename InputSequence, unsigned int lbegin, unsigned int lend, unsigned int rbegin, unsigned int rend>
        struct SwapNodeIndices
        {
            typedef typename boost::mpl::begin<InputSequence>::type Begin;
            typedef Begin Iter0;
            typedef typename boost::mpl::advance_c< Begin, lbegin >::type Iter1;
            typedef typename boost::mpl::advance_c< Begin, lend>::type Iter2;
            typedef typename boost::mpl::advance_c<Begin, rbegin>::type Iter3;
            typedef typename boost::mpl::advance_c<Begin, rend>::type Iter4;
            typedef typename boost::mpl::end<InputSequence>::type Iter5;
            
            typedef typename boost::mpl::iterator_range<Iter0, Iter1>::type Range0;
            typedef typename boost::mpl::iterator_range<Iter1, Iter2>::type Range1;
            typedef typename boost::mpl::iterator_range<Iter2, Iter3>::type Range2;
            typedef typename boost::mpl::iterator_range<Iter3, Iter4>::type Range3;
            typedef typename boost::mpl::iterator_range<Iter4, Iter5>::type Range4;

            typedef typename boost::mpl::joint_view<Range0, Range3>::type View0;
            typedef typename boost::mpl::joint_view<Range2, Range1>::type View1;
            typedef typename boost::mpl::joint_view<View1, Range4>::type View2;

            typedef typename boost::mpl::joint_view<View0, View2>::type Type;
            
        };

        // Swaps two nodes in the tree.
        template<typename NodeType, typename TargetNode, typename ReplacementNode, unsigned int TargetIndex, unsigned int StartingIndex = 0, typename enabled=void>
        struct SwapNodes;
        
        template<typename NodeType, typename TargetNode, typename ReplacementNode, unsigned int TargetIndex, unsigned int StartingIndex>
        struct SwapNodes<NodeType, TargetNode, ReplacementNode, TargetIndex, StartingIndex, 
            typename boost::enable_if
            <
                boost::mpl::greater<boost::mpl::int_<StartingIndex>, boost::mpl::int_<TargetIndex> >
            >::type>
        {
            typedef NodeType Type;
        };

        // If the index is <= StartingIndex, and the node type doesn't match, then we need to recurse.
        template<typename T, typename TargetNode, typename ReplacementNode, unsigned int TargetIndex, unsigned int StartingIndex>
        struct SwapNodes<Node<T, void, void>, TargetNode, ReplacementNode, TargetIndex, StartingIndex, 
            typename boost::enable_if
            <
                boost::mpl::or_
                <
                    boost::mpl::less<boost::mpl::int_<StartingIndex>, boost::mpl::int_<TargetIndex> >,
                    boost::mpl::and_
                    <
                        boost::mpl::less_equal<boost::mpl::int_<StartingIndex>, boost::mpl::int_<TargetIndex> >,
                        boost::mpl::not_<boost::is_same<Node<T, void, void>, TargetNode> >
                    >
                >
            >::type>
        {
            typedef Node<T, void, void> NodeType;
            typedef NodeType Type;
        };

        template<typename T, typename TOp, typename TargetNode, typename ReplacementNode, unsigned int TargetIndex, unsigned int StartingIndex>
        struct SwapNodes<Node<T, TOp, void>, TargetNode, ReplacementNode, TargetIndex, StartingIndex, 
            typename boost::enable_if
            <
                boost::mpl::or_
                <
                    boost::mpl::less<boost::mpl::int_<StartingIndex>, boost::mpl::int_<TargetIndex> >,
                    boost::mpl::and_
                    <
                        boost::mpl::less_equal<boost::mpl::int_<StartingIndex>, boost::mpl::int_<TargetIndex> >,
                        boost::mpl::not_<boost::is_same<Node<T, TOp, void>, TargetNode> >
                    >
                >
            >::type>
        {
            typedef typename SwapNodes<T, TargetNode, ReplacementNode, TargetIndex, StartingIndex>::Type TransformedChildType;
            typedef Node<TransformedChildType, TOp, void> Type;
        };

        template<typename T1, typename TOp, typename T2, typename TargetNode, typename ReplacementNode, unsigned int TargetIndex, unsigned int StartingIndex>
        struct SwapNodes<Node<T1, TOp, T2>, TargetNode, ReplacementNode, TargetIndex, StartingIndex, 
            typename boost::enable_if
            <
                boost::mpl::and_
                <
                    IsBinaryNode<Node<T1, TOp, T2> >,
                    boost::mpl::or_
                    <
                        boost::mpl::less<boost::mpl::int_<StartingIndex>, boost::mpl::int_<TargetIndex> >,
                        boost::mpl::and_
                        <
                            boost::mpl::less_equal<boost::mpl::int_<StartingIndex>, boost::mpl::int_<TargetIndex> >,
                            boost::mpl::not_<boost::is_same<Node<T1, TOp, T2>, TargetNode> >
                        >
                    >
                >
            >::type>
        {
            typedef typename SwapNodes<T1, TargetNode, ReplacementNode, TargetIndex, StartingIndex>::Type TransformedT1;
            typedef typename SwapNodes<T2, TargetNode, ReplacementNode, TargetIndex, StartingIndex + T1::TotalCount>::Type TransformedT2;
            typedef Node<TransformedT1, TOp, TransformedT2> Type;
        };

        // If the index == starting index and the node types match, then we are good.
        template<typename NodeType, typename ReplacementNode, unsigned int TargetIndex>
        struct SwapNodes<NodeType, NodeType, ReplacementNode, TargetIndex, TargetIndex>
        {
            typedef ReplacementNode Type;
        };
    }

    // When done, this metafunction moves the node that requires the most temporaries to the first 
    // position in the AC cluster.  This is a toplevel metafunction that performs this function 
    // for all AC clusters in the tree.
    //
    // This metafunction assumes that the cluster is already in cannonical form.
    template<typename NodeType, typename IndicesType, unsigned int StartIndex, typename enabled=void>
    struct SwapACNodesIfNeeded;

    template<typename NodeType, typename IndicesType, unsigned int StartIndex, typename TestOpType, typename enabled=void>
    struct AdvanceToNextOperator;

    template<typename T, typename IndicesType, unsigned int StartIndex, typename TestOpType>
    struct AdvanceToNextOperator<Node<T, void, void>, IndicesType, StartIndex, TestOpType>
    {
        typedef Node<T, void, void> NodeType;
        typedef NodeType TransformedNodeType;
        typedef IndicesType TransformedIndicesType;
    };

    template<typename T, typename TOp, typename IndicesType, unsigned int StartIndex, typename TestOpType>
    struct AdvanceToNextOperator<Node<T, TOp, void>, IndicesType, StartIndex, TestOpType>
    {
        typedef Node<T, TOp, void> NodeType;
        typedef typename AdvanceToNextOperator<T, IndicesType, StartIndex, TestOpType>::TransformedNodeType TransformedChildType;
        typedef typename AdvanceToNextOperator<T, IndicesType, StartIndex, TestOpType>::TransformedIndicesType TransformedChildIndices;
        
        typedef Node<TransformedChildType, TOp, void> TransformedNodeType;
        typedef TransformedChildIndices TransformedIndicesType;
    };

    // If same, recurse
    template<typename LhsType, typename RhsType, typename IndicesType, unsigned int StartIndex, typename TestOpType>
    struct AdvanceToNextOperator<Node<LhsType, TestOpType, RhsType>, IndicesType, StartIndex, TestOpType>
    {
        typedef typename AdvanceToNextOperator<LhsType, IndicesType, StartIndex, TestOpType>::TransformedNodeType TransformedLhsNodeType;
        typedef typename AdvanceToNextOperator<LhsType, IndicesType, StartIndex, TestOpType>::TransformedIndicesType Indices0;

        static const int RhsIndex = LhsType::TotalCount + StartIndex;
        typedef typename AdvanceToNextOperator<RhsType, Indices0, RhsIndex, TestOpType>::TransformedNodeType TransformedRhsNodeType;
        typedef typename AdvanceToNextOperator<RhsType, Indices0, RhsIndex, TestOpType>::TransformedIndicesType Indices1;

        typedef Node<TransformedLhsNodeType, TestOpType, TransformedRhsNodeType> TransformedNodeType;
        typedef Indices1 TransformedIndicesType;

    };

    // If different, apply the Swap algorithm.
    template<typename LhsType, typename OpType, typename RhsType, typename IndicesType, unsigned int StartIndex, typename TestOpType>
    struct AdvanceToNextOperator<Node<LhsType, OpType, RhsType>, IndicesType, StartIndex, TestOpType>
    {
        typedef Node<LhsType, OpType, RhsType> ThisType;

        typedef typename SwapACNodesIfNeeded<ThisType, IndicesType, StartIndex>::TransformedNodeType TransformedNodeType;
        typedef typename SwapACNodesIfNeeded<ThisType, IndicesType, StartIndex>::TransformedIndicesType TransformedIndicesType;
    };

    

    template<typename T, typename IndicesType, unsigned int StartIndex>
    struct SwapACNodesIfNeeded<Node<T, void, void>, IndicesType, StartIndex>
    {
        typedef Node<T, void, void> NodeType;
        typedef NodeType TransformedNodeType;
        typedef IndicesType TransformedIndicesType;
    };

    template<typename T, typename TOp, typename IndicesType, unsigned int StartIndex>
    struct SwapACNodesIfNeeded<Node<T, TOp, void>, IndicesType, StartIndex>
    {
        typedef Node<T, TOp, void> NodeType;
        typedef typename SwapACNodesIfNeeded<T, IndicesType, StartIndex>::TransformedNodeType TransformedChildType;
        typedef typename SwapACNodesIfNeeded<T, IndicesType, StartIndex>::TransformedIndicesType TransformedChildIndices;
        
        typedef Node<TransformedChildType, TOp, void> TransformedNodeType;
        typedef TransformedChildIndices TransformedIndicesType;
    };

    namespace impl
    {
        template<typename NodeType, typename LhsNode, typename RhsNode, typename IndicesType, unsigned int LhsIndex, unsigned int RhsIndex, unsigned int StartingIndex, typename enabled=void>
        struct SwapLeftAndRightIfNeeded
        {
            typedef NodeType TransformedNodeType;
            typedef IndicesType TransformedIndicesType;
        };


        template<typename NodeType, typename LhsNode, typename RhsNode, typename IndicesType, unsigned int LhsIndex, unsigned int RhsIndex, unsigned int StartingIndex>
        struct SwapLeftAndRightIfNeeded<NodeType, LhsNode, RhsNode, IndicesType, LhsIndex, RhsIndex, StartingIndex,
            typename boost::enable_if
            <
                boost::mpl::and_
                <
                    boost::mpl::not_<boost::is_same< boost::mpl::int_<LhsIndex>, boost::mpl::int_<RhsIndex> > >,
                    IsConstantNode<LhsNode>,
                    boost::mpl::not_<IsConstantNode<RhsNode> >
                >
            >::type>
        {
            typedef typename SwapNodes<NodeType, RhsNode, LhsNode, RhsIndex, StartingIndex>::Type Type0;
            typedef typename SwapNodes<Type0, LhsNode, RhsNode, LhsIndex, StartingIndex>::Type TransformedNodeType;

            typedef typename SwapNodeIndices<IndicesType, LhsIndex, LhsIndex + LhsNode::TotalCount, RhsIndex, RhsIndex + RhsNode::TotalCount>::Type TransformedIndicesType;
        };

    }


    template<typename L1, typename OpType, typename RhsType, typename IndicesType, unsigned int StartIndex>
    struct SwapACNodesIfNeeded<Node<Node<L1, void, void>, OpType, RhsType>, IndicesType, StartIndex,
        typename boost::enable_if
        <
            IsBinaryNode<Node<Node<L1, void, void>, OpType, RhsType> >
        >::type>
    {
        static const int RhsIndex = StartIndex + 1;
        typedef typename SwapACNodesIfNeeded<RhsType, IndicesType, RhsIndex>::TransformedNodeType TransformedRhsType;
        typedef typename SwapACNodesIfNeeded<RhsType, IndicesType, RhsIndex>::TransformedIndicesType TransformedIndicesType;

        typedef Node<L1, void, void> LhsType;
        typedef Node<LhsType, OpType, RhsType> TransformedNodeType;
        
    };

    template<typename L1, typename LOp, typename OpType, typename RhsType, typename IndicesType, unsigned int StartIndex>
    struct SwapACNodesIfNeeded<Node<Node<L1, LOp, void>, OpType, RhsType>, IndicesType, StartIndex,
        typename boost::enable_if
        <
            IsBinaryNode<Node<Node<L1, LOp, void>, OpType, RhsType> >
        >::type>
    {
        typedef Node<L1, LOp, void> LhsType;
        static const int RhsIndex = StartIndex + LhsType::TotalCount;

        typedef typename SwapACNodesIfNeeded<RhsType, IndicesType, RhsIndex>::TransformedNodeType TransformedRhsType;
        typedef typename SwapACNodesIfNeeded<RhsType, IndicesType, RhsIndex>::TransformedIndicesType TransformedRhsIndices;

        typedef typename SwapACNodesIfNeeded<L1, TransformedRhsIndices, StartIndex>::TransformedNodeType TransformedLhsType;
        typedef typename SwapACNodesIfNeeded<L1, TransformedRhsIndices, StartIndex>::TransformedIndicesType TransformedIndicesType;

        typedef Node<Node<TransformedLhsType, LOp>, OpType, TransformedRhsType> TransformedNodeType;
        
    };

    template<typename L1, typename LOp, typename L2, typename OpType, typename RhsType, typename IndicesType, unsigned int StartIndex>
    struct SwapACNodesIfNeeded<Node<Node<L1, LOp, L2>, OpType, RhsType>, IndicesType, StartIndex,
        typename boost::enable_if
        <
            boost::mpl::and_
            <
                IsBinaryNode<Node<L1, LOp, L2> >,
                CommutativeTraits< typename Node<L1, LOp, L2>::ResultType, OpType, typename RhsType::ResultType>,
                AssociativeTraits< typename L1::ResultType, LOp, typename L2::ResultType, OpType>
            >
        >::type>
    {
        typedef Node<L1, LOp, L2> LhsType;
        typedef Node<LhsType, OpType, RhsType> ThisType;

        typedef typename FindFurthestLeftNode<ThisType, OpType>::Type LhsNode;
        // The index for the furthest lhs node will always be the start index of the entire tree.
        static const unsigned int LhsIndex = StartIndex;

        typedef typename FindNodeWithMostTemporaries<ThisType,OpType>::Type MostTemps;
        static const unsigned int RhsIndex = FindNodeWithMostTemporaries<ThisType, OpType, StartIndex>::Index;

        // The only way a temporary is avoided here is if the rhs is a non-constant node and the 
        // lhs is a constant node.  Anything else just re-arranges the tree for no gain.
        //
        // Either way, once the swap is or is not made, this ac cluster is complete, so we can 
        // traverse the tree until we find a different operator before starting up again.  We could
        // just traverse with this metafunction, but we would end up swapping a lot of internal nodes to no overall gain.
        typedef typename impl::SwapLeftAndRightIfNeeded<ThisType, LhsNode, MostTemps, IndicesType, LhsIndex, RhsIndex, StartIndex>::TransformedNodeType SwappedNodeType;
        typedef typename impl::SwapLeftAndRightIfNeeded<ThisType, LhsNode, MostTemps, IndicesType, LhsIndex, RhsIndex, StartIndex>::TransformedIndicesType SwappedIndicesType;

        typedef typename AdvanceToNextOperator<SwappedNodeType, SwappedIndicesType, StartIndex, OpType>::TransformedNodeType TransformedNodeType;
        typedef typename AdvanceToNextOperator<SwappedNodeType, SwappedIndicesType, StartIndex, OpType>::TransformedIndicesType TransformedIndicesType;

    };

    template<typename L1, typename LOp, typename L2, typename OpType, typename RhsType, typename IndicesType, unsigned int StartIndex>
    struct SwapACNodesIfNeeded<Node<Node<L1, LOp, L2>, OpType, RhsType>, IndicesType, StartIndex,
        typename boost::disable_if
        <
            boost::mpl::and_
            <
                IsBinaryNode<Node<L1, LOp, L2> >,
                CommutativeTraits< typename Node<L1, LOp, L2>::ResultType, OpType, typename RhsType::ResultType>,
                AssociativeTraits< typename L1::ResultType, LOp, typename L2::ResultType, OpType>
            >
        >::type>
    {
        typedef Node<L1, LOp, L2> LhsType;
        typedef typename SwapACNodesIfNeeded<LhsType, IndicesType, StartIndex>::TransformedNodeType TransformedLhsType;
        typedef typename SwapACNodesIfNeeded<LhsType, IndicesType, StartIndex>::TransformedIndicesType Indices0;

        static const unsigned int RhsStartIndex = StartIndex + LhsType::TotalCount;
        typedef typename SwapACNodesIfNeeded<RhsType, Indices0, RhsStartIndex>::TransformedNodeType TransformedRhsType;
        typedef typename SwapACNodesIfNeeded<RhsType, Indices0, RhsStartIndex>::TransformedIndicesType TransformedIndicesType;
        typedef Node<TransformedLhsType, OpType, TransformedRhsType> TransformedNodeType;
    };


}

#endif
#endif
