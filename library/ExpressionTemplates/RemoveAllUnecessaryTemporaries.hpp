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

#include <ExpressionTemplates/AssociativeTraits.hpp>
#include <ExpressionTemplates/CommutativeTraits.hpp>
#include <ExpressionTemplates/CommutativeTransform.hpp>
#include <ExpressionTemplates/AssociativeTransform.hpp>

namespace expt
{
    // Input is a node in an AssociativeCommutative cluster, and we are determining if we should swap nodes.
    template<typename NodeType, typename IndicesType, unsigned int StartIndex, typename enabled=void>
    struct SortAssociativeCommutativeClusters
    {
        static const unsigned int SpecializationId = 0;

        typedef NodeType TransformedNodeType;
        typedef IndicesType TransformedIndicesType;
    };

    template<typename L1, typename LOp, typename L2, typename OpType, typename R1, typename ROp, typename R2,
             typename IndicesType, unsigned int StartIndex>
    struct SortAssociativeCommutativeClusters< expt::Node< expt::Node<L1, LOp, L2>, OpType, expt::Node<R1, ROp, R2> >, IndicesType, StartIndex,
        typename boost::enable_if
         <
            boost::mpl::and_
            <
                 boost::mpl::and_
                 <
                     boost::is_same<LOp, OpType>,
                     boost::mpl::not_<boost::is_same<ROp, void> >,
                     boost::mpl::not_<boost::is_same<ROp, OpType> >,
                     boost::mpl::less<boost::mpl::int_<expt::TemporaryCount<L2>::Value>, boost::mpl::int_<expt::TemporaryCount<expt::Node<R1, ROp, R2> >::Value+1> >
                 >,
                 boost::mpl::and_
                 <
                    CommutativeTraits< typename Node<L1, LOp, L2>::ResultType, OpType, typename Node<R1, ROp, R2>::ResultType>,
                    AssociativeTraits< typename L1::ResultType, LOp, typename L2::ResultType, OpType>
                 >
            >
         >::type>
    {
        static const unsigned int SpecializationId = 1;

        // Example tree: A+B+CD
        typedef expt::Node<expt::Node<L1, LOp, L2>, OpType, expt::Node<R1, ROp, R2> > Tree0;

        // Inverse associative so we have A+(B+CD)
        typedef typename InverseAssociativeTransform<Tree0>::TransformedNodeType Tree1;
        typedef typename Tree1::Right RightTree1;
        typedef typename Tree1::Left LeftTree1;
        static const int RightStart1 = StartIndex + LeftTree1::TotalCount;

        // Commutative transform on the right so we have A+(CD+B)
        typedef typename expt::CommutativeTransform<RightTree1, IndicesType, RightStart1>::TransformedNodeType RightTree2;
        typedef typename expt::CommutativeTransform<RightTree1, IndicesType, RightStart1>::TransformedIndicesType Indices2;
        typedef typename expt::Node<LeftTree1, typename Tree1::OpType, RightTree2> Tree2;


        // Associative transform to get A + CD + B
        typedef typename AssociativeTransform<Tree2>::TransformedNodeType Tree3;
        typedef Indices2 Indices3;

        // Recurse down the left side.
        typedef typename SortAssociativeCommutativeClusters<typename Tree3::Left, Indices3, StartIndex>::TransformedNodeType LeftTree4;
        typedef typename SortAssociativeCommutativeClusters<typename Tree3::Left, Indices3, StartIndex>::TransformedIndicesType Indices4;

        typedef expt::Node<LeftTree4, OpType, typename Tree3::Right> TransformedNodeType;
        typedef Indices4 TransformedIndicesType;
    };

    template<typename L1, typename LOp, typename L2, typename OpType, typename R1, typename ROp, typename R2,
             typename IndicesType, unsigned int StartIndex>
    struct SortAssociativeCommutativeClusters< expt::Node< expt::Node<L1, LOp, L2>, OpType, expt::Node<R1, ROp, R2> >, IndicesType, StartIndex,
        typename boost::enable_if
         <
             boost::mpl::and_
             <
                CommutativeTraits< typename Node<L1, LOp, L2>::ResultType, OpType, typename Node<R1, ROp, R2>::ResultType>,
                boost::is_same<LOp, void>,
                boost::mpl::not_<boost::is_same<ROp, void> >
             >
         >::type>
    {
        static const unsigned int SpecializationId = 2;
        typedef expt::Node<expt::Node<L1, LOp, L2>, OpType, expt::Node<R1, ROp, R2> > Tree0Type;

        typedef typename CommutativeTransform<Tree0Type, IndicesType, StartIndex>::TransformedNodeType Tree1Type;
        typedef typename CommutativeTransform<Tree0Type, IndicesType, StartIndex>::TransformedIndicesType Indices1Type;

        typedef Tree1Type TransformedNodeType;
        typedef Indices1Type TransformedIndicesType;
    };

    // Updates the parse tree so there are not unecessary temporaries.
    template<typename NodeType, typename IndicesType, unsigned int StartIndex = 0, typename enabled = void>
    struct RemoveUnecessaryTemporariesInternal;

    // Constant node (Base Case 0) - don't need to do anything.
    template<typename T, typename IndicesType, unsigned int StartIndex>
    struct RemoveUnecessaryTemporariesInternal<expt::Node<T, void, void>, IndicesType, StartIndex>
    {
        static const unsigned int SpecializationId = 1;
        typedef expt::Node<T, void, void> TransformedNodeType;
        typedef IndicesType TransformedIndicesType;
    };

    // Unary node (Base Case 1) - don't need to do anything.
    template<typename T, typename UnaryOp, typename IndicesType, unsigned int StartIndex>
    struct RemoveUnecessaryTemporariesInternal<expt::Node<T, UnaryOp, void>, IndicesType, StartIndex>
    {
        static const unsigned int SpecializationId = 2;
        typedef expt::Node<T, UnaryOp, void> TransformedNodeType;
        typedef IndicesType TransformedIndicesType;
    };

    template<typename NodeType, typename IndicesType, unsigned int StartIndex, typename enabled = void>
    struct PerformCommutativeTransformIfNeeded
    {
        typedef NodeType TransformedNodeType;
        typedef IndicesType TransformedIndicesType;
    };
    
    template<typename LhsType, typename OpType, typename RhsType, typename IndicesType, unsigned int StartIndex>
    struct PerformCommutativeTransformIfNeeded< expt::Node<LhsType, OpType, RhsType>, IndicesType, StartIndex,
    typename boost::enable_if
    <
    boost::mpl::and_
    <
    expt::CommutativeTraits<typename LhsType::ResultType, OpType, typename RhsType::ResultType>,
    boost::is_same<typename LhsType::OpType, void>,
    boost::mpl::not_<boost::is_same<typename RhsType::OpType, void> >
    //boost::mpl::less<boost::mpl::int_<expt::TemporaryCount<expt::Node<RhsType, OpType, LhsType> >::Value>, boost::mpl::int_<expt::TemporaryCount<expt::Node<LhsType, OpType, RhsType> >::Value> >
    >
    >::type>
    {
        typedef expt::Node<LhsType, OpType, RhsType> NodeType;
        typedef typename expt::CommutativeTransform<NodeType, IndicesType, StartIndex>::TransformedNodeType TransformedNodeType;
        typedef typename expt::CommutativeTransform<NodeType, IndicesType, StartIndex>::TransformedIndicesType TransformedIndicesType;
    };
    
//    template<typename NodeType, typename IndicesType, unsigned int StartIndex, typename enabled = void>
//    struct PerformAssociativerCommutativeTransformIfNeeded
//    {
//        typedef NodeType TransformedNodeType;
//        typedef IndicesType TransformedIndicesType;
//    };
//    
//    struct PerformAssociativerCommutativeTransformIfNeeded
//    {
//        typedef NodeType TransformedNodeType;
//        typedef IndicesType TransformedIndicesType;
//    };


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

    template<typename NodeType, typename IndicesType, unsigned int StartIndex = 0, typename enabled = void>
    struct RemoveUnecessaryTemporaries
    {
        typedef typename RemoveUnecessaryTemporariesInternal<NodeType, IndicesType, StartIndex>::TransformedNodeType Tree0Type;
        typedef typename RemoveUnecessaryTemporariesInternal<NodeType, IndicesType, StartIndex>::TransformedIndicesType Indices0Type;

        //typedef typename SortAssociativeCommutativeClusters<Tree0Type, Indices0Type, StartIndex>::TransformedNodeType Tree1Type;
        //typedef typename SortAssociativeCommutativeClusters<Tree0Type, Indices0Type, StartIndex>::TransformedIndicesType Indices1Type;

        typedef Tree0Type TransformedNodeType;
        typedef Indices0Type TransformedIndicesType;
    };


}

namespace Nektar
{
    //template<typename NodeType, typename IndicesType, unsigned int StartIndex = 0>
    //class RemoveUnecessaryTemporaries;

    //template<typename NodeType, typename IndicesType, unsigned int StartIndex>
    //class Default
    //{
    //    private:
    //        template<typename InnerNodeType, typename InnerIndicesType, unsigned int InnerStartIndex>
    //        struct Apply;

    //        template<typename T, typename InnerIndicesType, unsigned int InnerStartIndex>
    //        struct Apply<Node<T, void, void>, InnerIndicesType, InnerStartIndex>
    //        {
    //            typedef Node<T, void, void> TransformedNodeType;
    //            typedef InnerIndicesType TransformedIndicesType;
    //        };

    //        template<typename T, typename Op, typename InnerIndicesType, unsigned int InnerStartIndex>
    //        struct Apply<Node<T, Op, void>, InnerIndicesType, InnerStartIndex>
    //        {
    //            // TODO - still need to finish the unary cases.
    //        };

    //        template<typename L, typename Op, typename R, typename InnerIndicesType, unsigned int InnerStartIndex>
    //        struct Apply<Node<L, Op, R>, InnerIndicesType, InnerStartIndex>
    //        {
    //            typedef typename RemoveUnecessaryTemporaries<L, InnerIndicesType, InnerStartIndex>::TransformedNodeType OptimizedLeft;
    //            typedef typename RemoveUnecessaryTemporaries<L, InnerIndicesType, InnerStartIndex>::TransformedIndicesType FirstOptimizedIndices;

    //            static const int RightStart = InnerStartIndex + L::TotalCount;
    //            typedef typename RemoveUnecessaryTemporaries<R, FirstOptimizedIndices, RightStart>::TransformedNodeType OptimizedRight;
    //            typedef typename RemoveUnecessaryTemporaries<R, FirstOptimizedIndices, RightStart>::TransformedIndicesType SecondOptimizedIndices;

    //            typedef Node<OptimizedLeft, Op, OptimizedRight> TransformedNodeType;
    //            typedef SecondOptimizedIndices TransformedIndicesType;
    //        };

    //    public:

    //        typedef typename Apply<NodeType, IndicesType, StartIndex>::TransformedNodeType TransformedNodeType;
    //        typedef typename Apply<NodeType, IndicesType, StartIndex>::TransformedIndicesType TransformedIndicesType;
    //};

    //// RemoveUnecessaryTemporaries is a parse tree, compile time optimizer.
    //template<typename NodeType, typename IndicesType, unsigned int StartIndex>
    //class RemoveUnecessaryTemporaries : public Default<NodeType, IndicesType, StartIndex>
    //{
    //};

    //// Constant-Constant specialization.  
    //template<typename LeftType, typename Op, typename RightType, typename IndicesType, unsigned int StartIndex>
    //class RemoveUnecessaryTemporaries<Node<Node<LeftType, void, void>, Op, Node<RightType, void, void> >, IndicesType, StartIndex > :
    //    public Default<Node<Node<LeftType, void, void>, Op, Node<RightType, void, void> >, IndicesType, StartIndex > {};

    //// Left anything, right constant.  Nothing do do.
    //template<typename L1, typename LOp, typename L2, typename Op, typename R, typename Indices, unsigned int StartIndex>
    //class RemoveUnecessaryTemporaries<Node<Node<L1, LOp, L2>, Op, Node<R, void, void> >, Indices, StartIndex > :
    //    public Default<Node<Node<L1, LOp, L2>, Op, Node<R, void, void> >, Indices, StartIndex > {};

    //// Case 1 - Left constant - Right anything.
    //template<typename LeftType, typename Op, typename RightArg1, typename RightOp, typename RightArg2, typename IndicesType, unsigned int StartIndex>
    //class RemoveUnecessaryTemporaries<Node<Node<LeftType, void, void>, Op, Node<RightArg1, RightOp, RightArg2> >, IndicesType, StartIndex >
    //{
    //    private:
    //        template<typename NodeType, typename enabled = void>
    //        struct Apply : public Default<NodeType, IndicesType, StartIndex> {};

    //        template<typename LeftNode, typename InnerOp, typename RightNode>
    //        struct Apply<Node<LeftNode, InnerOp, RightNode>,
    //                typename boost::enable_if
    //                    <
    //                        CommutativeTraits<typename LeftNode::ResultType, InnerOp, typename RightNode::ResultType>
    //                    >::type>
    //        {
    //            typedef Node<LeftNode, InnerOp, RightNode> InnerNodeType;
    //            typedef CommutativeTransform<InnerNodeType, IndicesType, StartIndex> Transform;

    //            typedef typename Transform::TransformedNodeType CommutativeNodeType;
    //            typedef typename Transform::TransformedIndicesType CommutativeIndices;

    //            typedef typename RemoveUnecessaryTemporaries<CommutativeNodeType, CommutativeIndices, StartIndex>::TransformedNodeType TransformedNodeType;
    //            typedef typename RemoveUnecessaryTemporaries<CommutativeNodeType, CommutativeIndices, StartIndex>::TransformedIndicesType TransformedIndicesType;
    //        };
    //    
    //    public:
    //        typedef Node<Node<LeftType, void, void>, Op, Node<RightArg1, RightOp, RightArg2> > NodeType;
    //        typedef typename Apply<NodeType>::TransformedNodeType TransformedNodeType;
    //        typedef typename Apply<NodeType>::TransformedIndicesType TransformedIndicesType;
    //};

    //// Case 2 - Left and Right nodes are arbitrary.
    //template<typename L1, typename LOp, typename L2, typename Op, typename R1, typename ROp, typename R2, typename Indices, unsigned int StartIndex>
    //struct RemoveUnecessaryTemporaries<Node<Node<L1, LOp, L2>, Op, Node<R1, ROp, R2> >, Indices, StartIndex >
    //{
    //    template<typename NodeType, typename InnerIndicesType, typename enabled = void>
    //    struct Apply : public Default<NodeType, InnerIndicesType, StartIndex> {};

    //    // case 2.1 - If the node is associative, then apply the associative transformation to get a new 
    //    // tree with a new root and re-apply to the new root.
    //    template<typename InnerLhsNode, typename InnerOpType, typename InnerR1, typename InnerROpType, typename InnerR2, typename InnerIndicesType>
    //    struct Apply<Node<InnerLhsNode, InnerOpType, Node<InnerR1, InnerROpType, InnerR2> >, InnerIndicesType,
    //            typename boost::enable_if
    //                <
    //                    TreeIsAssociative<Node<InnerLhsNode, InnerOpType, Node<InnerR1, InnerROpType, InnerR2> > >
    //                >::type>
    //    {
    //        typedef Node<InnerLhsNode, InnerOpType, Node<InnerR1, InnerROpType, InnerR2> > NodeType;
    //        typedef typename AssociativeTransform<NodeType>::TransformedNodeType NewTree;

    //        typedef typename RemoveUnecessaryTemporaries<NewTree, InnerIndicesType, StartIndex>::TransformedNodeType TransformedNodeType;
    //        typedef typename RemoveUnecessaryTemporaries<NewTree, InnerIndicesType, StartIndex>::TransformedIndicesType TransformedIndicesType;
    //    };

    //    // case 2.2 - If not associative, but commutative, and the swapped tree is associative, 
    //    // do a commutative, followed by an associative.
    //    template<typename InnerL1, typename InnerLOp, typename InnerL2, typename InnerOp, typename InnerR1, typename InnerROp, typename InnerR2, typename InnerIndicesType>
    //    struct Apply<Node<Node<InnerL1, InnerLOp, InnerL2>, InnerOp, Node<InnerR1, InnerROp, InnerR2> >, InnerIndicesType,
    //        typename boost::enable_if
    //        <
    //            boost::mpl::and_
    //            <
    //                boost::mpl::not_<TreeIsAssociative<Node<Node<InnerL1, InnerLOp, InnerL2>, InnerOp, Node<InnerR1, InnerROp, InnerR2> > > >,
    //                CommutativeTraits<typename Node<InnerL1, InnerLOp, InnerL2>::ResultType, InnerOp, typename Node<InnerR1, InnerROp, InnerR2>::ResultType>,
    //                TreeIsAssociative<Node<Node<InnerL1, InnerLOp, InnerL2>, InnerOp, Node<InnerR1, InnerROp, InnerR2> > >
    //            >
    //        >::type >
    //    {
    //        typedef Node<Node<InnerL1, InnerLOp, InnerL2>, InnerOp, Node<InnerR1, InnerROp, InnerR2> > InnerNodeType;
    //        typedef typename CommutativeTransform<InnerNodeType, Indices, StartIndex>::TransformedNodeType CommutativeTransformedNodeType;
    //        typedef typename CommutativeTransform<InnerNodeType, Indices, StartIndex>::TransformedIndicesType CommutativeIndicesType;

    //        typedef typename AssociativeTransform<CommutativeTransformedNodeType>::TransformedNodeType AssociativeTransformNodeType;

    //        typedef typename RemoveUnecessaryTemporaries<AssociativeTransformNodeType, CommutativeIndicesType, StartIndex>::TransformedNodeType TransformedNodeType;
    //        typedef typename RemoveUnecessaryTemporaries<AssociativeTransformNodeType, CommutativeIndicesType, StartIndex>::TransformedIndicesType TransformedIndicesType;
    //    };

    //    // Step 1 - Optimize both side independently.
    //    typedef Node<L1, LOp, L2> LhsNodeType;
    //    typedef Node<R1, ROp, R2> RhsNodeType;
    //    
    //    typedef typename RemoveUnecessaryTemporaries<LhsNodeType, Indices, StartIndex>::TransformedNodeType TransformedLhsNodeType;
    //    typedef typename RemoveUnecessaryTemporaries<LhsNodeType, Indices, StartIndex>::TransformedIndicesType TransformedIndicesFirstPass;

    //    static const int RightStart = StartIndex + LhsNodeType::TotalCount;
    //    typedef typename RemoveUnecessaryTemporaries<RhsNodeType, TransformedIndicesFirstPass, RightStart>::TransformedNodeType TransformedRhsNodeType;
    //    typedef typename RemoveUnecessaryTemporaries<RhsNodeType, TransformedIndicesFirstPass, RightStart>::TransformedIndicesType TransformedIndicesSecondPass;

    //    // Step 2 - Check if we can avoid a temporary through a full associative tree transform, and do 
    //    // it if we can.
    //    typedef Node<TransformedLhsNodeType, Op, TransformedRhsNodeType > NodeType;
    //    typedef typename Apply<NodeType, TransformedIndicesSecondPass>::TransformedNodeType TransformedNodeType;
    //    typedef typename Apply<NodeType, TransformedIndicesSecondPass>::TransformedIndicesType TransformedIndicesType;
    //};

}

#endif //NEKTAR_EXPRESSION_TEMPLATES_REMOVE_ALL_UNECESSARY_TEMPORARIES_H


































//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////
////
//// The MIT License
////
//// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//// Department of Aeronautics, Imperial College London (UK), and Scientific
//// Computing and Imaging Institute, University of Utah (USA).
////
//// License for the specific language governing rights and limitations under
//// Permission is hereby granted, free of charge, to any person obtaining a
//// copy of this software and associated documentation files (the "Software"),
//// to deal in the Software without restriction, including without limitation
//// the rights to use, copy, modify, merge, publish, distribute, sublicense,
//// and/or sell copies of the Software, and to permit persons to whom the
//// Software is furnished to do so, subject to the following conditions:
////
//// The above copyright notice and this permission notice shall be included
//// in all copies or substantial portions of the Software.
////
//// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//// DEALINGS IN THE SOFTWARE.
////
/////////////////////////////////////////////////////////////////////////////////
//
//#ifndef NEKTAR_EXPRESSION_TEMPLATES_REMOVE_ALL_UNECESSARY_TEMPORARIES_H
//#define NEKTAR_EXPRESSION_TEMPLATES_REMOVE_ALL_UNECESSARY_TEMPORARIES_H
//
//#include <ExpressionTemplates/AssociativeTraits.hpp>
//#include <ExpressionTemplates/CommutativeTraits.hpp>
//#include <ExpressionTemplates/CommutativeTransform.hpp>
//#include <ExpressionTemplates/AssociativeTransform.hpp>
//
//namespace Nektar
//{
//    // Input is a node in an AssociativeCommutative cluster, and we are determining if we should swap nodes.
//    template<typename NodeType, typename OpType, typename IndicesType, unsigned int StartIndex, typename enabled=void>
//    struct SortAssociativeCommutativeCluster
//    {
//        static const unsigned int Id = 0;
//
//        typedef NodeType TransformedNodeType;
//        typedef IndicesType TransformedIndicesType;
//        
//    };
//
//    template<typename L1, typename LOp, typename L2, typename OpType, typename R1, typename ROp, typename R2,
//             typename AssociativeCommutativeOpType, typename IndicesType, unsigned int StartIndex>
//    struct SortAssociativeCommutativeCluster< expt::Node< expt::Node<L1, LOp, L2>, OpType, expt::Node<R1, ROp, R2> >, AssociativeCommutativeOpType, IndicesType, StartIndex,
//        typename boost::enable_if
//         <
//             //boost::mpl::and_
//             //<
//                 boost::mpl::and_
//                 <
//                     boost::is_same<OpType, AssociativeCommutativeOpType>,
//                     boost::is_same<LOp, void>,
//                     boost::mpl::not_<boost::is_same<ROp, void> >
//                 >
//             //>
//         >::type>
//    {
//        static const unsigned int Id = 1;
//
//        typedef expt::Node<expt::Node<L1, LOp, L2>, OpType, expt::Node<R1, ROp, R2> > Tree0;
//        typedef typename expt::CommutativeTransform<Tree0, IndicesType, StartIndex>::TransformedNodeType TransformedNodeType;
//        typedef typename expt::CommutativeTransform<Tree0, IndicesType, StartIndex>::TransformedIndicesType TransformedIndicesType;
//    };
//
//    template<typename L1, typename LOp, typename L2, typename OpType, typename R1, typename ROp, typename R2,
//             typename AssociativeCommutativeOpType, typename IndicesType, unsigned int StartIndex>
//    struct SortAssociativeCommutativeCluster< expt::Node< expt::Node<L1, LOp, L2>, OpType, expt::Node<R1, ROp, R2> >, AssociativeCommutativeOpType, IndicesType, StartIndex,
//        typename boost::enable_if
//         <
//             //boost::mpl::and_
//             //<
//                 boost::mpl::and_
//                 <
//                     boost::is_same<OpType, AssociativeCommutativeOpType>,
//                     boost::is_same<LOp, AssociativeCommutativeOpType>,
//                     boost::mpl::not_<boost::is_same<ROp, void> >,
//                     boost::mpl::less<boost::mpl::int_<expt::TemporaryCount<L2>::Value>, boost::mpl::int_<expt::TemporaryCount<expt::Node<R1, ROp, R2> >::Value+1> >
//                 >
//             //>
//         >::type>
//    {
//        static const unsigned int Id = 2;
//
//        // Example tree: A+B+CD
//        typedef expt::Node<expt::Node<L1, LOp, L2>, OpType, expt::Node<R1, ROp, R2> > Tree0;
//
//        // Inverse associative so we have A+(B+CD)
//        typedef typename InverseAssociativeTransform<Tree0>::TransformedNodeType Tree1;
//        typedef typename Tree1::Right RightTree1;
//        typedef typename Tree1::Left LeftTree1;
//        static const int RightStart1 = StartIndex + LeftTree1::TotalCount;
//
//        // Commutative transform on the right so we have A+(CD+B)
//        typedef typename expt::CommutativeTransform<RightTree1, IndicesType, RightStart1>::TransformedNodeType RightTree2;
//        typedef typename expt::CommutativeTransform<RightTree1, IndicesType, RightStart1>::TransformedIndicesType Indices2;
//        typedef typename expt::Node<LeftTree1, typename Tree1::OpType, RightTree2> Tree2;
//
//
//        // Associative transform to get A + CD + B
//        typedef typename AssociativeTransform<Tree2>::TransformedNodeType Tree3;
//        typedef Indices2 Indices3;
//
//        // Recurse down the left side.
//        typedef typename SortAssociativeCommutativeCluster<typename Tree3::Left, AssociativeCommutativeOpType, Indices3, StartIndex>::TransformedNodeType LeftTree4;
//        typedef typename SortAssociativeCommutativeCluster<typename Tree3::Left, AssociativeCommutativeOpType, Indices3, StartIndex>::TransformedIndicesType Indices4;
//
//        typedef expt::Node<LeftTree4, OpType, typename Tree3::Right> TransformedNodeType;
//        typedef Indices4 TransformedIndicesType;
//    };
//
//    // Updates the parse tree so there are not unecessary temporaries.
//    template<typename NodeType, typename IndicesType, unsigned int StartIndex = 0, typename enabled = void>
//    struct RemoveUnecessaryTemporaries
//    {
//        static const unsigned int SpecializationId = 0;
//
//        typedef typename NodeType::Left Left0Type;
//        typedef typename NodeType::Right Right0Type;
//        typedef typename NodeType::OpType OpType;
//        typedef typename NodeType Tree0Type;
//        typedef IndicesType Indices0Type;
//
//
//        static const unsigned int RightStartIndex = StartIndex + Left0Type::TotalCount;
//        typedef typename RemoveUnecessaryTemporaries<Right0Type, Indices0Type, RightStartIndex>::TransformedNodeType Right1Type;
//        typedef typename RemoveUnecessaryTemporaries<Right0Type, Indices0Type, RightStartIndex>::TransformedIndicesType Indices1Type;
//        typedef Left0Type Left1Type;
//        typedef Node<Left1Type, OpType, Right1Type> Tree1Type;
//
//        typedef typename RemoveUnecessaryTemporaries<Left1Type, Indices1Type, StartIndex>::TransformedNodeType Left2Type;
//        typedef typename RemoveUnecessaryTemporaries<Left1Type, Indices1Type, StartIndex>::TransformedIndicesType Indices2Type;
//        typedef Right1Type Right2Type;
//        typedef Node<Left2Type, OpType, Right2Type> Tree2Type;
//
//        typedef Tree2Type TransformedNodeType;
//        typedef Indices2Type TransformedIndicesType;
//    };
//
//    // Constant node (Base Case 0) - don't need to do anything.
//    template<typename T, typename IndicesType, unsigned int StartIndex>
//    struct RemoveUnecessaryTemporaries<expt::Node<T, void, void>, IndicesType, StartIndex>
//    {
//        static const unsigned int SpecializationId = 1;
//        typedef expt::Node<T, void, void> TransformedNodeType;
//        typedef IndicesType TransformedIndicesType;
//    };
//
//    // Unary node (Base Case 1) - don't need to do anything.
//    template<typename T, typename UnaryOp, typename IndicesType, unsigned int StartIndex>
//    struct RemoveUnecessaryTemporaries<expt::Node<T, UnaryOp, void>, IndicesType, StartIndex>
//    {
//        static const unsigned int SpecializationId = 2;
//        typedef expt::Node<T, UnaryOp, void> TransformedNodeType;
//        typedef IndicesType TransformedIndicesType;
//    };
//
//    // Commutative-only node.
//    template<typename L1, typename LOp, typename L2, typename OpType, typename R1, typename ROp, typename R2,
//             typename IndicesType, unsigned int StartIndex>
//    struct RemoveUnecessaryTemporaries<Node<Node<L1, LOp, L2>, OpType, Node<R1, ROp, R2> >, IndicesType, StartIndex,
//             typename boost::enable_if
//             <
//                boost::mpl::and_
//                <
//                    CommutativeTraits<typename Node<L1, LOp, L2>::ResultType, OpType, typename Node<R1, ROp, R2>::ResultType>,
//                    boost::mpl::not_<AssociativeTraits<typename Node<L1, LOp, L2>::ResultType, OpType, typename R1::ResultType, ROp> >
//                >
//             >::type>
//    {
//        static const unsigned int SpecializationId = 3;
//
//        template<typename NodeType, typename IndicesType, unsigned int StartIndex, typename enabled = void>
//        struct PerformCommutativeTransformIfNeeded
//        {
//            typedef NodeType TransformedNodeType;
//            typedef IndicesType TransformedIndicesType;
//        };
//
//        template<typename LhsType, typename OpType, typename RhsType, typename IndicesType, unsigned int StartIndex>
//        struct PerformCommutativeTransformIfNeeded< expt::Node<LhsType, OpType, RhsType>, IndicesType, StartIndex,
//            typename boost::enable_if
//            <
//                boost::mpl::and_
//                <
//                    boost::is_same<typename LhsType::OpType, void>,
//                    boost::mpl::not_<boost::is_same<typename RhsType::OpType, void> >
//                >
//            >::type>
//        {
//            typedef expt::Node<LhsType, OpType, RhsType> NodeType;
//            typedef typename expt::CommutativeTransform<NodeType, IndicesType, StartIndex>::TransformedNodeType TransformedNodeType;
//            typedef typename expt::CommutativeTransform<NodeType, IndicesType, StartIndex>::TransformedIndicesType TransformedIndicesType;
//        };
//
//        typedef typename Node<L1, LOp, L2> Left0Type;
//        typedef typename Node<R1, ROp, R2> Right0Type;
//        typedef typename Node<Left0Type, OpType, Right0Type> Tree0Type;
//        typedef IndicesType Indices0Type;
//
//        typedef typename RemoveUnecessaryTemporaries<Left0Type, Indices0Type, StartIndex>::TransformedNodeType Left1Type;
//        typedef typename RemoveUnecessaryTemporaries<Left0Type, Indices0Type, StartIndex>::TransformedIndicesType Indices1Type;
//        typedef Right0Type Right1Type;
//
//        static const unsigned int RightStartIndex = StartIndex + Left1Type::TotalCount;
//        typedef typename RemoveUnecessaryTemporaries<Right1Type, Indices1Type, RightStartIndex>::TransformedNodeType Right2Type;
//        typedef typename RemoveUnecessaryTemporaries<Right1Type, Indices1Type, RightStartIndex>::TransformedIndicesType Indices2Type;
//        typedef Left1Type Left2Type;
//        typedef Node<Left2Type, OpType, Right2Type> Tree2Type;
//
//        typedef typename PerformCommutativeTransformIfNeeded<Tree2Type, Indices2Type, StartIndex>::TransformedNodeType Tree3Type;
//        typedef typename PerformCommutativeTransformIfNeeded<Tree2Type, Indices2Type, StartIndex>::TransformedIndicesType Indices3Type;
//
//        typedef Tree3Type TransformedNodeType;
//        typedef Indices3Type TransformedIndicesType;
//    };
//
//    // Associative Only node.
//    template<typename L1, typename LOp, typename L2, typename OpType, typename R1, typename ROp, typename R2,
//             typename IndicesType, unsigned int StartIndex>
//    struct RemoveUnecessaryTemporaries<Node<Node<L1, LOp, L2>, OpType, Node<R1, ROp, R2> >, IndicesType, StartIndex,
//             typename boost::enable_if
//             <
//                boost::mpl::and_
//                <
//                    boost::mpl::not_<CommutativeTraits<typename Node<L1, LOp, L2>::ResultType, OpType, typename Node<R1, ROp, R2>::ResultType> >,
//                    AssociativeTraits<typename Node<L1, LOp, L2>::ResultType, OpType, typename R1::ResultType, ROp>
//                >
//             >::type>
//    {
//        static const unsigned int SpecializationId = 4;
//
//        typedef typename Node<L1, LOp, L2> Left0Type;
//        typedef typename Node<R1, ROp, R2> Right0Type;
//        typedef typename Node<Left0Type, OpType, Right0Type> Tree0Type;
//        typedef IndicesType Indices0Type;
//        static const unsigned int Right0StartIndex = StartIndex + Left0Type::TotalCount;
//
//        typedef typename RemoveUnecessaryTemporaries<Right0Type, Indices0Type, Right0StartIndex>::TransformedNodeType Right1Type;
//        typedef typename RemoveUnecessaryTemporaries<Right0Type, Indices0Type, Right0StartIndex>::TransformedIndicesType Indices1Type;
//        typedef Left0Type Left1Type;
//        typedef Node<Left1Type, OpType, Right1Type> Tree1Type;
//
//        typedef typename AssociativeTransform<Tree1Type>::TransformedNodeType Tree2Type;
//        typedef typename Tree2Type::Left Left2Type;
//        typedef typename Tree2Type::Right Right2Type;
//        typedef Indices1Type Indices2Type;
//
//        typedef typename RemoveUnecessaryTemporaries<Left2Type, Indices2Type, StartIndex>::TransformedNodeType Left3Type;
//        typedef typename RemoveUnecessaryTemporaries<Left2Type, Indices2Type, StartIndex>::TransformedIndicesType Indices3Type;
//        typedef typename Right2Type Right3Type;
//        typedef typename Node<Left3Type, OpType, Right3Type> Tree3Type;
//
//        typedef Tree3Type TransformedNodeType;
//        typedef Indices3Type TransformedIndicesType;
//    };
//
//    // Associative and commutative Only node.
//    template<typename L1, typename LOp, typename L2, typename OpType, typename R1, typename ROp, typename R2,
//             typename IndicesType, unsigned int StartIndex>
//    struct RemoveUnecessaryTemporaries<Node<Node<L1, LOp, L2>, OpType, Node<R1, ROp, R2> >, IndicesType, StartIndex,
//             typename boost::enable_if
//             <
//                boost::mpl::and_
//                <
//                    CommutativeTraits<typename Node<L1, LOp, L2>::ResultType, OpType, typename Node<R1, ROp, R2>::ResultType>,
//                    AssociativeTraits<typename Node<L1, LOp, L2>::ResultType, OpType, typename R1::ResultType, ROp>
//                >
//             >::type>
//    {
//        static const unsigned int SpecializationId = 5;
//
//        typedef typename Node<L1, LOp, L2> Left0Type;
//        typedef typename Node<R1, ROp, R2> Right0Type;
//        typedef typename Node<Left0Type, OpType, Right0Type> Tree0Type;
//        typedef IndicesType Indices0Type;
//        static const unsigned int Right0StartIndex = StartIndex + Left0Type::TotalCount;
//
//        typedef typename RemoveUnecessaryTemporaries<Right0Type, Indices0Type, Right0StartIndex>::TransformedNodeType Right1Type;
//        typedef typename RemoveUnecessaryTemporaries<Right0Type, Indices0Type, Right0StartIndex>::TransformedIndicesType Indices1Type;
//        typedef Left0Type Left1Type;
//        typedef Node<Left1Type, OpType, Right1Type> Tree1Type;
//
//        typedef typename AssociativeTransform<Tree1Type>::TransformedNodeType Tree2Type;
//        typedef typename Tree2Type::Left Left2Type;
//        typedef typename Tree2Type::Right Right2Type;
//        typedef Indices1Type Indices2Type;
//
//        typedef typename RemoveUnecessaryTemporaries<Left2Type, Indices2Type, StartIndex>::TransformedNodeType Left3Type;
//        typedef typename RemoveUnecessaryTemporaries<Left2Type, Indices2Type, StartIndex>::TransformedIndicesType Indices3Type;
//        typedef typename Right2Type Right3Type;
//        typedef typename Node<Left3Type, OpType, Right3Type> Tree3Type;
//
//        typedef Tree3Type TransformedNodeType;
//        typedef Indices3Type TransformedIndicesType;
//    };
//
//
//
//    //template<typename L, typename OpType, typename R, typename IndicesType, unsigned int StartIndex>
//    //struct RemoveUnecessaryTemporaries<expt::Node<L, OpType, R>, IndicesType, StartIndex>
//    //{
//
//    //    template<typename NodeType, typename IndicesType, unsigned int StartIndex, typename enabled = void>
//    //    struct PerformCommutativeTransformIfNeeded
//    //    {
//    //        typedef NodeType TransformedNodeType;
//    //        typedef IndicesType TransformedIndicesType;
//    //    };
//
//    //    template<typename LhsType, typename OpType, typename RhsType, typename IndicesType, unsigned int StartIndex>
//    //    struct PerformCommutativeTransformIfNeeded< expt::Node<LhsType, OpType, RhsType>, IndicesType, StartIndex,
//    //        typename boost::enable_if
//    //        <
//    //            boost::mpl::and_
//    //            <
//    //                expt::CommutativeTraits<typename LhsType::ResultType, OpType, typename RhsType::ResultType>,
//    //                boost::mpl::less<boost::mpl::int_<expt::TemporaryCount<expt::Node<RhsType, OpType, LhsType> >::Value>, boost::mpl::int_<expt::TemporaryCount<expt::Node<LhsType, OpType, RhsType> >::Value> >
//    //            >
//    //        >::type>
//    //    {
//    //        typedef expt::Node<LhsType, OpType, RhsType> NodeType;
//    //        typedef typename expt::CommutativeTransform<NodeType, IndicesType, StartIndex>::TransformedNodeType TransformedNodeType;
//    //        typedef typename expt::CommutativeTransform<NodeType, IndicesType, StartIndex>::TransformedIndicesType TransformedIndicesType;
//    //    };
//
//    //    static const int RightStart = StartIndex + L::TotalCount;
//    //    
//    //    // Optimize the right child.
//    //    typedef typename RemoveUnecessaryTemporaries<R, IndicesType, RightStart>::TransformedNodeType RightNode0Type;
//    //    typedef typename RemoveUnecessaryTemporaries<R, IndicesType, RightStart>::TransformedIndicesType Indices0Type;
//    //    typedef typename expt::Node<L, OpType, RightNode0Type> Tree0Type;
//
//    //    // If Op and ROp are associative, apply the associative transformation.
//    //    // Note that the specializations for the AssociativeTransform object only apply when the right 
//    //    // child is a binary node with the appropriate operator, so this transformation will not mistakenly
//    //    // turn a constant right child into a binary right child.
//    //    typedef typename AssociativeTransform<Tree0Type>::TransformedNodeType Tree1Type;
//    //    typedef typename Tree1Type::Left LeftNode1Type;
//    //    typedef Indices0Type Indices1Type;
//
//    //    // Optimize the left child.
//    //    typedef typename RemoveUnecessaryTemporaries<LeftNode1Type, Indices1Type, StartIndex>::TransformedNodeType LeftNode2Type;
//    //    typedef typename RemoveUnecessaryTemporaries<LeftNode1Type, Indices1Type, StartIndex>::TransformedIndicesType Indices2Type;
//    //    typedef expt::Node<LeftNode2Type, OpType, typename Tree1Type::Right> Tree2Type;
//
//    //    // Apply commutative if Op is commutative and the number of temps for the rhs is greater 
//    //    // than the number for the left.
//    //    //typedef typename PerformCommutativeTransformIfNeeded<Tree2Type, Indices2Type, StartIndex>::TransformedNodeType TransformedNodeType;
//    //    //typedef typename PerformCommutativeTransformIfNeeded<Tree2Type, Indices2Type, StartIndex>::TransformedIndicesType TransformedIndicesType;
//
//    //    // Sort the tree so that the leftmost node of an associative cluster requires the most temporaries of any 
//    //    // descendant of the cluster.
//    //    typedef typename SortAssociativeCommutativeCluster<Tree2Type, OpType, Indices2Type, StartIndex>::TransformedNodeType Tree3Type;
//    //    typedef typename SortAssociativeCommutativeCluster<Tree2Type, OpType, Indices2Type, StartIndex>::TransformedIndicesType Indices3Type;
//
//    //    typedef Tree3Type TransformedNodeType;
//    //    typedef Indices3Type TransformedIndicesType;
//    //    
//
//    //};
//
//}
//
//namespace Nektar
//{
//    //template<typename NodeType, typename IndicesType, unsigned int StartIndex = 0>
//    //class RemoveUnecessaryTemporaries;
//
//    //template<typename NodeType, typename IndicesType, unsigned int StartIndex>
//    //class Default
//    //{
//    //    private:
//    //        template<typename InnerNodeType, typename InnerIndicesType, unsigned int InnerStartIndex>
//    //        struct Apply;
//
//    //        template<typename T, typename InnerIndicesType, unsigned int InnerStartIndex>
//    //        struct Apply<Node<T, void, void>, InnerIndicesType, InnerStartIndex>
//    //        {
//    //            typedef Node<T, void, void> TransformedNodeType;
//    //            typedef InnerIndicesType TransformedIndicesType;
//    //        };
//
//    //        template<typename T, typename Op, typename InnerIndicesType, unsigned int InnerStartIndex>
//    //        struct Apply<Node<T, Op, void>, InnerIndicesType, InnerStartIndex>
//    //        {
//    //            // TODO - still need to finish the unary cases.
//    //        };
//
//    //        template<typename L, typename Op, typename R, typename InnerIndicesType, unsigned int InnerStartIndex>
//    //        struct Apply<Node<L, Op, R>, InnerIndicesType, InnerStartIndex>
//    //        {
//    //            typedef typename RemoveUnecessaryTemporaries<L, InnerIndicesType, InnerStartIndex>::TransformedNodeType OptimizedLeft;
//    //            typedef typename RemoveUnecessaryTemporaries<L, InnerIndicesType, InnerStartIndex>::TransformedIndicesType FirstOptimizedIndices;
//
//    //            static const int RightStart = InnerStartIndex + L::TotalCount;
//    //            typedef typename RemoveUnecessaryTemporaries<R, FirstOptimizedIndices, RightStart>::TransformedNodeType OptimizedRight;
//    //            typedef typename RemoveUnecessaryTemporaries<R, FirstOptimizedIndices, RightStart>::TransformedIndicesType SecondOptimizedIndices;
//
//    //            typedef Node<OptimizedLeft, Op, OptimizedRight> TransformedNodeType;
//    //            typedef SecondOptimizedIndices TransformedIndicesType;
//    //        };
//
//    //    public:
//
//    //        typedef typename Apply<NodeType, IndicesType, StartIndex>::TransformedNodeType TransformedNodeType;
//    //        typedef typename Apply<NodeType, IndicesType, StartIndex>::TransformedIndicesType TransformedIndicesType;
//    //};
//
//    //// RemoveUnecessaryTemporaries is a parse tree, compile time optimizer.
//    //template<typename NodeType, typename IndicesType, unsigned int StartIndex>
//    //class RemoveUnecessaryTemporaries : public Default<NodeType, IndicesType, StartIndex>
//    //{
//    //};
//
//    //// Constant-Constant specialization.  
//    //template<typename LeftType, typename Op, typename RightType, typename IndicesType, unsigned int StartIndex>
//    //class RemoveUnecessaryTemporaries<Node<Node<LeftType, void, void>, Op, Node<RightType, void, void> >, IndicesType, StartIndex > :
//    //    public Default<Node<Node<LeftType, void, void>, Op, Node<RightType, void, void> >, IndicesType, StartIndex > {};
//
//    //// Left anything, right constant.  Nothing do do.
//    //template<typename L1, typename LOp, typename L2, typename Op, typename R, typename Indices, unsigned int StartIndex>
//    //class RemoveUnecessaryTemporaries<Node<Node<L1, LOp, L2>, Op, Node<R, void, void> >, Indices, StartIndex > :
//    //    public Default<Node<Node<L1, LOp, L2>, Op, Node<R, void, void> >, Indices, StartIndex > {};
//
//    //// Case 1 - Left constant - Right anything.
//    //template<typename LeftType, typename Op, typename RightArg1, typename RightOp, typename RightArg2, typename IndicesType, unsigned int StartIndex>
//    //class RemoveUnecessaryTemporaries<Node<Node<LeftType, void, void>, Op, Node<RightArg1, RightOp, RightArg2> >, IndicesType, StartIndex >
//    //{
//    //    private:
//    //        template<typename NodeType, typename enabled = void>
//    //        struct Apply : public Default<NodeType, IndicesType, StartIndex> {};
//
//    //        template<typename LeftNode, typename InnerOp, typename RightNode>
//    //        struct Apply<Node<LeftNode, InnerOp, RightNode>,
//    //                typename boost::enable_if
//    //                    <
//    //                        CommutativeTraits<typename LeftNode::ResultType, InnerOp, typename RightNode::ResultType>
//    //                    >::type>
//    //        {
//    //            typedef Node<LeftNode, InnerOp, RightNode> InnerNodeType;
//    //            typedef CommutativeTransform<InnerNodeType, IndicesType, StartIndex> Transform;
//
//    //            typedef typename Transform::TransformedNodeType CommutativeNodeType;
//    //            typedef typename Transform::TransformedIndicesType CommutativeIndices;
//
//    //            typedef typename RemoveUnecessaryTemporaries<CommutativeNodeType, CommutativeIndices, StartIndex>::TransformedNodeType TransformedNodeType;
//    //            typedef typename RemoveUnecessaryTemporaries<CommutativeNodeType, CommutativeIndices, StartIndex>::TransformedIndicesType TransformedIndicesType;
//    //        };
//    //    
//    //    public:
//    //        typedef Node<Node<LeftType, void, void>, Op, Node<RightArg1, RightOp, RightArg2> > NodeType;
//    //        typedef typename Apply<NodeType>::TransformedNodeType TransformedNodeType;
//    //        typedef typename Apply<NodeType>::TransformedIndicesType TransformedIndicesType;
//    //};
//
//    //// Case 2 - Left and Right nodes are arbitrary.
//    //template<typename L1, typename LOp, typename L2, typename Op, typename R1, typename ROp, typename R2, typename Indices, unsigned int StartIndex>
//    //struct RemoveUnecessaryTemporaries<Node<Node<L1, LOp, L2>, Op, Node<R1, ROp, R2> >, Indices, StartIndex >
//    //{
//    //    template<typename NodeType, typename InnerIndicesType, typename enabled = void>
//    //    struct Apply : public Default<NodeType, InnerIndicesType, StartIndex> {};
//
//    //    // case 2.1 - If the node is associative, then apply the associative transformation to get a new 
//    //    // tree with a new root and re-apply to the new root.
//    //    template<typename InnerLhsNode, typename InnerOpType, typename InnerR1, typename InnerROpType, typename InnerR2, typename InnerIndicesType>
//    //    struct Apply<Node<InnerLhsNode, InnerOpType, Node<InnerR1, InnerROpType, InnerR2> >, InnerIndicesType,
//    //            typename boost::enable_if
//    //                <
//    //                    TreeIsAssociative<Node<InnerLhsNode, InnerOpType, Node<InnerR1, InnerROpType, InnerR2> > >
//    //                >::type>
//    //    {
//    //        typedef Node<InnerLhsNode, InnerOpType, Node<InnerR1, InnerROpType, InnerR2> > NodeType;
//    //        typedef typename AssociativeTransform<NodeType>::TransformedNodeType NewTree;
//
//    //        typedef typename RemoveUnecessaryTemporaries<NewTree, InnerIndicesType, StartIndex>::TransformedNodeType TransformedNodeType;
//    //        typedef typename RemoveUnecessaryTemporaries<NewTree, InnerIndicesType, StartIndex>::TransformedIndicesType TransformedIndicesType;
//    //    };
//
//    //    // case 2.2 - If not associative, but commutative, and the swapped tree is associative, 
//    //    // do a commutative, followed by an associative.
//    //    template<typename InnerL1, typename InnerLOp, typename InnerL2, typename InnerOp, typename InnerR1, typename InnerROp, typename InnerR2, typename InnerIndicesType>
//    //    struct Apply<Node<Node<InnerL1, InnerLOp, InnerL2>, InnerOp, Node<InnerR1, InnerROp, InnerR2> >, InnerIndicesType,
//    //        typename boost::enable_if
//    //        <
//    //            boost::mpl::and_
//    //            <
//    //                boost::mpl::not_<TreeIsAssociative<Node<Node<InnerL1, InnerLOp, InnerL2>, InnerOp, Node<InnerR1, InnerROp, InnerR2> > > >,
//    //                CommutativeTraits<typename Node<InnerL1, InnerLOp, InnerL2>::ResultType, InnerOp, typename Node<InnerR1, InnerROp, InnerR2>::ResultType>,
//    //                TreeIsAssociative<Node<Node<InnerL1, InnerLOp, InnerL2>, InnerOp, Node<InnerR1, InnerROp, InnerR2> > >
//    //            >
//    //        >::type >
//    //    {
//    //        typedef Node<Node<InnerL1, InnerLOp, InnerL2>, InnerOp, Node<InnerR1, InnerROp, InnerR2> > InnerNodeType;
//    //        typedef typename CommutativeTransform<InnerNodeType, Indices, StartIndex>::TransformedNodeType CommutativeTransformedNodeType;
//    //        typedef typename CommutativeTransform<InnerNodeType, Indices, StartIndex>::TransformedIndicesType CommutativeIndicesType;
//
//    //        typedef typename AssociativeTransform<CommutativeTransformedNodeType>::TransformedNodeType AssociativeTransformNodeType;
//
//    //        typedef typename RemoveUnecessaryTemporaries<AssociativeTransformNodeType, CommutativeIndicesType, StartIndex>::TransformedNodeType TransformedNodeType;
//    //        typedef typename RemoveUnecessaryTemporaries<AssociativeTransformNodeType, CommutativeIndicesType, StartIndex>::TransformedIndicesType TransformedIndicesType;
//    //    };
//
//    //    // Step 1 - Optimize both side independently.
//    //    typedef Node<L1, LOp, L2> LhsNodeType;
//    //    typedef Node<R1, ROp, R2> RhsNodeType;
//    //    
//    //    typedef typename RemoveUnecessaryTemporaries<LhsNodeType, Indices, StartIndex>::TransformedNodeType TransformedLhsNodeType;
//    //    typedef typename RemoveUnecessaryTemporaries<LhsNodeType, Indices, StartIndex>::TransformedIndicesType TransformedIndicesFirstPass;
//
//    //    static const int RightStart = StartIndex + LhsNodeType::TotalCount;
//    //    typedef typename RemoveUnecessaryTemporaries<RhsNodeType, TransformedIndicesFirstPass, RightStart>::TransformedNodeType TransformedRhsNodeType;
//    //    typedef typename RemoveUnecessaryTemporaries<RhsNodeType, TransformedIndicesFirstPass, RightStart>::TransformedIndicesType TransformedIndicesSecondPass;
//
//    //    // Step 2 - Check if we can avoid a temporary through a full associative tree transform, and do 
//    //    // it if we can.
//    //    typedef Node<TransformedLhsNodeType, Op, TransformedRhsNodeType > NodeType;
//    //    typedef typename Apply<NodeType, TransformedIndicesSecondPass>::TransformedNodeType TransformedNodeType;
//    //    typedef typename Apply<NodeType, TransformedIndicesSecondPass>::TransformedIndicesType TransformedIndicesType;
//    //};
//
//}
//
//#endif //NEKTAR_EXPRESSION_TEMPLATES_REMOVE_ALL_UNECESSARY_TEMPORARIES_H
