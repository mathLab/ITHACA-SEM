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

#ifndef EXPRESSION_TEMPLATES_PUSH_DOWN_UNARY_NODES_HPP
#define EXPRESSION_TEMPLATES_PUSH_DOWN_UNARY_NODES_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <ExpressionTemplates/InvertNode.hpp>

namespace expt
{
    // Pushes unary operations as far down the tree as possible.  This 
    // makes helps maximize the size of the clusters.
    //
    // For example A - (B-C) can be transformed to 
    // A + (-B + C) through the use of this class, creating a two node 
    // cluster of the + operator.

    template<typename NodeType>
    struct PushDownUnaryNodes;

    template<typename DataType>
    struct PushDownUnaryNodes<Node<DataType, void, void> >
    {
        typedef Node<DataType, void, void> Type;
    };

    namespace impl
    {
        template<typename ChildType, typename OpType, typename enabled = void>
        struct PushDownUnaryNodesUnaryImpl;

        template<typename OpType, typename OriginalType, typename TransformedType, typename RecursiveType, typename enabled=void>
        struct UnaryHelper;

        // If the types are the same, then the inversion added a unary node.
        // So recurse the child, then add the node.
        template<typename OpType, typename OriginalType, typename TransformedType, typename RecursiveType>
        struct UnaryHelper<OpType, OriginalType, TransformedType, RecursiveType, 
            typename boost::enable_if<boost::is_same<OriginalType, TransformedType> >::type>
        {
            typedef typename PushDownUnaryNodes<RecursiveType>::Type T0;
            typedef Node<T0, OpType> Type;
        };

        // If the types are different, just recurse as normal.
        template<typename OpType, typename OriginalType, typename TransformedType, typename RecursiveType>
        struct UnaryHelper<OpType, OriginalType, TransformedType, RecursiveType, 
            typename boost::enable_if<boost::mpl::not_<boost::is_same<OriginalType, TransformedType> > >::type>
        {
            typedef typename PushDownUnaryNodes<RecursiveType>::Type Type;
        };

        template<typename ChildType, typename OpType>
        struct PushDownUnaryNodesUnaryImpl<ChildType, OpType,
            typename boost::enable_if
            <
                boost::mpl::and_
                <
                    UnaryInverseOperatorTypeTraits<OpType>,
                    boost::mpl::not_<IsConstantNode<ChildType> >
                >
            >::type>
        {
            typedef typename InvertNode<OpType, ChildType>::Type TransformedChildType;
            typedef typename InvertNode<OpType, ChildType>::RecursionType RecursionType;

            typedef typename UnaryHelper<OpType, Node<ChildType, OpType>, TransformedChildType, RecursionType>::Type Type;
        };

        template<typename ChildType, typename OpType>
        struct PushDownUnaryNodesUnaryImpl<ChildType, OpType,
            typename boost::enable_if
            <
                boost::mpl::and_
                <
                    UnaryInverseOperatorTypeTraits<OpType>,
                    IsConstantNode<ChildType>
                >
            >::type>
        {
            typedef typename InvertNode<OpType, ChildType>::Type Type;
        };

        template<typename ChildType, typename OpType>
        struct PushDownUnaryNodesUnaryImpl<ChildType, OpType,
            typename boost::enable_if
            <
                boost::mpl::not_<UnaryInverseOperatorTypeTraits<OpType> >
            >::type>
        {
            typedef typename PushDownUnaryNodes<ChildType>::Type TransformedChildType;
            typedef Node<TransformedChildType, OpType> Type;
        };
    }

    template<typename ChildType, typename OpType>
    struct PushDownUnaryNodes<Node<ChildType, OpType, void> >
    {
        typedef typename impl::PushDownUnaryNodesUnaryImpl<ChildType, OpType>::Type Type;
    };

    template<typename LhsType, typename OpType, typename RhsType>
    struct PushDownUnaryNodes<Node<LhsType, OpType, RhsType> >
    {
        typedef typename PushDownUnaryNodes<LhsType>::Type TransformedLhsType;
        typedef typename PushDownUnaryNodes<RhsType>::Type TransformedRhsType;
        typedef Node<TransformedLhsType, OpType, TransformedRhsType> Type;
    };


}

#endif
#endif
