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

#ifndef EXPRESSION_TEMPLATES_INVERT_NODE_HPP
#define EXPRESSION_TEMPLATES_INVERT_NODE_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <ExpressionTemplates/Node.hpp>
#include <ExpressionTemplates/Operators.hpp>
#include <ExpressionTemplates/InverseOperatorTypeTraits.hpp>

namespace expt
{
    // Inverts a node with the given unary type.
    // Each specialization defines a RecursionType, which is used for algorithms 
    // that recursively traverse a tree.  Since some specialization add a unary 
    // node, the RecursionType is used to prevent infinite loops.
    template<typename InverseType, typename NodeType, typename enabled = void>
    struct InvertNode;

    // A -> -A
    template<typename InverseType, typename DataType>
    struct InvertNode<InverseType, Node<DataType, void, void> >
    {
        typedef Node<DataType, void, void> LeafType;
        typedef Node<LeafType, InverseType, void> Type;
        typedef LeafType RecursionType;
    };

    ////////////////////////////////////////////////
    // Unary Node Inversion
    ////////////////////////////////////////////////
    namespace impl
    {
        template<typename InverseType, typename ChildNodeType, typename UnaryOpType, typename enabled = void>
        struct InvertUnaryNode;

        // -A -> A
        // If the unary node corresponds to the InverseTypes unary node, then
        // they cancel out.
        template<typename InverseType, typename ChildNodeType, typename UnaryOpType>
        struct InvertUnaryNode<InverseType, ChildNodeType, UnaryOpType, 
            typename boost::enable_if
            <
                boost::is_same<InverseType, UnaryOpType>
            >::type>
        {
            typedef ChildNodeType Type;
            typedef Type RecursionType;
        };

        // 1/A -> - 1/A
        // If the unary operation does not relate to the InverseType, then this 
        // creates a new unary node.
        template<typename InverseType, typename ChildNodeType, typename UnaryOpType>
        struct InvertUnaryNode<InverseType, ChildNodeType, UnaryOpType, 
            typename boost::enable_if
            <
                boost::mpl::not_<boost::is_same<InverseType, UnaryOpType> >
            >::type>
        {
            typedef Node<ChildNodeType, UnaryOpType, void> NewChildType;
            typedef Node<NewChildType, InverseType, void> Type;
            typedef NewChildType RecursionType;
        };
    }

    template<typename InverseType, typename ChildNodeType, typename UnaryOpType>
    struct InvertNode<InverseType, Node<ChildNodeType, UnaryOpType, void> >
    {
        typedef typename impl::InvertUnaryNode<InverseType, ChildNodeType, UnaryOpType>::Type Type;
        typedef typename impl::InvertUnaryNode<InverseType, ChildNodeType, UnaryOpType>::RecursionType RecursionType;
    };

    ////////////////////////////////////////////////
    // Binary Node Inversion
    ////////////////////////////////////////////////
    namespace impl
    {
        template<typename InverseType, typename LhsNodeType, typename OpType, typename RhsNodeType, typename enabled = void>
        struct InvertBinaryNode;

        // A+B -> (-A)-B
        template<typename InverseType, typename LhsNodeType, typename OpType, typename RhsNodeType>
        struct InvertBinaryNode<InverseType, LhsNodeType, OpType, RhsNodeType,
            typename boost::enable_if
            <
                boost::mpl::not_<boost::is_same<typename UnaryInverseOperatorTypeTraits<InverseType>::template Invert<OpType>::Type, void> >
            >::type>
        {
            typedef typename UnaryInverseOperatorTypeTraits<InverseType>::template Invert<OpType>::Type TransformedOpType;
            typedef typename InvertNode<InverseType, LhsNodeType>::Type TransformedLhsType;
            typedef Node<TransformedLhsType, TransformedOpType, RhsNodeType> Type;
            typedef Type RecursionType;
        };

        // A*B -> -(A*B)
        template<typename InverseType, typename LhsNodeType, typename OpType, typename RhsNodeType>
        struct InvertBinaryNode<InverseType, LhsNodeType, OpType, RhsNodeType,
            typename boost::enable_if
            <
                boost::is_same<typename UnaryInverseOperatorTypeTraits<InverseType>::template Invert<OpType>::Type, void>
            >::type>
        {
            typedef Node<LhsNodeType, OpType, RhsNodeType> NewChildNodeType;
            typedef Node<NewChildNodeType, InverseType> Type;
            typedef NewChildNodeType RecursionType;
        };
    }

    template<typename InverseType, typename LhsNodeType, typename OpType, typename RhsNodeType>
    struct InvertNode<InverseType, Node<LhsNodeType, OpType, RhsNodeType> >
    {
        typedef typename impl::InvertBinaryNode<InverseType, LhsNodeType, OpType, RhsNodeType>::Type Type;
        typedef typename impl::InvertBinaryNode<InverseType, LhsNodeType, OpType, RhsNodeType>::RecursionType RecursionType;
    };
}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //EXPRESSION_TEMPLATES_INVERT_NODE_HPP


