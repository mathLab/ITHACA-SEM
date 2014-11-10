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

#ifndef EXPRESSION_TEMPLATES_FORWARD_ASSOCIATIVE_TRANSFORM_HPP
#define EXPRESSION_TEMPLATES_FORWARD_ASSOCIATIVE_TRANSFORM_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <ExpressionTemplates/Node.hpp>
#include <ExpressionTemplates/Operators.hpp>
#include <ExpressionTemplates/AssociativeTraits.hpp>

namespace expt
{
    /// \brief Performas a forward inverse transform.
    ///
    /// A forward inverse transform takes operators that are not 
    /// commutative or associative, and replaces them with their
    /// inverses, if their inverse is commutative or associative.
    /// This allows for better optimizations.  For example, the 
    /// expression A - (B-C) will require a temporary, since -
    /// is not commutative nor associative.  However, after rewriting
    /// as A + (-B) + C, then restructuring, and finally applying 
    /// the backward transform, we get A - B + C, which can be evaluated
    /// with no temporaries.
    template<typename NodeType>
    struct ForwardInverseTransform;

    /////////////////////////////////////////////////////////
    // Constant Node
    /////////////////////////////////////////////////////////
    template<typename DataType>
    struct ForwardInverseTransform<Node<DataType, void, void> >
    {
        typedef Node<DataType, void, void> Type;
    };

    /////////////////////////////////////////////////////////
    // Unary Node
    /////////////////////////////////////////////////////////
    template<typename ChildNodeType, typename UnaryOperatorType>
    struct ForwardInverseTransform<Node<ChildNodeType, UnaryOperatorType, void> >
    {
        typedef typename ForwardInverseTransform<ChildNodeType>::Type TransformedChildType;
        typedef Node<TransformedChildType, UnaryOperatorType, void> Type;
    };

    /////////////////////////////////////////////////////////
    // Binary Node
    /////////////////////////////////////////////////////////

    /// Even if the transform is mathematically valid, there are cases where the implementation
    /// does not provide the required unary operation, so A-B will compile, but A + (-B) will 
    /// not.  Specializations of this class will enable the forward transform to skip 
    /// these scenarios.
    template<typename OpType, typename DataType, typename enabled = void>
    struct HasUnaryOp : public boost::true_type {};

    namespace impl
    {
        template<typename OpType, typename DataType, typename enabled = void>
        struct OperatorCanBeInverted : public boost::false_type {};

        template<typename OpType, typename DataType>
        struct OperatorCanBeInverted<OpType, DataType, 
            typename boost::enable_if
            <
                boost::mpl::and_
                <
                    HasUnaryOp<typename InverseOperatorTypeTraits<OpType>::UnaryInversionType, DataType>,
                    boost::mpl::not_<AssociativeTraits<DataType, OpType, DataType, OpType> >,
                    boost::mpl::not_<CommutativeTraits<DataType, OpType, DataType> >,
                    boost::mpl::or_
                    <
                        AssociativeTraits<DataType, typename InverseOperatorTypeTraits<OpType>::InverseOpType, DataType, typename InverseOperatorTypeTraits<OpType>::InverseOpType>,
                        CommutativeTraits<DataType, typename InverseOperatorTypeTraits<OpType>::InverseOpType, DataType>
                    >
                >
            >::type> : public boost::true_type {};

        template<typename LhsType, typename OpType, typename RhsType, typename enabled=void>
        struct BinaryForwardAssociativeTransform;

        template<typename LhsType, typename OpType, typename RhsType>
        struct BinaryForwardAssociativeTransform<LhsType, OpType, RhsType,
            typename boost::enable_if
            <
                OperatorCanBeInverted<OpType, typename RhsType::ResultType>
            >::type>
        {
            typedef typename InverseOperatorTypeTraits<OpType>::InverseOpType InvertedOpType;
            typedef typename InvertNode<typename InverseOperatorTypeTraits<OpType>::UnaryInversionType, RhsType>::Type InvertedRhsType;
            typedef typename ForwardInverseTransform<InvertedRhsType>::Type TransformedRhsType;
            typedef typename ForwardInverseTransform<LhsType>::Type TransformedLhsType;

            typedef Node<TransformedLhsType, InvertedOpType, TransformedRhsType> Type;
        };

        template<typename LhsType, typename OpType, typename RhsType>
        struct BinaryForwardAssociativeTransform<LhsType, OpType, RhsType,
            typename boost::disable_if
            <
                OperatorCanBeInverted<OpType, typename RhsType::ResultType>
            >::type>
        {
            typedef typename ForwardInverseTransform<LhsType>::Type TransformedLhsType;
            typedef typename ForwardInverseTransform<RhsType>::Type TransformedRhsType;
            typedef Node<TransformedLhsType, OpType, TransformedRhsType> Type;
        };
    }

    template<typename LhsType, typename OpType, typename R1, typename ROp, typename R2>
    struct ForwardInverseTransform<Node<LhsType, OpType, Node<R1, ROp, R2> > >
    {
        typedef typename impl::BinaryForwardAssociativeTransform<LhsType, OpType, Node<R1, ROp, R2> >::Type Type;
    };
}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //EXPRESSION_TEMPLATES_FORWARD_ASSOCIATIVE_TRANSFORM_HPP
