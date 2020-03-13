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

#ifndef EXPRESSION_TEMPLATES_BACKWARD_ASSOCIATIVE_TRANSFORM_HPP
#define EXPRESSION_TEMPLATES_BACKWARD_ASSOCIATIVE_TRANSFORM_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

namespace expt
{
    /// \brief Reverses a forward transform.
    ///
    /// Given a tree of type A + (-B), where +/- are arbitrary inverse operators,
    /// and A and B are arbitrary expressions, this transform produces A - B.
    /// This transform does nothing if these conditions are not met.
    template<typename NodeType, typename IndicesType, unsigned int StartIndex, typename enbled = void>
    struct BackwardInverseTransform;

    ////////////////////////
    // Constant
    ////////////////////////
    template<typename DataType, typename IndicesType, unsigned int StartIndex>
    struct BackwardInverseTransform<Node<DataType, void, void>, IndicesType, StartIndex >
    {
        typedef Node<DataType, void, void> TransformedNodeType;
        typedef IndicesType TransformedIndicesType;
    };

    /////////////////////////////////////////////////////////
    // Unary Node
    /////////////////////////////////////////////////////////
    template<typename ChildNodeType, typename UnaryOperatorType, typename IndicesType, unsigned int StartIndex>
    struct BackwardInverseTransform<Node<ChildNodeType, UnaryOperatorType, void>, IndicesType, StartIndex >
    {
        typedef typename BackwardInverseTransform<ChildNodeType, IndicesType, StartIndex>::TransformedNodeType TransformedChildType;
        typedef typename BackwardInverseTransform<ChildNodeType, IndicesType, StartIndex>::TransformedIndicesType TransformedIndicesType;
        typedef Node<TransformedChildType, UnaryOperatorType, void> TransformedNodeType;
    };

    /////////////////////////////////////////////////////////
    // Binary Node
    /////////////////////////////////////////////////////////
    namespace impl
    {
        template<typename LhsType, typename OpType, typename RhsType, typename IndicesType, unsigned int StartIndex, typename enabled = void>
        struct BinaryBackwardInverseTransform
        {
            static const unsigned int RhsIndex = StartIndex + RhsType::TotalCount;

            typedef typename BackwardInverseTransform<LhsType, IndicesType, StartIndex>::TransformedNodeType TransformedLhsType;
            typedef typename BackwardInverseTransform<LhsType, IndicesType, StartIndex>::TransformedIndicesType Indices0;

            typedef typename BackwardInverseTransform<RhsType, Indices0, RhsIndex>::TransformedNodeType TransformedRhsType;
            typedef typename BackwardInverseTransform<RhsType, Indices0, RhsIndex>::TransformedIndicesType TransformedIndicesType;

            typedef Node<TransformedLhsType, OpType, TransformedRhsType> TransformedNodeType;
            
        };

        template<typename LhsType, typename OpType, typename RhsType, typename enabled = void>
        struct IsFixupNode : public boost::false_type {};

        template<typename L, typename LOp, typename OpType, typename RhsType>
        struct IsFixupNode<Node<L, LOp, void>, OpType, RhsType,
            typename boost::enable_if
            <
                boost::mpl::and_
                <
                    boost::is_same<LOp, typename InverseOperatorTypeTraits<OpType>::UnaryInversionType>,
                    CommutativeTraits<typename Node<L, LOp, void>::ResultType, OpType, typename RhsType::ResultType>,
                    IsConstantNode<L>
                >
            >::type> : public boost::true_type {};

        template<typename LhsType, typename OpType, typename R1, typename ROp, typename IndicesType, unsigned int StartIndex>
        struct BinaryBackwardInverseTransform<LhsType, OpType, Node<R1, ROp, void>, IndicesType, StartIndex,
            typename boost::enable_if
            <
                boost::mpl::and_
                <
                    boost::is_same<ROp, typename InverseOperatorTypeTraits<OpType>::UnaryInversionType>,
                    IsUnaryNode<Node<R1, ROp, void> >,
                    boost::mpl::not_<IsFixupNode<LhsType, OpType, Node<R1, ROp, void> > >
                >
            >::type>
        {
            static const unsigned int RhsIndex = StartIndex + R1::TotalCount;
            typedef typename InverseOperatorTypeTraits<OpType>::InverseOpType InvertedOpType;

            typedef typename BackwardInverseTransform<LhsType, IndicesType, StartIndex>::TransformedNodeType TransformedLhsType;
            typedef typename BackwardInverseTransform<LhsType, IndicesType, StartIndex>::TransformedIndicesType Indices0;

            typedef typename BackwardInverseTransform<R1, Indices0, RhsIndex>::TransformedNodeType TransformedRhsType;
            typedef typename BackwardInverseTransform<R1, Indices0, RhsIndex>::TransformedIndicesType TransformedIndicesType;

            typedef Node<TransformedLhsType, InvertedOpType, TransformedRhsType> TransformedNodeType;
        }; 

        // This case corresponds to something like
        // (-A) + BC.  The forward inverse transform creates negative nodes to maximize the 
        // size of clusters, but can create this situation where the negated node is moved to 
        // the far left side, even though optimally the tree should be BC-A.
        // This fixes that up.
        template<typename L, typename LOp, typename OpType, typename RhsType, typename IndicesType, unsigned int StartIndex>
        struct BinaryBackwardInverseTransform<Node<L, LOp, void>, OpType, RhsType, IndicesType, StartIndex,
            typename boost::enable_if
            <
                IsFixupNode<Node<L, LOp, void>, OpType, RhsType>
            >::type>
        {

            static const unsigned int RhsIndex = StartIndex + RhsType::TotalCount;
            typedef typename InverseOperatorTypeTraits<OpType>::InverseOpType InvertedOpType;
            typedef Node<L, LOp, void> LhsType;

            typedef typename BackwardInverseTransform<RhsType, IndicesType, RhsIndex>::TransformedNodeType TransformedRhsType;
            typedef typename BackwardInverseTransform<RhsType, IndicesType, RhsIndex>::TransformedIndicesType Indices0;

            typedef Node<LhsType, OpType, TransformedRhsType> NodeType0;
            typedef typename CommutativeTransform<NodeType0, Indices0, StartIndex>::TransformedIndicesType TransformedIndicesType;

            typedef Node<TransformedRhsType, InvertedOpType, L> TransformedNodeType;
        };
    }

    template<typename LhsType, typename OpType, typename RhsType,  typename IndicesType, unsigned int StartIndex>
    struct BackwardInverseTransform<Node<LhsType, OpType, RhsType>, IndicesType, StartIndex >
    {
        typedef typename impl::BinaryBackwardInverseTransform<LhsType, OpType, RhsType, IndicesType, StartIndex >::TransformedNodeType TransformedNodeType;
        typedef typename impl::BinaryBackwardInverseTransform<LhsType, OpType, RhsType,IndicesType, StartIndex>::TransformedIndicesType TransformedIndicesType;
    };

}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //EXPRESSION_TEMPLATES_BACKWARD_ASSOCIATIVE_TRANSFORM_HPP