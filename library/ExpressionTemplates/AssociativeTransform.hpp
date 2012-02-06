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

#ifndef EXPRESSION_TEMPLATES_ASSOCIATIVE_TRANSFORM_HPP
#define EXPRESSION_TEMPLATES_ASSOCIATIVE_TRANSFORM_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <ExpressionTemplates/Operators.hpp>
#include <ExpressionTemplates/AssociativeTraits.hpp>
#include <ExpressionTemplates/Node.hpp>

#include <boost/utility/enable_if.hpp>

namespace expt
{
    /// \brief Transforms the given tree using an associative transform if appropriate.
    ///
    /// Given a tree of the form A Op (B Op C), if Op is associative, this metafunction
    /// transforms it to (A Op B) Op C.
    template<typename NodeType, typename enabled = void>
    struct AssociativeTransform
    {
        typedef NodeType TransformedNodeType;
    };

    template<typename LhsNodeType, typename OpType, typename R1, typename ROp, typename R2>
    struct AssociativeTransform<Node<LhsNodeType, OpType, Node<R1, ROp, R2> >,
        typename boost::enable_if
        <
            AssociativeTraits
            <
                typename LhsNodeType::ResultType,
                OpType,
                typename R1::ResultType,
                ROp
            >
        >::type>
    {
        typedef Node<Node<LhsNodeType, OpType, R1>, ROp, R2> TransformedNodeType;
    };

    /// \brief Transforms the given tree using an associative transform if appropriate.
    ///
    /// Given a tree of the form (A Op B) Op C, if Op is associative, this metafunction
    /// transforms it to A Op (B Op C).
    template<typename NodeType, typename enabled = void>
    struct InverseAssociativeTransform
    {
        typedef NodeType TransformedNodeType;
    };

    template<typename L1, typename LOp, typename L2, typename OpType, typename RhsNodeType>
    struct InverseAssociativeTransform<Node<Node<L1, LOp, L2>, OpType, RhsNodeType >,
        typename boost::enable_if
        <
            AssociativeTraits
            <
                typename L1::ResultType,
                LOp,
                typename L2::ResultType,
                OpType
            >
        >::type>
    {
        typedef Node<L1, LOp, Node<L2, OpType, RhsNodeType> > TransformedNodeType;
    };
}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //EXPRESSION_TEMPLATES_ASSOCIATIVE_TRANSFORM_HPP

