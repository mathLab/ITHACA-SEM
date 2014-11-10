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

#ifndef EXPRESSION_TEMPLATES_PERFORM_COMMUTATIVE_TRANSFORM_IF_NEEDED_HPP
#define EXPRESSION_TEMPLATES_PERFORM_COMMUTATIVE_TRANSFORM_IF_NEEDED_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <ExpressionTemplates/CommutativeTraits.hpp>
#include <ExpressionTemplates/CommutativeTransform.hpp>
#include <ExpressionTemplates/Node.hpp>
#include <boost/type_traits.hpp>

namespace expt
{
    // Performs a commutative transform only if the transformation will remove a temporary.
    template<typename NodeType, typename IndicesType, unsigned int StartIndex, typename enabled = void>
    struct PerformCommutativeTransformIfNeeded
    {
        typedef NodeType TransformedNodeType;
        typedef IndicesType TransformedIndicesType;
    };
    
    // Expressions such as A + (BC), where a commutative transform will remove a temporary.
    template<typename LhsType, typename OpType, typename RhsType, typename IndicesType, unsigned int StartIndex>
    struct PerformCommutativeTransformIfNeeded< expt::Node<LhsType, OpType, RhsType>, IndicesType, StartIndex,
        typename boost::enable_if
        <
            boost::mpl::and_
            <
                expt::CommutativeTraits<typename LhsType::ResultType, OpType, typename RhsType::ResultType>,
                boost::is_same<typename LhsType::OpType, void>,
                boost::mpl::not_<boost::is_same<typename RhsType::OpType, void> >
            >
        >::type>
    {
        typedef expt::Node<LhsType, OpType, RhsType> NodeType;
        typedef typename expt::CommutativeTransform<NodeType, IndicesType, StartIndex>::TransformedNodeType TransformedNodeType;
        typedef typename expt::CommutativeTransform<NodeType, IndicesType, StartIndex>::TransformedIndicesType TransformedIndicesType;
    };

}

#endif
#endif

