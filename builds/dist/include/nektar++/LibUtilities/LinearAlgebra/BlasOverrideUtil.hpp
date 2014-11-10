///////////////////////////////////////////////////////////////////////////////
//
// For more information, please see: http://www.nektar.info
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
// Description: Defines the global functions needed for matrix operations.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BLAS_OVERRIDE_UTIL_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BLAS_OVERRIDE_UTIL_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

#include <ExpressionTemplates/Node.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

namespace expt
{
    template<typename NodeType, template <typename> class LhsTest, typename TestOpType,
        template <typename> class RhsTest, typename enabled=void>
    struct TestBinaryNode : public boost::false_type {};

    template<typename L, typename OpType, typename R, template <typename> class LhsTest,
        typename TestOpType, template <typename> class RhsTest>
    struct TestBinaryNode<Node<L, OpType, R>, LhsTest, TestOpType, RhsTest,
        typename boost::enable_if
        <
            boost::mpl::and_
            <
                boost::is_same<OpType, TestOpType>,
                LhsTest<typename L::ResultType>,
                RhsTest<typename R::ResultType>
            >
        >::type> : public boost::true_type {};

    template<typename NodeType, template <typename> class T1Type, typename TestOp1Type, template <typename> class T2Type,
        typename TestOp2Type, template <typename> class T3Type, typename enabled=void>
    struct Test3ArgumentAssociativeNode : public boost::false_type {};

    template<typename A1, typename Op1Type, typename A2, typename Op2Type, typename A3,
             template <typename> class T1Type, typename TestOp1Type, template <typename> class T2Type, typename TestOp2Type,
             template <typename> class T3Type>
    struct Test3ArgumentAssociativeNode<Node<Node<A1, Op1Type, A2>, Op2Type, A3>, T1Type, TestOp1Type, T2Type, TestOp2Type, T3Type,
        typename boost::enable_if
        <
            boost::mpl::and_
            <
                boost::is_same<Op1Type, TestOp1Type>,
                boost::is_same<Op2Type, TestOp2Type>,
                T1Type<typename A1::ResultType>,
                T2Type<typename A2::ResultType>,
                T3Type<typename A3::ResultType>
            >
        >::type> : public boost::true_type {};

}

#endif
#endif


