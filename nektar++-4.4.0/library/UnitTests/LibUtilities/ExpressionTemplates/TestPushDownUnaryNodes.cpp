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
///////////////////////////////////////////////////////////////////////////////

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

using namespace expt;

namespace Nektar
{
    namespace UnitTests
    {
        BOOST_AUTO_TEST_CASE(TestConstant)
        {
            typedef Node<NekVector<double> > LeafNode;
            typedef PushDownUnaryNodes<LeafNode>::Type ResultType;
            BOOST_MPL_ASSERT(( boost::is_same<ResultType, LeafNode> ));
        }

        BOOST_AUTO_TEST_CASE(TestOneLevelUnary)
        {
            typedef Node<NekVector<double> > LeafNode;
            typedef Node<LeafNode, NegateOp> NegateNode;
            typedef Node<LeafNode, InvertOp> InvertNode;

            typedef PushDownUnaryNodes<NegateNode>::Type NegateResult;
            typedef PushDownUnaryNodes<InvertNode>::Type InvertResult;

            BOOST_MPL_ASSERT(( boost::is_same<NegateNode, NegateResult> ));
            BOOST_MPL_ASSERT(( boost::is_same<InvertNode, InvertResult> ));

        }

        BOOST_AUTO_TEST_CASE(TestTwoLevelUnary)
        {
            typedef Node<NekVector<double> > LeafNode;
            typedef Node<LeafNode, AddOp, LeafNode> AddNode;
            typedef Node<LeafNode, MultiplyOp, LeafNode> MultiplyNode;
            typedef Node<AddNode, NegateOp> Test1Type;
            typedef Node<MultiplyNode, NegateOp> Test2Type;

            typedef PushDownUnaryNodes<Test1Type>::Type Test1Result;

            typedef Node<Node<LeafNode, NegateOp>, SubtractOp, LeafNode> ExpectedTest1Result;
            BOOST_MPL_ASSERT(( boost::is_same<Test1Result, ExpectedTest1Result> ));

            typedef PushDownUnaryNodes<Test2Type>::Type Test2Result;
            typedef Node<MultiplyNode, NegateOp> ExpectedTest2Result;
            BOOST_MPL_ASSERT(( boost::is_same<Test2Result, ExpectedTest2Result> ));

        }
    }
}
