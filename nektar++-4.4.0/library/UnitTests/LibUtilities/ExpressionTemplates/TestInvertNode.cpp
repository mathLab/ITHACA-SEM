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
#include <ExpressionTemplates/InvertNode.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

using namespace expt;

namespace Nektar
{
    namespace UnitTests
    {
        BOOST_AUTO_TEST_CASE(TestConstantTransformation)
        {
            typedef Node<NekVector<double> > LeafNode;
            typedef InvertNode<NegateOp, LeafNode>::Type Test1;

            typedef Node<LeafNode, NegateOp> ExpectedTest1;

            BOOST_MPL_ASSERT(( boost::is_same<Test1, ExpectedTest1> ));

            typedef InvertNode<InvertOp, LeafNode>::Type Test2;
            typedef Node<LeafNode, InvertOp> ExpectedTest2;

            BOOST_MPL_ASSERT(( boost::is_same<Test2, ExpectedTest2> ));
        }

        BOOST_AUTO_TEST_CASE(TestUnaryTransformation)
        {
            typedef Node<NekVector<double> > LeafNode;
            typedef Node<LeafNode, NegateOp> NegateNode;

            typedef InvertNode<NegateOp, NegateNode>::Type Test1;
            BOOST_MPL_ASSERT(( boost::is_same<Test1, LeafNode>));

            typedef Node<LeafNode, InvertOp> InvertNodeType;
            typedef InvertNode<NegateOp, InvertNodeType>::Type Test2;
            typedef Node<InvertNodeType, NegateOp> ExpectedTest2;
            BOOST_MPL_ASSERT(( boost::is_same<Test2, ExpectedTest2> ));

        }

        BOOST_AUTO_TEST_CASE(TestBinaryInversion)
        {
            typedef Node<NekVector<double> > LeafNode;
            typedef Node<LeafNode, AddOp, LeafNode> Test1;
            typedef Node<LeafNode, SubtractOp, LeafNode> Test2;
            typedef Node<LeafNode, MultiplyOp, LeafNode> Test3;

            typedef InvertNode<NegateOp, Test1>::Type Test1Result;
            typedef InvertNode<NegateOp, Test2>::Type Test2Result;
            typedef InvertNode<NegateOp, Test3>::Type Test3Result;

            typedef Node<Node<LeafNode, NegateOp>, SubtractOp, LeafNode> ExpectedTest1Result;
            typedef Node<Node<LeafNode, NegateOp>, AddOp, LeafNode> ExpectedTest2Result;
            typedef Node<Test3, NegateOp> ExpectedTest3Result;

            BOOST_MPL_ASSERT(( boost::is_same<expt::impl::InvertBinaryNode<NegateOp, LeafNode, AddOp, LeafNode>::TransformedOpType, SubtractOp > ));
            BOOST_MPL_ASSERT(( boost::is_same<expt::impl::InvertBinaryNode<NegateOp, LeafNode, AddOp, LeafNode>::TransformedLhsType, Node<LeafNode, NegateOp> > ));
            BOOST_MPL_ASSERT(( boost::is_same<Test1Result, ExpectedTest1Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test2Result, ExpectedTest2Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test3Result, ExpectedTest3Result> ));
         
            typedef InvertNode<NegateOp, Test1Result>::Type InverseTest1Result;
            typedef InvertNode<NegateOp, Test2Result>::Type InverseTest2Result;
            typedef InvertNode<NegateOp, Test3Result>::Type InverseTest3Result;

            BOOST_MPL_ASSERT(( boost::is_same<Test1, InverseTest1Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test2, InverseTest2Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test3, InverseTest3Result> ));
        }

        BOOST_AUTO_TEST_CASE(TestRecursion)
        {
            typedef Node<NekVector<double> > LeafNode;
            typedef Node<LeafNode, AddOp, LeafNode> AddNode;
            typedef Node<LeafNode, SubtractOp, LeafNode> SubtractNode;
            typedef Node<AddNode, AddOp, SubtractNode> TestType;

            typedef InvertNode<NegateOp, TestType>::Type ResultType;

            typedef Node<LeafNode, NegateOp> NegateNode;
            typedef Node<NegateNode, SubtractOp, LeafNode> ExpectedLhsNode;
            typedef Node<ExpectedLhsNode, SubtractOp, SubtractNode> ExpectedResultType;

            BOOST_MPL_ASSERT(( boost::is_same<ResultType, ExpectedResultType> ));

        }
    }
}
