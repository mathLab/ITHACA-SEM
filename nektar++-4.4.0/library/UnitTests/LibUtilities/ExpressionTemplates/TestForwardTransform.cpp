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
        BOOST_AUTO_TEST_CASE(TestConstantForwardTransform)
        {
            typedef Node<NekVector<double> > LeafNode;
            typedef ForwardInverseTransform<LeafNode>::Type Type;
            BOOST_MPL_ASSERT(( boost::is_same<LeafNode, Type> ));
        }

        BOOST_AUTO_TEST_CASE(TestUnaryForwardTransform)
        {
            typedef Node<NekVector<double> > LeafNode;
            typedef Node<LeafNode, NegateOp> NegateNod;
            typedef ForwardInverseTransform<NegateNod>::Type Type;
            BOOST_MPL_ASSERT(( boost::is_same<NegateNod, Type> ));
        }

        BOOST_AUTO_TEST_CASE(TestBinaryForwardTransformConstantRhs)
        {
            typedef Node<NekVector<double> > LeafNode;
            typedef Node<LeafNode, AddOp, LeafNode> AddNode;
            typedef Node<LeafNode, SubtractOp, LeafNode> SubtractNode;
            typedef Node<LeafNode, MultiplyOp, LeafNode> MultiplyNode;
            typedef Node<LeafNode, DivideOp, LeafNode> DivideNode;

            typedef ForwardInverseTransform<AddNode>::Type AddNodeResult;
            BOOST_MPL_ASSERT(( boost::is_same<AddNode, AddNodeResult> ));

            typedef ForwardInverseTransform<SubtractNode>::Type SubtractNodeResult;
            typedef Node<LeafNode, AddOp, Node<LeafNode, NegateOp> > ExpectedSubtractNodeResult;
            BOOST_MPL_ASSERT(( boost::is_same<SubtractNodeResult, ExpectedSubtractNodeResult> ));

            typedef ForwardInverseTransform<MultiplyNode>::Type MultiplyNodeResult;
            BOOST_MPL_ASSERT(( boost::is_same<MultiplyNodeResult, MultiplyNode> ));

            typedef ForwardInverseTransform<DivideNode>::Type DivideNodeResult;
            typedef Node<LeafNode, MultiplyOp, Node<LeafNode, InvertOp> > ExpectedDivideNodeResult;
            BOOST_MPL_ASSERT(( boost::is_same<ExpectedDivideNodeResult, DivideNodeResult> ));
        }

        BOOST_AUTO_TEST_CASE(TestNoInitalBackTransformRationale)
        {
            typedef Node<NekVector<double> > LeafNode;
            typedef Node<LeafNode, SubtractOp, Node<LeafNode, NegateOp> > TestNode;

            typedef ForwardInverseTransform<TestNode>::Type TestResult;
            typedef Node<LeafNode, AddOp, LeafNode> ExpectedTestResult;
            BOOST_MPL_ASSERT(( boost::is_same<TestResult, ExpectedTestResult> ));
        }

        BOOST_AUTO_TEST_CASE(TestBinaryForwardTransformUnaryRhs)
        {
            typedef Node<NekVector<double> > LeafNode;
            typedef Node<LeafNode, NegateOp> NegateNode;
            typedef Node<LeafNode, InvertOp> InvertNode;

            typedef Node<LeafNode, AddOp, NegateNode> Test1;
            typedef Node<LeafNode, SubtractOp, NegateNode> Test2;
            typedef Node<LeafNode, MultiplyOp, NegateNode> Test3;
            typedef Node<LeafNode, DivideOp, NegateNode> Test4;
            typedef Node<LeafNode, AddOp, InvertNode> Test5;
            typedef Node<LeafNode, SubtractOp, InvertNode> Test6;
            typedef Node<LeafNode, MultiplyOp, InvertNode> Test7;
            typedef Node<LeafNode, DivideOp, InvertNode> Test8;

            typedef ForwardInverseTransform<Test1>::Type Test1Result;
            typedef ForwardInverseTransform<Test2>::Type Test2Result;
            typedef ForwardInverseTransform<Test3>::Type Test3Result;
            typedef ForwardInverseTransform<Test4>::Type Test4Result;
            typedef ForwardInverseTransform<Test5>::Type Test5Result;
            typedef ForwardInverseTransform<Test6>::Type Test6Result;
            typedef ForwardInverseTransform<Test7>::Type Test7Result;
            typedef ForwardInverseTransform<Test8>::Type Test8Result;

            typedef Test1 ExpectedTest1Result;
            typedef Node<LeafNode, AddOp, LeafNode> ExpectedTest2Result;
            typedef Test3 ExpectedTest3Result;
            typedef Node<LeafNode, MultiplyOp, Node<Node<LeafNode, NegateOp>, InvertOp> > ExpectedTest4Result;
            typedef Test5 ExpectedTest5Result;
            typedef Node<LeafNode, AddOp, Node<Node<LeafNode, InvertOp>, NegateOp> > ExpectedTest6Result;
            typedef Test7 ExpectedTest7Result;
            typedef Node<LeafNode, MultiplyOp, LeafNode> ExpectedTest8Result;

            BOOST_MPL_ASSERT(( boost::is_same<Test1Result, ExpectedTest1Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test2Result, ExpectedTest2Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test3Result, ExpectedTest3Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test4Result, ExpectedTest4Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test5Result, ExpectedTest5Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test6Result, ExpectedTest6Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test7Result, ExpectedTest7Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test8Result, ExpectedTest8Result> ));
        }

        BOOST_AUTO_TEST_CASE(TestBinaryForwardTransformBinaryRhsAdditionRootOp)
        {
            typedef Node<NekVector<double> > LeafNode;
            typedef Node<LeafNode, AddOp, LeafNode> AddNode;
            typedef Node<LeafNode, SubtractOp, LeafNode> SubtractNode;
            typedef Node<LeafNode, MultiplyOp, LeafNode> MultiplyNode;
            typedef Node<LeafNode, DivideOp, LeafNode> DivideNode;

            typedef Node<LeafNode, AddOp, AddNode> Test1;
            typedef Node<LeafNode, AddOp, SubtractNode> Test2;
            typedef Node<LeafNode, AddOp, MultiplyNode> Test3;

            typedef ForwardInverseTransform<Test1>::Type Test1Result;
            typedef ForwardInverseTransform<Test2>::Type Test2Result;
            typedef ForwardInverseTransform<Test3>::Type Test3Result;

            typedef Test1 ExpectedTest1Result;
            typedef Node<LeafNode, AddOp, Node<LeafNode, AddOp, Node<LeafNode, NegateOp> > > ExpectedTest2Result;
            typedef Test3 ExpectedTest3Result;
            typedef Node<LeafNode, AddOp, Node<LeafNode, MultiplyOp, Node<LeafNode, InvertOp> > > ExpectedTest4Result;

            BOOST_MPL_ASSERT(( boost::is_same<Test1Result, ExpectedTest1Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test2Result, ExpectedTest2Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test3Result, ExpectedTest3Result> ));
        }

        BOOST_AUTO_TEST_CASE(TestBinaryForwardTransformBalancedTwoLevelBinary)
        {
            typedef Node<NekVector<double> > LeafNode;
            typedef Node<LeafNode, AddOp, LeafNode> AddNode;
            typedef Node<LeafNode, SubtractOp, LeafNode> SubtractNode;
            typedef Node<LeafNode, MultiplyOp, LeafNode> MultiplyNode;
            typedef Node<LeafNode, DivideOp, LeafNode> DivideNode;
            typedef Node<LeafNode, NegateOp> NegateNode;

            typedef Node<AddNode, SubtractOp, AddNode> Test1;
            typedef Node<SubtractNode, SubtractOp, SubtractNode> Test2;

            typedef ForwardInverseTransform<Test1>::Type Test1Result;
            typedef ForwardInverseTransform<Test2>::Type Test2Result;


            typedef Node<AddNode, AddOp, Node<NegateNode, AddOp, NegateNode> > ExpectedTest1Result;
            typedef Node<LeafNode, AddOp, NegateNode> ForwardTransformOfSubractNode;
            typedef Node<ForwardTransformOfSubractNode, AddOp, Node<NegateNode, AddOp, LeafNode> > ExpectedTest2Result;

            BOOST_MPL_ASSERT(( boost::is_same<Test1Result, ExpectedTest1Result> ));
            BOOST_MPL_ASSERT(( boost::is_same<Test2Result, ExpectedTest2Result> ));
        }

        BOOST_AUTO_TEST_CASE(TestMMinusMM)
        {
            double buf_0[] = {1.0, 2.0, 3.0, 4.0};
            double buf_1[] = {5.0, 6.0, 7.0, 8.0};
            double buf_2[] = {9.0, 10.0, 11.0, 12.0};

            NekMatrix<double> m0(2,2,buf_0);
            NekMatrix<double> m1(2,2,buf_1);
            NekMatrix<double> m2(2,2,buf_2);

            NekMatrix<double> result = m0 - m1*m2;

            double expected_result_buf[] = {1.0 - (5*9+7*10), 2.0 - (6*9+8*10), 3.0 - (5*11+7*12), 4.0 - (6*11+8*12)};
            NekMatrix<double> expected_result(2,2,expected_result_buf);

            BOOST_CHECK_EQUAL(result, expected_result);

        }

        BOOST_AUTO_TEST_CASE(TestForwardAssociative0)
        {
            double m_buf_0[] = {1.0, 2.0, 3.0, 4.0};
            double m_buf_1[] = {5.0, 6.0, 7.0, 8.0};
            double v_buf_0[] = {9.0, 10.0};
            double v_buf_1[] = {11.0, 12.0};
            double v_buf_2[] = {13.0, 14.0};

            NekMatrix<double> m0(2,2,m_buf_0);
            NekMatrix<double> m1(2,2,m_buf_1);
            NekVector<double> v0(2,v_buf_0);
            NekVector<double> v1(2,v_buf_1);
            NekVector<double> v2(2,v_buf_2);

            NekVector<double> result = v0+m0*v1 - m1*v2;

            double expected_result_buf[] = {-107, -110};
            NekVector<double> expected_result(2, expected_result_buf);

            BOOST_CHECK_EQUAL(expected_result, result);

        }

        BOOST_AUTO_TEST_CASE(TestDisabledBlockMatrixInversion)
        {
            typedef NekMatrix<double, StandardMatrixTag> Standard;
            typedef NekMatrix<Standard, BlockMatrixTag> Block;

            BOOST_MPL_ASSERT(( boost::mpl::not_<IsBlockMatrix<Standard> >));
            BOOST_MPL_ASSERT((IsBlockMatrix<Block>));

            BOOST_MPL_ASSERT(( boost::mpl::not_<HasUnaryOp<NegateOp, Block> >));
            BOOST_MPL_ASSERT((HasUnaryOp<NegateOp, Standard>));

            BOOST_MPL_ASSERT(( expt::impl::OperatorCanBeInverted<SubtractOp, Standard> ));
            BOOST_MPL_ASSERT(( boost::mpl::not_<expt::impl::OperatorCanBeInverted<SubtractOp, Block> > ));
        }
    }
}
