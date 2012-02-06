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
        BOOST_AUTO_TEST_CASE(TestLeafForwardTransform)
        {
        }

        BOOST_AUTO_TEST_CASE(TestOneLevelUnaryForwardTransform)
        {
        }

        BOOST_AUTO_TEST_CASE(TestOneLevelBinaryForwardTransform)
        {
            //typedef expt::Node<NekVector<double> > LeafNode;
            //typedef expt::Node<LeafNode, expt::AddOp, LeafNode> AddNode;
            //typedef expt::Node<LeafNode, expt::SubtractOp, LeafNode> SubtractNode;
            //typedef expt::Node<LeafNode, expt::MultiplyOp, LeafNode> MultiplyNode;
            //typedef expt::Node<LeafNode, expt::DivideOp, LeafNode> DivideNode;

            //BOOST_STATIC_ASSERT( expt::ForwardTransform<AddNode>::Id == 5);
            //typedef expt::ForwardTransform<AddNode>::Type TransformedAddNode;
            //typedef expt::ForwardTransform<SubtractNode>::Type TransformedSubtractNode;
            //typedef expt::ForwardTransform<MultiplyNode>::Type TransformedMultiplyNode;
            //typedef expt::ForwardTransform<DivideNode>::Type TransformedDivideNode;

            //BOOST_MPL_ASSERT(( boost::is_same<AddNode, TransformedAddNode> ));
            //BOOST_MPL_ASSERT(( boost::is_same<SubtractNode, TransformedSubtractNode> ));
            //BOOST_MPL_ASSERT(( boost::is_same<MultiplyNode, TransformedMultiplyNode> ));
            //BOOST_MPL_ASSERT(( boost::is_same<DivideNode, TransformedDivideNode> ));
        }

        BOOST_AUTO_TEST_CASE(TestLeftHeavyTwoLevelBinaryTransform)
        {
            //typedef expt::Node<NekVector<double> > L;
            //typedef expt::Node<L, expt::AddOp, L> AddNode;
            //typedef expt::Node<L, expt::SubtractOp, L> SubtractNode;
            //typedef expt::Node<L, expt::MultiplyOp, L> MultiplyNode;
            //typedef expt::Node<L, expt::DivideOp, L> DivideNode;
            //typedef expt::Node<L, expt::NegateOp> NegateNode;

            //typedef expt::Node<AddNode, expt::AddOp, L> Test1;
            //typedef expt::Node<SubtractNode, expt::AddOp, L> Test2;
            //typedef expt::Node<MultiplyNode, expt::AddOp, L> Test3;
            //typedef expt::Node<DivideNode, expt::AddOp, L> Test4;
            //typedef expt::Node<NegateNode, expt::AddOp, L> Test5;

            //typedef expt::ForwardTransform<Test1>::Type Test1Result;
            //typedef expt::ForwardTransform<Test2>::Type Test2Result;
            //typedef expt::ForwardTransform<Test3>::Type Test3Result;
            //typedef expt::ForwardTransform<Test4>::Type Test4Result;
            //typedef expt::ForwardTransform<Test5>::Type Test5Result;

            //typedef expt::Node<L, expt::AddOp, NegateNode> NegateFwdTransform;
            //typedef expt::Node<L, expt::InvertOp> InverseNode;
            //typedef expt::Node<L, expt::MultiplyOp, InverseNode> DivideFwdTransform;

            //BOOST_MPL_ASSERT(( boost::is_same<Test1Result, Test1> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Test2Result, Test2> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Test3Result, Test3> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Test4Result, Test4> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Test5Result, Test5> ));
            //            
            //typedef expt::Node<AddNode, expt::SubtractOp, L> Test6;
            //typedef expt::Node<SubtractNode, expt::SubtractOp, L> Test7;
            //typedef expt::Node<MultiplyNode, expt::SubtractOp, L> Test8;
            //typedef expt::Node<DivideNode, expt::SubtractOp, L> Test9;
            //typedef expt::Node<NegateNode, expt::SubtractOp, L> Test10;

            //typedef expt::ForwardTransform<Test6>::Type Test6Result;
            //typedef expt::ForwardTransform<Test7>::Type Test7Result;
            //typedef expt::ForwardTransform<Test8>::Type Test8Result;
            //typedef expt::ForwardTransform<Test9>::Type Test9Result;
            //typedef expt::ForwardTransform<Test10>::Type Test10Result;

            //BOOST_MPL_ASSERT(( boost::is_same<Test6Result, Test6> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Test7Result, Test7> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Test8Result, Test8> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Test9Result, Test9> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Test10Result, Test10> ));
        }

        BOOST_AUTO_TEST_CASE(TestRightHeavyTwoLevelBinaryTransform)
        {
            //typedef expt::Node<NekVector<double> > L;
            //typedef expt::Node<L, expt::AddOp, L> AddNode;
            //typedef expt::Node<L, expt::SubtractOp, L> SubtractNode;
            //typedef expt::Node<L, expt::MultiplyOp, L> MultiplyNode;
            //typedef expt::Node<L, expt::DivideOp, L> DivideNode;
            //typedef expt::Node<L, expt::NegateOp> NegateNode;

            //typedef expt::Node<L, expt::AddOp, AddNode> Test1;
            //typedef expt::Node<L, expt::AddOp, SubtractNode> Test2;
            //typedef expt::Node<L, expt::AddOp, MultiplyNode> Test3;
            //typedef expt::Node<L, expt::AddOp, DivideNode> Test4;
            //typedef expt::Node<L, expt::AddOp, NegateNode> Test5;

            //typedef expt::ForwardTransform<Test1>::Type Test1Result;
            //typedef expt::ForwardTransform<Test2>::Type Test2Result;
            //typedef expt::ForwardTransform<Test3>::Type Test3Result;
            //typedef expt::ForwardTransform<Test4>::Type Test4Result;
            //typedef expt::ForwardTransform<Test5>::Type Test5Result;

            //typedef expt::Node<L, expt::AddOp, NegateNode> NegateFwdTransform;
            //typedef expt::Node<L, expt::InvertOp> InverseNode;
            //typedef expt::Node<L, expt::MultiplyOp, InverseNode> DivideFwdTransform;

            //typedef Test1 ExpectedTest1Result;
            //typedef expt::Node<L, expt::AddOp, NegateFwdTransform> ExpectedTest2Result;
            //typedef Test3 ExpectedTest3Result;
            //typedef expt::Node<L, expt::AddOp, DivideFwdTransform> ExpectedTest4Result;
            //typedef Test5 ExpectedTest5Result;

            //BOOST_MPL_ASSERT(( boost::is_same<Test1Result, ExpectedTest1Result> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Test2Result, ExpectedTest2Result> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Test3Result, ExpectedTest3Result> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Test4Result, ExpectedTest4Result> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Test5Result, ExpectedTest5Result> ));

            //typedef expt::Node<L, expt::SubtractOp, AddNode> Test6;
            //typedef expt::Node<L, expt::SubtractOp, SubtractNode> Test7;
            //typedef expt::Node<L, expt::SubtractOp, MultiplyNode> Test8;
            //typedef expt::Node<L, expt::SubtractOp, DivideNode> Test9;
            //typedef expt::Node<L, expt::SubtractOp, NegateNode> Test10;

            //typedef expt::ForwardTransform<Test6>::Type Test6Result;
            //typedef expt::ForwardTransform<Test7>::Type Test7Result;
            //typedef expt::ForwardTransform<Test8>::Type Test8Result;
            //typedef expt::ForwardTransform<Test9>::Type Test9Result;
            //typedef expt::ForwardTransform<Test10>::Type Test10Result;

            //typedef expt::Node<L, expt::AddOp, expt::Node<NegateNode, expt::AddOp, NegateNode> > ExpectedTest6Result;
            //typedef expt::Node<L, expt::AddOp, expt::Node<NegateNode, expt::AddOp, L> > ExpectedTest7Result;
            //typedef expt::Node<L, expt::AddOp, expt::Node<MultiplyNode, expt::NegateOp> > ExpectedTest8Result;
            //typedef expt::Node<L, expt::AddOp, expt::Node<expt::Node<L, expt::MultiplyOp, InverseNode>, expt::NegateOp> > ExpectedTest9Result;
            //typedef expt::Node<L, expt::AddOp, expt::Node<NegateNode, expt::NegateOp> > ExpectedTest10Result;

            //BOOST_MPL_ASSERT(( boost::is_same<Test6Result, ExpectedTest6Result> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Test7Result, ExpectedTest7Result> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Test8Result, ExpectedTest8Result> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Test9Result, ExpectedTest9Result> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Test10Result, ExpectedTest10Result> ));
        }

        BOOST_AUTO_TEST_CASE(TestBalancedTwoLevelBinaryTransform)
        {
        }

        BOOST_AUTO_TEST_CASE(TestBackwardTransform)
        {
            //typedef expt::Node<NekVector<double> > L;
            //typedef expt::Node<L, expt::AddOp, L> AddNode;
            //typedef expt::Node<L, expt::SubtractOp, L> SubtractNode;
            //typedef expt::Node<L, expt::MultiplyOp, L> MultiplyNode;
            //typedef expt::Node<L, expt::DivideOp, L> DivideNode;
            //typedef expt::Node<L, expt::NegateOp> NegateNode;
            //typedef expt::Node<L, expt::InvertOp> InvertNode;

            //typedef expt::Node<L, expt::AddOp, NegateNode> Test1;
            //typedef expt::Node<L, expt::AddOp, InvertNode> Test2;
            //typedef expt::Node<L, expt::SubtractOp, NegateNode> Test3;
            //typedef expt::Node<L, expt::SubtractOp, InvertNode> Test4;
            //typedef expt::Node<L, expt::AddOp, Test1> Test5;

            //typedef expt::BackwardsTransform<Test1>::Type Test1Result;
            //typedef expt::BackwardsTransform<Test2>::Type Test2Result;
            //typedef expt::BackwardsTransform<Test3>::Type Test3Result;
            //typedef expt::BackwardsTransform<Test4>::Type Test4Result;
            //typedef expt::BackwardsTransform<Test5>::Type Test5Result;
            //BOOST_STATIC_ASSERT(( expt::BackwardsTransform<Test5>::Id == 10));
            //BOOST_MPL_ASSERT(( boost::is_same<expt::BackwardsTransform<Test5>::TransformedLhsType, L> ));
            //BOOST_MPL_ASSERT(( boost::is_same<expt::BackwardsTransform<Test5>::TransformedRhsType, Test1Result> ));

            //typedef SubtractNode ExpectedTest1Result;
            //typedef Test2 ExpectedTest2Result;
            //typedef AddNode ExpectedTest3Result;
            //typedef Test4 ExpectedTest4Result;
            //typedef expt::Node<L, expt::AddOp, SubtractNode> ExpectedTest5Result;

            //BOOST_MPL_ASSERT(( boost::is_same<Test1Result, ExpectedTest1Result> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Test2Result, ExpectedTest2Result> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Test3Result, ExpectedTest3Result> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Test4Result, ExpectedTest4Result> ));
            //BOOST_MPL_ASSERT(( boost::is_same<Test5Result, ExpectedTest5Result> ));

        }

        BOOST_AUTO_TEST_CASE(Test0)
        {
            double buf0[] = {1.0, 2.0};
            //double buf1[] = {3.0, 4.0};
            //double m_buf[] = {5.0, 6.0, 7.0, 8.0};

            //NekVector<double> v0(2, buf0);
            //NekVector<double> v1(2, buf1);
            //NekMatrix<double> m(2, 2, m_buf);

            //typedef expt::Node<NekVector<double> > VectorLeaf;
            //typedef expt::Node<NekMatrix<double> >  MatrixLeaf;
            //typedef expt::Node<VectorLeaf, expt::SubtractOp, expt::Node<MatrixLeaf, expt::MultiplyOp, VectorLeaf> > ExpectedType;

            //ExpectedType e = v0 - m*v1;

            //typedef expt::RemoveUnecessaryTemporaries<ExpectedType, ExpectedType::Indices, 0>::TransformedNodeType ResultExpNode;

            //typedef expt::Node<expt::Node<MatrixLeaf, expt::MultiplyOp, VectorLeaf>, expt::NegateOp> ExpectedLhs;
            //typedef expt::Node<ExpectedLhs, expt::AddOp, VectorLeaf> ExpectedTransformedType;

            //BOOST_MPL_ASSERT(( boost::is_same<ExpectedTransformedType, ResultExpNode> ));

            //NekVector<double> result = v0 - m*v1;
            //double expected_result_buf[] = {1.0 - (15+28), 2.0 - (18+32)};
            //NekVector<double> expected_result(2, expected_result_buf);
            //BOOST_CHECK_EQUAL(result, expected_result);

        }

        BOOST_AUTO_TEST_CASE(Test1)
        {
            //double m_buf[] = {5.0, 6.0, 7.0, 8.0};
            //double s_buf[] = {5.0, 6.0, 7.0, 8.0};

            //NekMatrix<NekDouble> m(2,2,m_buf);
            //boost::shared_ptr<NekMatrix<NekDouble> > inner(new NekMatrix<NekDouble>(2,2,s_buf));
            //DNekScalMat scaled(2.0, inner);

            //NekMatrix<NekDouble> result = m - scaled;
            //typedef expt::Node<NekMatrix<NekDouble> > MatrixNode;
            //typedef expt::Node<DNekScalMat> ScaledMatrixNode;

            //typedef expt::Node<MatrixNode, expt::SubtractOp, ScaledMatrixNode> InitialTree;
            //InitialTree e = m-scaled;

            //typedef expt::BackwardsTransform<InitialTree>::Type Step1Tree;
            //typedef expt::Node<MatrixNode, expt::AddOp, expt::Node<ScaledMatrixNode, expt::NegateOp> > ExpectedStep1Tree;
            //BOOST_MPL_ASSERT(( boost::is_same<Step1Tree, InitialTree> ));
        }

    }
}
