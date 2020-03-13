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
        BOOST_AUTO_TEST_CASE(ExpressionFromIncNavierStokes)
        {
            typedef Node<NekVector<double> > V;
            typedef Node<NekMatrix<double> > M;
            //typedef Node<V, NegateOp> Lhs;
            typedef V Lhs;
            typedef Node<M, MultiplyOp, V> Rhs;

            typedef Node<Lhs, SubtractOp, Rhs> TestNode;
            BOOST_MPL_ASSERT(( expt::impl::AlphaAXParameterAccess<Rhs, Rhs::Indices, 0> ));
            BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<V, SubtractOp, Node<M, MultiplyOp, V>, TestNode::Indices, 0> ));

            

            double buf_0[] = {1.0, 2.0};
            NekVector<double> v0(2, buf_0);

            double buf_1[] = {5.0, 6.0};
            NekVector<double> v1(2, buf_1);

            double m_buf[] = {9.0, 10.0, 11.0, 12.0};
            NekMatrix<double> m(2,2,m_buf);

            //NekVector<double> result = v0 - m*v1;
            NekVector<double> test1 = m*v1 + v0;
            double expected_buf1[] = {112, 124};
            NekVector<double> expected1(2, expected_buf1);
            BOOST_CHECK_EQUAL( test1, expected1);

            NekVector<double> test2 = m*v1 - v0;
            double expected_buf2[] = {110, 120};
            NekVector<double> expected2(2, expected_buf2);
            BOOST_CHECK_EQUAL( test2, expected2);

            NekVector<double> test3 = v0 - m*v1;
            typedef Node<Node<Rhs, NegateOp>, AddOp, V> ExpectedOptimizedTree;
            typedef Node<V, SubtractOp, Rhs> ExpectedT1;

            BOOST_MPL_ASSERT(( boost::is_same<expt::RemoveUnecessaryTemporaries<TestNode>::T1, ExpectedT1> ));
            BOOST_MPL_ASSERT(( boost::is_same<RemoveUnecessaryTemporaries<TestNode>::TransformedNodeType, ExpectedOptimizedTree> ));
            BOOST_MPL_ASSERT(( expt::impl::AlphaAXParameterAccess<Node<Rhs, NegateOp>, Node<Rhs, NegateOp>::Indices, 0> ));
            BOOST_MPL_ASSERT(( expt::impl::BetaYParameterAccess<V, V::Indices,0> ));
            BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<Node<Rhs, NegateOp>, AddOp, V, ExpectedOptimizedTree::Indices, 0> ));
            double expected_buf3[] = {-110, -120};
            NekVector<double> expected3(2, expected_buf3);
            BOOST_CHECK_EQUAL( test3, expected3);

            //double expected_result_buf[] = {1.0 - (45+66), 2.0 - (50 + 72)};
            //NekVector<double> expected_result(2, expected_result_buf);
            //BOOST_CHECK_EQUAL(expected_result, result);
        }

        BOOST_AUTO_TEST_CASE(ExpressionFromIncNavierStokesTypesOnly)
        {
            typedef Node<NekVector<double> > V;
            typedef Node<NekMatrix<double> > M;
            //typedef Node<V, NegateOp> Lhs;
            typedef V Lhs;
            typedef Node<M, MultiplyOp, V> Rhs;

            typedef Node<Lhs, SubtractOp, Rhs> TestNode;
            BOOST_MPL_ASSERT(( expt::impl::AlphaAXParameterAccess<Rhs, Rhs::Indices, 0> ));
            BOOST_MPL_ASSERT(( BinaryBinaryEvaluateNodeOverride<V, SubtractOp, Node<M, MultiplyOp, V>, TestNode::Indices, 0> ));

            typedef Node<Node<Rhs, NegateOp>, AddOp, V> ExpectedOptimizedTree;
            //typedef Node<V, AddOp, Node<Rhs, NegateOp> > ExpectedOptimizedTree;
            typedef Node<V, SubtractOp, Rhs> ExpectedT1;

            BOOST_MPL_ASSERT(( boost::is_same<expt::RemoveUnecessaryTemporaries<TestNode>::T1, ExpectedT1> ));
            BOOST_MPL_ASSERT(( boost::is_same<RemoveUnecessaryTemporaries<TestNode>::TransformedNodeType, ExpectedOptimizedTree> ));

            typedef RemoveUnecessaryTemporaries<TestNode>::TransformedIndicesType IndicesType;

            BOOST_STATIC_ASSERT(( boost::mpl::at_c<IndicesType, 0>::type::value == 1));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<IndicesType, 1>::type::value == 2));
            BOOST_STATIC_ASSERT(( boost::mpl::at_c<IndicesType, 2>::type::value == 0));
        }
    }
}

