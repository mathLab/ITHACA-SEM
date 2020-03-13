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

using namespace expt;
namespace Nektar
{
    namespace UnitTests
    {
        BOOST_AUTO_TEST_CASE(TestVectorNegation)
        {
            typedef expt::Node<expt::Node<NekVector<double> >, expt::NegateOp > UnaryTree;

            double buf[] = {1.0, 2.0, 3.0, 4.0, 5.0};
            NekVector<double> v0(5, buf);

            //UnaryTree tree;
            UnaryTree tree = -v0;
            NekVector<double> v1 = -v0;

            double expected_result_buf[] = {-1.0, -2.0, -3.0, -4.0, -5.0};
            NekVector<double> expected_result(5, expected_result_buf);

            BOOST_CHECK_EQUAL(expected_result, v1);
        }

        BOOST_AUTO_TEST_CASE(TestUnaryInTree)
        {

            typedef expt::Node<expt::Node<NekVector<double> >, expt::NegateOp > UnaryTree;

            double buf0[] = {1.0, 2.0, 3.0, 4.0, 5.0};
            double buf1[] = {6.0, 7.0, 8.0, 9.0, 10.0};
            double buf2[] = {11.0, 12.0, 13.0, 14.0, 15.0};
            double buf3[] = {16.0, 17.0, 18.0, 19.0, 20.0};
            double buf4[] = {21.0, 22.0, 23.0, 24.0, 25.0};

            NekVector<double> v0(5, buf0);
            NekVector<double> v1(5, buf1);
            NekVector<double> v2(5, buf2);
            NekVector<double> v3(5, buf3);
            NekVector<double> v4(5, buf4);

            NekVector<double> result = (v0 + v1) - ( v2 + (-(v3+v4)));

            double expected_result_buf[5];
            for(unsigned int i = 0; i < 5; ++i)
            {
                expected_result_buf[i] = (buf0[i] + buf1[i]) - (buf2[i] + (-(buf3[i]+buf4[i])));
            }

            NekVector<double> expected_result(5, expected_result_buf);

            BOOST_CHECK_EQUAL(expected_result, result);
        }

        BOOST_AUTO_TEST_CASE(TestUnaryLhsEvaluatorConstLhs)
        {
            // Vector
            //{
            //    double buf0[] = {1.0, 2.0, 3.0, 4.0, 5.0};
            //    double buf1[] = {6.0, 7.0, 8.0, 9.0, 10.0};

            //    NekVector<double> v0(5, buf0);
            //    NekVector<double> v1(5, buf1);

            //    NekVector<double> v = v0*v1;

            //    double expectedResultBuf1[5];
            //    for(int i = 0; i < 5; ++i)
            //    {
            //        expectedResultBuf1[i] = buf0[i]*buf1[i];
            //    }
            //    NekVector<double> expectedResult1(5, expectedResultBuf1);
            //    BOOST_CHECK_EQUAL( expectedResult1, v);

            //    typedef Node<NekVector<double> > V;
            //    typedef Node<Node<V, NegateOp>, MultiplyOp, V> ExpressionType;
            //    ExpressionType tree = (-v0)*v1;
            //    NekVector<double> test2 = tree;
            //    double expectedResultBuf2[5];
            //    for(int i = 0; i < 5; ++i)
            //    {
            //        expectedResultBuf2[i] = -buf0[i]*buf1[i];
            //    }
            //    NekVector<double> expectedResult2(5, expectedResultBuf2);
            //    BOOST_CHECK_EQUAL( expectedResult2, test2);
            //}
        }
    }
}
