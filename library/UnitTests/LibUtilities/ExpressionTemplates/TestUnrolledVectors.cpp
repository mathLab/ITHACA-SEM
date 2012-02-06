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
#ifndef NEKTAR_USE_EXPRESSION_TEMPLATES
#define NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //NEKTAR_USE_EXPRESSION_TEMPLATES

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <UnitTests/LibUtilities/ExpressionTemplates/CountedObjectExpression.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <ExpressionTemplates/Node.hpp>
#include <ExpressionTemplates/RemoveAllUnecessaryTemporaries.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits.hpp>
#include <LibUtilities/LinearAlgebra/MatrixSize.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>

namespace Nektar
{
    namespace UnitTests
    {
        BOOST_AUTO_TEST_CASE(Test2VectorUnrolled)
        {
            typedef NekVector<NekDouble> Vector;

            double buf1[] = {1.0, 2.0, 3.0, 4.0};
            double buf2[] = {5.0, 6.0, 7.0, 8.0};
            
            Vector v1(4, buf1);
            Vector v2(4, buf2);
            
            Vector result = v1+v2;

            double expectedResultBuf[] = {1.0 + 5.0,
                2.0 + 6.0, 3.0 + 7.0, 4.0 + 8.0};
            Vector expectedResult(4, expectedResultBuf);
            
            BOOST_CHECK_EQUAL(result, expectedResult);
        }

        BOOST_AUTO_TEST_CASE(Test3VectorUnrolled)
        {
            typedef NekVector<NekDouble> Vector;

            double buf1[] = {1.0, 2.0, 3.0, 4.0};
            double buf2[] = {5.0, 6.0, 7.0, 8.0};
            double buf3[] = {9.0, 10.0, 11.0, 12.0};

            Vector v1(4, buf1);
            Vector v2(4, buf2);
            Vector v3(4, buf3);

            Vector result = v1+v2+v3;

            double expectedResultBuf[] = {1.0 + 5.0 + 9.0,
                2.0 + 6.0 + 10.0, 3.0 + 7.0 + 11.0, 4.0 + 8.0 + 12.0};
            Vector expectedResult(4, expectedResultBuf);
            
            BOOST_CHECK_EQUAL(result, expectedResult);
        }

            BOOST_AUTO_TEST_CASE(Test4VectorUnrolled)
        {
            typedef NekVector<NekDouble> Vector;

            double buf1[] = {1.0, 2.0, 3.0, 4.0};
            double buf2[] = {5.0, 6.0, 7.0, 8.0};
            double buf3[] = {9.0, 10.0, 11.0, 12.0};
            double buf4[] = {13.0, 14.0, 15.0, 16.0};

            Vector v1(4, buf1);
            Vector v2(4, buf2);
            Vector v3(4, buf3);
            Vector v4(4, buf4);

            Vector result = v1+v2+v3+v4;

            double expectedResultBuf[] = {1.0 + 5.0 + 9.0+13.0,
                2.0 + 6.0 + 10.0+14.0, 3.0 + 7.0 + 11.0+15.0, 4.0 + 8.0 + 12.0+16.0};
            Vector expectedResult(4, expectedResultBuf);
            
            BOOST_CHECK_EQUAL(result, expectedResult);
        }

        //BOOST_AUTO_TEST_CASE(Test3VectorUnrolledWithSubtraction)
        //{
        //    typedef NekVector<NekDouble> Vector;

        //    double buf1[] = {1.0, 2.0, 3.0, 4.0};
        //    double buf2[] = {5.0, 6.0, 7.0, 8.0};
        //    double buf3[] = {9.0, 10.0, 11.0, 12.0};

        //    Vector v1(4, buf1);
        //    Vector v2(4, buf2);
        //    Vector v3(4, buf3);

        //    Vector result = v1+v2-v3;

        //    double expectedResultBuf[] = {1.0 + 5.0 - 9.0,
        //        2.0 + 6.0 - 10.0, 3.0 + 7.0 - 11.0, 4.0 + 8.0 - 12.0};
        //    Vector expectedResult(4, expectedResultBuf);
        //    
        //    BOOST_CHECK_EQUAL(result, expectedResult);
        //}

        BOOST_AUTO_TEST_CASE(TestVectorUnrolling)
        {
            //typedef expt::Node< NekVector<double> > ConstantNode;
            //BOOST_MPL_ASSERT(( expt::NodeCanUnroll<ConstantNode> ));

            //typedef expt::Node< double > ConstantDoubleNode;
            //BOOST_MPL_ASSERT(( boost::mpl::not_<expt::NodeCanUnroll<ConstantDoubleNode> > ));

            //typedef expt::Node<ConstantNode, expt::AddOp, ConstantNode> AddNode;
            //BOOST_MPL_ASSERT(( expt::NodeCanUnroll<AddNode> ));

            //typedef expt::Node<ConstantNode, expt::SubtractOp, ConstantNode> SubtractNode;
            //BOOST_MPL_ASSERT(( expt::NodeCanUnroll<SubtractNode> ));

            //typedef expt::Node<AddNode, expt::SubtractOp, AddNode> TwoLevelNode;
            //BOOST_MPL_ASSERT(( expt::NodeCanUnroll<TwoLevelNode> ));

            //typedef expt::Node<ConstantNode, expt::MultiplyOp, ConstantNode> MultiplyNode;
            //BOOST_MPL_ASSERT(( boost::mpl::not_<expt::NodeCanUnroll<MultiplyNode> > ));
        }
    }
}
