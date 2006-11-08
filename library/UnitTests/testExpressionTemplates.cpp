///////////////////////////////////////////////////////////////////////////////
//
// File: testExpressionTemplates.cpp
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#include <UnitTests/testExpressionTemplates.h>

#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>

#include <LibUtilities/LinearAlgebra/NekPoint.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>

#include <iostream>

using std::cout;
using std::endl;

namespace Nektar
{
    namespace UnitTests
    {
        using namespace Nektar;

        void testConstantExpressions()
        {
            using namespace Nektar;
            using namespace Nektar::expt;

            {
                Expression<ConstantExpressionPolicy<int> > e1(7);
                BOOST_CHECK_EQUAL(*e1, 7);
                BOOST_CHECK_EQUAL(e1.GetValue(), 7);

                Expression<ConstantExpressionPolicy<int> > e2(e1);
                BOOST_CHECK_EQUAL(*e2, *e1);
            }

            typedef NekPoint<unsigned int, 3> Point;
            {
                Point p(1,2,3);
                Expression<ConstantExpressionPolicy<Point> > e1(p);

                Point p1(e1);
                BOOST_CHECK_EQUAL(p, p1);

                Point p2 = e1;
                BOOST_CHECK_EQUAL(p, p2);

                Point p3(9, 10, 11);
                BOOST_CHECK(p != p3);
                p3 = e1;
                BOOST_CHECK_EQUAL(p, p3);
            }

            {
                // TODO - Find a way to prevent temporaries (meaning that the parameter to
                // this call is temporary and that could cause problems).
                Nektar::expt::Expression<ConstantExpressionPolicy<Point> > e2(Point(1,2,3));
            }
        }

        void testUnaryExpressions()
        {
            using namespace Nektar;

            NekPoint<double, 3> p(1,2,3);
            NekPoint<double, 3> p1(-(-p));

            BOOST_CHECK_EQUAL(p, p1);

            NekMatrix<double, eFull> m(3,3);
            NekMatrix<double, eFull> m1(-m);
        }

        void testNekMatrixMetadata()
        {
            using namespace Nektar;
            using namespace Nektar::expt;

            // Constant
            NekMatrix<double> m(3,3);
            Expression<ConstantExpressionPolicy<NekMatrix<double> > > m_exp(m);
            BOOST_CHECK_EQUAL(m.GetRows(), m_exp.GetMetadata().Rows);
            BOOST_CHECK_EQUAL(m.GetColumns(), m_exp.GetMetadata().Columns);

            // Unary
            NekMatrix<double> m1(3,3);
            Expression<UnaryExpressionPolicy<Expression<ConstantExpressionPolicy<NekMatrix<double> > >, NegateOp> > m1_exp = 
                Expression<UnaryExpressionPolicy<Expression<ConstantExpressionPolicy<NekMatrix<double> > >, NegateOp> >(
                Expression<ConstantExpressionPolicy<NekMatrix<double> > >(m1));
            BOOST_CHECK_EQUAL(m1.GetRows(), m1_exp.GetMetadata().Rows);
            BOOST_CHECK_EQUAL(m1.GetColumns(), m1_exp.GetMetadata().Columns);

            // Binary
            // TODO
        }

        void testBinaryExpressions()
        {
            using namespace Nektar;
            using namespace Nektar::expt;

            unsigned int buf1[] = {1, 2, 3, 4};
            unsigned int buf2[] = {3, 4, 5, 6};

            // The following crashes because of the temporary.  What to do here?
            //Nektar::ConstantExpression<Nektar::LibUtilities::NekMatrix<unsigned int> > m1(Nektar::LibUtilities::NekMatrix<unsigned int>(2, 2, buf1));
            //Nektar::ConstantExpression<Nektar::LibUtilities::NekMatrix<unsigned int> > m2(Nektar::LibUtilities::NekMatrix<unsigned int>(2, 2, buf2));

            NekMatrix<unsigned int> lhs(2,2,buf1);
            NekMatrix<unsigned int> rhs(2,2,buf2);
            Expression<ConstantExpressionPolicy<Nektar::NekMatrix<unsigned int> > > m1(lhs);
            Expression<ConstantExpressionPolicy<Nektar::NekMatrix<unsigned int> > > m2(rhs);


            Expression<BinaryExpressionPolicy<
                Expression<ConstantExpressionPolicy<NekMatrix<unsigned int> > >,
                Expression<ConstantExpressionPolicy<NekMatrix<unsigned int> > >,
                AddOp > > bexp(m1, m2);

            unsigned int result_buf[] = {4, 6, 8, 10};
            NekMatrix<unsigned int> result(2, 2, result_buf);
            NekMatrix<unsigned int> evaluated_result(bexp);
            BOOST_CHECK_EQUAL(result, evaluated_result);

            {
                // Nested additions.
                NekMatrix<unsigned int> m1(2,2,buf1);
                NekMatrix<unsigned int> m2(2,2,buf2);
                NekMatrix<unsigned int> m3(m1);
                NekMatrix<unsigned int> m4(m2);

                NekMatrix<unsigned int> m5 = m1+m2+m3+m4;

                unsigned int buf3[] = { 8, 12, 16, 20 };
                NekMatrix<unsigned int> nested_add_result(2,2,buf3);

                BOOST_CHECK_EQUAL(nested_add_result, m5);
            }
        }
        
        void testNekMatrixMultiplication()
        {
            double m1_buf[] = {-85, -55, -37, -35, 97, 50, 79, 56, 49};
            double m2_buf[] = {63, 57, -59, 45, -8, -93, 92, 43, -62};
            double m3_buf[] = {77, 66, 54, -5, 99, -61, -50, -12, -18};

            NekMatrix<double> m1(3, 3, m1_buf);
            NekMatrix<double> m2(3, 3, m2_buf);
            NekMatrix<double> m3(3, 3, m3_buf);
            
            NekMatrix<double> result = m1*m2*m3;
            
            double result_buf[] = {-1456238, -1484136, -464512, 1026425, 505353, 583929, 1538925, 1557252, 504714};
            NekMatrix<double> expectedResult(3, 3, result_buf);
            
            double epsilon = 1e-11;
            for(unsigned int i = 0; i < 3; ++i)
            {
                for(unsigned int j = 0; j < 3; ++j)
                {
                    BOOST_CHECK_CLOSE(result(i,j), expectedResult(i,j), epsilon);
                }
            }
            
            NekMatrix<double> result1 = m1*m2*m3*m1*m2*m3*m1*m2*m3;
            double result_buf1[] = {223791291531519928.0, -139146145309301688.0, 241968403742002232.0, -81497861322837100.0, 109613922100149537.0, -116433655760219405.0, -233781330982473300.0, 141216567102193860.0, -250757429804037708.0};
            NekMatrix<double> expectedResult1(3,3,result_buf1);
            for(unsigned int i = 0; i < 3; ++i)
            {
                for(unsigned int j = 0; j < 3; ++j)
                {
                    BOOST_CHECK_CLOSE(result1(i,j), expectedResult1(i,j), epsilon);
                }
            }
        }
        
        void testNekMatrixSomewhatComplicatedExpression()
        {
            {
                double m1_buf[] = {-85, -55, -37, -35, 97, 50, 79, 56, 49};
                double m2_buf[] = {63, 57, -59, 45, -8, -93, 92, 43, -62};
                double m3_buf[] = {77, 66, 54, -5, 99, -61, -50, -12, -18};
    
                NekMatrix<double> m1(3, 3, m1_buf);
                NekMatrix<double> m2(3, 3, m2_buf);
                NekMatrix<double> m3(3, 3, m3_buf);
                
                NekMatrix<double> result = (m1*m2) + m3;
                
                double result_buf[] = {-11157, -5930, 12478, 6755, -522, -10117, 11955, 6150, -12925};
                NekMatrix<double> expectedResult(3,3,result_buf);
                double epsilon = 1e-11;
                for(unsigned int i = 0; i < 3; ++i)
                {
                    for(unsigned int j = 0; j < 3; ++j)
                    {
                        BOOST_CHECK_CLOSE(result(i,j), expectedResult(i,j), epsilon);
                    }
                }
                
                NekMatrix<double> result1 = m3 + (m1*m2);
                for(unsigned int i = 0; i < 3; ++i)
                {
                    for(unsigned int j = 0; j < 3; ++j)
                    {
                        BOOST_CHECK_CLOSE(result1(i,j), expectedResult(i,j), epsilon);
                    }
                }
            }
        }
        
        void testNekMatrixComplicatedExpression()
        {
            {
                double m1_buf[] = {-85, -55, -37, -35, 97, 50, 79, 56, 49};
                double m2_buf[] = {63, 57, -59, 45, -8, -93, 92, 43, -62};
                double m3_buf[] = {31, -26, -62, 1, -47, -91, -47, -61, 41};
    
                NekMatrix<double> m1(3, 3, m1_buf);
                NekMatrix<double> m2(3, 3, m2_buf);
                NekMatrix<double> m3(3, 3, m3_buf);
                
                NekMatrix<double> result = (((m3-m1)*m3) * (m1-m2)) + m1*m2*m3;
                
                double result_buf[] = {-1279130, -1162366, 243990, -1663904, 1197403, 2021293, 547959, 1802365, 1422677};
                NekMatrix<double> expectedResult(3,3,result_buf);
                
                double epsilon = 1e-11;
                for(unsigned int i = 0; i < 3; ++i)
                {
                    for(unsigned int j = 0; j < 3; ++j)
                    {
                        BOOST_CHECK_CLOSE(result(i,j), expectedResult(i,j), epsilon);
                    }
                }
            }
        }
        
     }
}

/**
    $Log: testExpressionTemplates.cpp,v $
    Revision 1.8  2006/11/06 17:10:04  bnelson
    *** empty log message ***

    Revision 1.7  2006/09/30 15:38:29  bnelson
    no message

    Revision 1.6  2006/09/15 02:01:16  bnelson
    no message

    Revision 1.5  2006/09/11 03:28:41  bnelson
    no message

    Revision 1.4  2006/08/28 02:40:51  bnelson
    *** empty log message ***

    Revision 1.3  2006/08/27 02:14:09  bnelson
    Added support for negating an expression.

    Revision 1.2  2006/08/25 01:37:34  bnelson
    no message

    Revision 1.1  2006/08/25 01:36:25  bnelson
    no message


**/
