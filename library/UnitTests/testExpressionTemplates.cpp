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

#include <LibUtilities/ExpressionTemplates/ConstantExpression.hpp>
#include <LibUtilities/LinearAlgebra/NekPoint.hpp>

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
            {
                ConstantExpression<int> e1(7);
                BOOST_CHECK_EQUAL(*e1, 7);
                BOOST_CHECK_EQUAL(e1.GetValue(), 7);

                ConstantExpression<int> e2(e1);
                BOOST_CHECK_EQUAL(*e2, *e1);
            }

            using namespace Nektar::LibUtilities;

            typedef NekPoint<unsigned int, 3> Point;
            {
                Point p(1,2,3);
                Nektar::ConstantExpression<Point> e1(p);

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
                Nektar::ConstantExpression<Point> e2(Point(1,2,3));    
            }
        }

        void testUnaryExpressions()
        {
            using namespace Nektar::LibUtilities;

            NekPoint<double, 3> p(1,2,3);
            //NekPoint<double, 3> p1(-(-p));

            //BOOST_CHECK_EQUAL(p, p1);
        }

     }
}

/**
    $Log: testExpressionTemplates.cpp,v $
    Revision 1.1  2006/08/25 01:36:25  bnelson
    no message


**/
