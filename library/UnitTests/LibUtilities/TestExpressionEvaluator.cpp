///////////////////////////////////////////////////////////////////////////////
//
// File: TestExpressionEvaluator.cpp
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

#include <LibUtilities/Interpreter/ExpressionEvaluator.h>
#include <UnitTests/util.h>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

#include <math.h>

namespace Nektar
{
    namespace ExpressionTests
    {
        BOOST_AUTO_TEST_CASE(TestIdentityExpression)
        {
            //LibUtilities::ExpressionEvaluator eval;
            //eval.DefineFunction("x", "x");
            //
            //double x1 = 7.932;
            //double r1 = eval.Evaluate(x1);
            //BOOST_CHECK_EQUAL(x1, r1);
        }
        
        BOOST_AUTO_TEST_CASE(TestExpressionMultiplication)
        {
            //LibUtilities::ExpressionEvaluator eval;
            //eval.DefineFunction("x", "x*x");
            //
            //double r1 = eval.Evaluate(2.0);
            //BOOST_CHECK_EQUAL(r1, 2.0*2.0);
            //
            //double p1 = 82.45;
            //double r2 = eval.Evaluate(p1);
            //BOOST_CHECK_EQUAL(r2, p1*p1);
        }

        BOOST_AUTO_TEST_CASE(TestSingleParameterExpression)
        {
            //LibUtilities::ExpressionEvaluator eval;
            //eval.DefineFunction("x", "sin(x)");
            //
            //double x1 = .23612;
            //double r1 = eval.Evaluate(x1);
            //BOOST_CHECK_EQUAL(r1, sin(x1));
            
        }
        
        BOOST_AUTO_TEST_CASE(TestSingleParameterExpressionWithTwoVariables)
        {
            //LibUtilities::ExpressionEvaluator eval;
            //eval.DefineFunction("x y", "sin(x*y)");
            //
            //double x1 = .23612;
            //double y1 = .19;
            //double r1 = eval.Evaluate(x1, y1);
            //BOOST_CHECK_EQUAL(r1, sin(x1*y1));
            
        }
        
        BOOST_AUTO_TEST_CASE(TestTwoParameterExpression)
        {
            //LibUtilities::ExpressionEvaluator eval;
            //eval.DefineFunction("x y", "x*y");
            //
            //double x1 = 6.2;
            //double y1 = 5.6;
            //double r1 = eval.Evaluate(x1, y1);
            //
            //BOOST_CHECK_EQUAL(r1, x1*y1);
        }
        
        BOOST_AUTO_TEST_CASE(TestFMod)
        {
            //LibUtilities::ExpressionEvaluator eval;
            //eval.DefineFunction("x y", "fmod(x, y)");
            //
            //double x1 = 5.6;
            //double y1 = 4.5;
            //double r1 = eval.Evaluate(x1, y1);
            //BOOST_CHECK_EQUAL(r1, fmod(x1, y1));
        }
        
//        BOOST_AUTO_TEST_CASE(TestTwoParameterFunctions)
//        {
//            LibUtilities::ExpressionEvaluator eval;
//            eval.DefineFunction("x y", "pow(x, y)");
//            
//            double x1 = 6.2;
//            double y1 = 5.6;
//            double r1 = eval.Evaluate(x1, y1);
//            
//            BOOST_CHECK_EQUAL(r1, pow(x1, y1));
//        }

    }
}

