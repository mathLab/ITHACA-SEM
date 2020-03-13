///////////////////////////////////////////////////////////////////////////////
//
// File: TestAliasing.cpp
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

#ifndef NEKTAR_USE_EXPRESSION_TEMPLATES
#define NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //NEKTAR_USE_EXPRESSION_TEMPLATES

#include <UnitTests/LibUtilities/ExpressionTemplates/CountedObjectExpression.h>
#include <UnitTests/LibUtilities/ExpressionTemplates/ExpressionTemplateObjects.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

using namespace expt;

namespace Nektar
{
    namespace UnitTests
    {
        BOOST_AUTO_TEST_CASE(TestObjectAliasing)
        {
            typedef Node<TestObjectA> ConstantNode;
            typedef Node<ConstantNode, AddOp, ConstantNode> NodeType;

            TestObjectA obj1(2.3);
            TestObjectA obj2(3.4);
            TestObjectA obj3(4.5);

            NodeType expr1 = obj1+obj2;
            NodeType expr2 = obj1+obj3;

            BOOST_CHECK(ExpressionEvaluator::ContainsAlias(expr1, obj1));
            BOOST_CHECK(ExpressionEvaluator::ContainsAlias(expr1, obj2));
            BOOST_CHECK(!ExpressionEvaluator::ContainsAlias(expr1, obj3));

            BOOST_CHECK(ExpressionEvaluator::ContainsAlias(expr2, obj1));
            BOOST_CHECK(!ExpressionEvaluator::ContainsAlias(expr2, obj2));
            BOOST_CHECK(ExpressionEvaluator::ContainsAlias(expr2, obj3));

            
        }

        BOOST_AUTO_TEST_CASE(TestAliasing)
        {
            double lhs_buf[] = {1, 2, 3, 4,
                                5, 6, 7, 8,
                                9, 10, 11, 12,
                                13, 14, 15, 16};
            
            NekMatrix<double> lhs(4, 4, lhs_buf);
            BOOST_AUTO(expression, lhs+lhs+lhs);
            BOOST_CHECK(ExpressionEvaluator::ContainsAlias(expression, lhs));

            lhs = lhs+lhs+lhs;
            
            double expected_result_buf[] = {3, 6, 9, 12,
                                            15, 18, 21, 24,
                                            27, 30, 33, 36,
                                            39, 42, 45, 48};
                                            
            NekMatrix<double> expected_result(4, 4, expected_result_buf);
            BOOST_CHECK_EQUAL(expected_result, lhs);
        }
    }
}
