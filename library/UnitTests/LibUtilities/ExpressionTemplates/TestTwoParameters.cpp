///////////////////////////////////////////////////////////////////////////////
//
// File: TestTwoParameters.cpp
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

#define BOOST_TEST_MODULE ExpressionTemplateUnitTests test

#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

namespace Nektar
{
    namespace TestTwoParameters
    {
        // Tests 1op Specialization #1
        BOOST_AUTO_TEST_CASE(TestAddingTwoMatrices)
        {
            double lhs_buf[] = {1, 2, 3, 4,
                                5, 6, 7, 8,
                                9, 10, 11, 12,
                                13, 14, 15, 16};
            double rhs_buf[] = {2, 4, 6, 8,
                                10, 12, 14, 16,
                                18, 20, 22, 24,
                                26, 28, 30, 32};

            NekMatrix<double> lhs(4, 4, lhs_buf);
            NekMatrix<double> rhs(4, 4, rhs_buf);
            Expression<BinaryExpressionPolicy<ConstantExpressionPolicy<NekMatrix<double> >, 
                                              ConstantExpressionPolicy<NekMatrix<double> >,
                                              AddOp> > expr = lhs + rhs;
            NekMatrix<double> result(expr);
        }

        BOOST_AUTO_TEST_CASE(Test1OpSpecialization2_WithNoOperatorChange)
        {
//             double lhs_buf[] = {1, 2, 3, 4,
//                                 5, 6, 7, 8,
//                                 9, 10, 11, 12,
//                                 13, 14, 15, 16};
//             double middle_buf[] = {4, 8, 16, 20,
//                                    24, 28, 36, 40,
//                                    44, 48, 56, 60,
//                                    64, 58, 76, 80};
//             double rhs_buf[] = {2, 4, 6, 8,
//                                 10, 12, 14, 16,
//                                 18, 20, 22, 24,
//                                 26, 28, 30, 32};
// 
//             NekMatrix<double> lhs(4, 4, lhs_buf);
//             NekMatrix<double> rhs(4, 4, rhs_buf);
//             NekMatrix<double> middle(4, 4, middle_buf);
//             
//             NekMatrix<double> result(lhs + (middle+rhs));
// 
//             double expected_result_buf[] = {7, 14, 25, 32,
//                                         39, 46, 57, 64,
//                                         71, 78, 89, 96,
//                                         103, 100, 121, 128};
//             NekMatrix<double> expected_result(4, 4, expected_result_buf);
//             BOOST_CHECK_EQUAL(expected_result, result);
        }

        // Tests aliasing.
        BOOST_AUTO_TEST_CASE(TestMatrixAliasing)
        {
        }
    }
}
