///////////////////////////////////////////////////////////////////////////////
//
// File: TestBinaryExpressionEvaluatorWithConstantConstantParameters.cpp
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

#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

namespace Nektar
{
    namespace TestTwoParameters
    {
        BOOST_AUTO_TEST_CASE(TestConstantConstantSpecialization1)
        {
//             double lhs_buf[] = {1, 2, 3, 4,
//                                 5, 6, 7, 8,
//                                 9, 10, 11, 12,
//                                 13, 14, 15, 16};
//             double rhs_buf[] = {2, 4, 6, 8,
//                                 10, 12, 14, 16,
//                                 18, 20, 22, 24,
//                                 26, 28, 30, 32};
// 
//             NekMatrix<double> lhs(4, 4, lhs_buf);
//             NekMatrix<double> rhs(4, 4, rhs_buf);
// 
//             Expression<ConstantExpressionPolicy<NekMatrix<double> > > lhs_exp(lhs);
//             Expression<ConstantExpressionPolicy<NekMatrix<double> > > rhs_exp(rhs);
//             NekMatrix<double> result(4,4);
//             Accumulator<NekMatrix<double> > a(result);
// 
//             BinaryExpressionEvaluator<ConstantExpressionPolicy<NekMatrix<double> >,
//                                       ConstantExpressionPolicy<NekMatrix<double> >,
//                                       NekMatrix<double>, AddOp<NekMatrix<double>, NekMatrix<double> >, 
//                                       BinaryNullOp>::Eval(lhs_exp, rhs_exp, a);
// 
//             
//             double expected_result_buf[] = {3, 6, 9, 12,
//                                             15, 18, 21, 24,
//                                             27, 30, 33, 36,
//                                             39, 42, 45, 48};
// 
//             NekMatrix<double> expected_result(4, 4, expected_result_buf);
//             BOOST_CHECK_EQUAL(expected_result, *a);
        }

        BOOST_AUTO_TEST_CASE(TestConstantConstantSpecialization2_AddAdd)
        {
//             typedef ConstantExpressionPolicy<NekMatrix<double> > C;
// 
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
//             Expression<C> middle_exp(middle);
//             Expression<C> rhs_exp(rhs);
//             Accumulator<NekMatrix<double> > a(lhs);
// 
//             // Dereferencing the accumulator indicates that it is initialized.
//             *a;
// 
//             BinaryExpressionEvaluator<ConstantExpressionPolicy<NekMatrix<double> >,
//                                       ConstantExpressionPolicy<NekMatrix<double> >,
//                                       NekMatrix<double>, AddOp, AddOp>::Eval(middle_exp, rhs_exp, a);
// 
//             double expected_result_buf[] = {7, 14, 25, 32,
//                                         39, 46, 57, 64,
//                                         71, 78, 89, 96,
//                                         103, 100, 121, 128};
//             NekMatrix<double> expected_result(4, 4, expected_result_buf);
//             BOOST_CHECK_EQUAL(expected_result, *a);
        }

        BOOST_AUTO_TEST_CASE(TestConstantConstantSpecialization2_AddSubtract)
        {
//             typedef ConstantExpressionPolicy<NekMatrix<double> > C;
// 
//             double lhs_buf[] = {1, 2, 3, 4,
//                                 5, 6, 7, 8,
//                                 9, 10, 11, 12,
//                                 13, 14, 15, 16};
//             double middle_buf[] = {4, 8, 16, 20,
//                                    24, 28, 36, 40,
//                                    44, 48, 56, 60,
//                                    64, 68, 76, 80};
//             double rhs_buf[] = {2, 4, 6, 8,
//                                 10, 12, 14, 16,
//                                 18, 20, 22, 24,
//                                 26, 28, 30, 32};
// 
//             NekMatrix<double> lhs(4, 4, lhs_buf);
//             NekMatrix<double> rhs(4, 4, rhs_buf);
//             NekMatrix<double> middle(4, 4, middle_buf);
//             
//             Expression<C> middle_exp(middle);
//             Expression<C> rhs_exp(rhs);
//             Accumulator<NekMatrix<double> > a(lhs);
// 
//             // Dereferencing the accumulator indicates that it is initialized.
//             *a;
// 
//             BinaryExpressionEvaluator<ConstantExpressionPolicy<NekMatrix<double> >,
//                                       ConstantExpressionPolicy<NekMatrix<double> >,
//                                       NekMatrix<double>, SubtractOp, AddOp>::Eval(middle_exp, rhs_exp, a);
// 
//             double expected_result_buf[] = {3, 6, 13, 16,
//                                         19, 22, 29, 32,
//                                         35, 38, 45, 48,
//                                         51, 54, 61, 64};
//             NekMatrix<double> expected_result(4, 4, expected_result_buf);
//             BOOST_CHECK_EQUAL(expected_result, *a);
        }

        BOOST_AUTO_TEST_CASE(TestConstantConstantSpecialization2_SubtractAdd)
        {
//             typedef ConstantExpressionPolicy<NekMatrix<double> > C;
// 
//             double lhs_buf[] = {1, 2, 3, 4,
//                                 5, 6, 7, 8,
//                                 9, 10, 11, 12,
//                                 13, 14, 15, 16};
//             double middle_buf[] = {4, 8, 16, 20,
//                                    24, 28, 36, 40,
//                                    44, 48, 56, 60,
//                                    64, 68, 76, 80};
//             double rhs_buf[] = {2, 4, 6, 8,
//                                 10, 12, 14, 16,
//                                 18, 20, 22, 24,
//                                 26, 28, 30, 32};
// 
//             NekMatrix<double> lhs(4, 4, lhs_buf);
//             NekMatrix<double> rhs(4, 4, rhs_buf);
//             NekMatrix<double> middle(4, 4, middle_buf);
//             
//             Expression<C> middle_exp(middle);
//             Expression<C> rhs_exp(rhs);
//             Accumulator<NekMatrix<double> > a(lhs);
// 
//             // Dereferencing the accumulator indicates that it is initialized.
//             *a;
// 
//             BinaryExpressionEvaluator<ConstantExpressionPolicy<NekMatrix<double> >,
//                                       ConstantExpressionPolicy<NekMatrix<double> >,
//                                       NekMatrix<double>, AddOp, SubtractOp>::Eval(middle_exp, rhs_exp, a);
// 
//             double expected_result_buf[] = {-5, -10, -19, -24,
//                                             -29, -34, -43, -48,
//                                             -53, -58, -67, -72,
//                                             -77, -82, -91, -96};
//             NekMatrix<double> expected_result(4, 4, expected_result_buf);
//             BOOST_CHECK_EQUAL(expected_result, *a);
        }
    }
}
