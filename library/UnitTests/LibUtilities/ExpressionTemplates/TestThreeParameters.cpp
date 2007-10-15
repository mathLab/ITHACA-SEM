///////////////////////////////////////////////////////////////////////////////
//
// File: TestThreeParameters.cpp
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
    namespace TestThreeParametersWithNekMatrix
    {
        class Foo
        {
            public:
                Foo() : Value(0) {}
                Foo(int v) : Value(v) {}
                int Value;
        };
        
        class Bar
        {
            public:
                Bar() : Value(0) {}
                Bar(int v) : Value(v) {}
                int Value;
        };
        
        void NekAdd(Bar& result, const Foo& lhs, const Bar& rhs)
        {
            result.Value = lhs.Value + rhs.Value;
        }
        
        Bar NekAdd(const Foo& lhs, const Bar& rhs)
        {
            Bar result;
            NekAdd(result, lhs, rhs);
            return result;
        }
        
        void NekAdd(Foo& result, const Bar& lhs, const Bar& rhs)
        {
            result.Value = lhs.Value + rhs.Value;
        }
        
        Foo NekAdd(const Bar& lhs, const Bar& rhs)
        {
            Foo result;
            NekAdd(result, lhs, rhs);
            return result;
        }
        
            
         Expression<BinaryExpressionPolicy<ConstantExpressionPolicy<Bar>, ConstantExpressionPolicy<Bar>,AddOp > > 
         operator+(const Bar& lhs, const Bar& rhs) 
         { 
             return CreateBinaryExpression<AddOp>(lhs, rhs); 
         }
         
         Expression<BinaryExpressionPolicy<ConstantExpressionPolicy<Foo>, ConstantExpressionPolicy<Bar>,AddOp > > 
         operator+(const Foo& lhs, const Bar& rhs) 
         { 
             return CreateBinaryExpression<AddOp>(lhs, rhs); 
         }
        
        BOOST_AUTO_TEST_CASE(TestLhsBinaryRhsConstantSpecialization1)
        {
            double lhs_buf[] = {1, 2, 3, 4,
                                5, 6, 7, 8,
                                9, 10, 11, 12,
                                13, 14, 15, 16};
            double middle_buf[] = {2, 4, 6, 8,
                                10, 12, 14, 16,
                                18, 20, 22, 24,
                                26, 28, 30, 32};

            double rhs_buf[] = { 3, 6, 9, 12,
                                 15, 18, 21, 24,
                                 27, 30, 33, 36,
                                 39, 42, 45, 48};

            NekMatrix<double> lhs(4, 4, lhs_buf);
            NekMatrix<double> rhs(4, 4, rhs_buf);
            NekMatrix<double> middle(4, 4, middle_buf);

            Expression<BinaryExpressionPolicy<
                       BinaryExpressionPolicy<ConstantExpressionPolicy<NekMatrix<double> >, 
                                              ConstantExpressionPolicy<NekMatrix<double> >, AddOp>,
                                              ConstantExpressionPolicy<NekMatrix<double> >,
                                              AddOp> > expr = (lhs + middle) + rhs;
            NekMatrix<double> result(expr);

            double expected_result_buf[] = {6, 12, 18, 24,
                                            30, 36, 42, 48,
                                            54, 60, 66, 72,
                                            78, 84, 90, 96};
            NekMatrix<double> expected_result(4, 4, expected_result_buf);
            BOOST_CHECK_EQUAL(expected_result, result);
        }
        
        
        BOOST_AUTO_TEST_CASE(TestLhsBinaryRhsConstantSpecialization2)
        {
            Foo lhs(1);
            Bar middle(2);
            Bar rhs(3);
            
            Foo result;
            
            ((lhs + middle) + rhs).Apply(result);
            
            BOOST_CHECK_EQUAL(result.Value, 6);
        }
        
        // Tests the creation of a single binary expression from the + 
        // operator.
        BOOST_AUTO_TEST_CASE(TestAddingThreeMatrices_LhsIsBinary_TopOperationIsCommutative)
        {
            double lhs_buf[] = {1, 2, 3, 4,
                                5, 6, 7, 8,
                                9, 10, 11, 12,
                                13, 14, 15, 16};
            double middle_buf[] = {2, 4, 6, 8,
                                10, 12, 14, 16,
                                18, 20, 22, 24,
                                26, 28, 30, 32};

            double rhs_buf[] = { 3, 6, 9, 12,
                                 15, 18, 21, 24,
                                 27, 30, 33, 36,
                                 39, 42, 45, 48};

            NekMatrix<double> lhs(4, 4, lhs_buf);
            NekMatrix<double> rhs(4, 4, rhs_buf);
            NekMatrix<double> middle(4, 4, middle_buf);

            Expression<BinaryExpressionPolicy<
                       BinaryExpressionPolicy<ConstantExpressionPolicy<NekMatrix<double> >, 
                                              ConstantExpressionPolicy<NekMatrix<double> >, AddOp>,
                                              ConstantExpressionPolicy<NekMatrix<double> >,
                                              AddOp> > expr = (lhs + middle) + rhs;
            NekMatrix<double> result(expr);

            double expected_result_buf[] = {6, 12, 18, 24,
                                            30, 36, 42, 48,
                                            54, 60, 66, 72,
                                            78, 84, 90, 96};
            NekMatrix<double> expected_result(4, 4, expected_result_buf);
            BOOST_CHECK_EQUAL(expected_result, result);
        }

        BOOST_AUTO_TEST_CASE(TestAddingThreeMatrices_RhsIsBinary_TopOperationIsCommutative)
        {
            double lhs_buf[] = {1, 2, 3, 4,
                                5, 6, 7, 8,
                                9, 10, 11, 12,
                                13, 14, 15, 16};
            double middle_buf[] = {2, 4, 6, 8,
                                10, 12, 14, 16,
                                18, 20, 22, 24,
                                26, 28, 30, 32};

            double rhs_buf[] = { 3, 6, 9, 12,
                                 15, 18, 21, 24,
                                 27, 30, 33, 36,
                                 39, 42, 45, 48};

            NekMatrix<double> lhs(4, 4, lhs_buf);
            NekMatrix<double> rhs(4, 4, rhs_buf);
            NekMatrix<double> middle(4, 4, middle_buf);

            Expression<BinaryExpressionPolicy<
                       ConstantExpressionPolicy<NekMatrix<double> >,
                       BinaryExpressionPolicy<ConstantExpressionPolicy<NekMatrix<double> >, 
                                              ConstantExpressionPolicy<NekMatrix<double> >, AddOp>,
                                              AddOp> > expr = lhs + (middle + rhs);
            NekMatrix<double> result(expr);

            double expected_result_buf[] = {6, 12, 18, 24,
                                            30, 36, 42, 48,
                                            54, 60, 66, 72,
                                            78, 84, 90, 96};
            NekMatrix<double> expected_result(4, 4, expected_result_buf);
            BOOST_CHECK_EQUAL(expected_result, result);
        }

        BOOST_AUTO_TEST_CASE(TestAddingThreeMatrices_LhsIsBinary_TopOperationIsNotCommutative)
        {
            double lhs_buf[] = {1, 2, 3, 4,
                                5, 6, 7, 8,
                                9, 10, 11, 12,
                                13, 14, 15, 16};
            double middle_buf[] = {2, 4, 6, 8,
                                10, 12, 14, 16,
                                18, 20, 22, 24,
                                26, 28, 30, 32};

            double rhs_buf[] = { 3, 6, 9, 12,
                                 15, 18, 21, 24,
                                 27, 30, 33, 36,
                                 39, 42, 45, 48};

            NekMatrix<double> lhs(4, 4, lhs_buf);
            NekMatrix<double> rhs(4, 4, rhs_buf);
            NekMatrix<double> middle(4, 4, middle_buf);

            Expression<BinaryExpressionPolicy<
                       BinaryExpressionPolicy<ConstantExpressionPolicy<NekMatrix<double> >, 
                                              ConstantExpressionPolicy<NekMatrix<double> >, AddOp>,
                                              ConstantExpressionPolicy<NekMatrix<double> >,
                                              SubtractOp> > expr = (lhs + middle) - rhs;
            NekMatrix<double> result(expr);

            double expected_result_buf[] = {0, 0, 0, 0,
                                            0, 0, 0, 0,
                                            0, 0, 0, 0,
                                            0, 0, 0, 0};
            NekMatrix<double> expected_result(4, 4, expected_result_buf);
            BOOST_CHECK_EQUAL(expected_result, result);
        }

        BOOST_AUTO_TEST_CASE(TestAddingThreeMatrices_RhsIsBinary_TopOperationIsNotCommutative)
        {
            double lhs_buf[] = {1, 2, 3, 4,
                                5, 6, 7, 8,
                                9, 10, 11, 12,
                                13, 14, 15, 16};
            double middle_buf[] = {2, 4, 6, 8,
                                10, 12, 14, 16,
                                18, 20, 22, 24,
                                26, 28, 30, 32};

            double rhs_buf[] = { 3, 6, 9, 12,
                                 15, 18, 21, 24,
                                 27, 30, 33, 36,
                                 39, 42, 45, 48};

            NekMatrix<double> lhs(4, 4, lhs_buf);
            NekMatrix<double> rhs(4, 4, rhs_buf);
            NekMatrix<double> middle(4, 4, middle_buf);

//             Expression<BinaryExpressionPolicy<
//                        ConstantExpressionPolicy<NekMatrix<double> >,
//                        BinaryExpressionPolicy<ConstantExpressionPolicy<NekMatrix<double> >, 
//                                               ConstantExpressionPolicy<NekMatrix<double> >, AddOp>,
//                                               SubtractOp> > expr = lhs - (middle + rhs);
//             NekMatrix<double> result(lhs - (middle+rhs));
// 
//             double expected_result_buf[] = {-4, -8, -12, -16,
//                                             -20, -24, -28, -32,
//                                             -36, -40, -44, -48,
//                                             -52, -56, -60, -64};
//             NekMatrix<double> expected_result(4, 4, expected_result_buf);
//             BOOST_CHECK_EQUAL(expected_result, result);
        }
    }
}
