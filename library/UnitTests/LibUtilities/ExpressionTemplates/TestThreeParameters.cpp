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
        
        GENERATE_ADDITION_OPERATOR(Foo, 0, Bar, 0);
        
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
        
        GENERATE_ADDITION_OPERATOR(Bar, 0, Bar, 0);
        
        class ConstantBinarySpecialization0_Parameter
        {
            public:
                /// \todo Doc note - all classes in expression templates need a default constructor.
                ConstantBinarySpecialization0_Parameter() : m_value(0) {}
                explicit ConstantBinarySpecialization0_Parameter(int v) : m_value(v) {}
                ConstantBinarySpecialization0_Parameter(const ConstantBinarySpecialization0_Parameter& rhs) : m_value(rhs.m_value) {}
                ConstantBinarySpecialization0_Parameter operator=(const  ConstantBinarySpecialization0_Parameter& rhs)
                {
                    m_value = rhs.m_value;
                    return *this;
                }
                
                int m_value;
        };
        
        void NekMultiply(ConstantBinarySpecialization0_Parameter& result,
                         const ConstantBinarySpecialization0_Parameter& lhs,
                         const ConstantBinarySpecialization0_Parameter& rhs)
        {
            result.m_value = lhs.m_value * rhs.m_value;
        }
        
        ConstantBinarySpecialization0_Parameter NekMultiply(const ConstantBinarySpecialization0_Parameter& lhs,
                         const ConstantBinarySpecialization0_Parameter& rhs)
        {
            ConstantBinarySpecialization0_Parameter result(0);
            NekMultiply(result, lhs, rhs);
            return result;
        }
                         
        class ConstantBinarySpecialization0_Result
        {
            public:
                ConstantBinarySpecialization0_Result() : m_value() {}
                ConstantBinarySpecialization0_Result(int v) : m_value(v) {}
                ConstantBinarySpecialization0_Result(const ConstantBinarySpecialization0_Result& rhs) : m_value(rhs.m_value) {}
                ConstantBinarySpecialization0_Result operator=(const  ConstantBinarySpecialization0_Result& rhs)
                {
                    m_value = rhs.m_value;
                    return *this;
                }
                
                int m_value;
  
        };
  
        void NekAdd(ConstantBinarySpecialization0_Result& result,
                    const ConstantBinarySpecialization0_Parameter& lhs,
                    const ConstantBinarySpecialization0_Parameter& rhs)
        {
            result.m_value = lhs.m_value + rhs.m_value;
        }
        
        ConstantBinarySpecialization0_Result NekAdd(const ConstantBinarySpecialization0_Parameter& lhs,
                         const ConstantBinarySpecialization0_Parameter& rhs)
        {
            ConstantBinarySpecialization0_Result result(0);
            NekAdd(result, lhs, rhs);
            return result;
        }
        
        GENERATE_ADDITION_OPERATOR(ConstantBinarySpecialization0_Parameter, 0, ConstantBinarySpecialization0_Parameter, 0);
        GENERATE_MULTIPLICATION_OPERATOR(ConstantBinarySpecialization0_Parameter, 0, ConstantBinarySpecialization0_Parameter, 0);
        
        BOOST_AUTO_TEST_CASE(TestConstantBinarySpecialization0_TwoTemporariesRequired)
        {
            // This test case requires a left side and a right side that 
            // don't evaluate to the result type individually, but do add to the
            // result type.
            ConstantBinarySpecialization0_Parameter p0(8);
            ConstantBinarySpecialization0_Parameter p1(2);
            ConstantBinarySpecialization0_Parameter p2(20);
            
            typedef BinaryExpressionPolicy<ConstantExpressionPolicy<ConstantBinarySpecialization0_Parameter>,
                                           MultiplyOp,
                                           ConstantExpressionPolicy<ConstantBinarySpecialization0_Parameter> > RhsExpressionType;
            typedef ConstantExpressionPolicy<ConstantBinarySpecialization0_Parameter> LhsExpressionType;
            
            //BOOST_STATIC_ASSERT(( BinaryExpressionEvaluator<LhsExpressionType, RhsExpressionType, ConstantBinarySpecialization0_Result, AddOp, BinaryNullOp>::ClassNum == 1 ));
             
//            Expression<BinaryExpressionPolicy<LhsExpressionType, AddOp, RhsExpressionType> > exp =
//                p0 + (p1*p2);
//
//            ConstantBinarySpecialization0_Result r;
//            Assign(r, exp);
//            BOOST_CHECK_EQUAL(r.m_value, 8+2*20);
        }
//          Expression<BinaryExpressionPolicy<ConstantExpressionPolicy<Bar>, ConstantExpressionPolicy<Bar>,AddOp > > 
//          operator+(const Bar& lhs, const Bar& rhs) 
//          { 
//              return CreateBinaryExpression<AddOp>(lhs, rhs); 
//          }
//          
//          Expression<BinaryExpressionPolicy<ConstantExpressionPolicy<Foo>, ConstantExpressionPolicy<Bar>,AddOp > > 
//          operator+(const Foo& lhs, const Bar& rhs) 
//          { 
//              return CreateBinaryExpression<AddOp>(lhs, rhs); 
//          }
        
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
                        BinaryExpressionPolicy<ConstantExpressionPolicy<NekMatrix<double> >, AddOp,
                                               ConstantExpressionPolicy<NekMatrix<double> > >, AddOp,
                                               ConstantExpressionPolicy<NekMatrix<double> >
                                               > > expr = (lhs + middle) + rhs;
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
//             Foo lhs(1);
//             Bar middle(2);
//             Bar rhs(3);
//             
//             Foo result;
//             
//             ((lhs + middle) + rhs).Evaluate(result);
//             
//             BOOST_CHECK_EQUAL(result.Value, 6);
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
                        BinaryExpressionPolicy<ConstantExpressionPolicy<NekMatrix<double> >, AddOp,
                                               ConstantExpressionPolicy<NekMatrix<double> > >, AddOp,
                                               ConstantExpressionPolicy<NekMatrix<double> >
                                               > > expr = (lhs + middle) + rhs;
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
                        ConstantExpressionPolicy<NekMatrix<double> >, AddOp,
                        BinaryExpressionPolicy<ConstantExpressionPolicy<NekMatrix<double> >, AddOp, 
                                               ConstantExpressionPolicy<NekMatrix<double> > >
                                               > > expr = lhs + (middle + rhs);
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
                        BinaryExpressionPolicy<ConstantExpressionPolicy<NekMatrix<double> >, AddOp,
                                               ConstantExpressionPolicy<NekMatrix<double> > >, SubtractOp,
                                               ConstantExpressionPolicy<NekMatrix<double> >
                                               > > expr = (lhs + middle) - rhs;
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

             Expression<BinaryExpressionPolicy<
                        ConstantExpressionPolicy<NekMatrix<double> >, SubtractOp,
                        BinaryExpressionPolicy<ConstantExpressionPolicy<NekMatrix<double> >, AddOp,
                                               ConstantExpressionPolicy<NekMatrix<double> > >
                                               > > expr = lhs - (middle + rhs);
             NekMatrix<double> result(lhs - (middle+rhs));
 
             double expected_result_buf[] = {-4, -8, -12, -16,
                                             -20, -24, -28, -32,
                                             -36, -40, -44, -48,
                                             -52, -56, -60, -64};
             NekMatrix<double> expected_result(4, 4, expected_result_buf);
             BOOST_CHECK_EQUAL(expected_result, result);
        }
        
        BOOST_AUTO_TEST_CASE(Test3ParameterExpressionWithoutParens)
        {
            double buf1[] = {1, 2, 3, 4};
            double buf2[] = {5, 6, 7, 9};
            double buf3[] = {10, 11, 12, 13};
            
            NekMatrix<double> m1(2, 2, buf1);
            NekMatrix<double> m2(2, 2, buf2);
            NekMatrix<double> m3(2, 2, buf3);
            
            double expected_result_buf[] = {buf1[0] - buf2[0] + buf3[0],
                buf1[1] - buf2[1] + buf3[1],
                buf1[2] - buf2[2] + buf3[2],
                buf1[3] - buf2[3] + buf3[3]};
            NekMatrix<double> expected_result(2, 2, expected_result_buf);
            NekMatrix<double> result = m1-m2+m3;
            BOOST_CHECK_EQUAL(expected_result, result);
        }
        
        BOOST_AUTO_TEST_CASE(Test4ParameterExpressionWithoutParens)
        {
            double buf1[] = {1, 2, 3, 4};
            double buf2[] = {5, 6, 7, 9};
            double buf3[] = {10, 11, 12, 13};
            double buf4[] = {14, 15, 16, 17};
            
            NekMatrix<double> m1(2, 2, buf1);
            NekMatrix<double> m2(2, 2, buf2);
            NekMatrix<double> m3(2, 2, buf3);
            NekMatrix<double> m4(2, 2, buf4);
            
            double expected_result_buf[] = {-112, -142, -132, -168};
            NekMatrix<double> expected_result(2, 2, expected_result_buf);
            NekMatrix<double> result = m1-m2*m3+m4;
            BOOST_CHECK_EQUAL(expected_result, result);
        }
    }
}
