///////////////////////////////////////////////////////////////////////////////
//
// File: TestStandardMatrix.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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

#define BOOST_TEST_MODULE LinearAlgebraUnitTests test
#include <boost/test/unit_test.hpp>

#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <functional>
#include <UnitTests/LibUtilities/LinearAlgebra/TestCombinationRunner.h>

namespace Nektar
{
    namespace StandardMatrixUnitTests
    {
        BOOST_AUTO_TEST_CASE(TestIterators)
        {
            double buf1[] = { 1, 4,
                              2, 5,
                              3, 6};

            NekMatrix<double> m1(2, 3, buf1);
            const NekMatrix<double>& m2 = m1;

            NekMatrix<double>::iterator m1_iter = m1.begin();
            NekMatrix<double>::const_iterator m2_iter = m2.begin();

            BOOST_CHECK_EQUAL(1, *m1_iter);
            BOOST_CHECK_EQUAL(1, *m2_iter);
            BOOST_CHECK_EQUAL(1, *m1_iter++);
            BOOST_CHECK_EQUAL(1, *m2_iter++);
            BOOST_CHECK_EQUAL(4, *m1_iter++);
            BOOST_CHECK_EQUAL(4, *m2_iter++);
            BOOST_CHECK_EQUAL(2, *m1_iter++);
            BOOST_CHECK_EQUAL(2, *m2_iter++);
            BOOST_CHECK_EQUAL(5, *m1_iter++);
            BOOST_CHECK_EQUAL(5, *m2_iter++);
            BOOST_CHECK_EQUAL(3, *m1_iter++);
            BOOST_CHECK_EQUAL(3, *m2_iter++);
            BOOST_CHECK_EQUAL(6, *m1_iter++);
            BOOST_CHECK_EQUAL(6, *m2_iter++);
            BOOST_CHECK(m1_iter == m1.end());
            BOOST_CHECK(m2_iter == m2.end());

            m1.Transpose();
            m1_iter = m1.begin();
            m2_iter = m2.begin();

            BOOST_CHECK_EQUAL(1, *m1_iter);
            BOOST_CHECK_EQUAL(1, *m2_iter);
            BOOST_CHECK_EQUAL(2, *(++m1_iter));
            BOOST_CHECK_EQUAL(2, *(++m2_iter));
            BOOST_CHECK_EQUAL(3, *(++m1_iter));
            BOOST_CHECK_EQUAL(3, *(++m2_iter));
            BOOST_CHECK_EQUAL(4, *(++m1_iter));
            BOOST_CHECK_EQUAL(4, *(++m2_iter));
            BOOST_CHECK_EQUAL(5, *(++m1_iter));
            BOOST_CHECK_EQUAL(5, *(++m2_iter));
            BOOST_CHECK_EQUAL(6, *(++m1_iter));
            BOOST_CHECK_EQUAL(6, *(++m2_iter));
            BOOST_CHECK(m1_iter != m1.end());
            BOOST_CHECK(m2_iter != m2.end());
            ++m1_iter;
            ++m2_iter;
            BOOST_CHECK(m1_iter == m1.end());
            BOOST_CHECK(m2_iter == m2.end());
        }
        
        BOOST_AUTO_TEST_CASE(TestWrappedCopyConstructor)
        {
            double buf1[] = { 1, 4,
                              2, 5,
                              3, 6};

            NekMatrix<double> m1(2, 3, buf1);
            NekMatrix<double> m2 = NekMatrix<double>::CreateWrapper(m1);
            
            BOOST_CHECK_EQUAL(1.0, m1(0,0));
            BOOST_CHECK_EQUAL(1.0, m2(0,0));
            
            m1(0,0) = -4.5;
            
            BOOST_CHECK_EQUAL(-4.5, m1(0,0));
            BOOST_CHECK_EQUAL(-4.5, m2(0,0));
            
            BOOST_CHECK_EQUAL(4.0, m1(1,0));
            BOOST_CHECK_EQUAL(4.0, m2(1,0));
            
            m2(1,0) = 6.7;
            
            BOOST_CHECK_EQUAL(6.7, m1(1,0));
            BOOST_CHECK_EQUAL(6.7, m2(1,0));

        }
    }

    namespace StandardMatrixOperationsUnitTests
    {
        BOOST_AUTO_TEST_CASE(TestAddition)
        {
            double lhs_buf[] = {  2.0, 4.0, 6.0, 8.0,
                                  10.0, 12.0, 14.0, 16.0,
                                  18.0, 20.0, 22.0, 24.0,
                                  26.0, 28.0, 30.0, 32.0};
            double rhs_buf[] = { 3.0, 6.0, 9.0, 12.0,
                                  15.0, 18.0, 21.0, 24.0,
                                  27.0, 30.0, 33.0, 36.0,
                                  39.0, 42.0, 45.0, 48.0};

            NekMatrix<double> lhs1(4, 4, lhs_buf);
            std::shared_ptr<NekMatrix<NekMatrix<double>, ScaledMatrixTag> > lhs2;
            std::shared_ptr<NekMatrix<NekMatrix<double>, BlockMatrixTag> > lhs3;

            NekMatrix<double> rhs1(4, 4, rhs_buf);
            std::shared_ptr<NekMatrix<NekMatrix<double>, ScaledMatrixTag> > rhs2;
            std::shared_ptr<NekMatrix<NekMatrix<double>, BlockMatrixTag> > rhs3;

            GenerateMatrices(lhs1, 2.0, 2, 2, lhs2, lhs3);
            GenerateMatrices(rhs1, 3.0, 2, 2, rhs2, rhs3);

            double result_buf[] = {5, 10, 15, 20,
                                   25, 30, 35, 40,
                                   45, 50, 55, 60,
                                   65, 70, 75, 80};
            NekMatrix<double> result(4, 4, result_buf);

            RunAllAddCombinations(lhs1, *lhs2, *lhs3, rhs1, *rhs2, *rhs3, result);

            NekMatrix<double> rhs_transposed = Transpose(rhs1);
            NekMatrix<double> result_with_one_operand_transposed = lhs1 + rhs_transposed;

            double expected_result_with_one_operand_transposed[] = 
                { 5, 19, 33, 47,
                  16, 30, 44, 58,
                  27, 41, 55, 69,
                  38, 52, 66, 80};
            BOOST_CHECK_EQUAL(NekMatrix<double>(4, 4, expected_result_with_one_operand_transposed),
                result_with_one_operand_transposed);
        }

        BOOST_AUTO_TEST_CASE(TestSubtraction)
        {
            double lhs_buf[] = {  2.0, 4.0, 6.0, 8.0,
                                  10.0, 12.0, 14.0, 16.0,
                                  18.0, 20.0, 22.0, 24.0,
                                  26.0, 28.0, 30.0, 32.0};
            double rhs_buf[] = { 3.0, 6.0, 9.0, 12.0,
                                  15.0, 18.0, 21.0, 24.0,
                                  27.0, 30.0, 33.0, 36.0,
                                  39.0, 42.0, 45.0, 48.0};

            NekMatrix<double> lhs1(4, 4, lhs_buf);
            std::shared_ptr<NekMatrix<NekMatrix<double>, ScaledMatrixTag> > lhs2;
            std::shared_ptr<NekMatrix<NekMatrix<double>, BlockMatrixTag> > lhs3;

            NekMatrix<double> rhs1(4, 4, rhs_buf);
            std::shared_ptr<NekMatrix<NekMatrix<double>, ScaledMatrixTag> > rhs2;
            std::shared_ptr<NekMatrix<NekMatrix<double>, BlockMatrixTag> > rhs3;

            GenerateMatrices(lhs1, 2.0, 2, 2, lhs2, lhs3);
            GenerateMatrices(rhs1, 3.0, 2, 2, rhs2, rhs3);

            double result_buf[] = {-1, -2, -3, -4,
                                   -5, -6, -7, -8,
                                   -9, -10, -11, -12,
                                   -13, -14, -15, -16};
            NekMatrix<double> result(4, 4, result_buf);
            
            RunAllSubCombinations(lhs1, *lhs2, *lhs3, rhs1, *rhs2, *rhs3, result);
        }

        BOOST_AUTO_TEST_CASE(TestThreeAdditions)
        {
            double lhs_buf[] = {  2.0, 4.0, 6.0, 8.0,
                                  10.0, 12.0, 14.0, 16.0,
                                  18.0, 20.0, 22.0, 24.0,
                                  26.0, 28.0, 30.0, 32.0};
            double middle_buf[] = {  2.0, 4.0, 6.0, 8.0,
                                    10.0, 12.0, 14.0, 16.0,
                                    18.0, 20.0, 22.0, 24.0,
                                    26.0, 28.0, 30.0, 32.0};
            double rhs_buf[] = { 3.0, 6.0, 9.0, 12.0,
                                  15.0, 18.0, 21.0, 24.0,
                                  27.0, 30.0, 33.0, 36.0,
                                  39.0, 42.0, 45.0, 48.0};

            NekMatrix<double> lhs(4, 4, lhs_buf);
            NekMatrix<double> middle(4, 4, middle_buf);
            NekMatrix<double> rhs(4, 4, rhs_buf);

            NekMatrix<double> result = lhs+middle+rhs;
        }

        BOOST_AUTO_TEST_CASE(TestTransposedMultiplication)
        {
            {
                double buf1[] = { 1, 4,
                                  2, 5,
                                  3, 6};
                double buf2[] = { 10, 11, 12,
                                  13, 14, 15};

                double transposed_buf1[] = {1, 2, 3,
                                            4, 5, 6};

                double transposed_buf2[] = {10, 13,
                                            11, 14,
                                            12, 15};

                double result_buf[] = {68, 167, 
                                       86, 212};

                double transposed_result_buf[] = {62, 85, 108,
                                                  67, 92, 117,
                                                  72, 99, 126};

                NekMatrix<double> m1(2, 3, buf1);
                NekMatrix<double> m2(3, 2, buf2);
                NekMatrix<double> transposed_m1(3, 2, transposed_buf1);
                NekMatrix<double> transposed_m2(2, 3, transposed_buf2);

                NekMatrix<double> expected_result(2, 2, result_buf);
                NekMatrix<double> expected_transposed_result(3, 3, transposed_result_buf);

                BOOST_CHECK_EQUAL(expected_result, m1*m2);
                BOOST_CHECK_EQUAL(expected_transposed_result, transposed_m1*transposed_m2);
                m1.Transpose();
                m2.Transpose();
                BOOST_CHECK_EQUAL(expected_transposed_result, m1*m2);
                m1.Transpose();
                m2.Transpose();
                BOOST_CHECK_EQUAL(expected_result, m1*m2);

                NekMatrix<double> tm1 = Transpose(m1);
                NekMatrix<double> tm2 = Transpose(m2);
                BOOST_CHECK_EQUAL(expected_result, m1*m2);
                BOOST_CHECK_EQUAL(expected_transposed_result, tm1*tm2);
            }
            
        }

        BOOST_AUTO_TEST_CASE(TestGlobalTransposeMethod)
        {
            double buf1[] = { 1, 4,
                                  2, 5,
                                  3, 6};
            double buf2[] = { 10, 11, 12,
                              13, 14, 15};

            double transposed_buf1[] = {1, 2, 3,
                                        4, 5, 6};

            double transposed_buf2[] = {10, 13,
                                        11, 14,
                                        12, 15};

            NekMatrix<double> m1(2, 3, buf1);
            NekMatrix<double> m2(3, 2, buf2);

            NekMatrix<double> tm1 = Transpose(m1);
            NekMatrix<double> tm2 = Transpose(m2);

            NekMatrix<double> expected_tm1(3, 2, transposed_buf1);
            NekMatrix<double> expected_tm2(2, 3, transposed_buf2);

            BOOST_CHECK_EQUAL(tm1, expected_tm1);
            BOOST_CHECK_EQUAL(tm2, expected_tm2);

            BOOST_CHECK_EQUAL(m1(0, 1), 2);
            BOOST_CHECK_EQUAL(tm1(1, 0), 2);

            m1(0, 1) = 89;
            BOOST_CHECK_EQUAL(m1(0, 1), 89);
            BOOST_CHECK_EQUAL(tm1(1, 0), 89);

            BOOST_CHECK_EQUAL(m1.GetTransposeFlag(), 'N');
            BOOST_CHECK_EQUAL(tm1.GetTransposeFlag(), 'T');
        }
        
        BOOST_AUTO_TEST_CASE(TestScalarMatrixMultiply)
        {
            double buf1[] = { 1, 4,
                                  2, 5,
                                  3, 6};
                                  
            NekMatrix<double> m1(2, 3, buf1);
            NekMatrix<double> m2 = m1*2.0;
            NekMatrix<double> m3 = 3.0*m1;
            
            double m2_expected_result_buf[] = { 2.0, 8.0,
                                          4.0, 10.0,
                                          6.0, 12.0 };
                                          
            double m3_expected_result_buf[] = { 3.0, 12.0,
                                          6.0, 15.0,
                                          9.0, 18.0 };
                                          
            NekMatrix<double> m2_expected_result(2, 3, m2_expected_result_buf);
            NekMatrix<double> m3_expected_result(2, 3, m3_expected_result_buf);
            
            BOOST_CHECK_EQUAL(m2_expected_result, m2);
            BOOST_CHECK_EQUAL(m3_expected_result, m3);
        }

    }

}


