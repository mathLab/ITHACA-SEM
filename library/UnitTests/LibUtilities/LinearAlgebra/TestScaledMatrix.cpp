///////////////////////////////////////////////////////////////////////////////
//
// File: TestScaledMatrix.cpp
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

#include <boost/test/unit_test.hpp>

#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <UnitTests/LibUtilities/LinearAlgebra/TestCombinationRunner.h>

namespace Nektar
{
    namespace ScaledMatrixUnitTests
    {
        typedef NekMatrix<double> InnerMatrix;
        typedef NekMatrix<InnerMatrix, FullMatrixTag, ScaledMatrixTag> SMat;
        
        typedef NekMatrix<unsigned int> IntInnerMatrix;
        typedef NekMatrix<IntInnerMatrix, FullMatrixTag, ScaledMatrixTag> IntSMat;

        BOOST_AUTO_TEST_CASE(TestDefaultConstructor)
        {
            SMat m;
            BOOST_CHECK_EQUAL(0, m.GetRows());
            BOOST_CHECK_EQUAL(0, m.GetColumns());
            BOOST_CHECK(boost::shared_ptr<InnerMatrix>() != m.GetOwnedMatrix());
            BOOST_CHECK_EQUAL('N', m.GetTransposeFlag());
            BOOST_CHECK_EQUAL(0, m.GetStorageSize());
            BOOST_CHECK_EQUAL(eFULL, m.GetStorageType());
            BOOST_CHECK_EQUAL(0.0, m.Scale());

        }

        BOOST_AUTO_TEST_CASE(TestConstructorWithInnerMatrix)
        {
            double buf[] = {1.0, 2.0, 3.0,
                            4.0, 5.0, 6.0};

            boost::shared_ptr<InnerMatrix> in(new InnerMatrix(3, 2, buf));
            SMat m(2.7, in);
            BOOST_CHECK_EQUAL(3, m.GetRows());
            BOOST_CHECK_EQUAL(2, m.GetColumns());
            BOOST_CHECK_EQUAL(in, m.GetOwnedMatrix());
            BOOST_CHECK_EQUAL('N', m.GetTransposeFlag());
            BOOST_CHECK_EQUAL(6, m.GetStorageSize());
            BOOST_CHECK_EQUAL(eFULL, m.GetStorageType());
            BOOST_CHECK_EQUAL(2.7, m.Scale());

        }

        BOOST_AUTO_TEST_CASE(TestElementAccess)
        {
            double buf[] = {1.0, 2.0, 3.0,
                            4.0, 5.0, 6.0};

            boost::shared_ptr<InnerMatrix> in1(new InnerMatrix(3, 2, buf));
            boost::shared_ptr<InnerMatrix> in2(new InnerMatrix(3, 2, buf));
            boost::shared_ptr<InnerMatrix> in3(new InnerMatrix(3, 2, buf));
            boost::shared_ptr<InnerMatrix> in4(new InnerMatrix(3, 2, buf));
            in3->Transpose();
            in4->Transpose();

            SMat m1(1.0, in1);
            SMat m2(2.0, in2);
            SMat m3(3.0, in3);
            SMat m4(4.0, in4);

            m2.Transpose();
            m4.Transpose();

            // m1 is N/N
            // m2 is N/T
            // m3 is T/N
            // m4 is T/T
            
            BOOST_CHECK_EQUAL(m1(0, 0), 1.0);
            BOOST_CHECK_EQUAL(m1(1, 0), 2.0);
            BOOST_CHECK_EQUAL(m1(2, 0), 3.0);
            BOOST_CHECK_EQUAL(m1(0, 1), 4.0);
            BOOST_CHECK_EQUAL(m1(1, 1), 5.0);
            BOOST_CHECK_EQUAL(m1(2, 1), 6.0);

            BOOST_CHECK_EQUAL(m4(0, 0), 4.0);
            BOOST_CHECK_EQUAL(m4(1, 0), 8.0);
            BOOST_CHECK_EQUAL(m4(2, 0), 12.0);
            BOOST_CHECK_EQUAL(m4(0, 1), 16.0);
            BOOST_CHECK_EQUAL(m4(1, 1), 20.0);
            BOOST_CHECK_EQUAL(m4(2, 1), 24.0);

            BOOST_CHECK_EQUAL(m2(0, 0), 2.0);
            BOOST_CHECK_EQUAL(m2(0, 1), 4.0);
            BOOST_CHECK_EQUAL(m2(0, 2), 6.0);
            BOOST_CHECK_EQUAL(m2(1, 0), 8.0);
            BOOST_CHECK_EQUAL(m2(1, 1), 10.0);
            BOOST_CHECK_EQUAL(m2(1, 2), 12.0);

            BOOST_CHECK_EQUAL(m3(0, 0), 3.0);
            BOOST_CHECK_EQUAL(m3(0, 1), 6.0);
            BOOST_CHECK_EQUAL(m3(0, 2), 9.0);
            BOOST_CHECK_EQUAL(m3(1, 0), 12.0);
            BOOST_CHECK_EQUAL(m3(1, 1), 15.0);
            BOOST_CHECK_EQUAL(m3(1, 2), 18.0);
        }


        BOOST_AUTO_TEST_CASE(TestIteration)
        {
            double buf[] = {1.0, 2.0, 3.0,
                            4.0, 5.0, 6.0};

            boost::shared_ptr<InnerMatrix> in1(new InnerMatrix(3, 2, buf));
            boost::shared_ptr<InnerMatrix> in2(new InnerMatrix(3, 2, buf));
            boost::shared_ptr<InnerMatrix> in3(new InnerMatrix(3, 2, buf));
            boost::shared_ptr<InnerMatrix> in4(new InnerMatrix(3, 2, buf));
            in3->Transpose();
            in4->Transpose();

            SMat m1(1.0, in1);
            SMat m2(2.0, in2);
            SMat m3(3.0, in3);
            SMat m4(4.0, in4);

            m2.Transpose();
            m4.Transpose();

            // m1 is N/N
            // m2 is N/T
            // m3 is T/N
            // m4 is T/T
            
            SMat::const_iterator i1 = m1.begin();
            SMat::const_iterator i2 = m2.begin();
            SMat::const_iterator i3 = m3.begin();
            SMat::const_iterator i4 = m4.begin();

            BOOST_CHECK(1.0 == *(i1++));
            BOOST_CHECK(2.0 == *(i1++));
            BOOST_CHECK(3.0 == *(i1++));
            BOOST_CHECK(4.0 == *(i1++));
            BOOST_CHECK(5.0 == *(i1++));
            BOOST_CHECK(6.0 == *(i1++));
            BOOST_CHECK(m1.end() == i1);

            BOOST_CHECK(2.0 == *(i2++));
            BOOST_CHECK(8.0 == *(i2++));
            BOOST_CHECK(4.0 == *(i2++));
            BOOST_CHECK(10.0 == *(i2++));
            BOOST_CHECK(6.0 == *(i2++));
            BOOST_CHECK(12.0 == *(i2++));
            BOOST_CHECK(m2.end() == i2);

            BOOST_CHECK(3.0 == *(i3++));
            BOOST_CHECK(12.0 == *(i3++));
            BOOST_CHECK(6.0 == *(i3++));
            BOOST_CHECK(15.0 == *(i3++));
            BOOST_CHECK(9.0 == *(i3++));
            BOOST_CHECK(18.0 == *(i3++));
            BOOST_CHECK(m3.end() == i3);

            BOOST_CHECK(4.0 == *(i4++));
            BOOST_CHECK(8.0 == *(i4++));
            BOOST_CHECK(12.0 == *(i4++));
            BOOST_CHECK(16.0 == *(i4++));
            BOOST_CHECK(20.0 == *(i4++));
            BOOST_CHECK(24.0 == *(i4++));
            BOOST_CHECK(m4.end() == i4);
        }

        BOOST_AUTO_TEST_CASE(TestNMatrixNMatrixMultiplication)
        {
            {
                double lhs_buf[] = {1.0, 2.0, 3.0,
                                    4.0, 5.0, 6.0};
                double rhs_buf[] = {7.0, 8.0,
                                    9.0, 10.0,
                                    11.0, 12.0};
                boost::shared_ptr<InnerMatrix> in1(new InnerMatrix(3, 2, lhs_buf));
                boost::shared_ptr<InnerMatrix> in2(new InnerMatrix(2, 3, rhs_buf));

                SMat m1(1.0, in1);
                SMat m2(2.0, in2);

                double expected_result_buf[] = {78.0, 108.0, 138.0,
                                          98.0, 136.0, 174.0,
                                          118.0, 164.0, 210.0 };
                InnerMatrix expected_result(3, 3, expected_result_buf);
                InnerMatrix result1 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result1);
                
                in1->Transpose();
                m1.Transpose();
                InnerMatrix result2 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result2);

                in2->Transpose();
                m2.Transpose();
                InnerMatrix result3 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result3);

                in1->Transpose();
                m1.Transpose();
                InnerMatrix result4 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result4);
            }

            {
                unsigned int lhs_buf[] = {1, 2, 3,
                                    4, 5, 6};
                unsigned int rhs_buf[] = {7, 8,
                                    9, 10,
                                    11, 12};
                boost::shared_ptr<IntInnerMatrix> in1(new IntInnerMatrix(3, 2, lhs_buf));
                boost::shared_ptr<IntInnerMatrix> in2(new IntInnerMatrix(2, 3, rhs_buf));

                IntSMat m1(1, in1);
                IntSMat m2(2, in2);

                unsigned int expected_result_buf[] = {78, 108, 138,
                                          98, 136, 174,
                                          118, 164, 210 };
                IntInnerMatrix expected_result(3, 3, expected_result_buf);
                IntInnerMatrix result1 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result1);
                
                in1->Transpose();
                m1.Transpose();
                IntInnerMatrix result2 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result2);

                in2->Transpose();
                m2.Transpose();
                IntInnerMatrix result3 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result3);

                in1->Transpose();
                m1.Transpose();
                IntInnerMatrix result4 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result4);
            }

        }

        BOOST_AUTO_TEST_CASE(TestTMatrixNMatrixMultiplication)
        {
            {
                double lhs_buf[] = {1.0, 2.0, 
                                    3.0, 4.0, 
                                    5.0, 6.0};
                double rhs_buf[] = {7.0, 8.0,
                                    9.0, 10.0,
                                    11.0, 12.0};
                boost::shared_ptr<InnerMatrix> in1(new InnerMatrix(2, 3, lhs_buf));
                boost::shared_ptr<InnerMatrix> in2(new InnerMatrix(2, 3, rhs_buf));

                SMat m1(1.0, in1);
                SMat m2(2.0, in2);

                double expected_result_buf[] = {46.0, 106.0, 166.0,
                                          58.0, 134.0, 210.0,
                                          70.0, 162.0, 254.0 };
                InnerMatrix expected_result(3, 3, expected_result_buf);

                in1->Transpose();
                InnerMatrix result1 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result1);

                in1->Transpose();
                m1.Transpose();
                InnerMatrix result2 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result2);

                in1->Transpose();
                m1.Transpose();
                in2->Transpose();
                m2.Transpose();
                InnerMatrix result3 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result3);

                in1->Transpose();
                m1.Transpose();
                InnerMatrix result4 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result4);
            }

            {
                unsigned int lhs_buf[] = {1, 2, 
                                    3, 4, 
                                    5, 6};
                unsigned int rhs_buf[] = {7, 8,
                                    9, 10,
                                    11, 12};
                boost::shared_ptr<IntInnerMatrix> in1(new IntInnerMatrix(2, 3, lhs_buf));
                boost::shared_ptr<IntInnerMatrix> in2(new IntInnerMatrix(2, 3, rhs_buf));

                IntSMat m1(1, in1);
                IntSMat m2(2, in2);

                unsigned int expected_result_buf[] = {46, 106, 166,
                                          58, 134, 210,
                                          70, 162, 254 };
                IntInnerMatrix expected_result(3, 3, expected_result_buf);

                in1->Transpose();
                IntInnerMatrix result1 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result1);

                in1->Transpose();
                m1.Transpose();
                IntInnerMatrix result2 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result2);

                in1->Transpose();
                m1.Transpose();
                in2->Transpose();
                m2.Transpose();
                IntInnerMatrix result3 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result3);

                in1->Transpose();
                m1.Transpose();
                IntInnerMatrix result4 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result4);
            }
        }

        BOOST_AUTO_TEST_CASE(TestNMatrixTMatrixMultiplication)
        {
            {
                double lhs_buf[] = {1.0, 2.0, 
                                    3.0, 4.0, 
                                    5.0, 6.0};
                double rhs_buf[] = {7.0, 8.0,
                                    9.0, 10.0,
                                    11.0, 12.0};
                boost::shared_ptr<InnerMatrix> in1(new InnerMatrix(2, 3, lhs_buf));
                boost::shared_ptr<InnerMatrix> in2(new InnerMatrix(2, 3, rhs_buf));

                SMat m1(1.0, in1);
                SMat m2(2.0, in2);

                double expected_result_buf[] = {178.0, 232.0,
                                                196.0, 256.0 };
                InnerMatrix expected_result(2, 2, expected_result_buf);

                in2->Transpose();
                InnerMatrix result1 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result1);

                in2->Transpose();
                m2.Transpose();
                InnerMatrix result2 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result2);

                in2->Transpose();
                m2.Transpose();
                in1->Transpose();
                m1.Transpose();
                InnerMatrix result3 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result3);

                in2->Transpose();
                m2.Transpose();
                InnerMatrix result4 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result4);
            }

            {
                unsigned int lhs_buf[] = {1, 2, 
                                    3, 4, 
                                    5, 6};
                unsigned int rhs_buf[] = {7, 8,
                                    9, 10,
                                    11, 12};
                boost::shared_ptr<IntInnerMatrix> in1(new IntInnerMatrix(2, 3, lhs_buf));
                boost::shared_ptr<IntInnerMatrix> in2(new IntInnerMatrix(2, 3, rhs_buf));

                IntSMat m1(1, in1);
                IntSMat m2(2, in2);

                unsigned int expected_result_buf[] = {178, 232,
                                                196, 256 };
                IntInnerMatrix expected_result(2, 2, expected_result_buf);

                in2->Transpose();
                IntInnerMatrix result1 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result1);

                in2->Transpose();
                m2.Transpose();
                IntInnerMatrix result2 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result2);

                in2->Transpose();
                m2.Transpose();
                in1->Transpose();
                m1.Transpose();
                IntInnerMatrix result3 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result3);

                in2->Transpose();
                m2.Transpose();
                IntInnerMatrix result4 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result4);
            }
        }

        BOOST_AUTO_TEST_CASE(TestTMatrixTMatrixMultiplication)
        {
            {
                double lhs_buf[] = {1.0, 2.0, 3.0,
                                    4.0, 5.0, 6.0};
                double rhs_buf[] = {7.0, 8.0,
                                    9.0, 10.0,
                                    11.0, 12.0};
                boost::shared_ptr<InnerMatrix> in1(new InnerMatrix(3, 2, lhs_buf));
                boost::shared_ptr<InnerMatrix> in2(new InnerMatrix(2, 3, rhs_buf));

                SMat m1(1.0, in1);
                SMat m2(2.0, in2);

                double expected_result_buf[] = {116.0, 278.0,
                                                128.0, 308.0};
                InnerMatrix expected_result(2, 2, expected_result_buf);

                in1->Transpose();
                in2->Transpose();
                InnerMatrix result1 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result1);
                
                in1->Transpose();
                m1.Transpose();
                InnerMatrix result2 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result2);

                in2->Transpose();
                m2.Transpose();
                InnerMatrix result3 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result3);

                in1->Transpose();
                m1.Transpose();
                InnerMatrix result4 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result4);
            }

            {
                unsigned int lhs_buf[] = {1, 2, 3,
                                    4, 5, 6};
                unsigned int rhs_buf[] = {7, 8,
                                    9, 10,
                                    11, 12};
                boost::shared_ptr<IntInnerMatrix> in1(new IntInnerMatrix(3, 2, lhs_buf));
                boost::shared_ptr<IntInnerMatrix> in2(new IntInnerMatrix(2, 3, rhs_buf));

                IntSMat m1(1, in1);
                IntSMat m2(2, in2);

                unsigned int expected_result_buf[] = {116, 278,
                                                128, 308};
                IntInnerMatrix expected_result(2, 2, expected_result_buf);

                in1->Transpose();
                in2->Transpose();
                IntInnerMatrix result1 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result1);
                
                in1->Transpose();
                m1.Transpose();
                IntInnerMatrix result2 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result2);

                in2->Transpose();
                m2.Transpose();
                IntInnerMatrix result3 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result3);

                in1->Transpose();
                m1.Transpose();
                IntInnerMatrix result4 = m1*m2;
                BOOST_CHECK_EQUAL(expected_result, result4);
            }
        }

        BOOST_AUTO_TEST_CASE(TestScaledNMatrixVectorMultiply)
        {
            {
                double lhs_buf[] = {1.0, 3.0, 5.0, 7.0,
                    2.0, 4.0, 6.0, 8.0};
                double rhs_buf[] = {1.0, 2.0};

                boost::shared_ptr<InnerMatrix> in1(new InnerMatrix(4, 2, lhs_buf));
                NekVector<double> rhs(2, rhs_buf);
                SMat m1(2.0, in1);

                double expected_result_buf[] = {10, 22, 34, 46};
                NekVector<double> expected_result(4, expected_result_buf);

                NekVector<double> result1 = m1*rhs;
                BOOST_CHECK_EQUAL(expected_result, result1);

                in1->Transpose();
                m1.Transpose();
                NekVector<double> result2 = m1*rhs;
                BOOST_CHECK_EQUAL(expected_result, result2);
            }

            {
                unsigned int lhs_buf[] = {1, 3, 5, 7,
                    2, 4, 6, 8};
                unsigned int rhs_buf[] = {1, 2};

                boost::shared_ptr<IntInnerMatrix> in1(new IntInnerMatrix(4, 2, lhs_buf));
                NekVector<unsigned int> rhs(2, rhs_buf);
                IntSMat m1(2, in1);

                unsigned int expected_result_buf[] = {10, 22, 34, 46};
                NekVector<unsigned int> expected_result(4, expected_result_buf);

                NekVector<unsigned int> result1 = m1*rhs;
                BOOST_CHECK_EQUAL(expected_result, result1);

                in1->Transpose();
                m1.Transpose();
                NekVector<unsigned int> result2 = m1*rhs;
                BOOST_CHECK_EQUAL(expected_result, result2);
            }
        }

        BOOST_AUTO_TEST_CASE(TestScaledTMatrixVectorMultiply)
        {
            {
                double lhs_buf[] = {1.0, 2.0,
                                    3.0, 4.0,
                                    5.0, 6.0,
                                    7.0, 8.0};
                double rhs_buf[] = {1.0, 2.0};

                boost::shared_ptr<InnerMatrix> in1(new InnerMatrix(2, 4, lhs_buf));
                NekVector<double> rhs(2, rhs_buf);
                SMat m1(2.0, in1);

                double expected_result_buf[] = {10, 22, 34, 46};
                NekVector<double> expected_result(4, expected_result_buf);

                in1->Transpose();
                NekVector<double> result1 = m1*rhs;
                BOOST_CHECK_EQUAL(expected_result, result1);

                in1->Transpose();
                m1.Transpose();
                NekVector<double> result2 = m1*rhs;
                BOOST_CHECK_EQUAL(expected_result, result2);
            }

            {
                unsigned int lhs_buf[] = {1, 2,
                                    3, 4,
                                    5, 6,
                                    7, 8};
                unsigned int rhs_buf[] = {1, 2};

                boost::shared_ptr<IntInnerMatrix> in1(new IntInnerMatrix(2, 4, lhs_buf));
                NekVector<unsigned int> rhs(2, rhs_buf);
                IntSMat m1(2, in1);

                unsigned int expected_result_buf[] = {10, 22, 34, 46};
                NekVector<unsigned int> expected_result(4, expected_result_buf);

                in1->Transpose();
                NekVector<unsigned int> result1 = m1*rhs;
                BOOST_CHECK_EQUAL(expected_result, result1);

                in1->Transpose();
                m1.Transpose();
                NekVector<unsigned int> result2 = m1*rhs;
                BOOST_CHECK_EQUAL(expected_result, result2);
            }
        }
    }
}


