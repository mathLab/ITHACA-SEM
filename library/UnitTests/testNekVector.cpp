///////////////////////////////////////////////////////////////////////////////
//
// File: testNekVector.cpp
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
// Description: Test code for NekVector
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/NekVector.hpp>
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

        class VectorTestClass
        {
            public:
                VectorTestClass() : m_dataValue(0) {}
                explicit VectorTestClass(int data) : m_dataValue(data) {}
                VectorTestClass(const VectorTestClass& in) = default;
                ~VectorTestClass(){}

                VectorTestClass operator=(const VectorTestClass& rhs)
                {
                    m_dataValue = rhs.m_dataValue;
                    return *this;
                }

                int value() const { return m_dataValue; }

            private:
                int m_dataValue;
        };

        bool operator==(const VectorTestClass& lhs, const VectorTestClass& rhs)
        {
            return lhs.value() == rhs.value();
        }

        bool operator!=(const VectorTestClass& lhs, const VectorTestClass& rhs)
        {
            return !(lhs == rhs);
        }

        struct TenD 
        {
            static const unsigned int Value = 10;
        };

        BOOST_AUTO_TEST_CASE(TestNekVectorConstruction)
        {
            // Test the constructors on a numeric type.
            {
                NekVector<double> p1(10, 2.7);
                BOOST_CHECK(p1.GetDimension() == 10);
                for(unsigned int i = 0; i < p1.GetDimension(); ++i)
                {
                    BOOST_CHECK(p1(i) == 2.7);
                    BOOST_CHECK(p1[i] == 2.7);
                }

                NekVector<double> p2(1,0.0);
                BOOST_CHECK(p2.GetDimension() == 1);
                BOOST_CHECK(p2(0) == 0.0);

                NekVector<double> p3(p1);
                BOOST_CHECK(p3.GetDimension() == 10);
                for(unsigned int i = 0; i < p3.GetDimension(); ++i)
                {
                    BOOST_CHECK(p3(i) == 2.7);
                    BOOST_CHECK(p3[i] == 2.7);
                }

                p2 = p3;
                BOOST_CHECK(p2.GetDimension() == 10);
                for(unsigned int i = 0; i < p2.GetDimension(); ++i)
                {
                    BOOST_CHECK_EQUAL(p2(i), 2.7);
                    BOOST_CHECK_EQUAL(p2[i], 2.7);
                }
            }
        }

        BOOST_AUTO_TEST_CASE(TestNekVectorOperators)
        {
            NekVector<double> v1(3, 1.0);
            v1(0) = 1.1;
            v1(1) = 1.2;
            v1(2) = 1.3;

            BOOST_CHECK(v1(0) == 1.1);
            BOOST_CHECK(v1(1) == 1.2);
            BOOST_CHECK(v1(2) == 1.3);

            v1[0] = 1.4;
            v1[1] = 1.5;
            v1[2] = 1.6;

            BOOST_CHECK(v1[0] == 1.4);
            BOOST_CHECK(v1[1] == 1.5);
            BOOST_CHECK(v1[2] == 1.6);

            v1.x() = 1.7;
            v1.y() = 1.8;
            v1.z() = 1.9;

            BOOST_CHECK(v1[0] == 1.7);
            BOOST_CHECK(v1[1] == 1.8);
            BOOST_CHECK(v1[2] == 1.9);

            BOOST_CHECK(v1.x() == 1.7);
            BOOST_CHECK(v1.y() == 1.8);
            BOOST_CHECK(v1.z() == 1.9);

            NekVector<double> v2(3, 1.0);
            v2.x() = 1.7;
            v2.y() = 1.8;
            v2.z() = 1.9;

            BOOST_CHECK(v1 == v2);
            BOOST_CHECK(!(v1 != v2));

            NekVector<double> v3(4, 2.0);
            BOOST_CHECK(v3 != v1);

            NekVector<double> v4(3, 0.0);
            BOOST_CHECK(v4 != v1);

        }

        BOOST_AUTO_TEST_CASE(TestNekVectorArithmetic)
        {

        }
        
        BOOST_AUTO_TEST_CASE(TestNorms)
        {
            double vals[] = {1,-2,3};
            NekVector<double> v(3, vals);
                
            double epsilon = 1e-11;
            BOOST_CHECK_EQUAL(v.L1Norm(), 1.0 + 2.0 + 3.0);
            BOOST_CHECK_CLOSE(v.L2Norm(), sqrt(1.0+4.0+9.0), epsilon);
            BOOST_CHECK_EQUAL(v.InfinityNorm(), 3);
        }

        BOOST_AUTO_TEST_CASE(TestMatrixVectorMultiply)
        {
            {
                //double matrix_buf[] = {1.0, 2.0, 3.0,
                //                    4.0, 5.0, 6.0,
                //                    7.0, 8.0, 9.0};
                double matrix_buf[] = {1.0, 4.0, 7.0,
                                       2.0, 5.0, 8.0,
                                       3.0, 6.0, 9.0};
                double vector_buf[] = {20.0, 30.0, 40.0};
                
                NekMatrix<double> m(3, 3, matrix_buf);
                NekVector<double> v(3, vector_buf);
                NekVector<double> result = m*v;
                
                double result_buf[] = {200.0, 470.0, 740.0};
                NekVector<double> expected_result(3, result_buf);
                BOOST_CHECK_EQUAL(result, expected_result);
                BOOST_CHECK_EQUAL(3u, result.GetDimension());
            }
            
            {
                //double matrix_buf[] = {1.0, 2.0, 3.0,
                //                    4.0, 5.0, 6.0,
                //                    7.0, 8.0, 9.0};
                double matrix_buf[] = {1.0, 4.0, 7.0,
                                       2.0, 5.0, 8.0,
                                       3.0, 6.0, 9.0};
                double vector_buf[] = {20.0, 30.0, 40.0};
                
                std::shared_ptr<NekMatrix<double> > m(new NekMatrix<double>(3, 3, matrix_buf));
                NekMatrix<NekMatrix<double>, ScaledMatrixTag> s(2.0, m);
                NekVector<double> v(3, vector_buf);
                NekVector<double> result = s*v;
                
                double result_buf[] = {400.0, 940.0, 1480.0};
                NekVector<double> expected_result(3, result_buf);
                BOOST_CHECK_EQUAL(result, expected_result);
                BOOST_CHECK_EQUAL(3u, result.GetDimension());
            }
            
            {
                //double m1_buf[] = {1.0, 2.0, 3.0, 4.0};
                //double m2_buf[] = {5.0, 6.0, 7.0, 8.0};
                //double m3_buf[] = {9.0, 10.0, 11.0, 12.0};
                double m1_buf[] = {1.0, 3.0, 2.0, 4.0};
                double m2_buf[] = {5.0, 7.0, 6.0, 8.0};
                double m3_buf[] = {9.0, 11.0, 10.0, 12.0};
                double vector_buf[] = {20.0, 30.0};
                
                std::shared_ptr<NekMatrix<double> > m1(new NekMatrix<double>(2, 2, m1_buf));
                std::shared_ptr<NekMatrix<double> > m2(new NekMatrix<double>(2, 2, m2_buf));
                std::shared_ptr<NekMatrix<double> > m3(new NekMatrix<double>(2, 2, m3_buf));
                
                NekMatrix<NekMatrix<double>, BlockMatrixTag> b(3, 1, 2, 2);
                b.SetBlock(0,0,m1);
                b.SetBlock(1,0,m2);
                b.SetBlock(2,0,m3);
                
                NekVector<double> v(2, vector_buf);
                NekVector<double> result = b*v;
                
                double result_buf[] = {80.0, 180.0, 280.0, 380.0, 480.0, 580.0};
                NekVector<double> expected_result(6, result_buf);
                BOOST_CHECK_EQUAL(result, expected_result);
                BOOST_CHECK_EQUAL(6u, result.GetDimension());
            }

            {
                double m_buf[] = {1, 4, 7,
                                   2, 5, 8,
                                   3, 6, 9};
                double v_buf[] = {10, 11, 12};
                double expected_result_buf[] = {68, 167, 266};
                double expected_transpose_result_buf[] = {138, 171, 204};

                NekMatrix<double> m(3, 3, m_buf);
                NekVector<double> v_variable(3, v_buf);
                NekVector<double> v_constant(3, v_buf);
                NekVector<double> expected_result(3, expected_result_buf);
                NekVector<double> expected_transpose_result(3, expected_transpose_result_buf);

                NekVector<double> constantResult = m*v_constant;
                NekVector<double> variableResult = m*v_variable;
                BOOST_CHECK_EQUAL(variableResult, expected_result);
                BOOST_CHECK_EQUAL(constantResult, expected_result);

                m.Transpose();
                constantResult = m*v_constant;
                variableResult = m*v_variable;
                BOOST_CHECK_EQUAL(variableResult, expected_transpose_result);
                BOOST_CHECK_EQUAL(constantResult, expected_transpose_result);
                
                m.Transpose();
                constantResult = m*v_constant;
                variableResult = m*v_variable;
                BOOST_CHECK_EQUAL(variableResult, expected_result);
                BOOST_CHECK_EQUAL(constantResult, expected_result);

                NekMatrix<double> transposed = Transpose(m);
                constantResult = transposed*v_constant;
                variableResult = transposed*v_variable;
                BOOST_CHECK_EQUAL(variableResult, expected_transpose_result);
                BOOST_CHECK_EQUAL(constantResult, expected_transpose_result);
            }

            {
                double m_buf[] = {1, 4, 
                                  2, 5,
                                  3, 6};
                double v_non_transposed_buf[] = {10, 11, 12};
                double v_transposed_buf[] = {20, 21};
                double expected_result_buf[] = {68, 167};
                double expected_transpose_result_buf[] = {104, 145, 186};

                NekMatrix<double> m(2, 3, m_buf);
                NekVector<double> v_non_transposed_variable(3, v_non_transposed_buf);
                NekVector<double> v_non_transposed_constant(3, v_non_transposed_buf);
                NekVector<double> v_transposed_variable(2, v_transposed_buf);
                NekVector<double> v_transposed_constant(2, v_transposed_buf);
                NekVector<double> expected_result(2, expected_result_buf);
                NekVector<double> expected_transpose_result(3, expected_transpose_result_buf);

                NekVector<double> variableResult = m*v_non_transposed_variable;
                NekVector<double> constantResult = m*v_non_transposed_constant;
                BOOST_CHECK_EQUAL(variableResult, expected_result);
                BOOST_CHECK_EQUAL(constantResult, expected_result);

                m.Transpose();
                variableResult = m*v_transposed_variable;
                constantResult = m*v_transposed_constant;
                BOOST_CHECK_EQUAL(variableResult, expected_transpose_result);
                BOOST_CHECK_EQUAL(constantResult, expected_transpose_result);
                
                m.Transpose();
                variableResult = m*v_non_transposed_variable;
                constantResult = m*v_non_transposed_constant;
                BOOST_CHECK_EQUAL(variableResult, expected_result);
                BOOST_CHECK_EQUAL(constantResult, expected_result);

                NekMatrix<double> transposed = Transpose(m);
                variableResult = transposed*v_transposed_variable;
                constantResult = transposed*v_transposed_constant;
                BOOST_CHECK_EQUAL(variableResult, expected_transpose_result);
                BOOST_CHECK_EQUAL(constantResult, expected_transpose_result);
            }
           
        }

        BOOST_AUTO_TEST_CASE(TestVectorConstructorsWithSizeArguments)
        {
            {
                double buf[] = {1.0, 2.0, 3.0, 4.0};
                Array<OneD, double> a(4, buf);
                NekVector<double> b(3, a);

                BOOST_CHECK_EQUAL(b.GetRows(), 3u);
                BOOST_CHECK_EQUAL(b.GetDimension(), 3u);
                BOOST_CHECK_EQUAL(b[0], 1.0);
                BOOST_CHECK_EQUAL(b[1], 2.0);
                BOOST_CHECK_EQUAL(b[2], 3.0);

                NekVector<double>::iterator iter = b.begin();
                BOOST_CHECK_EQUAL(*iter, 1.0);
                ++iter;
                BOOST_CHECK_EQUAL(*iter, 2.0);
                ++iter;
                BOOST_CHECK_EQUAL(*iter, 3.0);
                ++iter;
                BOOST_CHECK(iter == b.end());
            }

            {
                double buf[] = {1.0, 2.0, 3.0, 4.0};
                Array<OneD, const double> a(4, buf);
                NekVector<double> b(3, a);

                BOOST_CHECK_EQUAL(b.GetRows(), 3u);
                BOOST_CHECK_EQUAL(b.GetDimension(), 3u);
                BOOST_CHECK_EQUAL(b[0], 1.0);
                BOOST_CHECK_EQUAL(b[1], 2.0);
                BOOST_CHECK_EQUAL(b[2], 3.0);

                NekVector<double>::iterator iter = b.begin();
                BOOST_CHECK_EQUAL(*iter, 1.0);
                ++iter;
                BOOST_CHECK_EQUAL(*iter, 2.0);
                ++iter;
                BOOST_CHECK_EQUAL(*iter, 3.0);
                ++iter;
                BOOST_CHECK(iter == b.end());
            }

            {
                double buf[] = {1.0, 2.0, 3.0, 4.0};
                Array<OneD, double> a(4, buf);
                NekVector<double> b(3, a, eWrapper);

                BOOST_CHECK_EQUAL(b.GetRows(), 3u);
                BOOST_CHECK_EQUAL(b.GetDimension(), 3u);
                BOOST_CHECK_EQUAL(b[0], 1.0);
                BOOST_CHECK_EQUAL(b[1], 2.0);
                BOOST_CHECK_EQUAL(b[2], 3.0);

                NekVector<double>::iterator iter = b.begin();
                BOOST_CHECK_EQUAL(*iter, 1.0);
                ++iter;
                BOOST_CHECK_EQUAL(*iter, 2.0);
                ++iter;
                BOOST_CHECK_EQUAL(*iter, 3.0);
                ++iter;
                BOOST_CHECK(iter == b.end());

                b[0] = 5.0;
                BOOST_CHECK_EQUAL(b[0], a[0]);
            }
        }

    }
}

