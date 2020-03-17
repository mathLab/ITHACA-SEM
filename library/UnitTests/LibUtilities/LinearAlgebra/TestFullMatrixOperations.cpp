///////////////////////////////////////////////////////////////////////////////
//
// File: TestFullMatrixOperations.cpp
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

#include <boost/test/unit_test.hpp>

#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <UnitTests/CountedObject.h>

#include <boost/core/ignore_unused.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

namespace Nektar
{
    namespace FullMatrixOperationsUnitTests
    {
        template<typename DataType, typename MatrixType>
        int foo(NekMatrix<DataType, MatrixType>& d)
        {
            boost::ignore_unused(d);
            return 1;
        }
        
        template<typename DataType>
        int foo(NekMatrix<DataType, BlockMatrixTag>& d)
        {
            boost::ignore_unused(d);
            return 2;
        }
            
        BOOST_AUTO_TEST_CASE(TestDoubleSquareFullVectorMultiplication)
        {
            NekMatrix<double> m1;
            NekMatrix<NekMatrix<double>, BlockMatrixTag> m2(2,2,2,2);
            foo(m1);
            foo(m2);
            double m_buf[] = {1, 2, 3,
                              4, 5, 6,
                              7, 8, 9};
            double v_buf[] = {4, 5, 6};
            
            NekMatrix<double> m(3, 3, m_buf);
            NekVector<double> v(3, v_buf);
            
            double expected_result_buf[] = {66, 81, 96};
            NekVector<double> expected_result(3, expected_result_buf);
            
            NekVector<double> result = m*v;
            
            BOOST_CHECK_EQUAL(expected_result, result);
        }
        
        BOOST_AUTO_TEST_CASE(TestDoubleSquareFullVectorMultiplicationWithAliasing)
        {
            double m_buf[] = {1, 2, 3,
                              4, 5, 6,
                              7, 8, 9};
            double v_buf[] = {4, 5, 6};
            
            NekMatrix<double> m(3, 3, m_buf);
            NekVector<double> v(3, v_buf);
            
            double expected_result_buf[] = {66, 81, 96};
            NekVector<double> expected_result(3, expected_result_buf);
            
            v = m*v;
            
            BOOST_CHECK_EQUAL(expected_result, v);
        }
        
        BOOST_AUTO_TEST_CASE(TestDoubleSquareFullVectorMultiplicationWithSharedArrayAliasing)
        {
            double m_buf[] = {1, 2, 3,
                              4, 5, 6,
                              7, 8, 9};
            double v_buf[] = {4, 5, 6};
            
            Array<OneD, double> vector_array(3, v_buf);
            
            NekMatrix<double> m(3, 3, m_buf);
            NekVector<double> v(3, vector_array, eWrapper);
            NekVector<double> result(3, vector_array, eWrapper);
            
            double expected_result_buf[] = {66, 81, 96};
            NekVector<double> expected_result(3, expected_result_buf);
            
            result = m*v;
            
            BOOST_CHECK_EQUAL(expected_result, result);
            BOOST_CHECK_EQUAL(vector_array[0], 66);
            BOOST_CHECK_EQUAL(vector_array[1], 81);
            BOOST_CHECK_EQUAL(vector_array[2], 96);
        }

        BOOST_AUTO_TEST_CASE(TestDoubleSquareFullVectorMultiplicationWithSharedArrayAliasingAndMatrixTranspose)
        {
            double m_buf[] = {1, 4, 7,
                              2, 5, 8,
                              3, 6, 9};
            double v_buf[] = {4, 5, 6};
            
            Array<OneD, double> vector_array(3, v_buf);
            
            NekMatrix<double> m(3, 3, m_buf);
            NekVector<double> v(3, vector_array, eWrapper);
            NekVector<double> result(3, vector_array, eWrapper);
            
            double expected_result_buf[] = {66, 81, 96};
            NekVector<double> expected_result(3, expected_result_buf);
            
            result = Transpose(m)*v;
            
            BOOST_CHECK_EQUAL(expected_result, result);
            BOOST_CHECK_EQUAL(vector_array[0], 66);
            BOOST_CHECK_EQUAL(vector_array[1], 81);
            BOOST_CHECK_EQUAL(vector_array[2], 96);
        }
        
        BOOST_AUTO_TEST_CASE(TestDoubleSquareFullVectorMultiplicationWithSharedArrayAliasingAndOverlap)
        {
            double m_buf[] = {1, 2, 3,
                              4, 5, 6,
                              7, 8, 9};
            double v_buf[] = {4, 5, 6, 0};
            
            Array<OneD, double> vector_array(4, v_buf);
            Array<OneD, double> tmp;

            NekMatrix<double> m(3, 3, m_buf);
            NekVector<double> v(3, vector_array, eWrapper);
            NekVector<double> result(3, tmp=vector_array+1, eWrapper);
            
            double expected_result_buf[] = {66, 81, 96};
            NekVector<double> expected_result(3, expected_result_buf);
            
            result = m*v;
            
            BOOST_CHECK_EQUAL(expected_result, result);
            BOOST_CHECK_EQUAL(vector_array[1], 66);
            BOOST_CHECK_EQUAL(vector_array[2], 81);
            BOOST_CHECK_EQUAL(vector_array[3], 96);
        }


        
        BOOST_AUTO_TEST_CASE(TestScaledDoubleSquareFullVectorMultiplication)
        {
            double m_buf[] = {1, 2, 3,
                              4, 5, 6,
                              7, 8, 9};
            double v_buf[] = {4, 5, 6};
            
            std::shared_ptr<NekMatrix<double> > inner(new NekMatrix<double>(3, 3, m_buf));
            NekMatrix<NekMatrix<double>, ScaledMatrixTag> m(7.0, inner);
            NekVector<double> v(3, v_buf);
            
            double expected_result_buf[] = {462, 567, 672};
            NekVector<double> expected_result(3, expected_result_buf);
            
            NekVector<double> result = m*v;
            
            BOOST_CHECK_EQUAL(expected_result, result);
        }


        BOOST_AUTO_TEST_CASE(TestThreeMatrixMultiplication)
        {
            double buf1[] = {1.0, 2.0,
                             3.0, 4.0};
            double buf2[] = {5.0, 6.0,
                             7.0, 8.0};
            double buf3[] = {9.0, 10.0,
                             11.0, 12.0};

            NekMatrix<double> m1(2, 2, buf1);
            NekMatrix<double> m2(2, 2, buf2);
            NekMatrix<double> m3(2, 2, buf3);

            double expected_result_buf[] = {517, 766,
                                            625, 926};
            NekMatrix<double> expected_result(2, 2, expected_result_buf);

            NekMatrix<double> result1 = m1*m2*m3;
            NekMatrix<double> result2 = m1*(m2*m3);
            BOOST_CHECK_EQUAL(expected_result, m1*m2*m3);

            BOOST_CHECK_EQUAL(expected_result, m1*(m2*m3));                                                    
        }

        BOOST_AUTO_TEST_CASE(TestThreeMatrixMultiplicationWithTranspose)
        {
            double buf1[] = {1.0, 2.0,
                             3.0, 4.0};
            double buf2[] = {5.0, 6.0,
                             7.0, 8.0};
            double buf3[] = {9.0, 11.0,
                             10.0, 12.0};

            NekMatrix<double> m1(2, 2, buf1);
            NekMatrix<double> m2(2, 2, buf2);
            NekMatrix<double> m3(2, 2, buf3);

            double expected_result_buf[] = {517, 766,
                                            625, 926};
            NekMatrix<double> expected_result(2, 2, expected_result_buf);

            BOOST_CHECK_EQUAL(expected_result, m1*m2*Transpose(m3));
        }

        
        BOOST_AUTO_TEST_CASE(TestThreeWrappedMatrixMultiplicationWithTranspose)
        {
            double buf1[] = {1.0, 2.0,
                             3.0, 4.0};
            double buf2[] = {5.0, 6.0,
                             7.0, 8.0};
            double buf3[] = {9.0, 11.0,
                             10.0, 12.0};
            double out_buf[] = {0.0, 0.0, 0.0, 0.0};

            Array<OneD, double> array_buf1(4, buf1);
            Array<OneD, double> array_buf2(4, buf2);
            Array<OneD, double> array_buf3(4, buf3);
            Array<OneD, double> array_out_buf(4, out_buf);

            NekMatrix<double> m1(2, 2, array_buf1, eWrapper);
            NekMatrix<double> m2(2, 2, array_buf2, eWrapper);
            NekMatrix<double> m3(2, 2, array_buf3, eWrapper);
            NekMatrix<double> result(2, 2, array_out_buf, eWrapper);

            double expected_result_buf[] = {517, 766,
                                            625, 926};
            NekMatrix<double> expected_result(2, 2, expected_result_buf);

            result = m1*m2*Transpose(m3);
            BOOST_CHECK_EQUAL(expected_result, result);
            BOOST_CHECK_EQUAL(array_out_buf[0], 517.0);
            BOOST_CHECK_EQUAL(array_out_buf[1], 766.0);
            BOOST_CHECK_EQUAL(array_out_buf[2], 625.0);
            BOOST_CHECK_EQUAL(array_out_buf[3], 926.0);
            BOOST_CHECK_EQUAL(result(0,0), 517.0);
            BOOST_CHECK_EQUAL(result(1,0), 766.0);
            BOOST_CHECK_EQUAL(result(0,1), 625.0);
            BOOST_CHECK_EQUAL(result(1,1), 926.0);
        }
    }
}
