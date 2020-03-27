///////////////////////////////////////////////////////////////////////////////
//
// File: TestTriangularMatrixOperations.cpp
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

#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekLinSys.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

namespace Nektar 
{
    namespace TriangularMatrixVectorMultiplicationUnitTests
    {
        BOOST_AUTO_TEST_CASE(TestUpperTriangularMatrixVectorMultiplication)
        {
            // [1 2 3 4]
            // [0 5 6 7]
            // [0 0 8 9]
            // [0 0 0 10]

            NekDouble a_buf[] = {1, 2, 5, 3, 6, 8, 4, 7, 9, 10};
            NekDouble x_buf[] = {10, 20, 30, 40};

            NekMatrix<NekDouble, StandardMatrixTag> m(4, 4, a_buf, eUPPER_TRIANGULAR);
            NekVector<NekDouble> x(4, x_buf);

            NekVector<NekDouble> result = m*x;

            NekDouble expected_result_buf[] = { 300, 560, 600, 400 };
            NekVector<NekDouble> expected_result(4, expected_result_buf);
            BOOST_CHECK_EQUAL(expected_result, result);
        }
                
        BOOST_AUTO_TEST_CASE(TestScaledUpperTriangularMatrixVectorMultiplication)
        {
            // [1 2 3 4]
            // [0 5 6 7]
            // [0 0 8 9]
            // [0 0 0 10]

            NekDouble a_buf[] = {1, 2, 5, 3, 6, 8, 4, 7, 9, 10};
            NekDouble x_buf[] = {10, 20, 30, 40};

            std::shared_ptr<NekMatrix<NekDouble, StandardMatrixTag> > m(
                new NekMatrix<NekDouble>(4, 4, a_buf, eUPPER_TRIANGULAR));
            NekMatrix<NekMatrix<NekDouble>, ScaledMatrixTag> scaled(2, m);
            NekVector<NekDouble> x(4, x_buf);

            NekVector<NekDouble> result = scaled*x;

            NekDouble expected_result_buf[] = { 600, 1120, 1200, 800 };
            NekVector<NekDouble> expected_result(4, expected_result_buf);
            BOOST_CHECK_EQUAL(expected_result, result);
        }

        BOOST_AUTO_TEST_CASE(TestLowerTriangularMatrixVectorMultiplication)
        {
            // [1 0 0 0]
            // [2 6 0 0]
            // [5 8 7 0]
            // [3 4 9 10]

            NekDouble a_buf[] = {1, 2, 5, 3, 6, 8, 4, 7, 9, 10};
            NekDouble x_buf[] = {10, 20, 30, 40};

            NekMatrix<NekDouble, StandardMatrixTag> m(4, 4, a_buf, eLOWER_TRIANGULAR);
            NekVector<NekDouble> x(4, x_buf);

            NekVector<NekDouble> result = m*x;

            NekDouble expected_result_buf[] = { 10, 140, 420, 780 };
            NekVector<NekDouble> expected_result(4, expected_result_buf);
            BOOST_CHECK_EQUAL(expected_result, result);
        }
        
        BOOST_AUTO_TEST_CASE(TestScaledLowerTriangularMatrixVectorMultiplication)
        {
            // [1 0 0 0]
            // [2 6 0 0]
            // [5 8 7 0]
            // [3 4 9 10]

            NekDouble a_buf[] = {1, 2, 5, 3, 6, 8, 4, 7, 9, 10};
            NekDouble x_buf[] = {10, 20, 30, 40};

            std::shared_ptr<NekMatrix<NekDouble, StandardMatrixTag> > m(
                new NekMatrix<NekDouble>(4, 4, a_buf, eLOWER_TRIANGULAR));
            NekMatrix<NekMatrix<NekDouble>, ScaledMatrixTag > scaled(3.0, m);
            NekVector<NekDouble> x(4, x_buf);

            NekVector<NekDouble> result = scaled*x;

            NekDouble expected_result_buf[] = { 30, 420, 1260, 2340 };
            NekVector<NekDouble> expected_result(4, expected_result_buf);
            BOOST_CHECK_EQUAL(expected_result, result);
        }


        BOOST_AUTO_TEST_CASE(TestUpperTriangularSolve)
        {
            // [1 2 3 4]
            // [0 5 6 7]
            // [0 0 8 9]
            // [0 0 0 10]

            NekDouble a_buf[] = {1, 2, 5, 3, 6, 8, 4, 7, 9, 10};
            NekMatrix<NekDouble, StandardMatrixTag> m(4, 4, a_buf, eUPPER_TRIANGULAR);

            NekDouble b_buf[] = {300, 560, 600, 400};
            NekVector<NekDouble> b(4, b_buf);

            LinearSystem sys(m);
            NekVector<NekDouble> x = sys.Solve(b);

            NekDouble expected_result_buf[] = {10, 20, 30, 40};
            NekVector<NekDouble> expected_result(4, expected_result_buf);

            BOOST_CHECK_EQUAL(expected_result, x);
        }

        BOOST_AUTO_TEST_CASE(TestUpperTriangularTransposeSolve)
        {
            // [1 2 3 4]
            // [0 5 6 7]
            // [0 0 8 9]
            // [0 0 0 10]

            NekDouble a_buf[] = {1, 2, 5, 3, 6, 8, 4, 7, 9, 10};
            NekMatrix<NekDouble, StandardMatrixTag> m(4, 4, a_buf, eUPPER_TRIANGULAR);

            NekDouble b_buf[] = {10, 120, 390, 850};
            NekVector<NekDouble> b(4, b_buf);

            LinearSystem sys(m);
            NekVector<NekDouble> x = sys.SolveTranspose(b);

            NekDouble expected_result_buf[] = {10, 20, 30, 40};
            NekVector<NekDouble> expected_result(4, expected_result_buf);

            BOOST_CHECK_EQUAL(expected_result, x);
        }

        BOOST_AUTO_TEST_CASE(TestLowerTriangularSolve)
        {
            // [1 2 3 4]
            // [0 5 6 7]
            // [0 0 8 9]
            // [0 0 0 10]

            NekDouble a_buf[] = {1, 2, 5, 3, 6, 8, 4, 7, 9, 10};
            NekMatrix<NekDouble, StandardMatrixTag> m(4, 4, a_buf, eLOWER_TRIANGULAR);

            NekDouble b_buf[] = {10, 140, 420, 780};
            NekVector<NekDouble> b(4, b_buf);

            LinearSystem sys(m);
            NekVector<NekDouble> x = sys.Solve(b);

            NekDouble expected_result_buf[] = {10, 20, 30, 40};
            NekVector<NekDouble> expected_result(4, expected_result_buf);

            BOOST_CHECK_EQUAL(expected_result, x);
        }

        BOOST_AUTO_TEST_CASE(TestLowerTriangularTransposeSolve)
        {
            // [1 2 3 4]
            // [0 5 6 7]
            // [0 0 8 9]
            // [0 0 0 10]

            NekDouble a_buf[] = {1, 2, 5, 3, 6, 8, 4, 7, 9, 10};
            NekMatrix<NekDouble, StandardMatrixTag> m(4, 4, a_buf, eLOWER_TRIANGULAR);

            NekDouble b_buf[] = {320, 520, 570, 400};
            NekVector<NekDouble> b(4, b_buf);

            LinearSystem sys(m);
            NekVector<NekDouble> x = sys.SolveTranspose(b);

            NekDouble expected_result_buf[] = {10, 20, 30, 40};
            NekVector<NekDouble> expected_result(4, expected_result_buf);

            BOOST_CHECK_EQUAL(expected_result, x);
        }
    }
}


