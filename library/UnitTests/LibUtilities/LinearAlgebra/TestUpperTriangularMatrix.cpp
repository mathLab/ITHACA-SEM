///////////////////////////////////////////////////////////////////////////////
//
// File: TestUpperTriangularMatrix.cpp
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

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

namespace Nektar
{
    namespace UpperTriangularMatrixUnitTests
    {
        typedef UpperTriangularMatrixFuncs Policy;

        BOOST_AUTO_TEST_CASE(TestMatrixVectorMultiplyUpper)
        {
            {
                double matrix_buf[] = {1, 2,
                                          3};
                NekMatrix<double> matrix(2,2,matrix_buf,eUPPER_TRIANGULAR);

                double vector_buf[] = {10, 11};
                NekVector<double> vector(2, vector_buf);

                NekVector<double> result = matrix*vector;

                double expected_buf[] = {32, 33};
                NekVector<double> expected_result(2, expected_buf);

                BOOST_CHECK_EQUAL(expected_result, result);
            }
        }

        BOOST_AUTO_TEST_CASE(Test3x3MatrixVectorMultiplyUpper)
        {
            {
                //double matrix_buf[] = {1, 2, 3,
                //                          4, 5,
                //                             6};
                double matrix_buf[] = {1,
                                       2, 4,
                                       3, 5, 6};
                NekMatrix<double> matrix(3,3,matrix_buf,eUPPER_TRIANGULAR);

                double vector_buf[] = {10, 11, 12};
                NekVector<double> vector(3, vector_buf);

                NekVector<double> result = matrix*vector;

                double expected_buf[] = {68, 104, 72};
                NekVector<double> expected_result(3, expected_buf);

                BOOST_CHECK_EQUAL(expected_result, result);
            }

        }
    }
}


