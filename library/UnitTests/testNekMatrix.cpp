///////////////////////////////////////////////////////////////////////////////
//
// File: testNekMatrix.cpp
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
// Description: Tests NekMatrix functionality.
//
///////////////////////////////////////////////////////////////////////////////

#include <UnitTests/testNekMatrix.h>
#include <LibUtilties/nekMatrix.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>

namespace Nektar
{
    namespace UnitTests
    {
        using namespace Nektar::LibUtilities;

        void testNekMatrixConstruction()
        {
            // Basic, dense matrix construction.
            {
                double buf[] = { 1.0, 2.0, 3.0,
                                4.0, 5.0, 6.0,
                                7.0, 8.0, 9.0,
                                10.0, 11.0, 12.0 };

                NekMatrix<double, 4, 3> static_matrix(buf);
                NekMatrix<double> dynamic_matrix(buf, 4, 3);

                BOOST_CHECK(static_matrix.rows() == 4);
                BOOST_CHECK(static_matrix.columns() == 3);
                BOOST_CHECK(dynamic_matrix.rows() == 4);
                BOOST_CHECK(dynamic_matrix.columns() == 3);

                for(unsigned int i = 0; i < 4; ++i)
                {
                    for(unsigned int j = 0; j < 3; ++j)
                    {
                        BOOST_CHECK(static_matrix(i,j) == buf[4*i + j]);
                        BOOST_CHECK(dynamic_matrix(i,j) == buf[4*i + j]);
                    }
                }
            }

            {
                NekMatrix<float, 7, 3> static_matrix(7.8);
                NekMatrix<float> dynamic_matrix(7.8, 7, 3);

                for(unsigned int i = 0; i < 7; ++i)
                {
                    for(unsigned int j = 0; j < 3; ++j)
                    {
                        BOOST_CHECK(static_matrix(i,j) == 7.8);
                        BOOST_CHECK(dynamic_matrix(i,j) == 7.8);
                    }
                }
            }

        }

        void testNekMatrixAccess()
        {
            // We need to be able to access any element in the matrix, and
            // assign into the matrix at any location.
            NekMatrix<unsigned int, 3, 3> static_matrix;

            // Test read access.
            for(unsigned int i = 0; i < 3; ++i)
            {
                for(unsigned int j = 0; j < 3; ++j)
                {
                    BOOST_CHECK(static_matrix(i,j) == 0.0);
                }
            }

            for(unsigned int i = 0; i < 3; ++i)
            {
                for(unsigned int j = 0; j < 3; ++j)
                {
                    static_matrix(i,j) = 10*i + j;
                }
            }

            for(unsigned int i = 0; i < 3; ++i)
            {
                for(unsigned int j = 0; j < 3; ++j)
                {
                    static_matrix(i,j) == 10*i + j;
                }
            }

            // Invalid access is an unrecoverable error.
            BOOST_CHECK_THROW(static_matrix(3,2), NekMatrix::OutOfBoundsError);
            BOOST_CHECK_THROW(static_matrix(2,3), NekMatrix::OutOfBoundsError);
            BOOST_CHECK_NO_THROW(static_matrix(2,2));

        }

        void testNekMatrixBasicMath()
        {
            // Addition tests.
            {
                double buf[] = {1.0, 2.0, 3.0,
                    4.0, 5.0, 6.0,
                    7.0, 8.0, 9.0 };

                NekMatrix<double, 3, 3> m1(buf);
                NekMatrix<double, 3, 3> m2(buf);
                NekMatrix<double, 3, 3> m3 = m1 + m2;

                for(unsigned int i = 0; i < 3; ++i)
                {
                    for(unsigned int j = 0; j < 3; ++j)
                    {
                        BOOST_CHECK(m3(i,j) == buf[3*i+j] + buf[3*i+j]);
                    }
                }

                NekMatrix<double> m4(buf, 3, 3);
                NekMatrix<double> m5(buf, 3, 3);
                NekMatrix<double> m6 = m4+m5;

                for(unsigned int i = 0; i < 3; ++i)
                {
                    for(unsigned int j = 0; j < 3; ++j)
                    {
                        BOOST_CHECK(m6(i,j) == buf[3*i+j] + buf[3*i+j]);
                    }
                }

                // Do a couple of tests that shouldn't compile.
                NekMatrix<double, 3, 3> m7(buf);
                NekMatrix<double, 2, 2> m8(buf);
                NekMatrix<double> m9 = m7 + m8; // This line should fail.
                NekMatrix<double, 3, 3> m10 = m7 + m8; // This line should fail.

                // Mixed mode.
                NekMatrix<double> m11 = m7 + m4;
                NekMatrix<double, 3, 3> m12 = m7 + m4;
                BOOST_CHECK(m11 == m12);
                BOOST_CHECK(m11 == m3);
            }

            // Multiply

            // Transpose

            // Determinant.

            // Invert/Check for singularity.

            // Eigenvalues/vectors

            // Condition number wrt various norms.

            // Various norm computations.

            // LU Decomposition?  More appropriate in LinAlg?
        }
    }
}


/**
    $Log: $
 **/

