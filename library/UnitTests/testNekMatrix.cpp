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
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

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

                NekMatrix<double> dynamic_matrix(4, 3, buf);

                BOOST_CHECK(dynamic_matrix.GetRows() == 4);
                BOOST_CHECK(dynamic_matrix.GetColumns() == 3);

                for(unsigned int i = 0; i < 4; ++i)
                {
                    for(unsigned int j = 0; j < 3; ++j)
                    {
                        BOOST_CHECK(dynamic_matrix(i,j) == buf[3*i + j]);
                    }
                }
            }

            {

                NekMatrix<float> dynamic_matrix(7, 3, 7.8);

                for(unsigned int i = 0; i < 7; ++i)
                {
                    for(unsigned int j = 0; j < 3; ++j)
                    {
                        BOOST_CHECK(dynamic_matrix(i,j) == 7.8f);
                    }
                }
            }

        }

        void testNekMatrixAccess()
        {
            // We need to be able to access any element in the matrix, and
            // assign into the matrix at any location.
            NekMatrix<unsigned int> static_matrix(3,3);

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
                    BOOST_CHECK(static_matrix(i,j) == 10*i + j);
                }
            }

            // Invalid access is an unrecoverable error.
            //BOOST_CHECK_THROW(static_matrix(3,2), OutOfBoundsError);
            //BOOST_CHECK_THROW(static_matrix(2,3), OutOfBoundsError);
            BOOST_CHECK_NO_THROW(static_matrix(2,2));

        }

        void testNekMatrixBasicMath()
        {
            // Addition tests.
            {
                double buf[] = {1.0, 2.0, 3.0,
                    4.0, 5.0, 6.0,
                    7.0, 8.0, 9.0 };

                NekMatrix<double> m1(3, 3, buf);
                NekMatrix<double> m2(3, 3, buf);
                NekMatrix<double> m3 = m1 + m2;

                for(unsigned int i = 0; i < 3; ++i)
                {
                    for(unsigned int j = 0; j < 3; ++j)
                    {
                        BOOST_CHECK(m3(i,j) == buf[3*i+j] + buf[3*i+j]);
                    }
                }

                NekMatrix<double> m4(3, 3, buf);
                NekMatrix<double> m5(3, 3, buf);
                NekMatrix<double> m6 = m4+m5;

                for(unsigned int i = 0; i < 3; ++i)
                {
                    for(unsigned int j = 0; j < 3; ++j)
                    {
                        BOOST_CHECK(m6(i,j) == buf[3*i+j] + buf[3*i+j]);
                    }
                }
            }
            // Multiply
            {
                unsigned int buf1[] = {1, 2, 3,
                                       4, 5, 6,
                                       7, 8, 9};
                unsigned int buf2[] = { 10, 11, 12, 14,
                                        15, 16, 17, 18,
                                        19, 20, 21, 22 };

                NekMatrix<unsigned int> lhs(3, 3, buf1);
                NekMatrix<unsigned int> rhs(3, 4, buf2);
                NekMatrix<unsigned int> result = lhs*rhs;

                BOOST_CHECK(result.GetRows() == 3);
                BOOST_CHECK(result.GetColumns() == 4);

                BOOST_CHECK(result(0,0) == 97);
                BOOST_CHECK(result(0,1) == 103);
                BOOST_CHECK(result(0,2) == 109);
                BOOST_CHECK(result(0,3) == 116);

                BOOST_CHECK(result(1,0) == 229);
                BOOST_CHECK(result(1,1) == 244);
                BOOST_CHECK(result(1,2) == 259);
                BOOST_CHECK(result(1,3) == 278);

                BOOST_CHECK(result(2,0) == 361);
                BOOST_CHECK(result(2,1) == 385);
                BOOST_CHECK(result(2,2) == 409);
                BOOST_CHECK(result(2,3) == 440);
            }

            {
                unsigned int buf1[] = {1, 2, 3,
                                       4, 5, 6,
                                       7, 8, 9};
                unsigned int buf2[] = { 1, 2, 3};

                NekMatrix<unsigned int> lhs(3, 3, buf1);
                NekVector<unsigned int> rhs(3, buf2);

                NekVector<unsigned int> result = lhs*rhs;

                BOOST_CHECK(result[0] == 14);
                BOOST_CHECK(result[1] == 32);
                BOOST_CHECK(result[2] == 50);
            }

            // Transpose

            // Determinant.

            // Invert/Check for singularity.

            // Eigenvalues/vectors

            // Condition number wrt various norms.

            // Various norm computations.

            // LU Decomposition?  More appropriate in LinAlg?
        }

        void testDiagonalMatrix()
        {
        }
    }
}


/**
    $Log: testNekMatrix.cpp,v $
    Revision 1.10  2006/05/31 23:24:47  bnelson
    Updated NekMatrix method names for the coding standard.

    Revision 1.9  2006/05/31 04:19:36  bnelson
    Removed a test for invalid access to a matrix.

    Revision 1.8  2006/05/29 03:40:48  bnelson
    Updated the tests to reflect the changed parameter order in the NekMatrix constructor.

    Revision 1.7  2006/05/25 02:54:54  bnelson
    Added Matrix/Vector multiplication test.

    Revision 1.6  2006/05/18 04:25:19  bnelson
    Added a multiplication test.

    Revision 1.5  2006/05/16 20:35:30  jfrazier
    Added the float literal specifier to make the unit test happy.

    Revision 1.4  2006/05/15 05:06:07  bnelson
    Added addition tests.

    Revision 1.3  2006/05/15 04:10:35  bnelson
    no message

    Revision 1.2  2006/05/14 21:33:58  bnelson
    *** empty log message ***

    Revision 1.1  2006/05/07 21:10:09  bnelson
    *** empty log message ***

 **/

