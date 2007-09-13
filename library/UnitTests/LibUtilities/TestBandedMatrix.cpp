///////////////////////////////////////////////////////////////////////////////
//
// File: TestBandedMatrix.cpp
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

#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

namespace Nektar
{
    namespace BandedMatrixPolicyUnitTests
    {

    }

    namespace BandedMatrixMultiplicationUnitTests
    {
        BOOST_AUTO_TEST_CASE(TestDirectBlasCall)
        {
            // [1 2 0 0 ]
            // [3 4 5 0 ]
            // [0 6 7 8 ]
            // [0 0 9 10]
            NekDouble x[] = {1, 2, 3, 4};
            NekDouble y[] = {0,0,0,0};

            NekDouble a[] = {0, 1, 2,
                             3, 4, 5,
                             6, 7, 8, 
                             9, 10, 0};
            Blas::Dgbmv('T', 4, 4, 1, 1, 1.0, a, 3, x, 1, 0.0, y, 1);

            NekDouble expected_result_buf[] = { 5, 26, 65, 67 };
            NekVector<NekDouble, 4> expected_result(expected_result_buf);
            NekVector<NekDouble, 4> result(y);
            BOOST_CHECK_EQUAL(expected_result, result);
        }
        
        BOOST_AUTO_TEST_CASE(TestSquareDirectBlasCall)
        {
            // [1  2  3  4 ]
            // [5  6  7  8 ]
            // [9  10 11 12]
            // [13 14 15 16]
            NekDouble x[] = {1, 2, 3, 4};
            NekDouble y[] = {0,0,0,0};

            NekDouble a[] = {0, 0, 0, 1, 2, 3, 4,
                             0, 0, 5, 6, 7, 8, 0, 
                             0, 9, 10, 11, 12, 0, 0,
                             13, 14, 15, 16, 0, 0, 0 };
            Blas::Dgbmv('T', 4, 4, 3, 3, 1.0, a, 7, x, 1, 0.0, y, 1);


            NekDouble expected_result_buf[] = { 30, 70, 110, 150 };
            NekVector<NekDouble, 4> expected_result(expected_result_buf);
            NekVector<NekDouble, 4> result(y);
            BOOST_CHECK_EQUAL(expected_result, result);
        }

        BOOST_AUTO_TEST_CASE(TestDifferentNumberOfLowerAndUpperDiagonalsDirectBlasCall)
        {
            // [1 2 11 0 ]
            // [3 4 5  12 ]
            // [0 6 7  8 ]
            // [0 0 9  10]
            NekDouble x[] = {1, 2, 3, 4};
            NekDouble y[] = {0,0,0,0};

            NekDouble a[] = {0, 1, 2, 11,
                             3, 4, 5, 12,
                             6, 7, 8, 0,
                             9, 10, 0, 0};

            Blas::Dgbmv('T', 4, 4, 2, 1, 1.0, a, 4, x, 1, 0.0, y, 1);

            NekDouble expected_result_buf[] = { 38, 74, 65, 67 };
            NekVector<NekDouble, 4> expected_result(expected_result_buf);
            NekVector<NekDouble, 4> result(y);
            BOOST_CHECK_EQUAL(expected_result, result);
        }

        BOOST_AUTO_TEST_CASE(TestLargerPackedMatrixThanNormal)
        {
            // [ 1 2 0 ]
            // [ 3 4 5 ]
            // [ 6 7 8 ]
            NekDouble x[] = {1, 2, 3};
            NekDouble y[] = {0, 0, 0};

            NekDouble a[] = {0, 1, 3, 6,
                             2, 4, 7, 0,
                             5, 8, 0, 0};

            Blas::Dgbmv('N', 3, 3, 2, 1, 1.0, a, 4, x, 1, 0.0, y, 1);

            NekDouble expected_result_buf[] = { 5, 26, 44 };
            NekVector<NekDouble, 3> expected_result(expected_result_buf);
            NekVector<NekDouble, 3> result(y);
            BOOST_CHECK_EQUAL(expected_result, result);
        }

        BOOST_AUTO_TEST_CASE(TestStandardMatrixVectorMultiply)
        {
            
            {
                typedef NekMatrix<NekDouble, BandedMatrixTag> MatrixType;

                // This is the matrix
                // [ 1 2 0 0 ]
                // [ 3 4 5 0 ]
                // [ 0 6 7 8 ]
                // [ 0 0 9 10]
                
                NekDouble buf[] = {0, 1, 3,
                                   2, 4, 6,
                                   5, 7, 9,
                                   8, 10, 0};
                              
    
                MatrixType m(4, 4, buf, MatrixType::PolicySpecificDataHolderType(1, 1));

                NekDouble vector_buf[] = {1.0, 2.0, 3.0, 4.0};
                NekVector<NekDouble> v(4, vector_buf);

                NekVector<NekDouble> result = m*v;

                NekDouble expected_result_buf[] = { 5, 26, 65, 67 };
                NekVector<NekDouble> expected_result(4, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }

            {
                typedef NekMatrix<unsigned int, BandedMatrixTag> MatrixType;

                // This is the matrix
                // [ 1 2 0 0 ]
                // [ 3 4 5 0 ]
                // [ 0 6 7 8 ]
                // [ 0 0 9 10]
                
                unsigned int buf[] = {0, 1, 3,
                                      2, 4, 6,
                                      5, 7, 9,
                                      8, 10, 0};
                              
    
                MatrixType m(4, 4, buf, MatrixType::PolicySpecificDataHolderType(1, 1));

                unsigned int vector_buf[] = {1, 2, 3, 4};
                NekVector<unsigned int> v(4, vector_buf);

                NekVector<unsigned int> result = m*v;

                unsigned int expected_result_buf[] = { 5, 26, 65, 67 };
                NekVector<unsigned int> expected_result(4, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }
        }
    }
}


