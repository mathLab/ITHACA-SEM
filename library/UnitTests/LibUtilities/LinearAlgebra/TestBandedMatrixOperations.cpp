///////////////////////////////////////////////////////////////////////////////
//
// File: TestBandedMatrixOperations.cpp
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
#include <UnitTests/LibUtilities/LinearAlgebra/TestCombinationRunner.h>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

namespace Nektar 
{
    namespace BandedMatrixVectorMultiplicationUnitTests
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

        BOOST_AUTO_TEST_CASE(TestUnequalNumbersOfSubAndSuperDiagonalsMatrixVectorMultiply)
        {
            {
                typedef NekMatrix<NekDouble, BandedMatrixTag> MatrixType;

                // This is the matrix
                // [  1  6 10   0  0 ]
                // [ 13  2  7  11  0 ]
                // [  0 14  3   8 12 ]
                // [  0  0 15   4  9 ]
                // [  0  0  0  16  5 ]
                
                NekDouble buf[] = {0, 0, 1, 13,
                                   0, 6, 2, 14,
                                   10, 7, 3, 15,
                                   11, 8, 4, 16,
                                   12, 9, 5, 0};
                              
    
                MatrixType m(5, 5, buf, MatrixType::PolicySpecificDataHolderType(1, 2));

                NekDouble vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<NekDouble> v(5, vector_buf);

                NekVector<NekDouble> result = m*v;

                NekDouble expected_result_buf[] = { 43, 82, 129, 106, 89 };
                NekVector<NekDouble> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }

            {
                typedef NekMatrix<unsigned int, BandedMatrixTag> MatrixType;

                // This is the matrix
                // [  1  6 10   0  0 ]
                // [ 13  2  7  11  0 ]
                // [  0 14  3   8 12 ]
                // [  0  0 15   4  9 ]
                // [  0  0  0  16  5 ]
                
                unsigned int buf[] = {0, 0, 1, 13,
                                      0, 6, 2, 14,
                                     10, 7, 3, 15,
                                     11, 8, 4, 16,
                                     12, 9, 5, 0};
                              
    
                MatrixType m(5, 5, buf, MatrixType::PolicySpecificDataHolderType(1, 2));

                unsigned int vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<unsigned int> v(5, vector_buf);

                NekVector<unsigned int> result = m*v;

                unsigned int expected_result_buf[] = { 43, 82, 129, 106, 89 };
                NekVector<unsigned int> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }
        }

        BOOST_AUTO_TEST_CASE(TestNoSubdiagonalsTwoSuperDiagonalsMatrixVectorMultiply)
        {
            {
                typedef NekMatrix<NekDouble, BandedMatrixTag> MatrixType;

                // This is the matrix
                // [  1  6 10   0  0 ]
                // [  0  2  7  11  0 ]
                // [  0  0  3   8 12 ]
                // [  0  0  0   4  9 ]
                // [  0  0  0   0  5 ]
                
                NekDouble buf[] = { 0, 0, 1,
                                    0, 6, 2,
                                    10, 7, 3,
                                    11, 8, 4,
                                    12, 9, 5};
                              
    
                MatrixType m(5, 5, buf, MatrixType::PolicySpecificDataHolderType(0, 2));

                NekDouble vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<NekDouble> v(5, vector_buf);

                NekVector<NekDouble> result = m*v;

                NekDouble expected_result_buf[] = { 43, 69, 101, 61, 25 };
                NekVector<NekDouble> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }

            {
                typedef NekMatrix<unsigned int, BandedMatrixTag> MatrixType;

                // This is the matrix
                // [  1  6 10   0  0 ]
                // [  0  2  7  11  0 ]
                // [  0  0  3   8 12 ]
                // [  0  0  0   4  9 ]
                // [  0  0  0   0  5 ]
                
                unsigned int buf[] = {0, 0, 1,
                                    0, 6, 2,
                                    10, 7, 3,
                                    11, 8, 4,
                                    12, 9, 5};
                              
    
                MatrixType m(5, 5, buf, MatrixType::PolicySpecificDataHolderType(0, 2));

                unsigned int vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<unsigned int> v(5, vector_buf);

                NekVector<unsigned int> result = m*v;

                unsigned int expected_result_buf[] = { 43, 69, 101, 61, 25 };
                NekVector<unsigned int> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }
        }

        BOOST_AUTO_TEST_CASE(TestNoSubdiagonalsThreeSuperDiagonalsMatrixVectorMultiply)
        {
            {
                typedef NekMatrix<NekDouble, BandedMatrixTag> MatrixType;

                // This is the matrix
                // [  1  6 10   2  0 ]
                // [  0  2  7  11  5 ]
                // [  0  0  3   8 12 ]
                // [  0  0  0   4  9 ]
                // [  0  0  0   0  5 ]
                
                NekDouble buf[] = { 0, 0, 0, 1,
                                    0, 0, 6, 2,
                                    0, 10, 7, 3,
                                    2, 11, 8, 4,
                                    5, 12, 9, 5};
                              
    
                MatrixType m(5, 5, buf, MatrixType::PolicySpecificDataHolderType(0, 3));

                NekDouble vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<NekDouble> v(5, vector_buf);

                NekVector<NekDouble> result = m*v;

                NekDouble expected_result_buf[] = { 51, 94, 101, 61, 25 };
                NekVector<NekDouble> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }

            {
                typedef NekMatrix<unsigned int, BandedMatrixTag> MatrixType;

                // This is the matrix
                // [  1  6 10   0  0 ]
                // [  0  2  7  11  0 ]
                // [  0  0  3   8 12 ]
                // [  0  0  0   4  9 ]
                // [  0  0  0   0  5 ]
                
                unsigned int buf[] = {0, 0, 0, 1,
                                    0, 0, 6, 2,
                                    0, 10, 7, 3,
                                    2, 11, 8, 4,
                                    5, 12, 9, 5};
                              
                              
    
                MatrixType m(5, 5, buf, MatrixType::PolicySpecificDataHolderType(0, 3));

                unsigned int vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<unsigned int> v(5, vector_buf);

                NekVector<unsigned int> result = m*v;

                unsigned int expected_result_buf[] = { 51, 94, 101, 61, 25 };
                NekVector<unsigned int> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }
        }

        BOOST_AUTO_TEST_CASE(TestNoSubdiagonalsFourSuperDiagonalsMatrixVectorMultiply)
        {
            {
                typedef NekMatrix<NekDouble, BandedMatrixTag> MatrixType;

                // This is the matrix
                // [  1  6 10   2  8 ]
                // [  0  2  7  11  5 ]
                // [  0  0  3   8 12 ]
                // [  0  0  0   4  9 ]
                // [  0  0  0   0  5 ]
                
                NekDouble buf[] = { 0, 0, 0, 0, 1,
                                    0, 0, 0, 6, 2,
                                    0, 0, 10, 7, 3,
                                    0, 2, 11, 8, 4,
                                    8, 5, 12, 9, 5};
                              
    
                MatrixType m(5, 5, buf, MatrixType::PolicySpecificDataHolderType(0, 4));

                NekDouble vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<NekDouble> v(5, vector_buf);

                NekVector<NekDouble> result = m*v;

                NekDouble expected_result_buf[] = { 91, 94, 101, 61, 25 };
                NekVector<NekDouble> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }

            {
                typedef NekMatrix<unsigned int, BandedMatrixTag> MatrixType;

                // This is the matrix
                // [  1  6 10   0  0 ]
                // [  0  2  7  11  0 ]
                // [  0  0  3   8 12 ]
                // [  0  0  0   4  9 ]
                // [  0  0  0   0  5 ]
                
                unsigned int buf[] = {0, 0, 0, 0, 1,
                                    0, 0, 0, 6, 2,
                                    0, 0, 10, 7, 3,
                                    0, 2, 11, 8, 4,
                                    8, 5, 12, 9, 5};
                              
                              
    
                MatrixType m(5, 5, buf, MatrixType::PolicySpecificDataHolderType(0, 4));

                unsigned int vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<unsigned int> v(5, vector_buf);

                NekVector<unsigned int> result = m*v;

                unsigned int expected_result_buf[] = { 91, 94, 101, 61, 25 };
                NekVector<unsigned int> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }
        }

        BOOST_AUTO_TEST_CASE(Test2Sub2SuperMatrixVectorMultiply)
        {
            {
                typedef NekMatrix<NekDouble, BandedMatrixTag> MatrixType;

                // This is the matrix
                // [  1  6 10   0  0 ]
                // [ 13  2  7  11  0 ]
                // [  7 14  3   8 12 ]
                // [  0  5 15   4  9 ]
                // [  0  0  4  16  5 ]
                
                NekDouble buf[] = {0, 0, 1, 13, 7,
                                   0, 6, 2, 14, 5,
                                   10, 7, 3, 15, 4,
                                   11, 8, 4, 16, 0,
                                   12, 9, 5, 0, 0};
                              
    
                MatrixType m(5, 5, buf, MatrixType::PolicySpecificDataHolderType(2, 2));

                NekDouble vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<NekDouble> v(5, vector_buf);

                NekVector<NekDouble> result = m*v;

                NekDouble expected_result_buf[] = { 43, 82, 136, 116, 101 };
                NekVector<NekDouble> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }

            {
                typedef NekMatrix<unsigned int, BandedMatrixTag> MatrixType;

                // This is the matrix
                // [  1  6 10   0  0 ]
                // [ 13  2  7  11  0 ]
                // [  7 14  3   8 12 ]
                // [  0  5 15   4  9 ]
                // [  0  0  4  16  5 ]
                
                unsigned int buf[] = {0, 0, 1, 13, 7,
                                   0, 6, 2, 14, 5,
                                   10, 7, 3, 15, 4,
                                   11, 8, 4, 16, 0,
                                   12, 9, 5, 0, 0};
                              
    
                MatrixType m(5, 5, buf, MatrixType::PolicySpecificDataHolderType(2, 2));

                unsigned int vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<unsigned int> v(5, vector_buf);

                NekVector<unsigned int> result = m*v;

                unsigned int expected_result_buf[] = { 43, 82, 136, 116, 101 };
                NekVector<unsigned int> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }
        }


        BOOST_AUTO_TEST_CASE(Test3Sub2SuperMatrixVectorMultiply)
        {
            {
                typedef NekMatrix<NekDouble, BandedMatrixTag> MatrixType;

                // This is the matrix
                // [  1  6 10   0  0 ]
                // [ 13  2  7  11  0 ]
                // [  7 14  3   8 12 ]
                // [ 20  5 15   4  9 ]
                // [  0 21  4  16  5 ]
                
                NekDouble buf[] = {0, 0, 1, 13, 7, 20,
                                   0, 6, 2, 14, 5, 21, 
                                   10, 7, 3, 15, 4, 0,
                                   11, 8, 4, 16, 0, 0,
                                   12, 9, 5, 0, 0, 0};
                              
    
                MatrixType m(5, 5, buf, MatrixType::PolicySpecificDataHolderType(3, 2));

                NekDouble vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<NekDouble> v(5, vector_buf);

                NekVector<NekDouble> result = m*v;

                NekDouble expected_result_buf[] = { 43, 82, 136, 136, 143 };
                NekVector<NekDouble> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }

            {
                typedef NekMatrix<unsigned int, BandedMatrixTag> MatrixType;

                // This is the matrix
                // [  1  6 10   0  0 ]
                // [ 13  2  7  11  0 ]
                // [  7 14  3   8 12 ]
                // [ 20  5 15   4  9 ]
                // [  0 21  4  16  5 ]
                
                unsigned int buf[] = {0, 0, 1, 13, 7, 20,
                                   0, 6, 2, 14, 5, 21, 
                                   10, 7, 3, 15, 4, 0,
                                   11, 8, 4, 16, 0, 0,
                                   12, 9, 5, 0, 0, 0};
                              
    
                MatrixType m(5, 5, buf, MatrixType::PolicySpecificDataHolderType(3, 2));

                unsigned int vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<unsigned int> v(5, vector_buf);

                NekVector<unsigned int> result = m*v;

                unsigned int expected_result_buf[] = { 43, 82, 136, 136, 143 };
                NekVector<unsigned int> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }
        }


        BOOST_AUTO_TEST_CASE(Test4Sub2SuperMatrixVectorMultiply)
        {
            {
                typedef NekMatrix<NekDouble, BandedMatrixTag> MatrixType;

                // This is the matrix
                // [  1  6 10   0  0 ]
                // [ 13  2  7  11  0 ]
                // [  7 14  3   8 12 ]
                // [ 20  5 15   4  9 ]
                // [ 30 21  4  16  5 ]
                
                NekDouble buf[] = {0, 0, 1, 13, 7, 20, 30,
                                   0, 6, 2, 14, 5, 21, 0,
                                   10, 7, 3, 15, 4, 0, 0,
                                   11, 8, 4, 16, 0, 0, 0,
                                   12, 9, 5, 0, 0, 0, 0};
                              
    
                MatrixType m(5, 5, buf, MatrixType::PolicySpecificDataHolderType(4, 2));

                NekDouble vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<NekDouble> v(5, vector_buf);

                NekVector<NekDouble> result = m*v;

                NekDouble expected_result_buf[] = { 43, 82, 136, 136, 173 };
                NekVector<NekDouble> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }

            {
                typedef NekMatrix<unsigned int, BandedMatrixTag> MatrixType;

                // This is the matrix
                // [  1  6 10   0  0 ]
                // [ 13  2  7  11  0 ]
                // [  7 14  3   8 12 ]
                // [ 20  5 15   4  9 ]
                // [ 30 21  4  16  5 ]
                
                unsigned int buf[] = {0, 0, 1, 13, 7, 20, 30,
                                   0, 6, 2, 14, 5, 21, 0,
                                   10, 7, 3, 15, 4, 0, 0,
                                   11, 8, 4, 16, 0, 0, 0,
                                   12, 9, 5, 0, 0, 0, 0};
                              
    
                MatrixType m(5, 5, buf, MatrixType::PolicySpecificDataHolderType(4, 2));

                unsigned int vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<unsigned int> v(5, vector_buf);

                NekVector<unsigned int> result = m*v;

                unsigned int expected_result_buf[] = { 43, 82, 136, 136, 173 };
                NekVector<unsigned int> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }
        }

        BOOST_AUTO_TEST_CASE(TestDiagonalOnlyMatrixVectorMultiply)
        {
            {
                typedef NekMatrix<NekDouble, BandedMatrixTag> MatrixType;

                // This is the matrix
                // [  5  0  0   0  0 ]
                // [  0  1  0   0  0 ]
                // [  0  0  9   0  0 ]
                // [  0  0  0   4  0 ]
                // [  0  0  0   0  2 ]
                
                NekDouble buf[] = {5, 1, 9, 4, 2};
                              
    
                MatrixType m(5, 5, buf, MatrixType::PolicySpecificDataHolderType(0, 0));

                NekDouble vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<NekDouble> v(5, vector_buf);

                NekVector<NekDouble> result = m*v;

                NekDouble expected_result_buf[] = { 5, 2, 27, 16, 10 };
                NekVector<NekDouble> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }

            {
                typedef NekMatrix<unsigned int, BandedMatrixTag> MatrixType;

                // This is the matrix
                // [  5  0  0   0  0 ]
                // [  0  1  0   0  0 ]
                // [  0  0  9   0  0 ]
                // [  0  0  0   4  0 ]
                // [  0  0  0   0  2 ]
                
                unsigned int buf[] = {5, 1, 9, 4, 2};
                              
    
                MatrixType m(5, 5, buf, MatrixType::PolicySpecificDataHolderType(0, 0));

                unsigned int vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<unsigned int> v(5, vector_buf);

                NekVector<unsigned int> result = m*v;

                unsigned int expected_result_buf[] = { 5, 2, 27, 16, 10 };
                NekVector<unsigned int> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }
        }
    }

    namespace BandedMatrixMatrixMultiplicationTests
    {
        BOOST_AUTO_TEST_CASE(TestMatrixMatrixMultiplication)
        {
            // [ 2  10   0  0 ]    [ 30 38  0  0 ]
            // [16   4  12  0 ]    [  0 32 40  0 ]
            // [22  18   6 14 ] *  [  0  0 34 42 ]
            // [ 0  24  20  8 ]    [  0  0  0 36 ]
            
            {
                typedef MatrixStoragePolicy<double, BandedMatrixTag> Policy;
                double lhs_buf[] = { 0, 2, 16, 22, 
                                     10, 4, 18, 24, 
                                     12, 6, 20, 0,
                                     14, 8, 0, 0 };
                double rhs_buf[] = { 0, 30,
                                     38, 32,
                                     40, 34,
                                     42, 36 };

                Policy::PolicySpecificDataHolderType lhsPolicyData(2, 1);
                NekMatrix<double, BandedMatrixTag> lhs1(4, 4, lhs_buf, lhsPolicyData);
                boost::shared_ptr<NekMatrix<NekMatrix<double, BandedMatrixTag, StandardMatrixTag>, BandedMatrixTag, ScaledMatrixTag> > lhs2;
                boost::shared_ptr<NekMatrix<NekMatrix<double>, BandedMatrixTag, BlockMatrixTag> > lhs3;

                Policy::PolicySpecificDataHolderType rhsPolicyData(0, 1);
                NekMatrix<double, BandedMatrixTag> rhs1(4, 4, rhs_buf, rhsPolicyData);
                boost::shared_ptr<NekMatrix<NekMatrix<double, BandedMatrixTag>, BandedMatrixTag, ScaledMatrixTag> > rhs2;
                boost::shared_ptr<NekMatrix<NekMatrix<double>, BandedMatrixTag, BlockMatrixTag> > rhs3;

                GenerateBandedMatrices(lhs1, 2.0, 2, 2, lhs2, lhs3);
                GenerateBandedMatrices(rhs1, 2.0, 2, 2, rhs2, rhs3);

                double result_buf[] = {60, 480, 660, 0,
                                       396, 736, 1412, 768,
                                       400, 568, 924, 1640,
                                       0, 504, 756, 1128};
                NekMatrix<double> result(4, 4, result_buf);
                
                RunAllTestCombinations(lhs1, *lhs2, *lhs3, rhs1, *rhs2, *rhs3, result, DoMultiplication());
            }

            {
            }
        }
    }
}


