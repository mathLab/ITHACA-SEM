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
#include <LibUtilities/LinearAlgebra/NekLinSys.hpp>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>

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
            NekVector<NekDouble> expected_result(4, expected_result_buf);
            NekVector<NekDouble> result(4, y);
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
            NekVector<NekDouble> expected_result(4, expected_result_buf);
            NekVector<NekDouble> result(4, y);
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
            NekVector<NekDouble> expected_result(4, expected_result_buf);
            NekVector<NekDouble> result(4, y);
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
            NekVector<NekDouble> expected_result(3, expected_result_buf);
            NekVector<NekDouble> result(3, y);
            BOOST_CHECK_EQUAL(expected_result, result);
        }

        BOOST_AUTO_TEST_CASE(TestStandardMatrixVectorMultiply)
        {
            
            {
                // This is the matrix
                // [ 1 2 0 0 ]
                // [ 3 4 5 0 ]
                // [ 0 6 7 8 ]
                // [ 0 0 9 10]
                
                NekDouble buf[] = {0, 1, 3,
                                   2, 4, 6,
                                   5, 7, 9,
                                   8, 10, 0};
                              
    
                std::shared_ptr<DenseMatrix> m(new DenseMatrix(4, 4, buf, eBANDED, 1, 1));

                std::shared_ptr<ScaledMatrix> scaled(new ScaledMatrix(2.0, m));
                
                NekDouble vector_buf[] = {1.0, 2.0, 3.0, 4.0};
                NekVector<NekDouble> v(4, vector_buf);

                NekVector<NekDouble> result = (*m)*v;
                NekVector<NekDouble> scaledResult = (*scaled)*v;

                NekDouble expected_result_buf[] = { 5, 26, 65, 67 };
                NekVector<NekDouble> expected_result(4, expected_result_buf);
                NekVector<NekDouble> expectedScaledResult = 2.0*expected_result;
                BOOST_CHECK_EQUAL(expected_result, result);

                BOOST_CHECK_EQUAL(expectedScaledResult, scaledResult);
            }
        }

        BOOST_AUTO_TEST_CASE(TestUnequalNumbersOfSubAndSuperDiagonalsMatrixVectorMultiply)
        {
            {
                typedef NekMatrix<NekDouble> MatrixType;

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
                              
    
                MatrixType m(5, 5, buf, eBANDED, 1, 2);

                NekDouble vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<NekDouble> v(5, vector_buf);

                NekVector<NekDouble> result = m*v;

                NekDouble expected_result_buf[] = { 43, 82, 129, 106, 89 };
                NekVector<NekDouble> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }

        }

        BOOST_AUTO_TEST_CASE(TestNoSubdiagonalsTwoSuperDiagonalsMatrixVectorMultiply)
        {
            {
                typedef NekMatrix<NekDouble> MatrixType;

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
                              
    
                MatrixType m(5, 5, buf, eBANDED, 0, 2);

                NekDouble vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<NekDouble> v(5, vector_buf);

                NekVector<NekDouble> result = m*v;

                NekDouble expected_result_buf[] = { 43, 69, 101, 61, 25 };
                NekVector<NekDouble> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }

        }

        BOOST_AUTO_TEST_CASE(TestNoSubdiagonalsThreeSuperDiagonalsMatrixVectorMultiply)
        {
            {
                typedef NekMatrix<NekDouble> MatrixType;

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
                              
    
                MatrixType m(5, 5, buf, eBANDED, 0, 3);

                NekDouble vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<NekDouble> v(5, vector_buf);

                NekVector<NekDouble> result = m*v;

                NekDouble expected_result_buf[] = { 51, 94, 101, 61, 25 };
                NekVector<NekDouble> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }

        }

        BOOST_AUTO_TEST_CASE(TestNoSubdiagonalsFourSuperDiagonalsMatrixVectorMultiply)
        {
            {
                typedef NekMatrix<NekDouble> MatrixType;

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
                              
    
                MatrixType m(5, 5, buf, eBANDED, 0, 4);

                NekDouble vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<NekDouble> v(5, vector_buf);

                NekVector<NekDouble> result = m*v;

                NekDouble expected_result_buf[] = { 91, 94, 101, 61, 25 };
                NekVector<NekDouble> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }

        }

        BOOST_AUTO_TEST_CASE(Test2Sub2SuperMatrixVectorMultiply)
        {
            {
                typedef NekMatrix<NekDouble> MatrixType;

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
                              
    
                MatrixType m(5, 5, buf, eBANDED, 2, 2);

                NekDouble vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<NekDouble> v(5, vector_buf);

                NekVector<NekDouble> result = m*v;

                NekDouble expected_result_buf[] = { 43, 82, 136, 116, 101 };
                NekVector<NekDouble> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }
        }


        BOOST_AUTO_TEST_CASE(Test3Sub2SuperMatrixVectorMultiply)
        {
            {
                typedef NekMatrix<NekDouble> MatrixType;

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
                              
    
                MatrixType m(5, 5, buf, eBANDED, 3, 2);

                NekDouble vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<NekDouble> v(5, vector_buf);

                NekVector<NekDouble> result = m*v;

                NekDouble expected_result_buf[] = { 43, 82, 136, 136, 143 };
                NekVector<NekDouble> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }

        }


        BOOST_AUTO_TEST_CASE(Test4Sub2SuperMatrixVectorMultiply)
        {
            {
                typedef NekMatrix<NekDouble> MatrixType;

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
                              
    
                MatrixType m(5, 5, buf, eBANDED, 4, 2);

                NekDouble vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<NekDouble> v(5, vector_buf);

                NekVector<NekDouble> result = m*v;

                NekDouble expected_result_buf[] = { 43, 82, 136, 136, 173 };
                NekVector<NekDouble> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }
        }

        BOOST_AUTO_TEST_CASE(TestDiagonalOnlyMatrixVectorMultiply)
        {
            {
                typedef NekMatrix<NekDouble> MatrixType;

                // This is the matrix
                // [  5  0  0   0  0 ]
                // [  0  1  0   0  0 ]
                // [  0  0  9   0  0 ]
                // [  0  0  0   4  0 ]
                // [  0  0  0   0  2 ]
                
                NekDouble buf[] = {5, 1, 9, 4, 2};
                              
    
                MatrixType m(5, 5, buf, eBANDED, 0, 0);

                NekDouble vector_buf[] = {1, 2, 3, 4, 5};
                NekVector<NekDouble> v(5, vector_buf);

                NekVector<NekDouble> result = m*v;

                NekDouble expected_result_buf[] = { 5, 2, 27, 16, 10 };
                NekVector<NekDouble> expected_result(5, expected_result_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }
        }
        
        BOOST_AUTO_TEST_CASE(TestBandedMatrixLinearSystemSolves)
        {
            typedef NekMatrix<double> MatrixType;

            {
                double buf[] = {0, 0, -.23, -6.98,
                                0, 2.54, 2.46, 2.56,
                                -3.66, -2.73, 2.46, -4.78,
                                -2.13, 4.07, -3.82, 0};
                              
    
                MatrixType m(4,4, buf, eBANDED, 1, 2);                
                LinearSystem l(m);

                double b_buf[] =  { 4.42, 27.13, -6.14, 10.5 };
                NekVector<double> b(4, b_buf);

                NekVector<double> result = l.Solve(b);
            }

            {
                

                // This is the matrix
                // [  1  6 10   0  0 ]
                // [ 13  2  7  11  0 ]
                // [  7 14  3   8 12 ]
                // [ 20  5 15   4  9 ]
                // [ 30 21  4  16  5 ]
                
                double buf[] = {0, 0, 1, 13, 7, 20, 30,
                                   0, 6, 2, 14, 5, 21, 0,
                                   10, 7, 3, 15, 4, 0, 0,
                                   11, 8, 4, 16, 0, 0, 0,
                                   12, 9, 5, 0, 0, 0, 0};
                              
    
                MatrixType m(5, 5, buf, eBANDED, 4, 2);

                double b_buf[] =  { 43, 82, 136, 136, 173 };
                NekVector<double> b(5, b_buf);
                
                LinearSystem l(m);
                NekVector<double> result(5, 0.0);
                result = l.Solve(b);
                
                double expected_result_buf[] = {1, 2, 3, 4, 5};
                NekVector<double> expected_result(5, expected_result_buf);

                BOOST_CHECK_CLOSE(expected_result[0], result[0], .00001);
                BOOST_CHECK_CLOSE(expected_result[1], result[1], .00001);
                BOOST_CHECK_CLOSE(expected_result[2], result[2], .00001);
                BOOST_CHECK_CLOSE(expected_result[3], result[3], .00001);
                BOOST_CHECK_CLOSE(expected_result[4], result[4], .00001);
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
                double lhs_buf[] = { 0, 2, 16, 22, 
                                     10, 4, 18, 24, 
                                     12, 6, 20, 0,
                                     14, 8, 0, 0 };
                double rhs_buf[] = { 0, 30,
                                     38, 32,
                                     40, 34,
                                     42, 36 };

                NekMatrix<double> lhs1(4, 4, lhs_buf, eBANDED, 2, 1);
                std::shared_ptr<NekMatrix<NekMatrix<double, StandardMatrixTag>, ScaledMatrixTag> > lhs2;
                std::shared_ptr<NekMatrix<NekMatrix<double>, BlockMatrixTag> > lhs3;

                NekMatrix<double> rhs1(4, 4, rhs_buf, eBANDED, 0, 1);
                std::shared_ptr<NekMatrix<NekMatrix<double>, ScaledMatrixTag> > rhs2;
                std::shared_ptr<NekMatrix<NekMatrix<double>, BlockMatrixTag> > rhs3;

                GenerateMatrices(lhs1, 2.0, 2, 2, lhs2, lhs3);
                GenerateMatrices(rhs1, 2.0, 2, 2, rhs2, rhs3);

                double result_buf[] = {60, 480, 660, 0,
                                       396, 736, 1412, 768,
                                       400, 568, 924, 1640,
                                       0, 504, 756, 1128};
                NekMatrix<double> result(4, 4, result_buf);
            }
        }
    }
}


