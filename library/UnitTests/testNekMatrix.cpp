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
#include <iostream>

namespace Nektar
{
    //template class NekMatrix<double, eFull>;
    namespace UnitTests
    {
       
        using namespace Nektar;

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

                NekMatrix<float> dynamic_matrix(7, 3, (float)7.8);

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
                    BOOST_CHECK(static_matrix(i,j) == 10*i + j);

                    const NekMatrix<unsigned int>& ref = static_matrix;

                    BOOST_CHECK(ref(i,j) == 10*i + j);
                    BOOST_CHECK(ref(i,j) == 10*i + j);

                }
            }

            // Invalid access is an unrecoverable error.
            //BOOST_CHECK_THROW(static_matrix(3,2), OutOfBoundsError);
            //BOOST_CHECK_THROW(static_matrix(2,3), OutOfBoundsError);
            BOOST_CHECK_NO_THROW(static_matrix(2,2));

        }

        void testNekMatrixBasicMath()
        {
//             // Addition tests.
//             {
//                 double buf[] = {1.0, 2.0, 3.0,
//                     4.0, 5.0, 6.0,
//                     7.0, 8.0, 9.0 };
// 
//                 NekMatrix<double> m1(3, 3, buf);
//                 NekMatrix<double> m2(3, 3, buf);
//                 NekMatrix<double> m3 = m1 + m2;
// 
//                 for(unsigned int i = 0; i < 3; ++i)
//                 {
//                     for(unsigned int j = 0; j < 3; ++j)
//                     {
//                         BOOST_CHECK(m3(i,j) == buf[3*i+j] + buf[3*i+j]);
//                     }
//                 }
// 
//                 NekMatrix<double> m4(3, 3, buf);
//                 NekMatrix<double> m5(3, 3, buf);
//                 NekMatrix<double> m6 = m4+m5;
// 
//                 for(unsigned int i = 0; i < 3; ++i)
//                 {
//                     for(unsigned int j = 0; j < 3; ++j)
//                     {
//                         BOOST_CHECK(m6(i,j) == buf[3*i+j] + buf[3*i+j]);
//                     }
//                 }
//             }
//             // Multiply
//             {
//                 unsigned int buf1[] = {1, 2, 3,
//                                        4, 5, 6,
//                                        7, 8, 9};
//                 unsigned int buf2[] = { 10, 11, 12, 14,
//                                         15, 16, 17, 18,
//                                         19, 20, 21, 22 };
// 
//                                        
//                 NekMatrix<unsigned int> lhs(3, 3, buf1);
//                 NekMatrix<unsigned int> rhs(3, 4, buf2);
//                 NekMatrix<unsigned int> result = lhs*rhs;
// 
//                 BOOST_CHECK(result.GetRows() == 3);
//                 BOOST_CHECK(result.GetColumns() == 4);
// 
//                 BOOST_CHECK(result(0,0) == 97);
//                 BOOST_CHECK(result(0,1) == 103);
//                 BOOST_CHECK(result(0,2) == 109);
//                 BOOST_CHECK(result(0,3) == 116);
// 
//                 BOOST_CHECK(result(1,0) == 229);
//                 BOOST_CHECK(result(1,1) == 244);
//                 BOOST_CHECK(result(1,2) == 259);
//                 BOOST_CHECK(result(1,3) == 278);
// 
//                 BOOST_CHECK(result(2,0) == 361);
//                 BOOST_CHECK(result(2,1) == 385);
//                 BOOST_CHECK(result(2,2) == 409);
//                 BOOST_CHECK(result(2,3) == 440);
//             }
// 
//             {
//                 double buf1[] = {1, 2, 3,
//                                 4, 5, 6,
//                                 7, 8, 9};
//                 double buf2[] = { 10, 11, 12,
//                     15, 16, 17,
//                     19, 20, 21 };
// 
//                                  
//                 NekMatrix<double> lhs(3, 3, buf1);
//                 NekMatrix<double> rhs(3, 3, buf2);
// 
//                 NekMatrix<double> result = lhs*rhs;
// 
//                 BOOST_CHECK(result.GetRows() == 3);
//                 BOOST_CHECK(result.GetColumns() == 3);
// 
//                 double epsilon = 1e-12;
//                 BOOST_CHECK_CLOSE(result(0,0), 97.0, epsilon);
//                 BOOST_CHECK_CLOSE(result(0,1), 103.0, epsilon);
//                 BOOST_CHECK_CLOSE(result(0,2), 109.0, epsilon);
// 
//                 BOOST_CHECK_CLOSE(result(1,0), 229.0, epsilon);
//                 BOOST_CHECK_CLOSE(result(1,1), 244.0, epsilon);
//                 BOOST_CHECK_CLOSE(result(1,2), 259.0, epsilon);
// 
//                 BOOST_CHECK_CLOSE(result(2,0), 361.0, epsilon);
//                 BOOST_CHECK_CLOSE(result(2,1), 385.0, epsilon);
//                 BOOST_CHECK_CLOSE(result(2,2), 409.0, epsilon);
//             }
//         
//             {
//                 double buf1[] = {1, 2, 3,
//                                         4, 5, 6,
//                                         7, 8, 9};
//                 double buf2[] = { 10, 11, 12, 14,
//                                     15, 16, 17, 18,
//                                     19, 20, 21, 22 };
// 
//                 NekMatrix<double> lhs(3, 3, buf1);
//                 NekMatrix<double> rhs(3, 4, buf2);
// 
//                 NekMatrix<double> result = lhs*rhs;
// 
//                 BOOST_CHECK(result.GetRows() == 3);
//                 BOOST_CHECK(result.GetColumns() == 4);
// 
//                 double epsilon = 1e-12;
//                 BOOST_CHECK_CLOSE(result(0,0), 97.0, epsilon);
//                 BOOST_CHECK_CLOSE(result(0,1), 103.0, epsilon);
//                 BOOST_CHECK_CLOSE(result(0,2), 109.0, epsilon);
//                 BOOST_CHECK_CLOSE(result(0,3), 116.0, epsilon);
// 
//                 BOOST_CHECK_CLOSE(result(1,0), 229.0, epsilon);
//                 BOOST_CHECK_CLOSE(result(1,1), 244.0, epsilon);
//                 BOOST_CHECK_CLOSE(result(1,2), 259.0, epsilon);
//                 BOOST_CHECK_CLOSE(result(1,3), 278.0, epsilon);
// 
//                 BOOST_CHECK_CLOSE(result(2,0), 361.0, epsilon);
//                 BOOST_CHECK_CLOSE(result(2,1), 385.0, epsilon);
//                 BOOST_CHECK_CLOSE(result(2,2), 409.0, epsilon);
//                 BOOST_CHECK_CLOSE(result(2,3), 440.0, epsilon);
//             }
//             
//             {
//                 unsigned int buf1[] = {1, 2, 3,
//                                        4, 5, 6,
//                                        7, 8, 9};
// 
//                                        
//                 unsigned int buf2[] = { 1, 2, 3};
// 
//                 NekMatrix<unsigned int> lhs(3, 3, buf1);
//                 NekVector<unsigned int> rhs(3, buf2);
// 
//                 NekVector<unsigned int> result = lhs*rhs;
// 
//                 BOOST_CHECK(result[0] == 14);
//                 BOOST_CHECK(result[1] == 32);
//                 BOOST_CHECK(result[2] == 50);
//             }
// 
//             // Negation
//             {
//                 int buf[] = {1, 2, 3,
//                                        4, 5, 6,
//                                        7, 8, 9};
// 
//                 int neg_buf[] = {-1, -2, -3,
//                                        -4, -5, -6,
//                                        -7, -8, -9};
// 
//                 NekMatrix<int> m1(3,3, buf);
//                 NekMatrix<int> m2 = -m1;
//                 NekMatrix<int> m3 = -m2;
// 
//                 BOOST_CHECK_EQUAL(m1, m3);
// 
//                 NekMatrix<int> negated(3,3,neg_buf);
//                 BOOST_CHECK_EQUAL(m2, negated);
// 
//                 NekMatrix<int> m4 = -(-m1);
//                 BOOST_CHECK_EQUAL(m4, m1);
// 
//                 NekMatrix<int> m5 = -(-(-(-(-(-(-(-(-m1))))))));
//                 BOOST_CHECK_EQUAL(m5, negated);
//             }
// 
//             // Transpose
// 
//             // Determinant.
// 
//             // Invert/Check for singularity.
// 
//             // Eigenvalues/vectors
// 
//             // Condition number wrt various norms.
// 
//             // Various norm computations.
// 
//             // LU Decomposition?  More appropriate in LinAlg?
        }

        void testDiagonalMatrix()
        {
        }

        void testNekMatrixFullDiagonalOperations()
        {
//             unsigned int fullValues[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
//             
//             NekMatrix<unsigned int> full(3, 3, fullValues);
// 
//             unsigned int diagonalValues[] = {6, 12, 5};
//             NekMatrix<unsigned int, eDiagonal> diag(3, diagonalValues);
// 
//             NekMatrix<unsigned int> result1 = full+diag;
//             NekMatrix<unsigned int> result2 = diag+full;
// 
//             BOOST_CHECK_EQUAL(result1, result2);
//             BOOST_CHECK_EQUAL(result1(0,0),7);
//             BOOST_CHECK_EQUAL(result1(0,1),2);
//             BOOST_CHECK_EQUAL(result1(0,2),3);
// 
//             BOOST_CHECK_EQUAL(result1(1,0),4);
//             BOOST_CHECK_EQUAL(result1(1,1),17);
//             BOOST_CHECK_EQUAL(result1(1,2),6);
// 
//             BOOST_CHECK_EQUAL(result1(2,0),7);
//             BOOST_CHECK_EQUAL(result1(2,1),8);
//             BOOST_CHECK_EQUAL(result1(2,2),14);
// 
         }
 
         void testUserManagedMatrixData()
         {    
//             {
//                 unsigned int matrixValues[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
//                 
//                 NekMatrix<unsigned int> m(3, 3, matrixValues, eWrapper);
// 
//                 BOOST_CHECK_EQUAL(m(0,0), 1);
//                 BOOST_CHECK_EQUAL(m(0,1), 2);
//                 BOOST_CHECK_EQUAL(m(0,2), 3);
// 
//                 BOOST_CHECK_EQUAL(m(1,0), 4);
//                 BOOST_CHECK_EQUAL(m(1,1), 5);
//                 BOOST_CHECK_EQUAL(m(1,2), 6);
// 
//                 BOOST_CHECK_EQUAL(m(2,0), 7);
//                 BOOST_CHECK_EQUAL(m(2,1), 8);
//                 BOOST_CHECK_EQUAL(m(2,2), 9);
// 
//                 m(0,0) = 18;
//                 BOOST_CHECK_EQUAL(m(0,0), 18);
//                 BOOST_CHECK_EQUAL(matrixValues[0], m(0,0));
// 
//                 NekMatrix<unsigned int> m1(m);
//                 m1(1,0) = 900;
//                 BOOST_CHECK_EQUAL(m1(1,0), 900);
//                 BOOST_CHECK_EQUAL(m(1,0), 4);
//                 BOOST_CHECK_EQUAL(matrixValues[3], 4);
// 
//                 NekMatrix<unsigned int> m2(3, 3);
//                 m2 = m1;
//                 m2(0,0) = 800;
//                 BOOST_CHECK_EQUAL(m1(0,0), 18);
//                 BOOST_CHECK_EQUAL(m2(0,0), 800);
//             }
// 
//             {
//                 unsigned int matrixValues[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
//                 NekMatrix<unsigned int> m(3, 3, matrixValues, eCopy);
// 
//                 m(0,0) = 800;
//                 BOOST_CHECK_EQUAL(m(0,0), 800);
//                 BOOST_CHECK_EQUAL(matrixValues[0], 1);
//             }
        }
        
        void testBlockDiagonalMatrices()
        {
/*            {
                typedef NekMatrix<double, eDiagonal, eBlock> BlockMatrix;
                typedef BlockMatrix::InnerMatrixType InnerMatrixType;
                
                // Create a 6x6 block diagonal matrix.  Each block is a 2x2 matrix.Blas
                BlockMatrix m(3, 2, 2);
                InnerMatrixType d1(2, 2, 7.0);
                InnerMatrixType d2(2, 2, 8.0);
                InnerMatrixType d3(2, 2, 9.0);
                InnerMatrixType d4(2, 2, 0.0);
                InnerMatrixType d5(3, 3, 9.0);
                
                m.GetBlock(0,0) = d1;
                m.GetBlock(1,1) = d2;
                m.GetBlock(2,2) = d3;
                
                BOOST_CHECK_EQUAL(m.GetBlock(0,0), d1);
                BOOST_CHECK_EQUAL(m.GetBlock(1,1), d2);
                BOOST_CHECK_EQUAL(m.GetBlock(2,2), d3);
                BOOST_CHECK(m.GetBlock(0,0) != d3);
                
                BOOST_CHECK_EQUAL(m.GetBlock(0,1), d4);
                BOOST_CHECK_EQUAL(m.GetBlock(0,2), d4);
                BOOST_CHECK_EQUAL(m.GetBlock(1,0), d4);
                BOOST_CHECK_EQUAL(m.GetBlock(1,2), d4);
                BOOST_CHECK_EQUAL(m.GetBlock(2,0), d4);
                BOOST_CHECK_EQUAL(m.GetBlock(2,1), d4);
                
                BOOST_CHECK(m.GetBlock(2,2) != d5);
                */
//            }
            
//             {
//                 typedef NekMatrix<double, eDiagonal, ePointerBlock> BlockMatrix;
//                 typedef BlockMatrix::InnerMatrixType InnerMatrixType;
//                 
//                 // Create a 6x6 block diagonal matrix.  Each block is a 2x2 matrix.Blas
//                 BlockMatrix m(3, 2, 2);
//                 boost::shared_ptr<InnerMatrixType> d1 = boost::shared_ptr<InnerMatrixType>(new InnerMatrixType(2, 2, 7.0));
//                 boost::shared_ptr<InnerMatrixType> d2 = boost::shared_ptr<InnerMatrixType>(new InnerMatrixType(2, 2, 8.0));
//                 boost::shared_ptr<InnerMatrixType> d3 = boost::shared_ptr<InnerMatrixType>(new InnerMatrixType(2, 2, 9.0));
//                 boost::shared_ptr<InnerMatrixType> d4 = boost::shared_ptr<InnerMatrixType>(new InnerMatrixType(2, 2, 0.0));
//                 boost::shared_ptr<InnerMatrixType> d5 = boost::shared_ptr<InnerMatrixType>(new InnerMatrixType(3, 3, 9.0));
// 
//                 m(0,0) = d1;
//                 m(1,1) = d2;
//                 m(2,2) = d3;
// 
//                 BOOST_CHECK(m(0,0) == d1);
//                 BOOST_CHECK(m(1,1) == d2);
//                 BOOST_CHECK(m(2,2) == d3);
//                 BOOST_CHECK(m(0,0) != d3);
// 
//                 BOOST_CHECK(m(0,1) == d4);
//                 BOOST_CHECK(m(0,2) == d4);
//                 BOOST_CHECK(m(1,0) == d4);
//                 BOOST_CHECK(m(1,2) == d4);
//                 BOOST_CHECK(m(2,0) == d4);
//                 BOOST_CHECK(m(2,1) == d4);
// 
//                 BOOST_CHECK(m(2,2) != d5);
//             }
            
//             // Block Matrix operators.
//             {
//                 typedef NekMatrix<int, eDiagonal, eBlock> BlockMatrix;
//                 typedef BlockMatrix::InnerMatrixType InnerMatrixType;
//                 
//                 BlockMatrix m1(2, 2, 2);
//                 const int m1_buf0[] = {-85, -55, 97, 50};
//                 const int m1_buf3[] = {57, -59, -93, 92};
//                 
//                 m1.GetBlock(0,0) = InnerMatrixType(2, 2, m1_buf0);
//                 m1.GetBlock(1,1) = InnerMatrixType(2, 2, m1_buf3);
//                 
//                 BOOST_CHECK_EQUAL(m1(0,0), -85);
//                 BOOST_CHECK_EQUAL(m1(0,1), -55);
//                 BOOST_CHECK_EQUAL(m1(0,2), 0);
//                 BOOST_CHECK_EQUAL(m1(0,3), 0);
//                 
//                 BOOST_CHECK_EQUAL(m1(1,0), 97);
//                 BOOST_CHECK_EQUAL(m1(1,1), 50);
//                 BOOST_CHECK_EQUAL(m1(1,2), 0);
//                 BOOST_CHECK_EQUAL(m1(1,3), 0);
//                 
//                 BOOST_CHECK_EQUAL(m1(2,0), 0);
//                 BOOST_CHECK_EQUAL(m1(2,1), 0);
//                 BOOST_CHECK_EQUAL(m1(2,2), 57);
//                 BOOST_CHECK_EQUAL(m1(2,3), -59);
//                 
//                 BOOST_CHECK_EQUAL(m1(3,0), 0);
//                 BOOST_CHECK_EQUAL(m1(3,1), 0);
//                 BOOST_CHECK_EQUAL(m1(3,2), -93);
//                 BOOST_CHECK_EQUAL(m1(3,3), 92);
//                 
//                 BlockMatrix m2(2, 2, 2);
//                 
//                 const int m2_buf0[] = {43, -62, 54, -5};
//                 const int m2_buf3[] = {-18, 31, 1, -47};
//                 
//                 m2.GetBlock(0,0) = InnerMatrixType(2, 2, m2_buf0);
//                 m2.GetBlock(1,1) = InnerMatrixType(2, 2, m2_buf3);
//                 
//                 BOOST_CHECK_EQUAL(m2(0,0), 43);
//                 BOOST_CHECK_EQUAL(m2(0,1), -62);
//                 BOOST_CHECK_EQUAL(m2(0,2), 0);
//                 BOOST_CHECK_EQUAL(m2(0,3), 0);
//                 
//                 BOOST_CHECK_EQUAL(m2(1,0), 54);
//                 BOOST_CHECK_EQUAL(m2(1,1), -5);
//                 BOOST_CHECK_EQUAL(m2(1,2), 0);
//                 BOOST_CHECK_EQUAL(m2(1,3), 0);
//                 
//                 BOOST_CHECK_EQUAL(m2(2,0), 0);
//                 BOOST_CHECK_EQUAL(m2(2,1), 0);
//                 BOOST_CHECK_EQUAL(m2(2,2), -18);
//                 BOOST_CHECK_EQUAL(m2(2,3), 31);
//                 
//                 BOOST_CHECK_EQUAL(m2(3,0), 0);
//                 BOOST_CHECK_EQUAL(m2(3,1), 0);
//                 BOOST_CHECK_EQUAL(m2(3,2), 1);
//                 BOOST_CHECK_EQUAL(m2(3,3), -47);
//                 
//                 BlockMatrix R = m1+m2;
//                 BOOST_CHECK_EQUAL(R(0,0), m1(0,0) + m2(0,0));
//                 BOOST_CHECK_EQUAL(R(0,1), m1(0,1) + m2(0,1));
//                 BOOST_CHECK_EQUAL(R(0,2), 0);
//                 BOOST_CHECK_EQUAL(R(0,3), 0);
//                 
//                 BOOST_CHECK_EQUAL(R(1,0), m1(1,0) + m2(1,0));
//                 BOOST_CHECK_EQUAL(R(1,1), m1(1,1) + m2(1,1));
//                 BOOST_CHECK_EQUAL(R(1,2), 0);
//                 BOOST_CHECK_EQUAL(R(1,3), 0);
//                 
//                 BOOST_CHECK_EQUAL(R(2,0), 0);
//                 BOOST_CHECK_EQUAL(R(2,1), 0);
//                 BOOST_CHECK_EQUAL(R(2,2), m1(2,2) + m2(2,2));
//                 BOOST_CHECK_EQUAL(R(2,3), m1(2,3) + m2(2,3));
//                 
//                 BOOST_CHECK_EQUAL(R(3,0), 0);
//                 BOOST_CHECK_EQUAL(R(3,1), 0);
//                 BOOST_CHECK_EQUAL(R(3,2), m1(3,2) + m2(3,2));
//                 BOOST_CHECK_EQUAL(R(3,3), m1(3,3) + m2(3,3));
//             }
//             

            
        }
        
        void testBlockDiagonalTimesEqual()
        {
/*            {
                typedef NekMatrix<int, eDiagonal, eBlock> BlockMatrix;
                typedef BlockMatrix::InnerMatrixType InnerMatrixType;
                
                BlockMatrix m1(4, 2, 2);
                const int m1_buf0[] = {-91, -47, 94, 83};
                const int m1_buf1[] = {49, 78, 80, 72};
                const int m1_buf2[] = {87, 79, 31, -34};
                const int m1_buf3[] = {9, 29, -17, -98};
                
                m1.GetBlock(0,0) = InnerMatrixType(2, 2, m1_buf0);
                m1.GetBlock(1,1) = InnerMatrixType(2, 2, m1_buf1);
                m1.GetBlock(2,2) = InnerMatrixType(2, 2, m1_buf2);
                m1.GetBlock(3,3) = InnerMatrixType(2, 2, m1_buf3);
                
                BlockMatrix m2(4, 2, 2);
                const int m2_buf0[] = {-36, 40, 4, -59};
                const int m2_buf1[] = {62, 11, 4, -11};
                const int m2_buf2[] = {-94, -68, -73, -91};
                const int m2_buf3[] = {-39, 8, 68, 45};
                
                m2.GetBlock(0,0) = InnerMatrixType(2, 2, m2_buf0);
                m2.GetBlock(1,1) = InnerMatrixType(2, 2, m2_buf1);
                m2.GetBlock(2,2) = InnerMatrixType(2, 2, m2_buf2);
                m2.GetBlock(3,3) = InnerMatrixType(2, 2, m2_buf3);
                
                m1 *= m2;
                
                BOOST_CHECK_EQUAL(m1(0,0), 3088);
                BOOST_CHECK_EQUAL(m1(0,1), -867);
                BOOST_CHECK_EQUAL(m1(1,0), -3052);
                BOOST_CHECK_EQUAL(m1(1,1), -1137);
                
                BOOST_CHECK_EQUAL(m1(2,2), 3350);
                BOOST_CHECK_EQUAL(m1(2,3), -319);
                BOOST_CHECK_EQUAL(m1(3,2), 5248);
                BOOST_CHECK_EQUAL(m1(3,3), 88);
                
                BOOST_CHECK_EQUAL(m1(4,4), -13945);
                BOOST_CHECK_EQUAL(m1(4,5), -13105);
                BOOST_CHECK_EQUAL(m1(5,4), -432);
                BOOST_CHECK_EQUAL(m1(5,5), 986);
                
                BOOST_CHECK_EQUAL(m1(6,6), 1621);
                BOOST_CHECK_EQUAL(m1(6,7), 1377);
                BOOST_CHECK_EQUAL(m1(7,6), -6001);
                BOOST_CHECK_EQUAL(m1(7,7), -4546);
                
                BOOST_CHECK_EQUAL(m1(0,2), 0);
                BOOST_CHECK_EQUAL(m1(0,3), 0);
                BOOST_CHECK_EQUAL(m1(0,4), 0);
                BOOST_CHECK_EQUAL(m1(0,5), 0);
                BOOST_CHECK_EQUAL(m1(0,6), 0);
                BOOST_CHECK_EQUAL(m1(0,7), 0);
                
                BOOST_CHECK_EQUAL(m1(1,2), 0);
                BOOST_CHECK_EQUAL(m1(1,3), 0);
                BOOST_CHECK_EQUAL(m1(1,4), 0);
                BOOST_CHECK_EQUAL(m1(1,5), 0);
                BOOST_CHECK_EQUAL(m1(1,6), 0);
                BOOST_CHECK_EQUAL(m1(1,7), 0);
                
                BOOST_CHECK_EQUAL(m1(2,0), 0);
                BOOST_CHECK_EQUAL(m1(2,1), 0);
                BOOST_CHECK_EQUAL(m1(2,4), 0);
                BOOST_CHECK_EQUAL(m1(2,5), 0);
                BOOST_CHECK_EQUAL(m1(2,6), 0);
                BOOST_CHECK_EQUAL(m1(2,7), 0);
                
                BOOST_CHECK_EQUAL(m1(3,0), 0);
                BOOST_CHECK_EQUAL(m1(3,1), 0);
                BOOST_CHECK_EQUAL(m1(3,4), 0);
                BOOST_CHECK_EQUAL(m1(3,5), 0);
                BOOST_CHECK_EQUAL(m1(3,6), 0);
                BOOST_CHECK_EQUAL(m1(3,7), 0);

                BOOST_CHECK_EQUAL(m1(4,0), 0);
                BOOST_CHECK_EQUAL(m1(4,1), 0);
                BOOST_CHECK_EQUAL(m1(4,2), 0);
                BOOST_CHECK_EQUAL(m1(4,3), 0);
                BOOST_CHECK_EQUAL(m1(4,6), 0);
                BOOST_CHECK_EQUAL(m1(4,7), 0);
                
                BOOST_CHECK_EQUAL(m1(5,0), 0);
                BOOST_CHECK_EQUAL(m1(5,1), 0);
                BOOST_CHECK_EQUAL(m1(5,2), 0);
                BOOST_CHECK_EQUAL(m1(5,3), 0);
                BOOST_CHECK_EQUAL(m1(5,6), 0);
                BOOST_CHECK_EQUAL(m1(5,7), 0);
                
                BOOST_CHECK_EQUAL(m1(6,0), 0);
                BOOST_CHECK_EQUAL(m1(6,1), 0);
                BOOST_CHECK_EQUAL(m1(6,2), 0);
                BOOST_CHECK_EQUAL(m1(6,3), 0);
                BOOST_CHECK_EQUAL(m1(6,4), 0);
                BOOST_CHECK_EQUAL(m1(6,5), 0);
                                
                BOOST_CHECK_EQUAL(m1(7,0), 0);
                BOOST_CHECK_EQUAL(m1(7,1), 0);
                BOOST_CHECK_EQUAL(m1(7,2), 0);
                BOOST_CHECK_EQUAL(m1(7,3), 0);
                BOOST_CHECK_EQUAL(m1(7,4), 0);
                BOOST_CHECK_EQUAL(m1(7,5), 0);

            }*/
        }
        

        
        void testBlockMatrices()
        {
//             typedef NekMatrix<int, eFull, eBlock> BlockMatrix;
//             
//             {
//                 unsigned int rows[] = {8};
//                 unsigned int columns[] = {8};
//                 BlockMatrix M1(1, 1, rows, columns);
//                 BlockMatrix M2(1, 1, rows, columns);
//                 
//                 int m1_buf[] = {11, 93, -14, -99, -67, 68, 45, 76, 6, 72, -28, -61, -59, 6, -87, 72, -46, -68, -42, -47, -32, 37, -93, -58, -90, -53, -69, -84, 46, 59, -56, -83, -91, 92, -93, 91, -54, 10, -77, -63, -90, 61, -3, -82, 16, -40, 21, -94, -98, 75, 39, 95, -68, 98, -36, -95, 8, 92, 8, -95, -18, 44, 66, -62};
//                 int m2_buf[] = {-32, 78, 39, 94, 68, -17, -98, -36, 40, 22, 5, -88, -43, -73, 25, 4, -59, 62, -55, 25, 9, 40, 61, 40, -78, 62, 11, 88, 1, 30, 81, -5, -28, 4, -11, 10, 57, -82, -48, -11, 38, -7, 58, -94, -68, 14, -35, -14, -9, -51, -73, -73, -91, 1, 5, -86, 43, -4, -50, 50, 67, -39, 8, -49};
//                 
//                 NekMatrix<int> m1(8, 8, m1_buf);
//                 NekMatrix<int> m2(8, 8, m2_buf);
//                  
//                 M1(0, 0) = m1;
//                 M2(0, 0) = m2;
//                 
//                 BlockMatrix R = M1*M2;
//             }
//             
//             
//             unsigned int rowSizes[] = {3, 3, 3};
//             unsigned int columnSizes[] = {3, 3, 3};
//             
//             BlockMatrix m1(3, 3, rowSizes, columnSizes);

        }
    }
}


/**
    $Log: testNekMatrix.cpp,v $
    Revision 1.19  2007/01/16 05:31:34  bnelson
    Major improvements for expression templates.

    Revision 1.18  2006/10/30 05:08:13  bnelson
    Added preliminary linear system and block matrix support.

    Revision 1.17  2006/10/02 01:20:48  bnelson
    Started working on adding BLAS and LAPACK

    Revision 1.16  2006/09/30 15:38:29  bnelson
    no message

    Revision 1.15  2006/09/11 03:28:41  bnelson
    no message

    Revision 1.14  2006/08/25 01:38:59  bnelson
    no message

    Revision 1.13  2006/08/25 01:36:25  bnelson
    no message

    Revision 1.12  2006/08/14 02:35:45  bnelson
    Added many LinearAlgebra tests

    Revision 1.11  2006/06/05 02:23:17  bnelson
    Updates for the reorganization of LibUtilities.

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

