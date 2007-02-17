///////////////////////////////////////////////////////////////////////////////
//
// File: testLinearSystem.cpp
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
// Description: Test code for NekVector
//
///////////////////////////////////////////////////////////////////////////////

#include <UnitTests/testLinearSystem.h>
#include <LibUtilities/LinearAlgebra/NekLinSys.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

namespace Nektar
{
    namespace LinearSystemUnitTests
    {
        void testMixedInputParameterTypes()
        {
            {
                unsigned int matrix_buf[] = { 10, 5, 2 };

                unsigned int result_buf[] = { 20, 50, 10 };

                boost::shared_ptr<NekMatrix<unsigned int, eDiagonal> >  A(new NekMatrix<unsigned int, eDiagonal>(3,matrix_buf));
                NekVector<unsigned int, 3> b(result_buf);

                LinearSystem<NekMatrix<unsigned int, eDiagonal> > linsys(A);

                NekVector<unsigned int, 3> result = linsys.Solve(b);
                
                unsigned int expected_result_buf[] = { 2, 10, 5 };
                NekVector<unsigned int, 3> expectedResult(expected_result_buf);

                BOOST_CHECK_EQUAL(result, expectedResult);
            }

            {
                unsigned int matrix_buf[] = { 10, 5, 2 };

                unsigned int result_buf[] = { 20, 50, 10 };

                boost::shared_ptr<NekMatrix<unsigned int, eDiagonal> >  A(new NekMatrix<unsigned int, eDiagonal>(3,matrix_buf));
                boost::shared_ptr<NekVector<unsigned int, 3> > b(new NekVector<unsigned int, 3>(result_buf));

                LinearSystem<NekMatrix<unsigned int, eDiagonal> > linsys(A);

                NekVector<unsigned int, 3> result = linsys.Solve(b);
                
                unsigned int expected_result_buf[] = { 2, 10, 5 };
                NekVector<unsigned int, 3> expectedResult(expected_result_buf);

                BOOST_CHECK_EQUAL(result, expectedResult);
            }

            {
                unsigned int matrix_buf[] = { 10, 5, 2 };

                unsigned int result_buf[] = { 20, 50, 10 };

                boost::shared_ptr<NekMatrix<unsigned int, eDiagonal> >  A(new NekMatrix<unsigned int, eDiagonal>(3,matrix_buf));
                boost::shared_ptr<NekVector<unsigned int, 3> > b(new NekVector<unsigned int, 3>(result_buf));
                boost::shared_ptr<NekVector<unsigned int, 3> > result(new NekVector<unsigned int, 3>());

                LinearSystem<NekMatrix<unsigned int, eDiagonal> > linsys(A);

                linsys.Solve(b,result);
                
                unsigned int expected_result_buf[] = { 2, 10, 5 };
                NekVector<unsigned int, 3> expectedResult(expected_result_buf);

                BOOST_CHECK_EQUAL(*result, expectedResult);
            }

            {
                unsigned int matrix_buf[] = { 10, 5, 2 };

                unsigned int result_buf[] = { 20, 50, 10 };

                boost::shared_ptr<NekMatrix<unsigned int, eDiagonal> >  A(new NekMatrix<unsigned int, eDiagonal>(3,matrix_buf));
                NekVector<unsigned int, 3> b(result_buf);
                NekVector<unsigned int, 3> result;

                LinearSystem<NekMatrix<unsigned int, eDiagonal> > linsys(A);

                linsys.Solve(b,result);
                
                unsigned int expected_result_buf[] = { 2, 10, 5 };
                NekVector<unsigned int, 3> expectedResult(expected_result_buf);

                BOOST_CHECK_EQUAL(result, expectedResult);
            }

            {
                unsigned int matrix_buf[] = { 10, 5, 2 };

                unsigned int result_buf[] = { 20, 50, 10 };

                boost::shared_ptr<NekMatrix<unsigned int, eDiagonal> >  A(new NekMatrix<unsigned int, eDiagonal>(3,matrix_buf));
                boost::shared_ptr<NekVector<unsigned int, 3> > b(new NekVector<unsigned int, 3>(result_buf));
                NekVector<unsigned int, 3> result;

                LinearSystem<NekMatrix<unsigned int, eDiagonal> > linsys(A);

                linsys.Solve(b,result);
                
                unsigned int expected_result_buf[] = { 2, 10, 5 };
                NekVector<unsigned int, 3> expectedResult(expected_result_buf);

                BOOST_CHECK_EQUAL(result, expectedResult);
            }

            {
                unsigned int matrix_buf[] = { 10, 5, 2 };

                unsigned int result_buf[] = { 20, 50, 10 };

                boost::shared_ptr<NekMatrix<unsigned int, eDiagonal> >  A(new NekMatrix<unsigned int, eDiagonal>(3,matrix_buf));
                NekVector<unsigned int, 3> b(result_buf);
                boost::shared_ptr<NekVector<unsigned int, 3> > result(new NekVector<unsigned int, 3>());

                LinearSystem<NekMatrix<unsigned int, eDiagonal> > linsys(A);

                linsys.Solve(b,result);
                
                unsigned int expected_result_buf[] = { 2, 10, 5 };
                NekVector<unsigned int, 3> expectedResult(expected_result_buf);

                BOOST_CHECK_EQUAL(*result, expectedResult);
            }
        }

        void testDiagonalSystem()
        {
            unsigned int matrix_buf[] = { 10, 5, 2 };

            unsigned int result_buf[] = { 20, 50, 10 };

            boost::shared_ptr<NekMatrix<unsigned int, eDiagonal> >  A(new NekMatrix<unsigned int, eDiagonal>(3,matrix_buf));
            boost::shared_ptr<NekVector<unsigned int, 3> > b(new NekVector<unsigned int, 3>(result_buf));

            LinearSystem<NekMatrix<unsigned int, eDiagonal> > linsys(A);

            NekVector<unsigned int, 3> result = linsys.Solve(b);
            
            unsigned int expected_result_buf[] = { 2, 10, 5 };
            NekVector<unsigned int, 3> expectedResult(expected_result_buf);

            BOOST_CHECK_EQUAL(result, expectedResult);
        }
        
        void testFullSystem()
        {
            {
                double matrix_buf[] = { 10, 0, 0,
                                        0, 5, 0, 
                                        0, 0, 2 };
    
                double result_buf[] = { 20, 50, 10 };
    
                boost::shared_ptr<NekMatrix<double, eFull> >  A(new NekMatrix<double, eFull>(3,3,matrix_buf));
                boost::shared_ptr<NekVector<double, 3> > b(new NekVector<double, 3>(result_buf));
    
                LinearSystem<NekMatrix<double, eFull> > linsys(A);
    
                NekVector<double, 3> result = linsys.Solve(b);
                
                unsigned int expected_result_buf[] = { 2, 10, 5 };
                NekVector<unsigned int, 3> expectedResult(expected_result_buf);
    
                double epsilon = 1e-14;
                
                BOOST_CHECK_CLOSE(result[0], 2.0, epsilon);
                BOOST_CHECK_CLOSE(result[1], 10.0, epsilon);
                BOOST_CHECK_CLOSE(result[2], 5.0, epsilon);
            }
            
            {
                double matrix_buf[] = {81, -5, 
                                        -28, 4};
                            
                double b_buf[] = {-941, 348};
                
                boost::shared_ptr<NekMatrix<double, eFull> > A(new NekMatrix<double, eFull>(2, 2, matrix_buf));
                boost::shared_ptr<NekVector<double> > b(new NekVector<double>(2, b_buf));
                LinearSystem<NekMatrix<double, eFull> > linsys(A);
            
                NekVector<double> result = linsys.Solve(b);
                double epsilon = 1e-11;
                BOOST_CHECK_CLOSE(result[0], -11.0, epsilon);
                BOOST_CHECK_CLOSE(result[1], 10.0, epsilon);
            }
            
            
            {
                // Larger matrix.
                double matrix_buf[] = {-85, -55, -37, -35, 97, 50, 79, 56, 49, 63, 
                    57, -59, 45, -8, -93, 92, 43, -62, 77, 66,
                    54, -5, 99, -61, -50, -12, -18, 31, -26, -62,
                    1, -47, -91, -47, -61, 41, -58, -90, 53, -1,
                    94, 83, -86, 23, -84, 19, -50, 88, -53, 85,
                    49, 78, 17, 72, -99, -85, -86, 30, 80, 72,
                    66, -29, -91, -53, -19, -47, 68, -72, -87, 79,
                    43, -66, -53, -61, -23, -37, 31, -34, -42, 88,
                    -76, -65, 25, 28, -61, -60, 9, 29, -66, -32,
                    78, 39, 94, 68, -17, -98, -36, 40, 22, 5 };
                            
                double b_buf[] = {12719, -3169, -16810, 7408, -14945, -6822, 10166, 7023, 8679, -11826};
                
                boost::shared_ptr<NekMatrix<double, eFull> > A(new NekMatrix<double, eFull>(10, 10, matrix_buf));
                boost::shared_ptr<NekVector<double> > b(new NekVector<double>(10, b_buf));
                LinearSystem<NekMatrix<double, eFull> > linsys(A);
                
                NekVector<double> result = linsys.Solve(b);
                double epsilon = 1e-11;
                BOOST_CHECK_CLOSE(result[0], -88.0, epsilon);
                BOOST_CHECK_CLOSE(result[1], -43.0, epsilon);
                BOOST_CHECK_CLOSE(result[2], -73.0, epsilon);
                BOOST_CHECK_CLOSE(result[3], 25.0, epsilon);
                BOOST_CHECK_CLOSE(result[4], 4.0, epsilon);
                BOOST_CHECK_CLOSE(result[5], -59.0, epsilon);
                BOOST_CHECK_CLOSE(result[6], 62.0, epsilon);
                BOOST_CHECK_CLOSE(result[7], -55.0, epsilon);
                BOOST_CHECK_CLOSE(result[8], 25.0, epsilon);
                BOOST_CHECK_CLOSE(result[9], 9.0, epsilon);
                    

            }
            
            {
                //boost::shared_ptr<NekMatrix<double, eFull, eBlock> > A(new NekMatrix<double, eFull, eBlock>(10, 10));
                //double b_buf[] = {12719, -3169, -16810, 7408, -14945, -6822, 10166, 7023, 8679, -11826};
                //boost::shared_ptr<NekVector<double> > b(new NekVector<double>(10, b_buf));
            }
        }
        
        void testSolvingBlockDiagonalMatrices()
        {
            //typedef NekMatrix<double, eDiagonal, eBlock> BlockMatrix;
            //typedef BlockMatrix::InnerMatrixType InnerMatrixType;
            //    
            //boost::shared_ptr<BlockMatrix> m1(new BlockMatrix(4, 2, 2));
            //const double m1_buf0[] = {-91, -47, 94, 83};
            //const double m1_buf1[] = {49, 78, 80, 72};
            //const double m1_buf2[] = {87, 79, 31, -34};
            //const double m1_buf3[] = {9, 29, -17, -98};
            //    
            //m1->GetBlock(0,0) = InnerMatrixType(2, 2, m1_buf0);
            //m1->GetBlock(1,1) = InnerMatrixType(2, 2, m1_buf1);
            //m1->GetBlock(2,2) = InnerMatrixType(2, 2, m1_buf2);
            //m1->GetBlock(3,3) = InnerMatrixType(2, 2, m1_buf3);
            //
            //double b_vals[] = {-7198, 7642, 1344, 3744, -9968, 115, -2469, 8424};
            //boost::shared_ptr<NekVector<double> > b(new NekVector<double>(8, b_vals));
            //
            //LinearSystem<BlockMatrix, NekVector<double> > linsys(m1, b);
            //
            //NekVector<double> result = linsys.Solve();
            //double epsilon = 1e-11;
            //BOOST_CHECK_CLOSE(result[0], 76.0, epsilon);
            //BOOST_CHECK_CLOSE(result[1], 6.0, epsilon);
            //BOOST_CHECK_CLOSE(result[2], 72.0, epsilon);
            //BOOST_CHECK_CLOSE(result[3], -28.0, epsilon);
            //BOOST_CHECK_CLOSE(result[4], -61.0, epsilon);
            //BOOST_CHECK_CLOSE(result[5], -59.0, epsilon);
            //BOOST_CHECK_CLOSE(result[6], 6.0, epsilon);
            //BOOST_CHECK_CLOSE(result[7], -87.0, epsilon);
        }
    }
}
