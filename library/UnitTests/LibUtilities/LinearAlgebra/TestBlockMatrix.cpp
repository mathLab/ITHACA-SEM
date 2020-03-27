///////////////////////////////////////////////////////////////////////////////
//
// File: TestBlockMatrix.cpp
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

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

namespace Nektar
{
    namespace BlockMatrixUnitTests
    {
        BOOST_AUTO_TEST_CASE(TestRowsAndColumnsPerBlockAccess)
        {
            NekMatrix<NekMatrix<double>, BlockMatrixTag> m(2, 3, 7, 10);
            BOOST_CHECK_EQUAL(7, m.GetNumberOfRowsInBlockRow(0));
            BOOST_CHECK_EQUAL(7, m.GetNumberOfRowsInBlockRow(1));
            BOOST_CHECK_EQUAL(10, m.GetNumberOfColumnsInBlockColumn(0));
            BOOST_CHECK_EQUAL(10, m.GetNumberOfColumnsInBlockColumn(1));
            BOOST_CHECK_EQUAL(10, m.GetNumberOfColumnsInBlockColumn(2));
        }

        BOOST_AUTO_TEST_CASE(TestMultiplication)
        {
            typedef NekMatrix<double> InnerType;
            typedef NekMatrix<InnerType, BlockMatrixTag> BlockType;

            double m1_buf[] = {1, 6,
                               2, 7,
                               3, 8};
            double m2_buf[] = {4, 9,
                               5, 10};
            double m3_buf[] = {11, 16,
                               12, 17,
                               13, 18};
            double m4_buf[] = {14, 19,
                               15, 20};

            std::shared_ptr<InnerType> m1(new InnerType(2, 3, m1_buf));
            std::shared_ptr<InnerType> m2(new InnerType(2, 2, m2_buf));
            std::shared_ptr<InnerType> m3(new InnerType(2, 3, m3_buf));
            std::shared_ptr<InnerType> m4(new InnerType(2, 2, m4_buf));

            unsigned int rowCounts[] = {2, 2};
            unsigned int colCounts[] = {3, 2};

            BlockType b(2, 2, rowCounts, colCounts);
            b.SetBlock(0, 0, m1);
            b.SetBlock(0, 1, m2);
            b.SetBlock(1, 0, m3);
            b.SetBlock(1, 1, m4);

            double rhs_buf[] = {10, 20, 30, 40, 50};
            NekVector<double> rhs(5, rhs_buf);

            NekVector<double> result = b*rhs;

            double expected_buf[] = {550, 1300, 2050, 2800};
            NekVector<double> expected_result(4, expected_buf);

            BOOST_CHECK_EQUAL(expected_result, result);
        }

        BOOST_AUTO_TEST_CASE(TestBlockMatrixTranspose)
        {
            typedef NekMatrix<double> InnerType;
            typedef NekMatrix<InnerType, BlockMatrixTag> BlockType;

            double m00_buf[] = {1, 2,
                                3, 4,
                                5, 6,
                                7, 8};
            double m01_buf[] = {9, 10,
                                11, 12,
                                13, 14,
                                15, 16,
                                17, 18};
            double m02_buf[] = {19, 20,
                                21, 22,
                                23, 24,
                                25, 26,
                                27, 28,
                                29, 30};
            double m10_buf[] = {31, 32, 33,
                                34, 35, 36,
                                37, 38, 39,
                                40, 41, 42};
            double m11_buf[] = {43, 44, 45,
                                46, 47, 48,
                                49, 50, 51,
                                52, 53, 54,
                                55, 56, 57};
            double m12_buf[] = { 58, 59, 60,
                                61, 62, 63,
                                64, 65, 66,
                                67, 68, 69,
                                70, 71, 72,
                                73, 74, 75};

            std::shared_ptr<InnerType> m00(new InnerType(2, 4, m00_buf));
            std::shared_ptr<InnerType> m01(new InnerType(2, 5, m01_buf));
            std::shared_ptr<InnerType> m02(new InnerType(2, 6, m02_buf));
            std::shared_ptr<InnerType> m10(new InnerType(3, 4, m10_buf));
            std::shared_ptr<InnerType> m11(new InnerType(3, 5, m11_buf));
            std::shared_ptr<InnerType> m12(new InnerType(3, 6, m12_buf));

            unsigned int rowCounts[] = {2, 3};
            unsigned int colCounts[] = {4, 5, 6};

            BlockType b(2, 3, rowCounts, colCounts);
            b.SetBlock(0, 0, m00);
            b.SetBlock(0, 1, m01);
            b.SetBlock(0, 2, m02);
            b.SetBlock(1, 0, m10);
            b.SetBlock(1, 1, m11);
            b.SetBlock(1, 2, m12);

            BOOST_CHECK_EQUAL(2, b.GetNumberOfBlockRows());
            BOOST_CHECK_EQUAL(3, b.GetNumberOfBlockColumns());
            BOOST_CHECK_EQUAL(2, b.GetNumberOfRowsInBlockRow(0));
            BOOST_CHECK_EQUAL(3, b.GetNumberOfRowsInBlockRow(1));
            BOOST_CHECK_EQUAL(4, b.GetNumberOfColumnsInBlockColumn(0));
            BOOST_CHECK_EQUAL(5, b.GetNumberOfColumnsInBlockColumn(1));
            BOOST_CHECK_EQUAL(6, b.GetNumberOfColumnsInBlockColumn(2));

            BOOST_CHECK_EQUAL(0, b.CalculateBlockIndex(0,0));
            BOOST_CHECK_EQUAL(1, b.CalculateBlockIndex(1,0));
            BOOST_CHECK_EQUAL(2, b.CalculateBlockIndex(0,1));
            BOOST_CHECK_EQUAL(3, b.CalculateBlockIndex(1,1));
            BOOST_CHECK_EQUAL(4, b.CalculateBlockIndex(0,2));
            BOOST_CHECK_EQUAL(5, b.CalculateBlockIndex(1,2));

            BOOST_CHECK_EQUAL(m00->GetRawPtr(), b.GetBlock(0,0)->GetRawPtr());
            BOOST_CHECK_EQUAL(m01->GetRawPtr(), b.GetBlock(0,1)->GetRawPtr());
            BOOST_CHECK_EQUAL(m02->GetRawPtr(), b.GetBlock(0,2)->GetRawPtr());
            BOOST_CHECK_EQUAL(m10->GetRawPtr(), b.GetBlock(1,0)->GetRawPtr());
            BOOST_CHECK_EQUAL(m11->GetRawPtr(), b.GetBlock(1,1)->GetRawPtr());
            BOOST_CHECK_EQUAL(m12->GetRawPtr(), b.GetBlock(1,2)->GetRawPtr());

            BOOST_CHECK_EQUAL('N', b.GetBlock(0,0)->GetTransposeFlag());
            BOOST_CHECK_EQUAL('N', b.GetBlock(0,1)->GetTransposeFlag());
            BOOST_CHECK_EQUAL('N', b.GetBlock(0,2)->GetTransposeFlag());
            BOOST_CHECK_EQUAL('N', b.GetBlock(1,0)->GetTransposeFlag());
            BOOST_CHECK_EQUAL('N', b.GetBlock(1,1)->GetTransposeFlag());
            BOOST_CHECK_EQUAL('N', b.GetBlock(1,2)->GetTransposeFlag());

            double v_buf[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
            NekVector<double> v(15, v_buf);
            NekVector<double> result = b*v;

            BOOST_CHECK_EQUAL(5, result.GetRows());
            BOOST_CHECK_EQUAL(2360, result[0]);
            BOOST_CHECK_EQUAL(2480, result[1]);
            BOOST_CHECK_EQUAL(7080, result[2]);
            BOOST_CHECK_EQUAL(7200, result[3]);
            BOOST_CHECK_EQUAL(7320, result[4]);



            b.Transpose();

            BOOST_CHECK_EQUAL(3, b.GetNumberOfBlockRows());
            BOOST_CHECK_EQUAL(2, b.GetNumberOfBlockColumns());

            BOOST_CHECK_EQUAL(2, b.GetNumberOfColumnsInBlockColumn(0));
            BOOST_CHECK_EQUAL(3, b.GetNumberOfColumnsInBlockColumn(1));
            BOOST_CHECK_EQUAL(4, b.GetNumberOfRowsInBlockRow(0));
            BOOST_CHECK_EQUAL(5, b.GetNumberOfRowsInBlockRow(1));
            BOOST_CHECK_EQUAL(6, b.GetNumberOfRowsInBlockRow(2));

            BOOST_CHECK_EQUAL(0, b.CalculateBlockIndex(0,0));
            BOOST_CHECK_EQUAL(1, b.CalculateBlockIndex(0,1));
            BOOST_CHECK_EQUAL(2, b.CalculateBlockIndex(1,0));
            BOOST_CHECK_EQUAL(3, b.CalculateBlockIndex(1,1));
            BOOST_CHECK_EQUAL(4, b.CalculateBlockIndex(2,0));
            BOOST_CHECK_EQUAL(5, b.CalculateBlockIndex(2,1));

            BOOST_CHECK_EQUAL(m00->GetRawPtr(), b.GetBlock(0,0)->GetRawPtr());
            BOOST_CHECK_EQUAL(m01->GetRawPtr(), b.GetBlock(1,0)->GetRawPtr());
            BOOST_CHECK_EQUAL(m02->GetRawPtr(), b.GetBlock(2,0)->GetRawPtr());
            BOOST_CHECK_EQUAL(m10->GetRawPtr(), b.GetBlock(0,1)->GetRawPtr());
            BOOST_CHECK_EQUAL(m11->GetRawPtr(), b.GetBlock(1,1)->GetRawPtr());
            BOOST_CHECK_EQUAL(m12->GetRawPtr(), b.GetBlock(2,1)->GetRawPtr());

            BOOST_CHECK_EQUAL('T', b.GetBlock(0,0)->GetTransposeFlag());
            BOOST_CHECK_EQUAL('T', b.GetBlock(1,0)->GetTransposeFlag());
            BOOST_CHECK_EQUAL('T', b.GetBlock(2,0)->GetTransposeFlag());
            BOOST_CHECK_EQUAL('T', b.GetBlock(0,1)->GetTransposeFlag());
            BOOST_CHECK_EQUAL('T', b.GetBlock(1,1)->GetTransposeFlag());
            BOOST_CHECK_EQUAL('T', b.GetBlock(2,1)->GetTransposeFlag());

            double v_transpose_buf[] = {1, 2, 3, 4, 5};
            NekVector<double> v_transpose(5, v_transpose_buf);
            NekVector<double> result_transpose = b*v_transpose;

            BOOST_CHECK_EQUAL(15, result_transpose.GetRows());
            BOOST_CHECK_EQUAL(391, result_transpose[0]);
            BOOST_CHECK_EQUAL(433, result_transpose[1]);
            BOOST_CHECK_EQUAL(475, result_transpose[2]);
            BOOST_CHECK_EQUAL(517, result_transpose[3]);
            BOOST_CHECK_EQUAL(559, result_transpose[4]);
            BOOST_CHECK_EQUAL(601, result_transpose[5]);
            BOOST_CHECK_EQUAL(643, result_transpose[6]);
            BOOST_CHECK_EQUAL(685, result_transpose[7]);
            BOOST_CHECK_EQUAL(727, result_transpose[8]);
            BOOST_CHECK_EQUAL(769, result_transpose[9]);
            BOOST_CHECK_EQUAL(811, result_transpose[10]);
            BOOST_CHECK_EQUAL(853, result_transpose[11]);
            BOOST_CHECK_EQUAL(895, result_transpose[12]);
            BOOST_CHECK_EQUAL(937, result_transpose[13]);
            BOOST_CHECK_EQUAL(979, result_transpose[14]);

            BlockType globalTranspose = Transpose(b);
            NekVector<double> globalResult = globalTranspose*v;
            BOOST_CHECK_EQUAL(5, globalResult.GetRows());
            BOOST_CHECK_EQUAL(2360, globalResult[0]);
            BOOST_CHECK_EQUAL(2480, globalResult[1]);
            BOOST_CHECK_EQUAL(7080, globalResult[2]);
            BOOST_CHECK_EQUAL(7200, globalResult[3]);
            BOOST_CHECK_EQUAL(7320, globalResult[4]);
            //
        }


        BOOST_AUTO_TEST_CASE(TestBlockScaledMatrix)
        {
            typedef NekMatrix<NekMatrix<double>, ScaledMatrixTag> InnerType;
            typedef NekMatrix<InnerType, BlockMatrixTag> BlockType;

            double m00_buf[] = {1, 2,
                                3, 4,
                                5, 6,
                                7, 8};
            double m01_buf[] = {9, 10,
                                11, 12,
                                13, 14,
                                15, 16,
                                17, 18};
            double m02_buf[] = {19, 20,
                                21, 22,
                                23, 24,
                                25, 26,
                                27, 28,
                                29, 30};
            double m10_buf[] = {31, 32, 33,
                                34, 35, 36,
                                37, 38, 39,
                                40, 41, 42};
            double m11_buf[] = {43, 44, 45,
                                46, 47, 48,
                                49, 50, 51,
                                52, 53, 54,
                                55, 56, 57};
            double m12_buf[] = { 58, 59, 60,
                                61, 62, 63,
                                64, 65, 66,
                                67, 68, 69,
                                70, 71, 72,
                                73, 74, 75};

            std::shared_ptr<NekMatrix<double> > m00(new NekMatrix<double>(2, 4, m00_buf));
            std::shared_ptr<NekMatrix<double> > m01(new NekMatrix<double>(2, 5, m01_buf));
            std::shared_ptr<NekMatrix<double> > m02(new NekMatrix<double>(2, 6, m02_buf));
            std::shared_ptr<NekMatrix<double> > m10(new NekMatrix<double>(3, 4, m10_buf));
            std::shared_ptr<NekMatrix<double> > m11(new NekMatrix<double>(3, 5, m11_buf));
            std::shared_ptr<NekMatrix<double> > m12(new NekMatrix<double>(3, 6, m12_buf));

            std::shared_ptr<InnerType> sm00(new InnerType(2.0, m00));
            std::shared_ptr<InnerType> sm01(new InnerType(2.0, m01));
            std::shared_ptr<InnerType> sm02(new InnerType(2.0, m02));
            std::shared_ptr<InnerType> sm10(new InnerType(2.0, m10));
            std::shared_ptr<InnerType> sm11(new InnerType(2.0, m11));
            std::shared_ptr<InnerType> sm12(new InnerType(2.0, m12));

            unsigned int rowCounts[] = {2, 3};
            unsigned int colCounts[] = {4, 5, 6};

            BlockType b(2, 3, rowCounts, colCounts);
            b.SetBlock(0, 0, sm00);
            b.SetBlock(0, 1, sm01);
            b.SetBlock(0, 2, sm02);
            b.SetBlock(1, 0, sm10);
            b.SetBlock(1, 1, sm11);
            b.SetBlock(1, 2, sm12);

            BOOST_CHECK_EQUAL(2, b.GetNumberOfBlockRows());
            BOOST_CHECK_EQUAL(3, b.GetNumberOfBlockColumns());
            BOOST_CHECK_EQUAL(2, b.GetNumberOfRowsInBlockRow(0));
            BOOST_CHECK_EQUAL(3, b.GetNumberOfRowsInBlockRow(1));
            BOOST_CHECK_EQUAL(4, b.GetNumberOfColumnsInBlockColumn(0));
            BOOST_CHECK_EQUAL(5, b.GetNumberOfColumnsInBlockColumn(1));
            BOOST_CHECK_EQUAL(6, b.GetNumberOfColumnsInBlockColumn(2));

            BOOST_CHECK_EQUAL(0, b.CalculateBlockIndex(0,0));
            BOOST_CHECK_EQUAL(1, b.CalculateBlockIndex(1,0));
            BOOST_CHECK_EQUAL(2, b.CalculateBlockIndex(0,1));
            BOOST_CHECK_EQUAL(3, b.CalculateBlockIndex(1,1));
            BOOST_CHECK_EQUAL(4, b.CalculateBlockIndex(0,2));
            BOOST_CHECK_EQUAL(5, b.CalculateBlockIndex(1,2));

            BOOST_CHECK_EQUAL(m00->GetRawPtr(), b.GetBlock(0,0)->GetRawPtr());
            BOOST_CHECK_EQUAL(m01->GetRawPtr(), b.GetBlock(0,1)->GetRawPtr());
            BOOST_CHECK_EQUAL(m02->GetRawPtr(), b.GetBlock(0,2)->GetRawPtr());
            BOOST_CHECK_EQUAL(m10->GetRawPtr(), b.GetBlock(1,0)->GetRawPtr());
            BOOST_CHECK_EQUAL(m11->GetRawPtr(), b.GetBlock(1,1)->GetRawPtr());
            BOOST_CHECK_EQUAL(m12->GetRawPtr(), b.GetBlock(1,2)->GetRawPtr());

            BOOST_CHECK_EQUAL('N', b.GetBlock(0,0)->GetTransposeFlag());
            BOOST_CHECK_EQUAL('N', b.GetBlock(0,1)->GetTransposeFlag());
            BOOST_CHECK_EQUAL('N', b.GetBlock(0,2)->GetTransposeFlag());
            BOOST_CHECK_EQUAL('N', b.GetBlock(1,0)->GetTransposeFlag());
            BOOST_CHECK_EQUAL('N', b.GetBlock(1,1)->GetTransposeFlag());
            BOOST_CHECK_EQUAL('N', b.GetBlock(1,2)->GetTransposeFlag());

            double v_buf[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
            NekVector<double> v(15, v_buf);
            NekVector<double> result = b*v;

            BOOST_CHECK_EQUAL(5, result.GetRows());
            BOOST_CHECK_EQUAL(2360*2, result[0]);
            BOOST_CHECK_EQUAL(2480*2, result[1]);
            BOOST_CHECK_EQUAL(7080*2, result[2]);
            BOOST_CHECK_EQUAL(7200*2, result[3]);
            BOOST_CHECK_EQUAL(7320*2, result[4]);



            b.Transpose();

            BOOST_CHECK_EQUAL(3, b.GetNumberOfBlockRows());
            BOOST_CHECK_EQUAL(2, b.GetNumberOfBlockColumns());

            BOOST_CHECK_EQUAL(2, b.GetNumberOfColumnsInBlockColumn(0));
            BOOST_CHECK_EQUAL(3, b.GetNumberOfColumnsInBlockColumn(1));
            BOOST_CHECK_EQUAL(4, b.GetNumberOfRowsInBlockRow(0));
            BOOST_CHECK_EQUAL(5, b.GetNumberOfRowsInBlockRow(1));
            BOOST_CHECK_EQUAL(6, b.GetNumberOfRowsInBlockRow(2));

            BOOST_CHECK_EQUAL(0, b.CalculateBlockIndex(0,0));
            BOOST_CHECK_EQUAL(1, b.CalculateBlockIndex(0,1));
            BOOST_CHECK_EQUAL(2, b.CalculateBlockIndex(1,0));
            BOOST_CHECK_EQUAL(3, b.CalculateBlockIndex(1,1));
            BOOST_CHECK_EQUAL(4, b.CalculateBlockIndex(2,0));
            BOOST_CHECK_EQUAL(5, b.CalculateBlockIndex(2,1));

            BOOST_CHECK_EQUAL(m00->GetRawPtr(), b.GetBlock(0,0)->GetRawPtr());
            BOOST_CHECK_EQUAL(m01->GetRawPtr(), b.GetBlock(1,0)->GetRawPtr());
            BOOST_CHECK_EQUAL(m02->GetRawPtr(), b.GetBlock(2,0)->GetRawPtr());
            BOOST_CHECK_EQUAL(m10->GetRawPtr(), b.GetBlock(0,1)->GetRawPtr());
            BOOST_CHECK_EQUAL(m11->GetRawPtr(), b.GetBlock(1,1)->GetRawPtr());
            BOOST_CHECK_EQUAL(m12->GetRawPtr(), b.GetBlock(2,1)->GetRawPtr());

            BOOST_CHECK_EQUAL('T', b.GetBlock(0,0)->GetTransposeFlag());
            BOOST_CHECK_EQUAL('T', b.GetBlock(1,0)->GetTransposeFlag());
            BOOST_CHECK_EQUAL('T', b.GetBlock(2,0)->GetTransposeFlag());
            BOOST_CHECK_EQUAL('T', b.GetBlock(0,1)->GetTransposeFlag());
            BOOST_CHECK_EQUAL('T', b.GetBlock(1,1)->GetTransposeFlag());
            BOOST_CHECK_EQUAL('T', b.GetBlock(2,1)->GetTransposeFlag());

            double v_transpose_buf[] = {1, 2, 3, 4, 5};
            NekVector<double> v_transpose(5, v_transpose_buf);
            NekVector<double> result_transpose = b*v_transpose;

            BOOST_CHECK_EQUAL(15, result_transpose.GetRows());
            BOOST_CHECK_EQUAL(391*2, result_transpose[0]);
            BOOST_CHECK_EQUAL(433*2, result_transpose[1]);
            BOOST_CHECK_EQUAL(475*2, result_transpose[2]);
            BOOST_CHECK_EQUAL(517*2, result_transpose[3]);
            BOOST_CHECK_EQUAL(559*2, result_transpose[4]);
            BOOST_CHECK_EQUAL(601*2, result_transpose[5]);
            BOOST_CHECK_EQUAL(643*2, result_transpose[6]);
            BOOST_CHECK_EQUAL(685*2, result_transpose[7]);
            BOOST_CHECK_EQUAL(727*2, result_transpose[8]);
            BOOST_CHECK_EQUAL(769*2, result_transpose[9]);
            BOOST_CHECK_EQUAL(811*2, result_transpose[10]);
            BOOST_CHECK_EQUAL(853*2, result_transpose[11]);
            BOOST_CHECK_EQUAL(895*2, result_transpose[12]);
            BOOST_CHECK_EQUAL(937*2, result_transpose[13]);
            BOOST_CHECK_EQUAL(979*2, result_transpose[14]);
        }

        BOOST_AUTO_TEST_CASE(TestBlockBlockMatrixTranspose)
        {
            typedef NekMatrix<double> Matrix;
            typedef NekMatrix<Matrix, BlockMatrixTag> InnerMatrix;
            typedef NekMatrix<InnerMatrix, BlockMatrixTag> BlockMatrix;

            double m_00_buf[] = { 1, 2,
                                  3, 4,
                                  5, 6};
            double m_01_buf[] = {7, 8,
                                 9, 10,
                                11, 12,
                                13, 14};
            double m_10_buf[] = {15, 16,
                                 17, 18,
                                 19, 20};
            double m_11_buf[] = {21, 22,
                                 23, 24,
                                 25, 26,
                                 27, 28};
            double m_20_buf[] = {29, 30,
                                 31, 32,
                                 33, 34};
            double m_21_buf[] = {35, 36,
                                 37, 38,
                                 39, 40,
                                 41, 42};
            std::shared_ptr<Matrix> m00(new Matrix(2, 3, m_00_buf));
            std::shared_ptr<Matrix> m01(new Matrix(2, 4, m_01_buf));
            std::shared_ptr<Matrix> m10(new Matrix(2, 3, m_10_buf));
            std::shared_ptr<Matrix> m11(new Matrix(2, 4, m_11_buf));
            std::shared_ptr<Matrix> m20(new Matrix(2, 3, m_20_buf));
            std::shared_ptr<Matrix> m21(new Matrix(2, 4, m_21_buf));

            unsigned int rowCounts[] = {2, 2, 2};
            unsigned int colCounts[] = {3, 4};

            std::vector<std::shared_ptr<InnerMatrix> > innerMatrices;

            for(int i = 0; i < 6; ++i)
            {
                std::shared_ptr<Matrix> t00(new Matrix(*m00));
                std::shared_ptr<Matrix> t01(new Matrix(*m01));
                std::shared_ptr<Matrix> t10(new Matrix(*m10));
                std::shared_ptr<Matrix> t11(new Matrix(*m11));
                std::shared_ptr<Matrix> t20(new Matrix(*m20));
                std::shared_ptr<Matrix> t21(new Matrix(*m21));

                (*t00) *= (i+1);
                (*t01) *= (i+1);
                (*t10) *= (i+1);
                (*t11) *= (i+1);
                (*t20) *= (i+1);
                (*t21) *= (i+1);

                std::shared_ptr<InnerMatrix> matrix(new InnerMatrix(3, 2, rowCounts, colCounts));
                matrix->SetBlock(0, 0, t00);
                matrix->SetBlock(0, 1, t01);
                matrix->SetBlock(1, 0, t10);
                matrix->SetBlock(1, 1, t11);
                matrix->SetBlock(2, 0, t20);
                matrix->SetBlock(2, 1, t21);

                innerMatrices.push_back(matrix);
            }

            unsigned int topLevelRowCounts[] = {6, 6};
            unsigned int topLevelColumnCounts[] = {7, 7, 7};
            BlockMatrix topLevelMatrix(2, 3, topLevelRowCounts, topLevelColumnCounts);
            topLevelMatrix.SetBlock(0, 0, innerMatrices[0]);
            topLevelMatrix.SetBlock(0, 1, innerMatrices[1]);
            topLevelMatrix.SetBlock(0, 2, innerMatrices[2]);
            topLevelMatrix.SetBlock(1, 0, innerMatrices[3]);
            topLevelMatrix.SetBlock(1, 1, innerMatrices[4]);
            topLevelMatrix.SetBlock(1, 2, innerMatrices[5]);

            double v_buf[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
            NekVector<double> v(21, v_buf);

            NekVector<double> result = topLevelMatrix*v;

            double expected_result_buf[] = {4256, 4816, 12096, 12656, 19936, 20496, 9611, 10864, 27153, 28406, 44695, 45948};
            NekVector<double> expected_result(12, expected_result_buf);

            BOOST_CHECK_EQUAL(expected_result, result);


            topLevelMatrix.Transpose();

            double transposed_v_buf[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
            NekVector<double> transposed_v(12, transposed_v_buf);
            NekVector<double> transposed_result = topLevelMatrix*transposed_v;

            double expected_transpose_result_buf[] = {4427, 4925, 5423, 5921, 6419, 6917, 7415, 5863, 6517,
                7171, 7825, 8479, 9133, 9787, 7299, 8109, 8919, 9729, 10539, 11349, 12159};
            NekVector<double> expected_transpose_result(21, expected_transpose_result_buf);

            BOOST_CHECK_EQUAL(expected_transpose_result, transposed_result);

        }

        BOOST_AUTO_TEST_CASE(TestBlockBlockScaledMatrixTranspose)
        {
            typedef NekMatrix<double> Matrix;
            typedef NekMatrix<NekMatrix<double>, ScaledMatrixTag>  ScaledMatrix;
            typedef NekMatrix<ScaledMatrix, BlockMatrixTag> InnerMatrix;
            typedef NekMatrix<InnerMatrix, BlockMatrixTag> BlockMatrix;

            double m_00_buf[] = { 1, 2,
                                  3, 4,
                                  5, 6};
            double m_01_buf[] = {7, 8,
                                 9, 10,
                                11, 12,
                                13, 14};
            double m_10_buf[] = {15, 16,
                                 17, 18,
                                 19, 20};
            double m_11_buf[] = {21, 22,
                                 23, 24,
                                 25, 26,
                                 27, 28};
            double m_20_buf[] = {29, 30,
                                 31, 32,
                                 33, 34};
            double m_21_buf[] = {35, 36,
                                 37, 38,
                                 39, 40,
                                 41, 42};
            std::shared_ptr<Matrix> m00(new Matrix(2, 3, m_00_buf));
            std::shared_ptr<Matrix> m01(new Matrix(2, 4, m_01_buf));
            std::shared_ptr<Matrix> m10(new Matrix(2, 3, m_10_buf));
            std::shared_ptr<Matrix> m11(new Matrix(2, 4, m_11_buf));
            std::shared_ptr<Matrix> m20(new Matrix(2, 3, m_20_buf));
            std::shared_ptr<Matrix> m21(new Matrix(2, 4, m_21_buf));

            unsigned int rowCounts[] = {2, 2, 2};
            unsigned int colCounts[] = {3, 4};

            std::vector<std::shared_ptr<InnerMatrix> > innerMatrices;

            for(int i = 0; i < 6; ++i)
            {
                std::shared_ptr<Matrix> t00(new Matrix(*m00));
                std::shared_ptr<Matrix> t01(new Matrix(*m01));
                std::shared_ptr<Matrix> t10(new Matrix(*m10));
                std::shared_ptr<Matrix> t11(new Matrix(*m11));
                std::shared_ptr<Matrix> t20(new Matrix(*m20));
                std::shared_ptr<Matrix> t21(new Matrix(*m21));

                std::shared_ptr<ScaledMatrix> s00(new ScaledMatrix(i+1, t00));
                std::shared_ptr<ScaledMatrix> s01(new ScaledMatrix(i+1, t01));
                std::shared_ptr<ScaledMatrix> s10(new ScaledMatrix(i+1, t10));
                std::shared_ptr<ScaledMatrix> s11(new ScaledMatrix(i+1, t11));
                std::shared_ptr<ScaledMatrix> s20(new ScaledMatrix(i+1, t20));
                std::shared_ptr<ScaledMatrix> s21(new ScaledMatrix(i+1, t21));

                std::shared_ptr<InnerMatrix> matrix(new InnerMatrix(3, 2, rowCounts, colCounts));
                matrix->SetBlock(0, 0, s00);
                matrix->SetBlock(0, 1, s01);
                matrix->SetBlock(1, 0, s10);
                matrix->SetBlock(1, 1, s11);
                matrix->SetBlock(2, 0, s20);
                matrix->SetBlock(2, 1, s21);

                innerMatrices.push_back(matrix);
            }

            unsigned int topLevelRowCounts[] = {6, 6};
            unsigned int topLevelColumnCounts[] = {7, 7, 7};
            BlockMatrix topLevelMatrix(2, 3, topLevelRowCounts, topLevelColumnCounts);
            topLevelMatrix.SetBlock(0, 0, innerMatrices[0]);
            topLevelMatrix.SetBlock(0, 1, innerMatrices[1]);
            topLevelMatrix.SetBlock(0, 2, innerMatrices[2]);
            topLevelMatrix.SetBlock(1, 0, innerMatrices[3]);
            topLevelMatrix.SetBlock(1, 1, innerMatrices[4]);
            topLevelMatrix.SetBlock(1, 2, innerMatrices[5]);

            double v_buf[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
            NekVector<double> v(21, v_buf);

            NekVector<double> result = topLevelMatrix*v;

            double expected_result_buf[] = {4256, 4816, 12096, 12656, 19936, 20496, 9611, 10864, 27153, 28406, 44695, 45948};
            NekVector<double> expected_result(12, expected_result_buf);

            BOOST_CHECK_EQUAL(expected_result, result);


            topLevelMatrix.Transpose();

            double transposed_v_buf[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
            NekVector<double> transposed_v(12, transposed_v_buf);
            NekVector<double> transposed_result = topLevelMatrix*transposed_v;

            double expected_transpose_result_buf[] = {4427, 4925, 5423, 5921, 6419, 6917, 7415, 5863, 6517,
                7171, 7825, 8479, 9133, 9787, 7299, 8109, 8919, 9729, 10539, 11349, 12159};
            NekVector<double> expected_transpose_result(21, expected_transpose_result_buf);

            BOOST_CHECK_EQUAL(expected_transpose_result, transposed_result);
        }

        BOOST_AUTO_TEST_CASE(TestMultiplicationBlock_1)
        {
            typedef NekMatrix<double> InnerType;
            typedef NekMatrix<InnerType, BlockMatrixTag> BlockType;

            double m_00_buf[] = {3, 1,
                                 5, 5};
            double m_01_buf[] = {2, 4,
                                 4, 1,
                                 1, 3};
            double m_02_buf[] = {3, 1,
                                 2, 2,
                                 4, 1,
                                 4, 2};

            double m_10_buf[] = {1, 3, 5,
                                 1, 4, 2};
            double m_11_buf[] = {4, 4, 1,
                                 1, 1, 4,
                                 5, 5, 3};
            double m_12_buf[] = {5, 2, 1,
                                 2, 3, 1,
                                 1, 3, 1,
                                 4, 1, 1};

            double m_20_buf[] = {3, 1, 4, 2,
                                 5, 1, 4, 2 };
            double m_21_buf[] = {4, 5, 2, 4,
                                 4, 2, 3, 5,
                                 2, 2, 5, 3};
            double m_22_buf[] = {3, 4, 3, 1,
                                 2, 1, 2, 4,
                                 4, 5, 2, 3,
                                 5, 4, 1, 1};

            double m_30_buf[] = {2, 2, 1,
                                 4, 5, 5};
            double m_31_buf[] = {1, 1, 3,
                                 2, 1, 2,
                                 3, 3, 2};
            double m_32_buf[] = {3, 1, 3,
                                 1, 2, 2,
                                 4, 2, 5,
                                 5, 1, 1};


            std::shared_ptr<InnerType> m00(new InnerType(2, 2, m_00_buf));
            std::shared_ptr<InnerType> m01(new InnerType(2, 3, m_01_buf));
            std::shared_ptr<InnerType> m02(new InnerType(2, 4, m_02_buf));

            std::shared_ptr<InnerType> m10(new InnerType(3, 2, m_10_buf));
            std::shared_ptr<InnerType> m11(new InnerType(3, 3, m_11_buf));
            std::shared_ptr<InnerType> m12(new InnerType(3, 4, m_12_buf));

            std::shared_ptr<InnerType> m20(new InnerType(4, 2, m_20_buf));
            std::shared_ptr<InnerType> m21(new InnerType(4, 3, m_21_buf));
            std::shared_ptr<InnerType> m22(new InnerType(4, 4, m_22_buf));

            std::shared_ptr<InnerType> m30(new InnerType(3, 2, m_30_buf));
            std::shared_ptr<InnerType> m31(new InnerType(3, 3, m_31_buf));
            std::shared_ptr<InnerType> m32(new InnerType(3, 4, m_32_buf));

            unsigned int rowCounts[] = {2, 3, 4, 3};
            unsigned int colCounts[] = {2, 3, 4};

            BlockType b(4, 3, rowCounts, colCounts);
            b.SetBlock(0, 0, m00);
            b.SetBlock(0, 1, m01);
            b.SetBlock(0, 2, m02);
            b.SetBlock(1, 0, m10);
            b.SetBlock(1, 1, m11);
            b.SetBlock(1, 2, m12);
            b.SetBlock(2, 0, m20);
            b.SetBlock(2, 1, m21);
            b.SetBlock(2, 2, m22);
            b.SetBlock(3, 0, m30);
            b.SetBlock(3, 1, m31);
            b.SetBlock(3, 2, m32);

            Array<OneD, unsigned int> arrayRowCounts(4, rowCounts);
            Array<OneD, unsigned int> arrayColCounts(3, colCounts);
            BlockType b_initializedWithArray(arrayRowCounts, arrayColCounts);
            b_initializedWithArray.SetBlock(0, 0, m00);
            b_initializedWithArray.SetBlock(0, 1, m01);
            b_initializedWithArray.SetBlock(0, 2, m02);
            b_initializedWithArray.SetBlock(1, 0, m10);
            b_initializedWithArray.SetBlock(1, 1, m11);
            b_initializedWithArray.SetBlock(1, 2, m12);
            b_initializedWithArray.SetBlock(2, 0, m20);
            b_initializedWithArray.SetBlock(2, 1, m21);
            b_initializedWithArray.SetBlock(2, 2, m22);
            b_initializedWithArray.SetBlock(3, 0, m30);
            b_initializedWithArray.SetBlock(3, 1, m31);
            b_initializedWithArray.SetBlock(3, 2, m32);

            double rhs_buf[] = {4, 2, 5, 5, 3, 1, 3, 2, 3};
            NekVector<double> rhs(9, rhs_buf);

            NekVector<double> result = b*rhs;
            NekVector<double> result1 = b_initializedWithArray*rhs;

            double expected_buf[] = {84, 63, 71, 80, 67, 100, 76, 80, 88, 69, 51, 67};
            NekVector<double> expected_result(12, expected_buf);

            BOOST_CHECK_EQUAL(expected_result, result);
            BOOST_CHECK_EQUAL(expected_result, result1);
        }

        BOOST_AUTO_TEST_CASE(TestDiagonalBlockMatrixMultiplication)
        {
            typedef NekMatrix<double> InnerType;
            typedef NekMatrix<InnerType, BlockMatrixTag> BlockType;

            double m_00_buf[] = {1, 3, 2, 4};
            double m_11_buf[] = {1, 4, 2, 5, 3, 6};
            double m_22_buf[] = {1, 3, 5, 2, 4, 6};


            std::shared_ptr<InnerType> m00(new InnerType(2, 2, m_00_buf));
            std::shared_ptr<InnerType> m11(new InnerType(2, 3, m_11_buf));
            std::shared_ptr<InnerType> m22(new InnerType(3, 2, m_22_buf));

            unsigned int rowCounts[] = {2, 2, 3};
            unsigned int colCounts[] = {2, 3, 2};

            BlockType b(3, 3, rowCounts, colCounts, eDIAGONAL);
            b.SetBlock(0, 0, m00);
            b.SetBlock(1, 1, m11);
            b.SetBlock(2, 2, m22);

            double rhs_buf[] = {1, 2, 3, 4, 5, 6, 7};
            NekVector<double> rhs(7, rhs_buf);

            NekVector<double> result = b*rhs;

            double expected_buf[] = {5, 11, 26, 62, 20, 46, 72};
            NekVector<double> expected_result(7, expected_buf);

            BOOST_CHECK_EQUAL(expected_result, result);
        }

        BOOST_AUTO_TEST_CASE(TestDiagonalBlockMatrixMultiplicationWith0SizeBlocks)
        {
            typedef NekMatrix<double> InnerType;
            typedef NekMatrix<InnerType, BlockMatrixTag> BlockType;

            double m_22_buf[] = {1, 2, 3, 4, 5, 6, 7, 8};

            std::shared_ptr<InnerType> m22(new InnerType(1, 8, m_22_buf));

            unsigned int rowCounts[] = {0, 0, 1};
            unsigned int colCounts[] = {6, 6, 8};

            BlockType b(3, 3, rowCounts, colCounts, eDIAGONAL);
            b.SetBlock(2, 2, m22);

            double rhs_buf[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8};
            NekVector<double> rhs(20, rhs_buf);

            NekVector<double> result = b*rhs;

            double expected_buf[] = {1*1 + 2*2 + 3*3 + 4*4 + 5*5 + 6*6 + 7*7 + 8*8};
            NekVector<double> expected_result(1, expected_buf);

            BOOST_CHECK_EQUAL(expected_result, result);

        }

        BOOST_AUTO_TEST_CASE(TestBlockMatrixMultiplicationWith0SizeBlocks)
        {
            typedef NekMatrix<double> InnerType;
            typedef NekMatrix<InnerType, BlockMatrixTag> BlockType;

            double m_22_buf[] = {1, 2, 3, 4, 5, 6, 7, 8};

            std::shared_ptr<InnerType> m22(new InnerType(1, 8, m_22_buf));

            unsigned int rowCounts[] = {0, 0, 1};
            unsigned int colCounts[] = {6, 6, 8};

            BlockType b(3, 3, rowCounts, colCounts, eFULL);
            b.SetBlock(2, 2, m22);

            double rhs_buf[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8};
            NekVector<double> rhs(20, rhs_buf);

            NekVector<double> result = b*rhs;

            double expected_buf[] = {1*1 + 2*2 + 3*3 + 4*4 + 5*5 + 6*6 + 7*7 + 8*8};
            NekVector<double> expected_result(1, expected_buf);

            BOOST_CHECK_EQUAL(expected_result, result);

        }


        BOOST_AUTO_TEST_CASE(TestBlockMatrixErrorFrom6_10)
        {

            //Array<OneD, NekDouble> f_bnd(m_BCinv->GetRows());
            //NekVector< NekDouble > F_bnd(f_bnd.size(), f_bnd, eWrapper);

            //Array<OneD, NekDouble> f_int(m_BCinv->GetColumns());
            //NekVector< NekDouble > F_int(f_int.size(),f_int, eWrapper);

            //Array<OneD, NekDouble > f_p(m_D_int->GetRows());
            //NekVector<  NekDouble > F_p(f_p.size(),f_p,eWrapper);

            //typedef NekMatrix<DNekScalBlkMat, BlockMatrixTag> BlkMatDNekScalBlkMat;
            //typedef std::shared_ptr<BlkMatDNekScalBlkMat>  BlkMatDNekScalBlkMatSharedPtr;

            //BlkMatDNekScalBlkMatSharedPtr      m_Btilde = MemoryManager<BlkMatDNekScalBlkMat>
            //    ::AllocateSharedPtr(nsize_int,nsize_bndry,blkmatStorage);
            //BlkMatDNekScalBlkMatSharedPtr   m_Cinv = MemoryManager<BlkMatDNekScalBlkMat>
            //    ::AllocateSharedPtr(nsize_int,nsize_int,eFull);

            typedef NekMatrix<double> M1;
            typedef NekMatrix<M1, ScaledMatrixTag> M2;
            typedef NekMatrix<M2, BlockMatrixTag> M3;
            typedef NekMatrix<M3, BlockMatrixTag> M4;

            std::vector< std::shared_ptr<M2> > scaledMatrices;
            for(unsigned int i = 0; i < 36; ++i)
            {
                double values[] = {i*4.0, i*4.0+1.0, i*4.0+2.0, i*4.0+3.0};
                std::shared_ptr<M1> matrix(new M1(2,2, values));
                std::shared_ptr<M2> scaled(new M2(1.0, matrix));
                scaledMatrices.push_back(scaled);
            }

            std::vector<std::shared_ptr<M3> > blockMatrices;
            for(unsigned int i = 0; i < 6; ++i)
            {
                std::shared_ptr<M3> blockMatrix(new M3(2,3,2,2));
                blockMatrix->SetBlock(0, 0, scaledMatrices[i*6]);
                blockMatrix->SetBlock(0, 1, scaledMatrices[i*6+1]);
                blockMatrix->SetBlock(0, 2, scaledMatrices[i*6+2]);
                blockMatrix->SetBlock(1, 0, scaledMatrices[i*6+3]);
                blockMatrix->SetBlock(1, 1, scaledMatrices[i*6+4]);
                blockMatrix->SetBlock(1, 2, scaledMatrices[i*6+5]);
                blockMatrices.push_back(blockMatrix);
            }

            std::shared_ptr<M4> matrix(new M4(3, 2, 4, 6));
            matrix->SetBlock(0,0, blockMatrices[0]);
            matrix->SetBlock(0,1, blockMatrices[1]);
            matrix->SetBlock(1,0, blockMatrices[2]);
            matrix->SetBlock(1,1, blockMatrices[3]);
            matrix->SetBlock(2,0, blockMatrices[4]);
            matrix->SetBlock(2,1, blockMatrices[5]);

            double b_buf[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
            NekVector<double> b(12, b_buf);

            {
                NekVector<double> result = (*matrix)*b;

                double expected_buf[] = {  1828,
                    1906,
                    2764,
                    2842,
                    5572,
                    5650,
                    6508,
                    6586,
                    9316,
                    9394,
                    10252,
                    10330 };
                NekVector<double> expected_result(12, expected_buf);
                BOOST_CHECK_EQUAL(result, expected_result);
            }

            {
                NekVector<double> result = b - (*matrix)*b;

                double expected_buf[] = {  -1827,
                    -1904,
                    -2761,
                    -2838,
                    -5567,
                    -5644,
                    -6501,
                    -6578,
                    -9307,
                    -9384,
                    -10241,
                    -10318 };

                NekVector<double> expected_result(12, expected_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }

            {
                NekVector<double> result = b;
                result = result - (*matrix)*b;

                double expected_buf[] = {  -1827,
                    -1904,
                    -2761,
                    -2838,
                    -5567,
                    -5644,
                    -6501,
                    -6578,
                    -9307,
                    -9384,
                    -10241,
                    -10318 };

                NekVector<double> expected_result(12, expected_buf);
                BOOST_CHECK_EQUAL(expected_result, result);
            }

//#if 0
//            F_int = (F_int - (*m_Btilde)*F_bnd);
//#else
//            for(i = 0; i < m_Btilde->GetRows(); ++i)
//            {
//                for(j = 0; j < m_Btilde->GetColumns(); ++j)
//                {
//                    F_int[i] -= (*m_Btilde)(i,j)*F_bnd[j];
//                }
//            }
//#endif
//
//            F_int = (F_int + Transpose(*m_D_int)*F_p);
//
//#if 0
//            F_int = (*m_Cinv)*F_int;
//#else
//            Array<OneD, NekDouble> ftmp(F_int.GetDimension(),0.0);
//            for(i = 0; i < m_Cinv->GetRows(); ++i)
//            {
//                for(j = 0; j < m_Cinv->GetColumns(); ++j)
//                {
//                    ftmp[i] += (*m_Cinv)(i,j)*F_int[j];
//                }
//            }
//            for(i = 0; i < ftmp.size(); ++i)
//            {
//                F_int[i] = ftmp[i];
//            }
//#endif

        }

    }
}
