///////////////////////////////////////////////////////////////////////////////
//
// File: TestScaledBlockMatrixOperations.cpp
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
    namespace ScaledBlockMatrixUnitTests
    {
        BOOST_AUTO_TEST_CASE(TestMatrixVectorMultiplication)
        {
            typedef NekMatrix<double> InnerType;
            typedef NekMatrix<InnerType, ScaledMatrixTag> ScaleType;

            double m1_buf[] = {1.0/2.0, 6.0/2.0,
                               2.0/2.0, 7.0/2.0,
                               3.0/2.0, 8.0/2.0};
            double m2_buf[] = {4.0/2.0, 9.0/2.0,
                               5.0/2.0, 10.0/2.0};
            double m3_buf[] = {11.0/2.0, 16.0/2.0,
                               12.0/2.0, 17.0/2.0,
                               13.0/2.0, 18/2.0};
            double m4_buf[] = {14.0/2.0, 19.0/2.0,
                               15.0/2.0, 20.0/2.0};

            std::shared_ptr<InnerType> in_m1(new InnerType(2, 3, m1_buf));
            std::shared_ptr<InnerType> in_m2(new InnerType(2, 2, m2_buf));
            std::shared_ptr<InnerType> in_m3(new InnerType(2, 3, m3_buf));
            std::shared_ptr<InnerType> in_m4(new InnerType(2, 2, m4_buf));

            std::shared_ptr<ScaleType> m1(new ScaleType(2.0, in_m1));
            std::shared_ptr<ScaleType> m2(new ScaleType(2.0, in_m2));
            std::shared_ptr<ScaleType> m3(new ScaleType(2.0, in_m3));
            std::shared_ptr<ScaleType> m4(new ScaleType(2.0, in_m4));
            
            unsigned int rowCounts[] = {2, 2};
            unsigned int colCounts[] = {3, 2};
            
            typedef NekMatrix<ScaleType, BlockMatrixTag> BlockType;

            std::shared_ptr<BlockType> b(new BlockType(2, 2, rowCounts, colCounts));
            b->SetBlock(0, 0, m1);
            b->SetBlock(0, 1, m2);
            b->SetBlock(1, 0, m3);
            b->SetBlock(1, 1, m4);
            
            double rhs_buf[] = {10, 20, 30, 40, 50};
            NekVector<double> rhs(5, rhs_buf);

            NekVector<double> result = (*b)*rhs;

            double expected_buf[] = {550, 1300, 2050, 2800};
            NekVector<double> expected_result(4, expected_buf);

            BOOST_CHECK_EQUAL(expected_result, result);
        }

        BOOST_AUTO_TEST_CASE(TestMultiplicationScaledBlock_1)
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
    }
}
