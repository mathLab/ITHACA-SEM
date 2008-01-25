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
            NekMatrix<NekMatrix<double>, FullMatrixTag, BlockMatrixTag> m(2, 3, 7, 10);
            BOOST_CHECK_EQUAL(7, m.GetNumberOfRowsInBlockRow(0));
            BOOST_CHECK_EQUAL(7, m.GetNumberOfRowsInBlockRow(1));
            BOOST_CHECK_EQUAL(10, m.GetNumberOfColumnsInBlockColumn(0));
            BOOST_CHECK_EQUAL(10, m.GetNumberOfColumnsInBlockColumn(1));
            BOOST_CHECK_EQUAL(10, m.GetNumberOfColumnsInBlockColumn(2));
        }

        BOOST_AUTO_TEST_CASE(TestMultiplication)
        {
            typedef NekMatrix<double> InnerType;
            typedef NekMatrix<InnerType, FullMatrixTag, BlockMatrixTag> BlockType;

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

            boost::shared_ptr<InnerType> m1(new InnerType(2, 3, m1_buf));
            boost::shared_ptr<InnerType> m2(new InnerType(2, 2, m2_buf));
            boost::shared_ptr<InnerType> m3(new InnerType(2, 3, m3_buf));
            boost::shared_ptr<InnerType> m4(new InnerType(2, 2, m4_buf));

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
    }
}
