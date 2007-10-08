///////////////////////////////////////////////////////////////////////////////
//
// File: TestStandardMatrix.cpp
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

#define BOOST_TEST_MODULE LinearAlgebraUnitTests test
#include <boost/test/unit_test.hpp>

#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <UnitTests/LibUtilities/LinearAlgebra/TestCombinationRunner.h>

namespace Nektar
{
    namespace StandardMatrixUnitTests
    {

    }

    namespace StandardMatrixOperationsUnitTests
    {
        BOOST_AUTO_TEST_CASE(TestAddition)
        {
            double lhs_buf[] = {  2.0, 4.0, 6.0, 8.0,
                                  10.0, 12.0, 14.0, 16.0,
                                  18.0, 20.0, 22.0, 24.0,
                                  26.0, 28.0, 30.0, 32.0};
            double rhs_buf[] = { 3.0, 6.0, 9.0, 12.0,
                                  15.0, 18.0, 21.0, 24.0,
                                  27.0, 30.0, 33.0, 36.0,
                                  39.0, 42.0, 45.0, 48.0};

            NekMatrix<double> lhs1(4, 4, lhs_buf);
            boost::shared_ptr<NekMatrix<NekMatrix<double>, FullMatrixTag, ScaledMatrixTag> > lhs2;
            boost::shared_ptr<NekMatrix<NekMatrix<double>, FullMatrixTag, BlockMatrixTag> > lhs3;

            NekMatrix<double> rhs1(4, 4, rhs_buf);
            boost::shared_ptr<NekMatrix<NekMatrix<double>, FullMatrixTag, ScaledMatrixTag> > rhs2;
            boost::shared_ptr<NekMatrix<NekMatrix<double>, FullMatrixTag, BlockMatrixTag> > rhs3;

            GenerateMatrices(lhs1, 2.0, 2, 2, lhs2, lhs3);
            GenerateMatrices(rhs1, 3.0, 2, 2, rhs2, rhs3);

            double result_buf[] = {5, 10, 15, 20,
                                   25, 30, 35, 40,
                                   45, 50, 55, 60,
                                   65, 70, 75, 80};
            NekMatrix<double> result(4, 4, result_buf);
            
            RunAllTestCombinations(lhs1, *lhs2, *lhs3, rhs1, *rhs2, *rhs3, result, DoAddition());
        }

        BOOST_AUTO_TEST_CASE(TestSubtraction)
        {
            double lhs_buf[] = {  2.0, 4.0, 6.0, 8.0,
                                  10.0, 12.0, 14.0, 16.0,
                                  18.0, 20.0, 22.0, 24.0,
                                  26.0, 28.0, 30.0, 32.0};
            double rhs_buf[] = { 3.0, 6.0, 9.0, 12.0,
                                  15.0, 18.0, 21.0, 24.0,
                                  27.0, 30.0, 33.0, 36.0,
                                  39.0, 42.0, 45.0, 48.0};

            NekMatrix<double> lhs1(4, 4, lhs_buf);
            boost::shared_ptr<NekMatrix<NekMatrix<double>, FullMatrixTag, ScaledMatrixTag> > lhs2;
            boost::shared_ptr<NekMatrix<NekMatrix<double>, FullMatrixTag, BlockMatrixTag> > lhs3;

            NekMatrix<double> rhs1(4, 4, rhs_buf);
            boost::shared_ptr<NekMatrix<NekMatrix<double>, FullMatrixTag, ScaledMatrixTag> > rhs2;
            boost::shared_ptr<NekMatrix<NekMatrix<double>, FullMatrixTag, BlockMatrixTag> > rhs3;

            GenerateMatrices(lhs1, 2.0, 2, 2, lhs2, lhs3);
            GenerateMatrices(rhs1, 3.0, 2, 2, rhs2, rhs3);

            double result_buf[] = {-1, -2, -3, -4,
                                   -5, -6, -7, -8,
                                   -9, -10, -11, -12,
                                   -13, -14, -15, -16};
            NekMatrix<double> result(4, 4, result_buf);
            
            RunAllTestCombinations(lhs1, *lhs2, *lhs3, rhs1, *rhs2, *rhs3, result, DoSubtraction());
        }

        BOOST_AUTO_TEST_CASE(TestThreeAdditions)
        {
            double lhs_buf[] = {  2.0, 4.0, 6.0, 8.0,
                                  10.0, 12.0, 14.0, 16.0,
                                  18.0, 20.0, 22.0, 24.0,
                                  26.0, 28.0, 30.0, 32.0};
            double middle_buf[] = {  2.0, 4.0, 6.0, 8.0,
                                    10.0, 12.0, 14.0, 16.0,
                                    18.0, 20.0, 22.0, 24.0,
                                    26.0, 28.0, 30.0, 32.0};
            double rhs_buf[] = { 3.0, 6.0, 9.0, 12.0,
                                  15.0, 18.0, 21.0, 24.0,
                                  27.0, 30.0, 33.0, 36.0,
                                  39.0, 42.0, 45.0, 48.0};

            NekMatrix<double> lhs(4, 4, lhs_buf);
            NekMatrix<double> middle(4, 4, middle_buf);
            NekMatrix<double> rhs(4, 4, rhs_buf);

            NekMatrix<double> result = lhs+middle+rhs;
        }
    }

}


