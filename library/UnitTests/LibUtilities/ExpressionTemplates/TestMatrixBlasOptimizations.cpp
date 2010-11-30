///////////////////////////////////////////////////////////////////////////////
//
// File: TestMatrixBlasOptimizations.cpp
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

#ifndef NEKTAR_USE_EXPRESSION_TEMPLATES
#define NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //NEKTAR_USE_EXPRESSION_TEMPLATES

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <LibUtilities/BasicUtils/OperatorGenerators.hpp>
#include <UnitTests/CountedObject.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <boost/timer.hpp>

namespace Nektar
{
    BOOST_AUTO_TEST_CASE(DgemmTest)
    {
        double a_buf[] = {1, 2, 3, 4};
        double b_buf[] = {4, 5, 6, 7};
        double c_buf[] = {8, 9, 10, 11};
        
        NekMatrix<double> a(2, 2, a_buf);
        NekMatrix<double> b(2, 2, b_buf);
        NekMatrix<double> c(2, 2, c_buf);
        
        NekMatrix<double> result = a*b + c;
        
        double expected_result_buf[] = {27, 37, 37, 51};
        NekMatrix<double> expected_result(2, 2, expected_result_buf);
        BOOST_CHECK_EQUAL(expected_result, result);
    }
    


    //    BOOST_AUTO_TEST_CASE(DgemmScaledMatrixTest)
    //{
    //    double a_buf[] = {1, 2, 3, 4};
    //    double b_buf[] = {4, 5, 6, 7};
    //    double c_buf[] = {8, 9, 10, 11};
    //    
    //    boost::shared_ptr<NekMatrix<double> > ia(new NekMatrix<double>(2, 2, a_buf));
    //    boost::shared_ptr<NekMatrix<double> > ib(new NekMatrix<double>(2, 2, b_buf));
    //    boost::shared_ptr<NekMatrix<double> > ic(new NekMatrix<double>(2, 2, c_buf));
    //    
    //    NekMatrix<NekMatrix<double>, ScaledMatrixTag> a(2.0, ia);
    //    NekMatrix<NekMatrix<double>, ScaledMatrixTag> b(3.0, ib);
    //    NekMatrix<NekMatrix<double>, ScaledMatrixTag> c(-6.0, ic);
    //    
    //    NekMatrix<double> result = a*b + c;
    //    
    //    double expected_result_buf[] = {66, 114, 102, 174};
    //    NekMatrix<double> expected_result(2, 2, expected_result_buf);
    //    BOOST_CHECK_EQUAL(expected_result, result);
    //}

    //
    //BOOST_AUTO_TEST_CASE(DgemmDiagonalTest)
    //{
    //    double a_buf[] = {1, 2, 3, 4};
    //    double b_buf[] = {4, 5, 6, 7};
    //    double c_buf[] = {8, 9, 10, 11};
    //    
    //    NekMatrix<double> a(4, 4, a_buf, eDIAGONAL);
    //    NekMatrix<double> b(4, 4, b_buf, eDIAGONAL);
    //    NekMatrix<double> c(4, 4, c_buf, eDIAGONAL);
    //    
    //    NekMatrix<double> result = a*b + c;
    //    
    //    
    //    double expected_result_buf[] = {12, 0, 0, 0,
    //                                    0, 19, 0, 0,
    //                                    0, 0, 28, 0,
    //                                    0, 0, 0, 39};
    //    NekMatrix<double> expected_result(4, 4, expected_result_buf);
    //    BOOST_CHECK_EQUAL(expected_result, result);
    //}
    //
    //BOOST_AUTO_TEST_CASE(TestMatrixMultiplyEqual01)
    //{
    //    double a_buf[] = {1, 2, 3, 4};
    //    double b_buf[] = {4, 5, 6, 7};
    //    double c_buf[] = {8, 9, 10, 11};
    //    
    //    NekMatrix<double> a(2, 2, a_buf);
    //    NekMatrix<double> b(2, 2, b_buf);
    //    NekMatrix<double> c(2, 2, c_buf);
    //    
    //    NekMatrix<double> result = a*b*c;

    //    double expected_result_buf[] = {395, 584, 487, 720};
    //    NekMatrix<double> expected_result(2, 2, expected_result_buf);
    //    BOOST_CHECK_EQUAL(expected_result, result);
    //}

    //BOOST_AUTO_TEST_CASE(TestMatrixMultiplyEqual02)
    //{
    //    double a_buf[] = {1, 2, 3, 4};
    //    double b_buf[] = {4, 5, 6, 7};
    //    double c_buf[] = {8, 9, 10, 11};
    //    double d_buf[] = {12, 13, 14, 15};
    //    
    //    NekMatrix<double> a(2, 2, a_buf);
    //    NekMatrix<double> b(2, 2, b_buf);
    //    NekMatrix<double> c(2, 2, c_buf);
    //    NekMatrix<double> d(2, 2, d_buf);
    //    
    //    NekMatrix<double> result = a*b*c*d;

    //    double expected_result_buf[] = {11071, 16368, 12835, 18976};
    //    NekMatrix<double> expected_result(2, 2, expected_result_buf);
    //    BOOST_CHECK_EQUAL(expected_result, result);
    //}

}

