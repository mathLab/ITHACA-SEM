///////////////////////////////////////////////////////////////////////////////
//
// File: TestFullMatrixOperations.cpp
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
    namespace FullMatrixOperationsUnitTests
    {
        BOOST_AUTO_TEST_CASE(TestDoubleSquareFullVectorMultiplication)
        {
            double m_buf[] = {1, 2, 3,
                              4, 5, 6,
                              7, 8, 9};
            double v_buf[] = {4, 5, 6};
            
            NekMatrix<double> m(3, 3, m_buf);
            NekVector<double> v(3, v_buf);
            
            double expected_result_buf[] = {66, 81, 96};
            NekVector<double> expected_result(3, expected_result_buf);
            
            NekVector<double> result = m*v;
            
            BOOST_CHECK_EQUAL(expected_result, result);
        }

        BOOST_AUTO_TEST_CASE(TestIntSquareFullVectorMultiplication)
        {
            int m_buf[] = {1, 2, 3,
                              4, 5, 6,
                              7, 8, 9};
            int v_buf[] = {4, 5, 6};
            
            NekMatrix<int> m(3, 3, m_buf);
            NekVector<int> v(3, v_buf);
            
            int expected_result_buf[] = {66, 81, 96};
            NekVector<int> expected_result(3, expected_result_buf);
            
            NekVector<int> result = m*v;
            
            BOOST_CHECK_EQUAL(expected_result, result);
        }
        
        BOOST_AUTO_TEST_CASE(TestScaledDoubleSquareFullVectorMultiplication)
        {
            double m_buf[] = {1, 2, 3,
                              4, 5, 6,
                              7, 8, 9};
            double v_buf[] = {4, 5, 6};
            
            boost::shared_ptr<NekMatrix<double> > inner(new NekMatrix<double>(3, 3, m_buf));
            NekMatrix<NekMatrix<double>, FullMatrixTag, ScaledMatrixTag> m(7.0, inner);
            NekVector<double> v(3, v_buf);
            
            double expected_result_buf[] = {462, 567, 672};
            NekVector<double> expected_result(3, expected_result_buf);
            
            NekVector<double> result = m*v;
            
            BOOST_CHECK_EQUAL(expected_result, result);
        }

        BOOST_AUTO_TEST_CASE(TestScaledIntSquareFullVectorMultiplication)
        {
            int m_buf[] = {1, 2, 3,
                              4, 5, 6,
                              7, 8, 9};
            int v_buf[] = {4, 5, 6};
            
            boost::shared_ptr<NekMatrix<int> > inner(new NekMatrix<int>(3, 3, m_buf));
            NekMatrix<NekMatrix<int>, FullMatrixTag, ScaledMatrixTag> m(7, inner);
            NekVector<int> v(3, v_buf);
            
            int expected_result_buf[] = {462, 567, 672};
            NekVector<int> expected_result(3, expected_result_buf);
            
            NekVector<int> result = m*v;
            
            BOOST_CHECK_EQUAL(expected_result, result);
        }

    }
}
