///////////////////////////////////////////////////////////////////////////////
//
// File: TestDiagonalMatrixOperations.cpp
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
    namespace DiagonalMatrixOperationsUnitTests
    {
        BOOST_AUTO_TEST_CASE(TestDoubleDiagonalVectorMultiplication)
        {
            double m_buf[] = {1, 2, 3};
            double v_buf[] = {4, 5, 6};
            
            NekMatrix<double> m(3, 3, m_buf, eDIAGONAL);
            NekVector<double> v(3, v_buf);
            
            double expected_result_buf[] = {4, 10, 18};
            NekVector<double> expected_result(3, expected_result_buf);
            
            NekVector<double> result = m*v;
            
            BOOST_CHECK_EQUAL(expected_result, result);
        }
        
        BOOST_AUTO_TEST_CASE(TestDoubleScaledDiagonalVectorMultiplication)
        {
            double m_buf[] = {1, 2, 3};
            double v_buf[] = {4, 5, 6};
            
            std::shared_ptr<NekMatrix<double> > inner(
                new NekMatrix<double>(3, 3, m_buf, eDIAGONAL));
            NekMatrix<NekMatrix<double>, ScaledMatrixTag> 
                m(5.0, inner);
                
            NekVector<double> v(3, v_buf);
            
            double expected_result_buf[] = {20, 50, 90};
            NekVector<double> expected_result(3, expected_result_buf);
            
            NekVector<double> result = m*v;
            
            BOOST_CHECK_EQUAL(expected_result, result);

        }
       
    }
}

