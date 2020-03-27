///////////////////////////////////////////////////////////////////////////////
//
// File: TestNekVector.cpp
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

#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

namespace Nektar 
{
    namespace VariableSizedNekVectorUnitTests
    {
        BOOST_AUTO_TEST_CASE(TestConstructorWithArrayAndUserSpecifiedSize)
        {        
            int nint=1, nbndry=1;
            Array<OneD, NekDouble> offset;

            NekDouble buf[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
            Array<OneD, NekDouble> out(10, buf);

            NekVector<NekDouble> Fint(nint, offset = out + nbndry);
            BOOST_CHECK_EQUAL(1u, Fint.GetDimension());
            BOOST_CHECK_EQUAL(2u, Fint[0]);

            NekVector<NekDouble> Vint(nint, offset = out + nbndry, eWrapper); 
            BOOST_CHECK_EQUAL(1u, Vint.GetDimension());
            BOOST_CHECK_EQUAL(2, Vint[0]);

            NekVector<NekDouble> testVector(4, offset = out + 3, eWrapper); 
            BOOST_CHECK_EQUAL(4u, testVector.GetDimension());
            BOOST_CHECK_EQUAL(4, testVector[0]);
            BOOST_CHECK_EQUAL(5, testVector[1]);
            BOOST_CHECK_EQUAL(6, testVector[2]);
            BOOST_CHECK_EQUAL(7, testVector[3]);

            testVector[0] = 9.9;
            BOOST_CHECK_EQUAL(9.9, testVector[0]);
            BOOST_CHECK_EQUAL(9.9, offset[0]);
            BOOST_CHECK_EQUAL(9.9, out[3]);
        }
        
        BOOST_AUTO_TEST_CASE(ItemsThatShouldNotCompile)
        {
            Array<OneD, const NekDouble> test(10, 0.0);
            //NekVector<NekDouble> wrapper(test, eWrapper);
        }
        
        BOOST_AUTO_TEST_CASE(TestVectorSubtraction)
        {
            double lhs_buf[] = { 1.0, 2.0, 3.0, 4.0 };
            double rhs_buf[] = {10.0, 11.0, 12.0, 13.0};
            
            NekVector<double> lhs(4, lhs_buf);
            NekVector<double> rhs(4, rhs_buf);
            
            NekVector<double> result = lhs - rhs;
            
            double expected_result_buf[] = {1.0-10.0, 2.0-11.0, 3.0-12.0, 4.0-13.0};
            NekVector<double> expected_result(4, expected_result_buf);
            
            BOOST_CHECK_EQUAL(expected_result, result);
        }

        BOOST_AUTO_TEST_CASE(TestVectorAddition)
        {
            double lhs_buf[] = { 1.0, 2.0, 3.0, 4.0 };
            double rhs_buf[] = {10.0, 11.0, 12.0, 13.0};

            NekVector<double> lhs(4, lhs_buf);
            NekVector<double> rhs(4, rhs_buf);

            NekVector<double> result = lhs + rhs;

            double expected_result_buf[] = {1.0+10.0, 2.0+11.0, 3.0+12.0, 4.0+13.0};
            NekVector<double> expected_result(4, expected_result_buf);

            BOOST_CHECK_EQUAL(expected_result, result);
        }
        
    }
}


