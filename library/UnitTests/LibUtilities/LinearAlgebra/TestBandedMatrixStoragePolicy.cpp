///////////////////////////////////////////////////////////////////////////////
//
// File: TestBandedMatrixStoragePolicy.cpp
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
#include <UnitTests/util.h>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

namespace Nektar
{
    namespace BandedMatrixStoragePolicyUnitTests
    {
        typedef BandedMatrixFuncs Policy;

        BOOST_AUTO_TEST_CASE(TestCalculateStorageSizeAndCalculateNumberOfRows)
        {
            UnitTests::RedirectCerrIfNeeded();
            BOOST_CHECK_EQUAL(16u, Policy::GetRequiredStorageSize(4, 4, 2, 1));
            BOOST_CHECK_EQUAL(4u, Policy::CalculateNumberOfRows(4, 2, 1));
            BOOST_CHECK_EQUAL(4u, Policy::CalculateNumberOfRows(3, 2, 1));
        }

        BOOST_AUTO_TEST_CASE(TestDiagonalOnlyCalculateIndex)
        {
            UnitTests::RedirectCerrIfNeeded();

            BOOST_CHECK_EQUAL(0, Policy::CalculateIndex(3, 3, 0, 0, 0, 0));
            BOOST_CHECK_EQUAL(1, Policy::CalculateIndex(3, 3, 1, 1, 0, 0));
            BOOST_CHECK_EQUAL(2, Policy::CalculateIndex(3, 3, 2, 2, 0, 0));
            BOOST_CHECK_EQUAL(0, Policy::CalculateIndex(1, 1, 0, 0, 0, 0));

            BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), Policy::CalculateIndex(3, 3, 0, 1, 0, 0)); 
            BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), Policy::CalculateIndex(3, 3, 0, 2, 0, 0)); 
            BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), Policy::CalculateIndex(3, 3, 1, 0, 0, 0)); 
            BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), Policy::CalculateIndex(3, 3, 1, 2, 0, 0)); 
            BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), Policy::CalculateIndex(3, 3, 2, 0, 0, 0)); 
            BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), Policy::CalculateIndex(3, 3, 2, 1, 0, 0)); 
        }

        BOOST_AUTO_TEST_CASE(TestSubDiagonalsOnlyCalculateIndex)
        {
            UnitTests::RedirectCerrIfNeeded();

            BOOST_CHECK_EQUAL(0, Policy::CalculateIndex(3, 3, 0, 0, 1, 0));
            BOOST_CHECK_EQUAL(2, Policy::CalculateIndex(3, 3, 1, 1, 1, 0));
            BOOST_CHECK_EQUAL(4, Policy::CalculateIndex(3, 3, 2, 2, 1, 0));
            BOOST_CHECK_EQUAL(1, Policy::CalculateIndex(3, 3, 1, 0, 1, 0));
            BOOST_CHECK_EQUAL(3, Policy::CalculateIndex(3, 3, 2, 1, 1, 0));

            BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), Policy::CalculateIndex(3, 3, 0, 1, 1, 0));
            BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), Policy::CalculateIndex(3, 3, 0, 2, 1, 0));
            BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), Policy::CalculateIndex(3, 3, 2, 0, 1, 0));
            BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), Policy::CalculateIndex(3, 3, 1, 2, 1, 0));
        }

        BOOST_AUTO_TEST_CASE(TestSuperDiagonalsOnlyCalculateIndex)
        {
            UnitTests::RedirectCerrIfNeeded();

            BOOST_CHECK_EQUAL(2, Policy::CalculateIndex(3, 3, 0, 0, 0, 2));
            BOOST_CHECK_EQUAL(5, Policy::CalculateIndex(3, 3, 1, 1, 0, 2));
            BOOST_CHECK_EQUAL(8, Policy::CalculateIndex(3, 3, 2, 2, 0, 2));
            BOOST_CHECK_EQUAL(4, Policy::CalculateIndex(3, 3, 0, 1, 0, 2));
            BOOST_CHECK_EQUAL(7, Policy::CalculateIndex(3, 3, 1, 2, 0, 2));
            BOOST_CHECK_EQUAL(6, Policy::CalculateIndex(3, 3, 0, 2, 0, 2));

            BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), Policy::CalculateIndex(3, 3, 1, 0, 0, 2));
            BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), Policy::CalculateIndex(3, 3, 2, 0, 0, 2));
            BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), Policy::CalculateIndex(3, 3, 2, 1, 0, 2));
        }

        BOOST_AUTO_TEST_CASE(TestSupAndSuperDiagonalsCalculateIndex)
        {
            UnitTests::RedirectCerrIfNeeded();


            BOOST_CHECK_EQUAL(2, Policy::CalculateIndex(4, 4, 0, 0, 1, 2));
            BOOST_CHECK_EQUAL(5, Policy::CalculateIndex(4, 4, 0, 1, 1, 2));
            BOOST_CHECK_EQUAL(8, Policy::CalculateIndex(4, 4, 0, 2, 1, 2));
            BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), Policy::CalculateIndex(4, 4, 0, 3, 1, 2));

            BOOST_CHECK_EQUAL(3, Policy::CalculateIndex(4, 4, 1, 0, 1, 2));
            BOOST_CHECK_EQUAL(6, Policy::CalculateIndex(4, 4, 1, 1, 1, 2));
            BOOST_CHECK_EQUAL(9, Policy::CalculateIndex(4, 4, 1, 2, 1, 2));
            BOOST_CHECK_EQUAL(12, Policy::CalculateIndex(4, 4, 1, 3, 1, 2));

            BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), Policy::CalculateIndex(4, 4, 2, 0, 1, 2));
            BOOST_CHECK_EQUAL(7, Policy::CalculateIndex(4, 4, 2, 1, 1, 2));
            BOOST_CHECK_EQUAL(10, Policy::CalculateIndex(4, 4, 2, 2, 1, 2));
            BOOST_CHECK_EQUAL(13, Policy::CalculateIndex(4, 4, 2, 3, 1, 2));

            BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), Policy::CalculateIndex(4, 4, 3, 0, 1, 2));
            BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), Policy::CalculateIndex(4, 4, 3, 1, 1, 2));
            BOOST_CHECK_EQUAL(11, Policy::CalculateIndex(4, 4, 3, 2, 1, 2));
            BOOST_CHECK_EQUAL(14, Policy::CalculateIndex(4, 4, 3, 3, 1, 2));
        }

        BOOST_AUTO_TEST_CASE(TestSetValue)
        {
            UnitTests::RedirectCerrIfNeeded();

            // [ 1 2 3 0 ]
            // [ 4 5 6 7 ]
            // [ 0 8 9 10 ]
            // [ 0 0 11 12 ]
            NekDouble buf[] = { 0, 0, 1, 4,
                                               0, 2, 5, 8,
                                               3, 6, 9, 11,
                                               7, 10, 12, 0 };
            NekMatrix<NekDouble> m(4, 4, buf, eBANDED, 1,2);

            BOOST_CHECK_EQUAL(1, m(0,0));
            m.SetValue(0, 0, 19);
            BOOST_CHECK_EQUAL(19, m(0,0));
            BOOST_CHECK_EQUAL(2, m(0,1));
            m.SetValue(0, 1, 20);
            BOOST_CHECK_EQUAL(20, m(0,1));
            BOOST_CHECK_EQUAL(3, m(0,2));
            m.SetValue(0, 2, 21);
            BOOST_CHECK_EQUAL(21, m(0,2));
            BOOST_CHECK_EQUAL(0, m(0,3));
            BOOST_CHECK_THROW(m.SetValue(0, 3, 21), ErrorUtil::NekError);

            BOOST_CHECK_EQUAL(4, m(1,0));
            m.SetValue(1, 0, 22);
            BOOST_CHECK_EQUAL(22, m(1,0));
            BOOST_CHECK_EQUAL(5, m(1,1));
            m.SetValue(1, 1, 23);
            BOOST_CHECK_EQUAL(23, m(1,1));
            BOOST_CHECK_EQUAL(6, m(1,2));
            m.SetValue(1, 2, 24);
            BOOST_CHECK_EQUAL(24, m(1,2));
            BOOST_CHECK_EQUAL(7, m(1,3));
            m.SetValue(1, 3, 25);
            BOOST_CHECK_EQUAL(25, m(1,3));

            BOOST_CHECK_EQUAL(0, m(2,0));
            BOOST_CHECK_THROW(m.SetValue(2, 0, 21), ErrorUtil::NekError);
            BOOST_CHECK_EQUAL(8, m(2,1));
            m.SetValue(2, 1, 26);
            BOOST_CHECK_EQUAL(26, m(2,1));
            BOOST_CHECK_EQUAL(9, m(2,2));
            m.SetValue(2, 2, 27);
            BOOST_CHECK_EQUAL(27, m(2,2));
            BOOST_CHECK_EQUAL(10, m(2,3));
            m.SetValue(2, 3, 28);
            BOOST_CHECK_EQUAL(28, m(2,3));

            BOOST_CHECK_EQUAL(0, m(3,0));
            BOOST_CHECK_THROW(m.SetValue(3, 0, 21), ErrorUtil::NekError);
            BOOST_CHECK_EQUAL(0, m(3,1));
            BOOST_CHECK_THROW(m.SetValue(3, 1, 21), ErrorUtil::NekError);
            BOOST_CHECK_EQUAL(11, m(3,2));
            m.SetValue(3, 2, 29);
            BOOST_CHECK_EQUAL(29, m(3,2));
            BOOST_CHECK_EQUAL(12, m(3,3));
            m.SetValue(3, 3, 30);
            BOOST_CHECK_EQUAL(30, m(3,3));
        }
    }
}




