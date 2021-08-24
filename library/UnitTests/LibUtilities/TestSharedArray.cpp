///////////////////////////////////////////////////////////////////////////////
//
// File: TestSharedArray.cpp
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

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <boost/core/ignore_unused.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

namespace Nektar
{
    namespace SharedArrayUnitTests
    {
        BOOST_AUTO_TEST_CASE(TestArrayConstructionFromConstantArray)
        {
            Array<OneD, const double> const_array_1(10, 7.0);
            Array<OneD, const double> const_array_2(10, 3.0);

            Array<OneD, double> array_1(const_array_1);
            Array<OneD, double> array_2(5, const_array_2);

            BOOST_CHECK_EQUAL(array_1.size(), const_array_1.size());
            BOOST_CHECK_EQUAL(array_2.size(), 5);

            array_1[2] = -1.0;
            array_2[2] = -1.0;

            BOOST_CHECK_EQUAL(array_1[2], -1.0);
            BOOST_CHECK_EQUAL(array_2[2], -1.0);

            BOOST_CHECK_EQUAL(const_array_1[2], 7.0);
            BOOST_CHECK_EQUAL(const_array_2[2], 3.0);
        }

        void CheckAddresses(Array<TwoD, double>::reference d, double* expectedAddress)
        {
            BOOST_CHECK_EQUAL(d.size(), 7);
            BOOST_CHECK_EQUAL(d.origin(), expectedAddress);
        }

        BOOST_AUTO_TEST_CASE(TestRowPointers)
        {
            Array<TwoD, double> array_1(10, 7, 0.0);
            CheckAddresses(array_1[0], array_1.data());
            CheckAddresses(array_1[1], array_1.data()+7);
        }
    }
}
