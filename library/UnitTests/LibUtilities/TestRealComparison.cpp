///////////////////////////////////////////////////////////////////////////////
//
// File: TestRealComparison.cpp
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
// Description: unit test for IsRealEqual()
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/RealComparison.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/unit_test.hpp>

namespace Nektar
{
namespace RealComparisonUnitTests
{
    BOOST_AUTO_TEST_CASE(TestRealComparisonSmall)
    {
        // Small difference, should return true
        double ad1 = 1.00001;
        double ad2 = ad1 * (1.0 + 2*std::numeric_limits<double>::epsilon());
        BOOST_CHECK(LibUtilities::IsRealEqual(ad1, ad2));

        double ad3 = 0.00001;
        double ad4 = ad3 * (1.0 + 2*std::numeric_limits<double>::epsilon());
        BOOST_CHECK(LibUtilities::IsRealEqual(ad3, ad4));
    }

    BOOST_AUTO_TEST_CASE(TestRealComparisonNotSmall)
    {
        // "Large" difference, should return false
        double ad1 = 1.00001;
        double ad2 = ad1 * (1.0 + 50*std::numeric_limits<double>::epsilon());
        BOOST_CHECK(!LibUtilities::IsRealEqual(ad1, ad2));

        double ad3 = 0.00001;
        double ad4 = ad3 * (1.0 + 50*std::numeric_limits<double>::epsilon());
        BOOST_CHECK(!LibUtilities::IsRealEqual(ad3, ad4));
    }

    BOOST_AUTO_TEST_CASE(TestRealComparisonSmallCustomTolerance)
    {
        // Small difference, custom tolerance should return true
        double ad1 = 1.00001;
        double ad2 = ad1 * (1.0 + std::numeric_limits<double>::epsilon());
        BOOST_CHECK(LibUtilities::IsRealEqual(ad1, ad2, 1));

        double ad3 = 0.00001;
        double ad4 = ad3 * (1.0 + std::numeric_limits<double>::epsilon());
        BOOST_CHECK(LibUtilities::IsRealEqual(ad3, ad4, 1));
    }

    BOOST_AUTO_TEST_CASE(TestRealComparisonNotSmallCustomTolerance)
    {
        // "Large" difference, custom tolerance should return true
        double ad1 = 1.00001;
        double ad2 = ad1 * (1.0 + 50*std::numeric_limits<double>::epsilon());
        BOOST_CHECK(LibUtilities::IsRealEqual(ad1, ad2, 100));
    }

    BOOST_AUTO_TEST_CASE(TestRealComparisonZero)
    {
        // Both zero, dist and ref value will be zero
        double ad1 = 0.0;
        double ad2 = 0.0;
        BOOST_CHECK(LibUtilities::IsRealEqual(ad1, ad2));
    }

    BOOST_AUTO_TEST_CASE(TestRealComparisonZeroNotZero)
    {
        // This should return false
        double ad1 = 0.0;
        double ad2 = 1.0;
        BOOST_CHECK(!LibUtilities::IsRealEqual(ad1, ad2));
    }

    BOOST_AUTO_TEST_CASE(TestRealComparisonRealLiteral)
    {
        // Large difference, should return false
        double ad1 = 1.00001;
        BOOST_CHECK(!LibUtilities::IsRealEqual(ad1, 1.000011));
    }

    BOOST_AUTO_TEST_CASE(TestRealComparisonLiteral)
    {
        // Large difference, should return false
        BOOST_CHECK(!LibUtilities::IsRealEqual(1.00001, 1.000011));
    }

    BOOST_AUTO_TEST_CASE(TestRealComparisonConstNotConst)
    {
        // This should return false
        const double ad1 = 0.0;
        double ad2 = 1.0;
        BOOST_CHECK(!LibUtilities::IsRealEqual(ad1, ad2));
    }

    BOOST_AUTO_TEST_CASE(TestRealComparisonNotConstConst)
    {
        // This should return false
        double ad1 = 0.0;
        const double ad2 = 1.0;
        BOOST_CHECK(!LibUtilities::IsRealEqual(ad1, ad2));
    }

    BOOST_AUTO_TEST_CASE(TestRealComparisonRef)
    {
        // This should return false
        double ad1 = 0.0;
        double& adr = ad1;
        BOOST_CHECK(!LibUtilities::IsRealEqual(adr, 1.0));
    }

    // // This should not compile because comparing to an int
    // BOOST_AUTO_TEST_CASE(TestRealComparisonRef)
    // {
    //     // This should return false
    //     double ad1 = 0.0;
    //     double& adr = ad1;
    //     BOOST_CHECK(!LibUtilities::IsRealEqual(adr, 1));
    // }

    // The precondition is tested only in debug mode
    #if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)

    BOOST_AUTO_TEST_CASE(TestRealComparisonThrow)
    {
        // This should throw
        double ad1 = 0.0;
        double ad2 = 0.0;
        BOOST_CHECK_THROW(LibUtilities::IsRealEqual(ad1, ad2, 0),
            std::runtime_error);
    }

    BOOST_AUTO_TEST_CASE(TestRealComparisonNoThrow)
    {
        // This should not throw
        double ad1 = 0.0;
        double ad2 = 0.0;
        BOOST_CHECK(LibUtilities::IsRealEqual(ad1, ad2, 1));
    }

    #endif

    // This should not compile because of the different float types
    // BOOST_AUTO_TEST_CASE(TestRealComparisonMixRealTypeThrow)
    // {
    //     double ad1 = 0.0;
    //     float ad2 = 0.0;
    //     BOOST_CHECK_THROW(LibUtilities::IsRealEqual(ad1, ad2, 0),
    //         std::runtime_error);
    // }
}
}
