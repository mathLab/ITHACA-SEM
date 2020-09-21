#include <LibUtilities/BasicUtils/CheckedCast.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

#include <iostream>

namespace Nektar
{
namespace LibUtilities
{
namespace CheckCastUnitTest
{

BOOST_AUTO_TEST_CASE(TestDoubleToInt)
{
    // expecting to convert
    {
        double adouble = std::numeric_limits<int>::max();
        int aint = checked_cast<int>(adouble);
        BOOST_CHECK_EQUAL(aint, adouble);
    }
    {
        double adouble = std::numeric_limits<int>::min();
        int aint = checked_cast<int>(adouble);
        BOOST_CHECK_EQUAL(aint, adouble);
    }

    // expecting to fail and throw
    try
    {
        double adouble = std::numeric_limits<int>::max()+1.0;
        int aint = checked_cast<int>(adouble);
        BOOST_CHECK_EQUAL(aint, adouble);
    }
    catch (std::runtime_error& e)
    {
        std::string errmss = e.what();
        BOOST_CHECK_EQUAL("Level 0 assertion violation", errmss.substr(0,27));
    }

    try
    {
        double adouble = std::numeric_limits<int>::min()-1.0;
        int aint = checked_cast<int>(adouble);
        BOOST_CHECK_EQUAL(aint, adouble);
    }
    catch (std::runtime_error& e)
    {
        std::string errmss = e.what();
        BOOST_CHECK_EQUAL("Level 0 assertion violation", errmss.substr(0,27));
    }

}



}
}
}