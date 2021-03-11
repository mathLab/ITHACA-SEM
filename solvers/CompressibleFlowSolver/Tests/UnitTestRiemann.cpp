#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>
namespace Nektar
{
namespace RiemannTests
{

    BOOST_AUTO_TEST_CASE(RRR)
    {

        BOOST_CHECK_CLOSE(1., 1., 1e-10);
    }

}
}