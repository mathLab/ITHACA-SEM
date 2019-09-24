#include <LibUtilities/Interpreter/Interpreter.h>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

#include <iostream>

namespace Nektar
{
namespace InterpreterUnitTests
{

BOOST_AUTO_TEST_CASE(TestConstant)
{
    LibUtilities::Interpreter interp;
    int func1 = interp.DefineFunction("x", "-2");
    Array<OneD, NekDouble> in(1, 2.0), out(1);

    interp.Evaluate(func1, { in }, out);
    BOOST_CHECK_EQUAL(out[0], -2);
}

BOOST_AUTO_TEST_CASE(TestPowOperator)
{
    LibUtilities::Interpreter interp;
    int func1 = interp.DefineFunction("x", "5*(-(2^x)^4)");
    int func2 = interp.DefineFunction("x", "-x^2");
    int func3 = interp.DefineFunction("x", "2^2^4");
    Array<OneD, NekDouble> in(1, 2.0), out(1);

    interp.Evaluate(func1, { in }, out);
    BOOST_CHECK_EQUAL(out[0], -1280);

    interp.Evaluate(func2, { in }, out);
    BOOST_CHECK_EQUAL(out[0], -4);

    interp.Evaluate(func3, { in }, out);
    BOOST_CHECK_EQUAL(out[0], 65536);
}

}
}
