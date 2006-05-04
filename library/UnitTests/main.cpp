////////////////////////////////////////////////////////////////////////////////
//
// main.cpp
// Blake Nelson
//
// Driver for the unit tests.
//
////////////////////////////////////////////////////////////////////////////////

#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/test/included/unit_test_framework.hpp>

using boost::unit_test_framework::test_suite;

#include <UnitTests/testNekPoint.h>
#include <UnitTests/testNekVector.h>

// The boost unit test framework provides the main function for us.
// All we need to do is provide a test suite.
test_suite* init_unit_test_suite( int, char* [] )
{
    test_suite* test= BOOST_TEST_SUITE( "Nektar++ Test Suite" );

    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekPointConstruction), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekPointArithmetic), 0);


    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekVectorConstruction), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekVectorOperators), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekVectorArithmetic), 0);

    return test;
}

/**
    $Log: main.cpp,v $
    Revision 1.2  2006/04/11 02:02:13  bnelson
    Added more tests.

    Revision 1.1  2006/01/31 14:06:23  bnelson
    Added the new UnitTest project.

**/

