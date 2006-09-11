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
#include <UnitTests/testBoostUtil.h>
#include <UnitTests/testNekMatrix.h>
#include <UnitTests/testExpressionTemplates.h>
//#include <UnitTests/testNekManager.h>

// The boost unit test framework provides the main function for us.
// All we need to do is provide a test suite.
test_suite* init_unit_test_suite( int, char* [] )
{
    test_suite* test= BOOST_TEST_SUITE( "Nektar++ Test Suite" );

    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekPointConstruction), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekPointArithmetic), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekPointDataAccess), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekPointMisc), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekPointAssignment), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekPointPointerManipulation), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekPointComparison), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekPointOperators), 0);


    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekVectorConstruction), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekVectorOperators), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekVectorArithmetic), 0);


    //test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testMakePtr), 0);


    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixConstruction), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixAccess), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixBasicMath), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixFullDiagonalOperations), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testUserManagedMatrixData), 0);


    // These tests were originally added because it appeared that a NekObjectFactory
    // would be needed instead of the LokiObject factory so that the factory would
    // play nice with the NekMemoryManager.  This may not be the case, but in case
    // it comes along in the near future I'll leave these statements in.
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNoParameterConstruction), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testSingleParameterConstruction), 0);

    // Unit tests for NekManager
    //test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekManager), 0);





    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testConstantExpressions), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testUnaryExpressions), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixMetadata), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testBinaryExpressions), 0);

    return test;
}

/**
    $Log: main.cpp,v $
    Revision 1.10  2006/08/28 02:40:50  bnelson
    *** empty log message ***

    Revision 1.9  2006/08/25 01:38:58  bnelson
    no message

    Revision 1.8  2006/08/25 01:36:25  bnelson
    no message

    Revision 1.7  2006/08/14 02:35:45  bnelson
    Added many LinearAlgebra tests

    Revision 1.6  2006/07/05 20:21:24  jfrazier
    Added NekManager test case.

    Revision 1.5  2006/06/05 02:23:17  bnelson
    Updates for the reorganization of LibUtilities.

    Revision 1.4  2006/05/25 02:54:53  bnelson
    Added Matrix/Vector multiplication test.

    Revision 1.3  2006/05/07 21:11:13  bnelson
    Added NekMatrix tests.

    Revision 1.2  2006/05/07 15:23:28  bnelson
    Added tests for boost utilities.

    Revision 1.1  2006/05/04 18:59:55  kirby
    *** empty log message ***

    Revision 1.2  2006/04/11 02:02:13  bnelson
    Added more tests.

    Revision 1.1  2006/01/31 14:06:23  bnelson
    Added the new UnitTest project.

**/

