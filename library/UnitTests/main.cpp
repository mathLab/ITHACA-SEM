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
#include <UnitTests/testLinearSystem.h>
#include <UnitTests/testNekLinAlgAlgorithms.h>
#include <UnitTests/testNekSharedArray.h>

#include <UnitTests/testNekManager.h>

#include <UnitTests/Memory/TestNekMemoryManager.h>
#include <UnitTests/StdRegions/testStdSegExp.h>
#include <UnitTests/testFoundation/testFoundation.h>
#include <UnitTests/testFoundation/testInterpolation.h>
#include <UnitTests/testFoundation/testDerivation.h>

// The boost unit test framework provides the main function for us.
// All we need to do is provide a test suite.
test_suite* init_unit_test_suite( int, char* [] )
{
    test_suite* test= BOOST_TEST_SUITE( "Nektar++ Test Suite" );

    
    
    // shared array
    test->add(BOOST_TEST_CASE(&Nektar::SharedArrayUnitTests::TestEmptyConstructor), 0);
    test->add(BOOST_TEST_CASE(&Nektar::SharedArrayUnitTests::TestUninitializedConstructor), 0);
    test->add(BOOST_TEST_CASE(&Nektar::SharedArrayUnitTests::TestParameterPopulation), 0);
    test->add(BOOST_TEST_CASE(&Nektar::SharedArrayUnitTests::TestSingleValueInitialization), 0);
    test->add(BOOST_TEST_CASE(&Nektar::SharedArrayUnitTests::TestPopulationFromCArray), 0);
    test->add(BOOST_TEST_CASE(&Nektar::SharedArrayUnitTests::TestCopyConstruction), 0);
    test->add(BOOST_TEST_CASE(&Nektar::SharedArrayUnitTests::Test1DAssignmentOperator), 0);
    test->add(BOOST_TEST_CASE(&Nektar::SharedArrayUnitTests::Test2DAssignmentOperator), 0);
    test->add(BOOST_TEST_CASE(&Nektar::SharedArrayUnitTests::TestOffsetAssignmentOperator), 0);
    test->add(BOOST_TEST_CASE(&Nektar::SharedArrayUnitTests::Test1DAccessOperator), 0);
    test->add(BOOST_TEST_CASE(&Nektar::SharedArrayUnitTests::Test2DAccessOperator), 0);
    

       //Test Derivation
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testDrivGaussGaussLegendre),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testDrivGaussRadauMLegendre),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testDrivGaussRadauPLegendre),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testDrivGaussLobattoLegendre),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testDrivGaussGaussChebyshev),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testDrivGaussRadauMChebyshev),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testDrivGaussRadauPChebyshev),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testDrivGaussLobattoChebyshev),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testDrivGaussRadauMAlpha0Beta1),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testDrivGaussRadauMAlpha0Beta2),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testDrivPolyEvenlySpaced),0);
//        test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testDrivFourierEvenlySpaced),0);


    //Test Interpolation
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testInterGaussGaussLegendre),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testInterGaussRadauMLegendre),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testInterGaussRadauPLegendre),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testInterGaussLobattoLegendre),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testInterGaussGaussChebyshev),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testInterGaussRadauMChebyshev),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testInterGaussRadauPChebyshev),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testInterGaussLobattoChebyshev),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testInterGaussRadauMAlpha0Beta1),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testInterGaussRadauMAlpha0Beta2),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testInterPolyEvenlySpaced),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testInterFourierEvenlySpaced),0);
   
	
	// Test Foundation
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testGaussGaussLegendre),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testGaussRadauMLegendre),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testGaussRadauPLegendre),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testGaussLobattoLegendre),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testGaussGaussChebyshev),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testGaussRadauMChebyshev),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testGaussRadauPChebyshev),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testGaussLobattoChebyshev),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testGaussRadauMAlpha0Beta1),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testGaussRadauMAlpha0Beta2),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testPolyEvenlySpaced),0);
    test->add(BOOST_TEST_CASE(&Nektar::foundationUnitTests::testFourierEvenlySpaced),0);
				
    // StdSegExp algorithms
     test->add(BOOST_TEST_CASE(&Nektar::StdSegExpUnitTests::testMassMatrix), 0);
    //test->add(BOOST_TEST_CASE(&Nektar::StdSegExpUnitTests::testLapMatrix), 0);
    test->add(BOOST_TEST_CASE(&Nektar::StdSegExpUnitTests::testIntegration), 0);
    test->add(BOOST_TEST_CASE(&Nektar::StdSegExpUnitTests::testDifferentiation), 0);
    test->add(BOOST_TEST_CASE(&Nektar::StdSegExpUnitTests::testIProductWRTBase), 0);
    test->add(BOOST_TEST_CASE(&Nektar::StdSegExpUnitTests::testFwdTrans), 0);
    test->add(BOOST_TEST_CASE(&Nektar::StdSegExpUnitTests::testBwdTrans), 0);
    test->add(BOOST_TEST_CASE(&Nektar::StdSegExpUnitTests::testPhysEvaluate), 0);
    test->add(BOOST_TEST_CASE(&Nektar::StdSegExpUnitTests::testNorms), 0);

    // Memory Manager
    test->add(BOOST_TEST_CASE(&Nektar::MemManagerUnitTests::testParameterizedConstructors), 0);
    test->add(BOOST_TEST_CASE(&Nektar::MemManagerUnitTests::testSmartPointerAllocation), 0);
    test->add(BOOST_TEST_CASE(&Nektar::MemManagerUnitTests::testArrayAllocation), 0);
    test->add(BOOST_TEST_CASE(&Nektar::MemManagerUnitTests::testSharedArrayAllocation), 0);

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
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNorms), 0);


    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testMakePtr), 0);


    //////////////////////////////
    // Matrix Tests
    //////////////////////////////
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixConstruction), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixAccess), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixBasicMath), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixFullDiagonalOperations), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testUserManagedMatrixData), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testBlockMatrices), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testBlockDiagonalMatrices), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testBlockDiagonalTimesEqual), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixTemp), 0);
    test->add(BOOST_TEST_CASE(&Nektar::MatrixUnitTests::TestFullNekMatrixGetValue), 0);
    test->add(BOOST_TEST_CASE(&Nektar::MatrixUnitTests::TestDiagonalMatrixGetValue), 0);
    test->add(BOOST_TEST_CASE(&Nektar::MatrixUnitTests::TestFullNekMatrixSetValue), 0);
    test->add(BOOST_TEST_CASE(&Nektar::MatrixUnitTests::TestDiagonalNekMatrixSetValue), 0);
    test->add(BOOST_TEST_CASE(&Nektar::MatrixUnitTests::TestFullFullMatrixAddition), 0);
    test->add(BOOST_TEST_CASE(&Nektar::MatrixUnitTests::TestFullDiagonalMatrixAddition), 0);
    test->add(BOOST_TEST_CASE(&Nektar::MatrixUnitTests::TestDiagonalDiagonalMatrixAddition), 0);
    
    test->add(BOOST_TEST_CASE(&Nektar::ScaledMatrixUnitTests::TestConstruction), 0);
    test->add(BOOST_TEST_CASE(&Nektar::ScaledMatrixUnitTests::TestElementAccess), 0);
    test->add(BOOST_TEST_CASE(&Nektar::ScaledMatrixUnitTests::TestGetNumElements), 0);
    test->add(BOOST_TEST_CASE(&Nektar::ScaledMatrixUnitTests::TestGetStorageType), 0);


    // These tests were originally added because it appeared that a NekObjectFactory
    // would be needed instead of the LokiObject factory so that the factory would
    // play nice with the NekMemoryManager.  This may not be the case, but in case
    // it comes along in the near future I'll leave these statements in.
    //test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNoParameterConstruction), 0);
    //test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testSingleParameterConstruction), 0);

    // Unit tests for NekManager
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekManager), 0);

    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testConstantExpressions), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testUnaryExpressions), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixMetadata), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testBinaryExpressions), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixMultiplication), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixSomewhatComplicatedExpression), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixComplicatedExpression), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testTemporaryGenerationFromSingleLevelBinaryExpressions), 0);
    //test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testExhaustiveSingleLevelBinaryExpressions), 0);
    //test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testExhaustive2OpBinaryExpressions), 0);


    /// Linear system tests.
    test->add(BOOST_TEST_CASE(&Nektar::LinearSystemUnitTests::testDiagonalSystem), 0);
    test->add(BOOST_TEST_CASE(&Nektar::LinearSystemUnitTests::testFullSystem), 0);
    test->add(BOOST_TEST_CASE(&Nektar::LinearSystemUnitTests::testSolvingBlockDiagonalMatrices), 0);
    test->add(BOOST_TEST_CASE(&Nektar::LinearSystemUnitTests::testMixedInputParameterTypes), 0);

    /// Linear algebra algorithsm.
    test->add(BOOST_TEST_CASE(&Nektar::NekLinAlgTests::testGramSchmidtOrthogonalizationBookExample), 0);
    

    
    return test;
}

/**
    $Log: main.cpp,v $
    Revision 1.34  2007/05/27 16:40:18  bnelson
    *** empty log message ***

    Revision 1.33  2007/05/27 14:40:56  bnelson
    Added 2D Access Operator tests.

    Revision 1.32  2007/05/25 23:46:02  ehan
    checking in derivation and interpolation

    Revision 1.31  2007/05/22 02:02:53  bnelson
    Changed Array::size to Array::num_elements.

    Fixed some compiler errors in assertions.

    Revision 1.30  2007/05/18 00:45:06  bnelson
    Added copy constructor tests.

    Revision 1.29  2007/05/15 05:19:54  bnelson
    Updated to use the new Array object.

    Revision 1.28  2007/03/29 19:42:02  bnelson
    *** empty log message ***

    Revision 1.27  2007/03/26 11:16:12  pvos
    made testStdRegions back working

    Revision 1.26  2007/03/22 05:31:09  ehan
    create testFoundation

    Revision 1.25  2007/03/16 12:09:49  pvos
    switched testStdSegExp and testNekMemoryManager back on

    Revision 1.24  2007/03/14 21:24:09  sherwin
    Update for working version of MultiRegions up to ExpList1D

    Revision 1.23  2007/03/14 11:59:45  pvos
    *** empty log message ***

    Revision 1.22  2007/03/08 17:09:03  pvos
    added StdSegExp tests

    Revision 1.21  2007/02/13 02:47:20  bnelson
    *** empty log message ***

    Revision 1.20  2007/01/29 01:37:16  bnelson
    *** empty log message ***

    Revision 1.19  2006/12/17 21:47:16  bnelson
    Added NekManager tests.

    Revision 1.18  2006/11/20 03:39:41  bnelson
    Added Gram-Schmidt tests

    Revision 1.17  2006/11/18 17:18:46  bnelson
    Added L1, L2, and Infinity norm tests.

    Revision 1.16  2006/11/11 01:32:52  bnelson
    *** empty log message ***

    Revision 1.15  2006/11/08 04:18:22  bnelson
    Added more expression template tests.

    Revision 1.14  2006/11/06 17:10:04  bnelson
    *** empty log message ***

    Revision 1.13  2006/10/30 05:08:13  bnelson
    Added preliminary linear system and block matrix support.

    Revision 1.12  2006/09/30 15:38:29  bnelson
    no message

    Revision 1.11  2006/09/11 03:28:41  bnelson
    no message

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

