////////////////////////////////////////////////////////////////////////////////
//
// main.cpp
// Sophia
// Blake Nelson
//
// Driver for the unit tests.
//
////////////////////////////////////////////////////////////////////////////////



#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/test/included/unit_test_framework.hpp>

using boost::unit_test_framework::test_suite;


#include <iostream>
using namespace std;

// FoundationsTest.cpp : Defines the entry point for the console application.
//

#include <LibUtilities/BasicUtils/NekManager.hpp>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/GaussPoints.h>
#include <LibUtilities/Foundations/PolyEPoints.h>
#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Foundations/NodalTriFekete.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace Nektar;
using namespace Nektar::LibUtilities;


namespace Nektar
{
    namespace MyNewUnitTests
    {        
        void simpleTest();
		//double polyFunc(double);
	}
}

// #include <UnitTests/testNekPoint.h>
// #include <UnitTests/testNekVector.h>
// #include <UnitTests/testBoostUtil.h>
// #include <UnitTests/testNekMatrix.h>
// #include <UnitTests/testExpressionTemplates.h>
// #include <UnitTests/testLinearSystem.h>
// #include <UnitTests/testNekLinAlgAlgorithms.h>
// 
// #include <UnitTests/testNekManager.h>
// 
// #include <UnitTests/Memory/TestNekMemoryManager.h>
// The boost unit test framework provides the main function for us.
// All we need to do is provide a test suite.
test_suite* init_unit_test_suite( int, char* [] )
{
    test_suite* test= BOOST_TEST_SUITE( "Nektar++ Test Suite" );
	
	//test->add(BOOST_TEST_CASE(&Nektar::MyNewUnitTests::simpleTest), 0);
	//test->add(BOOST_TEST_CASE(&Nektar::MyNewUnitTests::polyFunc(2)), 0);

//     // Memory Manager
//     test->add(BOOST_TEST_CASE(&Nektar::MemManagerUnitTests::testParameterizedConstructors), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::MemManagerUnitTests::testSmartPointerAllocation), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::MemManagerUnitTests::testArrayAllocation), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::MemManagerUnitTests::testSharedArrayAllocation), 0);
// 
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekPointConstruction), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekPointArithmetic), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekPointDataAccess), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekPointMisc), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekPointAssignment), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekPointPointerManipulation), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekPointComparison), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekPointOperators), 0);
// 
// 
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekVectorConstruction), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekVectorOperators), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekVectorArithmetic), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNorms), 0);
// 
// 
//     //test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testMakePtr), 0);
// 
// 
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixConstruction), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixAccess), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixBasicMath), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixFullDiagonalOperations), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testUserManagedMatrixData), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testBlockMatrices), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testBlockDiagonalMatrices), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testBlockDiagonalTimesEqual), 0);
//     
// 
//     // These tests were originally added because it appeared that a NekObjectFactory
//     // would be needed instead of the LokiObject factory so that the factory would
//     // play nice with the NekMemoryManager.  This may not be the case, but in case
//     // it comes along in the near future I'll leave these statements in.
// //     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNoParameterConstruction), 0);
// //     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testSingleParameterConstruction), 0);
// 
//     // Unit tests for NekManager
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekManager), 0);
// 
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testConstantExpressions), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testUnaryExpressions), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixMetadata), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testBinaryExpressions), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixMultiplication), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixSomewhatComplicatedExpression), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testNekMatrixComplicatedExpression), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testTemporaryGenerationFromSingleLevelBinaryExpressions), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testExhaustiveSingleLevelBinaryExpressions), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::UnitTests::testExhaustive2OpBinaryExpressions), 0);
// 
// 
//     /// Linear system tests.
//     test->add(BOOST_TEST_CASE(&Nektar::LinearSystemUnitTests::testDiagonalSystem), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::LinearSystemUnitTests::testFullSystem), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::LinearSystemUnitTests::testSolvingBlockDiagonalMatrices), 0);
//     test->add(BOOST_TEST_CASE(&Nektar::LinearSystemUnitTests::testMixedInputParameterTypes), 0);
// 
//     /// Linear algebra algorithsm.
//     test->add(BOOST_TEST_CASE(&Nektar::NekLinAlgTests::testGramSchmidtOrthogonalizationBookExample), 0);
//     
    return test;
}

namespace Nektar
{
    namespace MyNewUnitTests
    {        
        void simpleTest()
        {
			int size = 10;
			
			//PointsKey gaussKey(size, eGaussLobattoLegendre);
			
			PointsKey gaussKey(size, eGaussGaussLegendre);			
	    	PointsKey polyKey(size, ePolyEvenlySpaced);
    		boost::shared_ptr<Points<double> > ptr = PointsManager()[gaussKey];
			boost::shared_ptr<NekMatrix<double> > mat = ptr->GetI(gaussKey);
			
			int m = mat->GetRows();
			int n = mat->GetColumns();
			cout << "(m, n) = (" << m << ", " << n << ")" <<  endl;
			
			
			BOOST_CHECK(mat->GetRows() == size); // Fails here
			BOOST_CHECK(mat->GetColumns() == size);
			
			
			
			Points<double> & points = *ptr;
			cout << "points.GetPointsDim()    = " << points.GetPointsDim() << endl;
			cout << "points.GetNumPoints()    = " << points.GetNumPoints() << endl;
			cout << "points.GetTotNumPoints() = " << points.GetTotNumPoints() << endl;
			double *z = points.GetZ();
			for( int i=0; i<size; ++i ) {
				cout << "z[i] = " << z[i] <<  endl;
			}
			cout << endl;
			double *w = points.GetW();
			for( int i=0; i<size; ++i ) {
				cout << "w[i] = " << w[i] <<  endl;
			}
			
			cout << "Happy Testing!" << endl;
			
        }
		
// 	}
// 	
// 	namespace MyNewUnitTests
// 	  {	
		
// 		double polyFunc(double x)
// 		{
// 		  return (3 + x * (4 - 5 * x));
// 		  {
// 		    const double *z, *w;
// 			double Area = 0.0;			
// 			//int numPoints = GetNumPoints();
// 		    PointsKey key(5, eGaussGaussLegendre);
// 			boost::shared_ptr<Points<double> > ptr = PointsManager()[key];
// 			Points<double> & points = *ptr;
// 			int numPoints = points.GetNumPoints();
// 			ptr->GetZW(z, w);
// 			for(int i =0; i<numPoints; ++i)
// 			{
// 			  Area += w[i]*polyFunc(z[i]);
// 			   {
// 			     BOOST_CHECK(Area == 20.567);
// 			   }
// 			}
// 		  }
// 		}
			
//         void testParameterizedConstructors()
//         {
//             CountedObject<int>::ClearCounters();
// 
//             CountedObject<int>* ob1 = MemoryManager::Allocate<CountedObject<int> >();
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberDefaultConstructed, 1);
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberOf1ParameterConstructions, 0);
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberOf2ParameterConstructions, 0);
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberOf3ParameterConstructions, 0);
// 
//             CountedObject<int>* ob2 = MemoryManager::Allocate<CountedObject<int> >(1);
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberDefaultConstructed, 1);
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberOf1ParameterConstructions, 1);
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberOf2ParameterConstructions, 0);
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberOf3ParameterConstructions, 0);
// 
//             CountedObject<int>* ob3 = MemoryManager::Allocate<CountedObject<int> >(1, 2);
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberDefaultConstructed, 1);
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberOf1ParameterConstructions, 1);
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberOf2ParameterConstructions, 1);
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberOf3ParameterConstructions, 0);
// 
//             CountedObject<int>* ob4 = MemoryManager::Allocate<CountedObject<int> >(1, 2, 3);
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberDefaultConstructed, 1);
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberOf1ParameterConstructions, 1);
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberOf2ParameterConstructions, 1);
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberOf3ParameterConstructions, 1);
//             BOOST_CHECK_EQUAL(ob4->value, 6);
// 
//             MemoryManager::Deallocate<CountedObject<int> >(ob1);
//             MemoryManager::Deallocate<CountedObject<int> >(ob2);
//             MemoryManager::Deallocate<CountedObject<int> >(ob3);
//             MemoryManager::Deallocate<CountedObject<int> >(ob4);
// 
//             BOOST_CHECK(ob1 == NULL);
//             BOOST_CHECK(ob2 == NULL);
//             BOOST_CHECK(ob3 == NULL);
//             BOOST_CHECK(ob4 == NULL);
// 
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberDefaultConstructed, 1);
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberOf1ParameterConstructions, 1);
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberOf2ParameterConstructions, 1);
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberOf3ParameterConstructions, 1);
//             BOOST_CHECK_EQUAL(CountedObject<int>::numberDestroyed, 4);
//         }
	}
}
// int main(int argc, char* argv[])
// {
//   
//     PointsKey gaussKey(10, eGaussLobattoLegendre);
//     PointsKey polyKey(10,ePolyEvenlySpaced);
// 
//     boost::shared_ptr<Points<double> > ptr = PointsManager()[gaussKey];
// 
//     boost::shared_ptr<NekMatrix<double> > mat = ptr->GetI(gaussKey);
//  
//     return 0;
// }

/**
    $Log: main.cpp,v $

    Revision 1.1  2007/2/26 04:05:23  ehan
    Added the NewUnitTest project.

**/

