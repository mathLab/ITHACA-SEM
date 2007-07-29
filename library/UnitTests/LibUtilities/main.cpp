///////////////////////////////////////////////////////////////////////////////
//
// File: main.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: Unit tests for LibUtilities
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/test/included/unit_test_framework.hpp>
#include <UnitTests/LibUtilities/TestMatrixStoragePolicies.h>
#include <UnitTests/LibUtilities/TestUpperTriangularMatrix.h>
#include <UnitTests/LibUtilities/TestNekMatrixOperations.h>

using boost::unit_test_framework::test_suite;

// The boost unit test framework provides the main function for us.
// All we need to do is provide a test suite.

// On Windows, to turn off memory leak detection, --detect_memory_leaks=0
test_suite* init_unit_test_suite( int, char* [] )
{   
    test_suite* test= BOOST_TEST_SUITE( "Lib Utilities Test Suite" );

    // Nektar::ptr tests.
    test->add(BOOST_TEST_CASE(&Nektar::UpperTriangularUnitTests::Test0ParameterInitialize), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UpperTriangularUnitTests::Test2ParameterInitialize), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UpperTriangularUnitTests::TestSingleValuePopulationInitialize), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UpperTriangularUnitTests::TestCArrayInitialization), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UpperTriangularUnitTests::TestArrayInitialization), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UpperTriangularUnitTests::TestConstGetValue), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UpperTriangularUnitTests::TestSetValue), 0);
    test->add(BOOST_TEST_CASE(&Nektar::UpperTriangularUnitTests::TestAdvance), 0);
    
    // Matrix operations
    test->add(BOOST_TEST_CASE(&Nektar::MatrixOperationTests::TestLhsFullRhsFull), 0);
    test->add(BOOST_TEST_CASE(&Nektar::MatrixOperationTests::TestLhsFullRhsDiagonal), 0);
    test->add(BOOST_TEST_CASE(&Nektar::MatrixOperationTests::TestComboExpression), 0);
    test->add(BOOST_TEST_CASE(&Nektar::MatrixOperationTests::TestLhsUpperTriangularRhsUpperTriangular), 0);
    test->add(BOOST_TEST_CASE(&Nektar::MatrixOperationTests::TestLhsLowerTriangularRhsLowerTriangular), 0);

    test->add(BOOST_TEST_CASE(&Nektar::MatrixSubtractionTests::TestLhsFullRhsFull), 0);

    return test;
}

