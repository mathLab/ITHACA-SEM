///////////////////////////////////////////////////////////////////////////////
//
// File: TestBandedMatrixStoragePolicy.cpp
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
// Description: 
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/test/unit_test.hpp>

#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <UnitTests/CountedObject.h>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

namespace Nektar
{
    namespace BandedMatrixStoragePolicyUnitTests
    {
        typedef MatrixStoragePolicy<NekDouble, BandedMatrixTag> Policy;
        typedef Policy::PolicySpecificDataHolderType DataHolderType;

        BOOST_AUTO_TEST_CASE(TestUninitializedDataConstructionWithDefaultDataHolder)
        {
            Array<OneD, NekDouble> result = Policy::Initialize(10, 10, DataHolderType());
            BOOST_CHECK_EQUAL(result.num_elements(), 190);

            BOOST_CHECK_THROW(Policy::Initialize(10, 9, DataHolderType()), ErrorUtil::NekError);
            BOOST_CHECK_THROW(Policy::Initialize(9, 10, DataHolderType()), ErrorUtil::NekError);
        }

        BOOST_AUTO_TEST_CASE(TestUninitializedDataConstructionWithUserDefinedNumberOfDiagonals)
        {
            DataHolderType d1(1,0);
            Array<OneD, NekDouble> result1 = Policy::Initialize(10, 10, d1);
            BOOST_CHECK_EQUAL(result1.num_elements(), 20);

            DataHolderType d2(0,1);
            Array<OneD, NekDouble> result2 = Policy::Initialize(10, 10, d2);
            BOOST_CHECK_EQUAL(result2.num_elements(), 20);

            DataHolderType d3(1,1);
            Array<OneD, NekDouble> result3 = Policy::Initialize(10, 10, d3);
            BOOST_CHECK_EQUAL(result3.num_elements(), 30);

            DataHolderType d4(3,1);
            Array<OneD, NekDouble> result4 = Policy::Initialize(10, 10, d4);
            BOOST_CHECK_EQUAL(result4.num_elements(), 50);

            DataHolderType d5(5,1);
            Array<OneD, NekDouble> result5 = Policy::Initialize(10, 10, d5);
            BOOST_CHECK_EQUAL(result5.num_elements(), 70);

            DataHolderType d6(5,6);
            Array<OneD, NekDouble> result6 = Policy::Initialize(10, 10, d6);
            BOOST_CHECK_EQUAL(result6.num_elements(), 120);

            DataHolderType d7(8,6);
            Array<OneD, NekDouble> result7 = Policy::Initialize(10, 10, d7);
            BOOST_CHECK_EQUAL(result7.num_elements(), 150);
        }

        BOOST_AUTO_TEST_CASE(TestSingleValuePopulationInitialize)
        {
            typedef MatrixStoragePolicy<CountedObject<NekDouble>, BandedMatrixTag> Policy;
            Policy::PolicySpecificDataHolderType policyData(2, 1);
            {
                CountedObject<NekDouble> initValue(7);
                CountedObject<NekDouble>::ClearCounters();
                Array<OneD, CountedObject<NekDouble> > result = Policy::Initialize(5, 5, 
                    initValue, policyData);
                BOOST_CHECK_EQUAL(result.num_elements(), 20);
                CountedObject<NekDouble>::check(0, 0, 0, 20, 0, 0);
                for(Array<OneD, CountedObject<NekDouble> >::iterator iter = result.begin(); iter != result.end(); ++iter)
                {
                    BOOST_CHECK_EQUAL(*iter, CountedObject<NekDouble>(7) );
                }
            }

            {
                BOOST_CHECK_THROW( (Policy::Initialize(5,6, CountedObject<NekDouble>(7), policyData)), ErrorUtil::NekError);
                BOOST_CHECK_THROW( (Policy::Initialize(6,5,CountedObject<NekDouble>(7), policyData)), ErrorUtil::NekError);
            }
        }

        BOOST_AUTO_TEST_CASE(TestCArrayInitialization)
        {
            DataHolderType policyData(2, 1);
            NekDouble buf[] = {0, 1, 3, 6,
                               2, 4, 7, 0,
                               5, 8, 0, 0};

            {
                Array<OneD, NekDouble> result = Policy::Initialize(3, 3, buf, policyData);
                BOOST_CHECK_EQUAL(result.num_elements(), 12);

                BOOST_CHECK_EQUAL(result[0], 0.0);
                BOOST_CHECK_EQUAL(result[1], 1.0);
                BOOST_CHECK_EQUAL(result[2], 3.0);
                BOOST_CHECK_EQUAL(result[3], 6.0);
                BOOST_CHECK_EQUAL(result[4], 2.0);
                BOOST_CHECK_EQUAL(result[5], 4.0);
                BOOST_CHECK_EQUAL(result[6], 7.0);
                BOOST_CHECK_EQUAL(result[7], 0.0);
                BOOST_CHECK_EQUAL(result[8], 5.0);
                BOOST_CHECK_EQUAL(result[9], 8.0);
                BOOST_CHECK_EQUAL(result[10], 0.0);
                BOOST_CHECK_EQUAL(result[11], 0.0);
            }

            
            {
                NekDouble buf[] = {0, 1, 3, 6,
                                   2, 4, 7, 0,
                                   5, 8, 0, 0};

                BOOST_CHECK_THROW(Policy::Initialize(3, 2, buf, policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::Initialize(2, 3, buf, policyData), ErrorUtil::NekError);
            }
        }

        BOOST_AUTO_TEST_CASE(TestArrayInitialization)
        {
            DataHolderType policyData(2, 1);
            {
                NekDouble buf[] = {0, 1, 3, 6,
                               2, 4, 7, 0,
                               5, 8, 0, 0};
                Array<OneD, NekDouble> array_buf(12, buf);

                Array<OneD, NekDouble> result = Policy::Initialize(3, 3, array_buf, policyData);
                BOOST_CHECK_EQUAL(result.num_elements(), 12);

                BOOST_CHECK_EQUAL(result[0], 0.0);
                BOOST_CHECK_EQUAL(result[1], 1.0);
                BOOST_CHECK_EQUAL(result[2], 3.0);
                BOOST_CHECK_EQUAL(result[3], 6.0);
                BOOST_CHECK_EQUAL(result[4], 2.0);
                BOOST_CHECK_EQUAL(result[5], 4.0);
                BOOST_CHECK_EQUAL(result[6], 7.0);
                BOOST_CHECK_EQUAL(result[7], 0.0);
                BOOST_CHECK_EQUAL(result[8], 5.0);
                BOOST_CHECK_EQUAL(result[9], 8.0);
                BOOST_CHECK_EQUAL(result[10], 0.0);
                BOOST_CHECK_EQUAL(result[11], 0.0);
            }

            {
                NekDouble buf[] = {0, 1, 3, 6,
                               2, 4, 7, 0,
                               5, 8, 0, 0, 10, 11, 12};
                Array<OneD, NekDouble> array_buf(15, buf);

                Array<OneD, NekDouble> result = Policy::Initialize(3, 3, array_buf, policyData);
                BOOST_CHECK_EQUAL(result.num_elements(), 12);

                BOOST_CHECK_EQUAL(result[0], 0.0);
                BOOST_CHECK_EQUAL(result[1], 1.0);
                BOOST_CHECK_EQUAL(result[2], 3.0);
                BOOST_CHECK_EQUAL(result[3], 6.0);
                BOOST_CHECK_EQUAL(result[4], 2.0);
                BOOST_CHECK_EQUAL(result[5], 4.0);
                BOOST_CHECK_EQUAL(result[6], 7.0);
                BOOST_CHECK_EQUAL(result[7], 0.0);
                BOOST_CHECK_EQUAL(result[8], 5.0);
                BOOST_CHECK_EQUAL(result[9], 8.0);
                BOOST_CHECK_EQUAL(result[10], 0.0);
                BOOST_CHECK_EQUAL(result[11], 0.0);
            }

            {
                NekDouble buf[] = {0, 1, 3, 6,
                                   2, 4, 7, 0,
                                   5, 8, 0, 0};
                Array<OneD, NekDouble> array_buf(12, buf);

                BOOST_CHECK_THROW(Policy::Initialize(4, 4, array_buf, policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::Initialize(3, 2, array_buf, policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::Initialize(2, 3, array_buf, policyData), ErrorUtil::NekError);
            }
        }


        BOOST_AUTO_TEST_CASE(TestCalculateStorageSizeAndCalculateNumberOfRows)
        {
            DataHolderType policyData(2, 1);
            BOOST_CHECK_EQUAL(16, Policy::CalculateStorageSize(4, 4, policyData));
            BOOST_CHECK_EQUAL(4, Policy::CalculateNumberOfRows(4, policyData));
            BOOST_CHECK_EQUAL(4, Policy::CalculateNumberOfRows(3, policyData));
        }

        BOOST_AUTO_TEST_CASE(TestDiagonalOnlyCalculateIndex)
        {
            DataHolderType policyData(0, 0);
            BOOST_CHECK_EQUAL(0, *Policy::CalculateIndex(3, 3, 0, 0, policyData));
            BOOST_CHECK_EQUAL(1, *Policy::CalculateIndex(3, 3, 1, 1, policyData));
            BOOST_CHECK_EQUAL(2, *Policy::CalculateIndex(3, 3, 2, 2, policyData));
            BOOST_CHECK_EQUAL(0, *Policy::CalculateIndex(1, 1, 0, 0, policyData));

            BOOST_CHECK(!Policy::CalculateIndex(3, 3, 0, 1, policyData)); 
            BOOST_CHECK(!Policy::CalculateIndex(3, 3, 0, 2, policyData)); 
            BOOST_CHECK(!Policy::CalculateIndex(3, 3, 1, 0, policyData)); 
            BOOST_CHECK(!Policy::CalculateIndex(3, 3, 1, 2, policyData)); 
            BOOST_CHECK(!Policy::CalculateIndex(3, 3, 2, 0, policyData)); 
            BOOST_CHECK(!Policy::CalculateIndex(3, 3, 2, 1, policyData)); 
        }

        BOOST_AUTO_TEST_CASE(TestSubDiagonalsOnlyCalculateIndex)
        {
            DataHolderType oneSubData(1, 0);
            BOOST_CHECK_EQUAL(0, *Policy::CalculateIndex(3, 3, 0, 0, oneSubData));
            BOOST_CHECK_EQUAL(2, *Policy::CalculateIndex(3, 3, 1, 1, oneSubData));
            BOOST_CHECK_EQUAL(4, *Policy::CalculateIndex(3, 3, 2, 2, oneSubData));
            BOOST_CHECK_EQUAL(1, *Policy::CalculateIndex(3, 3, 1, 0, oneSubData));
            BOOST_CHECK_EQUAL(3, *Policy::CalculateIndex(3, 3, 2, 1, oneSubData));

            BOOST_CHECK(!Policy::CalculateIndex(3, 3, 0, 1, oneSubData));
            BOOST_CHECK(!Policy::CalculateIndex(3, 3, 0, 2, oneSubData));
            BOOST_CHECK(!Policy::CalculateIndex(3, 3, 2, 0, oneSubData));
            BOOST_CHECK(!Policy::CalculateIndex(3, 3, 1, 2, oneSubData));
        }

        BOOST_AUTO_TEST_CASE(TestSuperDiagonalsOnlyCalculateIndex)
        {
            DataHolderType oneSubData(0, 2);
            BOOST_CHECK_EQUAL(2, *Policy::CalculateIndex(3, 3, 0, 0, oneSubData));
            BOOST_CHECK_EQUAL(5, *Policy::CalculateIndex(3, 3, 1, 1, oneSubData));
            BOOST_CHECK_EQUAL(8, *Policy::CalculateIndex(3, 3, 2, 2, oneSubData));
            BOOST_CHECK_EQUAL(4, *Policy::CalculateIndex(3, 3, 0, 1, oneSubData));
            BOOST_CHECK_EQUAL(7, *Policy::CalculateIndex(3, 3, 1, 2, oneSubData));
            BOOST_CHECK_EQUAL(6, *Policy::CalculateIndex(3, 3, 0, 2, oneSubData));

            BOOST_CHECK(!Policy::CalculateIndex(3, 3, 1, 0, oneSubData));
            BOOST_CHECK(!Policy::CalculateIndex(3, 3, 2, 0, oneSubData));
            BOOST_CHECK(!Policy::CalculateIndex(3, 3, 2, 1, oneSubData));
        }

        BOOST_AUTO_TEST_CASE(TestSupAndSuperDiagonalsCalculateIndex)
        {
            DataHolderType oneSubData(1, 2);

            BOOST_CHECK_EQUAL(2, *Policy::CalculateIndex(4, 4, 0, 0, oneSubData));
            BOOST_CHECK_EQUAL(5, *Policy::CalculateIndex(4, 4, 0, 1, oneSubData));
            BOOST_CHECK_EQUAL(8, *Policy::CalculateIndex(4, 4, 0, 2, oneSubData));
            BOOST_CHECK(!Policy::CalculateIndex(4, 4, 0, 3, oneSubData));

            BOOST_CHECK_EQUAL(3, *Policy::CalculateIndex(4, 4, 1, 0, oneSubData));
            BOOST_CHECK_EQUAL(6, *Policy::CalculateIndex(4, 4, 1, 1, oneSubData));
            BOOST_CHECK_EQUAL(9, *Policy::CalculateIndex(4, 4, 1, 2, oneSubData));
            BOOST_CHECK_EQUAL(12, *Policy::CalculateIndex(4, 4, 1, 3, oneSubData));

            BOOST_CHECK(!Policy::CalculateIndex(4, 4, 2, 0, oneSubData));
            BOOST_CHECK_EQUAL(7, *Policy::CalculateIndex(4, 4, 2, 1, oneSubData));
            BOOST_CHECK_EQUAL(10, *Policy::CalculateIndex(4, 4, 2, 2, oneSubData));
            BOOST_CHECK_EQUAL(13, *Policy::CalculateIndex(4, 4, 2, 3, oneSubData));

            BOOST_CHECK(!Policy::CalculateIndex(4, 4, 3, 0, oneSubData));
            BOOST_CHECK(!Policy::CalculateIndex(4, 4, 3, 1, oneSubData));
            BOOST_CHECK_EQUAL(11, *Policy::CalculateIndex(4, 4, 3, 2, oneSubData));
            BOOST_CHECK_EQUAL(14, *Policy::CalculateIndex(4, 4, 3, 3, oneSubData));
        }


        BOOST_AUTO_TEST_CASE(TestGetValue)
        {
            typedef MatrixStoragePolicy<CountedObject<NekDouble>, BandedMatrixTag> Policy;
            Policy::PolicySpecificDataHolderType policyData(1, 2);

            // [ 1 2 3 0 ]
            // [ 4 5 6 7 ]
            // [ 0 8 9 10 ]
            // [ 0 0 11 12 ]
            CountedObject<NekDouble> buf[] = { CountedObject<NekDouble>(0), CountedObject<NekDouble>(0), CountedObject<NekDouble>(1), CountedObject<NekDouble>(4),
                                               CountedObject<NekDouble>(0), CountedObject<NekDouble>(2), CountedObject<NekDouble>(5), CountedObject<NekDouble>(8),
                                               CountedObject<NekDouble>(3), CountedObject<NekDouble>(6), CountedObject<NekDouble>(9), CountedObject<NekDouble>(11),
                                               CountedObject<NekDouble>(7), CountedObject<NekDouble>(10), CountedObject<NekDouble>(12), CountedObject<NekDouble>(0) };
            NekMatrix<CountedObject<NekDouble>, BandedMatrixTag> m(4, 4, buf, policyData);

            CountedObject<NekDouble>::ClearCounters();

            BOOST_CHECK_EQUAL(1, m(0,0).GetValue());
            BOOST_CHECK_EQUAL(2, m(0,1).GetValue());
            BOOST_CHECK_EQUAL(3, m(0,2).GetValue());
            BOOST_CHECK_EQUAL(0, m(0,3).GetValue());

            BOOST_CHECK_EQUAL(4, m(1,0).GetValue());
            BOOST_CHECK_EQUAL(5, m(1,1).GetValue());
            BOOST_CHECK_EQUAL(6, m(1,2).GetValue());
            BOOST_CHECK_EQUAL(7, m(1,3).GetValue());

            BOOST_CHECK_EQUAL(0, m(2,0).GetValue());
            BOOST_CHECK_EQUAL(8, m(2,1).GetValue());
            BOOST_CHECK_EQUAL(9, m(2,2).GetValue());
            BOOST_CHECK_EQUAL(10, m(2,3).GetValue());

            BOOST_CHECK_EQUAL(0, m(3,0).GetValue());
            BOOST_CHECK_EQUAL(0, m(3,1).GetValue());
            BOOST_CHECK_EQUAL(11, m(3,2).GetValue());
            BOOST_CHECK_EQUAL(12, m(3,3).GetValue());

            CountedObject<NekDouble>::check(0, 0, 0, 0, 0, 0);

        }
    }
}




