///////////////////////////////////////////////////////////////////////////////
//
// File: TestDiagonalMatrixStoragePolicy.cpp
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
    namespace DiagonalMatrixStoragePolicyUnitTests
    {
        typedef MatrixStoragePolicy<CountedObject<double>, DiagonalMatrixTag> Policy;

        BOOST_AUTO_TEST_CASE(Test0ParameterInitialize)
        {
            Array<OneD, CountedObject<double> > result = Policy::Initialize();
            BOOST_CHECK_EQUAL(result.num_elements(), 0);
            BOOST_CHECK_EQUAL(result.GetOffset(), 0);
            BOOST_CHECK(result.data() != 0);
        }

        BOOST_AUTO_TEST_CASE(Test2ParameterInitialize)
        {
            Policy::PolicySpecificDataHolderType policyData;
            CountedObject<double>::ClearCounters();

            {
                Array<OneD, CountedObject<double> > result = Policy::Initialize(5, 5, policyData);
                BOOST_CHECK_EQUAL(result.num_elements(), 5);
                CountedObject<double>::check(5, 0, 0, 0, 0, 0);
            }

            {
                BOOST_CHECK_THROW(Policy::Initialize(5, 6, policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::Initialize(6, 5, policyData), ErrorUtil::NekError);
            }
        }

        BOOST_AUTO_TEST_CASE(TestSingleValuePopulationInitialize)
        {
            Policy::PolicySpecificDataHolderType policyData;
            CountedObject<double> initValue(1);

            {    
                CountedObject<double>::ClearCounters();
                Array<OneD, CountedObject<double> > result = Policy::Initialize(5, 5, initValue, policyData);
                BOOST_CHECK_EQUAL(result.num_elements(), 5);
                for(Array<OneD, CountedObject<double> >::iterator iter = result.begin(); iter != result.end(); ++iter)
                {
                    BOOST_CHECK_EQUAL(*iter, initValue);
                }
                CountedObject<double>::check(0, 0, 0, 5, 0, 0);
            }

            {
                BOOST_CHECK_THROW(Policy::Initialize(5, 6, initValue, policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::Initialize(6, 5, initValue, policyData), ErrorUtil::NekError);
            }
        }

        BOOST_AUTO_TEST_CASE(TestCArrayInitialization)
        {
            Policy::PolicySpecificDataHolderType policyData;
            CountedObject<double> buf[] = {CountedObject<double>(1), CountedObject<double>(2), CountedObject<double>(3) };

            {
                CountedObject<double>::ClearCounters();
                Array<OneD, CountedObject<double> > result = Policy::Initialize(3, 3, buf, policyData);
                BOOST_CHECK_EQUAL(result.num_elements(), 3);
                CountedObject<double>::check(0, 0, 0, 3, 0, 0);

                BOOST_CHECK_EQUAL(result[0], CountedObject<double>(1));
                BOOST_CHECK_EQUAL(result[1], CountedObject<double>(2));
                BOOST_CHECK_EQUAL(result[2], CountedObject<double>(3));
            }

            {
                BOOST_CHECK_THROW(Policy::Initialize(3, 4, buf, policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::Initialize(4, 3, buf, policyData), ErrorUtil::NekError);
            }
        }

        BOOST_AUTO_TEST_CASE(TestArrayInitialization)
        {
            Policy::PolicySpecificDataHolderType policyData;
            CountedObject<double> buf[] = {CountedObject<double>(1), CountedObject<double>(2), CountedObject<double>(3) };
            Array<OneD, CountedObject<double> > array_buf(3, buf);
            
            {
                CountedObject<double>::ClearCounters();
                Array<OneD, CountedObject<double> > result = Policy::Initialize(3, 3, array_buf, policyData);
                BOOST_CHECK_EQUAL(result.num_elements(), 3);
                CountedObject<double>::check(3, 0, 0, 0, 0, 3);

                BOOST_CHECK_EQUAL(result[0], CountedObject<double>(1));
                BOOST_CHECK_EQUAL(result[1], CountedObject<double>(2));
                BOOST_CHECK_EQUAL(result[2], CountedObject<double>(3));
            }

            {
                CountedObject<double>::ClearCounters();
                Array<OneD, CountedObject<double> > result = Policy::Initialize(2, 2, array_buf, policyData);
                BOOST_CHECK_EQUAL(result.num_elements(), 2);
                CountedObject<double>::check(2, 0, 0, 0, 0, 2);

                BOOST_CHECK_EQUAL(result[0], CountedObject<double>(1));
                BOOST_CHECK_EQUAL(result[1], CountedObject<double>(2));
            }

            {
                BOOST_CHECK_THROW(Policy::Initialize(4, 4, array_buf, policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::Initialize(4, 3, array_buf, policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::Initialize(3, 4, array_buf, policyData), ErrorUtil::NekError);
            }
        }

        BOOST_AUTO_TEST_CASE(TestGetValue)
        {
            typedef MatrixStoragePolicy<double, DiagonalMatrixTag> Policy;
            Policy::PolicySpecificDataHolderType policyData;
            NekDouble buf[] = {1.0, 2.0, 3.0};
            Array<OneD, NekDouble> data = Policy::Initialize(3, 3, buf, policyData);
            ConstArray<OneD, NekDouble>& cdata = data;
            BOOST_CHECK_EQUAL(data.num_elements(), 3);
            {
                BOOST_CHECK_EQUAL(1.0, Policy::GetValue(3, 3, 0, 0, data, 'N', policyData));
                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 0, 1, data, 'N', policyData));
                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 0, 2, data, 'N', policyData));

                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 1, 0, data, 'N', policyData));
                BOOST_CHECK_EQUAL(2.0, Policy::GetValue(3, 3, 1, 1, data, 'N', policyData));
                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 1, 2, data, 'N', policyData));

                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 2, 0, data, 'N', policyData));
                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 2, 1, data, 'N', policyData));
                BOOST_CHECK_EQUAL(3.0, Policy::GetValue(3, 3, 2, 2, data, 'N', policyData));
                

                BOOST_CHECK_EQUAL(1.0, Policy::GetValue(3, 3, 0, 0, cdata, 'N', policyData));
                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 0, 1, cdata, 'N', policyData));
                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 0, 2, cdata, 'N', policyData));

                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 1, 0, cdata, 'N', policyData));
                BOOST_CHECK_EQUAL(2.0, Policy::GetValue(3, 3, 1, 1, cdata, 'N', policyData));
                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 1, 2, cdata, 'N', policyData));

                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 2, 0, cdata, 'N', policyData));
                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 2, 1, cdata, 'N', policyData));
                BOOST_CHECK_EQUAL(3.0, Policy::GetValue(3, 3, 2, 2, cdata, 'N', policyData));
            }

            #if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)
            {
                BOOST_CHECK_THROW(Policy::GetValue(3,3,3,4,data, 'N', policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::GetValue(3,3,4,3,data, 'N', policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::GetValue(3,3,4,4,data, 'N', policyData), ErrorUtil::NekError);
            }
            #endif
        }

        BOOST_AUTO_TEST_CASE(TestSetValue)
        {
            typedef MatrixStoragePolicy<double, DiagonalMatrixTag> Policy;
            Policy::PolicySpecificDataHolderType policyData;
            NekDouble buf[] = {1.0, 2.0, 3.0};
            Array<OneD, NekDouble> data = Policy::Initialize(3, 3, buf, policyData);
            BOOST_CHECK_EQUAL(data.num_elements(), 3);

            {
                BOOST_CHECK_EQUAL(1.0, Policy::GetValue(3, 3, 0, 0, data, 'N', policyData));
                Policy::SetValue(3, 3, 0, 0, data, 100.0, 'N', policyData);
                BOOST_CHECK_EQUAL(100.0, Policy::GetValue(3, 3, 0, 0, data, 'N', policyData));

                BOOST_CHECK_EQUAL(2.0, Policy::GetValue(3, 3, 1, 1, data, 'N', policyData));
                Policy::SetValue(3, 3, 1, 1, data, 101.0, 'N', policyData);
                BOOST_CHECK_EQUAL(101.0, Policy::GetValue(3, 3, 1, 1, data, 'N', policyData));

                BOOST_CHECK_EQUAL(3.0, Policy::GetValue(3, 3, 2, 2, data, 'N', policyData));
                Policy::SetValue(3, 3, 2, 2, data, 103.0, 'N', policyData);
                BOOST_CHECK_EQUAL(103.0, Policy::GetValue(3, 3, 2, 2, data, 'N', policyData));
            }

            {
                BOOST_CHECK_THROW(Policy::SetValue(3, 3, 0, 1, data, 1.0, 'N', policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::SetValue(3, 3, 0, 2, data, 1.0, 'N', policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::SetValue(3, 3, 1, 0, data, 1.0, 'N', policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::SetValue(3, 3, 1, 2, data, 1.0, 'N', policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::SetValue(3, 3, 2, 0, data, 1.0, 'N', policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::SetValue(3, 3, 2, 1, data, 1.0, 'N', policyData), ErrorUtil::NekError);
            }

            #if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)
            {
                BOOST_CHECK_THROW(Policy::SetValue(3,3,3,3,data, 1.0, 'N', policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::SetValue(3,3,4,3,data, 1.0, 'N', policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::SetValue(3,3,4,4,data, 1.0, 'N', policyData), ErrorUtil::NekError);
            }
            #endif
        }
                    
        BOOST_AUTO_TEST_CASE(TestAdvance)
        {
            typedef MatrixStoragePolicy<double, DiagonalMatrixTag> Policy;
            Policy::PolicySpecificDataHolderType policyData;

            {
                NekDouble buf[] = {1.0, 2.0, 3.0, 4.0 };
                Array<OneD, NekDouble> data = Policy::Initialize(4, 4, buf, policyData);

                unsigned int curRow = 0; 
                unsigned int curColumn = 0;
                boost::tie(curRow, curColumn) = Policy::Advance(4, 4, curRow, curColumn, 'N', policyData);
                BOOST_CHECK_EQUAL(1, curRow);
                BOOST_CHECK_EQUAL(1, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(4, 4, curRow, curColumn, 'N', policyData);
                BOOST_CHECK_EQUAL(2, curRow);
                BOOST_CHECK_EQUAL(2, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(4, 4, curRow, curColumn, 'N', policyData);
                BOOST_CHECK_EQUAL(3, curRow);
                BOOST_CHECK_EQUAL(3, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(2, 2, curRow, curColumn, 'N', policyData);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curRow);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curColumn);
            }

            {
                NekDouble buf[] = {1.0};
                Array<OneD, NekDouble> data = Policy::Initialize(1, 1, buf, policyData);

                unsigned int curRow = 0; 
                unsigned int curColumn = 0;
                boost::tie(curRow, curColumn) = Policy::Advance(1, 1, curRow, curColumn, 'N', policyData);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curRow);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curColumn);
            }
        }

        BOOST_AUTO_TEST_CASE(TestOneByOneMatrixInversion)
        {
            typedef MatrixStoragePolicy<double, FullMatrixTag> Policy;
            Policy::PolicySpecificDataHolderType policyData;

            NekDouble buf[] = {8.0};
            Array<OneD, NekDouble> data = Policy::Initialize(1, 1, buf, policyData);

            Policy::Invert(1, 1, data, 'N', policyData);
            BOOST_CHECK_EQUAL(1.0/8.0, Policy::GetValue(1, 1, 0, 0, data, 'N', policyData));
            Policy::Invert(1, 1, data, 'N', policyData);
            BOOST_CHECK_EQUAL(8.0, Policy::GetValue(1, 1, 0, 0, data, 'N', policyData));
        }

        BOOST_AUTO_TEST_CASE(TestTwoByTwoMatrixInversion)
        {
            typedef MatrixStoragePolicy<double, DiagonalMatrixTag> Policy;
            Policy::PolicySpecificDataHolderType policyData;

            NekDouble buf[] = {8.0, 10.0, 17.0, 2.0};
            Array<OneD, NekDouble> data = Policy::Initialize(4, 4, buf, policyData);
            Policy::Invert(4, 4, data, 'N', policyData);

            BOOST_CHECK_CLOSE(1.0/8.0, Policy::GetValue(4, 4, 0, 0, data, 'N', policyData), .0001);
            BOOST_CHECK_CLOSE(1.0/10.0, Policy::GetValue(4, 4, 1, 1, data, 'N', policyData), .0001);
            BOOST_CHECK_CLOSE(1.0/17.0, Policy::GetValue(4, 4, 2, 2, data, 'N', policyData), .0001);
            BOOST_CHECK_CLOSE(1.0/2.0, Policy::GetValue(4, 4, 3, 3, data, 'N', policyData), .0001);
        }
    }
}


