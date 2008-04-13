///////////////////////////////////////////////////////////////////////////////
//
// File: TestUpperTriangularMatrixStoragePolicy.cpp
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
    namespace UpperTriangularUnitTests
    {
        typedef MatrixStoragePolicy<NekDouble, UpperTriangularMatrixTag> Policy;

        BOOST_AUTO_TEST_CASE(Test0ParameterInitialize)
        {
            Array<OneD, NekDouble> result = Policy::Initialize();
            BOOST_CHECK_EQUAL(result.num_elements(), 0);
            BOOST_CHECK_EQUAL(result.GetOffset(), 0);
            BOOST_CHECK(result.data() == 0);
        }

        BOOST_AUTO_TEST_CASE(Test2ParameterInitialize)
        {
            typedef MatrixStoragePolicy<CountedObject<NekDouble>, UpperTriangularMatrixTag> Policy;
            Policy::PolicySpecificDataHolderType policyData;
            {
                CountedObject<NekDouble>::ClearCounters();
                Array<OneD, CountedObject<NekDouble> > result = Policy::Initialize(5, 5, policyData);
                BOOST_CHECK_EQUAL(result.num_elements(), 15);
                CountedObject<NekDouble>::Check(15, 0, 0, 0, 0, 0);
            }

            {
                BOOST_CHECK_THROW( (Policy::Initialize(5,6, policyData)), ErrorUtil::NekError);
                BOOST_CHECK_THROW( (Policy::Initialize(6,5, policyData)), ErrorUtil::NekError);
            }
        }

        BOOST_AUTO_TEST_CASE(TestSingleValuePopulationInitialize)
        {
            typedef MatrixStoragePolicy<CountedObject<NekDouble>, UpperTriangularMatrixTag> Policy;
            Policy::PolicySpecificDataHolderType policyData;
            {
                CountedObject<NekDouble> initValue(7);
                CountedObject<NekDouble>::ClearCounters();
                Array<OneD, CountedObject<NekDouble> > result = Policy::Initialize(5, 5, 
                    initValue, policyData);
                BOOST_CHECK_EQUAL(result.num_elements(), 15);
                CountedObject<NekDouble>::Check(0, 0, 0, 15, 0, 0);
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
            typedef MatrixStoragePolicy<CountedObject<NekDouble>, UpperTriangularMatrixTag> Policy;
            Policy::PolicySpecificDataHolderType policyData;
            CountedObject<NekDouble> buf[] = {CountedObject<NekDouble>(1), 
                               CountedObject<NekDouble>(2), CountedObject<NekDouble>(3),
                               CountedObject<NekDouble>(5), CountedObject<NekDouble>(6), CountedObject<NekDouble>(9)};

            {
                CountedObject<NekDouble>::ClearCounters();
                Array<OneD, CountedObject<NekDouble> > result = Policy::Initialize(3, 3, buf, policyData);
                BOOST_CHECK_EQUAL(result.num_elements(), 6);
                CountedObject<NekDouble>::Check(0, 0, 0, 6, 0, 0);

                BOOST_CHECK_EQUAL(result[0], CountedObject<NekDouble>(1));
                BOOST_CHECK_EQUAL(result[1], CountedObject<NekDouble>(2));
                BOOST_CHECK_EQUAL(result[2], CountedObject<NekDouble>(3));
                BOOST_CHECK_EQUAL(result[3], CountedObject<NekDouble>(5));
                BOOST_CHECK_EQUAL(result[4], CountedObject<NekDouble>(6));
                BOOST_CHECK_EQUAL(result[5], CountedObject<NekDouble>(9));
            }

            {

                BOOST_CHECK_THROW( (Policy::Initialize(5,6,buf, policyData)), ErrorUtil::NekError);
                BOOST_CHECK_THROW( (Policy::Initialize(6,5,buf, policyData)), ErrorUtil::NekError);
            }
        }

        BOOST_AUTO_TEST_CASE(TestArrayInitialization)
        {
            typedef MatrixStoragePolicy<CountedObject<NekDouble>, UpperTriangularMatrixTag> Policy;
            Policy::PolicySpecificDataHolderType policyData;
            CountedObject<NekDouble> buf[] = {CountedObject<NekDouble>(1), 
                               CountedObject<NekDouble>(2), CountedObject<NekDouble>(3),
                               CountedObject<NekDouble>(5), CountedObject<NekDouble>(6), CountedObject<NekDouble>(9)};

            Array<OneD, CountedObject<NekDouble> > array_buf(6, buf);

            {
                CountedObject<NekDouble>::ClearCounters();
                Array<OneD, CountedObject<NekDouble> > result = Policy::Initialize(3, 3, array_buf, policyData);
                BOOST_CHECK_EQUAL(result.num_elements(), 6);
                CountedObject<NekDouble>::Check(0, 0, 0, 6, 0, 0);

                BOOST_CHECK_EQUAL(result[0], CountedObject<NekDouble>(1));
                BOOST_CHECK_EQUAL(result[1], CountedObject<NekDouble>(2));
                BOOST_CHECK_EQUAL(result[2], CountedObject<NekDouble>(3));
                BOOST_CHECK_EQUAL(result[3], CountedObject<NekDouble>(5));
                BOOST_CHECK_EQUAL(result[4], CountedObject<NekDouble>(6));
                BOOST_CHECK_EQUAL(result[5], CountedObject<NekDouble>(9));
            }

            {

                BOOST_CHECK_THROW( (Policy::Initialize(5,6,buf, policyData)), ErrorUtil::NekError);
                BOOST_CHECK_THROW( (Policy::Initialize(6,5,buf, policyData)), ErrorUtil::NekError);
            }

            {
                CountedObject<NekDouble> buf[] = {CountedObject<NekDouble>(1), 
                               CountedObject<NekDouble>(2), CountedObject<NekDouble>(3),
                               CountedObject<NekDouble>(5), CountedObject<NekDouble>(6), CountedObject<NekDouble>(9),
                               CountedObject<NekDouble>(10) };
                Array<OneD, CountedObject<NekDouble> > array_buf(7, buf);
                Array<OneD, CountedObject<NekDouble> > result;
                CountedObject<NekDouble>::ClearCounters();
                BOOST_CHECK_NO_THROW(result = Policy::Initialize(3, 3, array_buf, policyData));
                BOOST_CHECK_EQUAL(6, result.num_elements());
                CountedObject<NekDouble>::Check(0, 0, 0, 6, 0, 0);
            }

            {
                CountedObject<NekDouble> buf[] = {CountedObject<NekDouble>(1), 
                               CountedObject<NekDouble>(2), CountedObject<NekDouble>(3),
                               CountedObject<NekDouble>(5), CountedObject<NekDouble>(6), CountedObject<NekDouble>(9),
                               CountedObject<NekDouble>(10) };
                Array<OneD, CountedObject<NekDouble> > array_buf(5, buf);
                BOOST_CHECK_THROW(Policy::Initialize(3, 3, array_buf, policyData), ErrorUtil::NekError);
            }
        }



        BOOST_AUTO_TEST_CASE(TestGetValue)
        {
            Policy::PolicySpecificDataHolderType policyData;
            NekDouble buf[] = {1.0, 
                               2.0, 3.0,
                               5.0, 6.0, 9.0};
            Array<OneD, NekDouble> array_buf(6, buf);
            Array<OneD, NekDouble> data = Policy::Initialize(3, 3, array_buf, policyData);
            Array<OneD, const NekDouble>& cdata = data;
            BOOST_CHECK_EQUAL(data.num_elements(), 6);

            {
                BOOST_CHECK_EQUAL(1.0, Policy::GetValue(3, 3, 0, 0, data, 'N', policyData));
                BOOST_CHECK_EQUAL(2.0, Policy::GetValue(3, 3, 0, 1, data, 'N', policyData));
                BOOST_CHECK_EQUAL(5.0, Policy::GetValue(3, 3, 0, 2, data, 'N', policyData));
                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 1, 0, data, 'N', policyData));
                BOOST_CHECK_EQUAL(3.0, Policy::GetValue(3, 3, 1, 1, data, 'N', policyData));
                BOOST_CHECK_EQUAL(6.0, Policy::GetValue(3, 3, 1, 2, data, 'N', policyData));
                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 2, 0, data, 'N', policyData));
                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 2, 1, data, 'N', policyData));
                BOOST_CHECK_EQUAL(9.0, Policy::GetValue(3, 3, 2, 2, data, 'N', policyData));

                BOOST_CHECK_EQUAL(1.0, Policy::GetValue(3, 3, 0, 0, cdata, 'N', policyData));
                BOOST_CHECK_EQUAL(2.0, Policy::GetValue(3, 3, 0, 1, cdata, 'N', policyData));
                BOOST_CHECK_EQUAL(5.0, Policy::GetValue(3, 3, 0, 2, cdata, 'N', policyData));
                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 1, 0, cdata, 'N', policyData));
                BOOST_CHECK_EQUAL(3.0, Policy::GetValue(3, 3, 1, 1, cdata, 'N', policyData));
                BOOST_CHECK_EQUAL(6.0, Policy::GetValue(3, 3, 1, 2, cdata, 'N', policyData));
                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 2, 0, cdata, 'N', policyData));
                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 2, 1, cdata, 'N', policyData));
                BOOST_CHECK_EQUAL(9.0, Policy::GetValue(3, 3, 2, 2, cdata, 'N', policyData));
            }

            #if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)
            {
                BOOST_CHECK_THROW(Policy::GetValue(3,4,3,3,data, 'N', policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::GetValue(3,3,4,3,data, 'N', policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::GetValue(3,3,3,4,data, 'N', policyData), ErrorUtil::NekError);
            }
            #endif
        }

        BOOST_AUTO_TEST_CASE(TestSetValue)
        {
            Policy::PolicySpecificDataHolderType policyData;

            NekDouble buf[] = {1.0, 
                               2.0, 3.0,
                               5.0, 6.0, 9.0};
            Array<OneD, NekDouble> array_buf(6, buf);
            Array<OneD, NekDouble> data = Policy::Initialize(3, 3, array_buf, policyData);
            BOOST_CHECK_EQUAL(data.num_elements(), 6);
            
            Policy::SetValue(3,3,0,0,data, 7.2, 'N', policyData);
            BOOST_CHECK_EQUAL(7.2, Policy::GetValue(3,3,0,0,data, 'N', policyData));
            BOOST_CHECK_EQUAL(1.0, array_buf[0]);
            BOOST_CHECK_EQUAL(7.2, data[0]);

            Policy::SetValue(3,3,0,1,data,10.0, 'N', policyData);
            BOOST_CHECK_EQUAL(10.0, Policy::GetValue(3,3,0,1,data, 'N', policyData));
            Policy::SetValue(3,3,0,2,data,20.0, 'N', policyData);
            BOOST_CHECK_EQUAL(20.0, Policy::GetValue(3,3,0,2,data, 'N', policyData));
            Policy::SetValue(3,3,1,1,data,30.0, 'N', policyData);
            BOOST_CHECK_EQUAL(30.0, Policy::GetValue(3,3,1,1,data, 'N', policyData));
            Policy::SetValue(3,3,1,2,data,40.0, 'N', policyData);
            BOOST_CHECK_EQUAL(40.0, Policy::GetValue(3,3,1,2,data, 'N', policyData));
            Policy::SetValue(3,3,2,2,data,50.0, 'N', policyData);
            BOOST_CHECK_EQUAL(50.0, Policy::GetValue(3,3,2,2,data, 'N', policyData));

            #if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)
                // Can't set in the lower triangle.
                BOOST_CHECK_THROW(Policy::SetValue(3,3,1,0,data,8.0, 'N', policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::SetValue(3,3,2,0,data,8.0, 'N', policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::SetValue(3,3,2,1,data,8.0, 'N', policyData), ErrorUtil::NekError);

                BOOST_CHECK_THROW(Policy::SetValue(3,4,3,3,data,8.0, 'N', policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::SetValue(3,3,4,3,data,8.0, 'N', policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::SetValue(3,3,3,4,data,8.0, 'N', policyData), ErrorUtil::NekError);
            #endif
        }
                    
        BOOST_AUTO_TEST_CASE(TestAdvance)
        {
            Policy::PolicySpecificDataHolderType policyData;

            {    
                NekDouble buf[] = {1.0, 
                                   2.0, 3.0,
                                   5.0, 6.0, 9.0};
                Array<OneD, NekDouble> array_buf(6, buf);
                Array<OneD, NekDouble> data = Policy::Initialize(3, 3, array_buf, policyData);
                BOOST_CHECK_EQUAL(data.num_elements(), 6);

                unsigned int curRow = 0; 
                unsigned int curColumn = 0;
                boost::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn, policyData);
                BOOST_CHECK_EQUAL(0, curRow);
                BOOST_CHECK_EQUAL(1, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn, policyData);
                BOOST_CHECK_EQUAL(1, curRow);
                BOOST_CHECK_EQUAL(1, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn, policyData);
                BOOST_CHECK_EQUAL(0, curRow);
                BOOST_CHECK_EQUAL(2, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn, policyData);
                BOOST_CHECK_EQUAL(1, curRow);
                BOOST_CHECK_EQUAL(2, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn, policyData);
                BOOST_CHECK_EQUAL(2, curRow);
                BOOST_CHECK_EQUAL(2, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn, policyData);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curRow);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curColumn);
            }

            {
                NekDouble buf[] = {1.0};
                Array<OneD, NekDouble> array_buf(1, buf);
                Array<OneD, NekDouble> data = Policy::Initialize(1, 1, array_buf, policyData);
                BOOST_CHECK_EQUAL(data.num_elements(), 1);

                unsigned int curRow = 0; 
                unsigned int curColumn = 0;
                boost::tie(curRow, curColumn) = Policy::Advance(1, 1, curRow, curColumn, policyData);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curRow);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curColumn);
            }

            {
                NekDouble buf[] = {1.0, 2.0,
                                        5.0};
                Array<OneD, NekDouble> array_buf(3, buf);
                Array<OneD, NekDouble> data = Policy::Initialize(2, 2, array_buf, policyData);
                BOOST_CHECK_EQUAL(data.num_elements(), 3);

                unsigned int curRow = 0; 
                unsigned int curColumn = 0;
                boost::tie(curRow, curColumn) = Policy::Advance(2, 2, curRow, curColumn, policyData);
                BOOST_CHECK_EQUAL(0, curRow);
                BOOST_CHECK_EQUAL(1, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(2, 2, curRow, curColumn, policyData);
                BOOST_CHECK_EQUAL(1, curRow);
                BOOST_CHECK_EQUAL(1, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(2, 2, curRow, curColumn, policyData);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curRow);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curColumn);
            }
        }
    }
}


