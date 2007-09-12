///////////////////////////////////////////////////////////////////////////////
//
// File: TestSymmetricMatrixStoragePolicy.cpp
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
    namespace SymmetricMatrixStoragePolicyUnitTests
    {
        typedef MatrixStoragePolicy<NekDouble, SymmetricMatrixTag> Policy;

        BOOST_AUTO_TEST_CASE(TestConstGetValue)
        {
            Policy::PolicySpecificDataHolderType policyData;
            NekDouble buf[] = {1.0, 2.0, 3.0,
                                    5.0, 6.0,
                                         9.0};
            Array<OneD, NekDouble> array_buf(6, buf);
            Array<OneD, NekDouble> data = Policy::Initialize(3, 3, array_buf, policyData);
            ConstArray<OneD, NekDouble>& cdata = data;

            BOOST_CHECK_EQUAL(data.num_elements(), 6);

            {
                BOOST_CHECK_EQUAL(1.0, Policy::GetValue(3, 3, 0, 0, data, policyData));
                BOOST_CHECK_EQUAL(2.0, Policy::GetValue(3, 3, 0, 1, data, policyData));
                BOOST_CHECK_EQUAL(5.0, Policy::GetValue(3, 3, 0, 2, data, policyData));
                BOOST_CHECK_EQUAL(2.0, Policy::GetValue(3, 3, 1, 0, data, policyData));
                BOOST_CHECK_EQUAL(3.0, Policy::GetValue(3, 3, 1, 1, data, policyData));
                BOOST_CHECK_EQUAL(6.0, Policy::GetValue(3, 3, 1, 2, data, policyData));
                BOOST_CHECK_EQUAL(5.0, Policy::GetValue(3, 3, 2, 0, data, policyData));
                BOOST_CHECK_EQUAL(6.0, Policy::GetValue(3, 3, 2, 1, data, policyData));
                BOOST_CHECK_EQUAL(9.0, Policy::GetValue(3, 3, 2, 2, data, policyData));
            }

            #if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)
            {
                BOOST_CHECK_THROW(Policy::GetValue(3,4,3,3,data, policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::GetValue(3,3,4,3,data, policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::GetValue(3,3,3,4,data, policyData), ErrorUtil::NekError);
            }
            #endif
        }

        BOOST_AUTO_TEST_CASE(TestSetValue)
        {
            Policy::PolicySpecificDataHolderType policyData;
            NekDouble buf[] = {1.0, 2.0, 3.0,
                                    5.0, 6.0,
                                         9.0};
            Array<OneD, NekDouble> array_buf(6, buf);
            Array<OneD, NekDouble> data = Policy::Initialize(3, 3, array_buf, policyData);
            BOOST_CHECK_EQUAL(data.num_elements(), 6);
            
            Policy::SetValue(3,3,0,0,data, 7.2, policyData);
            BOOST_CHECK_EQUAL(7.2, Policy::GetValue(3,3,0,0,data, policyData));
            BOOST_CHECK_EQUAL(1.0, array_buf[0]);
            BOOST_CHECK_EQUAL(7.2, data[0]);

            Policy::SetValue(3,3,0,1,data,10.0, policyData);
            BOOST_CHECK_EQUAL(10.0, Policy::GetValue(3,3,0,1,data, policyData));
            BOOST_CHECK_EQUAL(10.0, Policy::GetValue(3,3,1,0,data, policyData));
            Policy::SetValue(3,3,0,2,data,20.0, policyData);
            BOOST_CHECK_EQUAL(20.0, Policy::GetValue(3,3,0,2,data, policyData));
            BOOST_CHECK_EQUAL(20.0, Policy::GetValue(3,3,2,0,data, policyData));
            Policy::SetValue(3,3,1,1,data,30.0, policyData);
            BOOST_CHECK_EQUAL(30.0, Policy::GetValue(3,3,1,1,data, policyData));
            Policy::SetValue(3,3,1,2,data,40.0, policyData);
            BOOST_CHECK_EQUAL(40.0, Policy::GetValue(3,3,1,2,data, policyData));
            BOOST_CHECK_EQUAL(40.0, Policy::GetValue(3,3,2,1,data, policyData));
            Policy::SetValue(3,3,2,2,data,50.0, policyData);
            BOOST_CHECK_EQUAL(50.0, Policy::GetValue(3,3,2,2,data, policyData));

            #if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)
                BOOST_CHECK_THROW(Policy::SetValue(3,3,4,3,data,8.0, policyData), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::SetValue(3,3,3,4,data,8.0, policyData), ErrorUtil::NekError);
            #endif
        }
                    
        BOOST_AUTO_TEST_CASE(TestAdvance)
        {
            Policy::PolicySpecificDataHolderType policyData;
            {
                NekDouble buf[] = {1.0, 2.0, 3.0,
                                        5.0, 6.0,
                                             9.0};
                Array<OneD, NekDouble> array_buf(6, buf);
                Array<OneD, NekDouble> data = Policy::Initialize(3, 3, array_buf, policyData);
                BOOST_CHECK_EQUAL(data.num_elements(), 6);

                unsigned int curRow = 0; 
                unsigned int curColumn = 0;
                boost::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn, policyData);
                BOOST_CHECK_EQUAL(1, curRow);
                BOOST_CHECK_EQUAL(0, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn, policyData);
                BOOST_CHECK_EQUAL(2, curRow);
                BOOST_CHECK_EQUAL(0, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn, policyData);
                BOOST_CHECK_EQUAL(0, curRow);
                BOOST_CHECK_EQUAL(1, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn, policyData);
                BOOST_CHECK_EQUAL(1, curRow);
                BOOST_CHECK_EQUAL(1, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn, policyData);
                BOOST_CHECK_EQUAL(2, curRow);
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
                BOOST_CHECK_EQUAL(1, curRow);
                BOOST_CHECK_EQUAL(0, curColumn);

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


