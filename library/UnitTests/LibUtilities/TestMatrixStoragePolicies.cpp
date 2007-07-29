///////////////////////////////////////////////////////////////////////////////
//
// File: TestMatrixStoragePolicies.cpp
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

#include <UnitTests/LibUtilities/TestMatrixStoragePolicies.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

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

        void Test0ParameterInitialize()
        {
            Array<OneD, NekDouble> result = Policy::Initialize();
            BOOST_CHECK_EQUAL(result.num_elements(), 0);
            BOOST_CHECK_EQUAL(result.GetOffset(), 0);
            BOOST_CHECK(result.data() != 0);
        }

        void Test2ParameterInitialize()
        {
            {
                Array<OneD, NekDouble> result = Policy::Initialize(5, 5);
                BOOST_CHECK_EQUAL(result.num_elements(), 15);
            }

            {
                BOOST_CHECK_THROW( (Policy::Initialize(5,6)), ErrorUtil::NekError);
                BOOST_CHECK_THROW( (Policy::Initialize(6,5)), ErrorUtil::NekError);
            }
        }

        void TestSingleValuePopulationInitialize()
        {
            {
                Array<OneD, NekDouble> result = Policy::Initialize(5, 5, 7.2);
                BOOST_CHECK_EQUAL(result.num_elements(), 15);
                for(Array<OneD, NekDouble>::iterator iter = result.begin(); iter != result.end(); ++iter)
                {
                    BOOST_CHECK_EQUAL(*iter, 7.2);
                }
            }

            {
                BOOST_CHECK_THROW( (Policy::Initialize(5,6,7.2)), ErrorUtil::NekError);
                BOOST_CHECK_THROW( (Policy::Initialize(6,5,7.2)), ErrorUtil::NekError);
            }
        }

        void TestCArrayInitialization()
        {
            NekDouble buf[] = {1.0, 2.0, 3.0,
                                    5.0, 6.0,
                                         9.0};

            {
                Array<OneD, NekDouble> result = Policy::Initialize(3, 3, buf);
                BOOST_CHECK_EQUAL(result.num_elements(), 6);

                BOOST_CHECK_EQUAL(result[0], 1.0);
                BOOST_CHECK_EQUAL(result[1], 2.0);
                BOOST_CHECK_EQUAL(result[2], 3.0);
                BOOST_CHECK_EQUAL(result[3], 5.0);
                BOOST_CHECK_EQUAL(result[4], 6.0);
                BOOST_CHECK_EQUAL(result[5], 9.0);
            }

            {

                BOOST_CHECK_THROW( (Policy::Initialize(5,6,buf)), ErrorUtil::NekError);
                BOOST_CHECK_THROW( (Policy::Initialize(6,5,buf)), ErrorUtil::NekError);
            }
        }

        void TestArrayInitialization()
        {
            NekDouble buf[] = {1.0, 2.0, 3.0,
                                    5.0, 6.0,
                                         9.0};
            Array<OneD, NekDouble> array_buf(6, buf);

            {
                Array<OneD, NekDouble> result = Policy::Initialize(3, 3, array_buf);
                BOOST_CHECK_EQUAL(result.num_elements(), 6);

                BOOST_CHECK_EQUAL(result[0], 1.0);
                BOOST_CHECK_EQUAL(result[1], 2.0);
                BOOST_CHECK_EQUAL(result[2], 3.0);
                BOOST_CHECK_EQUAL(result[3], 5.0);
                BOOST_CHECK_EQUAL(result[4], 6.0);
                BOOST_CHECK_EQUAL(result[5], 9.0);
            }

            {

                BOOST_CHECK_THROW( (Policy::Initialize(5,6,buf)), ErrorUtil::NekError);
                BOOST_CHECK_THROW( (Policy::Initialize(6,5,buf)), ErrorUtil::NekError);
            }

            {
                NekDouble buf[] = {1.0, 2.0, 3.0,
                                    5.0, 6.0,
                                         9.0, 10.0};
                Array<OneD, NekDouble> array_buf(7, buf);
                BOOST_CHECK_NO_THROW(Policy::Initialize(3, 3, array_buf));
            }

            {
                NekDouble buf[] = {1.0, 2.0, 3.0,
                                    5.0, 6.0,
                                         9.0, 10.0};
                Array<OneD, NekDouble> array_buf(5, buf);
                BOOST_CHECK_THROW(Policy::Initialize(3, 3, array_buf), ErrorUtil::NekError);
            }
        }



        void TestConstGetValue()
        {
            NekDouble buf[] = {1.0, 2.0, 3.0,
                                    5.0, 6.0,
                                         9.0};
            Array<OneD, NekDouble> array_buf(6, buf);
            Array<OneD, NekDouble> data = Policy::Initialize(3, 3, array_buf);
            BOOST_CHECK_EQUAL(data.num_elements(), 6);

            {
                BOOST_CHECK_EQUAL(1.0, Policy::GetValue(3, 3, 0, 0, data));
                BOOST_CHECK_EQUAL(2.0, Policy::GetValue(3, 3, 0, 1, data));
                BOOST_CHECK_EQUAL(3.0, Policy::GetValue(3, 3, 0, 2, data));
                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 1, 0, data));
                BOOST_CHECK_EQUAL(5.0, Policy::GetValue(3, 3, 1, 1, data));
                BOOST_CHECK_EQUAL(6.0, Policy::GetValue(3, 3, 1, 2, data));
                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 2, 0, data));
                BOOST_CHECK_EQUAL(0.0, Policy::GetValue(3, 3, 2, 1, data));
                BOOST_CHECK_EQUAL(9.0, Policy::GetValue(3, 3, 2, 2, data));
            }

            #if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)
            {
                BOOST_CHECK_THROW(Policy::GetValue(3,4,3,3,data), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::GetValue(3,3,4,3,data), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::GetValue(3,3,3,4,data), ErrorUtil::NekError);
            }
            #endif
        }

        void TestSetValue()
        {
            NekDouble buf[] = {1.0, 2.0, 3.0,
                                    5.0, 6.0,
                                         9.0};
            Array<OneD, NekDouble> array_buf(6, buf);
            Array<OneD, NekDouble> data = Policy::Initialize(3, 3, array_buf);
            BOOST_CHECK_EQUAL(data.num_elements(), 6);
            
            Policy::SetValue(3,3,0,0,data, 7.2);
            BOOST_CHECK_EQUAL(7.2, Policy::GetValue(3,3,0,0,data));
            BOOST_CHECK_EQUAL(1.0, array_buf[0]);
            BOOST_CHECK_EQUAL(7.2, data[0]);

            Policy::SetValue(3,3,0,1,data,10.0);
            BOOST_CHECK_EQUAL(10.0, Policy::GetValue(3,3,0,1,data));
            Policy::SetValue(3,3,0,2,data,20.0);
            BOOST_CHECK_EQUAL(20.0, Policy::GetValue(3,3,0,2,data));
            Policy::SetValue(3,3,1,1,data,30.0);
            BOOST_CHECK_EQUAL(30.0, Policy::GetValue(3,3,1,1,data));
            Policy::SetValue(3,3,1,2,data,40.0);
            BOOST_CHECK_EQUAL(40.0, Policy::GetValue(3,3,1,2,data));
            Policy::SetValue(3,3,2,2,data,50.0);
            BOOST_CHECK_EQUAL(50.0, Policy::GetValue(3,3,2,2,data));

            // Can't set in the lower triangle.
            BOOST_CHECK_THROW(Policy::SetValue(3,3,1,0,data,8.0), ErrorUtil::NekError);
            BOOST_CHECK_THROW(Policy::SetValue(3,3,2,0,data,8.0), ErrorUtil::NekError);
            BOOST_CHECK_THROW(Policy::SetValue(3,3,2,1,data,8.0), ErrorUtil::NekError);

            #if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)
            {
                BOOST_CHECK_THROW(Policy::SetValue(3,4,3,3,data,8.0), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::SetValue(3,3,4,3,data,8.0), ErrorUtil::NekError);
                BOOST_CHECK_THROW(Policy::SetValue(3,3,3,4,data,8.0), ErrorUtil::NekError);
            }
            #endif
        }
        void TestAdvance()
        {
            {
                NekDouble buf[] = {1.0, 2.0, 3.0,
                                        5.0, 6.0,
                                             9.0};
                Array<OneD, NekDouble> array_buf(6, buf);
                Array<OneD, NekDouble> data = Policy::Initialize(3, 3, array_buf);
                BOOST_CHECK_EQUAL(data.num_elements(), 6);

                unsigned int curRow = 0; 
                unsigned int curColumn = 0;
                boost::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn);
                BOOST_CHECK_EQUAL(0, curRow);
                BOOST_CHECK_EQUAL(1, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn);
                BOOST_CHECK_EQUAL(0, curRow);
                BOOST_CHECK_EQUAL(2, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn);
                BOOST_CHECK_EQUAL(1, curRow);
                BOOST_CHECK_EQUAL(1, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn);
                BOOST_CHECK_EQUAL(1, curRow);
                BOOST_CHECK_EQUAL(2, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn);
                BOOST_CHECK_EQUAL(2, curRow);
                BOOST_CHECK_EQUAL(2, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(3, 3, curRow, curColumn);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curRow);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curColumn);
            }

            {
                NekDouble buf[] = {1.0};
                Array<OneD, NekDouble> array_buf(1, buf);
                Array<OneD, NekDouble> data = Policy::Initialize(1, 1, array_buf);
                BOOST_CHECK_EQUAL(data.num_elements(), 1);

                unsigned int curRow = 0; 
                unsigned int curColumn = 0;
                boost::tie(curRow, curColumn) = Policy::Advance(1, 1, curRow, curColumn);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curRow);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curColumn);
            }

            {
                NekDouble buf[] = {1.0, 2.0,
                                        5.0};
                Array<OneD, NekDouble> array_buf(3, buf);
                Array<OneD, NekDouble> data = Policy::Initialize(2, 2, array_buf);
                BOOST_CHECK_EQUAL(data.num_elements(), 3);

                unsigned int curRow = 0; 
                unsigned int curColumn = 0;
                boost::tie(curRow, curColumn) = Policy::Advance(2, 2, curRow, curColumn);
                BOOST_CHECK_EQUAL(0, curRow);
                BOOST_CHECK_EQUAL(1, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(2, 2, curRow, curColumn);
                BOOST_CHECK_EQUAL(1, curRow);
                BOOST_CHECK_EQUAL(1, curColumn);

                boost::tie(curRow, curColumn) = Policy::Advance(2, 2, curRow, curColumn);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curRow);
                BOOST_CHECK_EQUAL(std::numeric_limits<unsigned int>::max(), curColumn);
            }
        }
    }
}


