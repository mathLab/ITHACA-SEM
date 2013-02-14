///////////////////////////////////////////////////////////////////////////////
//
// File: testBoostUtil.cpp
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
// Description: Tests the boost utility functions.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/BoostUtil.hpp>

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>
#include <LibUtilities/BasicUtils/ConsistentObjectAccess.hpp>

namespace Nektar
{
    namespace UnitTests
    {
        class TestClass
        {
            public:
                TestClass()
                {
                    ++constructionCount;
                }

                ~TestClass()
                {
                    ++destructionCount;
                }

                static unsigned int constructionCount;
                static unsigned int destructionCount;
        };

        unsigned int TestClass::constructionCount = 0;
        unsigned int TestClass::destructionCount = 0;

        template<typename T>
        class FakeClass {};
        
        
        BOOST_AUTO_TEST_CASE(testMakePtr)
        {
            boost::shared_ptr<FakeClass<int> > a(new FakeClass<int>());
            //FakeClass<int>& b = ConsistentObjectAccess<boost::shared_ptr<FakeClass<int> > >::reference(a);
            
            {
                boost::shared_ptr<TestClass> p = MakePtr(new TestClass());
                BOOST_CHECK(TestClass::constructionCount == 1);
                BOOST_CHECK(TestClass::destructionCount == 0);
            }

            BOOST_CHECK(TestClass::constructionCount == 1);
            BOOST_CHECK(TestClass::destructionCount == 1);

        }
    }
}


