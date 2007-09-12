/////////////////////////////////////////////////////////////////////////////////
////
//// File: testBoostUtil.cpp
////
//// For more information, please see: http://www.nektar.info
////
//// The MIT License
////
//// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//// Department of Aeronautics, Imperial College London (UK), and Scientific
//// Computing and Imaging Institute, University of Utah (USA).
////
//// License for the specific language governing rights and limitations under
//// Permission is hereby granted, free of charge, to any person obtaining a
//// copy of this software and associated documentation files (the "Software"),
//// to deal in the Software without restriction, including without limitation
//// the rights to use, copy, modify, merge, publish, distribute, sublicense,
//// and/or sell copies of the Software, and to permit persons to whom the
//// Software is furnished to do so, subject to the following conditions:
////
//// The above copyright notice and this permission notice shall be included
//// in all copies or substantial portions of the Software.
////
//// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//// DEALINGS IN THE SOFTWARE.
////
//// Description: Tests the boost utility functions.
////
/////////////////////////////////////////////////////////////////////////////////
//
//#include <UnitTests/testBoostUtil.h>
//#include <LibUtilities/BasicUtils/BoostUtil.hpp>
//
//#include <boost/test/auto_unit_test.hpp>
//#include <boost/test/test_case_template.hpp>
//#include <boost/test/floating_point_comparison.hpp>
//#include <boost/test/unit_test.hpp>
//#include <boost/progress.hpp>
//
//namespace Nektar
//{
//    namespace UnitTests
//    {
//        class TestClass
//        {
//            public:
//                TestClass()
//                {
//                    ++constructionCount;
//                }
//
//                ~TestClass()
//                {
//                    ++destructionCount;
//                }
//
//                static unsigned int constructionCount;
//                static unsigned int destructionCount;
//        };
//
//        unsigned int TestClass::constructionCount = 0;
//        unsigned int TestClass::destructionCount = 0;
//
//        void testMakePtr()
//        {
//            {
//                boost::shared_ptr<TestClass> p = MakePtr(new TestClass());
//                BOOST_CHECK(TestClass::constructionCount == 1);
//                BOOST_CHECK(TestClass::destructionCount == 0);
//            }
//
//            BOOST_CHECK(TestClass::constructionCount == 1);
//            BOOST_CHECK(TestClass::destructionCount == 1);
//
//        }
//    }
//}
//
//
///**
//    $Log: testBoostUtil.cpp,v $
//    Revision 1.6  2007/07/22 23:04:28  bnelson
//    Backed out Nektar::ptr.
//
//    Revision 1.5  2007/07/20 02:24:41  bnelson
//    Replaced boost::shared_ptr with Nektar::ptr
//
//    Revision 1.4  2007/06/13 22:00:48  bnelson
//    *** empty log message ***
//
//    Revision 1.3  2006/07/05 20:21:24  jfrazier
//    Added NekManager test case.
//
//    Revision 1.2  2006/05/07 21:10:34  bnelson
//    Added the log
//
// **/
//
