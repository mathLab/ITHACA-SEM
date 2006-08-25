/////////////////////////////////////////////////////////////////////////////////
////
//// File: testNekManager.cpp
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
//// Description: Test code for NekManager
////
/////////////////////////////////////////////////////////////////////////////////
//
//#include <LibUtilities/BasicUtils/NekManager.hpp>
//
//#include <boost/test/auto_unit_test.hpp>
//#include <boost/test/test_case_template.hpp>
//#include <boost/test/floating_point_comparison.hpp>
//#include <boost/test/unit_test.hpp>
//#include <boost/progress.hpp>
//
//#include <iostream>
//
//using std::cout;
//using std::endl;
//
//// Test stuff.
//double MyCreator(const int &key)
//{
//    return key + 2.5;
//};
//
//class Temp
//{
//public:
//    Temp() {};
//
//    static double create(const int &key)
//    {
//        return key*1.025;
//    }
//};
//
//namespace Nektar
//{
//    namespace UnitTests
//    {
//        void testNekManager()
//        {
//			typedef NekManager<int, double> ManagerType;
//
//			ManagerType manager;
//			int key = 10;
//
//			// Registering a C function
//			manager.RegisterCreator(key, (ManagerType::CreateFunc)MyCreator);
//
//			// Registering a static class method
//			manager.RegisterCreator(20, (ManagerType::CreateFunc)Temp::create);
//
//			double value = manager.GetValue(key);
//			BOOST_CHECK(value == 12.5);
//
//			value = manager.GetValue(20);
//			BOOST_CHECK(value == 20.5);
//
//			return;
//		}
//    }
//}
//
///**
//    $Log: testNekManager.cpp,v $
//    Revision 1.1  2006/07/05 20:22:05  jfrazier
//    Initial checkin.
//
//**/
