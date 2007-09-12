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
//#include <boost/shared_ptr.hpp>
//
//#include <iostream>
//
//using std::cout;
//using std::endl;
//
//
//
//namespace Nektar
//{
//   namespace UnitTests
//   {
//        class DoubleWrapper
//        {
//            public:
//                DoubleWrapper(double v) : val(v) {}
//                double val;
//        };
//
//        
//        // Test stuff.
//        boost::shared_ptr<DoubleWrapper> MyCreator(const int &key)
//        {
//            return boost::shared_ptr<DoubleWrapper>(new DoubleWrapper(key + 2.5));;
//        };
//        
//        class Temp
//        {
//            public:
//                Temp() {};
//                
//                static boost::shared_ptr<DoubleWrapper> create(const int &key)
//                {
//                    return boost::shared_ptr<DoubleWrapper>(new DoubleWrapper(key*1.025));
//                }
//        };
//        
//        class GlobalCreator
//        {
//            public:
//                boost::shared_ptr<DoubleWrapper> operator()(const int& key)
//                {
//                    return boost::shared_ptr<DoubleWrapper>(new DoubleWrapper(key));
//                }
//        };
//
//       void testNekManager()
//       {
//           typedef LibUtilities::NekManager<int, DoubleWrapper> ManagerType;
//
//            /// TODO - See why I can't do GlobalCreator() directly in the constructor call.
//            GlobalCreator c;
//            ManagerType manager(c);
//            int key = 10;
//
//            // Registering a C function
//            manager.RegisterCreator(key, MyCreator);
//
//            // Registering a static class method
//            manager.RegisterCreator(20, Temp::create);
//            
//            boost::shared_ptr<DoubleWrapper> value = manager[key];
//            BOOST_CHECK(value->val == 12.5);
//
//            value = manager[20];
//            BOOST_CHECK(value->val == 20.5);
//            
//            manager[17] = boost::shared_ptr<DoubleWrapper>(new DoubleWrapper(-2.0));
//            BOOST_CHECK_EQUAL(manager[17]->val, -2.0);
//        }
//        
//
////         typedef boost::shared_ptr<NekMatrix<double, FullMatrixTag, eBlock> > MatrixType;
////         
////         MatrixType create(int k)
////         {
////             return MatrixType(new NekMatrix<double, DiagonalMatrixTag>(k, k, 2, 2));
////         }
////         
////         void testNekMatrixManager()
////         {
////             NekManager<int, NekMatrix<double, DiagonalMatrixTag> > manager(create);
////             
////             manager[10] = boost::shared_ptr<NekMatrix<double, DiagonalMatrixTag> >(new NekMatrix<double, DiagonalMatrixTag>(3, 3));
////            
////            if( manager[11] )
////            {
////            }
////         }
//   }
//}
//
///**
//   $Log: testNekManager.cpp,v $
//   Revision 1.9  2007/07/22 23:04:28  bnelson
//   Backed out Nektar::ptr.
//
//   Revision 1.8  2007/07/20 02:24:41  bnelson
//   Replaced boost::shared_ptr with Nektar::ptr
//
//   Revision 1.7  2007/06/10 23:45:59  bnelson
//   Matrix updates.
//
//   Revision 1.6  2007/01/18 18:44:45  bnelson
//   Updates to compile on Visual Studio 2005.
//
//   Revision 1.5  2007/01/16 17:17:58  bnelson
//   Updated to use shared pointers in the NekManager.
//
//   Revision 1.4  2006/12/17 21:45:25  bnelson
//   Added tests for the global create function.
//
//   Revision 1.3  2006/12/10 21:36:08  bnelson
//   *** empty log message ***
//
//   Revision 1.2  2006/08/25 01:38:58  bnelson
//   no message
//
//   Revision 1.1  2006/07/05 20:22:05  jfrazier
//   Initial checkin.
//
//**/
