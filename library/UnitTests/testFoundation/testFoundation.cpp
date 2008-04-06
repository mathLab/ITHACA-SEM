/////////////////////////////////////////////////////////////////////////////////
////
//// File: testFoundation.cpp
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
//// Description: Test Foundation(Integration)
//// Author: Sophia Han
////
/////////////////////////////////////////////////////////////////////////////////
//#undef NEKTAR_MEMORY_POOL_ENABLED    
//#include <boost/test/unit_test.hpp>
//#include <boost/test/test_tools.hpp>
//#include <boost/test/floating_point_comparison.hpp>
//
//#include <iostream>
//#include <limits>
// using namespace std;
//    
//#include <UnitTests/testFoundation/testFoundation.h>
//#include <LibUtilities/BasicUtils/NekManager.hpp>
//#include <LibUtilities/Foundations/Points.h>
//#include <LibUtilities/Foundations/GaussPoints.h>
//#include <LibUtilities/Foundations/PolyEPoints.h>
//#include <LibUtilities/Foundations/Basis.h>
//#include <LibUtilities/Foundations/NodalTriFekete.h>
//#include <LibUtilities/Foundations/ManagerAccess.h>
//
//
//using namespace Nektar;
//using namespace Nektar::LibUtilities;
//
// namespace Nektar{
//    
//        namespace foundationUnitTests{
//            static long double polyFunc(long double x){
//                return ((1.0/20.0*x*x - 1.0/6.0)*x*x + 3.0/4.0)*x + 1.0;
//            }
//            static long double polyFunc2(long double x){
//               return  (((3.0*x*x - 5.0)*x + 1.0)*x - 2.0)*x + 3.0;
//            }
//
//            static long double polyFunc3(long double x){
//               return  (((33.0*x*x + 27.0)*x*x - 7.0)*x*x + 5.0)*x*x*x*x  + 3.0;
//            }
//
//            
//            void testGaussGaussLegendre(){
//
//                PointsType type = eGaussGaussLegendre;
//    
//                long double exact = 20.0/3.0;
//                for(int nPts = 4; nPts<=20; ++nPts){
//                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 2);
//                    
//                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
//
//                    const Array<OneD, const NekDouble> &z = points->GetZ();
//                    const Array<OneD, const NekDouble> &w = points->GetW();
//    
//                    int numPoints = points->GetNumPoints();
//                    long double numericIntegral = 0.0;
//        
//                    for(int j = 0; j < numPoints; ++j) {
//                        numericIntegral += w[j] * polyFunc2(z[j]);
//                    //   cout << "w["<<j<<"] = " << w[j] << ", summand["<<j<<"] = " << w[j] * polyFunc(z[j]) << endl;
//                    }
//                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
//                    long double relativeError = (exact - numericIntegral)/exact;
//          //          cout << "nPts: " << nPts << ",     relativeError = " << relativeError
//          //              << ",      relError / nPts = " << relativeError/(nPts + 2) << ",       epsilon = " << epsilon << endl;
//                }
//        
//         //   cout<<"End of Test GaussGaussLegendre()" << endl;
//         //   cout<<"" << endl;
//            }
//    
//            
//            
//            
//            void testGaussRadauMLegendre(){
//                
//                PointsType type = eGaussRadauMLegendre;
//                  long double exact = 20.0/3.0;
//    
//                for(int nPts = 4; nPts<=20; ++nPts){
//                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 2);
//                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
//
//                    const Array<OneD, const NekDouble> &z = points->GetZ();
//                    const Array<OneD, const NekDouble> &w = points->GetW();
//
//                    int numPoints = points->GetNumPoints();
//                    long double numericIntegral = 0.0;
//        
//                    for(int j = 0; j < numPoints; ++j) {
//                        numericIntegral += w[j] * polyFunc2(z[j]);
//                    }
//    
//                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
//                    long double relativeError = (exact - numericIntegral)/exact;
////                     cout << "nPts: " << nPts << ",     relativeError = " << relativeError
////                     << ",      relError / nPts = " << relativeError/(nPts + 2) << ",       epsilon = " << epsilon << endl;
//                }
////                 cout<<"End of Test GaussRadauMLegendre()" << endl;
////                 cout<<"" << endl;
//            }
//    
//            void testGaussRadauPLegendre(){
//                  
//                PointsType type = eGaussRadauPLegendre;
//                   long double exact = 2.0;
//                for(int nPts = 4; nPts<=20; ++nPts){
//                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 2);
//                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
//
//                    const Array<OneD, const NekDouble> &z = points->GetZ();
//                    const Array<OneD, const NekDouble> &w = points->GetW();
//
//                    int numPoints = points->GetNumPoints();
//                    long double numericIntegral = 0.0;
//        
//                    for(int j = 0; j < numPoints; ++j) {
//                        numericIntegral += w[j] * polyFunc(z[j]);
//                    }
//   
//                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
//                    long double relativeError = (exact - numericIntegral)/exact;
////                     cout << "nPts: " << nPts << ",     relativeError = " << relativeError
////                     << ",      relError / nPts = " << relativeError/(nPts + 2) << ",       epsilon = " << epsilon << endl;
//                }
////                 cout<<"End of Test GaussRadauPLegendre()" << endl;
////                 cout<<"" << endl;
//            }
//    
//            void testGaussLobattoLegendre(){
//            
//                PointsType type = eGaussLobattoLegendre;
//                   long double exact = 2.0;
//    
//                for(int nPts = 4; nPts<=20; ++nPts){
//                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 2);
//                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
//
//                    const Array<OneD, const NekDouble> &z = points->GetZ();
//                    const Array<OneD, const NekDouble> &w = points->GetW();
//
//                    int numPoints = points->GetNumPoints();
//                    long double numericIntegral = 0.0;
//        
//                    for(int j = 0; j < numPoints; ++j) {
//                        numericIntegral += w[j] * polyFunc(z[j]);
//                    }
//    
//                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
//                    long double relativeError = (exact - numericIntegral)/exact;
////                     cout << "nPts: " << nPts << ",     relativeError = " << relativeError
////                     << ",      relError / nPts = " << relativeError/(nPts + 2) << ",       epsilon = " << epsilon << endl;    
//                }
////                 cout<<"End of Test eGaussLobattoLegendre()" << endl;
////                 cout<<"" << endl;
//            }
//
//               void testGaussGaussChebyshev(){
//            
//                PointsType type = eGaussGaussChebyshev;
//                 long double exact = 7.0/2.0*M_PI;
//    
//                for(int nPts = 4; nPts<=20; ++nPts){
//                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 2);
//                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
//
//                    const Array<OneD, const NekDouble> &z = points->GetZ();
//                    const Array<OneD, const NekDouble> &w = points->GetW();
//
//                    int numPoints = points->GetNumPoints();
//                    long double numericIntegral = 0.0;
//        
//                    for(int j = 0; j < numPoints; ++j) {
//                        numericIntegral += w[j] * polyFunc2(z[j]);//*sqrt(1.0 - z[j]*z[j]);
//                    }
//    
//                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
//                    long double relativeError = (exact - numericIntegral)/exact;
//    /*                cout << "nPts: " << nPts << ",     relativeError = " << relativeError
//                    << ",      relError / nPts = " << relativeError/(nPts + 2) << ",       epsilon = " << epsilon << endl;   */ 
//                }
////                 cout<<"End of Test eGaussGaussChebyshev()" << endl;
////                 cout<<"" << endl;
//            }
//
//            void testGaussRadauMChebyshev(){
//            
//                PointsType type = eGaussRadauMChebyshev;
//                 long double exact = 7.0/2.0*M_PI;
//
//                for(int nPts = 4; nPts<=20; ++nPts){
//                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 2);
//                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
//
//                    const Array<OneD, const NekDouble> &z = points->GetZ();
//                    const Array<OneD, const NekDouble> &w = points->GetW();
//                
//                    int numPoints = points->GetNumPoints();
//                    long double numericIntegral = 0.0;
//        
//                    for(int j = 0; j < numPoints; ++j) {
//                        numericIntegral += w[j] * polyFunc2(z[j]);//*sqrt(1.0 - z[j]*z[j]);
//                    }
//    
//                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
//                    long double relativeError = (exact - numericIntegral)/exact;
////                     cout << "nPts: " << nPts << ",     relativeError = " << relativeError
////                     << ",      relError / nPts = " << relativeError/(nPts + 2) << ",       epsilon = " << epsilon << endl;    
//                }
////                 cout<<"End of Test eGaussRadauMChebyshev()" << endl;
////                 cout<<"" << endl;
//            }
//
//            void testGaussRadauPChebyshev(){
//            
//                PointsType type = eGaussRadauPChebyshev;
//                    long double exact = M_PI;
//    
//                for(int nPts = 4; nPts<=20; ++nPts){
//                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 2)*(nPts + 2);
//                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
//
//                    const Array<OneD, const NekDouble> &z = points->GetZ();
//                    const Array<OneD, const NekDouble> &w = points->GetW();
//
//                    int numPoints = points->GetNumPoints();
//                    long double numericIntegral = 0.0;
//        
//                    for(int j = 0; j < numPoints; ++j) {
//                        numericIntegral += w[j] * polyFunc(z[j]);//*sqrt(1.0 - z[j]*z[j]);
//                    }
//    
//                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
//                    long double relativeError = (exact - numericIntegral)/exact;
////                     cout << "nPts: " << nPts << ",     relativeError = " << relativeError
////                     << ",      relError / nPts = " << relativeError/((nPts + 2)*(nPts + 2)) << ",       epsilon = " << epsilon << endl;
//                }
////                 cout<<"End of Test eGaussRadauPChebyshev()" << endl;
////                 cout<<"" << endl;
//            }
//
//             void testGaussLobattoChebyshev(){
//            
//                PointsType type = eGaussLobattoChebyshev;
//                 long double exact = 7.0/2.0*M_PI;
//    
//                for(int nPts = 4; nPts<=20; ++nPts){
//                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 2)*(nPts + 2);
//                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
//
//                    const Array<OneD, const NekDouble> &z = points->GetZ();
//                    const Array<OneD, const NekDouble> &w = points->GetW();
//                
//                    int numPoints = points->GetNumPoints();
//                    long double numericIntegral = 0.0;
//        
//                    for(int j = 0; j < numPoints; ++j) {
//                        numericIntegral += w[j] * polyFunc2(z[j]);//*sqrt(1.0 - z[j]*z[j]);
//                    }
//    
//                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
//                    long double relativeError = (exact - numericIntegral)/exact;
////                     cout << "nPts: " << nPts << ",     relativeError = " << relativeError
////                     << ",      relError / nPts = " << relativeError/((nPts + 2)*(nPts + 2)) << ",       epsilon = " << epsilon << endl;
//                }
////                 cout<<"End of Test eGaussLobattoChebyshev()" << endl;
////                 cout<<"" << endl;
//            }
//
//             void testGaussRadauMAlpha0Beta1(){
//            
//                PointsType type = eGaussRadauMAlpha0Beta1;
//                   long double exact = 88.0/21.0;
//    
//                for(int nPts = 4; nPts<=20; ++nPts){
//                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 4);
//                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
//
//                    const Array<OneD, const NekDouble> &z = points->GetZ();
//                    const Array<OneD, const NekDouble> &w = points->GetW();
//                     
//                    int numPoints = points->GetNumPoints();
//                    long double numericIntegral = 0.0;
//        
//                    for(int j = 0; j < numPoints; ++j) {
//                        numericIntegral += w[j] * polyFunc2(z[j]);
//                    }
//    
//                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
//                    long double relativeError = (exact - numericIntegral)/exact;
////                     cout << "nPts: " << nPts << ",     relativeError = " << relativeError
////                     << ",      relError / nPts = " << relativeError/(nPts + 4) << ",       epsilon = " << epsilon << endl;
//                }
////                 cout<<"End of Test eGaussRadauMAlpha0Beta1()" << endl;
////                 cout<<"" << endl;
//            }                                          
//            void testGaussRadauMAlpha0Beta2(){
//            
//                PointsType type = eGaussRadauMAlpha0Beta2;
//                   long double exact = 144.0/35.0;
//    
//                for(int nPts = 4; nPts<=20; ++nPts){
//                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (2*(nPts + 2));
//                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
//
//                    const Array<OneD, const NekDouble> &z = points->GetZ();
//                    const Array<OneD, const NekDouble> &w = points->GetW();
//
//                    int numPoints = points->GetNumPoints();
//                    long double numericIntegral = 0.0;
//        
//                    for(int j = 0; j < numPoints; ++j) {
//                        numericIntegral += w[j] * polyFunc2(z[j]);
//                    }
//    
//                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
//                    long double relativeError = (exact - numericIntegral)/exact;
////                     cout << "nPts: " << nPts << ",     relativeError = " << relativeError
////                     << ",      relError / nPts = " << relativeError/(2*(nPts + 2)) << ",       epsilon = " << epsilon << endl;    
//                }
////                 cout<<"End of Test eGaussRadauMAlpha0Beta2()" << endl;
////                 cout<<"" << endl;
//            }
//            
//            void testPolyEvenlySpaced(){
//            
//                PointsType type = ePolyEvenlySpaced;
//                   long double exact = 20.0/3.0;
//    
//                for(int nPts = 3; nPts<=15; ++nPts){
//                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 2);
//                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
//
//                    const Array<OneD, const NekDouble> &z = points->GetZ();
//                    const Array<OneD, const NekDouble> &w = points->GetW();
//
//                    int numPoints = points->GetNumPoints();
//                    long double numericIntegral = 0.0;
//        
//                    for(int j = 0; j < numPoints; ++j) {
//                        numericIntegral += w[j] * polyFunc2(z[j]);
//                    }
//    
//                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
//                    long double relativeError = (exact - numericIntegral)/exact;
//         /*           cout << "nPts: " << nPts << ",     relativeError = " << relativeError
//                    << ",      relError / nPts = " << relativeError/(nPts + 2) << ",       epsilon = " << epsilon << endl; */   
//                }
////                 cout<<"End of Test ePolyEvenlySpaced()" << endl;
////                 cout<<"" << endl;
//            }
//
//            static long double fourierFunc(long double x, int N){
//                long double z = M_PI*(x + 1.0);
//                return (cos(N/2*z) + sin((N/2 - 2)*z))/M_PI;
//            }
//
//            // For fourierFunc, this should evaluate to zero
//            long double TrapezoidalRule(int N) {
//                long double a = -1.0, b = 1.0;
//                long double (*f)(long double, int) = fourierFunc;
//                long double y = (f(a, N) + f(b, N)) / 2.0;
//                for( int k = 1; k <= N-1; ++k ) {
//                    y += f(a + k*(b-a)/N, N);
//                }
//                y *= 2.0 / N;
//                
//                return y;
//            }
//            
//            
//            void testFourierEvenlySpaced(){
//            
//                PointsType type = eFourierEvenlySpaced;
//                long double exact = 0;
//                long double (*f)(long double, int) = fourierFunc;
//               
//                for(int nPts = 2; nPts<=50; nPts += 2){
//                    long double epsilon = nPts* numeric_limits<double>::epsilon();
//                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
//
//                    const Array<OneD, const NekDouble> &z = points->GetZ();
//                    const Array<OneD, const NekDouble> &w = points->GetW();
//                
//                    int numPoints = points->GetNumPoints();
//                    long double numericIntegral = 0.0;
//
//                    for(int j = 0; j < numPoints; ++j) {
//                        numericIntegral += w[j] * f(z[j], numPoints);
//                    }
////                     cout << "Numeric integral = " << numericIntegral << ", exact = " << exact
////                          << ",  trapezoidal rule = " << TrapezoidalRule(numPoints) << endl;
//
//                    long double absoluteError = exact - numericIntegral;
//                    BOOST_CHECK( fabs(absoluteError) < epsilon );
//    /*
//                    cout << "nPts: " << nPts << ",     absoluteError = " << absoluteError
//                         << ",       epsilon = " << epsilon << endl;*/
//                }
////                 cout<<"End of Test eFourierEvenlySpaced()" << endl;
////                 cout<<"" << endl;
//            }                       
//        }
//    }
//    
