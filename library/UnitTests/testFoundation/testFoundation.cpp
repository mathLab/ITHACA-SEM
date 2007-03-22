///////////////////////////////////////////////////////////////////////////////
//
// File: testFoundation.cpp
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
// Description: Test code for Foundation
// Author: Sophia Han
//
///////////////////////////////////////////////////////////////////////////////
    
    #include <boost/test/unit_test.hpp>
    #include <boost/test/test_tools.hpp>
    #include <boost/test/floating_point_comparison.hpp>
    
    
    #include <iostream>
    #include <limits>
    using namespace std;
    
	#include <UnitTests/testFoundation/testFoundation.h>
    #include <LibUtilities/BasicUtils/NekManager.hpp>
    #include <LibUtilities/Foundations/Points.h>
    #include <LibUtilities/Foundations/GaussPoints.h>
    #include <LibUtilities/Foundations/PolyEPoints.h>
    #include <LibUtilities/Foundations/Basis.h>
    #include <LibUtilities/Foundations/NodalTriFekete.h>
    #include <LibUtilities/Foundations/ManagerAccess.h>
    
    using namespace Nektar;
    using namespace Nektar::LibUtilities;
    
    namespace Nektar{
    
        namespace foundationUnitTests{
            long double polyFunc(long double x){
                return ((1.0/20.0*x*x - 1.0/6.0)*x*x + 3.0/4.0)*x + 1.0;
            }
            long double polyFunc2(long double x){
               return  (((3.0*x*x - 5.0)*x + 1.0)*x - 2.0)*x + 3.0;
            }

            long double polyFunc3(long double x){
               return  (((33.0*x*x + 27.0)*x*x - 7.0)*x*x + 5.0)*x*x*x*x  + 3.0;
            }

            
            void testGaussGaussLegendre(){
            
                PointsType type = eGaussGaussLegendre;
    
                long double exact = 20.0/3.0;
                //      long double exact = 2.0;
                //      long double exact = 18.0;
    
                for(int nPts = 4; nPts<=20; ++nPts){
                    //long double epsilon = numeric_limits<double>::epsilon() + numeric_limits<long double>::epsilon()/2.0 * nPts;
                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 2);
                    
                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
                    const double *z, *w;
                    points->GetZW(z, w);
    
                    int numPoints = points->GetNumPoints();
                    long double numericIntegral = 0.0;
        
                    for(int j = 0; j < numPoints; ++j) {
                        numericIntegral += w[j] * polyFunc2(z[j]);
                    //   cout << "w["<<j<<"] = " << w[j] << ", summand["<<j<<"] = " << w[j] * polyFunc(z[j]) << endl;
                    }
                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
                    long double relativeError = (exact - numericIntegral)/exact;
                    cout << "nPts: " << nPts << ",     relativeError = " << relativeError
                        << ",      relError / nPts = " << relativeError/(nPts + 2) << ",       epsilon = " << epsilon << endl;
                }
    
    //             cout << "The difference between 1 and the smallest "
    //                 << "value greater than 1\n for float objects is: "
    //                 << numeric_limits<float>::epsilon()
    //                 << endl;
    //             cout << "The difference between 1 and the smallest "
    //                 << "value greater than 1\n for double objects is: "
    //                 << numeric_limits<double>::epsilon()
    //                 << endl;
    //             cout << "The difference between 1 and the smallest "
    //                 << "value greater than 1\n for long double objects is: "
    //                 << numeric_limits<long double>::epsilon()
    //                 << endl;
    // 
    //             cout << "Half of epsilon is: " << numeric_limits<double>::epsilon() / 2.0 << endl;
    
            cout<<"End of Test GaussGaussLegendre()" << endl;
            cout<<"" << endl;
            }
    
            
            
            
            void testGaussRadauMLegendre(){
                
                PointsType type = eGaussRadauMLegendre;
                  long double exact = 20.0/3.0;
                  //long double exact = 2.0;
    
                for(int nPts = 4; nPts<=20; ++nPts){
                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 2);
                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
                    const double *z, *w;
                    points->GetZW(z, w);
                
                    int numPoints = points->GetNumPoints();
                    long double numericIntegral = 0.0;
        
                    for(int j = 0; j < numPoints; ++j) {
                        numericIntegral += w[j] * polyFunc2(z[j]);
                    }
    
                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
                    long double relativeError = (exact - numericIntegral)/exact;
                    cout << "nPts: " << nPts << ",     relativeError = " << relativeError
                    << ",      relError / nPts = " << relativeError/(nPts + 2) << ",       epsilon = " << epsilon << endl;
                }
                cout<<"End of Test GaussRadauMLegendre()" << endl;
                cout<<"" << endl;
            }
    
            void testGaussRadauPLegendre(){
                  
                PointsType type = eGaussRadauPLegendre;
                // long double exact = 20.0/3.0;
                   long double exact = 2.0;
                for(int nPts = 4; nPts<=20; ++nPts){
                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 2);
                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
                    const double *z, *w;
                    points->GetZW(z, w);
                
                    int numPoints = points->GetNumPoints();
                    long double numericIntegral = 0.0;
        
                    for(int j = 0; j < numPoints; ++j) {
                        numericIntegral += w[j] * polyFunc(z[j]);
                    }
   
                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
                    long double relativeError = (exact - numericIntegral)/exact;
                    cout << "nPts: " << nPts << ",     relativeError = " << relativeError
                    << ",      relError / nPts = " << relativeError/(nPts + 2) << ",       epsilon = " << epsilon << endl;
                }
                cout<<"End of Test GaussRadauPLegendre()" << endl;
                cout<<"" << endl;
            }
    
            void testGaussLobattoLegendre(){
            
                PointsType type = eGaussLobattoLegendre;
                   // long double exact = 20.0/3.0;
                   long double exact = 2.0;
    
                for(int nPts = 4; nPts<=20; ++nPts){
                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 2);
                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
                    const double *z, *w;
                    points->GetZW(z, w);
                
                    int numPoints = points->GetNumPoints();
                    long double numericIntegral = 0.0;
        
                    for(int j = 0; j < numPoints; ++j) {
                        numericIntegral += w[j] * polyFunc(z[j]);
                    }
    
                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
                    long double relativeError = (exact - numericIntegral)/exact;
                    cout << "nPts: " << nPts << ",     relativeError = " << relativeError
                    << ",      relError / nPts = " << relativeError/(nPts + 2) << ",       epsilon = " << epsilon << endl;    
                }
                cout<<"End of Test eGaussLobattoLegendre()" << endl;
                cout<<"" << endl;
            }

               void testGaussGaussChebyshev(){
            
                PointsType type = eGaussGaussChebyshev;
                 long double exact = 7.0/2.0*M_PI;
               // long double exact = M_PI;
    
                for(int nPts = 4; nPts<=20; ++nPts){
                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 2);
                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
                    const double *z, *w;
                    points->GetZW(z, w);
                
                    int numPoints = points->GetNumPoints();
                    long double numericIntegral = 0.0;
        
                    for(int j = 0; j < numPoints; ++j) {
                        numericIntegral += w[j] * polyFunc2(z[j]);//*sqrt(1.0 - z[j]*z[j]);
                    }
    
                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
                    long double relativeError = (exact - numericIntegral)/exact;
                    cout << "nPts: " << nPts << ",     relativeError = " << relativeError
                    << ",      relError / nPts = " << relativeError/(nPts + 2) << ",       epsilon = " << epsilon << endl;    
                }
                cout<<"End of Test eGaussGaussChebyshev()" << endl;
                cout<<"" << endl;
            }

            void testGaussRadauMChebyshev(){
            
                PointsType type = eGaussRadauMChebyshev;
                 long double exact = 7.0/2.0*M_PI;
               //    long double exact = M_PI;
    
                for(int nPts = 4; nPts<=20; ++nPts){
                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 2);
                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
                    const double *z, *w;
                    points->GetZW(z, w);
                
                    int numPoints = points->GetNumPoints();
                    long double numericIntegral = 0.0;
        
                    for(int j = 0; j < numPoints; ++j) {
                        numericIntegral += w[j] * polyFunc2(z[j]);//*sqrt(1.0 - z[j]*z[j]);
                    }
    
                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
                    long double relativeError = (exact - numericIntegral)/exact;
                    cout << "nPts: " << nPts << ",     relativeError = " << relativeError
                    << ",      relError / nPts = " << relativeError/(nPts + 2) << ",       epsilon = " << epsilon << endl;    
                }
                cout<<"End of Test eGaussRadauMChebyshev()" << endl;
                cout<<"" << endl;
            }

            void testGaussRadauPChebyshev(){
            
                PointsType type = eGaussRadauPChebyshev;

                    long double exact = M_PI;
                    // long double exact = 4657.0/256.0*M_PI;
                    // long double exact = 7.0/2.0*M_PI;
    
                for(int nPts = 4; nPts<=20; ++nPts){
                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 2)*(nPts + 2);
                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
                    const double *z, *w;
                    points->GetZW(z, w);
                
                    int numPoints = points->GetNumPoints();
                    long double numericIntegral = 0.0;
        
                    for(int j = 0; j < numPoints; ++j) {
                        numericIntegral += w[j] * polyFunc(z[j]);//*sqrt(1.0 - z[j]*z[j]);
                    }
    
                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
                    long double relativeError = (exact - numericIntegral)/exact;
                    cout << "nPts: " << nPts << ",     relativeError = " << relativeError
                    << ",      relError / nPts = " << relativeError/((nPts + 2)*(nPts + 2)) << ",       epsilon = " << epsilon << endl;
                }
                cout<<"End of Test eGaussRadauPChebyshev()" << endl;
                cout<<"" << endl;
            }

             void testGaussLobattoChebyshev(){
            
                PointsType type = eGaussLobattoChebyshev;
                 long double exact = 7.0/2.0*M_PI;
                 //    long double exact = M_PI;
                 //    long double exact = 2717.0/256.0*M_PI;
    
                for(int nPts = 4; nPts<=20; ++nPts){
                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 2)*(nPts + 2);
                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
                    const double *z, *w;
                    points->GetZW(z, w);
                
                    int numPoints = points->GetNumPoints();
                    long double numericIntegral = 0.0;
        
                    for(int j = 0; j < numPoints; ++j) {
                        numericIntegral += w[j] * polyFunc2(z[j]);//*sqrt(1.0 - z[j]*z[j]);
                    }
    
                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
                    long double relativeError = (exact - numericIntegral)/exact;
                    cout << "nPts: " << nPts << ",     relativeError = " << relativeError
                    << ",      relError / nPts = " << relativeError/((nPts + 2)*(nPts + 2)) << ",       epsilon = " << epsilon << endl;
                }
                cout<<"End of Test eGaussLobattoChebyshev()" << endl;
                cout<<"" << endl;
            }

             void testGaussRadauMAlpha0Beta1(){
            
                PointsType type = eGaussRadauMAlpha0Beta1;
                   long double exact = 88.0/21.0;
    
                for(int nPts = 4; nPts<=20; ++nPts){
                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 4);
                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
                    const double *z, *w;
                    points->GetZW(z, w);
                
                    int numPoints = points->GetNumPoints();
                    long double numericIntegral = 0.0;
        
                    for(int j = 0; j < numPoints; ++j) {
                        numericIntegral += w[j] * polyFunc2(z[j]);
                    }
    
                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
                    long double relativeError = (exact - numericIntegral)/exact;
                    cout << "nPts: " << nPts << ",     relativeError = " << relativeError
                    << ",      relError / nPts = " << relativeError/(nPts + 4) << ",       epsilon = " << epsilon << endl;
                }
                cout<<"End of Test eGaussRadauMAlpha0Beta1()" << endl;
                cout<<"" << endl;
            }                                          
            void testGaussRadauMAlpha0Beta2(){
            
                PointsType type = eGaussRadauMAlpha0Beta2;
                   long double exact = 144.0/35.0;
    
                for(int nPts = 4; nPts<=20; ++nPts){
                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (2*(nPts + 2));
                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
                    const double *z, *w;
                    points->GetZW(z, w);
                
                    int numPoints = points->GetNumPoints();
                    long double numericIntegral = 0.0;
        
                    for(int j = 0; j < numPoints; ++j) {
                        numericIntegral += w[j] * polyFunc2(z[j]);
                    }
    
                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
                    long double relativeError = (exact - numericIntegral)/exact;
                    cout << "nPts: " << nPts << ",     relativeError = " << relativeError
                    << ",      relError / nPts = " << relativeError/(2*(nPts + 2)) << ",       epsilon = " << epsilon << endl;    
                }
                cout<<"End of Test eGaussRadauMAlpha0Beta2()" << endl;
                cout<<"" << endl;
            }
            
            void testPolyEvenlySpaced(){
            
                PointsType type = ePolyEvenlySpaced;
                   long double exact = 20.0/3.0;
                   //  long double exact = 2.0;
    
                for(int nPts = 3; nPts<=15; ++nPts){
                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 2);
                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
                    const double *z, *w;
                    points->GetZW(z, w);
                
                    int numPoints = points->GetNumPoints();
                    long double numericIntegral = 0.0;
        
                    for(int j = 0; j < numPoints; ++j) {
                        numericIntegral += w[j] * polyFunc2(z[j]);
                    }
    
                    BOOST_CHECK_CLOSE( numericIntegral, exact, 100.0*epsilon);
                    long double relativeError = (exact - numericIntegral)/exact;
                    cout << "nPts: " << nPts << ",     relativeError = " << relativeError
                    << ",      relError / nPts = " << relativeError/(nPts + 2) << ",       epsilon = " << epsilon << endl;    
                }
                cout<<"End of Test ePolyEvenlySpaced()" << endl;
                cout<<"" << endl;
            }

             long double fourierFunc(int N, long double x){
               return cos(N/2*x) + sin((N/2 - 2)*x);
            }
            void testFourierEvenlySpaced(){
            
                PointsType type = eFourierEvenlySpaced;
               long double exact[] = {0.0, 0.0, 2*sin(1), 0.0, sin(2), 0.0, (2*sin(3))/3, 0.0, sin(4)/2, 0.0, (2*sin(5))/5};
               
                for(int nPts = 2; nPts<=10; nPts += 2){
                    long double epsilon = numeric_limits<double>::epsilon()/2.0 * (nPts + 2);
                    const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
                    const double *z, *w;
                    points->GetZW(z, w);
                
                    int numPoints = points->GetNumPoints();
                    long double numericIntegral = 0.0;

                      for(int j = 0; j < numPoints; ++j) {
                         numericIntegral += w[j] * fourierFunc(nPts, z[j]);

                        cout << "w["<<j<<"] = " << w[j] << ", z["<<j<<"] = " << z[j] << endl;
                    }
    
                    BOOST_CHECK_CLOSE( numericIntegral, exact[nPts], 100.0*epsilon);
                    long double relativeError = (exact[nPts] - numericIntegral)/exact[nPts];
                    cout << "nPts: " << nPts << ",     relativeError = " << relativeError
                    << ",      relError / nPts = " << relativeError/(nPts + 2) << ",       epsilon = " << epsilon << endl;    
                }
                cout<<"End of Test eFourierEvenlySpaced()" << endl;
                cout<<"" << endl;
            }                                              
        }
    }
    
    
