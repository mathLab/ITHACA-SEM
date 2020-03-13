/////////////////////////////////////////////////////////////////////////////////
////
//// File: testInterpolation.cpp
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
//// Description: Test Interpolation
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
//using namespace std;
//
//#include <UnitTests/testFoundation/testInterpolation.h>
//#include <LibUtilities/BasicUtils/NekManager.hpp>
//#include <LibUtilities/Foundations/Points.h>
//#include <LibUtilities/Foundations/GaussPoints.h>
//#include <LibUtilities/Foundations/PolyEPoints.h>
//#include <LibUtilities/Foundations/Basis.h>
//#include <LibUtilities/Foundations/NodalTriFekete.h>
//#include <LibUtilities/Foundations/ManagerAccess.h>
//
//using namespace Nektar;
//using namespace Nektar::LibUtilities;
//
//namespace Nektar {
// namespace foundationUnitTests{
//
//static long double polyFunc(long double x) {
//    return  (((3.0*x*x - 5.0)*x + 1.0)*x - 2.0)*x + 3.0;
//}
//static long double polyFunc2(long double x) {
//    return  (((3.0*x*x - 5.0)*x + 1.0)*x - 2.0)*x + 3.0;
//}
//long double polyFunc3(long double x) {
//    return  (((33.0*x*x + 27.0)*x*x - 7.0)*x*x + 5.0)*x*x*x*x  + 3.0;
//}
//
//void testPointsInterpolation(PointsType, long double );
//void testFourierInterpolation(PointsType, long double );
//
//void testInterGaussGaussLegendre(){
//     testPointsInterpolation(eGaussGaussLegendre, 1.1e-12);
//    // cout<<"End of Interpolation Test GaussGaussLegendre()" << endl << endl;
//}
//void testInterGaussRadauMLegendre(){
//     testPointsInterpolation(eGaussRadauMLegendre, 1.1e-11);
//    // cout<<"End of Interpolation Test GaussRadauMLegendre()" <<endl << endl;
//}
//void testInterGaussRadauPLegendre(){
//     testPointsInterpolation(eGaussRadauPLegendre, 1.1e-14);
//   //  cout<<"End of Interpolation Test GaussRadauPLegendre()"<<endl << endl;
//}
// void testInterGaussLobattoLegendre(){
//     testPointsInterpolation(eGaussLobattoLegendre, 1.1e-14);
//    // cout<<"End of Interpolation Test GaussLobattoLegendre()"<<endl << endl;
//}
//void testInterGaussGaussChebyshev(){
//     testPointsInterpolation(eGaussGaussChebyshev, 1.1e-12);
//    // cout<<"End of Interpolation Test GaussGaussChebyshev()"<< endl << endl;
//}
//void testInterGaussRadauMChebyshev(){
//     testPointsInterpolation(eGaussRadauMChebyshev, 1.1e-11);
//    // cout<<"End of Interpolation Test GaussRadauMChebyshev()"<< endl << endl;
//}
//void testInterGaussRadauPChebyshev(){
//     testPointsInterpolation(eGaussRadauPChebyshev, 1.1e-14);
//    // cout<<"End of Interpolation Test GaussRadauPChebyshev()" << endl << endl;
//}
// void testInterGaussLobattoChebyshev(){
//     testPointsInterpolation(eGaussLobattoChebyshev, 1.1e-14);
//    // cout<<"End of Interpolation Test GaussLobattoChebyshev()" << endl << endl;
//}
//void testInterGaussRadauMAlpha0Beta1(){
//     testPointsInterpolation(eGaussRadauMAlpha0Beta1, 1.1e-12);
//    // cout<<"End of Interpolation Test GaussRadauMAlpha0Beta1()"<< endl << endl;
//}
//void testInterGaussRadauMAlpha0Beta2(){
//     testPointsInterpolation(eGaussRadauMAlpha0Beta2, 1.1e-12);
//    // cout<<"End of Interpolation Test GaussRadauMAlpha0Beta2()"<< endl << endl;
//}
//void testInterPolyEvenlySpaced(){
//     testPointsInterpolation(ePolyEvenlySpaced, 1.1e-13);
//    // cout<<"End of Interpolation Test PolyEvenlySpaced()"<< endl << endl;
//}
//
//void testInterFourierEvenlySpaced() {
//    testFourierInterpolation(eFourierEvenlySpaced, 1.1e-14);
//   // cout<<"End of Interpolation Test FourierEvenlySpaced()"<< endl << endl;
//}
//
//long double LagrangeInterpolation(int j, long double x, const Array<OneD, const NekDouble> &z, int N) {
//    long double product = 1.0;
//    for(int i=0; i<N; ++i) {
//        if(i != j) {
//            product *= (x-z[i])/(z[j]-z[i]);
//        }
//    }
//    return product;
//}
//
//// fourier function using Trapezoidal rule
//static long double fourierFunc(long double x, int N) {
//    long double z = M_PI*(x + 1.0);
//    return (cos(N/2.0*z) + sin((N/2.0 - 2.0)*z))/M_PI;
//}
//
//
////long double fourierInterpolationFunction(int j, long double x, SharedArray<const NekDouble> z, int N){
//long double fourierInterpolationFunction(int j, long double x, const Array<OneD, const NekDouble> &z, int N){
//    long double t = (long double)M_PI*(x - z[j])/2;
//    if( fabs(t) < 1e-12 ) {
//        return 1.0;
//    }
//    return sin(t*N) / (N * tan(t));
//}
//
//
//void testFourierInterpolation(PointsType type, long double epsilon) {
////     long double (*f)(long double x, int N) = fourierFunc;
////     const long double eps = epsilon;
////     for(int nPts = 4; nPts<=10; nPts += 2) {
////        const boost::shared_ptr<Points<NekDouble> > points = PointsManager()[PointsKey(nPts, type)];
//// 
////        const Array<OneD, const NekDouble> &z = points->GetZ();
////           
////        // Get the interpolation matrix I
////        int nNodes = nPts * 2;
////        const boost::shared_ptr<Points<NekDouble> > nodes = PointsManager()[PointsKey(nNodes, type)];
////        const Array<OneD, const NekDouble> &zNodePtr = nodes->GetZ();
//// 
////        const boost::shared_ptr<NekMatrix<NekDouble> > Iptr = points->GetI(nNodes, zNodePtr);
////        const NekMatrix<NekDouble> & I = *Iptr;
//// 
////        //int numPoints = points->GetNumPoints();
////        epsilon = eps*nPts;
////        // Check the derivative at each of the i points
////        for(int i = 0; i < I.GetRows(); ++i) {
////            long double exact = f(zNodePtr[i], nPts);
////            long double numericInterpolation = 0.0;
//// 
////            // Multiply the derivative matrix with the sample point vector to get the derivative
////             for(int j = 0; j < I.GetColumns(); ++j) {
////                numericInterpolation += I(i,j) * f(z[j], nPts);
////            }
//// 
////            // Compute the relative error; deals appropiately with the case when x' = 0
////            long double relativeError = (exact - numericInterpolation)/exact;
////            if( fabs(exact) < numeric_limits<double>::epsilon() ) {
////                relativeError = exact - numericInterpolation;
////                BOOST_CHECK( fabs(relativeError) < epsilon );
////            } else {
////                BOOST_CHECK_CLOSE(exact, numericInterpolation, 100.0*epsilon);
////            }
//// 
////            // Output the diagnostics if the verbose flag is set
////            if( bool isVerbose = true ) {
////                cout << "nPts: " << nPts << ",     relativeError = " << relativeError
////                //<< ",      relError / nPts = " << relativeError/(nPts + 2) << ",       epsilon = " << epsilon
////                <<",       exact = " << exact <<",       numerical = " << numericInterpolation<< endl; //}
//// 
//// 
////                long double expectedInterpolation = 0.0;
////                for(int j = 0; j<I.GetColumns(); ++j) {
////                    long double F = fourierInterpolationFunction(j, zNodePtr[i], z, nPts);
////                    cout << "h[" << j << "](x[" << i << "]) = " << F << endl;
////                    expectedInterpolation += F * f(z[j], nPts);
////                }
////                cout << "expectedInterpolation = " << expectedInterpolation << endl;
////                    cout << "I = " << endl;
////                    for(int k=0; k<I.GetRows(); ++k) {
////                        for(int j =0; j<I.GetColumns(); ++j) {
////                            printf("% 5.3f  ", I(k, j));
////                        }
////                        cout << endl;
////                    }
//// 
////                long double L1_interpolationMatrixError = 0;
////                double Linf = 0;
////                 long double RMS = 0, mean = 0;
//// 
////                cout << "Expected I = " << endl;
////                for(int k=0; k<I.GetRows(); ++k) {
////                    for(int j=0; j<I.GetColumns(); ++j) {
////                        long double F = fourierInterpolationFunction(j, zNodePtr[k], z, nPts);
////                        long double error = I(k,j) - F;
////                        Linf = max(Linf, static_cast<double>(fabs(error)));
////                        L1_interpolationMatrixError += fabs(error);
////                        RMS += error*error;
////                        mean += error;
////                        
////                        printf("% 5.3Lf  ", F);
////                    }
////                    cout << endl;
////                }
////                int N = I.GetRows() * I.GetColumns();
////                RMS = sqrt(RMS)/N;
////                mean /= N;
////                long double variance = fabs(RMS*RMS - mean*mean);
////                cout << "Interpolation Matrix Errors:" << endl;
////                cout << "L1       = " << L1_interpolationMatrixError << endl;
////                cout << "Linf     = " << Linf << endl;
////                cout << "RMS      = " << RMS << endl;
////                cout << "mean     = " << mean << endl;
////                cout << "stdDev   = " << sqrt(variance) << endl;
//// 
////           //     cout << "//////////////////////////////////////////////////////" << endl;
////            }//
////        }
////     }
//}//
//
//void testPointsInterpolation(PointsType type, long double epsilon) {
////     long double (*f)(long double x) = polyFunc;
////     const long double eps = epsilon;
////      for(int nPts = 8; nPts<=15; ++nPts) {
////        const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
////        const Array<OneD, const NekDouble> & z = points->GetZ();
////        
////        
////        // Get the interpolation matrix I
////        int nNodes = nPts * 2;
////        const boost::shared_ptr<Points<double> > nodes = PointsManager()[PointsKey(nNodes, type)];
////         const Array<OneD, const NekDouble> & zNodePtr = nodes->GetZ() ;
//// 
////        const Points<NekDouble>::MatrixSharedPtrType Iptr = points->GetI(nNodes, zNodePtr);
////        const NekMatrix<NekDouble> & I = *Iptr;
//// 
////        int numPoints = points->GetNumPoints();
////        epsilon = eps*numPoints;
////        
////        // Check the derivative at each of the i points
////        for(int i = 0; i < I.GetRows(); ++i) {
////            long double exact = f(zNodePtr[i]);
////            long double numericInterpolation = 0.0;
//// 
////            // Multiply the derivative matrix with the sample point vector to get the derivative
////             for(int j = 0; j < I.GetColumns(); ++j) {
////                numericInterpolation += I(i,j) * f(z[j]);
////            }
//// 
////            // Compute the relative error; deals appropiately with the case when x' = 0
////            long double relativeError = (exact - numericInterpolation)/exact;
////            if( fabs(exact) < numeric_limits<double>::epsilon() ) {
////                relativeError = exact - numericInterpolation;
////                BOOST_CHECK( fabs(relativeError) < epsilon );
////            } else {
////                BOOST_CHECK_CLOSE(exact, numericInterpolation, 100.0*epsilon);
////            }
//// 
////            // Output the diagnostics if the verbose flag is set
////            if( bool isVerbose = true ) {
////                cout << "nPts: " << nPts << ",     relativeError = " << relativeError
////                << ",      relError / nPts = " << relativeError/(nPts + 2) << ",       epsilon = " << epsilon
////                <<",       exact = " << exact <<",       numerical = " << numericInterpolation<< endl; //}
//// 
//// 
////                long double expectedInterpolation = 0.0;
////                for(int j = 0; j<I.GetColumns(); ++j) {
////                    cout << "h[" << j << "](x[" << i << "]) = " << LagrangeInterpolation(j, zNodePtr[i], z, numPoints) << endl;
////                    expectedInterpolation += LagrangeInterpolation(j, zNodePtr[i], z, nPts) * f(z[j]);
////                }
////                cout << "expectedInterpolation = " << expectedInterpolation << endl;
//// 
////                    cout << "I = " << endl;
////                    for(int k=0; k<I.GetRows(); ++k) {
////                        for(int j =0; j<I.GetColumns(); ++j) {
////                            printf("% 5.3f  ", I(k, j));
////                        }
////                        cout << endl;
////                    }
//// 
////                long double L1_interpolationMatrixError = 0;
////                 double Linf = 0;
////                 long double RMS = 0, mean = 0;
//// 
////                cout << "Expected I = " << endl;
////                for(int k=0; k<I.GetRows(); ++k) {
////                    for(int j=0; j<I.GetColumns(); ++j) {
////                        long double L = LagrangeInterpolation(j, zNodePtr[k], z, nPts);
////                        long double error = I(k,j) - L;
////                        Linf = max(Linf, static_cast<double>(fabs(error)));
////                        L1_interpolationMatrixError += fabs(error);
////                        RMS += error*error;
////                        mean += error;
//// 
////                        printf("% 5.3Lf  ", L);
////                    }
////                    cout << endl;
////                }
////                int N = I.GetRows() * I.GetColumns();
////                RMS = sqrt(RMS)/N;
////                mean /= N;
////                long double variance = fabs(RMS*RMS - mean*mean);
////                cout << "Interpolation Matrix Errors:" << endl;
////                cout << "L1       = " << L1_interpolationMatrixError << endl;
////                cout << "Linf     = " << Linf << endl;
////                cout << "RMS      = " << RMS << endl;
////                cout << "mean     = " << mean << endl;
////                cout << "stdDev   = " << sqrt(variance) << endl;
////            }//
////        }
//    }
////   }//
//  }
//}
