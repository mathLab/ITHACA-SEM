/////////////////////////////////////////////////////////////////////////////////
////
//// File: testDerivation.cpp
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
//// Description: Testing for Derivation
//// Author: Sophia Han
////
/////////////////////////////////////////////////////////////////////////////////
//#undef NEKTAR_MEMORY_POOL_ENABLED 
//#include <boost/test/unit_test.hpp>
//#include <boost/test/test_tools.hpp>
//#include <boost/test/floating_point_comparison.hpp>
//
//
//#include <iostream>
//#include <limits>
//using namespace std;
//
//
//#include <UnitTests/testFoundation/testDerivation.h>
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
//
//namespace Nektar {
//
//  namespace foundationUnitTests{
//
//long double polyFunc(long double x) {
//    return  (((3.0*x*x - 5.0)*x + 1.0)*x - 2.0)*x + 3.0;
//}
//
//long double derivativeFunc(long double x) {
//    return ((15.0*x*x - 15.0)*x   + 2.0)*x - 2.0;
//}
//
//static long double polyFunc2(long double x) {
//    return  (((33.0*x*x + 27.0)*x*x - 7.0)*x*x + 5.0)*x*x*x*x  + 3.0;
//}
//
//long double derivativeFunc2(long double x) {
//    return (((330*x*x + 216)*x*x - 42)*x*x + 20)*x*x*x;
//}
//
//static long double polyFunc3(long double x) {
//    return x*x;
//}
//long double derivativeFunc3(long double x) {
//    return 2.0*x;
//}
//
//long double polyFunc4(long double x) {
//    return x*x*x*x;
//}
//long double derivativeFunc4(long double x) {
//    return 4.0*x*x*x;
//}
//
//void testDerivation(PointsType type, long double epsilon);
//void testDerivationForPolyES(PointsType type, long double epsilon);
//void testFourierDerivation(PointsType type, long double epsilon);
//
//void testDrivGaussGaussLegendre() {
//    testDerivation(eGaussGaussLegendre, 4e-14);
//   // cout<<"End of Derivation Test GaussGaussLegendre()" << endl << endl;
//}
//void testDrivGaussRadauMLegendre() {
//    testDerivation(eGaussRadauMLegendre, 9e-14);
//   // cout<<"End of Derivation Test GaussRadauMLegendre()" << endl << endl;
//}
//
//void testDrivGaussRadauPLegendre() {
//    testDerivation(eGaussRadauPLegendre, 9e-14);
//    //cout<<"End of Derivation Test GaussRadauPLegendre()" << endl << endl;
//}
//
//void testDrivGaussLobattoLegendre() {
//    testDerivation(eGaussLobattoLegendre, 9e-14);
//   // cout<<"End of Derivation Test GaussLobattoLegendre()" << endl << endl;
//}
//
//void testDrivGaussGaussChebyshev() {
//    testDerivation(eGaussGaussChebyshev, 1e-13);
//   // cout<<"End of Derivation Test eGaussGaussChebyshev()" << endl << endl;
//}
//
//void testDrivGaussRadauMChebyshev() {
//    testDerivation(eGaussRadauMChebyshev, 1e-13);
//   // cout<<"End of Derivation Test GaussRadauMChebyshev()" << endl << endl;
//}
//
//void testDrivGaussRadauPChebyshev() {
//    testDerivation(eGaussRadauPChebyshev, 1e-13);
//    //cout<<"End of Derivation Test GaussRadauPChebyshev()" << endl << endl;
//}
//
//void testDrivGaussLobattoChebyshev() {
//    testDerivation(eGaussLobattoChebyshev, 1e-13);
//    //cout<<"End of Derivation Test GaussLobattoChebyshev()" << endl << endl;
//}
//
//void testDrivGaussRadauMAlpha0Beta1() {
//    testDerivation(eGaussRadauMAlpha0Beta1, 9e-14);
//   // cout<<"End of Derivation Test GaussRadauMAlpha0Beta1()" << endl << endl;
//}
//
//void testDrivGaussRadauMAlpha0Beta2() {
//    testDerivation(eGaussRadauMAlpha0Beta2, 1e-13);
//   // cout<<"End of Derivation Test GaussRadauMAlpha0Beta2()" << endl << endl;
//}
//
//void testDrivPolyEvenlySpaced() {
//    testDerivation(ePolyEvenlySpaced, 1e-13);
//   // cout<<"End of Derivation Test PolyEvenlySpaced()" << endl << endl;
//}
//
//void testDrivFourierEvenlySpaced() {
//    testFourierDerivation(eFourierEvenlySpaced, 1e-10);
//    cout<<"End of Derivation Test PolyEvenlySpaced()" << endl << endl;
//}
//
//long double derivativeLegrangePoly(const Array<OneD, const NekDouble> &x, int N, int i, int j) {
//    long double derivative = 0.0;
//    if( i == j ) {
//        long double sum = 0.0;
//        for( int n = 0; n < N; ++n ) {
//            if( n != j ) {
//                sum += 1.0/(x[j] - x[n]);
//            }
//        }
//        derivative = sum;
//    } else {
//        long double product = 1.0;
//        for( int n = 0; n < N; ++n ) {
//            if( n != i  &&  n != j ) {
//                product *= (x[i] - x[n])/(x[j] - x[n]);
//            }
//        }
//        derivative = product / (x[j] - x[i]);
//    }
//    return derivative;
//}
//
//// IDFT of the derivative of the Fourier transform of an arbitrary function
//long double derivativeMatrixFourier(const Array<OneD, const NekDouble> &x, int N, int i, int j) {
//    long double derivative = 0.0;
//    if( i != j ) {
//        int k = i - j;
//        long double x = (long double)2.0 * M_PI * k / N;
//        long double sign = (k&1) ? -1.0 : 1.0;
//        derivative = M_PI * 0.5 * sign / tan(x/2.0);
//    } else {
//        derivative = 0.0;
//    }
//    return derivative;
//}
//
//// fourier function using Trapezoidal rule
//static long double fourierFunc(long double x, int N) {
//    long double z = M_PI*(x + 1.0);
//    return (cos(N/2.0*z) + sin((N/2.0 - 2.0)*z))*M_PI;
//}
// long double derivativeFourierFunc(long double x, int N) {
//    long double z = M_PI*(x + 1.0);
//    long double a = (N-4.0)/2.0;
//    long double b = N/2.0;
//    return M_PI*M_PI*(a*cos(a*z) - b*sin(b*z));
//}
//
//void testFourierDerivation(PointsType type, long double epsilon) {
//     long double (*f)(long double x, int N) = fourierFunc;
//     long double (*d)(long double x, int N) = derivativeFourierFunc;
//
//    const long double eps = epsilon;
//
//    for(int nPts = 4; nPts<=20; nPts += 2) {
//        const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
//        const Points<NekDouble>::MatrixSharedPtrType Dptr = points->GetD();
//        const NekMatrix<NekDouble> & D = *Dptr;
//        const Array<OneD, const NekDouble> &z = points->GetZ();
//
//               
//        int numPoints = points->GetNumPoints();
//        
//        // Check the derivative at each of the i points
//        for(int i = 0; i < numPoints; ++i) {
//            long double exact = d(z[i], numPoints);
//            long double numericDerivative = 0.0;
//
//            // Multiply the derivative matrix with the sample point vector to get the derivative
//            for(int j = 0; j < numPoints; ++j) {
//                numericDerivative += D(i,j) * f(z[j], numPoints);
//            }
//
//            // Compute the relative error; deals appropiately with the case when x' = 0
//            long double relativeError = (exact - numericDerivative)/exact;
//            if( fabs(exact) < numeric_limits<double>::epsilon() ) {
//                relativeError = exact - numericDerivative;
//                BOOST_CHECK( fabs(relativeError) < epsilon );
//            } else {
//                BOOST_CHECK_CLOSE(exact, numericDerivative, 100.0*epsilon);
//            }
//
//            // Output the diagnostics if the verbose flag is set
//            if( bool isVerbose = true ) {
////                 cout << "nPts: " << nPts << ",     relativeError = " << relativeError
////                 << ",      relError / nPts = " << relativeError/(nPts + 2) << ",       epsilon = " << epsilon
////                 <<",       exact = " << exact <<",       numerical = " << numericDerivative<< endl; 
//
//                long double expectedDerivative = 0.0;              
//                for(int j = 0; j < numPoints; ++j) {
//                     expectedDerivative += derivativeMatrixFourier(z, numPoints, i,j) * f(z[j], numPoints);
//                }
//          //      cout << "expectedFourierDerivative = " << expectedDerivative << endl;
//
//                if( i == numPoints - 1 ) {
//               //     cout << "D = " << endl;
//                    for( int k = 0; k < numPoints; ++k ) {
//                        for( int j = 0; j < numPoints; ++j ) {
//                      //      printf("% 5.3f  ", D(k,j));
//                        }
//                        cout << endl;
//                    }
//                //    cout << "Expected D = " << endl;
//                    for( int k = 0; k < numPoints; ++k ) {
//                        for( int j = 0; j < numPoints; ++j ) {
//                 //            printf("% 10.3Lf", derivativeMatrixFourier(z, numPoints, k, j));
//                        }
//                        cout << endl;
//                    }
//                  //  cout << "y  = [";
//                    for( int j = 0; j < numPoints; ++j ) {
//                //        printf("% 4.3Lf", f(z[j], numPoints));
//                        if( j < numPoints-1 ) {
//               //             cout << ", ";
//                        }
//                    }
//                  //  cout << "]" << endl;
//                  //  cout << "y' = [";
//                    for( int j = 0; j < numPoints; ++j ) {
//                    //    printf("% 4.3Lf", d(z[j], numPoints));
//                        if( j < numPoints-1 ) {
//                      //      cout << ", ";
//                        }
//                    }
//                    //cout << "]" << endl;
//                }
//            }//
//        }
//    }
//} //
//
//
//void testDerivation(PointsType type, long double epsilon) {
//
//    long double (*f)(long double x) = polyFunc3;
//    long double (*d)(long double x) = derivativeFunc3;
//
//    const long double eps = epsilon;
//
//    for(int nPts = 4; nPts<=10; ++nPts) {
//        const boost::shared_ptr<Points<double> > points = PointsManager()[PointsKey(nPts, type)];
//        const Points<double>::MatrixSharedPtrType Dptr = points->GetD();
//        const NekMatrix<double> & D = *Dptr;
//        const Array<OneD, const NekDouble> &z = points->GetZ();
//         
//        int numPoints = points->GetNumPoints();
//        epsilon = eps*numPoints;
//
//        // Check the derivative at each of the i points
//        for(int i = 0; i < numPoints; ++i) {
//            long double exact = d(z[i]);
//            long double numericDerivative = 0.0;
//
//            // Multiply the derivative matrix with the sample point vector to get the derivative
//            for(int j = 0; j < numPoints; ++j) {
//                numericDerivative += D(i,j) * f(z[j]);
//            }
//
//
//            // Compute the relative error; deals appropiately with the case when x' = 0
//            long double relativeError = (exact - numericDerivative)/exact;
//            if( fabs(exact) < numeric_limits<double>::epsilon() ) {
//                relativeError = exact - numericDerivative;
//                BOOST_CHECK( fabs(relativeError) < epsilon );
//            } else {
//                BOOST_CHECK_CLOSE(exact, numericDerivative, 100.0*epsilon);
//            }
//
//            // Output the diagnostics if the verbose flag is set
//            if( bool isVerbose = true ) {
//                //cout << "nPts: " << nPts << ",     relativeError = " << relativeError
////                 << ",      relError / nPts = " << relativeError/(nPts + 2) << ",       epsilon = " << epsilon
////                 <<",       exact = " << exact <<",       numerical = " << numericDerivative<< endl; 
//
//                long double expectedDerivative = 0.0;
//                for(int j = 0; j < numPoints; ++j) {
//                     expectedDerivative += derivativeLegrangePoly(z, numPoints, i,j) * f(z[j]);
//                }
//            //    cout << "expectedLegrangeDerivative = " << expectedDerivative << endl;
//
//                if( i == numPoints - 1 ) {
//              //      cout << "D = " << endl;
//                    for( int k = 0; k < numPoints; ++k ) {
//                        for( int j = 0; j < numPoints; ++j ) {
//                //            printf("% 5.3f  ", D(k,j));
//                        }
//                        cout << endl;
//                    }
//                    if( bool isEvenlySpaced = true ) {
//                  //      cout << "Expected D = " << endl;
//                        for( int k = 0; k < numPoints; ++k ) {
//                            for( int j = 0; j < numPoints; ++j ) {
//                    //             printf("% 10.3Lf", derivativeLegrangePoly(z, numPoints, k, j));
//                            }
//                            cout << endl;
//                        }
//                    }
////                    cout << "y  = [";
//                    for( int j = 0; j < numPoints; ++j ) {
//  //                      printf("% 4.3Lf", f(z[j]));
//                        if( j < numPoints-1 ) {
//    //                        cout << ", ";
//                        }
//                    }
//      //              cout << "]" << endl;
//      //              cout << "y' = [";
//                    for( int j = 0; j < numPoints; ++j ) {
//        //                printf("% 4.3Lf", d(z[j]));
//                        if( j < numPoints-1 ) {
//          //                  cout << ", ";
//                        }
//                    }
//             //       cout << "]" << endl;
//                }
//            }//
//        }
//    }
//} //
//
//}
//}
//
//
