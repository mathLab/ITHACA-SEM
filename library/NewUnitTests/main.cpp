////////////////////////////////////////////////////////////////////////////////
//
// main.cpp
// Sophia

// Main class for the LibUtilities::Points unit tests.
//
////////////////////////////////////////////////////////////////////////////////



#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/test/included/unit_test_framework.hpp>
#include <boost/test/floating_point_comparison.hpp>

using boost::unit_test_framework::test_suite;


#include <iostream>
using namespace std;



#include <LibUtilities/BasicUtils/NekManager.hpp>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/GaussPoints.h>
#include <LibUtilities/Foundations/PolyEPoints.h>
#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Foundations/NodalTriFekete.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace Nektar;
using namespace Nektar::LibUtilities;


namespace Nektar
{
    namespace MyNewUnitTests
    {
        void simpleTest();
        void testPolyFunc();
    }
}


// The boost unit test framework provides the main function for us.
// All we need to do is provide a test suite.
test_suite* init_unit_test_suite( int, char* [] )
{
    test_suite* test= BOOST_TEST_SUITE( "Nektar++ Test Suite" );

    test->add(BOOST_TEST_CASE(&Nektar::MyNewUnitTests::simpleTest), 0);
    test->add(BOOST_TEST_CASE(&Nektar::MyNewUnitTests::testPolyFunc),0);

    return test;
}

namespace Nektar
{
    namespace MyNewUnitTests
    {
        void simpleTest()
        {
            int size = 10;

            //PointsKey gaussKey(size, eGaussLobattoLegendre);

            PointsKey gaussKey(size, eGaussGaussLegendre);
            PointsKey polyKey(size, ePolyEvenlySpaced);
            boost::shared_ptr<Points<NekDouble> > ptr = PointsManager()[gaussKey];
            boost::shared_ptr<NekMatrix<NekDouble> > mat = ptr->GetI(gaussKey);

            int m = mat->GetRows();
            int n = mat->GetColumns();
            cout << "(m, n) = (" << m << ", " << n << ")" <<  endl;


            BOOST_CHECK(mat->GetRows() == size); // Fails here
            BOOST_CHECK(mat->GetColumns() == size);

            Points<double> & points = *ptr;
            cout << "points.GetPointsDim()    = " << points.GetPointsDim() << endl;
            cout << "points.GetNumPoints()    = " << points.GetNumPoints() << endl;
            cout << "points.GetTotNumPoints() = " << points.GetTotNumPoints() << endl;
            const double *z = points.GetZ()->data();
            for( int i=0; i<size; ++i )
            {
                cout << "z[i] = " << z[i] <<  endl;
            }
            cout << endl;
            const double *w = points.GetW()->data();
            for( int i=0; i<size; ++i )
            {
                cout << "w[i] = " << w[i] <<  endl;
            }

            cout << "Happy Testing!" << endl;

        }

    }

    namespace TestUtilities
    {

        vector<double> generatePolynomial(int degree) {
            //             double a[] = {
            //                 1.0, 1.4, -1.5, 1.2, -1.3, 1.5, -0.1, 1.4, -3.2, 2.4, -1.0, 1.6, -1.3, 4.5,
            //                 1.3, 1.9, 1.6, 1.3, 1.4, 1.7, 1.9, 0.3, 1.6, // 23
            //                 1.2, 0.5, 1.0, 0.3, 0.5, 0.7, 0.4, 0.6 }; // 31
            double a[] = {
                1, 1,  3,  3, 5, 5,  7, 7, 9,  9,
                11,11,13, 13,15,15, 17,17,19, 19,
                21,21,23, 23,25,25, 27,27,29, 29,
                31,31,33, 33,35,35, 37,37,39, 39,
                41,41,43, 43,45,45, 47,47,49, 49,
                51,51,53, 53,55,55, 57,57,59, 59,
                61
            };

            vector<double> coefficients(a, a + degree + 1);
            return coefficients;
        }

        /// Taylor's expansion for sqrt(1-x^2)
        vector<double> generateChebyshevPolynomial(int degree) {
            vector<double> coefficients;

            // Set up the basis cases
            coefficients.push_back(1);
            if( degree > 0 ) {
                coefficients.push_back(0);
            }

            // Recursively define the coefficients
            for(int n=2; n <= degree; ++n ) {
                coefficients.push_back( (1.0 - 3.0/n) * coefficients[n-2] );
            }

            return coefficients;
        }

        // Evaluates at x the polynomial that is given by the coefficents
        long double evaluatePolynomial(long double x, const vector<double> & coefficients) {
            int N = coefficients.size();
            long double y = coefficients[N-1];
            for( int i = N-2; i >= 0; --i ) {
                y = coefficients[i] + x*y;                
            }
            return y;
        }

        // Integrates the polynomial from [-1,1]
        long double integrate(const vector<double> & coefficients) {
            int M = coefficients.size(), N = M - 1;
            long double integral = 0;
            for( int n = 0; n <= N/2; ++n ) {
                integral += coefficients[2*n] / (2.0*n + 1.0);
                //                 if( N == 4 ) {
                //                     cout << "N/2 = " << N/2 << ", 2*n = " << 2*n << ", coefficients[2*n] = " << coefficients[2*n] << ", (2*n + 1) = " << (2*n + 1) << ", coefficients[2*n] / (2*n + 1) = " << coefficients[2*n] / (2*n + 1) << endl;
                //                 }
            }
            return  2.0 * integral;
        }

        void printPolynomial(const vector<double> & a) {
            cout << a[0];
            if( a.size() > 1 ) {
                cout << " + " << a[1] << "x";
            }
            for( int n = 2; n < a.size(); ++n ) {
                cout << " + " << a[n] << "x^" << n;
            }
        }
    }
    namespace MyNewUnitTests
    {
        void testPolynomialOnWeights(const boost::shared_ptr<Points<double> > points, const vector<double> & polynomialCoefficients, bool isVerbose = false) {

            ConstNekDouble1DSharedArray zt, wt;
            points->GetZW(zt, wt);
            const double *z = zt->data(), *w = wt->data();
            int numPoints = points->GetNumPoints();
            int startDegree = 0;
            if( points->GetPointsType() == eGaussGaussChebyshev ) {
                startDegree = polynomialCoefficients.size() - 1;
            }

            for(int i = startDegree; i < polynomialCoefficients.size(); ++i) {
                vector<double> a( polynomialCoefficients.begin(), polynomialCoefficients.begin() + i + 1 );
                //long double analyticIntegral = round(TestUtilities::integrate(a));
                long double analyticIntegral = TestUtilities::integrate(a);
                long double numericIntegral = 0;
                for(int j = 0; j < numPoints; ++j) {
                    long double summand = w[j] * TestUtilities::evaluatePolynomial( z[j], a );
                    if( points->GetPointsType() == eGaussGaussChebyshev ) {
                        summand *= sqrt(1.0 - z[j]*z[j]);
                    }
                    numericIntegral += summand;
                    //if( isVerbose ) cout << "w["<<j<<"] = " << w[j] << ", z["<<j<<"] = " << z[j] << endl;
                }
                if( isVerbose ) {
                    cout << "Points = " << numPoints << ", Polynomial degree = " << i << ", nCoef = " << a.size();
                    printf(", Integral = %1.20f, Computed area = %1.20f", (double)analyticIntegral, (double)numericIntegral);
                    cout << ", Error = " << analyticIntegral - numericIntegral;
                    if( polynomialCoefficients.size() < 10 ) {cout << ", P_" << i << " = "; TestUtilities::printPolynomial(a);}
                    cout << ", z[0] = " << z[0] << ", w[0] = " << w[0] << ", z[n-1] = " << z[numPoints-1] << ", w[n-1] = " << w[numPoints-1];
                    cout << endl;
                }

                //                 BOOST_CHECK_CLOSE( numericIntegral, analyticIntegral, 1e-16 );
                BOOST_CHECK_CLOSE( numericIntegral, analyticIntegral, 1e-12 );
            }

            //             if( isVerbose ) {
            //                 cout << " End of the polynomial Unit Test" << endl;
            //             }
        }


        const int MaxNumberOfPoints = 15;

        void testPolyFunc()
        {
            vector<double> coefficients;
            PointsType type;
            const char * cp = "";

            type = eGaussGaussLegendre;
            cp = "Testing eGaussGaussLegendre";
            cout << cp << endl;
            BOOST_CHECKPOINT(cp);
            for( int i = 1; i < MaxNumberOfPoints; ++i ) {
                int nPts = i, n = nPts - 1, degree = 2*n + 1;
                coefficients = TestUtilities::generatePolynomial(degree);
                testPolynomialOnWeights( PointsManager()[PointsKey(nPts, type)], coefficients );
            }
            cout << "End of testing eGaussGaussLegendre" <<endl;
            cout << " " << endl;

            type = eGaussLobattoLegendre;
            cp = "Testing eGaussLobattoLegendre";
            cout << cp << endl;
            BOOST_CHECKPOINT(cp);
            for( int i = 2; i < MaxNumberOfPoints; ++i ) {
                int nPts = i, n = nPts - 1, degree = 2*n - 1;
                coefficients = TestUtilities::generatePolynomial(degree);
                testPolynomialOnWeights( PointsManager()[PointsKey(nPts, type)], coefficients );
            }
            cout << "End of testing eGaussLobattoLegendre" <<endl;
            cout << " " << endl;

            type = eGaussRadauMLegendre;
            cp = "Testing eGaussRadauMLegendre";
            cout << cp << endl;
            BOOST_CHECKPOINT(cp);
            for( int i = 1; i< MaxNumberOfPoints; ++i) {
                int nPts = i, n = nPts - 1, degree = 2*n - 1;
                coefficients = TestUtilities::generatePolynomial(degree);
                testPolynomialOnWeights( PointsManager()[PointsKey(nPts, type)], coefficients );
            }
            cout << "End of testing eGaussRadauMLegendre" << endl;
            cout << " " << endl;

            type = eGaussRadauPLegendre;
            cp = "Testing eGaussRadauPLegendre";
            cout << cp << endl;
            BOOST_CHECKPOINT(cp);
            for( int i=1; i<MaxNumberOfPoints; ++i){
                int nPts = i, n = nPts - 1, degree = 2*n;
                coefficients = TestUtilities::generatePolynomial(degree);
                testPolynomialOnWeights( PointsManager()[PointsKey(nPts, type)], coefficients );
            }
            cout << "End of testing eGaussRadauPLegendre" << endl;
            cout << " " <<endl;


            type = eGaussGaussChebyshev;
            cp = "Testing eGaussGaussChebyshev";
            cout << cp << endl;
            BOOST_CHECKPOINT(cp);
            //for( int i=1; i<MaxNumberOfPoints; ++i){
            {
                int i = 90;
                //int nPts = i, n = nPts - 1, degree = max(0, 2*n - 3);
                int nPts = i, n = nPts - 1, degree = 500000;
                coefficients = TestUtilities::generateChebyshevPolynomial(degree);
                testPolynomialOnWeights( PointsManager()[PointsKey(nPts, type)], coefficients, true );
            }
            cout << "End of testing eGaussGaussChebyshev" << endl;
            cout << " " << endl;



        }
    }

}

//    enum PointsType
// 00075         {
// 00076             eNoPointsType,
// 00077             eGaussGaussLegendre,     
// 00078             eGaussRadauMLegendre,    
// 00079             eGaussRadauPLegendre,    
// 00080             eGaussLobattoLegendre,   
// 00081             eGaussGaussChebyshev,    
// 00082             eGaussRadauMChebyshev,   
// 00083             eGaussRadauPChebyshev,   
// 00084             eGaussLobattoChebyshev,  
// 00085             eGaussRadauMAlpha0Beta1, 
// 00086             eGaussRadauMAlpha0Beta2, 
// 00087             ePolyEvenlySpaced,       
// 00088             eFourierEvenlySpaced,    
// 00089             eNodalTriElec,           
// 00090             eNodalTriFekete,         
// 00091             eNodalTetElec,           
// 00092             SIZE_PointsType          
// 00093         };

/**
$Log: main.cpp,v $
Revision 1.7  2007/03/09 03:07:39  ehan
*** empty log message ***

Revision 1.5  2007/03/06 20:29:23  ehan
test passed eGaussRadauMLegendre, eGaussRadauPLegendre,         eGaussGaussLegendre, eGaussLobattoLegendre

Revision 1.4  2007/03/02 01:10:44  ehan
fixed some bugs

Revision 1.2  2007/03/01 14:09:14  ehan
*** empty log message ***

Revision 1.1  2007/03/01 04:36:20  ehan
Kdevelop for UnitTest


Revision 1.1  2007/2/26 04:05:23  ehan
Added the NewUnitTest project.

**/
