////////////////////////////////////////////////////////////////////////////////
//
// main.cpp
// Sophia
// Blake Nelson
//
// Driver for the unit tests.
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
            boost::shared_ptr<Points<double> > ptr = PointsManager()[gaussKey];
            boost::shared_ptr<NekMatrix<double> > mat = ptr->GetI(gaussKey);

            int m = mat->GetRows();
            int n = mat->GetColumns();
            cout << "(m, n) = (" << m << ", " << n << ")" <<  endl;


            BOOST_CHECK(mat->GetRows() == size); // Fails here
            BOOST_CHECK(mat->GetColumns() == size);



            Points<double> & points = *ptr;
            cout << "points.GetPointsDim()    = " << points.GetPointsDim() << endl;
            cout << "points.GetNumPoints()    = " << points.GetNumPoints() << endl;
            cout << "points.GetTotNumPoints() = " << points.GetTotNumPoints() << endl;
            double *z = points.GetZ();
            for( int i=0; i<size; ++i )
            {
                cout << "z[i] = " << z[i] <<  endl;
            }
            cout << endl;
            double *w = points.GetW();
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
//                 -1.3, 1.4, -1.5, 1.2, -1.3, 1.5, -0.1, 1.4, -3.2, 2.4, -1.0, 1.6, -1.3, 4.5,
//                 1.3, 1.9, 1.6, 1.3, 1.4, 1.7, 1.9, 0.3, 1.6, // 23
//                 1.2, 0.5, 1.0, 0.3, 0.5, 0.7, 0.4, 0.6 }; // 31
            double a[] = {
                1, 1,  3,  3, 5, 5,  7, 7, 9,  9,
                11,11,13, 13,15,15, 17,17,19, 19,
                21,21,23, 23,25,25, 27,27,29, 29,
                31
            };

            vector<double> coefficients(a, a + degree + 1);
            return coefficients;
        }
        
        // Evaluates at x the polynomial that is given by the coefficents
        double evaluatePolynomial(double x, const vector<double> & coefficients) {
            int N = coefficients.size();
            double y = coefficients[N];
            for( int i = N-1; i >= 0; --i ) {
                y = coefficients[i] + x*y;                
            }
            return y;
        }
        
        // Integrates the polynomial from [-1,1]
        double integrate(const vector<double> & coefficients) {
            int N = coefficients.size();
            double integral = 0;
            for( int i = 0; i <= (N-1)/2; ++i ) {
                integral += coefficients[2*i] / (2*i + 1);
                if( N == 8 ) {
                    cout << "N/2 = " << N/2 << ", 2*i = " << 2*i << ", coefficients[2*i] = " << coefficients[2*i] << ", (2*i + 1) = " << (2*i + 1) << ", coefficients[2*i] / (2*i + 1) = " << coefficients[2*i] / (2*i + 1) << endl;
                }
            }
            return  2.0 * integral;
        }
        
        
    }
        namespace MyNewUnitTests
    {
        void testPolynomialOnWeights(const boost::shared_ptr<Points<double> > points, const vector<double> & polynomialCoefficients) {
            const double *z, *w;
            points->GetZW(z, w);
            int numPoints = points->GetNumPoints();
            
            for(int i = 0; i < polynomialCoefficients.size(); ++i) {
                vector<double> a(polynomialCoefficients.begin(), polynomialCoefficients.begin() + i + 1);
                double analyticIntegral = TestUtilities::integrate(a);
                double numericIntegral = 0;
                for(int j = 0; j < numPoints; ++j) {
                    numericIntegral += w[j] * TestUtilities::evaluatePolynomial( z[j], a ); 
                }
                cout << "Points = " << numPoints << ", Polynomial degree = " << i;
                cout << ", Integral = " << analyticIntegral << ", Computed area = " << numericIntegral << endl;
                
                BOOST_CHECK_CLOSE( numericIntegral, analyticIntegral, 1e-11 );
            }
        }
//         const int MaxNumberOfPoints = 15;
        const int MaxNumberOfPoints = 31;
        void testPolyFunc()
        {
//             int P[] = {2,4,8,12};
//             int N[] = {2,3,5,7};

            vector<double> coefficients;
            PointsType type;

            BOOST_CHECKPOINT("Testing eGaussGaussLegendre");
            type = eGaussGaussLegendre;
            for( int i = 1; i < MaxNumberOfPoints; ++i ) {
                int nPts = i;
              //  int degree = 2*i - 2;                
               int degree = 2*i - 1;
//                 int degree = 2*int(i/2);
                coefficients = TestUtilities::generatePolynomial(degree);
                testPolynomialOnWeights( PointsManager()[PointsKey(nPts, type)], coefficients );
            }

            BOOST_CHECKPOINT("Testing eGaussLobattoLegendre");
            type = eGaussLobattoLegendre;
            for( int i = 1; i < MaxNumberOfPoints; ++i ) {
                int nPts = i, degree = max(1, 2*i - 3);
                coefficients = TestUtilities::generatePolynomial(degree);
                testPolynomialOnWeights( PointsManager()[PointsKey(nPts, type)], coefficients );
            }
        }
    }

}



/**
    $Log: main.cpp,v $
    Revision 1.2  2007/03/01 14:09:14  ehan
    *** empty log message ***

    Revision 1.1  2007/03/01 04:36:20  ehan
    Kdevelop for UnitTest
 
 
    Revision 1.1  2007/2/26 04:05:23  ehan
    Added the NewUnitTest project.
 
**/
