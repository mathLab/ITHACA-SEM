#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>

using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/ManagerAccess.h>  // for PointsManager, etc

using namespace Nektar;
using namespace Nektar::LibUtilities;


// Quintic polynomial
long double polyFunc(long double x)
{
    return  (((3.0*x*x - 5.0)*x + 1.0)*x - 2.0)*x + 3.0;
}

// Derivative of the Quintic polynomial
long double derivativePolyFunc(long double x)
{
    return ((15.0*x*x - 15.0)*x   + 2.0)*x - 2.0;
}

// A Fourier function that integrates to 0 with the Trapezoidal rule
long double fourierFunc(long double x, int N)
{
    long double z = M_PI*(x + 1.0);
    return (cos(N/2.0*z) + sin((N/2.0 - 2.0)*z))/M_PI;
}

// Derivative of the Fourier function above
long double derivativeFourierFunc(long double x, int N)
{
    long double z = M_PI*(x + 1.0);
    long double a = (N-4.0)/2.0;
    long double b = N/2.0;
    return M_PI*M_PI*(a*cos(a*z) - b*sin(b*z));
}

long double func(long double x, int N, PointsType type)
{
    long double y = 0;
    if( type == eFourierEvenlySpaced)
    {
        y = fourierFunc(x,N);
    } else
    {
        y = polyFunc(x);
    }
    return y;
}

long double derivativeFunction(long double x, int N, PointsType type)
{
    long double yDerivative = 0;
    if(type == eFourierEvenlySpaced)
    {
         yDerivative = derivativeFourierFunc(x, N);
    }
    else
    {
         yDerivative = derivativePolyFunc(x);
    }
    return yDerivative;
}

long double integrationFunction(int nPts, PointsType type)
{
    boost::ignore_unused(nPts);

   long double integral = 0;
   switch(type)
   {
       case  eFourierEvenlySpaced:
          integral = 0;
       break;
       default:
          integral = (long double)20.0/3.0;
   }
   return integral;
}

long double integrandWeightFunction(long double x, PointsType type)
{
    long double weight = 1;

    switch(type)
    {
        case eGaussGaussChebyshev:
        case eGaussRadauMChebyshev:
        case eGaussRadauPChebyshev:
        case eGaussLobattoChebyshev:
            weight = 1.0 / sqrt(1.0 - x*x);
        break;

        case eGaussRadauMAlpha0Beta1:
            // weight = 1.0 + x; // ?
        break;

        case eGaussRadauMAlpha0Beta2:
            // weight = (1.0 + x)*(1.0 + x); // ?
        break;


        default:
            weight = 1.0;
    }
    return weight;
}
// This routine projects a polynomial or trigonmetric functions which
// has energy in all modes of the expansions and report an error.

int main(int argc, char *argv[])
{

    // Argument check: Display a help message if the count is wrong
    if(argc != 3)
    {
        cerr << "Usage: FoundationDemo Points1D-Type nPts" << endl;

        cerr << "Where type is an integer value which dictates the basis as:" << endl;
        for(int i=0; i<SIZE_PointsType; ++i) {
            cerr << setw(30) << kPointsTypeStr[i] << " = " << i << endl;
        }

        cerr << "Note type = 18, 19, and 20 are for two dimensional bases" << endl;
        cerr << "\nExample: FoundationDemo 5 6" << endl;
        cerr << "\n\t Tests GaussGaussChebyshev on 6 points" << endl;

        return 1; // Aborts main() function
    }

    // Read in the type for the points from the caller
    PointsType pointsType = (PointsType) atoi(argv[1]);
    if(pointsType == eNoPointsType)
    {
        cerr << "pointsType = " << pointsType << endl;
        cerr << "kPointsTypeStr["<<pointsType<<"] = " << kPointsTypeStr[pointsType] << endl;
        ErrorUtil::Error(ErrorUtil::efatal,__FILE__, __LINE__,
                  "No Points Type requested",0);
    }

    // Read in the number of points at which the function will be evaluated
    int nPts = atoi(argv[2]);

    // Show the set up to the user
    cout << "Points type:               " << kPointsTypeStr[pointsType] << endl;
    cout << "Number of points:          " << nPts << endl;

    // Display the example test function to the user
    if( pointsType == eFourierEvenlySpaced )
    {
        cout << "Trigonometric function:    ";
        cout << "f(x) = (cos(N/2.0*z) + sin((N/2.0 - 2.0)*z))/M_PI, where z = M_PI*(x + 1.0)" << endl;
    }
    else
    {
        cout << "Quintic polynomial:        ";
        cout << "p(x) = (((3.0*x*x - 5.0)*x + 1.0)*x - 2.0)*x + 3.0" << endl;
    }

      // Obtain a reference to a Points object via an appropriate PointsKey object
    PointsKey key(nPts, pointsType);

    PointsSharedPtr points = PointsManager()[key];
    //std::shared_ptr<Points<NekDouble> > points = PointsManager()[key];
    //const ptr<Points<NekDouble> > points = PointsManager()[key];


    // Get the abscissas and their matching quadrature weights
   //    SharedArray<const NekDouble> z, w;
    Array<OneD, const NekDouble> z, w;
    points->GetZW(z,w);

    // Evaluate the example function at the z[i] points in the interval [-1,1].
    // This generates the data samples, which we store in the y vector and use later
    // during interpolation/integration/differentiation
    //    SharedArray<NekDouble> y(nPts);
    Array<OneD, NekDouble> y(nPts);
    for(int i = 0; i < nPts; ++i)
    {
        y[i] = func( z[i], nPts, pointsType );
    }



    // /////////////////////////////////////////////////////////
    //                     Interpolation                      //
    // /////////////////////////////////////////////////////////

    // Generate a list of interpolating nodes
    int nNodes = 2*nPts; // Number of interpolating nodes
    std::shared_ptr<Points<NekDouble> > nodes = PointsManager()[PointsKey(nNodes, pointsType)];
    Array<OneD, const NekDouble> zNode = nodes->GetZ();

    // Get the interpolation matrix I
    // Note that I is 2N rows by N columns
    const Points<NekDouble>::MatrixSharedPtrType Iptr = points->GetI(nNodes,zNode);
    const NekMatrix<NekDouble> & I = *Iptr;

    // Interpolate the data values in the y vector using the interpolation matrix I
    Array<OneD, NekDouble> u(I.GetRows());
    for(int i = 0; i < int(I.GetRows()); ++i)
    {
        u[i] = 0;
        for(int j = 0; j < int(I.GetColumns()); ++j)
        {
            u[i] += I(i,j) * y[j];
        }
    }

    // Display the original samples
    cout << setprecision(3);
    cout << "\nOriginal data: \nx      = ";
    for(int i = 0; i < nPts; ++i)
    {
        cout << setw(6) << z[i] << " ";
    }
    cout << "\ny      = ";
    for(int i = 0; i < nPts; ++i)
    {
        cout << setw(6) << y[i] << " ";
    }

    // Display the interpolated data
    cout << "\n\n\n              **** Interpolation ****";
    cout << "\n\nResults of interpolation with " << kPointsTypeStr[pointsType] << ":";
    cout << "\ninterp = ";
    for(int i = 0; i < nNodes; ++i)
    {
        cout << setw(6) << u[i] << " ";
    }

    // Determine the exact solutions
    cout << "\nexact  = ";
    for(int i = 0; i < nNodes; ++i)
    {
        cout << setw(6) << func(zNode[i], nPts, pointsType) << " ";
    }

    // Display the pointwise error
    cout << setprecision(1);
    cout << "\nerror  = ";
    long double Linf = 0, RMS = 0;
    for(int i = 0; i < int(I.GetRows()); ++i)
    {
        //long double exact = function(zNode[i], nNodes, pointsType);
        long double exact = func(zNode[i], nPts, pointsType);
        long double error = exact - u[i];
        Linf = max(Linf, fabs(error));
        RMS += error*error;
        long double epsilon = 1e-2;
        if( fabs(exact) > epsilon )
        {
            error /= exact;
        }
        cout << setw(6) << error << " ";
    }
    RMS = sqrt(RMS) / int(I.GetRows());
    cout << setprecision(6);
    cout << "\nLinf   = " << setw(6) << Linf;
    cout << "\nRMS    = " << setw(6) << RMS << endl;

    // Show the interpolation matrix
    cout << "\nI = " << endl;
    for(int i = 0; i < int(I.GetRows()); ++i)
    {
        cout << "     ";
        for(int j = 0; j < int(I.GetColumns()); ++j)
        {
            printf("% 5.3f  ", I(i,j));
        }
        cout << "" << endl;
    }


    ///////////////////////////////////////////////////////////
    //                 Galkerin Projection                   //
    ///////////////////////////////////////////////////////////

    // Generate a list of projection nodes
    Array<OneD, NekDouble> yNode(zNode.size());

    for(int i = 0; i < nNodes; ++i)
    {
        yNode[i] = func( zNode[i], nNodes, pointsType );
    }

    PointsKey key1(nNodes, pointsType);

    // Note that I is 2N rows by N columns
    const Points<NekDouble>::MatrixSharedPtrType GPptr = points->GetGalerkinProjection(key1);
    const NekMatrix<NekDouble> & GP = *GPptr;

    // Project the data values in the yNode vector using the projection matrix GP
    for(int i = 0; i < int(GP.GetRows()); ++i)
    {
        u[i] = 0;
        for(int j = 0; j < int(GP.GetColumns()); ++j)
        {
            u[i] += GP(i,j) * yNode[j];
        }
    }

    // Display the original samples
    cout << setprecision(3);
    cout << "\n\n\n              **** Galerkin Project ****";
    cout << "\n\nResults of Galkerin Project with " << kPointsTypeStr[pointsType] << ":";

    cout << "\nOriginal data: \nx      = ";
    for(int i = 0; i < nPts; ++i)
    {
        cout << setw(6) << z[i] << " ";
    }
    cout << "\ny      = ";
    for(int i = 0; i < nPts; ++i)
    {
        cout << setw(6) << y[i] << " ";
    }

    cout << "\nproject = ";
    for(int i = 0; i < nPts; ++i)
    {
        cout << setw(6) << u[i] << " ";
    }



    // Display the pointwise error
    cout << setprecision(1);
    cout << "\nerror  = ";
    Linf = 0;
    RMS = 0;
    for(int i = 0; i < int(GP.GetRows()); ++i)
    {
        long double exact = func(z[i], nPts, pointsType);
        long double error = exact - u[i];
        Linf = max(Linf, fabs(error));
        RMS += error*error;
        long double epsilon = 1e-2;
        if( fabs(exact) > epsilon )
        {
            error /= exact;
        }
        cout << setw(6) << error << " ";
    }
    RMS = sqrt(RMS) / int(GP.GetRows());
    cout << setprecision(6);
    cout << "\nLinf   = " << setw(6) << Linf;
    cout << "\nRMS    = " << setw(6) << RMS << endl;

    // Show the projection matrix
    cout << "\nI = " << endl;
    for(int i = 0; i < int(GP.GetRows()); ++i)
    {
        cout << "     ";
        for(int j = 0; j < int(GP.GetColumns()); ++j)
        {
            printf("% 5.3f  ", GP(i,j));
        }
        cout << "" << endl;
    }


      // /////////////////////////////////////////////////////////
     //                     Derivation                         //
    // /////////////////////////////////////////////////////////


    // Get the derivation matrix D
    // Note that, unlike I, D is N rows by N columns
    const Points<NekDouble>::MatrixSharedPtrType Dptr = points->GetD();
    const NekMatrix<NekDouble> & D = *Dptr;


    // Differentiate the data values in the y vector using the derivative matrix D
    Array<OneD, NekDouble> v(nPts);
    for(int i = 0; i < int(D.GetRows()); ++i)
    {
        v[i] = 0;
        for(int j = 0; j < int(D.GetColumns()); ++j)
        {
            v[i] += D(i,j) * y[j];
        }
    }


    // Display the derivative approximations
    cout << "\n\n\n              **** Differentiation ****" << setprecision(3);
    cout << "\n\nResults of differentiation with " << kPointsTypeStr[pointsType] << ":";
    cout << "\nderived = ";
    for(int i = 0; i < nPts; ++i)
    {
        cout << setw(6) << v[i] << " ";
    }

    // Determine the exact solutions
    cout << "\nexact   = ";
    for(int i = 0; i < nPts; ++i)
    {
        cout << setw(6) << derivativeFunction(z[i], nPts, pointsType) << " ";
    }

    // Display the pointwise error
    cout << setprecision(1);
    cout << "\nerror   = ";
    Linf = 0, RMS = 0;
    for(int i = 0; i < nPts; ++i)
    {
        long double exact = derivativeFunction(z[i], nPts, pointsType);
        long double error = exact - v[i];
        Linf = max(Linf, fabs(error));
        RMS += error*error;

        long double epsilon = 1e-2;
        if( fabs(exact) > epsilon )
        {
            error /= exact;
        }
        cout << setw(6) << error << " ";
    }

    // Display the global error
    RMS = sqrt(RMS) / I.GetRows();
    cout << setprecision(6);
    cout << "\nLinf    = " << setw(6) << Linf;
    cout << "\nRMS     = " << setw(6) << RMS << endl;


    // Show the derivation matrix
    cout << "\nD = " << endl;
    for(int i = 0; i < int(D.GetRows()); ++i)
    {
        cout << "     ";
        for(int j = 0; j < int(D.GetColumns()); ++j)
        {
            printf("% 5.3f  ", D(i,j));
        }
        cout << "" << endl;
    }




      // /////////////////////////////////////////////////////////
     //                     Integration                        //
    // /////////////////////////////////////////////////////////


    // Integrate the function by approximating the integral as a weighted
    // linear combination of function evaluations at the z[i] points: y[i] = f(z[i])
    long double numericalIntegration = 0.0;
    for(int i = 0; i < nPts; ++i)
    {
        numericalIntegration += w[i] * y[i] / integrandWeightFunction(z[i], pointsType);
    }


    // Display the integral approximation
    cout << "\n\n\n              **** Integration ****" << setprecision(6);
    cout << "\n\nResults of integration with " << kPointsTypeStr[pointsType] << ":";
    cout << "\nnumerical integral  = " << setw(12) << numericalIntegration;

    // Determine the exact solutions
    cout << "\nexact               = " << setw(12) << integrationFunction(nPts, pointsType);

    // Display the error
    cout << "\nerror               = ";
    {
        long double exact = integrationFunction(nPts, pointsType);
        long double error = exact - numericalIntegration;
        long double epsilon = 1e-2;
        if( fabs(exact) > epsilon )
        {
            error /= exact;
        }
        cout << setw(12) << error;
    }


    // Show the weights
    cout << "\nquadrature weights =      " << setprecision(3);
    for(int i = 0; i < nPts; ++i)
    {
        cout << setw(7) << w[i] << " ";
    }
    cout << endl;

}

