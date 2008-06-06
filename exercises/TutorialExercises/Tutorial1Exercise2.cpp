#include <cstdio>
#include <cstdlib>

#include <StdRegions/StdQuadExp.h>
#include <LibUtilities/Polylib/Polylib.h>
using namespace Nektar;

int main(int argc, char *argv[])
{
    cout << "======================================================" <<endl;
    cout << "===================== TUTORIAL 1 =====================" <<endl;
    cout << "======================================================" <<endl;

    cout << "--------------------------------------------------------------" <<endl;
    cout << "| EXERCISE 2: INTEGRATION ON A 2 DIMENSIONAL STANDARD REGION |" <<endl;
    cout << "--------------------------------------------------------------" <<endl;
    
    cout << "Assignment(a): Integrate the function f(xi1,xi2) = xi1^12 * xi2^14 on the standard" << endl; 
    cout << "quadrilateral using Gaussian quadrature" << endl;
    {
        // Specify the number of quadrature points in both directions
        int nQuadPointsDir1 = 7;
        int nQuadPointsDir2 = 7;

        // Specify the type of quadrature points in both directions 
        LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
        LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussLobattoLegendre;

        // Declare variables (of type Array) to hold the quadrature zeros and weights
        Array<OneD, NekDouble> quadZerosDir1(nQuadPointsDir1);
        Array<OneD, NekDouble> quadWeightsDir1(nQuadPointsDir1);
        Array<OneD, NekDouble> quadZerosDir2(nQuadPointsDir2);
        Array<OneD, NekDouble> quadWeightsDir2(nQuadPointsDir2);

        // Calculate the GLL-quadrature zeros and weights in both directions. This is done in 2 steps.
        // Step 1: Declare the PointsKeys which uniquely defines the quadrature points

        //==> Write your code here <==


        // Step 2: Using this key, the quadrature zeros and weights can now be retrieved through
        // the PointsManager

        //==> Write your code here <==


        // Now you have the quadrature zeros and weight, apply the Gaussian quadrature technique
        // to integrate the function f(xi1,xi2) = xi1^12 * xi2^14 on the standard quadrilateral.
        // To do so, write a (double) loop which performs the summation.
        //
        // Store the solution in the variable 'result'
        NekDouble result = 0.0;

        //==> Write your code here <==


        // Display the output
        cout << "Q1 = " << nQuadPointsDir1 << ", Q2 = " << nQuadPointsDir2 ;
        cout << ": Error = "<< fabs(result-4.0/195.0) << endl;
    }
    cout << endl;

    cout <<  "Assignment(b): Integrate the function f(xi1,xi2) = xi1^12 * xi2^14 on the standard" << endl;
    cout << "triangle using Gaussian quadrature" << endl;
    {
        // Specify the number of quadrature points in both directions
        int nQuadPointsDir1 = 7;
        int nQuadPointsDir2 = 7;

        // Specify the type of quadrature points in both directions 
        // Use Gauss-Lobatto-Legendre points in both directions
        LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
        LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussLobattoLegendre;

        // Declare variables (of type Array) to hold the quadrature zeros and weights
        Array<OneD, NekDouble> quadZerosDir1(nQuadPointsDir1);
        Array<OneD, NekDouble> quadWeightsDir1(nQuadPointsDir1);
        Array<OneD, NekDouble> quadZerosDir2(nQuadPointsDir2);
        Array<OneD, NekDouble> quadWeightsDir2(nQuadPointsDir2);

        // Calculate the GLL-quadrature zeros and weights in both directions. This is done in 2 steps.
        // Step 1: Declare the PointsKeys which uniquely defines the quadrature points

        //==> Write your code here <==


        // Step 2: Using this key, the quadrature zeros and weights can now be retrieved through
        // the PointsManager

        //==> Write your code here <==


        // Now we have the quadrature zeros and weight, so we can apply the Gaussian quadrature
        // to integrate the function f(xi1,xi2) = xi1^12 * xi2^14 on the standard triangle.
        // To do so, write a (double) loop which performs the summation.
        //
        // However, note the following:
        // - the quadrature zeros are expressed in the collapsed coordinates (eta1,eta2)
        // - for the evaluation of the integrand, we need the reference coordinates (xi1,xi2) 
        //   of the quadrature zeros.
        //
        // Therefore in your implementation, make sure to use the transformation between both 
        // coordinate systems before evaluation the integrand.
        //
        // The transformation of coordinates also introduces a Jacobian in the integral.
        // Make sure to properly deal with this.
        //
        // Store the solution in the variable 'result'
        NekDouble result = 0.0;

        //==> Write your code here <==


        // Display the output
        cout << "Q1 = " << nQuadPointsDir1 << ", Q2 = " << nQuadPointsDir2 ;
        cout << ": Error = "<< fabs(result-2.0/195.0) << endl;
    }
    cout << endl;   

    cout <<  "Assignment(c): Integrate the function f(xi1,xi2) = xi1^12 * xi2^14 on the standard" << endl;
    cout << "triangle using Gaussian quadrature with Gauss-Radau-Jacobi points" << endl;
    {
        // Repeat the previous assignment, but with the following changes:
        //
        // (1) Use Gauss-Radau-Jacobi(alpha=1,beta=0) quadrature points in 
        //     the xi2-direction. This has the following advantages:
        //     - The jacobian appears as weight term in the quadrature rule
        //       so there is no need to calculate it.
        //     - The Radau distribution of points does not included the 
        //       singular value at the collapsed vertex. (Although this does
        //       not affect the integration, it does affect numerical 
        //       differentiation)
        //     In Nektar++, these points correspond to the PointsType:
        //     LibUtilities::eGaussRadauMAlpha1Beta0
        //
        // (2) Set the quadrature order Q1 and Q2 to the minimal value
        //     required for exact integration. Calculate these value 
        //     analytically and check if the error is indeed equal to zero.

        //==> Write your code here <==

    }
}

