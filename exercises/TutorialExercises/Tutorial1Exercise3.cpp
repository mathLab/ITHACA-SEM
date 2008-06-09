#include <cstdio>
#include <cstdlib>

#include <StdRegions/StdQuadExp.h>
#include <LibUtilities/Polylib/Polylib.h>

int main(int argc, char *argv[])
{
    cout << "======================================================" <<endl;
    cout << "===================== TUTORIAL 1 =====================" <<endl;
    cout << "======================================================" <<endl;

    cout << "------------------------------------------------------------------" <<endl;
    cout << "| EXERCISE 3: INTERPOLATION ON A 1 DIMENSIONAL STANDARD REGION   |" <<endl;
    cout << "------------------------------------------------------------------" <<endl;
    
    cout << "Assignment (a): Interpolate the function f(xi) = xi^12 evaluated at the Q=13 Gauss-Lobatto-Legendre zeros" << endl;
    cout << "to Q=10 equidistant points in the standard interval xi=[-1,1]." << endl;

    // Specify the number of Gauss-Lobatto-Legendre quadrature points and the number of equidistant
    // interpolation points
    int nQuadPoints = 13;
    int nInterpolPoints = 10;
        
    // Specify the type of quadrature points and interpolation points. This is done using the proper
    // Nektar++ syntax.
    LibUtilities::PointsType quadPointsType     = LibUtilities::eGaussLobattoLegendre;
    LibUtilities::PointsType interpolPointsType = LibUtilities::ePolyEvenlySpaced;
        
    // Declare a variable (of type Array) to hold the quadrature zeros and interpolation zeros.
    Array<OneD, NekDouble> quadZeros(nQuadPoints);
    Array<OneD, NekDouble> interpolZeros(nInterpolPoints);
        
    // Calculate the quadrature zeros and interpolation points. This is done in 2 steps.
    // Step 1: Declare a PointsKey which uniquely defines the quadrature points
    const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);
    const LibUtilities::PointsKey interpolPointsKey(nInterpolPoints, interpolPointsType);

    // Step 2: Using this key, the quadrature and interpolation zeros can now be retrieved through
    // the PointsManager
    quadZeros       = (LibUtilities::PointsManager()[quadPointsKey])->GetZ();
    interpolZeros   = (LibUtilities::PointsManager()[interpolPointsKey])->GetZ();


    // The next step is to calculate the interpolation matrix. 
    // First we will declare a variable I which will represent the interpolation matrix. 
    // This is an object of the class NekMatrix and can be allocated as follows:
    NekMatrix<NekDouble> I(nInterpolPoints,nQuadPoints);

    // Rather than filling this matrix yourself, it can be retrieved through the PointsManager 
    // by the following call.
    I = *((LibUtilities::PointsManager()[quadPointsKey])->GetI(interpolPointsKey));
    // The elements of the differentiation matrix can now be accessed through the call I(i,j)


    // Now you have the interpolation matrix, numerically evaluate the value of 
    // of the function f(xi) = xi^12 at the equidistant interpolation points within the standard segment.
    // To do so, write a loop which performs the necessary summation.
    //
    // Store the solution at every interpolation point in the Array 'result'
    Array<OneD, NekDouble> result(nInterpolPoints,0.0);

    //==> Write your code here <==


    // Display the output
    for(int i = 0; i < nInterpolPoints; i++)
    { 
        cout << "Interpolation point " << i << ": Error = ";
        cout << fabs(result[i] - pow(interpolZeros[i],12)) << endl;
    }


    // Now you have the values of the function f at Q=10 equidistant point,
    // calculate the quadrature weights associated with these points and numerically
    // evaluate the integral of f on the standard interval.
    // Verify that for Q=10 equidistant points, the integral is not approximated exactly.
    // Look back to exercise 1(a) to check how many Gauss-Lobatto-Legendre points
    // were required to exactly integrate the function f.

    //==> Write your code here <==
}

