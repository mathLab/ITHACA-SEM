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


    cout << "------------------------------------------------------------------" <<endl;
    cout << "| EXERCISE 2: DIFFERENTIATION ON A 1 DIMENSIONAL STANDARD REGION |" <<endl;
    cout << "------------------------------------------------------------------" <<endl;
    
    cout << "Assignment (a): Numerically evaluate the derivative of the function " << endl;
    cout << "f(xi) = xi^7 on the standard segment xi=[-1,1]" << endl;

    // Specify the number of quadrature points
    int nQuadPoints = 7;
        
    // Specify the type of quadrature points. This is done using the proper
    // Nektar++ syntax.
    LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;
        
    // Declare a variable (of type Array) to hold the quadrature zeros
    Array<OneD, NekDouble> quadZeros(nQuadPoints);
        
    // Calculate the quadrature zeros. This is done in 2 steps.
    // Step 1: Declare a PointsKey which uniquely defines the quadrature points
    const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);

    // Step 2: Using this key, the quadrature zeros can now be retrieved through
    // the PointsManager
    quadZeros   = (LibUtilities::PointsManager()[quadPointsKey])->GetZ();

    // The next step is to calculate the differentiation matrix. 
    // First we will declare a variable D which will represent the differentiation matrix. 
    // This is an object of the class NekMatrix and can be allocated as follows:
    NekMatrix<NekDouble> D(nQuadPoints,nQuadPoints);

    // Rather than filling this matrix yourself, it can be retrieved through the PointsManager 
    // by the following call.
    D = *((LibUtilities::PointsManager()[quadPointsKey])->GetD());
    // The elements of the differentiation matrix can now be accessed through the call D(i,j)


    // Now you have the differentiation matrix, numerically evaluate the derivative
    // of the function f(xi) = xi^7 at the quadrature points of the standard segment.
    // To do so, write a (double) loop which performs the necessary summation.
    //
    // Store the solution at every quadrature point in the Array 'result'
    Array<OneD, NekDouble> result(nQuadPoints,0.0);

    //==> Write your code here <==


    // Display the output
    for(int i = 0; i < nQuadPoints; i++)
    { 
        cout << "Quadrature point " << i << ": Error = ";
        cout << fabs(result[i] - 7*pow(quadZeros[i],6)) << endl;
    }

    // Now evaluate the derivative for Q = 7,8,9 quadrature points and verify
    // how many quadrature points are required to get an exact answer to numerical precision.
}

