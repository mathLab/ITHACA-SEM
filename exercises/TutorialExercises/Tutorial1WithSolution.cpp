#include <cstdio>
#include <cstdlib>

#include <StdRegions/StdQuadExp.h>
#include <LibUtilities/Polylib/Polylib.h>
using namespace Nektar;

# define WITHSOLUTION 1

int main(int argc, char *argv[])
{
    cout << "======================================================" <<endl;
    cout << "===================== TUTORIAL 1 =====================" <<endl;
    cout << "======================================================" <<endl;

    cout << "--------------------------------------------------------------" <<endl;
    cout << "| EXERCISE 1: INTEGRATION ON A 1 DIMENSIONAL STANDARD REGION |" <<endl;
    cout << "--------------------------------------------------------------" <<endl;
    
    cout << "Assignment (a): Integrate the function f(xi) = xi^12 on the standard" << endl;
    cout << "segment xi=[-1,1] using Gaussian quadrature" << endl;
    {
        // Specify the number of quadrature points
        int nQuadPoints = 4;
        
        // Specify the type of quadrature points. This is done using the proper
        // Nektar++ syntax.
        LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;
        
        // Declare variables (of type Array) to hold the quadrature zeros and weights
        Array<OneD, NekDouble> quadZeros(nQuadPoints);
        Array<OneD, NekDouble> quadWeights(nQuadPoints);
        
        // Calculate the quadrature zeros and weights. This is done in 2 steps.
        // Step 1: Declare a PointsKey which uniquely defines the quadrature points
        const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);

        // Step 2: Using this key, the quadrature zeros and weights can now be retrieved through
        // the PointsManager
        quadZeros   = (LibUtilities::PointsManager()[quadPointsKey])->GetZ();
        quadWeights = (LibUtilities::PointsManager()[quadPointsKey])->GetW();

        // Now you have the quadrature zeros and weight, apply the Gaussian quadrature technique
        // to integrate the function f(xi) = xi^12 on the standard segment xi=[-1,1].
        // To do so, write a loop which performs the summation.
        //
        // In C++, a loop is implemented as:
        //     for(i = min; i < max; i++)
        //     {
        //         // do something  
        //     }
        // 
        // The function f(xi) can be evaluated using the command 'pow(xi,12)' of the 
        // cmath standard library
        //
        // Store the solution in the variable 'result'
        NekDouble result = 0.0;

#if WITHSOLUTION
        for(int i = 0; i < nQuadPoints; i++)
        {
            result += pow(quadZeros[i],12) * quadWeights[i];
        }
#endif
        // Display the output
        NekDouble exactResult = 2.0/13.0;
        cout << "Q = " << nQuadPoints << ": Error = "<< fabs(result-exactResult) << endl;


        // Now evaluate the integral for a quadrature order of Q = Q_max
        // where Q_max is the number of quadrature points required for an exact evaluation of 
        // the integral (calculate this value analytically).
        // Check that the error should then be zero (up to numerical precision).
    }   
    cout << endl;

    cout << "Assignment (b): Integrate the function f(xi) = cos(xi) on the standard" << endl;
    cout << "segment xi=[-1,1] using Gaussian quadrature for 2<Q<9" << endl;
    {
        // Your code can be based on the previous exercise.
        //
        // The function f(xi) can be evaluated using the command 'cos(xi)' of the 
        // cmath standard library.
        //
        // Display the error for every value of Q and notice that the error 
        // converges exponentially (i.e. the error drops an order of magnitude with increasing Q)
#if WITHSOLUTION
        for(int nQuadPoints = 2; nQuadPoints< 9; nQuadPoints++)
        {        
            // Specify the type of quadrature points. This is done using the proper
            // Nektar++ syntax.
            LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;
            
            // Declare variables (of type Array) to hold the quadrature zeros and weights
            Array<OneD, NekDouble> quadZeros(nQuadPoints);
            Array<OneD, NekDouble> quadWeights(nQuadPoints);
            
            // Calculate the quadrature zeros and weights. This is done in 2 steps.
            // Step 1: Declare a PointsKey which uniquely defines the quadrature points
            const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);
            
            // Step 2: Using this key, the quadrature zeros and weights can now be retrieved through
            // the PointsManager
            quadZeros   = (LibUtilities::PointsManager()[quadPointsKey])->GetZ();
            quadWeights = (LibUtilities::PointsManager()[quadPointsKey])->GetW();

            // Now you have the quadrature zeros and weight, apply the Gaussian quadrature technique
            // to integrate the function f(xi) = xi^12 on the standard segment xi=[-1,1].
            // To do so, write a loop which performs the summation.
            //
            // In C++, a loop is implemented as:
            //     for(i = min; i < max; i++)
            //     {
            //         // do something  
            //     }
            // 
            // The function f(xi) can be evaluated using the command 'sin(M_PI*xi)' of the 
            // cmath standard library
            //
            // Store the solution in the variable 'result'
            NekDouble result = 0.0;
            
            for(int i = 0; i < nQuadPoints; i++)
            {
                result += cos(quadZeros[i]) * quadWeights[i];
            }

            // Display the output
            NekDouble exactResult = (sin(1.0)-sin(-1.0));
            cout << "Q = " << nQuadPoints << ": Error = "<< fabs(result-exactResult) << endl;
        }
#endif

    }   
    cout << endl;

    cout << "------------------------------------------------------------------" <<endl;
    cout << "| EXERCISE 2: DIFFERENTIATION ON A 1 DIMENSIONAL STANDARD REGION |" <<endl;
    cout << "------------------------------------------------------------------" <<endl;
    
    cout << "Assignment (a): Numerically evaluate the derivative of the function " << endl;
    cout << "f(xi) = xi^7 on the standard segment xi=[-1,1]" << endl;
    {
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

#if WITHSOLUTION
        for(int i = 0; i < nQuadPoints; i++)
        { 
            for(int j = 0; j < nQuadPoints; j++)
            {
                result[i] += D(i,j) * pow(quadZeros[j],7);
            }
        }
#endif
        // Display the output
        for(int i = 0; i < nQuadPoints; i++)
        { 
            cout << "Quadrature point " << i << ": Error = ";
            cout << fabs(result[i] - 7*pow(quadZeros[i],6)) << endl;
        }

        // Now evaluate the derivative for Q = 7,8,9 quadrature points and verify
        // how many quadrature points are required to get an exact answer to numerical precision.
    }  
    cout << endl;

    cout << "------------------------------------------------------------------" <<endl;
    cout << "| EXERCISE 3: INTERPOLATION ON A 1 DIMENSIONAL STANDARD REGION   |" <<endl;
    cout << "------------------------------------------------------------------" <<endl;
    
    cout << "Assignment (a): Interpolate the function f(xi) = xi^12 evaluated at the Q=13 Gauss-Lobatto-Legendre zeros" << endl;
    cout << "to Q=10 equidistant points in the standard interval xi=[-1,1]." << endl;
    {
        // Specify the number of Gauss-Lobatto-Legendre quadrature points and the number of equidistant
        // interpolation points
        int nQuadPoints = 13;
        int nInterpolPoints = 10;
        
        // Specify the type of quadrature points and interpolation points. This is done using the proper
        // Nektar++ syntax.
        LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;
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

#if WITHSOLUTION
        for(int i = 0; i < nInterpolPoints; i++)
        { 
            for(int j = 0; j < nQuadPoints; j++)
            {
                result[i] += I(i,j) * pow(quadZeros[j],12);
            }
        }
#endif
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
#if WITHSOLUTION
        // Declare variables (of type Array) to hold the quadrature weights
        Array<OneD, NekDouble> interpolWeights(nInterpolPoints);
        
        // The quadrature weights can now be obtained through
        // the PointsManager
        interpolWeights = (LibUtilities::PointsManager()[interpolPointsKey])->GetW();
        
        NekDouble InterpolPointsResult = 0.0;  
        for(int i = 0; i < nInterpolPoints; i++)
        {
            InterpolPointsResult += pow(interpolZeros[i],12) * interpolWeights[i];
        }
        
        // Display the output
        NekDouble exactResult = 2.0/13.0;
        cout << "Q = " << nInterpolPoints << " equidistant quadrature points: Integration error = ";
        cout << fabs(InterpolPointsResult-exactResult) << endl;
#endif
    }
    cout << endl;
}

