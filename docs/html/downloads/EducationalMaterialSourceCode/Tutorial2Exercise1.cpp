#include <cstdio>
#include <cstdlib>

#include <StdRegions/StdQuadExp.h>
#include <LibUtilities/Polylib/Polylib.h>
using namespace Nektar;

int main(int argc, char *argv[])
{
    cout << "======================================================" <<endl;
    cout << "===================== TUTORIAL 2 =====================" <<endl;
    cout << "======================================================" <<endl;

    cout << "-----------------------------------------------------------" <<endl;
    cout << "| EXERCISE 1: INTEGRATION ON 2 DIMENSIONAL ELEMENTS       |" <<endl;
    cout << "-----------------------------------------------------------" <<endl;
    
    cout << "Assignment(a): Integrate the function f(xi1,xi2) = xi1^12 * xi2^14 on the standard" << endl; 
    cout << "quadrilateral element using Gaussian quadrature" << endl;
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
        const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
        const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);

        // Step 2: Using this key, the quadrature zeros and weights can now be retrieved through
        // the PointsManager
        quadZerosDir1   = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetZ();
        quadWeightsDir1 = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetW();

        quadZerosDir2   = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetZ();
        quadWeightsDir2 = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetW();

        // Now you have the quadrature zeros and weight, apply the Gaussian quadrature technique
        // to integrate the function f(xi1,xi2) = xi1^12 * xi2^14 on the standard quadrilateral.
        // To do so, write a (double) loop which performs the summation.
        //
        // Store the solution in the variable 'result'
        NekDouble result = 0.0;

        //==> Write your code here <==


        // Display the output
        NekDouble exactResult = 4.0/195.0;
        cout << "Q1 = " << nQuadPointsDir1 << ", Q2 = " << nQuadPointsDir2 ;
        cout << ": Error = "<< fabs(result-exactResult) << endl;
    }
    cout << endl;

    cout << "Assignment(b): Integrate the function f(x1,x2) = x1^12 * x2^14 on a local" << endl; 
    cout << "quadrilateral element using Gaussian quadrature" << endl;
    {
        // The local (straight-sided) quadrilateral element has the following vertices:
        // - Vertex A: (x1_A,x2_A) = (0,-1)
        // - Vertex B: (x1_A,x2_A) = (1,-1)
        // - Vertex C: (x1_A,x2_A) = (1,1)
        // - Vertex D: (x1_A,x2_A) = (0,0)
        NekDouble x1_A =  0.0;
        NekDouble x2_A = -1.0;
        NekDouble x1_B =  1.0;
        NekDouble x2_B = -1.0;
        NekDouble x1_C =  1.0;
        NekDouble x2_C =  1.0;
        NekDouble x1_D =  0.0;
        NekDouble x2_D =  0.0;       

        // Integrate the function f(x1,x2) = x1^12 * x2^14 on a local quadrilateral element (defined
        // above). Use 7th order Gauss-Lobatto-Legendre quadrature in both direction.
        //
        // Your code can be based on the previous exercise. However, as we are
        // calculating the integral of a function defined on a local element rather than on 
        // a reference element, we have to take into account the geometry of the element.
        // Therefore, the implementation should be altered in two ways:
        // (1) The quadrature zeros should be transformed to local coordinates to evaluate the
        //     integrand f(x1,x2)
        // (2) Take into account the Jacobian of the transformation between local and reference
        //     coordinates when evaluating the integral. (Evaluate the expression for the Jacobian
        //     analytically rather than using numerical differentiation)
        //
        // Store the solution in the variable 'result'
        NekDouble result = 0.0;

        //==> Write your code here <==



        // Display the output
        NekDouble exactResult = 1.0/195.0+1.0/420.0;
        cout << "Error = "<< fabs(result-exactResult) << endl;
    }
    cout << endl;
}

