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
    cout << "| EXERCISE 1: INTEGRATION ON A 2 DIMENSIONAL LOCAL REGION |" <<endl;
    cout << "-----------------------------------------------------------" <<endl;

    cout << "Assignment(a): Integrate the function f(x1,x2) = x1^12 * x2^14 on a local" << endl; 
    cout << "quadrilateral element using Gaussian quadrature" << endl;

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
    // Your code can be based on Exercise 2(a) of the first tutorial. However, as we are
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
    cout << "Q1 = " << ", Q2 = ";
    cout << ": Error = "<< fabs(result-(1.0/195.0+1.0/420.0)) << endl;
}

