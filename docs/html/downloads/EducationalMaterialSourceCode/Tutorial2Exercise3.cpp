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

    cout << "--------------------------------------------------------------" <<endl;
    cout << "| EXERCISE 3: ELEMENTAL PROJECTION PROBLEM                   |" <<endl;
    cout << "--------------------------------------------------------------" <<endl;
    
    cout << "Assignment(a): Project the function f(x1,x2)= x1^6 * x2^6 defined on a local" << endl;
    cout << "quadrilateral element onto the spectral/hp element expansion defined on this element." << endl;
   
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

    // As outlined in the manual, this problem corresponds to solving the
    // following linear system: 
    //
    //     Mu=F
    // 
    // where: - M is the mass matrix
    //        - F is the vector who's i-th entry contains the inner product
    //          of the function f(x1,x2) with the i-th expansion mode.
    //
    // The solution to this problem (i.e. the value of the resulting spectral/hp
    // element expansion evaluated at the quadrature points) can be obtained in four
    // steps:
    // - Step 1: construct the mass matrix M
    //   This can be implemented based on Exercise 3(a) of the first tutorial. However, as we 
    //   are now integrating on a local element, the Jacobian of the coordinate transformation
    //   should be taken into account when applying Gaussian quadrature. 
    //
    // - Step 2: construct the right hand side vector F
    //   This also involves the integration on a local element.     
    //
    // - Step 3: solve the linear system to obtain the coefficients of the expansion.
    //
    // - Step 4: do a backward transformation of the expansion coefficients to obtain the values
    //   of the expansion evaluated at the quadrature points

    // Use the following properties:
    // - 6th order modified C0 expansion
    // - 8th order gaussian quadrature


    // SET UP ALL INFORMATION REQUIRED FOR THE GAUSSIAN QUADRATURE
    // Specify the number of quadrature points in both directions
    int nQuadPointsDir1 = 8;
    int nQuadPointsDir2 = 8;

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

    // SET UP ALL INFORMATION OF THE ORTHOGONAL BASIS FUNCTIONS
    // Specify the number of expansion modes in both directions
    int nModesDir1 = 7;
    int nModesDir2 = 7;

    // Specify the type of expansion in both directions
    LibUtilities::BasisType basisTypeDir1 = LibUtilities::eModified_A;
    LibUtilities::BasisType basisTypeDir2 = LibUtilities::eModified_A;

    // Calculate the values of the basis functions evaluated at the quadrature points
    // in both directions. This is done in 2 steps.
    // Step 1: Declare the BasisKey's which uniquely define the (one-dimensional) expansion bases.
    // As Nektar++ only implements discrete expansion bases (i.e. evaluated at the quadrature points), 
    // a BasisKey requires information about the quadrature (i.e. a PointsKey).
    const LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,nModesDir1,quadPointsKeyDir1);
    const LibUtilities::BasisKey basisKeyDir2(basisTypeDir2,nModesDir2,quadPointsKeyDir2);

    // Step 2: Using this key, the values of the expansion functions evaluated at the quadrature points
    // can now be retrieved through the BasisManager. 
    //   - first declare the variables
    Array<OneD, NekDouble> basisFunctionsDir1(nModesDir1*nQuadPointsDir1);
    Array<OneD, NekDouble> basisFunctionsDir2(nModesDir2*nQuadPointsDir2);
    //   - call the manager
    basisFunctionsDir1 = LibUtilities::BasisManager()[basisKeyDir1]->GetBdata();
    basisFunctionsDir2 = LibUtilities::BasisManager()[basisKeyDir2]->GetBdata();

    // The discrete basis functions are stored vector-wise rather than matrix-wise. This is
    // done in the following way:
    // 
    // psi_i(xi_j) = B[i*Q+j]
    //
    // where: - B is basisFunctionsDir1 or basisFunctionsDir2
    //        - Q is nQuadPointsDir1 or nQuadPointsDir2
    //
    // This corresponds to storing the matrix B[m][n] = psi_n(xi_m) in column major format


    // STEP 1: CONSTRUCT THE MASS MATRIX
    // First we will declare a variable M which will represent the mass matrix. This is
    // an object of the class NekMatrix and can be allocated as follows:
    int nTotModes = nModesDir1*nModesDir2;
    NekMatrix<NekDouble> M(nTotModes,nTotModes,0.0);
    // The elements of the mass matrix can now be accessed by the call M(i,j)

    // Now we have the all the ingredients required to fill the mass matrix.
    // In order to do so, write a series of loops and use Gaussian quadrature
    // to calculate each element. Don't forget to include the jacobian of the transformation

    //==> Write your code here <==


    // STEP 2: CONSTRUCT THE RIGHT HAND SIDE VECTOR
    // First we will declare a variable F which will represent the rhs vector. This is
    // an object of the class NekVector and can be allocated as follows:
    NekVector<NekDouble> F(nTotModes,0.0);
    // The elements of the rhs vector can now be accessed by the call F[i]

    // In order to fill the rhs vector, write a series of loops and use Gaussian quadrature
    // to calculate each element. Don't forget to 
    //   - transform the quadrature zeros to local coordinates when evaluating the function f(x1,x2)
    //   - include the jacobian of the transformation in the Gaussian quadrature

    //==> Write your code here <==


    // STEP 3: SOLVE THE LINEAR SYSTEM
    // Solve the linear system for the coefficients of the expansion
    NekVector<NekDouble> u_coeff(nTotModes,0.0);
    //!! (Uncomment this line) !!//   M.Invert();
    //!! (Uncomment this line) !!//   u_coeff = M*F;

    // STEP 4: BACKWARD TRANSFORMATION OF THE COEFFICIENTS
    // Do a backward transformation of the expansion coefficients to obtain the values
    // of the expansion evaluated at the quadrature points
    int nTotQuadPoints = nQuadPointsDir1*nQuadPointsDir2;
    Array<OneD, NekDouble> u_phys(nTotQuadPoints,0.0);

    //==> Write your code here <==


    // Display the output
    // Verify that the spectral/hp expansion approximates the function f(x1,x2)    
    NekDouble x1;
    NekDouble x2;   
    int g;
    int h;
    int z; 
    for(g = 0; g < nQuadPointsDir1; g++)
    {
        for(h = 0; h < nQuadPointsDir2; h++)
        {  
            z = g*nQuadPointsDir2+h;
                
            x1 = x1_A * 0.25 * (1-quadZerosDir1[g]) * (1-quadZerosDir2[h]) + 
                x1_B * 0.25 * (1+quadZerosDir1[g]) * (1-quadZerosDir2[h]) + 
                x1_D * 0.25 * (1-quadZerosDir1[g]) * (1+quadZerosDir2[h]) + 
                x1_C * 0.25 * (1+quadZerosDir1[g]) * (1+quadZerosDir2[h]); 
                
            x2 = x2_A * 0.25 * (1-quadZerosDir1[g]) * (1-quadZerosDir2[h]) + 
                x2_B * 0.25 * (1+quadZerosDir1[g]) * (1-quadZerosDir2[h]) + 
                x2_D * 0.25 * (1-quadZerosDir1[g]) * (1+quadZerosDir2[h]) + 
                x2_C * 0.25 * (1+quadZerosDir1[g]) * (1+quadZerosDir2[h]); 
                
            cout << "Quadrature point " << z << ": Error = ";
            cout << fabs(u_phys[z] - pow(x1,6)*pow(x2,6)) << endl;
        }
    }

    // Although a 6th order spectral/hp expansion is used to approximate a 
    // 6th order polynomial, the solution is not exact. Try to reason why this is
    // in this case.
}

