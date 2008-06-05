#include <cstdio>
#include <cstdlib>
#include <sys/time.h>

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
        cout << "Q = " << nQuadPoints << ": Error = "<< fabs(result-2.0/13.0) << endl;


        // Modify the implementation such that the integral is evaluated for the following range
        // of quadrature orders: Q = 3,...,Q_max
        // where Q_max is the number of quadrature points required for an exact evaluation of 
        // the integral (calculate this value analytically).
        // Display the error for every value of Q and notice the following two things:
        // - the error converges exponentially (the error drops an order of magnitude with increasing Q)
        // - if Q = Q_max, the error is zero (up to machine accuracy)
    }   
    cout << endl;

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
#if WITHSOLUTION
        const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
        const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);
#endif

        // Step 2: Using this key, the quadrature zeros and weights can now be retrieved through
        // the PointsManager
#if WITHSOLUTION
        quadZerosDir1   = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetZ();
        quadWeightsDir1 = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetW();

        quadZerosDir2   = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetZ();
        quadWeightsDir2 = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetW();
#endif

        // Now you have the quadrature zeros and weight, apply the Gaussian quadrature technique
        // to integrate the function f(xi1,xi2) = xi1^12 * xi2^14 on the standard quadrilateral.
        // To do so, write a (double) loop which performs the summation.
        //
        // Store the solution in the variable 'result'
        NekDouble result = 0.0;

#if WITHSOLUTION
        int i,j;
        for(i = 0; i < nQuadPointsDir1; i++)
        { 
            for(j = 0; j < nQuadPointsDir2; j++)
            {
                result += pow(quadZerosDir1[i],12) *  pow(quadZerosDir2[j],14) * 
                    quadWeightsDir1[i] * quadWeightsDir2[j];
            }
        }
#endif
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
#if WITHSOLUTION
        const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
        const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);
#endif

        // Step 2: Using this key, the quadrature zeros and weights can now be retrieved through
        // the PointsManager
#if WITHSOLUTION
        quadZerosDir1   = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetZ();
        quadWeightsDir1 = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetW();

        quadZerosDir2   = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetZ();
        quadWeightsDir2 = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetW();
#endif

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

#if WITHSOLUTION
        NekDouble xi1;
        NekDouble xi2;
        NekDouble jacobian;
        int i,j;
        for(i = 0; i < nQuadPointsDir1; i++)
        { 
            for(j = 0; j < nQuadPointsDir2; j++)
            {
                xi1      = 0.5 * (1 + quadZerosDir1[i]) * (1 - quadZerosDir2[j]) - 1;
                xi2      = quadZerosDir2[j];
                jacobian = 0.5 * (1-quadZerosDir2[j]);
                result += pow(xi1,12) *  pow(xi2,14) * jacobian *
                    quadWeightsDir1[i] * quadWeightsDir2[j];
            }
        }
#endif
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

#if WITHSOLUTION
        // Specify the number of quadrature points in both directions
        int nQuadPointsDir1 = 8;
        int nQuadPointsDir2 = 8;

        // Specify the type of quadrature points in both directions 
        // Use Gauss-Lobatto-Legendre points in both directions
        LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
        LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussRadauMAlpha1Beta0;

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
        // Store the solution in the variable 'result'
        NekDouble result = 0.0;

        NekDouble xi1;
        NekDouble xi2;
        int i,j;
        for(i = 0; i < nQuadPointsDir1; i++)
        { 
            for(j = 0; j < nQuadPointsDir2; j++)
            {
                xi1      = 0.5 * (1 + quadZerosDir1[i]) * (1 - quadZerosDir2[j]) - 1;
                xi2      = quadZerosDir2[j];
                result += pow(xi1,12) *  pow(xi2,14) *
                    quadWeightsDir1[i] * 0.5 * quadWeightsDir2[j];
            }
        }

        // Display the output
        cout << "Q1 = " << nQuadPointsDir1 << ", Q2 = " << nQuadPointsDir2 ;
        cout << ": Error = "<< fabs(result-2.0/195.0) << endl;
#endif
    }
    cout << endl;

    cout << "--------------------------------------------------------------" <<endl;
    cout << "| EXERCISE 3: GENERATION OF THE MASS MATRIX                  |" <<endl;
    cout << "--------------------------------------------------------------" <<endl;
    
    cout << "Assignment(a): Generate the mass matrix for the 7th order orthogonal spectral/hp basis" << endl;
    cout << "on the quadrilateral standard element using 9th order Gaussian quadrature" << endl;
    {      
        // SET UP ALL INFORMATION REQUIRED FOR THE GAUSSIAN QUADRATURE
        // Specify the number of quadrature points in both directions
        int nQuadPointsDir1 = 9;
        int nQuadPointsDir2 = 9;

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
        LibUtilities::BasisType basisTypeDir1 = LibUtilities::eOrtho_A;
        LibUtilities::BasisType basisTypeDir2 = LibUtilities::eOrtho_A;

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


        // SET UP THE MASS MATRIX
        // First we will declare a variable M which will represent the mass matrix. This is
        // an object of the class NekMatrix and can be allocated as follows:
        int nTotModes = nModesDir1*nModesDir2;
        NekMatrix<NekDouble> M(nTotModes,nTotModes,0.0);
        // The elements of the mass matrix can now be accessed by the call M(i,j)

        // Now we have the all the ingredients required to fill the mass matrix.
        // In order to do so, write a series of loops and use Gaussian quadrature
        // to calculate each element.

#if WITHSOLUTION
        int p,q,r,s,i,j;
        int n,m;
        for(p = 0; p < nModesDir1; p++)
        {
            for(q = 0; q < nModesDir2; q++)
            {
                for(r = 0; r < nModesDir1; r++)
                {
                    for(s = 0; s < nModesDir2; s++)
                    {
                        for(i = 0; i < nQuadPointsDir1; i++)
                        {
                            for(j = 0; j < nQuadPointsDir2; j++)
                            {
                                n = p*nModesDir2+q;
                                m = r*nModesDir2+s;

                                M(n,m) +=  
                                    basisFunctionsDir1[p*nQuadPointsDir1+i] *
                                    basisFunctionsDir2[q*nQuadPointsDir2+j] *
                                    basisFunctionsDir1[r*nQuadPointsDir1+i] *
                                    basisFunctionsDir2[s*nQuadPointsDir2+j] *
                                    quadWeightsDir1[i] * quadWeightsDir2[j];
                            }
                        }
                    }
                }
            }
        }
#endif

        // Display the output
        cout << "Mass Matrix structure: " << endl;    
        cout<< "(Due to the ortogonality of the basis functions, the matrix should be diagonal)" <<endl;               
        for(i = 0; i < nTotModes; i++)
        {
            for(j = 0; j < nTotModes; j++)
            {
                if(fabs(M(i,j))<1e-12)
                {
                    cout<<"-";
                }
                else
                {
                    cout<<"o";
                }   
                    
            }
            cout<<endl;
        } 
    }
    cout << endl;
}

