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
    cout << "| EXERCISE 2: GENERATION OF THE MASS MATRIX                  |" <<endl;
    cout << "--------------------------------------------------------------" <<endl;
    
    cout << "Assignment(a): Generate the mass matrix for the 7th order orthogonal spectral/hp basis" << endl;
    cout << "on a local standard element using 9th order Gaussian quadrature" << endl;
   
    // The local (straight-sided) quadrilateral element has the following vertices:
    // - Vertex A: (x1_A,x2_A) = (0,-1)
    // - Vertex B: (x1_A,x2_A) = (1,-1)
    // - Vertex C: (x1_A,x2_A) = (1,1)
    // - Vertex D: (x1_A,x2_A) = (0,1)
    NekDouble x1_A =  0.0;
    NekDouble x2_A = -1.0;
    NekDouble x1_B =  1.0;
    NekDouble x2_B = -1.0;
    NekDouble x1_C =  1.0;
    NekDouble x2_C =  1.0;
    NekDouble x1_D =  0.0;
    NekDouble x2_D =  1.0;      

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
    // psi_i(xi_j) = basisFunctionsDir[i*nQuadPointsDir+j]
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

    //==> Write your code here <==


    // Display the output
    cout << "Mass Matrix structure: " << endl;    
    cout<< "(Due to the ortogonality of the basis functions, the matrix should be diagonal)" <<endl;  
    int g;
    int h;             
    for(g = 0; g < nTotModes; g++)
    {
        for(h = 0; h < nTotModes; h++)
        {
            if(fabs(M(g,h))<1e-12)
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

