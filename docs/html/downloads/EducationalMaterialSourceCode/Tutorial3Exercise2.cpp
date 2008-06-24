#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ContField2D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    cout << "======================================================" <<endl;
    cout << "===================== TUTORIAL 3 =====================" <<endl;
    cout << "======================================================" <<endl;

    cout << "Assignment(a): Calculate the gradient of the function" << endl; 
    cout << "f(x1,x2) = x1^7 * x2^7 on a local quadrilateral element" << endl;

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

    // Specify the number of quadrature points in both directions
    int nQuadPointsDir1 = 8;
    int nQuadPointsDir2 = 8;

    // Specify the type of quadrature points in both directions 
    LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
    LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussLobattoLegendre;

    // Declare variables (of type Array) to hold the quadrature zeros in both directions
    Array<OneD, NekDouble> quadZerosDir1(nQuadPointsDir1);
    Array<OneD, NekDouble> quadZerosDir2(nQuadPointsDir2);

    // Calculate the GLL-quadrature zeros in both directions. This is done in 2 steps.
    // Step 1: Declare the PointsKeys which uniquely defines the quadrature points
    const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
    const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);

    // Step 2: Using this key, the quadrature zeros can now be retrieved through
    // the PointsManager
    quadZerosDir1   = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetZ();
    quadZerosDir2   = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetZ();

    // The next step is to calculate the differentiation matrices. 
    // First we will declare a variables Ddir1 and Ddir2 which will represent the differentiation 
    // matrices in both directions. 
    // These matrices are an object of the class NekMatrix and can be allocated as follows:
    NekMatrix<NekDouble> Ddir1(nQuadPointsDir1,nQuadPointsDir1);
    NekMatrix<NekDouble> Ddir2(nQuadPointsDir2,nQuadPointsDir2);

    // Rather than filling this matrix yourself, it can be retrieved through the PointsManager 
    // by the following call.
    Ddir1 = *((LibUtilities::PointsManager()[quadPointsKeyDir1])->GetD());
    Ddir2 = *((LibUtilities::PointsManager()[quadPointsKeyDir2])->GetD());
    // The elements of the differentiation matrix can now be accessed through the call Ddir(i,j)

    // The gradient can now be calculated using the following steps:
    // - Step 1: Calculate the derivatives of f with respect to the reference coordinates xi1 and xi2
    // - Step 2: Calculate the metric terms of the transformation bewteen reference and local coordinates
    // - Step 3: apply the chain rule to calculate the gradient of f with respect to the local coordinates

    // - Step 1: Calculate the derivatives of f with respect to the reference coordinates xi1 and xi2
    // First allocate two arrays which will store the values of df/dxi at all the quadrature points:
    int nTotQuadPoints = nQuadPointsDir1*nQuadPointsDir2;
    Array<OneD, NekDouble> df_dxi1(nTotQuadPoints,0.0);
    Array<OneD, NekDouble> df_dxi2(nTotQuadPoints,0.0);
    // The index m of those arrays can be related to the indices i and j of the one-dimensional quadrature 
    // points in the following way: m = i*nQuadPointsDir2+j

    // Fill the arrays df_dxi1 and df_dxi2

    //==> Write your code here <==


    // - Step 2: Calculate the metric terms of the transformation bewteen reference and local coordinates
    // First allocate four arrays which will store the values of dxi/dx at all the quadrature points:
    Array<OneD, NekDouble> dxi1_dx1(nTotQuadPoints,0.0);
    Array<OneD, NekDouble> dxi1_dx2(nTotQuadPoints,0.0);
    Array<OneD, NekDouble> dxi2_dx1(nTotQuadPoints,0.0);
    Array<OneD, NekDouble> dxi2_dx2(nTotQuadPoints,0.0);

    // fill these arrays

    //==> Write your code here <==


    // - Step 3: apply the chain rule to calculate the gradient of f with respect to the local coordinates
    // First allocate two arrays which will store the values of df/dx at all the quadrature points:
    Array<OneD, NekDouble> df_dx1(nTotQuadPoints,0.0);
    Array<OneD, NekDouble> df_dx2(nTotQuadPoints,0.0);

    // Fill these arrays

    //==> Write your code here <==


    // Display the output
    NekDouble exactResultDir1; 
    NekDouble exactResultDir2; 
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

            exactResultDir1 = 7*pow(x1,6)*pow(x2,7);
            exactResultDir2 = 7*pow(x1,7)*pow(x2,6);

            cout << "Quadrature point " << z << ": Error = ";
            cout << fabs(df_dx1[z] - exactResultDir1);
            cout << " (df_dx1), ";
            cout << fabs(df_dx2[z] - exactResultDir2);
            cout << " (df_dx2)" << endl;
        }
    }
}

