#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ContField2D.h>

using namespace Nektar;

# define WITHSOLUTION 1

int main(int argc, char *argv[])
{
    cout << "======================================================" <<endl;
    cout << "===================== TUTORIAL 3 =====================" <<endl;
    cout << "======================================================" <<endl;

    cout << "-----------------------------------------------------------" <<endl;
    cout << "| EXERCISE 1:  A 2D HELMHOLTZ SOLVER                       |" <<endl;
    cout << "-----------------------------------------------------------" <<endl;

    cout << "Assignment(a): Solve the helmholtz equation with forcing function f(x1,x2) = " << endl;
    cout << "-(lambda+2*pi^2)*cos(pi*x1)*cos(pi*x2) on the mesh defined in the input" << endl;
    cout << "file QuadMesh.xml using Nektar++" << endl << endl;
    {
        // Information about the mesh and the boundary conditions are contained
        // in the input file QuadMesh.cpp
        // The following configuration of the problem is stored in this file
        //  - Domain of the problem: the bi-unit square xi1=[0,1], xi2=[0,2]
        //  - The mesh: 4 identical quadrilateral (square) elements
        //  - The expansion: 4th order C0 continuous modified spectral/hp expansion
        //  - Forcing function: f(x1,x2) = -(lambda+2*pi^2)*cos(pi*x1)*cos(pi*x2)
        //  - Boundary Conditions: Dirichlet boundary conditions g(xi1,xi2) = cos(pi*x1)*cos(pi*x2) 
        //    on the entire domain boundary
        //  - Exact solution: u(xi1,xi2) = cos(pi*x1)*cos(pi*x2) 
        //
        // This exercise contains of two parts:
        //
        // PART 1: Complete the Helmholtz solver
        // Solve the exercise using the proper Nektar++ classes and functions.
        // The major part is already implemented. However, some essential Nektar++ calls 
        // have been left out.
        // The task is to, using the Doxygen documentation of the Nektar++ code, find out
        // how the code should be properly completed.
        // Once this is finished, the solution should have been written to the file
        // Helmholtz2DSolution.pos. Open this file using the program Gmsh to see a plot
        // of the solution. 
        //
        // PART 2: Convergence study
        // Now you have a working Helmholtz solver, run it for the following configurations:
        //
        //  -  1 element , P = 2 (NUMMODES= P+1 = 3);
        // h-refinement:
        //  -  4 elements, P = 2;
        //  -  9 elements, P = 2;
        //  - 16 elements, P = 2;
        // P-refinement:
        // -  1 element , P = 4;
        // -  1 element , P = 6;
        // -  1 element , P = 8;
        // 
        // To do so, you will have to modify (or create a new) inputfile for every case.
        // Plot the approximation error (calculated in the L2 norm) for every case
        // in function of the number of degrees of freedom in a semi-log plot
        // If your implementation is correct, you should observe two different
        // convergence, depending on the strategy (i.e. h-refinement versus P-refinement)


        // The first step is to read in the input file to obtain: 
        // - the mesh
        // - the boundary conditions
        string fileName = "QuadMesh.xml";

        // First read the mesh. All relevant mesh information will be contained 
        // in an object of the class MeshGraph2D.
        SpatialDomains::MeshGraph2D mesh; 
        // Both the geometry and the expansion information should be read
        mesh.ReadGeometry(fileName);
        mesh.ReadExpansions(fileName);

        // Also read the boundary conditions. This information will be stored in
        // an object of the class BoundaryConditions
        SpatialDomains::BoundaryConditions boundaryConds(&mesh); 
        boundaryConds.Read(fileName);

        // The next step is to declare an object of a suitable Nektar++ class
        // which will represent the spectral/hp element expansion of the given problem.
        // For the problem considered, the suitable Nektar++ class can be selected as follows:
        //
        // The spectral/hp expansion
        // - involves multiple elements => a class of the MultiRegions library
        // - is two dimensional => a *2D* class
        // - is C0 continuous => a *Cont* class
        // - involves boundary conditions => a *Field* class
        //
        // => hence, the appropriate class to represent the spectral/hp expansion
        // for this problem is ContField2D
        //
        // Now, to construct the spectral/hp expansion, define an object from the 
        // class ContField2D. This object can be constructed based on the input mesh
        // and the boundary conditions.
        MultiRegions::ContField2DSharedPtr spectralHpExp =
            MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(mesh,boundaryConds);

        // Before actually starting to solve the Helmholtz problem, we need to evaluate the 
        // forcing function at the quadrature points of the (discrete) spectral/hp expansion
        int nTotQuadPoints = spectralHpExp->GetTotPoints();
        // Get the (local) coordinates of all the quadrature points
        Array<OneD,NekDouble> x1(nTotQuadPoints,0.0);
        Array<OneD,NekDouble> x2(nTotQuadPoints,0.0);
        Array<OneD,NekDouble> x3(nTotQuadPoints,0.0);
        spectralHpExp->GetCoords(x1,x2);

        // Declare an array which will hold the values of the forcing function evaluated at all
        // the quadrature points
        Array<OneD, NekDouble> forcingFunction(nTotQuadPoints);
        // The equation for the forcing function was defined in the input file, and has been stored
        // in the object of the BoundaryCondition class.
        // Extract this equation
        SpatialDomains::ConstForcingFunctionShPtr forcingFunctionEquation 
            = boundaryConds.GetForcingFunction(boundaryConds.GetVariable(0));
        // Evaluate it at the quadrature points
        for(int i = 0; i < nTotQuadPoints; ++i)
        {
            forcingFunction[i] = forcingFunctionEquation->Evaluate(x1[i],x2[i],x3[i]);
        }

        // We will now store the values of the forcing function as the values of a second
        // spectral/hp expansion evaluated at the quadrature points. To do so, 
        // a new multi-elemental spectral/hp expansion is first created (based on the original one 
        // such that they share the same mesh).
        MultiRegions::ContField2DSharedPtr forcingExp =
            MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(*spectralHpExp);
        // Next, the Array forcingFunction is copied into the data member m_phys of the new object
        // forcingExp. This can be done through the following function of the class
        // ContField2D:
        forcingExp->SetPhys(forcingFunction);

        // Now using the doxygen documentation, find the proper Nektar++ function
        // which solves the helmholtz problem.
        // Complete the following line of code with this function and use the correct
        // function arguments:
#if WITHSOLUTION
        NekDouble lambda =boundaryConds.GetParameter("Lambda");
        spectralHpExp->HelmSolve(forcingExp->GetPhys(),spectralHpExp->UpdateCoeffs(),lambda);
#endif

        // The function above solves the linear system corresponding to the Helmholtz problem
        // for the expansion coefficients u.
        // However, in order to calculate the value of the spectral/hp approximation at the 
        // quadrature points, we need a second function. 
        // Insert this function here:
#if WITHSOLUTION
        spectralHpExp->BwdTrans(spectralHpExp->GetCoeffs(),spectralHpExp->UpdatePhys());
#endif

        // The last step is to calculate the L2 error of the approximation. Therefore, we should also
        // know the exact solution. The equation of the exact solution was defined in the input file.
        // First, we will evaluate the exact solution at the quadrature points. This can be done in a similar 
        // fashion as for the forcing function:
        Array<OneD, NekDouble> exactSolution(nTotQuadPoints);
        SpatialDomains::ConstForcingFunctionShPtr exactSolutionEquation 
            = boundaryConds.GetExactSolution(boundaryConds.GetVariable(0));
        for(int i = 0; i < nTotQuadPoints; ++i)
        {
            exactSolution[i] = exactSolutionEquation->Evaluate(x1[i],x2[i],x3[i]);
        }       
        // We will now store the values of the exact solution as the values of another
        // spectral/hp expansion evaluated at the quadrature points. This can be analogue
        // as for the forcing function.
        MultiRegions::ContField2DSharedPtr exactSolutionExp =
            MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(*spectralHpExp);
        exactSolutionExp->SetPhys(exactSolution);

        // Now, call the right function to evaluate the L2 error of the spectral/hp approximation
        NekDouble error;
#if WITHSOLUTION
        error = spectralHpExp->L2(exactSolutionExp->GetPhys());
#endif

        // Display the output              
        cout << "Error = " << error << endl;
        cout << "Number of degrees of freedom = " << spectralHpExp->GetContNcoeffs() << endl << endl;

        // In addition, we will use one the Nektar++ output formats to visualise the output.
        // The solution is written to the file Helmholtz2DSolution.pos
        // which can be opened using the program Gmsh. Given the values of the coefficients of 
        // the expansion, this program then plots the expansion in a high-order fashion.  
        ofstream outfile("Helmholtz2DSolution.pos");      
        spectralHpExp->WriteToFile(outfile,eGmsh);
        outfile.close();   
        cout << "To see a plot of the solution, open the file Helmholtz2DSolution.pos" <<endl;
        cout << "using the program Gmsh."<<endl;
        cout << endl;    
    }   

    cout << "---------------------------------------------------------------" <<endl;
    cout << "| EXERCISE 2: DIFFERENTIATION ON A 2 DIMENSIONAL LOCAL REGION |" <<endl;
    cout << "---------------------------------------------------------------" <<endl;

    cout << "Assignment(a): Calculate the gradient of the function" << endl; 
    cout << "f(x1,x2) = x1^7 * x2^7 on a local quadrilateral element" << endl;
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
#if WITHSOLUTION
        int i,j,p,q,m,n;

        NekDouble x1;
        NekDouble x2;
        for(i = 0; i < nQuadPointsDir1; i++)
        {
            for(j = 0; j < nQuadPointsDir2; j++)
            {   
                m = i*nQuadPointsDir2+j;

                for(p = 0; p < nQuadPointsDir1; p++)
                {
                    x1 = x1_A * 0.25 * (1-quadZerosDir1[p]) * (1-quadZerosDir2[j]) + 
                        x1_B * 0.25 * (1+quadZerosDir1[p]) * (1-quadZerosDir2[j]) + 
                        x1_D * 0.25 * (1-quadZerosDir1[p]) * (1+quadZerosDir2[j]) + 
                        x1_C * 0.25 * (1+quadZerosDir1[p]) * (1+quadZerosDir2[j]); 

                    x2 = x2_A * 0.25 * (1-quadZerosDir1[p]) * (1-quadZerosDir2[j]) + 
                        x2_B * 0.25 * (1+quadZerosDir1[p]) * (1-quadZerosDir2[j]) + 
                        x2_D * 0.25 * (1-quadZerosDir1[p]) * (1+quadZerosDir2[j]) + 
                        x2_C * 0.25 * (1+quadZerosDir1[p]) * (1+quadZerosDir2[j]);

                    df_dxi1[m] += Ddir1(i,p)*pow(x1,7)*pow(x2,7);
                }

                for(q = 0; q < nQuadPointsDir2; q++)
                {
                    x1 = x1_A * 0.25 * (1-quadZerosDir1[i]) * (1-quadZerosDir2[q]) + 
                        x1_B * 0.25 * (1+quadZerosDir1[i]) * (1-quadZerosDir2[q]) + 
                        x1_D * 0.25 * (1-quadZerosDir1[i]) * (1+quadZerosDir2[q]) + 
                        x1_C * 0.25 * (1+quadZerosDir1[i]) * (1+quadZerosDir2[q]); 
                    
                    x2 = x2_A * 0.25 * (1-quadZerosDir1[i]) * (1-quadZerosDir2[q]) + 
                        x2_B * 0.25 * (1+quadZerosDir1[i]) * (1-quadZerosDir2[q]) + 
                        x2_D * 0.25 * (1-quadZerosDir1[i]) * (1+quadZerosDir2[q]) + 
                        x2_C * 0.25 * (1+quadZerosDir1[i]) * (1+quadZerosDir2[q]); 
                    
                    df_dxi2[m] += Ddir2(j,q)*pow(x1,7)*pow(x2,7);
                }
            }
        }
#endif

        // - Step 2: Calculate the metric terms of the transformation bewteen reference and local coordinates
        // First allocate four arrays which will store the values of dxi/dx at all the quadrature points:
        Array<OneD, NekDouble> dxi1_dx1(nTotQuadPoints,0.0);
        Array<OneD, NekDouble> dxi1_dx2(nTotQuadPoints,0.0);
        Array<OneD, NekDouble> dxi2_dx1(nTotQuadPoints,0.0);
        Array<OneD, NekDouble> dxi2_dx2(nTotQuadPoints,0.0);

        // fill these arrays
#if WITHSOLUTION
        NekDouble dx1_dxi1;
        NekDouble dx1_dxi2;
        NekDouble dx2_dxi1;
        NekDouble dx2_dxi2;
        NekDouble jacobian;
        for(i = 0; i < nQuadPointsDir1; i++)
        {
            for(j = 0; j < nQuadPointsDir2; j++)
            {
                m = i*nQuadPointsDir2+j;

                dx1_dxi1 = 0.25 * (1-quadZerosDir2[j]) * (x1_B-x1_A) + 0.25 * (1+quadZerosDir2[j]) * (x1_C-x1_D);
                dx2_dxi2 = 0.25 * (1-quadZerosDir1[i]) * (x2_D-x2_A) + 0.25 * (1+quadZerosDir1[i]) * (x2_C-x2_B);
                dx1_dxi2 = 0.25 * (1-quadZerosDir1[i]) * (x1_D-x1_A) + 0.25 * (1+quadZerosDir1[i]) * (x1_C-x1_B);
                dx2_dxi1 = 0.25 * (1-quadZerosDir2[j]) * (x2_B-x2_A) + 0.25 * (1+quadZerosDir2[j]) * (x2_C-x2_D);

                jacobian = dx1_dxi1*dx2_dxi2 - dx1_dxi2*dx2_dxi1;
                dxi1_dx1[m] =  dx2_dxi2/jacobian;
                dxi1_dx2[m] = -dx1_dxi2/jacobian;
                dxi2_dx1[m] = -dx2_dxi1/jacobian;
                dxi2_dx2[m] =  dx1_dxi1/jacobian;
            }
        }
#endif

        // - Step 3: apply the chain rule to calculate the gradient of f with respect to the local coordinates
        // First allocate two arrays which will store the values of df/dx at all the quadrature points:
        Array<OneD, NekDouble> df_dx1(nTotQuadPoints,0.0);
        Array<OneD, NekDouble> df_dx2(nTotQuadPoints,0.0);

        // Fill these arrays
#if WITHSOLUTION
        for(m = 0; m < nTotQuadPoints; m++)
        {
            df_dx1[m] = dxi1_dx1[m] * df_dxi1[m] + dxi2_dx1[m] * df_dxi2[m];
            df_dx2[m] = dxi1_dx2[m] * df_dxi1[m] + dxi2_dx2[m] * df_dxi2[m];
        }
#endif

        // Display the output
        NekDouble exactResultDir1; 
        NekDouble exactResultDir2;   
        for(i = 0; i < nQuadPointsDir1; i++)
        {
            for(j = 0; j < nQuadPointsDir2; j++)
            {  
                m = i*nQuadPointsDir2+j;
                
                x1 = x1_A * 0.25 * (1-quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                    x1_B * 0.25 * (1+quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                    x1_D * 0.25 * (1-quadZerosDir1[i]) * (1+quadZerosDir2[j]) + 
                    x1_C * 0.25 * (1+quadZerosDir1[i]) * (1+quadZerosDir2[j]); 
                
                x2 = x2_A * 0.25 * (1-quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                    x2_B * 0.25 * (1+quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                    x2_D * 0.25 * (1-quadZerosDir1[i]) * (1+quadZerosDir2[j]) + 
                    x2_C * 0.25 * (1+quadZerosDir1[i]) * (1+quadZerosDir2[j]); 

                exactResultDir1 = 7*pow(x1,6)*pow(x2,7);
                exactResultDir2 = 7*pow(x1,7)*pow(x2,6);

                cout << "Quadrature point " << m << ": Error = ";
                cout << fabs(df_dx1[m] - exactResultDir1);
                cout << " (df_dx1), ";
                cout << fabs(df_dx2[m] - exactResultDir2);
                cout << " (df_dx2)" << endl;
            }
        }
    }
    cout << endl;
}

