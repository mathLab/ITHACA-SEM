#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ContField2D.h>

using namespace Nektar;

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
    cout << "file: QuadMesh.xml using Nektar++" << endl << endl;

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
    NekDouble lambda = boundaryConds.GetParameter("Lambda");

    //==> Write your code here <==


    // The function above solves the linear system corresponding to the Helmholtz problem
    // for the expansion coefficients u.
    // However, in order to calculate the value of the spectral/hp approximation at the 
    // quadrature points, we need a second function. 
    // Insert this function here:

    //==> Write your code here <==


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

    //==> Write your code here <==


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

