#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ContField2D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    cout << "=====================================================" <<endl;
    cout << "================ EXERCISES CHAPTER 4 ================" <<endl;
    cout << "=====================================================" <<endl;

    cout << "-------------------" <<endl;
    cout << "--- EXERCISE 5 ----" <<endl;
    cout << "-------------------" <<endl;

    // We will not solve this exercise without using the Nektar++ library

    {
        for(int p = 4; p < 11; p++)
        {
            // The mesh and boundary condition information are
            // contained in the input file Chapter4Exercise4_Quad_Pp.xml
            // The first step is to read in this file
            stringstream fileName;
            fileName << "Chapter4Exercise4_Quad_P" << p << ".xml";

            string fileNameString = fileName.str();

            SpatialDomains::MeshGraph2D mesh; 
            // Both the geometry and th expansion information should be read
            mesh.ReadGeometry(fileNameString);
            mesh.ReadExpansions(fileNameString);
            // Also read the boundary conditions
            SpatialDomains::BoundaryConditions boundaryConds(&mesh); 
            boundaryConds.Read(fileNameString);

            // Construct an object from the class ContField2D.
            // This is the class which represents a multi-elemental
            // continuous spectral/hp expansion with boundary conditions.
            // This object can be constructed based on the input mesh
            // and the boundary conditions.
            MultiRegions::ContField2DSharedPtr multiElementExp =
                MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(mesh,boundaryConds);

            // Evaluate the forcing function at the quadrature points
            int nTotQuadPoints = multiElementExp->GetTotPoints();
            Array<OneD,NekDouble> x1(nTotQuadPoints,0.0);
            Array<OneD,NekDouble> x2(nTotQuadPoints,0.0);
            Array<OneD,NekDouble> x3(nTotQuadPoints,0.0);
            multiElementExp->GetCoords(x1,x2);

            Array<OneD, NekDouble> forcingFunction(nTotQuadPoints);
            // Read the forcing function equation as defined in the input file
            SpatialDomains::ConstForcingFunctionShPtr forcingFunctionEquation 
                = boundaryConds.GetForcingFunction(boundaryConds.GetVariable(0));
        
            for(int i = 0; i < nTotQuadPoints; ++i)
            {
                forcingFunction[i] = forcingFunctionEquation->Evaluate(x1[i],x2[i],x3[i]);
            }

            // Store the forcing function as the physical values of an
            // object of the class ContExpList2D
            MultiRegions::ContField2DSharedPtr forcingExp =
                MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(*multiElementExp);
            forcingExp->SetPhys(forcingFunction);

            // Do the projection to obtain the coefficients of the expansion
            // The result is stored in the data member m_contCoeffs of the ContExpList2D 
            // object multiElementExp.
            multiElementExp->FwdTrans(forcingExp->GetPhys(),multiElementExp->UpdateCoeffs());

            // Perform a backward transformation to obtain the solution at the quadrature points
            // The result is stored in the data member m_phys of the ContExpList2D 
            // object multiElementExp.
            multiElementExp->BwdTrans(multiElementExp->GetCoeffs(),multiElementExp->UpdatePhys());

            // Calculate the error
            NekDouble error = multiElementExp->L2(forcingExp->GetPhys());

            // Display the output              
            cout << "P = " << p << " => Error = " << error << endl;
            // In addition, we will use one the Nektar++ output formats to visualise the output.
            // The solution is written to the files Chapter4Exercise5_Quad_Pi.pos
            // which can be opened using the program Gmsh. Given the values of the coefficients of 
            // the expansion, this program then plots the expansion in a high-order fashion.  
            stringstream outfileName;
            outfileName << "Chapter4Exercise5_Quad_P" << p << ".pos";
            ofstream outfile((outfileName.str()).data());      
            multiElementExp->WriteToFile(outfile,eGmsh);
        }
        cout << "To see a plot of the solution, open the files Chapter4Exercise5_Quad_Pi.pos" <<endl;
        cout << "using the program Gmsh."<<endl;
        cout << endl;       
    }


}

