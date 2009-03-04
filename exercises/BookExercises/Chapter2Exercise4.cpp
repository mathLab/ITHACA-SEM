#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ContField1D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
	cout << "=====================================================" <<endl;
    cout << "================ EXERCISES CHAPTER 2 ================" <<endl;
    cout << "=====================================================" <<endl;

    cout << "-------------------" <<endl;
    cout << "-- EXERCISE 4 a --" <<endl;
    cout << "-------------------" <<endl;

	// we will not solve this exercise without using the Nektar++ library
	cout << "METHOD 1:  Nektar++" << endl << endl;
	{
		int i;
		// The mesh is contained in the input file Chapter1Exercise4a.xml
        // The first step is to read in this file
        stringstream fileName;
        fileName << "Chapter2Exercise4a.xml";


		string fileNameString = fileName.str();

		SpatialDomains::MeshGraph1D mesh; 

		// Both the geometry and the expansion information should be read
        mesh.ReadGeometry(fileNameString);
        mesh.ReadExpansions(fileNameString);

		// Construct an object from the class ContExpList1D.
        // This is the class which represents a multi-elemental
        // This object can be constructed based on the input mesh
        MultiRegions::ContExpList1DSharedPtr multiElementExp =
            MemoryManager<MultiRegions::ContExpList1D>::AllocateSharedPtr(mesh);

		// Evaluate the forcing function at the quadrature points
        int nTotQuadPoints = multiElementExp->GetTotPoints();
        Array<OneD,NekDouble> x(nTotQuadPoints);
        multiElementExp->GetCoords(x);

		// Evaluate the forcing function
		Array<OneD, NekDouble> forcingFunction(nTotQuadPoints);
		for(int i = 0; i < nTotQuadPoints; ++i)
        {
           forcingFunction[i] = sin(x[i]);
		}

		// Store the forcing function as the physical values of an
        // object of the class ContExpList1D
        MultiRegions::ContExpList1DSharedPtr forcingExp =
            MemoryManager<MultiRegions::ContExpList1D>::AllocateSharedPtr(*multiElementExp);
        forcingExp->SetPhys(forcingFunction);

		// Do the projection to obtain the coefficients of the expansion
        // The result is stored in the data member m_contCoeffs of the ContExpList1D 
        // object multiElementExp.
        multiElementExp->FwdTrans(forcingExp->GetPhys(),multiElementExp->UpdateCoeffs());
		
		// Perform a backward transformation to obtain the solution at the quadrature points
        // The result is stored in the data member m_phys of the ContExpList1D 
        // object multiElementExp.
        multiElementExp->BwdTrans(multiElementExp->GetCoeffs(),multiElementExp->UpdatePhys());
		
		for(i = 0; i < nTotQuadPoints; i++)
        {               
            cout << "Quadrature point " << i << ":       Result = ";
            cout << (multiElementExp->GetPhys())[i]<< endl;
            cout << "                     Exact Result = "<< sin(x[i]) << endl; 
			cout << endl;
        }

	}

	cout << "-------------------" <<endl;
    cout << "-- EXERCISE 4 b --" <<endl;
    cout << "-------------------" <<endl;

	cout << "METHOD 1:  Nektar++" << endl << endl;
	{
		int i;
		// The mesh is contained in the input file Chapter1Exercise4b.xml
        // The first step is to read in this file
        stringstream fileName;
        fileName << "Chapter2Exercise4b.xml";

		string fileNameString = fileName.str();

		SpatialDomains::MeshGraph1D mesh; 

		// Both the geometry and the expansion information should be read
        mesh.ReadGeometry(fileNameString);
        mesh.ReadExpansions(fileNameString);

		// Also read the boundary conditions
        SpatialDomains::BoundaryConditions boundaryConds(&mesh); 
        boundaryConds.Read(fileNameString);

		// Construct an object from the class ContField1D.
        // This is the class which represents a multi-elemental
        // continuous spectral/hp expansion with boundary conditions.
        // This object can be constructed based on the input mesh
        // and the boundary conditions.
        MultiRegions::ContField1DSharedPtr multiElementExp =
            MemoryManager<MultiRegions::ContField1D>::AllocateSharedPtr(mesh,boundaryConds);

		// Evaluate the forcing function at the quadrature points
        int nTotQuadPoints = multiElementExp->GetTotPoints();
        Array<OneD,NekDouble> x(nTotQuadPoints);
        multiElementExp->GetCoords(x);

		Array<OneD, NekDouble> forcingFunction(nTotQuadPoints);

		// Read the forcing function equation as defined in the input file
        SpatialDomains::ConstForcingFunctionShPtr forcingFunctionEquation 
			  = boundaryConds.GetForcingFunction(boundaryConds.GetVariable(0));

		for(i = 0; i < nTotQuadPoints; i++)
        {
            forcingFunction[i] = forcingFunctionEquation->Evaluate(x[i]);
		}

		// Store the forcing function as the physical values of an
        // object of the class ContField1D
        MultiRegions::ContField1DSharedPtr forcingExp =
            MemoryManager<MultiRegions::ContField1D>::AllocateSharedPtr(*multiElementExp);
        forcingExp->SetPhys(forcingFunction);

		// Do the projection to obtain the coefficients of the expansion
        // The result is stored in the data member m_contCoeffs of the ContExpList1D 
        // object multiElementExp.
        multiElementExp->FwdTrans(forcingExp->GetPhys(),multiElementExp->UpdateCoeffs());
		
		// Perform a backward transformation to obtain the solution at the quadrature points
        // The result is stored in the data member m_phys of the ContField1D 
        // object multiElementExp.
        multiElementExp->BwdTrans(multiElementExp->GetCoeffs(),multiElementExp->UpdatePhys());
		
		// Getting the exact solution
		SpatialDomains::ConstExactSolutionShPtr exSol =
			boundaryConds.GetExactSolution(boundaryConds.GetVariable(0));

		// Evaluate the exact solution at the quarature points
		Array<OneD, NekDouble> exactSolution(nTotQuadPoints);
		for(i = 0; i < nTotQuadPoints; i++)
        {
            exactSolution[i] = exSol->Evaluate(x[i]);
        }
		
		for(i = 0; i < nTotQuadPoints; i++)
        {               
            cout << "Quadrature point " << i << ":       Result = ";
            cout << (multiElementExp->GetPhys())[i]<< endl;
            cout << "                     Exact Result = "<< exactSolution[i] << endl; 
			cout << endl;
        }
	
	}

        return 0;
}
