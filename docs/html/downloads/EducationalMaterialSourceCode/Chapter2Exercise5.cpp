#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ContField1D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    cout << "-------------------" <<endl;
    cout << "-- EXERCISE 3(a) --" <<endl;
    cout << "-------------------" <<endl;

    // We will not solve the exercise without using the Nektar++ library
    cout << "METHOD 1: Nektar++" << endl << endl;
    {
        int     i;
		 
        // The mesh is contained in the input file hemlholtz1D.xml
        // The first step is to read in this file
        stringstream fileName;
        fileName << "helmholtz1D.xml";

        string fileNameString = fileName.str();

        SpatialDomains::MeshGraph1D mesh; 

        // Both the geometry and the expansion information should be read
        mesh.ReadGeometry(fileNameString);
        mesh.ReadExpansions(fileNameString);
		
        // Also read the boundary conditions
        SpatialDomains::BoundaryConditions boundaryConds(&mesh); 
        boundaryConds.Read(fileNameString);

        // Print summary of solution details
        NekDouble  lambda = boundaryConds.GetParameter("Lambda");
        const SpatialDomains::CompositeVector domain = (mesh.GetDomain());
        cout << "Solving 1D Helmholtz:"  << endl; 
        cout << "         Lambda     : " << lambda << endl; 

        const SpatialDomains::ExpansionVector &expansions = mesh.GetExpansions();
        for(i = 0; i < domain.size(); ++i)
        {
            LibUtilities::BasisKey bkey = expansions[0]->m_BasisKeyVector[0];
            cout << "      Composite " << i << "   : " << endl;
            cout << "         Expansion  : " << LibUtilities::BasisTypeMap[bkey.GetBasisType()] << endl;
            cout << "         No. modes  : " << bkey.GetNumModes() << endl;
        }
        cout << endl;
		
        // Construct an object from the class ContField1D.
        // This is the class which represents a multi-elemental
        // continuous spectral/hp expansion with boundary conditions.
        // This object can be constructed based on the input mesh
        // and the boundary conditions.
        // Define Expansion 
        MultiRegions::ContField1DSharedPtr multiElementExp = 
            MemoryManager<MultiRegions::ContField1D>::AllocateSharedPtr(mesh,boundaryConds);
		
        // Set up coordinates of mesh for Forcing function evaluation
        int coordim = multiElementExp->GetCoordim(0);
        int nTotQuadPoints      = multiElementExp->GetTotPoints();
	    
        Array<OneD,NekDouble>  xc0,xc1,xc2; 
        xc0 = Array<OneD,NekDouble>(nTotQuadPoints);
        xc1 = Array<OneD,NekDouble>(nTotQuadPoints);
        xc2 = Array<OneD,NekDouble>(nTotQuadPoints);
	    
        switch(coordim)
        {
        case 1:
            multiElementExp->GetCoords(xc0);
            Vmath::Zero(nTotQuadPoints,&xc1[0],1);
            Vmath::Zero(nTotQuadPoints,&xc2[0],1);
            break;
        case 2:
            multiElementExp->GetCoords(xc0,xc1);
            Vmath::Zero(nTotQuadPoints,&xc2[0],1);
            break;
        case 3:
            multiElementExp->GetCoords(xc0,xc1,xc2);
            break;
        }
		
        // Define forcing function for first variable defined in file 
        Array<OneD,NekDouble>  fce(nTotQuadPoints); 
        SpatialDomains::ConstForcingFunctionShPtr ffunc 
            = boundaryConds.GetForcingFunction(boundaryConds.GetVariable(0));
        for(i = 0; i < nTotQuadPoints; ++i)
        {
            fce[i] = ffunc->Evaluate(xc0[i],xc1[i],xc2[i]);
        }
        //----------------------------------------------

        //----------------------------------------------
        // Setup expansion containing the  forcing function
        MultiRegions::ContField1DSharedPtr Fce = MemoryManager<MultiRegions::ContField1D>::AllocateSharedPtr(*multiElementExp);
        Fce->SetPhys(fce);
        //----------------------------------------------
	  
        //----------------------------------------------
        // Helmholtz solution taking physical forcing 
        multiElementExp->HelmSolve(Fce->GetPhys(),multiElementExp->UpdateCoeffs(),lambda);
        //----------------------------------------------
	    
        //----------------------------------------------
        // Backward Transform Solution to get solved values at 
        multiElementExp->BwdTrans(multiElementExp->GetCoeffs(),multiElementExp->UpdatePhys());
        //----------------------------------------------
	    
        //----------------------------------------------
        // Write solution 
        ofstream outfile("HelmholtzFile1D.dat");
        multiElementExp->WriteToFile(outfile);
        //----------------------------------------------
	    
        //----------------------------------------------
        // See if there is an exact solution, if so 
        // evaluate and plot errors
        SpatialDomains::ConstExactSolutionShPtr ex_sol =
            boundaryConds.GetExactSolution(boundaryConds.GetVariable(0));

        if(ex_sol)
        {
            //----------------------------------------------
            // evaluate exact solution 
            for(i = 0; i < nTotQuadPoints; ++i)
            {
                fce[i] = ex_sol->Evaluate(xc0[i],xc1[i],xc2[i]);
            }
            //----------------------------------------------

            //--------------------------------------------
            // Calculate L_inf error 
            Fce->SetPhys(fce);
            cout << "L infinity error: " << multiElementExp->Linf(Fce->GetPhys()) << endl;
            cout << "L 2 error:        " << multiElementExp->L2  (Fce->GetPhys()) << endl;
            //--------------------------------------------        
        }
        //----------------------------------------------
    }

    return 0;

}

