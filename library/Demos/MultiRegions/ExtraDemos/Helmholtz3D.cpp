#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ContField3D.h>

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    MultiRegions::ContField3DSharedPtr Exp, Fce;
    MultiRegions::ExpListSharedPtr DerExp1, DerExp2, DerExp3;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 
    NekDouble  lambda;

    if(argc != 3)
    {
        fprintf(stderr,"Usage: Helmholtz3D  meshfile boundaryfile\n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[1]);
    SpatialDomains::MeshGraph3D graph3D;
    graph3D.ReadGeometry(meshfile);
    graph3D.ReadExpansions(meshfile);
    //----------------------------------------------

    //----------------------------------------------
    // read the problem parameters from input file
    string bcfile(argv[2]);
    SpatialDomains::BoundaryConditions bcs(&graph3D);
    bcs.Read(bcfile);
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    lambda = bcs.GetParameter("Lambda");
    cout << "Solving 3D Helmholtz:"  << endl;
    cout << "         Lambda     : " << lambda << endl; 
    cout << endl;
    //----------------------------------------------
   
    //----------------------------------------------
    // Define Expansion 
    Exp = MemoryManager<MultiRegions::ContField3D>::AllocateSharedPtr(graph3D,bcs);
    //----------------------------------------------

    
    //----------------------------------------------
    // Set up coordinates of mesh for Forcing function evaluation
    coordim = Exp->GetCoordim(0);
    nq      = Exp->GetTotPoints();
    
    xc0 = Array<OneD,NekDouble>(nq,0.0);
    xc1 = Array<OneD,NekDouble>(nq,0.0);
    xc2 = Array<OneD,NekDouble>(nq,0.0);
    
    switch(coordim)
    {
    case 1:
        Exp->GetCoords(xc0);
        break;
    case 2:
        Exp->GetCoords(xc0,xc1);
        break;
    case 3:
        Exp->GetCoords(xc0,xc1,xc2);
        break;
    }
    //----------------------------------------------
    
    //----------------------------------------------
    // Define forcing function for first variable defined in file 
    fce = Array<OneD,NekDouble>(nq);
    SpatialDomains::ConstForcingFunctionShPtr ffunc = bcs.GetForcingFunction(bcs.GetVariable(0));

    ffunc->Evaluate(xc0,xc1,xc2,fce);
    //----------------------------------------------

    //----------------------------------------------
    // Setup expansion containing the  forcing function
    Fce = MemoryManager<MultiRegions::ContField3D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //----------------------------------------------
  
    //----------------------------------------------
    // Helmholtz solution taking physical forcing 
    Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), lambda);
    //----------------------------------------------
    
    //----------------------------------------------
    // Backward Transform Solution to get solved values at 
    Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys());
    //----------------------------------------------
    
    //----------------------------------------------
    // Write solution 
    ofstream outfile("HelmholtzFile3D.pos");
    Exp->WriteToFile(outfile,eGmsh);
    //----------------------------------------------
    
    //----------------------------------------------
    // See if there is an exact solution, if so 
    // evaluate and plot errors
    SpatialDomains::ConstExactSolutionShPtr ex_sol = bcs.GetExactSolution(bcs.GetVariable(0));

    if(ex_sol)
    {
        //----------------------------------------------
        // evaluate exact solution 
        ex_sol->Evaluate(xc0,xc1,xc2,fce);

        //----------------------------------------------

        //--------------------------------------------
        // Calculate L_inf error 
        Fce->SetPhys(fce);
        Fce->SetPhysState(true);


        cout << "L infinity error: " << Exp->Linf(Fce->GetPhys()) << endl;
        cout << "L 2 error:        " << Exp->L2  (Fce->GetPhys()) << endl;
        //--------------------------------------------        
    }
    //----------------------------------------------        
        return 0;
}

