#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ContField1D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    MultiRegions::ContField1DSharedPtr Exp,Fce;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 
    NekDouble  lambda;

    if(argc != 3)
    {
        fprintf(stderr,"Usage: Helmholtz1D  meshfile boundaryfile\n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[1]);
    SpatialDomains::MeshGraph1D graph1D; 
    graph1D.ReadGeometry(meshfile);
    graph1D.ReadExpansions(meshfile);
    //----------------------------------------------

    //----------------------------------------------
    // read the problem parameters from input file
    string bcfile(argv[2]);
    SpatialDomains::BoundaryConditions bcs(&graph1D); 
    bcs.Read(bcfile);
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    lambda = bcs.GetParameter("Lambda");
    const SpatialDomains::CompositeVector domain = (graph1D.GetDomain());
    cout << "Solving 1D Helmholtz:"  << endl; 
    cout << "         Lambda     : " << lambda << endl; 
    for(i = 0; i < domain.size(); ++i)
    {
        LibUtilities::BasisKey bkey = graph1D.GetBasisKey(domain[i],0);
        cout << "      Composite " << i << "   : " << endl;
        cout << "         Expansion  : " << LibUtilities::BasisTypeMap[bkey.GetBasisType()] << endl;
        cout << "         No. modes  : " << bkey.GetNumModes() << endl;
    }
    cout << endl;
    //----------------------------------------------
   
    //----------------------------------------------
    // Define Expansion 
    Exp = MemoryManager<MultiRegions::ContField1D>::
        AllocateSharedPtr(graph1D,bcs);
    //----------------------------------------------
    
    //----------------------------------------------
    // Set up coordinates of mesh for Forcing function evaluation
    coordim = Exp->GetCoordim(0);
    nq      = Exp->GetPointsTot();
    
    xc0 = Array<OneD,NekDouble>(nq);
    xc1 = Array<OneD,NekDouble>(nq);
    xc2 = Array<OneD,NekDouble>(nq);
    
    switch(coordim)
    {
    case 1:
        Exp->GetCoords(xc0);
        Vmath::Zero(nq,&xc1[0],1);
        Vmath::Zero(nq,&xc2[0],1);
        break;
    case 2:
        Exp->GetCoords(xc0,xc1);
        Vmath::Zero(nq,&xc2[0],1);
        break;
    case 3:
        Exp->GetCoords(xc0,xc1,xc2);
        break;
    }
    //----------------------------------------------
    
    //----------------------------------------------
    // Define forcing function for first variable defined in file 
    fce = Array<OneD,NekDouble>(nq);
    SpatialDomains::ConstForcingFunctionShPtr ffunc 
        = bcs.GetForcingFunction(bcs.GetVariable(0));
    for(i = 0; i < nq; ++i)
    {
        fce[i] = ffunc->Evaluate(xc0[i],xc1[i],xc2[i]);
    }
    //----------------------------------------------

    //----------------------------------------------
    // Setup expansion containing the  forcing function
    Fce = MemoryManager<MultiRegions::ContField1D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //----------------------------------------------
  
    //----------------------------------------------
    // Helmholtz solution taking physical forcing 
    Exp->HelmSolve(*Fce, lambda);
    //----------------------------------------------
    
    //----------------------------------------------
    // Backward Transform Solution to get solved values at 
    Exp->BwdTrans(*Exp);
    //----------------------------------------------
    
    //----------------------------------------------
    // Write solution 
    ofstream outfile("HelmholtzFile1D.dat");
    Exp->WriteToFile(outfile);
    //----------------------------------------------
    
    //----------------------------------------------
    // See if there is an exact solution, if so 
    // evaluate and plot errors
    SpatialDomains::ConstExactSolutionShPtr ex_sol =
        bcs.GetExactSolution(bcs.GetVariable(0));

    if(ex_sol)
    {
        //----------------------------------------------
        // evaluate exact solution 
        for(i = 0; i < nq; ++i)
        {
            fce[i] = ex_sol->Evaluate(xc0[i],xc1[i],xc2[i]);
        }
        //----------------------------------------------

        //--------------------------------------------
        // Calculate L_inf error 
        Fce->SetPhys(fce);
        cout << "L infinity error: " << Exp->Linf(*Fce) << endl;
        cout << "L 2 error:        " << Exp->L2  (*Fce) << endl;
        //--------------------------------------------        
    }
    //----------------------------------------------

    return 0;
}

