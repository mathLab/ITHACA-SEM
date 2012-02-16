#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
//#include <SpatialDomains/MeshPartition.h>
#include <MultiRegions/DisContField3D.h>

#define TIMING

#ifdef TIMING
#include <time.h>
#define Timing(s) \
 fprintf(stdout,"%s Took %g seconds\n",s,(clock()-st)/cps); \
 st = clock();
#else
#define Timing(s) \
 /* Nothing */
#endif

using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    LibUtilities::CommSharedPtr vComm = vSession->GetComm();
    string meshfile(argv[1]);

    MultiRegions::DisContField3DSharedPtr Exp, Fce;
    MultiRegions::ExpListSharedPtr DerExp1, DerExp2, DerExp3;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 
    NekDouble  lambda;
    double   st, cps = (double)CLOCKS_PER_SEC;

    if(argc != 3)
    {
        fprintf(stderr,"Usage: Helmholtz3D  meshfile boundaryfile\n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graph3D = MemoryManager<SpatialDomains::MeshGraph3D>::AllocateSharedPtr(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    lambda = vSession->GetParameter("Lambda");
    const SpatialDomains::ExpansionVector &expansions = graph3D->GetExpansions();
    cout << "Solving 3D Helmholtz:"  << endl;
    cout << "         Lambda     : " << lambda << endl; 
#if 0 
    for(i = 0; i < expansions.size(); ++i)
    {
        LibUtilities::BasisKey bkey = graph2D.GetBasisKey(expansions[i],0);
        cout << "      Element " << i << "   : " << endl;
        cout << "         Expansion  : " << LibUtilities::BasisTypeMap[bkey.GetBasisType()] << endl;
        cout << "         No. modes  : " << bkey.GetNumModes() << endl;
    }
    cout << endl;
#endif
    //----------------------------------------------
   
    //----------------------------------------------
    // Define Expansion 
    Exp = MemoryManager<MultiRegions::DisContField3D>::
        AllocateSharedPtr(vSession,graph3D,vSession->GetVariable(0));
    //----------------------------------------------
    Timing("Read files and define exp ..");
    
    //----------------------------------------------
    // Set up coordinates of mesh for Forcing function evaluation
    coordim = Exp->GetCoordim(0);
    nq      = Exp->GetPointsTot();
    
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
    LibUtilities::EquationSharedPtr ffunc = vSession->GetFunction("Forcing", 0);

    ffunc->Evaluate3Array(xc0, xc1, xc2, fce);

    //----------------------------------------------


    //----------------------------------------------
    // Setup expansion containing the  forcing function
    Fce = MemoryManager<MultiRegions::DisContField3D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //----------------------------------------------
    Timing("Define forcing ..");
  
    //----------------------------------------------
    // Helmholtz solution taking physical forcing 
    Exp->HelmSolve(*Fce, lambda);
    //----------------------------------------------
    
    Timing("Helmholtz Solve ..");

#if 0
    for(i = 0; i < 100; ++i)
    {
        Exp->HelmSolve(*Fce, lambda);
    }
    
    Timing("100 Helmholtz Solves:... ");
#endif 

    //----------------------------------------------
    // Backward Transform Solution to get solved values at 
    Exp->BwdTrans(*Exp);
    //----------------------------------------------
    Timing("Backard Transform ..");
    
    //----------------------------------------------
    // Write solution 
    ofstream outfile("UDGHelmholtzFile3D.dat");
    Exp->WriteToFile(outfile);
    //----------------------------------------------
    
    //----------------------------------------------
    // See if there is an exact solution, if so 
    // evaluate and plot errors
    LibUtilities::EquationSharedPtr ex_sol =
                            vSession->GetFunction("ExactSolution", 0);

    if(ex_sol)
    {
        //----------------------------------------------
        // evaluate exact solution 

       ex_sol->Evaluate3Array(xc0, xc1, xc2,  fce);

        //----------------------------------------------

        //--------------------------------------------
        // Calculate L_inf error 
        Fce->SetPhys(fce);
        Fce->SetPhysState(true);


        cout << "L infinity error: " << Exp->Linf(*Fce) << endl;
        cout << "L 2 error:        " << Exp->L2  (*Fce) << endl;
        //--------------------------------------------        
    }
    
    Timing("Output ..");
    //----------------------------------------------        
    return 0;
}

