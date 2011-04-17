
#include <cstdio>
#include <cstdlib>

#include <MultiRegions/DisContField3DHomogeneous1D.h>

//#define TIMING
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
    MultiRegions::DisContField3DHomogeneous1DSharedPtr Exp,Fce;
    MultiRegions::ExpListSharedPtr DerExp1,DerExp2,DerExp3;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce;
    Array<OneD,NekDouble>  xc0,xc1,xc2;
    NekDouble  lambda,lz;
    NekDouble   st, cps = (double)CLOCKS_PER_SEC;

    if(argc != 2)
    {
        fprintf(stderr,"Usage: Helmholtz2D  meshfile\n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[1]);
    SpatialDomains::MeshGraph2D graph2D;
    graph2D.ReadGeometry(meshfile);
    graph2D.ReadExpansions(meshfile);
    //----------------------------------------------

    //----------------------------------------------
    // read the problem parameters from input file
    SpatialDomains::BoundaryConditions bcs(&graph2D);
    bcs.Read(meshfile);
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    lz      = bcs.GetParameter("Lz");
	bool useFFT = false;
    int nplanes = 8;
    const LibUtilities::PointsKey Pkey(nplanes,LibUtilities::eFourierEvenlySpaced);
    const LibUtilities::BasisKey Bkey(LibUtilities::eFourier,nplanes,Pkey);
    Exp = MemoryManager<MultiRegions::DisContField3DHomogeneous1D>::
        AllocateSharedPtr(Bkey,lz,useFFT,graph2D,bcs);
    //----------------------------------------------
    Timing("Read files and define exp ..");



    //----------------------------------------------
    // Print summary of solution details
    lambda  = bcs.GetParameter("Lambda");

    const SpatialDomains::ExpansionVector &expansions = graph2D.GetExpansions();
    LibUtilities::BasisKey bkey0 = expansions[0]->m_basisKeyVector[0];
    cout << "Solving 3D Helmholtz (Homogeneous in z-direction):"  << endl;
    cout << "         Lambda         : " << lambda << endl;
    cout << "         Lz             : " << lz << endl;
    cout << "         No. modes      : " << bkey0.GetNumModes() << endl;
    cout << "         No. hom. modes : " << Bkey.GetNumModes() << endl;
    cout << endl;
    //----------------------------------------------

    //----------------------------------------------
    // Set up coordinates of mesh for Forcing function evaluation
    coordim = Exp->GetCoordim(0);
    nq      = Exp->GetTotPoints();

    xc0 = Array<OneD,NekDouble>(nq,0.0);
    xc1 = Array<OneD,NekDouble>(nq,0.0);
    xc2 = Array<OneD,NekDouble>(nq,0.0);

    Exp->GetCoords(xc0,xc1,xc2);
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
    Fce = MemoryManager<MultiRegions::DisContField3DHomogeneous1D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //----------------------------------------------
    Timing("Define forcing ..");

    //----------------------------------------------
    // Helmholtz solution taking physical forcing
    Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), lambda);
    //----------------------------------------------

    Timing("Helmholtz Solve ..");

#ifdef TIMING
    for(i = 0; i < 100; ++i)
    {
        Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), lambda);
    }

    Timing("100 Helmholtz Solves:... ");
#endif

    //-----------------------------------------------
    // Backward Transform Solution to get solved values at
    Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys());
    //-----------------------------------------------
    Timing("Backard Transform ..");

    //-----------------------------------------------
    // Write solution to file
    string   out = meshfile.substr(0, meshfile.find_last_of(".")) + ".fld";
    std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef
                                                = Exp->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    for(i = 0; i < FieldDef.size(); ++i)
    {
        FieldDef[i]->m_fields.push_back("u");
        Exp->AppendFieldData(FieldDef[i], FieldData[i]);
    }
    graph2D.Write(out, FieldDef, FieldData);

    //-----------------------------------------------

    //-----------------------------------------------
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
        // Calculate error
        Fce->SetPhys(fce);
        Fce->SetPhysState(true);

        cout << "L infinity error:  " << Exp->Linf(Fce->GetPhys()) << endl;
        cout << "L 2 error  :       " << Exp->L2  (Fce->GetPhys()) << endl;
        //--------------------------------------------
    }

    Timing("Output ..");
    //----------------------------------------------
    return 0;
}

