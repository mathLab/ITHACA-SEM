#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ContField2D.h>

using namespace Nektar;

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

int main(int argc, char *argv[])
{
    MultiRegions::ContField2DSharedPtr Exp,Fce;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce;
    Array<OneD,NekDouble>  xc0,xc1,xc2;
    NekDouble  lambda;
    NekDouble   st, cps = (double)CLOCKS_PER_SEC;

    if( (argc != 2) && (argc != 4) )
    {
        fprintf(stderr,"Usage: Helmholtz2D meshfile "
                    " [GlobalOptimizationFile ElementalOptimizationFile]\n");
        exit(1);
    }

    //----------------------------------------------
    // Load the ELEMENTAL optimization parameters if they
    // have been given as argument
    // (use ElementalOptimisationParameters.xml as an example)
    if( argc == 4 )
    {
        string eloptfile(argv[3]);
        NekOptimize::LoadElementalOptimizationParameters(eloptfile);
    }
    //----------------------------------------------

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
    // Print summary of solution details
    lambda = bcs.GetParameter("Lambda");
    const SpatialDomains::ExpansionVector &expansions = graph2D.GetExpansions();
    LibUtilities::BasisKey bkey0 = expansions[0]->m_BasisKeyVector[0];
    cout << "Solving 2D Helmholtz:"  << endl;
    cout << "         Lambda     : " << lambda << endl;
    cout << "         No. modes  : " << bkey0.GetNumModes() << endl;
    cout << endl;
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    int bc_val = 0;
    // Use this definition in the constructor below if you want a
    // standard static condensation solver. Note default is mulit level. 
    MultiRegions::GlobalSysSolnType staticcond = MultiRegions::eDirectStaticCond;
    // Use this definition in the constructor below if you want a
    // full direct matrix  solver. Note default is mulit level. 
    MultiRegions::GlobalSysSolnType fulldirect = MultiRegions::eDirectFullMatrix;
    // Use this definition in the constructor below if you want a
    // multi level static condensation inversion. 
    MultiRegions::GlobalSysSolnType multilevel = MultiRegions::eDirectMultiLevelStaticCond;

    Exp = MemoryManager<MultiRegions::ContField2D>::
        AllocateSharedPtr(graph2D,bcs,bc_val,multilevel);
    //----------------------------------------------

    //----------------------------------------------
    // Load the global optimization parameters
    // as specified in the (optional) input file
    // (use helmholtz2D.xml as an example)
    if( argc == 4 )
    {
        string globoptfile(argv[2]);
        Exp->ReadGlobalOptimizationParameters(globoptfile);
    }
    //----------------------------------------------

    Timing("Read files and define exp ..");

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
    SpatialDomains::ConstForcingFunctionShPtr ffunc
        = bcs.GetForcingFunction(bcs.GetVariable(0));
    for(i = 0; i < nq; ++i)
    {
        fce[i] = ffunc->Evaluate(xc0[i],xc1[i],xc2[i]);
    }
    //----------------------------------------------

    //----------------------------------------------
    // Setup expansion containing the  forcing function
    Fce = MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //----------------------------------------------
    Timing("Define forcing ..");

    //----------------------------------------------
    // Helmholtz solution taking physical forcing
    Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateContCoeffs(), lambda, true);
    //----------------------------------------------
    Timing("Helmholtz Solve ..");

#ifdef TIMING
    for(i = 0; i < 1000; ++i)
    {
        Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateContCoeffs(), lambda, true);
    }

    Timing("1000 Helmholtz Solves:... ");
#endif

    //----------------------------------------------
    // Backward Transform Solution to get solved values
    Exp->BwdTrans(Exp->GetContCoeffs(), Exp->UpdatePhys(), true);
    //----------------------------------------------

    //-----------------------------------------------
    // Write solution to file
    string   out(strtok(argv[1],"."));
    string   endfile(".fld");
    out += endfile;
    std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef 
                                                = Exp->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    Exp->GlobalToLocal(Exp->GetContCoeffs(),Exp->UpdateCoeffs());
    for(i = 0; i < FieldDef.size(); ++i)
    {
        FieldDef[i]->m_Fields.push_back("u");
        Exp->AppendFieldData(FieldDef[i], FieldData[i]);
    }
    graph2D.Write(out, FieldDef, FieldData);

    //-----------------------------------------------

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
        Fce->SetPhysState(true);


        cout << "L infinity error: " << Exp->Linf(Fce->GetPhys()) << endl;
        cout << "L 2 error:        " << Exp->L2  (Fce->GetPhys()) << endl;
        cout << "H 1 error:        " << Exp->H1  (Fce->GetPhys()) << endl;
        //--------------------------------------------
    }
    //----------------------------------------------

    return 0;
}

