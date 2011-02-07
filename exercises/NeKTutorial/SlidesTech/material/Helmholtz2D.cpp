#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ContField2D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    MultiRegions::ContField2DSharedPtr      Exp, Fce;
    int                                     i, nq,  coordim;
    Array<OneD,NekDouble>                   fce;
    Array<OneD,NekDouble>                   xc0, xc1, xc2;
    NekDouble                               lambda;

    if(argc != 2 )
    {
        fprintf(stderr,"Usage: Helmholtz2D meshfile\n");
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
    // Print summary of solution details
    lambda = bcs.GetParameter("Lambda");
    const SpatialDomains::ExpansionVector &expansions = graph2D.GetExpansions();
    LibUtilities::BasisKey bkey0 = expansions[0]->m_basisKeyVector[0];
    cout << "Solving 2D Helmholtz: "  << endl;
    cout << "         Lambda     : " << lambda << endl;
    cout << "         No. modes  : " << bkey0.GetNumModes() << endl;
    cout << endl;
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    Exp = MemoryManager<MultiRegions::ContField2D>::
        AllocateSharedPtr(graph2D,bcs);
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
    case 2:
        Exp->GetCoords(xc0,xc1);
        break;
    case 3:
        Exp->GetCoords(xc0,xc1,xc2);
        break;
    default:
        ASSERTL0(false,"Coordim not valid");
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
    // Setup expansion containing the forcing function
    Fce = MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //----------------------------------------------

    //----------------------------------------------
    // Helmholtz solution taking physical forcing
    Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateContCoeffs(), lambda, true);
    //----------------------------------------------

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
        FieldDef[i]->m_fields.push_back("u");
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

