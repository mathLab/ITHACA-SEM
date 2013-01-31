#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/DisContField3D.h>
#include <SpatialDomains/MeshGraph3D.h>

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
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    LibUtilities::CommSharedPtr vComm = vSession->GetComm();
    string meshfile(argv[1]);

    MultiRegions::DisContField3DSharedPtr Exp, Fce;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 
    StdRegions::ConstFactorMap factors;
    double   st, cps = (double)CLOCKS_PER_SEC;

    if(argc != 2)
    {
        fprintf(stderr,"Usage: Helmholtz3D  meshfile\n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graph3D = 
        MemoryManager<SpatialDomains::MeshGraph3D>::AllocateSharedPtr(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    factors[StdRegions::eFactorLambda] = vSession->GetParameter("Lambda");
    factors[StdRegions::eFactorTau] = 1.0;
    const SpatialDomains::ExpansionMap &expansions = graph3D->GetExpansions();
    LibUtilities::BasisKey bkey0
                            = expansions.begin()->second->m_basisKeyVector[0];

    if (vComm->GetRank() == 0)
    {
        cout << "Solving 3D Helmholtz:"  << endl;
        cout << "         Lambda     : " << factors[StdRegions::eFactorLambda] << endl;
        cout << "         No. modes  : " << bkey0.GetNumModes() << endl;
        cout << endl;
    }
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
    LibUtilities::EquationSharedPtr ffunc = vSession->GetFunction("Forcing", 0);

    ffunc->Evaluate(xc0, xc1, xc2, fce);

    //----------------------------------------------


    //----------------------------------------------
    // Setup expansion containing the  forcing function
    Fce = MemoryManager<MultiRegions::DisContField3D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //----------------------------------------------
    Timing("Define forcing ..");
  
    //----------------------------------------------
    // Helmholtz solution taking physical forcing 
    Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), NullFlagList, factors);
    //----------------------------------------------
    
    Timing("Helmholtz Solve ..");

#if 0
    for(i = 0; i < 100; ++i)
    {
        Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), NullFlagList, factors);
    }
    
    Timing("100 Helmholtz Solves:... ");
#endif 

    //----------------------------------------------
    // Backward Transform Solution to get solved values at 
    Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys());
    //----------------------------------------------
    Timing("Backward Transform ..");
    
    //-----------------------------------------------
    // Write solution to file
    string out = vSession->GetSessionName();
    if (vComm->GetSize() > 1)
    {
        out += "_P" + boost::lexical_cast<string>(vComm->GetRank());
    }
    out += ".fld";
    std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef
                                                = Exp->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    for(i = 0; i < FieldDef.size(); ++i)
    {
        FieldDef[i]->m_fields.push_back("u");
        Exp->AppendFieldData(FieldDef[i], FieldData[i]);
    }
    graph3D->Write(out, FieldDef, FieldData);
    //-----------------------------------------------
    
    //----------------------------------------------
    // See if there is an exact solution, if so 
    // evaluate and plot errors
    LibUtilities::EquationSharedPtr ex_sol =
        vSession->GetFunction("ExactSolution", 0);
    
    if(ex_sol)
    {
        //----------------------------------------------
        // evaluate exact solution 
        ex_sol->Evaluate(xc0, xc1, xc2,  fce);

        //----------------------------------------------

        //--------------------------------------------
        // Calculate L_inf error 
        Fce->SetPhys(fce);
        Fce->SetPhysState(true);

        NekDouble vLinfError = Exp->Linf(Fce->GetPhys());
        NekDouble vL2Error   = Exp->L2  (Fce->GetPhys());
        NekDouble vH1Error   = Exp->H1  (Fce->GetPhys());

        if (vComm->GetRank() == 0)
        {
            cout << "L infinity error: " << vLinfError << endl;
            cout << "L 2 error       : " << vL2Error   << endl;
            cout << "H 1 error       : " << vH1Error   << endl;
        }
        //--------------------------------------------        
    }
    
    Timing("Output ..");

    //----------------------------------------------        
    
    vSession->Finalise();
    
    return 0;
}

