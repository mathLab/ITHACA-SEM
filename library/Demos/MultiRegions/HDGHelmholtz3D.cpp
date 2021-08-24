#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/DisContField.h>
#include <SpatialDomains/MeshGraph.h>

//#define TIMING
#ifdef TIMING
#include <time.h>
#define Timing(s) \
 fprintf(stdout,"%s Took %g seconds\n",s,(clock()-st)/(double)CLOCKS_PER_SEC); \
 st = clock();
#else
#define Timing(s) \
 /* Nothing */
#endif

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    LibUtilities::CommSharedPtr vComm = vSession->GetComm();

    MultiRegions::DisContFieldSharedPtr Exp, Fce;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 
    StdRegions::ConstFactorMap factors;

    if(argc < 2)
    {
        fprintf(stderr,"Usage: HDGHelmholtz3D  meshfile [solntype]\n");
        exit(1);
    }

    LibUtilities::FieldIOSharedPtr fld = LibUtilities::FieldIO::CreateDefault(vSession);

    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graph3D = 
        SpatialDomains::MeshGraph::Read(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    factors[StdRegions::eFactorLambda] = vSession->GetParameter("Lambda");
    factors[StdRegions::eFactorTau] = 1.0;
    const SpatialDomains::ExpansionInfoMap &expansions = graph3D->GetExpansionInfo();
    LibUtilities::BasisKey bkey0
                            = expansions.begin()->second->m_basisKeyVector[0];

    if (vComm->GetRank() == 0)
    {
            cout << "Solving 3D Helmholtz:"  << endl;
            cout << "  - Communication: " 
                 << vSession->GetComm()->GetType() << " (" 
                 << vSession->GetComm()->GetSize() 
                 << " processes)" << endl;
            cout << "  - Solver type  : " 
                 << vSession->GetSolverInfo("GlobalSysSoln") << endl;
            cout << "  - Lambda       : " 
                 << factors[StdRegions::eFactorLambda] << endl;
            cout << "  - No. modes    : " 
                 << bkey0.GetNumModes() << endl;
            cout << endl;
    }
    //----------------------------------------------
   
    //----------------------------------------------
    // Define Expansion 
    Exp = MemoryManager<MultiRegions::DisContField>::
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
    Fce = MemoryManager<MultiRegions::DisContField>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //----------------------------------------------
    Timing("Define forcing ..");
  
    //----------------------------------------------
    // Helmholtz solution taking physical forcing 
    Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(),  factors);
    //----------------------------------------------
    
    Timing("Helmholtz Solve ..");

    //----------------------------------------------
    // Backward Transform Solution to get solved values at 
    Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys());
    //----------------------------------------------
    Timing("Backward Transform ..");
    
    //-----------------------------------------------
    // Write solution to file
    string out = vSession->GetSessionName() + ".fld";
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                                                = Exp->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    for(i = 0; i < FieldDef.size(); ++i)
    {
        FieldDef[i]->m_fields.push_back("u");
        Exp->AppendFieldData(FieldDef[i], FieldData[i]);
    }
    fld->Write(out, FieldDef, FieldData);
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

        NekDouble vLinfError = Exp->Linf(Exp->GetPhys(), Fce->GetPhys());
        NekDouble vL2Error   = Exp->L2  (Exp->GetPhys(), Fce->GetPhys());
        NekDouble vH1Error   = Exp->H1  (Exp->GetPhys(), Fce->GetPhys());

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

