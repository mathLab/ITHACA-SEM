#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/DisContField1D.h>
#include <SpatialDomains/MeshGraph.h>

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    LibUtilities::CommSharedPtr vComm = vSession->GetComm();

    MultiRegions::DisContField1DSharedPtr Exp,Fce;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce;
    Array<OneD,NekDouble>  xc0,xc1,xc2;
    StdRegions::ConstFactorMap factors;

    if(argc != 2)
    {
        fprintf(stderr,"Usage: Helmholtz1D  meshfile\n");
        exit(1);
    }

    LibUtilities::FieldIOSharedPtr fld = LibUtilities::FieldIO::CreateDefault(vSession);

    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graph1D = SpatialDomains::MeshGraph::Read(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    factors[StdRegions::eFactorLambda] = vSession->GetParameter("Lambda");
    factors[StdRegions::eFactorTau] = 1.0;
    cout << "Solving 1D Helmholtz:"  << endl;
    cout << "         Lambda     : " << factors[StdRegions::eFactorLambda] << endl;
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    Exp = MemoryManager<MultiRegions::DisContField1D>::
        AllocateSharedPtr(vSession,graph1D,vSession->GetVariable(0));
    //----------------------------------------------

    //----------------------------------------------
    // Set up coordinates of mesh for Forcing function evaluation
    coordim = Exp->GetCoordim(0);
    nq      = Exp->GetTotPoints();

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
    LibUtilities::EquationSharedPtr ffunc = vSession->GetFunction("Forcing", 0);

    ffunc->Evaluate(xc0, xc1, xc2, fce);

    //----------------------------------------------

    //----------------------------------------------
    // Setup expansion containing the  forcing function
    Fce = MemoryManager<MultiRegions::DisContField1D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //----------------------------------------------

    //----------------------------------------------
    // Helmholtz solution taking physical forcing
    Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), factors);
    //----------------------------------------------

    //----------------------------------------------
    // Backward Transform Solution to get solved values at
    Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys());
    //----------------------------------------------

    //----------------------------------------------
    // Write solution
    string   out(strtok(argv[1],"."));
    string   endfile(".fld");
    out += endfile;
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
        = Exp->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
    for(i = 0; i < FieldDef.size(); ++i)
    {
        FieldDef[i]->m_fields.push_back("u");
        Exp->AppendFieldData(FieldDef[i], FieldData[i]);
    }
    fld->Write(out, FieldDef, FieldData);
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

        ex_sol->Evaluate(xc0, xc1, xc2, fce);

        //----------------------------------------------

        //--------------------------------------------
        // Calculate L_inf error
        Fce->SetPhys(fce);
        cout << "L infinity error: " << Exp->Linf(Exp->GetPhys(), Fce->GetPhys()) << endl;
        cout << "L 2 error:        " << Exp->L2  (Exp->GetPhys(), Fce->GetPhys()) << endl;
        //--------------------------------------------
    }
    //----------------------------------------------

    vComm->Finalise();

    return 0;
}

