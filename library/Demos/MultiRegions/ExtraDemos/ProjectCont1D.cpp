#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <MultiRegions/ContField1D.h>
#include <SpatialDomains/MeshGraph.h>

using namespace std;
using namespace Nektar;

// compile using Builds/Demos/StdRegions -> make DEBUG=1 ProjectCont1D

// This routine projects a which has energy in all mdoes of the
// expansions and report an error.

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    MultiRegions::ContField1DSharedPtr Exp,Sol;

    int     i,j;
    int     order, nq;
    int     coordim;
    Array<OneD,NekDouble> sol;
    Array<OneD,NekDouble>  xc0,xc1,xc2;

    if(argc != 2)
    {
        fprintf(stderr,"Usage: ProjectCont1D mesh \n");
        exit(1);
    }

    //----------------------------------------------
    // read the problem parameters from input file
    SpatialDomains::MeshGraphSharedPtr graph1D = SpatialDomains::MeshGraph::Read(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    const SpatialDomains::ExpansionMap &expansions = graph1D->GetExpansions();
    LibUtilities::BasisKey bkey0
                            = expansions.begin()->second->m_basisKeyVector[0];
    int nmodes = bkey0.GetNumModes(); 
    cout << "Solving 1D Continuous Projection"  << endl; 
    cout << "    Expansion  : (" << LibUtilities::BasisTypeMap[bkey0.GetBasisType()] <<")" << endl;
    cout << "    No. modes  : " << nmodes << endl;
    cout << endl;
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    Exp = MemoryManager<MultiRegions::ContField1D>
                ::AllocateSharedPtr(vSession,graph1D,vSession->GetVariable(0));
    //----------------------------------------------

    //----------------------------------------------
    // Define solution to be projected
    coordim = Exp->GetCoordim(0);
    nq      = Exp->GetTotPoints();
    order   = Exp->GetExp(0)->GetNcoeffs();

    // define coordinates and solution
    sol = Array<OneD,NekDouble>(nq);

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

    for(i = 0; i < nq; ++i)
    {
        sol[i] = 0.0;
        for(j = 0; j < order; ++j)
        {
            sol[i] += pow(xc0[i],j);
            sol[i] += pow(xc1[i],j);
            sol[i] += pow(xc2[i],j);
        }
    }
    //----------------------------------------------

    //----------------------------------------------
    // Setup Temporary expansion and put in solution
    Sol = MemoryManager<MultiRegions::ContField1D>
                                ::AllocateSharedPtr(*Exp);
    Sol->SetPhys(sol);
    //----------------------------------------------

    //---------------------------------------------
    // Project onto Expansion
    Exp->FwdTrans(Sol->GetPhys(), Exp->UpdateCoeffs());
    //---------------------------------------------

    //-------------------------------------------
    // Backward Transform Solution to get projected values
    Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys());
    //-------------------------------------------

    //--------------------------------------------
    // Calculate L_inf error
    cout << "L infinity error: " << Exp->Linf(Sol->GetPhys()) << endl;
    cout << "L 2 error:        " << Exp->L2  (Sol->GetPhys()) << endl;
    //--------------------------------------------

    vSession->Finalise();

    return 0;
}

