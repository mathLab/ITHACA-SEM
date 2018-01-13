#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/ExpList1D.h>
#include <SpatialDomains/MeshGraph.h>

using namespace std;
using namespace Nektar;

// compile using Builds/Demos/StdRegions -> make DEBUG=1 ProjectLoc1D

// This routine projects a polynomial which has energy in all mdoes of
// the expansions and report an error.

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    MultiRegions::ExpList1DSharedPtr Exp,Sol;
    int i,j;
    int nq;
    int coordim;
    Array<OneD, NekDouble> sol; 
    Array<OneD, NekDouble> xc0,xc1,xc2; 

    // read in mesh
    SpatialDomains::MeshGraphSharedPtr graph1D = SpatialDomains::MeshGraph::Read(vSession);

    // Define Expansion
    const SpatialDomains::ExpansionMap &expansions = graph1D->GetExpansions();
    LibUtilities::BasisKey bkey0 = expansions.begin()->second->m_basisKeyVector[0];
    int nmodes = bkey0.GetNumModes();

    Exp = MemoryManager<MultiRegions::ExpList1D>::AllocateSharedPtr(vSession,bkey0,graph1D);

    //----------------------------------------------
    // Define solution to be projected 
    coordim = Exp->GetCoordim(0);
    nq      = Exp->GetTotPoints();

    // define coordinates and solution
    sol = Array<OneD, NekDouble>(nq);

    xc0 = Array<OneD, NekDouble>(nq);
    xc1 = Array<OneD, NekDouble>(nq);
    xc2 = Array<OneD, NekDouble>(nq);

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
        for(j = 0; j < nmodes; ++j)
        {
            sol[i] += pow(xc0[i],j);
            sol[i] += pow(xc1[i],j);
            sol[i] += pow(xc2[i],j);
        }
    }

    //---------------------------------------------
    // Set up ExpList1D containing the solution 
    Sol = MemoryManager<MultiRegions::ExpList1D>::AllocateSharedPtr(*Exp);
    Sol->SetPhys(sol);
    //---------------------------------------------

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
    if (vSession->GetComm()->GetRank() == 0)
    {
        cout << "L infinity error: " << Exp->Linf(Sol->GetPhys()) << endl;
        cout << "L 2 error:        " << Exp->L2  (Sol->GetPhys()) << endl;
    }
    //--------------------------------------------

    vSession->Finalise();

    return 0;
}
