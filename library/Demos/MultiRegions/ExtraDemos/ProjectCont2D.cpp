#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <MultiRegions/ContField2D.h>
#include <SpatialDomains/MeshGraph.h>

using namespace std;
using namespace Nektar;

// This routine projects a polynomial which has energy in all mdoes of
// the expansions and report an error.

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    MultiRegions::ContField2DSharedPtr Exp,Fce;
    int     i, j, nq,  coordim;
    Array<OneD,NekDouble>  fce; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 
    
    if(argc < 2)
    {
        fprintf(stderr,"Usage: ProjectCont2D  meshfile \n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graph2D = SpatialDomains::MeshGraph::Read(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    const SpatialDomains::ExpansionMap &expansions = graph2D->GetExpansions();
    LibUtilities::BasisKey bkey0
                            = expansions.begin()->second->m_basisKeyVector[0];
    LibUtilities::BasisKey bkey1
                            = expansions.begin()->second->m_basisKeyVector[1];
    int nmodes = bkey0.GetNumModes(); 
    cout << "Solving 2D Continuous Projection"  << endl; 
    cout << "    Expansion  : (" << LibUtilities::BasisTypeMap[bkey0.GetBasisType()] <<","<< LibUtilities::BasisTypeMap[bkey1.GetBasisType()]  << ")" << endl;
    cout << "    No. modes  : " << nmodes << endl;
    cout << endl;
    //----------------------------------------------
   
    //----------------------------------------------
    // Define Expansion 
    Exp = MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(vSession,graph2D,vSession->GetVariable(0));
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
    // Define forcing function
    fce = Array<OneD,NekDouble>(nq);  
 
    for(i = 0; i < nq; ++i)
    {
        fce[i] = 0.0;
        for(j = 0; j < nmodes; ++j)
        {
            fce[i] += pow(xc0[i],j);
            fce[i] += pow(xc1[i],j);
            fce[i] += pow(xc2[i],j);
        }
    }
    
    //---------------------------------------------
    // Set up ContField2DD containing the solution 
    Fce = MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //---------------------------------------------

    //---------------------------------------------
    // Project onto Expansion 
    Exp->FwdTrans(Fce->GetPhys(), Exp->UpdateCoeffs());
    //---------------------------------------------
    
    //-------------------------------------------
    // Backward Transform Solution to get projected values
    Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys());
    //-------------------------------------------  

    //-------------------------------------------
    // Calculate error
    Array<OneD,NekDouble> error(nq,0.0);
    Vmath::Vsub(nq, Exp->GetPhys(), 1, Fce->GetPhys(), 1, error, 1);
    //-------------------------------------------  

    //--------------------------------------------
    // Calculate L_inf and L_2 of error
    NekDouble vL2Error   = Exp->L2  (error);
    NekDouble vLinfError = Exp->Linf(error);
    if (vSession->GetComm()->GetRank() == 0)
    {
        cout << "L infinity error: " << vLinfError << endl;
        cout << "L 2 error       : " << vL2Error << endl;
    }  
    //--------------------------------------------

    vSession->Finalise();

    return 0;
}
