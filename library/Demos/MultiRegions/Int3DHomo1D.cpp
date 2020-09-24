#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <SpatialDomains/MeshGraph.h>

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession = LibUtilities::SessionReader::CreateInstance(argc, argv);

    LibUtilities::CommSharedPtr vComm = vSession->GetComm();
    
    MultiRegions::ContField3DHomogeneous1DSharedPtr Exp_u, Exp_v, Exp_w;
    
    StdRegions::ConstFactorMap factors;
    FlagList flags;

    if( (argc != 2) && (argc != 3))
    {
        fprintf(stderr,"Usage: Int3DHomo2D meshfile [SysSolnType]   \n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graph2D = SpatialDomains::MeshGraph::Read(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    int nzpoints;
    NekDouble lz;
    int FFT;

    vSession->LoadParameter("HomModesZ", nzpoints);
    vSession->LoadParameter("LZ",        lz);
    vSession->LoadParameter("USEFFT",    FFT);
    
    bool useFFT = false;
    bool deal = false;
    if(FFT==1){useFFT = true;}
    
    const LibUtilities::PointsKey PkeyZ(nzpoints,LibUtilities::eFourierEvenlySpaced);
    const LibUtilities::BasisKey  BkeyZ(LibUtilities::eFourier,nzpoints,PkeyZ);
    
    Exp_u = MemoryManager<MultiRegions::ContField3DHomogeneous1D>::AllocateSharedPtr(vSession,BkeyZ,lz,useFFT,deal,graph2D,vSession->GetVariable(0));
    Exp_v = MemoryManager<MultiRegions::ContField3DHomogeneous1D>::AllocateSharedPtr(vSession,BkeyZ,lz,useFFT,deal,graph2D,vSession->GetVariable(1));
    Exp_w = MemoryManager<MultiRegions::ContField3DHomogeneous1D>::AllocateSharedPtr(vSession,BkeyZ,lz,useFFT,deal,graph2D,vSession->GetVariable(2));
    
    //----------------------------------------------
    
    //----------------------------------------------
    // Print summary of solution details
    flags.set(eUseGlobal, false);
    
    const SpatialDomains::ExpansionInfoMap &expansions = graph2D->GetExpansionInfo();
    
    LibUtilities::BasisKey bkey0 = expansions.begin()->second->m_basisKeyVector[0];
    
    cout << "Calculating Integral (Homogeneous in z-plane):"  << endl;
    cout << "         Lz              : " << lz << endl;
    cout << "         N.modes         : " << bkey0.GetNumModes() << endl;
    cout << "         N.Z homo modes  : " << BkeyZ.GetNumModes() << endl;
    cout << endl;
    //----------------------------------------------
    
    //----------------------------------------------
    // Set up coordinates of mesh for Forcing function evaluation
    int nq  = Exp_u->GetTotPoints();
    
    Array<OneD,NekDouble>  xc0,xc1,xc2;
    
    xc0 = Array<OneD,NekDouble>(nq,0.0);
    xc1 = Array<OneD,NekDouble>(nq,0.0);
    xc2 = Array<OneD,NekDouble>(nq,0.0);
    
    Exp_u->GetCoords(xc0,xc1,xc2);
    //----------------------------------------------
    // Define fields
    LibUtilities::EquationSharedPtr ffunc_u = vSession->GetFunction("InitialCondition", 0);
    LibUtilities::EquationSharedPtr ffunc_v = vSession->GetFunction("InitialCondition", 1);
    LibUtilities::EquationSharedPtr ffunc_w = vSession->GetFunction("InitialCondition", 2);    
    
    ffunc_u->Evaluate(xc0,xc1,xc2,Exp_u->UpdatePhys());
    ffunc_v->Evaluate(xc0,xc1,xc2,Exp_v->UpdatePhys());
    ffunc_w->Evaluate(xc0,xc1,xc2,Exp_w->UpdatePhys());
    
    //----------------------------------------------
    //Calculate integral and print it as L inf error
    
    NekDouble intU = Exp_u->Integral(Exp_u->GetPhys());
    
    cout << "L infinity error (variable intU): " << intU << endl;
        
    NekDouble intV = Exp_v->Integral(Exp_v->GetPhys());
    
    cout << "L infinity error (variable intV): " << intV << endl;
        
    NekDouble intW = Exp_w->Integral(Exp_w->GetPhys());
    
    cout << "L infinity error (variable intW): " << intW << endl;
        
    return 0;
}
