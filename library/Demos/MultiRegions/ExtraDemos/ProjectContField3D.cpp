#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ContField3D.h>
#include <SpatialDomains/MeshGraph3D.h>

using namespace Nektar;

// This routine projects a polynomial which has energy in all mdoes of
// the expansions and report an error.

int main(int argc, char *argv[])
{ 
    LibUtilities::SessionReaderSharedPtr vSession
                    = LibUtilities::SessionReader::CreateInstance(argc, argv);

    MultiRegions::ContField3DSharedPtr Exp,Fce;
    int nq, coordim;
    Array<OneD,NekDouble>  fce; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 
    std::string meshfile(vSession->GetFilename());
    
    if(argc != 3)
    {
        fprintf(stderr,"Usage: ProjectContField3D  meshfile \n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graph3D = MemoryManager<SpatialDomains::MeshGraph3D>::AllocateSharedPtr(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    const SpatialDomains::ExpansionMap &expansions = graph3D->GetExpansions();
    LibUtilities::BasisKey bkey
                            = expansions.begin()->second->m_basisKeyVector[0];
    cout << "Solving 3D C0 continuous Projection (with boundary conditions)"  << endl; 
    cout << "    No. modes  : " << bkey.GetNumModes() << endl;
    cout << endl;
    //----------------------------------------------
   
    //----------------------------------------------
    // Define Expansion 
    Exp = MemoryManager<MultiRegions::ContField3D>::AllocateSharedPtr(vSession,graph3D,vSession->GetVariable(0));
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
    LibUtilities::EquationSharedPtr ffunc = vSession->GetFunction("Forcing",0);

    ffunc->Evaluate(xc0,xc1,xc2,fce);
    //----------------------------------------------
    
    //---------------------------------------------
    // Set up ExpList containing the solution 
    Fce = MemoryManager<MultiRegions::ContField3D>::AllocateSharedPtr(*Exp);
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

    //----------------------------------------------
    // Write solution 
    ofstream outfile2("ProjectContFieldFile3D.dat");
    Exp->WriteToFile(outfile2,eTecplot);
    outfile2.close();
    //----------------------------------------------
    
    //--------------------------------------------
    // Calculate L_inf error 
    cout << "L infinity error: " << Exp->Linf(Fce->GetPhys()) << endl;
    cout << "L 2 error:        " << Exp->L2  (Fce->GetPhys()) << endl;
    //--------------------------------------------

    vSession->Finalise();

    return 0;
}
