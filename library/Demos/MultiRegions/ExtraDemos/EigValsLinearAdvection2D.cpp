#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SpatialDomains/MeshPartition.h>
#include <MultiRegions/ContField2D.h>

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    string meshfile(vSession->GetFilename());

    MultiRegions::ContField2DSharedPtr Exp;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 
    NekDouble  ax,ay;

    if(argc != 3)
    {
        fprintf(stderr,"Usage: EigValsLinearAdvection2D  meshfile boundaryfile\n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraph2D graph2D; 

    graph2D.ReadGeometry(meshfile);
    graph2D.ReadExpansions(meshfile);
    //----------------------------------------------

    //----------------------------------------------
    // read the problem parameters from input file
    string bcfile(argv[argc-1]);
    SpatialDomains::BoundaryConditions bcs(&graph2D); 
    bcs.Read(bcfile);
    //----------------------------------------------

    //----------------------------------------------
    // Get Advection Velocity
    ax = bcs.GetParameter("Advection_x");
    ay = bcs.GetParameter("Advection_y");
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    const SpatialDomains::ExpansionMap &expansions = graph2D.GetExpansions();
    LibUtilities::BasisKey bkey0
                            = expansions.begin()->second->m_basisKeyVector[0];
    LibUtilities::BasisKey bkey1
                            = expansions.begin()->second->m_basisKeyVector[1];
    cout << "Calc. LinearAdvection E-vals:"  << endl; 
    cout << "             Advection_x    : " << ax << endl; 
    cout << "             Advection_y    : " << ay << endl; 
    cout << "             Expansion      : (" << SpatialDomains::kExpansionTypeStr[bkey0.GetBasisType()] <<","<< SpatialDomains::kExpansionTypeStr[bkey1.GetBasisType()]  << ")" << endl;
    cout << "             No. modes      : " << bkey0.GetNumModes() << endl;
    cout << endl;
    //----------------------------------------------
   
    //----------------------------------------------
    // Define Expansion 
    Exp = MemoryManager<MultiRegions::ContField2D>::
        AllocateSharedPtr(vSession,graph2D,bcs);
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
    // Evaluate Eigenvalues/vector of Linear Advection Op.
    int nmodes = Exp->GetNcoeffs();
    Array<OneD, NekDouble> Real(nmodes), Imag(nmodes), Evecs(nmodes*nmodes);
    
    Exp->LinearAdvectionEigs(ax,ay,Real,Imag,Evecs);
    //----------------------------------------------
    
    cout << "Eignvalues: [Real, Imag]" << endl;
    for(i = 0; i < Real.num_elements(); ++i)
    {
        cout << Real[i] <<", "<< Imag[i] << endl;
    }
    
    return 0;
}

