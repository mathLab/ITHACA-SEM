#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <SpatialDomains/MeshPartition.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList3D.h>

using namespace Nektar;

// This routine projects a polynomial which has energy in all mdoes of
// the expansions and report an error.

int main(int argc, char *argv[])
{
    LibUtilities::CommSharedPtr vComm
        = LibUtilities::GetCommFactory().CreateInstance("ParallelMPI",argc,argv);

    string meshfile(argv[1]);

    if (vComm->GetSize() > 1)
    {
        if (vComm->GetRank() == 0)
        {
            LibUtilities::SessionReaderSharedPtr vSession = MemoryManager<LibUtilities::SessionReader>::AllocateSharedPtr(meshfile);
            SpatialDomains::MeshPartitionSharedPtr vPartitioner = MemoryManager<SpatialDomains::MeshPartition>::AllocateSharedPtr(vSession);
            vPartitioner->PartitionMesh(vComm->GetSize());
            vPartitioner->WritePartitions(vSession, meshfile);
        }

        vComm->Block();

        meshfile = meshfile + "." + boost::lexical_cast<std::string>(vComm->GetRank());
    }

    MultiRegions::ExpList3DSharedPtr Exp,Fce;
    int     i, j, nq,  coordim;
    Array<OneD,NekDouble>  fce, tmp, tmp2;
    Array<OneD,NekDouble>  xc0,xc1,xc2;
    NekDouble  lambda;

    if(argc != 2)
    {
        fprintf(stderr,"Usage: ProjectLoc3D  meshfile \n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
//    string meshfile(argv[1]);
    SpatialDomains::MeshGraph3D graph3D;
    graph3D.ReadGeometry(meshfile);
    graph3D.ReadExpansions(meshfile);
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    const SpatialDomains::ExpansionMap &expansions = graph3D.GetExpansions();
    LibUtilities::BasisKey bkey = expansions.begin()->second->m_basisKeyVector[0];
    int nmodes = bkey.GetNumModes();
    if (vComm->GetRank() == 0)
    {
        cout << "Solving 3D Local Projection"  << endl;
        cout << "    No. modes  : " << nmodes << endl;
        cout << endl;
    }
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    Exp = MemoryManager<MultiRegions::ExpList3D>::AllocateSharedPtr(vComm,graph3D);
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
    tmp = Array<OneD,NekDouble>(nq);
    tmp2 = Array<OneD,NekDouble>(nq);
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
    // Set up ExpList1D containing the solution
    Fce = MemoryManager<MultiRegions::ExpList3D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //---------------------------------------------

    //---------------------------------------------
    // Project onto Expansion
    Exp->FwdTrans(Fce->GetPhys(), Exp->UpdateCoeffs());
    //---------------------------------------------

    //---------------------------------------------
    // Transform back to Physical Space
    Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys());
    //----------------------------------------------

    //--------------------------------------------
    // Calculate L_inf error
    if (vComm->GetRank() == 0)
    {
        cout << "Rank 0:" << endl;
        cout << "L infinity error: " << Exp->Linf(Fce->GetPhys()) << endl;
        cout << "L 2 error:        " << Exp->L2  (Fce->GetPhys()) << endl;
        for (unsigned int i = 1; i < vComm->GetSize(); ++i)
        {
            Array<OneD, NekDouble> data(2);
            vComm->Recv(i, data);
            cout << "Rank " << i << ":" << endl;
            cout << "L infinity error: " << data[0] << endl;
            cout << "L 2 error:        " << data[1] << endl;
        }
    }
    else
    {
        Array<OneD, NekDouble> data(2);
        data[0] = Exp->Linf(Fce->GetPhys());
        data[1] = Exp->L2  (Fce->GetPhys());
        vComm->Send(0,data);
    }
    //--------------------------------------------

    return 0;
}
