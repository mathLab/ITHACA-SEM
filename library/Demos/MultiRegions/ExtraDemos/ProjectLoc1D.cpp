#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <SpatialDomains/MeshPartition.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList1D.h>

using namespace Nektar;

// compile using Builds/Demos/StdRegions -> make DEBUG=1 ProjectLoc1D

// This routine projects a polynomial which has energy in all mdoes of
// the expansions and report an error.

int main(int argc, char *argv[])
{
    LibUtilities::CommSharedPtr vComm
        = LibUtilities::GetCommFactory().CreateInstance("ParallelMPI",argc,argv);

    string meshfile(argv[4]);

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

    MultiRegions::ExpList1DSharedPtr Exp,Sol;
    int i,j,k;
    int     order, nq;
    int     coordim;
    LibUtilities::PointsType Qtype;
    LibUtilities::BasisType  btype;  
    Array<OneD, NekDouble> sol; 
    Array<OneD, NekDouble> xc0,xc1,xc2; 
    
    if(argc != 5)
    {
    fprintf(stderr,"Usage: ProjectLoc1D Type order nq  mesh \n");
    
    fprintf(stderr,"Where type is an integer value which "
        "dictates the basis as:\n");
    fprintf(stderr,"\t Ortho_A    = 1\n");
    fprintf(stderr,"\t Modified_A = 4\n");
    fprintf(stderr,"\t Lagrange   = 8\n");
    fprintf(stderr,"\t Legendre   = 9\n"); 
    fprintf(stderr,"\t Chebyshev  = 10\n");
    
    fprintf(stderr,"Note type = 1,2,4,5 are for higher dimensional basis\n");
    
    exit(1);
    }
    
    btype =   (LibUtilities::BasisType) atoi(argv[1]);
    
    // Check to see that only 1D Expansions are used
    if((btype == LibUtilities::eOrtho_B)||(btype == LibUtilities::eOrtho_B)||
       (btype == LibUtilities::eModified_B)||(btype == LibUtilities::eModified_C))
    NEKERROR(ErrorUtil::efatal,
             "This basis is for 2 or 3D expansions");
    
    // Do not use Fourier expansion
    if(btype == LibUtilities::eFourier)
    {
    NEKERROR(ErrorUtil::efatal,
             "Demo not set up for Fourier Expanison");
    }
    
    order  =   atoi(argv[2]);
    nq     =   atoi(argv[3]);
    
    if(btype != LibUtilities::eFourier)
    {
    Qtype = LibUtilities::eGaussLobattoLegendre; 
    }
    else
    {
    Qtype = LibUtilities::eFourierEvenlySpaced;
    }
  
    // read in mesh
    SpatialDomains::MeshGraph1D graph1D;
    graph1D.ReadGeometry(meshfile);
    graph1D.ReadExpansions(meshfile);
    
    // Define Expansion
    const LibUtilities::PointsKey Pkey(nq,Qtype);
    const LibUtilities::BasisKey Bkey(btype,order,Pkey);
    Exp = MemoryManager<MultiRegions::ExpList1D>::AllocateSharedPtr(vComm,Bkey,graph1D);
    
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
    for(j = 0; j < order; ++j)
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
    // Write solution 
    ofstream outfile("ProjectLocFile1D.dat");
    Exp->WriteToFile(outfile);
    //-------------------------------------------
    
    //--------------------------------------------
    // Calculate L_inf error
    if (vComm->GetRank() == 0)
    {
        cout << "Rank 0:" << endl;
        cout << "L infinity error: " << Exp->Linf(Sol->GetPhys()) << endl;
        cout << "L 2 error:        " << Exp->L2  (Sol->GetPhys()) << endl;
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
        data[0] = Exp->Linf(Sol->GetPhys());
        data[1] = Exp->L2  (Sol->GetPhys());
        vComm->Send(0,data);
    }
    //--------------------------------------------

    return 0;
}
