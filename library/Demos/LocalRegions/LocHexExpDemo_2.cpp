#include <cstdio>
#include <cstdlib>

#include <LocalRegions/HexExp.h>
#include <SpatialDomains/MeshGraph3D.h>

using namespace Nektar;


NekDouble Hex_sol(NekDouble x, NekDouble y, NekDouble z, int order1, int order2, int order3);

int main(int argc, char *argv[])
{
    if(argc != 5)
    {
        cerr << "usage: LocHexExpDemo MeshFile nummodes0 nummodes1 nummodes2" << endl;
        exit(1);
    }

    int i;

    string in(argv[1]);
    int nummodes0 = atoi(argv[2]);
    int nummodes1 = atoi(argv[3]);
    int nummodes2 = atoi(argv[4]);

    int nquad0 = nummodes0 + 1;
    int nquad1 = nummodes1 + 1;
    int nquad2 = nummodes2 + 1;
    int ntotquad = nquad0*nquad1*nquad2;

    SpatialDomains::MeshGraph3D graph3D; 
    graph3D.ReadGeometry(in);
    
    LibUtilities::BasisType  bType = LibUtilities::eModified_A;
    LibUtilities::PointsType qtype = LibUtilities::eGaussLobattoLegendre;
    
    const LibUtilities::PointsKey   pointsKey0( nquad0, qtype );
    const LibUtilities::PointsKey   pointsKey1( nquad1, qtype );
    const LibUtilities::PointsKey   pointsKey2( nquad2, qtype );
    
    const LibUtilities::BasisKey    basisKey0( bType, nummodes0, pointsKey0 );
    const LibUtilities::BasisKey    basisKey1( bType, nummodes1, pointsKey1 );
    const LibUtilities::BasisKey    basisKey2( bType, nummodes2, pointsKey2 );

    SpatialDomains::HexGeomSharedPtr geom;
    if(!(geom = boost::dynamic_pointer_cast<SpatialDomains::HexGeom>(graph3D.GetCompositeItem(0,0))))
    {
        cerr << "Could not find HexGeom in input file" << endl;
        exit(1);
    }

    LocalRegions::HexExpSharedPtr E = MemoryManager<LocalRegions::HexExp>::
        AllocateSharedPtr(basisKey0,basisKey1,basisKey2,geom);
    
   
    Array<OneD,NekDouble> x(ntotquad);
    Array<OneD,NekDouble> y(ntotquad);
    Array<OneD,NekDouble> z(ntotquad);        
    E->GetCoords(x,y,z); 
    
    Array<OneD,NekDouble> sol(ntotquad); 
    //----------------------------------------------
    // Define solution to be projected
    for(i = 0; i < ntotquad; i++)
    {
        sol[i]  = Hex_sol(x[i],y[i],z[i],nummodes0,nummodes1,nummodes2);
    }
    //----------------------------------------------
           
    //---------------------------------------------
    // Project onto Expansion 
    E->FwdTrans( sol, E->UpdateCoeffs() );
    //---------------------------------------------
    
    //-------------------------------------------
    // Backward Transform Solution to get projected values
    E->BwdTrans( E->GetCoeffs(), E->UpdatePhys() );
    //-------------------------------------------  
    
    //--------------------------------------------
    // Calculate L_p error 
    cout << "L infinity error: " << E->Linf(sol) << endl;
    cout << "L 2 error:        " << E->L2  (sol) << endl;
    //--------------------------------------------
    return 0;
}

NekDouble Hex_sol(NekDouble x, NekDouble y, NekDouble z, int nummodes1, int nummodes2, int nummodes3)
{
    NekDouble sol = 0;
 
    for(int i = 0; i < nummodes1; ++i) 
    {
        for(int j = 0; j < nummodes2; ++j) 
        {
            for(int k = 0; k < nummodes3; ++k) 
            {
                sol += pow(x,i) * pow(y,j) * pow(z,k);
            }
        }
    }

    return sol;
}

