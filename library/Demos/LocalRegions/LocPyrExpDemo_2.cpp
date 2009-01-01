#include <cstdio>
#include <cstdlib>

#include <LocalRegions/PyrExp.h>
#include <SpatialDomains/MeshGraph3D.h>

using namespace Nektar;


NekDouble Pyr_sol(NekDouble x, NekDouble y, NekDouble z, int order1, int order2, int order3);
NekDouble Pyramid_sol2(NekDouble x, NekDouble y, NekDouble z, int P, int Q, int R);

int main(int argc, char *argv[])
{
    if(argc != 5)
    {
        cerr << "usage: LocPyrExpDemo MeshFile nummodes0 nummodes1 nummodes2" << endl;
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
    
    LibUtilities::BasisType  bTypeX = LibUtilities::eModified_A;
    LibUtilities::BasisType  bTypeY = LibUtilities::eModified_A;
    LibUtilities::BasisType  bTypeZ = LibUtilities::eModified_C;
    
    LibUtilities::PointsType qtypeX = LibUtilities::eGaussLobattoLegendre;
    LibUtilities::PointsType qtypeY = LibUtilities::eGaussLobattoLegendre;
    LibUtilities::PointsType qtypeZ = LibUtilities::eGaussRadauMAlpha2Beta0;
    
    const LibUtilities::PointsKey   pointsKey0( nquad0, qtypeX );
    const LibUtilities::PointsKey   pointsKey1( nquad1, qtypeY );
    const LibUtilities::PointsKey   pointsKey2( nquad2, qtypeZ );
    
    const LibUtilities::BasisKey    basisKey0( bTypeX, nummodes0, pointsKey0 );
    const LibUtilities::BasisKey    basisKey1( bTypeY, nummodes1, pointsKey1 );
    const LibUtilities::BasisKey    basisKey2( bTypeZ, nummodes2, pointsKey2 );

    SpatialDomains::PyrGeomSharedPtr geom;
    if(!(geom = boost::dynamic_pointer_cast<SpatialDomains::PyrGeom>(graph3D.GetCompositeItem(0,0))))
    {
        cerr << "Could not find PyrGeom in input file" << endl;
        exit(1);
    }

    LocalRegions::PyrExpSharedPtr E = MemoryManager<LocalRegions::PyrExp>::
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
        sol[i]  = Pyr_sol(x[i],y[i],z[i],nummodes0,nummodes1,nummodes2);
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

NekDouble Pyr_sol(NekDouble x, NekDouble y, NekDouble z, int nummodes1, int nummodes2, int nummodes3)
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

NekDouble Pyramid_sol2(NekDouble x, NekDouble y, NekDouble z, int P, int Q, int R) {
    NekDouble sol = 0;

    for(int p = 0; p <= P; ++p) {
        for(int q = 0; q <= Q ; ++q) {
            for(int r = 0; r <= R - p - q; ++r) {
                sol += pow(x,p) * pow(y,q) * pow(z,r);
            }
        }
    }

    return sol;
}