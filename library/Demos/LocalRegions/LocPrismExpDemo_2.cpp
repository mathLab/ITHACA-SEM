#include <cstdio>
#include <cstdlib>

#include <LocalRegions/PrismExp.h>
#include <SpatialDomains/MeshGraph3D.h>

using namespace Nektar;


NekDouble Prism_sol(NekDouble x, NekDouble y, NekDouble z, int order1, int order2, int order3);
NekDouble Prism_sol2(NekDouble x, NekDouble y, NekDouble z, int P, int Q, int R,
                    LibUtilities::BasisType bType_x,
                    LibUtilities::BasisType bType_y,
                    LibUtilities::BasisType bType_z );

int main(int argc, char *argv[])
{
    if(argc != 5)
    {
        cerr << "usage: LocPrismExpDemo MeshFile nummodes0 nummodes1 nummodes2" << endl;
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
    LibUtilities::BasisType  bTypeZ = LibUtilities::eModified_B;
    
    LibUtilities::PointsType qtypeX = LibUtilities::eGaussLobattoLegendre;
    LibUtilities::PointsType qtypeY = LibUtilities::eGaussLobattoLegendre;
    LibUtilities::PointsType qtypeZ = LibUtilities::eGaussRadauMAlpha1Beta0;
    
    const LibUtilities::PointsKey   pointsKey0( nquad0, qtypeX );
    const LibUtilities::PointsKey   pointsKey1( nquad1, qtypeY );
    const LibUtilities::PointsKey   pointsKey2( nquad2, qtypeZ );
    
    const LibUtilities::BasisKey    basisKey0( bTypeX, nummodes0, pointsKey0 );
    const LibUtilities::BasisKey    basisKey1( bTypeY, nummodes1, pointsKey1 );
    const LibUtilities::BasisKey    basisKey2( bTypeZ, nummodes2, pointsKey2 );

    SpatialDomains::PrismGeomSharedPtr geom;
    if(!(geom = boost::dynamic_pointer_cast<SpatialDomains::PrismGeom>(graph3D.GetCompositeItem(0,0))))
    {
        cerr << "Could not find PrismGeom in input file" << endl;
        exit(1);
    }

    LocalRegions::PrismExpSharedPtr E = MemoryManager<LocalRegions::PrismExp>::
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
        sol[i]  = Prism_sol(x[i],y[i],z[i],nummodes0,nummodes1,nummodes2);
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

NekDouble Prism_sol(NekDouble x, NekDouble y, NekDouble z, int nummodes1, int nummodes2, int nummodes3)
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

NekDouble Prism_sol2(NekDouble x, NekDouble y, NekDouble z, int P, int Q, int R,
                    LibUtilities::BasisType bType_x,
                    LibUtilities::BasisType bType_y,
                    LibUtilities::BasisType bType_z )
{

    NekDouble sol = 0;
 
    // case 1 -- Common case
    if( (bType_x != LibUtilities::eFourier) && (bType_y != LibUtilities::eFourier) && (bType_z != LibUtilities::eFourier)  )
    {
        for(int p = 0; p <= P; ++p) {
            for(int q = 0; q <= Q; ++q) {
                for(int r = 0; r <= R - p; ++r) {
                    sol += pow(x,p) * pow(y,q) * pow(z,r);
                }
            }
        }
    } else  // case 2
    if((bType_x != LibUtilities::eFourier) && (bType_y != LibUtilities::eFourier) && (bType_z == LibUtilities::eFourier))
    {
        for(int i = 0; i <= P; ++i) {
            for(int j = 0; j <= Q; ++j) {
                for(int k = 0; k <= R/2; ++k) {
                    sol += pow(x,i) * pow(y,j) * sin(M_PI*k*z) + pow(x,i) * pow(y,j) * cos(M_PI*k*z);
                }
            }
        }
    }else // case 3
    if((bType_x == LibUtilities::eFourier) && (bType_y != LibUtilities::eFourier) && (bType_z != LibUtilities::eFourier))
    {
        for(int i = 0; i <= P/2; ++i) {
            for(int j = 0; j <= Q; ++j) {
                for(int k = 0; k <= R - i; ++k) {
                    sol += sin(M_PI*i*x)* pow(y,j) * pow(z,k) + cos(M_PI*i*x)* pow(y,j) * pow(z,k);                      
                }
            }
        }
    }
    else // case 4
    if((bType_x != LibUtilities::eFourier) && (bType_y == LibUtilities::eFourier) && (bType_z != LibUtilities::eFourier))
    {
        for(int i = 0; i <= P; ++i) {
            for(int j = 0; j <= Q/2; ++j) {
                for(int k = 0; k <= R - i; ++k) {
                    sol += pow(x,i)*sin(M_PI*j*y)*pow(z,k) + pow(x,i)*cos(M_PI*j*y)*pow(z,k);                      
                }
            }
        }
    }
    
    return sol;
}


