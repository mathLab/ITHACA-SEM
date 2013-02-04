#include <StdRegions/StdExpansion3D.h>
#include <LocalRegions/PrismExp.h>

#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Basis.h>
#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/QuadGeom.h>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>


#include <algorithm>
#include <iostream> 
#include <cstdlib>   
#include <cmath>
#include <iomanip>     
#include <iosfwd>
 

using namespace std;     


using namespace Nektar;  


NekDouble Prism_sol(NekDouble x, NekDouble y, NekDouble z, int order1, int order2, int order3,
                  LibUtilities::BasisType bType_x, LibUtilities::BasisType bType_y, LibUtilities::BasisType bType_z);

// using namespace boost;
using namespace Nektar::LibUtilities;  
using namespace Nektar::LocalRegions;
using namespace Nektar::StdRegions;  
using namespace Nektar::SpatialDomains;


int main(int argc, char *argv[])
 {
    if( argc != 10 ) {
        cerr << "Usage: PrismDemo Type_x Type_y Type_z numModes_x numModes_y numModes_z Qx Qy Qz" << endl;
        cerr << "Where type is an integer value which dictates the basis as:" << endl;
        cerr << "\t Ortho_A  = 1\n";
        cerr << "\t Ortho_B  = 2\n";
        cerr << "\t Ortho_C  = 3\n";
        cerr << "\t Modified_A = 4\n";
        cerr << "\t Modified_B = 5\n";
        cerr << "\t Modified_C = 6\n";
        cerr << "\t Fourier    = 7\n";
        cerr << "\t Lagrange   = 8\n";
        cerr << "\t Legendre   = 9\n";
        cerr << "\t Chebyshev  = 10\n";
        cerr << "\t Nodal Tet (Electro) = 13    (3D Nodal Electrostatic Points on a Tetrahedron)\n";
        cerr << "\n\n" << "Example: " << argv[0] << " 4 4 5 3 3 3 5 5 5" << endl;
        cerr << endl;
        
        exit(1);
    }

    StdRegions::ExpansionType regionShape = StdRegions::ePrism;
    int bType_x_val = atoi(argv[1]);
    int bType_y_val = atoi(argv[2]);
    int bType_z_val = atoi(argv[3]);
    
    LibUtilities::BasisType   bType_x = static_cast<LibUtilities::BasisType>( bType_x_val );
    LibUtilities::BasisType   bType_y = static_cast<LibUtilities::BasisType>( bType_y_val );
    LibUtilities::BasisType   bType_z = static_cast<LibUtilities::BasisType>( bType_z_val );
    
    if( (bType_x_val == 13) || (bType_y_val == 13) || (bType_z_val == 13) )
    {
        bType_x =   LibUtilities::eOrtho_A;
        bType_y =   LibUtilities::eOrtho_B;
        bType_z =   LibUtilities::eOrtho_C;
    }

    // Check to see that correct Expansions are used
    if( regionShape == StdRegions::ePrism) 
    {
        if( (bType_x == LibUtilities::eOrtho_B) || (bType_x == LibUtilities::eModified_B) ) {
            NEKERROR(ErrorUtil::efatal, "Basis 1 cannot be of type Ortho_B or Modified_B");
        }
        if( (bType_x == LibUtilities::eOrtho_C) || (bType_x == LibUtilities::eModified_C) ) {
            NEKERROR(ErrorUtil::efatal, "Basis 1 cannot be of type Ortho_C or Modified_C");
        }
        if( (bType_y == LibUtilities::eOrtho_C) || (bType_y == LibUtilities::eModified_C) ) {
            NEKERROR(ErrorUtil::efatal, "Basis 2 cannot be of type Ortho_C or Modified_C");
        }
    }
      
    int xModes   = atoi(argv[4]);
    int yModes   = atoi(argv[5]);
    int zModes   = atoi(argv[6]);
    int Qx = atoi(argv[7]);
    int Qy = atoi(argv[8]);
    int Qz = atoi(argv[9]);
    int P = xModes - 1, Q = yModes - 1, R = zModes - 1;
    const int three = 3;
    
    Array<OneD, NekDouble> solution( Qx * Qy * Qz );

    LibUtilities::PointsType    Qtype_x = eGaussLobattoLegendre;
    LibUtilities::PointsType    Qtype_y = eGaussLobattoLegendre;
    LibUtilities::PointsType    Qtype_z = eGaussRadauMAlpha1Beta0;
        
        
    //-----------------------------------------------
    // Define a 3D expansion based on basis definition    
    StdRegions::StdExpansion3D *lpr = 0;
    
    if( regionShape == StdRegions::ePrism ) 
    {
        // //////////////////////////////////////////////////////
        // Set up Prism vertex coordinates

        const int nVerts = 6;
        const double point[][3] = {
            {0,0,0}, {1,0,0}, {1,1,0}, 
            {0,1,0}, {0,0,1}, {0,1,1}
        };

       // ////////////////////////////////////////////////////////////////////////////////////
       // Populate the list of verts
       // VertexComponent (const int coordim, const int vid, double x, double y, double z)
       VertexComponentSharedPtr verts[6];
       for(int i=0; i < nVerts; ++i){
         verts[i] =  MemoryManager<VertexComponent>::
         AllocateSharedPtr( three, i, point[i][0], point[i][1], point[i][2] );
       }


       // /////////////////////////////////////////////////////////////////////
       // Set up Prism Edges
       // SegGeom (int id, const int coordim), EdgeComponent(id, coordim)

        const int nEdges = 9;
        const int vertexConnectivity[][2] = {
            {0,1}, {1,2}, {3,2}, {0,3}, {0,4}, 
            {1,4}, {2,5}, {3,5}, {4,5}
        };


        // Populate the list of edges
        SegGeomSharedPtr edges[nEdges]; 
        for(int i=0; i < nEdges; ++i){
            VertexComponentSharedPtr vertsArray[2];
            for(int j=0; j<2; ++j){
                vertsArray[j] = verts[vertexConnectivity[i][j]];
            }
            edges[i] = MemoryManager<SegGeom>::AllocateSharedPtr(i, three, vertsArray);
        }
 
        // //////////////////////////////////////////////////////////////
        // Set up Prism faces
        const int nFaces = 5;
        //quad-edge connectivity base-face0, vertical-quadface2, vertical-quadface4
        const int quadEdgeConnectivity[][4] = { {0,1,2,3}, {1,6,8,5}, {3,7,8,4} }; 
        const bool   isQuadEdgeFlipped[][4] = { {0,0,1,1}, {0,0,1,1}, {0,0,1,1} };
        // QuadId ordered as 0, 1, 2, otherwise return false
        const int                  quadId[] = { 0,-1,1,-1,2 }; 

         //triangle-edge connectivity side-triface-1, side triface-3 
        const int  triEdgeConnectivity[][3] = { {0,5,4}, {2,6,7} };
        const bool    isTriEdgeFlipped[][3] = { {0,0,1}, {0,0,1} };
        // TriId ordered as 0, 1, otherwise return false
        const int                   triId[] = { -1,0,-1,1,-1 }; 

        // Populate the list of faces  
        Geometry2DSharedPtr faces[nFaces]; 
        for(int f = 0; f < nFaces; ++f){
            if(f == 1 || f == 3) {
                int i = triId[f];
                SegGeomSharedPtr edgeArray[3];
                Orientation eorientArray[3];
                for(int j = 0; j < 3; ++j){
                    edgeArray[j] = edges[triEdgeConnectivity[i][j]];
                    eorientArray[j] = isTriEdgeFlipped[i][j] ? eBackwards : eForwards;
                }
                faces[f] = MemoryManager<TriGeom>::AllocateSharedPtr(f, edgeArray, eorientArray);
            }            
            else {
                int i = quadId[f];
                SegGeomSharedPtr edgeArray[4];
                Orientation eorientArray[4]; 
                for(int j=0; j < 4; ++j){
                    edgeArray[j] = edges[quadEdgeConnectivity[i][j]];
                    eorientArray[j] = isQuadEdgeFlipped[i][j] ? eBackwards : eForwards;
                }
                faces[f] = MemoryManager<QuadGeom>::AllocateSharedPtr(f, edgeArray, eorientArray);
            }
        } 
        
         const LibUtilities::PointsKey   pointsKey_x( Qx, Qtype_x );    
         const LibUtilities::PointsKey   pointsKey_y( Qy, Qtype_y ); 
         const LibUtilities::PointsKey   pointsKey_z( Qz, Qtype_z );

         const LibUtilities::BasisKey    basisKey_x( bType_x, xModes, pointsKey_x );
         const LibUtilities::BasisKey    basisKey_y( bType_y, yModes, pointsKey_y ); 
         const LibUtilities::BasisKey    basisKey_z( bType_z, zModes, pointsKey_z );

         Array<OneD, StdRegions::StdExpansion3DSharedPtr> xMap(3);
         for(int i=0; i < 3; ++i){
            xMap[i] = MemoryManager<StdRegions::StdPrismExp>::AllocateSharedPtr(basisKey_x, basisKey_y, basisKey_z);

         }

         SpatialDomains::PrismGeomSharedPtr geom = MemoryManager<SpatialDomains::PrismGeom>::AllocateSharedPtr(faces);
         geom->SetOwnData();
 
   
        if( bType_x_val < 10 ) {  
            lpr = new LocalRegions::PrismExp( basisKey_x, basisKey_y, basisKey_z, geom );
        } else {
            cerr << "Implement the NodalTetExp!!!!!!" << endl;
            //lpr = new StdRegions::StdNodalTetExp( basisKey_x, basisKey_y, basisKey_z, NodalType );
            exit(1);
        }
    
        Array<OneD,NekDouble> x = Array<OneD,NekDouble>( Qx * Qy * Qz );
        Array<OneD,NekDouble> y = Array<OneD,NekDouble>( Qx * Qy * Qz );
        Array<OneD,NekDouble> z = Array<OneD,NekDouble>( Qx * Qy * Qz );
        
        lpr->GetCoords(x,y,z); 
    
        //----------------------------------------------
        // Define solution to be projected
        for(int n = 0; n < Qx * Qy * Qz; ++n) { 
            solution[n]  = Prism_sol( x[n], y[n], z[n], P, Q, R, bType_x, bType_y, bType_z );
        } 
        //----------------------------------------------
    }
           
    //---------------------------------------------
    // Project onto Expansion 
    lpr->FwdTrans( solution, lpr->UpdateCoeffs() );
    //---------------------------------------------
    
    //-------------------------------------------
    // Backward Transform Solution to get projected values
    lpr->BwdTrans( lpr->GetCoeffs(), lpr->UpdatePhys() );
    //-------------------------------------------  
    
    //--------------------------------------------
    // Calculate L_p error 
    cout << "L infinity error: " << lpr->Linf(solution) << endl;
    cout << "L 2 error:        " << lpr->L2  (solution) << endl;
    //--------------------------------------------
    
    //-------------------------------------------
    // Evaulate solution at mid point  and print error
    Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);
    t[0] = 0.5;
    t[1] = 0.5;
    t[2] = 0.2;
    
    if( regionShape == StdRegions::ePrism ) {
        solution[0] = Prism_sol( t[0], t[1], t[2], P, Q, R, bType_x, bType_y, bType_z );
    }
    
    NekDouble numericSolution = lpr->PhysEvaluate(t);
    cout << "Interpolation difference from actual solution at x = ( " << 
        t[0] << ", " << t[1] << ", " << t[2] << " ): " << numericSolution - solution[0] << endl;
    //-------------------------------------------
    
    return 0;
}


NekDouble Prism_sol(NekDouble x, NekDouble y, NekDouble z, int P, int Q, int R,
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

