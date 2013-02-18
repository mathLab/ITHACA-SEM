#include <StdRegions/StdExpansion3D.h>
#include <LocalRegions/PyrExp.h>

#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/QuadGeom.h>


#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <iosfwd>


using namespace std;


using namespace Nektar;


NekDouble Pyr_sol(NekDouble x, NekDouble y, NekDouble z, int order1, int order2, int order3);

// using namespace boost;
using namespace Nektar::LibUtilities;
using namespace Nektar::LocalRegions;
using namespace Nektar::StdRegions;
using namespace Nektar::SpatialDomains;


int main(int argc, char *argv[])
 {
    if( argc != 10 ) {
        cerr << "Usage: PyramidDemo Type_x Type_y Type_z numModes_x numModes_y numModes_z Qx Qy Qz" << endl;
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
        cerr << "\n\n" << "Example: " << argv[0] << " 4 4 4 3 3 3 5 5 5" << endl;
        cerr << endl;
        
        exit(1);
    }

    StdRegions::ExpansionType regionShape = StdRegions::ePyramid;
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
    if( regionShape == StdRegions::ePyramid ) 
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
    LibUtilities::PointsType    Qtype_z = eGaussRadauMAlpha2Beta0;
        
    //-----------------------------------------------
    // Define a 3D expansion based on basis definition
    
    StdRegions::StdExpansion3D *lpe = 0;
    
    if( regionShape == StdRegions::ePyramid ) 
    {
        // //////////////////////////////////////////////////////
        // Set up Prism vertex coordinates

        const int nVerts = 5;
        const double point[][3] = {
            {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {0,0,1}
        };

       // ////////////////////////////////////////////////////////////////////////////////////
       // Populate the list of verts
       // VertexComponent (const int coordim, const int vid, double x, double y, double z)
       VertexComponentSharedPtr verts[5];
       for(int i=0; i < nVerts; ++i){
         verts[i] =  MemoryManager<VertexComponent>::
         AllocateSharedPtr( three, i, point[i][0], point[i][1], point[i][2] );
       }


       // /////////////////////////////////////////////////////////////////////
       // Set up Pyramid Edges
       // SegGeom (int id, const int coordim), EdgeComponent(id, coordim)

        const int nEdges = 8;
        const int vertexConnectivity[][2] = {
            {0,1},{1,2},{2,3},{0,3},{0,4},
            {1,4},{2,4},{3,4}
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
        const int quadEdgeConnectivity[][4] = { {0,1,2,3} }; 
        const bool   isQuadEdgeFlipped[][4] = { {0,0,1,1} };
        // QuadId ordered as 0, 1, 2, otherwise return false
        const int                  quadId[] = { 0,-1,-1,-1,-1 };

         //triangle-edge connectivity side-triface-1, side triface-3 
        const int  triEdgeConnectivity[][3] = { {0,5,4}, {1,6,5}, {2,7,6}, {3,7,4} };
        const bool    isTriEdgeFlipped[][3] = { {0,0,1}, {0,0,1}, {1,0,1}, {0,0,1} };  
        // TriId ordered as 0, 1, otherwise return false
        const int                   triId[] = { -1,0,1,2,3 };

        // Populate the list of faces  
        Geometry2DSharedPtr faces[nFaces]; 
        for(int f = 0; f < nFaces; ++f){
            cout << "f = " << f << endl;
            if(f==0){
                int i = quadId[f];
                SegGeomSharedPtr edgeArray[4];
                Orientation eorientArray[4]; 
                for(int j=0; j < 4; ++j){
                    edgeArray[j] = edges[quadEdgeConnectivity[i][j]];
                    eorientArray[j] = isQuadEdgeFlipped[i][j] ? eBackwards : eForwards;
                }
                faces[f] = MemoryManager<QuadGeom>::AllocateSharedPtr(f, edgeArray, eorientArray);
                cout << "quad faces[" << f << "] = " << faces[f] << endl;
            }
            else {
                int i = triId[f];
                SegGeomSharedPtr edgeArray[3];
                Orientation eorientArray[3];
                for(int j=0; j < 3; ++j){
                    edgeArray[j] = edges[triEdgeConnectivity[i][j]];
                    eorientArray[j] = isTriEdgeFlipped[i][j] ? eBackwards : eForwards;
                }
                faces[f] = MemoryManager<TriGeom>::AllocateSharedPtr(f, edgeArray, eorientArray);
                cout << "tri faces[" << f << "] = " << faces[f] << endl;

            }

            cout << "faces[" << f << "] = " << faces[f] << endl;
            for(int k = 0; k < faces[f]->GetNumEdges(); ++k) {
                cout << "faces[" << f << "]->GetEid(" << k << ") = " << faces[f]->GetEid(k) << endl;
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
            xMap[i] = MemoryManager<StdRegions::StdPyrExp>::AllocateSharedPtr(basisKey_x, basisKey_y, basisKey_z);

         }

         SpatialDomains::PyrGeomSharedPtr geom =
         MemoryManager<SpatialDomains::PyrGeom>::AllocateSharedPtr(faces);
         geom->SetOwnData();


        if( bType_x_val < 10 ) {
            lpe = new LocalRegions::PyrExp( basisKey_x, basisKey_y, basisKey_z, geom );
        } else {
            cerr << "Implement the NodalTetExp!!!!!!" << endl;
            //lpe = new StdRegions::StdNodalTetExp( basisKey_x, basisKey_y, basisKey_z, NodalType );
            exit(1);
        }
    
        Array<OneD,NekDouble> x = Array<OneD,NekDouble>( Qx * Qy * Qz );
        Array<OneD,NekDouble> y = Array<OneD,NekDouble>( Qx * Qy * Qz );
        Array<OneD,NekDouble> z = Array<OneD,NekDouble>( Qx * Qy * Qz );
        
        lpe->GetCoords(x,y,z); 
    
        //----------------------------------------------
        // Define solution to be projected
        for(int n = 0; n < Qx * Qy * Qz; ++n) {
            solution[n]  = Pyr_sol( x[n], y[n], z[n], P, Q, R );
        }
        //----------------------------------------------
    }
           
    //---------------------------------------------
    // Project onto Expansion 
    lpe->FwdTrans( solution, lpe->UpdateCoeffs() );
    //---------------------------------------------
    
    //-------------------------------------------
    // Backward Transform Solution to get projected values
    lpe->BwdTrans( lpe->GetCoeffs(), lpe->UpdatePhys() );
    //-------------------------------------------  
    
    //--------------------------------------------
    // Calculate L_p error 
    cout << "L infinity error: " << lpe->Linf(solution) << endl;
    cout << "L 2 error:        " << lpe->L2  (solution) << endl;
    //--------------------------------------------
    
    //-------------------------------------------
    // Evaulate solution at x = y = z = 0  and print error
    Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);
     t[0] = -0.39;
     t[1] = -0.25;
     t[2] =  0.5;
    
    if( regionShape == StdRegions::ePyramid ) {
        solution[0] = Pyr_sol( t[0], t[1], t[2], P, Q, R );
    }
 
    NekDouble numericSolution = lpe->PhysEvaluate(t);
    cout << "Interpolation difference from actual solution at x = ( " << 
        t[0] << ", " << t[1] << ", " << t[2] << " ): " << numericSolution - solution[0] << endl;
    //-------------------------------------------
    
    return 0;
}

NekDouble Pyr_sol(NekDouble x, NekDouble y, NekDouble z, int P, int Q, int R) {
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

