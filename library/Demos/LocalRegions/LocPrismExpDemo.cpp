#include "StdRegions/StdExpansion3D.h"
#include "LocalRegions/PrismExp.h"
#include "LocalRegions/LocalRegions.hpp"

#include "LibUtilities/Foundations/Foundations.hpp"
#include "LibUtilities/Foundations/Basis.h"


#include <algorithm>
#include <iostream>
#include <limits>
#include <cstdio>
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


int main(int argc, char *argv[])
 {
    if( argc != 19 ) {  // arg[0]  arg[1]  arg[2]  arg[3]   arg[4]     arg[5]     arg[6]  arg[7] arg[8] arg[9] arg[10] arg[11] arg[12]
        cerr << "Usage: LocPrismDemo Type_x  Type_y  Type_z numModes_x numModes_y numModes_z    Qx    Qy     Qz     x1     y1    z1"
                       //arg[13] arg[14] arg[15] arg[16] arg[17] arg[18] arg[19] arg[20] arg[21] arg[22] arg[23] arg[24] ...arg[27]
                      "x2     y2      z2     x3      y3     z3      x4     y4      z4    x5   y5   z5  x6  y6  z6" << endl;
        
        cerr << "Where type is an interger value which dictates the basis as:" << endl;
        for(int i=0; i<SIZE_PointsType; ++i)
        {
            cerr << setw(30) << PointsTypeMap[i] << " =" << i << endl;
        }
//                   NoPointsType =0
//             GaussGaussLegendre =1
//            GaussRadauMLegendre =2
//            GaussRadauPLegendre =3
//           GaussLobattoLegendre =4
//            GaussGaussChebyshev =5
//           GaussRadauMChebyshev =6
//           GaussRadauPChebyshev =7
//          GaussLobattoChebyshev =8
//         GaussRadauMAlpha0Beta1 =9
//         GaussRadauMAlpha0Beta2 =10
//         GaussRadauMAlpha1Beta0 =11
//         GaussRadauMAlpha2Beta0 =12
//               PolyEvenlySpaced =13
//            FourierEvenlySpaced =14
//                   NodalTriElec =15
//                 NodalTriFekete =16
//           NodalTriEvenlySpaced =17
//           NodalTetEvenlySpaced =18
//                   NodalTetElec =19
        cerr << "\t Nodal Tet (Electro) = 19    (3D Nodal Electrostatic Points on a Tetrahedron)\n";
//         cerr << "\n\n" << "Example: " << argv[0] << " 1 1 3 2 2 2 5 5 5 0.2 0.7 0.5 0.1 0.3 0.6 0.2 0.7 0.9 0.3 0.5 0.1 0.2 0.5 0.7" << endl;
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
    LibUtilities::PointsType  NodalType = LibUtilities::eNoPointsType;
    
    if( (bType_x_val == 19) || (bType_y_val == 19) || (bType_z_val == 19) )
    {
        bType_x =   LibUtilities::eOrtho_A;
        bType_y =   LibUtilities::eOrtho_B;
        bType_z =   LibUtilities::eOrtho_C;
        
        NodalType = LibUtilities::eNodalTetElec;
    }

    // Check to see that correct Expansions are used
    if( regionShape == StdRegions::ePrism ) 
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
    
    Array<OneD, NekDouble> solution( Qx * Qy * Qz );

    LibUtilities::PointsType    Qtype_x = eGaussLobattoLegendre;
    LibUtilities::PointsType    Qtype_y = eGaussRadauMAlpha1Beta0;
    LibUtilities::PointsType    Qtype_z = eGaussRadauMAlpha2Beta0;
        
    //-----------------------------------------------
    // Define a 3D expansion based on basis definition
    
    StdRegions::StdExpansion3D *lpe;
    
    if( regionShape == StdRegions::ePrism ) 
    {
          Array<OneD, NekDouble> coords(12);
          coords[0]    =   atof(argv[10]);
          coords[1]    =   atof(argv[11]);
          coords[2]    =   atof(argv[12]);
          coords[3]    =   atof(argv[13]);
          coords[4]    =   atof(argv[14]);
          coords[5]    =   atof(argv[15]);
          coords[6]    =   atof(argv[16]);
          coords[7]    =   atof(argv[17]);
          coords[8]    =   atof(argv[18]);
          coords[9]    =   atof(argv[19]);
          coords[10]   =   atof(argv[20]);
          coords[11]   =   atof(argv[21]);
          coords[12]   =   atof(argv[22]);
          coords[13]   =   atof(argv[23]);
          coords[14]   =   atof(argv[24]);
          coords[15]   =   atof(argv[25]);
          coords[16]   =   atof(argv[26]);
          coords[17]   =   atof(argv[27]);

          // Set up Prism vertex coordinates
          const int zero  = 0;
          const int one   = 1;
          const int two   = 2;
          const int three = 3;
          const int four  = 4;
          const int five  = 5;
          const int six   = 6;
          const int seven = 7;
          const int eight = 8;
          
          SpatialDomains::VertexComponentSharedPtr verts[6];
          //    VertexComponent (const int coordim, const int vid, double x, double y, double z)
          verts[0] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(three,zero, coords[0],coords[1], coords[2]);
          verts[1] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(three,one,  coords[3],coords[4], coords[5]);
          verts[2] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(three,two,  coords[6],coords[7], coords[8]);
          verts[3] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(three,three,coords[9],coords[10],coords[11]);
          verts[4] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(three,four,coords[12],coords[13],coords[14]);
          verts[5] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(three,four,coords[15],coords[16],coords[17]);

          // Set up Prism Edges
          //SegGeom (int id, const int coordim), EdgeComponent(id, coordim)
          SpatialDomains::SegGeomSharedPtr edges[9];
          edges[0] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(zero, three);
          edges[1] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(one, three);
          edges[2] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(two, three);
          edges[3] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(three, three);
          edges[4] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(four, three);
          edges[5] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(five, three);
          edges[6] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(six, three);
          edges[7] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(seven, three);
          edges[8] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(eight, three);
      

          StdRegions::EdgeOrientation edgeDir = StdRegions::eForwards;
          StdRegions::EdgeOrientation eorient[9];
          eorient[0] = edgeDir; 
          eorient[1] = edgeDir; 
          eorient[2] = edgeDir;
          eorient[3] = edgeDir;
          eorient[4] = edgeDir;
          eorient[5] = edgeDir;
          eorient[6] = edgeDir;
          eorient[7] = edgeDir;
          eorient[8] = edgeDir;

          // Set up Prism faces
          SpatialDomains::QuadGeomSharedPtr qfaces[3];
          qfaces[0] = MemoryManager<SpatialDomains::QuadGeom>::AllocateSharedPtr(zero, three);
          qfaces[1] = MemoryManager<SpatialDomains::QuadGeom>::AllocateSharedPtr(one, three);
          qfaces[2] = MemoryManager<SpatialDomains::QuadGeom>::AllocateSharedPtr(two, three);
          SpatialDomains::TriGeomSharedPtr tfaces[2];
          tfaces[0] = MemoryManager<SpatialDomains::TriGeom>::AllocateSharedPtr(one, three);
          tfaces[1] = MemoryManager<SpatialDomains::TriGeom>::AllocateSharedPtr(two, three);
          tfaces[2] = MemoryManager<SpatialDomains::TriGeom>::AllocateSharedPtr(three, three);
          tfaces[3] = MemoryManager<SpatialDomains::TriGeom>::AllocateSharedPtr(four, three);


          //TODO:  must check the face direction
          StdRegions::FaceOrientation faceDir = StdRegions::eDir1FwdDir1_Dir2FwdDir2;
          StdRegions::FaceOrientation forient[5];
          forient[0] = faceDir;
          forient[1] = faceDir;
          forient[2] = faceDir;
          forient[3] = faceDir;
          forient[4] = faceDir;

         SpatialDomains::PrismGeomSharedPtr geom = MemoryManager<SpatialDomains::PrismGeom>::AllocateSharedPtr(verts,edges,tfaces,
                                                                                                      qfaces, eorient, forient);
         geom->SetOwnData();
         const LibUtilities::PointsKey   pointsKey_x( Qx, Qtype_x );
         const LibUtilities::PointsKey   pointsKey_y( Qy, Qtype_y );
         const LibUtilities::PointsKey   pointsKey_z( Qz, Qtype_z );

         const LibUtilities::BasisKey    basisKey_x( bType_x, xModes, pointsKey_x );
         const LibUtilities::BasisKey    basisKey_y( bType_y, yModes, pointsKey_y );
         const LibUtilities::BasisKey    basisKey_z( bType_z, zModes, pointsKey_z );

        if( bType_x_val < 15 ) {
            lpe = new LocalRegions::PrismExp( basisKey_x, basisKey_y, basisKey_z, geom );
        } else {
            cerr << "Implement the NodalTetExp!!!!!!" << endl;
            //lpe = new StdRegions::StdNodalTetExp( basisKey_x, basisKey_y, basisKey_z, NodalType );
            exit(1);
        }
    
        Array<OneD,NekDouble> x = Array<OneD,NekDouble>( Qx * Qy * Qz );
        Array<OneD,NekDouble> y = Array<OneD,NekDouble>( Qx * Qy * Qz );
        Array<OneD,NekDouble> z = Array<OneD,NekDouble>( Qx * Qy * Qz );
        
        lpe->GetCoords(x,y,z); //TODO: must test PrismExp::GetCoords()
    
        //----------------------------------------------
        // Define solution to be projected
        for(int n = 0; n < Qx * Qy * Qz; ++n) {
            solution[n]  = Prism_sol( x[n], y[n], z[n], P, Q, R, bType_x, bType_y, bType_z );
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
    
    if( regionShape == StdRegions::ePrism ) {
        solution[0] = Prism_sol( t[0], t[1], t[2], P, Q, R, bType_x, bType_y, bType_z );
    }
 
    NekDouble numericSolution = lpe->PhysEvaluate(t);
    cout << "Interpolation difference from actual solution at x = ( " << 
        t[0] << ", " << t[1] << ", " << t[2] << " ): " << numericSolution - solution[0] << endl;
    //-------------------------------------------
    
    return 0;
}

NekDouble Prism_sol(NekDouble x, NekDouble y, NekDouble z, int P, int Q, int R, LibUtilities::BasisType bType_x,
                    LibUtilities::BasisType bType_y, LibUtilities::BasisType bType_z ) {
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

