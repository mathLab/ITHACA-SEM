#include <StdRegions/StdExpansion3D.h>
#include <LocalRegions/HexExp.h>

#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include <SpatialDomains/MeshComponents.h>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <iosfwd>

using namespace std;
using namespace boost;
using namespace Nektar;
using namespace Nektar::LibUtilities;
using namespace Nektar::LocalRegions;
using namespace Nektar::StdRegions;
using namespace Nektar::SpatialDomains;

/// Defines the test solution which excites all modes in the Hex.
NekDouble Hex_sol(NekDouble x, NekDouble y, NekDouble z,
                  int order1, int order2, int order3,
                  LibUtilities::BasisType btype1,
                  LibUtilities::BasisType btype2,
                  LibUtilities::BasisType btype3);

/// Specialised pow function returning 1.0 for exponents < 0.
static double  pow_loc(const double val, const int i)
{
  return (i < 0)? 1.0: pow(val,i);
}


int main(int argc, char *argv[])
{
    if( argc != 10 ) {
        cerr << "Usage: HexDemo Type_x Type_y Type_z numModes_x numModes_y "
                "numModes_z Qx Qy Qz" << endl;
        cerr << "Where type is an interger value which dictates the basis as:"
             << endl;
        cerr << "\t Ortho_A    = 1\n";
        cerr << "\t Modified_A = 4\n";
        cerr << "\t Fourier    = 7\n";
        cerr << "\t Lagrange   = 8\n";
        cerr << "\t Legendre   = 9\n";
        cerr << "\t Chebyshev  = 10\n";
        cerr << "\n\n" << "Example: " << argv[0] << " 4 4 4 3 3 3 5 5 5"
             << endl;
        cerr << endl;

        exit(1);
    }

    // Set up the region shape and expansion types
    int bType_x_val = atoi(argv[1]);
    int bType_y_val = atoi(argv[2]);
    int bType_z_val = atoi(argv[3]);

    BasisType   bType_x   = static_cast<BasisType>( bType_x_val );
    BasisType   bType_y   = static_cast<BasisType>( bType_y_val );
    BasisType   bType_z   = static_cast<BasisType>( bType_z_val );

    if( (bType_x_val == 13) || (bType_y_val == 13) || (bType_z_val == 13) )
    {
        bType_x =   LibUtilities::eOrtho_A;
        bType_y =   LibUtilities::eOrtho_B;
        bType_z =   LibUtilities::eOrtho_C;
    }

    if( (bType_x == eOrtho_B) || (bType_x == eModified_B) ) {
        NEKERROR(ErrorUtil::efatal,
                 "Basis 1 cannot be of type Ortho_B or Modified_B");
    }
    if( (bType_x == eOrtho_C) || (bType_x == eModified_C) ) {
        NEKERROR(ErrorUtil::efatal,
                 "Basis 1 cannot be of type Ortho_C or Modified_C");
    }
    if( (bType_y == eOrtho_C) || (bType_y == eModified_C) ) {
        NEKERROR(ErrorUtil::efatal,
                 "Basis 2 cannot be of type Ortho_C or Modified_C");
    }

    // Set up the number of quadrature points, order of bases, etc
    int xModes   = atoi(argv[4]);
    int yModes   = atoi(argv[5]);
    int zModes   = atoi(argv[6]);
    int Qx       = atoi(argv[7]);
    int Qy       = atoi(argv[8]);
    int Qz       = atoi(argv[9]);
    int P        = xModes - 1;
    int Q        = yModes - 1;
    int R        = zModes - 1;
    const int three = 3;

    Array<OneD, NekDouble> solution( Qx * Qy * Qz );

    LibUtilities::PointsType    Qtype_x, Qtype_y, Qtype_z;

    Array<OneD,NekDouble> x = Array<OneD,NekDouble>( Qx * Qy * Qz );
    Array<OneD,NekDouble> y = Array<OneD,NekDouble>( Qx * Qy * Qz );
    Array<OneD,NekDouble> z = Array<OneD,NekDouble>( Qx * Qy * Qz );

    if(bType_x != LibUtilities::eFourier)
    {
        Qtype_x = LibUtilities::eGaussLobattoLegendre;
    }
    else
    {
        Qtype_x = LibUtilities::eFourierEvenlySpaced;
    }

    if(bType_y != LibUtilities::eFourier)
    {
        Qtype_y = LibUtilities::eGaussLobattoLegendre;
    }
    else
    {
        Qtype_y = LibUtilities::eFourierEvenlySpaced;
    }

    if(bType_z != LibUtilities::eFourier)
    {
        Qtype_z = LibUtilities::eGaussLobattoLegendre;
    }
    else
    {
        Qtype_z = LibUtilities::eFourierEvenlySpaced;
    }


    //-----------------------------------------------
    // Define a 3D expansion based on basis definition
    StdRegions::StdExpansion3D *lhe = 0;

    // ////////////////////////////////////////////////////////////////
    // Set up Hexahedron vertex coordinates
    // VertexComponent (const int coordim, const int vid, double x,
    //    double y, double z)

    const int nVerts = 8;
    const double point[][3] = {
        {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
        {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}
    };

    // Populate the list of verts
    VertexComponentSharedPtr verts[8];
    for( int i = 0; i < nVerts; ++i ) {
        verts[i] = MemoryManager<VertexComponent>
                            ::AllocateSharedPtr(three,  i,   point[i][0],
                                                point[i][1], point[i][2]);
    }

    // ////////////////////////////////////////////////////////////////
    // Set up Hexahedron Edges
    // SegGeom (int id, const int coordim), EdgeComponent(id, coordim)
    const int nEdges = 12;
    const int vertexConnectivity[][2] = {
        {0,1}, {1,2}, {2,3}, {0,3}, {0,4}, {1,5},
        {2,6}, {3,7}, {4,5}, {5,6}, {6,7}, {4,7}
    };

    // Populate the list of edges
    SegGeomSharedPtr edges[nEdges];
    for( int i = 0; i < nEdges; ++i ) {
        VertexComponentSharedPtr vertsArray[2];
        for( int j = 0; j < 2; ++j ) {
            vertsArray[j] = verts[vertexConnectivity[i][j]];
        }
        edges[i] = MemoryManager<SegGeom>::
            AllocateSharedPtr( i, three, vertsArray);
    }

    // ////////////////////////////////////////////////////////////////
    // Set up Hexahedron faces
    const int nFaces = 6;
    const int edgeConnectivity[][4] = {
          {0,1,2,3}, {0,5,8,4}, {1,6,9,5},
          {2,7,10,6}, {3,7,11,4}, {8,9,10,11}
    };
    const bool isEdgeFlipped[][4] = {
          {0,0,0,1}, {0,0,1,1}, {0,0,1,1},
          {0,0,1,1}, {0,0,1,1}, {0,0,0,1}
    };

    // Populate the list of faces
    QuadGeomSharedPtr faces[nFaces];
    for( int i = 0; i < nFaces; ++i ) {
        SegGeomSharedPtr edgeArray[4];
        Orientation eorientArray[4];
        for( int j = 0; j < 4; ++j ) {
            edgeArray[j]    = edges[edgeConnectivity[i][j]];
            eorientArray[j] = isEdgeFlipped[i][j] ? eBackwards : eForwards;
        }
        faces[i] = MemoryManager<QuadGeom>::AllocateSharedPtr(i, edgeArray,
                                                              eorientArray);
    }


    const LibUtilities::PointsKey   pkey1( Qx, Qtype_x );
    const LibUtilities::PointsKey   pkey2( Qy, Qtype_y );
    const LibUtilities::PointsKey   pkey3( Qz, Qtype_z );

    const LibUtilities::BasisKey    bkey1( bType_x, xModes, pkey1 );
    const LibUtilities::BasisKey    bkey2( bType_y, yModes, pkey2 );
    const LibUtilities::BasisKey    bkey3( bType_z, zModes, pkey3 );

    Array<OneD, StdRegions::StdExpansion3DSharedPtr> xMap(3);
    for(int i = 0; i < 3; ++i) {
        xMap[i] = MemoryManager<StdRegions::StdHexExp>
                                ::AllocateSharedPtr(bkey1, bkey2, bkey3);
    }

    SpatialDomains::HexGeomSharedPtr geom
            = MemoryManager<SpatialDomains::HexGeom>::AllocateSharedPtr(faces);
    geom->SetOwnData();


    if( bType_x_val < 10 ) {
        lhe = new LocalRegions::HexExp( bkey1, bkey2, bkey3, geom );
    }

    lhe->GetCoords(x,y,z);

    //----------------------------------------------
    // Define solution to be projected
    for(int n = 0; n < Qx * Qy * Qz; ++n) {
        solution[n]  = Hex_sol( x[n], y[n], z[n], P, Q, R,
                                bType_x, bType_y, bType_z );
    }
    //----------------------------------------------

    //---------------------------------------------
    // Project onto Expansion
    lhe->FwdTrans( solution, lhe->UpdateCoeffs() );
    //---------------------------------------------

    //-------------------------------------------
    // Backward Transform Solution to get projected values
    lhe->BwdTrans( lhe->GetCoeffs(), lhe->UpdatePhys() );
    //-------------------------------------------

    //--------------------------------------------
    // Calculate L_p error
    cout << "L infinity error: " << lhe->Linf(solution) << endl;
    cout << "L 2 error:        " << lhe->L2  (solution) << endl;
    //--------------------------------------------

    //-------------------------------------------
    // Evaulate solution at x = y = z = 0  and print error
    Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);
    t[0] =  0.5;
    t[1] =  0.5;
    t[2] =  0.5;

    NekDouble numericSolution = lhe->PhysEvaluate(t);

    solution[0] = Hex_sol( t[0], t[1], t[2], P, Q, R,
                           bType_x, bType_y, bType_z );

    cout << "Solution = " << solution[0] << endl;
    cout << "Numeric  = " << numericSolution << endl;
    cout << "Interpolation difference from actual solution at x = ( "
         << t[0] << ", " << t[1] << ", " << t[2] << " ): "
         << numericSolution - solution[0] << endl;
    //-------------------------------------------

    return 0;
}



NekDouble Hex_sol(NekDouble x, NekDouble y, NekDouble z,
                  int order1, int order2, int order3,
                  LibUtilities::BasisType btype1,
                  LibUtilities::BasisType btype2,
                  LibUtilities::BasisType btype3 ) {
    int i,j,k;
    NekDouble sol = 0.0;

    int Nx = (btype1 == LibUtilities::eFourier ? order1/2 : order1);
    int Ny = (btype2 == LibUtilities::eFourier ? order2/2 : order2);
    int Nz = (btype3 == LibUtilities::eFourier ? order3/2 : order3);
    bool Fx = (btype1 == LibUtilities::eFourier);
    bool Fy = (btype2 == LibUtilities::eFourier);
    bool Fz = (btype3 == LibUtilities::eFourier);
    NekDouble a;

    for (i = 0; i < Nx; ++i)
    {
        for (j = 0; j < Ny; ++j)
        {
            for (k = 0; k < Nz; ++k)
            {
                a  = (Fx ? sin(M_PI*i*x) + cos(M_PI*i*x) : pow_loc(x,i));
                a *= (Fy ? sin(M_PI*j*y) + cos(M_PI*j*y) : pow_loc(y,j));
                a *= (Fz ? sin(M_PI*k*z) + cos(M_PI*k*z) : pow_loc(z,k));
                sol += a;
            }
        }
    }
    return sol;
}

