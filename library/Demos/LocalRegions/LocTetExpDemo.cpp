#include <StdRegions/StdExpansion3D.h>
#include <LocalRegions/TetExp.h>

#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/QuadGeom.h>

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <iosfwd>

// using namespace boost;

using namespace Nektar;
using namespace Nektar::LibUtilities;
using namespace Nektar::LocalRegions;
using namespace Nektar::StdRegions;
using namespace Nektar::SpatialDomains;
using namespace std;

NekDouble Tet_sol(NekDouble x, NekDouble y, NekDouble z, int P, int Q, int R);

int main(int argc, char *argv[])
{
    if( argc != 10 ) {
        cerr << "Usage: TetDemo Type_x Type_y Type_z numModes_x numModes_y "
                "numModes_z Qx Qy Qz" << endl;
        cerr << "Where type is an integer value which dictates the basis as:"
             << endl;
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
        cerr << "\t Nodal Tet (Electro) = 13    (3D Nodal Electrostatic Points"
                " on a Tetrahedron)\n";
        cerr << "\n\n" << "Example: " << argv[0] << " 4 4 4 3 3 3 5 5 5"
             << endl << endl;
        exit(1);
    }

    StdRegions::ExpansionType regionShape = StdRegions::eTetrahedron;
    int bType_x_val = atoi(argv[1]);
    int bType_y_val = atoi(argv[2]);
    int bType_z_val = atoi(argv[3]);

    BasisType   bType_x = static_cast<BasisType>( bType_x_val );
    BasisType   bType_y = static_cast<BasisType>( bType_y_val );
    BasisType   bType_z = static_cast<BasisType>( bType_z_val );

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
    LibUtilities::PointsType    Qtype_y = eGaussRadauMAlpha1Beta0;
    LibUtilities::PointsType    Qtype_z = eGaussRadauMAlpha2Beta0;

    //-----------------------------------------------
    // Define a 3D expansion based on basis definition

    StdRegions::StdExpansion3D *lte = 0;

    const int nVerts = 4;
    const double point[][3] = {
        {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}
    };

    // //////////////////////////////////////////////////////////////////////
    // Populate the list of verts
    // VertexComponent (const int coordim, const int vid, double x, double y,
    //   double z)
    VertexComponentSharedPtr verts[4];
    for(int i=0; i < nVerts; ++i){
        verts[i] =  MemoryManager<VertexComponent>::
        AllocateSharedPtr( three, i, point[i][0], point[i][1], point[i][2] );
    }

    // /////////////////////////////////////////////////////////////////////
    // Set up Tetrahedron Edges
    // SegGeom (int id, const int coordim), EdgeComponent(id, coordim)
    const int nEdges = 6;
    const int vertexConnectivity[][2] = {
        {0,1},{1,2},{0,2},{0,3},{1,3},{2,3}
    };

    // Populate the list of edges
    SegGeomSharedPtr edges[nEdges];
    for(int i=0; i < nEdges; ++i){
        VertexComponentSharedPtr vertsArray[2];
        for(int j=0; j<2; ++j){
            vertsArray[j] = verts[vertexConnectivity[i][j]];
        }
        edges[i] = MemoryManager<SegGeom>
                                    ::AllocateSharedPtr(i, three, vertsArray);
    }

    // //////////////////////////////////////////////////////////////////////
    // Set up Tetrahedron faces
    const int nFaces = 4;
    const int edgeConnectivity[][3] = {
        {0,1,2}, {0,4,3}, {1,5,4}, {2,5,3}
    };
    const bool   isEdgeFlipped[][3] = {
        {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}
    };

    // Populate the list of faces
    TriGeomSharedPtr faces[nFaces];
    for(int i=0; i < nFaces; ++i)
    {
        SegGeomSharedPtr edgeArray[3];
        Orientation eorientArray[3];
        for(int j=0; j < 3; ++j)
        {
            edgeArray[j] = edges[edgeConnectivity[i][j]];
            eorientArray[j] = isEdgeFlipped[i][j] ? eBackwards : eForwards;
        }
        faces[i] = MemoryManager<TriGeom>
                            ::AllocateSharedPtr(i, edgeArray, eorientArray);
    }

    const LibUtilities::PointsKey   pointsKey_x( Qx, Qtype_x );
    const LibUtilities::PointsKey   pointsKey_y( Qy, Qtype_y );
    const LibUtilities::PointsKey   pointsKey_z( Qz, Qtype_z );
    const LibUtilities::BasisKey    basisKey_x( bType_x, xModes, pointsKey_x );
    const LibUtilities::BasisKey    basisKey_y( bType_y, yModes, pointsKey_y );
    const LibUtilities::BasisKey    basisKey_z( bType_z, zModes, pointsKey_z );

    Array<OneD, StdRegions::StdExpansion3DSharedPtr> xMap(3);
    for(int i=0; i < 3; ++i)
    {
        xMap[i] = MemoryManager<StdRegions::StdTetExp>
                        ::AllocateSharedPtr(basisKey_x, basisKey_y, basisKey_z);
    }

    SpatialDomains::TetGeomSharedPtr geom =
    MemoryManager<SpatialDomains::TetGeom>::AllocateSharedPtr(faces);
    geom->SetOwnData();

    if( bType_x_val < 10 )
    {
        lte = new LocalRegions::TetExp( basisKey_x, basisKey_y,
                                        basisKey_z, geom );
    }
    else
    {
        cerr << "Implement the NodalTetExp!!!!!!" << endl;
        //lte = new StdRegions::StdNodalTetExp( basisKey_x, basisKey_y,
        //                                      basisKey_z, NodalType );
        exit(1);
    }

    Array<OneD,NekDouble> x = Array<OneD,NekDouble>( Qx * Qy * Qz );
    Array<OneD,NekDouble> y = Array<OneD,NekDouble>( Qx * Qy * Qz );
    Array<OneD,NekDouble> z = Array<OneD,NekDouble>( Qx * Qy * Qz );

    lte->GetCoords(x,y,z);

    //----------------------------------------------
    // Define solution to be projected
    for(int n = 0; n < Qx * Qy * Qz; ++n) {
        solution[n]  = Tet_sol( x[n], y[n], z[n], P, Q, R );
    }
    //----------------------------------------------

    //---------------------------------------------
    // Project onto Expansion
    lte->FwdTrans( solution, lte->UpdateCoeffs() );
    //---------------------------------------------

    //-------------------------------------------
    // Backward Transform Solution to get projected values
    lte->BwdTrans( lte->GetCoeffs(), lte->UpdatePhys() );
    //-------------------------------------------

    //--------------------------------------------
    // Calculate L_p error
    cout << "L infinity error: " << lte->Linf(solution) << endl;
    cout << "L 2 error:        " << lte->L2  (solution) << endl;
    //--------------------------------------------

    //-------------------------------------------
    // Evaulate solution at x = y = z = 0  and print error
    Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);

    if( regionShape == StdRegions::eTetrahedron ) {
        solution[0] = Tet_sol( t[0], t[1], t[2], P, Q, R );
    }

    t[0] = -0.0;
    t[1] = -1.0;
    t[2] = -1.0;
    NekDouble numericSolution = lte->PhysEvaluate(t);
    cout << "Numeric Solution: " << numericSolution << endl;
    cout << "Actual Solution:  " << solution[0] << endl;
    cout << "Interpolation difference from actual solution at x = ( "
         << t[0] << ", " << t[1] << ", " << t[2] << " ): "
         << numericSolution - solution[0] << endl;
    //-------------------------------------------

    return 0;
}


NekDouble Tet_sol(NekDouble x, NekDouble y, NekDouble z, int P, int Q, int R) {
    NekDouble sol = 0;

    for(int i = 0; i <= P; ++i) {
        for(int j = 0; j <= Q - j; ++j) {
            for(int k = 0; k <= R - i - j; ++k) {
                sol += pow(x,i) * pow(y,j) * pow(z,k);
            }
        }
    }

    return sol;
}

