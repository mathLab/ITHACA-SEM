
#include <StdRegions/StdTetExp.h>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Basis.h>


#include <iostream>
#include <cstdlib>
#include <cmath>
using namespace std;


using namespace Nektar;


NekDouble Tet_sol(NekDouble x, NekDouble y, NekDouble z, int order1, int order2, int order3);

// using namespace boost;
using namespace Nektar::LibUtilities;
using namespace Nektar::StdRegions;


int main(int argc, char *argv[])
 {
    if( argc != 10 ) {
        cerr << "Usage: StdTetDemo Type_x Type_y Type_z numModes_x numModes_y numModes_z Qx Qy Qz" << endl;
        cerr << "Where type is an interger value which dictates the basis as:" << endl;
        cerr << "\t Ortho_A    = 1\n";
        cerr << "\t Ortho_B    = 2\n";
        cerr << "\t Ortho_C    = 3\n";
        cerr << "\t Modified_A = 4\n";
        cerr << "\t Modified_B = 5\n";
        cerr << "\t Modified_C = 6\n";
        cerr << "\t Nodal Tet (Electro) = 7    (3D Nodal Electrostatic Points on a Tetrahedron)\n";
        cerr << "\n\n" << "Example: " << argv[0] << " 1 2 3 2 2 2 5 5 5" << endl;
        cerr << endl;

        exit(1);
    }

    StdRegions::ExpansionType regionShape = StdRegions::eTetrahedron;
    int bType_x_val = atoi(argv[1]);
    int bType_y_val = atoi(argv[2]);
    int bType_z_val = atoi(argv[3]);
    
    LibUtilities::BasisType   bType_x = static_cast<LibUtilities::BasisType>( bType_x_val );
    LibUtilities::BasisType   bType_y = static_cast<LibUtilities::BasisType>( bType_y_val );
    LibUtilities::BasisType   bType_z = static_cast<LibUtilities::BasisType>( bType_z_val );
    
    if( (bType_x_val == 7) || (bType_y_val == 7) || (bType_z_val == 7) )
    {
        bType_x =   LibUtilities::eOrtho_A;
        bType_y =   LibUtilities::eOrtho_B;
        bType_z =   LibUtilities::eOrtho_C;
    }

    // Check to see that correct Expansions are used
    if( regionShape == StdRegions::eTetrahedron ) 
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
    
    StdRegions::StdExpansion *ste;
    
    if( regionShape == StdRegions::eTetrahedron ) 
    {
        const LibUtilities::PointsKey   pointsKey_x( Qx, Qtype_x );
        const LibUtilities::PointsKey   pointsKey_y( Qy, Qtype_y );
        const LibUtilities::PointsKey   pointsKey_z( Qz, Qtype_z );
        
        const LibUtilities::BasisKey    basisKey_x( bType_x, xModes, pointsKey_x );
        const LibUtilities::BasisKey    basisKey_y( bType_y, yModes, pointsKey_y );
        const LibUtilities::BasisKey    basisKey_z( bType_z, zModes, pointsKey_z );
        
        if( bType_x_val < 7 ) 
        {
            ste = new StdTetExp( basisKey_x, basisKey_y, basisKey_z );
        }
        else {
            cerr << "Implement the next line!!!!!!" << endl;
            //ste = new StdRegions::StdNodalTetExp( basisKey_x, basisKey_y, basisKey_z, NodalType );
            exit(1);
        }
        
        Array<OneD,NekDouble> x = Array<OneD,NekDouble>( Qx * Qy * Qz );
        Array<OneD,NekDouble> y = Array<OneD,NekDouble>( Qx * Qy * Qz );
        Array<OneD,NekDouble> z = Array<OneD,NekDouble>( Qx * Qy * Qz );
        
        ste->GetCoords(x,y,z);
        
        //----------------------------------------------
        // Define solution to be projected
        for(int n = 0; n < Qx * Qy * Qz; ++n) {
            solution[n]  = Tet_sol( x[n], y[n], z[n], P, Q, R );
        }
        //----------------------------------------------
    }
    
    //---------------------------------------------
    // Project onto Expansion 
    ste->FwdTrans( solution, ste->UpdateCoeffs() );
    //---------------------------------------------
    
    //-------------------------------------------
    // Backward Transform Solution to get projected values
    ste->BwdTrans( ste->GetCoeffs(), ste->UpdatePhys() );
    //-------------------------------------------  
    
    //--------------------------------------------
    // Calculate L_p error 
    cout << "L infinity error: " << ste->Linf(solution) << endl;
    cout << "L 2 error:        " << ste->L2  (solution) << endl;
    //--------------------------------------------
    
    //-------------------------------------------
    // Evaulate solution at x = y = z = 0  and print error
    Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);
    t[0] = -0.9;
    t[1] = -0.75;
    t[2] = -0.85;
    
    solution[0] = Tet_sol( t[0], t[1], t[2], P, Q, R ); 
 
    NekDouble numericSolution = ste->PhysEvaluate(t);
    cout << "Interpolation error at x = (" << 
        t[0] << ", " << t[1] << ", " << t[2] << " ): " << numericSolution - solution[0] << endl;
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

