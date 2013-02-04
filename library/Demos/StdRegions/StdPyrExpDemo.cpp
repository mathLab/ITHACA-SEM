
#include <StdRegions/StdPyrExp.h>
#include <StdRegions/StdRegions.hpp>
#include <LibUtilities/Foundations/Foundations.hpp>

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
using namespace std;
using namespace Nektar;

                    
NekDouble Pyramid_sol(NekDouble x, NekDouble y, NekDouble z, int P, int Q, int R);                 
                  
using namespace Nektar::LibUtilities;
using namespace Nektar::StdRegions;


int main(int argc, char *argv[]) {
    if( argc != 10 ) {
        cerr << "Usage: Pyramid Demo Type_x Type_y Type_z numModes_x numModes_y numModes_z Qx Qy Qz" << endl;
        cerr << "Where type is an interger value which dictates the basis as:" << endl;
        cerr << "\t Ortho_A    = 1\n";
        cerr << "\t Ortho_B    = 2\n";
        cerr << "\t Ortho_C    = 3\n";
        cerr << "\t Modified_A = 4\n";
        cerr << "\t Modified_B = 5\n";
        cerr << "\t Modified_C = 6\n";
        cerr << "\t Fourier    = 7\n";
        cerr << "\t Lagrange   = 8\n";
        cerr << "\t Legendre   = 9\n";
        cerr << "\t Chebyshev  = 10\n";
        cerr << "\t Nodal Tet (Electro) = 13    (3D Nodal Electrostatic Points on a Tetrahedron)\n";
        cerr << "\n\n" << "Example: " << argv[0] << " 1 1 3 2 2 2 5 5 5" << endl;
        cerr << endl;

        exit(1);
    }

    StdRegions::ExpansionType regionShape = ePyramid;    
    
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
        if( (bType_z == LibUtilities::eOrtho_A) || (bType_z == LibUtilities::eModified_A) ) {
            NEKERROR(ErrorUtil::efatal, "Basis 3 cannot be of type Ortho_A or Modified_A");
        }
        if( (bType_z == LibUtilities::eOrtho_B) || (bType_z == LibUtilities::eModified_B) ) {
            NEKERROR(ErrorUtil::efatal, "Basis 3 cannot be of type Ortho_B or Modified_B");
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
    Array<OneD, NekDouble> diff_solution( Qx * Qy * Qz );

    Array<OneD,NekDouble> x = Array<OneD,NekDouble>( Qx * Qy * Qz );
    Array<OneD,NekDouble> y = Array<OneD,NekDouble>( Qx * Qy * Qz );
    Array<OneD,NekDouble> z = Array<OneD,NekDouble>( Qx * Qy * Qz );
    
    LibUtilities::PointsType    Qtype_x = eGaussLobattoLegendre;
    LibUtilities::PointsType    Qtype_y = eGaussLobattoLegendre;
    LibUtilities::PointsType    Qtype_z = eGaussRadauMAlpha2Beta0;
    
        
    //-----------------------------------------------
    // Define a 3D expansion based on basis definition
    
    StdRegions::StdExpansion *sPyrE;
    
    if( regionShape == StdRegions::ePyramid ) 
    { 
        const LibUtilities::PointsKey   pointsKey_x( Qx, Qtype_x );
        const LibUtilities::PointsKey   pointsKey_y( Qy, Qtype_y );
        const LibUtilities::PointsKey   pointsKey_z( Qz, Qtype_z );

        const LibUtilities::BasisKey    basisKey_x( bType_x, xModes, pointsKey_x );
        const LibUtilities::BasisKey    basisKey_y( bType_y, yModes, pointsKey_y );
        const LibUtilities::BasisKey    basisKey_z( bType_z, zModes, pointsKey_z );
    
        if( bType_x_val < 10 ) {
            sPyrE = new StdRegions::StdPyrExp( basisKey_x, basisKey_y, basisKey_z );
        } else {
            cerr << "Implement the next line!!!!!!" << endl;
            //sPyrE = new StdRegions::StdNodalTetExp( basisKey_x, basisKey_y, basisKey_z, NodalType );
            exit(1);
        } 
            
        sPyrE->GetCoords(x,y,z);
    
        //----------------------------------------------
        // Define solution to be projected
        for(int n = 0; n < Qx * Qy * Qz; ++n) {
            solution[n]  = Pyramid_sol( x[n], y[n], z[n], P, Q, R );
            //cout << "Pyramid_solution["<<n<<"] = " << solution[n] << ",  x = " << x[n] << ", y = " << y[n] << ", z = " << z[n] << endl;
                        
        }
        //----------------------------------------------
                
    }
       
             
    //---------------------------------------------
    // Project onto Expansion 
   sPyrE->FwdTrans( solution, sPyrE->UpdateCoeffs() );
    //---------------------------------------------
    
    //-------------------------------------------
    // Backward Transform Solution to get projected values
   // sPyrE->SetPhys( sPyrE->BwdTrans( sPyrE->GetCoeffs() ) );
   sPyrE->BwdTrans( sPyrE->GetCoeffs(), sPyrE->UpdatePhys() );
    //-------------------------------------------  
    
    //--------------------------------------------
    // Calculate L_p error 
   cout << "\n*****************************************************\n " << endl;
   cout << "L infinity error: " << sPyrE->Linf(solution) << endl;
   cout << "L 2 error:        " << sPyrE->L2  (solution) << endl;
   cout << "\n*****************************************************\n " << endl;
    //--------------------------------------------
        
    //-------------------------------------------
    // Evaulate solution at x = y = z = 0  and print error
    Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);

     t[0] = -0.39;
     t[1] = -0.25;
     t[2] = 0.5;
    
    if( regionShape == StdRegions::ePyramid ) {
        solution[0] = Pyramid_sol( t[0], t[1], t[2], P, Q, R); 
    }
    

    NekDouble numericSolution = sPyrE->PhysEvaluate(t);
    cout << "\n*****************************************************\n " << endl;
    cout << "Interpolation difference from actual solution at x = ( " << 
        t[0] << ", " << t[1] << ", " << t[2] << " ): " << numericSolution - solution[0] << endl;
    //-------------------------------------------
    cout << "numericSolution = " << numericSolution << endl;
    cout << "exact solution = " << solution[0] << endl;
    cout << "\n*****************************************************\n " << endl;
    
     // Testing the physical evaluate(u_phys): projection on to the polynomial space given by the prismatic basis function
    // The result of output should converge to the interpolation solution 
    //Array<OneD, const NekDouble> const& u_phys = sPyrE->GetPhys();
    Array<OneD, NekDouble> pyramid_solution( Qx * Qy * Qz, 0.0 );
    cout << setprecision(4);
    for(int n = 0; n < Qx * Qy * Qz; ++n) {
        pyramid_solution[n]  = Pyramid_sol( x[n], y[n], z[n], P, Q, R);   
//         cout << "pyramid_solution[" << setw(2)<<n<<"] = " << setw(8) << pyramid_solution[n] << ",  u_phys["<< setw(2) << n<<"] = " << setw(8) << u_phys[n]  << 
//             ",               x = " << setw(8) << x[n] << ", y = " << setw(8) << y[n] << ", z = " << setw(8) << z[n] << endl;
//         
        //cout << "         u_p["<<n<<"] = " <<      u_p[n] << endl;
    }
    cout << setprecision(9);
    
    
    return 0;
}



NekDouble Pyramid_sol(NekDouble x, NekDouble y, NekDouble z, int P, int Q, int R) {
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

