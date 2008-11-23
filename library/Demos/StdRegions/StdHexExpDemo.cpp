#include <algorithm>
#include <iostream>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cmath>
using namespace std;

#include <StdRegions/StdExpUtil.h>
#include <StdRegions/StdHexExp.h>
#include "StdRegions/StdExpansion3D.h"
#include "StdRegions/StdRegions.hpp"
#include "LibUtilities/Foundations/Foundations.hpp"
#include "LibUtilities/Foundations/Basis.h"


using namespace Nektar;

NekDouble Hex_sol(NekDouble x, NekDouble y, NekDouble z, int order1, int order2, int order3,
                  LibUtilities::BasisType bType_x, LibUtilities::BasisType bType_y, LibUtilities::BasisType bType_z);

using namespace Nektar::LibUtilities;
using namespace Nektar::StdRegions;

// This routine projects a polynomial or trigonmetric functions which 
// has energy in all mdoes of the expansions and reports and error

int main(int argc, char *argv[]) {
    if( argc != 10 ) {
        cerr << "Usage: HexDemo Type_x Type_y Type_z numModes_x numModes_y numModes_z Qx Qy Qz" << endl;
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
        cerr << "\n\n" << "Example: " << argv[0] << " 1 2 3 3 3 3 5 5 5" << endl;
        cerr << endl;

        exit(1);
    }

    StdRegions::ExpansionType regionShape = eHexahedron;    
    
    int bType_x_val = atoi(argv[1]);
    int bType_y_val = atoi(argv[2]);
    int bType_z_val = atoi(argv[3]);
    
    LibUtilities::BasisType   bType_x = static_cast<LibUtilities::BasisType>( bType_x_val );
    LibUtilities::BasisType   bType_y = static_cast<LibUtilities::BasisType>( bType_y_val );
    LibUtilities::BasisType   bType_z = static_cast<LibUtilities::BasisType>( bType_z_val );

    // Check to see that correct Expansions are used
    if( regionShape == StdRegions::eHexahedron ) 
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
    Array<OneD, NekDouble> diff_solution( Qx * Qy * Qz );
    
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
    
    StdRegions::StdExpansion *she;
    
    if( regionShape == StdRegions::eHexahedron ) 
    { 
        const LibUtilities::PointsKey   pointsKey_x( Qx, Qtype_x );
        const LibUtilities::PointsKey   pointsKey_y( Qy, Qtype_y );
        const LibUtilities::PointsKey   pointsKey_z( Qz, Qtype_z );

        const LibUtilities::BasisKey    basisKey_x( bType_x, xModes, pointsKey_x );
        const LibUtilities::BasisKey    basisKey_y( bType_y, yModes, pointsKey_y );
        const LibUtilities::BasisKey    basisKey_z( bType_z, zModes, pointsKey_z );
    
        she = new StdRegions::StdHexExp( basisKey_x, basisKey_y, basisKey_z );
            
        she->GetCoords(x,y,z);
    
        //----------------------------------------------
        // Define solution to be projected
        for(int n = 0; n < Qx * Qy * Qz; ++n) {
            solution[n]  = Hex_sol( x[n], y[n], z[n], P, Q, R, bType_x, bType_y, bType_z );
//            cout << "Hex_solution["<<n<<"] = " << solution[n] << ",  x = " << x[n] << ", y = " << y[n] << ", z = " << z[n] << endl;
                        
        }
        //----------------------------------------------
                
    }
       
                     
    //---------------------------------------------
    // Project onto Expansion 
   she->FwdTrans( solution, she->UpdateCoeffs() );
    //---------------------------------------------
    
    //-------------------------------------------
    // Backward Transform Solution to get projected values
   // ste->SetPhys( ste->BwdTrans( ste->GetCoeffs() ) );
   she->BwdTrans( she->GetCoeffs(), she->UpdatePhys() );
    //-------------------------------------------  
    
    
    //--------------------------------------------
    // Calculate L_p error 
   cout << "\n ********************************************************** \n" << endl; 
   cout << "L infinity error: " << she->Linf(solution) << endl;
   cout << "L 2 error:        " << she->L2  (solution) << endl;
   cout << "\n ********************************************************** \n" << endl; 
    //--------------------------------------------
                        
    //-------------------------------------------
    // Evaulate solution at x = y = z = 0  and print error
    Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);
//     t[0] = -0.9;
//     t[1] = -0.75;
//     t[2] = -0.85;

     t[0] = -0.39;
     t[1] = -0.25;
     t[2] = 0.5;
    
    if( regionShape == StdRegions::eHexahedron ) {
        solution[0] = Hex_sol( t[0], t[1], t[2], P, Q, R, bType_x, bType_y, bType_z ); 
    }

     
        
    cout << "\n ********************************************************** \n" << endl;    
    NekDouble numericSolution = she->PhysEvaluate(t);
    cout << "Interpolation difference from actual solution at x = ( " << 
        t[0] << ", " << t[1] << ", " << t[2] << " ): " << numericSolution - solution[0] << endl;
    //-------------------------------------------
    cout << "numericSolution = " << numericSolution << endl;
    cout << "exact solution = " << solution[0] << endl;        
        
    cout << "\n ********************************************************** \n" << endl;    
    
    
    
    // To test this:  Turn it off the comment.
    // Testing the physical evaluate(u_phys): projection on to the 
    // polynomial space given by the Hexahedral basis function
    // The result of output should converge to the interpolation solution 
    
    // Array<OneD, const NekDouble> const& u_phys = she->GetPhys();
    // Array<OneD, NekDouble> hex_solution( Qx * Qy * Qz, 0.0 );
    // for(int n = 0; n < Qx * Qy * Qz; ++n) {
    //   hex_solution[n]  = Hex_sol( x[n], y[n], z[n], P, Q, R, bType_x, bType_y, bType_z );   
    //   cout << "hex_solution["<<n<<"] = " << hex_solution[n] << ",  u_phys["<<n<<"] = " << u_phys[n]  << 
    //   ",               x = " << x[n] << ", y = " << y[n] << ", z = " << z[n] << endl;
    //   //cout << "         u_p["<<n<<"] = " <<      u_p[n] << endl;
    // }
    
        
    return 0;
}  


NekDouble Hex_sol(NekDouble x, NekDouble y, NekDouble z, int P, int Q, int R, LibUtilities::BasisType bType_x,
                  LibUtilities::BasisType bType_y, LibUtilities::BasisType bType_z ) {
    NekDouble sol = 0;
 
    // case 1
    if( (bType_x != LibUtilities::eFourier) && (bType_y != LibUtilities::eFourier) && (bType_z != LibUtilities::eFourier)  )
    {
        for(int i = 0; i <= P; ++i) {
            for(int j = 0; j <= Q; ++j) {
                for(int k = 0; k <= R; ++k) {
                    sol += pow(x,i) * pow(y,j) * pow(z,k);
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
    if((bType_x != LibUtilities::eFourier) && (bType_y == LibUtilities::eFourier) && (bType_z == LibUtilities::eFourier))
    {
        for(int i = 0; i <= P; ++i) {
            for(int j = 0; j <= Q/2; ++j) {
                for(int k = 0; k <= R/2; ++k) {
                    sol += pow(x,i)*sin(M_PI*j*y)* sin(M_PI*k*z) + pow(x,i)*sin(M_PI*j*y)* cos(M_PI*k*z)
                         + pow(x,i)*cos(M_PI*j*y)* sin(M_PI*k*z) + pow(x,i)*cos(M_PI*j*y)* cos(M_PI*k*z);
                }
            }
        }
    }else // case 4
    if((bType_x == LibUtilities::eFourier) && (bType_y != LibUtilities::eFourier) && (bType_z != LibUtilities::eFourier))
    {
        for(int i = 0; i <= P/2; ++i) {
            for(int j = 0; j <= Q; ++j) {
                for(int k = 0; k <= R; ++k) {
                    sol += sin(M_PI*i*x)* pow(y,j) * pow(z,k) + cos(M_PI*i*x)* pow(y,j) * pow(z,k);                      
                }
            }
        }
    }else  // case 5
    if((bType_x == LibUtilities::eFourier) && (bType_y == LibUtilities::eFourier) && (bType_z != LibUtilities::eFourier))
    {
        for(int i = 0; i <= P/2; ++i) {
            for(int j = 0; j <= Q/2; ++j) {
                for(int k = 0; k <= R; ++k) {
                    sol += sin(M_PI*i*x)*sin(M_PI*j*y)* pow(z,k) + sin(M_PI*i*x)*cos(M_PI*j*y)* pow(z,k)
                         + cos(M_PI*i*x)*sin(M_PI*j*y)* pow(z,k) + cos(M_PI*i*x)*cos(M_PI*j*y)* pow(z,k);
                }
            } 
        }
    }else  // case 6
    if((bType_x == LibUtilities::eFourier) && (bType_y == LibUtilities::eFourier) && (bType_z == LibUtilities::eFourier))
    {
        for(int i = 0; i <= P/2; ++i) {
            for(int j = 0; j <= Q/2 ; ++j) {
                for(int k = 0; k <= R/2 ; ++k) {
                    sol += sin(M_PI*i*x)*sin(M_PI*j*y)*sin(M_PI*k*z) + sin(M_PI*i*x)*sin(M_PI*j*y)*cos(M_PI*k*z)
                         + sin(M_PI*i*x)*cos(M_PI*j*y)*sin(M_PI*k*z) + cos(M_PI*i*x)*sin(M_PI*j*y)*sin(M_PI*k*z)
                         + cos(M_PI*i*x)*cos(M_PI*j*y)*sin(M_PI*k*z) + cos(M_PI*i*x)*sin(M_PI*j*y)*cos(M_PI*k*z)
                         + sin(M_PI*i*x)*cos(M_PI*j*y)*cos(M_PI*k*z) + cos(M_PI*i*x)*cos(M_PI*j*y)*cos(M_PI*k*z);
                }
            }
        }
    }
    
    return sol;
}


