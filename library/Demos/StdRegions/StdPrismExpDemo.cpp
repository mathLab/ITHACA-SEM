
#include <StdRegions/StdRegions.hpp>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <StdRegions/StdPrismExp.h>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>


#include <iostream>     
#include <cstdlib>
#include <cmath>    
#include <iomanip>  
using namespace std;
using namespace Nektar;


NekDouble Prism_sol(NekDouble x, NekDouble y, NekDouble z, int order1, int order2, int order3,
                  LibUtilities::BasisType bType_x, LibUtilities::BasisType bType_y, LibUtilities::BasisType bType_z);

using namespace Nektar::LibUtilities;
using namespace Nektar::StdRegions;


int main(int argc, char *argv[]) {
    if( argc != 10 ) {
        cerr << "Usage: StdPrismExpDemo Type_x Type_y Type_z numModes_x numModes_y numModes_z Qx Qy Qz" << endl;
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
        cerr << "\n\n" << "Example: " << argv[0] << " 1 1 2 3 3 3 5 5 5" << endl;
        cerr << endl;

        exit(1);
    }

    LibUtilities::ShapeType regionShape = LibUtilities::ePrism;    
    
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
    if( regionShape == LibUtilities::ePrism ) 
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
      Qtype_z = LibUtilities::eGaussRadauMAlpha1Beta0; 
  }
  else
  {
      Qtype_z = LibUtilities::eFourierEvenlySpaced;
  }
        
    //-----------------------------------------------
    // Define a 3D expansion based on basis definition
    
    StdRegions::StdExpansion *spe;
    
    if( regionShape == LibUtilities::ePrism ) 
    { 
        const LibUtilities::PointsKey   pointsKey_x( Qx, Qtype_x );
        const LibUtilities::PointsKey   pointsKey_y( Qy, Qtype_y );
        const LibUtilities::PointsKey   pointsKey_z( Qz, Qtype_z );

        const LibUtilities::BasisKey    basisKey_x( bType_x, xModes, pointsKey_x );
        const LibUtilities::BasisKey    basisKey_y( bType_y, yModes, pointsKey_y );
        const LibUtilities::BasisKey    basisKey_z( bType_z, zModes, pointsKey_z );
    
        if( bType_x_val < 10 ) {
            spe = new StdRegions::StdPrismExp( basisKey_x, basisKey_y, basisKey_z );
        } else {
            cerr << "Implement the next line!!!!!!" << endl;
            //spe = new StdRegions::StdNodalTetExp( basisKey_x, basisKey_y, basisKey_z, NodalType );
            exit(1);
        } 
            
        spe->GetCoords(x,y,z);
    
        //----------------------------------------------
        // Define solution to be projected
        for(int n = 0; n < Qx * Qy * Qz; ++n) {
            solution[n]  = Prism_sol( x[n], y[n], z[n], P, Q, R, bType_x, bType_y, bType_z );
            //cout << "Prism_solution["<<n<<"] = " << solution[n] << ",  x = " << x[n] << ", y = " << y[n] << ", z = " << z[n] << endl;
                        
        }
        //----------------------------------------------
                
    }
                       
    //---------------------------------------------
    // Project onto Expansion 
   spe->FwdTrans( solution, spe->UpdateCoeffs() );
    //---------------------------------------------
    
    //-------------------------------------------
    // Backward Transform Solution to get projected values
   spe->BwdTrans( spe->GetCoeffs(), spe->UpdatePhys() );
    //-------------------------------------------  
    
    //--------------------------------------------
    // Calculate L_p error 
   cout << "\n*****************************************************\n " << endl;
   cout << "L infinity error: " << spe->Linf(solution) << endl;
   cout << "L 2 error:        " << spe->L2  (solution) << endl;
   cout << "\n*****************************************************\n " << endl;
    //--------------------------------------------
        
    //-------------------------------------------
    // Evaulate solution at x = y = z = 0  and print error
    Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);

    t[0] = -0.39;
    t[1] = -0.25;
    t[2] = 0.5;
    
    if(regionShape == LibUtilities::ePrism ) {
        solution[0] = Prism_sol( t[0], t[1], t[2], P, Q, R, bType_x, bType_y, bType_z ); 
    }
    

    NekDouble numericSolution = spe->PhysEvaluate(t);
    cout << "\n*****************************************************\n " << endl;
    cout << "Interpolation difference from actual solution at x = ( " << 
        t[0] << ", " << t[1] << ", " << t[2] << " ): " << numericSolution - solution[0] << endl;
    //-------------------------------------------
    cout << "numericSolution = " << numericSolution << endl;
    cout << "exact solution = " << solution[0] << endl;
    cout << "\n*****************************************************\n " << endl;
    
     // Testing the physical evaluate(u_phys): projection on to the polynomial space given by the prismatic basis function
    // The result of output should converge to the interpolation solution 
    //Array<OneD, const NekDouble> const& u_phys = spe->GetPhys();
    Array<OneD, NekDouble> prism_solution( Qx * Qy * Qz, 0.0 );
    cout << setprecision(4);
    for(int n = 0; n < Qx * Qy * Qz; ++n) {
        prism_solution[n]  = Prism_sol( x[n], y[n], z[n], P, Q, R, bType_x, bType_y, bType_z );   
/*        cout << "prism_solution[" << setw(2)<<n<<"] = " << setw(8) << prism_solution[n] << ",  u_phys["<< setw(2) << n<<"] = " << setw(8) << u_phys[n]  << 
            ",               x = " << setw(8) << x[n] << ", y = " << setw(8) << y[n] << ", z = " << setw(8) << z[n] << endl;*/
        
        //cout << "         u_p["<<n<<"] = " <<      u_p[n] << endl;
    }
    cout << setprecision(9);
    
    
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
