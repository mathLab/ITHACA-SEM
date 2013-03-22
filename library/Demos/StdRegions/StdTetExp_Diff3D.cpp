
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdExpansion3D.h>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Basis.h>


#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;
using namespace Nektar;
using namespace boost;
using namespace Nektar::LibUtilities;
using namespace Nektar::StdRegions;


NekDouble Tet_sol(NekDouble x, NekDouble y, NekDouble z, int order1, int order2, int order3);
NekDouble Tet_Diff_Sol(NekDouble x, NekDouble y, NekDouble z, int P, int Q, int R, int dir);

// modification to deal with exact solution. Return 1 if integer < 0
static double  pow_loc(const double val, const int i)
{
  return (i < 0)? 1.0: pow(val,i);
}



// This routine projects a polynomial or trigonmetric functions which 
// has energy in all mdoes of the expansions and reports and error
int main(int argc, char *argv[]);

int main(int argc, char *argv[]) {
    if( argc != 10 ) {
        cerr << "Usage: StdTetExp_Diff3D Demo Type_x Type_y Type_z numModes_x numModes_y numModes_z Qx Qy Qz" << endl;
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
        cerr << "\n\n" << "Example: " << argv[0] << " 1 2 3 2 2 2 5 5 5" << endl;
        cerr << "\n\n" << "Example: " << argv[0] << " 1 2 3 5 5 5 8 8 8" << endl;
        cerr << endl;

        exit(1);
    }

    LibUtilities::ShapeType regionShape = LibUtilities::eTetrahedron;    
    
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
    if( regionShape == LibUtilities::eTetrahedron ) 
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
    
    //double *dx, *dy, *dz;
    
    Array<OneD, NekDouble> dx(Qx*Qy*Qz,0.0);
    Array<OneD, NekDouble> dy(Qx*Qy*Qz,0.0);
    Array<OneD, NekDouble> dz(Qx*Qy*Qz,0.0);
    
    Array<OneD, NekDouble> solution( Qx * Qy * Qz, 0.0 );
    Array<OneD, NekDouble> diff_solution_x( Qx * Qy * Qz, 0.0 ), diff_solution_y( Qx * Qy * Qz, 0.0 ), diff_solution_z( Qx * Qy * Qz, 0.0 );
    Array<OneD, NekDouble> derivatives( Qx * Qy * Qz, 0.0 );    

    LibUtilities::PointsType    Qtype_x = eGaussLobattoLegendre;
    LibUtilities::PointsType    Qtype_y = eGaussRadauMAlpha1Beta0;
    LibUtilities::PointsType    Qtype_z = eGaussRadauMAlpha2Beta0;
    
    Array<OneD,NekDouble> x = Array<OneD,NekDouble>( Qx * Qy * Qz, 0.0 );
    Array<OneD,NekDouble> y = Array<OneD,NekDouble>( Qx * Qy * Qz, 0.0 );
    Array<OneD,NekDouble> z = Array<OneD,NekDouble>( Qx * Qy * Qz, 0.0 );

    //-----------------------------------------------
    // Define a 3D expansion based on basis definition
    
    StdRegions::StdExpansion3D *ste;
    
    if( regionShape == LibUtilities::eTetrahedron ) 
    { 
        const LibUtilities::PointsKey   pointsKey_x( Qx, Qtype_x );
        const LibUtilities::PointsKey   pointsKey_y( Qy, Qtype_y );
        const LibUtilities::PointsKey   pointsKey_z( Qz, Qtype_z );

        const LibUtilities::BasisKey    basisKey_x( bType_x, xModes, pointsKey_x );
        const LibUtilities::BasisKey    basisKey_y( bType_y, yModes, pointsKey_y );
        const LibUtilities::BasisKey    basisKey_z( bType_z, zModes, pointsKey_z );
    
        if( bType_x_val < 10 ) 
        {
            ste = new StdRegions::StdTetExp( basisKey_x, basisKey_y, basisKey_z );
        } else 
        {
            cerr << "Implement the next line!!!!!!" << endl;
            //ste = new StdRegions::StdNodalTetExp( basisKey_x, basisKey_y, basisKey_z, NodalType );
            exit(1);
        } 
            
        ste->GetCoords(x,y,z);
    
        //----------------------------------------------
        // Define solution to be projected
        for(int n = 0; n < Qx * Qy * Qz; ++n) {
            solution[n]  = Tet_sol( x[n], y[n], z[n], P, Q, R );
            cout << "tet_solution["<<n<<"] = " << solution[n] << ",  x = " << x[n] << ", y = " << y[n] << ", z = " << z[n] << endl;
        }
        //----------------------------------------------
        
        
        //----------------------------------------------
        // Define the derivative solution
         for(int n = 0; n < Qx*Qy*Qz; ++n){
            diff_solution_x[n] = Tet_Diff_Sol( x[n], y[n], z[n], P, Q, R, 1);
            diff_solution_y[n] = Tet_Diff_Sol( x[n], y[n], z[n], P, Q, R, 2);
            diff_solution_z[n] = Tet_Diff_Sol( x[n], y[n], z[n], P, Q, R, 3);
        }
        //---------------------------------------------
    }
     
    //---------------------------------------------
    // Project onto Expansion 
    ste->FwdTrans( solution, ste->UpdateCoeffs()   );
    //---------------------------------------------
    
    
    //-------------------------------------------
    // Backward Transform Solution to get projected values
    ste->BwdTrans( ste->GetCoeffs( ), ste->UpdatePhys() );
    //---------------------------------------------
    
    
    //--------------------------------------------
    // Write solution 
//     FILE *outfile = fopen("ProjectFile3D.dat","w");
//     ste->WriteToFile(outfile);
//     fclose(outfile);
    //-------------------------------------------   
       
       
       
    //--------------------------------------------
     // Calculate L_p error 
    cout << "L infinity error: " << ste->Linf(solution) << endl;
    cout << "L 2 error:        " << ste->L2  (solution) << endl; 
     //--------------------------------------------  
    

    //--------------------------------------------
    // Taking the physical derivative and putting them into dx, dy, dz.
    ste->PhysDeriv( ste->GetPhys(), dx, dy, dz );        
    //--------------------------------------------     
         
    
    double error_x = 0, error_y=0, error_z=0;
    
    for( int n = 0; n < Qx*Qy*Qz; ++n ) {
    
        error_x += fabs(diff_solution_x[n] - dx[n]);
        error_y += fabs(diff_solution_y[n] - dy[n]);
        error_z += fabs(diff_solution_z[n] - dz[n]);
        cout << "diff_solution_x[n] = " << diff_solution_x[n] << ",    dx[n] = " << dx[n] << ",     " <<
                "diff_solution_y[n] = " << diff_solution_y[n] << ",    dy[n] = " << dy[n] << ",     " <<
                "diff_solution_z[n] = " << diff_solution_z[n] << ",    dz[n] = " << dz[n] << ",     " <<
        endl;
        
        
    }
    
    cout << "\n ******************************************** " << endl;
    cout << "L 1 error of derivatives X =  " <<  error_x << endl;
    cout << "L 1 error of derivatives Y =  " <<  error_y << endl;
    cout << "L 1 error of derivatives Z =  " <<  error_z << endl;
    cout << "******************************************** \n" << endl;
    //-------------------------------------------
    // Evaulate solution at x = y = z = 0  and print error
    Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);
    t[0] = -0.9; 
    t[1] = -0.75;
    t[2] = -0.85;
    
    if( regionShape ==  LibUtilities::eTetrahedron ) {
        solution[0] = Tet_sol( t[0], t[1], t[2], P, Q, R );        
    }
    
    NekDouble numericSolution = ste->PhysEvaluate(t);
    cout << "Interpolation difference from actual solution at x = ( " << 
        t[0] << ", " << t[1] << ", " << t[2] << " ): " << numericSolution - solution[0] << endl;
    //-------------------------------------------
      
   
    // Testing the physical evaluate(u_phys): projection on to the polynomial space given by the Tetrahedral basis function
    // The result of output should converge to the interpolation solution 
//     Array<OneD, const NekDouble> const& u_phys = ste->GetPhys();
//     Array<OneD, NekDouble> higherDegreeSolution( Qx * Qy * Qz, 0.0 );
//     for(int n = 0; n < Qx * Qy * Qz; ++n) {
//         higherDegreeSolution[n]  = Tet_sol( x[n], y[n], z[n], P, Q, R ); 
//         cout << "higherDegreeSolution["<<n<<"] = " << higherDegreeSolution[n] << ",  u_phys["<<n<<"] = " << u_phys[n]  << 
//             ",               x = " << x[n] << ", y = " << y[n] << ", z = " << z[n] << endl;
//         //cout << "         u_p["<<n<<"] = " <<      u_p[n] << endl;
//     }
//         
//     
    return 0;
}




NekDouble Tet_sol(NekDouble x, NekDouble y, NekDouble z, int P, int Q, int R) {
    NekDouble sol = 0;

    for(int i = 0; i <= P; ++i) {
        for(int j = 0; j <= Q - i; ++j) {
            for(int k = 0; k <= R - i - j; ++k) {
                sol += pow_loc(x,i) * pow_loc(y,j) * pow_loc(z,k);
            }
        }
    }

    return sol;
}

NekDouble Tet_Diff_Sol(NekDouble x, NekDouble y, NekDouble z, int P, int Q, int R, int dir){
    NekDouble sol = 0;
    for(int i = 0; i<= P; ++i){
        for(int j = 0; j<= Q-i; ++j){
            for(int k=0; k <= R-i-j; ++k){
                if(dir == 1 && i > 0) { sol += i * pow_loc(x,i-1)* pow_loc(y,j)   * pow_loc(z,k); }
                if(dir == 2 && j > 0) { sol += j * pow_loc(x,i)  * pow_loc(y,j-1) * pow_loc(z,k); }
                if(dir == 3 && k > 0) { sol += k * pow_loc(x,i)  * pow_loc(y,j)   * pow_loc(z,k-1); }
            }
        }
    }
    return sol;
}

