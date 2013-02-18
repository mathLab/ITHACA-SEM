
#include <StdRegions/StdPyrExp.h>
#include <StdRegions/StdRegions.hpp>
#include <LibUtilities/Foundations/Foundations.hpp>

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;
using namespace Nektar;



                                 
NekDouble Pyramid_sol(NekDouble x, NekDouble y, NekDouble z, int P, int Q, int R);                                 
NekDouble Pyramid_Diff_Sol(NekDouble x, NekDouble y, NekDouble z, int P, int Q, int R, int dir);                                

// modification to deal with exact solution. Return 1 if integer < 0
static double  pow_loc(const double val, const int i)
{
  return (i < 0)? 1.0: pow(val,i);
}

using namespace Nektar::LibUtilities;
using namespace Nektar::StdRegions;


int main(int argc, char *argv[]) {
    if( argc != 10 ) {
        cerr << "Usage: Pyramid_Diff3D Demo Type_x Type_y Type_z numModes_x numModes_y numModes_z Qx Qy Qz" << endl;
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
    Array<OneD, NekDouble> diff_solution_x( Qx * Qy * Qz, 0.0 ), diff_solution_y( Qx * Qy * Qz, 0.0 ), diff_solution_z( Qx * Qy * Qz, 0.0 );
    
    Array<OneD, NekDouble> dx(Qx*Qy*Qz,0.0);
    Array<OneD, NekDouble> dy(Qx*Qy*Qz,0.0);
    Array<OneD, NekDouble> dz(Qx*Qy*Qz,0.0);
    
    Array<OneD,NekDouble> x = Array<OneD,NekDouble>( Qx * Qy * Qz );
    Array<OneD,NekDouble> y = Array<OneD,NekDouble>( Qx * Qy * Qz );
    Array<OneD,NekDouble> z = Array<OneD,NekDouble>( Qx * Qy * Qz );
    
    LibUtilities::PointsType    Qtype_x = eGaussLobattoLegendre;
    LibUtilities::PointsType    Qtype_y = eGaussLobattoLegendre;
    LibUtilities::PointsType    Qtype_z = eGaussRadauMAlpha2Beta0;

    //-----------------------------------------------
    // Define a 3D expansion based on basis definition
    
    StdRegions::StdExpansion3D *sPyrE;
    
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
            //cout << "Prism_solution["<<n<<"] = " << solution[n] << ",  x = " << x[n] << ", y = " << y[n] << ", z = " << z[n] << endl;
                        
        }
        //----------------------------------------------
        
       //----------------------------------------------
        // Define the derivative solution
         for(int n = 0; n < Qx*Qy*Qz; ++n){
            diff_solution_x[n] =  Pyramid_Diff_Sol( x[n], y[n], z[n], P, Q, R, 1);
            diff_solution_y[n] =  Pyramid_Diff_Sol( x[n], y[n], z[n], P, Q, R, 2);
            diff_solution_z[n] =  Pyramid_Diff_Sol( x[n], y[n], z[n], P, Q, R, 3);
        }
        //---------------------------------------------
                
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
        
        
 //--------------------------------------------
    // Taking the physical derivative and putting them into dx, dy, dz.
    sPyrE->PhysDeriv( sPyrE->GetPhys(), dx, dy, dz );        
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
        
        
    //-------------------------------------------
    // Evaulate solution at x = y = z = 0  and print error
    Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);
    
    t[0] = -0.19;
    t[1] = -0.5;
    t[2] = 0.25;


 if( regionShape == StdRegions::ePyramid ) {
        diff_solution_x[0] = Pyramid_Diff_Sol( t[0], t[1], t[2], P, Q, R, 1 );  
    }
    
 
    
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


NekDouble Pyramid_Diff_Sol(NekDouble x, NekDouble y, NekDouble z, int P, int Q, int R, int dir){
    NekDouble sol = 0;
    for(int i = 0; i<= P; ++i){
        for(int j = 0; j<= Q; ++j){
            for(int k=0; k <= R - i - j; ++k){
                if(dir == 1 && i > 0) { sol += i * pow_loc(x,i-1)* pow_loc(y,j)   * pow_loc(z,k); }
                if(dir == 2 && j > 0) { sol += j * pow_loc(x,i)  * pow_loc(y,j-1) * pow_loc(z,k); }
                if(dir == 3 && k > 0) { sol += k * pow_loc(x,i)  * pow_loc(y,j)   * pow_loc(z,k-1); }
            }
        }
    }
    return sol;
}


