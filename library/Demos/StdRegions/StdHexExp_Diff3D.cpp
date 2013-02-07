#include <iostream>
#include <cstdlib>
#include <cmath>
using namespace std;

#include <StdRegions/StdHexExp.h>
#include <StdRegions/StdRegions.hpp>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

using namespace Nektar;

using namespace Nektar::LibUtilities;
using namespace Nektar::StdRegions;

/// Defines a solution which excites all modes in a Hex expansion.
NekDouble Hex_sol(NekDouble x, NekDouble y, NekDouble z,
                    int order1, int order2, int order3,
                    LibUtilities::BasisType bType_x,
                    LibUtilities::BasisType bType_y,
                    LibUtilities::BasisType bType_z);

/// Derivative of solution on Hex expansion.
NekDouble Hex_Diff_Sol(NekDouble x, NekDouble y, NekDouble z,
                    int P, int Q, int R,
                    LibUtilities::BasisType bType_x,
                    LibUtilities::BasisType bType_y,
                    LibUtilities::BasisType bType_z, int dir);

/// modification to deal with exact solution. Return 1 if integer < 0
static double  pow_loc(const double val, const int i)
{
  return (i < 0)? 1.0: pow(val,i);
}


/// This routine projects a polynomial or trigonmetric functions which
/// has energy in all mdoes of the expansions and reports and error
int main(int argc, char *argv[]) {
    if( argc != 10 ) {
        cerr << "Usage: HexDemo Type_x Type_y Type_z numModes_x numModes_y "
                "numModes_z Qx Qy Qz" << endl;
        cerr << "Where type is an interger value which dictates the basis as:"
             << endl;
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
        cerr << "\n\n" << "Example: " << argv[0] << " 1 2 3 3 3 3 5 5 5"
             << endl << endl;
        exit(1);
    }

    int bType_x_val = atoi(argv[1]);
    int bType_y_val = atoi(argv[2]);
    int bType_z_val = atoi(argv[3]);

    BasisType   bType_x = static_cast<BasisType>( bType_x_val );
    BasisType   bType_y = static_cast<BasisType>( bType_y_val );
    BasisType   bType_z = static_cast<BasisType>( bType_z_val );

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

    int xModes  = atoi(argv[4]);
    int yModes  = atoi(argv[5]);
    int zModes  = atoi(argv[6]);
    int Qx      = atoi(argv[7]);
    int Qy      = atoi(argv[8]);
    int Qz      = atoi(argv[9]);
    int P       = xModes - 1;
    int Q       = yModes - 1;
    int R       = zModes - 1;

    Array<OneD, NekDouble> solution( Qx * Qy * Qz );
    Array<OneD, NekDouble> dx(Qx*Qy*Qz,0.0);
    Array<OneD, NekDouble> dy(Qx*Qy*Qz,0.0);
    Array<OneD, NekDouble> dz(Qx*Qy*Qz,0.0);
    Array<OneD, NekDouble> diff_solution_x( Qx * Qy * Qz, 0.0 );
    Array<OneD, NekDouble> diff_solution_y( Qx * Qy * Qz, 0.0 );
    Array<OneD, NekDouble> diff_solution_z( Qx * Qy * Qz, 0.0 );
    Array<OneD, NekDouble> derivatives( Qx * Qy * Qz, 0.0 );

    LibUtilities::PointsType    Qtype_x, Qtype_y, Qtype_z;
    Array<OneD, NekDouble> x( Qx * Qy * Qz );
    Array<OneD, NekDouble> y( Qx * Qy * Qz );
    Array<OneD, NekDouble> z( Qx * Qy * Qz );

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
        solution[n]  = Hex_sol( x[n], y[n], z[n], P, Q, R,
                                bType_x, bType_y, bType_z );
    }
    //----------------------------------------------

    //----------------------------------------------
    // Define the derivative solution
    for(int n = 0; n < Qx*Qy*Qz; ++n){
        diff_solution_x[n] = Hex_Diff_Sol( x[n], y[n], z[n], P, Q, R,
                                           bType_x, bType_y, bType_z, 1);
        diff_solution_y[n] = Hex_Diff_Sol( x[n], y[n], z[n], P, Q, R,
                                           bType_x, bType_y, bType_z, 2);
        diff_solution_z[n] = Hex_Diff_Sol( x[n], y[n], z[n], P, Q, R,
                                           bType_x, bType_y, bType_z, 3);
    }
    //---------------------------------------------

    //---------------------------------------------
    // Project onto Expansion
    she->FwdTrans( solution, she->UpdateCoeffs() );
    //---------------------------------------------

    //-------------------------------------------
    // Backward Transform Solution to get projected values
    she->BwdTrans( she->GetCoeffs(), she->UpdatePhys() );
    //-------------------------------------------

    //--------------------------------------------
    // Calculate L_p error
    cout << "L infinity error: " << she->Linf(solution) << endl;
    cout << "L 2 error:        " << she->L2  (solution) << endl;
    //--------------------------------------------


    //-------------------------------------------
    // Evaulate solution at x = y = z = 0  and print error
    Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);

    t[0] = -0.39;
    t[1] = -0.25;
    t[2] = 0.5;
    solution[0] = Hex_sol( t[0], t[1], t[2], P, Q, R,
                           bType_x, bType_y, bType_z );

    NekDouble numericSolution = she->PhysEvaluate(t);
    cout << "Interpolation difference from actual solution at x = ( "
         << t[0] << ", " << t[1] << ", " << t[2] << " ): "
         << numericSolution - solution[0] << endl;
    cout << "numericSolution = " << numericSolution << endl;
    cout << "exact solution = " << solution[0] << endl;
    //-------------------------------------------

    //--------------------------------------------
    // Taking the physical derivative and putting them into dx, dy, dz.
    she->PhysDeriv( she->GetPhys(), dx, dy, dz );
    //--------------------------------------------

    double error_x = 0, error_y=0, error_z=0;

    for( int n = 0; n < Qx*Qy*Qz; ++n )
    {
        error_x += fabs(diff_solution_x[n] - dx[n]);
        error_y += fabs(diff_solution_y[n] - dy[n]);
        error_z += fabs(diff_solution_z[n] - dz[n]);
    }

    cout << "L 1 error of derivatives X =  " <<  error_x << endl;
    cout << "L 1 error of derivatives Y =  " <<  error_y << endl;
    cout << "L 1 error of derivatives Z =  " <<  error_z << endl;

    return 0;
}


NekDouble Hex_sol(NekDouble x, NekDouble y, NekDouble z,
                    int order1, int order2, int order3,
                    LibUtilities::BasisType btype1,
                    LibUtilities::BasisType btype2,
                    LibUtilities::BasisType btype3 )
{
    int i,j,k;
    NekDouble sol = 0.0;

    int  Nx = (btype1 == LibUtilities::eFourier ? order1/2 : order1);
    int  Ny = (btype2 == LibUtilities::eFourier ? order2/2 : order2);
    int  Nz = (btype3 == LibUtilities::eFourier ? order3/2 : order3);
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


NekDouble Hex_Diff_Sol(NekDouble x, NekDouble y, NekDouble z,
                    int order1, int order2, int order3,
                    LibUtilities::BasisType btype1,
                    LibUtilities::BasisType btype2,
                    LibUtilities::BasisType btype3, int dir)
{
    int i,j,k;
    NekDouble sol = 0.0;

    int  Nx = (btype1 == LibUtilities::eFourier ? order1/2 : order1);
    int  Ny = (btype2 == LibUtilities::eFourier ? order2/2 : order2);
    int  Nz = (btype3 == LibUtilities::eFourier ? order3/2 : order3);
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
                a  = i*(Fx ? M_PI*(cos(M_PI*i*x) - sin(M_PI*i*x))
                           : pow_loc(x,i-1));
                a *= (Fy ? (sin(M_PI*j*y) + cos(M_PI*j*y)) : pow_loc(y,j));
                a *= (Fz ? (sin(M_PI*k*z) + cos(M_PI*k*z)) : pow_loc(z,k));
                sol += a;
                a  = (Fx ? (sin(M_PI*i*x) + cos(M_PI*i*x)) : pow_loc(x,i));
                a *= j*(Fy ? M_PI*(cos(M_PI*j*y) - sin(M_PI*j*y))
                           : pow_loc(y,j-1));
                a *= (Fz ? (sin(M_PI*k*z) + cos(M_PI*k*z)) : pow_loc(z,k));
                sol += a;
                a  = (Fx ? (sin(M_PI*i*x) + cos(M_PI*i*x)) : pow_loc(x,i));
                a *= (Fy ? (sin(M_PI*j*y) + cos(M_PI*j*y)) : pow_loc(y,j));
                a *= k*(Fz ? M_PI*(cos(M_PI*k*z) - sin(M_PI*k*z))
                           : pow_loc(z,k-1));
                sol += a;
            }
        }
    }
    return sol;
}

