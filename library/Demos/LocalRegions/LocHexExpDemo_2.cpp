#include <cstdio>
#include <cstdlib>

#include <LocalRegions/HexExp.h>
#include <SpatialDomains/MeshGraph3D.h>

using namespace Nektar;

/// modification to deal with exact solution. Return 1 if integer < 0
static double  pow_loc(const double val, const int i)
{
  return (i < 0)? 1.0: pow(val,i);
}

/// Solution to excite all modes in Hex expansion.
NekDouble Hex_sol(  NekDouble x, NekDouble y, NekDouble z,
                    int order1, int order2, int order3,
                    LibUtilities::BasisType btype1,
                    LibUtilities::BasisType btype2,
                    LibUtilities::BasisType btype3);

/// Derivative of solution in Hex expansion.
NekDouble Hex_Diff_Sol( NekDouble x, NekDouble y, NekDouble z,
                    int P, int Q, int R,
                    LibUtilities::BasisType btype1,
                    LibUtilities::BasisType btype2,
                    LibUtilities::BasisType btype3, int dir);


int main(int argc, char *argv[])
{
    if(argc != 8)
    {
        cerr << "usage: LocHexExpDemo MeshFile nummodes0 nummodes1 nummodes2 "
                "nx ny nz" << endl;
        exit(1);
    }

    int i;

    string in(argv[1]);
    int nummodes0   = atoi(argv[2]);
    int nummodes1   = atoi(argv[3]);
    int nummodes2   = atoi(argv[4]);
    int nquad0      = atoi(argv[5]);
    int nquad1      = atoi(argv[6]);
    int nquad2      = atoi(argv[7]);
    int P           = nummodes0 - 1;
    int Q           = nummodes1 - 1;
    int R           = nummodes2 - 1;

    int ntotquad = nquad0*nquad1*nquad2;

    SpatialDomains::MeshGraph3D graph3D;
    graph3D.ReadGeometry(in);

    LibUtilities::BasisType  bType = LibUtilities::eModified_A;
    LibUtilities::PointsType qtype = LibUtilities::eGaussLobattoLegendre;

    const LibUtilities::PointsKey   pointsKey0( nquad0, qtype );
    const LibUtilities::PointsKey   pointsKey1( nquad1, qtype );
    const LibUtilities::PointsKey   pointsKey2( nquad2, qtype );

    const LibUtilities::BasisKey    basisKey0( bType, nummodes0, pointsKey0 );
    const LibUtilities::BasisKey    basisKey1( bType, nummodes1, pointsKey1 );
    const LibUtilities::BasisKey    basisKey2( bType, nummodes2, pointsKey2 );

    SpatialDomains::HexGeomSharedPtr geom;
    if(!(geom = boost::dynamic_pointer_cast<SpatialDomains::HexGeom>(
                                                graph3D.GetCompositeItem(0,0))))
    {
        cerr << "Could not find HexGeom in input file" << endl;
        exit(1);
    }

    LocalRegions::HexExpSharedPtr E = MemoryManager<LocalRegions::HexExp>
                        ::AllocateSharedPtr(basisKey0,basisKey1,basisKey2,geom);

    Array<OneD, NekDouble> sol( ntotquad );
    Array<OneD, NekDouble> dx( ntotquad, 0.0 );
    Array<OneD, NekDouble> dy( ntotquad, 0.0 );
    Array<OneD, NekDouble> dz( ntotquad, 0.0 );
    Array<OneD, NekDouble> diff_solution_x( ntotquad, 0.0 );
    Array<OneD, NekDouble> diff_solution_y( ntotquad, 0.0 );
    Array<OneD, NekDouble> diff_solution_z( ntotquad, 0.0 );
    Array<OneD, NekDouble> derivatives( ntotquad, 0.0 );
    Array<OneD, NekDouble> x( ntotquad );
    Array<OneD, NekDouble> y( ntotquad );
    Array<OneD, NekDouble> z( ntotquad );

    E->GetCoords(x,y,z);

    //----------------------------------------------
    // Define solution to be projected
    for(i = 0; i < ntotquad; i++)
    {
        sol[i]  = Hex_sol(x[i],y[i],z[i],P,Q,R, bType, bType, bType);
    }
    //----------------------------------------------

    //----------------------------------------------
    // Define the derivative solution
     for(int n = 0; n < ntotquad; ++n){
        diff_solution_x[n] = Hex_Diff_Sol( x[n], y[n], z[n], P,Q,R,
                                           bType, bType, bType, 1);
        diff_solution_y[n] = Hex_Diff_Sol( x[n], y[n], z[n], P,Q,R,
                                           bType, bType, bType, 2);
        diff_solution_z[n] = Hex_Diff_Sol( x[n], y[n], z[n], P,Q,R,
                                           bType, bType, bType, 3);
    }
    //---------------------------------------------

    //---------------------------------------------
    // Project onto Expansion
    E->FwdTrans( sol, E->UpdateCoeffs() );
    //---------------------------------------------

    //-------------------------------------------
    // Backward Transform Solution to get projected values
    E->BwdTrans( E->GetCoeffs(), E->UpdatePhys() );
    //-------------------------------------------

    //--------------------------------------------
    // Calculate L_p error
    cout << "L infinity error: " << E->Linf(sol) << endl;
    cout << "L 2 error:        " << E->L2  (sol) << endl;
    //--------------------------------------------

    //-------------------------------------------
    // Evaulate solution at x = y = z = 0  and print error
    Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);

    t[0] = -0.39;
    t[1] = -0.25;
    t[2] = 0.5;
    NekDouble exact_sol = 0.0;

    exact_sol = Hex_sol( t[0], t[1], t[2], nummodes0, nummodes1, nummodes2,
                         bType, bType, bType);

    NekDouble numericSolution = E->PhysEvaluate(t);
    cout << "Numeric solution = " << numericSolution << endl;
    cout << "Exact solution = " << exact_sol << endl;
    cout << "Difference at x = ( "
         << t[0] << ", " << t[1] << ", " << t[2] << " ): "
         << numericSolution - exact_sol << endl;
    //-------------------------------------------

    E->PhysDeriv( sol, dx, dy, dz);

    double error_x = 0, error_y=0, error_z=0;

    for( int n = 0; n < ntotquad; ++n ) {
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
                    LibUtilities::BasisType btype3)
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
                    LibUtilities::BasisType btype3, int dir){
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
                switch (dir)
                {
                case 1:
                    a  = i*(Fx ? M_PI*(cos(M_PI*i*x) - sin(M_PI*i*x))
                               : pow_loc(x,i-1));
                    a *= (Fy ? (sin(M_PI*j*y) + cos(M_PI*j*y)) : pow_loc(y,j));
                    a *= (Fz ? (sin(M_PI*k*z) + cos(M_PI*k*z)) : pow_loc(z,k));
                    sol += a;
                    break;
                case 2:
                    a  = (Fx ? (sin(M_PI*i*x) + cos(M_PI*i*x)) : pow_loc(x,i));
                    a *= j*(Fy ? M_PI*(cos(M_PI*j*y) - sin(M_PI*j*y))
                               : pow_loc(y,j-1));
                    a *= (Fz ? (sin(M_PI*k*z) + cos(M_PI*k*z)) : pow_loc(z,k));
                    sol += a;
                    break;
                case 3:
                    a  = (Fx ? (sin(M_PI*i*x) + cos(M_PI*i*x)) : pow_loc(x,i));
                    a *= (Fy ? (sin(M_PI*j*y) + cos(M_PI*j*y)) : pow_loc(y,j));
                    a *= k*(Fz ? M_PI*(cos(M_PI*k*z) - sin(M_PI*k*z))
                               : pow_loc(z,k-1));
                    sol += a;
                    break;
                }
            }
        }
    }
    return sol;
}

