#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <StdRegions/StdHexExp.h>
#include <StdRegions/StdPrismExp.h>
#include <StdRegions/StdTetExp.h>

#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

using namespace Nektar;
using namespace Nektar::LibUtilities;
using namespace Nektar::StdRegions;

/// Defines a solution which excites all modes in a Tet expansion.
NekDouble Tet_sol(NekDouble x, NekDouble y, NekDouble z,
                  int order1, int order2, int order3);

/// Derivative of solution on Tet expansion.
NekDouble Tet_Dsol(NekDouble x, NekDouble y, NekDouble z,
                   int order1, int order2, int order3);

/// Defines a solution which excites all modes in a prismatic expansion.
NekDouble Prism_sol(NekDouble x, NekDouble y, NekDouble z,
                    int order1, int order2, int order3);

/// Derivative of solution on a prismatic expansion.
NekDouble Prism_Dsol(NekDouble x, NekDouble y, NekDouble z,
                     int order1, int order2, int order3);

/// Defines a solution which excites all modes in a Hex expansion.
NekDouble Hex_sol(NekDouble x, NekDouble y, NekDouble z,
                  int order1, int order2, int order3,
                  LibUtilities::BasisType btype1,
                  LibUtilities::BasisType btype2,
                  LibUtilities::BasisType btype3);

/// Derivative of solution on a Hex expansion.
NekDouble Hex_Dsol(NekDouble x, NekDouble y, NekDouble z,
                   int order1, int order2, int order3,
                   LibUtilities::BasisType btype1,
                   LibUtilities::BasisType btype2,
                   LibUtilities::BasisType btype3);

/// modification to deal with exact solution. Return 1 if integer < 0
static double  pow_loc(const double val, const int i)
{
    return (i < 0)? 1.0: pow(val,i);
}


/// This routine projects a polynomial or trigonmetric functions which
/// has energy in all mdoes of the expansions and reports the error
int main(int argc, char *argv[]){
    int           i;

    int           order1,order2,order3, nq1,nq2,nq3;
    LibUtilities::PointsType    Qtype1,Qtype2,Qtype3;
    LibUtilities::BasisType     btype1 =   LibUtilities::eOrtho_A;
    LibUtilities::BasisType     btype2 =   LibUtilities::eOrtho_B;
    LibUtilities::BasisType     btype3 =   LibUtilities::eOrtho_C;
    LibUtilities::ShapeType     regionshape;
    StdRegions::StdExpansion *E;
    Array<OneD, NekDouble> x, y, z, sol, dx, dy, dz;

    if(argc != 11)
    {
        fprintf(stderr,"Usage: StdProject2D RegionShape Type1 Type2 Type3 "
                       "order1 order2 order3 nq1 nq2 nq3 \n");

        fprintf(stderr,"Where RegionShape is an integer value which "
                       "dictates the region shape:\n");
        fprintf(stderr,"\t Tetrahedron   = 5\n");
        fprintf(stderr,"\t Prism         = 7\n");
        fprintf(stderr,"\t Hexahedron    = 8\n");

        fprintf(stderr,"Where type is an integer value which "
                       "dictates the basis as:\n");

        fprintf(stderr,"\t Ortho_A    = 1\n");
        fprintf(stderr,"\t Ortho_B    = 2\n");
        fprintf(stderr,"\t Ortho_C    = 3\n");
        fprintf(stderr,"\t Modified_A = 4\n");
        fprintf(stderr,"\t Modified_B = 5\n");
        fprintf(stderr,"\t Modified_C = 6\n");
        fprintf(stderr,"\t Fourier    = 7\n");
        fprintf(stderr,"\t Lagrange   = 8\n");
        fprintf(stderr,"\t Legendre   = 10\n");
        fprintf(stderr,"\t Chebyshev  = 11\n");
        fprintf(stderr,"\t Nodal tri (Electro) = 12\n");
        fprintf(stderr,"\t Nodal tri (Fekete)  = 13\n");

        exit(1);
    }

    regionshape = (LibUtilities::ShapeType) atoi(argv[1]);

    // Check to see if 3D region
    if (regionshape != LibUtilities::eTetrahedron &&
        regionshape != LibUtilities::ePrism       &&
        regionshape != LibUtilities::eHexahedron)
    {
        NEKERROR(ErrorUtil::efatal,"This shape is not a 3D region");
    }

    int btype1_val = atoi(argv[2]);
    int btype2_val = atoi(argv[3]);
    int btype3_val = atoi(argv[4]);
    if(( btype1_val <= 11)&&( btype2_val <= 11))
    {
        btype1 =   (LibUtilities::BasisType) btype1_val;
        btype2 =   (LibUtilities::BasisType) btype2_val;
        btype3 =   (LibUtilities::BasisType) btype3_val;
    }
    else if(( btype1_val >=12)&&(btype2_val <= 13))
    {
        btype1 =   LibUtilities::eOrtho_A;
        btype2 =   LibUtilities::eOrtho_B;
        btype3 =   LibUtilities::eOrtho_C;
    }

    // Check to see that correct Expansions are used
    switch(regionshape)
    {
    case LibUtilities::eTetrahedron:
        if((btype1 == eOrtho_B) || (btype1 == eOrtho_C)
           || (btype1 == eModified_B) || (btype1 == eModified_C))
        {
            NEKERROR(ErrorUtil::efatal, "Basis 1 cannot be of type Ortho_B, "
                     "Ortho_C, Modified_B or Modified_C");
        }
        if((btype2 == eOrtho_A) || (btype2 == eOrtho_C)
           || (btype2 == eModified_A) || (btype2 == eModified_C))
        {
            NEKERROR(ErrorUtil::efatal, "Basis 2 cannot be of type Ortho_A, "
                     "Ortho_C, Modified_A or Modified_C");
        }
        if((btype3 == eOrtho_A) || (btype3 == eOrtho_B)
           || (btype3 == eModified_A) || (btype3 == eModified_B))
        {
            NEKERROR(ErrorUtil::efatal, "Basis 3 cannot be of type Ortho_A, "
                     "Ortho_B, Modified_A or Modified_B");
        }
        break;
    case LibUtilities::ePrism:
        if((btype1 == eOrtho_B) || (btype1 == eOrtho_C)
           || (btype1 == eModified_B) || (btype1 == eModified_C))
        {
            NEKERROR(ErrorUtil::efatal, "Basis 1 cannot be of type Ortho_B, "
                     "Ortho_C, Modified_B or Modified_C");
        }
        if((btype2 == eOrtho_B) || (btype2 == eOrtho_C)
           || (btype2 == eModified_B) || (btype2 == eModified_C))
        {
            NEKERROR(ErrorUtil::efatal, "Basis 2 cannot be of type Ortho_B, "
                     "Ortho_C, Modified_B or Modified_C");
        }
        if((btype3 == eOrtho_A) || (btype3 == eOrtho_C)
           || (btype3 == eModified_A) || (btype3 == eModified_C))
        {
            NEKERROR(ErrorUtil::efatal, "Basis 3 cannot be of type Ortho_A, "
                     "Ortho_C, Modified_A or Modified_C");
        }
        break;
    case LibUtilities::eHexahedron:
        if((btype1 == eOrtho_B) || (btype1 == eOrtho_C)
           || (btype1 == eModified_B) || (btype1 == eModified_C))
        {
            NEKERROR(ErrorUtil::efatal, "Basis 1 is for 2 or 3D expansions");
        }
        if((btype2 == eOrtho_B) || (btype2 == eOrtho_C)
           || (btype2 == eModified_B) || (btype2 == eModified_C))
        {
            NEKERROR(ErrorUtil::efatal, "Basis 2 is for 2 or 3D expansions");
        }
        if((btype3 == eOrtho_B) || (btype3 == eOrtho_C)
           || (btype3 == eModified_B) || (btype3 == eModified_C))
        {
            NEKERROR(ErrorUtil::efatal, "Basis 3 is for 2 or 3D expansions");
        }
        break;
    default:
        ASSERTL0(false, "Not a 3D expansion.");
        break;
    }
    
    order1 =   atoi(argv[5]);
    order2 =   atoi(argv[6]);
    order3 =   atoi(argv[7]);
    nq1    =   atoi(argv[8]);
    nq2    =   atoi(argv[9]);
    nq3    =   atoi(argv[10]);

    sol = Array<OneD, NekDouble>(nq1*nq2*nq3);
    dx  = Array<OneD, NekDouble>(nq1*nq2*nq3);
    dy  = Array<OneD, NekDouble>(nq1*nq2*nq3);
    dz  = Array<OneD, NekDouble>(nq1*nq2*nq3);
    x   = Array<OneD, NekDouble>(nq1*nq2*nq3);
    y   = Array<OneD, NekDouble>(nq1*nq2*nq3);
    z   = Array<OneD, NekDouble>(nq1*nq2*nq3);

    if(btype1 != LibUtilities::eFourier)
    {
        Qtype1 = LibUtilities::eGaussLobattoLegendre;
    }
    else
    {
        Qtype1 = LibUtilities::eFourierEvenlySpaced;
    }

    if(btype2 != LibUtilities::eFourier)
    {
        if (regionshape == LibUtilities::eTetrahedron) {
            Qtype2 = LibUtilities::eGaussRadauMAlpha1Beta0;
        }
        else
        {
            Qtype2 = LibUtilities::eGaussLobattoLegendre;
        }
    }
    else
    {
        Qtype2 = LibUtilities::eFourierEvenlySpaced;
    }

    if(btype3 != LibUtilities::eFourier)
    {
        if (regionshape == LibUtilities::eTetrahedron) 
        {
            Qtype3 = LibUtilities::eGaussRadauMAlpha2Beta0;
        }
        else if (regionshape == LibUtilities::ePrism) 
        {
            Qtype3 = LibUtilities::eGaussRadauMAlpha1Beta0;
        }
        else
        {
            Qtype3 = LibUtilities::eGaussLobattoLegendre;
        }
    }
    else
    {
        Qtype3 = LibUtilities::eFourierEvenlySpaced;
    }

    //-----------------------------------------------
    // Define a 3D expansion based on basis definition

    switch(regionshape)
    {
    case LibUtilities::eTetrahedron:
        {
            const LibUtilities::PointsKey Pkey1(nq1,Qtype1);
            const LibUtilities::PointsKey Pkey2(nq2,Qtype2);
            const LibUtilities::PointsKey Pkey3(nq3,Qtype3);
            const LibUtilities::BasisKey  Bkey1(btype1,order1,Pkey1);
            const LibUtilities::BasisKey  Bkey2(btype2,order2,Pkey2);
            const LibUtilities::BasisKey  Bkey3(btype3,order3,Pkey3);
            
            E = new StdRegions::StdTetExp(Bkey1,Bkey2,Bkey3);
            E->GetCoords(x,y,z);

            //----------------------------------------------
            // Define solution to be projected
            for(i = 0; i < nq1*nq2*nq3; ++i)
            {
                sol[i]  = Tet_sol(x[i],y[i],z[i],order1,order2,order3);
            }
            //----------------------------------------------
        }
        break;
    case LibUtilities::ePrism:
        {
            const LibUtilities::PointsKey Pkey1(nq1,Qtype1);
            const LibUtilities::PointsKey Pkey2(nq2,Qtype2);
            const LibUtilities::PointsKey Pkey3(nq3,Qtype3);
            const LibUtilities::BasisKey  Bkey1(btype1,order1,Pkey1);
            const LibUtilities::BasisKey  Bkey2(btype2,order2,Pkey2);
            const LibUtilities::BasisKey  Bkey3(btype3,order3,Pkey3);
            
            E = new StdRegions::StdPrismExp(Bkey1,Bkey2,Bkey3);
            E->GetCoords(x,y,z);

            //----------------------------------------------
            // Define solution to be projected
            for(i = 0; i < nq1*nq2*nq3; ++i)
            {
                sol[i]  = Prism_sol(x[i],y[i],z[i],order1,order2,order3);
            }
            //----------------------------------------------
        }
        break;
    case LibUtilities::eHexahedron:
        {
            const LibUtilities::PointsKey Pkey1(nq1,Qtype1);
            const LibUtilities::PointsKey Pkey2(nq2,Qtype2);
            const LibUtilities::PointsKey Pkey3(nq3,Qtype3);
            const LibUtilities::BasisKey Bkey1(btype1,order1,Pkey1);
            const LibUtilities::BasisKey Bkey2(btype2,order2,Pkey2);
            const LibUtilities::BasisKey Bkey3(btype3,order3,Pkey3);

            E = new StdRegions::StdHexExp(Bkey1,Bkey2,Bkey3);
            E->GetCoords(x,y,z);

            //----------------------------------------------
            // Define solution to be projected
            for(i = 0; i < nq1*nq2*nq3; ++i)
            {
                sol[i]  = Hex_sol(x[i], y[i], z[i], order1, order2, order3,
                                  btype1, btype2, btype3);
            }
            //---------------------------------------------
        }
        break;
    default:
            ASSERTL0(false, "Not a 3D expansion.");
            break;
    }

    //---------------------------------------------
    // Evaluate derivative of solution, add together and put in sol
    E->PhysDeriv(sol,dx,dy,dz);
    Vmath::Vadd(nq1*nq2*nq3,dx,1,dy,1,sol,1);
    Vmath::Vadd(nq1*nq2*nq3,dz,1,sol,1,sol,1);
    //---------------------------------------------

    //---------------------------------------------
    // Project onto Expansion
    E->FwdTrans(sol,E->UpdateCoeffs());
    //---------------------------------------------

    //-------------------------------------------
    // Backward Transform Solution to get projected values
    E->BwdTrans(E->GetCoeffs(),E->UpdatePhys());
    //-------------------------------------------

    //----------------------------------------------
    // Define exact solution of differential
    switch(regionshape)
    {
    case LibUtilities::eTetrahedron:
        {
            for(i = 0; i < nq1*nq2*nq3; ++i)
            {
                sol[i] = Tet_Dsol(x[i],y[i],z[i],order1,order2,order3);
            }
        }
        break;
    case LibUtilities::ePrism:
        {
            for(i = 0; i < nq1*nq2*nq3; ++i)
            {
                sol[i] = Prism_Dsol(x[i],y[i],z[i],order1,order2,order3);
            }
        }
        break;
    case LibUtilities::eHexahedron:
        {
            for(i = 0; i < nq1*nq2*nq3; ++i)
            {
                sol[i] = Hex_Dsol(x[i], y[i], z[i], order1, order2, order3,
                                  btype1, btype2, btype3);
            }
        }
        break;
    default:
        ASSERTL0(false, "Not a 3D expansion.");
        break;
    }
    
    //--------------------------------------------
    // Calculate L_inf error
    cout << "L infinity error: " << E->Linf(sol) << endl;
    cout << "L 2 error:        " << E->L2  (sol) << endl;
    //--------------------------------------------

    return 0;
}


NekDouble Tet_sol(NekDouble x, NekDouble y, NekDouble z,
                  int order1, int order2, int order3)
{
    int    l,k,m;
    NekDouble sol = 0;

    for(k = 0; k < order1; ++k)
    {
        for(l = 0; l < order2-k; ++l)
        {
            for(m = 0; m < order3-k-l; ++m)
            {
                sol += pow_loc(x,k)*pow_loc(y,l)*pow_loc(z,m);
            }
        }
    }

    return sol;
}

NekDouble Prism_sol(NekDouble x, NekDouble y, NekDouble z,
                    int order1, int order2, int order3)
{
    int    l,k,m;
    NekDouble sol = 0;

    for(k = 0; k < order1; ++k)
    {
        for(l = 0; l < order2; ++l)
        {
            for(m = 0; m < order3-k; ++m)
            {
                sol += pow_loc(x,k)*pow_loc(y,l)*pow_loc(z,m);
            }
        }
    }

    return sol;
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


NekDouble Tet_Dsol(NekDouble x, NekDouble y, NekDouble z,
                   int order1, int order2, int order3)
{
    int    l,k,m;
    NekDouble sol = 0;

    for(k = 0; k < order1; ++k)
    {
        for(l = 0; l < order2-k; ++l)
        {
            for(m = 0; m < order3-k-l; ++m)
            {
                sol += k*pow_loc(x,k-1)*pow_loc(y,l)*pow_loc(z,m)
                    + pow_loc(x,k)*l*pow_loc(y,l-1)*pow_loc(z,m)
                    + pow_loc(x,k)*pow_loc(y,l)*m*pow_loc(z,m-1);
            }
        }
    }
    return sol;
}

NekDouble Prism_Dsol(NekDouble x, NekDouble y, NekDouble z,
                     int order1, int order2, int order3)
{
    int    l,k,m;
    NekDouble sol = 0;

    for(k = 0; k < order1; ++k)
    {
        for(l = 0; l < order2; ++l)
        {
            for(m = 0; m < order3-k; ++m)
            {
                sol += k*pow_loc(x,k-1)*pow_loc(y,l)*pow_loc(z,m)
                    + pow_loc(x,k)*l*pow_loc(y,l-1)*pow_loc(z,m)
                    + pow_loc(x,k)*pow_loc(y,l)*m*pow_loc(z,m-1);
            }
        }
    }
    return sol;
}

NekDouble Hex_Dsol(NekDouble x, NekDouble y, NekDouble z,
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

