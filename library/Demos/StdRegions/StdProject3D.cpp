#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <StdRegions/StdHexExp.h>
#include <StdRegions/StdPrismExp.h>
#include <StdRegions/StdNodalPrismExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdNodalTetExp.h>

#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

using namespace Nektar;
using namespace Nektar::LibUtilities;
using namespace Nektar::StdRegions;

/// Defines a solution which excites all modes in a StdTet expansion.
NekDouble Tet_sol(NekDouble x, NekDouble y, NekDouble z,
                  int order1, int order2, int order3);

/// Defines a solution which excites all modes in a StdPrism expansion.
NekDouble Prism_sol(NekDouble x, NekDouble y, NekDouble z,
                    int order1, int order2, int order3);

/// Defines a solution which excites all modes in a StdHex expansion.
NekDouble Hex_sol(NekDouble x, NekDouble y, NekDouble z,
                  int order1, int order2, int order3,
                  LibUtilities::BasisType btype1, LibUtilities::BasisType btype2,
                  LibUtilities::BasisType btype3);

/// This routine projects a polynomial or trigonmetric functions which
/// has energy in all mdoes of the expansions and reports and error
int main(int argc, char *argv[]){
    int           i;

    int           order1,order2,order3, nq1,nq2,nq3;
    LibUtilities::PointsType    Qtype1,Qtype2,Qtype3;
    LibUtilities::BasisType     btype1,btype2,btype3;
    LibUtilities::PointsType     NodalType;
    StdRegions::ExpansionType    regionshape;
    StdRegions::StdExpansion *E;
    Array<OneD, NekDouble> sol;

    if(argc != 11)
    {
        fprintf(stderr,"Usage: StdProject2D RegionShape Type1 Type2 Type3 "
                       "order1 order2 order3 nq1 nq2 nq3 \n");
        fprintf(stderr,"Where RegionShape is an integer value which "
                       "dictates the region shape:\n");
        fprintf(stderr,"\t Tetrahedron   = 4\n");
        fprintf(stderr,"\t Prism         = 6\n");
        fprintf(stderr,"\t Hexahedron    = 7\n");

        fprintf(stderr,"Where type is an integer value which "
                       "dictates the basis as:\n");

        fprintf(stderr,"\t Ortho_A             =  1\n");
        fprintf(stderr,"\t Ortho_B             =  2\n");
        fprintf(stderr,"\t Ortho_C             =  3\n");
        fprintf(stderr,"\t Modified_A          =  4\n");
        fprintf(stderr,"\t Modified_B          =  5\n");
        fprintf(stderr,"\t Modified_C          =  6\n");
        fprintf(stderr,"\t Fourier             =  7\n");
        fprintf(stderr,"\t Lagrange            =  8\n");
        fprintf(stderr,"\t Legendre            = 10\n");
        fprintf(stderr,"\t Chebyshev           = 11\n");
        fprintf(stderr,"\t Nodal tri (Electro) = 12\n");
        fprintf(stderr,"\t Nodal tri (Fekete)  = 13\n");
        fprintf(stderr,"\t Nodal tet (Electro) = 14\n");
        fprintf(stderr,"\t Nodal tet (Even)    = 15\n");
        fprintf(stderr,"\t Nodal prism (Even)  = 16\n");

        exit(1);
    }

    regionshape = (StdRegions::ExpansionType) atoi(argv[1]);

    // Check to see if 3D region
    if ((regionshape != StdRegions::eTetrahedron)
        && (regionshape != StdRegions::ePrism)
        && (regionshape != StdRegions::eHexahedron))
    {
        NEKERROR(ErrorUtil::efatal,"This shape is not a 3D region");
    }

    int btype1_val = atoi(argv[2]);
    int btype2_val = atoi(argv[3]);
    int btype3_val = atoi(argv[4]);
    
    if (btype1_val <= 11 && btype2_val <= 11)
    {
        btype1 =   (LibUtilities::BasisType) btype1_val;
        btype2 =   (LibUtilities::BasisType) btype2_val;
        btype3 =   (LibUtilities::BasisType) btype3_val;
    }
    else if(btype1_val >=12 && btype2_val <= 16)
    {
        if (regionshape == StdRegions::eTetrahedron)
        {
            btype1 = LibUtilities::eOrtho_A;
            btype2 = LibUtilities::eOrtho_B;
            btype3 = LibUtilities::eOrtho_C;
        }
        else if (regionshape == StdRegions::ePrism)
        {
            btype1 = LibUtilities::eOrtho_A;
            btype2 = LibUtilities::eOrtho_A;
            btype3 = LibUtilities::eOrtho_B;
        }
        
        if(btype1_val == 12)
        {
            NodalType = LibUtilities::eNodalTriElec;
        }
        else if (btype1_val == 13)
        {
            NodalType = LibUtilities::eNodalTriFekete;
        }
        else if (btype1_val == 14)
        {
            NodalType = LibUtilities::eNodalTetElec;
        }
        else if (btype1_val == 15)
        {
            NodalType = LibUtilities::eNodalTetEvenlySpaced;
        }
        else
        {
            NodalType = LibUtilities::eNodalPrismEvenlySpaced;
        }
    }

    // Check to see that correct Expansions are used
    switch(regionshape)
    {
        case StdRegions::eTetrahedron:
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
        case StdRegions::ePrism:
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
        case StdRegions::eHexahedron:
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
        if (regionshape == StdRegions::eTetrahedron) 
        {
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
        if (regionshape == StdRegions::eTetrahedron) 
        {
            Qtype3 = LibUtilities::eGaussRadauMAlpha2Beta0;
        }
        else if (regionshape == StdRegions::ePrism)
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
    Array<OneD,NekDouble> x(nq1*nq2*nq3);
    Array<OneD,NekDouble> y(nq1*nq2*nq3);
    Array<OneD,NekDouble> z(nq1*nq2*nq3);

    switch(regionshape)
    {
        case StdRegions::eTetrahedron:
        {
            const LibUtilities::PointsKey Pkey1(nq1,Qtype1);
            const LibUtilities::PointsKey Pkey2(nq2,Qtype2);
            const LibUtilities::PointsKey Pkey3(nq3,Qtype3);
            const LibUtilities::BasisKey  Bkey1(btype1,order1,Pkey1);
            const LibUtilities::BasisKey  Bkey2(btype2,order2,Pkey2);
            const LibUtilities::BasisKey  Bkey3(btype3,order3,Pkey3);

            if(btype1_val >= 10)
            {
                E = new StdRegions::StdNodalTetExp(Bkey1,Bkey2,Bkey3,NodalType);
            }
            else
            {
                E = new StdRegions::StdTetExp(Bkey1,Bkey2,Bkey3);
            }
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
        case StdRegions::ePrism:
        {
            const LibUtilities::PointsKey Pkey1(nq1,Qtype1);
            const LibUtilities::PointsKey Pkey2(nq2,Qtype2);
            const LibUtilities::PointsKey Pkey3(nq3,Qtype3);
            const LibUtilities::BasisKey  Bkey1(btype1,order1,Pkey1);
            const LibUtilities::BasisKey  Bkey2(btype2,order2,Pkey2);
            const LibUtilities::BasisKey  Bkey3(btype3,order3,Pkey3);

            if (btype1_val >= 10)
            {
                E = new StdRegions::StdNodalPrismExp(Bkey1,Bkey2,Bkey3,NodalType);
            }
            else
            {
                E = new StdRegions::StdPrismExp(Bkey1,Bkey2,Bkey3);
            }
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
        case StdRegions::eHexahedron:
        {
            const LibUtilities::PointsKey Pkey1(nq1,Qtype1);
            const LibUtilities::PointsKey Pkey2(nq2,Qtype2);
            const LibUtilities::PointsKey Pkey3(nq3,Qtype3);
            const LibUtilities::BasisKey Bkey1(btype1,order1,Pkey1);
            const LibUtilities::BasisKey Bkey2(btype2,order2,Pkey2);
            const LibUtilities::BasisKey Bkey3(btype3,order3,Pkey3);
            E = new StdRegions::StdHexExp(Bkey1,Bkey2,Bkey3);

            //----------------------------------------------
            // Define solution to be projected
            Array<OneD,NekDouble> x = Array<OneD,NekDouble>(nq1*nq2*nq3);
            Array<OneD,NekDouble> y = Array<OneD,NekDouble>(nq1*nq2*nq3);
            Array<OneD,NekDouble> z = Array<OneD,NekDouble>(nq1*nq2*nq3);
            E->GetCoords(x,y,z);

            for(i = 0; i < nq1*nq2*nq3; ++i)
            {
                sol[i]  = Hex_sol(x[i],y[i],z[i],order1,order2,order3,
                                  btype1,btype2,btype3);
            }
            //---------------------------------------------
        }
        break;
        default:
            ASSERTL0(false, "Not a 3D expansion.");
            break;
    }

    //---------------------------------------------
    // Project onto Expansion
    E->FwdTrans(sol,E->UpdateCoeffs());
    //---------------------------------------------

    //-------------------------------------------
    // Backward Transform Solution to get projected values
    E->BwdTrans(E->GetCoeffs(),E->UpdatePhys());
    //-------------------------------------------

    //--------------------------------------------
    // Calculate L_inf error
    cout << "L infinity error: " << E->Linf(sol) << endl;
    cout << "L 2 error:        " << E->L2  (sol) << endl;
    //--------------------------------------------

    //-------------------------------------------
    // Evaulate solution at x = y = z = 0  and print error
    Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);
    t[0] = -0.5;
    t[1] = -0.25;
    t[2] = -0.3;

    if(regionshape == StdRegions::eTetrahedron)
    {
        sol[0] = Tet_sol(t[0], t[1], t[2], order1, order2, order3);
    }
    else if (regionshape == StdRegions::ePrism)
    {
        sol[0] = Prism_sol(t[0], t[1], t[2], order1, order2, order3);
    }
    else
    {
        sol[0] = Hex_sol(t[0], t[1], t[2], order1, order2, order3,
                         btype1, btype2, btype3);
    }

    NekDouble nsol = E->PhysEvaluate(x);
    cout << "error at x = (" << t[0] << "," << t[1] << "," << t[2] << "): "
         << nsol - sol[0] << endl;
    //-------------------------------------------

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
                sol += pow(x,k)*pow(y,l)*pow(z,m);
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
                sol += pow(x,k)*pow(y,l)*pow(z,m);
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
                a  = (Fx ? sin(M_PI*i*x) + cos(M_PI*i*x) : pow(x,i));
                a *= (Fy ? sin(M_PI*j*y) + cos(M_PI*j*y) : pow(y,j));
                a *= (Fz ? sin(M_PI*k*z) + cos(M_PI*k*z) : pow(z,k));
                sol += a;
            }
        }
    }
    return sol;
}

