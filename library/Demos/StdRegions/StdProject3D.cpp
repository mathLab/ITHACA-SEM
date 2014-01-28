#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <boost/format.hpp>
#include <StdRegions/StdHexExp.h>
#include <StdRegions/StdPrismExp.h>
#include <StdRegions/StdNodalPrismExp.h>
#include <StdRegions/StdPyrExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdNodalTetExp.h>

#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

using namespace Nektar;
using namespace Nektar::LibUtilities;
using namespace Nektar::StdRegions;

void printSolution(StdRegions::StdExpansion *F,
                   std::string name,
                   Array<OneD, NekDouble> &x,
                   Array<OneD, NekDouble> &y,
                   Array<OneD, NekDouble> &z,
                   Array<OneD, NekDouble> &phys)
{
    int nPts = F->GetTotPoints();
    ofstream outf(name);
    int numBlocks = 1;
    for (int j = 0; j < 3; ++j)
    {
        numBlocks *= F->GetNumPoints(j)-1;
    }

    outf << "VARIABLES = x, y, z, mnpp" << endl;
    outf << "Zone, N=" << nPts << ", E="
         << numBlocks << ", F=FEBlock, ET=BRICK" << endl;
    
    const int np0 = F->GetNumPoints(0);
    const int np1 = F->GetNumPoints(1);
    const int np2 = F->GetNumPoints(2);
    const int np01 = np0*np1;

    for (int j = 0; j < nPts; ++j)
    {
        outf << x[j] << endl;
    }

    for (int j = 0; j < nPts; ++j)
    {
        outf << y[j] << endl;
    }

    for (int j = 0; j < nPts; ++j)
    {
        outf << z[j] << endl;
    }

    for (int j = 0; j < nPts; ++j)
    {
        outf << phys[j] << endl;
    }

    for(int j = 1; j < np2; ++j)
    {
        for(int k = 1; k < np1; ++k)
        {
            for(int l = 1; l < np0; ++l)
            {
                outf << (j-1)*np01 + (k-1)*np0 + l   << " ";
                outf << (j-1)*np01 + (k-1)*np0 + l+1 << " ";
                outf << (j-1)*np01 +  k   *np0 + l+1 << " ";
                outf << (j-1)*np01 +  k   *np0 + l   << " ";
                outf <<  j   *np01 + (k-1)*np0 + l   << " ";
                outf <<  j   *np01 + (k-1)*np0 + l+1 << " ";
                outf <<  j   *np01 +  k   *np0 + l+1 << " ";
                outf <<  j   *np01 +  k   *np0 + l   << endl;
            }
        }
    }
    
    outf.close();
}

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
    LibUtilities::BasisType     btype1 =   LibUtilities::eOrtho_A;
    LibUtilities::BasisType     btype2 =   LibUtilities::eOrtho_B;
    LibUtilities::BasisType     btype3 =   LibUtilities::eOrtho_C;
    LibUtilities::PointsType    NodalType = LibUtilities::eNodalTriElec;
    LibUtilities::ShapeType     regionshape;
    StdRegions::StdExpansion *E = NULL;
    Array<OneD, NekDouble> sol;

    if(argc != 11)
    {
        fprintf(stderr,"Usage: StdProject2D RegionShape Type1 Type2 Type3 "
                       "order1 order2 order3 nq1 nq2 nq3 \n");
        fprintf(stderr,"Where RegionShape is an integer value which "
                       "dictates the region shape:\n");
        fprintf(stderr,"\t Tetrahedron   = 5\n");
        fprintf(stderr,"\t Pyramid       = 6\n");
        fprintf(stderr,"\t Prism         = 7\n");
        fprintf(stderr,"\t Hexahedron    = 8\n");

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
        fprintf(stderr,"\t Gauss Lagrange      =  9\n");
        fprintf(stderr,"\t Legendre            = 10\n");
        fprintf(stderr,"\t Chebyshev           = 11\n");
        fprintf(stderr,"\t Nodal tri (Electro) = 12\n");
        fprintf(stderr,"\t Nodal tri (Fekete)  = 13\n");
        fprintf(stderr,"\t Nodal tet (Electro) = 14\n");
        fprintf(stderr,"\t Nodal tet (Even)    = 15\n");
        fprintf(stderr,"\t Nodal prism (Even)  = 16\n");

        exit(1);
    }

    regionshape = (LibUtilities::ShapeType) atoi(argv[1]);

    // Check to see if 3D region
    if ((regionshape != LibUtilities::eTetrahedron)
        && (regionshape != LibUtilities::ePyramid)
        && (regionshape != LibUtilities::ePrism)
        && (regionshape != LibUtilities::eHexahedron))
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
        if (regionshape == LibUtilities::eTetrahedron)
        {
            btype1 = LibUtilities::eOrtho_A;
            btype2 = LibUtilities::eOrtho_B;
            btype3 = LibUtilities::eOrtho_C;
        }
        else if (regionshape == LibUtilities::ePrism)
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
        case LibUtilities::ePyramid:
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

    vector<StdPyrExp::triple> pyrIdx;
    
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
        if (regionshape == LibUtilities::eTetrahedron) 
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
        if (regionshape == LibUtilities::eTetrahedron ||
            regionshape == LibUtilities::ePyramid) 
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
    Array<OneD,NekDouble> x(nq1*nq2*nq3);
    Array<OneD,NekDouble> y(nq1*nq2*nq3);
    Array<OneD,NekDouble> z(nq1*nq2*nq3);

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
    case LibUtilities::ePyramid:
        {
            const LibUtilities::PointsKey Pkey1(nq1,Qtype1);
            const LibUtilities::PointsKey Pkey2(nq2,Qtype2);
            const LibUtilities::PointsKey Pkey3(nq3,Qtype3);
            const LibUtilities::BasisKey  Bkey1(btype1,order1,Pkey1);
            const LibUtilities::BasisKey  Bkey2(btype2,order2,Pkey2);
            const LibUtilities::BasisKey  Bkey3(btype3,order3,Pkey3);

            StdRegions::StdPyrExp *F = new StdRegions::StdPyrExp(Bkey1,Bkey2,Bkey3);
            E = F;
            pyrIdx = F->GetMap();
            vector<int> &rmap = F->GetRMap();
            E->GetCoords(x,y,z);

            //----------------------------------------------
            // Define solution to be projected
            for(i = 0; i < nq1*nq2*nq3; ++i)
            {
                sol[i]  = Tet_sol(x[i],y[i],z[i],order1,order2,order3);
            }

#if 0
            int nCoeffs = F->GetNcoeffs();
            int nPts = F->GetTotPoints();
            const Array<OneD, const NekDouble> &bx = F->GetBase()[0]->GetBdata();
            const Array<OneD, const NekDouble> &by = F->GetBase()[1]->GetBdata();
            const Array<OneD, const NekDouble> &bz = F->GetBase()[2]->GetBdata();

            for (int cnt = 0; cnt < nCoeffs; ++cnt)
            {
                StdPyrExp::triple &idx = pyrIdx[cnt];
                const int p = boost::get<0>(idx);
                const int q = boost::get<1>(idx);
                const int r = boost::get<2>(idx);
                if (r == 0 && p == 4 && q == 4)
                {
                    Array<OneD, NekDouble> asd(nPts);
                    for (int k = 0; k < nq3; ++k)
                    {
                        cout << z[k*nq1*nq2] << " " << bz[k + nq3*F->GetTetMode(3,2,1)] << endl;
                        /*
                        for (int j = 0; j < nq2; ++j)
                        {
                            for (int i = 0; i < nq1; ++i)
                            {
                                int s = i + nq1*(j + nq2*k);
                                asd[s] =
                                    bx[i + nq1*p]*
                                    by[j + nq2*q]*
                                    bz[k + nq3*rmap[cnt]]*
                                    bz[k + nq3*F->GetTetMode(q-1,0,0)];
                            }
                        }
                        */
                    }
                    printSolution(F, "quadmode.dat", x, y, z, asd);
                    exit(0);
                }
            }
            
#endif
            
#if 0
            int nCoeffs = F->GetNcoeffs();
            int nPts = F->GetTotPoints();
            Array<OneD, NekDouble> blah(nCoeffs), blah2(nPts);
            for (i = 0; i < nCoeffs; ++i)
            {
                StdPyrExp::triple &idx = pyrIdx[i];
                const int p = boost::get<0>(idx);
                const int q = boost::get<1>(idx);
                const int r = boost::get<2>(idx);

                Vmath::Zero(nCoeffs, blah, 1);
                blah[i] = 1.0;
                F->BwdTrans(blah, blah2);

                boost::format pad("mode-%03d.dat");
                pad % i;

                printSolution(F, pad.str(), x, y, z, blah2);
            }
            exit(0);
#endif

#if 0
            int nCoeffs = F->GetNcoeffs();
            int nPts = F->GetTotPoints();
            Array<OneD, NekDouble> blah(nCoeffs), blah2(nPts);
            for (i = 0; i < nCoeffs; ++i)
            {
                Vmath::Zero(nCoeffs, blah, 1);
                blah[i] = 1.0;
                F->BwdTrans(blah, blah2);
                F->FwdTrans(blah2, blah);
                int zeroCount = 0;
                int oneCount = 0;
                for (int j = 0; j < nCoeffs; ++j)
                {
                    if (fabs(blah[j]) < 1e-8)
                    {
                        zeroCount++;
                    }
                    if (fabs(blah[j]-1.0) < 1e-8)
                    {
                        oneCount++;
                    }
                }

                cout << i << " " << zeroCount << " " << oneCount << endl;
            }
            exit(0);
#endif

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
    case LibUtilities::eHexahedron:
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

    
    Array<OneD, NekDouble> phys (nq1*nq2*nq3);
    Array<OneD, NekDouble> coeffs(order1*order2*order3);

    //---------------------------------------------
    // Project onto Expansion
    E->FwdTrans(sol,coeffs);
    //---------------------------------------------

    //-------------------------------------------
    // Backward Transform Solution to get projected values
    E->BwdTrans(coeffs,phys);
    //-------------------------------------------

#if 0
    int nPts = E->GetTotPoints();
    Array<OneD, NekDouble> errArr(nPts);
    Vmath::Vsub(nPts, phys, 1, sol, 1, errArr, 1);
    printSolution(E, "out.dat", x, y, z, errArr);
    printSolution(E, "sol.dat", x, y, z, phys);
#endif

    //--------------------------------------------
    // Calculate L_inf error
    cout << "L infinity error: " << E->Linf(phys,sol) << endl;
    cout << "L 2 error:        " << E->L2  (phys,sol) << endl;
    //--------------------------------------------

    //-------------------------------------------
    // Evaulate solution at x = y = z = 0  and print error
    Array<OneD, NekDouble> t = Array<OneD, NekDouble>(3);
    t[0] = -0.5;
    t[1] = -0.25;
    t[2] = -0.3;

    if(regionshape == LibUtilities::eTetrahedron)
    {
        sol[0] = Tet_sol(t[0], t[1], t[2], order1, order2, order3);
    }
    else if (regionshape == LibUtilities::ePyramid)
    {
        sol[0] = Tet_sol(t[0], t[1], t[2], order1, order2, order3);
    }
    else if (regionshape == LibUtilities::ePrism)
    {
        sol[0] = Prism_sol(t[0], t[1], t[2], order1, order2, order3);
    }
    else
    {
        sol[0] = Hex_sol(t[0], t[1], t[2], order1, order2, order3,
                         btype1, btype2, btype3);
    }

    NekDouble nsol = E->PhysEvaluate(x,phys);
    cout << "error at x = (" << t[0] << "," << t[1] << "," << t[2] << "): "
         << nsol - sol[0] << endl;
    //-------------------------------------------

    return 0;
}


NekDouble Tet_sol(NekDouble x, NekDouble y, NekDouble z,
                  int order1, int order2, int order3)
{
    int    l,k,m;
    NekDouble sol = 0.0;

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
                sol += pow(z,m);
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

