#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdNodalTriExp.h>

#include <LibUtilities/Foundations/Foundations.hpp>

using namespace Nektar;

NekDouble Tri_sol(NekDouble x, NekDouble y, int order1, int order2);
NekDouble Quad_sol(NekDouble x, NekDouble y, int order1, int order2,
                   LibUtilities::BasisType btype1, LibUtilities::BasisType btype2);


// This routine projects a polynomial or trigonmetric functions which
// has energy in all mdoes of the expansions and reports and error

int main(int argc, char *argv[]){
    int           i;

    int           order1,order2, nq1,nq2;
    LibUtilities::PointsType    Qtype1,Qtype2;
    LibUtilities::BasisType     btype1 =   LibUtilities::eOrtho_A;
    LibUtilities::BasisType     btype2 =   LibUtilities::eOrtho_B;
    LibUtilities::PointsType    NodalType = LibUtilities::eNodalTriElec;
    LibUtilities::ShapeType     regionshape;
    StdRegions::StdExpansion *E;
    Array<OneD, NekDouble> sol;

    if(argc != 8)
    {
        fprintf(stderr,"Usage: StdProject2D RegionShape Type1 Type2 order1 "
                "order2  nq1 nq2  \n");

        fprintf(stderr,"Where RegionShape is an integer value which "
                "dictates the region shape:\n");
        fprintf(stderr,"\t Triangle      = 3\n");
        fprintf(stderr,"\t Quadrilateral = 4\n");


        fprintf(stderr,"Where type is an integer value which "
                "dictates the basis as:\n");

	  
        fprintf(stderr,"\t Ortho_A    = 1\n");
        fprintf(stderr,"\t Ortho_B    = 2\n");
        fprintf(stderr,"\t Modified_A = 4\n");
        fprintf(stderr,"\t Modified_B = 5\n");
        fprintf(stderr,"\t Fourier    = 7\n");
        fprintf(stderr,"\t Lagrange   = 8\n");
        fprintf(stderr,"\t Legendre   = 9\n");
        fprintf(stderr,"\t Chebyshev  = 10\n");
        fprintf(stderr,"\t FourierSingleMode  = 11\n");
        fprintf(stderr,"\t Nodal tri (Electro) = 12\n");
        fprintf(stderr,"\t Nodal tri (Fekete)  = 13\n");


        fprintf(stderr,"Note type = 3,6 are for three-dimensional basis\n");

        exit(1);
    }

    regionshape = (LibUtilities::ShapeType) atoi(argv[1]);

    // Check to see if 2D region
    if((regionshape != LibUtilities::eTriangle)&&(regionshape != LibUtilities::eQuadrilateral))
    {
        NEKERROR(ErrorUtil::efatal,"This shape is not a 2D region");
    }

    int btype1_val = atoi(argv[2]);
    int btype2_val = atoi(argv[3]);

    if(( btype1_val <= 11)&&( btype2_val <= 11))
    {
        btype1 =   (LibUtilities::BasisType) btype1_val;
        btype2 =   (LibUtilities::BasisType) btype2_val;
    }
    else if(( btype1_val >=12)&&(btype2_val <= 13))
    {
        btype1 =   LibUtilities::eOrtho_A;
        btype2 =   LibUtilities::eOrtho_B;

        if(btype1_val == 12)
        {
            NodalType = LibUtilities::eNodalTriElec;
        }
        else
        {
            NodalType = LibUtilities::eNodalTriFekete;
        }

    }

    // Check to see that correct Expansions are used
    switch(regionshape)
    {
    case LibUtilities::eTriangle:
        if((btype1 == LibUtilities::eOrtho_B)||(btype1 == LibUtilities::eModified_B))
        {
            NEKERROR(ErrorUtil::efatal,
                     "Basis 1 cannot be of type Ortho_B or Modified_B");
            }
        
        break;
    case LibUtilities::eQuadrilateral:
        if((btype1 == LibUtilities::eOrtho_B)||(btype1 == LibUtilities::eOrtho_C)||
           (btype1 == LibUtilities::eModified_B)||(btype1 == LibUtilities::eModified_C))
        {
            NEKERROR(ErrorUtil::efatal,
                     "Basis 1 is for 2 or 3D expansions");
        }
        
        if((btype2 == LibUtilities::eOrtho_B)||(btype2 == LibUtilities::eOrtho_C)||
           (btype2 == LibUtilities::eModified_B)||(btype2 == LibUtilities::eModified_C))
        {
            NEKERROR(ErrorUtil::efatal,
                     "Basis 2 is for 2 or 3D expansions");
        }
        break;
    default:
        ASSERTL0(false, "Not a 2D expansion.");
        break;
    }
    
    order1 =   atoi(argv[4]);
    order2 =   atoi(argv[5]);
    nq1    =   atoi(argv[6]);
    nq2    =   atoi(argv[7]);
    
    sol = Array<OneD, NekDouble>(nq1*nq2);
    
    if(btype1== LibUtilities::eFourier)
    {
        Qtype1 = LibUtilities::eFourierEvenlySpaced;
    }
    else if(btype1== LibUtilities::eFourierSingleMode)
    {
        Qtype1 = LibUtilities::eFourierSingleModeSpaced;	
    }
    else 
    {
        Qtype1 = LibUtilities::eGaussLobattoLegendre;
    }

    if(btype2== LibUtilities::eFourier)
    {
        Qtype2 = LibUtilities::eFourierEvenlySpaced;
		
    }
    else if(btype2== LibUtilities::eFourierSingleMode)
    {
        Qtype2 = LibUtilities::eFourierSingleModeSpaced;	
        
    }
    else
    {
        if (regionshape == LibUtilities::eTriangle) 
        {
            Qtype2 = LibUtilities::eGaussRadauMAlpha1Beta0;
        }
        else
        {
            Qtype2 = LibUtilities::eGaussLobattoLegendre;
        }
    }

    //-----------------------------------------------
    // Define a 2D expansion based on basis definition
    
    switch(regionshape)
    {
        case LibUtilities::eTriangle:
        {
            const LibUtilities::PointsKey Pkey1(nq1,Qtype1);
            const LibUtilities::PointsKey Pkey2(nq2,Qtype2);
            const LibUtilities::BasisKey  Bkey1(btype1,order1,Pkey1);
            const LibUtilities::BasisKey  Bkey2(btype2,order2,Pkey2);

            if(btype1_val >= 10)
            {
                E = new StdRegions::StdNodalTriExp(Bkey1,Bkey2,NodalType);
            }
            else
            {
                E = new StdRegions::StdTriExp(Bkey1,Bkey2);
            }

            Array<OneD,NekDouble> x = Array<OneD,NekDouble>(nq1*nq2);
            Array<OneD,NekDouble> y = Array<OneD,NekDouble>(nq1*nq2);
            E->GetCoords(x,y);

            //----------------------------------------------
            // Define solution to be projected
            for(i = 0; i < nq1*nq2; ++i)
            {
                sol[i]  = Tri_sol(x[i],y[i],order1,order2);
            }
            //----------------------------------------------
        }
        break;

        case LibUtilities::eQuadrilateral:
        {
            const LibUtilities::PointsKey Pkey1(nq1,Qtype1);
            const LibUtilities::PointsKey Pkey2(nq2,Qtype2);
            const LibUtilities::BasisKey Bkey1(btype1,order1,Pkey1);
            const LibUtilities::BasisKey Bkey2(btype2,order2,Pkey2);
            E = new StdRegions::StdQuadExp(Bkey1,Bkey2);

            //----------------------------------------------
            // Define solution to be projected

            Array<OneD, NekDouble> x = Array<OneD, NekDouble>(nq1*nq2);
            Array<OneD, NekDouble> y = Array<OneD, NekDouble>(nq1*nq2);
            E->GetCoords(x,y);

            for(i = 0; i < nq1*nq2; ++i)
            {
                sol[i]  = Quad_sol(x[i],y[i],order1,order2,btype1,btype2);
            }
            //---------------------------------------------
        }
        break;
        
        default:
            ASSERTL0(false, "Not a 2D expansion.");
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
	
    // Evaulate solution at x = y =0  and print error
    Array<OneD, NekDouble> x = Array<OneD, NekDouble>(2);
    x[0] = 0;
    x[1] = -0.25;

    if(regionshape == LibUtilities::eTriangle)
    {
        sol[0] = Tri_sol(x[0],x[1],order1,order2);
    }
    else
    {
        sol[0] = Quad_sol(x[0],x[1],order1,order2,btype1,btype2);
    }

    NekDouble nsol = E->PhysEvaluate(x);
    cout << "error at x = (" <<x[0] <<","<<x[1] <<"): " << nsol - sol[0] << endl;
    //-------------------------------------------

    return 0;
}

NekDouble Tri_sol(NekDouble x, NekDouble y, int order1, int order2){
    int    l,k;
    NekDouble sol = 0;

    for(k = 0; k < order1; ++k)
    {
        for(l = 0; l < order2-k; ++l)
        {
            sol += pow(x,k)*pow(y,l);
        }
    }

    return sol;
}

NekDouble Quad_sol(NekDouble x, NekDouble y, int order1, int order2,
                   LibUtilities::BasisType btype1,
                   LibUtilities::BasisType btype2)
{
    int k,l;
    NekDouble sol = 0.0;
    if(btype1 == LibUtilities::eFourier)
    {
        if(btype2 == LibUtilities::eFourier)
        {	
            for(k = 0; k < order1/2; ++k)
            {
                for(l = 0; l < order2/2; ++l)
                {
                    sol += sin(M_PI*k*x)*sin(M_PI*l*y)
                        + sin(M_PI*k*x)*cos(M_PI*l*y)
                        + cos(M_PI*k*x)*sin(M_PI*l*y)
                        + cos(M_PI*k*x)*cos(M_PI*l*y);
                }
            }
        }
        else if(btype2 == LibUtilities::eFourierSingleMode)
        {
            for(k = 0; k < order1/2; ++k)
            {
                sol += sin(M_PI*k*x)*sin(M_PI*y)
                    + sin(M_PI*k*x)*cos(M_PI*y)
                    + cos(M_PI*k*x)*sin(M_PI*y)
                    + cos(M_PI*k*x)*cos(M_PI*y);
            }			
        }
        else
        {
            for(k = 0; k < order1/2; ++k)
            {
                for(l = 0; l < order2; ++l)
                {
                    sol += sin(M_PI*k*x)*pow(y,l) + cos(M_PI*k*x)*pow(y,l);
                }
            }
        }
    }
    else if(btype1 == LibUtilities::eFourierSingleMode)
    {
        if(btype2 == LibUtilities::eFourier)
        {	
            for(l = 0; l < order2/2; ++l)
            {
                sol += sin(M_PI*x)*sin(M_PI*l*y)
                    + sin(M_PI*x)*cos(M_PI*l*y)
                    + cos(M_PI*x)*sin(M_PI*l*y)
                    + cos(M_PI*x)*cos(M_PI*l*y);
            }
        }
        else if(btype2 == LibUtilities::eFourierSingleMode)
        {
            sol += sin(M_PI*x)*sin(M_PI*y)
                + sin(M_PI*x)*cos(M_PI*y)
                + cos(M_PI*x)*sin(M_PI*y)
                + cos(M_PI*x)*cos(M_PI*y);
        }			
        else
        {
            for(l = 0; l < order2; ++l)
            {
                sol += sin(M_PI*x)*pow(y,l) + cos(M_PI*x)*pow(y,l);
            }
        }
		
    }
    //case for Non Fourier expansion
    else
    {		
        if(btype2 == LibUtilities::eFourier)
        {	
            for(k = 0; k < order1; ++k)
            {	
                for(l = 0; l < order2/2; ++l)
                {
                    sol += pow(x,k)*sin(M_PI*l*y) + pow(x,k)*cos(M_PI*l*y);
                }
            }
        }
        else if(btype2 == LibUtilities::eFourierSingleMode)
        {
            for(k = 0; k < order1; ++k)
            {	

                sol += pow(x,k)*sin(M_PI*y) + pow(x,k)*cos(M_PI*y);
            }
        }
        else 
        {
            for(k = 0; k < order1; ++k)
            {
                for(l = 0; l < order2; ++l)
                {
                    sol += pow(x,k)*pow(y,l);
                }
            }
        }
    }
	
    return sol;
}
