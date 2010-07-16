#include <cstdio>
#include <cstdlib>

#include <LocalRegions/LocalRegions.hpp>
#include <LocalRegions/SegExp.h>
#include <SpatialDomains/SegGeom.h>

using namespace Nektar;

static double solutionpoly(double x, int order);
static double solutionfourier(double x, int order, double a, double b);

static double deriv_solutionpoly(double x, int order);
static double deriv_solutionfourier(double x, int order, double a, double b);

// This routine projects a polynomial or trigonmetric functions which 
// has energy in all mdoes of the expansions and report an error.

int main(int argc, char *argv[])
{
    int i;
    int order, nq;
    Array<OneD, const NekDouble> z,w;

    NekDouble x[2];
    LibUtilities::PointsType Qtype;
    LibUtilities::BasisType  btype;
    StdRegions::StdExpansion1D  *E;

    if(argc != 6)
    {
        fprintf(stderr,"Usage: Project1D Type order nq  x0 x1 \n");

        fprintf(stderr,"Where type is an integer value which "
            "dictates the basis as:\n");
        fprintf(stderr,"\t Ortho_A    = 1\n");
        fprintf(stderr,"\t Modified_A = 4\n");
        fprintf(stderr,"\t Fourier    = 7\n");
        fprintf(stderr,"\t Lagrange   = 8\n");
        fprintf(stderr,"\t Legendre   = 9\n"); 
        fprintf(stderr,"\t Chebyshev  = 10\n");

        fprintf(stderr,"Note type = 1,2,4,5 are for higher dimensional basis\n");

        exit(1);
    }

    btype =   (LibUtilities::BasisType) atoi(argv[1]);

    // Check to see that only 1D Expansions are used
    if((btype == LibUtilities::eOrtho_B)||(btype == LibUtilities::eOrtho_B)||
        (btype == LibUtilities::eModified_B)||(btype == LibUtilities::eModified_C))
    {
        NEKERROR(ErrorUtil::efatal,
            "This basis is for 2 or 3D expansions");
    }

    order =   atoi(argv[2]);
    nq    =   atoi(argv[3]);
    x[0]  =   atof(argv[4]);
    x[1]  =   atof(argv[5]);

    Array<OneD,NekDouble> sol(nq);

    if(btype != LibUtilities::eFourier)
    {
        Qtype = LibUtilities::eGaussLobattoLegendre;
    }
    else
    {
        Qtype = LibUtilities::eFourierEvenlySpaced;
    }

    //-----------------------------------------------
    // Define a segment expansion based on basis definition
    SpatialDomains::VertexComponentSharedPtr  vert1(new SpatialDomains::VertexComponent(1,0,x[0],0,0));
    SpatialDomains::VertexComponentSharedPtr  vert2(new SpatialDomains::VertexComponent(1,0,x[1],0,0));
    SpatialDomains::SegGeomSharedPtr geom(new SpatialDomains::SegGeom(0,vert1,vert2));
    geom->SetOwnData();

    const LibUtilities::PointsKey Pkey(nq,Qtype);
    const LibUtilities::BasisKey Bkey(btype,order,Pkey);
    E = new LocalRegions::SegExp(Bkey,geom);
    //-----------------------------------------------

    //----------------------------------------------
    // Define solution to be projected 
    Array<OneD,NekDouble> xc(nq);

    E->GetCoords(xc);

    if(btype != LibUtilities::eFourier)
    {
        for(i = 0; i < nq; ++i)
        {
            sol[i] = solutionpoly(xc[i],order);
        }
    }
    else
    {
        for(i = 0; i < nq; ++i)
        {
            sol[i] = solutionfourier(xc[i],order,x[0],x[1]);
        }
    }
    //---------------------------------------------

    //--------------------------------------------
    // Take the numerical derivative of the solutiion
    E->PhysDeriv(sol,sol);
    //---------------------------------------------

    //---------------------------------------------
    // Project onto Expansion 
    E->FwdTrans(sol,E->UpdateCoeffs());
    //---------------------------------------------
    
    //-------------------------------------------
    // Backward Transform Solution to get projected values
    E->BwdTrans(E->GetCoeffs(),E->UpdatePhys());
    //-------------------------------------------  

    //--------------------------------------------
    // Write solution 
    ofstream outfile("ProjectFile1D_Diff.dat");
    E->WriteToFile(outfile,eTecplot);
    outfile.close();
    //-------------------------------------------

    //-------------------------------------------------
    // Define derivative of the solution
    if(btype != LibUtilities::eFourier)
    {
        for(i = 0; i < nq; ++i)
        {
            sol[i] = deriv_solutionpoly(xc[i],order);
        }
    }
    else
    {
        for(i = 0; i < nq; ++i)
        {
            sol[i] = deriv_solutionfourier(xc[i],order,x[0],x[1]);
        }
    }
    //---------------------------------------------------

    //--------------------------------------------
    // Calculate L_inf error 
    cout << "L infinity error: " << E->Linf(sol) << endl;
    cout << "L 2 error:        " << E->L2  (sol) << endl;
    //--------------------------------------------

    return 0;
}

static double solutionpoly(double x, int order)
{
    int j;
    double sol;

    sol = 0.0;
    for(j = 0; j < order; ++j)
    {
        sol += pow(x,j);
    }

    return sol;
}

static double solutionfourier(double x, int order, double a, double b)
{
    int j;
    double sol;
    double xx = 2.0*(x-a)/(b-a) - 1.0;

    sol = 1.0;
    for(j = 0; j < order/2-1; ++j)
    {
        sol += sin(j*M_PI*xx) + cos(j*M_PI*xx);
    }

    return sol;
}


static double deriv_solutionpoly(double x, int order)
{
    int j;
    double sol;

    sol = 0.0;
    for(j = 1; j < order; ++j)
    {
        sol += j*pow(x,j-1);
    }

    return sol;
}

static double deriv_solutionfourier(double x, int order, double a, double b)
{
    int j;
    double sol;
    double xx = 2.0*(x-a)/(b-a) - 1.0;

    sol = 0.0;
    for(j = 0; j < order/2-1; ++j)
    {
        sol += j*M_PI*(cos(j*M_PI*xx) - sin(j*M_PI*xx));
    }

    return sol;
}
