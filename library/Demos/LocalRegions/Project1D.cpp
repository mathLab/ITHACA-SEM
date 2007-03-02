#include <cstdio>
#include <cstdlib>

#include <LocalRegions/LocalRegions.hpp>
#include <LocalRegions/SegExp.h>
#include <SpatialDomains/SegGeom.h>

using namespace Nektar;
  
static double solution(double x, int order, LibUtilities::BasisType btype);

// This routine projects a polynomial or trigonmetric functions which 
// has energy in all mdoes of the expansions and report an error.

main(int argc, char *argv[])
{
  int i,j;
  int order, nq;
  const double *z,*w;
  double *sol,*phys,L2_err;
  double x[2];
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
    ErrorUtil::Error(ErrorUtil::efatal,__FILE__,__LINE__,
		     "This basis is for 2 or 3D expansions");
  }
  
  order =   atoi(argv[2]);
  nq    =   atoi(argv[3]);
  x[0]  =   atof(argv[4]);
  x[1]  =   atof(argv[5]);

  BstShrDArray stmp = GetDoubleTmpSpace(nq);
  sol = stmp.get();
  
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
  E->SetMinfo(E->GenMinfo());
  
  //-----------------------------------------------
  
  //----------------------------------------------
  // Define solution to be projected 
  BstShrDArray tmp = GetDoubleTmpSpace(nq);
  double  *xc =  tmp.get();
  E->GetCoords(&xc);
  for(i = 0; i < nq; ++i)
  {
      sol[i] = solution(xc[i],order,btype);
  }

  //---------------------------------------------

  //---------------------------------------------
  // Project onto Expansion 
  E->FwdTrans(sol);
  //---------------------------------------------

  //-------------------------------------------
  // Backward Transform Solution to get projected values
  E->BwdTrans(&E->GetPhys()[0]);
  //-------------------------------------------  

  //--------------------------------------------
  // Write solution 
  FILE *outfile = fopen("ProjectFile1D.dat","w");
  E->WriteToFile(outfile);
  fclose(outfile);
  //-------------------------------------------

  //--------------------------------------------
  // Calculate L_inf error 
  cout << "L infinity error: " << E->Linf(sol) << endl;
  cout << "L 2 error:        " << E->L2  (sol) << endl;
  //--------------------------------------------

  //-------------------------------------------
  // Evaulate solution at mid point and print error
  sol[0] = solution(0.5*(x[1]+x[0]),order,btype);
  double lcoord = 0;
  E->GetCoord(&lcoord,xc);
  double nsol = E->Evaluate(xc);
  cout << "error at (xi = 0) x = " << xc[0] << " : " << nsol - sol[0] << endl;

  //-------------------------------------------
}

static double solution(double x, int order, LibUtilities::BasisType btype)
{
    int j;
    double sol;
    
    if(btype != LibUtilities::eFourier)
    {
	sol = 0.0;
	for(j = 0; j < order; ++j)
	{
	    sol += pow(x,j);
	}
    }
    else
    {
	sol = 0.0;
	for(j = 0; j < order/2; ++j)
	{
	    sol += sin(j*M_PI*x) + cos(j*M_PI*0.5*x);
	}
    }
}
