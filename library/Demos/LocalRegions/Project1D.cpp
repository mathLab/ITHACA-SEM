#include <cstdio>
#include <cstdlib>

#include <LocalRegions/LocalRegions.hpp>
#include <LocalRegions/SegExp.h>
#include <SpatialDomains/SegGeom.h>

using namespace Nektar;
using namespace StdRegions; 
using namespace LocalRegions; 
using namespace SpatialDomains; 
using namespace std;
  
// compile using Builds/Demos/StdRegions -> make DEBUG=1 ProjectS1D

// This routine projects a polynomial or trigonmetric functions which 
// has energy in all mdoes of the expansions and report an error.

main(int argc, char *argv[]){
  int i,j;
  int order, nq;
  const double *z,*w;
  double *sol,*phys,L2_err;
  double x[2];
  PointsType Qtype;
  BasisType  btype;
  StdExpansion1D  *E;
  
  if(argc != 6){
    fprintf(stderr,"Usage: Project1D Type order nq  x0 x1 \n");

    fprintf(stderr,"Where type is an integer value which "
	    "dictates the basis as:\n");
    fprintf(stderr,"\t Ortho_A    = 0\n");
    fprintf(stderr,"\t Modified_A = 3\n");
    fprintf(stderr,"\t Fourier    = 6\n");
    fprintf(stderr,"\t Lagrange   = 7\n");
    fprintf(stderr,"\t Legendre   = 8\n"); 
    fprintf(stderr,"\t Chebyshev  = 9\n");
 
    fprintf(stderr,"Note type = 1,2,4,5 are for higher dimensional basis\n");

    exit(1);
  }
  
  btype =   (BasisType) atoi(argv[1]);
  
  // Check to see that only 1D Expansions are used
  if((btype == eOrtho_B)||(btype == eOrtho_B)||
     (btype == eModified_B)||(btype == eModified_C))
  {
    ErrorUtil::Error(ErrorUtil::efatal,__FILE__,__LINE__,
		     "This basis is for 2 or 3D expansions");
  }
  
  order =   atoi(argv[2]);
  nq    =   atoi(argv[3]);
  x[0]  =   atof(argv[4]);
  x[1]  =   atof(argv[5]);

  sol = new double [nq];
  
  if(btype != eFourier)
  {
    Qtype = eLobatto; 
  }
  else
  {
    Qtype = eFourierEvenSp;
  }
  
  //-----------------------------------------------
  // Define a segment expansion based on basis definition
  const BasisKey b1(btype,order, Qtype, nq,0,0);

  VertexComponentSharedPtr  vert1(new VertexComponent(1,0,x[0],0,0));
  VertexComponentSharedPtr  vert2(new VertexComponent(1,0,x[1],0,0));
  SegGeomSharedPtr geom(new SegGeom(0,vert1,vert2));
  geom->SetOwnData();

  E = new SegExp(b1,geom);
  E->SetGeoFac(E->GenGeoFac());
  
  //-----------------------------------------------
  
  //----------------------------------------------
  // Define solution to be projected 
  E->GetZW(0,z,w); 

  double *xc = new double [nq];    
  E->GetCoords(&xc);
  
  if(btype != eFourier)
  {
    for(i = 0; i < nq; ++i)
    {
      sol[i] = 0.0;
      for(j = 0; j < order; ++j)
      {
	sol[i] += pow(xc[i],j);
      }
    }
  }
  else
  {
    for(i = 0; i < nq; ++i)
    {
      sol[i] = 0.0;
      for(j = 0; j < order/2; ++j)
      {
	sol[i] += sin(j*M_PI*xc[i]/(x[1]-x[0])) +
	  cos(j*M_PI*xc[i]/(x[1]-x[0]));
      }
    }
  }
  //---------------------------------------------

  //---------------------------------------------
  // Project onto Expansion 
  E->FwdTrans(sol);
  //---------------------------------------------

  //-------------------------------------------
  // Backward Transform Solution to get projected values
  E->BwdTrans(E->GetPhys());
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
  //xc[0] = 0.5*(x[1]+x[0]);
  double lcoord = 0;
  
  E->GetCoord(&lcoord,xc);
  
  if(btype != eFourier)
  {
    sol[0] = 0.0;
    for(j = 0; j < order; ++j)
    {
      sol[0] += pow(xc[0],j);
    }
  }
  else
  {
    sol[0] = 0.0;
    for(j = 0; j < order/2; ++j)
    {
      sol[0] += sin(j*M_PI*xc[0]/(x[1]-x[0])) +
	cos(j*M_PI*0.5*(x[1]+x[0])/(x[1]-x[0]));
    }
  }

  double nsol = E->Evaluate(xc);
  cout << "error at (xi = 0) x = " << xc[0] << " : " << nsol - sol[0] << endl;

  //-------------------------------------------

  delete[] sol; 
}
