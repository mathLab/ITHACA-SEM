#include <cstdio>
#include <cstdlib>
#include "StdRegions/StdSegExp.h"

#include "StdRegions/StdRegions.hpp"

using namespace Nektar;
using namespace StdRegions; 
using namespace std;

// compile using Builds/Demos/StdRegions -> make DEBUG=1 Project1D

// This routine projects a polynomial or trigonmetric functions which 
// has energy in all mdoes of the expansions and report an error.

main(int argc, char *argv[])
{
  int i,j;
  int order, nq;
  const double *z,*w;
  double *sol,*phys,L2_err;
  PointsType Qtype;
  BasisType  btype;
  StdExpansion1D  *E;
  
  if(argc != 4)
  {
    fprintf(stderr,"Usage: Project1D Type order nq \n");

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
  E = new StdSegExp(b1);
  //-----------------------------------------------
  
  //----------------------------------------------
  // Define solution to be projected
  E->GetZW(0,z,w);
  if(btype != eFourier)
  {
    for(i = 0; i < nq; ++i)
    {
      sol[i] = 0.0;
      for(j = 0; j < order; ++j)
      {
	sol[i] += pow(z[i],j);
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
	sol[i] += sin(j*M_PI*z[i]) + cos(j*M_PI*z[i]);
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
  // Calculate L_inf error 
  cout << "L infinity error: " << E->Linf(sol) << endl;
  cout << "L 2 error:        " << E->L2  (sol) << endl;
  //--------------------------------------------

  //-------------------------------------------
  // Evaulate solution at mid point and print error
  if(btype != eFourier)
  {
    sol[0] = 1;
  }
  else
  {
    sol[0] =  order/2;
  }

  double x = 0;
  double nsol = E->Evaluate(&x);
  cout << "error at x = 0: " << nsol - sol[0] << endl;
  //-------------------------------------------

  delete[] sol; 
}
