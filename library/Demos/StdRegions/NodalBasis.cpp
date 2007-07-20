#include <cstdio>
#include <cstdlib>
#include "StdRegions/StdExpansion2D.h"
#include "StdRegions/StdNodalTriExp.h"
#include "StdRegions/NodalBasisManager.h"

#include "StdRegions/StdRegions.hpp"

using namespace Nektar;
using namespace StdRegions; 
using namespace  std;


// compile using Builds/Demos/StdRegions -> make DEBUG=1
// PROG=NodalBasis //Gosse: I don't think "PROG=" has to be there

// This routine projects a polynomial or trigonmetric functions which 
// has energy in all mdoes of the expansions and report an error.

int main(int argc, char *argv[])
{
  int i,j;
  int order;
  NodalBasisType btype;
  NodalBasisManager B;
  

  if(argc != 3)
  {
    fprintf(stderr,"Usage: NodalBasis Type order\n");

    fprintf(stderr,"Where type is an integer value which "
        "dictates the basis as:\n");
    fprintf(stderr,"\t NodalTriElec    = 0\n");
    fprintf(stderr,"\t NodalTriFekete  = 1\n");
    fprintf(stderr,"\t NodalTetElec    = 2\n");
    exit(1);
  }
  
  btype =   (NodalBasisType) atoi(argv[1]);
  order =   atoi(argv[2]);
 
  int npts;
  const double *x,*y,*z;
  
  //----------------------------------------------
  // Output 1D basis using only basis information. 

  fprintf(stdout,"VARIABLES = x y\n");

  for(i = 2; i <= order; ++i)
  {
    npts = B.GetNodePoints(btype,i,x,y,z);
  
    fprintf(stdout,"ZONE T = \"Order %d\" I=%d\n",i,npts);

    for(j = 0; j < npts; ++j)
    {
      fprintf(stdout, "%9.6lf %9.6lf \n",x[j],y[j]);
    }
  }
  //-----------------------------------------------
}
