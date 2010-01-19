#include <cstdio>
#include <cstdlib>
#include "StdRegions/BasisManager.h"

#include "StdRegions/StdRegions.hpp"

using namespace Nektar;
using namespace StdRegions; 
using namespace  std;

// compile using Builds/Demos/StdRegions -> make DEBUG=1 Basis1D

// This routine generates a tecplot output file of a 1D basis

main(int argc, char *argv[])
{
  int i,j;
  int order, nq;
  const double *z,*w;
  const double *b;
  BasisType btype;
  BasisManager B;
  

  if(argc != 4)
  {
    fprintf(stderr,"Usage: Basis1D Type order nq \n");
    
    fprintf(stderr,"Where type is an integer value which "
        "dictates the basis as:\n");
    fprintf(stderr,"\t Ortho_A    = 0\n");
    fprintf(stderr,"\t Ortho_B    = 1\n");
    fprintf(stderr,"\t Ortho_C    = 2\n");
    fprintf(stderr,"\t Modified_A = 3\n");
    fprintf(stderr,"\t Modified_B = 4\n");
    fprintf(stderr,"\t Modified_C = 5\n");
    fprintf(stderr,"\t Fourier    = 6\n");
    fprintf(stderr,"\t Lagrange   = 7\n");
    fprintf(stderr,"\t Legendre   = 8\n"); 
    fprintf(stderr,"\t Chebyshev  = 9\n");
    
    fprintf(stderr,"Note type = 0,1,2,4,5 are for higher dimensional basis\n");

    exit(1);
  }
  
  btype =   (BasisType) atoi(argv[1]);
  order =   atoi(argv[2]);
  nq    =   atoi(argv[3]);

  //-----------------------------------------------
  // Define points at which basis is to be evaluated
  B.GetZW(StdRegions::eLobatto,nq,z,w,0,0);
  //-----------------------------------------------

  
  //----------------------------------------------
  // Generate basis 
  b = B.GetBasisArray(btype, order, eLobatto, nq, 0.0, 0.0);
  //----------------------------------------------


  //----------------------------------------------
  // Output 1D basis using only basis information. 

  fprintf(stdout,"VARIABLES = z phi(z)\n");

  for(i = 0; i < order; ++i)
  {
    fprintf(stdout,"ZONE T = \"Mode %d\" I=%d\n",i+1,nq);
    for(j = 0; j < nq; ++j)
    {
      fprintf(stdout, "%9.6lf %9.6lf \n",z[j],b[i*nq+j]);
    }
  }
  //-----------------------------------------------
    
}
