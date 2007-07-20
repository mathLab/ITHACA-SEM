#include <cstdio>
#include <cstdlib>
#include "StdRegions/StdSegExp.h"

#include "StdRegions/StdRegions.hpp"

using namespace Nektar; 
using namespace StdRegions; 
using namespace std;

// compile using Builds/Demos/StdRegions -> make DEBUG=1 PROG=MMatrix1D

main(int argc, char *argv[]){
  int i,j;
  int order, nq;
  const double *z,*w;
  PointsType Qtype;
  BasisType  btype;
  StdSegExp  *E;
  
  if(argc != 4)
  {
    fprintf(stderr,"Usage: MMatrix1D Type order nq \n");

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
    ErrorUtil::Error(ErrorUtil::efatal,"MMatrix1D",
             "This basis is for 2 or 3D expansions");
  }
  
  order =   atoi(argv[2]);
  nq    =   atoi(argv[3]);

  if(btype != eFourier)
  {
    Qtype = eLobatto; 
  }
  else
  {
    Qtype = eFourierEvenSp;
  }
  
 
  BasisKey BK(btype,order,Qtype,nq,0,0);

  //-----------------------------------------------
  // Define a segment expansion based on basis definition
  E = new StdSegExp (BK);
  //-----------------------------------------------
  
  // ---------------------------------------------
  // Generate mass matrix based upon basis definition and dump to stdout
  (E->GetMassMatrix())->DumpMatrix(stdout);
  //--------------------------------------------

  // ---------------------------------------------
  // Show Matrix structure 
  (E->GetMassMatrix())->ShowMatrixStructure(stdout);
  //--------------------------------------------
}
