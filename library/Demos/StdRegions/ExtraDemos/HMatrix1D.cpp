#include <cstdio>
#include <cstdlib>
#include "StdRegions.h"
#include "StdSegExp.h"

using namespace StdRegions; 
using namespace std;

/* 
   g++ -g -c HMatrix1D.cpp  -I ../../../include
// Linux default
   g++ -g -o HMatrix1D HMatrix1D.o -L../../  -lStdRegions -lPolylib -lm  -lblas -lg2c
// MaxOSX default
   g++ -g -o HMatrix1D HMatrix1D.o -L../../  -lStdRegions -lPolylib -lm  -framework Accelerate 
*/

int main(int argc, char *argv[]){
  int i,j;
  int order, nq;
  const double *z,*w;
  double *hmat;
  PointType Qtype;
  BasisType btype;
  StdSegExp  *E;
  
  if(argc != 4){
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
  if((btype == Ortho_B)   ||(btype == Ortho_B)||
     (btype == Modified_B)||(btype == Modified_C))
    ErrorUtil::Error(fatal,"HMatrix1D",
             "This basis is for 2 or 3D expansions");
  
  order =   atoi(argv[2]);
  nq    =   atoi(argv[3]);

  if(btype != Fourier)
    Qtype = Lobatto; 
  else
    Qtype = FourierEvenSp;


  BasisKey BK(btype,order,Qtype,nq,0,0);

  //-----------------------------------------------
  // Define a segment expansion based on basis definition
  E = new StdSegExp(BK);
  //-----------------------------------------------
  
  // ---------------------------------------------
  // Generate mass matrix based upon basis definition and put in mmat
  hmat = (E->GetLapMatrix())->get_matrix();
  //--------------------------------------------

  //----------------------------------------------
  // Dump mass matrix 

  fprintf(stdout,"VARIABLES = hmat\n");

  fprintf(stdout,"ZONE T = \"Helmholtz Matrix order %d\" I=%d J=%d\n",order,order,
        order);
  for(i = 0; i < order; ++i){
    for(j = 0; j < order; ++j)
      fprintf(stdout, "%9.6lf ",hmat[i*order+j]);
    fputc('\n',stdout);
  }
  //-----------------------------------------------
 return 0;   
}
