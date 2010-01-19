#include <cstdio>
#include <cstdlib>
#include "StdRegions/StdExpansion.h"
#include "StdRegions/StdHexExp.h"
#include "StdRegions/StdTetExp.h"
#include "StdRegions/StdPrismExp.h"
#include "StdRegions/StdPyrExp.h"

#include "StdRegions/StdRegions.hpp"

using namespace Nektar;
using namespace StdRegions; 
using namespace std;

// compile using Builds/Demos/StdRegions -> make DEBUG=1
// PROG=MMatrix3D //Gosse: should be: make DEBUG=1 MMatrix3D

main(int argc, char *argv[])
{
  int           i,j;
  int           order1,order2,order3,nq1,nq2,nq3;
  PointsType     Qtype1,Qtype2,Qtype3;
  BasisType     btype1,btype2,btype3;
  ShapeType     regionshape;
  StdExpansion   *E;
  
  if(argc != 11)
  {
    fprintf(stderr,"Usage: MMatrix3D RegionShape Type1 Type2 Type3 order1 "
        "order2  order3 nq1 nq2 nq3 \n");

    fprintf(stderr,"Where RegionShape is an integer value which "
        "dictates the region shape:\n");
    fprintf(stderr,"\t Tetrahedron      = 4\n"); 
    fprintf(stderr,"\t Pyramid          = 5\n");
    fprintf(stderr,"\t Prism            = 6\n");
    fprintf(stderr,"\t Hexahedron       = 7\n");
    
    
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
    
    exit(1);
  }
  
  regionshape = (ShapeType) atoi(argv[1]);
  
  // Check to see if 3D region 
  if((regionshape != eTetrahedron)&&(regionshape != eHexahedron)&&
     (regionshape != ePrism)&&(regionshape != ePyramid))
  {
    ErrorUtil::Error(ErrorUtil::efatal,"MMatrix3D","This shape is "
             "not a 3D region");
  }

  btype1 =   (BasisType) atoi(argv[2]);
  btype2 =   (BasisType) atoi(argv[3]);
  btype3 =   (BasisType) atoi(argv[4]);                                       
  
  // Check to see that correct Expansions are used

  order1 =   atoi(argv[5]);
  order2 =   atoi(argv[6]);
  order3 =   atoi(argv[7]);

  nq1    =   atoi(argv[8]);
  nq2    =   atoi(argv[9]);
  nq3    =   atoi(argv[9]);

  if(btype1 != eFourier)
  {
    Qtype1 = eLobatto; 
  }
  else
  {
    Qtype1 = eFourierEvenSp;
  }

  if(btype2 != eFourier)
  {
    Qtype2 = eLobatto; 
  }
  else
  {
    Qtype2 = eFourierEvenSp;
  }

  if(btype3 != eFourier)
  {
    Qtype3 = eLobatto; 
  }
  else
  {
    Qtype3 = eFourierEvenSp;
  }
  
  
  //-----------------------------------------------
  // Define a segment expansion based on basis definition

  switch(regionshape)
  {
  case eTetrahedron:
    fprintf(stderr,"Tetrahedron is not yet set up\n");
    exit(1);
    break;

  case ePyramid:
    fprintf(stderr,"Pyramid is not yet set up\n");
    exit(1);
    break;

  case ePrism:
    fprintf(stderr,"Prism is not yet set up\n");
    exit(1);
    break;

  case eHexahedron:
  {
    const BasisKey b1(btype1,order1,Qtype1,nq1,0,0);
    const BasisKey b2(btype2,order2,Qtype2,nq2,0,0);
    const BasisKey b3(btype3,order3,Qtype3,nq3,0,0);
    E = new StdHexExp (b1,b2,b3);
    break;
  }
  }

  
  // ---------------------------------------------
  // Generate mass matrix based upon basis definition and dump to stdout
  (E->GetMassMatrix())->DumpMatrix(stdout);
  //--------------------------------------------


  // ---------------------------------------------
  // Show Matrix structure 
  (E->GetMassMatrix())->ShowMatrixStructure(stdout);
  //--------------------------------------------
}
