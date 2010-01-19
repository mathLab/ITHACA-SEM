#include <cstdio>
#include <cstdlib>
#include "StdRegions/StdExpansion.h"
#include "StdRegions/StdQuadExp.h"
#include "StdRegions/StdTriExp.h"

#include "StdRegions/StdRegions.hpp"

using namespace Nektar;
using namespace StdRegions; 
using namespace std;

// compile using Builds/Demos/StdRegions -> make DEBUG=1
// PROG=MMatrix2D  //Gosse: "PROG=" should be removed

main(int argc, char *argv[])
{
  int           i,j;
  int           order1,order2, nq1,nq2;
  PointsType     Qtype1,Qtype2;
  BasisType     btype1,btype2;
  ShapeType     regionshape;
  StdExpansion   *E;
  
  if(argc != 8)
  {
    fprintf(stderr,"Usage: MMatrix2D RegionShape Type1 Type2 order1 "
        "order2  nq1 nq2  \n");

    fprintf(stderr,"Where RegionShape is an integer value which "
        "dictates the region shape:\n");
    fprintf(stderr,"\t Triangle      = 1\n"); //Gosse: should be:
                          //Triangle =2
    fprintf(stderr,"\t Quadrilateral = 2\n"); //Gosse: should be:
                          //Quadrilateral = 3
    
    
    fprintf(stderr,"Where type is an integer value which "
        "dictates the basis as:\n");

    fprintf(stderr,"\t Ortho_A    = 0\n");
    fprintf(stderr,"\t Ortho_B    = 1\n");
    fprintf(stderr,"\t Modified_A = 3\n");
    fprintf(stderr,"\t Modified_B = 4\n");
    fprintf(stderr,"\t Fourier    = 6\n");
    fprintf(stderr,"\t Lagrange   = 7\n");
    fprintf(stderr,"\t Legendre   = 8\n"); 
    fprintf(stderr,"\t Chebyshev  = 9\n");
 
    fprintf(stderr,"Note type = 2,5 are for three-dimensional basis\n");

    exit(1);
  }
  
  regionshape = (ShapeType) atoi(argv[1]);
  
  // Check to see if 2D region 
  if((regionshape != eTriangle)&&(regionshape != eQuadrilateral))
  {
    ErrorUtil::Error(ErrorUtil::efatal,"MMatrix2D","This shape is not a 2D region");
  }

  btype1 = (BasisType) atoi(argv[2]);
  btype2 = (BasisType) atoi(argv[3]);

  // Check to see that correct Expansions are used
  switch(regionshape)
  {
  case eTriangle:
    if((btype1 == eOrtho_B)||(btype1 == eModified_B))
    {
      ErrorUtil::Error(ErrorUtil::efatal,"MMatrix2D",
               "Basis 1 cannot be of type Ortho_B or Modified_B");
    }

    if((btype2 != eOrtho_B)&&(btype2 != eModified_B))
    {
      ErrorUtil::Error(ErrorUtil::efatal,"MMatrix2D",
               "Basis 2 must be of type Ortho_B or Modified_B");
    }
    break;
  case eQuadrilateral:
    if((btype1 == eOrtho_B)||(btype1 == eOrtho_B)||
       (btype1 == eModified_B)||(btype1 == eModified_C))
    {
      ErrorUtil::Error(ErrorUtil::efatal,"MMatrix2D",
               "Basis 1 is for 2 or 3D expansions");
    }

    if((btype2 == eOrtho_B)||(btype2 == eOrtho_B)||
       (btype2 == eModified_B)||(btype2 == eModified_C))
     {
      ErrorUtil::Error(ErrorUtil::efatal,"MMatrix2D",
             "Basis 2 is for 2 or 3D expansions");
     }
    break;
  }
  
  order1 =   atoi(argv[4]);
  order2 =   atoi(argv[5]);
  nq1    =   atoi(argv[6]);
  nq2    =   atoi(argv[7]);

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
  
  //-----------------------------------------------
  // Define a segment expansion based on basis definition

  switch(regionshape)
  {
  case eTriangle:
  {
    const BasisKey b1(btype1,order1,Qtype1,nq1,0,0);
    const BasisKey b2(btype2,order2,Qtype2,nq2,1,0);
    E = new StdTriExp (b1,b2);
    break;
  }
  case eQuadrilateral:
  {
    const BasisKey b1(btype1,order1,Qtype1,nq1,0,0);
    const BasisKey b2(btype2,order2,Qtype2,nq2,0,0);
    E = new StdQuadExp (b1,b2);
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
