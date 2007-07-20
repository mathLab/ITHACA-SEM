#include <cstdio>
#include <cstdlib>
#include <fstream>
#include "StdRegions/StdExpansion2D.h"
#include "StdRegions/StdQuadExp.h"
#include "StdRegions/StdTriExp.h"

#include "StdRegions/StdRegions.hpp"

using namespace Nektar;
using namespace StdRegions; 
using namespace std;

double Tri_sol(double x, double y, int order1, int order2);
double Quad_sol(double x, double y, int order1, int order2, 
        BasisType btype1, BasisType btype2);

// compile using Builds/Demos/StdRegions -> make DEBUG=1 AliasingProject2D

// This routine projects a polynomial or trigonmetric functions which 
// has energy in all mdoes of the expansions and reports and error

main(int argc, char *argv[])
{
  int            i,j,k,l;
  const          double *z1,*z2,*w;
  int            order1,order2, nq1,nq2;
  PointsType     Qtype1,Qtype2;
  BasisType      btype1,btype2;
  ShapeType      regionshape;
  StdExpansion2D *E;
  double         *sol;
  
  if(argc != 8)
  {
    fprintf(stderr,"Usage: AliasingProject2D RegionShape Type1 Type2 exporder1"
        "exporder2  nq1 nq2  \n");

    fprintf(stderr,"Where RegionShape is an integer value which "
        "dictates the region shape:\n");
    fprintf(stderr,"\t Triangle      = 2\n");
    fprintf(stderr,"\t Quadrilateral = 3\n");
    
    
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
    ErrorUtil::Error(ErrorUtil::efatal,"AliasingProject2D",
             "This shape is not a 2D region");
  }

  btype1 =   (BasisType) atoi(argv[2]);
  btype2 =   (BasisType) atoi(argv[3]);

  // Check to see that correct Expansions are used
  switch(regionshape)
  {
  case eTriangle:
    if((btype1 == eOrtho_B)||(btype1 == eModified_B))
    {
      ErrorUtil::Error(ErrorUtil::efatal,"Project2D",
               "Basis 1 cannot be of type Ortho_B or Modified_B");
    }

    if((btype2 != eOrtho_B)&&(btype2 != eModified_B))
    {
      ErrorUtil::Error(ErrorUtil::efatal,"Project2D",
               "Basis 2 must be of type Ortho_B or Modified_B");
    }
    break;
  case eQuadrilateral:
    if((btype1 == eOrtho_B)||(btype1 == eOrtho_B)||
       (btype1 == eModified_B)||(btype1 == eModified_C))
    {
      ErrorUtil::Error(ErrorUtil::efatal,"Project2D",
             "Basis 1 is for 2 or 3D expansions");
    }

    if((btype2 == eOrtho_B)||(btype2 == eOrtho_B)||
       (btype2 == eModified_B)||(btype2 == eModified_C))
    {
      ErrorUtil::Error(ErrorUtil::efatal,"Project2D",
             "Basis 2 is for 2 or 3D expansions");
    }
    break;
  }
  
  order1 = atoi(argv[4]);
  order2 = atoi(argv[5]);
  nq1    = atoi(argv[6]);
  nq2    = atoi(argv[7]);

  sol    = new double [nq1*nq2];

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
  case eTriangle:{
    const BasisKey b1(btype1,order1,Qtype1,nq1,0,0);
    const BasisKey b2(btype2,order2,eRadauM,nq2,1,0);

    E = new StdTriExp (b1,b2);
    
    E->GetZW(0,z1,w);
    E->GetZW(1,z2,w);

    //----------------------------------------------
    // Define solution to be projected
    double x,y;
    for(i = 0; i < nq1; ++i)
    {
      for(j = 0; j < nq2; ++j)
      {
    x = (1+z1[i])*(1-z2[j])/2-1.0;
    y = z2[j];
    sol[i+nq1*j]  = Tri_sol(x,y,order1,order2);
    sol[i+nq1*j] *= sol[i+nq1*j]; 
      }
    }
    //----------------------------------------------
  }
    break;
  case eQuadrilateral:
  {
    const BasisKey b1(btype1,order1,Qtype1,nq1,0,0);
    const BasisKey b2(btype2,order2,Qtype2,nq2,0,0);
    E = new StdQuadExp (b1,b2);

    //----------------------------------------------
    // Define solution to be projected
    
    E->GetZW(0,z1,w);
    E->GetZW(1,z2,w);
    
    for(i = 0; i < nq1; ++i)
    {
      for(j = 0; j < nq2; ++j)
      {
    sol[i*nq1 +j]  = Quad_sol(z1[i],z2[j],order1,order2,btype1,btype2);
    sol[i*nq1 +j] *= sol[i*nq1 +j];
      }
    }
    //---------------------------------------------
  }
    break;
  }
  

  //---------------------------------------------
  // Project onto Expansion 
  E->FwdTrans(sol);
  //---------------------------------------------

  // --------------------------------------------
  // Dump coefficients to file 
  ofstream out("AliasingModes.m");
  E->WriteCoeffsToFile(out);
  //--------------------------------------------

  //-------------------------------------------
  // Backward Transform Solution to get projected values
  E->BwdTrans(E->GetPhys());
  //-------------------------------------------  

  // --------------------------------------------
  // Dump physical points to file 
  ofstream out1("AliasingModes.dat");
  E->WriteToFile(out1);
  //--------------------------------------------

  //--------------------------------------------
  // Calculate L_inf error 
  cout << "L infinity error: " << E->Linf(sol) << endl;
  cout << "L 2 error:        " << E->L2  (sol) << endl;
  //--------------------------------------------

  delete[] sol;
}

double Tri_sol(double x, double y, int order1, int order2)
{
  int    l,k;
  double sol = 0;
  
  for(k = 0; k < order1; ++k)
  {
    for(l = 0; l < order2-k; ++l)
    {
      sol += pow(x,k)*pow(y,l);
    }
  }
  
  return sol;
}

double Quad_sol(double x, double y, int order1, int order2, BasisType btype1,
        BasisType btype2)
{

  int k,l;
  double sol = 0;

  if(btype1 != eFourier)
  {
    if(btype2 != eFourier)
    {
      for(k = 0; k < order1; ++k)
      {
    for(l = 0; l < order2; ++l)
    {
      sol += pow(x,k)*pow(y,l);
    }
      }
    }
    else
    {
      for(k = 0; k < order1; ++k)
      {
    for(l = 0; l < order2/2; ++l)
      {
      sol += pow(x,k)*sin(M_PI*l*y) + pow(x,k)*cos(M_PI*l*y);
      }
      }
    }
  }
  else
  {
    if(btype2 != eFourier){
      for(k = 0; k < order1/2; ++k)
      {
    for(l = 0; l < order2; ++l)
    {
      sol += sin(M_PI*k*x)*pow(y,l) + cos(M_PI*k*x)*pow(y,l);
    }
      }
    }
    else
    {
      for(k = 0; k < order1/2; ++k)
      {
    for(l = 0; l < order2/2; ++l)
    {
      sol += sin(M_PI*k*x)*sin(M_PI*l*y)
        + sin(M_PI*k*x)*cos(M_PI*l*y)
        + cos(M_PI*k*x)*sin(M_PI*l*y)
        + cos(M_PI*k*x)*cos(M_PI*l*y);
    }
      }
    }
  }
  
  return sol;
}
