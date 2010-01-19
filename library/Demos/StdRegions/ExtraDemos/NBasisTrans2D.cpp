#include <cstdio>
#include <cstdlib>
#include "StdRegions/StdExpansion2D.h"
#include "StdRegions/StdNodalTriExp.h"
#include "StdRegions/NodalBasisManager.h"

#include "StdRegions/StdRegions.hpp"

using namespace Nektar;
using namespace StdRegions; 
using namespace std;

// compile using Builds/Demos/StdRegions -> make DEBUG=1 NBasisTrans2D

// This routine generates a basis transformation from a nodal expansion to 
// the orthogonal expansion

int main(int argc, char *argv[])
{

  int           i,j,k,l;
  const         double *z1,*z2,*w;
  int           order, nq1,nq2;
  StdExpansion2D  *E;
  double        *sol,*coeffs;

  if(argc != 3)
  {
    fprintf(stderr,"Usage: NBasisTrans2D PointDis order \n");

    fprintf(stderr,"Where PointDist is an integer value which "
        "dictates the nodal distrbituion:\n");
    fprintf(stderr,"\t Electrostatic points  = 0\n");
    fprintf(stderr,"\t Fekete points         = 1\n");
    
    exit(1);
  }
  
  order =   atoi(argv[argc-1]);
  
  // generate Nodal triangle expansion
  NodalBasisType btype = (NodalBasisType) atoi(argv[argc-2]);

  const BasisKey b1(eOrtho_A,order,eLobatto,order+1,0,0);
  const BasisKey b2(eOrtho_B,order,eRadauM, order,1,0);
  // const BasisKey b1(Modified_A,order,Lobatto,order+1,0,0);
  // const BasisKey b2(Modified_B,order,RadauM, order,1,0);

  E = new StdNodalTriExp (b1,b2,btype);

  nq1 = E->GetPointsOrder(0);
  nq2 = E->GetPointsOrder(1);
  sol = new double [nq1*nq2];
  
  coeffs = E->GetCoeffs();
  E->GetZW(0,z1,w);
  E->GetZW(1,z2,w);
  
  //----------------------------------------------
  // Define solution at  Nodal points and quadrature points 

  double x,y;
  for(i = 0; i < nq1; ++i)
  {
    for(j = 0; j < nq2; ++j)
    {
      x = (1+z1[i])*(1-z2[j])/2-1.0;
      y = z2[j];
      sol[i+nq1*j] = 0.0;
      // activate all modes
      for(k = 0; k < order; ++k)
      {
    for(l = 0; l < order-k; ++l)
    {
      sol[i+nq1*j] += pow(x,k)*pow(y,l);
    }
      }
    }
  }

  const double *xp, *yp;
  int npts = E->GetNodalPoints(xp,yp);

  for(i = 0; i < npts; ++i)
  {
    coeffs[i] = 0.0;
    // activate all modes
    for(k = 0; k < order; ++k)
    {
    for(l = 0; l < order-k; ++l)
    {
        coeffs[i] += pow(xp[i],k)*pow(yp[i],l);
    }
    }
  }

  
  //-------------------------------------------
  // Backward Transform Solution to get projected values
  E->BwdTrans(E->GetPhys());
  //-------------------------------------------  

  //--------------------------------------------
  // Calculate L_inf error 
  cout << "L infinity error: " << E->Linf(sol) << endl;
  cout << "L 2 error:        " << E->L2  (sol) << endl;
  //--------------------------------------------

  delete[] sol;
  return 0;
}
