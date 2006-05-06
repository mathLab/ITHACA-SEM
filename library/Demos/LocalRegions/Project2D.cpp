#include <cstdio>
#include <cstdlib>

#include <LocalRegions/LocalRegions.hpp>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/TriExp.h>

#include <SpatialDomains/QuadGeom.h>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/MeshComponents.h>

using namespace Nektar;
using namespace StdRegions; 
using namespace LocalRegions; 
using namespace SpatialDomains; 
using namespace std;

static double quadsol(double x1,double x2,int order1, int order2, 
		      BasisType btype1, BasisType btype2);

static double trisol(double x1,double x2,int order1, int order2);

// compile using Builds/Demos/StdRegions -> make DEBUG=1 Project2D

// This routine projects a polynomial or trigonmetric functions which 
// has energy in all mdoes of the expansions and reports and error

int main(int argc, char *argv[])
{
  int           i,j,k,l;
  const         double *z1,*z2,*w;
  int           order1,order2, nq1,nq2;
  PointsType    Qtype1,Qtype2;
  BasisType     btype1,btype2;
  ShapeType     regionshape;
  StdExpansion2D *E;
  double        *sol, coords[8];
  EdgeOrientation edgeDir = eForwards; 
  
  if((argc != 16)&&(argc != 14))
  {
    fprintf(stderr,"Usage: Project2D RegionShape Type1 Type2 order1 "	   
	    "order2  nq1 nq2 x1, y1, x2, y2,.... \n");
    
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

    fprintf(stderr,"The last series of values are the coordinates\n");
    exit(1);
  }
  
  regionshape = (ShapeType) atoi(argv[1]);
  
  // Check to see if 2D region 
  if((regionshape != eTriangle)&&(regionshape != eQuadrilateral))
  {
    ErrorUtil::Error(ErrorUtil::efatal,__FILE__,__LINE__,
		     "This shape is not a 2D region");
  }

  btype1 = (BasisType) atoi(argv[2]);
  btype2 = (BasisType) atoi(argv[3]);

  // Check to see that correct Expansions are used
  switch(regionshape)
  {
  case eTriangle:
    if((btype1 == eOrtho_B)||(btype1 == eModified_B))
    {
      ErrorUtil::Error(ErrorUtil::efatal,__FILE__,__LINE__,
		       "Basis 1 cannot be of type Ortho_B or Modified_B");
    }
    
    if((btype2 != eOrtho_B)&&(btype2 != eModified_B))
    {
      ErrorUtil::Error(ErrorUtil::efatal,__FILE__,__LINE__,
		       "Basis 2 must be of type Ortho_B or Modified_B");
    }
    break;
  case eQuadrilateral:

    if((btype1 == eOrtho_B)||(btype1 == eOrtho_B)||
       (btype1 == eModified_B)||(btype1 == eModified_C))
    {
      ErrorUtil::Error(ErrorUtil::efatal,__FILE__,__LINE__,
		     "Basis 1 is for 2 or 3D expansions");
    }

    if((btype2 == eOrtho_B)||(btype2 == eOrtho_B)||
       (btype2 == eModified_B)||(btype2 == eModified_C))
    {
      ErrorUtil::Error(ErrorUtil::efatal,__FILE__,__LINE__,
		       "Basis 2 is for 2 or 3D expansions");
    }
    break;
  }
  
  order1 =   atoi(argv[4]);
  order2 =   atoi(argv[5]);
  nq1    =   atoi(argv[6]);
  nq2    =   atoi(argv[7]);

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
  case eTriangle:
  {
    coords[0]    =   atof(argv[8]);
    coords[1]    =   atof(argv[9]);
    coords[2]    =   atof(argv[10]);
    coords[3]    =   atof(argv[11]);
    coords[4]    =   atof(argv[12]);
    coords[5]    =   atof(argv[13]);

    // Set up coordinates
    VertexComponent *vert[3];
    vert[0] = new VertexComponent(2,0,coords);
    vert[1] = new VertexComponent(2,1,coords+2);
    vert[2] = new VertexComponent(2,2,coords+4);

    // Set up Edges
    EdgeComponent *edge[3];
    edge[0] = new EdgeComponent(2);
    edge[1] = new EdgeComponent(2);
    edge[2] = new EdgeComponent(2);

    EdgeOrientation eorient[3];
    
    eorient[0] = edgeDir; 
    eorient[1] = edgeDir; 
    eorient[2] = edgeDir; 
      
    TriGeom *geom = new TriGeom(vert,edge,eorient);
    geom->SetOwnData();

    const BasisKey b1(btype1,order1,Qtype1,nq1,0,0);
    const BasisKey b2(btype2,order2,Qtype2,nq2,1,0);

    E = new TriExp(b1,b2,geom);
    E->SetGeoFac(E->GenGeoFac());
    
    E->GetZW(0,z1,w);
    E->GetZW(1,z2,w);

    //----------------------------------------------
    // Define solution to be projected

    double *xc[2];
    xc[0] = new double [2*nq1*nq2];  
    xc[1] = xc[0] + nq1*nq2;
    E->GetCoords(xc);
    
    for(i = 0; i < nq1; ++i)
    {
      for(j = 0; j < nq2; ++j)
      {
	sol[i+nq1*j] = trisol(xc[0][j*nq1+i],xc[1][j*nq1+i],
			       order1,order2);
      }
    }
    delete[] xc[0];
    //----------------------------------------------    
  }
  break;
  case eQuadrilateral:
  {
    // Gather coordinates
    coords[0]    =   atof(argv[8]);
    coords[1]    =   atof(argv[9]);
    coords[2]    =   atof(argv[10]);
    coords[3]    =   atof(argv[11]);
    coords[4]    =   atof(argv[12]);
    coords[5]    =   atof(argv[13]);
    coords[6]    =   atof(argv[14]);
    coords[7]    =   atof(argv[15]);

    // Set up coordinates
    VertexComponent *vert[4];
    vert[0] = new VertexComponent(2,0,coords);
    vert[1] = new VertexComponent(2,1,coords+2);
    vert[2] = new VertexComponent(2,2,coords+4);
    vert[3] = new VertexComponent(2,3,coords+6);

    // Set up Edges
    EdgeComponent *edge[4];
    edge[0] = new EdgeComponent(2);
    edge[1] = new EdgeComponent(2);
    edge[2] = new EdgeComponent(2);
    edge[3] = new EdgeComponent(2);

    EdgeOrientation eorient[4];
    
    eorient[0] = edgeDir; 
    eorient[1] = edgeDir; 
    eorient[2] = edgeDir; 
    eorient[3] = edgeDir; 
      
    QuadGeom *geom = new QuadGeom(vert,edge,eorient);
    geom->SetOwnData();

    // Define basis
    const BasisKey b1(btype1,order1,Qtype1,nq1,0,0);
    const BasisKey b2(btype2,order2,Qtype2,nq2,0,0);

    E = new QuadExp(b1,b2,geom);
    E->SetGeoFac(E->GenGeoFac());

    //----------------------------------------------
    // Define solution to be projected
    
    E->GetZW(0,z1,w);
    E->GetZW(1,z2,w);

    double *xc[2];
    xc[0] = new double [2*nq1*nq2];  
    xc[1] = xc[0] + nq1*nq2;
    E->GetCoords(xc);
    
    for(i = 0; i < nq1; ++i)
    {
      for(j = 0; j < nq2; ++j)
      {
	sol[i+nq1*j] = quadsol(xc[0][j*nq1+i],xc[1][j*nq1+i],
			       order1,order2,btype1,btype2);
      }
    }

    delete[] xc[0];
    //---------------------------------------------
  }

  break;
  }

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
  FILE *outfile = fopen("ProjectFile2D.dat","w");
  E->WriteToFile(outfile);
  fclose(outfile);
  //-------------------------------------------

  //--------------------------------------------
  // Calculate L_inf error 
  cout << "L infinity error: " << E->Linf(sol) << endl;
  cout << "L 2 error:        " << E->L2  (sol) << endl;
  //--------------------------------------------

  //-------------------------------------------
  // Evaulate solution at x = y =0  and print error
  
  double x[2];
  x[0] = (coords[0] + coords[2])*0.5;
  x[1] = (coords[1] + coords[5])*0.5;
  double nsol = E->Evaluate(x);
  
  if(regionshape == eTriangle)
  {
    sol[0] = trisol(x[0],x[1],order1,order2);
  }
  else
  {
    sol[0] = quadsol(x[0],x[1],order1,order2,btype1,btype2);
  }
  
  cout << "error at (x,y) = " << x[0] << "," << x[1] << " : " 
       << nsol - sol[0] << endl;

  //-------------------------------------------

  delete[] sol;
  return 0;
}

static double quadsol(double x1,double x2,int order1, int order2, 
		      BasisType btype1, BasisType btype2)
{
  int k,l;
  double sol = 0.0; 
  
  if(btype1 != eFourier)
  {
    if(btype2 != eFourier)
    {
      for(k = 0; k < order1; ++k)
      {
	for(l = 0; l < order2; ++l)
	{
	  sol += pow(x1,k)*pow(x2,l);
	}
      }
    }
    else
    {
      for(k = 0; k < order1; ++k)
      {
	for(l = 0; l < order2/2; ++l)
	{
	  sol += pow(x1,k)*sin(M_PI*l*x2) 
		  + pow(x1,k)*cos(M_PI*l*x2);
	}
      }
    }
  }
  else
  {
    if(btype2 != eFourier)
    {
      for(k = 0; k < order1/2; ++k)
      {
	for(l = 0; l < order2; ++l)
	{
	  sol += sin(M_PI*k*x1)*pow(x2,l) 
	    + cos(M_PI*k*x1)*pow(x2,l);
	}
      }
    }
    else
    {
      for(k = 0; k < order1/2; ++k)
      {
	for(l = 0; l < order2/2; ++l)
	{
	  sol += sin(M_PI*k*x1)*sin(M_PI*l*x2)
	    + sin(M_PI*k*x1)*cos(M_PI*l*x2)
	    + cos(M_PI*k*x1)*sin(M_PI*l*x2)
	    + cos(M_PI*k*x1)*cos(M_PI*l*x2);
	}
      }
    }
  }

  return sol;
  
}

static double trisol(double x1,double x2,int order1, int order2)
{
  int k,l;
  double sol = 0.0; 

  // activate all modes
  for(k = 0; k < order1; ++k)
  {
    for(l = 0; l < order2-k; ++l)
    {
      sol += pow(x1,k)*pow(x2,l);
    }
  }

  return sol;
}
