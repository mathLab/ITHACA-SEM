#include <cstdio>
#include <cstdlib>

#include <MultiRegions/MultiRegions.h>
#include <MultiRegions/ExpList1D.h>

using namespace StdRegions;
using namespace SpatialDomains; 
using namespace MultiRegions;
using namespace std;

// compile using Builds/Demos/StdRegions -> make DEBUG=1 ProjectLoc1D

// This routine projects a polynomial which has energy in all mdoes of
// the expansions and report an error.

main(int argc, char *argv[]){
  ExpList1D  *Exp;
  int i,j,k;
  int order, nq;
  double *sol;
  int coordim;
  double  **xc;
  char    *infile;
  PointType Qtype;
  BasisType btype;  

  if(argc != 5){
    fprintf(stderr,"Usage: ProjectLoc1D Type order nq  mesh \n");

    fprintf(stderr,"Where type is an integer value which "
	    "dictates the basis as:\n");
    fprintf(stderr,"\t Ortho_A    = 0\n");
    fprintf(stderr,"\t Modified_A = 3\n");
    fprintf(stderr,"\t Lagrange   = 7\n");
    fprintf(stderr,"\t Legendre   = 8\n"); 
    fprintf(stderr,"\t Chebyshev  = 9\n");
 
    fprintf(stderr,"Note type = 1,2,4,5 are for higher dimensional basis\n");

    exit(1);
  }
  
  btype =   (BasisType) atoi(argv[1]);
  
  // Check to see that only 1D Expansions are used
  if((btype == Ortho_B)||(btype == Ortho_B)||
     (btype == Modified_B)||(btype == Modified_C))
    ErrorUtil::Error(fatal,"ProjectLoc1D",
		     "This basis is for 2 or 3D expansions");

  // Do not use Fourier expansion
  if(btype == Fourier)
    ErrorUtil::Error(fatal,"ProjectLoc1D",
		     "Demo not set up for Fourier Expanison");
  
  order  =   atoi(argv[2]);
  nq     =   atoi(argv[3]);
  infile =   argv[4];

  if(btype != Fourier)
    Qtype = Lobatto; 
  else
    Qtype = FourierEvenSp;
  
  // read in mesh
  ifstream in(infile);
  MeshGraph1D graph1D(in);

  // Define Expansion
  const BasisKey B(btype,order, Qtype, nq,0,0);
  Exp = new ExpList1D (B,graph1D);
  
  //----------------------------------------------
  // Define solution to be projected 
  coordim = Exp->get_coordim(0);
  nq      = Exp->get_nquad();

  // define coordinates and solution
  sol   = new double  [nq];
  xc    = new double *[coordim];
  xc[0] = new double  [coordim*nq];

  for(i = 1; i < coordim; ++i)
    xc[i] = xc[i-1] + nq;

  Exp->GetCoords(xc);
  
  for(i = 0; i < nq; ++i){
    sol[i] = 0.0;
    for(j = 0; j < order; ++j)
      for(k = 0; k < coordim; ++k)
	sol[i] += pow(xc[k][i],j);
  }

  //---------------------------------------------
  // Project onto Expansion 
  Exp->FwdTrans(sol);
  //---------------------------------------------

  //-------------------------------------------
  // Backward Transform Solution to get projected values
  Exp->BwdTrans(Exp->get_phys());
  //-------------------------------------------  

  //--------------------------------------------
  // Write solution 
  ofstream outfile("ProjectLocFile1D.dat");
  Exp->WriteToFile(outfile);
  //-------------------------------------------

  //--------------------------------------------
  // Calculate L_inf error 
  cout << "L infinity error: " << Exp->Linf(sol) << endl;
  cout << "L 2 error:        " << Exp->L2  (sol) << endl;
  //--------------------------------------------

  delete[] sol; 
  delete[] xc[0];
  delete[] xc;
}
