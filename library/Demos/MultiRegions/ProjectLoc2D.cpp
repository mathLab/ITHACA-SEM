#include <cstdio>
#include <cstdlib>

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList2D.h>

using namespace Nektar;
using namespace StdRegions;
using namespace SpatialDomains; 
using namespace MultiRegions;
using namespace std;

// compile using Builds/Demos/StdRegions -> make DEBUG=1 ProjectLoc2D

// This routine projects a polynomial which has energy in all mdoes of
// the expansions and report an error.

int main(int argc, char *argv[]){
    ExpList2D  *Exp;
    int        i,j,k;
    int        nq;
    double     *sol;
    int        coordim;
    double     **xc;
    int        Tritype, Quadtype, order, Tri_input;
    NodalBasisType Tri_Nb;
    BasisType  Tri_btype1, Tri_btype2, Quad_btype;  
    
    if(argc != 5)
    {
	fprintf(stderr,"Usage: ProjectLoc2D  Tri_Type  "  
		"Quad_Type order mesh.xml\n");
	
	fprintf(stderr,"Where Type is an integer value which "
		"dictates the basis type as:\n");
	
	fprintf(stderr,"Triangle options:\n");
	fprintf(stderr,"\t Ortho_A, Ortho_B       = 1\n");
	fprintf(stderr,"\t Modified_A, Modified_B = 2\n");
	fprintf(stderr,"\t Nodal Elec             = 3\n");
	fprintf(stderr,"\t Nodal Fekete           = 4\n");

	fprintf(stderr,"Quadrilateral options:\n");
	fprintf(stderr,"\t Ortho_A, Ortho_A       = 5\n");
	fprintf(stderr,"\t Modified_A, Modified_A = 6\n");
	fprintf(stderr,"\t Lagrange, Lagrange     = 7\n");
	fprintf(stderr,"\t Chebyshev, Chebychev   = 8\n");
	
	exit(1);
    }
    
    if((Tri_input = atoi(argv[1]))  < 3)
    {
	Tritype  = (StdRegions::BasisType) atoi(argv[1]);
	Tri_Nb   = (StdRegions::NodalBasisType) NULL;
    }
    else
    {
	Tritype  = (StdRegions::BasisType) 1;
	if(Tri_input == 3)
	{
	    Tri_Nb   = (StdRegions::NodalBasisType) StdRegions::eNodalTriElec;
	}
	else
	{
	    Tri_Nb   = (StdRegions::NodalBasisType) StdRegions::eNodalTriFekete;
	}
    }

    Quadtype  = (BasisType) atoi(argv[2]);
    order     = atoi(argv[3]);
    
    if((Tritype < 1)||(Tritype > 4))
    {
      ErrorUtil::Error(ErrorUtil::efatal,__FILE__,__LINE__,
		       "Illegal option for Tri_Type\n");
    }

    if((Quadtype < 5)||(Quadtype > 8))
    {
      ErrorUtil::Error(ErrorUtil::efatal,__FILE__,__LINE__,
		       "Illegal option for Quad_Type\n");
    }
    // read in mesh
    string  in(argv[argc-1]);
    MeshGraph2D graph2D;
    graph2D.Read(in);
    
    switch(Tritype){
    case 1: case 3: case 4:
	Tri_btype1 = eOrtho_A;
	Tri_btype2 = eOrtho_B;

	break;
    case 2:
	Tri_btype1 = eModified_A;
	Tri_btype2 = eModified_B;
	break;
    }
    
    switch(Quadtype){
    case 5:
	Quad_btype = eOrtho_A;
	break;
    case 6:
	Quad_btype = eModified_A;
	break;
    case 7:
	Quad_btype = eGLL_Lagrange;
	break;
    case 8:
	Quad_btype = eChebyshev;
	break;
    }

    const BasisKey T_Ba(Tri_btype1, order,  eLobatto, order+1, 0,0);
    const BasisKey T_Bb(Tri_btype2, order,  eLobatto, order+1, 1,0);
    const BasisKey Q_Ba(Quad_btype ,order,  eLobatto, order+1,0,0);
    
    Exp = new ExpList2D (T_Ba,T_Bb,Q_Ba,Q_Ba,graph2D, Tri_Nb);
    
    //----------------------------------------------
    // Define solution to be projected 
    coordim = Exp->GetCoordim(0);
    nq      = Exp->GetPointsTot();

    // Define coordinates and solution
    sol   = new double  [nq];
    xc    = new double *[coordim];
    xc[0] = new double  [coordim*nq];
    
    for(i = 1; i < coordim; ++i)
    {
	xc[i] = xc[i-1] + nq;
    }
    
    Exp->GetCoords(xc);
    
    for(i = 0; i < nq; ++i)
    {
	sol[i] = 0.0;
	for(j = 0; j < order; ++j)
	{
	    for(k = 0; k < coordim; ++k)
	    {
		sol[i] += pow(xc[k][i],j);
	    }
	}
    }
    
    //---------------------------------------------
    // Project onto Expansion 
    Exp->FwdTrans(sol);
    //---------------------------------------------
    
    //-------------------------------------------
    // Backward Transform Solution to get projected values
    Exp->BwdTrans(Exp->GetPhys());
    //-------------------------------------------  
    
    //--------------------------------------------
    // Write solution 
    ofstream outfile("ProjectLocFile2D.dat");
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
	return 0;
}
