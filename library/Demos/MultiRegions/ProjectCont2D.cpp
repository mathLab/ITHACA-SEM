#include <cstdio>
#include <cstdlib>

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ContExpList2D.h>

using namespace Nektar;
using namespace StdRegions;
using namespace SpatialDomains; 
using namespace MultiRegions;
using namespace std;

// compile using Builds/Demos/StdRegions -> make DEBUG=1 ProjectLoc2D

// This routine projects a polynomial which has energy in all mdoes of
// the expansions and report an error.

int main(int argc, char *argv[]){
    ContExpList2D  *Exp;
    int        i,j,k;
    int        nq;
    double     *sol;
    int        coordim;
    double     **xc;
    int        Tritype, Quadtype, order;
    NodalBasisType Tri_Nb;
    BasisType  Tri_btype1, Tri_btype2, Quad_btype;  
    
    if(argc != 5)
    {
	fprintf(stderr,"Usage: ProjectLoc2D  Tri_Type  "  
		"Quad_Type order mesh.xml\n");
	
	fprintf(stderr,"Where Type is an integer value which "
		"dictates the basis type as:\n");
	
	fprintf(stderr,"Triangle options:\n");
	fprintf(stderr,"\t Modified_A, Modified_B = 1\n");
	fprintf(stderr,"\t Nodal Elec             = 2\n");
	fprintf(stderr,"\t Nodal Fekete           = 3\n");

	fprintf(stderr,"Quadrilateral options:\n");
	fprintf(stderr,"\t Modified_A, Modified_A = 4\n");
	fprintf(stderr,"\t Lagrange, Lagrange     = 5\n");
	
	exit(1);
    }
    
    Tritype   = (BasisType) atoi(argv[1]);
    Quadtype  = (BasisType) atoi(argv[2]);
    order     = atoi(argv[3]);
    
    if((Tritype < 1)||(Tritype > 3))
    {
	ErrorUtil::Error(ErrorUtil::efatal,__FILE__,__LINE__,
			 "Illegal option for Tri_Type\n");
    }
    
    if((Quadtype < 4)||(Quadtype > 5))
    {
	ErrorUtil::Error(ErrorUtil::efatal,__FILE__,__LINE__,
			 "Illegal option for Quad_Type\n");
    }

    // read in mesh
    string  in(argv[argc-1]);
    MeshGraph2D graph2D;
    graph2D.Read(in);
    
    switch(Tritype){
    case 1: 
	Tri_btype1 = eModified_A;
	Tri_btype2 = eModified_B;
	Tri_Nb     = (StdRegions::NodalBasisType) NULL;
	break;
    case 2:
	Tri_btype1 = (StdRegions::BasisType) StdRegions::eOrtho_A;
	Tri_btype2 = (StdRegions::BasisType) StdRegions::eOrtho_B;
	Tri_Nb     = (StdRegions::NodalBasisType) StdRegions::eNodalTriElec;
	break;
    case 3:
	Tri_btype1 = (StdRegions::BasisType) StdRegions::eOrtho_A;
	Tri_btype2 = (StdRegions::BasisType) StdRegions::eOrtho_B;
	Tri_Nb     = (StdRegions::NodalBasisType) StdRegions::eNodalTriFekete;
	break;
    }
    
    switch(Quadtype){
    case 4:
	Quad_btype = eModified_A;
	break;
    case 5:
	Quad_btype = eGLL_Lagrange;
	break;
    }

    const BasisKey T_Ba(Tri_btype1, order, eLobatto, order+1, 0,0);
    const BasisKey T_Bb(Tri_btype2, order, eLobatto, order+1, 1,0);
    const BasisKey Q_Ba(Quad_btype ,order, eLobatto, order+1,0,0);
    
    Exp = new ContExpList2D (T_Ba,T_Bb,Tri_Nb,Q_Ba,Q_Ba,graph2D);
    
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
    ofstream outfile("ProjectContFile2D.dat");
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
