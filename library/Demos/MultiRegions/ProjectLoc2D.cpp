#include <cstdio>
#include <cstdlib>

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList2D.h>

using namespace Nektar;

// This routine projects a polynomial which has energy in all mdoes of
// the expansions and report an error.

int main(int argc, char *argv[])
{
    MultiRegions::ExpList2DSharedPtr Exp,Sol;
    int        i,j,k;
    int        order, nq;
    int        coordim;
    char    *infile;
    int        Tritype, Quadtype, Tri_input;
    LibUtilities::PointsType Tri_Nb, Qtype;
    LibUtilities::BasisType  Tri_btype1, Tri_btype2, Quad_btype;
    Array<OneD, NekDouble> sol; 
    Array<OneD, NekDouble> xc0,xc1,xc2; 
    
    if(argc != 6)
    {
	fprintf(stderr,"Usage: ProjectLoc2D  Tri_Type  "  
		"Quad_Type order nq mesh\n");
	
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
	Tritype  = (LibUtilities::BasisType) atoi(argv[1]);
	Tri_Nb   = LibUtilities::SIZE_PointsType;
    }
    else
    {
	Tritype  = (LibUtilities::BasisType) 1;
	if(Tri_input == 3)
	{
	    Tri_Nb   = (LibUtilities::PointsType) LibUtilities::eNodalTriElec;
	}
	else
	{
	    Tri_Nb   = (LibUtilities::PointsType) LibUtilities::eNodalTriFekete;
	}
    }

    Quadtype  = (LibUtilities::BasisType) atoi(argv[2]);
    order     = atoi(argv[3]);
    nq        = atoi(argv[4]);
    infile    = argv[5];

    Qtype = LibUtilities::eGaussLobattoLegendre; 
    
    if((Tritype < 1)||(Tritype > 4))
    {
      NEKERROR(ErrorUtil::efatal,
		       "Illegal option for Tri_Type\n");
    }

    if((Quadtype < 5)||(Quadtype > 8))
    {
      NEKERROR(ErrorUtil::efatal,
		       "Illegal option for Quad_Type\n");
    }
    // read in mesh
    string  in(infile);
    SpatialDomains::MeshGraph2D graph2D;
    graph2D.Read(in);
    
    switch(Tritype){
    case 1: case 3: case 4:
	Tri_btype1 = LibUtilities::eOrtho_A;
	Tri_btype2 = LibUtilities::eOrtho_B;

	break;
    case 2:
	Tri_btype1 = LibUtilities::eModified_A;
	Tri_btype2 = LibUtilities::eModified_B;
	break;
    }
    
    switch(Quadtype){
    case 5:
	Quad_btype = LibUtilities::eOrtho_A;
	break;
    case 6:
	Quad_btype = LibUtilities::eModified_A;
	break;
    case 7:
	Quad_btype = LibUtilities::eGLL_Lagrange;
	break;
    case 8:
	Quad_btype = LibUtilities::eChebyshev;
	break;
    }

    // Define Expansion
    const LibUtilities::PointsKey PKey(nq,Qtype);
    const LibUtilities::BasisKey T_Ba(Tri_btype1,order,PKey);
    const LibUtilities::BasisKey T_Bb(Tri_btype2,order,PKey);
    const LibUtilities::BasisKey Q_Ba(Quad_btype,order,PKey);

    Exp = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(T_Ba,T_Bb,Q_Ba,Q_Ba,graph2D, Tri_Nb);    
    
    //----------------------------------------------
    // Define solution to be projected 
    coordim = Exp->GetCoordim(0);
    nq      = Exp->GetPointsTot();

    // Define coordinates and solution
    sol = Array<OneD, NekDouble>(nq);

    xc0 = Array<OneD, NekDouble>(nq);
    xc1 = Array<OneD, NekDouble>(nq);
    xc2 = Array<OneD, NekDouble>(nq);
    
    switch(coordim)
    {
    case 1:
        Exp->GetCoords(xc0);
        Vmath::Zero(nq,&xc1[0],1);
        Vmath::Zero(nq,&xc2[0],1);
        break;
    case 2:
        Exp->GetCoords(xc0,xc1);
        Vmath::Zero(nq,&xc2[0],1);
        break;
    case 3:
        Exp->GetCoords(xc0,xc1,xc2);
        break;
    }
    
    for(i = 0; i < nq; ++i)
    {
	sol[i] = 0.0;
	for(j = 0; j < order; ++j)
	{
            sol[i] += pow(xc0[i],j);
            sol[i] += pow(xc1[i],j);
            sol[i] += pow(xc2[i],j);
	}
    }
   
    
    //---------------------------------------------
    // Set up ExpList1D containing the solution 
    Sol = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(*Exp);
    Sol->SetPhys(sol);
    //---------------------------------------------
 
    //---------------------------------------------
    // Project onto Expansion 
    Exp->FwdTrans(*Sol);
    //---------------------------------------------
    
    //-------------------------------------------
    // Backward Transform Solution to get projected values
    Exp->BwdTrans(*Exp);
    //-------------------------------------------  
    
    //--------------------------------------------
    // Write solution 
    ofstream outfile("ProjectLocFile2D.dat");
    Exp->WriteToFile(outfile);
    //-------------------------------------------
    
    //--------------------------------------------
    // Calculate L_inf error 
    cout << "L infinity error: " << Exp->Linf(*Sol) << endl;
    cout << "L 2 error:        " << Exp->L2  (*Sol) << endl;
    //--------------------------------------------

    return 0;
}
