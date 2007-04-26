#include <cstdio>
#include <cstdlib>

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList1D.h>

using namespace Nektar;

// compile using Builds/Demos/StdRegions -> make DEBUG=1 ProjectLoc1D

// This routine projects a polynomial which has energy in all mdoes of
// the expansions and report an error.

int main(int argc, char *argv[])
{
    MultiRegions::ExpList1D  *Exp;
    int i,j,k;
    int     order, nq;
    int     coordim;
    char    *infile;
    LibUtilities::PointsType Qtype;
    LibUtilities::BasisType  btype;  
    NekDoubleSharedArray sol; 
    NekDoubleSharedArray xc0,xc1,xc2; 
    
    if(argc != 5)
    {
	fprintf(stderr,"Usage: ProjectLoc1D Type order nq  mesh \n");
	
	fprintf(stderr,"Where type is an integer value which "
		"dictates the basis as:\n");
	fprintf(stderr,"\t Ortho_A    = 1\n");
	fprintf(stderr,"\t Modified_A = 4\n");
	fprintf(stderr,"\t Lagrange   = 8\n");
	fprintf(stderr,"\t Legendre   = 9\n"); 
	fprintf(stderr,"\t Chebyshev  = 10\n");
	
	fprintf(stderr,"Note type = 1,2,4,5 are for higher dimensional basis\n");
	
	exit(1);
    }
    
    btype =   (LibUtilities::BasisType) atoi(argv[1]);
    
    // Check to see that only 1D Expansions are used
    if((btype == LibUtilities::eOrtho_B)||(btype == LibUtilities::eOrtho_B)||
       (btype == LibUtilities::eModified_B)||(btype == LibUtilities::eModified_C))
	ErrorUtil::Error(ErrorUtil::efatal,__FILE__,__LINE__,
			 "This basis is for 2 or 3D expansions");
    
    // Do not use Fourier expansion
    if(btype == LibUtilities::eFourier)
    {
	ErrorUtil::Error(ErrorUtil::efatal,__FILE__,__LINE__,
			 "Demo not set up for Fourier Expanison");
    }
    
    order  =   atoi(argv[2]);
    nq     =   atoi(argv[3]);
    infile =   argv[4];
    
    if(btype != LibUtilities::eFourier)
    {
	Qtype = LibUtilities::eGaussLobattoLegendre; 
    }
    else
    {
	Qtype = LibUtilities::eFourierEvenlySpaced;
    }
  
    // read in mesh
    string  in(infile);
    SpatialDomains::MeshGraph1D graph1D;
    graph1D.Read(in);
    
    // Define Expansion
    const LibUtilities::PointsKey Pkey(nq,Qtype);
    const LibUtilities::BasisKey Bkey(btype,order,Pkey);
    Exp = new MultiRegions::ExpList1D (Bkey,graph1D);
    
    //----------------------------------------------
    // Define solution to be projected 
    coordim = Exp->GetCoordim(0);
    nq      = Exp->GetPointsTot();
    
    // define coordinates and solution
    sol = GetDoubleTmpSpace(nq);

    xc0 = GetDoubleTmpSpace(nq);
    xc1 = GetDoubleTmpSpace(nq);
    xc2 = GetDoubleTmpSpace(nq);

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
    // Project onto Expansion 
    Exp->FwdTrans(sol,*Exp);
    //---------------------------------------------
    
    //-------------------------------------------
    // Backward Transform Solution to get projected values
    Exp->BwdTrans(*Exp);
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
}
