#include <cstdio>
#include <cstdlib>

//#include <MultiRegions/ContExpList1D.h>
#include <MultiRegions/ContField1D.h>

using namespace Nektar;

// This routine projects a which has energy in all mdoes of the
// expansions and report an error.

int main(int argc, char *argv[])
{
    //    MultiRegions::ContExpList1DSharedPtr Exp,Sol;
    MultiRegions::ContField1DSharedPtr Exp,Sol;
    int     i,j,k;
    int     order, nq;
    int     coordim;
    char    *infile;
    LibUtilities::PointsType Qtype;
    LibUtilities::BasisType  btype;  
    Array<OneD,NekDouble>  sol,fce; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 
    NekDouble  lambda = 1.0;

    if(argc != 5)
    {
        fprintf(stderr,"Usage: ProjectCont1D Type order nq  mesh \n");
        
        fprintf(stderr,"Where type is an integer value which "
                "dictates the basis as:\n");
        fprintf(stderr,"\t Modified_A = 4\n");
        fprintf(stderr,"\t GLL Lagrange   = 8\n");
        
        fprintf(stderr,"Note type = 1,2,4,5 are for higher dimensional basis\n");
        
        exit(1);
    }
    
    btype =   (LibUtilities::BasisType) atoi(argv[1]);
    
    // Check to see that only continuous 1D Expansions are used
    if((btype != LibUtilities::eModified_A)&&(btype != LibUtilities::eGLL_Lagrange))
    {
        NEKERROR(ErrorUtil::efatal,
                         "This basis is only for 1D Modified_A or GLL_Lagrange expansions");
    }
    
    // Do not use Fourier expansion
    if(btype == LibUtilities::eFourier)
    {
        NEKERROR(ErrorUtil::efatal,
                         "Demo not set up for Fourier Expanison");
    }
    
    order  =   atoi(argv[2]);
    nq     =   atoi(argv[3]);
    infile =   argv[4];
    
    Qtype = LibUtilities::eGaussLobattoLegendre; 
    
    // read in mesh
    string in(infile);
    SpatialDomains::MeshGraph1D graph1D; 
    graph1D.Read(in);

    SpatialDomains::BoundaryConditions bcs(&graph1D); 
    bcs.Read(in);

    
    // Define Expansion
    const LibUtilities::PointsKey Pkey(nq,Qtype);
    const LibUtilities::BasisKey Bkey(btype,order,Pkey);

    Exp.reset(new MultiRegions::ContField1D(Bkey,graph1D.GetComposite(0),bcs));
    
    //Exp = MemoryManager<MultiRegions::ContExpList1D>::AllocateSharedPtr(Bkey,graph1D.GetComposite(0),bcs);

    //----------------------------------------------
    // Define solution to be projected 
    coordim = Exp->GetCoordim(0);
    nq      = Exp->GetPointsTot();
    
    // define coordinates and solution
    sol = Array<OneD,NekDouble>(nq);
    fce = Array<OneD,NekDouble>(nq);
    
    xc0 = Array<OneD,NekDouble>(nq);
    xc1 = Array<OneD,NekDouble>(nq);
    xc2 = Array<OneD,NekDouble>(nq);
    
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
    
#if 0 // constant solution 
    for(i = 0; i < nq; ++i)
    {
        sol[i] = 2.0; 
        
        fce[i] = -2.0; 
    }
    Exp->SetBoundaryCondition(0,sol[0]);
    Exp->SetBoundaryCondition(1,sol[nq-1]);
#else
#if 1
#if 1 // non-zero dirichlet solution 
    for(i = 0; i < nq; ++i)
    {
        sol[i] = cos(M_PI*xc0[i]);

        fce[i] = -(1.0+M_PI*M_PI)*cos(M_PI*xc0[i]);
    }
    Exp->SetBoundaryCondition(0,sol[0]);
    Exp->SetBoundaryCondition(1,sol[nq-1]);
#else // zero Dirichlet solution 
    for(i = 0; i < nq; ++i)
    {
        sol[i] = sin(M_PI*(xc0[i] - xc0[0])/(xc0[nq-1]-xc0[0]));

        fce[i] = -(1.0+pow(M_PI/(xc0[nq-1]-xc0[0]),2))*sin(M_PI*(xc0[i] - xc0[0])/(xc0[nq-1]-xc0[0]));
    }
    Exp->SetBoundaryCondition(0,0.0);
    Exp->SetBoundaryCondition(1,0.0);
#endif
#else // zero neumann solution
    for(i = 0; i < nq; ++i)
    {
        sol[i] = cos(M_PI*(xc0[i] - xc0[0])/(xc0[nq-1]-xc0[0]));

        fce[i] = -(1.0+pow(M_PI/(xc0[nq-1]-xc0[0]),2))*cos(M_PI*(xc0[i] - xc0[0])/(xc0[nq-1]-xc0[0]));
    }
    Exp->SetBoundaryCondition(0,0.0);
    Exp->SetBoundaryCondition(1,0.0);
#endif
#endif

    //----------------------------------------------
    // Setup Temporary expansion and put in solution
    //Sol = MemoryManager<MultiRegions::ContExpList1D>::AllocateSharedPtr(*Exp);
    Sol = MemoryManager<MultiRegions::ContField1D>::AllocateSharedPtr(*Exp);
    Sol->SetPhys(fce);
    //----------------------------------------------
  
    //---------------------------------------------
    // Helmholtz solution taking physical forcing 
    Exp->HelmSolve(*Sol, lambda);
    //---------------------------------------------
    
    //-------------------------------------------
    // Backward Transform Solution to get projected values
    Exp->BwdTrans(*Exp);
    //-------------------------------------------  
    
    //--------------------------------------------
    // Write solution 
    ofstream outfile("HelmholtzFile1D.dat");
    Exp->WriteToFile(outfile);
    //-------------------------------------------
    
    //--------------------------------------------
    // Calculate L_inf error 
    Sol->SetPhys(sol);
    cout << "L infinity error: " << Exp->Linf(*Sol) << endl;
    cout << "L 2 error:        " << Exp->L2  (*Sol) << endl;
    //--------------------------------------------
    
    return 0;
}

