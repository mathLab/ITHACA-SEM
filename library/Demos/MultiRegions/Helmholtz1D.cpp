#include <cstdio>
#include <cstdlib>

//#include <MultiRegions/ContExpList1D.h>
#include <MultiRegions/ContField1D.h>

using namespace Nektar;

// This routine projects a which has energy in all mdoes of the
// expansions and report an error.

int main(int argc, char *argv[])
{
    MultiRegions::ContField1DSharedPtr Exp,Fce;
    int     i,j,k;
    int     order, nq;
    int     coordim;
    char    *infile;
    LibUtilities::PointsType Qtype;
    LibUtilities::BasisType  btype;  
    Array<OneD,NekDouble>  sol,fce; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 
    NekDouble  lambda;

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

    lambda = bcs.GetParameter("Lambda");
    cout << "Solving Helmholtz problem with Lambda = " << lambda << endl;
    
    // Define Expansion
    const LibUtilities::PointsKey Pkey(nq,Qtype);
    const LibUtilities::BasisKey Bkey(btype,order,Pkey);
    SpatialDomains::Composite domain = graph1D.GetDomain();

    Exp = MemoryManager<MultiRegions::ContField1D>::AllocateSharedPtr(Bkey,domain,bcs);

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
    
    // Get Forcing function for variable 0
    SpatialDomains::ConstForcingFunctionShPtr ffunc 
        = bcs.GetForcingFunction(bcs.GetVariable(0));

    for(i = 0; i < nq; ++i)
    {
        fce[i] = ffunc->Evaluate(xc0[i],xc1[i],xc2[i]);
    }

    //----------------------------------------------
    // Setup Temporary expansion and put in solution
    //Sol = MemoryManager<MultiRegions::ContExpList1D>::AllocateSharedPtr(*Exp);
    Fce = MemoryManager<MultiRegions::ContField1D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //----------------------------------------------
  
    //---------------------------------------------
    // Helmholtz solution taking physical forcing 
    Exp->HelmSolve(*Fce, lambda);
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
    
    SpatialDomains::ConstExactSolutionShPtr ex_sol =
        bcs.GetExactSolution(bcs.GetVariable(0));

    if(ex_sol)
    {
        for(i = 0; i < nq; ++i)
        {
            fce[i] = ex_sol->Evaluate(xc0[i],xc1[i],xc2[i]);
        }

        //--------------------------------------------
        // Calculate L_inf error 
        Fce->SetPhys(fce);
        cout << "L infinity error: " << Exp->Linf(*Fce) << endl;
        cout << "L 2 error:        " << Exp->L2  (*Fce) << endl;
        //--------------------------------------------
        
    }
    return 0;
}

