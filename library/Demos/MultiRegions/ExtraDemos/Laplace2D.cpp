#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ContField2D.h>

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    MultiRegions::ContField2DSharedPtr Exp,Fce;
    int     i, j, nq,  coordim;
    Array<OneD,NekDouble>  fce; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 
    NekDouble  lambda;

    if(argc != 3)
    {
        fprintf(stderr,"Usage: Laplace2D  meshfile boundaryfile\n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[1]);
    SpatialDomains::MeshGraph2D graph2D; 
    graph2D.ReadGeometry(meshfile);
    graph2D.ReadExpansions(meshfile);
    //----------------------------------------------

    //----------------------------------------------
    // read the problem parameters from input file
    string bcfile(argv[2]);
    SpatialDomains::BoundaryConditions bcs(&graph2D); 
    bcs.Read(bcfile);
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    const SpatialDomains::ExpansionVector &expansions = graph2D.GetExpansions();
    LibUtilities::BasisKey bkey = graph2D.GetBasisKey(expansions[0],0);
    cout << "Solving 2D Laplace:"  << endl; 
    cout << "         Expansion  : " << SpatialDomains::kExpansionTypeStr[expansions[0]->m_ExpansionType] << endl;
    cout << "         No. modes  : " << (int) expansions[0]->m_NumModesEqn.Evaluate() << endl;
    cout << endl;
    //----------------------------------------------
   
    //----------------------------------------------
    // Define Expansion 
    Exp = MemoryManager<MultiRegions::ContField2D>::
        AllocateSharedPtr(graph2D,bcs);
    //----------------------------------------------

    
    //----------------------------------------------
    // Set up coordinates of mesh for Forcing function evaluation
    coordim = Exp->GetCoordim(0);
    nq      = Exp->GetTotPoints();
    
    xc0 = Array<OneD,NekDouble>(nq,0.0);
    xc1 = Array<OneD,NekDouble>(nq,0.0);
    xc2 = Array<OneD,NekDouble>(nq,0.0);
    
    switch(coordim)
    {
    case 1:
        Exp->GetCoords(xc0);
        break;
    case 2:
        Exp->GetCoords(xc0,xc1);
        break;
    case 3:
        Exp->GetCoords(xc0,xc1,xc2);
        break;
    }
    //----------------------------------------------
    
    //----------------------------------------------
    // Define forcing function for first variable defined in file 
    fce = Array<OneD,NekDouble>(nq);
    SpatialDomains::ConstForcingFunctionShPtr ffunc 
        = bcs.GetForcingFunction(bcs.GetVariable(0));

    ffunc->Evaluate(xc0,xc1,xc2,fce);
    //----------------------------------------------

    //----------------------------------------------
    // Setup expansion containing the  forcing function
    Fce = MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //----------------------------------------------

    //----------------------------------------------
    // Set up the variable laplacian coefficients
    bool useVariableCoeffs = false;
    if(bcs.CheckForParameter("UseVariableCoefficients"))
    {
        if((int) bcs.GetParameter("UseVariableCoefficients") != 0)
        {
            useVariableCoeffs = true;
        }
    }

    Array<OneD, Array<OneD,NekDouble> > lapcoeff(3);
    if(useVariableCoeffs)
    {
        std::string lapcoeffstr[3] = {"LaplacianCoefficient_00",
                                      "LaplacianCoefficient_01",
                                      "LaplacianCoefficient_11"};                
        for(i = 0 ; i < 3; i++)
        {
            lapcoeff[i] = Array<OneD,NekDouble>(nq);
            
            SpatialDomains::ConstUserDefinedEqnShPtr cfunc = bcs.GetUserDefinedEqn(lapcoeffstr[i]);
            
            cfunc->Evaluate(xc0,xc1,xc2,lapcoeff[i]);
        }
    }
    //----------------------------------------------

    //----------------------------------------------
    // Laplacian solution taking physical forcing 
    if(useVariableCoeffs)
    {
        Exp->LaplaceSolve(*Fce,lapcoeff);
    }
    else
    {
        Exp->LaplaceSolve(*Fce);
    }
    //----------------------------------------------
    
    //----------------------------------------------
    // Backward Transform Solution to get solved values at 
    Exp->BwdTrans(*Exp);
    //----------------------------------------------
    
    //----------------------------------------------
    // Write solution 
    ofstream outfile("LaplaceFile2D.pos");
    Exp->WriteToFile(outfile,eGmsh);
    outfile.close();

    ofstream outfile2("LaplaceFile2D.dat");
    Exp->WriteToFile(outfile2,eTecplot);
    outfile2.close();
    //----------------------------------------------
    
    //----------------------------------------------
    // See if there is an exact solution, if so 
    // evaluate and plot errors
    SpatialDomains::ConstExactSolutionShPtr ex_sol =
        bcs.GetExactSolution(bcs.GetVariable(0));

    if(ex_sol)
    {
        //----------------------------------------------
        // evaluate exact solution 
        ex_sol->Evaluate(xc0,xc1,xc2,fce);
        //----------------------------------------------

        //--------------------------------------------
        // Calculate L_inf error 
        Fce->SetPhys(fce);
        Fce->SetPhysState(true);


        cout << "L infinity error: " << Exp->Linf(*Fce) << endl;
        cout << "L 2 error:        " << Exp->L2  (*Fce) << endl;
        //--------------------------------------------        
    }
    //----------------------------------------------        
        return 0;
}

