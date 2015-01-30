#include <sstream>
#include <time.h>

#include <sys/time.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <SpatialDomains/MeshPartition.h>
#include <MultiRegions/ContField2D.h>


using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession;
    LibUtilities::CommSharedPtr vComm;
    MultiRegions::ContField2DSharedPtr Exp,Fce,Sol;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce,sol; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 
    NekDouble  lambda;
    MultiRegions::GlobalSysSolnType SolnType = MultiRegions::eDirectMultiLevelStaticCond;
    string vCommModule("Serial");

    if(argc != 5)
    {
        fprintf(stderr,"Usage: TimingCGHelmSolve2D Type NumElements NumModes OptimisationLevel\n");
        exit(1);
    }

    int Type        = atoi(argv[1]);
    int NumElements = atoi(argv[2]);
    int NumModes    = atoi(argv[3]);
    int optLevel    = atoi(argv[4]);
    // optLevel = 0  -->  elemental and matrix free
    // optLevel = 2  -->  block matrix operations 
    // optLevel = 3  -->  global matrix operation 
    // optLevel = 4  -->  optimal implementation strategies (optimal evaluation for every different evaluation)

    //---------------------
    // Do a fix in case optLevel 3 is chosen. If the number of elements
    // and the polynomial order is chosen too high, change OptLevel to 2.
    // Otherwise, it may cause memory overflow
    int optLevelVal = optLevel;
    if( optLevel==3 )
    {
        if( (Type<3) && (NumElements>200) && (NumModes>7) )
        {
            return 0;
        }
        if( (Type==3) && (NumElements>14) && (NumModes>10) )
        {
            return 0;
        }
    }
    //----------------------------------------------

    //----------------------------------------------
    // Retrieve the necessary input files
    stringstream MeshFileSS;
    stringstream BCfileSS;
    stringstream ExpansionsFileSS;
    
    int noBis = 1;
    
    switch(Type)
    {
    case 1:
        {
            MeshFileSS << "/Users/ssherw/HDG/Meshes/RegularQuadMeshes/";
            MeshFileSS << "UnitSquare_RegularQuadMesh_" << NumElements << "Elements.xml";
            BCfileSS << "/Users/ssherw/HDG/Meshes/Conditions/UnitSquare_DirichletBoundaryConditions.xml";
        }
        break;
    case 2:
        {
            MeshFileSS << "/Users/ssherw/HDG/Meshes/DeformedQuadMeshes/";
            MeshFileSS << "UnitSquare_DeformedQuadMesh_" << NumElements << "Elements.xml";
            BCfileSS << "/Users/ssherw/HDG/Meshes/Conditions/UnitSquare_DirichletBoundaryConditions.xml";
        }
        break;
    case 3:
        {
            MeshFileSS << "/Users/ssherw/HDG/Meshes/RegularTriMeshes/";
            MeshFileSS << "UnitSquare_RegularTriMesh_h_1_" << NumElements << ".xml";
            BCfileSS << "/Users/ssherw/HDG/Meshes/Conditions/UnitSquare_DirichletBoundaryConditions.xml";
        }
        break;
    case 4:
        {
            MeshFileSS << "/Users/ssherw/HDG/Meshes/Kirby/kirby_10K.xml";

            BCfileSS << "/Users/ssherw/HDG/Meshes/Kirby/kirby_10K.xml";

            noBis = 0;
        }
        break;
    default:
        {
            cerr << "Type should be equal to one of the following values: "<< endl;
            cerr << "  1: Regular Quads" << endl;
            cerr << "  2: Deformed Quads" << endl;
            cerr << "  3: Regular Tris" << endl;
            cerr << "  4: Kirby mesh" << endl;
            exit(1);
        }
    }


    ExpansionsFileSS << "/Users/ssherw/HDG/Meshes/Expansions/NektarExpansionsNummodes";
    ExpansionsFileSS << NumModes << ".xml";

    string meshfile      = MeshFileSS.str();
    string expansionfile = ExpansionsFileSS.str();
    string bcfile        = BCfileSS.str();

    vSession = MemoryManager<LibUtilities::SessionReader>::AllocateSharedPtr(meshfile);


    string globoptfile;

    switch(optLevel)
    {
    case 0:
        {
            stringstream GlobOptFileSS;
            GlobOptFileSS << "/Users/ssherw/HDG/Meshes/Optimisation/NoGlobalMat.xml";
            globoptfile = GlobOptFileSS.str();
        }
        break;
    case 2:
        {
            stringstream GlobOptFileSS;
            GlobOptFileSS << "/Users/ssherw/HDG/Meshes/Optimisation/DoBlockMat.xml";
            globoptfile = GlobOptFileSS.str();
        }
        break;
    case 3:
        {
            stringstream GlobOptFileSS;
            GlobOptFileSS << "/Users/ssherw/HDG/Meshes/Optimisation/DoGlobalMat.xml";
            globoptfile = GlobOptFileSS.str();
        }
        break;
    case 4:
        {


            switch(Type)
            {
            case 1:
                {
                    stringstream GlobOptFileSS;
                    GlobOptFileSS << "/Users/ssherw/HDG/Meshes/Optimisation/";
                    GlobOptFileSS << "UnitSquare_RegularQuadMesh_" << NumElements << "Elements_" << NumModes << "Modes_GlobOpt.xml";
                    globoptfile = GlobOptFileSS.str();
                }
                break;
            case 3:
                {
                    stringstream GlobOptFileSS;
                    GlobOptFileSS << "/Users/ssherw/HDG/Meshes/Optimisation/";
                    GlobOptFileSS << "UnitSquare_RegularTriMesh_h_1_" << NumElements << "_" << NumModes << "Modes_GlobOpt.xml";
                    globoptfile = GlobOptFileSS.str();
                }
                break;
            default:
                {
                    cerr << "No optimal strategy defined for this case "<< endl;
                    exit(1);
                }
            }
            
        }
        break;
    default:
        {
            ASSERTL0(false,"Unrecognised optimisation level");
        }
    }
    //----------------------------------------------


    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraph2D graph2D; 
    graph2D.ReadGeometry(meshfile);
    graph2D.ReadExpansions(expansionfile);
    //----------------------------------------------

    //----------------------------------------------
    // read the problem parameters from input file
    SpatialDomains::BoundaryConditions bcs(&graph2D); 
    bcs.Read(bcfile);
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    lambda = bcs.GetParameter("Lambda");
    //----------------------------------------------
   
    //----------------------------------------------
    // Define Expansion 
    Exp = MemoryManager<MultiRegions::ContField2D>::
        AllocateSharedPtr(vSession,graph2D,bcs);
    //----------------------------------------------
    //    NumElements = Exp->GetExpSize();

    //----------------------------------------------
    // load global optimisation parameters
    Exp->ReadGlobalOptimizationParameters(globoptfile);
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
    // Setup expansion containing the  forcing function
    Fce = MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //----------------------------------------------
  
    //----------------------------------------------
    // Helmholtz solution taking physical forcing 
    Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateContCoeffs(),lambda,true);
    //----------------------------------------------

    //----------------------------------------------
    // Backward Transform Solution to get solved values at 
    Exp->BwdTrans(Exp->GetContCoeffs(), Exp->UpdatePhys(),true);
    //----------------------------------------------

    //----------------------------------------------
    // See if there is an exact solution, if so 
    // evaluate and plot errors
    SpatialDomains::ConstExactSolutionShPtr ex_sol =
        bcs.GetExactSolution(bcs.GetVariable(0));

    //----------------------------------------------
    // evaluate exact solution 
    sol = Array<OneD,NekDouble>(nq);
    ex_sol->Evaluate(xc0,xc1,xc2,sol);
    //----------------------------------------------
    
    //--------------------------------------------
    // Calculate L_inf error 
    Sol = MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(*Exp);
    Sol->SetPhys(sol);
    Sol->SetPhysState(true);   
    
    NekDouble L2Error    = Exp->L2  (Sol->GetPhys());
    NekDouble H1Error    = Exp->H1  (Sol->GetPhys());
    NekDouble LinfError  = Exp->Linf(Sol->GetPhys()); 
    //--------------------------------------------        
    // alternative error calculation
    NekDouble L2ErrorBis;
    NekDouble H1ErrorBis;
    NekDouble LinfErrorBis; 
    if(noBis)
    {
        const LibUtilities::PointsKey PkeyT1(30,LibUtilities::eGaussLobattoLegendre);
        const LibUtilities::PointsKey PkeyT2(30,LibUtilities::eGaussRadauMAlpha1Beta0);
        const LibUtilities::PointsKey PkeyQ1(30,LibUtilities::eGaussLobattoLegendre);
        const LibUtilities::PointsKey PkeyQ2(30,LibUtilities::eGaussLobattoLegendre);
        LibUtilities::BasisKeyVector  BkeyT;
        LibUtilities::BasisKeyVector  BkeyQ;
        BkeyT.push_back(LibUtilities::BasisKey(LibUtilities::eModified_A,
                                               NumModes, PkeyT1));
        BkeyT.push_back(LibUtilities::BasisKey(LibUtilities::eModified_B,
                                               NumModes, PkeyT2));
        BkeyQ.push_back(LibUtilities::BasisKey(LibUtilities::eModified_A,
                                               NumModes, PkeyQ1));
        BkeyQ.push_back(LibUtilities::BasisKey(LibUtilities::eModified_A,
                                               NumModes, PkeyQ2));

        graph2D.SetBasisKey(SpatialDomains::eTriangle, BkeyT);
        graph2D.SetBasisKey(SpatialDomains::eQuadrilateral, BkeyQ);

        MultiRegions::ExpList2DSharedPtr ErrorExp = 
            MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(vComm,
                                                                      graph2D);
    
    
        int ErrorCoordim = ErrorExp->GetCoordim(0);
        int ErrorNq      = ErrorExp->GetTotPoints();
        
        Array<OneD,NekDouble> ErrorXc0(ErrorNq,0.0);
        Array<OneD,NekDouble> ErrorXc1(ErrorNq,0.0);
        Array<OneD,NekDouble> ErrorXc2(ErrorNq,0.0);
        
        switch(ErrorCoordim)
        {
        case 1:
            ErrorExp->GetCoords(ErrorXc0);
            break;
        case 2:
            ErrorExp->GetCoords(ErrorXc0,ErrorXc1);
            break;
        case 3:
            ErrorExp->GetCoords(ErrorXc0,ErrorXc1,ErrorXc2);
            break;
        }

        // evaluate exact solution 
        Array<OneD,NekDouble> ErrorSol(ErrorNq);
        ex_sol->Evaluate(ErrorXc0,ErrorXc1,ErrorXc2,ErrorSol);

        // calcualte spectral/hp approximation on the quad points of this new
        // expansion basis
        Exp->GlobalToLocal(Exp->GetContCoeffs(),ErrorExp->UpdateCoeffs());
        ErrorExp->BwdTrans_IterPerExp(ErrorExp->GetCoeffs(),ErrorExp->UpdatePhys());

        L2ErrorBis    = ErrorExp->L2  (ErrorSol);
        H1ErrorBis    = ErrorExp->H1  (ErrorSol);
        LinfErrorBis  = ErrorExp->Linf(ErrorSol); 
    }
    else
    {
        L2ErrorBis    = L2Error;
        H1ErrorBis    = H1Error;
        LinfErrorBis  = LinfError; 
    }
    //----------------------------------------------       

    //----------------------------------------------
    // Do the timings
    timeval timer1, timer2;
    NekDouble time1, time2;
    NekDouble exeTime;

    // We first do a single run in order to estimate the number of calls 
    // we are going to make
    gettimeofday(&timer1, NULL);
    Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateContCoeffs(),lambda,true);
    Exp->BwdTrans (Exp->GetContCoeffs(),Exp->UpdatePhys(),true);
    gettimeofday(&timer2, NULL);
    time1 = timer1.tv_sec*1000000.0+(timer1.tv_usec);
    time2 = timer2.tv_sec*1000000.0+(timer2.tv_usec);
    exeTime = (time2-time1);

    int NumCalls = (int) ceil(1.0e6/exeTime);
    if(NumCalls < 1)
    {
        NumCalls = 1;
    }

    gettimeofday(&timer1, NULL);
    for(i = 0; i < NumCalls; ++i)
    {
        Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateContCoeffs(),lambda,true);
        Exp->BwdTrans (Exp->GetContCoeffs(),Exp->UpdatePhys(),true);
    }
    gettimeofday(&timer2, NULL);

    time1 = timer1.tv_sec*1000000.0+(timer1.tv_usec);
    time2 = timer2.tv_sec*1000000.0+(timer2.tv_usec);
    exeTime = (time2-time1);

    int nLocCoeffs     = Exp->GetNcoeffs();
    int nGlobCoeffs    = Exp->GetContNcoeffs();
    int nLocBndCoeffs  = Exp->GetLocalToGlobalMap()->GetNumLocalBndCoeffs();
    int nGlobBndCoeffs = Exp->GetLocalToGlobalMap()->GetNumGlobalBndCoeffs();
    int nLocDirCoeffs  = Exp->GetLocalToGlobalMap()->GetNumLocalDirBndCoeffs();
    int nGlobDirCoeffs = Exp->GetLocalToGlobalMap()->GetNumGlobalDirBndCoeffs();
    MultiRegions::GlobalMatrixKey key(StdRegions::eHelmholtz,lambda,Exp->GetLocalToGlobalMap());
    int nnz            = Exp->GetGlobalMatrixNnz(key);

    ofstream outfile("TimingCGHelmSolve2D_NektarppResult.dat");
    outfile.precision(0);
    outfile << setw(10) << Type << " ";
    outfile << setw(10) << NumElements << " ";
    outfile << setw(10) << NumModes << " ";
    outfile << setw(10) << NumCalls << " ";
    outfile << setw(10) << fixed << noshowpoint << exeTime << " ";
    outfile << setw(10) << fixed << noshowpoint << ((NekDouble) (exeTime/((NekDouble)NumCalls))) << " ";
    outfile.precision(7);
    outfile << setw(15) << scientific << noshowpoint << L2Error << " ";
    outfile << setw(15) << scientific << noshowpoint << L2ErrorBis << " ";
    outfile << setw(15) << scientific << noshowpoint << LinfError << " ";
    outfile << setw(15) << scientific << noshowpoint << LinfErrorBis << " ";
    outfile << setw(15) << scientific << noshowpoint << H1Error << " ";
    outfile << setw(15) << scientific << noshowpoint << H1ErrorBis << " ";
    outfile << setw(10) << nLocCoeffs  << " ";
    outfile << setw(10) << nGlobCoeffs << " ";
    outfile << setw(10) << nLocBndCoeffs  << " ";
    outfile << setw(10) << nGlobBndCoeffs << " ";
    outfile << setw(10) << nLocDirCoeffs  << " ";
    outfile << setw(10) << nGlobDirCoeffs << " ";
    outfile << setw(10) << nnz << " ";
    outfile << setw(10) << optLevel << " ";
    outfile << endl;

    outfile.close();
    //----------------------------------------------
    
    return 0;
}

