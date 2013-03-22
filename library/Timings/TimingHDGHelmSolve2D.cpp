#include <sstream>
#include <time.h>
#include <sys/time.h>
#include <iomanip>

#include <boost/filesystem/path.hpp>
#include <SpatialDomains/MeshGraph2D.h>
#include <MultiRegions/DisContField2D.h>

#ifdef NEKTAR_USING_CHUD
#include <CHUD/CHUD.h> 
//#define SHARK
#endif 

using namespace Nektar;

std::string PortablePath(const boost::filesystem::path& path);

int main(int argc, char *argv[])
{
    MultiRegions::DisContField2DSharedPtr Exp,Fce,Sol;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce,sol; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 
    NekDouble  lambda;
    vector<string> vFilenames;

    if(argc != 5)
    {
        fprintf(stderr,"Usage: TimingHDGHelmSolve2D Type MeshSize NumModes OptimisationLevel\n");
        fprintf(stderr,"    where: - Type is one of the following:\n");
        fprintf(stderr,"                  1: Regular  Quadrilaterals \n");
        fprintf(stderr,"                  2: Deformed Quadrilaterals (may not be supported) \n");
        fprintf(stderr,"                  3: Regular  Triangles \n");
        fprintf(stderr,"    where: - MeshSize is 1/h \n");
        fprintf(stderr,"    where: - NumModes is the number of 1D modes of the expansion \n");
        fprintf(stderr,"    where: - OptimisationLevel is one of the following:\n");
        fprintf(stderr,"                  0: Use elemental sum-factorisation evaluation \n");
        fprintf(stderr,"                  2: Use elemental matrix evaluation using blockmatrices \n");
        fprintf(stderr,"                  3: Use global matrix evaluation \n");
        fprintf(stderr,"                  4: Use optimal evaluation (this option requires optimisation-files being set-up) \n");
        exit(1);
    }

    boost::filesystem::path basePath(BASE_PATH);

     int Type        = atoi(argv[1]);
     int MeshSize    = atoi(argv[2]);
     int NumModes    = atoi(argv[3]);
     int optLevel    = atoi(argv[4]);

     //----------------------------------------------
     // Retrieve the necessary input files
     stringstream MeshFileName;
     stringstream MeshFileDirectory;
     stringstream BCfileName;
     stringstream ExpansionsFileName;
     stringstream GlobOptFileName;

     switch(Type)
    {
    case 1:
        {
            MeshFileDirectory << "RegularQuadMeshes";
            MeshFileName << "UnitSquare_RegularQuadMesh_h_1_" << MeshSize << ".xml";
        }
        break;
    case 2:
        {
            MeshFileDirectory << "DeformedQuadMeshes";
            MeshFileName << "UnitSquare_DeformedQuadMesh_h_1_" << MeshSize << ".xml";
        }
        break;
    case 3:
        {
            MeshFileDirectory << "RegularTriMeshes";
            MeshFileName << "UnitSquare_RegularTriMesh_h_1_" << MeshSize << ".xml";
        }
        break;
    default:
        {
            cerr << "Type should be equal to one of the following values: "<< endl;
            cerr << "  1: Regular Quads" << endl;
            cerr << "  2: Deformed Quads" << endl;
            cerr << "  3: Regular Tris" << endl;
            exit(1);
        }
    }

    BCfileName << "UnitSquare_DirichletBoundaryConditions.xml";
    ExpansionsFileName << "NektarExpansionsNummodes" << NumModes << ".xml";

    switch(optLevel)
    {
    case 0:
        {
            GlobOptFileName << "NoGlobalMat.xml";
        }
        break;
    case 2:
        {
            GlobOptFileName << "DoBlockMat.xml";
        }
        break;
    case 3:
        {
            GlobOptFileName << "DoGlobalMat.xml";
        }
        break;
    case 4:
        {
            ASSERTL0(false,"Optimisation level not set up");            
        }
        break;
    default:
        {
            ASSERTL0(false,"Unrecognised optimisation level");
        }
    }


    boost::filesystem::path MeshFilePath = basePath / 
        boost::filesystem::path("InputFiles") /
        boost::filesystem::path("Geometry") /
        boost::filesystem::path(MeshFileDirectory.str()) /  
        boost::filesystem::path(MeshFileName.str());
    vFilenames.push_back(PortablePath(MeshFilePath));

    boost::filesystem::path BCfilePath = basePath / 
        boost::filesystem::path("InputFiles") /
        boost::filesystem::path("Conditions") /
        boost::filesystem::path(BCfileName.str());
    vFilenames.push_back(PortablePath(BCfilePath));

    boost::filesystem::path ExpansionsFilePath = basePath / 
        boost::filesystem::path("InputFiles") /
        boost::filesystem::path("Expansions") /
        boost::filesystem::path(ExpansionsFileName.str());
    vFilenames.push_back(PortablePath(ExpansionsFilePath));

    boost::filesystem::path GlobOptFilePath = basePath / 
        boost::filesystem::path("InputFiles") /
        boost::filesystem::path("Optimisation") /
        boost::filesystem::path(GlobOptFileName.str());
    vFilenames.push_back(PortablePath(GlobOptFilePath));

    //----------------------------------------------

    LibUtilities::SessionReaderSharedPtr vSession
                = LibUtilities::SessionReader::CreateInstance(argc, argv, vFilenames);

    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graph2D = MemoryManager<SpatialDomains::MeshGraph2D>::AllocateSharedPtr(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    lambda = vSession->GetParameter("Lambda");
    //----------------------------------------------
   
    //----------------------------------------------
    // Define Expansion 
    Exp = MemoryManager<MultiRegions::DisContField2D>::
        AllocateSharedPtr(vSession,graph2D,vSession->GetVariable(0));
    //----------------------------------------------
    int NumElements = Exp->GetExpSize();
    
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
    LibUtilities::EquationSharedPtr ffunc = vSession->GetFunction("Forcing",0);
    for(i = 0; i < nq; ++i)
    {
        fce[i] = ffunc->Evaluate(xc0[i],xc1[i],xc2[i]);
    }
    //----------------------------------------------

    //----------------------------------------------
    // Setup expansion containing the  forcing function
    Fce = MemoryManager<MultiRegions::DisContField2D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //----------------------------------------------
  
    //----------------------------------------------
    // Helmholtz solution taking physical forcing
    FlagList flags;
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorLambda] = lambda;
    Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(),flags,factors);
    //----------------------------------------------

    //----------------------------------------------
    // Backward Transform Solution to get solved values at 
    Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys());
    //----------------------------------------------

    //----------------------------------------------
    // See if there is an exact solution, if so 
    // evaluate and plot errors
    LibUtilities::EquationSharedPtr ex_sol = vSession->GetFunction("ExactSolution",0);

    //----------------------------------------------
    // evaluate exact solution 
    sol = Array<OneD,NekDouble>(nq);
    for(i = 0; i < nq; ++i)
    {
        sol[i] = ex_sol->Evaluate(xc0[i],xc1[i],xc2[i]);
    }
    //----------------------------------------------
    
    //--------------------------------------------
    // Calculate L_inf error 
    Sol = MemoryManager<MultiRegions::DisContField2D>::AllocateSharedPtr(*Exp);
    Sol->SetPhys(sol);
    Sol->SetPhysState(true);   
    
    NekDouble L2Error    = Exp->L2  (Sol->GetPhys());
    NekDouble LinfError  = Exp->Linf(Sol->GetPhys()); 
    //--------------------------------------------        
    // alternative error calculation
    const LibUtilities::PointsKey PkeyT1(30,LibUtilities::eGaussLobattoLegendre);
    const LibUtilities::PointsKey PkeyT2(30,LibUtilities::eGaussRadauMAlpha1Beta0);
    const LibUtilities::PointsKey PkeyQ1(30,LibUtilities::eGaussLobattoLegendre);
    const LibUtilities::PointsKey PkeyQ2(30,LibUtilities::eGaussLobattoLegendre);
    const LibUtilities::BasisKey  BkeyT1(LibUtilities::eModified_A,NumModes,PkeyT1);
    const LibUtilities::BasisKey  BkeyT2(LibUtilities::eModified_B,NumModes,PkeyT2);
    const LibUtilities::BasisKey  BkeyQ1(LibUtilities::eModified_A,NumModes,PkeyQ1);
    const LibUtilities::BasisKey  BkeyQ2(LibUtilities::eModified_A,NumModes,PkeyQ2);
    
    
    MultiRegions::ExpList2DSharedPtr ErrorExp = 
        MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(vSession,BkeyT1,BkeyT2,BkeyQ1,BkeyQ2,graph2D);
    
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
    for(i = 0; i < ErrorNq; ++i)
    {
        ErrorSol[i] = ex_sol->Evaluate(ErrorXc0[i],ErrorXc1[i],ErrorXc2[i]);
    }
    
    // calcualte spectral/hp approximation on the quad points of this new
    // expansion basis
    ErrorExp->BwdTrans_IterPerExp(Exp->GetCoeffs(),ErrorExp->UpdatePhys());
    
    NekDouble L2ErrorBis    = ErrorExp->L2  (ErrorSol);
    NekDouble LinfErrorBis  = ErrorExp->Linf(ErrorSol); 
    //--------------------------------------------     
#if 0
    cout << "L infinity error: " << LinfErrorBis << endl;
    cout << "L 2 error:        " << L2ErrorBis   << endl;
#endif 
    //----------------------------------------------       

    //----------------------------------------------
    // Do the timings
    timeval timer1, timer2;
    NekDouble time1, time2;
    NekDouble exeTime;

    // We first do a single run in order to estimate the number of calls 
    // we are going to make
    gettimeofday(&timer1, NULL);
    Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(),flags,factors);
    Exp->BwdTrans (Exp->GetCoeffs(),Exp->UpdatePhys());
    gettimeofday(&timer2, NULL);
    time1 = timer1.tv_sec*1000000.0+(timer1.tv_usec);
    time2 = timer2.tv_sec*1000000.0+(timer2.tv_usec);
    exeTime = (time2-time1);

    int NumCalls = (int) ceil(1.0e6/exeTime);
    if(NumCalls < 1)
    {
        NumCalls = 1;
    }

#ifdef SHARK
    NumCalls *= 20;

    chudInitialize();
    chudSetErrorLogFile(stderr);
    chudUmarkPID(getpid(), TRUE);	
    chudAcquireRemoteAccess();
    chudStartRemotePerfMonitor("TimingHDGHelmSolve2D");
#endif

    gettimeofday(&timer1, NULL);
    for(i = 0; i < NumCalls; ++i)
    {
        Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(),flags,factors);
        Exp->BwdTrans (Exp->GetCoeffs(),Exp->UpdatePhys());
    }
    gettimeofday(&timer2, NULL);

#ifdef SHARK
    chudStopRemotePerfMonitor();
    chudReleaseRemoteAccess();
    chudCleanup();
#endif


    time1 = timer1.tv_sec*1000000.0+(timer1.tv_usec);
    time2 = timer2.tv_sec*1000000.0+(timer2.tv_usec);
    exeTime = (time2-time1);

    int nLocCoeffs     = Exp->GetNcoeffs();
    int nGlobCoeffs    = 0;
    int nLocBndCoeffs  = 0;
    int nGlobBndCoeffs = 0;
    int nLocDirCoeffs  = 0;
    int nGlobDirCoeffs = 0;
    int nnz            = 0;

    ofstream outfile("TimingHDGHelmSolve2D.dat");
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

std::string PortablePath(const boost::filesystem::path& path)
{
    boost::filesystem::path temp = path;
    #if BOOST_VERSION > 104200
    temp.make_preferred();
    return temp.string();
    #else
    return temp.file_string();
    #endif

}
