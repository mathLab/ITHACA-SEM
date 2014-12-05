#include <sstream>
#include <LibUtilities/BasicUtils/Timer.h>
#include <iomanip>

#include <boost/filesystem/path.hpp>
#include <SpatialDomains/MeshGraph3D.h>
#include <MultiRegions/DisContField3D.h>

#ifdef NEKTAR_USING_CHUD
#include <CHUD/CHUD.h> 
//#define SHARK
#endif 

using namespace Nektar;

std::string PortablePath(const boost::filesystem::path& path);

int main(int argc, char *argv[])
{
    MultiRegions::DisContField3DSharedPtr Exp,Fce,Sol;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce,sol; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 
    NekDouble  lambda;
    vector<string> vFilenames;
    //defining timing variables
    Timer timer;
    NekDouble exeTime, fullTime, ppTime;

    if(argc < 6)//< is used to be able to submit "verbose" option
    {
        fprintf(stderr,"Usage: TimingHDGHelmSolve3D Type MeshSize NumModes OptimisationLevel\n");
        fprintf(stderr,"    where: - Type is one of the following:\n");
        fprintf(stderr,"                  1: Regular  Hexahedrals \n");
        fprintf(stderr,"                  2: Deformed Hexahedrals (may not be supported) \n");
        fprintf(stderr,"                  3: Regular  Tetrahedrals \n");
        fprintf(stderr,"    where: - MeshSize is 1/h \n");
        fprintf(stderr,"    where: - NumModes is the number of 1D modes of the expansion \n");
        fprintf(stderr,"    where: - OptimisationLevel is one of the following:\n");
        fprintf(stderr,"                  0: Use elemental sum-factorisation evaluation \n");
        fprintf(stderr,"                  2: Use elemental matrix evaluation using blockmatrices \n");
        fprintf(stderr,"                  3: Use global matrix evaluation \n");
        fprintf(stderr,"                  4: Use optimal evaluation (this option requires optimisation-files being set-up) \n");
        fprintf(stderr,"    where: - LinSysSolver is one of the following:\n");
        fprintf(stderr,"                  0: Use DirectStaticCond Solver \n");
        fprintf(stderr,"                  1: Use DirectMultilevelStaticCond Solver\n");
        fprintf(stderr,"                  2: Use IterativeStaticCond Solver \n");
        fprintf(stderr,"                  3: Use IterativeMultilevelStaticCond Solver\n");
        exit(1);
    }

    boost::filesystem::path basePath(BASE_PATH);

    int Type        = atoi(argv[1]);
    std::string TypeStr;
    int MeshSize    = atoi(argv[2]);
    int NumModes    = atoi(argv[3]);
    int optLevel    = atoi(argv[4]);
    std::string optLevelStr;
    int SolverType  = atoi(argv[5]);
    std::string SolverTypeStr;

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
                MeshFileDirectory << "RegularHexMeshes";
                MeshFileName << "UnitCube_RegularHexMesh_h_1_" << MeshSize << ".xml";
                TypeStr = "RHex";
            }
            break;
        case 2:
            {
                MeshFileDirectory << "DeformedHexMeshes";
                MeshFileName << "UnitCube_DeformedHexMesh_h_1_" << MeshSize << ".xml";
                TypeStr = "DHex";
            }
            break;
        case 3:
            {
                MeshFileDirectory << "RegularTetMeshes";
                MeshFileName << "UnitCube_RegularTetMesh_h_1_" << MeshSize << ".xml";
                TypeStr = "RTet";
            }
            break;
        default:
            {
                cerr << "Type should be equal to one of the following values: "<< endl;
                cerr << "  1: Regular Hexes" << endl;
                cerr << "  2: Deformed Hexes" << endl;
                cerr << "  3: Regular Tets" << endl;
                exit(1);
            }
    }

    BCfileName << "UnitCube_DirichletBoundaryConditions.xml";
    ExpansionsFileName << "NektarExpansionsNummodes" << NumModes << ".xml";

    switch(optLevel)
    {
        case 0:
            {
                GlobOptFileName << "NoGlobalMat.xml";
                optLevelStr = "SumFac";
            }
            break;
        case 2:
            {
                GlobOptFileName << "DoBlockMat.xml";
                optLevelStr = "BlkMat";
            }
            break;
        case 3:
            {
                GlobOptFileName << "DoGlobalMat.xml";
                optLevelStr = "GlbMat";
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

    //////////////////////////////////////////////////////////////////////////////////////
    // solution and error computation
    //////////////////////////////////////////////////////////////////////////////////////

    //----------------------------------------------

    LibUtilities::SessionReaderSharedPtr vSession
        = LibUtilities::SessionReader::CreateInstance(argc, argv, vFilenames);

    switch(SolverType)
    {
        case 0:
            {
                vSession->SetSolverInfo("GlobalSysSoln", "DirectStaticCond");
                SolverTypeStr = "DSC";
            }
            break;
        case 1:
            {
                vSession->SetSolverInfo("GlobalSysSoln", "DirectMultiLevelStaticCond");
                SolverTypeStr = "DMSC";
            }
            break;
        case 2:
            {
                vSession->SetSolverInfo("GlobalSysSoln", "IterativeStaticCond");
                SolverTypeStr = "ISC";
            }
            break;
        case 3:
            {
                vSession->SetSolverInfo("GlobalSysSoln", "IterativeMultiLevelStaticCond");
                SolverTypeStr = "IMSC";
            }
            break;
        default:
            {
                ASSERTL0(false,"Unrecognised system solver");
            }
    }

    //timing the whole solve including mesh loading
    timer.Start();
    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graph3D = MemoryManager<SpatialDomains::MeshGraph3D>::AllocateSharedPtr(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    lambda = vSession->GetParameter("Lambda");
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion 
    Exp = MemoryManager<MultiRegions::DisContField3D>::
        AllocateSharedPtr(vSession,graph3D,vSession->GetVariable(0));
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
    ffunc->Evaluate(xc0,xc1,xc2,fce);
    //----------------------------------------------

    //----------------------------------------------
    // Setup expansion containing the  forcing function
    Fce = MemoryManager<MultiRegions::DisContField3D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //----------------------------------------------

    //----------------------------------------------
    // Helmholtz solution taking physical forcing
    FlagList flags;
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorLambda] = lambda;
    factors[StdRegions::eFactorTau] = 1.0;
    Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(),flags,factors);
    //----------------------------------------------

    //----------------------------------------------
    // Backward Transform Solution to get solved values at 
    Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys());
    //----------------------------------------------
    //end of full solve timing
    timer.Stop();
    fullTime = timer.TimePerTest(1);

    //----------------------------------------------
    // See if there is an exact solution, if so 
    // evaluate and plot errors
    LibUtilities::EquationSharedPtr ex_sol = vSession->GetFunction("ExactSolution",0);

    //----------------------------------------------
    // evaluate exact solution 
    sol = Array<OneD,NekDouble>(nq);
    ex_sol->Evaluate(xc0,xc1,xc2,sol);

    NekDouble L2Error    = Exp->L2  (Exp->GetPhys(), sol);
    NekDouble LinfError  = Exp->Linf(Exp->GetPhys(), sol); 

    //////////////////////////////////////////////////////////////////////////////////////
    // Alternative error computation (finer sampling) 
    //////////////////////////////////////////////////////////////////////////////////////

    const LibUtilities::PointsKey PkeyT1(30,LibUtilities::eGaussLobattoLegendre);
    const LibUtilities::PointsKey PkeyT2(30,LibUtilities::eGaussRadauMAlpha1Beta0);
    const LibUtilities::PointsKey PkeyT3(30,LibUtilities::eGaussRadauMAlpha2Beta0);
    const LibUtilities::PointsKey PkeyQ1(30,LibUtilities::eGaussLobattoLegendre);
    const LibUtilities::PointsKey PkeyQ2(30,LibUtilities::eGaussLobattoLegendre);
    const LibUtilities::PointsKey PkeyQ3(30,LibUtilities::eGaussLobattoLegendre);
    const LibUtilities::BasisKey  BkeyT1(LibUtilities::eModified_A,NumModes,PkeyT1);
    const LibUtilities::BasisKey  BkeyT2(LibUtilities::eModified_B,NumModes,PkeyT2);
    const LibUtilities::BasisKey  BkeyT3(LibUtilities::eModified_C,NumModes,PkeyT3);
    const LibUtilities::BasisKey  BkeyQ1(LibUtilities::eModified_A,NumModes,PkeyQ1);
    const LibUtilities::BasisKey  BkeyQ2(LibUtilities::eModified_A,NumModes,PkeyQ2);
    const LibUtilities::BasisKey  BkeyQ3(LibUtilities::eModified_A,NumModes,PkeyQ3);


    MultiRegions::ExpList3DSharedPtr ErrorExp = 
        MemoryManager<MultiRegions::ExpList3D>::AllocateSharedPtr(vSession,BkeyT1,BkeyT2,BkeyT3,BkeyQ1,BkeyQ2,BkeyQ3,graph3D);

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
    ErrorExp->BwdTrans_IterPerExp(Exp->GetCoeffs(),ErrorExp->UpdatePhys());

    NekDouble L2ErrorBis    = ErrorExp->L2  (ErrorExp->GetPhys(), ErrorSol);
    NekDouble LinfErrorBis  = ErrorExp->Linf(ErrorExp->GetPhys(), ErrorSol); 

    //////////////////////////////////////////////////////////////////////////////////////
    // postprocessing and error computation
    //////////////////////////////////////////////////////////////////////////////////////
    //timing postprocessing
    timer.Start();

    int num_points = NumModes + 1;
    //Tetrahedron
    const LibUtilities::PointsKey ppPkeyT1(num_points+1,LibUtilities::eGaussLobattoLegendre);
    const LibUtilities::PointsKey ppPkeyT2(num_points,LibUtilities::eGaussRadauMAlpha1Beta0);//need to doublecheck this one
    const LibUtilities::PointsKey ppPkeyT3(num_points,LibUtilities::eGaussRadauMAlpha2Beta0);//need to doublecheck this one
    LibUtilities::BasisKeyVector  BkeyT;
    BkeyT.push_back(LibUtilities::BasisKey(LibUtilities::eModified_A, NumModes+1, ppPkeyT1));
    BkeyT.push_back(LibUtilities::BasisKey(LibUtilities::eModified_B, NumModes+1, ppPkeyT2));
    BkeyT.push_back(LibUtilities::BasisKey(LibUtilities::eModified_C, NumModes+1, ppPkeyT3));
    //Hexahedron
    const LibUtilities::PointsKey ppPkeyH(num_points+1,LibUtilities::eGaussLobattoLegendre);
    LibUtilities::BasisKeyVector  BkeyH;
    BkeyH.push_back(LibUtilities::BasisKey(LibUtilities::eModified_A, NumModes+1, ppPkeyH));
    BkeyH.push_back(LibUtilities::BasisKey(LibUtilities::eModified_A, NumModes+1, ppPkeyH));
    BkeyH.push_back(LibUtilities::BasisKey(LibUtilities::eModified_A, NumModes+1, ppPkeyH));


    graph3D->SetBasisKey(LibUtilities::eTetrahedron, BkeyT);
    graph3D->SetBasisKey(LibUtilities::eHexahedron, BkeyH);

    MultiRegions::DisContField3DSharedPtr PostProc = 
        MemoryManager<MultiRegions::DisContField3D>::AllocateSharedPtr(vSession,graph3D,vSession->GetVariable(0));

    // calcualte spectral/hp approximation on the quad points of this new
    // expansion basis
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef 
        = Exp->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
    std::string fieldstr = "u";

    for(i = 0; i < FieldDef.size(); ++i)
    {
        FieldDef[i]->m_fields.push_back(fieldstr);
        Exp->AppendFieldData(FieldDef[i], FieldData[i]);
        PostProc->ExtractDataToCoeffs(FieldDef[i],FieldData[i],fieldstr,PostProc->UpdateCoeffs());
    }

    // Interpolation of trace 
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> TraceDef 
        = Exp->GetTrace()->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > TraceData(TraceDef.size());
    for(i = 0; i < TraceDef.size(); ++i)
    {
        TraceDef[i]->m_fields.push_back(fieldstr);
        Exp->GetTrace()->AppendFieldData(TraceDef[i], TraceData[i]);
        PostProc->GetTrace()->ExtractDataToCoeffs(TraceDef[i],TraceData[i],fieldstr,PostProc->GetTrace()->UpdateCoeffs());
    }

    PostProc->BwdTrans_IterPerExp(PostProc->GetCoeffs(),PostProc->UpdatePhys());

    PostProc->EvaluateHDGPostProcessing(PostProc->UpdateCoeffs());
    PostProc->BwdTrans_IterPerExp(PostProc->GetCoeffs(),PostProc->UpdatePhys());
    //end postprocessing timing
    timer.Stop();
    ppTime = timer.TimePerTest(1);

    //computing postprocessing error
    int ppCoordim = PostProc->GetCoordim(0);
    int ppNq      = PostProc->GetTotPoints();

    Array<OneD,NekDouble> ppXc0(ppNq,0.0);
    Array<OneD,NekDouble> ppXc1(ppNq,0.0);
    Array<OneD,NekDouble> ppXc2(ppNq,0.0);

    switch(ppCoordim)
    {
        case 1:
            PostProc->GetCoords(ppXc0);
            break;
        case 2:
            PostProc->GetCoords(ppXc0,ppXc1);
            break;
        case 3:
            PostProc->GetCoords(ppXc0,ppXc1,ppXc2);
            break;
    }


    // evaluate exact solution 
    Array<OneD,NekDouble> ppSol(ppNq);
    ex_sol->Evaluate(ppXc0,ppXc1,ppXc2,ppSol);

    NekDouble L2ErrorPostProc = PostProc->L2(PostProc->GetPhys(), ppSol);
    NekDouble LinfErrorPostProc = PostProc->Linf(PostProc->GetPhys(), ppSol); 

    //////////////////////////////////////////////////////////////////////////////////////
    // Timing section 
    //////////////////////////////////////////////////////////////////////////////////////

    //--------------------------------------------     
#if 0
    cout << "L infinity error: " << LinfErrorBis << endl;
    cout << "L 2 error:        " << L2ErrorBis   << endl;
#endif 

    // We first do a single run in order to estimate the number of calls 
    // we are going to make
    timer.Start();
    Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(),flags,factors);
    Exp->BwdTrans (Exp->GetCoeffs(),Exp->UpdatePhys());
    timer.Stop();
    exeTime = timer.TimePerTest(1);

    int NumCalls = (int) ceil(1.0/exeTime);
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
    chudStartRemotePerfMonitor("TimingHDGHelmSolve3D");
#endif

    timer.Start();
    for(i = 0; i < NumCalls; ++i)
    {
        Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(),flags,factors);
        Exp->BwdTrans (Exp->GetCoeffs(),Exp->UpdatePhys());
    }
    timer.Stop();
    exeTime = timer.TimePerTest(1);

#ifdef SHARK
    chudStopRemotePerfMonitor();
    chudReleaseRemoteAccess();
    chudCleanup();
#endif


    int nLocCoeffs     = Exp->GetTraceMap()->GetNumLocalCoeffs();
    int nGlobCoeffs    = Exp->GetTraceMap()->GetNumGlobalCoeffs();
    int nLocBndCoeffs  = Exp->GetTraceMap()->GetNumLocalBndCoeffs();
    int nGlobBndCoeffs = Exp->GetTraceMap()->GetNumGlobalBndCoeffs();
    int nLocDirCoeffs  = Exp->GetTraceMap()->GetNumLocalDirBndCoeffs();
    int nGlobDirCoeffs = Exp->GetTraceMap()->GetNumGlobalDirBndCoeffs();
    int nGlobBandwidth = Exp->GetTraceMap()->GetBndSystemBandWidth();
    int nGlobBndRank = nGlobBndCoeffs - nGlobDirCoeffs;
    int nnz = 0;

    ofstream outfile("TimingHDGHelmSolve3D.dat");
    outfile.precision(0);
    outfile << setw(10) << SolverTypeStr << " ";
    outfile << setw(10) << optLevelStr << " ";
    outfile << setw(10) << TypeStr << " ";
    outfile << setw(10) << NumElements << " ";
    outfile << setw(10) << NumModes << " ";
    outfile.precision(7);
    outfile << setw(10) << fixed << exeTime << " ";
    outfile.precision(0);
    outfile << setw(10) << NumCalls << " ";
    outfile.precision(7);
    outfile << setw(10) << fixed << ((NekDouble) (exeTime/((NekDouble)NumCalls))) << " ";
    outfile << setw(15) << scientific << noshowpoint << L2Error << " ";
    outfile << setw(15) << scientific << noshowpoint << L2ErrorPostProc << " ";
    outfile << setw(15) << scientific << noshowpoint << L2ErrorBis << " ";
    outfile << setw(15) << scientific << noshowpoint << LinfError << " ";
    outfile << setw(15) << scientific << noshowpoint << LinfErrorPostProc << " ";
    outfile << setw(15) << scientific << noshowpoint << LinfErrorBis << " ";
    outfile << setw(10) << nLocCoeffs  << " ";
    outfile << setw(10) << nGlobCoeffs << " ";
    outfile << setw(10) << nLocBndCoeffs  << " ";
    outfile << setw(10) << nGlobBndCoeffs << " ";
    outfile << setw(10) << nLocDirCoeffs  << " ";
    outfile << setw(10) << nGlobDirCoeffs << " ";
    outfile << setw(10) << nGlobBndRank << " ";
    outfile << setw(10) << nGlobBandwidth << " ";
    outfile << setw(10) << nnz << " ";
    outfile << setw(10) << fixed << fullTime << " ";
    outfile << setw(10) << fixed << ppTime << " ";
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
