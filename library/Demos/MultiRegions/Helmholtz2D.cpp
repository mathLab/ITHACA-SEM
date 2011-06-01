#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <SpatialDomains/MeshPartition.h>
#include <MultiRegions/ContField2D.h>

using namespace Nektar;

//#define TIMING
#ifdef TIMING
#include <time.h>
#define Timing(s) \
 fprintf(stdout,"%s Took %g seconds\n",s,(clock()-st)/cps); \
 st = clock();
#else
#define Timing(s) \
 /* Nothing */
#endif

int NoCaseStringCompare(const string & s1, const string& s2);

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession;
    LibUtilities::CommSharedPtr vComm;
    MultiRegions::ContField2DSharedPtr Exp,Fce;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce;
    Array<OneD,NekDouble>  xc0,xc1,xc2;
    NekDouble  lambda;
    NekDouble    cps = (double)CLOCKS_PER_SEC;
    MultiRegions::GlobalSysSolnType SolnType = MultiRegions::eDirectMultiLevelStaticCond;
    string meshfile(argv[1]);
    string vCommModule("Serial");

    vSession = MemoryManager<LibUtilities::SessionReader>::AllocateSharedPtr(meshfile);

    if (vSession->DefinesSolverInfo("Communication"))
    {
        vCommModule = vSession->GetSolverInfo("Communication");
    }
    else if (LibUtilities::GetCommFactory().ModuleExists("ParallelMPI"))
    {
        vCommModule = "ParallelMPI";
    }

    vComm = LibUtilities::GetCommFactory().CreateInstance(vCommModule,argc,argv);

    if (vComm->GetSize() > 1)
    {
        if (vComm->GetRank() == 0)
        {
            SpatialDomains::MeshPartitionSharedPtr vPartitioner = MemoryManager<SpatialDomains::MeshPartition>::AllocateSharedPtr(vSession);
            vPartitioner->PartitionMesh(vComm->GetSize());
            vPartitioner->WritePartitions(vSession, meshfile);
        }

        vComm->Block();

        meshfile = meshfile + "." + boost::lexical_cast<std::string>(vComm->GetRank());

        // Force use of Iterative solver for parallel execution
        SolnType = MultiRegions::eIterativeFull;
    }

    if( (argc != 2) && (argc != 3) && (argc != 4))
    {
        fprintf(stderr,"Usage: Helmholtz2D meshfile [SysSolnType]   or   \n");
        fprintf(stderr,"Usage: Helmholtz2D meshfile [SysSolnType]"
                " [GlobalOptimizationFile ElementalOptimizationFile]\n");
        exit(1);
    }

    //----------------------------------------------
    // Load the solver type so we can test full solve, static
    // condensation and the default multi-level statis condensation.
    if( argc >= 3 )
    {
        if(!NoCaseStringCompare(argv[2],"MultiLevelStaticCond"))
        {
            SolnType = MultiRegions::eDirectMultiLevelStaticCond;
            cout << "Solution Type: MultiLevel Static Condensation" << endl;
        }
        else if(!NoCaseStringCompare(argv[2],"StaticCond"))
        {
            SolnType = MultiRegions::eDirectStaticCond;
            cout << "Solution Type: Static Condensation" << endl;
        }
        else if(!NoCaseStringCompare(argv[2],"FullMatrix"))
        {
            SolnType = MultiRegions::eDirectFullMatrix;
            cout << "Solution Type: Full Matrix" << endl;
        }
        else if(!NoCaseStringCompare(argv[2],"IterativeCG"))
        {
            SolnType = MultiRegions::eIterativeFull;
            cout << "Solution Type: Iterative Full Matrix" << endl;
        }
        else
        {
            cerr << "SolnType not recognised" <<endl;
            exit(1);
        }

    }
    //----------------------------------------------

    try
    {
        //----------------------------------------------
        // Load the ELEMENTAL optimization parameters if they
        // have been given as argument
        // (use ElementalOptimisationParameters.xml as an example)
        if( argc == 5 )
        {
            string eloptfile(argv[4]);
            NekOptimize::LoadElementalOptimizationParameters(eloptfile);
        }
        //----------------------------------------------

        //----------------------------------------------
        // Read in mesh from input file
        SpatialDomains::MeshGraph2D graph2D;
        graph2D.ReadGeometry(meshfile);
        graph2D.ReadExpansions(meshfile);
        //----------------------------------------------

        //----------------------------------------------
        // read the problem parameters from input file
        SpatialDomains::BoundaryConditions bcs(&graph2D);
        bcs.Read(meshfile);
        //----------------------------------------------

        //----------------------------------------------
        // Print summary of solution details
        lambda = bcs.GetParameter("Lambda");
        const SpatialDomains::ExpansionMap &expansions = graph2D.GetExpansions();
        LibUtilities::BasisKey bkey0 = expansions.begin()->second->m_basisKeyVector[0];
        cout << "Solving 2D Helmholtz: " << endl;
        cout << "         Communication: " << vCommModule << endl;
        cout << "         Solver type  : " << MultiRegions::GlobalSysSolnTypeMap[SolnType] << endl;
        cout << "         Lambda       : " << lambda << endl;
        cout << "         No. modes    : " << bkey0.GetNumModes() << endl;
        cout << endl;
        //----------------------------------------------

        //----------------------------------------------
        // Define Expansion
        int bc_loc = 0;
        Exp = MemoryManager<MultiRegions::ContField2D>::
            AllocateSharedPtr(vComm,graph2D,bcs,bc_loc,SolnType);
        Exp->ReadGlobalOptimizationParameters(meshfile);
        //----------------------------------------------

        Timing("Read files and define exp ..");

        //----------------------------------------------
        // Set up coordinates of mesh for Forcing function evaluation
        coordim = Exp->GetCoordim(0);
        nq      = Exp->GetTotPoints();

        xc0 = Array<OneD,NekDouble>(nq,0.0);
        xc1 = Array<OneD,NekDouble>(nq,0.0);
        xc2 = Array<OneD,NekDouble>(nq,0.0);

        switch(coordim)
        {
        case 2:
            Exp->GetCoords(xc0,xc1);
            break;
        case 3:
            Exp->GetCoords(xc0,xc1,xc2);
            break;
        default:
            ASSERTL0(false,"Coordim not valid");
            break;
        }
        //----------------------------------------------

        //----------------------------------------------
        // Define forcing function for first variable defined in file
        fce = Array<OneD,NekDouble>(nq);
        SpatialDomains::ConstForcingFunctionShPtr ffunc
            = bcs.GetForcingFunction(bcs.GetVariable(0));
        for(i = 0; i < nq; ++i)
        {
            fce[i] = ffunc->Evaluate(xc0[i],xc1[i],xc2[i]);
        }
        //----------------------------------------------

        //----------------------------------------------
        // Setup expansion containing the  forcing function
        Fce = MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(*Exp);
        Fce->SetPhys(fce);
        //----------------------------------------------
        Timing("Define forcing ..");

        //----------------------------------------------
        // Helmholtz solution taking physical forcing
        Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateContCoeffs(), lambda, true);
        //----------------------------------------------
        Timing("Helmholtz Solve ..");

#ifdef TIMING
        for(i = 0; i < 1000; ++i)
        {
            Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateContCoeffs(), lambda, true);
        }

        Timing("1000 Helmholtz Solves:... ");
#endif

        //----------------------------------------------
        // Backward Transform Solution to get solved values
        Exp->BwdTrans(Exp->GetContCoeffs(), Exp->UpdatePhys(), true);
        //----------------------------------------------

        //-----------------------------------------------
        // Write solution to file
        string   out(strtok(argv[1],"."));
        string   endfile(".fld");
        out += endfile;
        if (vComm->GetSize() > 1)
        {
            out += "." + boost::lexical_cast<string>(vComm->GetRank());
        }
        std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef
                                                    = Exp->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

        Exp->GlobalToLocal(Exp->GetContCoeffs(),Exp->UpdateCoeffs());
        for(i = 0; i < FieldDef.size(); ++i)
        {
            FieldDef[i]->m_fields.push_back("u");
            Exp->AppendFieldData(FieldDef[i], FieldData[i]);
        }
        graph2D.Write(out, FieldDef, FieldData);

        //-----------------------------------------------

        //----------------------------------------------
        // See if there is an exact solution, if so
        // evaluate and plot errors
        SpatialDomains::ConstExactSolutionShPtr ex_sol =
            bcs.GetExactSolution(bcs.GetVariable(0));


        if(ex_sol)
        {
            //----------------------------------------------
            // evaluate exact solution
            for(i = 0; i < nq; ++i)
            {
                fce[i] = ex_sol->Evaluate(xc0[i],xc1[i],xc2[i]);
            }
            Fce->SetPhys(fce);
            Fce->SetPhysState(true);
            //--------------------------------------------

            //--------------------------------------------
            // Calculate errors
            NekDouble vLinfError = Exp->Linf(Fce->GetPhys());
            NekDouble vL2Error   = Exp->L2(Fce->GetPhys());
            NekDouble vH1Error   = Exp->H1(Fce->GetPhys());
            if (vComm->GetRank() == 0)
            {
                cout << "L infinity error: " << vLinfError << endl;
                cout << "L 2 error:        " << vL2Error << endl;
                cout << "H 1 error:        " << vH1Error << endl;
            }
            //--------------------------------------------

        }
        //----------------------------------------------
    }
    catch (const std::runtime_error& e)
    {
        cout << "Caught an error" << endl;
        return 1;
    }

    vComm->Finalise();

    return 0;
}



/**
 * Performs a case-insensitive string comparison (from web).
 * @param   s1          First string to compare.
 * @param   s2          Second string to compare.
 * @returns             0 if the strings match.
 */
int NoCaseStringCompare(const string & s1, const string& s2)
{
    string::const_iterator it1=s1.begin();
    string::const_iterator it2=s2.begin();

    //stop when either string's end has been reached
    while ( (it1!=s1.end()) && (it2!=s2.end()) )
    {
        if(::toupper(*it1) != ::toupper(*it2)) //letters differ?
        {
            // return -1 to indicate smaller than, 1 otherwise
            return (::toupper(*it1)  < ::toupper(*it2)) ? -1 : 1;
        }

        //proceed to the next character in each string
        ++it1;
        ++it2;
    }

    size_t size1=s1.size();
    size_t size2=s2.size();// cache lengths

    //return -1,0 or 1 according to strings' lengths
    if (size1==size2)
    {
        return 0;
    }

    return (size1 < size2) ? -1 : 1;
}

