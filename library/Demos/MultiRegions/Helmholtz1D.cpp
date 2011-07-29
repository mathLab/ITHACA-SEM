#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <SpatialDomains/MeshPartition.h>
#include <MultiRegions/ContField1D.h>

using namespace Nektar;

int NoCaseStringCompare(const string & s1, const string& s2);


int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession;
    LibUtilities::CommSharedPtr vComm;
    MultiRegions::ContField1DSharedPtr Exp,Fce;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce;
    Array<OneD,NekDouble>  xc0,xc1,xc2;
    NekDouble  lambda;
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

        // Force Iterative solver for parallel execution
        SolnType = MultiRegions::eIterativeFull;
    }


    if( (argc != 2) && (argc != 3) && (argc != 4))
    {
        fprintf(stderr,"Usage: Helmholtz1D  meshfile \n");
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
        // Read in mesh from input file
        SpatialDomains::MeshGraph1D graph1D;
        graph1D.ReadGeometry(meshfile);
        graph1D.ReadExpansions(meshfile);
        //----------------------------------------------

        //----------------------------------------------
        // read the problem parameters from input file
        SpatialDomains::BoundaryConditions bcs(vSession, &graph1D);
        //----------------------------------------------

        //----------------------------------------------
        // Print summary of solution details
        lambda = vSession->GetParameter("Lambda");
        const SpatialDomains::CompositeMap domain = (graph1D.GetDomain());
        cout << "Solving 1D Helmholtz: "  << endl;
        cout << "       Communication: " << vCommModule << endl;
        cout << "       Solver type  : " << MultiRegions::GlobalSysSolnTypeMap[SolnType] << endl;
        cout << "       Lambda       : " << lambda << endl;
        //----------------------------------------------

        //----------------------------------------------
        // Define Expansion
        int bc_loc = 0;
        Exp = MemoryManager<MultiRegions::ContField1D>::
            AllocateSharedPtr(vComm,graph1D,bcs,bc_loc,SolnType);
        //----------------------------------------------

        //----------------------------------------------
        // Set up coordinates of mesh for Forcing function evaluation
        coordim = Exp->GetCoordim(0);
        nq      = Exp->GetTotPoints();

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
        //----------------------------------------------

        //----------------------------------------------
        // Define forcing function for first variable defined in file
        fce = Array<OneD,NekDouble>(nq);
        LibUtilities::EquationSharedPtr ffunc
                                        = vSession->GetFunction("Forcing", 0);
        for(i = 0; i < nq; ++i)
        {
            fce[i] = ffunc->Evaluate(xc0[i],xc1[i],xc2[i]);
        }
        //----------------------------------------------

        //----------------------------------------------
        // Setup expansion containing the  forcing function
        Fce = MemoryManager<MultiRegions::ContField1D>::AllocateSharedPtr(*Exp);
        Fce->SetPhys(fce);
        //----------------------------------------------

        //----------------------------------------------
        // Helmholtz solution taking physical forcing
        Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), lambda);
        //----------------------------------------------

        //----------------------------------------------
        // Backward Transform Solution to get solved values at
        Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys());
        //----------------------------------------------

        //----------------------------------------------
        // Write solution
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
        for(i = 0; i < FieldDef.size(); ++i)
        {
            FieldDef[i]->m_fields.push_back("u");
            Exp->AppendFieldData(FieldDef[i], FieldData[i]);
        }
        graph1D.Write(out, FieldDef, FieldData);
        //----------------------------------------------

        //----------------------------------------------
        // See if there is an exact solution, if so
        // evaluate and plot errors
        LibUtilities::EquationSharedPtr ex_sol
                                = vSession->GetFunction("ExactSolution", 0);


        if(ex_sol)
        {
            //----------------------------------------------
            // evaluate exact solution
            for(i = 0; i < nq; ++i)
            {
                fce[i] = ex_sol->Evaluate(xc0[i],xc1[i],xc2[i]);
            }
            Fce->SetPhys(fce);
            //----------------------------------------------

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
        cerr << "Caught exception." << endl;
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


