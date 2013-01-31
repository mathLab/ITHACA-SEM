#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/ContField2D.h>
#include <SpatialDomains/MeshGraph2D.h>

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
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    MultiRegions::ContField2DSharedPtr Exp,Fce;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce;
    Array<OneD,NekDouble>  xc0,xc1,xc2;
    StdRegions::ConstFactorMap factors;
    StdRegions::VarCoeffMap varcoeffs;
    FlagList flags;
    NekDouble    cps = (double)CLOCKS_PER_SEC;
    NekDouble    st;

    if( (argc != 2) && (argc != 3) && (argc != 4))
    {
        fprintf(stderr,"Usage: Helmholtz2D meshfile [SysSolnType]   or   \n");
        exit(1);
    }

    try
    {
        //----------------------------------------------
        // Read in mesh from input file
        SpatialDomains::MeshGraphSharedPtr graph2D = MemoryManager<SpatialDomains::MeshGraph2D>::AllocateSharedPtr(vSession);
        //----------------------------------------------

        //----------------------------------------------
        // Print summary of solution details
        flags.set(eUseGlobal, true);
        factors[StdRegions::eFactorLambda] = vSession->GetParameter("Lambda");
        const SpatialDomains::ExpansionMap &expansions = graph2D->GetExpansions();
        LibUtilities::BasisKey bkey0 = expansions.begin()->second->m_basisKeyVector[0];
        cout << "Solving 2D Helmholtz: " << endl;
        cout << "         Communication: " << vSession->GetComm()->GetType() << endl;
        cout << "         Solver type  : " << vSession->GetSolverInfo("GlobalSysSoln") << endl;
        cout << "         Lambda       : " << factors[StdRegions::eFactorLambda] << endl;
        cout << "         No. modes    : " << bkey0.GetNumModes() << endl;
        cout << endl;
        //----------------------------------------------

        //----------------------------------------------
        // Define Expansion
        Exp = MemoryManager<MultiRegions::ContField2D>::
            AllocateSharedPtr(vSession,graph2D,vSession->GetVariable(0));
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
        // Set up variable coefficients if defined
        if (vSession->DefinesFunction("d00"))
        {
            Array<OneD, NekDouble> d00(nq,0.0);
            LibUtilities::EquationSharedPtr d00func = vSession->GetFunction("d00",0);
            d00func->Evaluate(xc0, xc1, xc2, d00);
            varcoeffs[StdRegions::eVarCoeffD00] = d00;
        }
        if (vSession->DefinesFunction("d11"))
        {
            Array<OneD, NekDouble> d11(nq,0.0);
            LibUtilities::EquationSharedPtr d11func = vSession->GetFunction("d11",0);
            d11func->Evaluate(xc0, xc1, xc2, d11);
            varcoeffs[StdRegions::eVarCoeffD11] = d11;
        }
        //----------------------------------------------

        //----------------------------------------------
        // Define forcing function for first variable defined in file
        fce = Array<OneD,NekDouble>(nq);
        LibUtilities::EquationSharedPtr ffunc = vSession->GetFunction("Forcing",0);
        ffunc->Evaluate(xc0, xc1, xc2, fce);

        //----------------------------------------------

        //----------------------------------------------
        // Setup expansion containing the  forcing function
        Fce = MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(*Exp);
        Fce->SetPhys(fce);
        //----------------------------------------------
        Timing("Define forcing ..");

        //---------------------------------------------- 
        //Helmholtz solution taking physical forcing after setting
        //initial condition to zero
        Vmath::Zero(Exp->GetNcoeffs(),Exp->UpdateCoeffs(),1);
        Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), flags, factors, varcoeffs);
        //----------------------------------------------
        Timing("Helmholtz Solve ..");

#ifdef TIMING
        for(i = 0; i < 1000; ++i)
        {
            Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), flags, factors, varcoeffs);
        }

        Timing("1000 Helmholtz Solves:... ");
#endif

        //----------------------------------------------
        // Backward Transform Solution to get solved values
        Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys(), MultiRegions::eGlobal);
        //----------------------------------------------

        //-----------------------------------------------
        // Write solution to file
        string out = vSession->GetSessionName();
        if (vSession->GetComm()->GetSize() > 1)
        {
            out += "_P" + boost::lexical_cast<string>(vSession->GetComm()->GetRank());
        }
        out += ".fld";
        std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef
                                                    = Exp->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

        Exp->GlobalToLocal(Exp->GetCoeffs(),Exp->UpdateCoeffs());
        for(i = 0; i < FieldDef.size(); ++i)
        {
            FieldDef[i]->m_fields.push_back("u");
            Exp->AppendFieldData(FieldDef[i], FieldData[i]);
        }
        graph2D->Write(out, FieldDef, FieldData);
        //-----------------------------------------------

        //----------------------------------------------
        // See if there is an exact solution, if so
        // evaluate and plot errors
        LibUtilities::EquationSharedPtr ex_sol = vSession->GetFunction("ExactSolution",0);

        if(ex_sol)
        {
            //----------------------------------------------
            // evaluate exact solution
            ex_sol->Evaluate(xc0, xc1, xc2, fce);

            Fce->SetPhys(fce);
            Fce->SetPhysState(true);
            //--------------------------------------------

            //--------------------------------------------
            // Calculate errors
            NekDouble vLinfError = Exp->Linf(Fce->GetPhys());
            NekDouble vL2Error   = Exp->L2(Fce->GetPhys());
            NekDouble vH1Error   = Exp->H1(Fce->GetPhys());
            if (vSession->GetComm()->GetRank() == 0)
            {
                cout << "L infinity error: " << vLinfError << endl;
                cout << "L 2 error:        " << vL2Error << endl;
                cout << "H 1 error:        " << vH1Error << endl;
            }
            //--------------------------------------------

        }
        //----------------------------------------------
    }
    catch (const std::runtime_error&)
    {
        cout << "Caught an error" << endl;
        return 1;
    }

    vSession->Finalise();

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

