#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/ContField1D.h>
#include <SpatialDomains/MeshGraph.h>

using namespace std;
using namespace Nektar;

int NoCaseStringCompare(const string & s1, const string& s2);


int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    LibUtilities::CommSharedPtr vComm = vSession->GetComm();
    MultiRegions::ContField1DSharedPtr Exp,Fce;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce;
    Array<OneD,NekDouble>  xc0,xc1,xc2;
    StdRegions::ConstFactorMap factors;

    if( (argc != 2) && (argc != 3) && (argc != 4))
    {
        fprintf(stderr,"Usage: Helmholtz1D  meshfile \n");
        exit(1);
    }

    try
    {
        LibUtilities::FieldIOSharedPtr fld =
            LibUtilities::FieldIO::CreateDefault(vSession);

        //----------------------------------------------
        // Read in mesh from input file
        SpatialDomains::MeshGraphSharedPtr graph1D =
            SpatialDomains::MeshGraph::Read(vSession);
        //----------------------------------------------

        //----------------------------------------------
        // Print summary of solution details
        factors[StdRegions::eFactorLambda] = vSession->GetParameter("Lambda");
        const SpatialDomains::ExpansionMap &expansions = graph1D->GetExpansions();
        LibUtilities::BasisKey bkey0 = expansions.begin()->second->m_basisKeyVector[0];

        if (vComm->GetRank() ==0)
        {
            cout << "Solving 1D Helmholtz: "  << endl;
            cout << "       Communication: " << vComm->GetType() << endl;
            cout << "       Solver type  : " << vSession->GetSolverInfo("GlobalSysSoln") << endl;
            cout << "       Lambda       : " << factors[StdRegions::eFactorLambda] << endl;
            cout << "       No. modes    : " << bkey0.GetNumModes() << endl;
        }
        //----------------------------------------------

        //----------------------------------------------
        // Define Expansion
        Exp = MemoryManager<MultiRegions::ContField1D>::
            AllocateSharedPtr(vSession,graph1D,vSession->GetVariable(0));
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

        ffunc->Evaluate(xc0,xc1,xc2, fce);

        //----------------------------------------------

        //----------------------------------------------
        // Setup expansion containing the  forcing function
        Fce = MemoryManager<MultiRegions::ContField1D>::AllocateSharedPtr(*Exp);
        Fce->SetPhys(fce);
        //----------------------------------------------

        //----------------------------------------------
        //Helmholtz solution taking physical forcing after setting
        //initial condition to zero
        Vmath::Zero(Exp->GetNcoeffs(),Exp->UpdateCoeffs(),1);
        Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), factors);
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
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
            = Exp->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
        for(i = 0; i < FieldDef.size(); ++i)
        {
            FieldDef[i]->m_fields.push_back("u");
            Exp->AppendFieldData(FieldDef[i], FieldData[i]);
        }
        fld->Write(out, FieldDef, FieldData);
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

            ex_sol->Evaluate(xc0,xc1,xc2, fce);

            Fce->SetPhys(fce);
            //----------------------------------------------

            //--------------------------------------------
            // Calculate errors
            NekDouble vLinfError = Exp->Linf(Exp->GetPhys(), Fce->GetPhys());
            NekDouble vL2Error   = Exp->L2(Exp->GetPhys(), Fce->GetPhys());
            NekDouble vH1Error   = Exp->H1(Exp->GetPhys(), Fce->GetPhys());
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
    catch (const std::runtime_error&)
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
/*
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

*/
