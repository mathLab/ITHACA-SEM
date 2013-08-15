#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <SpatialDomains/MeshGraph2D.h>

using namespace Nektar;

int NoCaseStringCompare(const string & s1, const string& s2);

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    LibUtilities::CommSharedPtr vComm = vSession->GetComm();
    string meshfile(argv[1]);

    MultiRegions::ContField3DHomogeneous1DSharedPtr Exp, Fce;
    int     i, nq;
    Array<OneD,NekDouble>  fce;
    Array<OneD,NekDouble>  xc0,xc1,xc2;
    StdRegions::ConstFactorMap factors;
    FlagList flags;

    if( (argc != 2) && (argc != 3))
    {
        fprintf(stderr,"Usage: Helmholtz3DHomo1D meshfile [SysSolnType]   \n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graph2D = MemoryManager<SpatialDomains::MeshGraph2D>::AllocateSharedPtr(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    int nplanes      = vSession->GetParameter("HomModesZ");
    NekDouble lz     = vSession->GetParameter("LZ");
	bool useFFT = false;
	bool deal = false;
    const LibUtilities::PointsKey Pkey(nplanes,LibUtilities::eFourierEvenlySpaced);
    const LibUtilities::BasisKey Bkey(LibUtilities::eFourier,nplanes,Pkey);
    Exp = MemoryManager<MultiRegions::ContField3DHomogeneous1D>
        ::AllocateSharedPtr(vSession,Bkey,lz,useFFT,deal,graph2D,vSession->GetVariable(0));
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    flags.set(eUseGlobal, true);
    factors[StdRegions::eFactorLambda] = vSession->GetParameter("Lambda");
    const SpatialDomains::ExpansionMap &expansions = graph2D->GetExpansions();
    LibUtilities::BasisKey bkey0 
        = expansions.begin()->second->m_basisKeyVector[0];

    cout << "Solving 3D Helmholtz (Homogeneous in z-direction):"  << endl;
    cout << "         Lambda          : " << factors[StdRegions::eFactorLambda] << endl;
    cout << "         Lz              : " << lz << endl;
    cout << "         No. modes       : " << bkey0.GetNumModes() << endl;
    cout << "         No. hom. modes  : " << Bkey.GetNumModes() << endl;
    cout << endl;
    //----------------------------------------------

    //----------------------------------------------
    // Set up coordinates of mesh for Forcing function evaluation
    nq  = Exp->GetTotPoints();
    xc0 = Array<OneD,NekDouble>(nq,0.0);
    xc1 = Array<OneD,NekDouble>(nq,0.0);
    xc2 = Array<OneD,NekDouble>(nq,0.0);

    Exp->GetCoords(xc0,xc1,xc2);
    //----------------------------------------------

    //----------------------------------------------
    // Define forcing function for first variable defined in file
    fce = Array<OneD,NekDouble>(nq);
    LibUtilities::EquationSharedPtr ffunc = vSession->GetFunction("Forcing", 0);

    ffunc->Evaluate(xc0, xc1, xc2, fce);

    //----------------------------------------------

    //----------------------------------------------
    // Setup expansion containing the  forcing function
    Fce = MemoryManager<MultiRegions::ContField3DHomogeneous1D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //----------------------------------------------

    //----------------------------------------------
    // Helmholtz solution taking physical forcing
    Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), flags, factors);
     //----------------------------------------------

    //----------------------------------------------
    // Backward Transform Solution to get solved values at
    Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys(),MultiRegions::eGlobal);
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

        ex_sol->Evaluate(xc0, xc1, xc2, fce);

        //----------------------------------------------

        //--------------------------------------------
        // Calculate L_inf error
        Fce->SetPhys(fce);
        Fce->SetPhysState(true);

        cout << "L infinity error: " << Exp->Linf(Fce->GetPhys()) << endl;
        cout << "L 2 error:        " << Exp->L2  (Fce->GetPhys()) << endl;
        //--------------------------------------------
    }
    //----------------------------------------------

    //-----------------------------------------------
    // Write solution to file
    string   out(strtok(argv[1],"."));
    string   endfile(".fld");
    out += endfile;
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;

    Exp->GetFieldDefinitions(FieldDef);

    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    Exp->GlobalToLocal();
    for(i = 0; i < FieldDef.size(); ++i)
    {
        FieldDef[i]->m_fields.push_back("u");
        Exp->AppendFieldData(FieldDef[i], FieldData[i]);
    }

    LibUtilities::Write(out, FieldDef, FieldData);
    //-----------------------------------------------

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
