#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/DisContField2D.h>
#include <SpatialDomains/MeshGraph.h>

//#define TIMING
#ifdef TIMING
#include <time.h>
#define Timing(s) \
 fprintf(stdout,"%s Took %g seconds\n",s,(clock()-st)/(double)CLOCKS_PER_SEC); \
 st = clock();
#else
#define Timing(s) \
 /* Nothing */
#endif

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    LibUtilities::CommSharedPtr vComm = vSession->GetComm();

    MultiRegions::DisContField2DSharedPtr Exp,Fce;
    MultiRegions::ExpListSharedPtr DerExp1,DerExp2;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce;
    Array<OneD,NekDouble>  xc0,xc1,xc2;
    StdRegions::ConstFactorMap factors;

    if(argc < 2)
    {
        fprintf(stderr,"Usage: PostProcHDG2D  meshfile [solntype]\n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graph2D = SpatialDomains::MeshGraph::Read(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    factors[StdRegions::eFactorLambda] = vSession->GetParameter("Lambda");
    factors[StdRegions::eFactorTau] = 1.0;
    const SpatialDomains::ExpansionMap &expansions = graph2D->GetExpansions();
    LibUtilities::BasisKey bkey0
                            = expansions.begin()->second->m_basisKeyVector[0];

	//WILL NEED ADJUSTMENT FOR TRIANGLE, PROBABLY
	int num_modes = bkey0.GetNumModes();
	int num_points = bkey0.GetNumPoints();

    if (vComm->GetRank() == 0)
    {
        cout << "Solving 2D Helmholtz: " << endl;
        cout << "         Communication: " << vSession->GetComm()->GetType() << endl;
        cout << "         Solver type  : " << vSession->GetSolverInfo("GlobalSysSoln") << endl;
        cout << "         Lambda       : " << factors[StdRegions::eFactorLambda] << endl;
        cout << "         No. modes    : " << bkey0.GetNumModes() << endl;
        cout << endl;
    }

    //----------------------------------------------
    // Define Expansion
    //----------------------------------------------
    Exp = MemoryManager<MultiRegions::DisContField2D>::
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
    LibUtilities::EquationSharedPtr ffunc
                                    = vSession->GetFunction("Forcing", 0);

    ffunc->Evaluate(xc0, xc1, xc2, fce);

    //----------------------------------------------


    //----------------------------------------------
    // Setup expansion containing the  forcing function
    Fce = MemoryManager<MultiRegions::DisContField2D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //----------------------------------------------
    Timing("Define forcing ..");

    //----------------------------------------------
    // Helmholtz solution taking physical forcing
    Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), factors);
    //----------------------------------------------

    Timing("Helmholtz Solve ..");

    //-----------------------------------------------
    // Backward Transform Solution to get solved values at
    Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys());
    //-----------------------------------------------
    Timing("Backward Transform ..");

    //-----------------------------------------------
    // Write solution to file
    //string out = vSession->GetSessionName() + ".fld";
    //std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
    //    = Exp->GetFieldDefinitions();
    //std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    //for(i = 0; i < FieldDef.size(); ++i)
    //{
    //    FieldDef[i]->m_fields.push_back("u");
    //    Exp->AppendFieldData(FieldDef[i], FieldData[i]);
    //}
    //LibUtilities::Write(out, FieldDef, FieldData);
    //--------------------------------------------        
    //-----------------------------------------------
    // See if there is an exact solution, if so
    // evaluate and plot errors
    LibUtilities::EquationSharedPtr ex_sol =
        vSession->GetFunction("ExactSolution", 0);

	//----------------------------------------------
	// evaluate exact solution
	ex_sol->Evaluate(xc0, xc1, xc2, fce);

	//----------------------------------------------

	//Triangles
	const LibUtilities::PointsKey PkeyT1(num_points+1,LibUtilities::eGaussLobattoLegendre);
	const LibUtilities::PointsKey PkeyT2(num_points,LibUtilities::eGaussRadauMAlpha1Beta0);//need to doublecheck this one
	LibUtilities::BasisKeyVector  BkeyT;
	BkeyT.push_back(LibUtilities::BasisKey(LibUtilities::eModified_A, num_modes+1, PkeyT1));
	BkeyT.push_back(LibUtilities::BasisKey(LibUtilities::eModified_B, num_modes+1, PkeyT2));
	//Quadrilaterals
	const LibUtilities::PointsKey PkeyQ(num_points+1,LibUtilities::eGaussLobattoLegendre);
	LibUtilities::BasisKeyVector  BkeyQ;
	BkeyQ.push_back(LibUtilities::BasisKey(LibUtilities::eModified_A, num_modes+1, PkeyQ));
	BkeyQ.push_back(LibUtilities::BasisKey(LibUtilities::eModified_A, num_modes+1, PkeyQ));


	graph2D->SetBasisKey(LibUtilities::eTriangle, BkeyT);
	graph2D->SetBasisKey(LibUtilities::eQuadrilateral, BkeyQ);

	MultiRegions::DisContField2DSharedPtr PostProc = 
		MemoryManager<MultiRegions::DisContField2D>::AllocateSharedPtr(vSession,graph2D,vSession->GetVariable(0));

	int ErrorCoordim = PostProc->GetCoordim(0);
	int ErrorNq      = PostProc->GetTotPoints();

	Array<OneD,NekDouble> ErrorXc0(ErrorNq,0.0);
	Array<OneD,NekDouble> ErrorXc1(ErrorNq,0.0);
	Array<OneD,NekDouble> ErrorXc2(ErrorNq,0.0);

	switch(ErrorCoordim)
	{
		case 1:
			PostProc->GetCoords(ErrorXc0);
			break;
		case 2:
			PostProc->GetCoords(ErrorXc0,ErrorXc1);
			break;
		case 3:
			PostProc->GetCoords(ErrorXc0,ErrorXc1,ErrorXc2);
			break;
	}
        
        
	// evaluate exact solution 
	Array<OneD,NekDouble> ppSol(ErrorNq);
	ex_sol->Evaluate(ErrorXc0,ErrorXc1,ErrorXc2,ppSol);

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
	
	NekDouble vLinfError = Exp->Linf(Exp->GetPhys(), fce);
	NekDouble vL2Error   = Exp->L2  (Exp->GetPhys(), fce);
	NekDouble L2ErrorPostProc = PostProc->L2(PostProc->GetPhys(), ppSol);
	NekDouble LinfErrorPostProc = PostProc->Linf(PostProc->GetPhys(), ppSol); 

	if (vSession->GetComm()->GetRank() == 0)
	{
		cout << "L infinity error : " << vLinfError << endl;
		cout << "L 2 error        : " << vL2Error   << endl;
		cout << "Postprocessed L infinity error : " << LinfErrorPostProc << endl;
		cout << "Postprocessed L 2 error        : " << L2ErrorPostProc   << endl;
	}

	vSession->Finalise();
    
    return 0;
}

