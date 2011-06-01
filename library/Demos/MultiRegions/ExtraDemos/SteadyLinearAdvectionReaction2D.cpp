#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <SpatialDomains/MeshPartition.h>
#include <MultiRegions/ContField2D.h>

using namespace Nektar;

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
    LibUtilities::CommSharedPtr vComm = LibUtilities::GetCommFactory().CreateInstance("ParallelMPI",argc,argv);

    string meshfile(argv[1]);

    if (vComm->GetSize() > 1)
    {
        if (vComm->GetRank() == 0)
        {
            LibUtilities::SessionReaderSharedPtr vSession = MemoryManager<LibUtilities::SessionReader>::AllocateSharedPtr(meshfile);
            SpatialDomains::MeshPartitionSharedPtr vPartitioner = MemoryManager<SpatialDomains::MeshPartition>::AllocateSharedPtr(vSession);
            vPartitioner->PartitionMesh(vComm->GetSize());
            vPartitioner->WritePartitions(vSession, meshfile);
        }

        vComm->Block();

        meshfile = meshfile + "." + boost::lexical_cast<std::string>(vComm->GetRank());
    }

    MultiRegions::ContField2DSharedPtr Exp,Fce;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 
    NekDouble  lambda;
    NekDouble  ax,ay;
    NekDouble   st, cps = (double)CLOCKS_PER_SEC;
    // default solution type multilevel statis condensation
    MultiRegions::GlobalSysSolnType SolnType = MultiRegions::eDirectMultiLevelStaticCond;

    if((argc != 2)&&(argc != 3))
    {
        fprintf(stderr,"Usage: SteadyLinearAdvectionReaction2D  meshfile [SysSolnType]\n");
        exit(1);
    }

    //----------------------------------------------
    // Load the solver type so we can test full solve, static
    // condensation and the default multi-level statis condensation.
    if( argc == 3 )
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
        else
        {
            cerr << "SolnType not recognised" <<endl;
            exit(1);
        }

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
    // Get Advection Velocity
    ax = bcs.GetParameter("Advection_x");
    ay = bcs.GetParameter("Advection_y");
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    lambda = bcs.GetParameter("Lambda");
    cout << "            Lambda         : " << lambda << endl;
    const SpatialDomains::ExpansionMap &expansions = graph2D.GetExpansions();
    LibUtilities::BasisKey bkey0 = expansions.at(0)->m_basisKeyVector[0];
    LibUtilities::BasisKey bkey1 = expansions.at(0)->m_basisKeyVector[1];
    cout << "Solving Steady 2D LinearAdvection :"  << endl; 
    cout << "            Advection_x    : " << ax << endl; 
    cout << "            Advection_y    : " << ay << endl; 
    cout << "            Expansion      : (" << LibUtilities::BasisTypeMap[bkey0.GetBasisType()] <<","<< LibUtilities::BasisTypeMap[bkey1.GetBasisType()]  << ")" << endl;
    cout << "            No. modes      : " << bkey0.GetNumModes() << endl;
    cout << endl;
    //----------------------------------------------
   
    //----------------------------------------------
    // Define Expansion 
    int bc_val = 0;
    Exp = MemoryManager<MultiRegions::ContField2D>::
        AllocateSharedPtr(vComm,graph2D,bcs,bc_val,SolnType);
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

    Array<OneD, Array< OneD, NekDouble> > Vel(2);
    Vel[0] = Array<OneD, NekDouble> (nq,ax);
    Vel[1] = Array<OneD, NekDouble> (nq,ay);
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
    Exp->LinearAdvectionReactionSolve(Vel, Fce->GetPhys(), Exp->UpdateContCoeffs(), lambda, true);
    //----------------------------------------------
    Timing("Linear Advection Solve ..");
    
    //----------------------------------------------
    // Backward Transform Solution to get solved values 
    Exp->BwdTrans(Exp->GetContCoeffs(), Exp->UpdatePhys(), true);
    //----------------------------------------------
    
    //----------------------------------------------
    // Write solution 
    ofstream outfile("LinearAdvectionFile2D.pos");
    Exp->WriteToFile(outfile,eGmsh);
    outfile.close();

    ofstream outfile2("LinearAdvectionFile2D.dat");
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
        for(i = 0; i < nq; ++i)
        {
            fce[i] = ex_sol->Evaluate(xc0[i],xc1[i],xc2[i]);
        }
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
