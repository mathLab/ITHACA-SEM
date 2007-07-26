#include <cstdio>
#include <cstdlib>
#include <iostream> 
#include "SpatialDomains/MeshGraph2D.h"
#include <SpatialDomains/BoundaryConditions.h>

#include <boost/shared_ptr.hpp>

using namespace Nektar;
using namespace SpatialDomains; 
using namespace std;

// compile using Builds/Demos/SpatialDomains -> make DEBUG=1 Graph1D

int main(int argc, char *argv[]){

    //if(argc != 2){
    //  cerr << "usage: Graph2D file" << endl;
    //  exit(1);
    //}

    //string in(argv[argc-1]);
    // If we all have the same relative structure, these should work for everyone.
#if 0 
    string in("../../../library/Demos/SpatialDomains/meshdef2D.xml");
    string bcfile("../../../library/Demos/SpatialDomains/BC1.xml");
#else 
    string in("C:/Data/PhD/Research/dev/Nektar++/library/Demos/SpatialDomains/meshdef2D.xml");
    string bcfile("c:/Data/PhD/Research/dev/Nektar++/library/Demos/SpatialDomains/BC1.xml");
#endif

    MeshGraph2D graph2D;
    BoundaryConditions bcs(&graph2D);

    graph2D.Read(in);
    bcs.Read(bcfile);

    ConstForcingFunctionShPtr ffunc  = bcs.GetForcingFunction("u");
    NekDouble val = ffunc->Evaluate(8.0);
    
    ConstInitialConditionShPtr ic = bcs.GetInitialCondition("u");
    val = ic->Evaluate();

    NekDouble tolerance = bcs.GetParameter("Tolerance");

    BoundaryConditionCollection &bConditions = bcs.GetBoundaryConditions();

    BoundaryConditionShPtr bcShPtr = (*bConditions[0])["u"];

    if (bcShPtr->GetBoundaryConditionType() == eDirichlet)
    {
        DirichletBCShPtr dirichletBCShPtr = boost::static_pointer_cast<DirichletBoundaryCondition>(bcShPtr);
        NekDouble val = dirichletBCShPtr->m_DirichletCondition.Evaluate();
    }

    ConstExactSolutionShPtr es = bcs.GetExactSolution("u");
    val = es->Evaluate(2.0, 4.0);
    
    return 0;
}
