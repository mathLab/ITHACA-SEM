#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <SpatialDomains/MeshGraph1D.h>
#include <SpatialDomains/Conditions.h>

#include <iostream>
#include <vector>
#include <string>

#include <boost/shared_ptr.hpp>

using namespace Nektar;
using namespace SpatialDomains;
using namespace std;

// compile using Builds/Demos/SpatialDomains -> make DEBUG=1 Graph1D

int main(int argc, char *argv[]){

    //if(argc != 2){
    //  cerr << "usage: Graph1D file" << endl;
    //  exit(1);
    //}
    LibUtilities::SessionReaderSharedPtr vSession;

    //string in(argv[argc-1]);
#ifdef PC
    string meshfile = "C:\\Data\\PhD\\Research\\dev\\Nektar++\\library\\Demos\\SpatialDomains\\meshdef1D.xml";
    string bcfile = "c:\\Data\\PhD\\Research\\dev\\Nektar++\\library\\Demos\\SpatialDomains\\BC1.xml";
#else
    string meshfile("../../../library/Demos/SpatialDomains/meshdef1D.xml");
    string bcfile("../../../library/Demos/SpatialDomains/BC1.xml");
#endif

    MeshGraph1D graph1D;
    BoundaryConditions bcs(&graph1D);

    graph1D.ReadGeometry(meshfile);
    graph1D.ReadExpansions(meshfile);
    bcs.Read(bcfile);

    BoundaryRegionCollection &boundaryRegions = bcs.GetBoundaryRegions();
    BoundaryConditionCollection &boundaryConditions = bcs.GetBoundaryConditions();

    // Region 1, v component
    BoundaryConditionShPtr bcShPtr((*boundaryConditions[1])["v"]);
    boost::shared_ptr<RobinBoundaryCondition> rbBC(boost::dynamic_pointer_cast<RobinBoundaryCondition>(bcShPtr));

    ConstForcingFunctionShPtr ffunc  = bcs.GetForcingFunction("u");
    ConstInitialConditionShPtr ic = bcs.GetInitialCondition("v");

    return 0;
}
