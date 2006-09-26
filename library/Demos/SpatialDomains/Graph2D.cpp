#include <cstdio>
#include <cstdlib>
#include <iostream> 
#include "SpatialDomains/MeshGraph2D.h"
#include "SpatialDomains/Domain.h"

using namespace Nektar;
using namespace SpatialDomains; 
using namespace std;

// compile using Builds/Demos/SpatialDomains -> make DEBUG=1 Graph1D

int main(int argc, char *argv[]){

    //if(argc != 2){
    //  cerr << "usage: Graph1D file" << endl;
    //  exit(1);
    //}

    //string in(argv[argc-1]);
    string in("C:\\Data\\PhD\\Research\\dev\\Nektar++\\libs\\Demos\\SpatialDomains\\meshdef2D.xml");
    string domainfile("C:\\Data\\PhD\\Research\\dev\\Nektar++\\libs\\Demos\\SpatialDomains\\domain.xml");
    MeshGraph2D graph2D;
    Domain domain2D(&graph2D);

    graph2D.Read(in);
    domain2D.Read(domainfile);

    CompositeVector compVector = domain2D.GetDomain();
    BoundaryVector boundVector = domain2D.GetBoundaries();
    
    return 0;
}
