#include <cstdio>
#include <cstdlib>
#include <iostream> 
#include "SpatialDomains/MeshGraph2D.h"
#include <SpatialDomains/BoundaryConditions.h>

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
    string in("C:/Data/PhD/Research/dev/Nektar++/library/Demos/SpatialDomains/meshdef2D.xml");
    string bcfile("c:/Data/PhD/Research/dev/Nektar++/library/Demos/SpatialDomains/BC1.xml");
    MeshGraph2D graph2D;
    BoundaryConditions bcs(&graph2D);

    graph2D.Read(in);
    bcs.Read(bcfile);
    
    return 0;
}
