#include <cstdio>
#include <cstdlib>
#include <iostream> 
#include <SpatialDomains/MeshGraph1D.h>
#include <SpatialDomains/BoundaryConditions.h>

#include <iostream>
#include <vector>
#include <string>

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
    string meshfile = "C:\\Data\\PhD\\Research\\dev\\Nektar++\\library\\Demos\\SpatialDomains\\meshdef1D.xml";
    string bcfile = "c:\\Data\\PhD\\Research\\dev\\Nektar++\\library\\Demos\\SpatialDomains\\BC1.xml";

    MeshGraph1D graph1D;
    BoundaryConditions bcs(&graph1D);

    graph1D.Read(meshfile);
    bcs.Read(bcfile);

    return 0;
}
