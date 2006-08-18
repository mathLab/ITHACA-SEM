#include <cstdio>
#include <cstdlib>
#include <iostream> 
#include "SpatialDomains/MeshGraph2D.h"

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
    MeshGraph2D graph2D;

    graph2D.Read(in);

    GeometrySharedPtr geom = graph2D.GetCompositeItem(2,2);

    return 0;
}
