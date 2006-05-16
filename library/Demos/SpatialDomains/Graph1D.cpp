#include <cstdio>
#include <cstdlib>
#include <iostream> 
#include "SpatialDomains/MeshGraph1D.h"

using namespace Nektar;
using namespace SpatialDomains; 
using namespace std;

// compile using Builds/Demos/SpatialDomains -> make DEBUG=1 Graph1D

int main(int argc, char *argv[]){

    if(argc != 2){
      cerr << "usage: Graph1D file" << endl;
      exit(1);
    }

    string in(argv[argc-1]);
    MeshGraph1D graph1D;

    graph1D.Read(in);

    return 0;
}
