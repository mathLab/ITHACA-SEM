#include <cstdio>
#include <cstdlib>
#include <math.h>

#include "StdRegions/StdExpansion2D.h"
#include "StdRegions/StdNodalTriExp.h"

#include "StdRegions/StdRegions.hpp"
#include "LibUtilities/Foundations/Foundations.hpp"

using namespace Nektar;
using namespace StdRegions; 
using namespace std;

#include "StdRegions/StdNodalTriExp.h"

// This routine projects a polynomial or trigonmetric functions which 
// has energy in all mdoes of the expansions and report an error.

int main(int argc, char *argv[])
{
    int i,j;
    int order;
    int npts; 
    
    LibUtilities::PointsType    Qtype1,Qtype2;
    LibUtilities::BasisType     btype1,btype2;
    LibUtilities::PointsType     NodalType;
    StdRegions::StdNodalTriExp *E;
    Array<OneD, const NekDouble> r,s;

    if(argc != 3)
    {
        fprintf(stderr,"Usage: NodalBasis Type order\n");
        
        fprintf(stderr,"Where type is an integer value which "
                "dictates the basis as:\n");
        fprintf(stderr,"\t NodalTriElec    = 0\n");
        fprintf(stderr,"\t NodalTriFekete  = 1\n");
        exit(1);
    }
    
    btype1 = LibUtilities::eOrtho_A;
    btype2 = LibUtilities::eOrtho_B;
    Qtype1 = LibUtilities::eGaussLobattoLegendre; 
    Qtype2 = LibUtilities::eGaussLobattoLegendre; 
    
    switch(atoi(argv[argc-2]))
    {
    case 0:
        NodalType = LibUtilities::eNodalTriElec;
        break;
    case 1:
        NodalType = LibUtilities::eNodalTriFekete;
        break;
    }

    order =   atoi(argv[argc-1]);
    ASSERTL0(order > 1,"Order must be larger than 1");

    const LibUtilities::PointsKey Pkey1(order+1,Qtype1);
    const LibUtilities::PointsKey Pkey2(order+1,Qtype2);
    const LibUtilities::BasisKey  Bkey1(btype1,order,Pkey1);
    const LibUtilities::BasisKey  Bkey2(btype2,order,Pkey2);

    E = new StdRegions::StdNodalTriExp(Bkey1,Bkey2,NodalType);
    
    E->GetNodalPoints(r,s);

    //----------------------------------------------
    // Output 1D basis using only basis information. 
    npts = order*(order+1)/2;
    fprintf(stdout,"VARIABLES = x y\n");
    fprintf(stdout,"ZONE T = \"Order %d\" I=%d\n",order,npts);
    
    for(i = 0; i  < npts; ++i)
    {
        fprintf(stdout, "%16.14lf %16.14lf \n",r[i],s[i]);
    }
    //-----------------------------------------------
}
