#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    Array<OneD,NekDouble>  fce; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 

    if(argc != 2)
    {
        fprintf(stderr,"Usage: XmlToTecplot  meshfile\n");
        exit(1);
    }

    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-1]);
    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(meshfile);
    //----------------------------------------------

    //----------------------------------------------        
    // Define Expansion 
    graphShPt->ReadExpansions(meshfile);
    int expdim   = graphShPt->GetMeshDimension();
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp(1);

    switch(expdim)
    {
    case 1:
        {
            MultiRegions::ExpList1DSharedPtr Exp1D;
            Exp1D = MemoryManager<MultiRegions::ExpList1D>::AllocateSharedPtr(vSession,graphShPt);
            Exp[0] = Exp1D;
        }
        break;
    case 2:
        {
            MultiRegions::ExpList2DSharedPtr Exp2D;
            Exp2D = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(vSession,graphShPt);
            Exp[0] =  Exp2D;
        }
        break;
    case 3:
        ASSERTL0(false,"3D not set up");
        break;
    default:
        ASSERTL0(false,"Expansion dimension not recognised");
        break;
    }

    //-----------------------------------------------
        
    //----------------------------------------------
    // Write solution  depending on #define
    string   outfile(strtok(argv[argc-1],"."));
    outfile +=  ".dat";
    ofstream outstrm(outfile.c_str());

    Exp[0]->WriteToFile(outstrm,eTecplot,"Dummy");
    outstrm.close();
    //----------------------------------------------
    return 0;
}

