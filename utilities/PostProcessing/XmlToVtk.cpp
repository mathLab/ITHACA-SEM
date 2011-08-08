#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    int i;

    if(argc != 2)
    {
        fprintf(stderr,"Usage: XmlToVtk  meshfile\n");
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
    int expdim  = graphShPt->GetMeshDimension();
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp(1);

    switch(expdim)
    {
    case 1:
        {
            MultiRegions::ExpList1DSharedPtr Exp1D;
            Exp1D = MemoryManager<MultiRegions::ExpList1D>
                                                    ::AllocateSharedPtr(vSession,graphShPt);
            Exp[0] = Exp1D;
        }
        break;
    case 2:
        {
            MultiRegions::ExpList2DSharedPtr Exp2D;
            Exp2D = MemoryManager<MultiRegions::ExpList2D>
                                                    ::AllocateSharedPtr(vSession,graphShPt);
            Exp[0] =  Exp2D;
        }
        break;
    case 3:
        {
            MultiRegions::ExpList3DSharedPtr Exp3D;
            Exp3D = MemoryManager<MultiRegions::ExpList3D>
                                                    ::AllocateSharedPtr(vSession,graphShPt);
            Exp[0] =  Exp3D;
        }
        break;
    default:
        ASSERTL0(false,"Expansion dimension not recognised");
        break;
    }
    //----------------------------------------------

    //----------------------------------------------
    // Write out VTK file.
    string   outname(strtok(argv[argc-1],"."));
    outname += ".vtu";
    ofstream outfile(outname.c_str());

    Exp[0]->WriteVtkHeader(outfile);
    // For each field write header and footer, since there is no field data.
    for(i = 0; i < Exp[0]->GetExpSize(); ++i)
    {
        Exp[0]->WriteVtkPieceHeader(outfile,i);
        Exp[0]->WriteVtkPieceFooter(outfile,i);
    }
    Exp[0]->WriteVtkFooter(outfile);
    //----------------------------------------------

    return 0;
}

