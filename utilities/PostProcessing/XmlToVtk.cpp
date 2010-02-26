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

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-1]);
    SpatialDomains::MeshGraph graph;
    SpatialDomains::MeshGraphSharedPtr graphShPt = graph.Read(meshfile);
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    int expdim  = graphShPt->GetMeshDimension();
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp(1);

    switch(expdim)
    {
    case 1:
        {
            SpatialDomains::MeshGraph1DSharedPtr mesh;

            if(!(mesh = boost::dynamic_pointer_cast<
                                    SpatialDomains::MeshGraph1D>(graphShPt)))
            {
                ASSERTL0(false,"Dynamics cast failed");
            }

            MultiRegions::ExpList1DSharedPtr Exp1D;
            Exp1D = MemoryManager<MultiRegions::ExpList1D>
                                                    ::AllocateSharedPtr(*mesh);
            Exp[0] = Exp1D;
        }
        break;
    case 2:
        {
            SpatialDomains::MeshGraph2DSharedPtr mesh;

            if(!(mesh = boost::dynamic_pointer_cast<
                                    SpatialDomains::MeshGraph2D>(graphShPt)))
            {
                ASSERTL0(false,"Dynamics cast failed");
            }

            MultiRegions::ExpList2DSharedPtr Exp2D;
            Exp2D = MemoryManager<MultiRegions::ExpList2D>
                                                    ::AllocateSharedPtr(*mesh);
            Exp[0] =  Exp2D;
        }
        break;
    case 3:
        {
            SpatialDomains::MeshGraph3DSharedPtr mesh;

            if(!(mesh = boost::dynamic_pointer_cast<
                                    SpatialDomains::MeshGraph3D>(graphShPt)))
            {
                ASSERTL0(false,"Dynamics cast failed");
            }

            MultiRegions::ExpList3DSharedPtr Exp3D;
            Exp3D = MemoryManager<MultiRegions::ExpList3D>
                                                    ::AllocateSharedPtr(*mesh);
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

