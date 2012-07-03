#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    Array<OneD,NekDouble>  fce; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 

    if(argc != 2)
    {
        fprintf(stderr,"Usage: XmlToTecplotWireFrame  meshfile\n");
        exit(1);
    }

    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-1]);
    SpatialDomains::MeshGraphSharedPtr mesh  = SpatialDomains::MeshGraph::Read(vSession);
    //----------------------------------------------

    //----------------------------------------------        
    // Define Expansion 
    int expdim   = mesh->GetMeshDimension();

    switch(expdim)
    {
    case 1:
        ASSERTL0(false,"3D not set up");
        break;
    case 2:
        {
            NekDouble x,y,z;
            string   outname(strtok(argv[argc-1],"."));
            outname += ".dat";
            FILE *fp = fopen(outname.c_str(),"w");
            
            SpatialDomains::TriGeomMap trigeom = mesh->GetAllTriGeoms();
            SpatialDomains::QuadGeomMap quadgeom = mesh->GetAllQuadGeoms();

            int nverts = mesh->GetNvertices();

            fprintf(fp,"Variables = x, y\n");
            fprintf(fp,"Zone,N=%d, E=%d,DATAPACKING=POINT,ZONETYPE=FEQUADRILATERAL\n",nverts,(int)(trigeom.size()+quadgeom.size()));
            
            for(int i = 0; i < nverts; ++i)
            {
                mesh->GetVertex(i)->GetCoords(x,y,z);
                fprintf(fp,"%lf %lf\n",x,y);
            }
            
            std::map<int,SpatialDomains::TriGeomSharedPtr>::iterator triIter;
            for(triIter = trigeom.begin(); triIter != trigeom.end(); ++triIter)
            {
                fprintf(fp,"%d %d %d %d\n",(triIter->second)->GetVid(0)+1,(triIter->second)->GetVid(1)+1,(triIter->second)->GetVid(2)+1,(triIter->second)->GetVid(2)+1);
            }
            
            std::map<int,SpatialDomains::QuadGeomSharedPtr>::iterator quadIter;
            for(quadIter = quadgeom.begin(); quadIter != quadgeom.end(); ++quadIter)
            {
                fprintf(fp,"%d %d %d %d\n",(quadIter->second)->GetVid(0)+1,(quadIter->second)->GetVid(1)+1,(quadIter->second)->GetVid(2)+1,(quadIter->second)->GetVid(3)+1);
            }
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
        
    return 0;
}

