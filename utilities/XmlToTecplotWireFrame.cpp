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
        ASSERTL0(false,"1D not set up");
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
        {
            int i;
            NekDouble x,y,z;
            string   outname(strtok(argv[argc-1],"."));
            outname += ".dat";
            FILE *fp = fopen(outname.c_str(),"w");
            
            SpatialDomains::TetGeomMap   tetgeom   = mesh->GetAllTetGeoms();
            SpatialDomains::PyrGeomMap   pyrgeom   = mesh->GetAllPyrGeoms();
            SpatialDomains::PrismGeomMap prismgeom = mesh->GetAllPrismGeoms();
            SpatialDomains::HexGeomMap   hexgeom   = mesh->GetAllHexGeoms();

            std::map<int,SpatialDomains::TetGeomSharedPtr>::iterator tetIter;
            std::map<int,SpatialDomains::PyrGeomSharedPtr>::iterator pyrIter;
            std::map<int,SpatialDomains::PrismGeomSharedPtr>::iterator prismIter;
            std::map<int,SpatialDomains::HexGeomSharedPtr>::iterator hexIter;

            int nverts = mesh->GetNvertices();

            map<int, int> vertid;
            vector<int> vertord;

            int cnt = 0;
            // set up list of all vertex ids
            for(tetIter = tetgeom.begin(); tetIter != tetgeom.end(); ++tetIter)
            {
                for(i = 0; i < 4; ++i)
                {
                    int vid = (tetIter->second)->GetVid(i);
                    if(vertid.count(vid) == 0)
                    {
                        vertord.push_back(vid);
                        vertid[vid] = cnt++;
                    }
                }
            }
            
            for(pyrIter = pyrgeom.begin(); pyrIter != pyrgeom.end(); ++pyrIter)
            {
                for(i = 0; i < 5; ++i)
                {
                    int vid = (pyrIter->second)->GetVid(i);
                    if(vertid.count(vid) == 0)
                    {
                        vertord.push_back(vid);
                        vertid[vid] = cnt++;
                    }
                }
            }
            

            for(prismIter = prismgeom.begin(); prismIter != prismgeom.end(); ++prismIter)
            {
                for(i = 0; i < 5; ++i)
                {
                    int vid = (prismIter->second)->GetVid(i);
                    if(vertid.count(vid) == 0)
                    {
                        vertord.push_back(vid);
                        vertid[vid] = cnt++;
                    }
                }
            }

            for(hexIter = hexgeom.begin(); hexIter != hexgeom.end(); ++hexIter)
            {
                for(i = 0; i < 6; ++i)
                {
                    int vid = (hexIter->second)->GetVid(i);
                    if(vertid.count(vid) == 0)
                    {
                        vertord.push_back(vid);
                        vertid[vid] = cnt++;
                    }
                }
            }
            
            ASSERTL0(cnt == nverts,"Vertex count did not match");

            fprintf(fp,"Variables = x, y, z\n");
            fprintf(fp,"Zone,N=%d, E=%d,DATAPACKING=POINT,ZONETYPE=FEBRICK\n",nverts,(int)(tetgeom.size() + pyrgeom.size() + prismgeom.size() + hexgeom.size()));
            
            Array<OneD, NekDouble> xc(nverts),yc(nverts),zc(nverts);
            
            map<int, int>::iterator viter;
            for(i = 0; i < nverts; ++i)
            {
                mesh->GetVertex(vertord[i])->GetCoords(x,y,z);
                fprintf(fp,"%lf %lf %lf \n",x,y,z);
                xc[i] = x;
                yc[i] = y;
                zc[i] = z;
            }
            
            cnt = 0;
            for(tetIter = tetgeom.begin(); tetIter != tetgeom.end(); ++tetIter)
            {

                fprintf(fp,"%d %d %d %d %d %d %d %d\n",
                        vertid[(tetIter->second)->GetVid(0)]+1,
                        vertid[(tetIter->second)->GetVid(1)]+1,
                        vertid[(tetIter->second)->GetVid(2)]+1,
                        vertid[(tetIter->second)->GetVid(2)]+1,
                        vertid[(tetIter->second)->GetVid(3)]+1,
                        vertid[(tetIter->second)->GetVid(3)]+1,
                        vertid[(tetIter->second)->GetVid(3)]+1,
                        vertid[(tetIter->second)->GetVid(3)]+1);
            }


            for(pyrIter = pyrgeom.begin(); pyrIter != pyrgeom.end(); ++pyrIter)
            {
                fprintf(fp,"%d %d %d %d %d %d %d %d\n",
                        vertid[(pyrIter->second)->GetVid(0)]+1,
                        vertid[(pyrIter->second)->GetVid(1)]+1,
                        vertid[(pyrIter->second)->GetVid(3)]+1,
                        vertid[(pyrIter->second)->GetVid(2)]+1,
                        vertid[(pyrIter->second)->GetVid(4)]+1,
                        vertid[(pyrIter->second)->GetVid(4)]+1,
                        vertid[(pyrIter->second)->GetVid(4)]+1,
                        vertid[(pyrIter->second)->GetVid(4)]+1);
            }

            for(prismIter = prismgeom.begin(); prismIter != prismgeom.end(); ++prismIter)
            {
                fprintf(fp,"%d %d %d %d %d %d %d %d\n",
                        vertid[(prismIter->second)->GetVid(0)]+1,
                        vertid[(prismIter->second)->GetVid(1)]+1,
                        vertid[(prismIter->second)->GetVid(3)]+1,
                        vertid[(prismIter->second)->GetVid(2)]+1,
                        vertid[(prismIter->second)->GetVid(4)]+1,
                        vertid[(prismIter->second)->GetVid(4)]+1,
                        vertid[(prismIter->second)->GetVid(5)]+1,
                        vertid[(prismIter->second)->GetVid(5)]+1);
            }

            for(hexIter = hexgeom.begin(); hexIter != hexgeom.end(); ++hexIter)
            {
                fprintf(fp,"%d %d %d %d %d %d %d %d\n",
                        vertid[(hexIter->second)->GetVid(0)]+1,
                        vertid[(hexIter->second)->GetVid(1)]+1,
                        vertid[(hexIter->second)->GetVid(3)]+1,
                        vertid[(hexIter->second)->GetVid(2)]+1,
                        vertid[(hexIter->second)->GetVid(4)]+1,
                        vertid[(hexIter->second)->GetVid(5)]+1,
                        vertid[(hexIter->second)->GetVid(7)]+1,
                        vertid[(hexIter->second)->GetVid(6)]+1);
            }

        }            

        break;
    default:
        ASSERTL0(false,"Expansion dimension not recognised");
        break;
    }

    //-----------------------------------------------
        
    return 0;
}

class Ord
{
public:
    double x;
    double y;
    double z;
};

