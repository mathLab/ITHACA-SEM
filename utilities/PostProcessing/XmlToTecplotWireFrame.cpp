#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>

using namespace Nektar;

void CheckTetRotation(Array<OneD, NekDouble> &xc, Array<OneD, NekDouble> &yc, 
                      Array<OneD, NekDouble> &xz, 
                      std::map<int,SpatialDomains::TetGeomSharedPtr>::iterator &tetIter,
                      int id);

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
            NekDouble x,y,z;
            string   outname(strtok(argv[argc-1],"."));
            outname += ".dat";
            FILE *fp = fopen(outname.c_str(),"w");
            
            SpatialDomains::TetGeomMap   tetgeom   = mesh->GetAllTetGeoms();
            SpatialDomains::PyrGeomMap   pyrgeom   = mesh->GetAllPyrGeoms();
            SpatialDomains::PrismGeomMap prismgeom = mesh->GetAllPrismGeoms();
            SpatialDomains::HexGeomMap   hexgeom   = mesh->GetAllHexGeoms();

            int nverts = mesh->GetNvertices();

            fprintf(fp,"Variables = x, y, z\n");
            fprintf(fp,"Zone,N=%d, E=%d,DATAPACKING=POINT,ZONETYPE=FEBRICK\n",nverts,(int)(tetgeom.size() + pyrgeom.size() + prismgeom.size() + hexgeom.size()));
            
            Array<OneD, NekDouble> xc(nverts),yc(nverts),zc(nverts);
            
            for(int i = 0; i < nverts; ++i)
            {
                mesh->GetVertex(i)->GetCoords(x,y,z);
                fprintf(fp,"%lf %lf %lf \n",x,y,z);
                xc[i] = x;
                yc[i] = y;
                zc[i] = z;
            }
            
            std::map<int,SpatialDomains::TetGeomSharedPtr>::iterator tetIter;
            int cnt = 0;
            for(tetIter = tetgeom.begin(); tetIter != tetgeom.end(); ++tetIter)
            {
                // check rotation and dump
                CheckTetRotation(xc,yc,zc,tetIter,cnt++);

                fprintf(fp,"%d %d %d %d %d %d %d %d\n",(tetIter->second)->GetVid(0)+1,
                        (tetIter->second)->GetVid(1)+1,(tetIter->second)->GetVid(2)+1,
                        (tetIter->second)->GetVid(3)+1,(tetIter->second)->GetVid(3)+1,
                        (tetIter->second)->GetVid(3)+1,(tetIter->second)->GetVid(3)+1,
                        (tetIter->second)->GetVid(3)+1);
            }


            std::map<int,SpatialDomains::PyrGeomSharedPtr>::iterator pyrIter;
            for(pyrIter = pyrgeom.begin(); pyrIter != pyrgeom.end(); ++pyrIter)
            {
                fprintf(fp,"%d %d %d %d %d %d %d %d\n",(pyrIter->second)->GetVid(0)+1,
                        (pyrIter->second)->GetVid(1)+1,(pyrIter->second)->GetVid(2)+1,
                        (pyrIter->second)->GetVid(3)+1,(pyrIter->second)->GetVid(4)+1,
                        (pyrIter->second)->GetVid(4)+1,(pyrIter->second)->GetVid(4)+1,
                        (pyrIter->second)->GetVid(5)+1);
            }


            std::map<int,SpatialDomains::PrismGeomSharedPtr>::iterator prismIter;
            for(prismIter = prismgeom.begin(); prismIter != prismgeom.end(); ++prismIter)
            {
                fprintf(fp,"%d %d %d %d %d %d %d %d\n",(prismIter->second)->GetVid(0)+1,
                        (prismIter->second)->GetVid(1)+1,(prismIter->second)->GetVid(2)+1,
                        (prismIter->second)->GetVid(3)+1,(prismIter->second)->GetVid(4)+1,
                        (prismIter->second)->GetVid(4)+1,(prismIter->second)->GetVid(5)+1,
                        (prismIter->second)->GetVid(5)+1);
            }

            std::map<int,SpatialDomains::HexGeomSharedPtr>::iterator hexIter;
            for(hexIter = hexgeom.begin(); hexIter != hexgeom.end(); ++hexIter)
            {
                fprintf(fp,"%d %d %d %d %d %d %d %d\n",(hexIter->second)->GetVid(0)+1,
                        (hexIter->second)->GetVid(1)+1,(hexIter->second)->GetVid(2)+1,
                        (hexIter->second)->GetVid(3)+1,(hexIter->second)->GetVid(4)+1,
                        (hexIter->second)->GetVid(5)+1,(hexIter->second)->GetVid(6)+1,
                        (hexIter->second)->GetVid(7)+1);
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

void CheckTetRotation(Array<OneD, NekDouble> &xc, Array<OneD, NekDouble> &yc, Array<OneD, NekDouble> &zc, std::map<int,SpatialDomains::TetGeomSharedPtr>::iterator &tetIter, int id)
{
    Ord       v[4];
    NekDouble abx,aby,abz; 
    
    v[0].x = xc[(tetIter->second)->GetVid(0)];
    v[0].y = yc[(tetIter->second)->GetVid(0)];
    v[0].z = zc[(tetIter->second)->GetVid(0)];

    v[1].x = xc[(tetIter->second)->GetVid(1)];
    v[1].y = yc[(tetIter->second)->GetVid(1)];
    v[1].z = zc[(tetIter->second)->GetVid(1)];

    v[2].x = xc[(tetIter->second)->GetVid(2)];
    v[2].y = yc[(tetIter->second)->GetVid(2)];
    v[2].z = zc[(tetIter->second)->GetVid(2)];

    v[3].x = xc[(tetIter->second)->GetVid(3)];
    v[3].y = yc[(tetIter->second)->GetVid(3)];
    v[3].z = zc[(tetIter->second)->GetVid(3)];
    
    // cross product of edge 0 and 2
    abx = (v[1].y-v[0].y)*(v[2].z-v[0].z) - 
        (v[1].z-v[0].z)*(v[2].y-v[0].y);
    aby = (v[1].z-v[0].z)*(v[2].x-v[0].x) -
        (v[1].x-v[0].x)*(v[2].z-v[0].z);
    abz = (v[1].x-v[0].x)*(v[2].y-v[0].y) -
        (v[1].y-v[0].y)*(v[2].x-v[0].x);

    // inner product of cross product with respect to edge 3 should be positive 
    if(((v[3].x-v[0].x)*abx + (v[3].y-v[0].y)*aby +
        (v[3].z-v[0].z)*abz)<0.0)
    {
        cerr << "ERROR: Element " << id + 1 << "is NOT counter-clockwise\n" << endl;
    }
    
}
