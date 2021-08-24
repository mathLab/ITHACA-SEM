#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>

using namespace std;
using namespace Nektar;

bool CheckTetRotation(Array<OneD, NekDouble> &xc, Array<OneD, NekDouble> &yc, 
                      Array<OneD, NekDouble> &xz, 
                      std::map<int,SpatialDomains::TetGeomSharedPtr>::iterator &tetIter,
                      int id);

int main(int argc, char *argv[])
{
    Array<OneD,NekDouble>  fce; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 

    if(argc != 2)
    {
        fprintf(stderr,"Usage: CheckXmlFile  meshfile.xml\n");
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

            Array<OneD, NekDouble> xc(nverts),yc(nverts),zc(nverts);

            for(int i = 0; i < nverts; ++i)
            {
                mesh->GetVertex(i)->GetCoords(x,y,z);
                xc[i] = x;
                yc[i] = y;
                zc[i] = z;
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
            
            SpatialDomains::TetGeomMap   tetgeom   = mesh->GetAllTetGeoms();
            SpatialDomains::PyrGeomMap   pyrgeom   = mesh->GetAllPyrGeoms();
            SpatialDomains::PrismGeomMap prismgeom = mesh->GetAllPrismGeoms();
            SpatialDomains::HexGeomMap   hexgeom   = mesh->GetAllHexGeoms();

            int nverts = mesh->GetNvertices();

            Array<OneD, NekDouble> xc(nverts),yc(nverts),zc(nverts);
            
            for(int i = 0; i < nverts; ++i)
            {
                // Check element roation 
                mesh->GetVertex(i)->GetCoords(x,y,z);
                xc[i] = x;
                yc[i] = y;
                zc[i] = z;
                
            }

            for (int i = 0; i < nverts; ++i)
            {
                for (int j = i+1; j < nverts; ++j)
                {
                    if ((xc[i]-xc[j])*(xc[i]-xc[j]) +
                        (yc[i]-yc[j])*(yc[i]-yc[j]) +
                        (zc[i]-zc[j])*(zc[i]-zc[j]) < 1e-10)
                    {
                        cout << "Duplicate vertices: " << i << " " << j << endl;
                    }
                }
            }
            
            std::map<int,SpatialDomains::TetGeomSharedPtr>::iterator tetIter;
            int cnt = 0;
            bool NoRotateIssues = true;
            bool NoOrientationIssues = true;
            for(tetIter = tetgeom.begin(); tetIter != tetgeom.end(); ++tetIter)
            {
                // check rotation and dump
                NoOrientationIssues = CheckTetRotation(xc,yc,zc,tetIter,cnt++);

                // Check face rotation 
                if((tetIter->second)->GetFace(0)->GetVid(2) != (tetIter->second)->GetVid(2))
                {
                    cout << "ERROR: Face " << tetIter->second->GetFid(0) << " (vert "<< (tetIter->second)->GetFace(0)->GetVid(2) << ") is not aligned with base vertex of Tet " << (tetIter->second)->GetGlobalID() << " (vert "<< (tetIter->second)->GetVid(2) <<")" << endl;
                    NoRotateIssues = false; 
                }

                for(int i = 1; i < 4; ++i)

                {
                    if((tetIter->second)->GetFace(i)->GetVid(2) != (tetIter->second)->GetVid(3))
                    {
                        cout << "ERROR: Face " << tetIter->second->GetFid(i) << " is not aligned with top Vertex of Tet " << (tetIter->second)->GetGlobalID() << endl;
                        NoRotateIssues = false; 
                    }
                }
                
            }
            if(NoOrientationIssues)
            {
                cout << "All Tet have correct ordering for anticlockwise rotation" << endl;
            }

            if(NoRotateIssues)
            {
                cout << "All Tet faces are correctly aligned" << endl;
            }


            std::map<int,SpatialDomains::PyrGeomSharedPtr>::iterator pyrIter;
            for(pyrIter = pyrgeom.begin(); pyrIter != pyrgeom.end(); ++pyrIter)
            {
                // Put pyramid checks in here 
            }


            // Put prism checks in here
            std::map<int,SpatialDomains::PrismGeomSharedPtr>::iterator Iter;
            NoRotateIssues = true;
            NoOrientationIssues = true;
            for(Iter = prismgeom.begin(); Iter != prismgeom.end(); ++Iter)
            {
                
                // Check face rotation 
                if((Iter->second)->GetFace(1)->GetVid(2) != (Iter->second)->GetVid(4))
                {
                    cout << "ERROR: Face " << Iter->second->GetFid(1) << " (vert "<< (Iter->second)->GetFace(1)->GetVid(2) << ") not aligned to face 1 singular vert of Prism " << (Iter->second)->GetGlobalID() << " (vert "<< (Iter->second)->GetVid(4) <<")" << endl;
                    NoRotateIssues = false; 
                }
                
                
                // Check face rotation 
                if((Iter->second)->GetFace(3)->GetVid(2) != (Iter->second)->GetVid(5))
                {
                    cout << "ERROR: Face " << Iter->second->GetFid(3) << " (vert "<< (Iter->second)->GetFace(3)->GetVid(2) << ") not aligned to face 3 singular vert of Prism " << (Iter->second)->GetGlobalID() << " (vert "<< (Iter->second)->GetVid(5) <<")" << endl;
                    NoRotateIssues = false; 
                }
                
            }
                
            if(NoRotateIssues)
            {
                cout << "All Prism Tri faces are correctly aligned" << endl;
            }
            
            std::map<int,SpatialDomains::HexGeomSharedPtr>::iterator hexIter;
            for(hexIter = hexgeom.begin(); hexIter != hexgeom.end(); ++hexIter)
            {
                // PUt Hex checks in here
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

bool CheckTetRotation(Array<OneD, NekDouble> &xc, Array<OneD, NekDouble> &yc, Array<OneD, NekDouble> &zc, std::map<int,SpatialDomains::TetGeomSharedPtr>::iterator &tetIter, int id)
{
    bool      RotationOK = true;
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
        RotationOK = false;
    }
    return RotationOK;
}
