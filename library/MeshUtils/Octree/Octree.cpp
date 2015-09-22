////////////////////////////////////////////////////////////////////////////////
//
//  File: octree.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <limits>

#include <LibUtilities/BasicUtils/Progressbar.hpp>

#include <MeshUtils/Octree/Octree.h>
#include <LibUtilities/CADSystem/CADSurf.h>

using namespace std;
namespace Nektar{
namespace MeshUtils {

NekDouble Octree::Query(Array<OneD, NekDouble> loc)
{
    //starting at master octant 0 move through succsesive octants which contain
    //the point loc until a leaf is found
    int n = 0;
    int quad;

    bool found=false;

    while(!found)
    {
        Array<OneD, NekDouble> octloc = OctantList[n]->GetLoc();

        if(loc[0]>=octloc[0] && loc[1]>octloc[1] && loc[2]>octloc[2])
        {
            quad=0;
        }
        else if(loc[0]>octloc[0] && loc[1]<=octloc[1] && loc[2]>octloc[2])
        {
            quad=1;
        }
        else if(loc[0]<=octloc[0] && loc[1]<=octloc[1] && loc[2]>octloc[2])
        {
            quad=2;
        }
        else if(loc[0]<=octloc[0] && loc[1]>octloc[1] && loc[2]>octloc[2])
        {
            quad=3;
        }
        else if(loc[0]>octloc[0] && loc[1]>octloc[1] && loc[2]<=octloc[2])
        {
            quad=4;
        }
        else if(loc[0]>octloc[0] && loc[1]<=octloc[1] && loc[2]<=octloc[2])
        {
            quad=5;
        }
        else if(loc[0]<=octloc[0] && loc[1]<=octloc[1] && loc[2]<=octloc[2])
        {
            quad=6;
        }
        else if(loc[0]<=octloc[0] && loc[1]>octloc[1] && loc[2]<=octloc[2])
        {
            quad=7;
        }
        else
        {
            ASSERTL0(false,"Cannot locate quadrant");
        }

        n=OctantList[n]->GetChild(quad);

        if(OctantList[n]->GetLeaf())
        {
            found=true;
        }
    }
    return OctantList[n]->GetDelta();
}

void Octree::Build(const NekDouble min, const NekDouble max,
                   const NekDouble eps)
{
    m_minDelta = min;
    m_maxDelta = max;
    m_eps = eps;

    BoundingBox = m_cad->GetBoundingBox();

    if(m_verbose)
        cout << endl << "Octree system" << endl;

    //build curvature samples
    CompileCuravturePointList();

    if(m_verbose)
        cout << "\tCurvature samples: " << m_cpList.size() << endl
        << "\tInitial subdivision: ";

    NekDouble maxdim = (BoundingBox[1]-BoundingBox[0])/2 >
                           (BoundingBox[3]-BoundingBox[2])/2 ?
                           (BoundingBox[1]-BoundingBox[0])/2 :
                           (BoundingBox[3]-BoundingBox[2])/2;
    maxdim = maxdim > (BoundingBox[5]-BoundingBox[4])/2 ?
                        maxdim : (BoundingBox[5]-BoundingBox[4])/2;

    //make master octant based on the bounding box of the domain
    OctantSharedPtr newOctant =
    MemoryManager<Octant>::AllocateSharedPtr
    ((BoundingBox[1]+BoundingBox[0])/2,
     (BoundingBox[3]+BoundingBox[2])/2,
     (BoundingBox[5]+BoundingBox[4])/2, maxdim, -1, 0, m_cpList);

    OctantList.push_back(newOctant);

    m_totNotDividing=0;

    if(OctantList[0]->GetDivide())
    {
        OctantList[0]->SetLeaf(false);
        InitialSubDivide(0);
    }

    if(m_verbose)
    {
        int ct=0;
        int maxLevel=0;

        for(int i = 0; i < OctantList.size(); i++)
        {
            if(OctantList[i]->GetLeaf())
                ct++;
            if(OctantList[i]->GetLevel()>maxLevel)
                maxLevel=OctantList[i]->GetLevel();
        }

        cout << "No. octant leaves: " << ct <<
        ", Max octree level: " << maxLevel << endl;
    }

    for(int i = 0; i < OctantList.size(); i++)
    {
        if(m_verbose)
        {
            LibUtilities::PrintProgressbar(i, OctantList.size(),
                                           "\tDetermining neigbours");
        }
        if(OctantList[i]->GetLeaf())
        {
            OctantList[i]->CreateNeighbourList(OctantList);
        }
    }
    if(m_verbose)
        cout << endl;

    SubDivideByLevel();


    if(m_verbose)
    {
        int ct=0;
        for(int i = 0; i < OctantList.size(); i++)
        {
            if(OctantList[i]->GetLeaf()){ct++;}
        }
        cout << "\tNew Stats: ";
        cout << "No. octant leaves: " << ct << endl;
        cout << "Smoothing across the geometry surface" << endl;
    }

    SmoothSurfaceOctants();

    if(m_verbose)
        cout << "\tPropagating spacing out to domain boundary" << endl;

    PropagateDomain();

    if(m_verbose)
        cout << "\tRecersively ensuring smoothness between all nodes" << endl;

    SmoothAllOctants();

    for(int i = 0; i < OctantList.size(); i++)
    {
        if(OctantList[i]->GetLeaf())
        {
            ASSERTL0(!(OctantList[i]->GetDelta()<0.0),
                     "Error in initial octree construction");
        }
    }

    if(m_verbose)
    {
        //int elem=CountElemt();
        //cout << endl<< "\tPredicted mesh: " << elem << " elements" << endl;
    }
}

int Octree::CountElemt()
{
    //by considering the volume of a tet evaluate the number of elements in the
    //mesh

    Array<OneD, NekDouble> box = m_cad->GetBoundingBox();
    NekDouble total=0.0;

    for(int i = 0 ; i < OctantList.size(); i++)
    {
        if(OctantList[i]->GetLeaf() && OctantList[i]->HasPoints())
        {
            NekDouble xmin, xmax, ymin, ymax, zmin, zmax;
            xmin = max(box[0], OctantList[i]->FX(-1.0));
            xmax = min(box[1], OctantList[i]->FX(+1.0));
            ymin = max(box[2], OctantList[i]->FY(-1.0));
            ymax = min(box[3], OctantList[i]->FY(+1.0));
            zmin = max(box[4], OctantList[i]->FZ(-1.0));
            zmax = min(box[5], OctantList[i]->FZ(+1.0));
            ASSERTL0(xmin < xmax, "error");
            ASSERTL0(ymin < ymax, "error");
            ASSERTL0(zmin < zmax, "error");
            NekDouble voloverlap = (xmax - xmin)*(ymax - ymin)*(zmax - zmin);

            NekDouble volumeTet = OctantList[i]->GetDelta()*
                                  OctantList[i]->GetDelta()*
                                  OctantList[i]->GetDelta()/6.0/sqrt(2.0);

            total += voloverlap/volumeTet;
        }
        else if(OctantList[i]->GetLeaf())
        {
            Array<OneD, NekDouble> loc = OctantList[i]->GetLoc();
            if(m_cad->InsideShape(loc))
            {
                NekDouble volumeTet = OctantList[i]->GetDelta()*
                OctantList[i]->GetDelta()*
                OctantList[i]->GetDelta()/6.0/sqrt(2.0);

                NekDouble volumeOct = OctantList[i]->DX()*OctantList[i]->DX()*
                OctantList[i]->DX()*8.0;

                total += volumeOct/volumeTet;
            }
        }
    }

    return int(total);
}

void Octree::SmoothAllOctants()
{
    //until no more changes occur smooth the mesh specification between all
    //octants
    int ct = 0;

    do
    {
        ct=0;
        for(int i = 0; i < OctantList.size(); i++)
        {
            if(OctantList[i]->GetLeaf())
            {
                vector<int> checkID;
                vector<int> nList = OctantList[i]->GetNeighbourList();

                for(int j = 0; j < nList.size(); j++)
                {
                    if(OctantList[nList[j]]->GetLeaf() == false) //this should not happen but does, this is a bit hacky but fixes it
                    {
                        continue;
                    }
                    if(OctantList[nList[j]]->GetDelta() <
                       OctantList[i]->GetDelta() &&
                       ddx(i, nList[j])>0.25)
                    {
                        checkID.push_back(nList[j]);
                    }
                }

                //for each neighbour listed in check_id, figure out the smoothed delta, and asign the miminum of these to nodes[i].GetDelta()
                if(checkID.size() > 0)
                {
                    NekDouble deltaSM = numeric_limits<double>::max();
                    for(int j = 0; j < checkID.size(); j++)
                    {
                        NekDouble r =
                          OctantList[i]->Distance(OctantList[checkID[j]]);

                        if(0.24*r +
                           OctantList[checkID[j]]->GetDelta() < deltaSM)
                        {
                            deltaSM = 0.24*r +
                            OctantList[checkID[j]]->GetDelta();
                        }
                    }
                    OctantList[i]->SetDelta(deltaSM);

                    ct+=1;
                }
            }
        }

    }while(ct>0);
}

void Octree::PropagateDomain()
{
    //until all octants know their delta specifcation and orientaion
    //look over all octants and if their neighours know either their orientation
    //or specifcation calculate one for this octant
    int ct=0;

    do
    {
        ct=0;
        for(int i = 0; i<OctantList.size(); i++)
        {
            if(OctantList[i]->GetLeaf() && !OctantList[i]->GetDeltaKnown())
            { //if it is leaf, has no points and delta has not been asigned
                ct+=1;
                vector<int> knownID;
                vector<int> nList = OctantList[i]->GetNeighbourList();

                for(int j = 0; j<nList.size(); j++)
                {
                    if(OctantList[nList[j]]->GetDeltaKnown())
                    {
                        knownID.push_back(nList[j]);
                    }
                }//create list of neighbours where delta is known.
                if(knownID.size() > 0)
                {
                    vector<NekDouble> deltaPrime;
                    for(int j = 0; j < knownID.size(); j++)
                    {
                        NekDouble r =
                         OctantList[i]->Distance(OctantList[knownID[j]]);

                        if(0.24*r +
                           OctantList[knownID[j]]->GetDelta() < m_maxDelta)
                        {
                            deltaPrime.push_back(0.24*r +
                                                 OctantList[knownID[j]]->
                                                 GetDelta());
                        }
                        else
                        {
                            deltaPrime.push_back(m_maxDelta);
                        }
                    }
                    NekDouble min=numeric_limits<double>::max();
                    for(int j = 0; j < deltaPrime.size(); j++)
                    {
                        if(deltaPrime[j] < min)
                        {
                            min=deltaPrime[j];
                        }
                    }
                    OctantList[i]->SetDelta(min);
                    ASSERTL0(!(min<m_minDelta),
                            "Delta assignment less than min delta");
                    deltaPrime.clear();
                }
                knownID.clear();
            }
        }

    }while(ct>0);

    for(int i = 0; i<OctantList.size(); i++)
    {
        if(OctantList[i]->GetLeaf())
        {
            ASSERTL0(OctantList[i]->GetDelta() > 0.0,
                            "leaf delta less than zero");
        }
    }
}

void Octree::SmoothSurfaceOctants()
{
    //for all the octants which are surface containing and know their delta
    //specification, look over all neighbours and ensure the specification
    //between them is smooth
    int ct = 0;

    do
    {
        ct=0;

        for(int i = 0; i < OctantList.size(); i++)
        {
            if(OctantList[i]->GetLeaf() && OctantList[i]->GetDeltaKnown())
            {
                vector<int> checkID;
                vector<int> nList = OctantList[i]->GetNeighbourList();

                for(int j = 0; j < nList.size(); j++)
                {
                    if(OctantList[nList[j]]->GetDeltaKnown()
                       && OctantList[nList[j]]->GetDelta()
                       < OctantList[i]->GetDelta() &&
                       ddx(i,nList[j]) > 0.075)
                    {
                        checkID.push_back(nList[j]);
                    }
                }

                //for each neighbour listed in check_id, figure out the smoothed
                // delta, and asign the miminum of these to nodes[i].GetDelta()
                if(checkID.size() > 0)
                {
                    NekDouble deltaSM = numeric_limits<double>::max();
                    for(int j = 0; j < checkID.size(); j++)
                    {
                        NekDouble r =
                            OctantList[i]->Distance(OctantList[checkID[j]]);

                        if(0.074*r +
                           OctantList[checkID[j]]->GetDelta() < deltaSM)
                        {
                            deltaSM = 0.074*r +
                            OctantList[checkID[j]]->GetDelta();
                        }
                    }
                    OctantList[i]->SetDelta(deltaSM);
                    ASSERTL0(!(deltaSM<m_minDelta),
                                    "Delta assignment less than min delta");
                    ct+=1;
                }
            }
        }

    }while(ct>0);
}

NekDouble Octree::ddx(int i, int j)
{
    NekDouble r = OctantList[i]->Distance(OctantList[j]);
    return fabs(OctantList[i]->GetDelta()-OctantList[j]->GetDelta())/r;
}

void Octree::SubDivideByLevel()
{
    //until all subdivision ceases, evaluate each octant and subdivide if
    //the neigbour levels are not smooth.
    int ct=0;
    int imax=0;

    do
    {
        ct=0;
        int run = OctantList.size();
        for(int j = 0; j < run; j++)
        {
            if(OctantList[j]->GetLeaf())
            {
                vector<int> nList = OctantList[j]->GetNeighbourList();

                for(int k = 0; k < nList.size(); k++)
                {
                    if(OctantList[nList[k]]->GetLevel() -
                       OctantList[j]->GetLevel() > 1)
                    {
                        ct+=1;
                        if(j>imax)
                            imax = j;
                        SubDivideLevel(j);
                        break;
                    }
                }
                if(ct>0)
                    break;

                nList.clear();
            }
        }

        if(m_verbose)
        {
            LibUtilities::PrintProgressbar(imax, OctantList.size(),
                                      "\tSubdivide by level");
        }

    }while(ct>0);
    if(m_verbose)
        cout <<endl;
}

void Octree::SubDivideLevel(int parent)
{
    //create 8 child octants in turn for octant parent
    //after creation, re-evaluate all nessercary neighbour lists

    OctantList[parent]->SetLeaf(false);

    Array<OneD, int> children(8);

    Array<OneD, NekDouble> parentloc = OctantList[parent]->GetLoc();

    for(int i = 0; i < 8; i++)
    {
        float pmx,pmy,pmz;
        if(i<4)
        {
            pmz=+1.0;
            if(i<2)
            {
                pmx=+1.0;
            }
            else
            {
                pmx=-1.0;
            }
            if(i==0||i==3)
            {
                pmy=+1.0;
            }
            else
            {
                pmy=-1.0;
            }
        }
        else
        {
            pmz=-1.0;
            if(i<6)
            {
                pmx=+1.0;
            }
            else
            {
                pmx=-1.0;
            }
            if(i==4||i==7)
            {
                pmy=+1.0;
            }
            else
            {
                pmy=-1.0;
            }
        }

        OctantSharedPtr newOctant =
        MemoryManager<Octant>::AllocateSharedPtr
        (parentloc[0] + pmx*OctantList[parent]->DX()/2,
         parentloc[1] + pmy*OctantList[parent]->DX()/2,
         parentloc[2] + pmz*OctantList[parent]->DX()/2,
         OctantList[parent]->DX()/2,  parent, OctantList[parent]->GetLevel()+1,
         OctantList[parent]->GetCPList());

        OctantList.push_back(newOctant);
        children[i]=OctantList.size()-1;
    }

    OctantList[parent]->SetChildren(children);

    for(int i = 0; i < 8; i++)
    {
        OctantList[children[i]]->CreateNeighbourList(OctantList);
    }

    //need to revaluate the neighbour list of all the neighbours of the parent
    vector<int> nList = OctantList[parent]->GetNeighbourList();
    for(int i = 0; i < nList.size(); i++)
    {
        OctantList[nList[i]]->CreateNeighbourList(OctantList);
    }
    nList.clear();
}


void Octree::InitialSubDivide(int parent)
{
    //in turn, create 8 child octants for octant parent
    //if that child also needs sub dividing, call this function recursively
    Array<OneD, int> children(8);

    Array<OneD, NekDouble> parentloc = OctantList[parent]->GetLoc();

    for(int i = 0; i < 8; i++)
    {
        float pmx,pmy,pmz;
        if(i<4)
        {
            pmz=+1.0;
            if(i<2)
            {
                pmx=+1.0;
            }
            else
            {
                pmx=-1.0;
            }
            if(i==0||i==3)
            {
                pmy=+1.0;
            }
            else
            {
                pmy=-1.0;
            }
        }
        else
        {
            pmz=-1.0;
            if(i<6)
            {
                pmx=+1.0;
            }
            else
            {
                pmx=-1.0;
            }
            if(i==4||i==7)
            {
                pmy=+1.0;
            }
            else
            {
                pmy=-1.0;
            }
        }

        OctantSharedPtr newOctant =
        MemoryManager<Octant>::AllocateSharedPtr
        (parentloc[0] + pmx*OctantList[parent]->DX()/2.0,
         parentloc[1] + pmy*OctantList[parent]->DX()/2.0,
         parentloc[2] + pmz*OctantList[parent]->DX()/2.0,
         OctantList[parent]->DX()/2.0, parent,
         OctantList[parent]->GetLevel() + 1,
         OctantList[parent]->GetCPList());

        OctantList.push_back(newOctant);
        children[i]=OctantList.size()-1;

        if(OctantList[children[i]]->GetDivide())
        {
            if(OctantList[children[i]]->DX()/2.0 > m_minDelta)
            {
                OctantList[children[i]]->SetLeaf(false);
                InitialSubDivide(children[i]);
            }
        }

    }

    OctantList[parent]->SetChildren(children);

}


void Octree::CompileCuravturePointList()
{

    for(int i = 1; i <= m_cad->GetNumSurf(); i++)
    {
        LibUtilities::CADSurfSharedPtr surf = m_cad->GetSurf(i);
        Array<OneD, NekDouble> bounds = surf->GetBounds();

        //to figure out the amount of curvature sampling to be conducted on
        //each parameter plane the surface is first sampled with a 40x40 grid
        //the real space lengths of this grid are analysed to find the largest
        //strecthing in the u and v directions
        //this stretching this then cosnidered with the mindelta user input
        //to find a number of sampling points in each direction which
        //enures that in the final octree each surface octant will have at least
        //one sample point within its volume.
        //the 40x40 grid is used to ensure each surface has a minimum of 40x40
        //samples.
        NekDouble du = (bounds[1]-bounds[0])/(40-1);
        NekDouble dv = (bounds[3]-bounds[2])/(40-1);

        NekDouble DeltaU = 0.0;
        NekDouble DeltaV = 0.0;

        Array<TwoD, Array<OneD, NekDouble> > samplepoints(40,40);

        for(int j = 0; j < 40; j++)
        {
            for(int k = 0; k < 40; k++)
            {
                Array<OneD, NekDouble> uv(2);
                uv[0] = k*du + bounds[0];
                uv[1] = j*dv + bounds[2];
                if(j==40-1) uv[1]=bounds[3];
                if(k==40-1) uv[0]=bounds[1];
                samplepoints[k][j] = surf->P(uv);
            }
        }

        for(int j = 0; j < 40-1; j++)
        {
            for(int k = 0; k < 40-1; k++)
            {
                NekDouble deltau = sqrt((samplepoints[k][j][0]-
                                         samplepoints[k+1][j][0])*
                                        (samplepoints[k][j][0]-
                                         samplepoints[k+1][j][0])+
                                        (samplepoints[k][j][1]-
                                         samplepoints[k+1][j][1])*
                                        (samplepoints[k][j][1]-
                                         samplepoints[k+1][j][1])+
                                        (samplepoints[k][j][2]-
                                         samplepoints[k+1][j][2])*
                                        (samplepoints[k][j][2]-
                                         samplepoints[k+1][j][2]));
                NekDouble deltav = sqrt((samplepoints[k][j][0]-
                                         samplepoints[k][j+1][0])*
                                        (samplepoints[k][j][0]-
                                         samplepoints[k][j+1][0])+
                                        (samplepoints[k][j][1]-
                                         samplepoints[k][j+1][1])*
                                        (samplepoints[k][j][1]-
                                         samplepoints[k][j+1][1])+
                                        (samplepoints[k][j][2]-
                                         samplepoints[k][j+1][2])*
                                        (samplepoints[k][j][2]-
                                         samplepoints[k][j+1][2]));

                if(deltau > DeltaU)
                    DeltaU = deltau;
                if(deltav > DeltaV)
                    DeltaV = deltav;
            }
        }

        //these are the acutal number of sample points in each parametric
        //direction
        int nu = ceil(DeltaU/m_minDelta)*40*1.2;
        int nv = ceil(DeltaV/m_minDelta)*40*1.2;

        for(int j = 0; j < nu; j++)
        {
            for(int k = 0; k < nv; k++)
            {
                Array<OneD, NekDouble> uv(2);
                uv[0] = (bounds[1]- bounds[0])/(nu-1)*j + bounds[0];
                uv[1] = (bounds[3]-bounds[2])/(nv-1)*k + bounds[2];

                //this prevents round off error at the end of the surface
                //may not be neseercary but works
                if(j==nu-1) uv[0]=bounds[1];
                if(k==nv-1) uv[1]=bounds[3];

                Array<OneD, NekDouble> N = surf->N(uv);

                //a zero normal occurs at a signularity, CurvaturePoint
                //cannot be sampled here
                if(N[0]==0 && N[1]==0 && N[2]==0)
                {
                    continue;
                }

                Array<OneD, NekDouble> r = surf->D2(uv);

                //metric and curvature tensors
                NekDouble E = r[3]*r[3] + r[4]*r[4] + r[5]*r[5];
                NekDouble F = r[3]*r[6] + r[4]*r[7] + r[5]*r[8];
                NekDouble G = r[6]*r[6] + r[7]*r[7] + r[8]*r[8];
                NekDouble e = N[0]*r[9] + N[1]*r[10] + N[2]*r[11];
                NekDouble f = N[0]*r[15] + N[1]*r[16] + N[2]*r[17];
                NekDouble g = N[0]*r[12] + N[1]*r[13] + N[2]*r[14];

                //if det is zero cannot invert matrix, R=0 so must skip
                if(E*G-F*F<1E-30)
                {
                    continue;
                }

                NekDouble K, H;

                K = (e*g-f*f)/(E*G-F*F);
                H = 0.5*(e*G-2*f*F+g*E)/(E*G-F*F);

                NekDouble kv[2];
                kv[0] = abs(H + sqrt(H*H-K));
                kv[1] = abs(H - sqrt(H*H-K));

                //create new point based on smallest R, flat surfaces have k=0
                //but still need a point for element estimation
                if(kv[0] != 0 || kv[1] != 0)
                {
                    CurvaturePointSharedPtr newCPoint =
                    MemoryManager<CurvaturePoint>::AllocateSharedPtr
                    (surf->P(uv),
                     1.0/(kv[0] > kv[1] ? kv[0] : kv[1]),
                     N);

                    m_cpList.push_back(newCPoint);
                }else
                {
                    CurvaturePointSharedPtr newCPoint =
                    MemoryManager<CurvaturePoint>::AllocateSharedPtr
                    (surf->P(uv), N);
                    m_cpList.push_back(newCPoint);
                }
            }
        }
    }

    for(int i = 0; i < m_cpList.size(); i++)
    {
        m_cpList[i]->Process(m_minDelta,m_maxDelta,m_eps);
    }
}

}
}
