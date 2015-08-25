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

#include <string>
#include <fstream>
#include <algorithm>
#include <limits>

#include <MeshUtils/Octree/Octree.h>
#include <LibUtilities/CADSystem/CADSurf.h>

using namespace std;
namespace Nektar{
namespace MeshUtils {

NekDouble Octree::Query(Array<OneD, NekDouble> loc)
{
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

    CompileCuravturePointList();

    if(m_verbose)
        cout << "\tCurvature samples: " << m_cpList.size() << endl
        << "\tInitial subdivision" << endl;

    NekDouble maxdim = (BoundingBox[1]-BoundingBox[0])/2 >
                           (BoundingBox[3]-BoundingBox[2])/2 ?
                           (BoundingBox[1]-BoundingBox[0])/2 :
                           (BoundingBox[3]-BoundingBox[2])/2;
    maxdim = maxdim > (BoundingBox[5]-BoundingBox[4])/2 ?
                        maxdim : (BoundingBox[5]-BoundingBox[4])/2;

    OctantSharedPtr newOctant =
    MemoryManager<Octant>::AllocateSharedPtr
    ((BoundingBox[1]+BoundingBox[0])/2,
     (BoundingBox[3]+BoundingBox[2])/2,
     (BoundingBox[5]+BoundingBox[4])/2, maxdim, -1, 0, m_cpList);

    OctantList.push_back(newOctant);
    //parent created.

    m_totNotDividing=0;

    if(OctantList[0]->Divide())
    {
        OctantList[0]->SetLeaf(false);
        InitialSubDivide(0);
    }

    int ct=0;
    int maxLevel=0;

    for(int i = 0; i < OctantList.size(); i++)
    {
        if(OctantList[i]->GetLeaf())
            ct++;
        if(OctantList[i]->GetLevel()>maxLevel)
            maxLevel=OctantList[i]->GetLevel();
    }

    if(m_verbose)
        cout << "\tNo. octant leaves: " << ct << endl <<
        "\tMax octree level: " << maxLevel << endl;

    if(m_verbose)
        cout << "\tPopulating initial neighbours list..." << endl;

    for(int i = 0; i < OctantList.size(); i++)
    {
        if(m_verbose)
        {
            int pos = 70*i/OctantList.size();
            cout << "\t[";
            for (int j = 0; j < 70; ++j) {
                if (j < pos) cout << "=";
                else if (j == pos) cout << ">";
                else cout << " ";
            }
            cout << "] " << int(float(pos)/(70-1)*100)<< " %\r";
            cout.flush();
        }
        if(OctantList[i]->GetLeaf())
        {
            OctantList[i]->CreateNeighbourList(OctantList);
        }
    }

    //begin smoothing

    //smooth levels first
    if(m_verbose)
        cout << endl << "\tSmoothing octant levels" << endl;

    SubDivideByLevel();

    ct=0;
    for(int i = 0; i < OctantList.size(); i++)
    {
        if(OctantList[i]->GetLeaf()){ct++;}
    }
    cout << "\tNew Stats" << endl;
    cout << "\tNo. octant leaves: " << ct << endl;

    cout << "\tSmoothing across the geometry surface" << endl;

    SmoothSurfaceOctants();

    cout << "\tPropagating spacing out to domain boundary" << endl;

    PropagateDomain();

    cout << "\tRecersively ensuring smoothness between all nodes" << endl;

    SmoothAllOctants();

    for(int i = 0; i < OctantList.size(); i++)
    {
        if(OctantList[i]->GetLeaf())
        {
            ASSERTL0(!(OctantList[i]->GetDelta()<m_minDelta),
                     "Error in initial octree construction");
        }
    }

    int elem=CountElemt();

    cout << endl<< "\tPredicted mesh: " << elem << " elements" << endl;

}

int Octree::CountElemt()
{
    NekDouble total=0.0;

    for(int i = 0 ; i < OctantList.size(); i++)
    {
        if(OctantList[i]->GetLeaf())
        {
            if(OctantList[i]->GetOrient() != 3 &&
               OctantList[i]->GetOrient() != -1)
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
                    ASSERTL0(!(deltaSM<m_minDelta),
                                "Delta assignment less than min delta");
                    ct+=1;
                }
            }
        }

    }while(ct>0);
}

void Octree::PropagateDomain()
{
    int ct=0;

    do
    {
        ct=0;
        for(int i = 0; i<OctantList.size(); i++)
        {
            if(OctantList[i]->GetLeaf() && !OctantList[i]->GetDeltaKnown())
            { //if it is leaf, has no points and delta has not been asigned

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
                    ct+=1;

                    deltaPrime.clear();
                }
                knownID.clear();
            }


            if(OctantList[i]->GetLeaf() && !OctantList[i]->GetOrientKnown())
            { //if the node does not know its location
                vector<int> knownID;
                vector<int> nList = OctantList[i]->GetNeighbourList();

                for(int j = 0; j < nList.size(); j++)
                {
                    if(OctantList[nList[j]]->GetOrientKnown())
                    {
                        knownID.push_back(nList[j]);
                    }
                }
                if(knownID.size() > 0)
                {
                    vector<int> idKnowsOrient;
                    for(int j = 0; j < knownID.size(); j++)
                    {
                        if(OctantList[knownID[j]]->GetOrientKnown())
                        {
                            idKnowsOrient.push_back(knownID[j]);
                        }
                    }
                    if(idKnowsOrient.size() > 0)
                    {

                        vector<int> isOrient2;
                        for(int j = 0; j < idKnowsOrient.size(); j++)
                        {
                            if(OctantList[idKnowsOrient[j]]->GetOrient()==2)
                            {
                                isOrient2.push_back(idKnowsOrient[j]);
                            }
                        }

                        if(isOrient2.size() == 0)
                        {
                            NekDouble dist=numeric_limits<double>::max();
                            int closestID;
                            for(int j = 0; j < idKnowsOrient.size(); j++)
                            {
                                NekDouble r = OctantList[i]->
                                    Distance(OctantList[idKnowsOrient[j]]);

                                if(r < dist)
                                {
                                    closestID=idKnowsOrient[j];
                                    dist = r;
                                }
                            }

                            OctantList[i]->SetOrient(OctantList[closestID]->
                                                    GetOrient());
                            ct+=1;
                            if(OctantList[closestID]->GetOrient()==2)
                            {
                                cout << "error in assignment" << endl;
                            }
                        }
                        else
                        {
                            NekDouble dist=numeric_limits<double>::max();
                            int closestID;
                            for(int j = 0; j < isOrient2.size(); j++)
                            {
                                NekDouble r = OctantList[i]->
                                    Distance(OctantList[isOrient2[j]]);

                                if(r < dist)
                                {
                                    closestID=isOrient2[j];
                                    dist = r;
                                }
                            }

                            CurvaturePointSharedPtr closestPoint;
                            dist = numeric_limits<double>::max();
                            vector<CurvaturePointSharedPtr> cu =
                                        OctantList[closestID]->GetCPList();
                            for(int j = 0; j < cu.size(); j++)
                            {
                                NekDouble r = OctantList[i]->CPDistance(cu[j]);

                                if(r < dist)
                                {
                                    closestPoint=cu[j];
                                    dist = r;
                                }
                            }

                            Array<OneD, NekDouble> r(3);
                            Array<OneD, NekDouble> ocloc =
                                                OctantList[i]->GetLoc();

                            r[0] =ocloc[0] - closestPoint->X();
                            r[1] =ocloc[1] - closestPoint->Y();
                            r[2] =ocloc[2] - closestPoint->Z();

                            Array<OneD, NekDouble> N(3);
                            closestPoint->GetNormal(N[0],N[1],N[2]);

                            NekDouble dot = r[0]*N[0]+r[1]*N[1]+r[2]*N[2];

                            if(dot <= 0.0)
                            {
                                OctantList[i]->SetOrient(1);
                            }else{
                                OctantList[i]->SetOrient(3);
                            }
                            ct+=1;

                        }
                    }

                }
                knownID.clear();
            }
        }

    }while(ct>0);

}

void Octree::SmoothSurfaceOctants()
{
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

                //for each neighbour listed in check_id, figure out the smoothed delta, and asign the miminum of these to nodes[i].GetDelta()
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
                    ASSERTL0(!(deltaSM<m_minDelta),"Delta assignment less than min delta");
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

        int pos = 70*imax/OctantList.size();
        cout << "\t[";
        for (int k = 0; k < 70; ++k) {
            if (k < pos) cout << "=";
            else if (k == pos) cout << ">";
            else cout << " ";
        }
        cout << "] " << int(float(pos)/(70-1)*100)<< " %\r";
        cout.flush();

    }while(ct>0);
    cout <<endl;
}

void Octree::SubDivideLevel(int parent)
{
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
    Array<OneD, int> children(8);
    //create 8 children for parent and check in turn.

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

        if(OctantList[children[i]]->Divide())
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
    NekDouble MaxDim = 0.0;
    if(BoundingBox[1]-BoundingBox[0]>MaxDim)
        MaxDim = BoundingBox[1]-BoundingBox[0];
    if(BoundingBox[3]-BoundingBox[2]>MaxDim)
        MaxDim = BoundingBox[3]-BoundingBox[2];
    if(BoundingBox[5]-BoundingBox[4]>MaxDim)
        MaxDim = BoundingBox[5]-BoundingBox[4];

    for(int i = 1; i <= m_cad->GetNumSurf(); i++)
    {
        LibUtilities::CADSurfSharedPtr surf = m_cad->GetSurf(i);
        Array<OneD, NekDouble> ParameterPlaneBounds = surf->GetBounds();

        NekDouble du = (ParameterPlaneBounds[1]-
                        ParameterPlaneBounds[0])/(40-1);
        NekDouble dv = (ParameterPlaneBounds[3]-
                        ParameterPlaneBounds[2])/(40-1);

        NekDouble DeltaU = 0.0;
        NekDouble DeltaV = 0.0;

        Array<TwoD, Array<OneD, NekDouble> > samplepoints(40,40);

        for(int j = 0; j < 40; j++)
        {
            for(int k = 0; k < 40; k++)
            {
                Array<OneD, NekDouble> uv(2);
                uv[0] = k*du + ParameterPlaneBounds[0];
                uv[1] = j*dv + ParameterPlaneBounds[2];
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

        int nu = ceil(DeltaU/m_minDelta)*40*1.5;
        int nv = ceil(DeltaV/m_minDelta)*40*1.5;

        for(int j = 0; j < nu; j++)
        {
            for(int k = 0; k < nv; k++)
            {
                Array<OneD, NekDouble> uv(2);
                uv[0] = (ParameterPlaneBounds[1]- ParameterPlaneBounds[0])
                                /(nu-1)*j + ParameterPlaneBounds[0];
                uv[1] = (ParameterPlaneBounds[3]-ParameterPlaneBounds[2])
                                /(nv-1)*k + ParameterPlaneBounds[2];
                if(j==nu-1)
                    uv[0]=ParameterPlaneBounds[1]; //These statements prevent floating point error at end of loop
                if(k==nv-1)
                    uv[1]=ParameterPlaneBounds[3];


                Array<OneD, NekDouble> N = surf->N(uv);

                if(N[0]==0 && N[1]==0 && N[2]==0)
                {
                    continue;
                }

                Array<OneD, NekDouble> r = surf->D2(uv);

                NekDouble E = r[3]*r[3] + r[4]*r[4] + r[5]*r[5];
                NekDouble F = r[3]*r[6] + r[4]*r[7] + r[5]*r[8];
                NekDouble G = r[6]*r[6] + r[7]*r[7] + r[8]*r[8];
                NekDouble e = N[0]*r[9] + N[1]*r[10] + N[2]*r[11];
                NekDouble f = N[0]*r[15] + N[1]*r[16] + N[2]*r[17];
                NekDouble g = N[0]*r[12] + N[1]*r[13] + N[2]*r[14];

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

                if(kv[0] != 0 || kv[1] != 0)
                {
                    CurvaturePointSharedPtr newCPoint =
                    MemoryManager<CurvaturePoint>::AllocateSharedPtr
                    (r[0],r[1],r[2],
                     1.0/(kv[0] > kv[1] ? kv[0] : kv[1]),
                     N[0],N[1],N[2]);

                    m_cpList.push_back(newCPoint);
                }else
                {
                    CurvaturePointSharedPtr newCPoint =
                    MemoryManager<CurvaturePoint>::AllocateSharedPtr
                    (r[0],r[1],r[2],
                     N[0],N[1],N[2]);
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
