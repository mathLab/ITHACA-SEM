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

            if(loc[0]>=OctantList[n]->X() && loc[1]>OctantList[n]->Y() && loc[2]>OctantList[n]->Z())
            {
                quad=0;
            }
            else if(loc[0]>OctantList[n]->X() && loc[1]<=OctantList[n]->Y() &&
                     loc[2]>OctantList[n]->Z())
            {
                quad=1;
            }
            else if(loc[0]<=OctantList[n]->X() && loc[1]<=OctantList[n]->Y() &&
                     loc[2]>OctantList[n]->Z())
            {
                quad=2;
            }
            else if(loc[0]<=OctantList[n]->X() && loc[1]>OctantList[n]->Y() &&
                     loc[2]>OctantList[n]->Z())
            {
                quad=3;
            }
            else if(loc[0]>OctantList[n]->X() && loc[1]>OctantList[n]->Y() &&
                     loc[2]<=OctantList[n]->Z()){
                quad=4;
            }
            else if(loc[0]>OctantList[n]->X() && loc[1]<=OctantList[n]->Y() &&
                     loc[2]<=OctantList[n]->Z())
            {
                quad=5;
            }
            else if(loc[0]<=OctantList[n]->X() && loc[1]<=OctantList[n]->Y() &&
                     loc[2]<=OctantList[n]->Z())
            {
                quad=6;
            }
            else if(loc[0]<=OctantList[n]->X() && loc[1]>OctantList[n]->Y() &&
                     loc[2]<=OctantList[n]->Z())
            {
                quad=7;
            }

            n=OctantList[n]->GetChild(quad);

            if(OctantList[n]->isLeaf())
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

        vector<int> dum;
        OctantSharedPtr newOctant =
        MemoryManager<Octant>::AllocateSharedPtr
        ((BoundingBox[1]+BoundingBox[0])/2,
         (BoundingBox[3]+BoundingBox[2])/2,
         (BoundingBox[5]+BoundingBox[4])/2,
         maxdim,
         maxdim,
         maxdim,
         -1, 0, m_cpList, dum);

        OctantList.push_back(newOctant);
        //parent created.

        m_totNotDividing=0;

        if(OctantList[0]->Divide())
        {
            OctantList[0]->LeafFalse();
            subdivide(0);
        }

        int ct=0;
        int maxLevel=0;

        for(int i = 0; i < OctantList.size(); i++)
        {
            if(OctantList[i]->isLeaf()){ct++;}
            if(OctantList[i]->GetLevel()>maxLevel){maxLevel=OctantList[i]->GetLevel();}
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
            if(OctantList[i]->isLeaf())
            {
                OctantList[i]->CreateNeighbourList(OctantList);
            }
        }

        //begin smoothing

        //smooth levels first
        if(m_verbose)
            cout << endl << "\tSmoothing octant levels" << endl;

        SmoothOctants();

        ct=0;
        for(int i = 0; i < OctantList.size(); i++)
        {
            if(OctantList[i]->isLeaf()){ct++;}
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
            if(OctantList[i]->isLeaf())
            {
                ASSERTL0(OctantList[i]->GetDelta()-m_minDelta>-0.0000001,
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
            if(OctantList[i]->isLeaf())
            {
                if(OctantList[i]->GetOrient() != 3 &&
                   OctantList[i]->GetOrient() != -1)
                {
                    NekDouble volumeTet = OctantList[i]->GetDelta()*
                    OctantList[i]->GetDelta()*
                    OctantList[i]->GetDelta()/6.0/sqrt(2.0);

                    NekDouble volumeOct = OctantList[i]->DX()*OctantList[i]->DY()*
                    OctantList[i]->DZ()*8.0;

                    if(OctantList[i]->GetDelta()-m_minDelta<-0.00001)
                    {
                        cout << "error " << OctantList[i]->GetDelta() << endl;
                    }

                    if(OctantList[i]->GetDelta() != OctantList[i]->GetDelta())
                    {
                        cout << "nan in math" <<endl;
                        exit(-1);
                    }
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
                if(OctantList[i]->isLeaf())
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
                        NekDouble deltaSM = 500.0;
                        for(int j = 0; j < checkID.size(); j++)
                        {
                            NekDouble r =
                            sqrt((OctantList[i]->X()-OctantList[checkID[j]]->X())*
                                 (OctantList[i]->X()-OctantList[checkID[j]]->X())+
                                 (OctantList[i]->Y()-OctantList[checkID[j]]->Y())*
                                 (OctantList[i]->Y()-OctantList[checkID[j]]->Y())+
                                 (OctantList[i]->Z()-OctantList[checkID[j]]->Z())*
                                 (OctantList[i]->Z()-OctantList[checkID[j]]->Z()));

                            if(0.24*r +
                               OctantList[checkID[j]]->GetDelta() < deltaSM)
                            {
                                deltaSM = 0.24*r +
                                OctantList[checkID[j]]->GetDelta();
                            }
                        }
                        OctantList[i]->SetDelta(deltaSM);
                        ASSERTL0(deltaSM>=m_minDelta,"Delta assignment less than min delta");
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
                if(OctantList[i]->isLeaf() && !OctantList[i]->isDeltaKnown())
                { //if it is leaf, has no points and delta has not been asigned

                    vector<int> knownID;
                    vector<int> nList = OctantList[i]->GetNeighbourList();

                    for(int j = 0; j<nList.size(); j++)
                    {
                        if(OctantList[nList[j]]->isDeltaKnown())
                        {
                            knownID.push_back(nList[j]);
                        }
                    }									//create list of neighbours where delta is known.
                    if(knownID.size() > 0)
                    {
                        vector<NekDouble> deltaPrime;
                        for(int j = 0; j < knownID.size(); j++)
                        {
                            NekDouble r =
                            sqrt((OctantList[i]->X()-OctantList[knownID[j]]->X())*
                                 (OctantList[i]->X()-OctantList[knownID[j]]->X())+
                                 (OctantList[i]->Y()-OctantList[knownID[j]]->Y())*
                                 (OctantList[i]->Y()-OctantList[knownID[j]]->Y())+
                                 (OctantList[i]->Z()-OctantList[knownID[j]]->Z())*
                                 (OctantList[i]->Z()-OctantList[knownID[j]]->Z()));

                            if(0.24*r +
                               OctantList[knownID[j]]->GetDelta() < m_maxDelta)
                            {
                                deltaPrime.push_back(0.24*r +
                                                     OctantList[knownID[j]]->
                                                     GetDelta());
                            }else{
                                deltaPrime.push_back(m_maxDelta);
                            }
                        }
                        NekDouble min=500.0;
                        for(int j = 0; j < deltaPrime.size(); j++)
                        {
                            if(deltaPrime[j] < min)
                            {
                                min=deltaPrime[j];
                            }
                        }
                        OctantList[i]->SetDelta(min);
                        ASSERTL0(min>=m_minDelta,"Delta assignment less than min delta");
                        ct+=1;

                        deltaPrime.clear();
                    }
                    knownID.clear();
                }


                if(OctantList[i]->isLeaf() && !OctantList[i]->isOrientKnown())
                { //if the node does not know its location
                    vector<int> knownID;
                    vector<int> nList = OctantList[i]->GetNeighbourList();

                    for(int j = 0; j < nList.size(); j++)
                    {
                        if(OctantList[nList[j]]->isOrientKnown())
                        {
                            knownID.push_back(nList[j]);
                        }
                    }
                    if(knownID.size() > 0)
                    {
                        vector<int> idKnowsOrient;
                        for(int j = 0; j < knownID.size(); j++)
                        {
                            if(OctantList[knownID[j]]->isOrientKnown())
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
                                NekDouble dist=10000.0;
                                int closestID;
                                for(int j = 0; j < idKnowsOrient.size(); j++)
                                {
                                    NekDouble r =
                                    sqrt((OctantList[i]->X()-OctantList[idKnowsOrient[j]]->X())*
                                         (OctantList[i]->X()-OctantList[idKnowsOrient[j]]->X())+
                                         (OctantList[i]->Y()-OctantList[idKnowsOrient[j]]->Y())*
                                         (OctantList[i]->Y()-OctantList[idKnowsOrient[j]]->Y())+
                                         (OctantList[i]->Z()-OctantList[idKnowsOrient[j]]->Z())*
                                         (OctantList[i]->Z()-OctantList[idKnowsOrient[j]]->Z()));

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
                                NekDouble dist=10000.0;
                                int closestID;
                                for(int j = 0; j < isOrient2.size(); j++)
                                {
                                    NekDouble r =
                                    sqrt((OctantList[i]->X()-OctantList[isOrient2[j]]->X())*
                                         (OctantList[i]->X()-OctantList[isOrient2[j]]->X())+
                                         (OctantList[i]->Y()-OctantList[isOrient2[j]]->Y())*
                                         (OctantList[i]->Y()-OctantList[isOrient2[j]]->Y())+
                                         (OctantList[i]->Z()-OctantList[isOrient2[j]]->Z())*
                                         (OctantList[i]->Z()-OctantList[isOrient2[j]]->Z()));

                                    if(r < dist)
                                    {
                                        closestID=isOrient2[j];
                                        dist = r;
                                    }
                                }

                                int closestPoint;
                                dist = 10000.0;
                                for(int j = 0; j <
                                    OctantList[closestID]->NumCurvePoint(); j++)
                                {
                                    int CPID = OctantList[closestID]->GetCPID(j);
                                    NekDouble r =
                                    sqrt((OctantList[i]->X()-m_cpList[CPID]->X())*
                                         (OctantList[i]->X()-m_cpList[CPID]->X())+
                                         (OctantList[i]->Y()-m_cpList[CPID]->Y())*
                                         (OctantList[i]->Y()-m_cpList[CPID]->Y())+
                                         (OctantList[i]->Z()-m_cpList[CPID]->Z())*
                                         (OctantList[i]->Z()-m_cpList[CPID]->Z()));

                                    if(r < dist)
                                    {
                                        closestPoint=j;
                                        dist = r;
                                    }
                                }

                                Array<OneD, NekDouble> r(3);
                                int CPID = OctantList[closestID]->GetCPID(closestPoint);
                                r[0] =OctantList[i]->X() - m_cpList[CPID]->X();
                                r[1] =OctantList[i]->Y() - m_cpList[CPID]->Y();
                                r[2] =OctantList[i]->Z() - m_cpList[CPID]->Z();

                                Array<OneD, NekDouble> N(3);
                                m_cpList[CPID]->GetNormal(N[0],N[1],N[2]);

                                NekDouble dot = r[0]*N[0]+r[1]*N[1]+r[2]*N[2];

                                if(dot <= 0.0)
                                {
                                    OctantList[i]->SetOrient(3);
                                }else{
                                    OctantList[i]->SetOrient(1);
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
                if(OctantList[i]->isLeaf() && OctantList[i]->isDeltaKnown())
                {
                    vector<int> checkID;
                    vector<int> nList = OctantList[i]->GetNeighbourList();

                    for(int j = 0; j < nList.size(); j++)
                    {
                        if(OctantList[nList[j]]->isDeltaKnown()
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
                        NekDouble deltaSM = 500.0;
                        for(int j = 0; j < checkID.size(); j++)
                        {
                            NekDouble r =
                            sqrt((OctantList[i]->X()-OctantList[checkID[j]]->X())*
                                 (OctantList[i]->X()-OctantList[checkID[j]]->X())+
                                 (OctantList[i]->Y()-OctantList[checkID[j]]->Y())*
                                 (OctantList[i]->Y()-OctantList[checkID[j]]->Y())+
                                 (OctantList[i]->Z()-OctantList[checkID[j]]->Z())*
                                 (OctantList[i]->Z()-OctantList[checkID[j]]->Z()));

                            if(0.074*r +
                               OctantList[checkID[j]]->GetDelta() < deltaSM)
                            {
                                deltaSM = 0.074*r +
                                OctantList[checkID[j]]->GetDelta();
                            }
                        }
                        OctantList[i]->SetDelta(deltaSM);
                        ASSERTL0(deltaSM>=m_minDelta,"Delta assignment less than min delta");
                        ct+=1;
                    }
                }
            }

        }while(ct>0);
    }

    NekDouble Octree::ddx(int i, int j)
    {
        NekDouble r = sqrt((OctantList[i]->X()-OctantList[j]->X())*
                           (OctantList[i]->X()-OctantList[j]->X())+
                           (OctantList[i]->Y()-OctantList[j]->Y())*
                           (OctantList[i]->Y()-OctantList[j]->Y())+
                           (OctantList[i]->Z()-OctantList[j]->Z())*
                           (OctantList[i]->Z()-OctantList[j]->Z()));
        return abs(OctantList[i]->GetDelta()-OctantList[j]->GetDelta())/r;
    }

    void Octree::SmoothOctants()
    {
        int ct=0;
        int imax=0;

        do
        {
            ct=0;
            int run = OctantList.size();
            for(int j = 0; j < run; j++)
            {
                if(OctantList[j]->isLeaf())
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
        OctantList[parent]->LeafFalse();

        Array<OneD, int> children(8);

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
            (OctantList[parent]->X()+pmx*OctantList[parent]->DX()/2,
             OctantList[parent]->Y()+pmy*OctantList[parent]->DY()/2,
             OctantList[parent]->Z()+pmz*OctantList[parent]->DZ()/2,
             OctantList[parent]->DX()/2,
             OctantList[parent]->DY()/2,
             OctantList[parent]->DZ()/2,
             parent, OctantList[parent]->GetLevel()+1,
             m_cpList, OctantList[parent]->GetCPList());

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


    void Octree::subdivide(int parent)
    {
        Array<OneD, int> children(8);
        //create 8 children for parent and check in turn.

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
            (OctantList[parent]->X()+pmx*OctantList[parent]->DX()/2,
             OctantList[parent]->Y()+pmy*OctantList[parent]->DY()/2,
             OctantList[parent]->Z()+pmz*OctantList[parent]->DZ()/2,
             OctantList[parent]->DX()/2,
             OctantList[parent]->DY()/2,
             OctantList[parent]->DZ()/2,
             parent, OctantList[parent]->GetLevel()+1,
             m_cpList, OctantList[parent]->GetCPList());

            OctantList.push_back(newOctant);
            children[i]=OctantList.size()-1;

            if(OctantList[children[i]]->Divide())
            {
                if(OctantList[children[i]]->DX() > m_minDelta &&
                   OctantList[children[i]]->DY() > m_minDelta &&
                   OctantList[children[i]]->DZ() > m_minDelta)
                {
                    OctantList[children[i]]->LeafFalse();
                    subdivide(children[i]);
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
