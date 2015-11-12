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
#include <MeshUtils/CADSystem/CADSurf.h>

using namespace std;
namespace Nektar
{
namespace MeshUtils
{

NekDouble Octree::Query(Array<OneD, NekDouble> loc)
{
    //starting at master octant 0 move through succsesive octants which contain
    //the point loc until a leaf is found
    OctantSharedPtr n = m_masteroct;
    int quad;

    bool found=false;

    while(!found)
    {
        Array<OneD, NekDouble> octloc = n->GetLoc();

        if(!(loc[0] < octloc[0]) &&
           !(loc[1] < octloc[1]) &&
           !(loc[2] < octloc[2]))
        {
            quad=0;
        }
        else if(!(loc[0] < octloc[0]) &&
                !(loc[1] > octloc[1]) &&
                !(loc[2] < octloc[2]))
        {
            quad=1;
        }
        else if(!(loc[0] > octloc[0]) &&
                !(loc[1] > octloc[1]) &&
                !(loc[2] < octloc[2]))
        {
            quad=2;
        }
        else if(!(loc[0] > octloc[0]) &&
                !(loc[1] < octloc[1]) &&
                !(loc[2] < octloc[2]))
        {
            quad=3;
        }
        else if(!(loc[0] < octloc[0]) &&
                !(loc[1] < octloc[1]) &&
                !(loc[2] > octloc[2]))
        {
            quad=4;
        }
        else if(!(loc[0] < octloc[0]) &&
                !(loc[1] > octloc[1]) &&
                !(loc[2] > octloc[2]))
        {
            quad=5;
        }
        else if(!(loc[0] > octloc[0]) &&
                !(loc[1] > octloc[1]) &&
                !(loc[2] > octloc[2]))
        {
            quad=6;
        }
        else if(!(loc[0] > octloc[0]) &&
                !(loc[1] < octloc[1]) &&
                !(loc[2] > octloc[2]))
        {
            quad=7;
        }
        else
        {
            ASSERTL0(false,"Cannot locate quadrant");
        }

        n = n->GetChild(quad);

        if(n->IsLeaf())
        {
            found=true;
        }
    }
    return n->GetDelta();
}

void Octree::Build()
{
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
    MemoryManager<Octant>::AllocateSharedPtr(0,
     (BoundingBox[1]+BoundingBox[0])/2,
     (BoundingBox[3]+BoundingBox[2])/2,
     (BoundingBox[5]+BoundingBox[4])/2, maxdim, m_cpList);

    m_masteroct = newOctant;
    Octants.insert(newOctant);

    m_totNotDividing=0;

    InitialSubDivide(m_masteroct);

    OctantSet::iterator it;

    int ct = 0;
    for(it = Octants.begin(); it != Octants.end(); it++, ct++)
    {
        OctantSharedPtr oct = *it;
        if(m_verbose)
        {
            LibUtilities::PrintProgressbar(ct, Octants.size(),
                                           "\tDetermining neigbours\t");
        }
        if(oct->IsLeaf())
        {
            AssignNeigbours(oct);
        }
    }
    if(m_verbose)
        cout << endl;

    //memory is good down to hear

    SubDivideByLevel();

    SmoothSurfaceOctants();

    if(m_verbose)
        cout << "\tPropagating spacing out to domain boundary" << endl;

    PropagateDomain();

    if(m_verbose)
        cout << "\tRecersively ensuring smoothness between all nodes" << endl;

    SmoothAllOctants();

    for(it = Octants.begin(); it != Octants.end(); it++)
    {
        if((*it)->IsLeaf())
        {
            ASSERTL0(!((*it)->GetDelta() < 0.0),
                     "Error in initial octree construction");
        }
    }

    if(m_verbose)
    {
        int elem = CountElemt();
        printf("\tPredicted mesh: %.2d elements\n", elem);
    }

    exit(-1);
}

void Octree::InitialSubDivide(OctantSharedPtr parent)
{
    //in turn, create 8 child octants for octant parent
    //if that child also needs sub dividing, call this function recursively
    Array<OneD, OctantSharedPtr> children(8);

    for(int i = 0; i < 8; i++)
    {
        //set up x,y,z ordering of the 8 octants
        Array<OneD, NekDouble> dir(3);
        if(i<4)
        {
            dir[2] = +1.0;
            if(i<2)
            {
                dir[0] = +1.0;
            }
            else
            {
                dir[0] = -1.0;
            }
            if(i==0||i==3)
            {
                dir[1] = +1.0;
            }
            else
            {
                dir[1] = -1.0;
            }
        }
        else
        {
            dir[2] = -1.0;
            if(i<6)
            {
                dir[0] = +1.0;
            }
            else
            {
                dir[0] = -1.0;
            }
            if(i==4||i==7)
            {
                dir[1] = +1.0;
            }
            else
            {
                dir[1] = -1.0;
            }
        }

        OctantSharedPtr newOctant = MemoryManager<Octant>::AllocateSharedPtr
                                                    (Octants.size(), parent, dir);

        pair<OctantSet::iterator, bool> test = Octants.insert(newOctant);
        ASSERTL0(test.second,"insertion should not fail");
        children[i] = newOctant;

        if(children[i]->GetDivide())
        {
            if(children[i]->DX() / 2.0 > m_minDelta)
            {
                InitialSubDivide(children[i]);
            }
        }

    }

    parent->SetChildren(children);
}

void Octree::SubDivideByLevel()
{
    //until all subdivision ceases, evaluate each octant and subdivide if
    //the neigbour levels are not smooth.
    bool brk;
    int j = 0;
    int imax=0;
    OctantSet::iterator it;

    do
    {
        brk = false;
        j = 0;
        OctantSharedPtr oct;
        Array<OneD, OctantSharedPtr> newoctants;
        for(it = Octants.begin(); it != Octants.end(); it++, j++)
        {
            if(j > imax)
                imax = j;

            oct = *it;
            if(oct->IsLeaf())
            {
                vector<OctantSharedPtr> nList = oct->GetNeighbourList();

                for(int k = 0; k < nList.size(); k++)
                {
                    if(nList[k]->GetLevel() - oct->GetLevel() > 1)
                    {
                        brk = true;
                        newoctants = SubDivideLevel(oct);
                        break;
                    }
                }
                if(brk)
                    break;
            }

            if(m_verbose)
            {
                LibUtilities::PrintProgressbar(imax, Octants.size(),
                                          "\tSubdivide by level\t");
            }
        }
        if(brk)
        {
            oct->SetChildren(newoctants);
            for(int i = 0; i < 8; i++)
            {
                pair<OctantSet::iterator, bool> test = Octants.insert(newoctants[i]);
                ASSERTL0(test.second,"insertion should not fail");
            }
            for(int i = 0; i < 8; i++)
            {
                AssignNeigbours(newoctants[i]);
            }
            vector<OctantSharedPtr> nList = oct->GetNeighbourList();
            for(int i = 0; i < nList.size(); i++)
            {
                AssignNeigbours(nList[i]);
            }
        }

    }while(brk);

    if(m_verbose)
        cout <<endl;
}

Array<OneD, OctantSharedPtr> Octree::SubDivideLevel(OctantSharedPtr parent)
{
    //create 8 child octants in turn for octant parent
    //after creation, re-evaluate all nessercary neighbour lists

    Array<OneD, OctantSharedPtr> children(8);

    for(int i = 0; i < 8; i++)
    {
        Array<OneD, NekDouble> dir(3);
        if(i<4)
        {
            dir[2] = +1.0;
            if(i<2)
            {
                dir[0] = +1.0;
            }
            else
            {
                dir[0] = -1.0;
            }
            if(i==0||i==3)
            {
                dir[1] = +1.0;
            }
            else
            {
                dir[1] = -1.0;
            }
        }
        else
        {
            dir[2] = -1.0;
            if(i<6)
            {
                dir[0] = +1.0;
            }
            else
            {
                dir[0] = -1.0;
            }
            if(i==4||i==7)
            {
                dir[1] = +1.0;
            }
            else
            {
                dir[1] = -1.0;
            }
        }

        OctantSharedPtr newOctant = MemoryManager<Octant>::AllocateSharedPtr
                                                    (Octants.size()+i, parent, dir);
        children[i] = newOctant;
    }

    return children;

    /*parent->SetChildren(children);

    for(int i = 0; i < 8; i++)
    {
        AssignNeigbours(children[i]);
    }

    //need to revaluate the neighbour list of all the neighbours of the parent
    vector<OctantSharedPtr> nList = parent->GetNeighbourList();
    for(int i = 0; i < nList.size(); i++)
    {
        AssignNeigbours(nList[i]);
    }*/
}

void Octree::SmoothSurfaceOctants()
{
    //for all the octants which are surface containing and know their delta
    //specification, look over all neighbours and ensure the specification
    //between them is smooth
    int ct = 0;
    OctantSet::iterator it;

    do
    {
        ct=0;
        for(it = Octants.begin(); it != Octants.end(); it++)
        {
            OctantSharedPtr oct = *it;

            if(oct->IsLeaf() && oct->IsDeltaKnown())
            {
                vector<OctantSharedPtr> check;
                vector<OctantSharedPtr> nList = oct->GetNeighbourList();

                for(int j = 0; j < nList.size(); j++)
                {
                    if(nList[j]->IsDeltaKnown() && nList[j]->GetDelta()
                                                    < oct->GetDelta()
                       && ddx(oct, nList[j]) > 0.075)
                    {
                        check.push_back(nList[j]);
                    }
                }

                //for each neighbour listed in check_id, figure out the smoothed
                // delta, and asign the miminum of these to nodes[i].GetDelta()
                if(check.size() > 0)
                {
                    NekDouble deltaSM = numeric_limits<double>::max();
                    for(int j = 0; j < check.size(); j++)
                    {
                        NekDouble r = oct->Distance(check[j]);

                        if(0.074 * r + check[j]->GetDelta() < deltaSM)
                        {
                            deltaSM = 0.074 * r + check[j]->GetDelta();
                        }
                    }
                    oct->SetDelta(deltaSM);
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
    //until all octants know their delta specifcation and orientaion
    //look over all octants and if their neighours know either their orientation
    //or specifcation calculate one for this octant
    int ct=0;
    OctantSet::iterator it;

    do
    {
        ct=0;
        for(it = Octants.begin(); it != Octants.end(); it++)
        {
            OctantSharedPtr oct = *it;

            if(oct->IsLeaf() && !oct->IsDeltaKnown())
            { //if it is leaf and delta has not been asigned
                ct+=1;
                vector<OctantSharedPtr> known;
                vector<OctantSharedPtr> nList = oct->GetNeighbourList();

                for(int j = 0; j < nList.size(); j++)
                {
                    if(nList[j]->IsDeltaKnown())
                    {
                        known.push_back(nList[j]);
                    }
                }//create list of neighbours where delta is known.

                if(known.size() > 0)
                {
                    vector<NekDouble> deltaPrime;
                    for(int j = 0; j < known.size(); j++)
                    {
                        NekDouble r = oct->Distance(known[j]);

                        if(0.24*r + known[j]->GetDelta() < m_maxDelta)
                        {
                            deltaPrime.push_back(0.24 * r + known[j]->GetDelta());
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
                    oct->SetDelta(min);
                    ASSERTL0(!(min<m_minDelta),
                            "Delta assignment less than min delta");
                    deltaPrime.clear();
                }
                known.clear();
            }

            //check to see if the octant knows its location if not try and find it
            if(oct->IsLeaf() && !oct->KnowsLocation())
            {
                ct+=1;
                vector<OctantSharedPtr> known;
                vector<OctantSharedPtr> nList = oct->GetNeighbourList();

                for(int j = 0; j < nList.size(); j++)
                {
                    if(nList[j]->KnowsLocation())
                    {
                        known.push_back(nList[j]);
                    }
                }

                if(known.size() > 0)
                {
                    vector<OctantSharedPtr> loc2; //need to know which octants are location 2
                    for(int j = 0; j < known.size(); j++)
                    {
                        if(known[j]->GetLocation() == 2)
                        {
                            loc2.push_back(known[j]);
                        }
                    }
                    if(loc2.size() > 0)
                    {
                        //has a neighbour with boundary section
                        //need to explicity calucalte orienation
                        //otherwise can just copy

                        //get list of curvature points in from all loc2 octants
                        vector<CurvaturePointSharedPtr> cplist;
                        for(int j = 0; j < loc2.size(); j++)
                        {
                            vector<CurvaturePointSharedPtr> l = loc2[j]->GetCPList();
                            cplist.insert(cplist.end(), l.begin(), l.end());
                        }
                        //of all these points find the closest
                        NekDouble dist = numeric_limits<double>::max();
                        CurvaturePointSharedPtr closest;
                        for(int j = 0; j < cplist.size(); j++)
                        {
                            if(oct->CPDistance(cplist[j]) < dist)
                            {
                                dist = oct->CPDistance(cplist[j]);
                                closest = cplist[j];
                            }
                        }

                        //we know the closest point and the octant
                        Array<OneD, NekDouble> octloc, cploc;
                        octloc = oct->GetLoc();
                        cploc = closest->GetLoc();
                        Array<OneD, NekDouble> uv(2);
                        int surf;
                        closest->GetCAD(surf, uv);
                        Array<OneD, NekDouble> SurfNorm = m_cad->GetSurf(surf)->N(uv);
                        Array<OneD, NekDouble> pdir(3);

                        pdir[0] = cploc[0] - octloc[0];
                        pdir[1] = cploc[1] - octloc[1];
                        pdir[2] = cploc[2] - octloc[2];

                        NekDouble dot = SurfNorm[0]*pdir[0] +
                                        SurfNorm[1]*pdir[1] +
                                        SurfNorm[2]*pdir[2];

                        if(dot < 0.0)
                        {
                            oct->SetLocation(1); //outside the domain
                        }
                        else
                        {
                            oct->SetLocation(3); //inside
                        }

                    }
                    else
                    {
                        int newloc = known[0]->GetLocation();
                        oct->SetLocation(newloc);
                        //doing this may cause issues, needs testing and fixing
                    }
                }
                known.clear();
            }
        }

    }while(ct>0);

    for(it = Octants.begin(); it != Octants.end(); it++)
    {
        if((*it)->IsLeaf())
        {
            ASSERTL0((*it)->GetDelta() > 0.0, "leaf delta less than zero");
        }
    }
}

void Octree::SmoothAllOctants()
{
    //until no more changes occur smooth the mesh specification between all
    //octants not particualrly strictly
    int ct = 0;
    OctantSet::iterator it;

    do
    {
        ct=0;
        for(it = Octants.begin(); it != Octants.end(); it++)
        {
            OctantSharedPtr oct = *it;

            if(oct->IsLeaf())
            {
                vector<OctantSharedPtr> check;
                vector<OctantSharedPtr> nList = oct->GetNeighbourList();

                for(int j = 0; j < nList.size(); j++)
                {
                    if(!nList[j]->IsLeaf()) //this should not happen but does, this is a bit hacky but fixes it
                    {
                        continue;
                    }
                    if(nList[j]->GetDelta() < oct->GetDelta() && ddx(oct, nList[j]) > 0.3)
                    {
                        check.push_back(nList[j]);
                    }
                }

                //for each neighbour listed in check, figure out the smoothed delta, and asign the miminum of these to nodes[i].GetDelta()
                if(check.size() > 0)
                {
                    NekDouble deltaSM = numeric_limits<double>::max();
                    for(int j = 0; j < check.size(); j++)
                    {
                        NekDouble r = oct->Distance(check[j]);

                        if(0.29 * r + check[j]->GetDelta() < deltaSM)
                        {
                            deltaSM = 0.29 * r + check[j]->GetDelta();
                        }
                    }
                    oct->SetDelta(deltaSM);
                    ct+=1;
                }
            }
        }

    }while(ct>0);
}

int Octree::CountElemt()
{
    //by considering the volume of a tet evaluate the number of elements in the
    //mesh

    NekDouble total = 0.0;
    NekDouble altot = 0.0;
    OctantSet::iterator it;

    int ctcor = 0;
    int ctinc = 0;
    int c2 = 0;

    for(it = Octants.begin(); it != Octants.end(); it++)
    {
        OctantSharedPtr oct = *it;
        if(oct->IsLeaf())
        {
            if(oct->GetLocation() == 1)
            {
                total += 8.0*oct->DX()*oct->DX()*oct->DX() /
                 (oct->GetDelta()*oct->GetDelta()*oct->GetDelta()/6.0/sqrt(2));

                Array<OneD, NekDouble> loc = oct->GetLoc();
                if(!m_cad->InsideShape(loc))
                {
                    ctinc++;
                }
                else
                {
                    ctcor++;
                }
            }
            if(oct->GetLocation() == 2)
                c2++;

            Array<OneD, NekDouble> loc = oct->GetLoc();
            if(m_cad->InsideShape(loc) || oct->GetLocation() == 2)
            {
                altot += 8.0*oct->DX()*oct->DX()*oct->DX() /
                 (oct->GetDelta()*oct->GetDelta()*oct->GetDelta()/6.0/sqrt(2));
            }
        }
    }

    cout << ctcor << " correct " << ctinc << " incorcect " << c2 << " in 2 loc" << endl;
    cout << "alternative total " << int(altot) << endl;

    return int(total);
}

void Octree::CompileCuravturePointList()
{
    for(int i = 1; i <= m_cad->GetNumSurf(); i++)
    {
        CADSurfSharedPtr surf = m_cad->GetSurf(i);
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
        int nu = ceil(DeltaU/m_minDelta)*40;
        int nv = ceil(DeltaV/m_minDelta)*40;

        for(int j = 0; j < nu; j++)
        {
            for(int k = 0; k < nv; k++)
            {
                Array<OneD, NekDouble> uv(2);
                uv[0] = (bounds[1]-bounds[0])/(nu-1)*j + bounds[0];
                uv[1] = (bounds[3]-bounds[2])/(nv-1)*k + bounds[2];

                //this prevents round off error at the end of the surface
                //may not be neseercary but works
                if(j==nu-1) uv[0]=bounds[1];
                if(k==nv-1) uv[1]=bounds[3];

                NekDouble C = surf->Curvature(uv);

                //create new point based on smallest R, flat surfaces have k=0
                //but still need a point for element estimation
                if(C != 0.0)
                {
                    bool minlimited = false;
                    NekDouble ideal;

                    NekDouble del = 2.0*(1.0/C)*sqrt(m_eps*(2.0-m_eps));

                    if(del>m_maxDelta)
                    {
                        del = m_maxDelta;
                    }
                    if(del<m_minDelta)
                    {
                        ideal = del;
                        del = m_minDelta;
                        minlimited = true;
                    }

                    if(minlimited)
                    {
                        CurvaturePointSharedPtr newCPoint =
                            MemoryManager<CurvaturePoint>::AllocateSharedPtr
                            (surf->GetId(), uv, surf->P(uv), del, ideal);

                        m_cpList.push_back(newCPoint);
                    }
                    else
                    {
                        CurvaturePointSharedPtr newCPoint =
                            MemoryManager<CurvaturePoint>::AllocateSharedPtr
                            (surf->GetId(), uv, surf->P(uv), del);

                        m_cpList.push_back(newCPoint);
                    }

                }else
                {
                    CurvaturePointSharedPtr newCPoint =
                        MemoryManager<CurvaturePoint>::AllocateSharedPtr
                        (surf->GetId(), uv, surf->P(uv));

                    m_cpList.push_back(newCPoint);
                }
            }
        }
    }
}

NekDouble Octree::ddx(OctantSharedPtr i, OctantSharedPtr j)
{
    return fabs(i->GetDelta() - j->GetDelta()) / i->Distance(j);
}

void Octree::AssignNeigbours(OctantSharedPtr const &o)
{
    //clear old list
    o->ClearNeigbourList();

    OctantSet::iterator it;
    vector<OctantSharedPtr> nlist;

    int eqhit = 0;

    //look over all octants and consider if they are neigbours
    for(it = Octants.begin(); it != Octants.end(); it++)
    {
        OctantSharedPtr oct = *it;

        if(oct->IsLeaf())
        {
            if(oct == o)
            {
                eqhit++;
                continue;
            }

            //work out the max distance between the two octants if they were
            //joined corner to corner
            NekDouble rmax = o->DiagonalDim() + oct->DiagonalDim();

            NekDouble ractual = o->Distance(oct);

            //if the actucal distance between them is greater that this
            //this octant should not be considered, massive speed up
            if(ractual > 1.1*rmax)
            {
                continue;
            }

            //check overlapping in 2-D for all three dimensions
            if(fabs(o->FX(-1) - oct->FX(+1)) < 1E-6 ||
               fabs(o->FX(+1) - oct->FX(-1)) < 1E-6 )
            {
                //check yz rects
                bool Cond1 = false;
                bool Cond2 = false;
                bool Cond3 = false;
                bool Cond4 = false;
                if(o->FY(+1) < oct->FY(-1))
                    Cond1 = true;
                if(o->FY(-1) > oct->FY(+1))
                    Cond1 = true;
                if(o->FZ(+1) < oct->FZ(-1))
                    Cond1 = true;
                if(o->FZ(-1) > oct->FZ(+1))
                    Cond1 = true;
                if(!Cond1 && !Cond2 && !Cond3 && !Cond4)
                {
                    nlist.push_back(oct);
                }
            }
            else if(fabs(o->FY(-1) - oct->FY(+1)) <1E-6 ||
                    fabs(o->FY(+1) - oct->FY(-1)) <1E-6)
            {
                //check xz rects
                bool Cond1 = false;
                bool Cond2 = false;
                bool Cond3 = false;
                bool Cond4 = false;
                if(o->FX(+1) < oct->FX(-1))
                    Cond1 = true;
                if(o->FX(-1) > oct->FX(+1))
                    Cond1 = true;
                if(o->FZ(+1) < oct->FZ(-1))
                    Cond1 = true;
                if(o->FZ(-1) > oct->FZ(+1))
                    Cond1 = true;
                if(!Cond1 && !Cond2 && !Cond3 && !Cond4)
                {
                    nlist.push_back(oct);
                }
            }
            else if(fabs(o->FZ(-1) - oct->FZ(+1)) <1E-6 ||
                    fabs(o->FZ(+1) - oct->FZ(-1)) <1E-6)
            {
                //check xy rects
                bool Cond1 = false;
                bool Cond2 = false;
                bool Cond3 = false;
                bool Cond4 = false;
                if(o->FX(+1) < oct->FX(-1))
                    Cond1 = true;
                if(o->FX(-1) > oct->FX(+1))
                    Cond1 = true;
                if(o->FY(+1) < oct->FY(-1))
                    Cond1 = true;
                if(o->FY(-1) > oct->FY(+1))
                    Cond1 = true;
                if(!Cond1 && !Cond2 && !Cond3 && !Cond4)
                {
                    nlist.push_back(oct);
                }
            }
        }
    }
    o->SetNeigbourList(nlist);
    ASSERTL0(nlist.size()>6,"not enough neigbours");
    ASSERTL0(eqhit == 1,"found itself more than once");
}

}
}
