////////////////////////////////////////////////////////////////////////////////
//
//  File: SurfaceMesh.cpp
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
//  Description: surfacemesh object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <limits>
#include <MeshUtils/SurfaceMeshing/FaceMesh.h>
#include <MeshUtils/ExtLibInterface/TriangleInterface.h>

#include <LocalRegions/MatrixKey.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;
namespace Nektar
{
namespace MeshUtils
{

void FaceMesh::Mesh()
{
    Stretching();

    OrientateCurves();

    int numPoints = 0;
    for(int i = 0; i < orderedLoops.size(); i++)
    {
        numPoints+=orderedLoops[i].size();
    }

    stringstream ss;
    ss << "3 points required for triangulation, " << numPoints << " in loop" << endl;
    ss << "curves: ";
    for(int i = 0; i < m_edgeloops.size(); i++)
    {
        for(int j = 0; j < m_edgeloops[i].edges.size(); j++)
        {
            ss << m_edgeloops[i].edges[j]->GetID() << " ";
        }
    }

    ASSERTL0(numPoints > 2, ss.str());

    //create interface to triangle thirdparty library
    TriangleInterfaceSharedPtr pplanemesh =
        MemoryManager<TriangleInterface>::AllocateSharedPtr();

    vector<Array<OneD, NekDouble> > centers;
    for(int i = 0; i < m_edgeloops.size(); i++)
    {
        centers.push_back(m_edgeloops[i].center);
    }

    pplanemesh->Assign(orderedLoops, centers, m_id, asr/pasr);

    pplanemesh->Mesh();

    pplanemesh->Extract(m_connec);

    bool repeat = true;
    int meshcounter = 1;

    while (repeat)
    {
        repeat = Validate();
        if(!repeat)
        {
            break;
        }
        m_connec.clear();
        pplanemesh->AssignStiener(m_stienerpoints);
        pplanemesh->Mesh();
        pplanemesh->Extract(m_connec);
        meshcounter++;
    }

    NekDouble area = numeric_limits<double>::max();
    //idiot check the elements
    for(int i = 0; i < m_connec.size(); i++)
    {
        Array<OneD, NekDouble> a = m_connec[i][0]->GetCADSurf(m_id);
        Array<OneD, NekDouble> b = m_connec[i][1]->GetCADSurf(m_id);
        Array<OneD, NekDouble> c = m_connec[i][2]->GetCADSurf(m_id);

        area = min(area, 0.5*(-b[0]*a[1] + c[0]*a[1] + a[0]*b[1] - c[0]*b[1] - a[0]*c[1] + b[0]*c[1]));

        if(area <= 0)
        {
            cout << area << endl;
            exit(-1);
        }
    }

    //make new elements and add to list from list of nodes and connectivity from triangle
    for(int i = 0; i < m_connec.size(); i++)
    {
        ElmtConfig conf(LibUtilities::eTriangle, 1, false, false);

        vector<int> tags;
        tags.push_back(m_id);
        ElementSharedPtr E = GetElementFactory().CreateInstance(
                                LibUtilities::eTriangle, conf, m_connec[i], tags);
        E->SetCADSurf(m_id);
        vector<EdgeSharedPtr> edgs = E->GetEdgeList();
        for(int j = 0; j < edgs.size(); j++)
        {
            edgs[j]->CADSurfID.push_back(m_id);
        }
        m_mesh->m_element[m_mesh->m_expDim].push_back(E);
    }
}

void FaceMesh::Report()
{
    int edgec = 0;
    for(int i = 0; i < m_edgeloops.size(); i++)
    {
        edgec+=m_edgeloops[i].edges.size();
    }
    int nump = 0;
    for(int i = 0; i < orderedLoops.size(); i++)
    {
        nump+=orderedLoops[i].size();
    }
    nump+=m_stienerpoints.size();

    cout << scientific << "\tPoints: " << nump << "\tTris: " << m_connec.size() << "\tCAD Edges: " << edgec <<  "\tLoops: " << orderedLoops.size() << endl;

}

void FaceMesh::Stretching()
{
    //define a sampling and calculate the aspect ratio of the paramter plane
    asr = 0.0;
    Array<OneD, NekDouble> bnds = m_cadsurf->GetBounds();
    pasr = (bnds[1] - bnds[0])/
           (bnds[3] - bnds[2]);

    Array<TwoD, Array<OneD,NekDouble> > stretch(40,40);

    NekDouble du = (bnds[1]-bnds[0])/(40-1);
    NekDouble dv = (bnds[3]-bnds[2])/(40-1);

    for(int i = 0; i < 40; i++)
    {
        for(int j = 0; j < 40; j++)
        {
            Array<OneD, NekDouble> uv(2);
            uv[0] = bnds[0] + i*du;
            uv[1] = bnds[2] + j*dv;
            if(i==40-1) uv[0]=bnds[1];
            if(j==40-1) uv[1]=bnds[3];
            stretch[i][j]=m_cadsurf->P(uv);
        }
    }

    int ct = 0;

    for(int i = 0; i < 40-1; i++)
    {
        for(int j = 0; j < 40-1; j++)
        {
            NekDouble ru = sqrt((stretch[i][j][0]-stretch[i+1][j][0])*
                                (stretch[i][j][0]-stretch[i+1][j][0])+
                                (stretch[i][j][1]-stretch[i+1][j][1])*
                                (stretch[i][j][1]-stretch[i+1][j][1])+
                                (stretch[i][j][2]-stretch[i+1][j][2])*
                                (stretch[i][j][2]-stretch[i+1][j][2]));
            NekDouble rv = sqrt((stretch[i][j][0]-stretch[i][j+1][0])*
                                (stretch[i][j][0]-stretch[i][j+1][0])+
                                (stretch[i][j][1]-stretch[i][j+1][1])*
                                (stretch[i][j][1]-stretch[i][j+1][1])+
                                (stretch[i][j][2]-stretch[i][j+1][2])*
                                (stretch[i][j][2]-stretch[i][j+1][2]));

            if(rv < 1E-8)
                continue;

            asr += ru/rv;
            ct++;
        }
    }

    asr/=ct;
}

bool FaceMesh::Validate()
{
    //check all edges in the current mesh for length against the octree
    //if the octree is not conformed to add a new point inside the triangle
    //if no new points are added meshing can stop
    int pointBefore = m_stienerpoints.size();
    for(int i = 0; i < m_connec.size(); i++)
    {
        Array<OneD, NekDouble> triDelta(3);

        Array<OneD, NekDouble> r(3);

        r[0]=m_connec[i][0]->Distance(m_connec[i][1]);
        r[1]=m_connec[i][1]->Distance(m_connec[i][2]);
        r[2]=m_connec[i][2]->Distance(m_connec[i][0]);

        triDelta[0] = m_octree->Query(m_connec[i][0]->GetLoc());
        triDelta[1] = m_octree->Query(m_connec[i][1]->GetLoc());
        triDelta[2] = m_octree->Query(m_connec[i][2]->GetLoc());

        int numValid = 0;

        if(r[0] < triDelta[0])
            numValid++;

        if(r[1] < triDelta[1])
            numValid++;

        if(r[2] < triDelta[2])
            numValid++;

        if(numValid != 3)
        {
            Array<OneD,NekDouble> ainfo,binfo,cinfo;
            ainfo = m_connec[i][0]->GetCADSurf(m_id);
            binfo = m_connec[i][1]->GetCADSurf(m_id);
            cinfo = m_connec[i][2]->GetCADSurf(m_id);

            Array<OneD, NekDouble> uvc(2);
            uvc[0] = (ainfo[0]+binfo[0]+cinfo[0])/3.0;
            uvc[1] = (ainfo[1]+binfo[1]+cinfo[1])/3.0;
            AddNewPoint(uvc);
        }
    }

    if(m_stienerpoints.size() == pointBefore)
    {
        return false;
    }
    else
    {
        return true;
    }
}

void FaceMesh::AddNewPoint(Array<OneD, NekDouble> uv)
{
    //adds a new point but checks that there are no other points nearby first
    Array<OneD, NekDouble> np = m_cadsurf->P(uv);
    NekDouble npDelta = m_octree->Query(np);

    NodeSharedPtr n = boost::shared_ptr<Node>(new Node(0,np[0],np[1],np[2]));

    bool add = true;

    for(int i = 0; i < orderedLoops.size(); i++)
    {
        for(int j = 0; j < orderedLoops[i].size(); j++)
        {
            NekDouble r = orderedLoops[i][j]->Distance(n);

            if(r<npDelta/1.414)
            {
                add = false;
                break;
            }
        }
    }

    if(add)
    {
        for(int i = 0; i < m_stienerpoints.size(); i++)
        {
            NekDouble r = m_stienerpoints[i]->Distance(n);

            if(r<npDelta/1.414)
            {
                add = false;
                break;
            }
        }
    }

    if(add)
    {
        n->SetCADSurf(m_id,uv);
        m_stienerpoints.push_back(n);
    }
}

void FaceMesh::OrientateCurves()
{
    //create list of bounding loop nodes
    for(int i = 0; i < m_edgeloops.size(); i++)
    {
        vector<NodeSharedPtr> cE;
        for(int j = 0; j < m_edgeloops[i].edges.size(); j++)
        {
            int cid = m_edgeloops[i].edges[j]->GetID();
            vector<NodeSharedPtr> edgePoints = m_curvemeshes[cid]->GetMeshPoints();

            int numPoints = m_curvemeshes[cid]->GetNumPoints();

            if(m_edgeloops[i].edgeo[j] == 0)
            {
                for(int k = 0; k < numPoints-1; k++)
                {
                    cE.push_back(edgePoints[k]);
                }
            }
            else
            {
                for(int k = numPoints-1; k >0; k--)
                {
                    cE.push_back(edgePoints[k]);
                }
            }
        }
        orderedLoops.push_back(cE);
    }

    //loops made need to orientate on which is biggest and define holes

    for(int i = 0; i < orderedLoops.size(); i++)
    {
        int half = int(orderedLoops[i].size()/2) - 1;

        NodeSharedPtr n1,n2,nh;

        n1 = orderedLoops[i][0];
        n2 = orderedLoops[i][1];
        nh = orderedLoops[i][half];

        Array<OneD,NekDouble> n1info,n2info,nhinfo;
        n1info = n1->GetCADSurf(m_id);
        n2info = n2->GetCADSurf(m_id);
        nhinfo = nh->GetCADSurf(m_id);

        NekDouble ua = (100.0*n1info[0]+
                        100.0*n2info[0]+
                        1.0* nhinfo[0])/201.0 ;
        NekDouble va = (100.0*n1info[1]+
                        100.0*n2info[1]+
                        1.0* nhinfo[1])/201.0 ;

        Array<OneD, NekDouble> tmp(2);
        tmp[0]=ua;
        tmp[1]=va;
        m_edgeloops[i].center=tmp;
    }

    for(int i = 0; i < orderedLoops.size(); i++)
    {
        NekDouble area=0.0;
        for(int j = 0; j < orderedLoops[i].size()-1; j++)
        {
            Array<OneD,NekDouble> n1info,n2info;
            n1info = orderedLoops[i][j]->GetCADSurf(m_id);
            n2info = orderedLoops[i][j+1]->GetCADSurf(m_id);

            area += -n2info[1]*(n2info[0]-n1info[0])
                    +n1info[0]*(n2info[1]-n1info[1]);
        }
        area*=0.5;
        m_edgeloops[i].area = area;
    }

    int ct=0;

    do
    {
        ct=0;
        for(int i = 0; i < m_edgeloops.size()-1; i++)
        {
            if(fabs(m_edgeloops[i].area)<fabs(m_edgeloops[i+1].area))
            {
                //swap
                vector<NodeSharedPtr> orderedlooptmp = orderedLoops[i];
                EdgeLoop edgeLoopstmp = m_edgeloops[i];

                orderedLoops[i]=orderedLoops[i+1];
                m_edgeloops[i]=m_edgeloops[i+1];

                orderedLoops[i+1]=orderedlooptmp;
                m_edgeloops[i+1]=edgeLoopstmp;

                ct+=1;
            }
        }

    }while(ct>0);

    if(m_edgeloops[0].area<0) //reverse the first uvLoop
    {
        vector<NodeSharedPtr> tmp = orderedLoops[0];
        reverse(tmp.begin(), tmp.end());
        orderedLoops[0]=tmp;
    }

    for(int i = 1; i < orderedLoops.size(); i++)
    {
        if(m_edgeloops[i].area>0) //reverse the loop
        {
            vector<NodeSharedPtr> tmp = orderedLoops[i];
            reverse(tmp.begin(), tmp.end());
            orderedLoops[i]=tmp;
        }
    }
}

}
}
