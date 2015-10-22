////////////////////////////////////////////////////////////////////////////////
//
//  File: SurfaceMeshing.cpp
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
//  Description: surfacemeshing object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <list>
#include <MeshUtils/SurfaceMeshing/SurfaceMesh.h>

#include <LibUtilities/BasicUtils/Progressbar.hpp>
#include <LocalRegions/MatrixKey.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;
namespace Nektar
{
namespace MeshUtils
{

/**
 * @todo write a optimistaion algorithm for edges on curves, not stricly
 * needed because OCC curves are not distorted, but should be included for
 * completness
 */
void SurfaceMeshing::HOSurf()
{
    if(m_verbose)
        cout << endl << "\tHigh-Order Surface meshing" << endl;

    map<int, MeshEdgeSharedPtr>::iterator eit;
    int counter = 0;

    LibUtilities::PointsKey ekey(m_order+1,
                                 LibUtilities::eGaussLobattoLegendre);
    Array<OneD, NekDouble> gll;

    LibUtilities::PointsManager()[ekey]->GetPoints(gll);

    //loop over all edges in the surface
    for(eit = Edges.begin(); eit != Edges.end(); eit++)
    {
        if(m_verbose)
        {
            LibUtilities::PrintProgressbar(counter,Edges.size(),
                                           "\t\tEdges");
        }
        counter++;

        MeshEdgeSharedPtr e = eit->second;
        Array<OneD, int> n = e->GetN();

        //check if edge in on a surface or curve
        if(e->GetCurve() != -1)
        {
            LibUtilities::CADCurveSharedPtr c = m_cad->GetCurve(
                        e->GetCurve());
            //edge is on curve and needs hoing that way
            NekDouble tb = Nodes[n[0]]->GetC(e->GetCurve());
            NekDouble te = Nodes[n[1]]->GetC(e->GetCurve());

            Array<OneD, NekDouble> ti(m_order+1);

            for(int i = 0; i < m_order+1; i++)
            {
                ti[i] = tb*(1.0 -  gll[i])/2.0 +
                        te*(1.0 +  gll[i])/2.0;
            }
            vector<int> Surfs = c->GetAdjSurf();

            ASSERTL0(Surfs.size() == 2, "Number of common surfs should be 2");

            vector<int> honodes(m_order-1);

            LibUtilities::CADSurfSharedPtr s1,s2;
            s1 = m_cad->GetSurf(Surfs[0]);
            s2 = m_cad->GetSurf(Surfs[1]);

            for(int i = 1; i < m_order +1 -1; i++)
            {
                Array<OneD, NekDouble> loc = c->P(ti[i]);
                MeshNodeSharedPtr nn = MemoryManager<MeshNode>::
                        AllocateSharedPtr(Nodes.size(),loc[0],
                                          loc[1],loc[2]);
                nn->SetCurve(c->GetID(),ti[i]);
                Array<OneD, NekDouble> uv = s1->locuv(loc);
                nn->SetSurf(Surfs[0],uv);
                uv = s2->locuv(loc);
                nn->SetSurf(Surfs[1],uv);
                honodes[i-1] = Nodes.size();
                Nodes[Nodes.size()] = nn;
            }

            e->SetHONodes(honodes);
        }
        else
        {
            //edge is on surface and needs 2d optimisation
            LibUtilities::CADSurfSharedPtr s = m_cad->GetSurf(e->GetSurf());
            Array<OneD, NekDouble> uvb,uve;
            uvb = Nodes[n[0]]->GetS(e->GetSurf());
            uve = Nodes[n[1]]->GetS(e->GetSurf());

            vector<int> honodes(m_order-1);
            vector<NekDouble> gllint(m_order-1);
            for(int i = 1; i < m_order+1 -1; i++)
            {
                Array<OneD, NekDouble> loc;
                Array<OneD, NekDouble> uv(2);
                uv[0] = uvb[0]*(1.0 - gll[i])/2.0 + uve[0]*(1.0 + gll[i])/2.0;
                uv[1] = uvb[1]*(1.0 - gll[i])/2.0 + uve[1]*(1.0 + gll[i])/2.0;
                loc = s->P(uv);
                MeshNodeSharedPtr nn = MemoryManager<MeshNode>::
                        AllocateSharedPtr(Nodes.size(),loc[0],
                                            loc[1],loc[2]);
                nn->SetSurf(e->GetSurf(),uv);
                honodes[i-1] = Nodes.size();
                Nodes[Nodes.size()] = nn;
                gllint[i-1] = gll[i];
            }

            //begin optimisation
            bool repeatoverallnodes = true;

            NekDouble tol = 1E-8;

            while(repeatoverallnodes)
            {
                int converged = 0;

                Array<OneD, NekDouble> uvi;

                Array<OneD, NekDouble> bounds = s->GetBounds();

                NekDouble za,zm,zb;

                for(int i = 0; i < honodes.size(); i++)
                {
                    vector<Array<OneD, NekDouble> > bcs;

                    if(i==0)
                    {
                        bcs.push_back(uvb);
                        za = -1.0;
                    }
                    else
                    {
                        bcs.push_back(Nodes[honodes[i-1]]->GetS(e->GetSurf()));
                        za = gllint[i-1];
                    }
                    if(i==honodes.size()-1)
                    {
                        bcs.push_back(uve);
                        zb = 1.0;
                    }
                    else
                    {
                        bcs.push_back(Nodes[honodes[i+1]]->GetS(e->GetSurf()));
                        zb = gllint[i+1];
                    }

                    zm = gllint[i];

                    vector<NekDouble> weights(2);
                    weights[0] = 1.0/(zb - zm);
                    weights[1] = 1.0 /(zm - za);

                    uvi = Nodes[honodes[i]]->GetS(e->GetSurf());

                    bool valid;
                    Array<OneD, NekDouble> df = EdgeGrad(uvi[0],uvi[1],bcs,
                                                         weights,e->GetSurf(),
                                                         valid);
                    if(!valid)
                    {
                        converged++;
                        continue;
                    }

                    NekDouble a,b;

                    Find1DBounds(a,b,uvi,df,bounds);

                    NekDouble fxi = EdgeF(uvi[0],uvi[1],bcs,weights,
                                          e->GetSurf(),valid);
                    NekDouble fx= fxi;

                    NekDouble xmin = BrentOpti(a,0,b,fx,tol,e->GetSurf(),
                                               uvi,df,bounds,bcs,weights,
                                               &SurfaceMeshing::EdgeF);

                    if(fabs(fx - fxi) < tol)
                    {
                        converged++;
                    }
                    else
                    {
                        uvi[0]+=xmin*df[0]; uvi[1]+=xmin*df[1];
                        Array<OneD, NekDouble> loc = s->P(uvi);
                        Nodes[honodes[i]]->Move(loc,uvi);
                    }
                }
                if(converged == honodes.size())
                {
                    repeatoverallnodes = false;
                }
            }
            e->SetHONodes(honodes);
        }
    }

    if(m_verbose)
        cout << endl;

    map<int, MeshTriSharedPtr>::iterator trit;

    //this section of code sets up the standard ho triangle and sorts out
    //node numbering for the optimsation scheme
    LibUtilities::PointsKey pkey(m_order+1,
                                 LibUtilities::eNodalTriFekete);
    Array<OneD, NekDouble> u,v;

    int TotNumPoints = LibUtilities::PointsManager()[pkey]->
                                                    GetTotNumPoints();
    int numInteriorPoints = (m_order-2)*(m_order-1)/2;

    LibUtilities::PointsManager()[pkey]->GetPoints(u,v);

    DNekMat c (3,3,1.0);
    c(0,0) = u[0];
    c(1,0) = v[0];
    c(2,0) = 1.0;
    c(0,1) = u[1];
    c(1,1) = v[1];
    c(2,1) = 1.0;
    c(0,2) = u[2];
    c(1,2) = v[2];
    c(2,2) = 1.0;
    c.Invert();

    DNekMat p (3,numInteriorPoints,1.0);
    for(int j = 0; j < numInteriorPoints; j++)
    {
        p(0,j) = u[TotNumPoints-numInteriorPoints+j];
        p(1,j) = v[TotNumPoints-numInteriorPoints+j];
        p(2,j) = 1.0;
    }

    map<pair<int,int>, int> nodeorder;
    pair<int, int> id;
    id.first = 1;
    id.second  = 1;
    nodeorder[id] = 0;
    id.first = m_order + 1;
    id.second  = 1;
    nodeorder[id] = 1;
    id.first = 1;
    id.second  = m_order + 1;
    nodeorder[id] = 2;

    for(int i = 1; i < m_order+1 -1; i++)
    {
        id.second = 1;
        id.first  = i+1;
        nodeorder[id] = i+2;
    }
    for(int i = 0; i < m_order-1; i++)
    {
        id.first = m_order - i;
        id.second = 2 + i;
        nodeorder[id] = m_order + 1 + 1 + i;
    }
    for(int i = 0; i < m_order-1; i++)
    {
        id.first = 1;
        id.second = m_order - i;
        nodeorder[id] = m_order+1 + m_order + i;
    }
    int i = 2;
    int j = 2;
    int limit = m_order - 1;
    for(int k = 0; k < numInteriorPoints; k++)
    {
        id.first = i;
        id.second = j;
        nodeorder[id] = 3*m_order + k;
        i++;
        if(i > limit)
        {
            limit--;
            j++;
            i=2;
        }
    }

    map<pair<int,int>, vector<NekDouble> > nodeweight;

    //calculates spring weights from the standard triangle
    for(j = 2; j <= m_order-1; j++)
    {
        for(i = 2; i <= m_order + 1 -j; i++)
        {
            NekDouble Zu,Zv,Z1u,Z1v,Z2u,Z2v,Z3u,Z3v;

            Z1u = 0.0; Z1v = 0.0;
            id.first = i-1; id.second = j;
            Z1u += 0.5*u[nodeorder[id]]; Z1v += 0.5*v[nodeorder[id]];
            id.first = i; id.second = j-1;
            Z1u += 0.5*u[nodeorder[id]]; Z1v += 0.5*v[nodeorder[id]];

            Z2u = 0.0; Z2v = 0.0;
            id.first = i+1; id.second = j-1;
            Z2u += 0.5*u[nodeorder[id]]; Z2v += 0.5*v[nodeorder[id]];
            id.first = i+1; id.second = j;
            Z2u += 0.5*u[nodeorder[id]]; Z2v += 0.5*v[nodeorder[id]];

            Z3u = 0.0; Z3v = 0.0;
            id.first = i; id.second = j+1;
            Z3u += 0.5*u[nodeorder[id]]; Z3v += 0.5*v[nodeorder[id]];
            id.first = i-1; id.second = j+1;
            Z3u += 0.5*u[nodeorder[id]]; Z3v += 0.5*v[nodeorder[id]];

            id.first = i; id.second = j;
            Zu = u[nodeorder[id]]; Zv = v[nodeorder[id]];

            vector<NekDouble> K(3);
            K[0] = (Zu - Z2u)*(Zv - Z3v) - (Zv - Z2v)*(Zu - Z3u);
            K[1] = -1.0*((Zu - Z1u)*(Zv - Z3v) - (Zv - Z1v)*(Zu - Z3u));
            K[2] = (Zu - Z1u)*(Zv - Z2v) - (Zv - Z1v)*(Zu - Z2u);

            nodeweight[id] = K;
        }
    }

    counter = 0;

    for(trit = Tris.begin(); trit != Tris.end(); trit++)
    {
        MeshTriSharedPtr t = trit->second;

        if(m_verbose)
        {
            LibUtilities::PrintProgressbar(counter,Tris.size(),
                                           "\t\tFaces");
        }
        counter++;

        Array<OneD, int> n = t->GetN();

        Array<OneD, NekDouble> uv1,uv2,uv3;
        uv1 = Nodes[n[0]]->GetS(t->Getcid());
        uv2 = Nodes[n[1]]->GetS(t->Getcid());
        uv3 = Nodes[n[2]]->GetS(t->Getcid());

        DNekMat a (3,3,1.0);
        a(0,0) = uv1[0];
        a(1,0) = uv1[1];
        a(2,0) = 1.0;
        a(0,1) = uv2[0];
        a(1,1) = uv2[1];
        a(2,1) = 1.0;
        a(0,2) = uv3[0];
        a(1,2) = uv3[1];
        a(2,2) = 1.0;

        DNekMat M = a*c;
        DNekMat result = M*p;

        vector<int> honodes(numInteriorPoints);
        for(int i = 0; i < numInteriorPoints; i++)
        {
            Array<OneD, NekDouble> loc;
            Array<OneD, NekDouble> uv(2);
            uv[0] = result(0,i);
            uv[1] = result(1,i);
            loc = m_cad->GetSurf(t->Getcid())->P(uv);
            MeshNodeSharedPtr nn = MemoryManager<MeshNode>::
                    AllocateSharedPtr(Nodes.size(),loc[0],loc[1],loc[2]);
            nn->SetSurf(t->Getcid(),uv);
            honodes[i] = Nodes.size();
            Nodes[Nodes.size()] = nn;

        }

        //construct a vector of all the uv coords of the triangle in nektar order to form boundary conditions
        vector<Array<OneD, NekDouble> > uvList;
        uvList.push_back(uv1); uvList.push_back(uv2); uvList.push_back(uv3);
        Array<OneD, int> e = t->GetE();
        for(int i = 0; i < 3; i++)
        {
            vector<int> hon = Edges[e[i]]->GetHONodes(n[i]);
            for(int j = 0; j < hon.size(); j++)
            {
                uvList.push_back(Nodes[hon[j]]->GetS(t->Getcid()));
            }
        }
        for(int j = 0; j < honodes.size(); j++)
        {
            uvList.push_back(Nodes[honodes[j]]->GetS(t->Getcid()));
        }

        LibUtilities::CADSurfSharedPtr s = m_cad->GetSurf(t->Getcid());

        bool repeatoverallnodes = true;

        NekDouble tol = 1E-8;

        int ct = 0;

        while(repeatoverallnodes)
        {
            ct++;
            int converged = 0;

            Array<OneD, NekDouble> uvi;

            Array<OneD, NekDouble> bounds = s->GetBounds();
            int node;
            int hocnt = 0;

            for(int j = 2; j <= m_order-1; j++)
            {
                for(int i = 2; i <= m_order + 1 -j; i++)
                {
                    vector<Array<OneD, NekDouble> > bcs;
                    id.first = i-1;
                    id.second = j;
                    node = nodeorder[id];
                    bcs.push_back(uvList[node]);

                    id.first = i;
                    id.second = j-1;
                    node = nodeorder[id];
                    bcs.push_back(uvList[node]);

                    id.first = i+1;
                    id.second = j-1;
                    node = nodeorder[id];
                    bcs.push_back(uvList[node]);

                    id.first = i+1;
                    id.second = j;
                    node = nodeorder[id];
                    bcs.push_back(uvList[node]);

                    id.first = i;
                    id.second = j+1;
                    node = nodeorder[id];
                    bcs.push_back(uvList[node]);

                    id.first = i-1;
                    id.second = j+1;
                    node = nodeorder[id];
                    bcs.push_back(uvList[node]);

                    id.first = i;
                    id.second = j;
                    node = nodeorder[id];
                    uvi = uvList[node];
                    vector<NekDouble> weights = nodeweight[id];

                    bool valid;
                    Array<OneD, NekDouble> df = FaceGrad(uvi[0],uvi[1],bcs,
                                                         weights,t->Getcid(),
                                                         valid);
                    if(!valid)
                    {
                        converged++;
                        continue;
                    }

                    NekDouble a,b;
                    Find1DBounds(a,b,uvi,df,bounds);

                    NekDouble fxi = FaceF(uvi[0],uvi[1],bcs,weights,t->Getcid(),valid);
                    NekDouble fx= fxi;

                    NekDouble xmin = BrentOpti(a,0,b,fx,tol,t->Getcid(),
                                               uvi,df,bounds,bcs,weights,
                                               &SurfaceMeshing::FaceF);

                    if(fabs(fx - fxi) < tol)
                    {
                        converged++;
                    }
                    else
                    {
                        uvi[0]+=xmin*df[0]; uvi[1]+=xmin*df[1];
                        Array<OneD, NekDouble> loc = s->P(uvi);
                        uvList[node] = uvi;
                        Nodes[honodes[hocnt]]->Move(loc,uvi);
                    }

                    hocnt++;
                }
            }
            if(converged == honodes.size())
            {
                repeatoverallnodes = false;
            }
            ASSERTL0(ct<10000,"falid to optismise face");
        }
        t->SetHONodes(honodes);
    }
}



}
}
