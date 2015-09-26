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
#include <MeshUtils/SurfaceMeshing/SurfaceMeshing.h>

#include <LibUtilities/BasicUtils/Progressbar.hpp>
#include <LocalRegions/MatrixKey.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;
namespace Nektar{
namespace MeshUtils {

void SurfaceMeshing::Mesh()
{
    if(m_verbose)
        cout << endl << "Surface meshing" << endl;

    vector<Array<OneD, NekDouble> > vertices = m_cad->GetVerts();
    for(int i = 0; i < vertices.size(); i++)
    {
        MeshNodeSharedPtr n = boost::shared_ptr<MeshNode>(
                          new MeshNode(Nodes.size(),vertices[i][0],
                          vertices[i][1],vertices[i][2]));
        Nodes[Nodes.size()] = n;
    }

    //linear mesh all curves
    for(int i = 1; i <= m_cad->GetNumCurve(); i++)
    {
        if(m_verbose)
        {
            LibUtilities::PrintProgressbar(i,m_cad->GetNumCurve(),
                                           "\tCurve meshing");
        }

        m_curvemeshes[i] =
            MemoryManager<CurveMesh>::AllocateSharedPtr(
                m_verbose, i, m_cad->GetCurve(i), m_octree);

        m_curvemeshes[i]->Mesh(Nodes, Edges);

    }

    //all nodes thus far exist on curves on sufaces but do not know about the surface
    map<int, MeshNodeSharedPtr>::iterator it;
    for(it = Nodes.begin(); it != Nodes.end(); it++)
    {
        Array<OneD, NekDouble> loc = it->second->GetLoc();
        list<int> l;
        map<int, NekDouble> curves = it->second->GetCurveMap();
        map<int, NekDouble>::iterator cit;
        for(cit = curves.begin(); cit != curves.end(); cit++)
        {
            vector<int> s = m_cad->GetCurve(cit->first)->GetAdjSurf();
            for(int i = 0; i < s.size(); i++)
            {
                l.push_back(s[i]);
            }
        }
        l.sort(); l.unique();
        for (list<int>::iterator lit=l.begin(); lit!=l.end(); ++lit)
        {
            Array<OneD, NekDouble> uv = m_cad->GetSurf(*lit)->locuv(loc);
            it->second->SetSurf(*lit,uv);
        }
    }

    //analyse for two node surfaces (not possible)
    for(int i = 1; i <= m_cad->GetNumSurf(); i++)
    {
        if(m_cad->GetSurf(i)->GetTwoC())
        {
            vector<vector<pair<int,int> > > es = m_cad->GetSurf(i)->GetEdges();
            vector<int> cs(2);
            for(int j = 0; j < 2; j++)
            {
                cs[j] = es[0][j].first;
            }
            if(m_curvemeshes[cs[0]]->GetNumPoints() == 2 && m_curvemeshes[cs[1]]->GetNumPoints() == 2)
            {
                vector<int> s1 = m_cad->GetCurve(cs[0])->GetAdjSurf();
                vector<int> s2 = m_cad->GetCurve(cs[1])->GetAdjSurf();

                int nn1 = m_curvemeshes[cs[0]]->SplitEdge(
                                m_curvemeshes[cs[0]]->GetFirstPoint(),
                                m_curvemeshes[cs[0]]->GetLastPoint(),
                                Nodes,Edges);

                int nn2 = m_curvemeshes[cs[1]]->SplitEdge(
                                m_curvemeshes[cs[1]]->GetFirstPoint(),
                                m_curvemeshes[cs[1]]->GetLastPoint(),
                                Nodes,Edges);
                //the new node need its surface information added
                Array<OneD, NekDouble> loc = Nodes[nn1]->GetLoc();
                for(int j = 0; j < s1.size(); j++)
                {
                    Array<OneD, NekDouble> uvn = m_cad->GetSurf(s1[j])->locuv(loc);
                    Nodes[nn1]->SetSurf(s1[j],uvn);
                }
                loc = Nodes[nn2]->GetLoc();
                for(int j = 0; j < s2.size(); j++)
                {
                    Array<OneD, NekDouble> uvn = m_cad->GetSurf(s2[j])->locuv(loc);
                    Nodes[nn2]->SetSurf(s2[j],uvn);
                }
            }
        }
    }

    // this is where high-order awareness can be inserted for curve edges
    map<int, MeshEdgeSharedPtr>::iterator eit;
    bool repeat = true;
    while(repeat)
    {
        repeat = false;
        for(eit = Edges.begin(); eit != Edges.end(); eit++)
        {
            Array<OneD, int> n = eit->second->GetN();
            int c = eit->second->GetCurve();
            ASSERTL0(c != -1, "edge not on curve");
            vector<int> s = m_cad->GetCurve(c)->GetAdjSurf();
            for(int i = 0; i < s.size(); i++)
            {
                Array<OneD, NekDouble> uv1,uv2;
                uv1 = Nodes[n[0]]->GetS(s[i]);
                uv2 = Nodes[n[1]]->GetS(s[i]);
                Array<OneD, NekDouble> N1,N2;
                N1 = m_cad->GetSurf(s[i])->N(uv1);
                N2 = m_cad->GetSurf(s[i])->N(uv2);
                NekDouble ang = N1[0]*N2[0] + N1[1]*N2[1] +N1[2]*N2[2];
                ang /= sqrt(N1[0]*N1[0] + N1[1]*N1[1] + N1[1]*N1[1]);
                ang /= sqrt(N2[0]*N2[0] + N2[1]*N2[1] + N2[1]*N2[1]);
                if(ang < 0.12 )
                {
                    int nn = m_curvemeshes[c]->SplitEdge(n[0],n[1],Nodes,Edges);
                    cout << "split edge" << endl;
                    repeat = true;
                    //the new node need its surface information added
                    Array<OneD, NekDouble> loc = Nodes[nn]->GetLoc();
                    for(int j = 0; j < s.size(); j++)
                    {
                        Array<OneD, NekDouble> uvn = m_cad->GetSurf(s[j])->locuv(loc);
                        Nodes[nn]->SetSurf(s[j],uvn);
                    }
                    continue;
                }
            }
        }
    }

    if(m_verbose)
    {
        cout << endl << "\tCurve mesh stats:" << endl;
        for(int i = 1; i <= m_cad->GetNumCurve(); i++)
        {
            cout << "\t\tCurve: " << i;
            m_curvemeshes[i]->Report();
        }
    }

    //linear mesh all surfaces
    for(int i = 1; i <= m_cad->GetNumSurf(); i++)
    {
        if(m_verbose)
        {
            LibUtilities::PrintProgressbar(i,m_cad->GetNumSurf(),
                                           "\tSurface meshing");
        }
        m_surfacemeshes[i] =
            MemoryManager<SurfaceMesh>::AllocateSharedPtr(i,m_verbose,
                m_cad->GetSurf(i), m_octree, m_curvemeshes);

        m_surfacemeshes[i]->Mesh(Nodes,Edges,Tris);

    }

    Validate();

    Optimise();

    repeat = true;
    while(repeat)
    {
        repeat = false;
        for(eit = Edges.begin(); eit != Edges.end(); eit++)
        {
            Array<OneD, int> n = eit->second->GetN();
            int c = eit->second->GetCurve();
            if(c != -1) continue;

            int s = eit->second->GetSurf();

            Array<OneD, NekDouble> uv1,uv2;
            uv1 = Nodes[n[0]]->GetS(s);
            uv2 = Nodes[n[1]]->GetS(s);
            Array<OneD, NekDouble> N1,N2;
            N1 = m_cad->GetSurf(s)->N(uv1);
            N2 = m_cad->GetSurf(s)->N(uv2);
            NekDouble ang = N1[0]*N2[0] + N1[1]*N2[1] +N1[2]*N2[2];
            ang /= sqrt(N1[0]*N1[0] + N1[1]*N1[1] + N1[1]*N1[1]);
            ang /= sqrt(N2[0]*N2[0] + N2[1]*N2[1] + N2[1]*N2[1]);
            if(ang < 0.05 )
            {
                cout << endl << "need to split edge on surface" << endl;

                vector<int> ts = eit->second->GetTri();
                ASSERTL0(ts.size() == 2, "wrong amount of edges");
                Array<OneD, int> t1n = Tris[ts[0]]->GetN();
                Array<OneD, int> t2n = Tris[ts[1]]->GetN();

                cout << n[0] << " " << n[1] << endl;
                for(int i = 0; i < 3; i++)
                    cout << t1n[i] << " ";
                cout << endl;
                for(int i = 0; i < 3; i++)
                    cout << t2n[i] << " ";
                cout << endl;

                int tn1 = -1 , tn2 = -1;
                for(int i = 0; i < 3; i++)
                {
                    if(t1n[i] != n[0] && t1n[i] != n[1])
                    tn1 = t1n[i];
                    break;
                }
                for(int i = 0; i < 3; i++)
                {
                    if(t2n[i] != n[0] && t2n[i] != n[1])
                    tn2 = t2n[i];
                    break;
                }

                cout << tn1 << " " << tn2 << endl;

                Nodes[n[1]]->RemoveTri(ts[0]); Nodes[n[1]]->RemoveTri(ts[1]);

                Nodes[n[1]]->RemoveEdge(eit->first);

                //find other edges
                int oe1, oe2, oe3, oe4;
                oe1 = Nodes[tn1]->EdgeInCommon(Nodes[n[0]]);
                oe2 = Nodes[n[0]]->EdgeInCommon(Nodes[tn2]);
                oe3 = Nodes[tn1]->EdgeInCommon(Nodes[n[1]]);
                oe4 = Nodes[n[1]]->EdgeInCommon(Nodes[tn2]);

                ASSERTL0(oe1 != -1 && oe2 != -1 && oe3 != -1 && oe4 != -1,
                        "faild to find outer edge");

                Array<OneD, NekDouble> uvn(2);
                uvn[0] = (uv1[0] + uv2[0])/2.0;
                uvn[1] = (uv1[1] + uv2[1])/2.0;

                //putting it in the mid point in the parameter plane is stupid this should be opitmised in some way to be the true center

                Array<OneD, NekDouble> loc = m_cad->GetSurf(s)->P(uvn);
                MeshNodeSharedPtr nsp = boost::shared_ptr<MeshNode>(
                                    new MeshNode(Nodes.size(),loc[0],loc[1],loc[2]));
                nsp->SetSurf(s, uvn);
                int nn = Nodes.size();
                Nodes[Nodes.size()] = nsp;

                int ne1, ne2, ne3;

                MeshEdgeSharedPtr e1 = MemoryManager<MeshEdge>::AllocateSharedPtr(Edges.size(), tn1, nn);
                ne1 = Edges.size();
                Edges[Edges.size()] = e1;
                MeshEdgeSharedPtr e2 = MemoryManager<MeshEdge>::AllocateSharedPtr(Edges.size(), nn, tn2);
                ne2 = Edges.size();
                Edges[Edges.size()] = e2;
                MeshEdgeSharedPtr e3 = MemoryManager<MeshEdge>::AllocateSharedPtr(Edges.size(), nn, n[1]);
                ne3 = Edges.size();
                Edges[Edges.size()] = e3;

                int nt3,nt4;

                MeshTriSharedPtr t3 = boost::shared_ptr<MeshTri>(
                                new MeshTri(Tris.size(),tn1,n[1],nn,oe3,
                                            ne3,ne1,s));
                nt3 = Tris.size();
                Tris[Tris.size()] = t3;

                MeshTriSharedPtr t4 = boost::shared_ptr<MeshTri>(
                                new MeshTri(Tris.size(),n[1],tn2,nn,oe4,
                                            ne2,ne3,s));
                nt4 = Tris.size();
                Tris[Tris.size()] = t4;

                Tris[ts[0]]->Swap(tn1,nn,n[0]);
                Tris[ts[0]]->ResetEdges(ne1,eit->first,oe1);

                Tris[ts[1]]->Swap(nn,tn2,n[0]);
                Tris[ts[1]]->ResetEdges(ne2,oe2,eit->first);

                eit->second->ModifyNodes(n[0],nn);

                Edges[oe3]->RemoveTri(ts[0]);
                Edges[oe3]->SetTri(nt3);

                Edges[oe4]->RemoveTri(ts[1]);
                Edges[oe4]->SetTri(nt4);

                Nodes[tn1]->SetEdge(ne1);
                Nodes[tn1]->SetTri(nt3);
                Nodes[tn2]->SetEdge(ne2);
                Nodes[tn2]->SetTri(nt4);
                Nodes[n[1]]->SetEdge(ne3);
                Nodes[n[1]]->SetTri(nt3);
                Nodes[n[1]]->SetTri(nt4);

                Nodes[nn]->SetEdge(eit->first);
                Nodes[nn]->SetEdge(ne1);
                Nodes[nn]->SetEdge(ne2);
                Nodes[nn]->SetEdge(ne3);
                Nodes[nn]->SetTri(ts[0]);
                Nodes[nn]->SetTri(ts[1]);
                Nodes[nn]->SetTri(nt3);
                Nodes[nn]->SetTri(nt4);

                Edges[ne1]->SetTri(ts[0]);
                Edges[ne1]->SetTri(nt3);
                Edges[ne2]->SetTri(ts[1]);
                Edges[ne2]->SetTri(nt4);
                Edges[ne3]->SetTri(nt4);
                Edges[ne3]->SetTri(nt3);

                Edges[ne1]->SetSurf(s);
                Edges[ne2]->SetSurf(s);
                Edges[ne3]->SetSurf(s);

                repeat = true;
                break;
            }
        }
    }


    if(m_verbose)
    {
        cout << endl << "\tSurface mesh stats:" << endl;
        for(int i = 1; i <= m_cad->GetNumSurf(); i++)
        {
            cout << "\t\tSurface: " << i;
            m_surfacemeshes[i]->Report();
        }
    }

    if(m_verbose)
    {
        int ns = Nodes.size();
        int es = Edges.size();
        int ts = Tris.size();
        int ep = ns - es + ts;
        cout << endl << "\tSurface mesh statistics" << endl;
        cout << "\t\tNodes: " << ns << endl;
        cout << "\t\tEdges: " << es << endl;
        cout << "\t\tTriangles " << ts << endl;
        cout << "\t\tEuler-Poincaré characteristic: " << ep << endl;
    }

}

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

            NekDouble dz = 2.0/m_order;

            Array<OneD, NekDouble> ti(m_order+1);

            for(int i = 0; i < m_order+1; i++)
            {
                ti[i] = tb*((1.0 - (-1.0 + dz*i))/2.0) +
                               te*((1.0 + (-1.0 + dz*i))/2.0);
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
            for(int i = 1; i < m_order+1 -1; i++)
            {
                Array<OneD, NekDouble> loc;
                Array<OneD, NekDouble> uv(2);
                uv[0] = uvb[0]+i*(uve[0]-uvb[0])/m_order;
                uv[1] = uvb[1]+i*(uve[1]-uvb[1])/m_order;
                loc = s->P(uv);
                MeshNodeSharedPtr nn = MemoryManager<MeshNode>::
                        AllocateSharedPtr(Nodes.size(),loc[0],
                                            loc[1],loc[2]);
                nn->SetSurf(e->GetSurf(),uv);
                honodes[i-1] = Nodes.size();
                Nodes[Nodes.size()] = nn;

            }

            //begin optimisation
            bool repeatoverallnodes = true;

            NekDouble tol = 1E-8;

            while(repeatoverallnodes)
            {
                int converged = 0;

                Array<OneD, NekDouble> uv1, uv2, uvx, uvi;

                Array<OneD, NekDouble> bounds = s->GetBounds();

                for(int i = 0; i < honodes.size(); i++)
                {
                    if(i==0)
                    {
                        uv1 = uvb;
                    }
                    else
                    {
                        uv1 = Nodes[honodes[i-1]]->GetS(e->GetSurf());
                    }
                    if(i==honodes.size()-1)
                    {
                        uv2 = uve;
                    }
                    else
                    {
                        uv2 = Nodes[honodes[i+1]]->GetS(e->GetSurf());
                    }

                    uvi = Nodes[honodes[i]]->GetS(e->GetSurf());

                    bool valid;
                    Array<OneD, NekDouble> df = EdgeGrad(uv1,uv2,uvi,e->GetSurf(),valid);
                    if(!valid)
                    {
                        converged++;
                        continue;
                    }

                    NekDouble a,b;

                    Find1DBounds(a,b,uvi,df,bounds);

                    //initial conditions
                    vector<Array<OneD, NekDouble> > bcs;
                    bcs.push_back(uv1); bcs.push_back(uv2);

                    NekDouble fxi = EdgeF(uvi[0],uvi[1],bcs,e->GetSurf());
                    NekDouble fx= fxi;

                    NekDouble xmin = BrentOpti(a,0,b,fx,tol,e->GetSurf(),
                                               uvi,df,bounds,bcs,
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
                                 LibUtilities::eNodalTriEvenlySpaced);
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

    map<pair<int,int>, Array<OneD, NekDouble> > nodeweight;

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

            Array<OneD, NekDouble> K(3);
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

        while(repeatoverallnodes)
        {
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
                    W = nodeweight[id];

                    bool valid;
                    Array<OneD, NekDouble> df = FaceGrad(uvi,bcs,t->Getcid(),
                                                         valid);
                    if(!valid)
                    {
                        converged++;
                        continue;
                    }

                    NekDouble a,b;
                    Find1DBounds(a,b,uvi,df,bounds);

                    NekDouble fxi = FaceF(uvi[0],uvi[1],bcs,t->Getcid());
                    NekDouble fx= fxi;

                    NekDouble xmin = BrentOpti(a,0,b,fx,tol,t->Getcid(),
                                               uvi,df,bounds,bcs,
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
        }
        t->SetHONodes(honodes);
    }
}

void SurfaceMeshing::Optimise()
{
    if(m_verbose)
        cout << endl << "\tOptimising linear mesh" << endl;

    //perform two runs of edge swapping based on node connection defect
    for(int q = 0; q <2; q++)
    {
        if(m_verbose)
            cout << "\t\t Edge swap defect run: " << q+1 << endl;

        map<int, MeshEdgeSharedPtr>::iterator it;
        for(it = Edges.begin(); it != Edges.end(); it++)
        {
            if(it->second->GetCurve() != -1)
                continue;

            vector<int> t = it->second->GetTri();
            ASSERTL0(t.size() == 2, "each edge should have two tri");
            Array<OneD, int> n = it->second->GetN();
            Array<OneD, int> nt = Tris[t[0]]->GetN();
            int nodecheck = 0;
            for(int i = 0; i < 3; i++)
            {
                for(int j = 0; j < 2; j++)
                {
                    if(n[j] == nt[i])
                    {
                        nodecheck++;
                    }
                }
            }
            ASSERTL0(nodecheck == 2, "edge and tri should have 2 n in commom");

            //identify node a,b,c,d of the swapping
            int A,B,C,D;
            if(nt[0] != n[0] && nt[0] != n[1])
            {
                C = nt[0];
                B = nt[1];
                A = nt[2];
            }
            else if(nt[1] != n[0] && nt[1] != n[1])
            {
                C = nt[1];
                B = nt[2];
                A = nt[0];
            }
            else if(nt[2] != n[0] && nt[2] != n[1])
            {
                C = nt[2];
                B = nt[0];
                A = nt[1];
            }

            nt = Tris[t[1]]->GetN();
            nodecheck = 0;
            for(int i = 0; i < 3; i++)
            {
                for(int j = 0; j < 2; j++)
                {
                    if(n[j] == nt[i])
                    {
                        nodecheck++;
                    }
                }
            }
            ASSERTL0(nodecheck == 2, "edge and tri should have 2 n in commom");

            if(nt[0] != n[0] && nt[0] != n[1])
            {
                D = nt[0];
            }
            else if(nt[1] != n[0] && nt[1] != n[1])
            {
                D = nt[1];
            }
            else if(nt[2] != n[0] && nt[2] != n[1])
            {
                D = nt[2];
            }

            //determine signed area of alternate config
            Array<OneD, NekDouble> ai,bi,ci,di;
            ai = Nodes[A]->GetS(it->second->GetSurf());
            bi = Nodes[B]->GetS(it->second->GetSurf());
            ci = Nodes[C]->GetS(it->second->GetSurf());
            di = Nodes[D]->GetS(it->second->GetSurf());

            NekDouble CDA, CBD;

            CDA = (ci[0]*di[1] + ci[1]*ai[0] + di[0]*ai[1]) -
                  (ai[0]*di[1] + ai[1]*ci[0] + di[0]*ci[1]);

            CBD = (ci[0]*bi[1] + ci[1]*di[0] + bi[0]*di[1]) -
                  (di[0]*bi[1] + di[1]*ci[0] + bi[0]*ci[1]);

            //if signed area of the swapping triangles is less than zero
            //that configuration is invalid and swap cannot be performed
            if(CDA < 0.0 || CBD < 0.0)
                continue;

            //determine high-order applicabilty of alternate config
            Array<OneD, NekDouble> Nc, Nd;
            Nc = m_cad->GetSurf(it->second->GetSurf())->N(ci);
            Nd = m_cad->GetSurf(it->second->GetSurf())->N(di);

            NekDouble dot = Nc[0]*Nd[0] + Nc[1]*Nd[1] + Nc[2]*Nd[2];
            if(dot < 0)
                continue;

            int nodedefectbefore = 0;
            nodedefectbefore += Nodes[A]->GetEdges().size() > 6 ?
                                Nodes[A]->GetEdges().size() - 6 :
                                6 - Nodes[A]->GetEdges().size();
            nodedefectbefore += Nodes[B]->GetEdges().size() > 6 ?
                                Nodes[B]->GetEdges().size() - 6 :
                                6 - Nodes[B]->GetEdges().size();
            nodedefectbefore += Nodes[C]->GetEdges().size() > 6 ?
                                Nodes[C]->GetEdges().size() - 6 :
                                6 - Nodes[C]->GetEdges().size();
            nodedefectbefore += Nodes[D]->GetEdges().size() > 6 ?
                                Nodes[D]->GetEdges().size() - 6 :
                                6 - Nodes[D]->GetEdges().size();



            int nodedefectafter = 0;
            nodedefectafter  += Nodes[A]->GetEdges().size() - 1 > 6 ?
                                Nodes[A]->GetEdges().size() - 1 - 6 :
                                6 - (Nodes[A]->GetEdges().size() - 1);
            nodedefectafter  += Nodes[B]->GetEdges().size() - 1 > 6 ?
                                Nodes[B]->GetEdges().size() - 1 - 6 :
                                6 - (Nodes[B]->GetEdges().size() - 1);
            nodedefectafter  += Nodes[C]->GetEdges().size() + 1 > 6 ?
                                Nodes[C]->GetEdges().size() + 1 - 6 :
                                6 - (Nodes[C]->GetEdges().size() + 1);
            nodedefectafter  += Nodes[D]->GetEdges().size() + 1 > 6 ?
                                Nodes[D]->GetEdges().size() + 1 - 6 :
                                6 - (Nodes[D]->GetEdges().size() + 1);

            //if the node defect of the two triangles will be imrpoved by the
            //swap perfrom the swap
            if(nodedefectafter < nodedefectbefore)
            {
                Edges[it->first]->Swap(C,D);
                Nodes[C]->SetEdge(it->first);
                Nodes[D]->SetEdge(it->first);
                Nodes[A]->RemoveEdge(it->first);
                Nodes[B]->RemoveEdge(it->first);
                Nodes[C]->SetTri(t[1]);
                Nodes[B]->RemoveTri(t[1]);
                Nodes[D]->SetTri(t[0]);
                Nodes[A]->RemoveTri(t[0]);

                int e1, e2, e3;
                e1 = Nodes[C]->EdgeInCommon(Nodes[B]);
                e2 = Nodes[B]->EdgeInCommon(Nodes[D]);
                Edges[e2]->RemoveTri(t[1]);
                Edges[e2]->SetTri(t[0]);
                e3 = it->first;
                ASSERTL0(e1 != -1 && e2 != -1 && e3 != -1,"no edge in common");
                Tris[t[0]]->Swap(C,B,D);
                Tris[t[0]]->ResetEdges(e1, e2, e3);

                Tris[t[1]]->Swap(C,D,A);
                e1 = it->first;
                e2 = Nodes[D]->EdgeInCommon(Nodes[A]);
                e3 = Nodes[A]->EdgeInCommon(Nodes[C]);
                Edges[e3]->RemoveTri(t[0]);
                Edges[e3]->SetTri(t[1]);
                ASSERTL0(e1 != -1 && e2 != -1 && e3 != -1,"no edge in common");
                Tris[t[1]]->ResetEdges(e1, e2 ,e3);
            }
        }
    }

    //perform 2 runs of edge swapping based on maximising smallest angle
    for(int q = 0; q <2; q++)
    {
        if(m_verbose)
            cout << "\t\t Edge swap angle run: " << q+1 << endl;

        map<int, MeshEdgeSharedPtr>::iterator it;
        for(it = Edges.begin(); it != Edges.end(); it++)
        {
            if(it->second->GetCurve() != -1)
                continue;

            vector<int> t = it->second->GetTri();
            ASSERTL0(t.size() == 2, "each edge should have two tri");
            Array<OneD, int> n = it->second->GetN();
            Array<OneD, int> nt = Tris[t[0]]->GetN();

            int A,B,C,D;
            if(nt[0] != n[0] && nt[0] != n[1])
            {
                C = nt[0];
                B = nt[1];
                A = nt[2];
            }
            else if(nt[1] != n[0] && nt[1] != n[1])
            {
                C = nt[1];
                B = nt[2];
                A = nt[0];
            }
            else if(nt[2] != n[0] && nt[2] != n[1])
            {
                C = nt[2];
                B = nt[0];
                A = nt[1];
            }

            nt = Tris[t[1]]->GetN();

            if(nt[0] != n[0] && nt[0] != n[1])
            {
                D = nt[0];
            }
            else if(nt[1] != n[0] && nt[1] != n[1])
            {
                D = nt[1];
            }
            else if(nt[2] != n[0] && nt[2] != n[1])
            {
                D = nt[2];
            }

            //determine signed area of alternate config
            Array<OneD, NekDouble> ai,bi,ci,di;
            ai = Nodes[A]->GetS(it->second->GetSurf());
            bi = Nodes[B]->GetS(it->second->GetSurf());
            ci = Nodes[C]->GetS(it->second->GetSurf());
            di = Nodes[D]->GetS(it->second->GetSurf());

            NekDouble CDA, CBD;

            CDA = (ci[0]*di[1] + ci[1]*ai[0] + di[0]*ai[1]) -
                  (ai[0]*di[1] + ai[1]*ci[0] + di[0]*ci[1]);

            CBD = (ci[0]*bi[1] + ci[1]*di[0] + bi[0]*di[1]) -
                  (di[0]*bi[1] + di[1]*ci[0] + bi[0]*ci[1]);

            if(CDA < 0.0 || CBD < 0.0)
                continue;

            //determine high-order applicabilty of alternate config
            Array<OneD, NekDouble> Nc, Nd;
            Nc = m_cad->GetSurf(it->second->GetSurf())->N(ci);
            Nd = m_cad->GetSurf(it->second->GetSurf())->N(di);

            NekDouble dot = Nc[0]*Nd[0] + Nc[1]*Nd[1] + Nc[2]*Nd[2];
            if(dot < 0)
                continue;

            NekDouble minangleb = Nodes[C]->Angle(Nodes[A],Nodes[B]);
            minangleb = min(minangleb,Nodes[B]->Angle(Nodes[C],Nodes[A]));
            minangleb = min(minangleb,Nodes[A]->Angle(Nodes[B],Nodes[C]));

            minangleb = min(minangleb,Nodes[B]->Angle(Nodes[A],Nodes[D]));
            minangleb = min(minangleb,Nodes[D]->Angle(Nodes[B],Nodes[A]));
            minangleb = min(minangleb,Nodes[A]->Angle(Nodes[D],Nodes[B]));

            NekDouble minanglea = Nodes[C]->Angle(Nodes[D],Nodes[B]);
            minanglea = min(minanglea,Nodes[B]->Angle(Nodes[C],Nodes[D]));
            minanglea = min(minanglea,Nodes[D]->Angle(Nodes[B],Nodes[C]));

            minanglea = min(minanglea,Nodes[C]->Angle(Nodes[A],Nodes[D]));
            minanglea = min(minanglea,Nodes[D]->Angle(Nodes[C],Nodes[A]));
            minanglea = min(minanglea,Nodes[A]->Angle(Nodes[D],Nodes[C]));

            if(minanglea > minangleb)
            {
                Edges[it->first]->Swap(C,D);
                Nodes[C]->SetEdge(it->first);
                Nodes[D]->SetEdge(it->first);
                Nodes[A]->RemoveEdge(it->first);
                Nodes[B]->RemoveEdge(it->first);
                Nodes[C]->SetTri(t[1]);
                Nodes[B]->RemoveTri(t[1]);
                Nodes[D]->SetTri(t[0]);
                Nodes[A]->RemoveTri(t[0]);


                int e1, e2, e3;
                e1 = Nodes[C]->EdgeInCommon(Nodes[B]);
                e2 = Nodes[B]->EdgeInCommon(Nodes[D]);
                Edges[e2]->RemoveTri(t[1]);
                Edges[e2]->SetTri(t[0]);
                e3 = it->first;
                ASSERTL0(e1 != -1 && e2 != -1 && e3 != -1,"no edge in common");
                Tris[t[0]]->Swap(C,B,D);
                Tris[t[0]]->ResetEdges(e1, e2, e3);

                Tris[t[1]]->Swap(C,D,A);
                e1 = it->first;
                e2 = Nodes[D]->EdgeInCommon(Nodes[A]);
                e3 = Nodes[A]->EdgeInCommon(Nodes[C]);
                Edges[e3]->RemoveTri(t[0]);
                Edges[e3]->SetTri(t[1]);
                ASSERTL0(e1 != -1 && e2 != -1 && e3 != -1,"no edge in common");
                Tris[t[1]]->ResetEdges(e1, e2 ,e3);

            }
        }
    }

    //perform 4 runs of elastic relaxation based on the octree
    for(int q = 0; q < 4; q++)
    {
        if(m_verbose)
            cout << "\t\t Elastic relaxation run: " << q+1 << endl;

        map<int, MeshNodeSharedPtr>::iterator it;
        for(it = Nodes.begin(); it!=Nodes.end(); it++)
        {
            if(it->second->IsOnACurve())
                continue;

            NekDouble d = m_octree->Query(it->second->GetLoc());

            map<int, Array<OneD, NekDouble> > surf = it->second->GetSurfMap();
            ASSERTL0(surf.size()==1,
                        "node should be interior and only be on one surface");

            map<int, Array<OneD, NekDouble> >::iterator sit = surf.begin();
            int surface = sit->first;
            Array<OneD, NekDouble> uvi = sit->second;

            vector<int> es = it->second->GetEdges();
            vector<int> connodes;
            for(int i = 0; i < es.size(); i++)
            {
                connodes.push_back(Edges[es[i]]->OtherNode(it->first));
            }

            vector<NekDouble> om;
            for(int i = 0; i < connodes.size(); i++)
            {
                om.push_back(it->second->Distance(Nodes[connodes[i]]) - d);
            }

            NekDouble u0=0.0,v0=0.0,fu=0.0,dfu=0.0,fv=0.0,dfv=0.0;
            for(int i = 0; i < connodes.size(); i++)
            {
                Array<OneD, NekDouble> uvj = Nodes[connodes[i]]->GetS(surface);
                u0+=uvj[0]/connodes.size();
                v0+=uvj[1]/connodes.size();
            }
            for(int i = 0; i < connodes.size();i++)
            {
                Array<OneD, NekDouble> uvj = Nodes[connodes[i]]->GetS(surface);
                NekDouble sqr = sqrt((uvj[0]-u0)*(uvj[0]-u0) +
                                     (uvj[1]-v0)*(uvj[1]-v0));
                fu+=om[i]*(uvj[0]-u0)/sqr;
                fv+=om[i]*(uvj[1]-v0)/sqr;
                dfu+=om[i]*sqr*(2*(uvj[0]-u0)*(u0-uvj[0])+
                                    (uvj[1]-v0)*(uvj[1]-v0));
                dfv+=om[i]*sqr*(2*(uvj[1]-v0)*(uvi[1]-v0)+
                                    (uvj[0]-u0)*(uvj[0]-u0));
            }
            Array<OneD, NekDouble> uv(2);
            Array<OneD, NekDouble> bounds = m_cad->GetSurf(surface)->GetBounds();
            uv[0] = u0-fu/dfu; uv[1] = v0-fv/dfv;
            if(!(uv[0] < bounds[0] ||
                       uv[0] > bounds[1] ||
                       uv[1] < bounds[2] ||
                       uv[1] > bounds[3]))
            {
                Array<OneD, NekDouble> l2 = m_cad->GetSurf(surface)->P(uv);
                Nodes[it->first]->Move(l2,uv);
            }
        }
    }
}

//validate the linear mesh from the EPC number caluclated from the CAD
//these should be the same as it is a geomtrical constant for each geometry
//this mesh is also valided that each egde is listed twice in the triangles
void SurfaceMeshing::Validate()
{
    if(m_verbose)
        cout << endl << "\tVerifying surface mesh" << endl;

    if(m_cad->GetEPC() != Nodes.size()-Edges.size()+Tris.size())
    {
        if(m_verbose)
            cout << "\t\tFailed" << endl;
        //ASSERTL0(false,"Euler-Poincaré characteristics do not match");
    }

    map<int,int> edgecheck;
    map<int, MeshEdgeSharedPtr>::iterator ite;
    for(ite=Edges.begin(); ite!=Edges.end(); ite++)
    {
        edgecheck[ite->first] = 0;
    }

    map<int, MeshTriSharedPtr>::iterator it;
    for(it=Tris.begin(); it!=Tris.end(); it++)
    {
        Array<OneD,int> ed = it->second->GetE();
        for(int i = 0; i < 3; i++)
        {
            edgecheck[ed[i]]++;
        }
    }

    map<int,int>::iterator ch;
    for(ch=edgecheck.begin(); ch!=edgecheck.end(); ch++)
    {
        if(ch->second != 2)
        {
            if(m_verbose)
                cout << "\t\tFailed" << endl;
            ASSERTL0(false,"edge not listed twice");
        }
    }
    if(m_verbose)
        cout << "\t\tPassed" << endl;

}

}
}
