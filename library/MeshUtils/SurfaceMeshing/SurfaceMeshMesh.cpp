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

void SurfaceMesh::Mesh()
{
    if(m_mesh->m_verbose)
        cout << endl << "Surface meshing" << endl;

    map<int, CADVertSharedPtr> verts = m_cad->GetVerts();
    map<int, CADVertSharedPtr>::iterator itv;
    for(itv = verts.begin(); itv != verts.end(); itv++)
    {
        Array<OneD, NekDouble> loc = itv->second->GetLoc();
        NodeSharedPtr n = boost::shared_ptr<Node>(
                          new Node(itv->first-1, loc[0], loc[1], loc[2]));
        m_mesh->m_meshnode.push_back(n);
    }

    //linear mesh all curves
    for(int i = 1; i <= m_cad->GetNumCurve(); i++)
    {
        if(m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(i,m_cad->GetNumCurve(),
                                           "\tCurve meshing");
        }

        m_curvemeshes[i] =
            MemoryManager<CurveMesh>::AllocateSharedPtr(i, m_mesh,
                            m_cad->GetCurve(i), m_octree);

        m_curvemeshes[i]->Mesh();

    }

    //all nodes thus far exist on curves on sufaces but do not know about the surface
    for(int i = 0; i < m_mesh->m_meshnode.size(); i++)
    {
        Array<OneD, NekDouble> loc = m_mesh->m_meshnode[i]->GetLoc();
        list<int> l;
        map<int, NekDouble> curves = m_mesh->m_meshnode[i]->GetCurveMap();
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

    /*
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
                    //modify octree
                    /*for(int j = 0; j < 2; j++)
                    {
                        NekDouble dist = Nodes[nn]->Distance(Nodes[n[j]]);
                        m_octree->Modify(Nodes[nn]->GetLoc(), dist*5.0);
                        m_octree->Modify(Nodes[n[j]]->GetLoc(), dist*5.0);
                    }*/
/*                    continue;
                }
            }
        }
    }

    /*m_octree->SmoothAllOctants();

    repeat = true;
    while(repeat)
    {
        repeat = false;
        for(eit = Edges.begin(); eit != Edges.end(); eit++)
        {
            Array<OneD, int> n = eit->second->GetN();
            int c = eit->second->GetCurve();
            ASSERTL0(c != -1, "edge not on curve");
            vector<int> s = m_cad->GetCurve(c)->GetAdjSurf();

            NekDouble dist = Nodes[n[0]]->Distance(Nodes[n[1]]);
            for(int i = 0; i < 2; i++)
            {
                if(dist > m_octree->Query(Nodes[n[i]]->GetLoc())*2.0)
                {
                    int nn = m_curvemeshes[c]->SplitEdge(n[0],n[1],Nodes,Edges);
                    cout << "edgesplit octree" << endl;
                    repeat = true;
                    Array<OneD, NekDouble> loc = Nodes[nn]->GetLoc();
                    for(int j = 0; j < s.size(); j++)
                    {
                        Array<OneD, NekDouble> uvn = m_cad->GetSurf(s[j])->locuv(loc);
                        Nodes[nn]->SetSurf(s[j],uvn);
                    }
                    break;
                }
            }
            if(repeat == true)
            {
                break;
            }
        }
    }*/
/*
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
            if(ang < 0.12 )
            {
                cout << "edge needs to split" << endl;
                vector<int> ts = eit->second->GetTri();
                ASSERTL0(ts.size() == 2, "wrong amount of edges");
                Array<OneD, int> t1n = Tris[ts[0]]->GetN();
                Array<OneD, int> t2n = Tris[ts[1]]->GetN();


                int tn1, tn2;
                tn1 = -1; tn2 = -1;
                for(int i = 0; i < 3; i++)
                {
                    if(t1n[i] == n[0])
                        continue;
                    if(t1n[i] == n[1])
                        continue;
                    tn1 = t1n[i];
                    break;
                }
                for(int i = 0; i < 3; i++)
                {
                    if(t2n[i] == n[0])
                        continue;
                    if(t2n[i] == n[1])
                        continue;
                    tn2 = t2n[i];
                    break;
                }

                ASSERTL0(tn1 != -1 && tn2 != -1, "failed to sort nodes");

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

    Optimise();

    if(m_verbose)
    {
        cout << endl << "\tSurface mesh stats:" << endl;
        for(int i = 1; i <= m_cad->GetNumSurf(); i++)
        {
            cout << "\t\tSurface: " << i;
            m_surfacemeshes[i]->Report();
        }
    }

    if(m_mesh->m_verbose)
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
*/
}

}
}
