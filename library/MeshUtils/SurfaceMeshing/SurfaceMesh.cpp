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
                          new Node(m_mesh->m_meshnode.size(), loc[0], loc[1], loc[2]));
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
        vector<int> cs = m_mesh->m_meshnode[i]->GetListCADCurve();
        for(int j = 0; j < cs.size(); j++)
        {
            vector<CADSurfSharedPtr> s = m_cad->GetCurve(cs[j])->GetAdjSurf();
            for(int k = 0; k < s.size(); k++)
            {
                l.push_back(s[k]->GetId());
            }
        }
        l.sort(); l.unique();
        for (list<int>::iterator lit=l.begin(); lit!=l.end(); ++lit)
        {
            Array<OneD, NekDouble> uv = m_cad->GetSurf(*lit)->locuv(loc);
            m_mesh->m_meshnode[i]->SetCADSurf(*lit,uv);
        }
    }

    //This needs improving to be more general.

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
    }*/

    /*
    //again this needs improving before being used again.

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
                    for(int j = 0; j < 2; j++)
                    {
                        NekDouble dist = Nodes[nn]->Distance(Nodes[n[j]]);
                        m_octree->Modify(Nodes[nn]->GetLoc(), dist*5.0);
                        m_octree->Modify(Nodes[n[j]]->GetLoc(), dist*5.0);
                    }
                    continue;
                }
            }
        }
    }*/

    if(m_mesh->m_verbose)
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
        if(m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(i,m_cad->GetNumSurf(),
                                           "\tSurface meshing");
        }
        m_facemeshes[i] =
            MemoryManager<FaceMesh>::AllocateSharedPtr(i,m_mesh,
                m_cad->GetSurf(i), m_octree, m_curvemeshes);

        m_facemeshes[i]->Mesh();

    }


}

//this mesh is valided that each egde is listed twice in the triangles
void SurfaceMesh::Validate()
{
    if(m_mesh->m_verbose)
        cout << endl << "\tVerifying surface mesh" << endl;

    EdgeSet::iterator it;

    for(it = m_mesh->m_edgeSet.begin(); it != m_mesh->m_edgeSet.end(); it++)
    {
        if((*it)->m_elLink.size() != 2)
        {
            if(m_mesh->m_verbose)
                cout << "\t\tFailed" << endl;
            ASSERTL0(false,"edge not listed twice");
        }
    }

    if(m_mesh->m_verbose)
        cout << "\t\tPassed" << endl;
}

void SurfaceMesh::Optimise()
{
    /*if(m_mesh->m_verbose)
        cout << endl << "\tOptimising linear mesh" << endl;

    //perform two runs of edge swapping based on node connection defect
    for(int q = 0; q <2; q++)
    {
        if(m_mesh->m_verbose)
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
    /*
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
    }*/
}

    /*
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

    if(m_mesh->m_verbose)
    {
        cout << endl << "\tSurface mesh stats:" << endl;
        for(int i = 1; i <= m_cad->GetNumSurf(); i++)
        {
            cout << "\t\tSurface: " << i;
            m_facemeshes[i]->Report();
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
