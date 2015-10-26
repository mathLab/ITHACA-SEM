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
    if(m_mesh->m_verbose)
        cout << endl << "\tOptimising linear mesh" << endl;

    //perform two runs of edge swapping based on node connection defect
    for(int q = 0; q <2; q++)
    {
        if(m_mesh->m_verbose)
            cout << "\t\t Edge swap defect run: " << q+1 << endl;

        int deletedEdges = 0;

        //loop over all the elements and count the connections to nodes to work out node defect
        NodeSet::iterator nit;
        for(nit = m_mesh->m_vertexSet.begin(); nit != m_mesh->m_vertexSet.end(); nit++)
        {
            (*nit)->m_elCount = 0;
        }
        EdgeSet::iterator it;
        for(it = m_mesh->m_edgeSet.begin(); it != m_mesh->m_edgeSet.end(); it++)
        {
            (*it)->m_n1->m_elCount++;
            (*it)->m_n2->m_elCount++;
        }

        for(it = m_mesh->m_edgeSet.begin(); it != m_mesh->m_edgeSet.end(); it++)
        {
            if((*it)->m_elLink.size() != 2)
                continue;

            ElementSharedPtr tri1 = (*it)->m_elLink[0].first;
            ElementSharedPtr tri2 = (*it)->m_elLink[1].first;

            //if the triangles are not on the same surface then the edge is on a curve
            //therefore cannot edge swap
            if(tri1->GetCADSurf() == tri2->GetCADSurf())
            {
                NodeSharedPtr n1 = (*it)->m_n1;
                NodeSharedPtr n2 = (*it)->m_n2;

                vector<NodeSharedPtr> nt = tri1->GetVertexList();

                //identify node a,b,c,d of the swapping
                NodeSharedPtr A,B,C,D;
                if(nt[0] != n1 && nt[0] != n2)
                {
                    C = nt[0];
                    B = nt[1];
                    A = nt[2];
                }
                else if(nt[1] != n1 && nt[1] != n2)
                {
                    C = nt[1];
                    B = nt[2];
                    A = nt[0];
                }
                else if(nt[2] != n1 && nt[2] != n2)
                {
                    C = nt[2];
                    B = nt[0];
                    A = nt[1];
                }

                nt = tri2->GetVertexList();

                if(nt[0] != n1 && nt[0] != n2)
                {
                    D = nt[0];
                }
                else if(nt[1] != n1 && nt[1] != n2)
                {
                    D = nt[1];
                }
                else if(nt[2] != n1 && nt[2] != n2)
                {
                    D = nt[2];
                }

                NekDouble CBA, BDA;

                //determine signed area of alternate config
                Array<OneD, NekDouble> ai,bi,ci,di;
                ai = A->GetCADSurf(tri1->GetCADSurf());
                bi = B->GetCADSurf(tri1->GetCADSurf());
                ci = C->GetCADSurf(tri1->GetCADSurf());
                di = D->GetCADSurf(tri1->GetCADSurf());

                CBA = (ci[0]*bi[1] + ci[1]*ai[0] + bi[0]*ai[1]) -
                      (ai[0]*bi[1] + ai[1]*ci[0] + bi[0]*ci[1]);

                BDA = (bi[0]*di[1] + bi[1]*ai[0] + di[0]*ai[1]) -
                      (ai[0]*di[1] + ai[1]*bi[0] + di[0]*bi[1]);

                //logic of the next step wont work if this is incorrect so assert it
                if(!(CBA > 0 && BDA > 0))
                    cout << CBA << " " << BDA << endl;
                ASSERTL0(CBA > 0 && BDA > 0, "inverted triangle")

                NekDouble CDA, CBD;

                CDA = (ci[0]*di[1] + ci[1]*ai[0] + di[0]*ai[1]) -
                      (ai[0]*di[1] + ai[1]*ci[0] + di[0]*ci[1]);

                CBD = (ci[0]*bi[1] + ci[1]*di[0] + bi[0]*di[1]) -
                      (di[0]*bi[1] + di[1]*ci[0] + bi[0]*ci[1]);

                //if signed area of the swapping triangles is less than zero
                //that configuration is invalid and swap cannot be performed
                if(CDA > 0.0 && CBD > 0.0)
                {
                    int nodedefectbefore = 0;
                    nodedefectbefore += A->m_elCount > 6 ? A->m_elCount - 6 :
                                                           6 - A->m_elCount;
                    nodedefectbefore += B->m_elCount > 6 ? B->m_elCount - 6 :
                                                           6 - B->m_elCount;
                    nodedefectbefore += C->m_elCount > 6 ? C->m_elCount - 6 :
                                                           6 - C->m_elCount;
                    nodedefectbefore += D->m_elCount > 6 ? D->m_elCount - 6 :
                                                           6 - D->m_elCount;

                    int nodedefectafter = 0;
                    nodedefectafter  += A->m_elCount - 1 > 6 ? A->m_elCount - 1 - 6 :
                                                               6 - (A->m_elCount - 1);
                    nodedefectafter  += B->m_elCount - 1 > 6 ? B->m_elCount - 1 - 6 :
                                                               6 - (B->m_elCount - 1);
                    nodedefectafter  += C->m_elCount + 1 > 6 ? C->m_elCount + 1 - 6 :
                                                               6 - (C->m_elCount + 1);
                    nodedefectafter  += D->m_elCount + 1 > 6 ? D->m_elCount + 1 - 6 :
                                                               6 - (D->m_elCount + 1);

                    //if the node defect of the two triangles will be imrpoved by the
                    //swap perfrom the swap
                    if(nodedefectafter < nodedefectbefore)
                    {
                        //make the 4 other edges
                        EdgeSharedPtr CA, AD, DB, BC, CAt, ADt, DBt, BCt;
                        CAt = boost::shared_ptr<Edge>(new Edge(C,A));
                        ADt = boost::shared_ptr<Edge>(new Edge(A,D));
                        DBt = boost::shared_ptr<Edge>(new Edge(D,B));
                        BCt = boost::shared_ptr<Edge>(new Edge(B,C));

                        for(int i = 0; i < 3; i++)
                        {
                            if(tri1->GetEdge(i) == CAt)
                            {
                                CA = tri1->GetEdge(i);
                            }
                            else if(tri1->GetEdge(i) == BCt)
                            {
                                BC = tri1->GetEdge(i);
                            }

                            if(tri2->GetEdge(i) == ADt)
                            {
                                AD = tri2->GetEdge(i);
                            }
                            else if(tri2->GetEdge(i) == DBt)
                            {
                                DB = tri2->GetEdge(i);
                            }
                        }

                        //now sort out links for the 4 edges surrounding the patch
                        vector<pair<ElementSharedPtr, int> > links;

                        links = CA->m_elLink;
                        CA->m_elLink.clear();
                        for(int i = 0; i < 2; i++)
                        {
                            if(links[i].first == tri1)
                                continue;
                            CA->m_elLink.push_back(links[i]);
                        }

                        links = BC->m_elLink;
                        BC->m_elLink.clear();
                        for(int i = 0; i < 2; i++)
                        {
                            if(links[i].first == tri1)
                                continue;
                            BC->m_elLink.push_back(links[i]);
                        }

                        links = AD->m_elLink;
                        AD->m_elLink.clear();
                        for(int i = 0; i < 2; i++)
                        {
                            if(links[i].first == tri2)
                                continue;
                            AD->m_elLink.push_back(links[i]);
                        }

                        links = DB->m_elLink;
                        DB->m_elLink.clear();
                        for(int i = 0; i < 2; i++)
                        {
                            if(links[i].first == tri2)
                                continue;
                            DB->m_elLink.push_back(links[i]);
                        }

                        vector<NodeSharedPtr> t1,t2;
                        t1.push_back(B); t1.push_back(D); t1.push_back(C);
                        t2.push_back(A); t2.push_back(C); t2.push_back(D);

                        ElmtConfig conf(LibUtilities::eTriangle,1,false,false);
                        vector<int> tags;
                        tags.push_back(tri1->GetCADSurf());

                        int id1 = tri1->GetId();
                        int id2 = tri2->GetId();

                        tri1 = GetElementFactory().
                                    CreateInstance(LibUtilities::eTriangle,
                                                   conf,t1,tags);
                        tri2 = GetElementFactory().
                                    CreateInstance(LibUtilities::eTriangle,
                                                   conf,t2,tags);
                        tri1->SetCADSurf(tags[0]);
                        tri2->SetCADSurf(tags[0]);
                        tri1->SetId(id1);
                        tri2->SetId(id2);

                        (*it)->m_elLink.clear(); //with no element links the code will know to avoid it

                        deletedEdges++;

                        for(int i = 0; i < 3; i++)
                        {
                            EdgeSharedPtr e = tri1->GetEdge(i);
                            pair<EdgeSet::iterator,bool> testIns;
                            testIns = m_mesh->m_edgeSet.insert(e);

                            if(testIns.second)
                            {
                                EdgeSharedPtr ed2 = *testIns.first;
                                ed2->m_id = (*it)->m_id;
                                ed2->m_elLink.push_back(
                                    pair<ElementSharedPtr,int>(tri1,i));
                            }
                            else
                            {
                                EdgeSharedPtr e2 = *testIns.first;
                                tri1->SetEdge(i, e2);
                                e2->m_elLink.push_back(
                                    pair<ElementSharedPtr,int>(tri1,i));
                            }
                        }

                        for(int i = 0; i < 3; i++)
                        {
                            EdgeSharedPtr e = tri2->GetEdge(i);
                            pair<EdgeSet::iterator,bool> testIns;
                            testIns = m_mesh->m_edgeSet.insert(e);

                            if(testIns.second)
                            {
                                EdgeSharedPtr ed2 = *testIns.first;
                                ed2->m_id = (*it)->m_id;
                                ed2->m_elLink.push_back(
                                    pair<ElementSharedPtr,int>(tri2,i));
                            }
                            else
                            {

                                EdgeSharedPtr e2 = *testIns.first;
                                tri2->SetEdge(i, e2);
                                e2->m_elLink.push_back(
                                    pair<ElementSharedPtr,int>(tri2,i));
                            }
                        }

                    }

                }

            }

        }

        if(m_mesh->m_verbose)
            cout << "\t\t\tEdges swapped: " << deletedEdges << endl;

    }

    //perform two runs of edge swapping based on node connection defect
    for(int q = 0; q <2; q++)
    {
        if(m_mesh->m_verbose)
            cout << "\t\t Edge swap angle run: " << q+1 << endl;

        int deletedEdges = 0;

        EdgeSet::iterator it;

        for(it = m_mesh->m_edgeSet.begin(); it != m_mesh->m_edgeSet.end(); it++)
        {
            if((*it)->m_elLink.size() != 2)
                continue;

            ElementSharedPtr tri1 = (*it)->m_elLink[0].first;
            ElementSharedPtr tri2 = (*it)->m_elLink[1].first;

            //if the triangles are not on the same surface then the edge is on a curve
            //therefore cannot edge swap
            if(tri1->GetCADSurf() == tri2->GetCADSurf())
            {
                NodeSharedPtr n1 = (*it)->m_n1;
                NodeSharedPtr n2 = (*it)->m_n2;

                vector<NodeSharedPtr> nt = tri1->GetVertexList();

                //identify node a,b,c,d of the swapping
                NodeSharedPtr A,B,C,D;
                if(nt[0] != n1 && nt[0] != n2)
                {
                    C = nt[0];
                    B = nt[1];
                    A = nt[2];
                }
                else if(nt[1] != n1 && nt[1] != n2)
                {
                    C = nt[1];
                    B = nt[2];
                    A = nt[0];
                }
                else if(nt[2] != n1 && nt[2] != n2)
                {
                    C = nt[2];
                    B = nt[0];
                    A = nt[1];
                }

                nt = tri2->GetVertexList();

                if(nt[0] != n1 && nt[0] != n2)
                {
                    D = nt[0];
                }
                else if(nt[1] != n1 && nt[1] != n2)
                {
                    D = nt[1];
                }
                else if(nt[2] != n1 && nt[2] != n2)
                {
                    D = nt[2];
                }

                NekDouble CBA, BDA;

                //determine signed area of alternate config
                Array<OneD, NekDouble> ai,bi,ci,di;
                ai = A->GetCADSurf(tri1->GetCADSurf());
                bi = B->GetCADSurf(tri1->GetCADSurf());
                ci = C->GetCADSurf(tri1->GetCADSurf());
                di = D->GetCADSurf(tri1->GetCADSurf());

                CBA = (ci[0]*bi[1] + ci[1]*ai[0] + bi[0]*ai[1]) -
                      (ai[0]*bi[1] + ai[1]*ci[0] + bi[0]*ci[1]);

                BDA = (bi[0]*di[1] + bi[1]*ai[0] + di[0]*ai[1]) -
                      (ai[0]*di[1] + ai[1]*bi[0] + di[0]*bi[1]);

                //logic of the next step wont work if this is incorrect so assert it
                if(!(CBA > 0 && BDA > 0))
                    cout << CBA << " " << BDA << endl;
                ASSERTL0(CBA > 0 && BDA > 0, "inverted triangle")

                NekDouble CDA, CBD;

                CDA = (ci[0]*di[1] + ci[1]*ai[0] + di[0]*ai[1]) -
                      (ai[0]*di[1] + ai[1]*ci[0] + di[0]*ci[1]);

                CBD = (ci[0]*bi[1] + ci[1]*di[0] + bi[0]*di[1]) -
                      (di[0]*bi[1] + di[1]*ci[0] + bi[0]*ci[1]);

                //if signed area of the swapping triangles is less than zero
                //that configuration is invalid and swap cannot be performed
                if(CDA > 0.0 && CBD > 0.0)
                {
                    NekDouble minangleb = C->Angle(A,B);
                    minangleb = min(minangleb,B->Angle(C,A));
                    minangleb = min(minangleb,A->Angle(B,C));

                    minangleb = min(minangleb,B->Angle(A,D));
                    minangleb = min(minangleb,D->Angle(B,A));
                    minangleb = min(minangleb,A->Angle(D,B));

                    NekDouble minanglea = C->Angle(D,B);
                    minanglea = min(minanglea,B->Angle(C,D));
                    minanglea = min(minanglea,D->Angle(B,C));

                    minanglea = min(minanglea,C->Angle(A,D));
                    minanglea = min(minanglea,D->Angle(C,A));
                    minanglea = min(minanglea,A->Angle(D,C));

                    //if the node defect of the two triangles will be imrpoved by the
                    //swap perfrom the swap
                    if(minanglea > minangleb)
                    {
                        //make the 4 other edges
                        EdgeSharedPtr CA, AD, DB, BC, CAt, ADt, DBt, BCt;
                        CAt = boost::shared_ptr<Edge>(new Edge(C,A));
                        ADt = boost::shared_ptr<Edge>(new Edge(A,D));
                        DBt = boost::shared_ptr<Edge>(new Edge(D,B));
                        BCt = boost::shared_ptr<Edge>(new Edge(B,C));

                        for(int i = 0; i < 3; i++)
                        {
                            if(tri1->GetEdge(i) == CAt)
                            {
                                CA = tri1->GetEdge(i);
                            }
                            else if(tri1->GetEdge(i) == BCt)
                            {
                                BC = tri1->GetEdge(i);
                            }

                            if(tri2->GetEdge(i) == ADt)
                            {
                                AD = tri2->GetEdge(i);
                            }
                            else if(tri2->GetEdge(i) == DBt)
                            {
                                DB = tri2->GetEdge(i);
                            }
                        }

                        //now sort out links for the 4 edges surrounding the patch
                        vector<pair<ElementSharedPtr, int> > links;

                        links = CA->m_elLink;
                        CA->m_elLink.clear();
                        for(int i = 0; i < 2; i++)
                        {
                            if(links[i].first == tri1)
                                continue;
                            CA->m_elLink.push_back(links[i]);
                        }

                        links = BC->m_elLink;
                        BC->m_elLink.clear();
                        for(int i = 0; i < 2; i++)
                        {
                            if(links[i].first == tri1)
                                continue;
                            BC->m_elLink.push_back(links[i]);
                        }

                        links = AD->m_elLink;
                        AD->m_elLink.clear();
                        for(int i = 0; i < 2; i++)
                        {
                            if(links[i].first == tri2)
                                continue;
                            AD->m_elLink.push_back(links[i]);
                        }

                        links = DB->m_elLink;
                        DB->m_elLink.clear();
                        for(int i = 0; i < 2; i++)
                        {
                            if(links[i].first == tri2)
                                continue;
                            DB->m_elLink.push_back(links[i]);
                        }

                        vector<NodeSharedPtr> t1,t2;
                        t1.push_back(B); t1.push_back(D); t1.push_back(C);
                        t2.push_back(A); t2.push_back(C); t2.push_back(D);

                        ElmtConfig conf(LibUtilities::eTriangle,1,false,false);
                        vector<int> tags;
                        tags.push_back(tri1->GetCADSurf());

                        int id1 = tri1->GetId();
                        int id2 = tri2->GetId();

                        tri1 = GetElementFactory().
                                    CreateInstance(LibUtilities::eTriangle,
                                                   conf,t1,tags);
                        tri2 = GetElementFactory().
                                    CreateInstance(LibUtilities::eTriangle,
                                                   conf,t2,tags);
                        tri1->SetCADSurf(tags[0]);
                        tri2->SetCADSurf(tags[0]);
                        tri1->SetId(id1);
                        tri2->SetId(id2);

                        (*it)->m_elLink.clear(); //with no element links the code will know to avoid it

                        deletedEdges++;

                        for(int i = 0; i < 3; i++)
                        {
                            EdgeSharedPtr e = tri1->GetEdge(i);
                            pair<EdgeSet::iterator,bool> testIns;
                            testIns = m_mesh->m_edgeSet.insert(e);

                            if(testIns.second)
                            {
                                EdgeSharedPtr ed2 = *testIns.first;
                                ed2->m_id = (*it)->m_id;
                                ed2->m_elLink.push_back(
                                    pair<ElementSharedPtr,int>(tri1,i));
                            }
                            else
                            {
                                EdgeSharedPtr e2 = *testIns.first;
                                tri1->SetEdge(i, e2);
                                e2->m_elLink.push_back(
                                    pair<ElementSharedPtr,int>(tri1,i));
                            }
                        }

                        for(int i = 0; i < 3; i++)
                        {
                            EdgeSharedPtr e = tri2->GetEdge(i);
                            pair<EdgeSet::iterator,bool> testIns;
                            testIns = m_mesh->m_edgeSet.insert(e);

                            if(testIns.second)
                            {
                                EdgeSharedPtr ed2 = *testIns.first;
                                ed2->m_id = (*it)->m_id;
                                ed2->m_elLink.push_back(
                                    pair<ElementSharedPtr,int>(tri2,i));
                            }
                            else
                            {

                                EdgeSharedPtr e2 = *testIns.first;
                                tri2->SetEdge(i, e2);
                                e2->m_elLink.push_back(
                                    pair<ElementSharedPtr,int>(tri2,i));
                            }
                        }

                    }

                }

            }

        }

        if(m_mesh->m_verbose)
            cout << "\t\t\tEdges swapped: " << deletedEdges << endl;

    }

    map<int, vector<NodeSharedPtr> > conectingNodes;
    NodeSet::iterator nit;
    for(nit = m_mesh->m_vertexSet.begin(); nit != m_mesh->m_vertexSet.end(); nit++)
    {
        conectingNodes[(*nit)->m_id] = vector<NodeSharedPtr>();
    }
    EdgeSet::iterator eit;
    for(eit = m_mesh->m_edgeSet.begin(); eit != m_mesh->m_edgeSet.end(); eit++)
    {
        if((*eit)->m_elLink.size() !=2)
            continue;

        conectingNodes[(*eit)->m_n1->m_id].push_back((*eit)->m_n2);
        conectingNodes[(*eit)->m_n2->m_id].push_back((*eit)->m_n1);
    }

    //perform 4 runs of elastic relaxation based on the octree
    for(int q = 0; q < 4; q++)
    {
        if(m_mesh->m_verbose)
            cout << "\t\t Elastic relaxation run: " << q+1 << endl;

        for(nit = m_mesh->m_vertexSet.begin(); nit != m_mesh->m_vertexSet.end(); nit++)
        {
            vector<int> c = (*nit)->GetListCADCurve();
            if(c.size()>0) //node is on curve so skip
                continue;

            NekDouble d = m_octree->Query((*nit)->GetLoc());

            vector<int> surfs = (*nit)->GetListCADSurf();
            ASSERTL0(surfs.size()==1, //idiot checking
                        "node should be interior and only be on one surface");

            Array<OneD, NekDouble> uvi = (*nit)->GetCADSurf(surfs[0]);

            vector<NodeSharedPtr> connodes;
            connodes = conectingNodes[(*nit)->m_id];

            vector<NekDouble> om;
            for(int i = 0; i < connodes.size(); i++)
            {
                om.push_back((*nit)->Distance(connodes[i]) - d);
            }

            NekDouble u0=0.0,v0=0.0,fu=0.0,dfu=0.0,fv=0.0,dfv=0.0;
            for(int i = 0; i < connodes.size(); i++)
            {
                Array<OneD, NekDouble> uvj = connodes[i]->GetCADSurf(surfs[0]);
                u0+=uvj[0]/connodes.size();
                v0+=uvj[1]/connodes.size();
            }
            for(int i = 0; i < connodes.size();i++)
            {
                Array<OneD, NekDouble> uvj = connodes[i]->GetCADSurf(surfs[0]);
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
            Array<OneD, NekDouble> bounds = m_cad->GetSurf(surfs[0])->GetBounds();
            uv[0] = u0-fu/dfu; uv[1] = v0-fv/dfv;
            if(!(uv[0] < bounds[0] ||
                       uv[0] > bounds[1] ||
                       uv[1] < bounds[2] ||
                       uv[1] > bounds[3]))
            {
                Array<OneD, NekDouble> l2 = m_cad->GetSurf(surfs[0])->P(uv);
                (*nit)->Move(l2,surfs[0],uv);
            }
        }
    }
}

void SurfaceMesh::Report()
{
    if(m_mesh->m_verbose)
    {
        cout << endl << "\tSurface mesh stats:" << endl;
        for(int i = 1; i <= m_cad->GetNumSurf(); i++)
        {
            cout << "\t\tSurface: " << i;
            m_facemeshes[i]->Report();
        }
    }

    if(m_mesh->m_verbose)
    {
        int ns = m_mesh->m_vertexSet.size();
        int es = m_mesh->m_edgeSet.size();
        int ts = m_mesh->m_element[2].size();
        int ep = ns - es + ts;
        cout << endl << "\tSurface mesh statistics" << endl;
        cout << "\t\tNodes: " << ns << endl;
        cout << "\t\tEdges: " << es << endl;
        cout << "\t\tTriangles " << ts << endl;
        cout << "\t\tEuler-Poincaré characteristic: " << ep << endl;
    }
}


}
}
