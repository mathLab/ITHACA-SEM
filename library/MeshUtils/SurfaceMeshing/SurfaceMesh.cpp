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

    if(m_mesh->m_verbose)
        cout << endl << "\tCurve meshing:" << endl << endl;

    m_mesh->m_numNodes = m_cad->GetNumVerts();

    //linear mesh all curves
    for(int i = 1; i <= m_cad->GetNumCurve(); i++)
    {
        if(m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(i,m_cad->GetNumCurve(),
                                           "Curve progress");
        }

        m_curvemeshes[i] =
            MemoryManager<CurveMesh>::AllocateSharedPtr(i, m_mesh,
                            m_cad->GetCurve(i), m_octree);

        m_curvemeshes[i]->Mesh();

    }

    if(m_mesh->m_verbose)
        cout << endl << "\tFace meshing:" << endl << endl;

    //linear mesh all surfaces
    for(int i = 1; i <= m_cad->GetNumSurf(); i++)
    {
        if(m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(i,m_cad->GetNumSurf(),
                                           "Face progress");
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

    //check elements have positive area in paramter plane after processing
    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        Array<OneD, NekDouble> a = m_mesh->m_element[2][i]->GetVertex(0)->GetCADSurf(m_mesh->m_element[2][i]->GetCADSurf());
        Array<OneD, NekDouble> b = m_mesh->m_element[2][i]->GetVertex(1)->GetCADSurf(m_mesh->m_element[2][i]->GetCADSurf());
        Array<OneD, NekDouble> c = m_mesh->m_element[2][i]->GetVertex(2)->GetCADSurf(m_mesh->m_element[2][i]->GetCADSurf());

        NekDouble area = 0.5*(-b[0]*a[1] + c[0]*a[1] + a[0]*b[1] - c[0]*b[1] - a[0]*c[1] + b[0]*c[1]);

        if(!(area > 0))
        {
            if(m_mesh->m_verbose)
                cout << "\t\tFailed" << endl;
            cout << i << " " << area << endl;
            ASSERTL0(false,"area negative");
        }
    }

    if(m_mesh->m_verbose)
        cout << "\t\tPassed" << endl;
}

/*void SurfaceMesh::HOAwareness()
{
    //this is where patches of the surface mesh are manually remeshed to ensure
    //that the surface mesh is suitable for high-ordering

    if(m_mesh->m_verbose)
        cout << endl << "\tIntroducing high-order awareness to linear mesh" << endl;

    bool repeat = true;
    int splitedges = 0;
    while(repeat)
    {
        repeat = false;
        int edgesStart = m_mesh->m_edgeSet.size();
        EdgeSet edges = m_mesh->m_edgeSet;
        m_mesh->m_edgeSet.clear();

        EdgeSet::iterator it;

        int edgesplitthisrun = 0;

        set<int> alteredtri;

        for(it = edges.begin(); it != edges.end(); it++)
        {
            EdgeSharedPtr e = *it;

            if(e->CADSurfID.size() == 2)
            {
                //edge is on cad so dont interfer
                m_mesh->m_edgeSet.insert(e);
                continue;
            }

            ElementSharedPtr tri1 = e->m_elLink[0].first;
            ElementSharedPtr tri2 = e->m_elLink[1].first;

            set<int>::iterator fd1 = alteredtri.find(tri1->GetId());
            set<int>::iterator fd2 = alteredtri.find(tri2->GetId());

            if(!(fd1 == alteredtri.end() && fd2 == alteredtri.end()))
            {
                m_mesh->m_edgeSet.insert(e);
                continue;
            }

            CADSurfSharedPtr s = m_cad->GetSurf(e->CADSurfID[0]);
            int surf = e->CADSurfID[0];

            NodeSharedPtr n1 = e->m_n1;
            NodeSharedPtr n2 = e->m_n2;

            Array<OneD,NekDouble> N1,N2;
            N1 = s->N(n1->GetCADSurf(surf));
            N2 = s->N(n2->GetCADSurf(surf));

            NekDouble dot = N1[0]*N2[0] + N1[1]*N2[1] + N1[2]*N2[2];
            if(acos(dot) > 3.142/2.0-0.1)
            {
                repeat = true;
                splitedges++;
                edgesplitthisrun++;

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

                Array<OneD, NekDouble> ai,bi,ci,di;
                ai = A->GetCADSurf(surf);
                bi = B->GetCADSurf(surf);
                ci = C->GetCADSurf(surf);
                di = D->GetCADSurf(surf);

                //make the 4 other edges
                EdgeSharedPtr CA, AD, DB, BC, CAt, ADt, DBt, BCt;
                CAt = boost::shared_ptr<Edge>(new Edge(C,A));
                ADt = boost::shared_ptr<Edge>(new Edge(A,D));
                DBt = boost::shared_ptr<Edge>(new Edge(D,B));
                BCt = boost::shared_ptr<Edge>(new Edge(B,C));

                vector<EdgeSharedPtr> es = tri1->GetEdgeList();
                for(int i = 0; i < 3; i++)
                {
                    if(es[i] == CAt)
                    {
                        CA = es[i];
                    }
                    if(es[i] == BCt)
                    {
                        BC = es[i];
                    }
                }
                es = tri2->GetEdgeList();
                for(int i = 0; i < 3; i++)
                {
                    if(es[i] == DBt)
                    {
                        DB = es[i];
                    }
                    if(es[i] == ADt)
                    {
                        AD = es[i];
                    }
                }

                //now sort out links for the 4 edges surrounding the patch
                vector<pair<ElementSharedPtr, int> > links;

                links = CA->m_elLink;
                CA->m_elLink.clear();
                for(int i = 0; i < 2; i++)
                {
                    if(links[i].first->GetId() == tri1->GetId())
                        continue;
                    CA->m_elLink.push_back(links[i]);
                }

                links = BC->m_elLink;
                BC->m_elLink.clear();
                for(int i = 0; i < 2; i++)
                {
                    if(links[i].first->GetId() == tri1->GetId())
                        continue;
                    BC->m_elLink.push_back(links[i]);
                }

                links = AD->m_elLink;
                AD->m_elLink.clear();
                for(int i = 0; i < 2; i++)
                {
                    if(links[i].first->GetId() == tri2->GetId())
                        continue;
                    AD->m_elLink.push_back(links[i]);
                }

                links = DB->m_elLink;
                DB->m_elLink.clear();
                for(int i = 0; i < 2; i++)
                {
                    if(links[i].first->GetId() == tri2->GetId())
                        continue;
                    DB->m_elLink.push_back(links[i]);
                }

                //nodes identified element links removed

                Array<OneD, NekDouble> uve(2);
                uve[0] = (ai[0] + bi[0] + ci[0] + di[0])/4.0;
                uve[1] = (ai[1] + bi[1] + ci[1] + di[1])/4.0;
                Array<OneD, NekDouble> loc = s->P(uve);

                NodeSharedPtr E = boost::shared_ptr<Node>(new Node(0,loc[0],loc[1],loc[2]));
                E->SetCADSurf(surf,uve);

                EdgeSharedPtr CE = boost::shared_ptr<Edge>(new Edge(C,E));
                EdgeSharedPtr AE = boost::shared_ptr<Edge>(new Edge(A,E));
                EdgeSharedPtr EB = boost::shared_ptr<Edge>(new Edge(E,B));
                EdgeSharedPtr ED = boost::shared_ptr<Edge>(new Edge(E,D));

                vector<NodeSharedPtr> t1,t2,t3,t4;
                t1.push_back(C); t1.push_back(E); t1.push_back(A);
                t2.push_back(E); t2.push_back(D); t2.push_back(A);
                t3.push_back(C); t3.push_back(B); t3.push_back(E);
                t4.push_back(B); t4.push_back(D); t4.push_back(E);

                ElmtConfig conf(LibUtilities::eTriangle,1,false,false);
                vector<int> tags;
                tags.push_back(tri1->GetCADSurf());

                int id1 = tri1->GetId();
                int id2 = tri2->GetId();

                ElementSharedPtr ntri1 = GetElementFactory().
                            CreateInstance(LibUtilities::eTriangle,
                                           conf,t1,tags);
                ElementSharedPtr ntri2 = GetElementFactory().
                            CreateInstance(LibUtilities::eTriangle,
                                           conf,t2,tags);
                ElementSharedPtr ntri3 = GetElementFactory().
                            CreateInstance(LibUtilities::eTriangle,
                                           conf,t3,tags);
                ElementSharedPtr ntri4 = GetElementFactory().
                            CreateInstance(LibUtilities::eTriangle,
                                           conf,t4,tags);

                ntri1->SetCADSurf(tags[0]);
                ntri2->SetCADSurf(tags[0]);
                ntri3->SetCADSurf(tags[0]);
                ntri4->SetCADSurf(tags[0]);
                ntri1->SetId(id1);
                ntri2->SetId(id2);
                ntri3->SetId(m_mesh->m_element[2].size());
                ntri4->SetId(m_mesh->m_element[2].size()+1);

                AE->CADSurfID.push_back(tri1->GetCADSurf());
                CE->CADSurfID.push_back(tri1->GetCADSurf());
                EB->CADSurfID.push_back(tri1->GetCADSurf());
                ED->CADSurfID.push_back(tri1->GetCADSurf());

                vector<EdgeSharedPtr> t1es = ntri1->GetEdgeList();
                for(int i = 0; i < 3; i++)
                {
                    if(t1es[i] == CA)
                    {
                        ntri1->SetEdge(i,CA);
                        CA->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri1,i));
                    }
                    else if(t1es[i] == CE)
                    {
                        ntri1->SetEdge(i,CE);
                        CE->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri1,i));
                    }
                    else if(t1es[i] == AE)
                    {
                        ntri1->SetEdge(i,AE);
                        AE->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri1,i));
                    }
                    else
                    {
                        ASSERTL0(false,"weird edge in new tri 1");
                    }
                }
                vector<EdgeSharedPtr> t2es = ntri2->GetEdgeList();
                for(int i = 0; i < 3; i++)
                {
                    if(t2es[i] == AD)
                    {
                        ntri2->SetEdge(i,AD);
                        AD->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri2,i));
                    }
                    else if(t2es[i] == ED)
                    {
                        ntri2->SetEdge(i,ED);
                        ED->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri2,i));
                    }
                    else if(t2es[i] == AE)
                    {
                        ntri2->SetEdge(i,AE);
                        AE->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri2,i));
                    }
                    else
                    {
                        ASSERTL0(false,"weird edge in new tri 2");
                    }
                }
                vector<EdgeSharedPtr> t3es = ntri3->GetEdgeList();
                for(int i = 0; i < 3; i++)
                {
                    if(t3es[i] == BC)
                    {
                        ntri3->SetEdge(i,BC);
                        BC->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri3,i));
                    }
                    else if(t3es[i] == CE)
                    {
                        ntri3->SetEdge(i,CE);
                        CE->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri3,i));
                    }
                    else if(t3es[i] == EB)
                    {
                        ntri3->SetEdge(i,EB);
                        EB->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri3,i));
                    }
                    else
                    {
                        ASSERTL0(false,"weird edge in new tri 3");
                    }
                }
                vector<EdgeSharedPtr> t4es = ntri4->GetEdgeList();
                for(int i = 0; i < 3; i++)
                {
                    if(t4es[i] == DB)
                    {
                        ntri4->SetEdge(i,DB);
                        DB->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri4,i));
                    }
                    else if(t4es[i] == ED)
                    {
                        ntri4->SetEdge(i,ED);
                        ED->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri4,i));
                    }
                    else if(t4es[i] == EB)
                    {
                        ntri4->SetEdge(i,EB);
                        EB->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri4,i));
                    }
                    else
                    {
                        ASSERTL0(false,"weird edge in new tri 1");
                    }
                }

                m_mesh->m_element[2][id1] = ntri1;
                m_mesh->m_element[2][id2] = ntri2;
                m_mesh->m_element[2].push_back(ntri3);
                m_mesh->m_element[2].push_back(ntri4);

                alteredtri.insert(id1);
                alteredtri.insert(id2);
                alteredtri.insert(ntri3->GetId());
                alteredtri.insert(ntri4->GetId());


                m_mesh->m_edgeSet.insert(AE);
                m_mesh->m_edgeSet.insert(CE);
                m_mesh->m_edgeSet.insert(ED);
                m_mesh->m_edgeSet.insert(EB);
            }
            else
            {
                m_mesh->m_edgeSet.insert(e);
            }
        }

        ASSERTL0(m_mesh->m_edgeSet.size() == edgesStart+3*edgesplitthisrun, "mismatch edge count");
    }

    if(m_mesh->m_verbose)
        cout << "\t\tEdges split: " << splitedges << endl;
}*/

void SurfaceMesh::Report()
{
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
