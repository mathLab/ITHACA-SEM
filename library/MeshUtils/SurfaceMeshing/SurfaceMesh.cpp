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

void SurfaceMesh::Optimise()
{
    if(m_mesh->m_verbose)
        cout << endl << "\tOptimising linear mesh" << endl;

    //perfrom edge swap based on node defect and then angle
    for(int q = 0; q < 8; q++)
    {
        if(m_mesh->m_verbose)
        {
            cout << "\t\t Edge swap ";
            if(q<4)
            {
                cout << "defect run: " << q+1;
            }
            else
            {
                cout << "angle run: " << q+1-4;
            }
        }

        int edgesStart = m_mesh->m_edgeSet.size();
        EdgeSet edges = m_mesh->m_edgeSet;
        m_mesh->m_edgeSet.clear();

        int swappedEdges = 0;

        set<int> alteredtri;

        EdgeSet::iterator it;

        for(it = edges.begin(); it != edges.end(); it++)
        {
            EdgeSharedPtr e = *it;

            ElementSharedPtr tri1 = e->m_elLink[0].first;
            ElementSharedPtr tri2 = e->m_elLink[1].first;

            if(tri1->GetCADSurf() != tri2->GetCADSurf()) //edge on curve cannot swap
            {
                m_mesh->m_edgeSet.insert(e);
                continue;
            }

            set<int>::iterator fd1 = alteredtri.find(tri1->GetId());
            set<int>::iterator fd2 = alteredtri.find(tri2->GetId());

            if(!(fd1 == alteredtri.end() && fd2 == alteredtri.end()))
            {
                m_mesh->m_edgeSet.insert(e);
                continue;
            }

            int surf = tri1->GetCADSurf();

            NodeSharedPtr n1 = e->m_n1;
            NodeSharedPtr n2 = e->m_n2;

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

            //uncomment this section to not allow swapping of curve adjacent elments
            /*if(A->GetListCADCurve().size() > 0 ||
               B->GetListCADCurve().size() > 0 ||
               C->GetListCADCurve().size() > 0 ||
               D->GetListCADCurve().size() > 0 )
            {
                m_mesh->m_edgeSet.insert(e);
                continue;
            }*/

            NekDouble CBA, BDA;

            //determine signed area of alternate config
            Array<OneD, NekDouble> ai,bi,ci,di;
            ai = A->GetCADSurf(surf);
            bi = B->GetCADSurf(surf);
            ci = C->GetCADSurf(surf);
            di = D->GetCADSurf(surf);

            CBA = 0.5*(-bi[0]*ci[1] + ai[0]*ci[1] + ci[0]*bi[1] - ai[0]*bi[1] -
                            ci[0]*ai[1] + bi[0]*ai[1]);

            BDA = 0.5*(-di[0]*bi[1] + ai[0]*bi[1] + bi[0]*di[1] - ai[0]*di[1] -
                            bi[0]*ai[1] + di[0]*ai[1]);

            //shouldnt have any negative area triangles in set
            ASSERTL0(CBA > 1E-8 && BDA > 1E-8, "inverted triangle inputted to swap routine");

            NekDouble CDA, CBD;

            CDA = 0.5*(-di[0]*ci[1] + ai[0]*ci[1] + ci[0]*di[1] - ai[0]*di[1] -
                            ci[0]*ai[1] + di[0]*ai[1]);

            CBD = 0.5*(-bi[0]*ci[1] + di[0]*ci[1] + ci[0]*bi[1] - di[0]*bi[1] -
                            ci[0]*di[1] + bi[0]*di[1]);

            //if signed area of the swapping triangles is less than zero
            //that configuration is invalid and swap cannot be performed
            if(!(CDA > 0.001 && CBD > 0.001))
            {
                m_mesh->m_edgeSet.insert(e);
                continue;
            }

            bool swap = false; //assume do not swap

            if(q<4)
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
                if(nodedefectafter < nodedefectbefore)
                {
                    swap = true;
                }

            }
            else
            {
                NekDouble minanglebefore = C->Angle(A,B);
                minanglebefore = min(minanglebefore, A->Angle(B,C));
                minanglebefore = min(minanglebefore, B->Angle(A,C));
                minanglebefore = min(minanglebefore, B->Angle(A,D));
                minanglebefore = min(minanglebefore, A->Angle(B,D));
                minanglebefore = min(minanglebefore, D->Angle(A,B));

                NekDouble minangleafter = C->Angle(B,D);
                minangleafter = min(minangleafter, D->Angle(B,C));
                minangleafter = min(minangleafter, B->Angle(C,D));
                minangleafter = min(minangleafter, C->Angle(A,D));
                minangleafter = min(minangleafter, A->Angle(C,D));
                minangleafter = min(minangleafter, D->Angle(A,C));

                if(minangleafter > minanglebefore)
                {
                    swap = true;
                }
            }

            if(swap)
            {
                A->m_elCount--;
                B->m_elCount--;
                C->m_elCount++;
                D->m_elCount++;

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

                EdgeSharedPtr newe = boost::shared_ptr<Edge>(new Edge(C,D));
                newe->m_id = e->m_id;

                vector<NodeSharedPtr> t1,t2;
                t1.push_back(B); t1.push_back(D); t1.push_back(C);
                t2.push_back(A); t2.push_back(C); t2.push_back(D);

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

                ntri1->SetCADSurf(tags[0]);
                ntri2->SetCADSurf(tags[0]);
                ntri1->SetId(id1);
                ntri2->SetId(id2);

                vector<EdgeSharedPtr> t1es = ntri1->GetEdgeList();
                for(int i = 0; i < 3; i++)
                {
                    if(t1es[i] == DB)
                    {
                        ntri1->SetEdge(i,DB);
                        DB->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri1,i));
                    }
                    else if(t1es[i] == BC)
                    {
                        ntri1->SetEdge(i,BC);
                        BC->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri1,i));
                    }
                    else if(t1es[i] == newe)
                    {
                        ntri1->SetEdge(i,newe);
                        newe->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri1,i));
                    }
                    else
                    {
                        ASSERTL0(false,"weird edge in new tri 1");
                    }
                }
                vector<EdgeSharedPtr> t2es = ntri2->GetEdgeList();
                for(int i = 0; i < 3; i++)
                {
                    if(t2es[i] == CA)
                    {
                        ntri2->SetEdge(i,CA);
                        CA->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri2,i));
                    }
                    else if(t2es[i] == AD)
                    {
                        ntri2->SetEdge(i,AD);
                        AD->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri2,i));
                    }
                    else if(t2es[i] == newe)
                    {
                        ntri2->SetEdge(i,newe);
                        newe->m_elLink.push_back(pair<ElementSharedPtr,int>(ntri2,i));
                    }
                    else
                    {
                        ASSERTL0(false,"weird edge in new tri 2");
                    }
                }

                newe->CADSurfID.push_back(tri1->GetCADSurf());
                m_mesh->m_edgeSet.insert(newe);

                m_mesh->m_element[2][id1] = ntri1;
                m_mesh->m_element[2][id2] = ntri2;

                alteredtri.insert(id1);
                alteredtri.insert(id2);

                swappedEdges++;

            }
            else
            {
                m_mesh->m_edgeSet.insert(e);
            }
        }

        ASSERTL0(m_mesh->m_edgeSet.size() == edgesStart, "mismatch edge count");

        if(m_mesh->m_verbose)
            cout << ".\tEdges swapped: " << swappedEdges << endl;
    }

    EdgeSet::iterator eit;
    for(eit = m_mesh->m_edgeSet.begin(); eit != m_mesh->m_edgeSet.end(); eit++)
    {
        if((*eit)->m_elLink.size() !=2)
            continue;

        (*eit)->m_n1->m_connectingNodes.push_back((*eit)->m_n2);
        (*eit)->m_n2->m_connectingNodes.push_back((*eit)->m_n1);
    }

    NodeSet::iterator nit;
    //perform 8 runs of elastic relaxation based on the octree
    for(int q = 0; q < 8; q++)
    {
        if(m_mesh->m_verbose)
            cout << "\t\t Elastic relaxation run: " << q+1 << endl;

        for(nit = m_mesh->m_vertexSet.begin(); nit != m_mesh->m_vertexSet.end(); nit++)
        {
            vector<int> c = (*nit)->GetListCADCurve();
            if(c.size()>0) //node is on curve so skip
                continue;

            vector<int> surfs = (*nit)->GetListCADSurf();
            ASSERTL0(surfs.size()==1, //idiot checking
                        "node should be interior and only be on one surface");

            int surf = surfs[0];
            CADSurfSharedPtr s = m_cad->GetSurf(surf);

            vector<NodeSharedPtr> connodes = (*nit)->m_connectingNodes;

            //figure out the convexity of the connodes system using graham scan
            NodeSharedPtr lowest;
            vector<NodeSharedPtr> orderedNodes;

            //find lowest v node
            lowest = connodes[0];
            for(int i = 1; i < connodes.size(); i++)
            {
                Array<OneD, NekDouble> uvlow = lowest->GetCADSurf(surf);
                Array<OneD, NekDouble> uvtest = connodes[i]->GetCADSurf(surf);
                if(uvtest[1] < uvlow[1])
                {
                    lowest = connodes[i];
                }
            }
            //build unordered list of others
            for(int i = 0; i < connodes.size(); i++)
            {
                if(connodes[i] == lowest)
                    continue;

                orderedNodes.push_back(connodes[i]);
            }

            Array<OneD, NekDouble> uvlow = lowest->GetCADSurf(surf);

            vector<NekDouble> angles;
            for(int i = 0; i < orderedNodes.size(); i++)
            {
                Array<OneD, NekDouble> uv = orderedNodes[i]->GetCADSurf(surf);
                Array<OneD, NekDouble> cs1(2), cs2(2), cn1(2);
                cs1[0]   = 1.0;
                cs1[1]   = 0.0;
                cs2[0]   = uv[0]-uvlow[0];
                cs2[1]   = uv[1]-uvlow[1];
                cn1[0]   = -cs1[1];
                cn1[1]   =  cs1[0];
                NekDouble an       = cn1[0]*cn1[0]+cn1[1]*cn1[1];
                an       = 1.0/sqrt(an);
                cn1[0]   = cn1[0]*an;
                cn1[1]   = cn1[1]*an;
                an       = cs1[0]*cs1[0]+cs1[1]*cs1[1];
                an       = 1.0/sqrt(an);
                cs1[0]   = cs1[0]*an;
                cs1[1]   = cs1[1]*an;
                an       = cs2[0]*cs2[0]+cs2[1]*cs2[1];
                an       = 1.0/sqrt(an);
                cs2[0]   = cs2[0]*an;
                cs2[1]   = cs2[1]*an;
                NekDouble cosw     = cs1[0]*cs2[0]+cs1[1]*cs2[1];
                NekDouble sinw     = cs2[0]*cn1[0]+cs2[1]*cn1[1];
                angles.push_back(atan2(sinw,cosw));
            }

            //sort the orderedNodes based on the angles
            bool repeat = true;
            while(repeat)
            {
                repeat = false;
                for(int i = 0; i < orderedNodes.size() -1; i++)
                {
                    if(angles[i+1] < angles[i])
                    {
                        NodeSharedPtr tmpn = orderedNodes[i];
                        NekDouble tmpa = angles[i];
                        orderedNodes[i] = orderedNodes[i+1];
                        angles[i] = angles[i+1];
                        angles[i+1] = tmpa;
                        orderedNodes[i+1] = tmpn;
                        repeat = true;
                    }
                }
            }

            bool concave = false;

            if(!(orderedNodes.size()>2))
                continue;

            Array<OneD, NekDouble> uva = orderedNodes[0]->GetCADSurf(surf);
            Array<OneD, NekDouble> uvb = orderedNodes[1]->GetCADSurf(surf);
            if((uvlow[1]-uva[1])*(uvb[0]-uva[0]) - (uvlow[0]-uva[0])*(uvb[1]-uva[1]) < 0)
            {
                concave = true;
            }
            else
            {
                for(int i = 1; i < orderedNodes.size() - 1; i++)
                {
                    Array<OneD, NekDouble> uvc = uva;
                    uva = uvb;
                    uvb = orderedNodes[i+1]->GetCADSurf(surf);
                    if((uvc[1]-uva[1])*(uvb[0]-uva[0]) - (uvc[0]-uva[0])*(uvb[1]-uva[1]) < 0)
                    {
                        concave = true;
                        break;
                    }
                }
                if(concave == false) //test last combo
                {
                    Array<OneD, NekDouble> uvc = uva;
                    uva = uvb;
                    uvb = uvlow;
                    if((uvc[1]-uva[1])*(uvb[0]-uva[0]) - (uvc[0]-uva[0])*(uvb[1]-uva[1]) < 0)
                    {
                        concave = true;
                    }
                }
            }

            if(concave)
            {
                continue;
            }

            Array<OneD, NekDouble> uv0(2);
            uv0[0]=0.0; uv0[1]=0.0;

            DNekMat f(2,1,0.0);
            DNekMat df(2,2,0.0);
            for(int i = 0; i < connodes.size(); i++)
            {
                Array<OneD, NekDouble> uvj = connodes[i]->GetCADSurf(surf);
                uv0[0]+=uvj[0]/connodes.size();
                uv0[1]+=uvj[1]/connodes.size();
            }

            Array<OneD, NekDouble> rui = m_cad->GetSurf(surf)->P(uv0);
            NekDouble d = m_octree->Query(rui);
            for(int i = 0; i < connodes.size();i++)
            {
                Array<OneD, NekDouble> uvj = connodes[i]->GetCADSurf(surf);
                Array<OneD, NekDouble> rj  = connodes[i]->GetLoc();

                NekDouble difR = sqrt((rui[0]-rj[0])*(rui[0]-rj[0]) +
                                      (rui[1]-rj[1])*(rui[1]-rj[1]) +
                                      (rui[2]-rj[2])*(rui[2]-rj[2]));

                NekDouble difU = sqrt((uvj[0]-uv0[0])*(uvj[0]-uv0[0]) +
                                      (uvj[1]-uv0[1])*(uvj[1]-uv0[1]));

                NekDouble A    = difR - d;

                NekDouble B    = (uvj[0]-uv0[0]) / difU;

                NekDouble C    = (uvj[1]-uv0[1]) / difU;

                Array<OneD, NekDouble> r = s->D1(uv0);

                NekDouble dAdu = ((r[0]-rj[0])*r[3] +
                                  (r[1]-rj[1])*r[4] +
                                  (r[2]-rj[2])*r[5]) / difR;

                NekDouble dAdv = ((r[0]-rj[0])*r[6] +
                                  (r[1]-rj[1])*r[7] +
                                  (r[2]-rj[2])*r[8]) / difR;

                NekDouble dBdu = (uvj[0]-uv0[0])*(uvj[0]-uv0[0])/
                                 difU/difU/difU - 1.0/difU;

                NekDouble dBdv = (uvj[0]-uv0[0])*(uvj[1]-uv0[1])/
                                 difU/difU/difU;

                NekDouble dCdv = (uvj[1]-uv0[1])*(uvj[1]-uv0[1])/
                                 difU/difU/difU - 1.0/difU;

                NekDouble dCdu = dBdv;

                f(0,0) += A*B;
                f(1,0) += A*C;

                df(0,0) += B*dAdu + A*dBdu;
                df(1,0) += C*dAdu + A*dCdu;
                df(0,1) += B*dAdv + A*dBdv;
                df(1,1) += C*dAdv + A*dCdv;
            }

            df.Invert();

            DNekMat ui = df*f;
            Array<OneD, NekDouble> uvn(2);
            uvn[0] = uv0[0];// - ui(0,0);
            uvn[1] = uv0[1];// - ui(1,0);

            Array<OneD, NekDouble> bounds = m_cad->GetSurf(surf)->GetBounds();
            Array<OneD, NekDouble> uvi = (*nit)->GetCADSurf(surf);

            if(!(uvn[0] < bounds[0] ||
                       uvn[0] > bounds[1] ||
                       uvn[1] < bounds[2] ||
                       uvn[1] > bounds[3]))
            {

                Array<OneD, NekDouble> l2 = m_cad->GetSurf(surf)->P(uvn);
                (*nit)->Move(l2,surf,uvn);
            }
        }
    }
}

void SurfaceMesh::HOAwareness()
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
