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
#include <NekMeshUtils/SurfaceMeshing/FaceMesh.h>
#include <NekMeshUtils/ExtLibInterface/TriangleInterface.h>

#include <LocalRegions/MatrixKey.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

void FaceMesh::Mesh()
{
    Stretching();

    OrientateCurves();

    if (m_makebl)
    {
        MakeBL();
    }

    int numPoints = 0;
    for (int i = 0; i < orderedLoops.size(); i++)
    {
        numPoints += orderedLoops[i].size();
    }

    stringstream ss;
    ss << "3 points required for triangulation, " << numPoints << " in loop"
       << endl;
    ss << "curves: ";
    for (int i = 0; i < m_edgeloops.size(); i++)
    {
        for (int j = 0; j < m_edgeloops[i].edges.size(); j++)
        {
            ss << m_edgeloops[i].edges[j]->GetId() << " ";
        }
    }

    ASSERTL0(numPoints > 2, ss.str());

    // create interface to triangle thirdparty library
    TriangleInterfaceSharedPtr pplanemesh =
        MemoryManager<TriangleInterface>::AllocateSharedPtr();

    vector<Array<OneD, NekDouble> > centers;
    for (int i = 0; i < m_edgeloops.size(); i++)
    {
        centers.push_back(m_edgeloops[i].center);
    }

    pplanemesh->Assign(orderedLoops, centers, m_id, m_str);

    pplanemesh->Mesh();

    pplanemesh->Extract(m_connec);

    bool repeat     = true;
    int meshcounter = 1;

    while (repeat)
    {
        repeat = Validate();
        if (!repeat)
        {
            break;
        }
        m_connec.clear();
        pplanemesh->AssignStiener(m_stienerpoints);
        pplanemesh->Mesh();
        pplanemesh->Extract(m_connec);
        meshcounter++;
    }

    BuildLocalMesh();

    OptimiseLocalMesh();

    // clear local element links
    EdgeSet::iterator eit;
    for (eit = m_localEdges.begin(); eit != m_localEdges.end(); eit++)
    {
        (*eit)->m_elLink.clear();
    }

    // make new elements and add to list from list of nodes and connectivity
    // from triangle
    for (int i = 0; i < m_localElements.size(); i++)
    {
        vector<EdgeSharedPtr> e = m_localElements[i]->GetEdgeList();
        for (int j = 0; j < e.size(); j++)
        {
            e[j]->m_elLink.clear();
        }
        m_mesh->m_element[m_mesh->m_expDim].push_back(m_localElements[i]);
    }

    cout << "\r                                                                "
            "                             ";
    cout << scientific << "\r\t\tFace " << m_id << endl
         << "\t\t\tNodes: " << m_localNodes.size() << endl
         << "\t\t\tEdges: " << m_localEdges.size() << endl
         << "\t\t\tTriangles: " << m_localElements.size() << endl
         << "\t\t\tLoops: " << m_edgeloops.size() << endl
         << "\t\t\tSTR: " << m_str << endl
         << endl;
}

void FaceMesh::MakeBL()
{
    for (int i = 1; i < orderedLoops.size(); i++) // dont do this to first loop
    {
        for (int j = 0; j < orderedLoops[i].size(); j++)
        {
            // for each of the nodes make a new node which exists off the
            // surface
            vector<pair<int, CADSurfSharedPtr> > surfs =
                orderedLoops[i][j]->GetCADSurfs();
            ASSERTL0(surfs.size() > 1,
                     "point does not have enough surfs to make quad");

            Array<OneD, NekDouble> AN(3, 0.0);
            // make a averaged normal ignoring surface mid
            for (int s = 0; s < surfs.size(); s++)
            {
                if (surfs[s].first == m_id)
                    continue; // does not contribute to norm

                Array<OneD, NekDouble> uv =
                    orderedLoops[i][j]->GetCADSurfInfo(surfs[s].first);
                Array<OneD, NekDouble> N = surfs[s].second->N(uv);
                for (int k = 0; k < 3; k++)
                    AN[k] += N[k];
            }

            // renormalise normal
            NekDouble mag = sqrt(AN[0] * AN[0] + AN[1] * AN[1] + AN[2] * AN[2]);
            for (int k = 0; k < 3; k++)
                AN[k] /= mag;

            Array<OneD, NekDouble> loc = orderedLoops[i][j]->GetLoc();
            Array<OneD, NekDouble> tp(3);
            for (int k = 0; k < 3; k++)
                tp[k] = m_bl * AN[k] + loc[k];
            // project tp onto to surface to get new point
            Array<OneD, NekDouble> uv(2);
            m_cadsurf->ProjectTo(tp, uv);
            NodeSharedPtr nn = boost::shared_ptr<Node>(
                new Node(m_mesh->m_numNodes++, tp[0], tp[1], tp[2]));
            nn->SetCADSurf(m_id, m_cadsurf, uv);

            blpairs.push_back(
                pair<NodeSharedPtr, NodeSharedPtr>(orderedLoops[i][j], nn));
            // place the new node into ordered loops for the surface
            // triangulation to work With
            orderedLoops[i][j] = nn;
        }
    }
}

void FaceMesh::OptimiseLocalMesh()
{
    DiagonalSwap();

    Smoothing();

    DiagonalSwap();

    Smoothing();
}

void FaceMesh::Smoothing()
{
    EdgeSet::iterator eit;
    NodeSet::iterator nit;

    map<int, vector<EdgeSharedPtr> > connectingedges;

    map<int, vector<ElementSharedPtr> > connectingelements;

    for (eit = m_localEdges.begin(); eit != m_localEdges.end(); eit++)
    {
        connectingedges[(*eit)->m_n1->m_id].push_back(*eit);
        connectingedges[(*eit)->m_n2->m_id].push_back(*eit);
    }

    for (int i = 0; i < m_localElements.size(); i++)
    {
        vector<NodeSharedPtr> v = m_localElements[i]->GetVertexList();
        for (int j = 0; j < 3; j++)
        {
            connectingelements[v[j]->m_id].push_back(m_localElements[i]);
        }
    }

    // perform 4 runs of elastic relaxation based on the octree
    for (int q = 0; q < 4; q++)
    {
        for (nit = m_localNodes.begin(); nit != m_localNodes.end(); nit++)
        {
            if ((*nit)->GetNumCadCurve() > 0) // node is on curve so skip
                continue;

            vector<NodeSharedPtr> connodes; // this can be real nodes or dummy
                                            // nodes depending on the system

            vector<EdgeSharedPtr> edges  = connectingedges[(*nit)->m_id];
            vector<ElementSharedPtr> els = connectingelements[(*nit)->m_id];

            bool perfrom = true;
            for (int i = 0; i < els.size(); i++)
            {
                if (els[i]->GetConf().m_e == LibUtilities::eQuadrilateral)
                {
                    perfrom = false;
                    break;
                }
            }

            if (!perfrom)
                continue;

            vector<NodeSharedPtr> nodesystem;
            vector<NekDouble> lamp;

            for (int i = 0; i < edges.size(); i++)
            {
                vector<NekDouble> lambda;

                NodeSharedPtr J;
                if (*nit == edges[i]->m_n1)
                    J = edges[i]->m_n2;
                else if (*nit == edges[i]->m_n2)
                    J = edges[i]->m_n1;
                else
                    ASSERTL0(false, "could not find node");

                Array<OneD, NekDouble> ui = (*nit)->GetCADSurfInfo(m_id);
                Array<OneD, NekDouble> uj = J->GetCADSurfInfo(m_id);

                for (int j = 0; j < els.size(); j++)
                {
                    vector<NodeSharedPtr> v = els[j]->GetVertexList();
                    if (v[0] == J || v[1] == J || v[2] == J)
                        continue; // elememt is adjacent to J therefore no
                                  // intersection on IJ

                    // need to find other edge
                    EdgeSharedPtr AtoB;
                    bool found               = false;
                    vector<EdgeSharedPtr> es = els[j]->GetEdgeList();
                    for (int k = 0; k < es.size(); k++)
                    {
                        if (!(es[k]->m_n1 == *nit || es[k]->m_n2 == *nit))
                        {
                            found = true;
                            AtoB  = es[k];
                            break;
                        }
                    }
                    ASSERTL0(found, "failed to find edge to test");

                    Array<OneD, NekDouble> A = AtoB->m_n1->GetCADSurfInfo(m_id);
                    Array<OneD, NekDouble> B = AtoB->m_n2->GetCADSurfInfo(m_id);

                    NekDouble lam = ((A[0] - uj[0]) * (B[1] - A[1]) -
                                     (A[1] - uj[1]) * (B[0] - A[0])) /
                                    ((ui[0] - uj[0]) * (B[1] - A[1]) -
                                     (ui[1] - uj[1]) * (B[0] - A[0]));

                    if (!(lam < 0) && !(lam > 1))
                        lambda.push_back(lam);
                }

                if (lambda.size() > 0)
                {
                    sort(lambda.begin(), lambda.end());
                    // make a new dummy node based on the system
                    Array<OneD, NekDouble> ud(2);
                    ud[0]                       = uj[0] + lambda[0] * (ui[0] - uj[0]);
                    ud[1]                       = uj[1] + lambda[0] * (ui[1] - uj[1]);
                    Array<OneD, NekDouble> locd = m_cadsurf->P(ud);
                    NodeSharedPtr dn = boost::shared_ptr<Node>(
                        new Node(0, locd[0], locd[1], locd[2]));
                    dn->SetCADSurf(m_id, m_cadsurf, ud);

                    nodesystem.push_back(dn);
                    lamp.push_back(lambda[0]);
                }
                else
                {
                    nodesystem.push_back(J);
                    lamp.push_back(1.0);
                }
            }

            Array<OneD, NekDouble> ui(2);
            ui[0] = 0.0;
            ui[1] = 0.0;

            for (int i = 0; i < nodesystem.size(); i++)
            {
                Array<OneD, NekDouble> uj = nodesystem[i]->GetCADSurfInfo(m_id);
                ui[0] += uj[0] / nodesystem.size();
                ui[1] += uj[1] / nodesystem.size();
            }

            Array<OneD, NekDouble> bounds = m_cadsurf->GetBounds();

            Array<OneD, NekDouble> uvn(2);

            uvn[0] = ui[0];
            uvn[1] = ui[1];

            if (!(uvn[0] < bounds[0] || uvn[0] > bounds[1] ||
                  uvn[1] < bounds[2] || uvn[1] > bounds[3]))
            {
                Array<OneD, NekDouble> l2 = m_cadsurf->P(uvn);
                (*nit)->Move(l2, m_id, uvn);
            }
        }
    }
}

void FaceMesh::DiagonalSwap()
{
    /// TODO fix this bit of code which figures out the node defect, doesnt work
    /// on quads or relfex angles
    map<int, int> idealConnec;
    map<int, int> actualConnec;
    map<int, vector<EdgeSharedPtr> > nodetoedge;
    // figure out ideal node count and actual node count
    EdgeSet::iterator eit;
    for (eit = m_localEdges.begin(); eit != m_localEdges.end(); eit++)
    {
        nodetoedge[(*eit)->m_n1->m_id].push_back(*eit);
        nodetoedge[(*eit)->m_n2->m_id].push_back(*eit);
    }
    NodeSet::iterator nit;
    for (nit = m_localNodes.begin(); nit != m_localNodes.end(); nit++)
    {
        if ((*nit)->GetNumCadCurve() == 0)
        {
            // node is interior
            idealConnec[(*nit)->m_id] = 6;
        }
        else
        {
            // need to identify the two other nodes on the boundary to find
            // interior angle
            vector<NodeSharedPtr> ns;
            vector<EdgeSharedPtr> e = nodetoedge[(*nit)->m_id];
            for (int i = 0; i < e.size(); i++)
            {
                if (!e[i]->onCurve)
                    continue; // the linking nodes are not going to exist on
                              // interior edges

                if (e[i]->m_n1 == (*nit))
                    ns.push_back(e[i]->m_n2);
                else
                    ns.push_back(e[i]->m_n1);
            }
            ASSERTL0(ns.size() == 2,
                     "failed to find 2 nodes in the angle system");

            idealConnec[(*nit)->m_id] =
                ceil((*nit)->Angle(ns[0], ns[1]) / 3.142 * 3) + 1;
        }
    }
    for (nit = m_localNodes.begin(); nit != m_localNodes.end(); nit++)
    {
        actualConnec[(*nit)->m_id] = nodetoedge[(*nit)->m_id].size();
    }

    // edgeswapping fun times
    // perfrom edge swap based on node defect and then angle
    for (int q = 0; q < 4; q++)
    {
        int edgesStart = m_localEdges.size();
        EdgeSet edges = m_localEdges;
        m_localEdges.clear();

        int swappedEdges = 0;

        EdgeSet::iterator it;

        for (it = edges.begin(); it != edges.end(); it++)
        {
            EdgeSharedPtr e = *it;

            if (e->m_elLink.size() != 2)
            {
                m_localEdges.insert(e);
                continue;
            }
            if (e->m_elLink[0].first->GetConf().m_e ==
                    LibUtilities::eQuadrilateral ||
                e->m_elLink[1].first->GetConf().m_e ==
                    LibUtilities::eQuadrilateral)
            {
                m_localEdges.insert(e);
                continue;
            }

            ElementSharedPtr tri1 = e->m_elLink[0].first;
            ElementSharedPtr tri2 = e->m_elLink[1].first;

            NodeSharedPtr n1 = e->m_n1;
            NodeSharedPtr n2 = e->m_n2;

            vector<NodeSharedPtr> nt = tri1->GetVertexList();

            // identify node a,b,c,d of the swapping
            NodeSharedPtr A, B, C, D;
            if (nt[0] != n1 && nt[0] != n2)
            {
                C = nt[0];
                B = nt[1];
                A = nt[2];
            }
            else if (nt[1] != n1 && nt[1] != n2)
            {
                C = nt[1];
                B = nt[2];
                A = nt[0];
            }
            else if (nt[2] != n1 && nt[2] != n2)
            {
                C = nt[2];
                B = nt[0];
                A = nt[1];
            }
            else
            {
                ASSERTL0(false, "failed to identify verticies in tri1");
            }

            nt = tri2->GetVertexList();

            if (nt[0] != n1 && nt[0] != n2)
            {
                D = nt[0];
            }
            else if (nt[1] != n1 && nt[1] != n2)
            {
                D = nt[1];
            }
            else if (nt[2] != n1 && nt[2] != n2)
            {
                D = nt[2];
            }
            else
            {
                ASSERTL0(false, "failed to identify verticies in tri2");
            }

            // determine signed area of alternate config
            Array<OneD, NekDouble> ai, bi, ci, di;
            ai = A->GetCADSurfInfo(m_id);
            bi = B->GetCADSurfInfo(m_id);
            ci = C->GetCADSurfInfo(m_id);
            di = D->GetCADSurfInfo(m_id);

            NekDouble CDA, CBD;

            CDA = 0.5 * (-di[0] * ci[1] + ai[0] * ci[1] + ci[0] * di[1] -
                         ai[0] * di[1] - ci[0] * ai[1] + di[0] * ai[1]);

            CBD = 0.5 * (-bi[0] * ci[1] + di[0] * ci[1] + ci[0] * bi[1] -
                         di[0] * bi[1] - ci[0] * di[1] + bi[0] * di[1]);

            // if signed area of the swapping triangles is less than zero
            // that configuration is invalid and swap cannot be performed
            if (!(CDA > 0.001 && CBD > 0.001))
            {
                m_localEdges.insert(e);
                continue;
            }

            bool swap = false; // assume do not swap

            if (q < 2)
            {
                int nodedefectbefore = 0;
                nodedefectbefore +=
                    abs(actualConnec[A->m_id] - idealConnec[A->m_id]);
                nodedefectbefore +=
                    abs(actualConnec[B->m_id] - idealConnec[B->m_id]);
                nodedefectbefore +=
                    abs(actualConnec[C->m_id] - idealConnec[C->m_id]);
                nodedefectbefore +=
                    abs(actualConnec[D->m_id] - idealConnec[D->m_id]);

                int nodedefectafter = 0;
                nodedefectafter +=
                    abs(actualConnec[A->m_id] - 1 - idealConnec[A->m_id]);
                nodedefectafter +=
                    abs(actualConnec[B->m_id] - 1 - idealConnec[B->m_id]);
                nodedefectafter +=
                    abs(actualConnec[C->m_id] + 1 - idealConnec[C->m_id]);
                nodedefectafter +=
                    abs(actualConnec[D->m_id] + 1 - idealConnec[D->m_id]);

                if (nodedefectafter < nodedefectbefore)
                {
                    swap = true;
                }
            }
            else
            {
                NekDouble minanglebefore = C->Angle(A, B);
                minanglebefore           = min(minanglebefore, A->Angle(B, C));
                minanglebefore           = min(minanglebefore, B->Angle(A, C));
                minanglebefore           = min(minanglebefore, B->Angle(A, D));
                minanglebefore           = min(minanglebefore, A->Angle(B, D));
                minanglebefore           = min(minanglebefore, D->Angle(A, B));

                NekDouble minangleafter = C->Angle(B, D);
                minangleafter           = min(minangleafter, D->Angle(B, C));
                minangleafter           = min(minangleafter, B->Angle(C, D));
                minangleafter           = min(minangleafter, C->Angle(A, D));
                minangleafter           = min(minangleafter, A->Angle(C, D));
                minangleafter           = min(minangleafter, D->Angle(A, C));

                if (minangleafter > minanglebefore)
                {
                    swap = true;
                }
            }

            if (swap)
            {
                actualConnec[A->m_id]--;
                actualConnec[B->m_id]--;
                actualConnec[C->m_id]++;
                actualConnec[D->m_id]++;

                // make the 4 other edges
                EdgeSharedPtr CA, AD, DB, BC, CAt, ADt, DBt, BCt;
                CAt = boost::shared_ptr<Edge>(new Edge(C, A));
                ADt = boost::shared_ptr<Edge>(new Edge(A, D));
                DBt = boost::shared_ptr<Edge>(new Edge(D, B));
                BCt = boost::shared_ptr<Edge>(new Edge(B, C));

                vector<EdgeSharedPtr> es = tri1->GetEdgeList();
                for (int i = 0; i < 3; i++)
                {
                    if (es[i] == CAt)
                    {
                        CA = es[i];
                    }
                    if (es[i] == BCt)
                    {
                        BC = es[i];
                    }
                }
                es = tri2->GetEdgeList();
                for (int i = 0; i < 3; i++)
                {
                    if (es[i] == DBt)
                    {
                        DB = es[i];
                    }
                    if (es[i] == ADt)
                    {
                        AD = es[i];
                    }
                }

                // now sort out links for the 4 edges surrounding the patch
                vector<pair<ElementSharedPtr, int> > links;

                links = CA->m_elLink;
                CA->m_elLink.clear();
                for (int i = 0; i < links.size(); i++)
                {
                    if (links[i].first->GetId() == tri1->GetId())
                        continue;
                    CA->m_elLink.push_back(links[i]);
                }

                links = BC->m_elLink;
                BC->m_elLink.clear();
                for (int i = 0; i < links.size(); i++)
                {
                    if (links[i].first->GetId() == tri1->GetId())
                        continue;
                    BC->m_elLink.push_back(links[i]);
                }

                links = AD->m_elLink;
                AD->m_elLink.clear();
                for (int i = 0; i < links.size(); i++)
                {
                    if (links[i].first->GetId() == tri2->GetId())
                        continue;
                    AD->m_elLink.push_back(links[i]);
                }

                links = DB->m_elLink;
                DB->m_elLink.clear();
                for (int i = 0; i < links.size(); i++)
                {
                    if (links[i].first->GetId() == tri2->GetId())
                        continue;
                    DB->m_elLink.push_back(links[i]);
                }

                EdgeSharedPtr newe = boost::shared_ptr<Edge>(new Edge(C, D));

                vector<NodeSharedPtr> t1, t2;
                t1.push_back(B);
                t1.push_back(D);
                t1.push_back(C);
                t2.push_back(A);
                t2.push_back(C);
                t2.push_back(D);

                ElmtConfig conf(LibUtilities::eTriangle, 1, false, false);
                vector<int> tags = tri1->GetTagList();

                int id1 = tri1->GetId();
                int id2 = tri2->GetId();

                ElementSharedPtr ntri1 = GetElementFactory().CreateInstance(
                    LibUtilities::eTriangle, conf, t1, tags);
                tags                   = tri2->GetTagList();
                ElementSharedPtr ntri2 = GetElementFactory().CreateInstance(
                    LibUtilities::eTriangle, conf, t2, tags);

                ntri1->SetId(id1);
                ntri2->SetId(id2);
                ntri1->CADSurfId = m_id;
                ntri2->CADSurfId = m_id;

                vector<EdgeSharedPtr> t1es = ntri1->GetEdgeList();
                for (int i = 0; i < 3; i++)
                {
                    if (t1es[i] == DB)
                    {
                        ntri1->SetEdge(i, DB);
                        DB->m_elLink.push_back(
                            pair<ElementSharedPtr, int>(ntri1, i));
                    }
                    else if (t1es[i] == BC)
                    {
                        ntri1->SetEdge(i, BC);
                        BC->m_elLink.push_back(
                            pair<ElementSharedPtr, int>(ntri1, i));
                    }
                    else if (t1es[i] == newe)
                    {
                        ntri1->SetEdge(i, newe);
                        newe->m_elLink.push_back(
                            pair<ElementSharedPtr, int>(ntri1, i));
                    }
                    else
                    {
                        ASSERTL0(false, "weird edge in new tri 1");
                    }
                }
                vector<EdgeSharedPtr> t2es = ntri2->GetEdgeList();
                for (int i = 0; i < 3; i++)
                {
                    if (t2es[i] == CA)
                    {
                        ntri2->SetEdge(i, CA);
                        CA->m_elLink.push_back(
                            pair<ElementSharedPtr, int>(ntri2, i));
                    }
                    else if (t2es[i] == AD)
                    {
                        ntri2->SetEdge(i, AD);
                        AD->m_elLink.push_back(
                            pair<ElementSharedPtr, int>(ntri2, i));
                    }
                    else if (t2es[i] == newe)
                    {
                        ntri2->SetEdge(i, newe);
                        newe->m_elLink.push_back(
                            pair<ElementSharedPtr, int>(ntri2, i));
                    }
                    else
                    {
                        ASSERTL0(false, "weird edge in new tri 2");
                    }
                }

                m_localEdges.insert(newe);

                m_localElements[id1] = ntri1;
                m_localElements[id2] = ntri2;

                swappedEdges++;
            }
            else
            {
                m_localEdges.insert(e);
            }
        }

        ASSERTL0(m_localEdges.size() == edgesStart, "mismatch edge count");
    }
}

void FaceMesh::BuildLocalMesh()
{
    /*************************
    // build a local set of nodes edges and elemenets for optimstaion prior to
    putting them into m_mesh
    */
    int tricomp = m_mesh->m_numcomp++;

    // first build quads is bl surface
    if (m_makebl)
    {
        int quadcomp = m_mesh->m_numcomp++;
        for (int i = 0; i < blpairs.size() - 1;
             i++) // this wont work for more than one sym plane loop
        {
            ElmtConfig conf(LibUtilities::eQuadrilateral, 1, false, false);
            vector<NodeSharedPtr> ns;
            ns.push_back(blpairs[i].first);
            ns.push_back(blpairs[i].second);
            ns.push_back(blpairs[i + 1].second);
            ns.push_back(blpairs[i + 1].first);
            vector<int> tags;
            tags.push_back(quadcomp);
            ElementSharedPtr E = GetElementFactory().CreateInstance(
                LibUtilities::eQuadrilateral, conf, ns, tags);
            E->CADSurfId = m_id;
            m_localElements.push_back(E);
        }
        {
            ElmtConfig conf(LibUtilities::eQuadrilateral, 1, false, false);
            vector<NodeSharedPtr> ns;
            ns.push_back(blpairs.back().first);
            ns.push_back(blpairs.back().second);
            ns.push_back(blpairs[0].second);
            ns.push_back(blpairs[0].first);
            vector<int> tags;
            tags.push_back(quadcomp);
            ElementSharedPtr E = GetElementFactory().CreateInstance(
                LibUtilities::eQuadrilateral, conf, ns, tags);
            E->CADSurfId = m_id;
            m_localElements.push_back(E);
        }
        for (int i = 0; i < m_localElements.size(); i++)
        {
            vector<NodeSharedPtr> nods = m_localElements[i]->GetVertexList();
            for (int j = 0; j < nods.size(); j++)
            {
                // nodes are already unique some will insert some wont
                m_localNodes.insert(nods[j]);
            }
            vector<EdgeSharedPtr> edgs = m_localElements[i]->GetEdgeList();
            for (int j = 0; j < edgs.size(); j++)
            {
                // look for edge in m_mesh edgeset from curves
                EdgeSet::iterator s = m_mesh->m_edgeSet.find(edgs[j]);
                if (!(s == m_mesh->m_edgeSet.end()))
                {
                    edgs[j] = *s;
                    m_localElements[i]->SetEdge(j, edgs[j]);
                }

                pair<EdgeSet::iterator, bool> test =
                    m_localEdges.insert(edgs[j]);

                if (test.second)
                {
                    (*test.first)
                        ->m_elLink.push_back(
                            pair<ElementSharedPtr, int>(m_localElements[i], j));
                }
                else
                {
                    m_localElements[i]->SetEdge(j, *test.first);
                    (*test.first)
                        ->m_elLink.push_back(
                            pair<ElementSharedPtr, int>(m_localElements[i], j));
                }
            }
            m_localElements[i]->SetId(i);
        }
    }

    for (int i = 0; i < m_connec.size(); i++)
    {
        ElmtConfig conf(LibUtilities::eTriangle, 1, false, false);

        vector<int> tags;
        tags.push_back(tricomp);
        ElementSharedPtr E = GetElementFactory().CreateInstance(
            LibUtilities::eTriangle, conf, m_connec[i], tags);
        E->CADSurfId = m_id;

        vector<NodeSharedPtr> nods = E->GetVertexList();
        for (int j = 0; j < nods.size(); j++)
        {
            // nodes are already unique some will insert some wont
            m_localNodes.insert(nods[j]);
        }
        vector<EdgeSharedPtr> edgs = E->GetEdgeList();
        for (int j = 0; j < edgs.size(); j++)
        {
            // look for edge in m_mesh edgeset from curves
            EdgeSet::iterator s = m_mesh->m_edgeSet.find(edgs[j]);
            if (!(s == m_mesh->m_edgeSet.end()))
            {
                edgs[j] = *s;
                E->SetEdge(j, edgs[j]);
            }

            pair<EdgeSet::iterator, bool> test = m_localEdges.insert(edgs[j]);

            if (test.second)
            {
                (*test.first)
                    ->m_elLink.push_back(pair<ElementSharedPtr, int>(E, j));
            }
            else
            {
                E->SetEdge(j, *test.first);
                (*test.first)
                    ->m_elLink.push_back(pair<ElementSharedPtr, int>(E, j));
            }
        }
        E->SetId(m_localElements.size());
        m_localElements.push_back(E);
    }
}

void FaceMesh::Stretching()
{
    // define a sampling and calculate the aspect ratio of the paramter plane
    m_str                       = 0.0;
    Array<OneD, NekDouble> bnds = m_cadsurf->GetBounds();

    NekDouble dxu = int(bnds[1] - bnds[0] < bnds[3] - bnds[2]
                            ? 40
                            : (bnds[1] - bnds[0]) / (bnds[3] - bnds[2]) * 40);
    NekDouble dxv = int(bnds[3] - bnds[2] < bnds[1] - bnds[0]
                            ? 40
                            : (bnds[3] - bnds[2]) / (bnds[1] - bnds[0]) * 40);

    NekDouble du = (bnds[1] - bnds[0]) / dxu;
    NekDouble dv = (bnds[3] - bnds[2]) / dxv;

    int ct = 0;

    for (int i = 0; i < dxu; i++)
    {
        for (int j = 0; j < dxv; j++)
        {
            Array<OneD, NekDouble> uv(2);
            uv[0] = bnds[0] + i * du;
            uv[1] = bnds[2] + j * dv;
            if (i == dxu - 1)
                uv[0] = bnds[1];
            if (j == dxv - 1)
                uv[1]                = bnds[3];
            Array<OneD, NekDouble> r = m_cadsurf->D1(uv);

            NekDouble ru = sqrt(r[3] * r[3] + r[4] * r[4] + r[5] * r[5]);
            NekDouble rv = sqrt(r[6] * r[6] + r[7] * r[7] + r[8] * r[8]);

            ru *= du;
            rv *= dv;

            if (rv < 1E-8)
                continue;

            m_str += ru / rv;
            ct++;
        }
    }

    m_str /= ct;
}

bool FaceMesh::Validate()
{
    // check all edges in the current mesh for length against the octree
    // if the octree is not conformed to add a new point inside the triangle
    // if no new points are added meshing can stop
    int pointBefore = m_stienerpoints.size();
    for (int i = 0; i < m_connec.size(); i++)
    {
        Array<OneD, NekDouble> triDelta(3);

        Array<OneD, NekDouble> r(3);

        r[0] = m_connec[i][0]->Distance(m_connec[i][1]);
        r[1] = m_connec[i][1]->Distance(m_connec[i][2]);
        r[2] = m_connec[i][2]->Distance(m_connec[i][0]);

        triDelta[0] = m_octree->Query(m_connec[i][0]->GetLoc());
        triDelta[1] = m_octree->Query(m_connec[i][1]->GetLoc());
        triDelta[2] = m_octree->Query(m_connec[i][2]->GetLoc());

        int numValid = 0;

        if (r[0] < triDelta[0] && r[2] < triDelta[0])
            numValid++;

        if (r[1] < triDelta[1] && r[0] < triDelta[1])
            numValid++;

        if (r[2] < triDelta[2] && r[1] < triDelta[2])
            numValid++;

        if (numValid != 3)
        {
            Array<OneD, NekDouble> ainfo, binfo, cinfo;
            ainfo = m_connec[i][0]->GetCADSurfInfo(m_id);
            binfo = m_connec[i][1]->GetCADSurfInfo(m_id);
            cinfo = m_connec[i][2]->GetCADSurfInfo(m_id);

            Array<OneD, NekDouble> uvc(2);
            uvc[0] = (ainfo[0] + binfo[0] + cinfo[0]) / 3.0;
            uvc[1] = (ainfo[1] + binfo[1] + cinfo[1]) / 3.0;
            AddNewPoint(uvc);
        }
    }

    if (m_stienerpoints.size() == pointBefore)
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
    // adds a new point but checks that there are no other points nearby first
    Array<OneD, NekDouble> np = m_cadsurf->P(uv);
    NekDouble npDelta         = m_octree->Query(np);

    NodeSharedPtr n = boost::shared_ptr<Node>(
        new Node(m_mesh->m_numNodes++, np[0], np[1], np[2]));

    bool add = true;

    for (int i = 0; i < orderedLoops.size(); i++)
    {
        for (int j = 0; j < orderedLoops[i].size(); j++)
        {
            NekDouble r = orderedLoops[i][j]->Distance(n);

            if (r < npDelta / 2.0)
            {
                add = false;
                break;
            }
        }
    }

    if (add)
    {
        for (int i = 0; i < m_stienerpoints.size(); i++)
        {
            NekDouble r = m_stienerpoints[i]->Distance(n);

            if (r < npDelta / 2.0)
            {
                add = false;
                break;
            }
        }
    }

    if (add)
    {
        n->SetCADSurf(m_id, m_cadsurf, uv);
        m_stienerpoints.push_back(n);
    }
}

void FaceMesh::OrientateCurves()
{
    // create list of bounding loop nodes
    for (int i = 0; i < m_edgeloops.size(); i++)
    {
        vector<NodeSharedPtr> cE;
        for (int j = 0; j < m_edgeloops[i].edges.size(); j++)
        {
            int cid = m_edgeloops[i].edges[j]->GetId();
            vector<NodeSharedPtr> edgePoints =
                m_curvemeshes[cid]->GetMeshPoints();

            int numPoints = m_curvemeshes[cid]->GetNumPoints();

            if (m_edgeloops[i].edgeo[j] == 0)
            {
                for (int k = 0; k < numPoints - 1; k++)
                {
                    cE.push_back(edgePoints[k]);
                }
            }
            else
            {
                for (int k = numPoints - 1; k > 0; k--)
                {
                    cE.push_back(edgePoints[k]);
                }
            }
        }
        orderedLoops.push_back(cE);
    }

    // loops made need to orientate on which is biggest and define holes
    for (int i = 0; i < orderedLoops.size(); i++)
    {
        NekDouble area = 0.0;
        for (int j = 0; j < orderedLoops[i].size() - 1; j++)
        {
            Array<OneD, NekDouble> n1info, n2info;
            n1info = orderedLoops[i][j]->GetCADSurfInfo(m_id);
            n2info = orderedLoops[i][j + 1]->GetCADSurfInfo(m_id);

            area += -n2info[1] * (n2info[0] - n1info[0]) +
                    n1info[0] * (n2info[1] - n1info[1]);
        }
        area *= 0.5;
        m_edgeloops[i].area = area;
    }

    int ct = 0;

    do
    {
        ct = 0;
        for (int i = 0; i < m_edgeloops.size() - 1; i++)
        {
            if (fabs(m_edgeloops[i].area) < fabs(m_edgeloops[i + 1].area))
            {
                // swap
                vector<NodeSharedPtr> orderedlooptmp = orderedLoops[i];
                EdgeLoop edgeLoopstmp                = m_edgeloops[i];

                orderedLoops[i] = orderedLoops[i + 1];
                m_edgeloops[i]  = m_edgeloops[i + 1];

                orderedLoops[i + 1] = orderedlooptmp;
                m_edgeloops[i + 1]  = edgeLoopstmp;

                ct += 1;
            }
        }

    } while (ct > 0);

    for (int i = 0; i < orderedLoops.size(); i++)
    {
        NodeSharedPtr n1, n2;

        n1 = orderedLoops[i][0];
        n2 = orderedLoops[i][1];

        Array<OneD, NekDouble> n1info, n2info;
        n1info = n1->GetCADSurfInfo(m_id);
        n2info = n2->GetCADSurfInfo(m_id);

        Array<OneD, NekDouble> N(2);
        NekDouble mag = sqrt((n1info[0] - n2info[0]) * (n1info[0] - n2info[0]) +
                             (n1info[1] - n2info[1]) * (n1info[1] - n2info[1]));
        ASSERTL0(mag > 1e-30, "infinity");
        N[0] = -1.0 * (n2info[1] - n1info[1]) / mag;
        N[1] = (n2info[0] - n1info[0]) / mag;

        Array<OneD, NekDouble> P(2);
        P[0] = (n1info[0] + n2info[0]) / 2.0 + 1e-6 * N[0];
        P[1] = (n1info[1] + n2info[1]) / 2.0 + 1e-6 * N[1];

        // now test to see if p is inside or outside the shape
        // vector to the right
        int intercepts = 0;
        for (int j = 0; j < orderedLoops[i].size() - 1; j++)
        {
            Array<OneD, NekDouble> nt1, nt2;
            nt1 = orderedLoops[i][j]->GetCADSurfInfo(m_id);
            nt2 = orderedLoops[i][j + 1]->GetCADSurfInfo(m_id);

            if (fabs(nt2[1] - nt1[1]) < 1e-30)
                continue;

            NekDouble lam = (P[1] - nt1[1]) / (nt2[1] - nt1[1]);
            NekDouble S   = nt1[0] - P[0] + (nt2[0] - nt1[0]) * lam;

            if (!(lam < 0) && !(lam > 1) && S > 0)
            {
                intercepts++;
            }
        }
        {
            Array<OneD, NekDouble> nt1, nt2;
            nt1 = orderedLoops[i].back()->GetCADSurfInfo(m_id);
            nt2 = orderedLoops[i][0]->GetCADSurfInfo(m_id);

            if (fabs(nt2[1] - nt1[1]) < 1e-30)
                continue;

            NekDouble lam = (P[1] - nt1[1]) / (nt2[1] - nt1[1]);
            NekDouble S   = nt1[0] - P[0] + (nt2[0] - nt1[0]) * lam;

            if (!(lam < 0) && !(lam > 1) && S > 0)
            {
                intercepts++;
            }
        }
        if (intercepts % 2 == 0)
        {
            P[0]       = (n1info[0] + n2info[0]) / 2.0 - 1e-6 * N[0];
            P[1]       = (n1info[1] + n2info[1]) / 2.0 - 1e-6 * N[1];
            intercepts = 0;
            for (int j = 0; j < orderedLoops[i].size() - 1; j++)
            {
                Array<OneD, NekDouble> nt1, nt2;
                nt1 = orderedLoops[i][j]->GetCADSurfInfo(m_id);
                nt2 = orderedLoops[i][j + 1]->GetCADSurfInfo(m_id);

                if (fabs(nt2[1] - nt1[1]) < 1e-30)
                    continue;

                NekDouble lam = (P[1] - nt1[1]) / (nt2[1] - nt1[1]);
                NekDouble S   = nt1[0] - P[0] + (nt2[0] - nt1[0]) * lam;

                if (!(lam < 0) && !(lam > 1) && S > 0)
                {
                    intercepts++;
                }
            }
            {
                Array<OneD, NekDouble> nt1, nt2;
                nt1 = orderedLoops[i].back()->GetCADSurfInfo(m_id);
                nt2 = orderedLoops[i][0]->GetCADSurfInfo(m_id);

                if (fabs(nt2[1] - nt1[1]) < 1e-30)
                    continue;

                NekDouble lam = (P[1] - nt1[1]) / (nt2[1] - nt1[1]);
                NekDouble S   = nt1[0] - P[0] + (nt2[0] - nt1[0]) * lam;

                if (!(lam < 0) && !(lam > 1) && S > 0)
                {
                    intercepts++;
                }
            }
            if (intercepts % 2 == 0)
            {
                cerr << "still failed to find point inside loop" << endl;
            }
        }

        m_edgeloops[i].center = P;
    }

    if (m_edgeloops[0].area < 0) // reverse the first uvLoop
    {
        vector<NodeSharedPtr> tmp = orderedLoops[0];
        reverse(tmp.begin(), tmp.end());
        orderedLoops[0] = tmp;
    }

    for (int i = 1; i < orderedLoops.size(); i++)
    {
        if (m_edgeloops[i].area > 0) // reverse the loop
        {
            vector<NodeSharedPtr> tmp = orderedLoops[i];
            reverse(tmp.begin(), tmp.end());
            orderedLoops[i] = tmp;
        }
    }
}
}
}
