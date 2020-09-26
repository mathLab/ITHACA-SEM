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

#include <NekMeshUtils/ExtLibInterface/TriangleInterface.h>
#include <NekMeshUtils/Octree/Octree.h>
#include <NekMeshUtils/SurfaceMeshing/FaceMesh.h>
#include <limits>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

bool FaceMesh::ValidateCurves()
{
    vector<int> curvesInSurface;
    for (int i = 0; i < m_edgeloops.size(); i++)
    {
        for (int j = 0; j < m_edgeloops[i]->edges.size(); j++)
        {
            curvesInSurface.push_back(m_edgeloops[i]->edges[j]->GetId());
        }
    }

    bool error = false;

    for (int i = 0; i < curvesInSurface.size(); i++)
    {
        vector<EdgeSharedPtr> es =
            m_curvemeshes[curvesInSurface[i]]->GetMeshEdges();

        for (int j = i; j < curvesInSurface.size(); j++)
        {
            if (i == j)
            {
                continue;
            }

            vector<EdgeSharedPtr> es2 =
                m_curvemeshes[curvesInSurface[j]]->GetMeshEdges();

            for (int l = 0; l < es.size(); l++)
            {
                Array<OneD, NekDouble> P1 = es[l]->m_n1->GetCADSurfInfo(m_id);
                Array<OneD, NekDouble> P2 = es[l]->m_n2->GetCADSurfInfo(m_id);
                for (int k = 0; k < es2.size(); k++)
                {
                    if (es[l]->m_n1 == es2[k]->m_n1 ||
                        es[l]->m_n1 == es2[k]->m_n2 ||
                        es[l]->m_n2 == es2[k]->m_n1 ||
                        es[l]->m_n2 == es2[k]->m_n2)
                    {
                        continue;
                    }

                    Array<OneD, NekDouble> P3 =
                        es2[k]->m_n1->GetCADSurfInfo(m_id);
                    Array<OneD, NekDouble> P4 =
                        es2[k]->m_n2->GetCADSurfInfo(m_id);

                    NekDouble den = (P4[0] - P3[0]) * (P2[1] - P1[1]) -
                                    (P2[0] - P1[0]) * (P4[1] - P3[1]);

                    if (fabs(den) < 1e-8)
                    {
                        continue;
                    }

                    NekDouble t = ((P1[0] - P3[0]) * (P4[1] - P3[1]) -
                                   (P4[0] - P3[0]) * (P1[1] - P3[1])) /
                                  den;
                    NekDouble u =
                        (P1[0] - P3[0] + t * (P2[0] - P1[0])) / (P4[0] - P3[0]);

                    if (t < 1.0 && t > 0.0 && u < 1.0 && u > 0.0)
                    {
                        Array<OneD, NekDouble> uv(2);
                        uv[0] = P1[0] + t * (P2[0] - P1[0]);
                        uv[1] = P1[1] + t * (P2[1] - P1[1]);
                        Array<OneD, NekDouble> loc = m_cadsurf->P(uv);
                        cout << endl
                             << "Curve mesh error at " << loc[0] << " "
                             << loc[1] << " " << loc[2] << " on face " << m_id
                             << endl;
                        error = true;
                    }
                }
            }
        }
    }
    return error;
}

void FaceMesh::ValidateLoops()
{
    OrientateCurves();

    for (int i = 0; i < orderedLoops.size(); i++)
    {
        int numPoints = orderedLoops[i].size();
        if(numPoints == 2)
        {
            //force a remesh of the curves
            for(int j = 0; j < m_edgeloops[i]->edges.size(); j++)
            {
                int cid = m_edgeloops[i]->edges[j]->GetId();
                m_curvemeshes[cid]->ReMesh();
            }
        }
    }
}

void FaceMesh::Mesh()
{
    Stretching();
    OrientateCurves();

    int numPoints = 0;
    for (int i = 0; i < orderedLoops.size(); i++)
    {
        numPoints += orderedLoops[i].size();
        for (int j = 0; j < orderedLoops[i].size(); j++)
        {
            m_inBoundary.insert(orderedLoops[i][j]);
        }
    }

    stringstream ss;
    ss << "3 points required for triangulation, " << numPoints << " in loop"
       << endl;
    ss << "curves: ";
    for (int i = 0; i < m_edgeloops.size(); i++)
    {
        for (int j = 0; j < m_edgeloops[i]->edges.size(); j++)
        {
            ss << m_edgeloops[i]->edges[j]->GetId() << " ";
        }
    }

    ASSERTL0(numPoints > 2, "number of verts in face is less than 3");

    // create interface to triangle thirdparty library
    TriangleInterfaceSharedPtr pplanemesh =
        MemoryManager<TriangleInterface>::AllocateSharedPtr();

    vector<Array<OneD, NekDouble> > centers;
    for (int i = 0; i < m_edgeloops.size(); i++)
    {
        centers.push_back(m_edgeloops[i]->center);
    }

    pplanemesh->Assign(orderedLoops, centers, m_id, m_str);

    pplanemesh->Mesh();

    pplanemesh->Extract(m_connec);

    bool repeat     = true;
    int meshcounter = 1;

    // continuously remesh until all triangles conform to the spacing in the
    // octree
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
        //    break;
    }

    // build a local version of the mesh (one set of triangles).  this is done
    // so edge connectivity infomration can be used for optimisation
    BuildLocalMesh();

    OptimiseLocalMesh();

    // make new elements and add to list from list of nodes and connectivity
    // from triangle removing unnesercary infomration from the elements
    for (int i = 0; i < m_localElements.size(); i++)
    {
        vector<EdgeSharedPtr> e = m_localElements[i]->GetEdgeList();
        for (int j = 0; j < e.size(); j++)
        {
            e[j]->m_elLink.clear();
        }
        m_mesh->m_element[2].push_back(m_localElements[i]);
    }

    if (m_mesh->m_verbose)
    {
        cout << "\r                               "
                "                                 "
                "                             ";
        cout << scientific << "\r\t\tFace " << m_id << endl
             << "\t\t\tNodes: " << m_localNodes.size() << endl
             << "\t\t\tEdges: " << m_localEdges.size() << endl
             << "\t\t\tTriangles: " << m_localElements.size() << endl
             << "\t\t\tLoops: " << m_edgeloops.size() << endl
             << "\t\t\tSTR: " << m_str << endl
             << endl;
    }
}

void FaceMesh::OptimiseLocalMesh()
{
    // each optimisation algorithm is based on the work in chapter 19
    DiagonalSwap();

    Smoothing();

    DiagonalSwap();

    Smoothing();
}

void FaceMesh::Smoothing()
{
    EdgeSet::iterator eit;
    NodeSet::iterator nit;

    Array<OneD, NekDouble> bounds = m_cadsurf->GetBounds();

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
            NodeSet::iterator f = m_inBoundary.find((*nit));

            // node is on curve so skip
            if (f != m_inBoundary.end())
            {
                continue;
            }

            // this can be real nodes or dummy nodes depending on the system
            vector<NodeSharedPtr> connodes;

            vector<EdgeSharedPtr> edges  = connectingedges[(*nit)->m_id];
            vector<ElementSharedPtr> els = connectingelements[(*nit)->m_id];

            vector<NodeSharedPtr> nodesystem;
            vector<NekDouble> lamp;

            for (int i = 0; i < edges.size(); i++)
            {
                vector<NekDouble> lambda;

                NodeSharedPtr J;
                if (*nit == edges[i]->m_n1)
                {
                    J = edges[i]->m_n2;
                }
                else if (*nit == edges[i]->m_n2)
                {
                    J = edges[i]->m_n1;
                }
                else
                {
                    ASSERTL0(false, "could not find node");
                }

                Array<OneD, NekDouble> ui = (*nit)->GetCADSurfInfo(m_id);
                Array<OneD, NekDouble> uj = J->GetCADSurfInfo(m_id);

                for (int j = 0; j < els.size(); j++)
                {
                    vector<NodeSharedPtr> v = els[j]->GetVertexList();

                    // elememt is adjacent to J therefore no intersection on IJ
                    if (v[0] == J || v[1] == J || v[2] == J)
                    {
                        continue;
                    }

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
                    ud[0] = uj[0] + lambda[0] * (ui[0] - uj[0]);
                    ud[1] = uj[1] + lambda[0] * (ui[1] - uj[1]);
                    Array<OneD, NekDouble> locd = m_cadsurf->P(ud);
                    NodeSharedPtr dn = std::shared_ptr<Node>(
                        new Node(0, locd[0], locd[1], locd[2]));
                    dn->SetCADSurf(m_cadsurf, ud);

                    nodesystem.push_back(dn);
                    lamp.push_back(lambda[0]);
                }
                else
                {
                    nodesystem.push_back(J);
                    lamp.push_back(1.0);
                }
            }

            Array<OneD, NekDouble> u0(2);
            u0[0] = 0.0;
            u0[1] = 0.0;

            for (int i = 0; i < nodesystem.size(); i++)
            {
                Array<OneD, NekDouble> uj = nodesystem[i]->GetCADSurfInfo(m_id);
                u0[0] += uj[0] / nodesystem.size();
                u0[1] += uj[1] / nodesystem.size();
            }

            /*Array<OneD, NekDouble> pu0 = m_cadsurf->P(u0);
            NekDouble di = m_mesh->m_octree->Query(pu0);
            Array<OneD, NekDouble> F(2, 0.0), dF(4, 0.0);
            for (int i = 0; i < nodesystem.size(); i++)
            {
                Array<OneD, NekDouble> rj = nodesystem[i]->GetLoc();
                Array<OneD, NekDouble> uj = nodesystem[i]->GetCADSurfInfo(m_id);
                NekDouble dj  = m_mesh->m_octree->Query(rj);
                NekDouble d   = (di + dj) / 2.0;
                NekDouble wij = sqrt((rj[0] - pu0[0]) * (rj[0] - pu0[0]) +
                                     (rj[1] - pu0[1]) * (rj[1] - pu0[1]) +
                                     (rj[2] - pu0[2]) * (rj[2] - pu0[2])) -
                                d;

                NekDouble umag = sqrt((uj[0] - u0[0]) * (uj[0] - u0[0]) +
                                      (uj[1] - u0[1]) * (uj[1] - u0[1]));

                F[0] += wij * (uj[0] - u0[0]) / umag;
                F[1] += wij * (uj[1] - u0[1]) / umag;

                Array<OneD, NekDouble> d1 = m_cadsurf->D1(u0);

                Array<OneD, NekDouble> dw(2, 0.0);
                dw[0] = -2.0 *
                        ((rj[0] - pu0[0]) * d1[3] + (rj[1] - pu0[1]) * d1[4] +
                         (rj[2] - pu0[2]) * d1[5]);
                dw[1] = -2.0 *
                        ((rj[0] - pu0[0]) * d1[6] + (rj[1] - pu0[1]) * d1[7] +
                         (rj[2] - pu0[2]) * d1[8]);

                Array<OneD, NekDouble> drhs(4, 0.0);
                drhs[0] = 2.0 * ((uj[0] - u0[0]) * (uj[0] - u0[0]) - umag) /
                          umag / umag;
                drhs[1] =
                    2.0 * ((uj[0] - u0[0]) * (uj[1] - u0[1])) / umag / umag;
                drhs[2] = drhs[1];
                drhs[3] = 2.0 * ((uj[1] - u0[1]) * (uj[1] - u0[1]) - umag) /
                          umag / umag;

                dF[0] += dw[0] * (uj[0] - u0[0]) / umag + wij * drhs[0];

                dF[1] += dw[0] * (uj[1] - u0[1]) / umag + wij * drhs[1];

                dF[2] += dw[1] * (uj[0] - u0[0]) / umag + wij * drhs[2];

                dF[3] += dw[1] * (uj[1] - u0[1]) / umag + wij * drhs[3];
            }

            NekDouble det = dF[0] * dF[3] - dF[1] * dF[2];
            NekDouble tmp = dF[0];
            dF[0]         = dF[3] / det;
            dF[3]         = tmp / det;
            dF[1] *= -1.0 / det;
            dF[2] *= -1.0 / det;

            u0[0] -= (dF[0] * F[0] + dF[2] * F[1]);
            u0[1] -= (dF[1] * F[0] + dF[3] * F[1]);*/

            bool inbounds = true;
            if (u0[0] < bounds[0])
            {
                inbounds = false;
            }
            else if (u0[0] > bounds[1])
            {
                inbounds = false;
            }
            else if (u0[1] < bounds[2])
            {
                inbounds = false;
            }
            else if (u0[1] > bounds[3])
            {
                inbounds = false;
            }

            if (!inbounds)
            {
                continue;
            }

            /*Array<OneD, NekDouble> FN(2, 0.0);
            pu0 = m_cadsurf->P(u0);
            di  = m_mesh->m_octree->Query(pu0);
            for (int i = 0; i < nodesystem.size(); i++)
            {
                Array<OneD, NekDouble> rj = nodesystem[i]->GetLoc();
                Array<OneD, NekDouble> uj = nodesystem[i]->GetCADSurfInfo(m_id);
                NekDouble dj  = m_mesh->m_octree->Query(rj);
                NekDouble d   = (di + dj) / 2.0;
                NekDouble wij = sqrt((rj[0] - pu0[0]) * (rj[0] - pu0[0]) +
                                     (rj[1] - pu0[1]) * (rj[1] - pu0[1]) +
                                     (rj[2] - pu0[2]) * (rj[2] - pu0[2])) - d;

                NekDouble umag = sqrt((uj[0] - u0[0]) * (uj[0] - u0[0]) +
                                      (uj[1] - u0[1]) * (uj[1] - u0[1]));

                FN[0] += wij * (uj[0] - u0[0]) / umag;
                FN[1] += wij * (uj[1] - u0[1]) / umag;
            }

            if (F[0] * F[0] + F[1] * F[1] < FN[0] * FN[0] + FN[1] * FN[1])
            {
                continue;
            }*/

            Array<OneD, NekDouble> l2 = m_cadsurf->P(u0);
            (*nit)->Move(l2, m_id, u0);
        }
    }
}

void FaceMesh::DiagonalSwap()
{
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
        // this routine is broken and needs looking at
        NodeSet::iterator f = m_inBoundary.find((*nit));
        if (f == m_inBoundary.end()) // node isnt on curve so skip
        {
            // node is interior
            idealConnec[(*nit)->m_id] = 6;
        }
        else
        {
            // need to identify the two other nodes on the boundary to find
            // interior angle
            /*vector<NodeSharedPtr> ns;
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
                ceil((*nit)->Angle(ns[0], ns[1]) / 3.142 * 3) + 1;*/
            idealConnec[(*nit)->m_id] = 4;
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
        EdgeSet edges  = m_localEdges;
        m_localEdges.clear();

        int swappedEdges = 0;

        EdgeSet::iterator it;

        for (it = edges.begin(); it != edges.end(); it++)
        {
            EdgeSharedPtr e = *it;

            NodeSet::iterator f1 = m_inBoundary.find((*it)->m_n1);
            NodeSet::iterator f2 = m_inBoundary.find((*it)->m_n2);
            if (f1 != m_inBoundary.end() && f2 != m_inBoundary.end())
            {
                m_localEdges.insert(e);
                continue;
            }

            ElementSharedPtr tri1 = e->m_elLink[0].first.lock();
            ElementSharedPtr tri2 = e->m_elLink[1].first.lock();

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
            // //cout<< A->GetNumCADSurf() << " " << B->GetNumCADSurf() << " "
            // //    << C->GetNumCADSurf() << " " << D->GetNumCADSurf() << endl;
            // ofstream file;
            // file.open("pts.3D");
            // file << "x y z value" << endl;
            // file << A->m_x << " " << A->m_y << " " << A->m_z << endl;
            // file << B->m_x << " " << B->m_y << " " << B->m_z << endl;
            // file << C->m_x << " " << C->m_y << " " << C->m_z << endl;
            // file << D->m_x << " " << D->m_y << " " << D->m_z << endl;
            // file.close();
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
                CAt = std::shared_ptr<Edge>(new Edge(C, A));
                ADt = std::shared_ptr<Edge>(new Edge(A, D));
                DBt = std::shared_ptr<Edge>(new Edge(D, B));
                BCt = std::shared_ptr<Edge>(new Edge(B, C));

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
                vector<pair<weak_ptr<Element>, int> > links;

                links = CA->m_elLink;
                CA->m_elLink.clear();
                for (int i = 0; i < links.size(); i++)
                {
                    if (links[i].first.lock()->GetId() == tri1->GetId())
                    {
                        continue;
                    }
                    CA->m_elLink.push_back(links[i]);
                }

                links = BC->m_elLink;
                BC->m_elLink.clear();
                for (int i = 0; i < links.size(); i++)
                {
                    if (links[i].first.lock()->GetId() == tri1->GetId())
                    {
                        continue;
                    }
                    BC->m_elLink.push_back(links[i]);
                }

                links = AD->m_elLink;
                AD->m_elLink.clear();
                for (int i = 0; i < links.size(); i++)
                {
                    if (links[i].first.lock()->GetId() == tri2->GetId())
                    {
                        continue;
                    }
                    AD->m_elLink.push_back(links[i]);
                }

                links = DB->m_elLink;
                DB->m_elLink.clear();
                for (int i = 0; i < links.size(); i++)
                {
                    if (links[i].first.lock()->GetId() == tri2->GetId())
                    {
                        continue;
                    }
                    DB->m_elLink.push_back(links[i]);
                }

                EdgeSharedPtr newe = std::shared_ptr<Edge>(new Edge(C, D));

                vector<NodeSharedPtr> t1, t2;
                t1.push_back(B);
                t1.push_back(D);
                t1.push_back(C);
                t2.push_back(A);
                t2.push_back(C);
                t2.push_back(D);

                ElmtConfig conf(LibUtilities::eTriangle, 1, false, false,
                                m_mesh->m_spaceDim != 3);
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
                ntri1->m_parentCAD = m_cadsurf;
                ntri2->m_parentCAD = m_cadsurf;

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

    for (int i = 0; i < m_connec.size(); i++)
    {
        ElmtConfig conf(LibUtilities::eTriangle, 1, false, false,
                        m_mesh->m_spaceDim != 3);

        vector<int> tags;
        tags.push_back(m_compId);
        ElementSharedPtr E = GetElementFactory().CreateInstance(
            LibUtilities::eTriangle, conf, m_connec[i], tags);
        E->m_parentCAD = m_cadsurf;

        vector<NodeSharedPtr> nods = E->GetVertexList();
        for (int j = 0; j < nods.size(); j++)
        {
            // nodes are already unique some will insert some wont
            m_localNodes.insert(nods[j]);
        }

        E->SetId(m_localElements.size());
        m_localElements.push_back(E);
    }

    for (int i = 0; i < m_localElements.size(); ++i)
    {
        for (int j = 0; j < m_localElements[i]->GetEdgeCount(); ++j)
        {
            pair<EdgeSet::iterator, bool> testIns;
            EdgeSharedPtr ed = m_localElements[i]->GetEdge(j);
            // look for edge in m_mesh edgeset from curves
            EdgeSet::iterator s = m_mesh->m_edgeSet.find(ed);
            if (!(s == m_mesh->m_edgeSet.end()))
            {
                ed = *s;
                m_localElements[i]->SetEdge(j, *s);
            }

            testIns = m_localEdges.insert(ed);

            if (testIns.second)
            {
                EdgeSharedPtr ed2 = *testIns.first;
                ed2->m_elLink.push_back(
                    pair<ElementSharedPtr, int>(m_localElements[i], j));
            }
            else
            {
                EdgeSharedPtr e2 = *(testIns.first);
                m_localElements[i]->SetEdge(j, e2);
                e2->m_elLink.push_back(
                    pair<ElementSharedPtr, int>(m_localElements[i], j));
            }
        }
    }
}

void FaceMesh::Stretching()
{
    // define a sampling and calculate the aspect ratio of the paramter plane
    m_str = 0.0;
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
            {
                uv[0] = bnds[1];
            }
            if (j == dxv - 1)
            {
                uv[1] = bnds[3];
            }
            Array<OneD, NekDouble> r = m_cadsurf->D1(uv);

            NekDouble ru = sqrt(r[3] * r[3] + r[4] * r[4] + r[5] * r[5]);
            NekDouble rv = sqrt(r[6] * r[6] + r[7] * r[7] + r[8] * r[8]);

            ru *= du;
            rv *= dv;

            if (rv < 1E-8)
            {
                continue;
            }

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
        Array<OneD, NekDouble> r(3), a(3);

        vector<Array<OneD, NekDouble> > info;

        for (int j = 0; j < 3; j++)
        {
            info.push_back(m_connec[i][j]->GetCADSurfInfo(m_id));
        }

        r[0] = m_connec[i][0]->Distance(m_connec[i][1]);
        r[1] = m_connec[i][1]->Distance(m_connec[i][2]);
        r[2] = m_connec[i][2]->Distance(m_connec[i][0]);

        a[0] = m_connec[i][0]->Angle(m_connec[i][1]->GetLoc(),
                                     m_connec[i][2]->GetLoc(),
                                     m_cadsurf->N(info[0]));
        a[1] = m_connec[i][1]->Angle(m_connec[i][2]->GetLoc(),
                                     m_connec[i][0]->GetLoc(),
                                     m_cadsurf->N(info[1]));
        a[2] = m_connec[i][2]->Angle(m_connec[i][0]->GetLoc(),
                                     m_connec[i][1]->GetLoc(),
                                     m_cadsurf->N(info[2]));

        NekDouble d1 = m_mesh->m_octree->Query(m_connec[i][0]->GetLoc());
        NekDouble d2 = m_mesh->m_octree->Query(m_connec[i][1]->GetLoc());
        NekDouble d3 = m_mesh->m_octree->Query(m_connec[i][2]->GetLoc());

        Array<OneD, NekDouble> uvc(2);
        uvc[0] = (info[0][0] + info[1][0] + info[2][0]) / 3.0;
        uvc[1] = (info[0][1] + info[1][1] + info[2][1]) / 3.0;

        Array<OneD, NekDouble> locc = m_cadsurf->P(uvc);
        NekDouble d4 = m_mesh->m_octree->Query(locc);

        NekDouble d = (d1 + d2 + d3 + d4) / 4.0;

        vector<bool> valid(3);
        valid[0] = r[0] < d * 1.5;
        valid[1] = r[1] < d * 1.5;
        valid[2] = r[2] < d * 1.5;

        vector<bool> angValid(3);
        angValid[0] = a[0] / M_PI * 180.0 > 20.0 && a[0] / M_PI * 180.0 < 120.0;
        angValid[1] = a[1] / M_PI * 180.0 > 20.0 && a[1] / M_PI * 180.0 < 120.0;
        angValid[2] = a[2] / M_PI * 180.0 > 20.0 && a[2] / M_PI * 180.0 < 120.0;

        int numValid    = 0;
        int numAngValid = 0;
        for (int j = 0; j < 3; j++)
        {
            if (valid[j])
            {
                numValid++;
            }
            if (angValid[j])
            {
                numAngValid++;
            }
        }

        // if numvalid is zero no work to be done
        /*if (numValid != 3)
        {
            AddNewPoint(uvc);
        }*/

        if (numValid != 3 || numAngValid != 3)
        {
            // break the bad edge
            /*int a=0, b=0;
            if(!valid[0])
            {
                a = 0;
                b = 1;
            }
            else if(!valid[1])
            {
                a = 1;
                b = 2;
            }
            else if(!valid[2])
            {
                a = 2;
                b = 0;
            }

            Array<OneD, NekDouble> uvn(2);
            uvn[0] = (info[a][0] + info[b][0]) / 2.0;
            uvn[1] = (info[a][1] + info[b][1]) / 2.0;*/
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
    NekDouble npDelta = m_mesh->m_octree->Query(np);

    NodeSharedPtr n = std::shared_ptr<Node>(
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
        n->SetCADSurf(m_cadsurf, uv);
        m_stienerpoints.push_back(n);
    }
}

void FaceMesh::OrientateCurves()
{
    // this could be a second run on orentate so clear some info
    orderedLoops.clear();

    // create list of bounding loop nodes
    for (int i = 0; i < m_edgeloops.size(); i++)
    {
        vector<NodeSharedPtr> cE;
        for (int j = 0; j < m_edgeloops[i]->edges.size(); j++)
        {
            int cid = m_edgeloops[i]->edges[j]->GetId();
            vector<NodeSharedPtr> edgePoints =
                m_curvemeshes[cid]->GetMeshPoints();

            int numPoints = m_curvemeshes[cid]->GetNumPoints();

            if (m_edgeloops[i]->edgeo[j] == CADOrientation::eForwards)
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
}
}
}
