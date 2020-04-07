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

#include <algorithm>
#include <list>

#include <NekMeshUtils/CADSystem/CADCurve.h>
#include <NekMeshUtils/CADSystem/CADSurf.h>
#include <NekMeshUtils/Optimisation/BGFS-B.h>
#include <NekMeshUtils/SurfaceMeshing/HOSurfaceMesh.h>
#include <NekMeshUtils/SurfaceMeshing/OptimiseFunctions.h>

#include <LibUtilities/BasicUtils/Progressbar.hpp>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LocalRegions/MatrixKey.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

ModuleKey HOSurfaceMesh::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "hosurface"), HOSurfaceMesh::create,
    "Generates a high-order surface mesh based on CAD");

HOSurfaceMesh::HOSurfaceMesh(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["opti"] =
        ConfigOption(true, "0", "Perform edge node optimisation.");
}

HOSurfaceMesh::~HOSurfaceMesh()
{
}

void HOSurfaceMesh::Process()
{
    if (m_mesh->m_verbose)
        cout << endl << "\tHigh-Order Surface meshing" << endl;

    LibUtilities::PointsKey ekey(m_mesh->m_nummode,
                                 LibUtilities::eGaussLobattoLegendre);
    Array<OneD, NekDouble> gll;

    LibUtilities::PointsManager()[ekey]->GetPoints(gll);

    LibUtilities::PointsKey pkey(m_mesh->m_nummode,
                                 LibUtilities::eNodalTriElec);

    Array<OneD, NekDouble> u, v;

    int nq = m_mesh->m_nummode;

    int np = nq * (nq + 1) / 2;

    int ni = (nq - 2) * (nq - 3) / 2;

    int npq = nq * nq;

    int niq = npq - 4 - 4 * (nq - 2);

    LibUtilities::PointsManager()[pkey]->GetPoints(u, v);

    bool qOpti = m_config["opti"].beenSet;

    // loop over all the faces in the surface mesh, check all three edges for
    // high order info, if nothing high-order the edge.

    EdgeSet surfaceEdges;
    EdgeSet completedEdges;

    for (int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        vector<EdgeSharedPtr> es = m_mesh->m_element[2][i]->GetEdgeList();
        for (int j = 0; j < es.size(); j++)
            surfaceEdges.insert(es[j]);
    }

    for (int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        if (m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(i, m_mesh->m_element[2].size(),
                                           "\t\tSurface elements");
        }

        if (!m_mesh->m_element[2][i]->m_parentCAD)
        {
            // no parent cad
            continue;
        }

        CADObjectSharedPtr o = m_mesh->m_element[2][i]->m_parentCAD;
        CADSurfSharedPtr s   = std::dynamic_pointer_cast<CADSurf>(o);
        int surf             = s->GetId();

        FaceSharedPtr f = m_mesh->m_element[2][i]->GetFaceLink();

        bool dumFace = false;

        if (!f)
        {
            //This uses a fake face to build the high-order info
            //in the case of 2D and manifold geometries without having to
            //rewrite the 3D code
            //important to note that face nodes need to be inserted into the
            //volume nodes of the surface element or they will be forgotton
            f = std::shared_ptr<Face>(new Face(
                m_mesh->m_element[2][i]->GetVertexList(),
                vector<NodeSharedPtr>(), m_mesh->m_element[2][i]->GetEdgeList(),
                LibUtilities::ePolyEvenlySpaced));

            dumFace = true;
        }

        f->m_parentCAD = s;

        vector<EdgeSharedPtr> edges = f->m_edgeList;
        for (int j = 0; j < edges.size(); j++)
        {
            EdgeSharedPtr e = edges[j];
            // test insert the edge into completedEdges
            // if the edge already exists move on
            // if not figure out its high-order information

            EdgeSet::iterator test = completedEdges.find(e);

            if (test != completedEdges.end())
            {
                continue;
            }

            // the edges in the element are different to those in the face
            // the cad information is stored in the element edges which are not
            // in the m_mesh->m_edgeSet groups
            // need to link them together and copy the cad information to be
            // able to identify how to make it high-order
            EdgeSet::iterator it = surfaceEdges.find(e);
            ASSERTL0(it != surfaceEdges.end(),
                     "could not find edge in surface");

            if ((*it)->m_parentCAD)
            {
                e->m_parentCAD = (*it)->m_parentCAD;
            }
            else
            {
                e->m_parentCAD = s;
            }

            vector<NodeSharedPtr> honodes(m_mesh->m_nummode - 2);

            if (e->m_parentCAD->GetType() == CADType::eCurve)
            {
                int cid = e->m_parentCAD->GetId();
                CADCurveSharedPtr c =
                    std::dynamic_pointer_cast<CADCurve>(e->m_parentCAD);
                NekDouble tb = e->m_n1->GetCADCurveInfo(cid);
                NekDouble te = e->m_n2->GetCADCurveInfo(cid);

                // distrobute points along curve as inital guess
                Array<OneD, NekDouble> ti(m_mesh->m_nummode);
                for (int k = 0; k < m_mesh->m_nummode; k++)
                {
                    ti[k] =
                        tb * (1.0 - gll[k]) / 2.0 + te * (1.0 + gll[k]) / 2.0;
                }

                if (qOpti)
                {
                    Array<OneD, NekDouble> xi(nq - 2);
                    for (int k = 1; k < nq - 1; k++)
                    {
                        xi[k - 1] = ti[k];
                    }

                    OptiEdgeSharedPtr opti =
                        MemoryManager<OptiEdge>::AllocateSharedPtr(ti, gll, c);

                    DNekMat B(nq - 2, nq - 2,
                              0.0); // approximate hessian (I to start)
                    for (int k = 0; k < nq - 2; k++)
                    {
                        B(k, k) = 1.0;
                    }
                    DNekMat H(nq - 2, nq - 2,
                              0.0); // approximate inverse hessian (I to start)
                    for (int k = 0; k < nq - 2; k++)
                    {
                        H(k, k) = 1.0;
                    }

                    DNekMat J = opti->dF(xi);

                    Array<OneD, NekDouble> bnds = c->GetBounds();

                    bool repeat = true;
                    int itct    = 0;
                    while (repeat)
                    {
                        NekDouble Norm = 0;
                        for (int k = 0; k < nq - 2; k++)
                        {
                            Norm += J(k, 0) * J(k, 0) / (bnds[1] - bnds[0]) /
                                    (bnds[1] - bnds[0]);
                        }
                        Norm = sqrt(Norm);

                        if (Norm < 1E-8)
                        {
                            repeat = false;
                            break;
                        }
                        if (itct > 1000)
                        {
                            cout << "failed to optimise on curve" << endl;
                            for (int k = 0; k < nq; k++)
                            {
                                ti[k] = tb * (1.0 - gll[k]) / 2.0 +
                                        te * (1.0 + gll[k]) / 2.0;
                            }
                            break;
                        }
                        itct++;

                        if (!BGFSUpdate(opti, J, B, H))
                        {
                            if (m_mesh->m_verbose)
                            {
                                cout << "BFGS reported no update, curve on "
                                     << c->GetId() << endl;
                            }
                            break;
                        }
                    }
                    // need to pull the solution out of opti
                    ti = opti->GetSolution();
                }
                vector<pair<weak_ptr<CADSurf>, CADOrientation::Orientation>> s =
                    c->GetAdjSurf();

                for (int k = 1; k < m_mesh->m_nummode - 1; k++)
                {
                    Array<OneD, NekDouble> loc = c->P(ti[k]);
                    NodeSharedPtr nn = std::shared_ptr<Node>(
                        new Node(0, loc[0], loc[1], loc[2]));

                    nn->SetCADCurve(c, ti[k]);
                    for (int m = 0; m < s.size(); m++)
                    {
                        Array<OneD, NekDouble> uv =
                            s[m].first.lock()->locuv(loc);
                        nn->SetCADSurf(s[m].first.lock(), uv);
                    }

                    honodes[k - 1] = nn;
                }
            }
            else
            {
                // edge is on surface and needs 2d optimisation
                Array<OneD, NekDouble> uvb, uve;
                uvb            = e->m_n1->GetCADSurfInfo(surf);
                uve            = e->m_n2->GetCADSurfInfo(surf);

                Array<OneD, NekDouble> l1 = e->m_n1->GetLoc();
                Array<OneD, NekDouble> l2 = e->m_n2->GetLoc();

                e->m_parentCAD = s;
                Array<OneD, Array<OneD, NekDouble>> uvi(nq);

                for (int k = 0; k < nq; k++)
                {
                    Array<OneD, NekDouble> uv(2);
                    uv[0] = uvb[0] * (1.0 - gll[k]) / 2.0 +
                            uve[0] * (1.0 + gll[k]) / 2.0;
                    uv[1] = uvb[1] * (1.0 - gll[k]) / 2.0 +
                            uve[1] * (1.0 + gll[k]) / 2.0;
                    uvi[k] = uv;
                }

                if (qOpti)
                {
                    Array<OneD, NekDouble> bnds = s->GetBounds();
                    Array<OneD, NekDouble> all(2 * nq);
                    for (int k = 0; k < nq; k++)
                    {
                        all[k * 2 + 0] = uvi[k][0];
                        all[k * 2 + 1] = uvi[k][1];
                    }

                    Array<OneD, NekDouble> xi(2 * (nq - 2));
                    for (int k = 1; k < nq - 1; k++)
                    {
                        xi[(k - 1) * 2 + 0] = all[k * 2 + 0];
                        xi[(k - 1) * 2 + 1] = all[k * 2 + 1];
                    }

                    OptiEdgeSharedPtr opti =
                        MemoryManager<OptiEdge>::AllocateSharedPtr(all, gll, s);

                    DNekMat B(2 * (nq - 2), 2 * (nq - 2),
                              0.0); // approximate hessian (I to start)
                    for (int k = 0; k < 2 * (nq - 2); k++)
                    {
                        B(k, k) = 1.0;
                    }
                    DNekMat H(2 * (nq - 2), 2 * (nq - 2),
                              0.0); // approximate inverse hessian (I to start)
                    for (int k = 0; k < 2 * (nq - 2); k++)
                    {
                        H(k, k) = 1.0;
                    }

                    DNekMat J = opti->dF(xi);

                    bool repeat = true;
                    int itct    = 0;
                    while (repeat)
                    {
                        NekDouble Norm = 0;
                        for (int k = 0; k < 2 * (nq - 2); k++)
                        {
                            if (k % 2 == 0)
                            {
                                Norm +=
                                    J(k, 0) * J(k, 0); // (bnds[1] - bnds[0]) /
                                                       //(bnds[1] - bnds[0]);
                            }
                            else
                            {
                                Norm +=
                                    J(k, 0) * J(k, 0); // (bnds[3] - bnds[2]) /
                                                       //(bnds[3] - bnds[2]);
                            }
                        }
                        Norm = sqrt(Norm);

                        if (Norm < 2E-2)
                        {
                            repeat = false;
                            break;
                        }

                        if (itct > 1000)
                        {
                            cout << "failed to optimise on edge " << Norm
                                 << endl;
                            for (int k = 0; k < nq; k++)
                            {
                                Array<OneD, NekDouble> uv(2);
                                uv[0] = uvb[0] * (1.0 - gll[k]) / 2.0 +
                                        uve[0] * (1.0 + gll[k]) / 2.0;
                                uv[1] = uvb[1] * (1.0 - gll[k]) / 2.0 +
                                        uve[1] * (1.0 + gll[k]) / 2.0;
                                uvi[k] = uv;
                            }
                            break;
                        }
                        itct++;

                        if (!BGFSUpdate(opti, J, B, H))
                        {
                            if (m_mesh->m_verbose)
                            {
                                cout << "BFGS reported no update, edge on "
                                     << surf << endl;
                            }
                            // exit(-1);
                            break;
                        }
                    }

                    all = opti->GetSolution();

                    // need to put all backinto uv
                    for (int k = 0; k < nq; k++)
                    {
                        uvi[k][0] = all[k * 2 + 0];
                        uvi[k][1] = all[k * 2 + 1];
                    }
                }

                for (int k = 1; k < nq - 1; k++)
                {
                    Array<OneD, NekDouble> loc = s->P(uvi[k]);
                    NodeSharedPtr nn = std::shared_ptr<Node>(
                        new Node(0, loc[0], loc[1], loc[2]));
                    nn->SetCADSurf(s, uvi[k]);
                    honodes[k - 1] = nn;
                }
            }

            e->m_edgeNodes = honodes;
            e->m_curveType = LibUtilities::eGaussLobattoLegendre;
            completedEdges.insert(e);
        }

        // just add the face interior nodes through interp and project
        vector<NodeSharedPtr> vertices = f->m_vertexList;

        SpatialDomains::GeometrySharedPtr geom = f->GetGeom(3);
        geom->FillGeom();
        StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();
        Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
        Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
        Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);

        Array<OneD, NekDouble> xc(xmap->GetTotPoints());
        Array<OneD, NekDouble> yc(xmap->GetTotPoints());
        Array<OneD, NekDouble> zc(xmap->GetTotPoints());

        xmap->BwdTrans(coeffs0,xc);
        xmap->BwdTrans(coeffs1,yc);
        xmap->BwdTrans(coeffs2,zc);

        if(vertices.size() == 3)
        {
            //build an array of all uvs
            vector<Array<OneD, NekDouble> > uvi;
            for(int j = np-ni; j < np; j++)
            {
                Array<OneD, NekDouble> xp(2);
                xp[0] = u[j];
                xp[1] = v[j];

                Array<OneD, NekDouble> loc(3);
                loc[0] = xmap->PhysEvaluate(xp, xc);
                loc[1] = xmap->PhysEvaluate(xp, yc);
                loc[2] = xmap->PhysEvaluate(xp, zc);

                Array<OneD, NekDouble> uv = s->locuv(loc);
                uvi.push_back(uv);
            }

            vector<NodeSharedPtr> honodes;
            for(int j = 0; j < ni; j++)
            {
                Array<OneD, NekDouble> loc;
                loc = s->P(uvi[j]);
                NodeSharedPtr nn = std::shared_ptr<Node>(new
                                                Node(0,loc[0],loc[1],loc[2]));
                nn->SetCADSurf(s, uvi[j]);
                honodes.push_back(nn);
            }

            f->m_faceNodes = honodes;
            f->m_curveType = LibUtilities::eNodalTriElec;
        }
        else if(vertices.size() == 4)
        {
            //build an array of all uvs
            vector<Array<OneD, NekDouble> > uvi;
            for(int i = 1; i < nq - 1; i++)
            {
                for(int j = 1; j < nq - 1; j++)
                {
                    Array<OneD, NekDouble> xp(2);
                    xp[0] = gll[j];
                    xp[1] = gll[i];

                    Array<OneD, NekDouble> loc(3);
                    loc[0] = xmap->PhysEvaluate(xp, xc);
                    loc[1] = xmap->PhysEvaluate(xp, yc);
                    loc[2] = xmap->PhysEvaluate(xp, zc);

                    Array<OneD, NekDouble> uv = s->locuv(loc);
                    uvi.push_back(uv);
                }
            }

            vector<NodeSharedPtr> honodes;
            for(int j = 0; j < niq; j++)
            {
                Array<OneD, NekDouble> loc;
                loc = s->P(uvi[j]);
                NodeSharedPtr nn = std::shared_ptr<Node>(new
                                                Node(0,loc[0],loc[1],loc[2]));
                nn->SetCADSurf(s, uvi[j]);
                honodes.push_back(nn);
            }

            f->m_faceNodes = honodes;
            f->m_curveType = LibUtilities::eGaussLobattoLegendre;
        }

        if(dumFace)
        {
            m_mesh->m_element[2][i]->SetVolumeNodes(f->m_faceNodes);
            m_mesh->m_element[2][i]->SetCurveType(f->m_curveType);
        }
    }

    if (m_mesh->m_verbose)
        cout << endl;
}
}
}
