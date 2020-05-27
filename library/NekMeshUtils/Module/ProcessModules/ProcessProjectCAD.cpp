////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJac.cpp
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
//  Description: Calculate Jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/MeshElements/Element.h>
#include "ProcessProjectCAD.h"

#include <NekMeshUtils/CADSystem/CADCurve.h>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/BasicUtils/Progressbar.hpp>

#include <boost/core/ignore_unused.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

const NekDouble prismU1[6] = {-1.0, 1.0, 1.0,-1.0,-1.0,-1.0};
const NekDouble prismV1[6] = {-1.0,-1.0, 1.0, 1.0,-1.0, 1.0};
const NekDouble prismW1[6] = {-1.0,-1.0,-1.0,-1.0, 1.0, 1.0};

ModuleKey ProcessProjectCAD::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "projectcad"),
    ProcessProjectCAD::create,
    "Projects mesh to CAD");

ProcessProjectCAD::ProcessProjectCAD(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["file"] =
        ConfigOption(false, "", "CAD file");
    m_config["order"] =
        ConfigOption(false, "4", "CAD file");
}

ProcessProjectCAD::~ProcessProjectCAD()
{
}

bool ProcessProjectCAD::findAndProject(bgi::rtree<boxI, bgi::quadratic<16> > &rtree,
                                       Array<OneD, NekDouble> in,
                                       int &surf)
{
    boost::ignore_unused(surf);

    point q(in[0], in[1], in[2]);
    vector<boxI> result;
    rtree.query(bgi::intersects(q), back_inserter(result));

    if(result.size() == 0)
    {
        //along a projecting edge the node is too far from any surface boxes
        //this is hardly surprising but rare, return false and linearise
        return false;
    }

    int minsurf = 0;
    NekDouble minDist = numeric_limits<double>::max();
    NekDouble dist;

    for (int j = 0; j < result.size(); j++)
    {
        m_mesh->m_cad->GetSurf(result[j].second)->locuv(in, dist);

        if (dist < minDist)
        {
            minDist = dist;
            minsurf = result[j].second;
        }
    }

    Array<OneD, NekDouble> uv = m_mesh->m_cad->GetSurf(minsurf)->locuv(in, dist);

    return true;
}

bool ProcessProjectCAD::IsNotValid(vector<ElementSharedPtr> &els)
{
    //short algebraic method to figure out the vailidy of elements
    for(int i = 0; i < els.size(); i++)
    {
        if(els[i]->GetShapeType() == LibUtilities::ePrism)
        {
            vector<NodeSharedPtr> ns = els[i]->GetVertexList();
            for(int j = 0; j < 6; j++)
            {
                NekDouble a2 = 0.5 * (1 + prismU1[j]);
                NekDouble b1 = 0.5 * (1 - prismV1[j]);
                NekDouble b2 = 0.5 * (1 + prismV1[j]);
                NekDouble c2 = 0.5 * (1 + prismW1[j]);
                NekDouble d  = 0.5 * (prismU1[j] + prismW1[j]);

                Array<OneD, NekDouble> jac(9,0.0);

                jac[0]  = -0.5 * b1 * ns[0]->m_x + 0.5 * b1 * ns[1]->m_x +
                           0.5 * b2 * ns[2]->m_x - 0.5 * b2 * ns[3]->m_x;
                jac[1]  = -0.5 * b1 * ns[0]->m_y + 0.5 * b1 * ns[1]->m_y +
                           0.5 * b2 * ns[2]->m_y - 0.5 * b2 * ns[3]->m_y;
                jac[2]  = -0.5 * b1 * ns[0]->m_z + 0.5 * b1 * ns[1]->m_z +
                           0.5 * b2 * ns[2]->m_z - 0.5 * b2 * ns[3]->m_z;

                jac[3]  = 0.5 * d  * ns[0]->m_x - 0.5 * a2 * ns[1]->m_x +
                          0.5 * a2 * ns[2]->m_x - 0.5 * d  * ns[3]->m_x -
                          0.5 * c2 * ns[4]->m_x + 0.5 * c2 * ns[5]->m_x;
                jac[4]  = 0.5 * d *  ns[0]->m_y - 0.5 * a2 * ns[1]->m_y +
                          0.5 * a2 * ns[2]->m_y - 0.5 * d *  ns[3]->m_y -
                          0.5 * c2 * ns[4]->m_y + 0.5 * c2 * ns[5]->m_y;
                jac[5]  = 0.5 * d *  ns[0]->m_z - 0.5 * a2 * ns[1]->m_z +
                          0.5 * a2 * ns[2]->m_z - 0.5 * d *  ns[3]->m_z -
                          0.5 * c2 * ns[4]->m_z + 0.5 * c2 * ns[5]->m_z;

                jac[6]  = -0.5 * b1 * ns[0]->m_x - 0.5 * b2 * ns[3]->m_x +
                           0.5 * b1 * ns[4]->m_x + 0.5 * b2 * ns[5]->m_x;
                jac[7]  = -0.5 * b1 * ns[0]->m_y - 0.5 * b2 * ns[3]->m_y +
                           0.5 * b1 * ns[4]->m_y + 0.5 * b2 * ns[5]->m_y;
                jac[8]  = -0.5 * b1 * ns[0]->m_z - 0.5 * b2 * ns[3]->m_z +
                           0.5 * b1 * ns[4]->m_z + 0.5 * b2 * ns[5]->m_z;

                NekDouble jc = jac[0] * (jac[4] * jac[8] - jac[5] * jac[7]) -
                               jac[3] * (jac[1] * jac[8] - jac[2] * jac[7]) +
                               jac[6] * (jac[1] * jac[5] - jac[2] * jac[4]);

                if (jc < NekConstants::kNekZeroTol)
                {
                    return true;
                }
            }
        }
        else if(els[i]->GetShapeType() == LibUtilities::eTetrahedron)
        {
            vector<NodeSharedPtr> ns = els[i]->GetVertexList();
            for(int j = 0; j < 4; j++)
            {
                Array<OneD, NekDouble> jac(9,0.0);

                jac[0] = 0.5*(ns[1]->m_x - ns[0]->m_x);
                jac[1] = 0.5*(ns[1]->m_y - ns[0]->m_y);
                jac[2] = 0.5*(ns[1]->m_z - ns[0]->m_z);
                jac[3] = 0.5*(ns[2]->m_x - ns[0]->m_x);
                jac[4] = 0.5*(ns[2]->m_y - ns[0]->m_y);
                jac[5] = 0.5*(ns[2]->m_z - ns[0]->m_z);
                jac[6] = 0.5*(ns[3]->m_x - ns[0]->m_x);
                jac[7] = 0.5*(ns[3]->m_y - ns[0]->m_y);
                jac[8] = 0.5*(ns[3]->m_z - ns[0]->m_z);

                NekDouble jc = jac[0] * (jac[4] * jac[8] - jac[5] * jac[7]) -
                               jac[3] * (jac[1] * jac[8] - jac[2] * jac[7]) +
                               jac[6] * (jac[1] * jac[5] - jac[2] * jac[4]);

                if (jc < NekConstants::kNekZeroTol)
                {
                    return true;
                }
            }
        }
        else
        {
            ASSERTL0(false, "not coded");
        }
    }

    return false;
}

void ProcessProjectCAD::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << "ProcessAssignCAD: Processing... " << endl;
    }

    cout << "ProcessAssignCAD: Warning: This module is designed for use with "
            "starccm+ meshes only, it also requires that the star mesh was  "
            "created in a certain way" << endl;

    if(!m_config["order"].beenSet)
    {
        cout << "Mesh order not set will assume 4" << endl;
    }

    m_mesh->m_nummode = m_config["order"].as<int>() + 1;

    //load instance of the CAD model
    ModuleSharedPtr module = GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "loadcad"), m_mesh);
    module->RegisterConfig("filename", m_config["file"].as<string>());
    module->SetDefaults();
    module->Process();

    //turn the bounding box of the CAD surfaces into a k-d tree
    vector<boxI> boxes;
    for (int i = 1; i <= m_mesh->m_cad->GetNumSurf(); i++)
    {
        LibUtilities::PrintProgressbar(i, m_mesh->m_cad->GetNumSurf(), "building surface bboxes", i-1);
        Array<OneD, NekDouble> bx = m_mesh->m_cad->GetSurf(i)->BoundingBox();
        boxes.push_back(make_pair(box(point(bx[0], bx[1], bx[2]), point(bx[3], bx[4], bx[5])), i));
    }
    cout << endl;

    cout << "building admin data structures" << endl;

    bgi::rtree<boxI, bgi::quadratic<16> > rtree(boxes);

    NodeSet surfNodes;

    //find unique nodes on the surface
    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[2][i];
        vector<NodeSharedPtr> ns = el->GetVertexList();
        for(int j = 0; j < ns.size(); j++)
        {
            surfNodes.insert(ns[j]);
        }
    }

    map<NodeSharedPtr, vector<ElementSharedPtr> > surfNodeToEl;

    //link surface nodes to their 3D element
    for(int i = 0; i < m_mesh->m_element[3].size(); i++)
    {
        if(m_mesh->m_element[3][i]->HasBoundaryLinks())
        {
            vector<NodeSharedPtr> ns = m_mesh->m_element[3][i]->GetVertexList();
            for(int j = 0; j < ns.size(); j++)
            {
                if(surfNodes.count(ns[j]) > 0)
                {
                    surfNodeToEl[ns[j]].push_back(m_mesh->m_element[3][i]);
                }
            }
        }
    }

    //link the surface node to a value for the shortest connecting edge to it
    map<NodeSharedPtr, NekDouble> minConEdge;
    for (int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[2][i];
        vector<NodeSharedPtr> ns = el->GetVertexList();
        NekDouble l1 = ns[0]->Distance(ns[1]);
        NekDouble l2 = ns[1]->Distance(ns[2]);
        NekDouble l3 = ns[2]->Distance(ns[0]);

        if (minConEdge.count(ns[0]))
        {
            NekDouble l = minConEdge[ns[0]];
            minConEdge[ns[0]] = min(l, min(l1, l3));
        }
        else
        {
            minConEdge.insert(make_pair(ns[0], min(l1, l3)));
        }
        if (minConEdge.count(ns[1]))
        {
            NekDouble l = minConEdge[ns[1]];
            minConEdge[ns[1]] = min(l, min(l1, l1));
        }
        else
        {
            minConEdge.insert(make_pair(ns[1], min(l1, l1)));
        }
        if (minConEdge.count(ns[2]))
        {
            NekDouble l = minConEdge[ns[2]];
            minConEdge[ns[2]] = min(l, min(l2, l3));
        }
        else
        {
            minConEdge.insert(make_pair(ns[2], min(l2, l3)));
        }
    }

    map<int, vector<int> > finds;

    cout << "searching tree" << endl;

    NekDouble maxNodeCor = 0;

    // this is a set of nodes which have a CAD failure
    // if touched in the HO stage they should be ignored and linearised
    NodeSet lockedNodes;

    //find nodes surface and parametric location
    int ct = 0;
    for (auto i = surfNodes.begin(); i != surfNodes.end(); i++, ct++)
    {
        LibUtilities::PrintProgressbar(ct, surfNodes.size(), "projecting verts", ct-1);

        point q((*i)->m_x, (*i)->m_y, (*i)->m_z);
        vector<boxI> result;
        rtree.query(bgi::intersects(q), back_inserter(result));

        ASSERTL0(result.size() > 0, "node thinks its not in any boxes");

        NekDouble tol = minConEdge[*i] * 0.5;
        tol = min(tol, 5e-4);
        tol = max(tol, 1e-6);

        vector<int> distId;
        vector<NekDouble> distList;

        for (int j = 0; j < result.size(); j++)
        {
            NekDouble dist;
            m_mesh->m_cad->GetSurf(result[j].second)->locuv((*i)->GetLoc(), dist);
            distList.push_back(dist);
            distId.push_back(result[j].second);
        }

        bool repeat = true;
        while (repeat)
        {
            repeat = false;
            for (int j = 0; j < distId.size() - 1; j++)
            {
                if (distList[j+1] < distList[j])
                {
                    repeat = true;
                    swap(distList[j+1], distList[j]);
                    swap(distId[j+1], distId[j]);
                }
            }
        }

        int pos = 0;
        for (int j = 0; j < distId.size(); j++)
        {
            if (distList[j] < tol)
            {
                pos++;
            }
        }

        distId.resize(pos);

        finds[pos].push_back(0);

        if (pos == 0)
        {
            lockedNodes.insert(*i);
            cout << endl << "WARNING: surface unknown " << distList[0] << " " << tol <<  endl;
        }
        else
        {
            NekDouble shift;
            bool st = false;
            for (int j = 0; j < distId.size(); j++)
            {
                if (distList[j] > tol)
                {
                    continue;
                }
                if (m_mesh->m_cad->GetSurf(distId[j])->IsPlanar())
                {
                    continue;
                }

                shift = distList[j];
                Array<OneD, NekDouble> uvt(2);
                NekDouble dist = 0;
                CADSurfSharedPtr s = m_mesh->m_cad->GetSurf(distId[j]);
                Array<OneD, NekDouble> l = (*i)->GetLoc();
                uvt = s->locuv(l, dist);

                NekDouble tmpX = (*i)->m_x;
                NekDouble tmpY = (*i)->m_y;
                NekDouble tmpZ = (*i)->m_z;

                (*i)->m_x = l[0];
                (*i)->m_y = l[1];
                (*i)->m_z = l[2];

                if(IsNotValid(surfNodeToEl[*i]))
                {
                    (*i)->m_x = tmpX;
                    (*i)->m_y = tmpY;
                    (*i)->m_z = tmpZ;
                    cout << "reset element projection" << endl;
                    break;
                }

                st = true;
                break;
            }

            if (!st)
            {
                lockedNodes.insert(*i);
                continue;
            }

            for (int j = 0; j < distId.size(); j++)
            {
                if (distList[j] > tol)
                {
                    continue;
                }
                if (m_mesh->m_cad->GetSurf(distId[j])->IsPlanar())
                {
                    continue;
                }

                CADSurfSharedPtr s = m_mesh->m_cad->GetSurf(distId[j]);
                Array<OneD, NekDouble> uv;
                NekDouble dist=0;
                Array<OneD, NekDouble> loc = (*i)->GetLoc();
                uv = s->locuv(loc, dist);
                (*i)->SetCADSurf(s, uv);
            }
            maxNodeCor = max(maxNodeCor, shift);
        }
    }

    cout << endl;

    cout << "max surface node correction " << maxNodeCor << endl;

    //make edges of surface mesh unique
    ClearElementLinks();
    EdgeSet surfEdges;
    vector<ElementSharedPtr> &elmt = m_mesh->m_element[2];
    map<int, int> surfIdToLoc;
    for (int i = 0; i < elmt.size(); i++)
    {
        surfIdToLoc.insert(make_pair(elmt[i]->GetId(), i));
        for (int j = 0; j < elmt[i]->GetEdgeCount(); ++j)
        {
            pair<EdgeSet::iterator,bool> testIns;
            EdgeSharedPtr ed = elmt[i]->GetEdge(j);
            testIns = surfEdges.insert(ed);

            if (testIns.second)
            {
                EdgeSharedPtr ed2 = *testIns.first;
                ed2->m_elLink.push_back(
                    pair<ElementSharedPtr,int>(elmt[i],j));
            }
            else
            {
                EdgeSharedPtr e2 = *(testIns.first);
                elmt[i]->SetEdge(j, e2);

                // Update edge to element map.
                e2->m_elLink.push_back( pair<ElementSharedPtr,int>(elmt[i],j));
            }
        }
    }

    int order = m_config["order"].as<int>();

    map<int, vector<int> > eds;

    LibUtilities::PointsKey ekey(order+1, LibUtilities::eGaussLobattoLegendre);
    Array<OneD, NekDouble> gll;
    LibUtilities::PointsManager()[ekey]->GetPoints(gll);

    //make surface edges high-order
    for (auto i = surfEdges.begin(); i != surfEdges.end(); i++)
    {
        if (lockedNodes.count((*i)->m_n1) || lockedNodes.count((*i)->m_n2))
        {
            continue;
        }
        vector<CADSurfSharedPtr> v1 = (*i)->m_n1->GetCADSurfs();
        vector<CADSurfSharedPtr> v2 = (*i)->m_n2->GetCADSurfs();

        vector<int> vi1, vi2;
        for (size_t j = 0; j < v1.size(); ++j)
        {
            vi1.push_back(v1[j]->GetId());
        }
        for (size_t j = 0; j < v2.size(); ++j)
        {
            vi2.push_back(v2[j]->GetId());
        }

        sort(vi1.begin(), vi1.end());
        sort(vi2.begin(), vi2.end());

        vector<int> cmn;
        set_intersection(vi1.begin(), vi1.end(), vi2.begin(), vi2.end(), back_inserter(cmn));
        eds[cmn.size()].push_back(0);

        (*i)->m_curveType = LibUtilities::eGaussLobattoLegendre;

        if (cmn.size() == 1 || cmn.size() == 2)
        {
            for (int j = 0; j < cmn.size(); j++)
            {
                if (m_mesh->m_cad->GetSurf(cmn[j])->IsPlanar())
                {
                    // if its planar dont care
                    continue;
                }

                Array<OneD, NekDouble> uvb = (*i)->m_n1->GetCADSurfInfo(cmn[j]);
                Array<OneD, NekDouble> uve = (*i)->m_n2->GetCADSurfInfo(cmn[j]);

                // can compare the loction of the projection to the
                // corresponding position of the straight sided edge
                // if the two differ by more than the length of the edge
                // something has gone wrong
                NekDouble len = (*i)->m_n1->Distance((*i)->m_n2);

                for (int k = 1; k < order+1 - 1; k++)
                {
                    Array<OneD, NekDouble> uv(2);
                    uv[0] = uvb[0] * (1.0 - gll[k]) / 2.0 + uve[0] * (1.0 + gll[k]) / 2.0;
                    uv[1] = uvb[1] * (1.0 - gll[k]) / 2.0 + uve[1] * (1.0 + gll[k]) / 2.0;
                    Array<OneD, NekDouble> loc;
                    loc              = m_mesh->m_cad->GetSurf(cmn[j])->P(uv);
                    Array<OneD, NekDouble> locT(3);
                    locT[0] = (*i)->m_n1->m_x * (1.0 - gll[k]) / 2.0 +
                    (*i)->m_n2->m_x * (1.0 + gll[k]) / 2.0;
                    locT[1] = (*i)->m_n1->m_y * (1.0 - gll[k]) / 2.0 +
                    (*i)->m_n2->m_y * (1.0 + gll[k]) / 2.0;
                    locT[2] = (*i)->m_n1->m_z * (1.0 - gll[k]) / 2.0 +
                    (*i)->m_n2->m_z * (1.0 + gll[k]) / 2.0;

                    NekDouble d = sqrt((locT[0] - loc[0]) * (locT[0] - loc[0]) +
                        (locT[1] - loc[1]) * (locT[1] - loc[1]) +
                        (locT[2] - loc[2]) * (locT[2] - loc[2]));

                    if (d > len)
                    {
                        (*i)->m_edgeNodes.clear();
                        break;
                    }

                    NodeSharedPtr nn = std::shared_ptr<Node>(
                        new Node(0, loc[0], loc[1], loc[2]));

                    (*i)->m_edgeNodes.push_back(nn);
                }

                if ((*i)->m_edgeNodes.size() != 0)
                {
                    // it suceeded on this surface so skip the other possibility
                    break;
                }
            }
        }
        else if (cmn.size() == 0)
        {
            // projection, if the projection requires more than two surfaces
            // including the edge nodes, then,  in theory projection shouldnt be
            // used

            set<int> sused;
            for (int k = 1; k < order+1 - 1; k++)
            {
                Array<OneD, NekDouble> locT(3);
                locT[0] = (*i)->m_n1->m_x * (1.0 - gll[k]) / 2.0 +
                          (*i)->m_n2->m_x * (1.0 + gll[k]) / 2.0;
                locT[1] = (*i)->m_n1->m_y * (1.0 - gll[k]) / 2.0 +
                          (*i)->m_n2->m_y * (1.0 + gll[k]) / 2.0;
                locT[2] = (*i)->m_n1->m_z * (1.0 - gll[k]) / 2.0 +
                          (*i)->m_n2->m_z * (1.0 + gll[k]) / 2.0;

                int s;
                if(!findAndProject(rtree, locT, s))
                {
                    (*i)->m_edgeNodes.clear();
                    break;
                }
                sused.insert(s);

                if (sused.size() > 2)
                {
                    (*i)->m_edgeNodes.clear();
                    break;
                }

                NodeSharedPtr nn = std::shared_ptr<Node>(
                    new Node(0, locT[0], locT[1], locT[2]));

                (*i)->m_edgeNodes.push_back(nn);
            }
        }
    }

}
}
}
