////////////////////////////////////////////////////////////////////////////////
//
//  File: Face.cpp
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
//  Description: Mesh face object.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/MeshElements/Face.h>

#include <NekMeshUtils/CADSystem/CADSurf.h>

#include <LibUtilities/Foundations/ManagerAccess.h>

namespace Nektar
{
namespace NekMeshUtils
{

void Face::GetCurvedNodes(
    std::vector<NodeSharedPtr> &nodeList)
{
    // Treat 2D point distributions differently to 3D.
    if (m_curveType == LibUtilities::eNodalTriFekete ||
        m_curveType == LibUtilities::eNodalTriEvenlySpaced ||
        m_curveType == LibUtilities::eNodalTriElec)
    {
        int n = m_edgeList[0]->GetNodeCount();
        int n2 = m_edgeList[1]->GetNodeCount();
        int n3 = m_edgeList[2]->GetNodeCount();

        bool same = (n == n2 ? (n2 == n3) : false);
        ASSERTL0(same, "Edges are not consistent");

        nodeList.insert(
            nodeList.end(), m_vertexList.begin(), m_vertexList.end());
        for (int k = 0; k < m_edgeList.size(); ++k)
        {
            nodeList.insert(nodeList.end(),
                            m_edgeList[k]->m_edgeNodes.begin(),
                            m_edgeList[k]->m_edgeNodes.end());
            if (m_edgeList[k]->m_n1 != m_vertexList[k])
            {
                // If edge orientation is reversed relative to node
                // ordering, we need to reverse order of nodes.
                std::reverse(nodeList.begin() + 3 + k * (n - 2),
                             nodeList.begin() + 3 + (k + 1) * (n - 2));
            }
        }
        nodeList.insert(
            nodeList.end(), m_faceNodes.begin(), m_faceNodes.end());
    }
    else
    {
        // Write out in 2D tensor product order.
        ASSERTL0(m_vertexList.size() == 4,
                 "Face nodes of tensor product only supported "
                 "for quadrilaterals.");

        int n = (int)sqrt((NekDouble)GetNodeCount());
        nodeList.resize(n * n);

        ASSERTL0(n * n == GetNodeCount(), "Wrong number of modes?");

        // Write vertices
        nodeList[0]         = m_vertexList[0];
        nodeList[n - 1]     = m_vertexList[1];
        nodeList[n * n - 1] = m_vertexList[2];
        nodeList[n * (n - 1)] = m_vertexList[3];

        // Write edge-interior
        int skips[4][2] = {
            {0, 1}, {n - 1, n}, {n * n - 1, -1}, {n * (n - 1), -n}};
        for (int i = 0; i < 4; ++i)
        {
            bool reverseEdge = m_edgeList[i]->m_n1 == m_vertexList[i];

            if (!reverseEdge)
            {
                for (int j = 1; j < n - 1; ++j)
                {
                    nodeList[skips[i][0] + j * skips[i][1]] =
                        m_edgeList[i]->m_edgeNodes[n - 2 - j];
                }
            }
            else
            {
                for (int j = 1; j < n - 1; ++j)
                {
                    nodeList[skips[i][0] + j * skips[i][1]] =
                        m_edgeList[i]->m_edgeNodes[j - 1];
                }
            }
        }

        // Write interior
        for (int i = 1; i < n - 1; ++i)
        {
            for (int j = 1; j < n - 1; ++j)
            {
                nodeList[i * n + j] =
                    m_faceNodes[(i - 1) * (n - 2) + (j - 1)];
            }
        }
    }
}


std::string Face::GetXmlCurveString()
{
    std::stringstream s;
    std::string str;
    std::vector<NodeSharedPtr> nodeList;

    // assemble listof nodes
    GetCurvedNodes(nodeList);

    // put them into a string
    for (int k = 0; k < nodeList.size(); ++k)
    {
        s << std::scientific << std::setprecision(8) << "    "
          << nodeList[k]->m_x << "  " << nodeList[k]->m_y << "  "
          << nodeList[k]->m_z << "    ";
    }

    return s.str();
}

void Face::MakeOrder(int                                order,
                     SpatialDomains::GeometrySharedPtr  geom,
                     LibUtilities::PointsType           pType,
                     int                                coordDim,
                     int                               &id)
{
    if (m_vertexList.size() == 3)
    {
        // Triangles of order < 3 have no interior volume points.
        if (order < 3)
        {
            m_faceNodes.clear();
            return;
        }

        int nPoints = order + 1;
        StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();

        Array<OneD, NekDouble> px, py;
        LibUtilities::PointsKey pKey(nPoints, pType);
        ASSERTL1(pKey.GetPointsDim() == 2, "Points distribution must be 2D");
        LibUtilities::PointsManager()[pKey]->GetPoints(px, py);

        Array<OneD, Array<OneD, NekDouble> > phys(coordDim);

        for (int i = 0; i < coordDim; ++i)
        {
            phys[i] = Array<OneD, NekDouble>(xmap->GetTotPoints());
            xmap->BwdTrans(geom->GetCoeffs(i), phys[i]);
        }

        const int nTriPts = nPoints * (nPoints + 1) / 2;
        const int nTriIntPts = (nPoints - 3) * (nPoints - 2) / 2;
        m_faceNodes.resize(nTriIntPts);

        for (int i = 3 + 3*(nPoints-2), cnt = 0; i < nTriPts; ++i, ++cnt)
        {
            Array<OneD, NekDouble> xp(2);
            xp[0] = px[i];
            xp[1] = py[i];

            Array<OneD, NekDouble> x(3, 0.0);
            for (int j = 0; j < coordDim; ++j)
            {
                x[j] = xmap->PhysEvaluate(xp, phys[j]);
            }

            m_faceNodes[cnt] = std::shared_ptr<Node>(
                new Node(id++, x[0], x[1], x[2]));
        }
        m_curveType = pType;
    }
    else if (m_vertexList.size() == 4)
    {
        // Quads of order < 2 have no interior volume points.
        if (order < 2)
        {
            m_faceNodes.clear();
            return;
        }

        int nPoints = order + 1;
        StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();

        Array<OneD, NekDouble> px;
        LibUtilities::PointsKey pKey(nPoints, pType);
        ASSERTL1(pKey.GetPointsDim() == 1, "Points distribution must be 1D");
        LibUtilities::PointsManager()[pKey]->GetPoints(px);

        Array<OneD, Array<OneD, NekDouble> > phys(coordDim);

        for (int i = 0; i < coordDim; ++i)
        {
            phys[i] = Array<OneD, NekDouble>(xmap->GetTotPoints());
            xmap->BwdTrans(geom->GetCoeffs(i), phys[i]);
        }

        int nQuadIntPts = (nPoints - 2) * (nPoints - 2);
        m_faceNodes.resize(nQuadIntPts);

        for (int i = 1, cnt = 0; i < nPoints-1; ++i)
        {
            for (int j = 1; j < nPoints-1; ++j, ++cnt)
            {
                Array<OneD, NekDouble> xp(2);
                xp[0] = px[j];
                xp[1] = px[i];

                Array<OneD, NekDouble> x(3, 0.0);
                for (int k = 0; k < coordDim; ++k)
                {
                    x[k] = xmap->PhysEvaluate(xp, phys[k]);
                }

                m_faceNodes[cnt] = std::shared_ptr<Node>(
                    new Node(id++, x[0], x[1], x[2]));
            }
        }

        m_curveType = pType;
    }
    else
    {
        ASSERTL0(false, "Unknown number of vertices");
    }

    if(m_parentCAD)
    {
        CADSurfSharedPtr s = std::dynamic_pointer_cast<CADSurf>(m_parentCAD);
        for(int i = 0; i < m_faceNodes.size(); i++)
        {
            Array<OneD, NekDouble> loc(3);
            loc[0] = m_faceNodes[i]->m_x;
            loc[1] = m_faceNodes[i]->m_y;
            loc[2] = m_faceNodes[i]->m_z;
            Array<OneD, NekDouble> uv = s->locuv(loc);
            loc = s->P(uv);
            m_faceNodes[i]->m_x = loc[0];
            m_faceNodes[i]->m_y = loc[1];
            m_faceNodes[i]->m_z = loc[2];
            m_faceNodes[i]->SetCADSurf(s,uv);
        }
    }
}

SpatialDomains::Geometry2DSharedPtr Face::GetGeom(int coordDim)
{
    int nEdge = m_edgeList.size();

    SpatialDomains::SegGeomSharedPtr edges[4];
    SpatialDomains::Geometry2DSharedPtr ret;
    StdRegions::Orientation edgeo[4];

    for (int i = 0; i < nEdge; ++i)
    {
        edges[i] = m_edgeList[i]->GetGeom(coordDim);
    }

    for (int i = 0; i < nEdge; ++i)
    {
        edgeo[i] = m_edgeList[i]->m_n1 == m_vertexList[i]
                       ? StdRegions::eForwards
                       : StdRegions::eBackwards;
    }

    if (m_faceNodes.size() > 0)
    {
        if (nEdge == 3)
        {
            SpatialDomains::CurveSharedPtr c =
                MemoryManager<SpatialDomains::Curve>::AllocateSharedPtr(
                    m_id, m_curveType);

            for (int j = 0; j < m_vertexList.size(); j++)
            {
                c->m_points.push_back(m_vertexList[j]->GetGeom(coordDim));
            }
            for (int j = 0; j < nEdge; j++)
            {
                std::vector<NodeSharedPtr> ed = m_edgeList[j]->m_edgeNodes;
                if (edgeo[j] == StdRegions::eBackwards)
                {
                    for (int k = ed.size() - 1; k >= 0; k--)
                    {
                        c->m_points.push_back(ed[k]->GetGeom(coordDim));
                    }
                }
                else
                {
                    for (int k = 0; k < ed.size(); k++)
                    {
                        c->m_points.push_back(ed[k]->GetGeom(coordDim));
                    }
                }
            }
            for (int j = 0; j < m_faceNodes.size(); j++)
            {
                c->m_points.push_back(m_faceNodes[j]->GetGeom(coordDim));
            }

            ret = MemoryManager<SpatialDomains::TriGeom>::AllocateSharedPtr(
                m_id, edges, c);
        }
        else
        {
            SpatialDomains::CurveSharedPtr c =
                MemoryManager<SpatialDomains::Curve>::AllocateSharedPtr(
                    m_id, m_curveType);

            ASSERTL0(m_vertexList.size() == 4,
                     "Face nodes of tensor product only supported "
                     "for quadrilaterals.");

            int n = (int)sqrt((NekDouble)GetNodeCount());
            std::vector<NodeSharedPtr> tmp(n * n);

            ASSERTL0(n * n == GetNodeCount(), "Wrong number of modes?");

            // Write vertices
            tmp[0]         = m_vertexList[0];
            tmp[n - 1]     = m_vertexList[1];
            tmp[n * n - 1] = m_vertexList[2];
            tmp[n * (n - 1)] = m_vertexList[3];

            // Write edge-interior
            int skips[4][2] = {
                {0, 1}, {n - 1, n}, {n * n - 1, -1}, {n * (n - 1), -n}};
            for (int i = 0; i < 4; ++i)
            {
                bool reverseEdge = edgeo[i] == StdRegions::eBackwards;

                if (reverseEdge)
                {
                    for (int j = 1; j < n - 1; ++j)
                    {
                        tmp[skips[i][0] + j * skips[i][1]] =
                            m_edgeList[i]->m_edgeNodes[n - 2 - j];
                    }
                }
                else
                {
                    for (int j = 1; j < n - 1; ++j)
                    {
                        tmp[skips[i][0] + j * skips[i][1]] =
                            m_edgeList[i]->m_edgeNodes[j - 1];
                    }
                }
            }

            // Write interior
            for (int i = 1; i < n - 1; ++i)
            {
                for (int j = 1; j < n - 1; ++j)
                {
                    tmp[i * n + j] =
                        m_faceNodes[(i - 1) * (n - 2) + (j - 1)];
                }
            }

            for (int k = 0; k < tmp.size(); ++k)
            {
                c->m_points.push_back(tmp[k]->GetGeom(coordDim));
            }

            ret =
                MemoryManager<SpatialDomains::QuadGeom>::AllocateSharedPtr(
                    m_id, edges, c);
        }
    }
    else
    {
        if (nEdge == 3)
        {
            ret = MemoryManager<SpatialDomains::TriGeom>::AllocateSharedPtr(
                m_id, edges);
        }
        else
        {
            ret =
                MemoryManager<SpatialDomains::QuadGeom>::AllocateSharedPtr(
                    m_id, edges);
        }
    }

    ret->Setup();
    return ret;
}

}
}
