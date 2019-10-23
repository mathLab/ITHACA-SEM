////////////////////////////////////////////////////////////////////////////////
//
//  File: Line.cpp
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
//  Description: Mesh line object.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/MeshElements/Line.h>

#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

LibUtilities::ShapeType Line::m_type =
    GetElementFactory().RegisterCreatorFunction(
        LibUtilities::eSegment, Line::create, "Line");

/**
 * @brief Create a line element.
 */
Line::Line(ElmtConfig pConf,
           vector<NodeSharedPtr> pNodeList,
           vector<int> pTagList)
    : Element(pConf, GetNumNodes(pConf), pNodeList.size())
{
    m_tag     = "S";
    m_dim     = 1;
    m_taglist = pTagList;
    int n     = m_conf.m_order - 1;

    // Add vertices
    for (int i = 0; i < 2; ++i)
    {
        m_vertex.push_back(pNodeList[i]);
    }

    if (m_conf.m_order > 1)
    {
        for (int j = 0; j < n; ++j)
        {
            m_volumeNodes.push_back(pNodeList[2 + j]);
        }
    }
}

SpatialDomains::GeometrySharedPtr Line::GetGeom(int coordDim)
{
    // Create edge vertices.
    SpatialDomains::PointGeomSharedPtr p[2];
    SpatialDomains::SegGeomSharedPtr ret;

    p[0] = m_vertex[0]->GetGeom(coordDim);
    p[1] = m_vertex[1]->GetGeom(coordDim);

    if (m_edge[0]->m_edgeNodes.size() > 0)
    {
        SpatialDomains::CurveSharedPtr c =
            MemoryManager<SpatialDomains::Curve>::AllocateSharedPtr(
                m_id, m_edge[0]->m_curveType);

        c->m_points.push_back(p[0]);
        for (int i = 0; i < m_edge[0]->m_edgeNodes.size(); ++i)
        {
            c->m_points.push_back(m_edge[0]->m_edgeNodes[i]->GetGeom(coordDim));
        }
        c->m_points.push_back(p[1]);

        ret = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(
            m_id, 2, p, c);
    }
    else
    {
        ret = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(
            m_id, 2, p);
    }

    ret->Setup();
    return ret;
}


void Line::GetCurvedNodes(std::vector<NodeSharedPtr> &nodeList) const
{
    nodeList.push_back(m_vertex[0]);
    for (int i = 0; i < m_volumeNodes.size(); ++i)
    {
        nodeList.push_back(m_volumeNodes[i]);
    }
    nodeList.push_back(m_vertex[1]);
}
void Line::MakeOrder(int                                order,
                     SpatialDomains::GeometrySharedPtr  geom,
                     LibUtilities::PointsType           pType,
                     int                                coordDim,
                     int                               &id,
                     bool                               justConfig)
{
    m_conf.m_order       = order;
    m_curveType          = pType;
    m_conf.m_volumeNodes = false;
    m_volumeNodes.clear();

    // Lines of order == 1 have no interior volume points.
    if (order == 1)
    {
        m_conf.m_faceNodes = false;
        return;
    }

    m_conf.m_faceNodes = true;

    if (justConfig)
    {
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
    m_volumeNodes.resize(nQuadIntPts);

    for (int i = 1, cnt = 0; i < nPoints-1; ++i)
    {
        Array<OneD, NekDouble> xp(1);
        xp[0] = px[i];

        Array<OneD, NekDouble> x(3, 0.0);
        for (int k = 0; k < coordDim; ++k)
        {
            x[k] = xmap->PhysEvaluate(xp, phys[k]);
        }

        m_volumeNodes[cnt] = std::shared_ptr<Node>(
            new Node(id++, x[0], x[1], x[2]));
    }
}

/**
 * @brief Return the number of nodes defining a line.
 */
unsigned int Line::GetNumNodes(ElmtConfig pConf)
{
    return pConf.m_order + 1;
}
}
}
