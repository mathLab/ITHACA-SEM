////////////////////////////////////////////////////////////////////////////////
//
//  File: Edge.h
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
//  Description: Mesh Edge.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESHUTILS_MESHELEMENTS_EDGE
#define NEKMESHUTILS_MESHELEMENTS_EDGE

#include <LibUtilities/BasicUtils/HashUtils.hpp>
#include <SpatialDomains/SegGeom.h>

#include <NekMeshUtils/NekMeshUtilsDeclspec.h>
#include <NekMeshUtils/MeshElements/Node.h>

namespace Nektar
{
namespace NekMeshUtils
{

class Element;
typedef std::shared_ptr<Element> ElementSharedPtr;

/**
 * @brief Represents an edge which joins two points.
 *
 * An edge is defined by two nodes (vertices) and, for high-order edges,
 * a set of control nodes defining the shape of the edge.
 */
class Edge
{
public:
    /// Creates a new edge.
    NEKMESHUTILS_EXPORT Edge(NodeSharedPtr              pVertex1,
                             NodeSharedPtr              pVertex2,
                             std::vector<NodeSharedPtr> pEdgeNodes,
                             LibUtilities::PointsType   pCurveType)
        : m_n1(pVertex1), m_n2(pVertex2), m_edgeNodes(pEdgeNodes),
          m_curveType(pCurveType), m_geom(){}

    /// Creates a new linear edge.
    NEKMESHUTILS_EXPORT Edge(NodeSharedPtr pVertex1, NodeSharedPtr pVertex2)
        : m_n1(pVertex1), m_n2(pVertex2), m_edgeNodes(), m_curveType(), m_geom()
    {}

    /// Copies an existing edge.
    NEKMESHUTILS_EXPORT Edge(const Edge &pSrc)
        : m_id(pSrc.m_id), m_n1(pSrc.m_n1), m_n2(pSrc.m_n2),
          m_edgeNodes(pSrc.m_edgeNodes), m_curveType(pSrc.m_curveType),
          m_elLink(pSrc.m_elLink), m_parentCAD(pSrc.m_parentCAD)
    {
    }

    NEKMESHUTILS_EXPORT ~Edge()
    {
    }

    /// Returns the total number of nodes defining the edge.
    NEKMESHUTILS_EXPORT unsigned int GetNodeCount() const
    {
        return m_edgeNodes.size() + 2;
    }

    NEKMESHUTILS_EXPORT void GetCurvedNodes(
        std::vector<NodeSharedPtr> &nodeList) const
    {
        nodeList.push_back(m_n1);
        for (int k = 0; k < m_edgeNodes.size(); ++k)
        {
            nodeList.push_back(m_edgeNodes[k]);
        }
        nodeList.push_back(m_n2);
    }

    /// Creates a Nektar++ string listing the coordinates of all the
    /// nodes.
    NEKMESHUTILS_EXPORT std::string GetXmlCurveString();

    /// Generate a SpatialDomains::SegGeom object for this edge.
    NEKMESHUTILS_EXPORT SpatialDomains::SegGeomSharedPtr GetGeom(int coordDim);

    /// Make this edge an order @p order edge. @see Element::MakeOrder.
    void MakeOrder(int                                order,
                   SpatialDomains::GeometrySharedPtr  geom,
                   LibUtilities::PointsType           edgeType,
                   int                                coordDim,
                   int                               &id);

    /// ID of edge.
    unsigned int m_id;
    /// First vertex node.
    NodeSharedPtr m_n1;
    /// Second vertex node.
    NodeSharedPtr m_n2;
    /// List of control nodes between the first and second vertices.
    std::vector<NodeSharedPtr> m_edgeNodes;
    /// Distributions of points along edge.
    LibUtilities::PointsType m_curveType;
    /// Element(s) which are linked to this edge.
    std::vector<std::pair<std::weak_ptr<Element>, int> > m_elLink;

    CADObjectSharedPtr m_parentCAD;

private:
    SpatialDomains::SegGeomSharedPtr m_geom;
};
/// Shared pointer to an edge.
typedef std::shared_ptr<Edge> EdgeSharedPtr;

NEKMESHUTILS_EXPORT bool operator==(EdgeSharedPtr const &p1,
                                    EdgeSharedPtr const &p2);
NEKMESHUTILS_EXPORT bool operator<(EdgeSharedPtr const &p1,
                                   EdgeSharedPtr const &p2);

/**
 * @brief Defines a hash function for edges.
 *
 * The hash of an edge is defined using the IDs of the two nodes which
 * define it. First the minimum ID is hashed, then the maximum
 * ID, which takes the two possible orientations into account.
 */
struct EdgeHash : std::unary_function<EdgeSharedPtr, std::size_t>
{
    std::size_t operator()(EdgeSharedPtr const &p) const
    {
        const unsigned int id1 = p->m_n1->m_id;
        const unsigned int id2 = p->m_n2->m_id;
        return id1 < id2 ? hash_combine(id1, id2) : hash_combine(id2, id1);
    }
};
typedef std::unordered_set<EdgeSharedPtr, EdgeHash> EdgeSet;
}
}

#endif
