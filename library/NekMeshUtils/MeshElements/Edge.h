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
//  Description: Mesh manipulation objects.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NekMeshUtils_MESHELEMENTS_EDGE
#define NekMeshUtils_MESHELEMENTS_EDGE

#include <iomanip>

#include <SpatialDomains/SegGeom.h>

#include <NekMeshUtils/NekMeshUtilsDeclspec.h>
#include <NekMeshUtils/MeshElements/Element.h>
#include <NekMeshUtils/MeshElements/Node.h>

namespace Nektar
{
namespace NekMeshUtils
{

class Element;
typedef boost::shared_ptr<Element> ElementSharedPtr;

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
          m_curveType(pCurveType), m_geom()
    {
#ifdef NEKTAR_USE_MESHGEN
        onCurve = false;
#endif
    }

    /// Creates a new linear edge.
    NEKMESHUTILS_EXPORT Edge(NodeSharedPtr pVertex1, NodeSharedPtr pVertex2)
        : m_n1(pVertex1), m_n2(pVertex2), m_edgeNodes(), m_curveType(), m_geom()
    {
#ifdef NEKTAR_USE_MESHGEN
        onCurve = false;
#endif
    }

    /// Copies an existing edge.
    NEKMESHUTILS_EXPORT Edge(const Edge &pSrc)
        : m_n1(pSrc.m_n1), m_n2(pSrc.m_n2), m_edgeNodes(pSrc.m_edgeNodes),
          m_curveType(pSrc.m_curveType), m_geom(pSrc.m_geom)
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
    NEKMESHUTILS_EXPORT std::string GetXmlCurveString() const
    {
        std::vector<NodeSharedPtr> nodeList;

        GetCurvedNodes(nodeList);

        std::stringstream s;
        std::string str;

        // put them into a string and return
        for (int k = 0; k < nodeList.size(); ++k)
        {
            s << std::scientific << std::setprecision(8) << "    "
              << nodeList[k]->m_x << "  " << nodeList[k]->m_y << "  "
              << nodeList[k]->m_z << "    ";
        }

        return s.str();
    }

    /// Generate a SpatialDomains::SegGeom object for this edge.
    NEKMESHUTILS_EXPORT SpatialDomains::SegGeomSharedPtr GetGeom(int coordDim)
    {
        // Create edge vertices.
        SpatialDomains::PointGeomSharedPtr p[2];
        SpatialDomains::SegGeomSharedPtr ret;

        p[0] = m_n1->GetGeom(coordDim);
        p[1] = m_n2->GetGeom(coordDim);

        // Create a curve if high-order information exists.
        if (m_edgeNodes.size() > 0)
        {
            SpatialDomains::CurveSharedPtr c =
                MemoryManager<SpatialDomains::Curve>::AllocateSharedPtr(
                    m_id, m_curveType);

            c->m_points.push_back(p[0]);
            for (int i = 0; i < m_edgeNodes.size(); ++i)
            {
                c->m_points.push_back(m_edgeNodes[i]->GetGeom(coordDim));
            }
            c->m_points.push_back(p[1]);

            ret = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(
                m_id, coordDim, p, c);
        }
        else
        {
            ret = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(
                m_id, coordDim, p);
        }

        return ret;
    }

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
    std::vector<std::pair<ElementSharedPtr, int> > m_elLink;

#ifdef NEKTAR_USE_MESHGEN
    bool onCurve;
    /// id of cad curve which edge lies on
    int CADCurveId;
    CADCurveSharedPtr CADCurve;
#endif

private:
    SpatialDomains::SegGeomSharedPtr m_geom;
};
/// Shared pointer to an edge.
typedef boost::shared_ptr<Edge> EdgeSharedPtr;

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
        std::size_t seed = 0;
        unsigned int id1 = p->m_n1->m_id;
        unsigned int id2 = p->m_n2->m_id;
        boost::hash_combine(seed, id1 < id2 ? id1 : id2);
        boost::hash_combine(seed, id2 < id1 ? id1 : id2);
        return seed;
    }
};
typedef boost::unordered_set<EdgeSharedPtr, EdgeHash> EdgeSet;
}
}

#endif
