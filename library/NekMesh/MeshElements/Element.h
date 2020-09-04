////////////////////////////////////////////////////////////////////////////////
//
//  File: Element.h
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
//  Description: Mesh element.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESH_MESHELEMENTS_ELEMENT
#define NEKMESH_MESHELEMENTS_ELEMENT

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/Foundations/PointsType.h>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>

#include <NekMesh/NekMeshDeclspec.h>
#include <NekMesh/MeshElements/Face.h>
#include <NekMesh/MeshElements/ElementConfig.h>

namespace Nektar
{
namespace NekMesh
{

/**
 * @brief Base class for element definitions.
 *
 * An element is defined by a list of vertices, edges and faces
 * (depending on the dimension of the problem). This base class
 * provides the underlying structure.
 */
class Element
{
public:
    NEKMESH_EXPORT Element(
        ElmtConfig pConf, unsigned int pNumNodes, unsigned int pGotNodes);

    NEKMESH_EXPORT virtual ~Element() = default;

    /// Returns the ID of the element (or associated edge or face for
    /// boundary elements).
    unsigned int GetId() const
    {
        if (m_faceLink.get() != 0)
            return m_faceLink->m_id;
        if (m_edgeLink.get() != 0)
            return m_edgeLink->m_id;
        return m_id;
    }
    /// Returns the expansion dimension of the element.
    unsigned int GetDim() const
    {
        return m_dim;
    }
    /// Returns the configuration of the element.
    ElmtConfig GetConf() const
    {
        return m_conf;
    }
    ///returns the shapetype
    LibUtilities::ShapeType GetShapeType() const
    {
        return m_conf.m_e;
    }
    /// Returns the tag which defines the element shape.
    std::string GetTag() const
    {
        if (m_faceLink.get() != 0)
            return "F";
        if (m_edgeLink.get() != 0)
            return "E";
        return m_tag;
    }
    /// Access a vertex node.
    NodeSharedPtr GetVertex(unsigned int i) const
    {
        return m_vertex[i];
    }
    /// Access an edge.
    EdgeSharedPtr GetEdge(unsigned int i) const
    {
        return m_edge[i];
    }
    /// Access a face.
    FaceSharedPtr GetFace(unsigned int i) const
    {
        return m_face[i];
    }
    /// Access the list of vertex nodes.
    std::vector<NodeSharedPtr> GetVertexList() const
    {
        return m_vertex;
    }
    /// Access the list of edges.
    std::vector<EdgeSharedPtr> GetEdgeList() const
    {
        return m_edge;
    }
    /// Access the list of faces.
    std::vector<FaceSharedPtr> GetFaceList() const
    {
        return m_face;
    }
    /// Access the list of volume nodes.
    std::vector<NodeSharedPtr> GetVolumeNodes() const
    {
        return m_volumeNodes;
    }
    void SetVolumeNodes(std::vector<NodeSharedPtr> &nodes)
    {
        m_volumeNodes = nodes;
    }
    LibUtilities::PointsType GetCurveType() const
    {
        return m_curveType;
    }
    void SetCurveType(LibUtilities::PointsType cT)
    {
        m_curveType = cT;
    }
    /// Returns the total number of nodes (vertices, edge nodes and
    /// face nodes and volume nodes).
    NEKMESH_EXPORT unsigned int GetNodeCount();
    /// Access the list of tags associated with this element.
    std::vector<int> GetTagList() const
    {
        return m_taglist;
    }
    /// Returns the number of vertices.
    unsigned int GetVertexCount() const
    {
        return m_vertex.size();
    }
    /// Returns the number of edges.
    unsigned int GetEdgeCount() const
    {
        return m_edge.size();
    }
    /// Returns the number of faces.
    unsigned int GetFaceCount() const
    {
        return m_face.size();
    }
    /// Change the ID of the element.
    void SetId(unsigned int p)
    {
        m_id = p;
    }
    /**
     * @brief Replace a vertex in the element.
     *
     * When a vertex is replaced, the element edges and faces are also
     * searched and the corresponding edge/face nodes are updated to
     * maintain consistency.
     *
     * @param  p        Index of the vertex to replace.
     * @param  pNew     New vertex.
     * @param  descend  If true, we loop over edges and faces and replace the
     *                  corresponding vertices with @p pNew.
     */
    NEKMESH_EXPORT void SetVertex(
        unsigned int p, NodeSharedPtr pNew, bool descend = true);
    /**
     * @brief Replace an edge in the element.
     *
     * When an edge is replaced, the element faces are also searched and
     * the corresponding face edges are updated to maintain consistency.
     *
     * @param  p        Index of the edge to replace.
     * @param  pNew     New edge.
     * @param  descend  If true, we loop over faces and replace the corresponding
     *                  face edge with @p pNew.
     */
    NEKMESH_EXPORT void SetEdge(
        unsigned int p, EdgeSharedPtr pNew, bool descend = true);
    /**
     * @brief Replace a face in the element.
     *
     * When a face is replaced, no other consistency checks are required.
     *
     * @param  p     Index of the face to replace.
     * @param  pNew  New face.
     */
    NEKMESH_EXPORT void SetFace(unsigned int p, FaceSharedPtr pNew);
    /// Set a correspondence between this element and an edge
    /// (2D boundary element).
    void SetEdgeLink(EdgeSharedPtr pLink)
    {
        m_edgeLink = pLink;
    }
    /// Get correspondence between this element and an edge.
    EdgeSharedPtr GetEdgeLink()
    {
        return m_edgeLink;
    }
    /// Set a correspondence between this element and a face
    /// (3D boundary element).
    void SetFaceLink(FaceSharedPtr pLink)
    {
        m_faceLink = pLink;
    }
    /// Get correspondence between this element and a face.
    FaceSharedPtr GetFaceLink()
    {
        return m_faceLink;
    }
    /// Set a correspondence between edge or face i and its
    /// representative boundary element m->element[expDim-1][j].
    void SetBoundaryLink(int i, int j)
    {
        m_boundaryLinks[i] = j;
    }
    /// Get the location of the boundary face/edge i for this element.
    int GetBoundaryLink(int i)
    {
        std::map<int, int>::iterator it = m_boundaryLinks.find(i);
        if (it == m_boundaryLinks.end())
        {
            return -1;
        }
        else
        {
            return it->second;
        }
    }
    /// Is this element connected to a boundary
    bool HasBoundaryLinks()
    {
        return m_boundaryLinks.size() > 0;
    }
    /// Set the list of tags associated with this element.
    void SetTagList(const std::vector<int> &tags)
    {
        m_taglist = tags;
    }
    /// Generate a list of vertices (1D), edges (2D), or faces (3D).
    NEKMESH_EXPORT virtual std::string GetXmlString();
    /// get list of volume interior nodes
    virtual void GetCurvedNodes(
        std::vector<NodeSharedPtr> &nodeList) const
    {
        boost::ignore_unused(nodeList);
        NEKERROR(ErrorUtil::efatal,
                 "This function should be implemented at a shape level.");
    }

    /// Generates a string listing the coordinates of all nodes
    /// associated with this element.
    NEKMESH_EXPORT std::string GetXmlCurveString();
    /// Generate a Nektar++ geometry object for this element.
    virtual SpatialDomains::GeometrySharedPtr GetGeom(
        int coordDim)
    {
        boost::ignore_unused(coordDim);
        NEKERROR(ErrorUtil::efatal,
                 "This function should be implemented at a shape level.");
        return std::shared_ptr<SpatialDomains::Geometry>();
    }

    /**
     * @brief Obtain the order of an element by looking at edges.
     */
    NEKMESH_EXPORT int GetMaxOrder();

    /**
     * @brief Determines whether an element is deformed by inspecting whether
     * there are any edge, face or volume interior nodes.
     */
    bool IsDeformed()
    {
        if (m_volumeNodes.size() > 0)
        {
            return true;
        }

        for (auto &edge : m_edge)
        {
            if (edge->m_edgeNodes.size() > 0)
            {
                return true;
            }
        }

        for (auto &face : m_face)
        {
            if (face->m_faceNodes.size() > 0)
            {
                return true;
            }
        }

        return false;
    }

    /**
     * @brief Returns the approximate bounding box of this element based on the
     * coordinates of all vertices, edges and faces of the element. Note that
     * this does not robustly take into account the curvature of the element.
     */
    std::pair<Node, Node> GetBoundingBox()
    {
#define SWAP_NODE(a)                                            \
        lower.m_x = std::min(lower.m_x, a->m_x);                  \
        lower.m_y = std::min(lower.m_y, a->m_y);                  \
        lower.m_z = std::min(lower.m_z, a->m_z);                  \
        upper.m_x = std::max(upper.m_x, a->m_x);                  \
        upper.m_y = std::max(upper.m_y, a->m_y);                  \
        upper.m_z = std::max(upper.m_z, a->m_z);

        Node lower(*m_vertex[0]), upper(*m_vertex[0]);

        for (int i = 1; i < m_vertex.size(); ++i)
        {
            SWAP_NODE(m_vertex[i])
        }
        for (auto &edge : m_edge)
        {
            for (auto &edgeNode : edge->m_edgeNodes)
            {
                SWAP_NODE(edgeNode);
            }
        }
        for (auto &face : m_face)
        {
            for (auto &faceNode : face->m_faceNodes)
            {
                SWAP_NODE(faceNode);
            }
        }

        return std::make_pair(lower, upper);
#undef SWAP_NODE
    }

    /**
     * @brief Insert interior (i.e. volume) points into this element to make the
     * geometry an order @p order representation.
     *
     * @param order          The desired polynomial order.
     * @param geom           The geometry object used to describe the curvature
     *                       mapping.
     * @param edgeType       The points distribution to use on the volume.
     * @param coordDim       The coordinate (i.e. space) dimension.
     * @param id             Counter which should be incremented to supply
     *                       consistent vertex IDs.
     * @param justConfig     If true, then the configuration Element::m_conf
     *                       will be updated but no nodes will be
     *                       generated. This is used when considering boundary
     *                       elements, which just require copying of face or
     *                       edge interior nodes.
     */
    virtual void MakeOrder(
        int                                order,
        SpatialDomains::GeometrySharedPtr  geom,
        LibUtilities::PointsType           edgeType,
        int                                coordDim,
        int                               &id,
        bool                               justConfig = false)
    {
        boost::ignore_unused(order, geom, edgeType, coordDim, id, justConfig);
        NEKERROR(ErrorUtil::efatal,
                 "This function should be implemented at a shape level.");
    }

    /**
     * @brief Get the edge orientation of @p edge with respect to the local
     * element, which lies at edge index @p edgeId.
     */
    virtual StdRegions::Orientation GetEdgeOrient(
        int edgeId, EdgeSharedPtr edge)
    {
        boost::ignore_unused(edgeId, edge);
        NEKERROR(ErrorUtil::efatal,
                 "This function should be implemented at a shape level.");
        return StdRegions::eNoOrientation;
    }

    /**
     * @brief Returns the local index of vertex @p j of face @p i.
     */
    virtual int GetFaceVertex(int i, int j)
    {
        boost::ignore_unused(i, j);
        NEKERROR(ErrorUtil::efatal,
                 "This function should be implemented at a shape level.");
        return 0;
    }

    template<class T>
    void Print(T &out)
    {
        int i, j;
        for (i = 0; i < m_vertex.size(); ++i)
        {
            out << m_vertex[i]->m_x << " " << m_vertex[i]->m_y << " "
                << m_vertex[i]->m_z << std::endl;
        }
        for (i = 0; i < m_edge.size(); ++i)
        {
            for (j = 0; j < m_edge[i]->m_edgeNodes.size(); ++j)
            {
                NodeSharedPtr n = m_edge[i]->m_edgeNodes[j];
                out << n->m_x << " " << n->m_y << " " << n->m_z
                    << std::endl;
            }
        }
        for (i = 0; i < m_face.size(); ++i)
        {
            for (j = 0; j < m_face[i]->m_faceNodes.size(); ++j)
            {
                NodeSharedPtr n = m_face[i]->m_faceNodes[j];
                out << n->m_x << " " << n->m_y << " " << n->m_z
                    << std::endl;
            }
        }
    }

    /**
     * @brief returns the normal to the element
     */
    virtual Array<OneD, NekDouble> Normal(bool inward = false)
    {
        boost::ignore_unused(inward);
        NEKERROR(ErrorUtil::efatal,
                 "This function should be implemented at a shape level.");
        return Array<OneD, NekDouble>();
    }

    CADObjectSharedPtr m_parentCAD;

protected:
    /// ID of the element.
    unsigned int                      m_id;
    /// Dimension of the element.
    unsigned int                      m_dim;
    /// Contains configuration of the element.
    ElmtConfig                        m_conf;
    /// Tag character describing the element.
    std::string                       m_tag;
    /// List of integers specifying properties of the element.
    std::vector<int>                  m_taglist;
    /// List of element vertex nodes.
    std::vector<NodeSharedPtr>        m_vertex;
    /// List of element edges.
    std::vector<EdgeSharedPtr>        m_edge;
    /// List of element faces.
    std::vector<FaceSharedPtr>        m_face;
    /// List of element volume nodes.
    std::vector<NodeSharedPtr>        m_volumeNodes;
    /// Volume curve type
    LibUtilities::PointsType          m_curveType;
    /// Pointer to the corresponding edge if element is a 2D boundary.
    EdgeSharedPtr                     m_edgeLink;
    /// Pointer to the corresponding face if element is a 3D boundary.
    FaceSharedPtr                     m_faceLink;
    /// Array mapping faces/edges to the location of the appropriate
    /// boundary elements in m->element.
    std::map<int, int>                m_boundaryLinks;
    /// Nektar++ geometry object for this element.
    SpatialDomains::GeometrySharedPtr m_geom;
};

typedef std::shared_ptr<Element> ElementSharedPtr;
/// Container for elements; key is expansion dimension, value is
/// vector of elements of that dimension.
typedef std::array<std::vector<ElementSharedPtr>, 4> ElementMap;
/// Element factory definition.
typedef LibUtilities::NekFactory<LibUtilities::ShapeType,
                                 Element,
                                 ElmtConfig,
                                 std::vector<NodeSharedPtr>,
                                 std::vector<int> > ElementFactory;

NEKMESH_EXPORT ElementFactory &GetElementFactory();

NEKMESH_EXPORT bool operator==(ElementSharedPtr const &e1,
                               ElementSharedPtr const &e2);

/// Define element ordering based on ID.
struct element_id_less_than
{
    typedef std::shared_ptr<Element> pT;
    bool operator()(const pT a, const pT b) const
    {
        // check for 0
        if (a.get() == 0)
        {
            // if b is also 0, then they are equal, hence a is not
            // less than b
            return b.get() != 0;
        }
        else if (b.get() == 0)
        {
            return false;
        }
        else
        {
            return a->GetId() < b->GetId();
        }
    }
};
}
}

#endif
