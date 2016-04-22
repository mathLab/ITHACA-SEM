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

#ifndef NekMeshUtils_MESHELEMENTS_ELEMENT
#define NekMeshUtils_MESHELEMENTS_ELEMENT

#include <iomanip>

#include <LibUtilities/Foundations/PointsType.h>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>

#include <NekMeshUtils/NekMeshUtilsDeclspec.h>
#include <NekMeshUtils/MeshElements/Face.h>

namespace Nektar
{
namespace NekMeshUtils
{
/**
 * @brief Basic information about an element.
 *
 * ElmtConfig contains four member variables which denote the
 * properties of an element when it is created.
 */
struct ElmtConfig
{
    ElmtConfig(LibUtilities::ShapeType  pE,
               unsigned int             pOrder,
               bool                     pFn,
               bool                     pVn,
               bool                     pReorient = true,
               LibUtilities::PointsType pECt = LibUtilities::ePolyEvenlySpaced,
               LibUtilities::PointsType pFCt = LibUtilities::ePolyEvenlySpaced)
        : m_e(pE), m_faceNodes(pFn), m_volumeNodes(pVn), m_order(pOrder),
          m_reorient(pReorient), m_edgeCurveType(pECt), m_faceCurveType(pFCt)
    {
    }

    ElmtConfig(ElmtConfig const &p)
        : m_e(p.m_e), m_faceNodes(p.m_faceNodes),
          m_volumeNodes(p.m_volumeNodes), m_order(p.m_order),
          m_reorient(p.m_reorient), m_edgeCurveType(p.m_edgeCurveType),
          m_faceCurveType(p.m_faceCurveType)
    {
    }

    ElmtConfig()
    {
    }

    /// Element type (e.g. triangle, quad, etc).
    LibUtilities::ShapeType m_e;
    /// Denotes whether the element contains face nodes. For 2D
    /// elements, if this is true then the element contains interior
    /// nodes.
    bool m_faceNodes;
    /// Denotes whether the element contains volume (i.e. interior)
    /// nodes. These are not supported by either the mesh converter or
    /// Nektar++ but are included for completeness and are required
    /// for some output modules (e.g. Gmsh).
    bool m_volumeNodes;
    /// Order of the element.
    unsigned int m_order;
    /// Denotes whether the element needs to be re-orientated for a
    /// spectral element framework.
    bool m_reorient;
    /// Distribution of points in edges.
    LibUtilities::PointsType m_edgeCurveType;
    /// Distribution of points in faces.
    LibUtilities::PointsType m_faceCurveType;
};

NEKMESHUTILS_EXPORT bool operator==(ElmtConfig const &c1, ElmtConfig const &c2);

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
    NEKMESHUTILS_EXPORT Element(
        ElmtConfig pConf, unsigned int pNumNodes, unsigned int pGotNodes);

    /// Returns the ID of the element (or associated edge or face for
    /// boundary elements).
    NEKMESHUTILS_EXPORT unsigned int GetId() const
    {
        if (m_faceLink.get() != 0)
            return m_faceLink->m_id;
        if (m_edgeLink.get() != 0)
            return m_edgeLink->m_id;
        return m_id;
    }
    /// Returns the expansion dimension of the element.
    NEKMESHUTILS_EXPORT unsigned int GetDim() const
    {
        return m_dim;
    }
    /// Returns the configuration of the element.
    NEKMESHUTILS_EXPORT ElmtConfig GetConf() const
    {
        return m_conf;
    }
    /// Returns the tag which defines the element shape.
    NEKMESHUTILS_EXPORT std::string GetTag() const
    {
        if (m_faceLink.get() != 0)
            return "F";
        if (m_edgeLink.get() != 0)
            return "E";
        return m_tag;
    }
    /// Access a vertex node.
    NEKMESHUTILS_EXPORT NodeSharedPtr GetVertex(unsigned int i) const
    {
        return m_vertex[i];
    }
    /// Access an edge.
    NEKMESHUTILS_EXPORT EdgeSharedPtr GetEdge(unsigned int i) const
    {
        return m_edge[i];
    }
    /// Access a face.
    NEKMESHUTILS_EXPORT FaceSharedPtr GetFace(unsigned int i) const
    {
        return m_face[i];
    }
    /// Access the list of vertex nodes.
    NEKMESHUTILS_EXPORT std::vector<NodeSharedPtr> GetVertexList() const
    {
        return m_vertex;
    }
    /// Access the list of edges.
    NEKMESHUTILS_EXPORT std::vector<EdgeSharedPtr> GetEdgeList() const
    {
        return m_edge;
    }
    /// Access the list of faces.
    NEKMESHUTILS_EXPORT std::vector<FaceSharedPtr> GetFaceList() const
    {
        return m_face;
    }
    /// Access the list of volume nodes.
    NEKMESHUTILS_EXPORT std::vector<NodeSharedPtr> GetVolumeNodes() const
    {
        return m_volumeNodes;
    }
    NEKMESHUTILS_EXPORT void SetVolumeNodes(std::vector<NodeSharedPtr> &nodes)
    {
        m_volumeNodes = nodes;
    }
    NEKMESHUTILS_EXPORT LibUtilities::PointsType GetCurveType() const
    {
        return m_curveType;
    }
    NEKMESHUTILS_EXPORT void SetCurveType(LibUtilities::PointsType cT)
    {
        m_curveType = cT;
    }
    /// Returns the total number of nodes (vertices, edge nodes and
    /// face nodes and volume nodes).
    NEKMESHUTILS_EXPORT unsigned int GetNodeCount() const
    {
        unsigned int n = m_volumeNodes.size();
        if (m_dim == 1)
        {
            n += 2;
        }
        else if (m_dim == 2)
        {
            for (int i = 0; i < m_edge.size(); ++i)
            {
                n += m_edge[i]->GetNodeCount();
            }
            n -= m_vertex.size();
        }
        else
        {
            std::cerr << "Not supported." << std::endl;
            exit(1);
        }
        return n;
    }
    /// Access the list of tags associated with this element.
    NEKMESHUTILS_EXPORT std::vector<int> GetTagList() const
    {
        return m_taglist;
    }
    /// Returns the number of vertices.
    NEKMESHUTILS_EXPORT unsigned int GetVertexCount() const
    {
        return m_vertex.size();
    }
    /// Returns the number of edges.
    NEKMESHUTILS_EXPORT unsigned int GetEdgeCount() const
    {
        return m_edge.size();
    }
    /// Returns the number of faces.
    NEKMESHUTILS_EXPORT unsigned int GetFaceCount() const
    {
        return m_face.size();
    }
    /// Change the ID of the element.
    NEKMESHUTILS_EXPORT void SetId(unsigned int p)
    {
        m_id = p;
    }
    /// Replace a vertex with another vertex object.
    NEKMESHUTILS_EXPORT void SetVertex(unsigned int p, NodeSharedPtr pNew);
    /// Replace an edge with another edge object.
    NEKMESHUTILS_EXPORT void SetEdge(unsigned int p, EdgeSharedPtr pNew);
    /// Replace a face with another face object.
    NEKMESHUTILS_EXPORT void SetFace(unsigned int p, FaceSharedPtr pNew);
    /// Set a correspondence between this element and an edge
    /// (2D boundary element).
    NEKMESHUTILS_EXPORT void SetEdgeLink(EdgeSharedPtr pLink)
    {
        m_edgeLink = pLink;
    }
    /// Get correspondence between this element and an edge.
    NEKMESHUTILS_EXPORT EdgeSharedPtr GetEdgeLink()
    {
        return m_edgeLink;
    }
    /// Set a correspondence between this element and a face
    /// (3D boundary element).
    NEKMESHUTILS_EXPORT void SetFaceLink(FaceSharedPtr pLink)
    {
        m_faceLink = pLink;
    }
    /// Get correspondence between this element and a face.
    NEKMESHUTILS_EXPORT FaceSharedPtr GetFaceLink()
    {
        return m_faceLink;
    }
    /// Set a correspondence between edge or face i and its
    /// representative boundary element m->element[expDim-1][j].
    NEKMESHUTILS_EXPORT void SetBoundaryLink(int i, int j)
    {
        m_boundaryLinks[i] = j;
    }
    /// Get the location of the boundary face/edge i for this element.
    NEKMESHUTILS_EXPORT int GetBoundaryLink(int i)
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
    /// Set the list of tags associated with this element.
    NEKMESHUTILS_EXPORT void SetTagList(const std::vector<int> &tags)
    {
        m_taglist = tags;
    }
    /// Generate a list of vertices (1D), edges (2D), or faces (3D).
    NEKMESHUTILS_EXPORT virtual std::string GetXmlString() const
    {
        std::stringstream s;
        switch (m_dim)
        {
            case 1:
                for (int j = 0; j < m_vertex.size(); ++j)
                {
                    s << std::setw(5) << m_vertex[j]->m_id << " ";
                }
                break;
            case 2:
                for (int j = 0; j < m_edge.size(); ++j)
                {
                    s << std::setw(5) << m_edge[j]->m_id << " ";
                }
                break;
            case 3:
                for (int j = 0; j < m_face.size(); ++j)
                {
                    s << std::setw(5) << m_face[j]->m_id << " ";
                }
                break;
        }
        return s.str();
    }

    NEKMESHUTILS_EXPORT void GetCurvedNodes(
        std::vector<NodeSharedPtr> &nodeList) const
    {
        // Node orderings are different for different elements.
        // Triangle
        if (m_vertex.size() == 2)
        {
            nodeList.push_back(m_vertex[0]);
            for (int i = 0; i < m_volumeNodes.size(); ++i)
            {
                nodeList.push_back(m_volumeNodes[i]);
            }
            nodeList.push_back(m_vertex[1]);
        }
        else if (m_vertex.size() == 3)
        {
            int n = m_edge[0]->GetNodeCount();
            nodeList.resize(n * (n + 1) / 2);

            // Populate nodelist
            std::copy(m_vertex.begin(), m_vertex.end(), nodeList.begin());
            for (int i = 0; i < 3; ++i)
            {
                std::copy(m_edge[i]->m_edgeNodes.begin(),
                          m_edge[i]->m_edgeNodes.end(),
                          nodeList.begin() + 3 + i * (n - 2));
                if (m_edge[i]->m_n1 != m_vertex[i])
                {
                    // If edge orientation is reversed relative to node
                    // ordering, we need to reverse order of nodes.
                    std::reverse(nodeList.begin() + 3 + i * (n - 2),
                                 nodeList.begin() + 3 + (i + 1) * (n - 2));
                }
            }

            // Copy volume nodes.
            std::copy(m_volumeNodes.begin(),
                      m_volumeNodes.end(),
                      nodeList.begin() + 3 * (n - 1));
        }
        // Quad
        else if (m_dim == 2 && m_vertex.size() == 4)
        {
            int n = m_edge[0]->GetNodeCount();
            nodeList.resize(n * n);

            // Write vertices
            nodeList[0]         = m_vertex[0];
            nodeList[n - 1]     = m_vertex[1];
            nodeList[n * n - 1] = m_vertex[2];
            nodeList[n * (n - 1)] = m_vertex[3];

            // Write edge-interior
            int skips[4][2] = {
                {0, 1}, {n - 1, n}, {n * n - 1, -1}, {n * (n - 1), -n}};
            for (int i = 0; i < 4; ++i)
            {
                bool reverseEdge = m_edge[i]->m_n1 == m_vertex[i];

                if (!reverseEdge)
                {
                    for (int j = 1; j < n - 1; ++j)
                    {
                        nodeList[skips[i][0] + j * skips[i][1]] =
                            m_edge[i]->m_edgeNodes[n - 2 - j];
                    }
                }
                else
                {
                    for (int j = 1; j < n - 1; ++j)
                    {
                        nodeList[skips[i][0] + j * skips[i][1]] =
                            m_edge[i]->m_edgeNodes[j - 1];
                    }
                }
            }

            // Write interior
            for (int i = 1; i < n - 1; ++i)
            {
                for (int j = 1; j < n - 1; ++j)
                {
                    nodeList[i * n + j] =
                        m_volumeNodes[(i - 1) * (n - 2) + (j - 1)];
                }
            }
        }
        else
        {
            std::cerr << "GetXmlCurveString for a " << m_vertex.size()
                      << "-vertex element is not yet implemented."
                      << std::endl;
        }
    }

    /// Generates a string listing the coordinates of all nodes
    /// associated with this element.
    NEKMESHUTILS_EXPORT std::string GetXmlCurveString() const
    {
        // Temporary node list for reordering
        std::vector<NodeSharedPtr> nodeList;

        GetCurvedNodes(nodeList);

        // Finally generate the XML string corresponding to our new
        // node reordering.
        std::stringstream s;
        std::string str;
        for (int k = 0; k < nodeList.size(); ++k)
        {
            s << std::scientific << std::setprecision(8) << "    "
              << nodeList[k]->m_x << "  " << nodeList[k]->m_y << "  "
              << nodeList[k]->m_z << "    ";
        }
        return s.str();
    }

    /// Generate a Nektar++ geometry object for this element.
    NEKMESHUTILS_EXPORT virtual SpatialDomains::GeometrySharedPtr GetGeom(
        int coordDim)
    {
        ASSERTL0(false,
                 "This function should be implemented at a shape level.");
        return boost::shared_ptr<SpatialDomains::Geometry>();
    }
    NEKMESHUTILS_EXPORT int GetMaxOrder();
    /// Complete this object.
    NEKMESHUTILS_EXPORT virtual void Complete(int order)
    {
        ASSERTL0(false,
                 "This function should be implemented at a shape level.");
    }
    NEKMESHUTILS_EXPORT void Print()
    {
        int i, j;
        for (i = 0; i < m_vertex.size(); ++i)
        {
            std::cout << m_vertex[i]->m_x << " " << m_vertex[i]->m_y << " "
                      << m_vertex[i]->m_z << std::endl;
        }
        for (i = 0; i < m_edge.size(); ++i)
        {
            for (j = 0; j < m_edge[i]->m_edgeNodes.size(); ++j)
            {
                NodeSharedPtr n = m_edge[i]->m_edgeNodes[j];
                std::cout << n->m_x << " " << n->m_y << " " << n->m_z
                          << std::endl;
            }
        }
        for (i = 0; i < m_face.size(); ++i)
        {
            for (j = 0; j < m_face[i]->m_faceNodes.size(); ++j)
            {
                NodeSharedPtr n = m_face[i]->m_faceNodes[j];
                std::cout << n->m_x << " " << n->m_y << " " << n->m_z
                          << std::endl;
            }
        }
    }

#ifdef NEKTAR_USE_MESHGEN
    int CADSurfId;
#endif

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

typedef boost::shared_ptr<Element> ElementSharedPtr;
/// Container for elements; key is expansion dimension, value is
/// vector of elements of that dimension.
typedef std::map<unsigned int, std::vector<ElementSharedPtr> > ElementMap;
/// Element factory definition.
typedef Nektar::LibUtilities::NekFactory<LibUtilities::ShapeType,
                                         Element,
                                         ElmtConfig,
                                         std::vector<NodeSharedPtr>,
                                         std::vector<int> > ElementFactory;

NEKMESHUTILS_EXPORT ElementFactory &GetElementFactory();

/// Define element ordering based on ID.
struct element_id_less_than
{
    typedef boost::shared_ptr<Element> pT;
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
