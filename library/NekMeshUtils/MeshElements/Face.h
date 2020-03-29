////////////////////////////////////////////////////////////////////////////////
//
//  File: Face.h
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

#ifndef NEKMESHUTILS_MESHELEMENTS_FACE
#define NEKMESHUTILS_MESHELEMENTS_FACE

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/HashUtils.hpp>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/QuadGeom.h>

#include <NekMeshUtils/NekMeshUtilsDeclspec.h>
#include <NekMeshUtils/MeshElements/Edge.h>
#include <NekMeshUtils/MeshElements/Node.h>

namespace Nektar
{
namespace NekMeshUtils
{

class Element;
typedef std::shared_ptr<Element> ElementSharedPtr;

/**
 * @brief Represents a face comprised of three or more edges.
 *
 * A face is defined by a list of vertices, a list of edges joining
 * these vertices, and a list of control nodes within the interior of
 * the face, defining the shape of the face.
 */
class Face
{
public:
    /// Create a new face.
    NEKMESHUTILS_EXPORT Face(std::vector<NodeSharedPtr> pVertexList,
                             std::vector<NodeSharedPtr> pFaceNodes,
                             std::vector<EdgeSharedPtr> pEdgeList,
                              LibUtilities::PointsType pCurveType)
                : m_vertexList(pVertexList), m_edgeList(pEdgeList),
                  m_faceNodes(pFaceNodes), m_curveType(pCurveType), m_geom()
    {
    }

    /// Copy an existing face.
    NEKMESHUTILS_EXPORT Face(const Face &pSrc)
        : m_id(pSrc.m_id), m_vertexList(pSrc.m_vertexList),
          m_edgeList(pSrc.m_edgeList), m_faceNodes(pSrc.m_faceNodes),
          m_curveType(pSrc.m_curveType), m_geom(pSrc.m_geom),
          m_parentCAD(pSrc.m_parentCAD)
    {
    }

    NEKMESHUTILS_EXPORT ~Face()
    {
    }

    /// Equality is defined by matching all vertices.
    NEKMESHUTILS_EXPORT bool operator==(Face &pSrc)
    {
        std::vector<NodeSharedPtr>::iterator it1;
        for (it1 = m_vertexList.begin(); it1 != m_vertexList.end(); ++it1)
        {
            if (find(pSrc.m_vertexList.begin(),
                     pSrc.m_vertexList.end(),
                     *it1) == pSrc.m_vertexList.end())
            {
                return false;
            }
        }
        return true;
    }

    /// Returns the total number of nodes (vertices, edge nodes and
    /// face nodes).
    NEKMESHUTILS_EXPORT unsigned int GetNodeCount() const
    {
        unsigned int n = m_faceNodes.size();
        for (int i = 0; i < m_edgeList.size(); ++i)
        {
            n += m_edgeList[i]->GetNodeCount();
        }
        n -= m_vertexList.size();
        return n;
    }

    /// Assemble a list of nodes on curved face
    NEKMESHUTILS_EXPORT void GetCurvedNodes(
        std::vector<NodeSharedPtr> &nodeList);

    /// Generates a string listing the coordinates of all nodes
    /// associated with this face.
    NEKMESHUTILS_EXPORT std::string GetXmlCurveString();

    /// Make this face an order @p order face. @see Element::MakeOrder.
    void MakeOrder(int                                order,
                   SpatialDomains::GeometrySharedPtr  geom,
                   LibUtilities::PointsType           pType,
                   int                                coordDim,
                   int                               &id);

    /// Generate either SpatialDomains::TriGeom or
    /// SpatialDomains::QuadGeom for this element.
    NEKMESHUTILS_EXPORT SpatialDomains::Geometry2DSharedPtr GetGeom(int coordDim);

    /// ID of the face.
    unsigned int                         m_id;
    /// List of vertex nodes.
    std::vector<NodeSharedPtr>           m_vertexList;
    /// List of corresponding edges.
    std::vector<EdgeSharedPtr>           m_edgeList;
    /// List of face-interior nodes defining the shape of the face.
    std::vector<NodeSharedPtr>           m_faceNodes;
    /// Distribution of points in this face.
    LibUtilities::PointsType             m_curveType;
    /// Element(s) which are linked to this face.
    std::vector<std::pair<std::weak_ptr<Element>, int> > m_elLink;
    /// Nektar++ representation of geometry
    SpatialDomains::Geometry2DSharedPtr  m_geom;

    CADObjectSharedPtr m_parentCAD;
};

typedef std::shared_ptr<Face> FaceSharedPtr;

NEKMESHUTILS_EXPORT bool operator==(FaceSharedPtr const &p1,
                                    FaceSharedPtr const &p2);
NEKMESHUTILS_EXPORT bool operator<(FaceSharedPtr const &p1,
                                   FaceSharedPtr const &p2);

struct FaceHash : std::unary_function<FaceSharedPtr, std::size_t>
{
    std::size_t operator()(FaceSharedPtr const &p) const
    {
        unsigned int nVert = p->m_vertexList.size();
        std::size_t seed = 0;
        std::vector<unsigned int> ids(nVert);

        for (int i = 0; i < nVert; ++i)
        {
            ids[i] = p->m_vertexList[i]->m_id;
        }

        std::sort(ids.begin(), ids.end());
        hash_range(seed, ids.begin(), ids.end());

        return seed;
    }
};
typedef std::unordered_set<FaceSharedPtr, FaceHash> FaceSet;

}
}

#endif
