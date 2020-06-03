////////////////////////////////////////////////////////////////////////////////
//
//  File: Hexahedron.h
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
//  Description: Mesh Hexahedral object.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NekMesh_MESHELEMENTS_HEX
#define NekMesh_MESHELEMENTS_HEX

#include <NekMesh/NekMeshDeclspec.h>
#include <NekMesh/MeshElements/Element.h>

namespace Nektar
{
namespace NekMesh
{
/**
 * @brief A 3-dimensional six-faced element.
 */
class Hexahedron : public Element
{
public:
    /// Creates an instance of this class
    NEKMESH_EXPORT static ElementSharedPtr create(
        ElmtConfig                 pConf,
        std::vector<NodeSharedPtr> pNodeList,
        std::vector<int>           pTagList)
    {
        return std::shared_ptr<Element>(
            new Hexahedron(pConf, pNodeList, pTagList));
    }
    /// Element type
    static LibUtilities::ShapeType m_type;

    /**
     * @brief Create a hexahedral element.
     */
    NEKMESH_EXPORT Hexahedron(ElmtConfig pConf,
                                   std::vector<NodeSharedPtr> pNodeList,
                                   std::vector<int> pTagList);
    /// copy hex element
    NEKMESH_EXPORT Hexahedron(const Hexahedron &pSrc);
    NEKMESH_EXPORT virtual ~Hexahedron()
    {
    }

    NEKMESH_EXPORT virtual SpatialDomains::GeometrySharedPtr GetGeom(
        int coordDim);
    NEKMESH_EXPORT virtual StdRegions::Orientation GetEdgeOrient(
        int edgeId, EdgeSharedPtr edge);
    NEKMESH_EXPORT virtual void MakeOrder(
        int                                order,
        SpatialDomains::GeometrySharedPtr  geom,
        LibUtilities::PointsType           pType,
        int                                coordDim,
        int                               &id,
        bool                               justConfig = false);

    NEKMESH_EXPORT static unsigned int GetNumNodes(ElmtConfig pConf);
    NEKMESH_EXPORT virtual int GetFaceVertex(int i, int j)
    {
        return m_faceIds[i][j];
    }

private:
    static int m_faceIds[6][4];
};
}
}

#endif
