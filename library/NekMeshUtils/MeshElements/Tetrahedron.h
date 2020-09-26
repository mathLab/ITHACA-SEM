////////////////////////////////////////////////////////////////////////////////
//
//  File: Tetrahedron.h
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
//  Description: Mesh tet object.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NekMeshUtils_MESHELEMENTS_TET
#define NekMeshUtils_MESHELEMENTS_TET

#include <NekMeshUtils/NekMeshUtilsDeclspec.h>
#include <NekMeshUtils/MeshElements/Element.h>

namespace Nektar
{
namespace NekMeshUtils
{

/**
 * @brief A 3-dimensional four-faced element.
 */
class Tetrahedron : public Element
{
public:
    /// Creates an instance of this class
    static ElementSharedPtr create(ElmtConfig pConf,
                                   std::vector<NodeSharedPtr> pNodeList,
                                   std::vector<int> pTagList)
    {
        return std::make_shared<Tetrahedron>(pConf, pNodeList, pTagList);
    }
    /// Element type
    static LibUtilities::ShapeType m_type;

    NEKMESHUTILS_EXPORT Tetrahedron(ElmtConfig pConf,
                                    std::vector<NodeSharedPtr> pNodeList,
                                    std::vector<int> pTagList);
    NEKMESHUTILS_EXPORT Tetrahedron(const Tetrahedron &pSrc);
    NEKMESHUTILS_EXPORT virtual ~Tetrahedron()
    {
    }

    NEKMESHUTILS_EXPORT virtual SpatialDomains::GeometrySharedPtr GetGeom(
        int coordDim);
    NEKMESHUTILS_EXPORT virtual void GetCurvedNodes(
        std::vector<NodeSharedPtr> &nodeList) const;
    NEKMESHUTILS_EXPORT virtual StdRegions::Orientation GetEdgeOrient(
        int edgeId, EdgeSharedPtr edge);
    NEKMESHUTILS_EXPORT virtual void MakeOrder(
        int                                order,
        SpatialDomains::GeometrySharedPtr  geom,
        LibUtilities::PointsType           pType,
        int                                coordDim,
        int                               &id,
        bool                               justConfig = false);

    NEKMESHUTILS_EXPORT static unsigned int GetNumNodes(ElmtConfig pConf);
    NEKMESHUTILS_EXPORT virtual int GetFaceVertex(int i, int j)
    {
        return m_faceVertMap[i][j];
    }

    int m_orientationMap[4];
    int m_origVertMap[4];

protected:
    void OrientTet();

private:
    static int m_edgeVertMap[6][2];
    static int m_faceVertMap[4][3];
    static int m_faceEdgeMap[4][3];
};
}
}

#endif
