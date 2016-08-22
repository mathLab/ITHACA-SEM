////////////////////////////////////////////////////////////////////////////////
//
//  File: Mesh.h
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

#ifndef NekMeshUtils_MESHELEMENTS_QUAD
#define NekMeshUtils_MESHELEMENTS_QUAD

#include <NekMeshUtils/NekMeshUtilsDeclspec.h>
#include <NekMeshUtils/MeshElements/Element.h>

namespace Nektar
{
namespace NekMeshUtils
{

/**
 * @brief A 2-dimensional four-sided element.
 */
class Quadrilateral : public Element
{
public:
    /// Creates an instance of this class
    static ElementSharedPtr create(ElmtConfig pConf,
                                   std::vector<NodeSharedPtr> pNodeList,
                                   std::vector<int> pTagList)
    {
        ElementSharedPtr e = boost::shared_ptr<Element>(
            new Quadrilateral(pConf, pNodeList, pTagList));
        std::vector<EdgeSharedPtr> m_edges = e->GetEdgeList();
        for (int i = 0; i < m_edges.size(); ++i)
        {
            m_edges[i]->m_elLink.push_back(std::pair<ElementSharedPtr, int>(e, i));
        }
        return e;
    }
    /// Element type
    static LibUtilities::ShapeType m_type;

    NEKMESHUTILS_EXPORT Quadrilateral(ElmtConfig pConf,
                                      std::vector<NodeSharedPtr> pNodeList,
                                      std::vector<int> pTagList);
    NEKMESHUTILS_EXPORT Quadrilateral(const Quadrilateral &pSrc);
    NEKMESHUTILS_EXPORT virtual ~Quadrilateral()
    {
    }

    NEKMESHUTILS_EXPORT virtual SpatialDomains::GeometrySharedPtr GetGeom(
        int coordDim);
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
};

/**
 * @brief A lightweight struct for dealing with high-order quadrilateral
 * alignment.
 */
template <typename T> struct HOQuadrilateral
{
    HOQuadrilateral(std::vector<int> pVertId, std::vector<T> pSurfVerts)
        : vertId(pVertId), surfVerts(pSurfVerts)
    {
    }

    HOQuadrilateral(std::vector<int> pVertId) : vertId(pVertId)
    {
    }

    /// The quadrilateral vertex IDs
    std::vector<int> vertId;

    /// The quadrilateral surface vertices -- templated so that this can either
    /// be nodes or IDs.
    std::vector<T> surfVerts;

    void ReverseX()
    {
        int np = (int)(sqrt((NekDouble)surfVerts.size()) + 0.5);
        for (int i = 0; i < np; ++i)
        {
            for (int j = 0; j < np/2; ++j)
            {
                swap(surfVerts[i*np + j], surfVerts[i*np + np-j-1]);
            }
        }
    }

    void ReverseY()
    {
        int np = (int)(sqrt((NekDouble)surfVerts.size()) + 0.5);
        // Reverse y direction
        for (int j = 0; j < np; ++j)
        {
            for (int i = 0; i < np/2; ++i)
            {
                swap(surfVerts[i*np + j], surfVerts[(np-i-1)*np + j]);
            }
        }
    }

    void Transpose()
    {
        int np = (int)(sqrt((NekDouble)surfVerts.size()) + 0.5);
        std::vector<T> tmp(surfVerts.size());

        for (int i = 0; i < np; ++i)
        {
            for (int j = 0; j < np; ++j)
            {
                tmp[i*np+j] = surfVerts[j*np+i];
            }
        }

        surfVerts = tmp;
    }

    /**
     * @brief Align this surface to a given vertex ID.
     */
    void Align(std::vector<int> vertId)
    {
        int vmap[4] = {-1, -1, -1, -1};

        // Determine which vertices map to vertId
        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 4; ++j)
            {
                if (this->vertId[j] == vertId[i])
                {
                    vmap[i] = j;
                    break;
                }
            }

            ASSERTL1(vmap[i] != -1,
                     "Could not determine mapping between vertex IDs");
        }

        StdRegions::Orientation orient = StdRegions::eNoOrientation;

        if (vmap[1] == (vmap[0]+1) % 4)
        {
            switch (vmap[0])
            {
                case 0:
                    orient = StdRegions::eDir1FwdDir1_Dir2FwdDir2;
                    break;
                case 1:
                    orient = StdRegions::eDir1BwdDir2_Dir2FwdDir1;
                    break;
                case 2:
                    orient = StdRegions::eDir1BwdDir1_Dir2BwdDir2;
                    break;
                case 3:
                    orient = StdRegions::eDir1FwdDir2_Dir2BwdDir1;
                    break;
            }
        }
        else
        {
            switch (vmap[0])
            {
                case 0:
                    orient = StdRegions::eDir1FwdDir2_Dir2FwdDir1;
                    break;
                case 1:
                    orient = StdRegions::eDir1BwdDir1_Dir2FwdDir2;
                    break;
                case 2:
                    orient = StdRegions::eDir1BwdDir2_Dir2BwdDir1;
                    break;
                case 3:
                    orient = StdRegions::eDir1FwdDir1_Dir2BwdDir2;
                    break;
            }
        }

        if (orient == StdRegions::eDir1FwdDir2_Dir2FwdDir1 ||
            orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1 ||
            orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1 ||
            orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
        {
            Transpose();
        }

        if (orient == StdRegions::eDir1BwdDir1_Dir2FwdDir2 ||
            orient == StdRegions::eDir1BwdDir1_Dir2BwdDir2 ||
            orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1 ||
            orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
        {
            ReverseX();
        }

        if (orient == StdRegions::eDir1FwdDir1_Dir2BwdDir2 ||
            orient == StdRegions::eDir1BwdDir1_Dir2BwdDir2 ||
            orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1 ||
            orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
        {
            ReverseY();
        }
    }
};

}
}

#endif
