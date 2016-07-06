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

#ifndef NekMeshUtils_MESHELEMENTS_TRIANGLE
#define NekMeshUtils_MESHELEMENTS_TRIANGLE

#include <NekMeshUtils/NekMeshUtilsDeclspec.h>
#include <NekMeshUtils/MeshElements/Element.h>

namespace Nektar
{
namespace NekMeshUtils
{

/**
 * @brief A 2-dimensional three-sided element.
 */
class Triangle : public Element
{
public:
    /// Creates an instance of this class
    static ElementSharedPtr create(ElmtConfig pConf,
                                   std::vector<NodeSharedPtr> pNodeList,
                                   std::vector<int> pTagList)
    {
        return boost::shared_ptr<Element>(
            new Triangle(pConf, pNodeList, pTagList));
    }
    /// Element type
    static LibUtilities::ShapeType m_type;

    NEKMESHUTILS_EXPORT Triangle(ElmtConfig pConf,
                                 std::vector<NodeSharedPtr> pNodeList,
                                 std::vector<int> pTagList);
    NEKMESHUTILS_EXPORT Triangle(const Triangle &pSrc);
    NEKMESHUTILS_EXPORT virtual ~Triangle()
    {
    }

    NEKMESHUTILS_EXPORT virtual SpatialDomains::GeometrySharedPtr GetGeom(
        int coordDim);
    NEKMESHUTILS_EXPORT virtual void Complete(int order);
    NEKMESHUTILS_EXPORT virtual void GetCurvedNodes(
        std::vector<NodeSharedPtr> &nodeList) const;
    NEKMESHUTILS_EXPORT static unsigned int GetNumNodes(ElmtConfig pConf);
};

typedef HOTriangle<NodeSharedPtr> HOSurf;
typedef boost::shared_ptr<HOSurf> HOSurfSharedPtr;

/**
 * Hash class for high-order surfaces.
 */
struct HOSurfHash : std::unary_function<HOSurfSharedPtr, std::size_t>
{
    /**
     * Calculate hash of a given high-order surface p by taking
     * successive hashes of the vertex IDs.
     */
    std::size_t operator()(HOSurfSharedPtr const &p) const
    {
        std::size_t seed     = 0;
        std::vector<int> ids = p->vertId;

        std::sort(ids.begin(), ids.end());
        for (int i = 0; i < ids.size(); ++i)
        {
            boost::hash_combine(seed, ids[i]);
        }
        return seed;
    }
};

NEKMESHUTILS_EXPORT bool operator==(HOSurfSharedPtr const &p1,
                                    HOSurfSharedPtr const &p2);

typedef boost::unordered_set<HOSurfSharedPtr, HOSurfHash> HOSurfSet;
}
}

#endif
