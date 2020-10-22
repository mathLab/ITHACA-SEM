////////////////////////////////////////////////////////////////////////////////
//
//  File: Point.h
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
//  Description: Mesh point object.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESH_MESHELEMENTS_POINT
#define NEKMESH_MESHELEMENTS_POINT

#include <NekMesh/NekMeshDeclspec.h>
#include <NekMesh/MeshElements/Element.h>

namespace Nektar
{
namespace NekMesh
{
/**
 * @brief A 0-dimensional vertex.
 */
class Point : public Element
{
public:
    /// Creates an instance of this class
    static ElementSharedPtr create(ElmtConfig pConf,
                                   std::vector<NodeSharedPtr> pNodeList,
                                   std::vector<int> pTagList)
    {
        return std::make_shared<Point>(pConf, pNodeList, pTagList);
    }
    /// Element type
    static LibUtilities::ShapeType m_type;

    NEKMESH_EXPORT Point(ElmtConfig pConf,
                              std::vector<NodeSharedPtr> pNodeList,
                              std::vector<int> pTagList);
    NEKMESH_EXPORT Point(const Point &pSrc);
    NEKMESH_EXPORT virtual ~Point()
    {
    }

    NEKMESH_EXPORT static unsigned int GetNumNodes(ElmtConfig pConf);
};
}
}

#endif
