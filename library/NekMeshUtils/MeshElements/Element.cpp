////////////////////////////////////////////////////////////////////////////////
//
//  File: Element.cpp
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
//  Description: Mesh element.
//
////////////////////////////////////////////////////////////////////////////////

#include <loki/Singleton.h>

#include <NekMeshUtils/MeshElements/Element.h>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

ElementFactory &GetElementFactory()
{
    typedef Loki::SingletonHolder<ElementFactory, Loki::CreateUsingNew,
                                  Loki::NoDestroy> Type;
    return Type::Instance();
}

Element::Element(ElmtConfig pConf, unsigned int pNumNodes,
                 unsigned int pGotNodes)
    : m_conf(pConf), m_curveType(LibUtilities::ePolyEvenlySpaced), m_geom()
{
    if (pNumNodes != pGotNodes)
    {
        cerr << "Number of modes mismatch for type " << pConf.m_e
             << "! Should be " << pNumNodes << " but got " << pGotNodes
             << " nodes." << endl;
        abort();
    }
}

/**
 * @brief Replace a vertex in the element.
 *
 * When a vertex is replaced, the element edges and faces are also
 * searched and the corresponding edge/face nodes are updated to
 * maintain consistency.
 *
 * @param  p     Index of the vertex to replace.
 * @param  pNew  New vertex.
 */
void Element::SetVertex(unsigned int p, NodeSharedPtr pNew)
{
    NodeSharedPtr vOld = m_vertex[p];
    m_vertex[p] = pNew;
    for (unsigned int i = 0; i < m_edge.size(); ++i)
    {
        if (m_edge[i]->m_n1 == vOld)
        {
            m_edge[i]->m_n1 = pNew;
        }
        else if (m_edge[i]->m_n2 == vOld)
        {
            m_edge[i]->m_n2 = pNew;
        }
    }
    for (unsigned int i = 0; i < m_face.size(); ++i)
    {
        // Replace vertices in faces
        for (unsigned int j = 0; j < m_face[i]->m_vertexList.size(); ++j)
        {
            if (m_face[i]->m_vertexList[j] == vOld)
            {
                m_face[i]->m_vertexList[j] = pNew;
            }
        }
        for (unsigned int j = 0; j < m_face[i]->m_edgeList.size(); ++j)
        {
            if (m_face[i]->m_edgeList[j]->m_n1 == vOld)
            {
                m_face[i]->m_edgeList[j]->m_n1 = pNew;
            }
            else if (m_face[i]->m_edgeList[j]->m_n2 == vOld)
            {
                m_face[i]->m_edgeList[j]->m_n2 = pNew;
            }
        }
    }
}

/**
 * @brief Replace an edge in the element.
 *
 * When an edge is replaced, the element faces are also searched and
 * the corresponding face edges are updated to maintain consistency.
 *
 * @param  p     Index of the edge to replace.
 * @param  pNew  New edge.
 */
void Element::SetEdge(unsigned int p, EdgeSharedPtr pNew)
{
    EdgeSharedPtr vOld = m_edge[p];
    m_edge[p] = pNew;
    for (unsigned int i = 0; i < m_face.size(); ++i)
    {
        for (unsigned int j = 0; j < m_face[i]->m_edgeList.size(); ++j)
        {
            if (m_face[i]->m_edgeList[j] == vOld)
            {
                m_face[i]->m_edgeList[j] = pNew;
            }
        }
    }
}

/**
 * @brief Replace a face in the element.
 *
 * When a face is replaced, no other consistency checks are required.
 *
 * @param  p     Index of the face to replace.
 * @param  pNew  New face.
 */
void Element::SetFace(unsigned int p, FaceSharedPtr pNew)
{
    m_face[p] = pNew;
}

/**
 * @brief Obtain the order of an element by looking at edges.
 */
int Element::GetMaxOrder()
{
    int i, ret = 1;

    for (i = 0; i < m_edge.size(); ++i)
    {
        int edgeOrder = m_edge[i]->GetNodeCount() - 1;
        if (edgeOrder > ret)
        {
            ret = edgeOrder;
        }
    }

    return ret;
}
}
}
