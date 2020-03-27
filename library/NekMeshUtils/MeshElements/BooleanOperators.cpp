////////////////////////////////////////////////////////////////////////////////
//
//  File: BooleanOperators.cpp
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
//  Description: Boolean operators for comparison of mesh objects.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/MeshElements/Triangle.h>
#include <NekMeshUtils/MeshElements/Mesh.h>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

/**
 * @brief Compares two element config structs
 */
bool operator==(ElmtConfig const &c1, ElmtConfig const &c2)
{
    return (c1.m_e == c2.m_e && c1.m_order == c2.m_order);
}

/**
 * @brief Compares two element shared pointers
 */
bool operator==(ElementSharedPtr const &e1, ElementSharedPtr const &e2)
{
    return e1->GetId() == e2->GetId();
}

/**
 * @brief Compares two %HOSurf objects referred to as shared pointers.
 *
 * Two %HOSurf objects are defined to be equal if they contain identical
 * vertex ids contained in HOSurf::vertId.
 */
bool operator==(HOSurfSharedPtr const &p1, HOSurfSharedPtr const &p2)
{
    if (p1->vertId.size() != p2->vertId.size())
    {
        return false;
    }

    vector<int> ids1 = p1->vertId;
    vector<int> ids2 = p2->vertId;
    sort(ids1.begin(), ids1.end());
    sort(ids2.begin(), ids2.end());

    for (int i = 0; i < ids1.size(); ++i)
    {
        if (ids1[i] != ids2[i])
            return false;
    }

    return true;
}

/**
 * @brief Test equality of two conditions - i.e. compare types, fields
 * and values but _not_ composite ids.
 */
bool operator==(ConditionSharedPtr const &c1, ConditionSharedPtr const &c2)
{
    int i, n = c1->type.size();

    if (n != c2->type.size())
    {
        return false;
    }

    for (i = 0; i < n; ++i)
    {
        if (c1->type[i] != c2->type[i])
        {
            return false;
        }

        if (c1->field[i] != c2->field[i] || c1->value[i] != c2->value[i])
        {
            return false;
        }
    }

    return true;
}

/**
 * @brief Defines equality between two #NodeSharedPtr objects.
 */
bool operator==(NodeSharedPtr const &p1, NodeSharedPtr const &p2)
{
    return *p1 == *p2;
}

/**
 * @brief Compares two nodes for inequality based on IDs
 */
bool operator!=(NodeSharedPtr const &p1, NodeSharedPtr const &p2)
{
    if (p1->m_id != p2->m_id)
    {
        return true;
    }
    else
    {
        return false;
    }
}

/**
 * @brief Defines ordering between two #NodeSharedPtr objects.
 */
bool operator< (NodeSharedPtr const &p1, NodeSharedPtr const &p2)
{
    return *p1 < *p2;
}

/**
 * @brief Print description of node to given ostream.
 */
std::ostream &operator<<(std::ostream &os, const NodeSharedPtr &n)
{
    os << n->m_x << " " << n->m_y << " " << n->m_z;
    return os;
}

/**
 * @brief Defines equality of two edges (equal if IDs of end nodes
 * match in either ordering).
 */
bool operator==(EdgeSharedPtr const &p1, EdgeSharedPtr const &p2)
{
    return ( ((*(p1->m_n1) == *(p2->m_n1)) && (*(p1->m_n2) == *(p2->m_n2)))
          || ((*(p1->m_n2) == *(p2->m_n1)) && (*(p1->m_n1) == *(p2->m_n2))));
}

/**
 * @brief Defines ordering between two edges (based on ID of edges).
 */
bool operator< (EdgeSharedPtr const &p1, EdgeSharedPtr const &p2)
{
    return p1->m_id < p2->m_id;
}

/**
 * @brief Defines equality of two faces (equal if IDs of vertices are
 * the same.)
 */
bool operator==(FaceSharedPtr const &p1, FaceSharedPtr const &p2)
{
    std::vector<NodeSharedPtr>::iterator it1;
    for (it1 = p1->m_vertexList.begin(); it1 != p1->m_vertexList.end(); ++it1)
    {
        if (find(p2->m_vertexList.begin(), p2->m_vertexList.end(), *it1)
            == p2->m_vertexList.end())
        {
            return false;
        }
    }
    return true;
}

/**
 * @brief Defines ordering between two faces (depending on ID of
 * faces).
 */
bool operator< (FaceSharedPtr const &p1, FaceSharedPtr const &p2)
{
    return p1->m_id < p2->m_id;
}

}
}
