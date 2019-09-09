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

#include <NekMeshUtils/MeshElements/Element.h>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

ElementFactory &GetElementFactory()
{
    static ElementFactory instance;
    return instance;
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

void Element::SetVertex(unsigned int p, NodeSharedPtr pNew, bool descend)
{
    NodeSharedPtr vOld = m_vertex[p];
    m_vertex[p] = pNew;

    if (!descend)
    {
        return;
    }

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

void Element::SetEdge(unsigned int p, EdgeSharedPtr pNew, bool descend)
{
    EdgeSharedPtr vOld = m_edge[p];
    m_edge[p] = pNew;

    if (!descend)
    {
        return;
    }

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

void Element::SetFace(unsigned int p, FaceSharedPtr pNew)
{
    m_face[p] = pNew;
}

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

unsigned int Element::GetNodeCount()
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
        for (int i = 0; i < m_face.size(); ++i)
        {
            n += m_face[i]->GetNodeCount();
        }
        for (int i = 0; i < m_edge.size(); ++i)
        {
            n -= m_edge[i]->GetNodeCount();
        }
        n += m_vertex.size();
        std::cerr << "Not supported." << std::endl;
        exit(1);
    }
    return n;
}

string Element::GetXmlString()
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

string Element::GetXmlCurveString()
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

}
}
