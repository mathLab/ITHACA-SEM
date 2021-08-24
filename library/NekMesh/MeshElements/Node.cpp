////////////////////////////////////////////////////////////////////////////////
//
//  File: Node.cpp
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
//  Description: Mesh node object.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMesh/MeshElements/Node.h>

namespace Nektar
{
namespace NekMesh
{

/// Define node equality based on coordinate with optional custom tolerance factor.
NEKMESH_EXPORT bool IsNodeEqual(const Node &n1, const Node &n2,
    const unsigned int fact)
{
    return LibUtilities::IsRealEqual(n1.m_x, n2.m_x, fact) &&
           LibUtilities::IsRealEqual(n1.m_y, n2.m_y, fact) &&
           LibUtilities::IsRealEqual(n1.m_z, n2.m_z, fact);
}

/// Define node equality based on coordinate with optional custom tolerance factor.
NEKMESH_EXPORT bool IsNodeClose(const Node &n1, const Node &n2,
    const NekDouble tol)
{
    return LibUtilities::IsRealClose(n1.m_x, n2.m_x, tol) &&
           LibUtilities::IsRealClose(n1.m_y, n2.m_y, tol) &&
           LibUtilities::IsRealClose(n1.m_z, n2.m_z, tol);
}

}
}
