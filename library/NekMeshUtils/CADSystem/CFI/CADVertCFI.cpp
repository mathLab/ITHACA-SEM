////////////////////////////////////////////////////////////////////////////////
//
//  File: CADVertCFI.cpp
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
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include "CADVertCFI.h"

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

std::string CADVertCFI::key = GetCADVertFactory().RegisterCreatorFunction(
        "cfi", CADVertCFI::create, "CAD vert cfi");

void CADVertCFI::Initialise(int i, cfi::Point* in, NekDouble s)
{
    m_id      = i;
    m_cfipoint = in;
    m_scal = s;

    cfi::Position pos = m_cfipoint->getGeometry();

    m_node = std::shared_ptr<Node>(
        new Node(i - 1, pos.x*m_scal, pos.y*m_scal, pos.z*m_scal));
    degen = false;
}

}
}
