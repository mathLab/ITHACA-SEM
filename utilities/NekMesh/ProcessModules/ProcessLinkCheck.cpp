////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJac.cpp
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
//  Description: Calculate Jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/MeshElements/Element.h>
#include "ProcessLinkCheck.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessLinkCheck::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "linkcheck"),
    ProcessLinkCheck::create,
    "Checks elemental links within elements.");

ProcessLinkCheck::ProcessLinkCheck(MeshSharedPtr m) : ProcessModule(m)
{

}

ProcessLinkCheck::~ProcessLinkCheck()
{
}

void ProcessLinkCheck::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << "ProcessLinkCheck: Checking links... " << endl;
    }

    //need to reset all links first to make sure there are no bugs!
    ClearElementLinks();
    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();

    int count = 0;

    if(m_mesh->m_expDim == 2)
    {
        for(auto &edge : m_mesh->m_edgeSet)
        {
            if(edge->m_elLink.size() != 2)
            {
                count++;
            }
        }
    }
    else
    {
        for(auto &face : m_mesh->m_faceSet)
        {
            if(face->m_elLink.size() != 2)
            {
                count++;
            }
        }
    }



    if (count - m_mesh->m_element[m_mesh->m_expDim-1].size() != 0)
    {
        cout << "Link Check Error: mesh contains incorrectly connected"
             << " entities and is not valid: "
             << count - m_mesh->m_element[m_mesh->m_expDim-1].size()
             << endl;
    }
}

}
}
