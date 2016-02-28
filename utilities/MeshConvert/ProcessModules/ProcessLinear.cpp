////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessLinear.cpp
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
//  Description: linearises mesh.
//
////////////////////////////////////////////////////////////////////////////////

#include "../MeshElements.h"
#include "ProcessLinear.h"

#include <vector>
using namespace std;

namespace Nektar
{
namespace Utilities
{
ModuleKey ProcessLinear::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "linearise"), ProcessLinear::create,
        "Linearises mesh.");

ProcessLinear::ProcessLinear(MeshSharedPtr m) : ProcessModule(m)
{

}

ProcessLinear::~ProcessLinear()
{

}

void ProcessLinear::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << "ProcessLinear: Linearising mesh... " << endl;
    }

    EdgeSet::iterator eit;
    for(eit = m_mesh->m_edgeSet.begin(); eit != m_mesh->m_edgeSet.end(); eit++)
    {
        (*eit)->m_edgeNodes.clear();
    }

    FaceSet::iterator fit;
    for(fit = m_mesh->m_faceSet.begin(); fit != m_mesh->m_faceSet.end(); fit++)
    {
        (*fit)->m_faceNodes.clear();
    }

    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        vector<NodeSharedPtr> empty;
        m_mesh->m_element[m_mesh->m_expDim][i]->SetVolumeNodes(empty);
    }

}

}
}
