////////////////////////////////////////////////////////////////////////////////
//
//  File: Outputstl.cpp
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
//  Description: STL surface writer.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/MeshElements/Element.h>

#include "OutputSTL.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey OutputSTL::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eOutputModule, "stl"), OutputSTL::create, "Writes STL file.");

OutputSTL::OutputSTL(MeshSharedPtr m) : OutputModule(m)
{
}

OutputSTL::~OutputSTL()
{
}

void OutputSTL::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << "Outputstl: Writing file..." << endl;
    }

    ASSERTL0(m_mesh->m_expDim == 2 || m_mesh->m_expDim == 3,
             "Only 2D or 3D meshes are supported.");

    OpenStream();

    m_mshFile << std::scientific << setprecision(8);

    if (m_mesh->m_expDim == 2)
    {
        vector<ElementSharedPtr> &el = m_mesh->m_element[2];

        m_mshFile << "solid comp:" << 0 << endl;

        for (int i = 0; i < el.size(); ++i)
        {
            vector<NodeSharedPtr> ns = el[i]->GetVertexList();

            Array<OneD, NekDouble> tmp(3, 0.0);
            tmp[0] = (ns[1]->m_y - ns[0]->m_y) * (ns[2]->m_z - ns[0]->m_z) -
                     (ns[1]->m_z - ns[0]->m_z) * (ns[2]->m_y - ns[0]->m_y);
            tmp[1] = (ns[1]->m_z - ns[0]->m_z) * (ns[2]->m_x - ns[0]->m_x) -
                     (ns[1]->m_x - ns[0]->m_x) * (ns[2]->m_z - ns[0]->m_z);
            tmp[2] = (ns[1]->m_x - ns[0]->m_x) * (ns[2]->m_y - ns[0]->m_y) -
                     (ns[1]->m_y - ns[0]->m_y) * (ns[2]->m_x - ns[0]->m_x);

            NekDouble mt = tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2];
            mt           = sqrt(mt);
            tmp[0] /= mt;
            tmp[1] /= mt;
            tmp[2] /= mt;

            m_mshFile << "facet normal " << tmp[0] << " " << tmp[1] << " "
                      << tmp[2] << endl;
            m_mshFile << "outer loop" << endl;
            for (int j = 0; j < ns.size(); j++)
            {
                m_mshFile << "vertex " << ns[j]->m_x << " " << ns[j]->m_y << " "
                          << ns[j]->m_z << endl;
            }
            m_mshFile << "endloop" << endl << "endfacet" << endl;
        }

        m_mshFile << "endsolid" << endl;
        return;
    }

    for (auto &it : m_mesh->m_composite)
    {
        if (it.second->m_tag != "F")
        {
            continue;
        }

        m_mshFile << "solid comp:" << it.second->m_id << endl;

        vector<ElementSharedPtr> el = it.second->m_items;

        for (int i = 0; i < el.size(); i++)
        {
            vector<NodeSharedPtr> ns = el[i]->GetVertexList();

            Array<OneD, NekDouble> tmp(3, 0.0);
            tmp[0] = (ns[1]->m_y - ns[0]->m_y) * (ns[2]->m_z - ns[0]->m_z) -
                     (ns[1]->m_z - ns[0]->m_z) * (ns[2]->m_y - ns[0]->m_y);
            tmp[1] = (ns[1]->m_z - ns[0]->m_z) * (ns[2]->m_x - ns[0]->m_x) -
                     (ns[1]->m_x - ns[0]->m_x) * (ns[2]->m_z - ns[0]->m_z);
            tmp[2] = (ns[1]->m_x - ns[0]->m_x) * (ns[2]->m_y - ns[0]->m_y) -
                     (ns[1]->m_y - ns[0]->m_y) * (ns[2]->m_x - ns[0]->m_x);

            NekDouble mt = tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2];
            mt           = sqrt(mt);
            tmp[0] /= mt;
            tmp[1] /= mt;
            tmp[2] /= mt;

            m_mshFile << "facet normal " << tmp[0] << " " << tmp[1] << " "
                      << endl;
            m_mshFile << "outer loop" << endl;
            for (int j = 0; j < ns.size(); j++)
            {
                m_mshFile << "vertex " << ns[j]->m_x << " " << ns[j]->m_y << " "
                          << ns[j]->m_z << endl;
            }
            m_mshFile << "endloop" << endl << "endfacet" << endl;
        }
        m_mshFile << "endsolid" << endl;
    }
}
}
}
