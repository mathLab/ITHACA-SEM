///////////////////////////////////////////////////////////////////////////////
//
//  File: OutputCADfix.cpp
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
//  Description: CADfix file format output.
//
///////////////////////////////////////////////////////////////////////////////

#include "OutputCADfix.h"

using namespace Nektar::NekMeshUtils;
using namespace Nektar::SpatialDomains;

namespace Nektar
{
namespace Utilities
{
ModuleKey OutputCADfix::className =
    GetModuleFactory().RegisterCreatorFunction(ModuleKey(eOutputModule, "fbm"),
                                               OutputCADfix::create,
                                               "Appends to a CADfix database file.");

OutputCADfix::OutputCADfix(MeshSharedPtr m) : OutputModule(m)
{
    m_config["from"] = ConfigOption(false, "", "FBM file to load");
}

OutputCADfix::~OutputCADfix()
{
}

void OutputCADfix::Process()
{
    ModuleSharedPtr module = GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "loadcad"), m_mesh);
    module->RegisterConfig("CFIMesh", "");
    if (m_mesh->m_verbose)
    {
        module->RegisterConfig("verbose", "");
    }

    // If no input file specified, load output file
    module->RegisterConfig("filename", m_config["from"].beenSet ? m_config["from"].as<string>() : m_config["outfile"].as<string>());

    module->SetDefaults();
    module->Process();
    
    if (m_mesh->m_verbose)
    {
        cout << "OutputCADfix: Writing file..." << endl;
    }

    m_mesh->m_expDim   = 3;
    m_mesh->m_spaceDim = 3;

    m_cad           = std::dynamic_pointer_cast<CADSystemCFI>(m_mesh->m_cad);
    m_model         = m_cad->GetCFIModel();
    NekDouble scal  = m_cad->GetScaling();

    // Force order 2
    m_mesh->MakeOrder(2, LibUtilities::ePolyEvenlySpaced);

    // Add edge- and face-interior nodes to vertex set.
    for (auto &eIt : m_mesh->m_edgeSet)
    {
        m_mesh->m_vertexSet.insert(eIt->m_edgeNodes.begin(),
                                   eIt->m_edgeNodes.end());
    }

    for (auto &fIt : m_mesh->m_faceSet)
    {
        m_mesh->m_vertexSet.insert(fIt->m_faceNodes.begin(),
                                   fIt->m_faceNodes.end());
    }

    // Do second pass over elements for volume nodes.
    for (int d = 1; d <= 3; ++d)
    {
        for (int i = 0; i < m_mesh->m_element[d].size(); ++i)
        {
            ElementSharedPtr e = m_mesh->m_element[d][i];
            vector<NodeSharedPtr> volList = e->GetVolumeNodes();
            m_mesh->m_vertexSet.insert(volList.begin(), volList.end());
        }
    }

    // map of new node numbers
    map<NodeSharedPtr, cfi::Node*> newMap;

    // Write out nodes
    for (auto &it : m_mesh->m_vertexSet)
    {
        newMap[it] = m_model->createOrphanFenode(0, it->m_x / scal, it->m_y / scal, it->m_z / scal);
    }

    // Write out elements
    for (auto &el : m_mesh->m_element[3])
    {
        vector<cfi::Node*> nodes;
        
        // Assuming it's a tet
        int type = CFI_SUBTYPE_TE10;

        for (auto &node : el->GetVertexList())
        {
            nodes.push_back(newMap[node]);
        }
        
        for (auto &edge : el->GetEdgeList())
        {
            nodes.push_back(newMap[edge->m_edgeNodes[0]]);
        }
        
        if (el->GetTag() != "A")
        {
            // Now assuming it's a prism
            type = CFI_SUBTYPE_PE18;

            for (auto &face : el->GetFaceList())
            {
                if (face->m_faceNodes.size())
                {
                    nodes.push_back(newMap[face->m_faceNodes[0]]);
                }
            }

            if (el->GetTag() != "R")
            {
                // It's definitely a hex now
                type = CFI_SUBTYPE_HE27;

                for (auto &node : el->GetVolumeNodes())
                {
                    nodes.push_back(newMap[node]);
                }
            }
        }

        m_model->createOrphanElement(0, cfi::EntitySubtype(type), nodes);
    }

    m_model->saveCopy(m_config["outfile"].as<string>());
}

}
}
