////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessExtractTetPrismInterface.cpp
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
//  Description: Extract the interface between prismatic and tetrahedral
//  elements.
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <SpatialDomains/MeshGraph.h>

#include <NekMeshUtils/MeshElements/Element.h>

#include "ProcessExtractTetPrismInterface.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessExtractTetPrismInterface::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "extracttetprism"),
        ProcessExtractTetPrismInterface::create,
        "Process elements to extract the faces between "
        "tets and prisms.");

ProcessExtractTetPrismInterface::ProcessExtractTetPrismInterface(
    MeshSharedPtr m)
    : ProcessModule(m)
{
}

ProcessExtractTetPrismInterface::~ProcessExtractTetPrismInterface()
{
}

void ProcessExtractTetPrismInterface::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << "ProcessExtractTetPrismInterface: Extracting tet-prism "
             << "interface... " << endl;
    }

    ASSERTL0(m_mesh->m_expDim == 3,
             "The prism-tet interface module"
             " only works for three-dimensional meshes.");

    // Copy 3D elements and 2D boundary elements and clear existing.
    vector<ElementSharedPtr> el    = m_mesh->m_element[m_mesh->m_expDim];
    vector<ElementSharedPtr> bndEl = m_mesh->m_element[m_mesh->m_expDim - 1];
    m_mesh->m_element[m_mesh->m_expDim].clear();
    m_mesh->m_element[m_mesh->m_expDim - 1].clear();

    // Extract prismatic elements.
    for (int i = 0; i < el.size(); ++i)
    {
        ElementSharedPtr elmt = el[i];

        if (elmt->GetConf().m_e == LibUtilities::ePrism)
        {
            m_mesh->m_element[m_mesh->m_expDim].push_back(elmt);
        }
    }

    ASSERTL0(m_mesh->m_element[m_mesh->m_expDim].size() > 0,
             "Mesh does not contain any prismatic elements!");

    // Extract boundary region already associated with prisms
    // (i.e. outer wall of the computational domain)
    for (auto &face : m_mesh->m_faceSet)
    {
        if (face->m_elLink.size() == 1)
        {
            ElementSharedPtr el = face->m_elLink[0].first.lock();

            if (el->GetConf().m_e != LibUtilities::eTetrahedron)
            {
                m_mesh->m_element[m_mesh->m_expDim - 1].push_back(
                    bndEl[el->GetBoundaryLink(face->m_elLink[0].second)]);
            }
        }
    }

    // Now extract prismatic faces that are not connected to any other
    // elements, which denotes the prism/tet boundary.
    for (auto &face : m_mesh->m_faceSet)
    {
        if (face->m_elLink.size() != 1)
        {
            ElementSharedPtr el1 = face->m_elLink[0].first.lock();
            ElementSharedPtr el2 = face->m_elLink[1].first.lock();

            if ((el1->GetConf().m_e == LibUtilities::ePrism &&
                 el2->GetConf().m_e == LibUtilities::eTetrahedron) ||
                (el2->GetConf().m_e == LibUtilities::ePrism &&
                 el1->GetConf().m_e == LibUtilities::eTetrahedron))
            {
                // Create a new linear triangle from face for boundary.
                vector<NodeSharedPtr> nodeList(3);
                vector<int> tags(1);
                tags[0] = m_mesh->m_composite.size();

                nodeList = face->m_vertexList;
                ElmtConfig conf(
                    LibUtilities::eTriangle, 1, false, false, false);
                ElementSharedPtr tri = GetElementFactory().CreateInstance(
                    LibUtilities::eTriangle, conf, nodeList, tags);

                m_mesh->m_element[m_mesh->m_expDim - 1].push_back(tri);
            }
        }
    }

    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();
}
}
}
