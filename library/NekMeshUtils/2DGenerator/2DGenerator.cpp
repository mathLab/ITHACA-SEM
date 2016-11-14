////////////////////////////////////////////////////////////////////////////////
//
//  File: Generator2D.cpp
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
//  Description: 2D generator object methods.
//
////////////////////////////////////////////////////////////////////////////////
#include <algorithm>

#include <NekMeshUtils/2DGenerator/2DGenerator.h>

#include <LibUtilities/BasicUtils/Progressbar.hpp>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

ModuleKey Generator2D::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "2dgenerator"), Generator2D::create,
    "Generates a 2D mesh");

Generator2D::Generator2D(MeshSharedPtr m) : ProcessModule(m)
{
}

Generator2D::~Generator2D()
{
}

void Generator2D::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << endl << "2D meshing" << endl;
        cout << endl << "\tCurve meshing:" << endl << endl;
    }

    m_mesh->m_numNodes = m_mesh->m_cad->GetNumVerts();

    // linear mesh all curves
    for (int i = 1; i <= m_mesh->m_cad->GetNumCurve(); i++)
    {
        if (m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(i, m_mesh->m_cad->GetNumCurve(),
                                           "Curve progress");
        }

        m_curvemeshes[i] =
            MemoryManager<CurveMesh>::AllocateSharedPtr(i, m_mesh);

        m_curvemeshes[i]->Mesh();
    }

    EdgeSet::iterator it;
    for (it = m_mesh->m_edgeSet.begin(); it != m_mesh->m_edgeSet.end(); it++)
    {
        vector<NodeSharedPtr> ns;
        ns.push_back((*it)->m_n1);
        ns.push_back((*it)->m_n2);

        // for each iterator create a LibUtilities::eSegement
        // push segment into m_mesh->m_element[1]
        // tag for the elements shoudl be the CAD number of the curves

        ElmtConfig conf(LibUtilities::eSegment, 1, false, false);

        vector<int> tags;
        tags.push_back((*it)->m_id);

        ElementSharedPtr E2 = GetElementFactory().CreateInstance(
            LibUtilities::eSegment, conf, ns, tags);
        
        E->CADSurfId = (*it)->m_id;

        vector<NodeSharedPtr> nods = E->GetVertexList();
        for (int j = 0; j < nods.size(); j++)
        {
            // nodes are already unique some will insert some wont
            m_localNodes.insert(nods[j]);
        }
        E->SetId(m_localElements.size());
        m_localElements.push_back(E);
    }

    for (int i = 0; i < m_localElements.size(); ++i)
    {
        for (int j = 0; j < m_localElements[i]->GetEdgeCount(); ++j)
        {
            pair<EdgeSet::iterator, bool> testIns;
            EdgeSharedPtr ed = m_localElements[i]->GetEdge(j);
            // look for edge in m_mesh edgeset from curves
            EdgeSet::iterator s = m_mesh->m_edgeSet.find(ed);
            if (!(s == m_mesh->m_edgeSet.end()))
            {
                ed = *s;
                m_localElements[i]->SetEdge(j, *s);
            }

            testIns = m_localEdges.insert(ed);

            if (testIns.second)
            {
                EdgeSharedPtr ed2 = *testIns.first;
                ed2->m_elLink.push_back(
                    pair<ElementSharedPtr, int>(m_localElements[i], j));
            }
            else
            {
                EdgeSharedPtr e2 = *(testIns.first);
                m_localElements[i]->SetEdge(j, e2);
                e2->m_elLink.push_back(
                    pair<ElementSharedPtr, int>(m_localElements[i], j));
            }
        }
    }

    // make new elements and add to list from list of nodes and connectivity
    // from triangle
    for (int i = 0; i < m_localElements.size(); i++)
    {
        vector<EdgeSharedPtr> e = m_localElements[i]->GetEdgeList();
        for (int j = 0; j < e.size(); j++)
        {
            e[j]->m_elLink.clear();
        }
        m_mesh->m_element[2].push_back(m_localElements[i]);
    }

    if (m_mesh->m_verbose)
    {
        cout << "\r                                                            "
                "    "
                "                             ";
        cout << scientific /*<< "\r\t\tFace " << m_id << endl*/
             << "\t\t\tNodes: " << m_localNodes.size() << endl
             << "\t\t\tEdges: " << m_localEdges.size() << endl
             << "\t\t\tTriangles: " << m_localElements.size() << endl
             // << "\t\t\tLoops: " << m_edgeloops.size() << endl
             // << "\t\t\tSTR: " << m_str << endl
             << endl;
    }

    if (m_mesh->m_verbose)
    {
        cout << endl << "\tFace meshing:" << endl << endl;
    }

    // linear mesh all surfaces
    for (int i = 1; i <= m_mesh->m_cad->GetNumSurf(); i++)
    {
        if (m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(i, m_mesh->m_cad->GetNumSurf(),
                                           "Face progress");
        }

        m_facemeshes[i] = MemoryManager<FaceMesh>::AllocateSharedPtr(
            i, m_mesh, m_curvemeshes, m_mesh->m_cad->GetNumSurf() > 100);

        m_facemeshes[i]->Mesh();
    }

    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();

    Report();
}

void Generator2D::Report()
{
    if (m_mesh->m_verbose)
    {
        int ns = m_mesh->m_vertexSet.size();
        int es = m_mesh->m_edgeSet.size();
        int ts = m_mesh->m_element[2].size();
        int ep = ns - es + ts;
        cout << endl << "\tSurface mesh statistics" << endl;
        cout << "\t\tNodes: " << ns << endl;
        cout << "\t\tEdges: " << es << endl;
        cout << "\t\tTriangles " << ts << endl;
        cout << "\t\tEuler-Poincaré characteristic: " << ep << endl;
    }
}
}
}
