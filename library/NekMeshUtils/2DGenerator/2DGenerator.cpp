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

#include <LibUtilities/BasicUtils/ParseUtils.hpp>
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
    m_config["blcurves"] =
        ConfigOption(false, "0", "Generate parallelograms on these curves");
    m_config["blthick"] =
        ConfigOption(false, "0", "Parallelogram layer thickness");
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
        tags.push_back((*it)->m_parentCAD->GetId());

        ElementSharedPtr E2 = GetElementFactory().CreateInstance(
            LibUtilities::eSegment, conf, ns, tags);

        m_mesh->m_element[1].push_back(E2);
    }

    if (m_config["blcurves"].beenSet)
    {
        MakeBL();
    }

    m_mesh->m_element[1].clear();

    /*

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

    */

    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();

    Report();
}

void Generator2D::MakeBL()
{
    if (m_mesh->m_verbose)
    {
        cout << endl << "\tBoundary layer meshing:" << endl << endl;
    }

    // identify the nodes which will become the boundary layer.

    vector<unsigned int> blCurves;
    ParseUtils::GenerateSeqVector(m_config["blcurves"].as<string>().c_str(),
                                  blCurves);
    NekDouble thickness = m_config["blthick"].as<NekDouble>();

    map<NodeSharedPtr, set<NodeSharedPtr> > nodesAdjacent;

    for (vector<unsigned int>::iterator it = blCurves.begin();
         it != blCurves.end(); ++it)
    {
        vector<NodeSharedPtr> localNodes = m_curvemeshes[*it]->GetMeshPoints();

        for (vector<NodeSharedPtr>::iterator it2 = localNodes.begin();
             it2 != localNodes.end(); ++it2)
        {
            if (it2 != localNodes.begin())
            {
                nodesAdjacent[*it2].insert(*(it2 - 1));
            }
            if (it2 != localNodes.end() - 1)
            {
                nodesAdjacent[*it2].insert(*(it2 + 1));
            }
        }
    }

    // on each node calculate a normal

    map<NodeSharedPtr, NodeSharedPtr> nodesNormal;

    for (map<NodeSharedPtr, set<NodeSharedPtr> >::iterator it =
             nodesAdjacent.begin();
         it != nodesAdjacent.end(); ++it)
    {
        size_t n = it->second.size();
        ASSERTL0(n == 2, "node not properly connected");

        Array<OneD, NekDouble> p1 = (*(it->second.begin()))->GetLoc();
        Array<OneD, NekDouble> p2 = it->first->GetLoc();
        Array<OneD, NekDouble> p3 = (*(it->second.rbegin()))->GetLoc();

        Node N12(0, p1[1] - p2[1], p2[0] - p1[0], 0);
        N12 /= sqrt(N12.abs2());

        Node N23(0, p2[1] - p3[1], p3[0] - p2[0], 0);
        N23 /= sqrt(N23.abs2());

        // NodeSharedPtr Nmean (new Node(N12 + N23));
        NodeSharedPtr Nmean = MemoryManager<Node>::AllocateSharedPtr(N12 + N23);
        *Nmean *= thickness / sqrt(Nmean->abs2());
        *Nmean += *(it->first);

        nodesNormal[it->first] = Nmean;
    }

    // create quadrilerals

    for (map<NodeSharedPtr, set<NodeSharedPtr> >::iterator it =
             nodesAdjacent.begin();
         it != nodesAdjacent.end(); ++it)
    {
        vector<NodeSharedPtr> ns;
        ns.push_back(it->first);
        ns.push_back(*(it->second.begin()));
        ns.push_back(nodesNormal[*(it->second.begin())]);
        ns.push_back(nodesNormal[it->first]);

        ElmtConfig conf(LibUtilities::eQuadrilateral, 1, false, false);

        vector<int> tags;
        tags.push_back(102);

        ElementSharedPtr E = GetElementFactory().CreateInstance(
            LibUtilities::eQuadrilateral, conf, ns, tags);

        m_mesh->m_element[2].push_back(E);
    }
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
