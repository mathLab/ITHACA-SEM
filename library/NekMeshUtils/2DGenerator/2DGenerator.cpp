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
        // we need to do the boundary layer generation in a face by face basis
        MakeBLPrep();

        // Im going to do a horrendous trick to get the edge orientaion.
        // Going to activate the first routine of facemeshing without actually
        // face meshing, this will orientate the edgeloop objects (hopefully);
        // which can be used by the makebl command to know the normal
        // orienation
        for (int i = 1; i <= m_mesh->m_cad->GetNumSurf(); i++)
        {
            m_facemeshes[i] = MemoryManager<FaceMesh>::AllocateSharedPtr(
                i, m_mesh, m_curvemeshes, m_mesh->m_cad->GetNumSurf() > 100);

            m_facemeshes[i]->OrientateCurves();
            MakeBL(i, m_facemeshes[i]->GetEdges());
        }
    }

    m_mesh->m_element[1].clear();

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

void Generator2D::MakeBLPrep()
{
    if (m_mesh->m_verbose)
    {
        cout << endl << "\tBoundary layer meshing:" << endl << endl;
    }

    // identify the nodes which will become the boundary layer.
    ParseUtils::GenerateSeqVector(m_config["blcurves"].as<string>().c_str(),
                                  m_blCurves);
    m_thickness = m_config["blthick"].as<NekDouble>();

    for (vector<unsigned>::iterator it = m_blCurves.begin();
         it != m_blCurves.end(); ++it)
    {
        vector<EdgeSharedPtr> localedges = m_curvemeshes[*it]->GetMeshEdges();
        for (int i = 0; i < localedges.size(); i++)
        {
            m_nodesToEdge[localedges[i]->m_n1].push_back(localedges[i]);
            m_nodesToEdge[localedges[i]->m_n2].push_back(localedges[i]);
        }
    }
}

void Generator2D::MakeBL(int i, std::vector<EdgeLoop> e)
{
    // on each node calculate a normal
    map<NodeSharedPtr, NodeSharedPtr> nodeNormals;

    for (map<NodeSharedPtr, vector<EdgeSharedPtr> >::iterator it =
             m_nodesToEdge.begin();
         it != m_nodesToEdge.end(); ++it)
    {
        size_t n = it->second.size();
        cout << n << endl;
        ASSERTL0(n == 2, "node not properly connected");

        Array<OneD, NekDouble> p1 = it->second[0]->m_n1->GetLoc();
        Array<OneD, NekDouble> p2 = it->first->GetLoc();
        Array<OneD, NekDouble> p3 = it->second[1]->m_n2->GetLoc();

        if (p1 == p2)
        {
            p1 = it->second[0]->m_n2->GetLoc();
        }

        if (p3 == p2)
        {
            p3 = it->second[1]->m_n1->GetLoc();
        }

        Node N12(0, p1[1] - p2[1], p2[0] - p1[0], 0);
        N12 /= sqrt(N12.abs2());

        Node N23(0, p2[1] - p3[1], p3[0] - p2[0], 0);
        N23 /= sqrt(N23.abs2());

        NodeSharedPtr Nmean = boost::shared_ptr<Node>(new Node(m_mesh->m_numNodes++,
                        (N12.m_x+N23.m_x)/2.0,
                        (N12.m_y+N23.m_y)/2.0,
                        (N12.m_z+N23.m_z)/2.0));
        //NodeSharedPtr Nmean = MemoryManager<Node>::AllocateSharedPtr(N12 + N23);
        *Nmean *= m_thickness / sqrt(Nmean->abs2());
        *Nmean += *(it->first);

        vector<pair<int, CADSurfSharedPtr> > surfs = it->first->GetCADSurfs();
        for (vector<pair<int, CADSurfSharedPtr> >::iterator it2 = surfs.begin();
             it2 != surfs.end(); ++it2)
        {
            Nmean->SetCADSurf(it2->first, it2->second,
                              it->first->GetCADSurfInfo(it2->first));
        }

        nodeNormals[it->first] = Nmean;
    }

    // create quadrilerals

    for (map<NodeSharedPtr, vector<EdgeSharedPtr> >::iterator it =
             m_nodesToEdge.begin();
         it != m_nodesToEdge.end(); ++it)
    {
        vector<NodeSharedPtr> ns;

        ns.push_back(it->first);

        if (it->second[0]->m_n1 == it->first)
        {
            ns.push_back(it->second[0]->m_n2);
            ns.push_back(nodeNormals[it->second[0]->m_n2]);
        }
        else
        {
            ns.push_back(it->second[0]->m_n1);
            ns.push_back(nodeNormals[it->second[0]->m_n1]);
        }

        ns.push_back(nodeNormals[it->first]);

        ElmtConfig conf(LibUtilities::eQuadrilateral, 1, false, false);

        vector<int> tags;
        tags.push_back(102);

        ElementSharedPtr E = GetElementFactory().CreateInstance(
            LibUtilities::eQuadrilateral, conf, ns, tags);

        m_mesh->m_element[2].push_back(E);
    }

    for (vector<unsigned>::iterator it = m_blCurves.begin();
         it != m_blCurves.end(); ++it)
    {
        vector<NodeSharedPtr> oldNodes = m_curvemeshes[*it]->GetMeshPoints();
        vector<NodeSharedPtr> newNodes;

        for (vector<NodeSharedPtr>::iterator it2 = oldNodes.begin();
             it2 != oldNodes.end(); ++it2)
        {
            newNodes.push_back(nodeNormals[*it2]);
        }

        m_curvemeshes[*it] =
            MemoryManager<CurveMesh>::AllocateSharedPtr(*it, m_mesh, newNodes);
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
