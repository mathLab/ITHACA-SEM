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

void Generator2D::MakeBL(int faceid, vector<EdgeLoop> e)
{
    map<NodeSharedPtr, NodeSharedPtr> nodeNormals;

    for (vector<EdgeLoop>::iterator lit = e.begin(); lit != e.end(); ++lit)
    {
        for (int i = 0; i < lit->edges.size(); ++i)
        {
            int id    = lit->edges[i]->GetId();
            int edgeo = lit->edgeo[i];

            bool isBl = false;
            for (vector<unsigned>::iterator it = m_blCurves.begin();
                 it != m_blCurves.end(); ++it)
            {
                if (*it == id)
                {
                    isBl = true;
                    break;
                }
            }

            if (!isBl)
            {
                continue;
            }

            vector<NodeSharedPtr> nodes = m_curvemeshes[id]->GetMeshPoints();

            // on each node calculate a normal

            for (int ni = 0; ni < nodes.size(); ++ni)
            {
                NodeSharedPtr node = nodes[ni];

                if (nodeNormals.count(node))
                {
                    continue;
                }

                vector<EdgeSharedPtr> edges = m_nodesToEdge[node];

                Array<OneD, NekDouble> p2 = node->GetLoc();

                Array<OneD, NekDouble> p1 = (node == edges[0]->m_n1)
                                                ? edges[0]->m_n2->GetLoc()
                                                : edges[0]->m_n1->GetLoc();
                Array<OneD, NekDouble> p3 = (node == edges[1]->m_n1)
                                                ? edges[1]->m_n2->GetLoc()
                                                : edges[1]->m_n1->GetLoc();

                Node N12(0, p1[1] - p2[1], p2[0] - p1[0], 0);
                N12 /= sqrt(N12.abs2());

                Node N23(0, p2[1] - p3[1], p3[0] - p2[0], 0);
                N23 /= sqrt(N23.abs2());

                NodeSharedPtr Nmean = boost::shared_ptr<Node>(
                    new Node(m_mesh->m_numNodes++, N12.m_x + N23.m_x,
                             N12.m_y + N23.m_y, N12.m_z + N23.m_z));
                *Nmean *= m_thickness / sqrt(Nmean->abs2());

                if (i == 0 && ni == (edgeo ? nodes.size() - 1 : 0))
                {
                    *Nmean *= -1;
                }

                *Nmean += *node;

                Nmean->SetCADSurf(faceid, node->GetCADSurf(faceid),
                                  node->GetCADSurfInfo(faceid));

                nodeNormals[node] = Nmean;
            }

            // create quadrilerals

            for (int ni = 0; ni < nodes.size() - 1; ++ni)
            {
                NodeSharedPtr node = nodes[ni];

                vector<NodeSharedPtr> ns;

                ns.push_back(node);

                EdgeSharedPtr edge;

                if (i == 0 && ni == 0)
                {
                    edge = m_nodesToEdge[node][edgeo ? 1 : 0];
                }
                else
                {
                    edge = m_nodesToEdge[node][edgeo ? 0 : 1];
                }

                if (edge->m_n1 == node)
                {
                    ns.push_back(edge->m_n2);
                    ns.push_back(nodeNormals[edge->m_n2]);
                }
                else
                {
                    ns.push_back(edge->m_n1);
                    ns.push_back(nodeNormals[edge->m_n1]);
                }

                ns.push_back(nodeNormals[node]);

                ElmtConfig conf(LibUtilities::eQuadrilateral, 1, false, false);

                vector<int> tags;
                tags.push_back(102);

                ElementSharedPtr E = GetElementFactory().CreateInstance(
                    LibUtilities::eQuadrilateral, conf, ns, tags);

                m_mesh->m_element[2].push_back(E);
            }

            // replace old curvemesh with offset

            vector<NodeSharedPtr> newNodes;

            for (vector<NodeSharedPtr>::iterator nit = nodes.begin();
                 nit != nodes.end(); ++nit)
            {
                newNodes.push_back(nodeNormals[*nit]);
            }

            m_curvemeshes[id] = MemoryManager<CurveMesh>::AllocateSharedPtr(
                id, m_mesh, newNodes);
        }
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
