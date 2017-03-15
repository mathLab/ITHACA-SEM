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
#include <math.h>

#include <NekMeshUtils/2DGenerator/2DGenerator.h>

#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/Progressbar.hpp>

#include <boost/algorithm/string.hpp>

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
        ConfigOption(false, "", "Generate parallelograms on these curves");
    m_config["blthick"] =
        ConfigOption(false, "0.0", "Parallelogram layer thickness");
    m_config["periodic"] =
        ConfigOption(false, "", "Set of pairs of periodic curves");
    m_config["bltadjust"] =
        ConfigOption(false, "2.0", "Boundary layer thickness adjustment");
    m_config["adjustblteverywhere"] =
        ConfigOption(true, "0", "Adjust thickness everywhere");
}

Generator2D::~Generator2D()
{
}

void Generator2D::Process()
{
    /*if (m_mesh->m_verbose)
    {
        cout << endl << "2D meshing" << endl;
        cout << endl << "\tCurve meshing:" << endl << endl;
    }

    m_mesh->m_numNodes = m_mesh->m_cad->GetNumVerts();

    set<unsigned> periodic;

    if (m_config["periodic"].beenSet)
    {
        m_periodicPairs.clear();

        // Build periodic curve pairs

        string s = m_config["periodic"].as<string>();
        vector<string> lines;

        boost::split(lines, s, boost::is_any_of(":"));

        for (vector<string>::iterator il = lines.begin(); il != lines.end();
             ++il)
        {
            vector<string> tmp;
            boost::split(tmp, *il, boost::is_any_of(","));

            ASSERTL0(tmp.size() == 2, "periodic pairs ill-defined");

            vector<unsigned> data(2);
            data[0] = boost::lexical_cast<unsigned>(tmp[0]);
            data[1] = boost::lexical_cast<unsigned>(tmp[1]);

            ASSERTL0(!periodic.count(data[0]), "curve already periodic");
            ASSERTL0(!periodic.count(data[1]), "curve already periodic");

            m_periodicPairs[data[0]] = data[1];
            periodic.insert(data[0]);
            periodic.insert(data[1]);
        }

        // Check compatibility

        for (map<unsigned, unsigned>::iterator it = m_periodicPairs.begin();
             it != m_periodicPairs.end(); ++it)
        {
            NekDouble L1 = m_mesh->m_cad->GetCurve(it->first)->GetTotLength();
            NekDouble L2 = m_mesh->m_cad->GetCurve(it->second)->GetTotLength();

            ASSERTL0(abs((L1 - L2) / L1) < 1.0e-3,
                     "periodic curves of different length");
        }
    }

    if (m_config["blcurves"].beenSet)
    {
        ParseUtils::GenerateSeqVector(m_config["blcurves"].as<string>().c_str(),
                                      m_blCurves);
        m_thickness_ID = m_thickness.DefineFunction(
            "x y z", m_config["blthick"].as<string>());
    }*/

    /*if (m_config["periodic"].beenSet)
    {
        // Override slave curves

        for (map<unsigned, unsigned>::iterator ip = m_periodicPairs.begin();
             ip != m_periodicPairs.end(); ++ip)
        {
            Array<OneD, NekDouble> A1 =
                m_curvemeshes[ip->first]->GetFirstPoint()->GetLoc();
            Array<OneD, NekDouble> A2 =
                m_curvemeshes[ip->first]->GetLastPoint()->GetLoc();
            Array<OneD, NekDouble> B1 =
                m_curvemeshes[ip->second]->GetFirstPoint()->GetLoc();
            Array<OneD, NekDouble> B2 =
                m_curvemeshes[ip->second]->GetLastPoint()->GetLoc();

            Array<OneD, NekDouble> T1(2);
            Array<OneD, NekDouble> T2(2);
            Array<OneD, NekDouble> dT(2);

            // Compute translation vector

            T1[0] = B1[0] - A1[0];
            T1[1] = B1[1] - A1[1];

            T2[0] = B2[0] - A2[0];
            T2[1] = B2[1] - A2[1];

            dT[0] = T2[0] - T1[0];
            dT[1] = T2[1] - T1[1];

            NekDouble dTmag = (dT[0] * dT[0] + dT[1] * dT[1]) /
                              (T1[0] * T1[0] + T1[1] * T1[1]);

            // Check if slave vector is reverse oriented

            bool reverse = false;

            if (dTmag > 1.0e-3)
            {
                reverse = true;

                T1[0] = B1[0] - A2[0];
                T1[1] = B1[1] - A2[1];

                T2[0] = B2[0] - A1[0];
                T2[1] = B2[1] - A1[1];

                dT[0] = T2[0] - T1[0];
                dT[1] = T2[1] - T1[1];

                dTmag = (dT[0] * dT[0] + dT[1] * dT[1]) /
                        (T1[0] * T1[0] + T1[1] * T1[1]);

                ASSERTL0(dTmag < 1.0e-3, "curve cannot be translated");
            }

            // Build vector of translated nodes

            vector<NodeSharedPtr> nodes =
                m_curvemeshes[ip->first]->GetMeshPoints();
            vector<NodeSharedPtr> nnodes;

            vector<pair<CADSurfSharedPtr, CADOrientation::Orientation> > surfs =
                m_curvemeshes[ip->second]->GetCADCurve()->GetAdjSurf();

            nnodes.push_back(m_curvemeshes[ip->second]->GetFirstPoint());

            for (vector<NodeSharedPtr>::iterator in = nodes.begin() + 1;
                 in != nodes.end() - 1; ++in)
            {
                Array<OneD, NekDouble> loc = (*in)->GetLoc();
                NodeSharedPtr nn = boost::shared_ptr<Node>(new Node(
                    m_mesh->m_numNodes++, loc[0] + T1[0], loc[1] + T1[1], 0.0));

                for (vector<pair<CADSurfSharedPtr,
                                 CADOrientation::Orientation> >::iterator is =
                         surfs.begin();
                     is != surfs.end(); ++is)
                {
                    nn->SetCADSurf(is->first->GetId(), is->first,
                                   is->first->locuv(nn->GetLoc()));
                }

                nn->SetCADCurve(ip->second,
                                m_curvemeshes[ip->second]->GetCADCurve(),
                                m_curvemeshes[ip->second]->GetCADCurve()->loct(
                                    nn->GetLoc()));

                nnodes.push_back(nn);
            }

            nnodes.push_back(m_curvemeshes[ip->second]->GetLastPoint());

            // Reverse internal nodes of the vector if necessary

            if (reverse)
            {
                std::reverse(++nnodes.begin(), --nnodes.end());
            }

            // Clean m_edgeSet and build new CurveMesh

            vector<EdgeSharedPtr> edges =
                m_curvemeshes[ip->second]->GetMeshEdges();
            for (vector<EdgeSharedPtr>::iterator ie = edges.begin();
                 ie != edges.end(); ++ie)
            {
                m_mesh->m_edgeSet.erase(*ie);
            }

            m_curvemeshes[ip->second] =
                MemoryManager<CurveMesh>::AllocateSharedPtr(ip->second, m_mesh,
                                                            nnodes, true);
        }

        if (m_mesh->m_verbose)
        {
            cout << "\t\tPeriodic boundary conditions" << endl;
            for (map<unsigned, unsigned>::iterator it = m_periodicPairs.begin();
                 it != m_periodicPairs.end(); ++it)
            {
                cout << "\t\t\tCurves " << it->first << " => " << it->second
                     << endl;
            }
            cout << endl;
        }
    }*/

    if (m_mesh->m_verbose)
    {
        cout << endl << "2D meshing" << endl;
        cout << endl << "\tCurve meshing:" << endl << endl;
    }
    m_mesh->m_numNodes = m_mesh->m_cad->GetNumVerts();
    m_thickness_ID =
        m_thickness.DefineFunction("x y z", m_config["blthick"].as<string>());
    ParseUtils::GenerateSeqVector(m_config["blcurves"].as<string>().c_str(),
                                  m_blCurves);

    // linear mesh all curves
    for (int i = 1; i <= m_mesh->m_cad->GetNumCurve(); i++)
    {
        if (m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(i, m_mesh->m_cad->GetNumCurve(),
                                           "Curve progress");
        }
        vector<unsigned int>::iterator f =
            find(m_blCurves.begin(), m_blCurves.end(), i);
        if (f == m_blCurves.end())
        {
            m_curvemeshes[i] =
                MemoryManager<CurveMesh>::AllocateSharedPtr(i, m_mesh);
        }
        else
        {
            m_curvemeshes[i] = MemoryManager<CurveMesh>::AllocateSharedPtr(
                i, m_mesh, m_config["blthick"].as<string>());
        }
        m_curvemeshes[i]->Mesh();
    }

    ////////////////////////////////////////

    if (m_config["blcurves"].beenSet)
    {
        // we need to do the boundary layer generation in a face by face basis
        MakeBLPrep();
        for (int i = 1; i <= m_mesh->m_cad->GetNumSurf(); i++)
        {
            MakeBL(i);
        }
    }


    /*if (m_mesh->m_verbose)
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
            i, m_mesh, m_curvemeshes, 99 + i);
        m_facemeshes[i]->Mesh();
    }*/


    ////////////////////////////////////

    /*EdgeSet::iterator it;
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
    }*/


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

void Generator2D::MakeBL(int faceid)
{
    map<int, Array<OneD, NekDouble> > edgeNormals;
    int eid = 0;
    for (vector<unsigned>::iterator it = m_blCurves.begin();
         it != m_blCurves.end(); ++it)
    {
        CADOrientation::Orientation edgeo =
            m_mesh->m_cad->GetCurve(*it)->GetOrienationWRT(faceid);
        vector<EdgeSharedPtr> es = m_curvemeshes[*it]->GetMeshEdges();
        // on each !!!EDGE!!! calculate a normal
        // always to the left unless edgeo is 1
        // normal must be done in the parametric space (and then projected back)
        // because of face orientation
        for (int j = 0; j < es.size(); j++)
        {
            es[j]->m_id = eid++;
            Array<OneD, NekDouble> p1, p2;
            p1 = es[j]->m_n1->GetCADSurfInfo(faceid);
            p2 = es[j]->m_n2->GetCADSurfInfo(faceid);
            if (edgeo == CADOrientation::eBackwards)
            {
                swap(p1, p2);
            }
            Array<OneD, NekDouble> n(2);
            n[0]          = p1[1] - p2[1];
            n[1]          = p2[0] - p1[0];
            NekDouble mag = sqrt(n[0] * n[0] + n[1] * n[1]);
            n[0] /= mag;
            n[1] /= mag;
            Array<OneD, NekDouble> np = p1;
            np[0] += n[0];
            np[1] += n[1];
            Array<OneD, NekDouble> loc  = es[j]->m_n1->GetLoc();
            Array<OneD, NekDouble> locp = m_mesh->m_cad->GetSurf(faceid)->P(np);
            n[0] = locp[0] - loc[0];
            n[1] = locp[1] - loc[1];
            mag  = sqrt(n[0] * n[0] + n[1] * n[1]);
            n[0] /= mag;
            n[1] /= mag;
            edgeNormals[es[j]->m_id] = n;
        }
    }

    bool adjust           = m_config["bltadjust"].beenSet;
    NekDouble divider     = m_config["bltadjust"].as<NekDouble>();
    bool adjustEverywhere = m_config["adjustblteverywhere"].beenSet;

    if (divider < 2.0)
    {
        WARNINGL1(false, "BndLayerAdjustment too low, corrected to 2.0");
        divider = 2.0;
    }

    map<NodeSharedPtr, NodeSharedPtr> nodeNormals;
    map<NodeSharedPtr, vector<EdgeSharedPtr> >::iterator it;
    for (it = m_nodesToEdge.begin(); it != m_nodesToEdge.end(); it++)
    {
        Array<OneD, NekDouble> n(3,0.0);
        ASSERTL0(it->second.size() == 2,
                 "wierdness, most likely bl_surfs are incorrect");
        Array<OneD, NekDouble> n1 = edgeNormals[it->second[0]->m_id];
        Array<OneD, NekDouble> n2 = edgeNormals[it->second[1]->m_id];
        n[0]          = (n1[0] + n2[0]) / 2.0;
        n[1]          = (n1[1] + n2[1]) / 2.0;
        NekDouble mag = sqrt(n[0] * n[0] + n[1] * n[1]);
        n[0] /= mag;
        n[1] /= mag;
        NekDouble t = m_thickness.Evaluate(m_thickness_ID, it->first->m_x,
                                           it->first->m_y, 0.0, 0.0);
        // Adjust thickness according to angle between normals
        if (adjust)
        {
            if (adjustEverywhere || it->first->GetNumCadCurve() > 1)
            {
                NekDouble angle = acos(n1[0] * n2[0] + n1[1] * n2[1]);
                angle           = (angle > M_PI) ? 2 * M_PI - angle : angle;
                t /= cos(angle / divider);
            }
        }

        n[0] = n[0] * t + it->first->m_x;
        n[1] = n[1] * t + it->first->m_y;
        NodeSharedPtr nn = boost::shared_ptr<Node>(
            new Node(m_mesh->m_numNodes++, n[0], n[1], 0.0));
        CADSurfSharedPtr s = m_mesh->m_cad->GetSurf(faceid);
        Array<OneD, NekDouble> uv = s->locuv(n);
        nn->SetCADSurf(faceid, s, uv);
        nodeNormals[it->first] = nn;
    }

    for (vector<unsigned>::iterator it = m_blCurves.begin();
         it != m_blCurves.end(); ++it)
    {
        CADOrientation::Orientation edgeo =
            m_mesh->m_cad->GetCurve(*it)->GetOrienationWRT(faceid);
        vector<NodeSharedPtr> ns = m_curvemeshes[*it]->GetMeshPoints();
        vector<NodeSharedPtr> newNs;
        for (int i = 0; i < ns.size(); i++)
        {
            newNs.push_back(nodeNormals[ns[i]]);
        }
        m_curvemeshes[*it] =
            MemoryManager<CurveMesh>::AllocateSharedPtr(*it, m_mesh, newNs);
        if (edgeo == CADOrientation::eBackwards)
        {
            reverse(ns.begin(), ns.end());
        }
        for (int i = 0; i < ns.size() - 1; ++i)
        {
            vector<NodeSharedPtr> qns;
            qns.push_back(ns[i]);
            qns.push_back(ns[i + 1]);
            qns.push_back(nodeNormals[ns[i + 1]]);
            qns.push_back(nodeNormals[ns[i]]);
            ElmtConfig conf(LibUtilities::eQuadrilateral, 1, false, false);
            vector<int> tags;
            tags.push_back(101);
            ElementSharedPtr E = GetElementFactory().CreateInstance(
                LibUtilities::eQuadrilateral, conf, qns, tags);
            E->m_parentCAD = m_mesh->m_cad->GetSurf(faceid);
            for (int j = 0; j < E->GetEdgeCount(); ++j)
            {
                pair<EdgeSet::iterator, bool> testIns;
                EdgeSharedPtr ed = E->GetEdge(j);
                // look for edge in m_mesh edgeset from curves
                EdgeSet::iterator s = m_mesh->m_edgeSet.find(ed);
                if (!(s == m_mesh->m_edgeSet.end()))
                {
                    ed = *s;
                    E->SetEdge(j, *s);
                }
            }
            m_mesh->m_element[2].push_back(E);
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
