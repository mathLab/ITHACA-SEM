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
#include <cmath>

#include <NekMeshUtils/2DGenerator/2DGenerator.h>

#include <NekMeshUtils/Octree/Octree.h>

#include <LibUtilities/BasicUtils/ParseUtils.h>
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
    m_config["smoothbl"] =
        ConfigOption(true, "0", "Smooth the BL normal directions to avoid "
                                "(nearly) intersecting normals");
    m_config["spaceoutbl"] = ConfigOption(
        false, "0.5", "Threshold to space out BL according to Delta");
    m_config["nospaceoutsurf"] =
        ConfigOption(false, "", "Surfaces where spacing out shouldn't be used");
}

Generator2D::~Generator2D()
{
}

void Generator2D::Process()
{
    // Check that cad is 2D
    Array<OneD, NekDouble> bndBox = m_mesh->m_cad->GetBoundingBox();
    ASSERTL0(fabs(bndBox[5] - bndBox[4]) < 1.0e-7, "CAD isn't 2D");

    if (m_mesh->m_verbose)
    {
        cout << endl << "2D meshing" << endl;
        cout << endl << "\tCurve meshing:" << endl << endl;
    }
    m_mesh->m_numNodes = m_mesh->m_cad->GetNumVerts();
    m_thickness_ID =
        m_thickness.DefineFunction("x y z", m_config["blthick"].as<string>());
    ParseUtils::GenerateSeqVector(m_config["blcurves"].as<string>(),
                                  m_blCurves);

    // find the ends of the BL curves
    if (m_config["blcurves"].beenSet)
    {
        FindBLEnds();
    }

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

            // Fheck if this curve is at an end of the BL
            // If so, define an offset for the second node, corresponding to the
            // BL thickness
            if (m_blends.count(i))
            {
                vector<CADVertSharedPtr> vertices =
                    m_mesh->m_cad->GetCurve(i)->GetVertex();
                Array<OneD, NekDouble> loc;
                NekDouble t;

                // offset needed at first node (or both)
                if (m_blends[i] == 0 || m_blends[i] == 2)
                {
                    loc = vertices[0]->GetLoc();
                    t   = m_thickness.Evaluate(m_thickness_ID, loc[0], loc[1],
                                             loc[2], 0.0);
                    m_curvemeshes[i]->SetOffset(0, t);
                }
                // offset needed at second node (or both)
                if (m_blends[i] == 1 || m_blends[i] == 2)
                {
                    loc = vertices[1]->GetLoc();
                    t   = m_thickness.Evaluate(m_thickness_ID, loc[0], loc[1],
                                             loc[2], 0.0);
                    m_curvemeshes[i]->SetOffset(1, t);
                }
            }
        }
        else
        {
            m_curvemeshes[i] = MemoryManager<CurveMesh>::AllocateSharedPtr(
                i, m_mesh, m_config["blthick"].as<string>());
        }
        m_curvemeshes[i]->Mesh();
    }

    ////////
    // consider periodic curves

    if (m_config["periodic"].beenSet)
    {
        PeriodicPrep();
        MakePeriodic();
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

        // If the BL doesn't form closed loops, we need to remove the outside
        // nodes from the curve meshes
        for (auto &ic : m_blends)
        {
            vector<NodeSharedPtr> nodes =
                m_curvemeshes[ic.first]->GetMeshPoints();

            if (ic.second == 0 || ic.second == 2)
            {
                nodes.erase(nodes.begin());
            }
            if (ic.second == 1 || ic.second == 2)
            {
                nodes.erase(nodes.end() - 1);
            }

            // Rebuild the curvemesh without the first node, the last node or
            // both
            m_curvemeshes[ic.first] =
                MemoryManager<CurveMesh>::AllocateSharedPtr(ic.first, m_mesh,
                                                            nodes);
        }
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
            i, m_mesh, m_curvemeshes, 99 + i);
        m_facemeshes[i]->Mesh();
    }

    ////////////////////////////////////

    EdgeSet::iterator it;
    for (auto &it : m_mesh->m_edgeSet)
    {
        vector<NodeSharedPtr> ns;
        ns.push_back(it->m_n1);
        ns.push_back(it->m_n2);
        // for each iterator create a LibUtilities::eSegement
        // push segment into m_mesh->m_element[1]
        // tag for the elements shoudl be the CAD number of the curves
        ElmtConfig conf(LibUtilities::eSegment, 1, false, false);
        vector<int> tags;
        tags.push_back(it->m_parentCAD->GetId());
        ElementSharedPtr E2 = GetElementFactory().CreateInstance(
            LibUtilities::eSegment, conf, ns, tags);
        m_mesh->m_element[1].push_back(E2);
    }

    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();
    Report();
}

void Generator2D::FindBLEnds()
{
    // Set of CAD vertices
    // Vertices of each curve are added to the set if not found and removed from
    // the set if found
    // This leaves us with a set of vertices that are at the end of BL open
    // loops
    set<CADVertSharedPtr> cadverts;

    for (auto &it : m_blCurves)
    {
        vector<CADVertSharedPtr> vertices =
            m_mesh->m_cad->GetCurve(it)->GetVertex();

        for (auto &iv : vertices)
        {
            set<CADVertSharedPtr>::iterator is = cadverts.find(iv);

            if (is != cadverts.end())
            {
                cadverts.erase(is);
            }
            else
            {
                cadverts.insert(iv);
            }
        }
    }

    // Build m_blends based on the previously constructed set of vertices
    // m_blends is a map of curve number (the curves right outside the BL open
    // loops) to the offset node number: 0, 1 or 2 (for both)
    for (int i = 1; i <= m_mesh->m_cad->GetNumCurve(); ++i)
    {
        if (find(m_blCurves.begin(), m_blCurves.end(), i) != m_blCurves.end())
        {
            continue;
        }

        vector<CADVertSharedPtr> vertices =
            m_mesh->m_cad->GetCurve(i)->GetVertex();

        for (int j = 0; j < 2; ++j)
        {
            if (!cadverts.count(vertices[j]))
            {
                continue;
            }

            if (m_blends.count(i))
            {
                m_blends[i] = 2;
            }
            else
            {
                m_blends[i] = j;
            }
        }
    }
}

void Generator2D::MakeBLPrep()
{
    if (m_mesh->m_verbose)
    {
        cout << endl << "\tBoundary layer meshing:" << endl << endl;
    }

    // identify the nodes and edges which will become the boundary layer.

    for (auto &it : m_blCurves)
    {
        vector<EdgeSharedPtr> localedges = m_curvemeshes[it]->GetMeshEdges();
        for (auto &ie : localedges)
        {
            m_blEdges.push_back(ie);
            m_nodesToEdge[ie->m_n1].push_back(ie);
            m_nodesToEdge[ie->m_n2].push_back(ie);
        }
    }
}

void Generator2D::MakeBL(int faceid)
{
    map<int, Array<OneD, NekDouble>> edgeNormals;
    int eid = 0;
    for (auto &it : m_blCurves)
    {
        CADOrientation::Orientation edgeo =
            m_mesh->m_cad->GetCurve(it)->GetOrienationWRT(faceid);
        vector<EdgeSharedPtr> es = m_curvemeshes[it]->GetMeshEdges();
        // on each !!!EDGE!!! calculate a normal
        // always to the left unless edgeo is 1
        // normal must be done in the parametric space (and then projected back)
        // because of face orientation
        for (auto &ie : es)
        {
            ie->m_id = eid++;
            Array<OneD, NekDouble> p1, p2;
            p1 = ie->m_n1->GetCADSurfInfo(faceid);
            p2 = ie->m_n2->GetCADSurfInfo(faceid);
            Array<OneD, NekDouble> n(2);
            n[0] = p1[1] - p2[1];
            n[1] = p2[0] - p1[0];
            if (edgeo == CADOrientation::eBackwards)
            {
                n[0] *= -1.0;
                n[1] *= -1.0;
            }
            NekDouble mag = sqrt(n[0] * n[0] + n[1] * n[1]);
            n[0] /= mag;
            n[1] /= mag;
            Array<OneD, NekDouble> np(2);
            np[0] = p1[0] + n[0];
            np[1] = p1[1] + n[1];
            Array<OneD, NekDouble> loc  = ie->m_n1->GetLoc();
            Array<OneD, NekDouble> locp = m_mesh->m_cad->GetSurf(faceid)->P(np);
            n[0] = locp[0] - loc[0];
            n[1] = locp[1] - loc[1];
            mag  = sqrt(n[0] * n[0] + n[1] * n[1]);
            n[0] /= mag;
            n[1] /= mag;
            edgeNormals[ie->m_id] = n;
        }
    }

    bool adjust           = m_config["bltadjust"].beenSet;
    NekDouble divider     = m_config["bltadjust"].as<NekDouble>();
    bool adjustEverywhere = m_config["adjustblteverywhere"].beenSet;
    bool smoothbl         = m_config["smoothbl"].beenSet;
    bool spaceoutbl       = m_config["spaceoutbl"].beenSet;
    NekDouble spaceoutthr = m_config["spaceoutbl"].as<NekDouble>();

    if (divider < 2.0)
    {
        WARNINGL0(false, "BndLayerAdjustment too low, corrected to 2.0");
        divider = 2.0;
    }

    map<NodeSharedPtr, NodeSharedPtr> nodeNormals;
    for (auto &it : m_nodesToEdge)
    {
        ASSERTL0(it.second.size() == 1 || it.second.size() == 2,
                 "weirdness, most likely bl_surfs are incorrect");

        // If node at the end of the BL open loop, the "normal node" isn't
        // constructed by computing a normal but found on the adjacent curve
        if (it.second.size() == 1)
        {
            vector<CADCurveSharedPtr> curves = it.first->GetCADCurves();

            vector<EdgeSharedPtr> edges =
                m_curvemeshes[curves[0]->GetId()]->GetMeshEdges();
            vector<EdgeSharedPtr>::iterator ie =
                find(edges.begin(), edges.end(), it.second[0]);
            int rightCurve =
                (ie == edges.end()) ? curves[0]->GetId() : curves[1]->GetId();

            vector<NodeSharedPtr> nodes =
                m_curvemeshes[rightCurve]->GetMeshPoints();
            nodeNormals[it.first] =
                (nodes[0] == it.first) ? nodes[1] : nodes[nodes.size() - 2];

            continue;
        }

        Array<OneD, NekDouble> n(3, 0.0);
        Array<OneD, NekDouble> n1 = edgeNormals[it.second[0]->m_id];
        Array<OneD, NekDouble> n2 = edgeNormals[it.second[1]->m_id];
        n[0]          = (n1[0] + n2[0]) / 2.0;
        n[1]          = (n1[1] + n2[1]) / 2.0;
        NekDouble mag = sqrt(n[0] * n[0] + n[1] * n[1]);
        n[0] /= mag;
        n[1] /= mag;
        NekDouble t = m_thickness.Evaluate(m_thickness_ID, it.first->m_x,
                                           it.first->m_y, 0.0, 0.0);
        // Adjust thickness according to angle between normals
        if (adjust)
        {
            if (adjustEverywhere || it.first->GetNumCadCurve() > 1)
            {
                NekDouble angle = acos(n1[0] * n2[0] + n1[1] * n2[1]);
                angle           = (angle > M_PI) ? 2 * M_PI - angle : angle;
                t /= cos(angle / divider);
            }
        }

        n[0]             = n[0] * t + it.first->m_x;
        n[1]             = n[1] * t + it.first->m_y;
        NodeSharedPtr nn = std::shared_ptr<Node>(
            new Node(m_mesh->m_numNodes++, n[0], n[1], 0.0));
        CADSurfSharedPtr s = m_mesh->m_cad->GetSurf(faceid);
        Array<OneD, NekDouble> uv = s->locuv(n);
        nn->SetCADSurf(s, uv);
        nodeNormals[it.first] = nn;
    }

    // Check for any intersecting boundary layer normals and smooth them if
    // needed
    if (smoothbl)
    {
        // Nodes that need normal smoothing and their unit normal
        map<NodeSharedPtr, vector<NodeSharedPtr>> unitNormals;
        // Nodes that need normal smoothing and their BL thickness
        map<NodeSharedPtr, NekDouble> dist;

        int count = 0;

        do
        {
            unitNormals.clear();
            dist.clear();

            for (const auto &it : m_blEdges)
            {
                // Line intersection based on
                // https://stackoverflow.com/a/565282/7241595
                NodeSharedPtr p = it->m_n1;
                NodeSharedPtr q = it->m_n2;

                Node r = *nodeNormals[p] - *p;
                Node s = *nodeNormals[q] - *q;

                NekDouble d = r.curl(s).m_z;

                // Should probably use tolerance to check parallelism
                if (d == 0)
                {
                    continue;
                }

                NekDouble t = (*q - *p).curl(s).m_z / d;
                NekDouble u = (*q - *p).curl(r).m_z / d;

                // Check for intersection of the infinite continuation of one
                // normal with the other. A tolerance of 0.5 times the length of
                // thenormalis used. Could maybe be decreased to a less
                // aggressive value.
                if ((-0.5 < t && t <= 1.5) || (-0.5 < u && u <= 1.5))
                {
                    dist[p] = sqrt(r.abs2());
                    dist[q] = sqrt(s.abs2());

                    NodeSharedPtr sum =
                        make_shared<Node>(r / dist[p] + s / dist[q]);

                    unitNormals[p].push_back(sum);
                    unitNormals[q].push_back(sum);
                }
            }

            // Smooth each normal one by one
            for (const auto &it : unitNormals)
            {
                Node avg(0, 0.0, 0.0, 0.0);

                for (const auto &i : it.second)
                {
                    avg += *i;
                }

                avg /= sqrt(avg.abs2());

                // Create new BL node with smoothed normal
                NodeSharedPtr nn = std::shared_ptr<Node>(
                    new Node(nodeNormals[it.first]->GetID(),
                             it.first->m_x + avg.m_x * dist[it.first],
                             it.first->m_y + avg.m_y * dist[it.first], 0.0));
                CADSurfSharedPtr s =
                    *nodeNormals[it.first]->GetCADSurfs().begin();
                Array<OneD, NekDouble> uv = s->locuv(nn->GetLoc());
                nn->SetCADSurf(s, uv);

                nodeNormals[it.first] = nn;
            }
        } while (unitNormals.size() && count++ < 50);

        if (m_mesh->m_verbose)
        {
            if (count < 50)
            {
                cout << "\t\tNormals smoothed in " << count << " iterations."
                     << endl;
            }
            else
            {
                cout << "\t\tNormals smoothed. Algorithm didn't converge after "
                     << count << " iterations." << endl;
            }
        }
    }

    // Space out the outter BL nodes to better fit the local required Delta
    if (spaceoutbl)
    {
        if (spaceoutthr < 0.0 || spaceoutthr > 1.0)
        {
            WARNINGL0(false, "The boundary layer space out threshold should be "
                             "between 0 and 1. It will now be adjusted to "
                             "0.5.");
            spaceoutthr = 0.5;
        }

        vector<unsigned int> nospaceoutsurf;
        ParseUtils::GenerateSeqVector(m_config["nospaceoutsurf"].as<string>(),
                                      nospaceoutsurf);

        // List of connected nodes at need spacing out
        vector<deque<NodeSharedPtr>> nodesToMove;

        int count = 0;

        // This will supposedly spread the number of nodes to be moved until
        // sufficient space is found
        do
        {
            nodesToMove.clear();

            // Find which nodes need to be spaced out
            for (const auto &ie : m_blEdges)
            {
                auto it = find(nospaceoutsurf.begin(), nospaceoutsurf.end(),
                               ie->m_parentCAD->GetId());
                if (it != nospaceoutsurf.end())
                {
                    continue;
                }

                NodeSharedPtr n1 = nodeNormals[ie->m_n1];
                NodeSharedPtr n2 = nodeNormals[ie->m_n2];

                NekDouble targetD =
                    m_mesh->m_octree->Query(((*n1 + *n2) / 2.0).GetLoc());
                NekDouble realD = sqrt((*n1 - *n2).abs2());

                // Add nodes if condition fulfilled
                if (realD < spaceoutthr * targetD)
                {
                    bool connected = false;

                    for (auto &il : nodesToMove)
                    {
                        if (il.front() == n1)
                        {
                            il.push_front(n2);
                            connected = true;
                            break;
                        }
                        if (il.front() == n2)
                        {
                            il.push_front(n1);
                            connected = true;
                            break;
                        }
                        if (il.back() == n1)
                        {
                            il.push_back(n2);
                            connected = true;
                            break;
                        }
                        if (il.back() == n2)
                        {
                            il.push_back(n1);
                            connected = true;
                            break;
                        }
                    }

                    // Create new set of connected nodes if necessary
                    if (!connected)
                    {
                        deque<NodeSharedPtr> newList;
                        newList.push_back(n1);
                        newList.push_back(n2);

                        nodesToMove.push_back(newList);
                    }
                }
            }

            for (int i = 0;; ++i)
            {
                // Reconnect sets of connected nodes together if need be. Done
                // once before ad once after expanding the set by one node to
                // find extra space.
                for (int i1 = 0; i1 < nodesToMove.size(); ++i1)
                {
                    NodeSharedPtr n11 = nodesToMove[i1].front();
                    NodeSharedPtr n12 = nodesToMove[i1].back();

                    for (int i2 = i1 + 1; i2 < nodesToMove.size(); ++i2)
                    {
                        NodeSharedPtr n21 = nodesToMove[i2].front();
                        NodeSharedPtr n22 = nodesToMove[i2].back();

                        if (n11 == n21 || n11 == n22 || n12 == n21 ||
                            n12 == n22)
                        {
                            if (n11 == n21 || n12 == n22)
                            {
                                reverse(nodesToMove[i2].begin(),
                                        nodesToMove[i2].end());
                                n21 = nodesToMove[i2].front();
                                n22 = nodesToMove[i2].back();
                            }

                            if (n11 == n22)
                            {
                                nodesToMove[i1].insert(nodesToMove[i1].begin(),
                                                       nodesToMove[i2].begin(),
                                                       nodesToMove[i2].end() -
                                                           1);
                            }
                            else
                            {
                                nodesToMove[i1].insert(nodesToMove[i1].end(),
                                                       nodesToMove[i2].begin() +
                                                           1,
                                                       nodesToMove[i2].end());
                            }

                            nodesToMove.erase(nodesToMove.begin() + i2);
                            continue;
                        }
                    }
                }

                if (i >= 1)
                {
                    break;
                }

                set<EdgeSharedPtr> addedEdges;

                // Expand each set of connected nodes by one node to allow for
                // extra space
                for (auto &il : nodesToMove)
                {
                    NodeSharedPtr n11 = *(il.begin() + 0);
                    NodeSharedPtr n12 = *(il.begin() + 1);

                    NodeSharedPtr n13 = *(il.rbegin() + 1);
                    NodeSharedPtr n14 = *(il.rbegin() + 0);

                    for (const auto &ie : m_blEdges)
                    {
                        auto it =
                            find(nospaceoutsurf.begin(), nospaceoutsurf.end(),
                                 ie->m_parentCAD->GetId());
                        if (addedEdges.count(ie) || it != nospaceoutsurf.end())
                        {
                            continue;
                        }

                        NodeSharedPtr n21 = nodeNormals[ie->m_n1];
                        NodeSharedPtr n22 = nodeNormals[ie->m_n2];

                        NodeSharedPtr frontPush, backPush;

                        if (n11)
                        {
                            if (n11 == n21 && n12 != n22)
                            {
                                frontPush = n22;
                            }
                            else if (n11 == n22 && n12 != n21)
                            {
                                frontPush = n21;
                            }

                            if (frontPush)
                            {
                                il.push_front(frontPush);
                                n11.reset();
                                addedEdges.insert(ie);
                            }
                        }
                        if (n14 && !frontPush)
                        {
                            if (n14 == n21 && n13 != n22)
                            {
                                backPush = n22;
                            }
                            if (n14 == n22 && n13 != n21)
                            {
                                backPush = n21;
                            }

                            if (backPush)
                            {
                                il.push_back(backPush);
                                n14.reset();
                                addedEdges.insert(ie);
                            }
                        }

                        if (!n11 && !n14)
                        {
                            break;
                        }
                    }
                }
            }

            // Actual spacing out of the nodes. Done by simple linear
            // interpolation between the 2 end nodes. Weights come from the
            // required Delta or each edge.
            for (const auto &il : nodesToMove)
            {
                NodeSharedPtr ni = il.front();
                NodeSharedPtr nf = il.back();

                vector<NekDouble> deltas;
                NekDouble total = 0.0;

                for (int i = 0; i < il.size() - 1; ++i)
                {
                    NodeSharedPtr n1 = il[i];
                    NodeSharedPtr n2 = il[i + 1];

                    deltas.push_back(
                        m_mesh->m_octree->Query(((*n1 + *n2) / 2.0).GetLoc()));
                    total += deltas.back();
                }

                for (auto &id : deltas)
                {
                    id /= total;
                }

                NekDouble runningTotal = 0.0;

                for (int i = 1; i < il.size() - 1; ++i)
                {
                    runningTotal += deltas[i - 1];
                    Array<OneD, NekDouble> loc =
                        (*ni * (1.0 - runningTotal) + *nf * runningTotal)
                            .GetLoc();

                    Array<OneD, NekDouble> uv =
                        m_mesh->m_cad->GetSurf(faceid)->locuv(loc);

                    il[i]->Move(loc, faceid, uv);
                }
            }
        } while (nodesToMove.size() && count++ < 50);

        if (m_mesh->m_verbose)
        {
            if (count < 50)
            {
                cout << "\t\tBL spaced out in " << count << " iterations."
                     << endl;
            }
            else
            {
                cout << "\t\tBL spaced out. Algorithm didn't converge after "
                     << count << " iterations." << endl;
            }
        }
    }

    for (auto &it : m_blCurves)
    {
        CADOrientation::Orientation edgeo =
            m_mesh->m_cad->GetCurve(it)->GetOrienationWRT(faceid);
        vector<NodeSharedPtr> ns = m_curvemeshes[it]->GetMeshPoints();
        vector<NodeSharedPtr> newNs;
        for (auto &in : ns)
        {
            newNs.push_back(nodeNormals[in]);
        }
        m_curvemeshes[it] =
            MemoryManager<CurveMesh>::AllocateSharedPtr(it, m_mesh, newNs);
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

void Generator2D::PeriodicPrep()
{
    m_periodicPairs.clear();
    set<unsigned> periodic;

    // Build periodic curve pairs
    string s = m_config["periodic"].as<string>();
    vector<string> lines;

    boost::split(lines, s, boost::is_any_of(":"));

    for (auto &il : lines)
    {
        vector<string> tmp;
        boost::split(tmp, il, boost::is_any_of(","));

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
}

void Generator2D::MakePeriodic()
{
    // Override slave curves

    for (auto &ip : m_periodicPairs)
    {
        m_curvemeshes[ip.second]->PeriodicOverwrite(m_curvemeshes[ip.first]);
    }

    if (m_mesh->m_verbose)
    {
        cout << "\t\tPeriodic boundary conditions" << endl;
        for (auto &it : m_periodicPairs)
        {
            cout << "\t\t\tCurves " << it.first << " => " << it.second << endl;
        }
        cout << endl;
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
