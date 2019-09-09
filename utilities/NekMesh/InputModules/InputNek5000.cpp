////////////////////////////////////////////////////////////////////////////////
//
//  File: InputNek5000.cpp
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
//  Description: Nektar file format converter.
//
////////////////////////////////////////////////////////////////////////////////

#include <map>
#include <vector>
#include <sstream>
#include <string>

#include <boost/algorithm/string.hpp>
#include <LibUtilities/BasicUtils/HashUtils.hpp>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <NekMeshUtils/MeshElements/Element.h>

#include "InputNek5000.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey InputNek5000::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eInputModule, "rea5000"),
    InputNek5000::create, "Reads Nektar rea file.");

InputNek5000::InputNek5000(MeshSharedPtr m) : InputModule(m)
{
}

InputNek5000::~InputNek5000()
{
}

/**
 * @brief Processes Nek5000 file format.
 *
 * Nek5000 sessions are defined by rea files, and contain sections defining a
 * DNS simulation in a specific order. The converter only reads mesh
 * information, curve information if it exists and boundary information. The
 * format is similar to the rea format supported by #InputNek, but the layout is
 * sufficiently different that this module is separate.
 */
void InputNek5000::Process()
{
    // Open the file stream.
    OpenStream();

    string line, word;
    int nParam, nElements, nCurves;
    int i, j, k, nodeCounter = 0;
    int nComposite           = 1;
    LibUtilities::ShapeType elType;
    double vertex[8][3];

    m_mesh->m_expDim   = 0;
    m_mesh->m_spaceDim = 0;

    if (m_mesh->m_verbose)
    {
        cout << "InputNek5000: Start reading file..." << endl;
    }

    // -- Read in parameters.

    // Ignore first 3 lines. 4th line contains number of parameters.
    for (i = 0; i < 4; ++i)
    {
        getline(m_mshFile, line);
    }

    stringstream s(line);
    s >> nParam;

    for (i = 0; i < nParam; ++i)
    {
        string tmp1, tmp2;
        getline(m_mshFile, line);
        s.str(line);
        s >> tmp1 >> tmp2;
    }

    // -- Read in passive scalars (ignore)
    getline(m_mshFile, line);
    s.clear();
    s.str(line);
    s >> j;
    for (i = 0; i < j; ++i)
    {
        getline(m_mshFile, line);
    }

    // -- Read in logical switches (ignore)
    getline(m_mshFile, line);
    s.clear();
    s.str(line);
    s >> j;
    for (i = 0; i < j; ++i)
    {
        getline(m_mshFile, line);
    }

    // -- Read in mesh data.

    // First hunt for MESH tag
    bool foundMesh = false;
    while (!m_mshFile.eof())
    {
        getline(m_mshFile, line);
        if (line.find("MESH") != string::npos)
        {
            foundMesh = true;
            break;
        }
    }

    if (!foundMesh)
    {
        cerr << "Couldn't find MESH tag inside file." << endl;
        abort();
    }

    // Now read in number of elements and space dimension.
    getline(m_mshFile, line);
    s.clear();
    s.str(line);
    s >> nElements >> m_mesh->m_expDim;
    m_mesh->m_spaceDim = m_mesh->m_expDim;

    // Set up field names.
    m_mesh->m_fields.push_back("u");
    m_mesh->m_fields.push_back("v");
    if (m_mesh->m_spaceDim > 2)
    {
        m_mesh->m_fields.push_back("w");
    }
    m_mesh->m_fields.push_back("p");

    // Loop over and create elements.
    for (i = 0; i < nElements; ++i)
    {
        int nNodes;
        getline(m_mshFile, line);

        if (m_mesh->m_expDim == 2)
        {
            // - quad: 2 lines with x-coords, y-coords
            elType = LibUtilities::eQuadrilateral;
            nNodes = 4;
            for (j = 0; j < 2; ++j)
            {
                getline(m_mshFile, line);
                s.clear();
                s.str(line);
                for (k = 0; k < 4; ++k)
                {
                    s >> vertex[k][j];
                }
            }
        }
        else
        {
            // - hex: 3 lines with x/y/z-coords for base 4 nodes, then 3 more
            //   for upper 4 nodes
            elType = LibUtilities::eHexahedron;
            nNodes = 8;
            for (j = 0; j < 6; ++j)
            {
                getline(m_mshFile, line);
                s.clear();
                s.str(line);
                int offset = j > 2 ? 4 : 0;
                for (k = 0; k < 4; ++k)
                {
                    s >> vertex[offset + k][j % 3];
                }
            }
        }

        // Nek5000 meshes do not contain a unique list of nodes, so this block
        // constructs a unique set so that elements can be created with unique
        // nodes.
        vector<NodeSharedPtr> nodeList(nNodes);
        for (k = 0; k < nNodes; ++k)
        {
            nodeList[k] = std::shared_ptr<Node>(
                new Node(
                    0, vertex[k][0], vertex[k][1], vertex[k][2]));
            auto testIns = m_mesh->m_vertexSet.insert(nodeList[k]);

            if (!testIns.second)
            {
                nodeList[k] = *(testIns.first);
            }
            else
            {
                nodeList[k]->m_id = nodeCounter++;
            }
        }

        vector<int> tags(1, 0);
        ElmtConfig conf(elType, 1, false, false);
        ElementSharedPtr E = GetElementFactory().CreateInstance(
            elType, conf, nodeList, tags);
        m_mesh->m_element[E->GetDim()].push_back(E);
    }

    // -- Read in curved data.
    getline(m_mshFile, line);
    if (line.find("CURVE") == string::npos)
    {
        cerr << "Cannot find curved side data." << endl;
        abort();
    }

    // Read number of curves.
    getline(m_mshFile, line);
    s.clear();
    s.str(line);
    s >> nCurves;

    Array<OneD, NekDouble> rp;
    int nq = 6;
    LibUtilities::PointsKey curveType(nq, LibUtilities::eGaussLobattoLegendre);
    LibUtilities::PointsManager()[curveType]->GetPoints(rp);

    // Map to reorder Nek5000 -> Nektar++ edge ordering. Nek5000 has the same
    // counter-clockwise ordering of edges/vertices; however the vertical
    // (i.e. t- or xi_3-direction) edges come last.
    int nek2nekedge[12] = {
        0, 1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7
    };

    // Map to reorder Nek5000 -> Nektar++ face ordering. Again we have the same
    // counter-clockwise ordering; however the 4 vertical faces of the hex are
    // first, followed by bottom face and then top face.
    int nek2nekface[6] = {
        1, 2, 3, 4, 0, 5
    };

    if (nCurves > 0)
    {
        for (i = 0; i < nCurves; ++i)
        {
            getline(m_mshFile, line);

            int elmt, side;
            NekDouble curveData[5];
            char curveType;

            if (nElements < 1000)
            {
                // side in first 3 characters, elmt in next 3
                s.str(line.substr(0, 3));
                s >> side;
                s.clear();
                s.str(line.substr(3, 3));
                s >> elmt;
                line = line.substr(6);
            }
            else if (nElements < 1000000)
            {
                // side in first 2 characters, elmt in next 6
                s.str(line.substr(0, 2));
                s >> side;
                s.clear();
                s.str(line.substr(2, 6));
                s >> elmt;
                line = line.substr(8);
            }
            else
            {
                // side in first 2 characters, elmt in next 12
                s.str(line.substr(0, 2));
                s >> side;
                s.clear();
                s.str(line.substr(2, 12));
                s >> elmt;
                line = line.substr(14);
            }

            s.clear();
            s.str(line);

            for (j = 0; j < 5; ++j)
            {
                s >> curveData[j];
            }
            s >> curveType;

            elmt--;
            side--;
            side = nek2nekedge[side];

            switch (curveType)
            {
                case 'C':
                {
                    // Apply circular curvature to edges. Nek5000 assumes that
                    // the curvature should be imposed in x-y planes and has no
                    // z-dependence. The following code is adapted from Semtex
                    // (src/mesh.C)
                    NekDouble radius = curveData[0];
                    int convexity = radius < 0 ? -1 : 1;
                    radius = fabs(radius);

                    ElementSharedPtr el =
                        m_mesh->m_element[m_mesh->m_expDim][elmt];
                    EdgeSharedPtr edge = el->GetEdge(side);
                    edge->m_curveType = LibUtilities::eGaussLobattoLegendre;

                    // Assume 2D projection
                    Node P1(*(edge->m_n1)), P2(*(edge->m_n2));

                    if (fabs(P1.m_z - P2.m_z) > 1e-8)
                    {
                        cout << "warning: detected non x-y edge." << endl;
                    }
                    P1.m_z = P2.m_z = 0.0;

                    Node unitNormal, link, centroid, centre;
                    Node midpoint = (P1 + P2)*0.5, dx = P2 - P1;
                    NekDouble l = sqrt(dx.abs2()), sign = 0.0, semiangle = 0.0;

                    unitNormal.m_x = -dx.m_y / l;
                    unitNormal.m_y = dx.m_x / l;

                    if (2.0 * radius < l)
                    {
                        cerr << "failure" << endl;
                    }
                    else
                    {
                        semiangle = asin (0.5 * l / radius);
                    }

                    // Calculate element centroid
                    vector<NodeSharedPtr> elNodes = el->GetVertexList();
                    int nNodes = elNodes.size();

                    for (int i = 0; i < nNodes; ++i)
                    {
                        // Assume 2D projection
                        Node tmp(*elNodes[i]);
                        tmp.m_z = 0.0;
                        centroid += tmp;
                    }

                    centroid /= (NekDouble)nNodes;
                    link      = centroid - midpoint;
                    sign      = link.dot(unitNormal);
                    sign      = convexity * sign / fabs(sign);
                    centre    = midpoint + unitNormal * (sign * cos(semiangle) *
                                                         radius);

                    NekDouble theta1, theta2, dtheta, phi;
                    theta1 = atan2 (P1.m_y - centre.m_y, P1.m_x - centre.m_x);
                    theta2 = atan2 (P2.m_y - centre.m_y, P2.m_x - centre.m_x);
                    dtheta = theta2 - theta1;

                    if (fabs(dtheta) > 2.0*semiangle + 1e-15)
                    {
                        dtheta += (dtheta < 0.0) ? 2.0*M_PI : -2.0*M_PI;
                    }

                    edge->m_edgeNodes.clear();

                    for (j = 1; j < nq-1; ++j) {
                        phi = theta1 + dtheta * 0.5 * (rp[j] + 1.0);
                        NodeSharedPtr asd(new Node(
                                              0,
                                              centre.m_x + radius * cos(phi),
                                              centre.m_y + radius * sin(phi),
                                              edge->m_n1->m_z));
                        edge->m_edgeNodes.push_back(asd);
                    }
                    break;
                }
                case 's':
                case 'S':
                case 'm':
                case 'M':
                    cerr << "Curve type '" << curveType << "' on side " << side
                         << " of element " << elmt << " is unsupported;"
                         << "will ignore." << endl;
                    break;
                default:
                    cerr << "Unknown curve type '" << curveType << "' on side "
                         << side << " of element " << elmt << "; will ignore."
                         << endl;
                    break;
            }
        }
    }

    // Read boundary conditions.
    getline(m_mshFile, line);
    getline(m_mshFile, line);
    if (line.find("BOUNDARY") == string::npos)
    {
        cerr << "Cannot find boundary conditions." << endl;
        abort();
    }

    int nSurfaces = 0;
    std::unordered_set<pair<int, int>, PairHash> periodicIn;
    int periodicInId = -1, periodicOutId = -1;

    // Boundary conditions: should be precisely nElements * nFaces lines to
    // read.
    int lineCnt = 0;
    int perIn = 0, perOut = 0;

    while (m_mshFile.good())
    {
        getline(m_mshFile, line);

        // Found a new section. We don't support anything in the rea file beyond
        // this point so we'll just quit.
        if (line.find("*") != string::npos)
        {
            break;
        }

        ConditionSharedPtr c = MemoryManager<Condition>::AllocateSharedPtr();
        char bcType;
        int elmt, side;
        NekDouble data[5];

        // type in chars 0-3
        s.clear();
        s.str(line.substr(0, 4));
        s >> bcType;

        // Some lines have no boundary condition entries
        if (s.fail())
        {
            lineCnt++;
            continue;
        }

        if (nElements < 1000)
        {
            // elmt in chars 4-6, side in next 3
            s.clear();
            s.str(line.substr(4, 3));
            s >> elmt;
            s.clear();
            s.str(line.substr(7, 3));
            s >> side;
            line = line.substr(10);
        }
        else if (nElements < 100000)
        {
            // elmt in chars 4-8, side in next 1
            s.clear();
            s.str(line.substr(4, 5));
            s >> elmt;
            s.clear();
            s.str(line.substr(9, 1));
            s >> side;
            line = line.substr(10);
        }
        else if (nElements < 1000000)
        {
            // elmt in chars 4-9, no side
            s.clear();
            s.str(line.substr(4, 6));
            s >> elmt;
            side = lineCnt % (2 * m_mesh->m_expDim);
            line = line.substr(9);
        }
        else
        {
            // elmt in chars 4-15, no side
            s.clear();
            s.str(line.substr(4, 12));
            s >> elmt;
            side = lineCnt % (2 * m_mesh->m_expDim);
            line = line.substr(15);
        }

        s.clear();
        s.str(line);

        for (i = 0; i < 5; ++i)
        {
            s >> data[i];
        }

        // Our ordering starts from 0, not 1.
        --elmt;
        --side;

        // Increment lines read
        lineCnt++;

        ElementSharedPtr el = m_mesh->m_element[m_mesh->m_spaceDim][elmt];

        std::string fields[] = { "u", "v", "w", "p" };

        switch (bcType)
        {
            case 'E':
                // Edge/face connectivity; ignore since we already have this, at
                // least for conformal meshes.
                continue;

            case 'W':
            {
                for (i = 0; i < m_mesh->m_fields.size() - 1; ++i)
                {
                    c->field.push_back(fields[i]);
                    c->value.push_back("0");
                    c->type.push_back(eDirichlet);
                }

                // Set high-order boundary condition for wall.
                c->field.push_back(fields[3]);
                c->value.push_back("0");
                c->type.push_back(eHOPCondition);
                break;
            }

            case 'P':
            {
                // Determine periodic element and face.
                int perElmt = (int)(data[0] + 0.5) - 1;
                int perFace = (int)(data[1] + 0.5) - 1;

                bool setup = false;
                if (periodicInId == -1)
                {
                    periodicInId = m_mesh->m_condition.size();
                    periodicOutId = m_mesh->m_condition.size()+1;
                    setup = true;
                }

                bool hasIn = periodicIn.find(make_pair(perElmt, perFace)) !=
                    periodicIn.end();

                if (hasIn)
                {
                    swap(periodicInId, periodicOutId);
                    perOut++;
                }
                else
                {
                    periodicIn.insert(make_pair(elmt, side));
                    perIn++;
                }

                std::string periodicInStr = "[" +
                    boost::lexical_cast<string>(periodicInId) + "]";
                std::string periodicOutStr = "[" +
                    boost::lexical_cast<string>(periodicOutId) + "]";

                for (i = 0; i < m_mesh->m_fields.size() - 1; ++i)
                {
                    c->field.push_back(fields[i]);
                    c->value.push_back(periodicOutStr);
                    c->type.push_back(ePeriodic);
                }

                c->field.push_back(fields[3]);
                c->value.push_back(periodicOutStr);
                c->type.push_back(ePeriodic);

                if (setup)
                {
                    ConditionSharedPtr c2 = MemoryManager<Condition>
                        ::AllocateSharedPtr();

                    c->m_composite.push_back(nComposite++);
                    c2->m_composite.push_back(nComposite++);

                    c2->field = c->field;
                    c2->type = c->type;
                    for (i = 0; i < c->type.size(); ++i)
                    {
                        c2->value.push_back(periodicInStr);
                    }

                    m_mesh->m_condition[periodicInId] = c;
                    m_mesh->m_condition[periodicOutId] = c2;
                }

                if (hasIn)
                {
                    swap(periodicInId, periodicOutId);
                }

                break;
            }

            default:
                continue;
        }

        int compTag, conditionId;
        ElementSharedPtr surfEl;

        // Create element for face (3D) or segment (2D).
        if (el->GetDim() == 3)
        {
            FaceSharedPtr f = el->GetFace(nek2nekface[side]);
            vector<NodeSharedPtr> nodeList;
            nodeList.insert(nodeList.begin(),
                            f->m_vertexList.begin(),
                            f->m_vertexList.end());

            vector<int> tags;
            ElmtConfig conf(
                LibUtilities::eQuadrilateral, 1, true, true, false,
                LibUtilities::eGaussLobattoLegendre);
            surfEl =
                GetElementFactory().CreateInstance(LibUtilities::eQuadrilateral,
                                                   conf, nodeList, tags);

            // Copy high-order surface information from edges.
            for (int i = 0; i < f->m_vertexList.size(); ++i)
            {
                surfEl->GetEdge(i)->m_edgeNodes = f->m_edgeList[i]->m_edgeNodes;
                surfEl->GetEdge(i)->m_curveType = f->m_edgeList[i]->m_curveType;
            }
        }
        else
        {
            EdgeSharedPtr f = el->GetEdge(side);

            vector<NodeSharedPtr> nodeList;
            nodeList.push_back(f->m_n1);
            nodeList.push_back(f->m_n2);

            vector<int> tags;

            ElmtConfig conf(
                LibUtilities::eSegment, 1, true, true, false,
                LibUtilities::eGaussLobattoLegendre);
            surfEl = GetElementFactory().CreateInstance(
                LibUtilities::eSegment, conf, nodeList, tags);
        }

        // Now attempt to find this boundary condition inside
        // m_mesh->condition. This is currently a linear search and should
        // probably be made faster!
        bool found = false;
        for (auto &it : m_mesh->m_condition)
        {
            if (c == it.second)
            {
                found = true;
                c = it.second;
                break;
            }
        }

        if (!found)
        {
            conditionId = m_mesh->m_condition.size();
            compTag = nComposite;
            c->m_composite.push_back(compTag);
            m_mesh->m_condition[conditionId] = c;
        }
        else
        {
            compTag = c->m_composite[0];
        }

        // Insert composite tag into element and insert element into
        // mesh.
        vector<int> existingTags = surfEl->GetTagList();

        existingTags.insert(existingTags.begin(), compTag);
        surfEl->SetTagList(existingTags);
        surfEl->SetId(nSurfaces);

        m_mesh->m_element[surfEl->GetDim()].push_back(surfEl);
        nSurfaces++;
    }

    if (lineCnt != nElements * (m_mesh->m_expDim * 2))
    {
        cerr << "Warning: boundary conditions may not have been correctly read "
             << "from Nek5000 input file." << endl;
    }

    if (perIn != perOut)
    {
        cerr << "Warning: number of periodic faces does not match." << endl;
    }

    m_mshFile.reset();

    // -- Process rest of mesh.
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();

    // -- Set periodic composites to not be reordered.
    if (periodicInId != -1)
    {
        m_mesh->m_composite[m_mesh->m_condition[periodicInId]
                            ->m_composite[0]]->m_reorder = false;
        m_mesh->m_composite[m_mesh->m_condition[periodicOutId]
                            ->m_composite[0]]->m_reorder = false;
    }
}

}
}
