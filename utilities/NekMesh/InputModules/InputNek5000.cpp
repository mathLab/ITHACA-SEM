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
//  Description: Nektar file format converter.
//
////////////////////////////////////////////////////////////////////////////////

#include <map>
#include <vector>
#include <sstream>
#include <string>

#include <boost/algorithm/string.hpp>
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
 * @brief Processes Nektar file format.
 *
 * Nektar sessions are defined by rea files, and contain sections defining a DNS
 * simulation in a specific order. The converter only reads mesh information,
 * curve information if it exists and boundary information.
 *
 * @param pFilename Filename of Nektar session file to read.
 */
void InputNek5000::Process()
{
    // Open the file stream.
    OpenStream();

    string line, word;
    int nParam, nElements, nCurves;
    int i, j, k, nodeCounter = 0;
    int nComposite           = 0;
    LibUtilities::ShapeType elType;
    double vertex[8][3];
    map<LibUtilities::ShapeType, int> domainComposite;
    map<LibUtilities::ShapeType, vector<vector<NodeSharedPtr> > > elNodes;
    map<LibUtilities::ShapeType, vector<int> > elIds;
    boost::unordered_map<int, int> elMap;

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
            // - hex: 4 lines with x/y/z-coords for base 4 nodes, then 4 more
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

        // Nektar meshes do not contain a unique list of nodes, so this
        // block constructs a unique set so that elements can be created
        // with unique nodes.
        vector<NodeSharedPtr> nodeList(nNodes);
        for (k = 0; k < nNodes; ++k)
        {
            nodeList[k] = boost::shared_ptr<Node>(
                new Node(
                    0, vertex[k][0], vertex[k][1], vertex[k][2]));
            pair<NodeSet::iterator, bool> testIns =
                m_mesh->m_vertexSet.insert(nodeList[k]);

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

    // Reorder edges
    // Nek5000 ordering: 0-3: r-direction, 4-7: s-direction, 8-11: t-direction
    int nek2nekedge[12] = {
        //0, 2, 8, 10, 3, 1, 11, 9, 4, 5, 7, 6
        0, 1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7
    };

    if (nCurves > 0)
    {
        for (i = 0; i < nCurves; ++i)
        {
            getline(m_mshFile, line);
            s.clear();
            s.str(line);

            int elmt, side;
            NekDouble curveData[5];
            char curveType;

            s >> side >> elmt;
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
                        cout << "nope " << fabs(P1.m_z - P2.m_z) << endl;
                    }
                    P1.m_z = P2.m_z = 0.0;

                    Node unitNormal, link, centroid, centre;
                    Node midpoint = (P1 + P2)*0.5, dx = P2 - P1;
                    NekDouble l = sqrt(dx.abs2()), sign = 0.0, semiangle = 0.0;

                    unitNormal.m_x = -dx.m_y / l;
                    unitNormal.m_y = dx.m_x / l;

                    if (2.0 * radius < l)
                    {
                        cerr << "fail, midpoint = (" << midpoint.m_x << ", " << midpoint.m_y
                             << "), radius = " << radius << ", l = "
                             << l << endl;
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
                    centre    = midpoint + unitNormal * (sign * cos(semiangle) * radius);

                    NekDouble theta1, theta2, dtheta, phi;
                    theta1 = atan2 (P1.m_y - centre.m_y, P1.m_x - centre.m_x);
                    theta2 = atan2 (P2.m_y - centre.m_y, P2.m_x - centre.m_x);
                    dtheta = theta2 - theta1;

                    if (fabs(dtheta) > 2.0*semiangle + 1e-15)
                    {
                        dtheta += (dtheta < 0.0) ? 2.0*M_PI : -2.0*M_PI;
                    }

                    edge->m_edgeNodes.clear();

                    //cout << edge->m_n1->m_x << " " << edge->m_n1->m_y << endl;

                    for (j = 1; j < nq-1; ++j) {
                        phi = theta1 + dtheta * 0.5 * (rp[j] + 1.0);
                        NodeSharedPtr asd(new Node(
                                              0,
                                              centre.m_x + radius * cos(phi),
                                              centre.m_y + radius * sin(phi),
                                              edge->m_n1->m_z));
                        edge->m_edgeNodes.push_back(asd);
                        //cout << asd->m_x << " " << asd->m_y << endl;
                    }

                    //cout << edge->m_n2->m_x << " " << edge->m_n2->m_y << endl;
                    break;
                }
                case 's':
                case 'S':
                case 'm':
                case 'M':
                    cerr << "Curve type '" << curveType << "' is unsupported."
                         << endl;
                    break;
                default:
                    cerr << "Unknown curve type '" << curveType << "'" << endl;
                    break;
            }
        }
    }

    m_mshFile.close();

    // -- Process rest of mesh.
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();
}

}
}
