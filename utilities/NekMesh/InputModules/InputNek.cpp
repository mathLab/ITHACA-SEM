////////////////////////////////////////////////////////////////////////////////
//
//  File: InputNek.cpp
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
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <NekMeshUtils/MeshElements/Element.h>
#include <NekMeshUtils/MeshElements/Prism.h>
#include <NekMeshUtils/MeshElements/Tetrahedron.h>

#include "InputNek.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey InputNek::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eInputModule, "rea"), InputNek::create, "Reads Nektar rea file.");

InputNek::InputNek(MeshSharedPtr m) : InputModule(m)
{
    m_config["scalar"] = ConfigOption(
        true, "0", "If defined then assume input rea is for scalar problem");
}

InputNek::~InputNek()
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
void InputNek::Process()
{
    // Open the file stream.
    OpenStream();

    string line, word;
    int nParam, nElements, nCurves;
    int i, j, k, nodeCounter = 0;
    int nComposite           = 0;
    LibUtilities::ShapeType elType;
    double vertex[3][8];
    map<LibUtilities::ShapeType, int> domainComposite;
    map<LibUtilities::ShapeType, vector<vector<NodeSharedPtr> > > elNodes;
    map<LibUtilities::ShapeType, vector<int> > elIds;
    std::unordered_map<int, int> elMap;
    vector<LibUtilities::ShapeType> elmOrder;

    bool scalar = m_config["scalar"].as<bool>();

    // Set up vector of processing orders.
    elmOrder.push_back(LibUtilities::eSegment);
    elmOrder.push_back(LibUtilities::eTriangle);
    elmOrder.push_back(LibUtilities::eQuadrilateral);
    elmOrder.push_back(LibUtilities::ePrism);
    elmOrder.push_back(LibUtilities::ePyramid);
    elmOrder.push_back(LibUtilities::eTetrahedron);
    elmOrder.push_back(LibUtilities::eHexahedron);

    m_mesh->m_expDim   = 0;
    m_mesh->m_spaceDim = 0;

    if (m_mesh->m_verbose)
    {
        cout << "InputNek: Start reading file..." << endl;
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

    if (!scalar)
    {
        m_mesh->m_fields.push_back("v");
        if (m_mesh->m_spaceDim > 2)
        {
            m_mesh->m_fields.push_back("w");
        }
        m_mesh->m_fields.push_back("p");
    }

    // Loop over and create elements.
    for (i = 0; i < nElements; ++i)
    {
        getline(m_mshFile, line);

        if (m_mesh->m_expDim == 2)
        {
            if (line.find("Qua") != string::npos ||
                line.find("qua") != string::npos)
            {
                elType = LibUtilities::eQuadrilateral;
            }
            else
            {
                // Default element type in 2D is triangle.
                elType = LibUtilities::eTriangle;
            }
        }
        else
        {
            if (line.find("Tet") != string::npos ||
                line.find("tet") != string::npos)
            {
                elType = LibUtilities::eTetrahedron;
            }
            else if (line.find("Hex") != string::npos ||
                     line.find("hex") != string::npos)
            {
                elType = LibUtilities::eHexahedron;
            }
            else if (line.find("Prism") != string::npos ||
                     line.find("prism") != string::npos)
            {
                elType = LibUtilities::ePrism;
            }
            else if (line.find("Pyr") != string::npos ||
                     line.find("pyr") != string::npos)
            {
                elType = LibUtilities::ePyramid;
            }
            else if (line.find("Qua") != string::npos ||
                     line.find("qua") != string::npos)
            {
                elType = LibUtilities::eQuadrilateral;
            }
            else
            {
                // Default element type in 2D is tetrahedron.
                elType = LibUtilities::eTetrahedron;
            }
        }

        // Read in number of vertices for element type.
        const int nNodes = GetNnodes(elType);

        for (j = 0; j < m_mesh->m_expDim; ++j)
        {
            getline(m_mshFile, line);
            s.clear();
            s.str(line);
            for (k = 0; k < nNodes; ++k)
            {
                s >> vertex[j][k];
            }
        }

        // Zero co-ordinates bigger than expansion dimension.
        for (j = m_mesh->m_expDim; j < 3; ++j)
        {
            for (k = 0; k < nNodes; ++k)
            {
                vertex[j][k] = 0.0;
            }
        }

        // Nektar meshes do not contain a unique list of nodes, so this
        // block constructs a unique set so that elements can be created
        // with unique nodes.
        vector<NodeSharedPtr> nodeList;
        for (k = 0; k < nNodes; ++k)
        {
            NodeSharedPtr n = std::shared_ptr<Node>(new Node(
                nodeCounter++, vertex[0][k], vertex[1][k], vertex[2][k]));
            nodeList.push_back(n);
        }

        elNodes[elType].push_back(nodeList);
        elIds[elType].push_back(i);
    }

    int reorderedId = 0;
    nodeCounter     = 0;

    for (i = 0; i < elmOrder.size(); ++i)
    {
        LibUtilities::ShapeType elType      = elmOrder[i];
        vector<vector<NodeSharedPtr> > &tmp = elNodes[elType];

        for (j = 0; j < tmp.size(); ++j)
        {
            vector<int> tags;
            auto compIt = domainComposite.find(elType);
            if (compIt == domainComposite.end())
            {
                tags.push_back(nComposite);
                domainComposite[elType] = nComposite;
                nComposite++;
            }
            else
            {
                tags.push_back(compIt->second);
            }

            elMap[elIds[elType][j]] = reorderedId++;

            vector<NodeSharedPtr> nodeList = tmp[j];

            for (k = 0; k < nodeList.size(); ++k)
            {
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

            // Create linear element
            ElmtConfig conf(elType, 1, false, false);
            ElementSharedPtr E = GetElementFactory().CreateInstance(
                elType, conf, nodeList, tags);
            m_mesh->m_element[E->GetDim()].push_back(E);
        }
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

    if (nCurves > 0)
    {
        string curveTag;

        for (i = 0; i < nCurves; ++i)
        {
            getline(m_mshFile, line);
            s.clear();
            s.str(line);
            s >> word;

            if (word == "File")
            {
                // Next line contains filename and curve tag.
                getline(m_mshFile, line);
                s.clear();
                s.str(line);
                s >> word >> curveTag;
                curveTags[curveTag] = make_pair(eFile, word);
            }
            else if (word == "Recon")
            {
                // Next line contains curve tag.
                getline(m_mshFile, line);
                s.clear();
                s.str(line);
                s >> word >> curveTag;
                curveTags[curveTag] = make_pair(eRecon, word);
            }
            else
            {
                cerr << "Unsupported curve type " << word << endl;
                abort();
            }
        }

        // Load high order surface information.
        LoadHOSurfaces();

        // Read in curve information. First line should contain number
        // of curved sides.
        getline(m_mshFile, line);

        if (line.find("side") == string::npos)
        {
            cerr << "Unable to read number of curved sides" << endl;
            abort();
        }

        int nCurvedSides;
        int faceId, elId;

        s.clear();
        s.str(line);
        s >> nCurvedSides;

        // Iterate over curved sides, and look up high-order surface
        // information in the HOSurfSet, then map this onto faces.
        for (i = 0; i < nCurvedSides; ++i)
        {
            getline(m_mshFile, line);
            s.clear();
            s.str(line);
            s >> faceId >> elId >> word;
            faceId--;
            elId                = elMap[elId - 1];
            ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][elId];

            int origFaceId = faceId;

            if (el->GetConf().m_e == LibUtilities::ePrism && faceId % 2 == 0)
            {
                std::shared_ptr<Prism> p =
                    std::dynamic_pointer_cast<Prism>(el);
                if (p->m_orientation == 1)
                {
                    faceId = (faceId + 2) % 6;
                }
                else if (p->m_orientation == 2)
                {
                    faceId = (faceId + 4) % 6;
                }
            }
            else if (el->GetConf().m_e == LibUtilities::eTetrahedron)
            {
                std::shared_ptr<Tetrahedron> t =
                    std::dynamic_pointer_cast<Tetrahedron>(el);
                faceId = t->m_orientationMap[faceId];
            }

            auto it = curveTags.find(word);
            if (it == curveTags.end())
            {
                cerr << "Unrecognised curve tag " << word << " in curved lines"
                     << endl;
                abort();
            }

            if (it->second.first == eRecon)
            {
                // Spherigon information: read in vertex normals.
                vector<NodeSharedPtr> &tmp = el->GetFace(faceId)->m_vertexList;
                vector<Node> n(tmp.size());

                int offset = 0;
                if (el->GetConf().m_e == LibUtilities::ePrism &&
                    faceId % 2 == 1)
                {
                    offset =
                        std::dynamic_pointer_cast<Prism>(el)->m_orientation;
                }

                // Read x/y/z coordinates.
                getline(m_mshFile, line);
                s.clear();
                s.str(line);
                for (j = 0; j < tmp.size(); ++j)
                {
                    s >> n[j].m_x;
                }

                getline(m_mshFile, line);
                s.clear();
                s.str(line);
                for (j = 0; j < tmp.size(); ++j)
                {
                    s >> n[j].m_y;
                }

                getline(m_mshFile, line);
                s.clear();
                s.str(line);
                for (j = 0; j < tmp.size(); ++j)
                {
                    s >> n[j].m_z;
                }

                for (j = 0; j < tmp.size(); ++j)
                {
                    int id = tmp[(j + offset) % tmp.size()]->m_id;
                    auto vIt = m_mesh->m_vertexNormals.find(id);

                    if (vIt == m_mesh->m_vertexNormals.end())
                    {
                        m_mesh->m_vertexNormals[id] = n[j];
                    }
                }

                // Add edge/face to list of faces to apply spherigons
                // to.
                m_mesh->m_spherigonSurfs.insert(make_pair(elId, faceId));
            }
            else if (it->second.first == eFile)
            {
                FaceSharedPtr f               = el->GetFace(faceId);
                static int tetFaceVerts[4][3] = {
                    {0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}};

                vector<int> vertId(3);
                s >> vertId[0] >> vertId[1] >> vertId[2];

                // Prisms and tets may have been rotated by OrientPrism
                // routine which reorders vertices. This block rotates
                // vertex ids accordingly.
                if (el->GetConf().m_e == LibUtilities::eTetrahedron)
                {
                    std::shared_ptr<Tetrahedron> tet =
                        std::static_pointer_cast<Tetrahedron>(el);
                    vector<int> tmpVertId = vertId;

                    for (j = 0; j < 3; ++j)
                    {
                        int v =
                            tet->GetVertex(tet->m_origVertMap
                                               [tetFaceVerts[origFaceId][j]])
                                ->m_id;

                        for (k = 0; k < 3; ++k)
                        {
                            int w = f->m_vertexList[k]->m_id;
                            if (v == w)
                            {
                                vertId[k] = tmpVertId[j];
                                break;
                            }
                        }
                    }
                }
                else if (el->GetConf().m_e == LibUtilities::ePrism)
                {
                    std::shared_ptr<Prism> pr =
                        std::static_pointer_cast<Prism>(el);
                    if (pr->m_orientation == 1)
                    {
                        swap(vertId[2], vertId[1]);
                        swap(vertId[0], vertId[1]);
                    }
                    else if (pr->m_orientation == 2)
                    {
                        swap(vertId[0], vertId[1]);
                        swap(vertId[2], vertId[1]);
                    }
                }

                HOSurfSharedPtr hs =
                    std::shared_ptr<HOSurf>(new HOSurf(vertId));
                // Find vertex combination in hoData.
                auto hoIt = hoData[word].find(hs);

                if (hoIt == hoData[word].end())
                {
                    cerr << "Unable to find high-order surface data "
                         << "for element id " << elId + 1 << endl;
                    abort();
                }

                // Depending on order of vertices in rea file, surface
                // information may need to be rotated or reflected.
                HOSurfSharedPtr surf = *hoIt;
                surf->Align(vertId);

                // Finally, add high order data to appropriate face.
                int Ntot = (*hoIt)->surfVerts.size();
                int N    = ((int)sqrt(8.0 * Ntot + 1.0) - 1) / 2;
                EdgeSharedPtr edge;

                // Apply high-order map to convert face data to Nektar++
                // ordering (vertices->edges->internal).
                vector<NodeSharedPtr> tmpVerts = (*hoIt)->surfVerts;
                for (j = 0; j < tmpVerts.size(); ++j)
                {
                    (*hoIt)->surfVerts[hoMap[j]] = tmpVerts[j];
                }

                for (j = 0; j < tmpVerts.size(); ++j)
                {
                    NodeSharedPtr a = (*hoIt)->surfVerts[j];
                }

                vector<int> faceVertIds(3);
                faceVertIds[0] = f->m_vertexList[0]->m_id;
                faceVertIds[1] = f->m_vertexList[1]->m_id;
                faceVertIds[2] = f->m_vertexList[2]->m_id;

                for (j = 0; j < f->m_edgeList.size(); ++j)
                {
                    edge = f->m_edgeList[j];

                    // Skip over edges which have already been
                    // populated.
                    if (edge->m_edgeNodes.size() > 0)
                    {
                        continue;
                    }

                    // edge->m_edgeNodes.clear();
                    edge->m_curveType = LibUtilities::eGaussLobattoLegendre;

                    for (int k = 0; k < N - 2; ++k)
                    {
                        edge->m_edgeNodes.push_back(
                            (*hoIt)->surfVerts[3 + j * (N - 2) + k]);
                    }

                    // Nodal triangle data is always
                    // counter-clockwise. Add this check to reorder
                    // where necessary.
                    if (edge->m_n1->m_id != faceVertIds[j])
                    {
                        reverse(edge->m_edgeNodes.begin(),
                                edge->m_edgeNodes.end());
                    }
                }

                // Insert interior face curvature.
                f->m_curveType = LibUtilities::eNodalTriElec;
                for (int j = 3 + 3 * (N - 2); j < Ntot; ++j)
                {
                    f->m_faceNodes.push_back((*hoIt)->surfVerts[j]);
                }
            }
        }
    }

    // -- Process fluid boundary conditions.

    // Define a fairly horrendous map: key is the condition ID, the
    // value is a vector of pairs of composites and element
    // types. Essentially this map takes conditions -> composites for
    // each element type.
    map<int, vector<pair<int, LibUtilities::ShapeType> > > surfaceCompMap;

    // Skip boundary conditions line.
    getline(m_mshFile, line);
    getline(m_mshFile, line);

    int nSurfaces = 0;

    while (true)
    {
        getline(m_mshFile, line);

        // Break out of loop at end of boundary conditions section.
        if (line.find("*") != string::npos || m_mshFile.eof() ||
            line.length() == 0)
        {
            break;
        }

        // Read boundary type, element ID and face ID.
        char bcType;
        int elId, faceId;
        s.clear();
        s.str(line);
        s >> bcType >> elId >> faceId;
        faceId--;
        elId = elMap[elId - 1];

        vector<string> vals;
        vector<ConditionType> type;
        ConditionSharedPtr c = MemoryManager<Condition>::AllocateSharedPtr();

        ElementSharedPtr elm = m_mesh->m_element[m_mesh->m_spaceDim][elId];

        // Ignore BCs for undefined edges/faces
        if ((elm->GetDim() == 2 && faceId >= elm->GetEdgeCount()) ||
            (elm->GetDim() == 3 && faceId >= elm->GetFaceCount()))
        {
            continue;
        }

        // First character on each line describes type of BC. Currently
        // only support V, W, and O. In this switch statement we
        // construct the quantities needed to search for the condition.
        switch (bcType)
        {
            // Wall boundary.
            case 'W':
            {
                if (scalar)
                {
                    vals.push_back("0");
                    type.push_back(eDirichlet);
                }
                else
                {
                    for (i = 0; i < m_mesh->m_fields.size() - 1; ++i)
                    {
                        vals.push_back("0");
                        type.push_back(eDirichlet);
                    }
                    // Set high-order boundary condition for wall.
                    vals.push_back("0");
                    type.push_back(eHOPCondition);
                }
                break;
            }

            // Velocity boundary condition (either constant or dependent
            // upon x,y,z).
            case 'V':
            case 'v':
            {
                if (scalar)
                {
                    getline(m_mshFile, line);
                    size_t p = line.find_first_of('=');
                    vals.push_back(
                        boost::algorithm::trim_copy(line.substr(p + 1)));
                    type.push_back(eDirichlet);
                }
                else
                {
                    for (i = 0; i < m_mesh->m_fields.size() - 1; ++i)
                    {
                        getline(m_mshFile, line);
                        size_t p = line.find_first_of('=');
                        vals.push_back(
                            boost::algorithm::trim_copy(line.substr(p + 1)));
                        type.push_back(eDirichlet);
                    }
                    // Set high-order boundary condition for Dirichlet
                    // condition.
                    vals.push_back("0");
                    type.push_back(eHOPCondition);
                }
                break;
            }

            // Natural outflow condition (default value = 0.0?)
            case 'O':
            {
                if (scalar)
                {
                    vals.push_back("0");
                    type.push_back(eNeumann);
                }
                else
                {
                    for (i = 0; i < m_mesh->m_fields.size(); ++i)
                    {
                        vals.push_back("0");
                        type.push_back(eNeumann);
                    }
                    // Set zero Dirichlet condition for outflow.
                    type[m_mesh->m_fields.size() - 1] = eDirichlet;
                }
                break;
            }

            // Ignore unsupported BCs
            case 'P':
            case 'E':
            case 'F':
            case 'f':
            case 'N':
                continue;
                break;

            default:
                cerr << "Unknown boundary condition type " << line[0] << endl;
                abort();
        }

        // Populate condition information.
        c->field = m_mesh->m_fields;
        c->type  = type;
        c->value = vals;

        // Now attempt to find this boundary condition inside
        // m_mesh->condition. This is currently a linear search and should
        // probably be made faster!
        bool found = false;
        auto it = m_mesh->m_condition.begin();
        for (; it != m_mesh->m_condition.end(); ++it)
        {
            if (c == it->second)
            {
                found = true;
                break;
            }
        }

        int compTag, conditionId;
        ElementSharedPtr surfEl;

        // Create element for face (3D) or segment (2D). At the moment
        // this is a bit of a hack since high-order nodes are not
        // copied, so some output modules (e.g. Gmsh) will not output
        // correctly.
        if (elm->GetDim() == 3)
        {
            // 3D elements may have been reoriented, so face IDs will
            // change.
            if (elm->GetConf().m_e == LibUtilities::ePrism && faceId % 2 == 0)
            {
                std::shared_ptr<Prism> p =
                    std::dynamic_pointer_cast<Prism>(elm);
                if (p->m_orientation == 1)
                {
                    faceId = (faceId + 2) % 6;
                }
                else if (p->m_orientation == 2)
                {
                    faceId = (faceId + 4) % 6;
                }
            }
            else if (elm->GetConf().m_e == LibUtilities::eTetrahedron)
            {
                std::shared_ptr<Tetrahedron> t =
                    std::dynamic_pointer_cast<Tetrahedron>(elm);
                faceId = t->m_orientationMap[faceId];
            }

            FaceSharedPtr f = elm->GetFace(faceId);
            bool tri        = f->m_vertexList.size() == 3;

            vector<NodeSharedPtr> nodeList;
            nodeList.insert(nodeList.begin(),
                            f->m_vertexList.begin(),
                            f->m_vertexList.end());

            vector<int> tags;

            LibUtilities::ShapeType seg =
                tri ? LibUtilities::eTriangle : LibUtilities::eQuadrilateral;
            ElmtConfig conf(
                seg, 1, true, true, false, LibUtilities::eGaussLobattoLegendre);
            surfEl =
                GetElementFactory().CreateInstance(seg, conf, nodeList, tags);

            // Copy high-order surface information from edges.
            for (int i = 0; i < f->m_vertexList.size(); ++i)
            {
                surfEl->GetEdge(i)->m_edgeNodes = f->m_edgeList[i]->m_edgeNodes;
                surfEl->GetEdge(i)->m_curveType = f->m_edgeList[i]->m_curveType;
            }
        }
        else if (faceId < elm->GetEdgeCount())
        {
            EdgeSharedPtr f = elm->GetEdge(faceId);

            vector<NodeSharedPtr> nodeList;
            nodeList.push_back(f->m_n1);
            nodeList.push_back(f->m_n2);

            vector<int> tags;

            ElmtConfig conf(LibUtilities::eSegment,
                            1,
                            true,
                            true,
                            false,
                            LibUtilities::eGaussLobattoLegendre);
            surfEl = GetElementFactory().CreateInstance(
                LibUtilities::eSegment, conf, nodeList, tags);
        }

        if (!surfEl)
        {
            continue;
        }

        LibUtilities::ShapeType surfElType = surfEl->GetConf().m_e;

        if (!found)
        {
            // If condition does not already exist, add to condition
            // list, create new composite tag and put inside
            // surfaceCompMap.
            conditionId = m_mesh->m_condition.size();
            compTag = nComposite;
            c->m_composite.push_back(compTag);
            m_mesh->m_condition[conditionId] = c;

            surfaceCompMap[conditionId].push_back(
                pair<int, LibUtilities::ShapeType>(nComposite, surfElType));

            nComposite++;
        }
        else
        {
            // Otherwise find existing composite inside surfaceCompMap.
            auto it2 = surfaceCompMap.find(it->first);

            found = false;
            if (it2 == surfaceCompMap.end())
            {
                // This should never happen!
                cerr << "Unable to find condition!" << endl;
                abort();
            }

            for (j = 0; j < it2->second.size(); ++j)
            {
                pair<int, LibUtilities::ShapeType> tmp = it2->second[j];
                if (tmp.second == surfElType)
                {
                    found   = true;
                    compTag = tmp.first;
                    break;
                }
            }

            // If no pairs were found, then this condition contains
            // multiple element types (i.e. both triangles and
            // quads). Create another composite for the new shape type
            // and insert into the map.
            if (!found)
            {
                it2->second.push_back(
                    pair<int, LibUtilities::ShapeType>(nComposite, surfElType));
                compTag = nComposite;
                m_mesh->m_condition[it->first]->m_composite.push_back(compTag);
                nComposite++;
            }

            conditionId = it->first;
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

    m_mshFile.reset();

    // -- Process rest of mesh.
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();
}

/**
 * Load high order surface information from hsf file.
 */
void InputNek::LoadHOSurfaces()
{
    int nodeId = m_mesh->GetNumEntities();

    for (auto &it : curveTags)
    {
        ifstream hsf;
        string line, fileName = it.second.second;
        size_t pos;
        int N, Nface, dot;

        if (it.second.first != eFile)
        {
            continue;
        }

        // Replace fro extension with hsf.
        dot      = fileName.find_last_of('.');
        fileName = fileName.substr(0, dot);
        fileName += ".hsf";

        // Open hsf file.
        hsf.open(fileName.c_str());
        if (!hsf.is_open())
        {
            cerr << "Could not open surface file " << fileName << endl;
            abort();
        }

        // Read in header line; determine element order, number of faces
        // from this line.
        getline(hsf, line);
        pos = line.find("=");
        if (pos == string::npos)
        {
            cerr << "hsf header error: cannot read number of "
                 << "nodal points." << endl;
            abort();
        }
        line = line.substr(pos + 1);
        stringstream ss(line);
        ss >> N;

        pos = line.find("=");
        if (pos == string::npos)
        {
            cerr << "hsf header error: cannot read number of "
                 << "faces." << endl;
            abort();
        }
        line = line.substr(pos + 1);
        ss.clear();
        ss.str(line);
        ss >> Nface;

        int Ntot = N * (N + 1) / 2;

        // Skip a line, then read in r,s positions inside the next
        // comments.
        Array<OneD, NekDouble> r(Ntot), s(Ntot);
        getline(hsf, line);

        for (int i = 0; i < 2; ++i)
        {
            string word;

            getline(hsf, line);
            ss.clear();
            ss.str(line);
            ss >> word;

            if (word != "#")
            {
                cerr << "hsf header error: cannot read in "
                     << "r/s points" << endl;
                abort();
            }

            for (int j = 0; j < Ntot; ++j)
            {
                ss >> (i == 0 ? r[j] : s[j]);
            }
        }

        // Generate electrostatic points so that re-mapping array can
        // be constructed.
        Array<OneD, NekDouble> rp(Ntot), sp(Ntot);
        LibUtilities::PointsKey elec(N, LibUtilities::eNodalTriElec);
        LibUtilities::PointsManager()[elec]->GetPoints(rp, sp);

        // Expensively construct remapping array nodemap. This will
        // map nodal ordering from hsf order to Nektar++ ordering
        // (i.e. vertices followed by edges followed by interior
        // points.)
        for (int i = 0; i < Ntot; ++i)
        {
            for (int j = 0; j < Ntot; ++j)
            {
                if (fabs(r[i] - rp[j]) < 1e-5 && fabs(s[i] - sp[j]) < 1e-5)
                {
                    hoMap[i] = j;
                    break;
                }
            }
        }

        // Skip variables line
        getline(hsf, line);

        // Read in nodal points for each face.
        map<int, vector<NodeSharedPtr> > faceMap;
        for (int i = 0; i < Nface; ++i)
        {
            getline(hsf, line);
            vector<NodeSharedPtr> faceNodes(Ntot);
            for (int j = 0; j < Ntot; ++j, ++nodeId)
            {
                double x, y, z;
                getline(hsf, line);
                ss.clear();
                ss.str(line);
                ss >> x >> y >> z;
                faceNodes[j] = NodeSharedPtr(new Node(nodeId, x, y, z));
            }
            // Skip over tecplot connectivity information.
            for (int j = 0; j < (N - 1) * (N - 1); ++j)
            {
                getline(hsf, line);
            }
            faceMap[i] = faceNodes;
        }

        // Finally, read in connectivity information to set up after
        // reading rea file.
        getline(hsf, line);
        for (int i = 0; i < Nface; ++i)
        {
            string tmp;
            int fid;
            vector<int> nodeIds(3);

            getline(hsf, line);
            ss.clear();
            ss.str(line);
            ss >> tmp >> fid >> nodeIds[0] >> nodeIds[1] >> nodeIds[2];

            if (tmp != "#")
            {
                cerr << "Unable to read hsf connectivity information." << endl;
                abort();
            }

            hoData[it.first].insert(
                HOSurfSharedPtr(new HOSurf(nodeIds, faceMap[i])));
        }

        hsf.close();
    }
}

/**
 * This routine aids the reading of Nektar files only; it returns the
 * number of nodes for a given entity typw.
 */
int InputNek::GetNnodes(LibUtilities::ShapeType InputNekEntity)
{
    int nNodes = 0;

    switch (InputNekEntity)
    {
        case LibUtilities::ePoint:
            nNodes = 1;
            break;
        case LibUtilities::eSegment:
            nNodes = 2;
            break;
        case LibUtilities::eTriangle:
            nNodes = 3;
            break;
        case LibUtilities::eQuadrilateral:
            nNodes = 4;
            break;
        case LibUtilities::eTetrahedron:
            nNodes = 4;
            break;
        case LibUtilities::ePyramid:
            nNodes = 5;
            break;
        case LibUtilities::ePrism:
            nNodes = 6;
            break;
        case LibUtilities::eHexahedron:
            nNodes = 8;
            break;
        default:
            cerr << "unknown Nektar element type" << endl;
    }

    return nNodes;
}
}
}
