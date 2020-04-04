////////////////////////////////////////////////////////////////////////////////
//
//  File: InputStarTec.cpp
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
//  Description: Tecplot file converter.
//
////////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>
#include <boost/algorithm/string.hpp>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <NekMeshUtils/MeshElements/Element.h>

#include "InputStarTec.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey InputTec::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eInputModule, "dat"),
    InputTec::create,
    "Reads Tecplot polyhedron ascii format converted from Star CCM (.dat).");

InputTec::InputTec(MeshSharedPtr m) : InputModule(m)
{
}

InputTec::~InputTec()
{
}

/**
 * Tecplot file Polyhedron format contains a list of nodes, a
 * node count per face, the node ids, Element ids that are on
 * the left of each face and Element ids which are on the
 * right of each face. There are then a series of zone of each
 * surface. In the case of a surface the number of nodes is
 * not provided indicating it is a 2D zone
 *
 * @param pFilename Filename of Tecplot file to read.
 */
void InputTec::Process()
{
    m_mesh->m_expDim   = 3;
    m_mesh->m_spaceDim = 3;

    if (m_mesh->m_verbose)
    {
        cout << "InputStarTec: Start reading file..." << endl;
    }

    string line, word;

    // Open the file stream.
    OpenStream();

    int nComposite = 0;

    // read first zone (Hopefully 3D)
    while (!m_mshFile.eof())
    {
        getline(m_mshFile, line);
        if (line.find("ZONE") != string::npos)
        {
            ReadZone(nComposite);
            break;
        }
    }

    // read remaining 2D zones
    while (!m_mshFile.eof())
    {
        if (line.find("ZONE") != string::npos)
        {
            ReadZone(nComposite);
        }
    }

    PrintSummary();
    m_mshFile.reset();

    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();
}

void InputTec::ReadZone(int &nComposite)
{
    int i;
    string line, tag;
    int nfaces, nnodes, nelements;
    int start, end;
    stringstream s;
    NekDouble value;
    static int zcnt = 1;

    // Read Zone Header
    nnodes = nfaces = nelements = 0;
    while (!m_mshFile.eof())
    {
        getline(m_mshFile, line);

        boost::to_upper(line);

        // cehck to see if readable data.
        if (sscanf(line.c_str(), "%lf", &value) == 1)
        {
            break;
        }

        if ((line.find("NODES") != string::npos) &&
            (line.find("TOTALNUMFACENODES") == string::npos))
        {
            s.clear();
            s.str(line);

            tag    = s.str();
            start  = tag.find("NODES=");
            end    = tag.find_first_of(',', start);
            nnodes = atoi(tag.substr(start + 6, end).c_str());
        }

        if ((line.find("FACES") != string::npos) &&
            (line.find("NUMCONNECTEDBOUNDARYFACES") == string::npos))
        {
            s.clear();
            s.str(line);

            tag    = s.str();
            start  = tag.find("FACES=");
            end    = tag.find_first_of(',', start);
            nfaces = atoi(tag.substr(start + 6, end).c_str());
        }

        if (line.find("ELEMENTS") != string::npos)
        {
            s.clear();
            s.str(line);

            tag       = s.str();
            start     = tag.find("ELEMENTS=");
            end       = tag.find_first_of(',', start);
            nelements = atoi(tag.substr(start + 9, end).c_str());
        }

        if (line.find("ZONETYPE") != string::npos)
        {
            s.clear();
            s.str(line);

            if ((line.find("FEPOLYGON") == string::npos) &&
                (line.find("FEPOLYHEDRON") == string::npos))
            {
                ASSERTL0(false,
                         "Routine only set up for FEPolygon or FEPolyhedron");
            }
        }
    }
    if (!nnodes) // No zone found
    {
        return;
    }

    cout << "Setting up zone " << zcnt++;

    int nodeCount = 3 * nnodes;
    vector<NekDouble> nodeLocs;

    while (nodeCount > 0 && !m_mshFile.eof())
    {
        s.clear();
        s.str(line);
        while (s >> value)
        {
            nodeLocs.push_back(value);
            nodeCount--;
        }
        if (nodeCount > 0)
        {
            getline(m_mshFile, line);
        }
    }

    ASSERTL0(nodeLocs.size() == 3*nnodes, "Unable to read correct number of "
             "nodes from Tecplot file");

    std::vector<NodeSharedPtr> Nodes;
    for (i = 0; i < nnodes; ++i)
    {
        Nodes.push_back(
            std::shared_ptr<Node>(
                new Node(i, nodeLocs[i], nodeLocs[i+nnodes],
                         nodeLocs[i+2*nnodes])));
    }

    // Read Node count per face
    getline(m_mshFile, line);
    if (line.find("node count per face") == string::npos)
    {
        if (line.find("face nodes") == string::npos)
        {
            getline(m_mshFile, line);
        }
    }

    s.clear();
    s.str(line);

    vector<int> Nodes_per_face;
    if (line.find("node count per face") != string::npos)
    {
        int nodes;
        for (i = 0; i < nfaces; ++i)
        {
            m_mshFile >> nodes;
            ASSERTL0(nodes <= 4,
                     "Can only handle meshes with "
                     "up to four nodes per face");
            Nodes_per_face.push_back(nodes);
        }
        // Read next line
        getline(m_mshFile, line);
    }

    // Read face nodes;
    if (line.find("face nodes") == string::npos)
    {
        getline(m_mshFile, line);
    }
    s.clear();
    s.str(line);

    vector<vector<int> > FaceNodes;

    if (line.find("face nodes") != string::npos)
    {

        for (i = 0; i < nfaces; ++i)
        {
            // check to see if Nodes_per_face is defined and
            // if not assume 2 nodes for 2D case
            int nodes = (Nodes_per_face.size()) ? Nodes_per_face[i] : 2;

            int nodeID;
            vector<int> Fnodes;
            for (int j = 0; j < nodes; ++j)
            {

                m_mshFile >> nodeID;

                Fnodes.push_back(nodeID - 1);
            }

            FaceNodes.push_back(Fnodes);
        }
    }
    else
    {
        ASSERTL0(false, "Failed to find face node section");
    }

    // Read left elements
    Array<OneD, vector<int> > ElementFaces(nelements);

    // check to see if next line contains left elements
    getline(m_mshFile, line);
    if (line.find("left elements") == string::npos)
    {
        getline(m_mshFile, line);
    }

    if (line.find("left elements") != string::npos)
    {
        int elmtID;

        for (i = 0; i < nfaces; ++i)
        {
            m_mshFile >> elmtID;

            if (elmtID > 0)
            {
                ElementFaces[elmtID - 1].push_back(i);
            }
        }
    }
    else
    {
        ASSERTL0(false, "Left element not found");
    }

    // check to see if next line contains right elements
    getline(m_mshFile, line);
    if (line.find("right elements") == string::npos)
    {
        getline(m_mshFile, line);
    }

    if (line.find("right elements") != string::npos)

    {
        int elmtID;

        for (i = 0; i < nfaces; ++i)
        {
            m_mshFile >> elmtID;

            if (elmtID > 0)
            {
                ElementFaces[elmtID - 1].push_back(i);
            }
        }

        // read to end of line
        getline(m_mshFile, line);
    }
    else
    {
        ASSERTL0(false, "Left element not found");
    }

    if (Nodes_per_face.size()) // 3D Zone
    {
        cout << " (3D) " << endl;

        // Reset node ordering so that all prism faces have
        // consistent numbering for singular vertex re-ordering
        ResetNodes(Nodes, ElementFaces, FaceNodes);

        m_mesh->m_node = Nodes;

        // create Prisms/Pyramids first
        for (i = 0; i < nelements; ++i)
        {
            if (ElementFaces[i].size() > 4)
            {
                GenElement3D(
                    Nodes, i, ElementFaces[i], FaceNodes, nComposite, true);
            }
        }

        nComposite++;

        // create Tets second
        for (i = 0; i < nelements; ++i)
        {
            if (ElementFaces[i].size() == 4)
            {
                GenElement3D(
                    Nodes, i, ElementFaces[i], FaceNodes, nComposite, true);
            }
        }
        nComposite++;

        ProcessVertices();
    }
    else // 2D Zone
    {
        cout << " (2D)" << endl;

        // find ids of VertNodes from m_mesh->m_vertexSet so that we can
        // identify
        for (i = 0; i < Nodes.size(); ++i)
        {
            auto it = m_mesh->m_vertexSet.find(Nodes[i]);

            if (it == m_mesh->m_vertexSet.end())
            {
                ASSERTL0(false, "Failed to find face vertex in 3D list");
            }
            else
            {
                Nodes[i] = *it;
            }
        }

        for (i = 0; i < nelements; ++i)
        {
            GenElement2D(Nodes, i, ElementFaces[i], FaceNodes, nComposite);
        }

        nComposite++;
    }
}

static void PrismLineFaces(int prismid,
                           map<int, int> &facelist,
                           vector<vector<int> > &FacesToPrisms,
                           vector<vector<int> > &PrismsToFaces,
                           vector<bool> &PrismDone);

void InputTec::ResetNodes(vector<NodeSharedPtr> &Vnodes,
                          Array<OneD, vector<int> > &ElementFaces,
                          vector<vector<int> > &FaceNodes)
{
    int i, j;
    Array<OneD, int> NodeReordering(Vnodes.size(), -1);
    int face1_map[3] = {0, 1, 4};
    int face3_map[3] = {3, 2, 5};
    int nodeid       = 0;
    map<int, bool> FacesRenumbered;

    // Determine Prism triangular face connectivity.
    vector<vector<int> > FaceToPrisms(FaceNodes.size());
    vector<vector<int> > PrismToFaces(ElementFaces.size());
    map<int, int> Prisms;

    // generate map of prism-faces to prisms and prism to
    // triangular-faces as well as ids of each prism.
    for (i = 0; i < ElementFaces.size(); ++i)
    {
        // Find Prism (and pyramids!).
        if (ElementFaces[i].size() == 5)
        {
            vector<int> LocTriFaces;
            // Find triangular faces
            for (j = 0; j < ElementFaces[i].size(); ++j)
            {
                if (FaceNodes[ElementFaces[i][j]].size() == 3)
                {
                    LocTriFaces.push_back(j);
                }
            }

            if (LocTriFaces.size() == 2) // prism otherwise a pyramid
            {
                Prisms[i] = i;

                PrismToFaces[i].push_back(ElementFaces[i][LocTriFaces[0]]);
                PrismToFaces[i].push_back(ElementFaces[i][LocTriFaces[1]]);

                FaceToPrisms[ElementFaces[i][LocTriFaces[0]]].push_back(i);
                FaceToPrisms[ElementFaces[i][LocTriFaces[1]]].push_back(i);
            }
        }
    }

    vector<bool> FacesDone(FaceNodes.size(), false);
    vector<bool> PrismDone(ElementFaces.size(), false);

    // For every prism find the list of prismatic elements
    // that represent an aligned block of cells. Then renumber
    // these blocks consecutativiesly
    for (auto &PrismIt : Prisms)
    {
        int elmtid = PrismIt.first;
        map<int, int> facelist;

        if (PrismDone[elmtid])
        {
            continue;
        }
        else
        {
            // Generate list of faces in list
            PrismLineFaces(
                elmtid, facelist, FaceToPrisms, PrismToFaces, PrismDone);

            // loop over faces and number vertices of associated prisms.
            for (auto &faceIt : facelist)
            {
                int faceid = faceIt.second;

                for (i = 0; i < FaceToPrisms[faceid].size(); ++i)
                {
                    int prismid = FaceToPrisms[faceid][i];

                    if ((FacesDone[PrismToFaces[prismid][0]] == true) &&
                        (FacesDone[PrismToFaces[prismid][1]] == true))
                    {
                        continue;
                    }

                    Array<OneD, int> Nodes =
                        SortFaceNodes(Vnodes, ElementFaces[prismid], FaceNodes);

                    if ((FacesDone[PrismToFaces[prismid][0]] == false) &&
                        (FacesDone[PrismToFaces[prismid][1]] == false))
                    {
                        // number all nodes consecutive since
                        // already correctly re-arranged.
                        for (i = 0; i < 3; ++i)
                        {
                            if (NodeReordering[Nodes[face1_map[i]]] == -1)
                            {
                                NodeReordering[Nodes[face1_map[i]]] = nodeid++;
                            }
                        }

                        for (i = 0; i < 3; ++i)
                        {
                            if (NodeReordering[Nodes[face3_map[i]]] == -1)
                            {
                                NodeReordering[Nodes[face3_map[i]]] = nodeid++;
                            }
                        }
                    }
                    else if ((FacesDone[PrismToFaces[prismid][0]] == false) &&
                             (FacesDone[PrismToFaces[prismid][1]] == true))
                    {
                        // find node of highest id
                        int max_id1, max_id2;

                        max_id1 = (NodeReordering[Nodes[face3_map[0]]] <
                                   NodeReordering[Nodes[face3_map[1]]])
                                      ? 1
                                      : 0;
                        max_id2 = (NodeReordering[Nodes[face3_map[max_id1]]] <
                                   NodeReordering[Nodes[face3_map[2]]])
                                      ? 2
                                      : max_id1;

                        // add numbering according to order of
                        int id0 = (max_id1 == 1) ? 0 : 1;

                        if (NodeReordering[Nodes[face1_map[id0]]] == -1)
                        {
                            NodeReordering[Nodes[face1_map[id0]]] = nodeid++;
                        }

                        if (NodeReordering[Nodes[face1_map[max_id1]]] == -1)
                        {
                            NodeReordering[Nodes[face1_map[max_id1]]] =
                                nodeid++;
                        }

                        if (NodeReordering[Nodes[face1_map[max_id2]]] == -1)
                        {
                            NodeReordering[Nodes[face1_map[max_id2]]] =
                                nodeid++;
                        }
                    }
                    else if ((FacesDone[PrismToFaces[prismid][0]] == true) &&
                             (FacesDone[PrismToFaces[prismid][1]] == false))
                    {
                        // find node of highest id
                        int max_id1, max_id2;

                        max_id1 = (NodeReordering[Nodes[face1_map[0]]] <
                                   NodeReordering[Nodes[face1_map[1]]])
                                      ? 1
                                      : 0;
                        max_id2 = (NodeReordering[Nodes[face1_map[max_id1]]] <
                                   NodeReordering[Nodes[face1_map[2]]])
                                      ? 2
                                      : max_id1;

                        // add numbering according to order of
                        int id0 = (max_id1 == 1) ? 0 : 1;

                        if (NodeReordering[Nodes[face3_map[id0]]] == -1)
                        {
                            NodeReordering[Nodes[face3_map[id0]]] = nodeid++;
                        }

                        if (NodeReordering[Nodes[face3_map[max_id1]]] == -1)
                        {
                            NodeReordering[Nodes[face3_map[max_id1]]] =
                                nodeid++;
                        }

                        if (NodeReordering[Nodes[face3_map[max_id2]]] == -1)
                        {
                            NodeReordering[Nodes[face3_map[max_id2]]] =
                                nodeid++;
                        }
                    }
                }
            }
        }
    }

    // fill in any unset nodes at from other shapes
    for (i = 0; i < NodeReordering.size(); ++i)
    {
        if (NodeReordering[i] == -1)
        {
            NodeReordering[i] = nodeid++;
        }
    }

    ASSERTL1(nodeid == NodeReordering.size(),
             "Have not renumbered all nodes");

    // Renumbering successfull so resort nodes and faceNodes;
    for (i = 0; i < FaceNodes.size(); ++i)
    {
        for (j = 0; j < FaceNodes[i].size(); ++j)
        {
            FaceNodes[i][j] = NodeReordering[FaceNodes[i][j]];
        }
    }

    vector<NodeSharedPtr> save(Vnodes);
    for (i = 0; i < Vnodes.size(); ++i)
    {
        Vnodes[NodeReordering[i]] = save[i];
        Vnodes[NodeReordering[i]]->SetID(NodeReordering[i]);
    }
}

static void PrismLineFaces(int prismid,
                           map<int, int> &facelist,
                           vector<vector<int> > &FaceToPrisms,
                           vector<vector<int> > &PrismToFaces,
                           vector<bool> &PrismDone)
{
    if (PrismDone[prismid] == false)
    {
        PrismDone[prismid] = true;

        // Add faces0
        int face       = PrismToFaces[prismid][0];
        facelist[face] = face;
        for (int i = 0; i < FaceToPrisms[face].size(); ++i)
        {
            PrismLineFaces(FaceToPrisms[face][i],
                           facelist,
                           FaceToPrisms,
                           PrismToFaces,
                           PrismDone);
        }

        // Add faces1
        face           = PrismToFaces[prismid][1];
        facelist[face] = face;
        for (int i = 0; i < FaceToPrisms[face].size(); ++i)
        {
            PrismLineFaces(FaceToPrisms[face][i],
                           facelist,
                           FaceToPrisms,
                           PrismToFaces,
                           PrismDone);
        }
    }
}

void InputTec::GenElement2D(vector<NodeSharedPtr> &VertNodes,
                            int i,
                            vector<int> &ElementFaces,
                            vector<vector<int> > &FaceNodes,
                            int nComposite)
{
    boost::ignore_unused(i);

    LibUtilities::ShapeType elType = (LibUtilities::ShapeType)0;
    // set up Node list

    if (ElementFaces.size() == 3)
    {
        elType = LibUtilities::eTriangle;
    }
    else if (ElementFaces.size() == 4)
    {
        elType = LibUtilities::eQuadrilateral;
    }
    else
    {
        NEKERROR(ErrorUtil::efatal,
                 "Not set up for elements which are not Tris or Quads");
    }

    // Create element tags
    vector<int> tags;
    tags.push_back(nComposite);

    // make unique node list
    vector<NodeSharedPtr> nodeList;
    Array<OneD, int> Nodes = SortEdgeNodes(VertNodes, ElementFaces, FaceNodes);
    for (int j = 0; j < Nodes.size(); ++j)
    {
        nodeList.push_back(VertNodes[Nodes[j]]);
    }

    // Create element
    ElmtConfig conf(elType, 1, true, true);
    ElementSharedPtr E =
        GetElementFactory().CreateInstance(elType, conf, nodeList, tags);

    m_mesh->m_element[E->GetDim()].push_back(E);
}

void InputTec::GenElement3D(vector<NodeSharedPtr> &VertNodes,
                            int i,
                            vector<int> &ElementFaces,
                            vector<vector<int> > &FaceNodes,
                            int nComposite,
                            bool DoOrient)
{
    boost::ignore_unused(i);

    LibUtilities::ShapeType elType = (LibUtilities::ShapeType)0;
    // set up Node list
    Array<OneD, int> Nodes = SortFaceNodes(VertNodes, ElementFaces, FaceNodes);
    int nnodes             = Nodes.size();
    map<LibUtilities::ShapeType, int> domainComposite;

    // Set Nodes  -- Not sure we need this so could
    // m_mesh->m_node = VertNodes;

    // element type
    if (nnodes == 4)
    {
        elType = LibUtilities::eTetrahedron;
    }
    else if (nnodes == 5)
    {
        elType = LibUtilities::ePyramid;
    }
    else if (nnodes == 6)
    {
        elType = LibUtilities::ePrism;
    }
    else
    {
        NEKERROR(ErrorUtil::efatal,
                 "Not set up for elements which are not Tets, "
                 "Prisms or Pyramids.");
    }

    // Create element tags
    vector<int> tags;
    tags.push_back(nComposite);

    // make unique node list
    vector<NodeSharedPtr> nodeList;
    for (int j = 0; j < Nodes.size(); ++j)
    {
        nodeList.push_back(VertNodes[Nodes[j]]);
    }

    // Create element
    if (elType != LibUtilities::ePyramid)
    {
        ElmtConfig conf(elType, 1, true, true, DoOrient);
        ElementSharedPtr E =
            GetElementFactory().CreateInstance(elType, conf, nodeList, tags);

        m_mesh->m_element[E->GetDim()].push_back(E);
    }
    else
    {
        cout << "Warning: Pyramid detected " << endl;
    }
}

Array<OneD, int> InputTec::SortEdgeNodes(vector<NodeSharedPtr> &Vnodes,
                                         vector<int> &ElementFaces,
                                         vector<vector<int> > &FaceNodes)
{
    int i, j;
    Array<OneD, int> returnval;

    if (ElementFaces.size() == 3) // Triangle
    {
        returnval = Array<OneD, int>(3);

        returnval[0] = FaceNodes[ElementFaces[0]][0];
        returnval[1] = FaceNodes[ElementFaces[0]][1];

        // Find third node index;
        for (i = 0; i < 2; ++i)
        {
            if ((FaceNodes[ElementFaces[1]][i] != returnval[0]) &&
                (FaceNodes[ElementFaces[1]][i] != returnval[1]))
            {
                returnval[2] = FaceNodes[ElementFaces[1]][i];
                break;
            }
        }
    }
    else if (ElementFaces.size() == 4) // quadrilateral
    {
        returnval = Array<OneD, int>(4);

        int indx0 = FaceNodes[ElementFaces[0]][0];
        int indx1 = FaceNodes[ElementFaces[0]][1];
        int indx2, indx3;

        indx2 = indx3 = -1;
        // Find third, fourth node index;
        for (j = 1; j < 4; ++j)
        {
            for (i = 0; i < 2; ++i)
            {
                if ((FaceNodes[ElementFaces[j]][i] != indx0) &&
                    (FaceNodes[ElementFaces[j]][i] != indx1))
                {
                    if (indx2 == -1)
                    {
                        indx2 = FaceNodes[ElementFaces[j]][i];
                    }
                    else if (indx2 != -1)
                    {
                        if (FaceNodes[ElementFaces[j]][i] != indx2)
                        {
                            indx3 = FaceNodes[ElementFaces[j]][i];
                        }
                    }
                }
            }
        }

        ASSERTL1((indx2 != -1) && (indx3 != -1),
                 "Failed to find vertex 3 or 4");

        // calculate 0-1,
        Node a = *(Vnodes[indx1]) - *(Vnodes[indx0]);
        // calculate 0-2,
        Node b      = *(Vnodes[indx2]) - *(Vnodes[indx0]);
        Node acurlb = a.curl(b);

        // calculate 2-1,
        Node c = *(Vnodes[indx1]) - *(Vnodes[indx2]);
        // calculate 3-2,
        Node d      = *(Vnodes[indx3]) - *(Vnodes[indx2]);
        Node acurld = a.curl(d);

        NekDouble acurlb_dot_acurld = acurlb.dot(acurld);
        if (acurlb_dot_acurld > 0.0)
        {
            returnval[0] = indx0;
            returnval[1] = indx1;
            returnval[2] = indx2;
            returnval[3] = indx3;
        }
        else
        {
            returnval[0] = indx0;
            returnval[1] = indx1;
            returnval[2] = indx3;
            returnval[3] = indx2;
        }
    }

    return returnval;
}

Array<OneD, int> InputTec::SortFaceNodes(vector<NodeSharedPtr> &Vnodes,
                                         vector<int> &ElementFaces,
                                         vector<vector<int> > &FaceNodes)
{

    int i, j;
    Array<OneD, int> returnval;

    if (ElementFaces.size() == 4) // Tetrahedron
    {
        ASSERTL1(FaceNodes[ElementFaces[0]].size() == 3,
                 "Face is not triangular");

        returnval = Array<OneD, int>(4);

        int indx0 = FaceNodes[ElementFaces[0]][0];
        int indx1 = FaceNodes[ElementFaces[0]][1];
        int indx2 = FaceNodes[ElementFaces[0]][2];
        int indx3 = -1;

        // calculate 0-1,
        Node a = *(Vnodes[indx1]) - *(Vnodes[indx0]);
        // calculate 0-2,
        Node b = *(Vnodes[indx2]) - *(Vnodes[indx0]);

        // Find fourth node index;
        ASSERTL1(FaceNodes[ElementFaces[1]].size() == 3,
                 "Face is not triangular");
        for (i = 0; i < 3; ++i)
        {

            if ((FaceNodes[ElementFaces[1]][i] != indx0) &&
                (FaceNodes[ElementFaces[1]][i] != indx1) &&
                (FaceNodes[ElementFaces[1]][i] != indx2))
            {
                indx3 = FaceNodes[ElementFaces[1]][i];
                break;
            }
        }

        // calculate 0-3,
        Node c      = *(Vnodes[indx3]) - *(Vnodes[indx0]);
        Node acurlb = a.curl(b);

        NekDouble acurlb_dotc = acurlb.dot(c);
        if (acurlb_dotc < 0.0)
        {
            returnval[0] = indx0;
            returnval[1] = indx1;
            returnval[2] = indx2;
            returnval[3] = indx3;
        }
        else
        {
            returnval[0] = indx1;
            returnval[1] = indx0;
            returnval[2] = indx2;
            returnval[3] = indx3;
        }
    }
    else if (ElementFaces.size() == 5) // prism or pyramid
    {
        int triface0, triface1;
        int quadface0, quadface1, quadface2;
        bool isPrism = true;

        // find ids of tri faces and first quad face
        triface0 = triface1 = -1;
        quadface0 = quadface1 = quadface2 = -1;
        for (i = 0; i < 5; ++i)
        {
            if (FaceNodes[ElementFaces[i]].size() == 3)
            {
                if (triface0 == -1)
                {
                    triface0 = i;
                }
                else if (triface1 == -1)
                {
                    triface1 = i;
                }
                else
                {
                    isPrism = false;
                }
            }

            if (FaceNodes[ElementFaces[i]].size() == 4)
            {
                if (quadface0 == -1)
                {
                    quadface0 = i;
                }
                else if (quadface1 == -1)
                {
                    quadface1 = i;
                }
                else if (quadface2 == -1)
                {
                    quadface2 = i;
                }
            }
        }

        if (isPrism) // Prism
        {
            returnval = Array<OneD, int>(6);
        }
        else // Pyramid
        {
            returnval = Array<OneD, int>(5);
        }

        // find matching nodes between triface0 and triquad0
        int indx0, indx1, indx2, indx3, indx4;

        indx0 = indx1 = indx2 = indx3 = indx4 = -1;
        // Loop over all quad nodes and if they match any
        // triangular nodes If they do set these to indx0 and
        // indx1 and if not set it to indx2, indx3

        for (i = 0; i < 4; ++i)
        {
            for (j = 0; j < 3; ++j)
            {
                if (FaceNodes[ElementFaces[triface0]][j] ==
                    FaceNodes[ElementFaces[quadface0]][i])
                {
                    break; // same node break
                }
            }

            if (j == 3) // Vertex not in quad face
            {
                if (indx2 == -1)
                {
                    indx2 = FaceNodes[ElementFaces[quadface0]][i];
                }
                else if (indx3 == -1)
                {
                    indx3 = FaceNodes[ElementFaces[quadface0]][i];
                }
                else
                {
                    ASSERTL0(
                        false,
                        "More than two vertices do not match triangular face");
                }
            }
            else // if found match then set indx0,indx1;
            {
                if (indx0 == -1)
                {
                    indx0 = FaceNodes[ElementFaces[quadface0]][i];
                }
                else
                {
                    indx1 = FaceNodes[ElementFaces[quadface0]][i];
                }
            }
        }

        // Finally check for top vertex
        for (int i = 0; i < 3; ++i)
        {
            if ((FaceNodes[ElementFaces[triface0]][i] != indx0) &&
                (FaceNodes[ElementFaces[triface0]][i] != indx1) &&
                (FaceNodes[ElementFaces[triface0]][i] != indx2))
            {
                indx4 = FaceNodes[ElementFaces[triface0]][i];
                break;
            }
        }

        // calculate 0-1,
        Node a = *(Vnodes[indx1]) - *(Vnodes[indx0]);
        // calculate 0-4,
        Node b = *(Vnodes[indx4]) - *(Vnodes[indx0]);
        // calculate 0-2,
        Node c      = *(Vnodes[indx2]) - *(Vnodes[indx0]);
        Node acurlb = a.curl(b);

        NekDouble acurlb_dotc = acurlb.dot(c);
        if (acurlb_dotc < 0.0)
        {
            returnval[0] = indx0;
            returnval[1] = indx1;
            returnval[4] = indx4;
        }
        else
        {
            returnval[0] = indx1;
            returnval[1] = indx0;
            returnval[4] = indx4;
        }

        // check to see if two vertices are shared between one of the other
        // faces
        // to define which is indx2 and indx3

        int cnt = 0;
        for (int i = 0; i < 4; ++i)
        {
            if ((FaceNodes[ElementFaces[quadface1]][i] == returnval[1]) ||
                (FaceNodes[ElementFaces[quadface1]][i] == indx2))
            {
                cnt++;
            }
        }

        if (cnt == 2) // have two matching vertices
        {
            returnval[2] = indx2;
            returnval[3] = indx3;
        }
        else
        {
            cnt = 0;
            for (int i = 0; i < 4; ++i)
            {
                if ((FaceNodes[ElementFaces[quadface2]][i] == returnval[1]) ||
                    (FaceNodes[ElementFaces[quadface2]][i] == indx2))
                {
                    cnt++;
                }
            }

            if (cnt != 2) // neither of the other faces has two matching nodes
                          // so reverse
            {
                returnval[2] = indx3;
                returnval[3] = indx2;
            }
            else // have two matching vertices
            {
                returnval[2] = indx2;
                returnval[3] = indx3;
            }
        }

        if (isPrism == true)
        {
            // finally need to find last vertex from second triangular face.
            for (int i = 0; i < 3; ++i)
            {
                if ((FaceNodes[ElementFaces[triface1]][i] != indx2) &&
                    (FaceNodes[ElementFaces[triface1]][i] != indx3) &&
                    (FaceNodes[ElementFaces[triface1]][i] != indx3))
                {
                    returnval[5] = FaceNodes[ElementFaces[triface1]][i];
                    break;
                }
            }
        }
    }
    else
    {
        ASSERTL0(false, "SortFaceNodes not set up for this number of faces");
    }

    return returnval;
}
}
}
