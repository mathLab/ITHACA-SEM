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
#include "InputStar.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
static char const kDefaultState[] = "default";
namespace Utilities
{

ModuleKey InputStar::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eInputModule, "ccm"),
    InputStar::create,
    "Reads mesh from Star CCM (.ccm).");

InputStar::InputStar(MeshSharedPtr m) : InputModule(m)
{
    m_config["writelabelsonly"] = ConfigOption(
        true,
        "0",
        "Just write out tags from star file for each surface/composite");
}

InputStar::~InputStar()
{
}

/**
 * Tecplot file Polyhedron format contains a list of nodes, a node count per
 * face, the node ids, Element ids that are on the left of each face and Element
 * ids which are on the right of each face. There are then a series of zone of
 * each surface. In the case of a surface the number of nodes is not provided
 * indicating it is a 2D zone.
 *
 * @param pFilename Filename of Tecplot file to read.
 */
void InputStar::Process()
{
    m_mesh->m_expDim   = 3;
    m_mesh->m_spaceDim = 3;

    if (m_mesh->m_verbose)
    {
        cout << "InputCCM: Start reading file..." << endl;
    }

    InitCCM();

    SetupElements();

    PrintSummary();

    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();
}

void InputStar::SetupElements(void)
{
    int i;
    string line, tag;
    stringstream s;
    streampos pos;
    int nComposite = 0;

    // Read in Nodes
    ReadNodes(m_mesh->m_node);

    // Get list of faces nodes and adjacents elements.
    unordered_map<int, vector<int> > FaceNodes;
    Array<OneD, vector<int> > ElementFaces;

    // Read interior faces and set up first part of Element
    // Faces and FaceNodes
    ReadInternalFaces(FaceNodes, ElementFaces);

    vector<vector<int> > BndElementFaces;
    vector<string> Facelabels;
    ReadBoundaryFaces(BndElementFaces, FaceNodes, ElementFaces, Facelabels);

    if (m_config["writelabelsonly"].beenSet)
    {
        nComposite = 2;
        // write boundary zones/composites
        for (i = 0; i < BndElementFaces.size(); ++i)
        {
            cout << " 2D Zone (composite = " << nComposite
                 << ", label = " << Facelabels[i] << ")" << endl;
            nComposite++;
        }
        exit(1);
    }

    // 3D Zone
    // Reset node ordering so that all prism faces have
    // consistent numbering for singular vertex re-ordering
    ResetNodes(m_mesh->m_node, ElementFaces, FaceNodes);

    // create Prisms/Pyramids first
    int nelements = ElementFaces.size();
    cout << " Generating 3D Zones: " << endl;
    int cnt = 0;
    for (i = 0; i < nelements; ++i)
    {

        if (ElementFaces[i].size() > 4)
        {
            GenElement3D(
                m_mesh->m_node, i, ElementFaces[i], FaceNodes, nComposite, true);
            ++cnt;
        }
    }
    cout <<" \t" << cnt << " Prisms" << endl;

    nComposite++;

    // create Tets second
    cnt = 0;
    for (i = 0; i < nelements; ++i)
    {
        if (ElementFaces[i].size() == 4)
        {
            GenElement3D(
                m_mesh->m_node, i, ElementFaces[i], FaceNodes, nComposite, true);
            ++cnt;
        }
    }
    cout <<"\t" << cnt << " Tets" << endl;
    nComposite++;

    // Insert vertices into map.
    for (auto &node : m_mesh->m_node)
    {
        m_mesh->m_vertexSet.insert(node);
    }

    // Add boundary zones/composites
    for (i = 0; i < BndElementFaces.size(); ++i)
    {
        cout << " Generating 2D Zone (composite = " << nComposite
             << ", label = " << Facelabels[i] << ")" << endl;

        for (int j = 0; j < BndElementFaces[i].size(); ++j)
        {
            auto it = FaceNodes.find(BndElementFaces[i][j]);
            if (it != FaceNodes.end())
            {
                GenElement2D(m_mesh->m_node, j, it->second, nComposite);
            }
            else
            {
                string msg = "Failed to find FaceNodes for Face ";
                msg += boost::lexical_cast<string>(BndElementFaces[i][j]);
                ASSERTL0(false, msg);
            }
        }

        m_mesh->m_faceLabels[nComposite] = Facelabels[i];
        nComposite++;
    }
}

static void PrismLineFaces(int prismid,
                           map<int, int> &facelist,
                           vector<vector<int> > &FacesToPrisms,
                           vector<vector<int> > &PrismsToFaces,
                           vector<bool> &PrismDone);

void InputStar::ResetNodes(vector<NodeSharedPtr> &Vnodes,
                           Array<OneD, vector<int> > &ElementFaces,
                           unordered_map<int, vector<int> > &FaceNodes)
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

    // Renumbering successfull so reset nodes and faceNodes;
    for (auto &it : FaceNodes)
    {
        for (j = 0; j < it.second.size(); ++j)
        {
            it.second[j] = NodeReordering[it.second[j]];
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

void InputStar::GenElement2D(vector<NodeSharedPtr> &VertNodes,
                             int i,
                             vector<int> &FaceNodes,
                             int nComposite)
{
    boost::ignore_unused(i);

    LibUtilities::ShapeType elType;

    if (FaceNodes.size() == 3)
    {
        elType = LibUtilities::eTriangle;
    }
    else if (FaceNodes.size() == 4)
    {
        elType = LibUtilities::eQuadrilateral;
    }
    else
    {
        ASSERTL0(false, "Not set up for elements which are not Tets or Prism");
    }

    // Create element tags
    vector<int> tags;
    tags.push_back(nComposite);

    // make unique node list
    vector<NodeSharedPtr> nodeList;
    Array<OneD, int> Nodes = SortEdgeNodes(VertNodes, FaceNodes);
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

void InputStar::GenElement3D(vector<NodeSharedPtr> &VertNodes,
                             int i,
                             vector<int> &ElementFaces,
                             unordered_map<int, vector<int> > &FaceNodes,
                             int nComposite,
                             bool DoOrient)
{
    boost::ignore_unused(i);

    LibUtilities::ShapeType elType;
    // set up Node list
    Array<OneD, int> Nodes = SortFaceNodes(VertNodes, ElementFaces, FaceNodes);
    int nnodes             = Nodes.size();
    map<LibUtilities::ShapeType, int> domainComposite;

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

        ASSERTL0(false, "Not set up for elements which are not Tets or Prism");
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

Array<OneD, int> InputStar::SortEdgeNodes(vector<NodeSharedPtr> &Vnodes,
                                          vector<int> &FaceNodes)
{
    Array<OneD, int> returnval;

    if (FaceNodes.size() == 3) // Triangle
    {
        returnval = Array<OneD, int>(3);

        returnval[0] = FaceNodes[0];
        returnval[1] = FaceNodes[1];
        returnval[2] = FaceNodes[2];
    }
    else if (FaceNodes.size() == 4) // quadrilateral
    {
        returnval = Array<OneD, int>(4);

        int indx0 = FaceNodes[0];
        int indx1 = FaceNodes[1];
        int indx2 = FaceNodes[2];
        int indx3 = FaceNodes[3];

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

Array<OneD, int> InputStar::SortFaceNodes(vector<NodeSharedPtr> &Vnodes,
                                          vector<int> &ElementFaces,
                                          unordered_map<int, vector<int> > &FaceNodes)
{

    int i, j;
    Array<OneD, int> returnval;

    if (ElementFaces.size() == 4) // Tetrahedron
    {
        ASSERTL1(FaceNodes[ElementFaces[0]].size() == 3,
                 "Face is not triangular");

        returnval = Array<OneD, int>(4);

        auto it = FaceNodes.find(ElementFaces[0]);
        int indx0 = it->second[0];
        int indx1 = it->second[1];
        int indx2 = it->second[2];
        int indx3 = -1;

        // calculate 0-1,
        Node a = *(Vnodes[indx1]) - *(Vnodes[indx0]);
        // calculate 0-2,
        Node b = *(Vnodes[indx2]) - *(Vnodes[indx0]);

        // Find fourth node index;
        ASSERTL1(FaceNodes[ElementFaces[1]].size() == 3,
                 "Face is not triangular");

        auto it2 = FaceNodes.find(ElementFaces[1]);
        for (i = 0; i < 3; ++i)
        {
            if ((it2->second[i] != indx0) && (it2->second[i] != indx1) &&
                (it2->second[i] != indx2))
            {
                indx3 = it2->second[i];
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
        int triface0, triface1, triface2, triface3;
        int quadface0, quadface1, quadface2;
        bool isPrism = true;

        // find ids of tri faces and first quad face
        triface0 = triface1 = triface2 = triface3 = -1;
        quadface0 = quadface1 = quadface2 = -1;
        for (i = 0; i < 5; ++i)
        {
            auto it = FaceNodes.find(ElementFaces[i]);
            if (it->second.size() == 3)
            {
                if (triface0 == -1)
                {
                    triface0 = i;
                }
                else if (triface1 == -1)
                {
                    triface1 = i;
                }
                else if (triface2 == -1)
                {
                    triface2 = i;
                    isPrism = false;
                }
                else if (triface3 == -1)
                {
                    triface3 = i;
                }
            }

            if (it->second.size() == 4)
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
            ASSERTL1(quadface0 != -1, "Quad face 0 not found");
            ASSERTL1(quadface1 != -1, "Quad face 1 not found");
            ASSERTL1(quadface2 != -1, "Quad face 2 not found");
            ASSERTL1(triface0 != -1, "Tri face 0 not found");
            ASSERTL1(triface1 != -1, "Tri face 1 not found");
        }
        else // Pyramid
        {
            ASSERTL1(quadface0 != -1, "Quad face 0 not found");
            ASSERTL1(triface0 != -1, "Tri face 0 not found");
            ASSERTL1(triface1 != -1, "Tri face 1 not found");
            ASSERTL1(triface2 != -1, "Tri face 2 not found");
            ASSERTL1(triface3 != -1, "Tri face 3 not found");
            ASSERTL0(false,"Pyramids still not sorted");
            returnval = Array<OneD, int>(5);
        }

        // find matching nodes between triface0 and triquad0
        int indx0, indx1, indx2, indx3, indx4;

        indx0 = indx1 = indx2 = indx3 = indx4 = -1;
        // Loop over all quad nodes and if they match any
        // triangular nodes If they do set these to indx0 and
        // indx1 and if not set it to indx2, indx3

        auto &triface0_vec = FaceNodes.find(ElementFaces[triface0])->second;
        auto &quadface0_vec = FaceNodes.find(ElementFaces[quadface0])->second;
        for (i = 0; i < 4; ++i)
        {
            for (j = 0; j < 3; ++j)
            {
                if (triface0_vec[j] == quadface0_vec[i])
                {
                    break; // same node break
                }
            }

            if (j == 3) // Vertex not in quad face
            {
                if (indx2 == -1)
                {
                    indx2 = quadface0_vec[i];
                }
                else if (indx3 == -1)
                {
                    indx3 = quadface0_vec[i];
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
                    indx0 = quadface0_vec[i];
                }
                else
                {
                    indx1 = quadface0_vec[i];
                }
            }
        }

        // Finally check for top vertex
        for (int i = 0; i < 3; ++i)
        {
            if (triface0_vec[i] != indx0 && triface0_vec[i] != indx1 &&
                triface0_vec[i] != indx2)
            {
                indx4 = triface0_vec[i];
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

        auto &quadface1_vec = FaceNodes.find(ElementFaces[quadface1])->second;
        auto &quadface2_vec = FaceNodes.find(ElementFaces[quadface2])->second;
        int cnt = 0;
        for (int i = 0; i < 4; ++i)
        {
            if (quadface1_vec[i] == returnval[1] || quadface1_vec[i] == indx2)
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
                if (quadface2_vec[i] == returnval[1] || quadface2_vec[i] == indx2)
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
            auto &triface1_vec = FaceNodes.find(ElementFaces[triface1])->second;
            for (int i = 0; i < 3; ++i)
            {
                if (triface1_vec[i] != indx2 && triface1_vec[i] != indx3)
                {
                    returnval[5] = triface1_vec[i];
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

// initialise and read ccm file to ccm structure
void InputStar::InitCCM(void)
{
    // Open ccm file for reading.
    CCMIOID root;
    // Open the file.  Because we did not initialize 'err' we
    // need to pass in NULL (which always means kCCMIONoErr)
    // and then assign the return value to 'err'.).
    string fname = m_config["infile"].as<string>();
    m_ccmErr     = CCMIOOpenFile(NULL, fname.c_str(), kCCMIORead, &root);
    ASSERTL0(m_ccmErr == kCCMIONoErr,"Error opening file");

    int i = 0;
    CCMIOID state, problem;

    // We are going to assume that we have a state with a
    // known name.  We could instead use CCMIONextEntity() to
    // walk through all the states in the file and present the
    // list to the user for selection.
    CCMIOGetState(&m_ccmErr, root, kDefaultState, &problem, &state);
    if (m_ccmErr != kCCMIONoErr)
    {
        cout << "No state named '" << kDefaultState << "'" << endl;
        exit(0);
    }

    // Find the first processor (i has previously been
    // initialized to 0) and read the mesh and solution
    // information.
    CCMIONextEntity(&m_ccmErr, state, kCCMIOProcessor, &i, &m_ccmProcessor);
    ASSERTL0(m_ccmErr == kCCMIONoErr,"Failed to find Next Entity");
}

void InputStar::ReadNodes(std::vector<NodeSharedPtr> &Nodes)
{
    CCMIOID mapID, vertices;
    CCMIOSize nVertices;
    int dims = 1;

    CCMIOReadProcessor(
        &m_ccmErr, m_ccmProcessor, &vertices, &m_ccmTopology, NULL, NULL);
    ASSERTL0(m_ccmErr == kCCMIONoErr,"Error Reading Processor");
    CCMIOEntitySize(&m_ccmErr, vertices, &nVertices, NULL);
    ASSERTL0(m_ccmErr == kCCMIONoErr,"Error Reading NextEntitySize in ReadNodes");

    // Read the vertices.  This involves reading both the vertex data and
    // the map, which maps the index into the data array with the ID number.
    // As we process the vertices we need to be sure to scale them by the
    // appropriate scaling factor.  The offset is just to show you can read
    // any chunk.  Normally this would be in a for loop.
    float scale;
    int nvert  = nVertices;
    vector<int> mapData;
    mapData.resize(nvert);
    vector<float> verts;
    verts.resize(3*nvert);

    for (int k = 0; k < nvert; ++k)
    {
        verts[3 * k] = verts[3 * k + 1] = verts[3 * k + 2] = 0.0;
        mapData[k] = 0;
    }

    CCMIOReadVerticesf(&m_ccmErr,
                       vertices,
                       &dims,
                       &scale,
                       &mapID,
                       &verts[0],
                       0,
                       nVertices);
    ASSERTL0(m_ccmErr == kCCMIONoErr,"Error Reading Vertices in ReadNodes");
    CCMIOReadMap(&m_ccmErr,
                 mapID,
                 &mapData[0],
                 0,
                 nVertices);
    ASSERTL0(m_ccmErr == kCCMIONoErr,"Error Reading Map in ReadNodes");

    for (int i = 0; i < nVertices; ++i)
    {
        Nodes.push_back(std::make_shared<Node>(
                            i, verts[3 * i], verts[3 * i + 1], verts[3 * i + 2]));
    }
}

void InputStar::ReadInternalFaces(unordered_map<int, vector<int> > &FacesNodes,
                                  Array<OneD, vector<int> > &ElementFaces)
{

    CCMIOID mapID, id;
    CCMIOSize nFaces, size;
    vector<int> faces, faceCells, mapData;

    // Read the internal faces.
    CCMIOGetEntity(&m_ccmErr, m_ccmTopology, kCCMIOInternalFaces, 0, &id);
    CCMIOEntitySize(&m_ccmErr, id, &nFaces, NULL);

    int nf = nFaces;
    mapData.resize(nf);
    faceCells.resize(2 * nf);

    CCMIOReadFaces(&m_ccmErr,
                   id,
                   kCCMIOInternalFaces,
                   NULL,
                   &size,
                   NULL,
                   kCCMIOStart,
                   kCCMIOEnd);
    faces.resize((size_t)size);
    CCMIOReadFaces(&m_ccmErr,
                   id,
                   kCCMIOInternalFaces,
                   &mapID,
                   NULL,
                   &faces[0],
                   kCCMIOStart,
                   kCCMIOEnd);
    CCMIOReadFaceCells(&m_ccmErr,
                       id,
                       kCCMIOInternalFaces,
                       &faceCells[0],
                       kCCMIOStart,
                       kCCMIOEnd);
    CCMIOReadMap(&m_ccmErr,
                 mapID,
                 &mapData[0],
                 kCCMIOStart,
                 kCCMIOEnd);

    // Add face nodes
    int cnt = 0;
    for (int i = 0; i < nf; ++i)
    {
        vector<int> Fnodes;
        int j;
        if (cnt < faces.size())
        {
            int nv = faces[cnt];
            ASSERTL0(nv <= 4,
                     "Can only handle meshes with "
                     "up to four nodes per face");

            for (j = 0; j < nv; ++j)
            {
                if (cnt + 1 + j < faces.size())
                {
                    Fnodes.push_back(faces[cnt + 1 + j] - 1);
                }
            }
            cnt += nv + 1;
        }
        FacesNodes[mapData[i] - 1] = Fnodes;
    }

    // find number of elements;
    int nelmt = 0;
    for (int i = 0; i < faceCells.size(); ++i)
    {
        nelmt = max(nelmt, faceCells[i]);
    }

    ElementFaces = Array<OneD, vector<int> >(nelmt);
    for (int i = 0; i < nf; ++i)
    {
        // left element
        if (faceCells[2 * i])
        {
            ElementFaces[faceCells[2 * i] - 1].push_back(mapData[i] - 1);
        }

        // right element
        if (faceCells[2 * i + 1])
        {
            ElementFaces[faceCells[2 * i + 1] - 1].push_back(mapData[i] - 1);
        }
    }
}

void InputStar::ReadBoundaryFaces(vector<vector<int> > &BndElementFaces,
                                  unordered_map<int, vector<int> > &FacesNodes,
                                  Array<OneD, vector<int> > &ElementFaces,
                                  vector<string> &Facelabels)
{
    // Read the boundary faces.
    int index = 0;
    CCMIOID mapID, id;
    CCMIOSize nFaces, size;
    vector<int> faces, faceCells, mapData;
    vector<string> facelabel;

    while (CCMIONextEntity(
               NULL, m_ccmTopology, kCCMIOBoundaryFaces, &index, &id) ==
           kCCMIONoErr)
    {
        int boundaryVal;

        CCMIOEntitySize(&m_ccmErr, id, &nFaces, NULL);
        CCMIOSize nf = nFaces;
        mapData.resize(nf);
        faceCells.resize(nf);
        CCMIOReadFaces(&m_ccmErr,
                       id,
                       kCCMIOBoundaryFaces,
                       NULL,
                       &size,
                       NULL,
                       kCCMIOStart,
                       kCCMIOEnd);

        faces.resize((size_t)size);
        CCMIOReadFaces(&m_ccmErr,
                       id,
                       kCCMIOBoundaryFaces,
                       &mapID,
                       NULL,
                       &faces[0],
                       kCCMIOStart,
                       kCCMIOEnd);
        CCMIOReadFaceCells(&m_ccmErr,
                           id,
                           kCCMIOBoundaryFaces,
                           &faceCells[0],
                           kCCMIOStart,
                           kCCMIOEnd);
        CCMIOReadMap(&m_ccmErr,
                     mapID,
                     &mapData[0],
                     kCCMIOStart,
                     kCCMIOEnd);

        CCMIOGetEntityIndex(&m_ccmErr, id, &boundaryVal);

        // check to see if we have a label for this boundary faces
        int size;
        char *name;
        if (CCMIOReadOptstr(NULL, id, "Label", &size, NULL) == kCCMIONoErr)
        {
            name = new char[size + 1];
            CCMIOReadOptstr(NULL, id, "Label", NULL, name);
            Facelabels.push_back(string(name));
        }
        else
        {
            Facelabels.push_back("Not known");
        }

        // Add face nodes
        int cnt = 0;
        for (int i = 0; i < nf; ++i)
        {
            vector<int> Fnodes;
            int j;
            if (cnt < faces.size())
            {
                int nv = faces[cnt];
                ASSERTL0(nv <= 4,
                         "Can only handle meshes with "
                         "up to four nodes per face");

                for (j = 0; j < nv; ++j)
                {
                    if (cnt + 1 + j < faces.size())
                    {
                        Fnodes.push_back(faces[cnt + 1 + j] - 1);
                    }
                }
                cnt += nv + 1;
            }
            FacesNodes[mapData[i] - 1] = Fnodes;
        }

        vector<int> BndFaces;
        for (int i = 0; i < nf; ++i)
        {
            if (faceCells[i])
            {
                ElementFaces[faceCells[i] - 1].push_back(mapData[i] - 1);
            }
            BndFaces.push_back(mapData[i] - 1);
        }
        BndElementFaces.push_back(BndFaces);
    }
}
}
}
