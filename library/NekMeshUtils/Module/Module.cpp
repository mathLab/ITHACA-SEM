////////////////////////////////////////////////////////////////////////////////
//
//  File: Module.cpp
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
//  Description: Abstract input/output modules.
//
////////////////////////////////////////////////////////////////////////////////

#include <boost/iostreams/filter/gzip.hpp>

#include "Module.h"

using namespace std;
namespace io = boost::iostreams;

namespace Nektar
{
namespace NekMeshUtils
{

/**
 * Returns an instance of the module factory, held as a singleton.
 */
ModuleFactory& GetModuleFactory()
{
    static ModuleFactory instance;
    return instance;
}

/**
 * Prints a given module key to a stream.
 */
std::ostream& operator<<(std::ostream& os, const ModuleKey& rhs)
{
    return os << ModuleTypeMap[rhs.first] << ": " << rhs.second;
}

InputModule::InputModule(MeshSharedPtr m) : Module(m)
{
    m_config["infile"] = ConfigOption(false, "", "Input filename.");
}

OutputModule::OutputModule(MeshSharedPtr m) : Module(m)
{
    m_config["outfile"] = ConfigOption(false, "", "Output filename.");
}

/**
 * @brief Open a file for input.
 */
void InputModule::OpenStream()
{
    string fname = m_config["infile"].as<string>();

    if (fname.size() > 3 && fname.substr(fname.size() - 3, 3) == ".gz")
    {
        m_mshFileStream.open(fname.c_str(), ios_base::in | ios_base::binary);
        m_mshFile.push(io::gzip_decompressor());
        m_mshFile.push(m_mshFileStream);
    }
    else
    {
        m_mshFileStream.open(fname.c_str());
        m_mshFile.push(m_mshFileStream);
    }

    if (!m_mshFile.good())
    {
        cerr << "Error opening file: " << fname << endl;
        abort();
    }
}

/**
 * @brief Open a file for output.
 */
void OutputModule::OpenStream()
{
    string fname = m_config["outfile"].as<string>();

    if (fname.size() > 3 && fname.substr(fname.size() - 3, 3) == ".gz")
    {
        m_mshFileStream.open(fname.c_str(), ios_base::out | ios_base::binary);
        m_mshFile.push(io::gzip_compressor());
        m_mshFile.push(m_mshFileStream);
    }
    else
    {
        m_mshFileStream.open(fname.c_str());
        m_mshFile.push(m_mshFileStream);
    }

    if (!m_mshFile.good())
    {
        cerr << "Error opening file: " << fname << endl;
        abort();
    }
}

/**
 * @brief Create a unique set of mesh vertices from elements stored in
 * Mesh::element.
 *
 * Each element is processed in turn and the vertices extracted and
 * inserted into #m_vertexSet, which at the end of the routine
 * contains all unique vertices in the mesh.
 */
void Module::ProcessVertices()
{
    vector<ElementSharedPtr> &elmt = m_mesh->m_element[m_mesh->m_expDim];

    m_mesh->m_vertexSet.clear();

    for (int i = 0, vid = 0; i < elmt.size(); ++i)
    {
        for (int j = 0; j < elmt[i]->GetVertexCount(); ++j)
        {
            pair<NodeSet::iterator,bool> testIns =
                m_mesh->m_vertexSet.insert(elmt[i]->GetVertex(j));

            if (testIns.second)
            {
                (*(testIns.first))->m_id = vid++;
            }
            else
            {
                elmt[i]->SetVertex(j,*testIns.first);
            }
        }
    }
}

/**
 * @brief Create a unique set of mesh edges from elements stored in
 * Mesh::element.
 *
 * All elements are first scanned and a list of unique, enumerated
 * edges produced in #m_edgeSet. Since each element generated its
 * edges independently, we must now ensure that each element only uses
 * edge objects from the #m_edgeSet set This ensures there are no
 * duplicate edge objects. Finally, we scan the list of elements for
 * 1-D boundary elements which correspond to an edge in
 * #m_edgeSet. For such elements, we set its edgeLink to reference the
 * corresponding edge in #m_edgeSet.
 *
 * This routine only proceeds if the expansion dimension is 2 or 3.
 */
void Module::ProcessEdges(bool ReprocessEdges)
{
    if (m_mesh->m_expDim < 2) return;

    if(ReprocessEdges)
    {
        vector<ElementSharedPtr> &elmt = m_mesh->m_element[m_mesh->m_expDim];

        m_mesh->m_edgeSet.clear();

        // Clear all edge links

        // Scan all elements and generate list of unique edges
        for (int i = 0, eid = 0; i < elmt.size(); ++i)
        {
            for (int j = 0; j < elmt[i]->GetEdgeCount(); ++j)
            {
                pair<EdgeSet::iterator,bool> testIns;
                EdgeSharedPtr ed = elmt[i]->GetEdge(j);
                testIns = m_mesh->m_edgeSet.insert(ed);

                if (testIns.second)
                {
                    EdgeSharedPtr ed2 = *testIns.first;
                    ed2->m_id = eid++;
                    ed2->m_elLink.push_back(
                        pair<ElementSharedPtr,int>(elmt[i],j));
                }
                else
                {
                    EdgeSharedPtr e2 = *(testIns.first);
                    elmt[i]->SetEdge(j, e2);
                    if (e2->m_edgeNodes.size() == 0 &&
                        ed->m_edgeNodes.size() > 0)
                    {
                        e2->m_curveType = ed->m_curveType;
                        e2->m_edgeNodes = ed->m_edgeNodes;

                        // Reverse nodes if appropriate.
                        if (e2->m_n1->m_id != ed->m_n1->m_id)
                        {
                            reverse(e2->m_edgeNodes.begin(),
                                    e2->m_edgeNodes.end());
                        }
                    }

                    if (ed->m_parentCAD)
                    {
                        e2->m_parentCAD = ed->m_parentCAD;
                    }

                    // Update edge to element map.
                    e2->m_elLink.push_back(
                        pair<ElementSharedPtr,int>(elmt[i],j));
                }
            }
        }
    }

    // Create links for 1D elements
    for (int i = 0; i < m_mesh->m_element[1].size(); ++i)
    {
        ElementSharedPtr elmt = m_mesh->m_element[1][i];
        NodeSharedPtr v0 = elmt->GetVertex(0);
        NodeSharedPtr v1 = elmt->GetVertex(1);
        vector<NodeSharedPtr> edgeNodes;
        EdgeSharedPtr E = std::shared_ptr<Edge>(
            new Edge(v0, v1, edgeNodes, elmt->GetConf().m_edgeCurveType));

        EdgeSet::iterator it = m_mesh->m_edgeSet.find(E);
        if (it == m_mesh->m_edgeSet.end())
        {
            cerr << "Cannot find corresponding element edge for "
                 << "1D element " << i << endl;
            abort();
        }
        elmt->SetEdgeLink(*it);

        // Update 2D element boundary map.
        pair<weak_ptr<Element>, int> eMap = (*it)->m_elLink.at(0);
        eMap.first.lock()->SetBoundaryLink(eMap.second, i);

        // Update vertices
        elmt->SetVertex(0, (*it)->m_n1, false);
        elmt->SetVertex(1, (*it)->m_n2, false);

        // Copy curvature to edge.
        if ((*it)->m_edgeNodes.size() > 0)
        {
            ElementSharedPtr edge = elmt;
            if (edge->GetVertex(0) == (*it)->m_n1)
            {
                edge->SetVolumeNodes((*it)->m_edgeNodes);
            }
            elmt->SetCurveType((*it)->m_curveType);
        }
    }
}


/**
 * @brief Create a unique set of mesh faces from elements stored in
 * Mesh::element.
 *
 * All elements are scanned and a unique list of enumerated faces is
 * produced in #m_faceSet. Since elements created their own faces
 * independently, we examine each element only uses face objects from
 * #m_faceSet. Duplicate faces of those in #m_face are replaced with
 * the corresponding entry in #m_faceSet. Finally, we scan the list of
 * elements for 2-D boundary faces which correspond to faces in
 * #m_faceSet. For such elements, we set its faceLink to reference the
 * corresponding face in #m_faceSet.
 *
 * This routine only proceeds if the expansion dimension is 3.
 */
void Module::ProcessFaces(bool ReprocessFaces)
{
    if (m_mesh->m_expDim < 3) return;

    if(ReprocessFaces)
    {
        vector<ElementSharedPtr> &elmt = m_mesh->m_element[m_mesh->m_expDim];

        m_mesh->m_faceSet.clear();

        // Scan all elements and generate list of unique faces
        for (int i = 0, fid = 0; i < elmt.size(); ++i)
        {
            for (int j = 0; j < elmt[i]->GetFaceCount(); ++j)
            {
                pair<FaceSet::iterator,bool> testIns;
                testIns = m_mesh->m_faceSet.insert(elmt[i]->GetFace(j));

                if (testIns.second)
                {
                    (*(testIns.first))->m_id = fid++;
                    (*(testIns.first))->m_elLink.push_back(
                        pair<ElementSharedPtr,int>(elmt[i],j));
                }
                else
                {
                    elmt[i]->SetFace(j,*testIns.first);
                    // Update face to element map.
                    (*(testIns.first))->m_elLink.push_back(
                        pair<ElementSharedPtr,int>(elmt[i],j));
                }
            }
        }
    }

    // Create links for 2D elements
    for (int i = 0; i < m_mesh->m_element[2].size(); ++i)
    {
        ElementSharedPtr elmt = m_mesh->m_element[2][i];
        vector<NodeSharedPtr> vertices = elmt->GetVertexList();
        vector<NodeSharedPtr> faceNodes;
        vector<EdgeSharedPtr> edgeList = elmt->GetEdgeList();
        FaceSharedPtr F = std::shared_ptr<Face>(
            new Face(vertices, faceNodes, edgeList,
                     elmt->GetConf().m_faceCurveType));

        FaceSet::iterator it = m_mesh->m_faceSet.find(F);
        if (it == m_mesh->m_faceSet.end())
        {
            cout << "Cannot find corresponding element face for 2D "
                 << "element " << i << endl;
            abort();
        }

        elmt->SetFaceLink(*it);

        // Set edges/vertices
        for (int j = 0; j < elmt->GetVertexCount(); ++j)
        {
            elmt->SetVertex(j, (*it)->m_vertexList[j], false);
            elmt->SetEdge(j, (*it)->m_edgeList[j], false);
        }

        EdgeSet tmp(edgeList.begin(),edgeList.end());

        for (int j = 0; j < elmt->GetEdgeCount(); ++j)
        {
            EdgeSharedPtr e = elmt->GetEdge(j);
            EdgeSet::iterator f = tmp.find(e);
            if(f != tmp.end())
            {
                if ((*f)->m_parentCAD)
                {
                    e->m_parentCAD = (*f)->m_parentCAD;
                }
            }
        }

        // Update 3D element boundary map.
        for (int j = 0; j < (*it)->m_elLink.size(); ++j)
        {
            pair<weak_ptr<Element>, int> eMap = (*it)->m_elLink.at(j);
            eMap.first.lock()->SetBoundaryLink(eMap.second, i);
        }

        // Copy face curvature
        if ((*it)->m_faceNodes.size() > 0)
        {
            elmt->SetVolumeNodes((*it)->m_faceNodes);
            elmt->SetCurveType((*it)->m_curveType);
        }
    }
}

/**
 * @brief Enumerate elements stored in Mesh::element.
 *
 * For all elements of equal dimension to the mesh dimension, we
 * enumerate sequentially. All other elements in the list should be of
 * lower dimension and have ID set by a corresponding edgeLink or
 * faceLink (as set in #ProcessEdges or #ProcessFaces).
 */
void Module::ProcessElements()
{
    int cnt = 0;
    for (int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); ++i)
    {
        m_mesh->m_element[m_mesh->m_expDim][i]->SetId(cnt++);
    }
}

/**
 * @brief Generate a list of composites (groups of elements) from tag
 * IDs stored in mesh vertices/edges/faces/elements.
 *
 * Each element is assigned to a composite ID by an input module. First
 * we scan the element list and generate a list of composite IDs. We
 * then generate the composite objects and populate them with a second
 * scan through the element list.
 */
void Module::ProcessComposites()
{
    m_mesh->m_composite.clear();

    // For each element, check to see if a composite has been
    // created. If not, create a new composite. Otherwise, add the
    // element to the composite.
    for (int d = 0; d <= m_mesh->m_expDim; ++d)
    {
        vector<ElementSharedPtr> &elmt = m_mesh->m_element[d];

        for (int i = 0; i < elmt.size(); ++i)
        {
            CompositeMap::iterator it;
            unsigned int tagid = elmt[i]->GetTagList()[0];

            it = m_mesh->m_composite.find(tagid);

            if (it == m_mesh->m_composite.end())
            {
                CompositeSharedPtr tmp = std::shared_ptr<Composite>(
                                                    new Composite());
                pair<CompositeMap::iterator, bool> testIns;
                tmp->m_id  = tagid;
                tmp->m_tag = elmt[i]->GetTag();
                if(m_mesh->m_faceLabels.count(tmp->m_id) != 0)
                {
                    tmp->m_label =  m_mesh->m_faceLabels[tmp->m_id];
                }

                testIns  = m_mesh->m_composite.insert(
                    pair<unsigned int, CompositeSharedPtr>(tagid,tmp));
                it       = testIns.first;
            }

            if (elmt[i]->GetTag() != it->second->m_tag)
            {
                cout << "Different types of elements in same composite!" << endl;
                cout << " -> Composite uses " << it->second->m_tag << endl;
                cout << " -> Element uses   " << elmt[i]->GetTag() << endl;
                cout << "Have you specified physical volumes and surfaces?" << endl;
            }
            it->second->m_items.push_back(elmt[i]);
        }
    }
}

/**
 * clear all element link information from mesh entities to be able to reprocess new mesh
 */
void Module::ClearElementLinks()
{
    EdgeSet::iterator eit;

    for(eit = m_mesh->m_edgeSet.begin(); eit != m_mesh->m_edgeSet.end(); eit++)
    {
        (*eit)->m_elLink.clear();
    }

    FaceSet::iterator fit;

    for(fit = m_mesh->m_faceSet.begin(); fit != m_mesh->m_faceSet.end(); fit++)
    {
        (*fit)->m_elLink.clear();
    }
}

/**
 * @brief Reorder node IDs so that prisms and tetrahedra are aligned
 * correctly.
 *
 * Orientation of prism lines (i.e. a large prism which has been split
 * into subprisms) cannot be guaranteed when elements are created
 * one-by-one, or when periodic boundary conditions are used. This
 * routine uses the following strategy:
 *
 *   - Destroy the existing node numbering.
 *   - Detect a line of prisms using the PrismLines routine.
 *   - For each line, renumber node IDs consistently so that highest ID
 *     per-element corresponds to the line of collapsed coordinate
 *     points.
 *   - Recreate each prism in the line using the new ordering, and apply
 *     the existing OrientPrism routine to orient nodes accordingly.
 *   - When all prism lines are processed, recreate all tetrahedra using
 *     the existing orientation code.
 *   - Finally renumber any other nodes (i.e. those belonging to
 *     hexahedra).
 *
 * The last step is to eliminate duplicate edges/faces and reenumerate.
 *
 * NOTE: This routine does not copy face-interior high-order information
 * yet!
 */
void Module::ReorderPrisms(PerMap &perFaces)
{
    // Loop over elements and extract any that are prisms.
    int i, j, k;

    if (m_mesh->m_expDim < 3)
    {
        return;
    }

    map<int, int>    lines;
    set<int>         prismsDone, tetsDone;
    PerMap::iterator pIt;

    // Compile list of prisms and tets.
    for (i = 0; i < m_mesh->m_element[3].size(); ++i)
    {
        ElementSharedPtr el = m_mesh->m_element[3][i];

        if (el->GetConf().m_e == LibUtilities::ePrism)
        {
            prismsDone.insert(i);
        }
        else if (el->GetConf().m_e == LibUtilities::eTetrahedron)
        {
            tetsDone.insert(i);
        }
    }

    // Destroy existing node numbering.
    NodeSet::iterator it;
    for (it = m_mesh->m_vertexSet.begin(); it != m_mesh->m_vertexSet.end(); ++it)
    {
        (*it)->m_id = -1;
    }

    // Counter for new node IDs.
    int nodeId = 0;
    int prismTris[2][3] = {{0,1,4}, {3,2,5}};

    // Warning flag for high-order curvature information.
    bool warnCurvature = false;

    // facesDone tracks face IDs inside prisms which have already been
    // aligned.
    std::unordered_set<int> facesDone;
    std::unordered_set<int>::iterator fIt[2], fIt2;

    // Loop over prisms until we've found all lines of prisms.
    while (prismsDone.size() > 0)
    {
        vector<ElementSharedPtr> line;

        // Call PrismLines to identify all prisms connected to
        // prismDone.begin() and place them in line[].
        PrismLines(*prismsDone.begin(), perFaces, prismsDone, line);

        // Loop over each prism, figure out which line of vertices
        // contains the vertex with highest ID.
        for (i = 0; i < line.size(); ++i)
        {
            // Copy tags and nodes from existing element.
            vector<int>           tags  = line[i]->GetTagList();
            vector<NodeSharedPtr> nodes = line[i]->GetVertexList();

            // See if either face of this prism has been renumbered
            // already.
            FaceSharedPtr f[2] = {
                line[i]->GetFace(1), line[i]->GetFace(3)
            };

            fIt[0] = facesDone.find(f[0]->m_id);
            fIt[1] = facesDone.find(f[1]->m_id);

            // See if either of these faces is periodic. If it is, then
            // assign ids accordingly.
            for (j = 0; j < 2; ++j)
            {
                pIt = perFaces.find(f[j]->m_id);

                if (pIt == perFaces.end())
                {
                    continue;
                }

                fIt2 = facesDone.find(pIt->second.first->m_id);

                if (fIt[j] == facesDone.end() &&
                    fIt2   != facesDone.end())
                {
                    fIt[j] = fIt2;
                }
            }

            if (fIt[0] != facesDone.end() &&
                fIt[1] != facesDone.end())
            {
                // Should not be the case that both faces have already
                // been renumbered.
                ASSERTL0(false, "Renumbering error!");
            }
            else if (fIt[0] == facesDone.end() &&
                     fIt[1] == facesDone.end())
            {
                // Renumber both faces.
                for (j = 0; j < 2; ++j)
                {
                    for (k = 0; k < 3; ++k)
                    {
                        NodeSharedPtr n = nodes[prismTris[j][k]];
                        if (n->m_id == -1)
                        {
                            n->m_id = nodeId++;
                        }
                    }
                }

                facesDone.insert(f[0]->m_id);
                facesDone.insert(f[1]->m_id);
            }
            else
            {
                // Renumber face. t identifies the face not yet
                // numbered, o identifies the other face.
                int t = fIt[0] == facesDone.end() ? 0 : 1;
                int o = (t+1) % 2;
                ASSERTL1(fIt[o] != facesDone.end(),"Renumbering error");

                // Determine which of the three vertices on the 'other'
                // face corresponds to the highest ID - this signifies
                // the singular point of the line of prisms.
                int tmp1[3] = {
                    nodes[prismTris[o][0]]->m_id,
                    nodes[prismTris[o][1]]->m_id,
                    nodes[prismTris[o][2]]->m_id
                };
                int tmp2[3] = {0,1,2};

                if (tmp1[0] > tmp1[1])
                {
                    swap(tmp1[0], tmp1[1]);
                    swap(tmp2[0], tmp2[1]);
                }

                if (tmp1[1] > tmp1[2])
                {
                    swap(tmp1[1], tmp1[2]);
                    swap(tmp2[1], tmp2[2]);
                }

                if (tmp1[0] > tmp1[2])
                {
                    swap(tmp1[0], tmp1[2]);
                    swap(tmp2[0], tmp2[2]);
                }

                // Renumber this face so that highest ID matches.
                for (j = 0; j < 3; ++j)
                {
                    NodeSharedPtr n = nodes[prismTris[t][tmp2[j]]];
                    if (n->m_id == -1)
                    {
                        n->m_id = nodeId++;
                    }
                }

                facesDone.insert(f[t]->m_id);
            }

            for (j = 0; j < 6; ++j)
            {
                ASSERTL1(nodes[j]->m_id != -1, "Renumbering error");
            }

            // Recreate prism with the new ordering.
            ElmtConfig conf(LibUtilities::ePrism, 1, false, false, true);
            ElementSharedPtr el = GetElementFactory().CreateInstance(
                LibUtilities::ePrism, conf, nodes, tags);

            // Now transfer high-order information back into
            // place. TODO: Face curvature.
            for (j = 0; j < 9; ++j)
            {
                EdgeSharedPtr e1 = line[i]->GetEdge(j);
                for (k = 0; k < 9; ++k)
                {
                    EdgeSharedPtr e2 = el->GetEdge(k);
                    if (e1->m_n1 == e2->m_n1 && e1->m_n2 == e2->m_n2)
                    {
                        e2->m_edgeNodes = e1->m_edgeNodes;
                    }
                    else if (e1->m_n1 == e2->m_n1 && e1->m_n2 == e2->m_n2)
                    {
                        e2->m_edgeNodes = e1->m_edgeNodes;
                        std::reverse(e2->m_edgeNodes.begin(),
                                     e2->m_edgeNodes.end());
                    }
                }
            }

            // Warn users that we're throwing away face curvature
            if (!warnCurvature)
            {
                for (j = 0; j < 5; ++j)
                {
                    if (line[i]->GetFace(j)->m_faceNodes.size() > 0)
                    {
                        warnCurvature = true;
                        break;
                    }
                }
            }

            // Replace old prism.
            m_mesh->m_element[3][line[i]->GetId()] = el;
        }
    }

    if (warnCurvature)
    {
        cerr << "[ReorderPrisms] WARNING: Face curvature detected in "
             << "some prisms; this will be ignored in further module "
             << "evaluations."
             << endl;
    }

    // Loop over periodic faces, enumerate vertices.
    for (pIt = perFaces.begin(); pIt != perFaces.end(); ++pIt)
    {
        FaceSharedPtr f2       = pIt->second.first;
        FaceSharedPtr f1       = perFaces[f2->m_id].first;
        vector<int>   perVerts = pIt->second.second;
        int           nVerts   = perVerts.size();

        // Number periodic vertices first.
        for (j = 0; j < nVerts; ++j)
        {
            NodeSharedPtr n1 = f1->m_vertexList[j];
            NodeSharedPtr n2 = f2->m_vertexList[perVerts[j]];

            if (n1->m_id == -1 && n2->m_id == -1)
            {
                n1->m_id = nodeId++;
                n2->m_id = nodeId++;
            }
            else if (n1->m_id != -1 && n2->m_id != -1)
            {
                continue;
            }
            else
            {
                ASSERTL0(false, "Periodic face renumbering error");
            }
        }
    }

    // Recreate tets.
    set<int>::iterator it2;
    for (it2 = tetsDone.begin(); it2 != tetsDone.end(); ++it2)
    {
        ElementSharedPtr el = m_mesh->m_element[3][*it2];
        vector<NodeSharedPtr> nodes = el->GetVertexList();
        vector<int> tags = el->GetTagList();

        for (i = 0; i < 4; ++i)
        {
            if (nodes[i]->m_id == -1)
            {
                nodes[i]->m_id = nodeId++;
            }
        }

        // Recreate tet.
        ElmtConfig conf(LibUtilities::eTetrahedron, 1, false, false, true);
        m_mesh->m_element[3][*it2] = GetElementFactory().CreateInstance(
            LibUtilities::eTetrahedron, conf, nodes, tags);
    }

    // Enumerate rest of vertices.
    for (it = m_mesh->m_vertexSet.begin(); it != m_mesh->m_vertexSet.end(); ++it)
    {
        if ((*it)->m_id == -1)
        {
            (*it)->m_id = nodeId++;
        }
    }

    ProcessEdges   ();
    ProcessFaces   ();
    ProcessElements();
}

void Module::PrismLines(int                       prism,
                        PerMap                   &perFaces,
                        set<int>                 &prismsDone,
                        vector<ElementSharedPtr> &line)
{
    int                i;
    set<int>::iterator it = prismsDone.find(prism);
    PerMap::iterator   it2;

    if (it == prismsDone.end())
    {
        return;
    }

    // Remove this prism from the list.
    prismsDone.erase(it);
    line.push_back(m_mesh->m_element[3][prism]);

    // Now find prisms connected to this one through a triangular face.
    for (i = 1; i <= 3; i += 2)
    {
        FaceSharedPtr f = m_mesh->m_element[3][prism]->GetFace(i);
        int nextId;

        // See if this face is periodic.
        it2 = perFaces.find(f->m_id);

        if (it2 != perFaces.end())
        {
            int id2 = it2->second.first->m_id;
            nextId  = it2->second.first->m_elLink[0].first.lock()->GetId();
            perFaces.erase(it2);
            perFaces.erase(id2);
            PrismLines(nextId, perFaces, prismsDone, line);
        }

        // Nothing else connected to this face.
        if (f->m_elLink.size() == 1)
        {
            continue;
        }

        nextId = f->m_elLink[0].first.lock()->GetId();
        if (nextId == m_mesh->m_element[3][prism]->GetId())
        {
            nextId = f->m_elLink[1].first.lock()->GetId();
        }

        PrismLines(nextId, perFaces, prismsDone, line);
    }
}

/**
 * @brief Register a configuration option with a module.
 */
void Module::RegisterConfig(string key, string val)
{
    map<string, ConfigOption>::iterator it = m_config.find(key);
    if (it == m_config.end())
    {
        cerr << "WARNING: Unrecognised config option " << key
             << ", proceeding anyway." << endl;
    }

    it->second.beenSet = true;

    if (it->second.isBool)
    {
        it->second.value = "1";
    }
    else
    {
        if(val.size() == 0)
        {
            it->second.value = it->second.defValue;
        }
        else
        {
            it->second.value = val;
        }
    }
}

/**
 * @brief Print out all configuration options for a module.
 */
void Module::PrintConfig()
{
    map<string, ConfigOption>::iterator it;

    if (m_config.size() == 0)
    {
        cerr << "No configuration options for this module." << endl;
        return;
    }

    for (it = m_config.begin(); it != m_config.end(); ++it)
    {
        cerr << setw(10) << it->first << ": " << it->second.desc
             << endl;
    }
}

/**
 * @brief Sets default configuration options for those which have not
 * been set.
 */
void Module::SetDefaults()
{
    map<string, ConfigOption>::iterator it;

    for (it = m_config.begin(); it != m_config.end(); ++it)
    {
        if (!it->second.beenSet)
        {
            it->second.value = it->second.defValue;
        }
    }
}

/**
 * @brief Print a brief summary of information.
 */
void InputModule::PrintSummary()
{
    // Compute the number of full-dimensional elements and boundary
    // elements.
    cerr << "Expansion dimension is " << m_mesh->m_expDim << endl;
    cerr << "Space dimension is " << m_mesh->m_spaceDim << endl;
    cerr << "Read " << m_mesh->m_node.size() << " nodes" << endl;
    cerr << "Read " << m_mesh->GetNumElements() << " "
         << m_mesh->m_expDim << "-D elements" << endl;
    cerr << "Read " << m_mesh->GetNumBndryElements()
         << " boundary elements" << endl;
}

}
}
