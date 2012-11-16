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
//  Description: Abstract input/output modules.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include "Module.h"

using namespace std;

namespace Nektar
{
    namespace Utilities
    {
        /**
         * Returns an instance of the module factory, held as a singleton.
         */
        ModuleFactory& GetModuleFactory()
        {
            typedef Loki::SingletonHolder<ModuleFactory,
                Loki::CreateUsingNew,
                Loki::NoDestroy > Type;
            return Type::Instance();
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
            config["infile"] = ConfigOption(false, "", "Input filename.");
        }
        
        OutputModule::OutputModule(MeshSharedPtr m) : Module(m)
        {
            config["outfile"] = ConfigOption(false, "", "Output filename.");
        }

        /**
         * @brief Open a file for input.
         */
        void InputModule::OpenStream()
        {
            string fname = config["infile"].as<string>();
            mshFile.open(fname.c_str());
            if (!mshFile.good())
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
            string fname = config["outfile"].as<string>();
            mshFile.open(fname.c_str());
            if (!mshFile.good())
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
            vector<ElementSharedPtr> &elmt = m->element[m->expDim];
            
            m->vertexSet.clear();
            
            for (int i = 0, vid = 0; i < elmt.size(); ++i)
            {
                for (int j = 0; j < elmt[i]->GetVertexCount(); ++j)
                {
                    pair<NodeSet::iterator,bool> testIns =
                        m->vertexSet.insert(elmt[i]->GetVertex(j));
                    if (testIns.second)
                    {
                        (*(testIns.first))->id = vid++;
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
        void Module::ProcessEdges()
        {
            if (m->expDim < 2) return;

            vector<ElementSharedPtr> &elmt = m->element[m->expDim];

            m->edgeSet.clear();
            
            // Scan all elements and generate list of unique edges
            for (int i = 0, eid = 0; i < elmt.size(); ++i)
            {
                for (int j = 0; j < elmt[i]->GetEdgeCount(); ++j)
                {
                    pair<EdgeSet::iterator,bool> testIns;
                    EdgeSharedPtr ed = elmt[i]->GetEdge(j);
                    testIns = m->edgeSet.insert(ed);
                    
                    if (testIns.second)
                    {
                        (*(testIns.first))->id = eid++;
                    }
                    else
                    {
                        EdgeSharedPtr e2 = *(testIns.first);
                        elmt[i]->SetEdge(j, e2);
                        if (e2->edgeNodes.size() == 0 && 
                            ed->edgeNodes.size() > 0)
                        {
                            e2->curveType = ed->curveType;
                            e2->edgeNodes = ed->edgeNodes;
                            
                            // Reverse nodes if appropriate.
                            if (e2->n1->id != ed->n1->id)
                            {
                                reverse(e2->edgeNodes.begin(),
                                        e2->edgeNodes.end());
                            }
                        }
                    }
                }
            }
            
            // Create links for 1D elements
            for (int i = 0; i < m->element[1].size(); ++i)
            {
                NodeSharedPtr v0 = m->element[1][i]->GetVertex(0);
                NodeSharedPtr v1 = m->element[1][i]->GetVertex(1);
                vector<NodeSharedPtr> edgeNodes;
                EdgeSharedPtr E = boost::shared_ptr<Edge>(
                    new Edge(v0, v1, edgeNodes, 
                             m->element[1][i]->GetConf().edgeCurveType));
                EdgeSet::iterator it = m->edgeSet.find(E);
                if (it == m->edgeSet.end())
                {
                    cerr << "Cannot find corresponding element face for "
                         << "1D element " << i << endl;
                    abort();
                }
                m->element[1][i]->SetEdgeLink(*it);
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
        void Module::ProcessFaces()
        {
            if (m->expDim < 3) return;

            vector<ElementSharedPtr> &elmt = m->element[m->expDim];
            
            m->faceSet.clear();
            
            // Scan all elements and generate list of unique faces
            for (int i = 0, fid = 0; i < elmt.size(); ++i)
            {
                for (int j = 0; j < elmt[i]->GetFaceCount(); ++j)
                {
                    pair<FaceSet::iterator,bool> testIns;
                    testIns = m->faceSet.insert(elmt[i]->GetFace(j));
                    
                    if (testIns.second)
                    {
                        (*(testIns.first))->id = fid++;
                    }
                    else
                    {
                        elmt[i]->SetFace(j,*testIns.first);
                        // Update face to element map.
                        (*(testIns.first))->elLink.push_back(
                            pair<ElementSharedPtr,int>(elmt[i],j));
                    }
                }
            }
            
            // Create links for 2D elements
            for (int i = 0; i < m->element[2].size(); ++i)
            {
                vector<NodeSharedPtr> vertices = m->element[2][i]->GetVertexList();
                vector<NodeSharedPtr> faceNodes;
                vector<EdgeSharedPtr> edgeList = m->element[2][i]->GetEdgeList();
                FaceSharedPtr F = boost::shared_ptr<Face>(
                    new Face(vertices, faceNodes, edgeList, 
                             m->element[2][i]->GetConf().faceCurveType));
                FaceSet::iterator it = m->faceSet.find(F);
                if (it == m->faceSet.end())
                {
                    cout << "Cannot find corresponding element face for 2D "
                         << "element " << i << endl;
                    abort();
                }
                m->element[2][i]->SetFaceLink(*it);
                
                // Update 3D element boundary map.
                ASSERTL0((*it)->elLink.size() != 0,
                         "Empty boundary map!");
                ASSERTL0((*it)->elLink.size() == 1,
                         "Too many elements in boundary map!");
                pair<ElementSharedPtr, int> eMap = (*it)->elLink.at(0);
                eMap.first->SetBoundaryLink(eMap.second, i);
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
            for (int i = 0; i < m->element[m->expDim].size(); ++i)
            {
                m->element[m->expDim][i]->SetId(cnt++);
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
            m->composite.clear();
            
            // For each element, check to see if a composite has been
            // created. If not, create a new composite. Otherwise, add the
            // element to the composite.
            for (int d = 0; d <= m->expDim; ++d)
            {
                vector<ElementSharedPtr> &elmt = m->element[d];
                
                for (int i = 0; i < elmt.size(); ++i)
                {
                    CompositeMap::iterator it;
                    unsigned int tagid = elmt[i]->GetTagList()[0];
                    
                    it = m->composite.find(tagid);
                    
                    if (it == m->composite.end())
                    {
                        CompositeSharedPtr tmp = boost::shared_ptr<Composite>(
                            new Composite);
                        pair<CompositeMap::iterator, bool> testIns;
                        tmp->id  = tagid;
                        tmp->tag = elmt[i]->GetTag();
                        testIns  = m->composite.insert(
                            pair<unsigned int, CompositeSharedPtr>(tagid,tmp));
                        it       = testIns.first;
                    }

                    if (elmt[i]->GetTag() != it->second->tag)
                    {
                        cout << "Different types of elements in same composite!" << endl;
                        cout << " -> Composite uses " << it->second->tag << endl;
                        cout << " -> Element uses   " << elmt[i]->GetTag() << endl;
                        cout << "Have you specified physical volumes and surfaces?" << endl;
                    }
                    it->second->items.push_back(elmt[i]);
                }
            }
        }
        
        /**
         * @brief Register a configuration option with a module.
         */
        void Module::RegisterConfig(string key, string val)
        {
            map<string, ConfigOption>::iterator it = config.find(key);
            if (it == config.end())
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
                it->second.value = val;
            }
        }
        
        /**
         * @brief Print out all configuration options for a module.
         */
        void Module::PrintConfig()
        {
            map<string, ConfigOption>::iterator it;
            
            if (config.size() == 0)
            {
                cerr << "No configuration options for this module." << endl;
                return;
            }
            
            for (it = config.begin(); it != config.end(); ++it)
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
            
            for (it = config.begin(); it != config.end(); ++it)
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
            cerr << "Expansion dimension is " << m->expDim << endl;
            cerr << "Space dimension is " << m->spaceDim << endl;
            cerr << "Read " << m->node.size() << " nodes" << endl;
            cerr << "Read " << m->GetNumElements() << " " 
                 << m->expDim << "-D elements" << endl;
            cerr << "Read " << m->GetNumBndryElements() 
                 << " boundary elements" << endl;
        }
    }
}
