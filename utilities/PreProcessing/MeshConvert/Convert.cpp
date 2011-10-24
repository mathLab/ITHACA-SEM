////////////////////////////////////////////////////////////////////////////////
//
//  File: Convert.cpp
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
//  Description: Mesh converter base class and XML writer.
//
////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <algorithm>
using namespace std;

#include "Convert.h"

namespace Nektar
{
    namespace Utilities
    {
        ConvertFactory& GetConvertFactory()
        {
            typedef Loki::SingletonHolder<ConvertFactory,
                Loki::CreateUsingNew,
                Loki::NoDestroy > Type;
            return Type::Instance();
        }

        /**
         * Each element is processed in turn and the vertices extracted and
         * inserted into #m_vertexSet, which at the end of the routine
         * contains all unique vertices in the mesh.
         */
        void Convert::ProcessVertices()
        {
            for (int i = 0, vid = 0; i < m_element.size(); ++i)
            {
                if (m_element[i]->GetDim() == m_expDim)
                {
                    for (int j = 0; j < m_element[i]->GetVertexCount(); ++j)
                    {
                        pair<NodeSet::iterator,bool> insTest =
                            m_vertexSet.insert(m_element[i]->GetVertex(j));
                        if (insTest.second)
                        {
                            (*(insTest.first))->id = vid++;
                        }
                    }
                }
            }
        }


        /**
         * This routine only proceeds if the expansion dimension is 2 or 3.
         *
         * All elements are first scanned and a list of unique, enumerated
         * edges produced in #m_edgeSet. Since each element generated its
         * edges independently, we must now ensure that each element only uses
         * edge objects from the #m_edgeSet set This ensures there are no
         * duplicate edge objects. Finally, we scan the list of elements for
         * 1-D boundary elements which correspond to an edge in
         * #m_edgeSet. For such elements, we set its edgeLink to reference the
         * corresponding edge in #m_edgeSet.
         */
        void Convert::ProcessEdges()
        {
            if (m_expDim < 2) return;

            // Scan all elements and generate list of unique edges
            for (int i = 0, eid = 0; i < m_element.size(); ++i)
            {
                if (m_element[i]->GetDim() == m_expDim)
                {
                    for (int j = 0; j < m_element[i]->GetEdgeCount(); ++j)
                    {
                        pair<EdgeSet::iterator,bool> testIns;
                        testIns = m_edgeSet.insert(m_element[i]->GetEdge(j));

                        if (testIns.second == false)
                        {
                            m_element[i]->SetEdge(j,*testIns.first);
                        }
                        else
                        {
                            (*(testIns.first))->id = eid++;
                        }
                    }
                }
            }
            // Create links for 1D elements
            for (int i = 0; i < m_element.size(); ++i)
            {
                if (m_element[i]->GetDim() == 1)
                {
                    NodeSharedPtr v0 = m_element[i]->GetVertex(0);
                    NodeSharedPtr v1 = m_element[i]->GetVertex(1);
                    vector<NodeSharedPtr> edgeNodes;
                    EdgeSharedPtr E = boost::shared_ptr<Edge>(new Edge(v0, v1, edgeNodes));
                    EdgeSet::iterator it = m_edgeSet.find(E);
                    if (it == m_edgeSet.end())
                    {
                        cout << "Cannot find corresponding element face for 1D element " << i << endl;
                        abort();
                    }
                    m_element[i]->SetEdgeLink(*it);
                }
            }
        }


        /**
         * This routine only proceeds if the expansion dimension is 3.
         *
         * All elements are scanned and a unique list of enumerated faces is
         * produced in #m_faceSet. Since elements created their own faces
         * independently, we examine each element only uses face objects from
         * #m_faceSet. Duplicate faces of those in #m_face are replaced with
         * the corresponding entry in #m_faceSet. Finally, we scan the list of
         * elements for 2-D boundary faces which correspond to faces in
         * #m_faceSet. For such elements, we set its faceLink to reference the
         * corresponding face in #m_faceSet.
         */
        void Convert::ProcessFaces()
        {
            if (m_expDim < 3) return;

            // Scan all elements and generate list of unique faces
            for (int i = 0, fid = 0; i < m_element.size(); ++i)
            {
                if (m_element[i]->GetDim() == m_expDim)
                {
                    for (int j = 0; j < m_element[i]->GetFaceCount(); ++j)
                    {
                        pair<FaceSet::iterator,bool> testIns;
                        testIns = m_faceSet.insert(m_element[i]->GetFace(j));
                        
                        if (testIns.second == false)
                        {
                            m_element[i]->SetFace(j,*testIns.first);
                        }
                        else
                        {
                            (*(testIns.first))->id = fid++;
                        }
                    }
                }
            }
            // Create links for 2D elements
            for (int i = 0; i < m_element.size(); ++i)
            {
                if (m_element[i]->GetDim() == 2)
                {
                    vector<NodeSharedPtr> vertices = m_element[i]->GetVertexList();
                    vector<NodeSharedPtr> faceNodes;
                    vector<EdgeSharedPtr> edgeList = m_element[i]->GetEdgeList();
                    FaceSharedPtr F = boost::shared_ptr<Face>(new Face(vertices, faceNodes, edgeList));
                    FaceSet::iterator it = m_faceSet.find(F);
                    if (it == m_faceSet.end())
                    {
                        cout << "Cannot find corresponding element face for 2D element " << i << endl;
                        abort();
                    }
                    m_element[i]->SetFaceLink(*it);
                }
            }
        }


        /**
         * For all elements of equal dimension to the mesh dimension, we
         * enumerate sequentially. All other elements in the list should be of
         * lower dimension and have ID set by a corresponding edgeLink or
         * faceLink (as set in #ProcessEdges or #ProcessFaces).
         */
        void Convert::ProcessElements()
        {
            int cnt = 0;
            for (int i = 0; i < m_element.size(); ++i)
            {
                if (m_element[i]->GetDim() == m_expDim)
                {
                    m_element[i]->SetId(cnt++);
                }
            }
        }


        /**
         * Each element is assigned to a composite ID by Gmsh. First we scan
         * the element list and generate a list of composite IDs. We then
         * generate the composite objects and populate them with a second scan
         * through the element list.
         */
        void Convert::ProcessComposites()
        {
            int p = 0;
            vector<ElementSharedPtr>::iterator it;
            list<int> compIdList;
            list<int>::iterator it2;

            // Generate sorted list of unique composite IDs to create.
            for (it = m_element.begin(); it != m_element.end(); ++it)
            {
                compIdList.push_back((*it)->GetTagList()[0]);
            }
            compIdList.sort();
            compIdList.unique();

            // Create a composite object for each unique composite ID.
            m_composite.resize(compIdList.size());
            it2 = compIdList.begin();
            for (int i = 0; i < m_composite.size(); ++i, ++it2)
            {
                m_composite[i] = boost::shared_ptr<Composite>(new Composite);
                m_composite[i]->id = *it2;
                m_composite[i]->tag = "";
            }

            // Populate composites with elements.
            for (int i = 0; i < m_element.size(); ++i)
            {
                int p = m_element[i]->GetTagList()[0];
                for (int j = 0; j < m_composite.size(); ++j)
                {
                    if (m_composite[j]->id == p)
                    {
                        if (m_composite[j]->tag == "")
                        {
                            // This is the first element in this composite
                            // so assign the composite tag.
                            m_composite[j]->tag = m_element[i]->GetTag();
                        }
                        else
                        {
                            // Otherwise, check element tag matches composite.
                            if (m_element[i]->GetTag() != m_composite[j]->tag)
                            {
                                cout << "Different types of elements in same composite!" << endl;
                                cout << " -> Composite uses " << m_composite[j]->tag << endl;
                                cout << " -> Element uses   " << m_element[i]->GetTag() << endl;
                                cout << "Have you specified physical volumes and surfaces?" << endl;
                            }
                        }
                        m_composite[j]->items.push_back(m_element[i]);
                        break;
                    }
                }
            }
        }


        /**
         * Writes the mesh information stored in the data structures out to an
         * XML file in the native Nektar++ XML file format.
         */
        void Convert::WriteXmlFile(const string pFilename)
        {
            TiXmlDocument doc;
            TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "utf-8", "");
            doc.LinkEndChild( decl );

            TiXmlElement * root = new TiXmlElement( "NEKTAR" );
            doc.LinkEndChild( root );

            // Begin <GEOMETRY> section
            TiXmlElement * geomTag = new TiXmlElement( "GEOMETRY" );
            geomTag->SetAttribute("DIM", m_expDim);
            geomTag->SetAttribute("SPACE", m_spaceDim);
            root->LinkEndChild( geomTag );

            WriteXmlNodes(geomTag);

            WriteXmlEdges(geomTag);

            WriteXmlFaces(geomTag);

            WriteXmlElements(geomTag);

            WriteXmlCurves(geomTag);

            WriteXmlComposites(geomTag);

            WriteXmlDomain(geomTag);

            WriteXmlExpansions(root);

            WriteXmlConditions(root);

            // Save the XML file.
            doc.SaveFile( pFilename );
        }

        void Convert::WriteXmlNodes(TiXmlElement * pRoot)
        {
            TiXmlElement* verTag = new TiXmlElement( "VERTEX" );
            NodeSet::iterator it;

            for (it = m_vertexSet.begin(); it != m_vertexSet.end(); ++it)
            {
                NodeSharedPtr n = *it;
                stringstream s;
                s << scientific << setprecision(3) 
                  << n->x << " " << n->y << " " << n->z;
                TiXmlElement * v = new TiXmlElement( "V" );
                v->SetAttribute("ID",n->id);
                v->LinkEndChild(new TiXmlText(s.str()));
                verTag->LinkEndChild(v);
            }
            pRoot->LinkEndChild(verTag);
        }

        void Convert::WriteXmlEdges(TiXmlElement * pRoot)
        {
            if (m_expDim >= 2)
            {
                int edgecnt = 0;
                TiXmlElement* verTag = new TiXmlElement( "EDGE" );
                EdgeSet::iterator it;
                
                for (it = m_edgeSet.begin(); it != m_edgeSet.end(); ++it)
                {
                    EdgeSharedPtr ed = *it;
                    stringstream s;

                    s << setw(5) << ed->n1->id << "  " << ed->n2->id << "   ";
                    TiXmlElement * e = new TiXmlElement( "E" );
                    e->SetAttribute("ID",ed->id);
                    e->LinkEndChild( new TiXmlText(s.str()) );
                    verTag->LinkEndChild(e);
                }
                pRoot->LinkEndChild( verTag );
            }
        }

        void Convert::WriteXmlFaces(TiXmlElement * pRoot)
        {
            if (m_expDim == 3)
            {
                TiXmlElement* verTag = new TiXmlElement( "FACE" );
                FaceSet::iterator it;

                for (it = m_faceSet.begin(); it != m_faceSet.end(); ++it)
                {
                    stringstream s;
                    FaceSharedPtr fa = *it;

                    for (int j = 0; j < fa->edgeList.size(); ++j)
                    {
                        s << setw(10) << fa->edgeList[j]->id;
                    }
                    TiXmlElement * f;
                    switch(fa->vertexList.size())
                    {
                        case 3:
                            f = new TiXmlElement("T");
                            break;
                        case 4:
                            f = new TiXmlElement("Q");
                            break;
                        default:
                            abort();
                    }
                    f->SetAttribute("ID", fa->id);
                    f->LinkEndChild( new TiXmlText(s.str()));
                    verTag->LinkEndChild(f);
                }
                pRoot->LinkEndChild( verTag );
            }
        }

        void Convert::WriteXmlElements(TiXmlElement * pRoot)
        {
            TiXmlElement* verTag = new TiXmlElement( "ELEMENT" );

            for(int i=0; i<m_element.size(); ++i)
            {
                if (m_element[i]->GetDim() == m_expDim)
                {
                    TiXmlElement *elm_tag = new TiXmlElement(m_element[i]->GetTag());
                    elm_tag->SetAttribute("ID", m_element[i]->GetId());
                    elm_tag->LinkEndChild( new TiXmlText(m_element[i]->GetXmlString()) );
                    verTag->LinkEndChild(elm_tag);
                }
            }
            pRoot->LinkEndChild( verTag );
        }

        void Convert::WriteXmlCurves(TiXmlElement * pRoot)
        {
            int edgecnt = 0;
            int facecnt = 0;

            bool curve = false;
            EdgeSet::iterator it;
            for (it = m_edgeSet.begin(); it != m_edgeSet.end(); ++it)
            {
                if ((*it)->edgeNodes.size() > 0) 
                {
                    curve = true;
                    break;
                }
            }
            if (!curve) return;

            TiXmlElement * curved = new TiXmlElement ("CURVED" );

            for (it = m_edgeSet.begin(); it != m_edgeSet.end(); ++it)
            {
                if ((*it)->edgeNodes.size() > 0)
                {
                    TiXmlElement * e = new TiXmlElement( "E" );
                    e->SetAttribute("ID", edgecnt++);
                    e->SetAttribute("EDGEID", (*it)->id);
                    e->SetAttribute("TYPE", "PolyEvenlySpaced");
                    e->SetAttribute("NUMPOINTS", (*it)->GetNodeCount());
                    TiXmlText * t0 = new TiXmlText((*it)->GetXmlCurveString());
                    e->LinkEndChild(t0);
                    curved->LinkEndChild(e);
                }
            }

            FaceSet::iterator it2;
            
            for (it2 = m_faceSet.begin(); it2 != m_faceSet.end(); ++it2)
            {
                if ((*it2)->faceNodes.size() > 0)
                {
                    TiXmlElement * f = new TiXmlElement( "F" );
                    f->SetAttribute("ID",facecnt++);
                    f->SetAttribute("FACEID",(*it2)->id);
                    f->SetAttribute("TYPE","PolyEvenlySpaced");
                    f->SetAttribute("NUMPOINTS",(*it2)->GetNodeCount());
                    TiXmlText * t0 = new TiXmlText((*it2)->GetXmlCurveString());
                    f->LinkEndChild(t0);
                    curved->LinkEndChild(f);
                }
            }
            pRoot->LinkEndChild( curved );
        }

        void Convert::WriteXmlComposites(TiXmlElement * pRoot)
        {
            TiXmlElement* verTag = new TiXmlElement("COMPOSITE");

            for(int i = 0; i < m_composite.size(); ++i){
                if (m_composite[i]->items.size() > 0) {
                    TiXmlElement *comp_tag = new TiXmlElement("C"); // Composite
                    TiXmlElement *elm_tag;

                    comp_tag->SetAttribute("ID", m_composite[i]->id);
                    comp_tag->LinkEndChild( new TiXmlText(m_composite[i]->GetXmlString()) );
                    verTag->LinkEndChild(comp_tag);
                }
                else
                {
                    cout << "Composite " << m_composite[i]->id << " contains nothing." << endl;
                }
            }

            pRoot->LinkEndChild( verTag );
        }

        void Convert::WriteXmlDomain(TiXmlElement * pRoot)
        {
            // Write the <DOMAIN> subsection.
            TiXmlElement * domain = new TiXmlElement ("DOMAIN" );
            std::string list;
            for (int i = 0; i < m_composite.size(); ++i)
            {
                if (m_composite[i]->items[0]->GetDim() == m_expDim)
                {
                    if (list.length() > 0)
                    {
                        list += ",";
                    }
                    list += boost::lexical_cast<std::string>(m_composite[i]->id);
                }
            }
            domain->LinkEndChild( new TiXmlText(" C[" + list + "] "));
            pRoot->LinkEndChild( domain );
        }

        void Convert::WriteXmlExpansions(TiXmlElement * pRoot)
        {
            // Write a default <EXPANSIONS> section.
            TiXmlElement * expansions = new TiXmlElement ("EXPANSIONS");
            for (int i = 0; i < m_composite.size(); ++i)
            {
                if (m_composite[i]->items[0]->GetDim() == m_expDim)
                {
                    TiXmlElement * exp = new TiXmlElement ( "E");
                    exp->SetAttribute("COMPOSITE", "C["
                        + boost::lexical_cast<std::string>(m_composite[i]->id)
                        + "]");
                    exp->SetAttribute("NUMMODES",7);
                    exp->SetAttribute("FIELDS","u");
                    exp->SetAttribute("TYPE","MODIFIED");
                    expansions->LinkEndChild(exp);
                }
            }
            pRoot->LinkEndChild(expansions);
        }
        
        void Convert::WriteXmlConditions(TiXmlElement * pRoot)
        {
            TiXmlElement * conditions = new TiXmlElement ("CONDITIONS");
            pRoot->LinkEndChild(conditions);
        }   
    }
}
