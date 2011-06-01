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

            // Write the <DOMAIN> subsection.
            TiXmlElement * domain = new TiXmlElement ("DOMAIN" );
            domain->LinkEndChild( new TiXmlText( " C[0] " ));
            geomTag->LinkEndChild( domain );

            // Write a default <EXPANSIONS> section.
            TiXmlElement * expansions = new TiXmlElement ("EXPANSIONS");
            TiXmlElement * exp1 = new TiXmlElement ( "E");
            exp1->SetAttribute("COMPOSITE","C[0]");
            exp1->SetAttribute("NUMMODES",8);
            exp1->SetAttribute("TYPE","MODIFIED");
            expansions->LinkEndChild(exp1);
            root->LinkEndChild(expansions);

            // Save the XML file.
            doc.SaveFile( pFilename );
        }

        void Convert::WriteXmlNodes(TiXmlElement * pRoot)
        {
            TiXmlElement* verTag = new TiXmlElement( "VERTEX" );

            sort(m_vertex.begin(), m_vertex.end(), shared_ptr_less_than<Node>());
            for( int i = 0; i < m_vertex.size(); ++i )
            {
                stringstream s;
                s << scientific << setprecision(3) <<  m_vertex[i]->x << " "
                  << m_vertex[i]->y << " " << m_vertex[i]->z;
                TiXmlElement * v = new TiXmlElement( "V" );
                v->SetAttribute("ID",m_vertex[i]->id);
                v->LinkEndChild( new TiXmlText(s.str()) );
                verTag->LinkEndChild(v);
            }
            pRoot->LinkEndChild( verTag );
        }

        void Convert::WriteXmlEdges(TiXmlElement * pRoot)
        {
            if (m_expDim >= 2)
            {
                int edgecnt = 0;
                int small, large;
                TiXmlElement* verTag = new TiXmlElement( "EDGE" );

                for( int i = 0; i < m_edge.size(); ++i )
                {
                    stringstream s;

                    s << setw(5) << m_edge[i]->n1->id << "  " << m_edge[i]->n2->id << "   ";
                    TiXmlElement * e = new TiXmlElement( "E" );
                    e->SetAttribute("ID",m_edge[i]->id);
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
                int small, large;
                TiXmlElement* verTag = new TiXmlElement( "FACE" );

                for( int i = 0; i < m_face.size(); ++i )
                {
                    stringstream s;

                    for ( int j = 0; j < m_face[i]->edgeList.size(); ++j)
                    {
                        s << setw(5) << m_face[i]->edgeList[j]->id;
                    }
                    TiXmlElement * f;
                    if (m_face[i]->vertexList.size() == 3)
                    {
                        f = new TiXmlElement( "T" );
                    }
                    else
                    {
                        f = new TiXmlElement( "Q" );
                    }
                    f->SetAttribute("ID", m_face[i]->id);
                    f->LinkEndChild( new TiXmlText(s.str()) );
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
            for (int i = 0; i < m_edge.size(); ++i)
            {
                if (m_edge[i]->edgeNodes.size() > 0) {
                    curve = true;
                }
            }
            if (!curve) return;

            TiXmlElement * curved = new TiXmlElement ("CURVED" );

            for( int i = 0; i < m_edge.size(); ++i )
            {
                if (m_edge[i]->edgeNodes.size() > 0)
                {
                    TiXmlElement * e = new TiXmlElement( "E" );
                    e->SetAttribute("ID",edgecnt++);
                    e->SetAttribute("EDGEID",m_edge[i]->id);
                    e->SetAttribute("TYPE","PolyEvenlySpaced");
                    e->SetAttribute("NUMPOINTS",m_edge[i]->GetNodeCount());
                    TiXmlText * t0 = new TiXmlText(m_edge[i]->GetXmlCurveString());
                    e->LinkEndChild(t0);
                    curved->LinkEndChild(e);
                }
            }

            for( int i = 0; i < m_face.size(); ++i)
            {
                if (m_face[i]->faceNodes.size() > 0)
                {
                    TiXmlElement * f = new TiXmlElement( "F" );
                    f->SetAttribute("ID",facecnt++);
                    f->SetAttribute("FACEID",m_face[i]->id);
                    f->SetAttribute("TYPE","PolyEvenlySpaced");
                    f->SetAttribute("NUMPOINTS",m_face[i]->GetNodeCount());
                    TiXmlText * t0 = new TiXmlText(m_face[i]->GetXmlCurveString());
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

        int Convert::FindNodeIndex(const NodeSharedPtr pSrc)
        {
            vector<NodeSharedPtr>::iterator it;
            it = find(m_node.begin(), m_node.end(), pSrc);
            return (it != m_node.end() ? it - m_node.begin() : -1);
        }

        int Convert::FindVertexIndex(const NodeSharedPtr pSrc)
        {
            vector<NodeSharedPtr>::iterator it;
            it = find(m_vertex.begin(), m_vertex.end(), pSrc);
            return (it != m_vertex.end() ? it - m_vertex.begin() : -1);
        }

        int Convert::FindEdgeIndex(const EdgeSharedPtr pSrc)
        {
            vector<EdgeSharedPtr>::iterator it;
            for (it = m_edge.begin(); it != m_edge.end(); ++it)
            {
                if (**it == *pSrc) return (it - m_edge.begin());
            }
            return -1;
        }

        int Convert::FindFaceIndex(const FaceSharedPtr pSrc)
        {
            vector<FaceSharedPtr>::iterator it;
            for (it = m_face.begin(); it != m_face.end(); ++it)
            {
                if (**it == *pSrc) return (it - m_face.begin());
            }
            return -1;
        }
    }
}
