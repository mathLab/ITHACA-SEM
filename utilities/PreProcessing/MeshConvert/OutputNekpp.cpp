////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputNekpp.cpp
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
//  Description: Nektar++ file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
using namespace std;

#include "MeshElements.h"
#include "OutputNekpp.h"

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey OutputNekpp::className = 
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eOutputModule, "xml"), OutputNekpp::create,
                "Writes a Nektar++ xml file.");

        OutputNekpp::OutputNekpp(MeshSharedPtr m) : OutputModule(m)
        {
            
        }

        OutputNekpp::~OutputNekpp()
        {

        }
        
        void OutputNekpp::Process()
        {
            if (m->verbose)
            {
                cout << "OutputNekpp: Writing file..." << endl;
            }

            TiXmlDocument doc;
            TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "utf-8", "");
            doc.LinkEndChild( decl );

            TiXmlElement * root = new TiXmlElement( "NEKTAR" );
            doc.LinkEndChild( root );

            // Begin <GEOMETRY> section
            TiXmlElement * geomTag = new TiXmlElement( "GEOMETRY" );
            geomTag->SetAttribute("DIM", m->expDim);
            geomTag->SetAttribute("SPACE", m->spaceDim);
            root->LinkEndChild( geomTag );

            WriteXmlNodes     (geomTag);
            WriteXmlEdges     (geomTag);
            WriteXmlFaces     (geomTag);
            WriteXmlElements  (geomTag);
            WriteXmlCurves    (geomTag);
            WriteXmlComposites(geomTag);
            WriteXmlDomain    (geomTag);
            WriteXmlExpansions(root);
            WriteXmlConditions(root);
            
            // Save the XML file.
            doc.SaveFile(config["outfile"].as<string>());
        }

        void OutputNekpp::WriteXmlNodes(TiXmlElement * pRoot)
        {
            TiXmlElement* verTag = new TiXmlElement( "VERTEX" );
            std::set<NodeSharedPtr>::iterator it;

            std::set<NodeSharedPtr> tmp(
                    m->vertexSet.begin(),
                    m->vertexSet.end());

            for (it = tmp.begin(); it != tmp.end(); ++it)
            {
                NodeSharedPtr n = *it;
                stringstream s;
                s << scientific << setprecision(8) 
                  << n->x << " " << n->y << " " << n->z;
                TiXmlElement * v = new TiXmlElement( "V" );
                v->SetAttribute("ID",n->id);
                v->LinkEndChild(new TiXmlText(s.str()));
                verTag->LinkEndChild(v);
            }
            pRoot->LinkEndChild(verTag);
        }

        void OutputNekpp::WriteXmlEdges(TiXmlElement * pRoot)
        {
            if (m->expDim >= 2)
            {
                int edgecnt = 0;
                TiXmlElement* verTag = new TiXmlElement( "EDGE" );
                std::set<EdgeSharedPtr>::iterator it;
                std::set<EdgeSharedPtr> tmp(
                        m->edgeSet.begin(),
                        m->edgeSet.end());
                
                for (it = tmp.begin(); it != tmp.end(); ++it)
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

        void OutputNekpp::WriteXmlFaces(TiXmlElement * pRoot)
        {
            if (m->expDim == 3)
            {
                TiXmlElement* verTag = new TiXmlElement( "FACE" );
                std::set<FaceSharedPtr>::iterator it;
                std::set<FaceSharedPtr> tmp(
                        m->faceSet.begin(),
                        m->faceSet.end());

                for (it = tmp.begin(); it != tmp.end(); ++it)
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

        void OutputNekpp::WriteXmlElements(TiXmlElement * pRoot)
        {
            TiXmlElement* verTag = new TiXmlElement( "ELEMENT" );
            vector<ElementSharedPtr> &elmt = m->element[m->expDim];

            for(int i = 0; i < elmt.size(); ++i)
            {
                TiXmlElement *elm_tag = new TiXmlElement(elmt[i]->GetTag());
                elm_tag->SetAttribute("ID", elmt[i]->GetId());
                elm_tag->LinkEndChild(new TiXmlText(elmt[i]->GetXmlString()));
                verTag->LinkEndChild(elm_tag);
            }
            pRoot->LinkEndChild(verTag);
        }

        void OutputNekpp::WriteXmlCurves(TiXmlElement * pRoot)
        {
            int edgecnt = 0;
            int facecnt = 0;

            bool curve = false;
            EdgeSet::iterator it;
            for (it = m->edgeSet.begin(); it != m->edgeSet.end(); ++it)
            {
                if ((*it)->edgeNodes.size() > 0) 
                {
                    curve = true;
                    break;
                }
            }
            if (!curve) return;

            TiXmlElement * curved = new TiXmlElement ("CURVED" );

            for (it = m->edgeSet.begin(); it != m->edgeSet.end(); ++it)
            {
                if ((*it)->edgeNodes.size() > 0)
                {
                    TiXmlElement * e = new TiXmlElement( "E" );
                    e->SetAttribute("ID",        edgecnt++);
                    e->SetAttribute("EDGEID",    (*it)->id);
                    e->SetAttribute("NUMPOINTS", (*it)->GetNodeCount());
                    e->SetAttribute("TYPE", 
                        LibUtilities::kPointsTypeStr[(*it)->curveType]);
                    TiXmlText * t0 = new TiXmlText((*it)->GetXmlCurveString());
                    e->LinkEndChild(t0);
                    curved->LinkEndChild(e);
                }
            }

            // 2D elements in 3-space, output face curvature information
            if (m->expDim == 2 && m->spaceDim == 3)
            {
                vector<ElementSharedPtr>::iterator it;
                for (it = m->element[m->expDim].begin(); it != m->element[m->expDim].end(); ++it)
                {
                    // Only generate face curve if there are volume nodes
                    if ((*it)->GetVolumeNodes().size() > 0)
                    {
                        TiXmlElement * e = new TiXmlElement( "F" );
                        e->SetAttribute("ID",        facecnt++);
                        e->SetAttribute("FACEID",    (*it)->GetId());
                        e->SetAttribute("NUMPOINTS", (*it)->GetNodeCount());

                        // Quad use PolyEvenlySpaced points, tri uses
                        // NodalTriEvenlySpaced points
                        if ((*it)->GetVertexCount() == 4)
                        {
                            e->SetAttribute("TYPE", "PolyEvenlySpaced");
                        }
                        else
                        {
                            e->SetAttribute("TYPE", "NodalTriEvenlySpaced");
                        }
                        TiXmlText * t0 = new TiXmlText((*it)->GetXmlCurveString());
                        e->LinkEndChild(t0);
                        curved->LinkEndChild(e);
                    }
                }
            }

            // This code is commented out until face nodes are fully supported
            // in Nektar++.
            /*
            FaceSet::iterator it2;
            for (it2 = m->faceSet.begin(); it2 != m->faceSet.end(); ++it2)
            {
                if ((*it2)->faceNodes.size() > 0)
                {
                    TiXmlElement * f = new TiXmlElement( "F" );
                    f->SetAttribute("ID",       facecnt++);
                    f->SetAttribute("FACEID",   (*it2)->id);
                    f->SetAttribute("NUMPOINTS",(*it2)->GetNodeCount());
                    f->SetAttribute("TYPE",
                        LibUtilities::kPointsTypeStr[(*it2)->curveType]);
                    TiXmlText * t0 = new TiXmlText((*it2)->GetXmlCurveString());
                    f->LinkEndChild(t0);
                    curved->LinkEndChild(f);
                }
            }
            */
            pRoot->LinkEndChild( curved );
        }

        void OutputNekpp::WriteXmlComposites(TiXmlElement * pRoot)
        {
            TiXmlElement* verTag = new TiXmlElement("COMPOSITE");
            CompositeMap::iterator it;
            ConditionMap::iterator it2;
            int j = 0;

            for (it = m->composite.begin(); it != m->composite.end(); ++it, ++j)
            {
                if (it->second->items.size() > 0) 
                {
                    TiXmlElement *comp_tag = new TiXmlElement("C"); // Composite
                    TiXmlElement *elm_tag;
                    bool doSort = true;
                    
                    // Ensure that this composite is not used for periodic BCs!
                    for (it2  = m->condition.begin(); 
                         it2 != m->condition.end(); ++it2)
                    {
                        ConditionSharedPtr c = it2->second;
                        
                        // Ignore non-periodic boundary conditions.
                        if (find(c->type.begin(), c->type.end(), ePeriodic) ==
                            c->type.end())
                        {
                            continue;
                        }

                        for (int i = 0; i < c->composite.size(); ++i)
                        {
                            if (c->composite[i] == j)
                            {
                                doSort = false;
                            }
                        }
                    }
                    
                    comp_tag->SetAttribute("ID", it->second->id);
                    comp_tag->LinkEndChild(
                        new TiXmlText(it->second->GetXmlString(doSort)));
                    verTag->LinkEndChild(comp_tag);
                }
                else
                {
                    cout << "Composite " << it->second->id << " "
                         << "contains nothing." << endl;
                }
            }

            pRoot->LinkEndChild(verTag);
        }

        void OutputNekpp::WriteXmlDomain(TiXmlElement * pRoot)
        {
            // Write the <DOMAIN> subsection.
            TiXmlElement * domain = new TiXmlElement ("DOMAIN" );
            std::string list;
            CompositeMap::iterator it;
            
            for (it = m->composite.begin(); it != m->composite.end(); ++it)
            {
                if (it->second->items[0]->GetDim() == m->expDim)
                {
                    if (list.length() > 0)
                    {
                        list += ",";
                    }
                    list += boost::lexical_cast<std::string>(it->second->id);
                }
            }
            domain->LinkEndChild( new TiXmlText(" C[" + list + "] "));
            pRoot->LinkEndChild( domain );
        }

        void OutputNekpp::WriteXmlExpansions(TiXmlElement * pRoot)
        {
            // Write a default <EXPANSIONS> section.
            TiXmlElement * expansions = new TiXmlElement ("EXPANSIONS");
            CompositeMap::iterator it;
            
            for (it = m->composite.begin(); it != m->composite.end(); ++it)
            {
                if (it->second->items[0]->GetDim() == m->expDim)
                {
                    TiXmlElement * exp = new TiXmlElement ( "E");
                    exp->SetAttribute("COMPOSITE", "C["
                        + boost::lexical_cast<std::string>(it->second->id)
                        + "]");
                    exp->SetAttribute("NUMMODES",7);
                    exp->SetAttribute("TYPE","MODIFIED");
                    
                    if (m->fields.size() == 0)
                    {
                        exp->SetAttribute("FIELDS","u");
                    }
                    else
                    {
                        string fstr;
                        for (int i = 0; i < m->fields.size(); ++i)
                        {
                            fstr += m->fields[i]+",";
                        }
                        fstr = fstr.substr(0,fstr.length()-1);
                        exp->SetAttribute("FIELDS", fstr);
                    }
                    
                    expansions->LinkEndChild(exp);
                }
            }
            pRoot->LinkEndChild(expansions);
        }
        
        void OutputNekpp::WriteXmlConditions(TiXmlElement * pRoot)
        {
            TiXmlElement *conditions = 
                new TiXmlElement("CONDITIONS");
            TiXmlElement *boundaryregions = 
                new TiXmlElement("BOUNDARYREGIONS");
            TiXmlElement *boundaryconditions = 
                new TiXmlElement("BOUNDARYCONDITIONS");
            TiXmlElement *variables = 
                new TiXmlElement("VARIABLES");
            ConditionMap::iterator it;
            
            for (it = m->condition.begin(); it != m->condition.end(); ++it)
            {
                ConditionSharedPtr c = it->second;
                string tmp;
                
                // First set up boundary regions.
                TiXmlElement *b = new TiXmlElement("B");
                b->SetAttribute("ID", boost::lexical_cast<string>(it->first));
                
                for (int i = 0; i < c->composite.size(); ++i)
                {
                    tmp += boost::lexical_cast<string>(c->composite[i]) + ",";
                }
                
                tmp = tmp.substr(0, tmp.length()-1);

                TiXmlText *t0 = new TiXmlText("C["+tmp+"]");
                b->LinkEndChild(t0);
                boundaryregions->LinkEndChild(b);
                
                TiXmlElement *region = new TiXmlElement("REGION");
                region->SetAttribute(
                    "REF", boost::lexical_cast<string>(it->first));
                
                for (int i = 0; i < c->type.size(); ++i)
                {
                    string tagId;
                    
                    switch(c->type[i])
                    {
                        case eDirichlet:    tagId = "D"; break;
                        case eNeumann:      tagId = "N"; break;
                        case ePeriodic:     tagId = "P"; break;
                        case eHOPCondition: tagId = "N"; break;
                        default:                         break;
                    }
                    
                    TiXmlElement *tag = new TiXmlElement(tagId);
                    tag->SetAttribute("VAR", c->field[i]);
                    tag->SetAttribute("VALUE", c->value[i]);
                    
                    if (c->type[i] == eHOPCondition)
                    {
                        tag->SetAttribute("USERDEFINEDTYPE", "H");
                    }
                    
                    region->LinkEndChild(tag);
                }
                
                boundaryconditions->LinkEndChild(region);
            }

            for (int i = 0; i < m->fields.size(); ++i)
            {
                TiXmlElement *v = new TiXmlElement("V");
                v->SetAttribute("ID", boost::lexical_cast<std::string>(i));
                TiXmlText *t0 = new TiXmlText(m->fields[i]);
                v->LinkEndChild(t0);
                variables->LinkEndChild(v);
            }
            
            if (m->fields.size() > 0)
            {
                conditions->LinkEndChild(variables);
            }
            
            if (m->condition.size() > 0)
            {
                conditions->LinkEndChild(boundaryregions);
                conditions->LinkEndChild(boundaryconditions);
            }
            
            pRoot->LinkEndChild(conditions);
        }   
    }
}
