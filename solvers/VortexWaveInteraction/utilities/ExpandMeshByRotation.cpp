///////////////////////////////////////////////////////////////////////////////
//
// File ExpandMeshByRotation.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Take a mesh and expand it by rotating it 180 degrees around x-width/2, 0
//
///////////////////////////////////////////////////////////////////////////////


#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SpatialDomains/MeshGraph.h>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <map>
#include <iomanip>
#include <tinyxml/tinyxml.h>

using namespace std;
using namespace Nektar;

void ExpandVertices(TiXmlElement* mesh, map<int,int> jointVerts, map<int,int> &fullVerts);

void ExpandEdges(TiXmlElement* mesh, map<int,int> &newVerts, 
                 map<int,int> jointEdges, map<int,int> &newEdges);

void ExpandElmts(TiXmlElement* mesh, map<int,int> &newEdges, int &nOrigElmts);

void ExpandComposites(TiXmlElement* mesh, map<int,int> fullEdges, int nOrigElmts);

int main(int argc, char *argv[])
{
    int i;
    
    if(argc !=3)
    {
        fprintf(stderr,"Usage: ./ExpandMeshByRotation meshfile  outfile\n");
        exit(1);
    }
    
    //------------------------------------------------------------
    // Create Session file which also reads meshfile. 
    LibUtilities::SessionReaderSharedPtr vSession
        = LibUtilities::SessionReader::CreateInstance(argc-1, argv);
    //-----------------------------------------------------------
    
    //-------------------------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-2]);
    SpatialDomains::MeshGraphSharedPtr graphShPt=SpatialDomains::MeshGraph::Read(vSession);
    SpatialDomains::Composite composite;
    composite = graphShPt->GetComposite(300);
    std::map<int,int> jointEdges, jointVerts, newVerts, newEdges; 
    int compsize = composite->size();
    for(i = 0; i < compsize; ++i)
    {
        SpatialDomains::Geometry1DSharedPtr tmp1 = boost::dynamic_pointer_cast<SpatialDomains::Geometry1D>((*composite)[i]);
        SpatialDomains::Geometry1DSharedPtr tmp2 = boost::dynamic_pointer_cast<SpatialDomains::Geometry1D>((*composite)[compsize-1-i]);
        jointEdges[tmp1->GetEid() ] = tmp2->GetEid();
        jointVerts[tmp1->GetVid(0)] = tmp2->GetVid(1);
        jointVerts[tmp1->GetVid(1)] = tmp2->GetVid(0);
    }


   //------------------------------------------------------------
    TiXmlDocument& meshdoc = vSession->GetDocument();

    TiXmlHandle docHandle(&meshdoc);
    TiXmlElement* mesh = docHandle.FirstChildElement("NEKTAR").FirstChildElement("GEOMETRY").Element();
    
    int nOrigElmts;
    //------------------------------------------------------------
    // Expand Mesh
    ExpandVertices(mesh,jointVerts,newVerts);

    ExpandEdges(mesh,newVerts,jointEdges,newEdges);
    
    ExpandElmts(mesh, newEdges,nOrigElmts);
    
    ExpandComposites( mesh, newEdges, nOrigElmts);

    meshdoc.SaveFile(argv[argc-1]);
    
}

void ExpandVertices(TiXmlElement* mesh, map<int,int> jointVerts, map<int,int> &newVerts)
{
    
    TiXmlElement* element = mesh->FirstChildElement("VERTEX");
    ASSERTL0(element, "Unable to find mesh VERTEX tag in file.");
    
    TiXmlElement *vertex = element->FirstChildElement("V");
    
    int indx;
    int nextVertexNumber = -1;
    int err;    /// Error value returned by TinyXML.
    
    vector<NekDouble> xpts,ypts,zpts;
    NekDouble xval, yval, zval;
        
    NekDouble yoffset = 0.0;
    while (vertex)
    {
        nextVertexNumber++;
        
        TiXmlAttribute *vertexAttr = vertex->FirstAttribute();
        std::string attrName(vertexAttr->Name());

        ASSERTL0(attrName == "ID", (std::string("Unknown attribute name: ") + attrName).c_str());
        
        err = vertexAttr->QueryIntValue(&indx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute ID.");
        
        // Now read body of vertex
        std::string vertexBodyStr;
        
        TiXmlNode *vertexBody = vertex->FirstChild();

        while (vertexBody)
        {
            // Accumulate all non-comment body data.
            if (vertexBody->Type() == TiXmlNode::TEXT)
            {
                vertexBodyStr += vertexBody->ToText()->Value();
                vertexBodyStr += " ";
            }
            
            vertexBody = vertexBody->NextSibling();
        }

        ASSERTL0(!vertexBodyStr.empty(), "Vertex definitions must contain vertex data.");

        // Get vertex data from the data string.
        std::istringstream vertexDataStrm(vertexBodyStr.c_str());

        try
        {
            while(!vertexDataStrm.fail())
            {
                vertexDataStrm >> xval >> yval >> zval;                
            }
            xpts.push_back(xval);
            ypts.push_back(yval +yoffset);
            zpts.push_back(zval);
        }
        catch(...)
        {
            ASSERTL0(false, "Unable to read VERTEX data.");
        }   
        vertex = vertex->NextSiblingElement("V");
    }

    // Add in newvertices
    int npts = xpts.size();
    NekDouble xmax = Vmath::Vmax(npts,&xpts[0],1);
    int cnt = npts;
    int cnt1 = npts;

    for(int i = 0; i < npts; ++i)
    {
        if(jointVerts.count(i) == 0)
        {
            stringstream s;
            xval = xmax - xpts[i];
            yval = -ypts[i];
            s << scientific << setprecision(8) 
              << xval << " " << yval << " " << zpts[i];
            TiXmlElement * v = new TiXmlElement( "V" );
            v->SetAttribute("ID",cnt1);
            v->LinkEndChild(new TiXmlText(s.str()));
            element->LinkEndChild(v);
            newVerts[cnt++] = cnt1++;
        }
        else
        {
            newVerts[cnt++] = jointVerts[i];
        }
    }       
}

void ExpandEdges(TiXmlElement* mesh, map<int,int> &newVerts, map<int,int> jointEdges, map<int,int> &newEdges)
{
    /// Look for elements in ELEMENT block.
    TiXmlElement* field = mesh->FirstChildElement("EDGE");
    
    ASSERTL0(field, "Unable to find EDGE tag in file.");
    
    /// All elements are of the form: "<E ID="#"> ... </E>", with
    /// ? being the element type.
    /// Read the ID field first.
    TiXmlElement *edge = field->FirstChildElement("E");
    
    /// Since all edge data is one big text block, we need to
    /// accumulate all TEXT data and then parse it.  This
    /// approach effectively skips all comments or other node
    /// types since we only care about the edge list.  We
    /// cannot handle missing edge numbers as we could with
    /// missing element numbers due to the text block format.
    std::string edgeStr;
    int indx;
    int nextEdgeNumber = -1;

    map<int,int> edgeVert0,edgeVert1;
    while(edge)
    {
        nextEdgeNumber++;
        
        int err = edge->QueryIntAttribute("ID",&indx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read edge attribute ID.");
        
        TiXmlNode *child = edge->FirstChild();
        edgeStr.clear();
        if (child->Type() == TiXmlNode::TEXT)
        {
            edgeStr += child->ToText()->ValueStr();
        }
        
        /// Now parse out the edges, three fields at a time.
        int vertex1, vertex2;
        std::istringstream edgeDataStrm(edgeStr.c_str());
        
        try
        {
            while (!edgeDataStrm.fail())
            {
                edgeDataStrm >> vertex1 >> vertex2;                
            }
            edgeVert0[indx] = vertex1;
            edgeVert1[indx] = vertex2;
        }
        catch(...)
        {
            NEKERROR(ErrorUtil::efatal, (std::string("Unable to read edge data: ") + edgeStr).c_str());
        }
        edge = edge->NextSiblingElement("E");
    }

    int nedges = edgeVert0.size();
    int cnt    = nedges;
    int cnt1   = nedges;
    int norigverts = newVerts.size();

    for(int i = 0; i < nedges; ++i)
    {
        if(jointEdges.count(i) == 0)
        {
            stringstream s;
            
            s << setw(5) << newVerts[edgeVert0[i]+norigverts] << "  "
              << newVerts[edgeVert1[i]+norigverts] << "   ";
            TiXmlElement * e = new TiXmlElement( "E" );
            e->SetAttribute("ID",cnt1);
            e->LinkEndChild( new TiXmlText(s.str()) );
            field->LinkEndChild(e);
            newEdges[cnt++] = cnt1++;
        }
        else
        {
            newEdges[cnt++] = jointEdges[i];
        }   
    }
}

void ExpandElmts(TiXmlElement* mesh, map<int,int> &newEdges, int &nelmts)
{

    TiXmlElement* field = mesh->FirstChildElement("ELEMENT");

    ASSERTL0(field, "Unable to find ELEMENT tag in file.");
    
    int nextElementNumber = -1;
    
    /// All elements are of the form: "<? ID="#"> ... </?>", with
    /// ? being the element type.
    
    TiXmlElement *element = field->FirstChildElement();

    map<int,vector<int> > ElmtEdges;
    
    while (element)
    {
        std::string elementType(element->ValueStr());
        
        ASSERTL0(elementType == "Q" || elementType == "T",
                 (std::string("Unknown 2D element type: ") + elementType).c_str());
        
        /// These should be ordered.
        nextElementNumber++;
        
        /// Read id attribute.
        int indx;
        int err = element->QueryIntAttribute("ID", &indx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read element attribute ID.");
        
        /// Read text element description.
        TiXmlNode* elementChild = element->FirstChild();
        std::string elementStr;
        while(elementChild)
        {
            if (elementChild->Type() == TiXmlNode::TEXT)
            {
                elementStr += elementChild->ToText()->ValueStr();
            }
            elementChild = elementChild->NextSibling();
        }

        ASSERTL0(!elementStr.empty(), "Unable to read element description body.");
        
        /// Parse out the element components corresponding to type of element.
        if (elementType == "T")
        {
            // Read three edge numbers
            int edge1, edge2, edge3;
            std::istringstream elementDataStrm(elementStr.c_str());
            
            try
            {
                elementDataStrm >> edge1;
                elementDataStrm >> edge2;
                elementDataStrm >> edge3;
                
                ASSERTL0(!elementDataStrm.fail(), (std::string("Unable to read element data for TRIANGLE: ") + elementStr).c_str());
                
                vector<int> edges;
                edges.push_back(edge1);
                edges.push_back(edge2);
                edges.push_back(edge3);

                ElmtEdges[indx] = edges;
            }
            catch(...)
            {
                NEKERROR(ErrorUtil::efatal, (std::string("Unable to read element data for TRIANGLE: ") + elementStr).c_str());
            }
        }
        else if (elementType == "Q")
        {
            // Read four edge numbers
            int edge1, edge2, edge3, edge4;
            std::istringstream elementDataStrm(elementStr.c_str());
            
            try
            {
                elementDataStrm >> edge1;
                elementDataStrm >> edge2;
                elementDataStrm >> edge3;
                elementDataStrm >> edge4;
                
                ASSERTL0(!elementDataStrm.fail(), (std::string("Unable to read element data for QUAD: ") + elementStr).c_str());

                vector<int> edges;
                edges.push_back(edge1);
                edges.push_back(edge2);
                edges.push_back(edge3);
                edges.push_back(edge4);

                ElmtEdges[indx] = edges;
                
            }
            catch(...)
            {
                NEKERROR(ErrorUtil::efatal,(std::string("Unable to read element data for QUAD: ") + elementStr).c_str());
            }
        }
        
        /// Keep looking
        element = element->NextSiblingElement();
    }
    
    nelmts = ElmtEdges.size();
    int nedges = newEdges.size();
    
    for(int i = 0; i < ElmtEdges.size(); ++i)
    {
        stringstream s;
        
        for(int j = 0; j < ElmtEdges[i].size() ; ++j)
        {
            s << setw(10) << newEdges[ElmtEdges[i][j]+nedges];
        }
        
        TiXmlElement * f;
        switch(ElmtEdges[i].size())
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
        f->SetAttribute("ID", i+nelmts);
        f->LinkEndChild( new TiXmlText(s.str()));
        field->LinkEndChild(f);
    }
}
    
// Generate a Nektar++ string describing the vector of ints 
string GetXmlString(char tag, vector<unsigned int> &ids)
{
    stringstream st;
    vector<unsigned int>::iterator it;
    bool range = false;
    int vId = ids[0];
    int prevId = vId;
    
    st << " " << tag << "[" << vId;
    
    for (it = ids.begin()+1; it != ids.end(); ++it){
        // store previous element ID and get current one
        prevId = vId;
        vId = (*it);
        
        // continue an already started range
        if (prevId > -1 && vId == prevId + 1)
        {
            range = true;
            // if this is the last element, it's the end of a range, so write
            if (*it == ids.back())
            {
                st << "-" << vId;
            }
            continue;
        }
        
        // terminate a range, if present
        if (range)
        {
            st << "-" << prevId;
            range = false;
        }
        
        // write what will be either a single entry or start of new range
        st << "," << vId;
    }
    // terminate
    st << "] ";
    return st.str();
}


void  ExpandComposites(TiXmlElement * mesh, map<int,int> newEdges, int nOrigElmts)
{
    TiXmlElement* field = mesh->FirstChildElement("COMPOSITE");
    ASSERTL0(field, "Unable to find COMPOSITE tag in file.");

    int nextCompositeNumber = -1;

    /// All elements are of the form: "<C ID = "N"> ... </C>".
    /// Read the ID field first.
    TiXmlElement *composite = field->FirstChildElement("C");

    while (composite)
    {
        nextCompositeNumber++;
        
        int indx;
        int err = composite->QueryIntAttribute("ID", &indx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute ID.");
        
        TiXmlNode* compositeChild = composite->FirstChild();
        // This is primarily to skip comments that may be present.
        // Comments appear as nodes just like elements.
        // We are specifically looking for text in the body
        // of the definition.
        while(compositeChild && compositeChild->Type() != TiXmlNode::TEXT)
        {
            compositeChild = compositeChild->NextSibling();
        }
        
        ASSERTL0(compositeChild, "Unable to read composite definition body.");
        std::string compositeStr = compositeChild->ToText()->ValueStr();
        
        /// Parse out the element components corresponding to type of element.
        std::istringstream compositeDataStrm(compositeStr.c_str());


        try
        {

            std::string compositeElementStr;
            compositeDataStrm >> compositeElementStr;
            
            std::istringstream tokenStream(compositeElementStr);
            char type;
            
            tokenStream >> type;
            
            // in what follows we are assuming there is only one block of data 
            std::string::size_type indxBeg = compositeElementStr.find_first_of('[') + 1;
            std::string::size_type indxEnd = compositeElementStr.find_last_of(']') - 1;
        
            ASSERTL0(indxBeg <= indxEnd, (std::string("Error reading index definition:") +  compositeElementStr).c_str());
            
            std::string indxStr = compositeElementStr.substr(indxBeg, indxEnd - indxBeg + 1);
            std::vector<unsigned int> seqVector;
            std::vector<unsigned int>::iterator seqIter;
            
            bool err = ParseUtils::GenerateSeqVector(indxStr.c_str(), seqVector);
            
            ASSERTL0(err, (std::string("Error reading composite elements: ") + indxStr).c_str());
            
            switch(type)
            {
            case 'E':   // Expand edges using newEdges map and known values
                {
                    int seqlen = seqVector.size();
                    int nedges = newEdges.size();
                    
                    map<unsigned int, unsigned int> seqMap;
                    
                    for(int i =0; i < seqlen; ++i) // set up a map of defined edges
                    {
                        seqMap[seqVector[i]] = 1;
                    }
                         
                    // if edge does not exist in composite add it to the list 
                    for(int i =0; i < seqlen; ++i)
                    {
                        if(seqMap.count(newEdges[seqVector[i]+nedges]) == 0)
                        {
                            seqVector.push_back(newEdges[seqVector[i]+nedges]);
                        }                        
                    }
                }
                break;
                
                case 'T':  case 'Q':  // Expand Triangles & Quads with new elements
                {
                    int seqlen = seqVector.size();
                    
                    for(int i = 0; i < seqlen; ++i)
                    {
                        seqVector.push_back(seqVector[i]+nOrigElmts);
                    }
                    
                    break;
                }
            default:
                NEKERROR(ErrorUtil::efatal, (std::string("Unrecognized composite token: ") + compositeElementStr).c_str());
            }
            
            // now redefine string in composite
            
            
            composite->ReplaceChild(compositeChild, TiXmlText(GetXmlString(type,seqVector)));
            

        }
        catch(...)
        {
            NEKERROR(ErrorUtil::efatal,
                     (std::string("Unable to read COMPOSITE data for composite: ") + compositeStr).c_str());
        }
        
        /// Keep looking
        composite = composite->NextSiblingElement("C");
    }
}


