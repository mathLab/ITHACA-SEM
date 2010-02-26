////////////////////////////////////////////////////////////////////////////////
//
//  File: History.cpp
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
//  Description: Stores history points.
//
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/HistoryPoints.h>
#include <string>

namespace Nektar
{
    namespace SpatialDomains
    {
        History::History(const MeshGraph *meshGraph) :
            m_MeshGraph(meshGraph)
        {    
        }
        
        void History::Read(std::string &infilename)
        {
            TiXmlDocument doc(infilename);
            bool loadOkay = doc.LoadFile();

            ASSERTL0(loadOkay, (std::string("Unable to load file: ") + 
                infilename).c_str());
            
            Read(doc);
        }


        void History::Read(TiXmlDocument &doc)
        {
            TiXmlHandle docHandle(&doc);

            TiXmlNode* node = NULL;
            TiXmlElement* history = NULL;

            /// Look for all data in HISTORY block.
            history = docHandle.FirstChildElement("NEKTAR").FirstChildElement("HISTORY").Element();

            if (!history) return;

            ReadHistoryPoints(history);
        }

        void History::ReadHistoryPoints(TiXmlElement *history)
        {
            if (!history) return;

            TiXmlElement *vertex = history->FirstChildElement("H");
        
            int indx;
            int nextVertexNumber = -1;
        
            while (vertex)
            {
                nextVertexNumber++;
        
                TiXmlAttribute *vertexAttr = vertex->FirstAttribute();
                std::string attrName(vertexAttr->Name());
        
                ASSERTL0(attrName == "ID", 
                         (std::string("Unknown attribute name: ") 
                                                        + attrName).c_str());
        
                int err = vertexAttr->QueryIntValue(&indx);
                ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute ID.");
                ASSERTL0(indx == nextVertexNumber, 
                         "Vertex IDs must begin with zero and be sequential.");
        
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
        
                ASSERTL0(!vertexBodyStr.empty(), 
                         "Vertex definitions must contain vertex data.");
        
                // Get vertex data from the data string.
                double xval, yval, zval;
                std::istringstream vertexDataStrm(vertexBodyStr.c_str());
                int dim = m_MeshGraph->GetSpaceDimension();
              
                try
                {
                    while(!vertexDataStrm.fail())
                    {
                        vertexDataStrm >> xval >> yval >> zval;
        
                        // Need to check it here because we may not be good after the read
                        // indicating that there was nothing to read.
                        if (!vertexDataStrm.fail())
                        {
                            VertexComponentSharedPtr vert(
                                        MemoryManager<VertexComponent>
                                                ::AllocateSharedPtr(dim, indx, 
                                                                    xval, yval, 
                                                                    zval));
                            m_HistoryPoints.push_back(vert);
                        }
                    }
                }
                catch(...)
                {
                    ASSERTL0(false, "Unable to read HISTORY data.");
                }
        
                vertex = vertex->NextSiblingElement("H");
            }
        }

        int History::GetNumHistoryPoints() const
        {
            return m_HistoryPoints.size();
        }
        
        VertexComponentSharedPtr History::GetHistoryPoint(int idx) 
                                                                        const
        {
            return m_HistoryPoints.at(idx);
        }


    }
}
