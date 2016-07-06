////////////////////////////////////////////////////////////////////////////////
//
//  File:  MeshGraph1D.cpp
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
//  Description:
//
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/MeshGraph1D.h>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <tinyxml.h>

using namespace std;

namespace Nektar
{
    namespace SpatialDomains
    {

        MeshGraph1D::MeshGraph1D()
        {
        }

        MeshGraph1D::MeshGraph1D(const LibUtilities::SessionReaderSharedPtr &pSession,
                                 const DomainRangeShPtr &rng)
            : MeshGraph(pSession,rng)
        {
            ReadGeometry  (pSession->GetDocument());
            ReadExpansions(pSession->GetDocument());
        }

        MeshGraph1D::~MeshGraph1D()
        {
        }

        // \brief Read segments (and general MeshGraph) given filename.
        void MeshGraph1D::ReadGeometry(const std::string &infilename)
        {
            TiXmlDocument doc(infilename);

            bool loadOkay = doc.LoadFile();
            std::stringstream errstr;
            errstr << "Unable to load file: " << infilename << "\n";
            errstr << doc.ErrorDesc() << " (Line " << doc.ErrorRow()
                   << ", Column " << doc.ErrorCol() << ")";
            ASSERTL0(loadOkay, errstr.str());

            ReadGeometry(doc);
        }

        // \brief Read segments (and general MeshGraph) given TiXmlDocument.
        void MeshGraph1D::ReadGeometry(TiXmlDocument &doc)
        {
             // Read mesh first
            MeshGraph::ReadGeometry(doc);
            TiXmlHandle docHandle(&doc);

            TiXmlElement* mesh = NULL;

            /// Look for all geometry related data in GEOMETRY block.
            mesh = docHandle.FirstChildElement("NEKTAR").FirstChildElement("GEOMETRY").Element();

            ASSERTL0(mesh, "Unable to find GEOMETRY tag in file.");

            ReadCurves(doc);
            ReadElements(doc);
            ReadComposites(doc);
            ReadDomain(doc);
        }

        void MeshGraph1D::ReadElements(TiXmlDocument &doc)
        {
            /// We know we have it since we made it this far.
            TiXmlHandle docHandle(&doc);
            TiXmlElement* mesh = docHandle.FirstChildElement("NEKTAR").FirstChildElement("GEOMETRY").Element();
            TiXmlElement* field = NULL;

            /// Look for elements in ELEMENT block.
            field = mesh->FirstChildElement("ELEMENT");

            ASSERTL0(field, "Unable to find ELEMENT tag in file.");

            /// All elements are of the form: "<S ID = n> ... </S>", with
            /// ? being the element type.

            TiXmlElement *segment = field->FirstChildElement("S");
            CurveMap::iterator it;

            while (segment)
            {
                string IsCompressed;
                segment->QueryStringAttribute("COMPRESSED",&IsCompressed); 
                
                if(IsCompressed.size()) 
                {
                    ASSERTL0(boost::iequals(IsCompressed,
                               LibUtilities::CompressData::GetCompressString()),
                             "Compressed formats do not match. Expected :"
                             + LibUtilities::CompressData::GetCompressString()
                             + " but got " + std::string(IsCompressed));

                    // Extract the face body
                    TiXmlNode* child = segment->FirstChild();
                    ASSERTL0(child, "Unable to extract the data from "
                             "the compressed face tag.");

                    std::string str;
                    if (child->Type() == TiXmlNode::TINYXML_TEXT)
                    {
                        str += child->ToText()->ValueStr();
                    }

                    int indx;

                    std::vector<LibUtilities::MeshEdge> data;
                    LibUtilities::CompressData::ZlibDecodeFromBase64Str(str,
                                                                        data);

                    for(int i = 0; i < data.size(); ++i)
                    {
                        indx = data[i].id;

                        /// See if this face has curves.
                        it = m_curvedEdges.find(indx);

                        PointGeomSharedPtr vertices[2] = {
                                GetVertex(data[i].v0), GetVertex(data[i].v1)};
                        SegGeomSharedPtr seg;

                        if (it == m_curvedEdges.end())
                        {
                            seg = MemoryManager<SegGeom>::AllocateSharedPtr(
                                            indx, m_spaceDimension, vertices);
                            seg->SetGlobalID(indx); // Set global mesh id
                        }
                        else
                        {
                            seg = MemoryManager<SegGeom>::AllocateSharedPtr(
                                            indx, m_spaceDimension,
                                            vertices, it->second);
                            seg->SetGlobalID(indx); //Set global mesh id
                        }
                        seg->SetGlobalID(indx);
                        m_segGeoms[indx] = seg;
                    }
                }
                else
                {
                    
                    int indx;
                    int err = segment->QueryIntAttribute("ID", &indx);
                    ASSERTL0(err == TIXML_SUCCESS, "Unable to read element attribute ID.");
                    
                    TiXmlNode* elementChild = segment->FirstChild();
                    while(elementChild && elementChild->Type() != TiXmlNode::TINYXML_TEXT)
                    {
                        elementChild = elementChild->NextSibling();
                    }
                    
                    ASSERTL0(elementChild, "Unable to read element description body.");
                    std::string elementStr = elementChild->ToText()->ValueStr();
                    
                    /// Parse out the element components corresponding to type of element.
                    /// Read two vertex numbers
                    int vertex1, vertex2;
                    std::istringstream elementDataStrm(elementStr.c_str());
                    
                    try
                    {
                        elementDataStrm >> vertex1;
                        elementDataStrm >> vertex2;
                        
                        ASSERTL0(!elementDataStrm.fail(), (std::string("Unable to read element data for SEGMENT: ") + elementStr).c_str());
                        
                        PointGeomSharedPtr vertices[2] = {GetVertex(vertex1), GetVertex(vertex2)};
                        SegGeomSharedPtr seg;
                        it = m_curvedEdges.find(indx);
                        
                        if (it == m_curvedEdges.end())
                        {
                            seg = MemoryManager<SegGeom>::AllocateSharedPtr(indx, m_spaceDimension, vertices);
                            seg->SetGlobalID(indx); // Set global mesh id
                        }
                        else
                        {
                            seg = MemoryManager<SegGeom>::AllocateSharedPtr(indx, m_spaceDimension, vertices, it->second);
                            seg->SetGlobalID(indx); //Set global mesh id
                        }
                        seg->SetGlobalID(indx);
                        m_segGeoms[indx] = seg;
                    }
                    catch(...)
                    {
                        NEKERROR(ErrorUtil::efatal,
                                 (std::string("Unable to read element data for segment: ") + elementStr).c_str());
                    }
                }
                /// Keep looking for additional segments
                segment = segment->NextSiblingElement("S");
            }
        }

        void MeshGraph1D::ReadComposites(TiXmlDocument &doc)
        {
            TiXmlHandle docHandle(&doc);

            /// We know we have it since we made it this far.
            TiXmlElement* mesh = docHandle.FirstChildElement("NEKTAR").FirstChildElement("GEOMETRY").Element();
            TiXmlElement* field = NULL;

            ASSERTL0(mesh, "Unable to find GEOMETRY tag in file.");

            /// Look for elements in ELEMENT block.
            field = mesh->FirstChildElement("COMPOSITE");

            ASSERTL0(field, "Unable to find COMPOSITE tag in file.");

            TiXmlElement *node = field->FirstChildElement("C");

            // Sequential counter for the composite numbers.
            int nextCompositeNumber = -1;

            while (node)
            {
                /// All elements are of the form: "<? ID="#"> ... </?>", with
                /// ? being the element type.

                nextCompositeNumber++;

                int indx;
                int err = node->QueryIntAttribute("ID", &indx);
                ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute ID.");
                //ASSERTL0(indx == nextCompositeNumber, "Composite IDs must begin with zero and be sequential.");

                TiXmlNode* compositeChild = node->FirstChild();
                // This is primarily to skip comments that may be present.
                // Comments appear as nodes just like elements.
                // We are specifically looking for text in the body
                // of the definition.
                while(compositeChild && compositeChild->Type() != TiXmlNode::TINYXML_TEXT)
                {
                    compositeChild = compositeChild->NextSibling();
                }

                ASSERTL0(compositeChild, "Unable to read composite definition body.");
                std::string compositeStr = compositeChild->ToText()->ValueStr();

                /// Parse out the element components corresponding to type of element.

                std::istringstream compositeDataStrm(compositeStr.c_str());

                try
                {
                    bool first = true;
                    std::string prevCompositeElementStr;

                    while (!compositeDataStrm.fail())
                    {
                        std::string compositeElementStr;
                        compositeDataStrm >> compositeElementStr;

                        if (!compositeDataStrm.fail())
                        {
                            if (first)
                            {
                                first = false;

                                Composite curVector = MemoryManager<std::vector<GeometrySharedPtr> >::AllocateSharedPtr();
                                m_meshComposites[indx] = curVector;
                            }

                            if (compositeElementStr.length() > 0)
                            {
                                ResolveGeomRef(prevCompositeElementStr, compositeElementStr, m_meshComposites[indx]);
                            }
                            prevCompositeElementStr = compositeElementStr;
                        }
                    }
                }
                catch(...)
                {
                    NEKERROR(ErrorUtil::efatal,
                        (std::string("Unable to read COMPOSITE data for composite: ") + compositeStr).c_str());
                }

                /// Keep looking for additional composite definitions.
                node = node->NextSiblingElement("C");
            }

            ASSERTL0(nextCompositeNumber >= 0, "At least one composite must be specified.");
        }


        // Take the string that is the composite reference and find the
        // pointer to the Geometry object corresponding to it.

        // Only allow segments to be grouped for 1D mesh.
        void MeshGraph1D::ResolveGeomRef(const std::string &prevToken, const std::string &token,
                Composite& composite)
        {
            try
            {
                std::istringstream tokenStream(token);
                std::istringstream prevTokenStream(prevToken);

                char type;
                char prevType;

                tokenStream >> type;

                std::string::size_type indxBeg = token.find_first_of('[') + 1;
                std::string::size_type indxEnd = token.find_last_of(']') - 1;

                ASSERTL0(indxBeg <= indxEnd, (std::string("Error reading index definition:") + token).c_str());

                std::string indxStr = token.substr(indxBeg, indxEnd - indxBeg + 1);

                typedef vector<unsigned int> SeqVectorType;
                SeqVectorType seqVector;

                if (!ParseUtils::GenerateSeqVector(indxStr.c_str(), seqVector))
                {
                    NEKERROR(ErrorUtil::efatal, (std::string("Ill-formed sequence definition: ") + indxStr).c_str());
                }

                prevTokenStream >> prevType;

                // All composites must be of the same dimension.
                bool validSequence = (prevToken.empty() ||         // No previous, then current is just fine.
                    (type == 'V' && prevType == 'V') ||
                    (type == 'S' && prevType == 'S'));

                ASSERTL0(validSequence, (std::string("Invalid combination of composite items: ")
                    + type + " and " + prevType + ".").c_str());

                switch(type)
                {
                case 'V':   // Vertex
                    for (SeqVectorType::iterator iter=seqVector.begin(); iter!=seqVector.end(); ++iter)
                    {
                        if (m_vertSet.find(*iter) == m_vertSet.end())
                        {
                            char errStr[16] = "";
                            ::sprintf(errStr, "%d", *iter);
                            NEKERROR(ErrorUtil::ewarning, (std::string("Unknown vertex index: ") + errStr).c_str());
                        }
                        else
                        {
                            composite->push_back(m_vertSet[*iter]);
                        }
                    }
                    break;

                case 'S':   // Segment
                    for (SeqVectorType::iterator iter=seqVector.begin(); iter!=seqVector.end(); ++iter)
                    {
                        if (m_segGeoms.find(*iter) == m_segGeoms.end())
                        {
                            char errStr[16] = "";
                            ::sprintf(errStr, "%d", *iter);
                            NEKERROR(ErrorUtil::ewarning, (std::string("Unknown segment index: ") + errStr).c_str());
                        }
                        else
                        {
                            composite->push_back(m_segGeoms[*iter]);
                        }
                    }
                    break;

                default:
                    NEKERROR(ErrorUtil::efatal, (std::string("Unrecognized composite token: ") + token).c_str());
                }
            }
            catch(...)
            {
                NEKERROR(ErrorUtil::efatal, (std::string("Problem processing composite token: ") + token).c_str());
            }

            return;
        }

    }; //end of namespace
}; //end of namespace
