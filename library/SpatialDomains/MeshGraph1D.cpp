////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/MeshGraph1D.cpp,v $
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
#include "pchSpatialDomains.h"

#include <SpatialDomains/MeshGraph1D.h>
#include <SpatialDomains/ParseUtils.hpp>

namespace Nektar
{
    namespace SpatialDomains
    {

        MeshGraph1D::MeshGraph1D()
        {
        }

        MeshGraph1D::~MeshGraph1D()
        {
        }

        // \brief Read segments (and general MeshGraph) given filename.
        void MeshGraph1D::Read(std::string &infilename)
        {
            SetFileName(infilename);
            TiXmlDocument doc(infilename);

            bool loadOkay = doc.LoadFile();

            ASSERTL0(loadOkay, (std::string("Unable to load file:") + infilename + ".").c_str());

            Read(doc);
        }

        // \brief Read segments (and general MeshGraph) given TiXmlDocument.
        void MeshGraph1D::Read(TiXmlDocument &doc)
        {
             // Read mesh first
            MeshGraph::Read(doc);
            TiXmlHandle docHandle(&doc);

            TiXmlNode* node = NULL;
            TiXmlElement* mesh = NULL;

            /// Look for all geometry related data in GEOMETRY block.
            mesh = docHandle.FirstChildElement("NEKTAR").FirstChildElement("GEOMETRY").Element();

            ASSERTL0(mesh, "Unable to find GEOMETRY tag in file.");

            ReadElements(doc);
            ReadComposites(doc);
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

            TiXmlAttribute *attr = field->FirstAttribute();
            int elementNumber = 0;
            int nextElementNumber = 0;
            int err = 0;

            /// All elements are of the form: "ID <?> ... </?>", with
            /// ? being the element type.

            /// Read the ID field first.
            TiXmlNode *child = field->FirstChild();

            /// Get first TEXT node, just in case any comments are present.
            /// If no TEXT appears before an ELEMENT, then we will assume
            /// a running number.
            int elementID = nextElementNumber++;

            while (child)
            {
                int type = child->Type();

                if (type == TiXmlNode::TEXT)
                {
                    /// Read the text...it should be our index number
                    nextElementNumber = ::atoi(child->ToText()->Value());
                }
                else if (type == TiXmlNode::ELEMENT)
                {
                    std::string elementType(child->ValueStr());

                    ASSERTL0(elementType == "S",
                        (std::string("Unknown 2D element type: ") + elementType).c_str());

                    elementID = nextElementNumber++;
                    /// Read text element description.

                    TiXmlNode* elementChild = child->FirstChild();
                    while(elementChild->Type() != TiXmlNode::TEXT)
                    {
                        elementChild = elementChild->NextSibling();
                    }

                    ASSERTL0(elementChild, "Unable to read element description body.");
                    std::string elementStr = elementChild->ToText()->ValueStr();

                    /// Parse out the element components corresponding to type of element.
                    if (elementType == "S")
                    {
                        // Read two vertex numbers
                        int vertex1, vertex2;
                        std::istrstream elementDataStrm(elementStr.c_str());

                        try
                        {
                            elementDataStrm >> vertex1;
                            elementDataStrm >> vertex2;

                            ASSERTL0(!elementDataStrm.fail(), (std::string("Unable to read element data for SEGMENT: ") + elementStr).c_str());

                            VertexComponentSharedPtr vertShPtr1 = GetVertex(vertex1);
                            VertexComponentSharedPtr vertShPtr2 = GetVertex(vertex2);

                            SegGeomSharedPtr seg = MemoryManager<SegGeom>::AllocateSharedPtr(elementID, 
                                GetVertex(vertex1), GetVertex(vertex2));
                            m_seggeoms.push_back(seg);                        
                        }
                        catch(...)
                        {
                            NEKERROR(ErrorUtil::efatal,
                                (std::string("Unable to read element data for TRIANGLE: ") + elementStr).c_str());
                        }
                    }
                }

                /// Keep looking
                child = child->NextSibling();
            }
        }

        //void MeshGraph1D::ReadSegments(TiXmlDocument &doc)
        //{
        //    // Read mesh first
        //    MeshGraph::Read(doc);

        //    TiXmlHandle docHandle(&doc);

        //    TiXmlNode* node = NULL;
        //    TiXmlElement *element = NULL;
        //    int err;    /// Error value returned by TinyXML.

        //    /// Look for segments in ELEMENT block.
        //    element = docHandle.FirstChildElement("NEKTAR").FirstChildElement("GEOMETRY").FirstChildElement("ELEMENT").Element();

        //    ASSERTL0(element, "Unable to find ELEMENT tag in file.");
        //    TiXmlAttribute *attr = element->FirstAttribute();
        //    int numSegments = -1;

        //    while (attr)
        //    {
        //        std::string attrName(attr->Name());
        //        if (attrName == "NUMBER")
        //        {
        //            err = attr->QueryIntValue(&numSegments);
        //            if (err)
        //            {
        //                numSegments = -1;
        //            }
        //        }
        //        // Get the next attribute.  Shouldn't be any more, but just in case.
        //        attr = attr->Next();
        //    }

        //    ASSERTL0(numSegments > 0, "Unable to read NUMBER attribute value.");
        //    int indx = -1;

        //    node = element->FirstChild();
        //    for (int i=0; i<numSegments; ++i)
        //    {
        //        /// First read index number
        //        /// Then read tagged element

        //        int nodeType = node->Type();
        //        if (nodeType == TiXmlNode::TEXT)
        //        {
        //            // Read index number
        //            indx = atoi(node->ToText()->Value());
        //            --i; // Don't count the index
        //        }
        //        else if (nodeType == TiXmlNode::COMMENT)
        //        {
        //            --numSegments; /// Leave it pointing where it is and skip it.
        //        }
        //        else if (nodeType == TiXmlNode::ELEMENT)
        //        {
        //            /// All segments are embedded in <S></S> tags.
        //            TiXmlElement* segment = element->FirstChildElement();
        //            if (std::string(segment->Value())!= "S")
        //            {
        //                NEKERROR(ErrorUtil::ewarning,(std::string("Element type not allowed: ") +
        //                    segment->Value()).c_str());
        //            }
        //            else
        //            {
        //                //// Get the entire data block then go through it one piece at a time.
        //                std::string segmentData;

        //                // Look through all immediate children and
        //                // accumulate all text (type TEXT).  This is just
        //                // in case it is broken up by other tagged
        //                // elements.  Using element->GetText() will not
        //                // work if text contains XML comments, for
        //                // example.
        //                TiXmlNode *child = node->FirstChild();
        //                while(child)
        //                {
        //                    if (child->Type() == TiXmlNode::TEXT)
        //                    {
        //                        segmentData += child->ToText()->Value();
        //                        segmentData += " "; // Don't jam together.
        //                    }
        //                    child = child->NextSibling();
        //                }

        //                // Get segment data from the data string.
        //                int vertex1, vertex2;
        //                std::istrstream vertexDataStrm(segmentData.c_str());

        //                try
        //                {
        //                    vertexDataStrm >> vertex1 >> vertex2;

        //                    ASSERTL0(!vertexDataStrm.fail(), 
        //                        (std::string("Unable to read segment data: ")
        //                        + std::string(segmentData) + ".").c_str());

        //                    SegGeomSharedPtr seg(new SegGeom(indx, 
        //                        GetVertex(vertex1), GetVertex(vertex2)));
        //                    m_seggeoms.push_back(seg);
        //                }
        //                catch(...)
        //                {
        //                    NEKERROR(ErrorUtil::efatal, (std::string("Unable to read segment data: "
        //                        + segmentData + ".").c_str()));
        //                }

        //                indx = -1;
        //            }
        //        }
        //        else
        //        {
        //            ASSERTL1(false, "Unknown node type.");
        //        }

        //        node = node->NextSibling();
        //    }
        //}

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

            TiXmlAttribute *attr = field->FirstAttribute();
            int nextCompositeNumber = 0;

            /// All elements are of the form: "ID <C ID="#"> ... </?>", with
            /// ? being the element type.

            /// Read the ID field first.
            TiXmlNode *child = field->FirstChild();

            while (child)
            {
                int type = child->Type();

                if (type == TiXmlNode::TEXT)
                {
                    /// Read the text...it should be our index number
                    int indx = ::atoi(child->ToText()->Value());
#pragma message("Fix this")
                    ASSERTL0(indx == nextCompositeNumber, "");
                    ++nextCompositeNumber;
                }
                else if (type == TiXmlNode::ELEMENT)// XML element, not mesh element
                {
                    std::string elementType(child->ValueStr());

                    ASSERTL0(elementType == "C",
                        (std::string("Unknown XML tag in COMPOSITE definition: ") + elementType).c_str());

                    /// Read text element description.

                    TiXmlNode* compositeChild = child->FirstChild();
                    // This is primarily to skip comments that may be present.
                    // Comments appear as nodes just like elements.
                    // We are specifically looking for text in the body
                    // of the definition.
                    while(compositeChild->Type() != TiXmlNode::TEXT)
                    {
                        compositeChild = compositeChild->NextSibling();
                    }

                    ASSERTL0(compositeChild, "Unable to read composite description body.");
                    std::string compositeStr = compositeChild->ToText()->ValueStr();

                    /// Parse out the element components corresponding to type of element.

                    std::istrstream compositeDataStrm(compositeStr.c_str());

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

                                    Composite curVector(new std::vector<GeometrySharedPtr>);
                                    m_MeshCompositeVector.push_back(curVector);
                                }

                                if (compositeElementStr.length() > 0)
                                {
                                    ResolveGeomRef(prevCompositeElementStr, compositeElementStr);
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
                }

                /// Keep looking
                child = child->NextSibling();
            }
        }
        // Take the string that is the composite reference and find the
        // pointer to the Geometry object corresponding to it.

        // Only allow segments to be grouped for 1D mesh.
        void MeshGraph1D::ResolveGeomRef(const std::string &prevToken, const std::string &token)
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
                        if (*iter >= m_vertset.size())
                        {
                            char errStr[16] = "";
                            ::sprintf(errStr, "%d", *iter);
                            NEKERROR(ErrorUtil::ewarning, (std::string("Unknown vertex index: ") + errStr).c_str());
                        }
                        else
                        {
                            m_MeshCompositeVector.back()->push_back(m_vertset[*iter]);
                        }
                    }
                    break;

                case 'S':   // Segment
                    for (SeqVectorType::iterator iter=seqVector.begin(); iter!=seqVector.end(); ++iter)
                    {
                        if (*iter >= m_seggeoms.size())
                        {
                            char errStr[16] = "";
                            ::sprintf(errStr, "%d", *iter);
                            NEKERROR(ErrorUtil::ewarning, (std::string("Unknown segment index: ") + errStr).c_str());
                        }
                        else
                        {
                            m_MeshCompositeVector.back()->push_back(m_seggeoms[*iter]);
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

        void MeshGraph1D::Write(std::string &outfilename)
        {
        }

    }; //end of namespace
}; //end of namespace

//
// $Log: MeshGraph1D.cpp,v $
// Revision 1.8  2007/03/14 21:24:08  sherwin
// Update for working version of MultiRegions up to ExpList1D
//
// Revision 1.7  2006/10/15 06:18:58  sherwin
// Moved NekPoint out of namespace LibUtilities
//
// Revision 1.6  2006/09/26 23:41:53  jfrazier
// Updated to account for highest level NEKTAR tag and changed the geometry tag to GEOMETRY.
//
// Revision 1.5  2006/06/01 14:15:30  sherwin
// Added typdef of boost wrappers and made GeoFac a boost shared pointer.
//
// Revision 1.4  2006/05/23 19:56:33  jfrazier
// These build and run, but the expansion pieces are commented out
// because they would not run.
//
// Revision 1.3  2006/05/16 22:28:31  sherwin
// Updates to add in FaceComponent call to constructors
//
// Revision 1.2  2006/05/16 20:12:59  jfrazier
// Minor fixes to correct bugs.
//
// Revision 1.1  2006/05/04 18:59:01  kirby
// *** empty log message ***
//
// Revision 1.17  2006/04/09 02:08:35  jfrazier
// Added precompiled header.
//
// Revision 1.16  2006/04/04 23:12:37  jfrazier
// More updates to readers.  Still work to do on MeshGraph2D to store tris and quads.
//
// Revision 1.15  2006/03/25 00:58:29  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.14  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.13  2006/03/12 07:42:03  sherwin
//
// Updated member names and StdRegions call. Still has not been compiled
//
// Revision 1.12  2006/02/26 21:19:43  bnelson
// Fixed a variety of compiler errors caused by updates to the coding standard.
//
// Revision 1.11  2006/02/19 01:37:33  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
