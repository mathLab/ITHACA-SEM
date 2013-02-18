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
#include <tinyxml/tinyxml.h>

namespace Nektar
{
    namespace SpatialDomains
    {

        MeshGraph1D::MeshGraph1D()
        {
        }

        MeshGraph1D::MeshGraph1D(const LibUtilities::SessionReaderSharedPtr &pSession)
            : MeshGraph(pSession)
        {
            ReadGeometry(pSession->GetDocument());
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

            int nextElementNumber = -1;

            /// All elements are of the form: "<S ID = n> ... </S>", with
            /// ? being the element type.

            TiXmlElement *segment = field->FirstChildElement("S");

            while (segment)
            {
                nextElementNumber++;

                int indx;
                int err = segment->QueryIntAttribute("ID", &indx);
                ASSERTL0(err == TIXML_SUCCESS, "Unable to read element attribute ID.");
//                ASSERTL0(indx == nextElementNumber, "Element IDs must begin with zero and be sequential.");

                TiXmlNode* elementChild = segment->FirstChild();
                while(elementChild && elementChild->Type() != TiXmlNode::TEXT)
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

                    VertexComponentSharedPtr v1 = GetVertex(vertex1);
                    VertexComponentSharedPtr v2 = GetVertex(vertex2);
                    SegGeomSharedPtr seg = MemoryManager<SegGeom>::AllocateSharedPtr(indx, v1,v2);
                    seg->SetGlobalID(indx);
                    m_segGeoms[indx] = seg;
                }
                catch(...)
                {
                    NEKERROR(ErrorUtil::efatal,
                        (std::string("Unable to read element data for segment: ") + elementStr).c_str());
                }

                /// Keep looking for additional segments
                segment = segment->NextSiblingElement("S");
            }

            ASSERTL0(nextElementNumber >= 0, "At least one element must be specified.");
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

//
// $Log: MeshGraph1D.cpp,v $
// Revision 1.17  2007/12/04 03:01:18  jfrazier
// Changed to stringstream.
//
// Revision 1.16  2007/09/20 22:25:06  jfrazier
// Added expansion information to meshgraph class.
//
// Revision 1.15  2007/07/28 05:44:27  sherwin
// Fixed for new MemoryManager call
//
// Revision 1.14  2007/07/26 01:38:32  jfrazier
// Cleanup of some attribute reading code.
//
// Revision 1.13  2007/07/24 16:52:09  jfrazier
// Added domain code.
//
// Revision 1.12  2007/07/23 16:54:30  jfrazier
// Change a dynamic allocation using new to memory manager allocate.
//
// Revision 1.11  2007/07/05 04:21:10  jfrazier
// Changed id format and propagated from 1d to 2d.
//
// Revision 1.10  2007/06/10 02:27:10  jfrazier
// Another checkin with an incremental completion of the boundary conditions reader.
//
// Revision 1.9  2007/06/07 23:55:24  jfrazier
// Intermediate revisions to add parsing for boundary conditions file.
//
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
