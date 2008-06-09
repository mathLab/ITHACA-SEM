////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/MeshGraph3D.cpp,v $
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

#include <SpatialDomains/MeshGraph3D.h>
#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/TetGeom.h>
#include <SpatialDomains/PrismGeom.h>
#include <SpatialDomains/PyrGeom.h>
#include <SpatialDomains/HexGeom.h>
#include <SpatialDomains/ParseUtils.hpp>

namespace Nektar
{
    namespace SpatialDomains
    {
        MeshGraph3D::MeshGraph3D()
        {
        }

        MeshGraph3D::~MeshGraph3D()
        {
        }

        void MeshGraph3D::ReadGeometry(std::string &infilename)
        {
            TiXmlDocument doc(infilename);
            bool loadOkay = doc.LoadFile();

            ASSERTL0(loadOkay, (std::string("Unable to load file: ") + 
                                infilename).c_str());

            ReadGeometry(doc);
        }

        // \brief Read segments (and general MeshGraph) given TiXmlDocument.
        void MeshGraph3D::ReadGeometry(TiXmlDocument &doc)
        {
            // Read mesh first
            MeshGraph::ReadGeometry(doc);
            TiXmlHandle docHandle(&doc);

            TiXmlNode* node = NULL;
            TiXmlElement* mesh = NULL;

            /// Look for all geometry related data in GEOMETRY block.
            mesh = docHandle.FirstChildElement("NEKTAR").FirstChildElement("GEOMETRY").Element();

            ASSERTL0(mesh, "Unable to find GEOMETRY tag in file.");

            ReadEdges(doc);
            ReadFaces(doc);
            ReadElements(doc);
            ReadComposites(doc);
            ReadDomain(doc);
        }

        void MeshGraph3D::ReadEdges(TiXmlDocument &doc)
        {
            /// We know we have it since we made it this far.
            TiXmlHandle docHandle(&doc);
            TiXmlElement* mesh = docHandle.FirstChildElement("NEKTAR").FirstChildElement("GEOMETRY").Element();
            TiXmlElement* field = NULL;

            /// Look for elements in ELEMENT block.
            field = mesh->FirstChildElement("EDGE");

            ASSERTL0(field, "Unable to find EDGE tag in file.");

            /// All elements are of the form: "<E ID="#"> ... </E>", with
            /// ? being the element type.
            /// Read the ID field first.
            TiXmlElement *edge = field->FirstChildElement("E");

            /// Since all edge data is one big text block, we need to accumulate
            /// all TEXT data and then parse it.  This approach effectively skips
            /// all comments or other node types since we only care about the
            /// edge list.  We cannot handle missing edge numbers as we could
            /// with missing element numbers due to the text block format.
            std::string edgeStr;
            int indx;
            int err = 0;
            int nextEdgeNumber = -1;

            while(edge)
            {
                nextEdgeNumber++;

                int err = edge->QueryIntAttribute("ID",&indx);
                ASSERTL0(err == TIXML_SUCCESS, "Unable to read edge attribute ID.");
                ASSERTL0(indx == nextEdgeNumber, "Edge IDs must begin with zero and be sequential.");

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

                        // Must check after the read because we may be at the end and not know it.
                        // If we are at the end we will add a duplicate of the last entry if we
                        // don't check here.
                        if (!edgeDataStrm.fail())
                        {
                            VertexComponentSharedPtr vertices[2] = {GetVertex(vertex1), GetVertex(vertex2)};
                            SegGeomSharedPtr edge(MemoryManager<SegGeom>::AllocateSharedPtr(indx, m_MeshDimension, vertices));

                            m_seggeoms.push_back(edge);
                        }
                    }
                }
                catch(...)
                {
                    NEKERROR(ErrorUtil::efatal, (std::string("Unable to read edge data: ") + edgeStr).c_str());
                }

                edge = edge->NextSiblingElement("E");
            }
        }

        void MeshGraph3D::ReadFaces(TiXmlDocument &doc)
        {
        }

        void MeshGraph3D::ReadElements(TiXmlDocument &doc)
        {
            /// We know we have it since we made it this far.
            TiXmlHandle docHandle(&doc);
            TiXmlElement* mesh = docHandle.FirstChildElement("NEKTAR").FirstChildElement("GEOMETRY").Element();
            TiXmlElement* field = NULL;

            /// Look for elements in ELEMENT block.
            field = mesh->FirstChildElement("ELEMENT");

            ASSERTL0(field, "Unable to find ELEMENT tag in file.");

            int nextElementNumber = -1;

            /// All elements are of the form: "<? ID="#"> ... </?>", with
            /// ? being the element type.

            TiXmlElement *element = field->FirstChildElement();

            while (element)
            {
                std::string elementType(element->ValueStr());

                //A - tet, P - pyramid, R - prism, H - hex
                ASSERTL0(elementType == "A" || elementType == "P" || elementType == "R" || elementType == "H",
                    (std::string("Unknown 3D element type: ") + elementType).c_str());

                /// These should be ordered.
                nextElementNumber++;

                /// Read id attribute.
                int indx;
                int err = element->QueryIntAttribute("ID", &indx);
                ASSERTL0(err == TIXML_SUCCESS, "Unable to read element attribute ID.");
                ASSERTL0(indx == nextElementNumber, "Element IDs must begin with zero and be sequential.");

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

                // Tetrahedral
                if (elementType == "A")
                {
                    // Read four face numbers
                    int face1, face2, face3, face4;
                    std::istringstream elementDataStrm(elementStr.c_str());

                    try
                    {
                        elementDataStrm >> face1;
                        elementDataStrm >> face2;
                        elementDataStrm >> face3;
                        elementDataStrm >> face4;

                        ASSERTL0(!elementDataStrm.fail(), (std::string("Unable to read element data for TETRAHEDRON: ") + elementStr).c_str());

                        /// Create a TriGeom to hold the new definition.
                        TriGeomSharedPtr faces[TetGeom::kNtfaces] = 
                        {
                            GetTriGeom(face1),
                            GetTriGeom(face2),
                            GetTriGeom(face3),
                            GetTriGeom(face4)
                        };

                        StdRegions::FaceOrientation faceorient[TetGeom::kNtfaces] = 
                        {
                            //TriGeom::GetFaceOrientation(*faces[0], *faces[1]),
                            //TriGeom::GetFaceOrientation(*faces[1], *faces[2]), 
                            //TriGeom::GetFaceOrientation(*faces[2], *faces[3])
                            //TriGeom::GetFaceOrientation(*faces[3], *faces[0])
                        };

                        TetGeomSharedPtr tetgeom(MemoryManager<TetGeom>::AllocateSharedPtr(faces, faceorient));
						tetgeom->SetGlobalID(indx);

                        m_tetgeoms.push_back(tetgeom);
                    }
                    catch(...)
                    {
                        NEKERROR(ErrorUtil::efatal,
                            (std::string("Unable to read element data for TETRAHEDRON: ") + elementStr).c_str());
                    }
                }
                // Pyramid
                else if (elementType == "P")
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

                        ASSERTL0(!elementDataStrm.fail(), (std::string("Unable to read element data for PRISM: ") + elementStr).c_str());

                        /// Create a QuadGeom to hold the new definition.
                        SegGeomSharedPtr edges[QuadGeom::kNedges] = 
                        {GetSegGeom(edge1),GetSegGeom(edge2),
                         GetSegGeom(edge3),GetSegGeom(edge4)};

                        StdRegions::EdgeOrientation edgeorient[QuadGeom::kNedges] =
                        {
                            SegGeom::GetEdgeOrientation(*edges[0], *edges[1]),
                            SegGeom::GetEdgeOrientation(*edges[1], *edges[2]),
                            SegGeom::GetEdgeOrientation(*edges[2], *edges[3]),
                            SegGeom::GetEdgeOrientation(*edges[3], *edges[0])
                        };

                        QuadGeomSharedPtr quadgeom(new QuadGeom(edges, edgeorient));
						quadgeom->SetGlobalID(indx);

                        m_quadgeoms.push_back(quadgeom);

                    }
                    catch(...)
                    {
                        NEKERROR(ErrorUtil::efatal,(std::string("Unable to read element data for PRISM: ") + elementStr).c_str());
                    }
                }
                // Prism
                else if (elementType == "R")
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

                        /// Create a QuadGeom to hold the new definition.
                        SegGeomSharedPtr edges[QuadGeom::kNedges] = 
                        {GetSegGeom(edge1),GetSegGeom(edge2),
                         GetSegGeom(edge3),GetSegGeom(edge4)};

                        StdRegions::EdgeOrientation edgeorient[QuadGeom::kNedges] =
                        {
                            SegGeom::GetEdgeOrientation(*edges[0], *edges[1]),
                            SegGeom::GetEdgeOrientation(*edges[1], *edges[2]),
                            SegGeom::GetEdgeOrientation(*edges[2], *edges[3]),
                            SegGeom::GetEdgeOrientation(*edges[3], *edges[0])
                        };

                        QuadGeomSharedPtr quadgeom(new QuadGeom(edges, edgeorient));
						quadgeom->SetGlobalID(indx);

                        m_quadgeoms.push_back(quadgeom);

                    }
                    catch(...)
                    {
                        NEKERROR(ErrorUtil::efatal,(std::string("Unable to read element data for QUAD: ") + elementStr).c_str());
                    }
                }
                // Hexahedral
                else if (elementType == "H")
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

                        /// Create a QuadGeom to hold the new definition.
                        SegGeomSharedPtr edges[QuadGeom::kNedges] = 
                        {GetSegGeom(edge1),GetSegGeom(edge2),
                         GetSegGeom(edge3),GetSegGeom(edge4)};

                        StdRegions::EdgeOrientation edgeorient[QuadGeom::kNedges] =
                        {
                            SegGeom::GetEdgeOrientation(*edges[0], *edges[1]),
                            SegGeom::GetEdgeOrientation(*edges[1], *edges[2]),
                            SegGeom::GetEdgeOrientation(*edges[2], *edges[3]),
                            SegGeom::GetEdgeOrientation(*edges[3], *edges[0])
                        };

                        QuadGeomSharedPtr quadgeom(new QuadGeom(edges, edgeorient));
						quadgeom->SetGlobalID(indx);

                        m_quadgeoms.push_back(quadgeom);

                    }
                    catch(...)
                    {
                        NEKERROR(ErrorUtil::efatal,(std::string("Unable to read element data for QUAD: ") + elementStr).c_str());
                    }
                }

                /// Keep looking
                element = element->NextSiblingElement();
            }
        }

        void MeshGraph3D::ReadComposites(TiXmlDocument &doc)
        {
            TiXmlHandle docHandle(&doc);

            /// We know we have it since we made it this far.
            TiXmlElement* mesh = docHandle.FirstChildElement("NEKTAR").FirstChildElement("GEOMETRY").Element();
            TiXmlElement* field = NULL;

            ASSERTL0(mesh, "Unable to find GEOMETRY tag in file.");

            /// Look for elements in ELEMENT block.
            field = mesh->FirstChildElement("COMPOSITE");

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
                ASSERTL0(indx == nextCompositeNumber, "Composite IDs must begin with zero and be sequential.");

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

                                Composite curVector(MemoryManager<GeometryVector>::AllocateSharedPtr());
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

                /// Keep looking
                composite = composite->NextSiblingElement("C");
            }
        }


        SegGeomSharedPtr MeshGraph3D::GetSegGeom(int eID)
        {
            SegGeomSharedPtr returnval;

            if (eID >= 0 && eID < int(m_seggeoms.size()))
            {
                returnval = m_seggeoms[eID];
            }

            return returnval;
        };

        TriGeomSharedPtr MeshGraph3D::GetTriGeom(int tID)
        {
            TriGeomSharedPtr returnval;

            if (tID >= 0 && tID < int(m_trigeoms.size()))
            {
                returnval = m_trigeoms[tID];
            }

            return returnval;
        };

        // Take the string that is the composite reference and find the
        // pointer to the Geometry object corresponding to it.

        // The only allowable combinations of previous and current items
        // are V (0D); E (1D); and T and Q (2D); A (Tet, 3D), P (Pyramid, 3D), R (Prism, 3D), H (Hex, 3D).
        // Only elements of the same dimension are allowed to be grouped.
        void MeshGraph3D::ResolveGeomRef(const std::string &prevToken, const std::string &token)
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

                std::vector<unsigned int> seqVector;
                std::vector<unsigned int>::iterator seqIter;

                bool err = ParseUtils::GenerateSeqVector(indxStr.c_str(), seqVector);

                ASSERTL0(err, (std::string("Error reading composite elements: ") + indxStr).c_str());

                prevTokenStream >> prevType;

                // All composites must be of the same dimension.  This map makes things clean to compare.
                map<char, int> typeMap;
                typeMap['V'] = 1; // Vertex
                typeMap['E'] = 1; // Edge
                typeMap['T'] = 2; // Triangle
                typeMap['Q'] = 2; // Quad
                typeMap['A'] = 3; // Tet
                typeMap['P'] = 3; // Pyramid
                typeMap['R'] = 3; // Prism
                typeMap['H'] = 3; // Hex

                // Make sure only geoms of the same dimension are combined.
                bool validSequence = (prevToken.empty() || (typeMap[type] == typeMap[prevType]));

                ASSERTL0(validSequence, (std::string("Invalid combination of composite items: ")
                    + type + " and " + prevType + ".").c_str()); 

                switch(type)
                {
                case 'V':   // Vertex
                    for (seqIter = seqVector.begin(); seqIter != seqVector.end(); ++seqIter)
                    {
                        if (*seqIter >= m_vertset.size())
                        {
                            char errStr[16] = "";
                            ::sprintf(errStr, "%d", *seqIter);
                            NEKERROR(ErrorUtil::ewarning, (std::string("Unknown vertex index: ") + errStr).c_str());
                        }
                        else
                        {
                            m_MeshCompositeVector.back()->push_back(m_vertset[*seqIter]);
                        }
                    }
                    break;

                case 'E':   // Edge
                    for (seqIter = seqVector.begin(); seqIter != seqVector.end(); ++seqIter)
                    {
                        if (*seqIter >= m_seggeoms.size())
                        {
                            char errStr[16] = "";
                            ::sprintf(errStr, "%d", *seqIter);
                            NEKERROR(ErrorUtil::ewarning, (std::string("Unknown edge index: ") + errStr).c_str());
                        }
                        else
                        {
                            m_MeshCompositeVector.back()->push_back(m_seggeoms[*seqIter]);
                        }
                    }
                    break;

                case 'T':   // Triangle
                    for (seqIter = seqVector.begin(); seqIter != seqVector.end(); ++seqIter)
                    {
                        if (*seqIter >= m_trigeoms.size())
                        {
                            char errStr[16] = "";
                            ::sprintf(errStr, "%d", *seqIter);
                            NEKERROR(ErrorUtil::ewarning, (std::string("Unknown triangle index: ") + errStr).c_str());
                        }
                        else
                        {
                            m_MeshCompositeVector.back()->push_back(m_trigeoms[*seqIter]);
                        }
                    }
                    break;

                case 'Q':   // Quad
                    for (seqIter = seqVector.begin(); seqIter != seqVector.end(); ++seqIter)
                    {
                        if (*seqIter >= m_quadgeoms.size())
                        {
                            char errStr[16] = "";
                            ::sprintf(errStr, "%d", *seqIter);
                            NEKERROR(ErrorUtil::ewarning, (std::string("Unknown quad index: ") + errStr).c_str());
                        }
                        else
                        {
                            m_MeshCompositeVector.back()->push_back(m_quadgeoms[*seqIter]);
                        }
                    }
                    break;

                // Tetrahedron
                case 'A':
                    for (seqIter = seqVector.begin(); seqIter != seqVector.end(); ++seqIter)
                    {
                        if (*seqIter >= m_tetgeoms.size())
                        {
                            char errStr[16] = "";
                            ::sprintf(errStr, "%d", *seqIter);
                            NEKERROR(ErrorUtil::ewarning, (std::string("Unknown tet index: ") + errStr).c_str());
                        }
                        else
                        {
                            m_MeshCompositeVector.back()->push_back(m_tetgeoms[*seqIter]);
                        }
                    }
                    break;

                // Pyramid
                case 'P':
                    for (seqIter = seqVector.begin(); seqIter != seqVector.end(); ++seqIter)
                    {
                        if (*seqIter >= m_pyrgeoms.size())
                        {
                            char errStr[16] = "";
                            ::sprintf(errStr, "%d", *seqIter);
                            NEKERROR(ErrorUtil::ewarning, (std::string("Unknown pyramid index: ") + errStr).c_str());
                        }
                        else
                        {
                            m_MeshCompositeVector.back()->push_back(m_pyrgeoms[*seqIter]);
                        }
                    }
                    break;

                // Prism
                case 'R':
                    for (seqIter = seqVector.begin(); seqIter != seqVector.end(); ++seqIter)
                    {
                        if (*seqIter >= m_prismgeoms.size())
                        {
                            char errStr[16] = "";
                            ::sprintf(errStr, "%d", *seqIter);
                            NEKERROR(ErrorUtil::ewarning, (std::string("Unknown prism index: ") + errStr).c_str());
                        }
                        else
                        {
                            m_MeshCompositeVector.back()->push_back(m_prismgeoms[*seqIter]);
                        }
                    }
                    break;

                // Hex
                case 'H':
                    for (seqIter = seqVector.begin(); seqIter != seqVector.end(); ++seqIter)
                    {
                        if (*seqIter >= m_hexgeoms.size())
                        {
                            char errStr[16] = "";
                            ::sprintf(errStr, "%d", *seqIter);
                            NEKERROR(ErrorUtil::ewarning, (std::string("Unknown hex index: ") + errStr).c_str());
                        }
                        else
                        {
                            m_MeshCompositeVector.back()->push_back(m_hexgeoms[*seqIter]);
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

        ElementEdgeVectorSharedPtr MeshGraph3D::GetElementsFromEdge(SegGeomSharedPtr edge)
        {
            // Search Tets, Pyramids, Prisms, and Hexes
            // Need to iterate through vectors because there may be multiple
            // occurrences.
            //ElementEdgeSharedPtr elementEdge;
            //TriGeomVector::iterator triIter;

            ElementEdgeVectorSharedPtr returnval = MemoryManager<ElementEdgeVector>::AllocateSharedPtr();

            //for(triIter = m_trigeoms.begin(); triIter != m_trigeoms.end(); ++triIter)
            //{
            //    int edgeNum;
            //    if ((edgeNum = (*triIter)->WhichEdge(edge)) > -1)
            //    {
            //        elementEdge = MemoryManager<ElementEdge>::AllocateSharedPtr();
            //        elementEdge->m_Element = *triIter;
            //        elementEdge->m_EdgeIndx = edgeNum;
            //        returnval->push_back(elementEdge);
            //    }
            //}

            //QuadGeomVector::iterator quadIter;

            //for(quadIter = m_quadgeoms.begin(); quadIter != m_quadgeoms.end(); ++quadIter)
            //{
            //    int edgeNum;
            //    if ((edgeNum = (*quadIter)->WhichEdge(edge)) > -1)
            //    {
            //        elementEdge = MemoryManager<ElementEdge>::AllocateSharedPtr();
            //        elementEdge->m_Element = *quadIter;
            //        elementEdge->m_EdgeIndx = edgeNum;
            //        returnval->push_back(elementEdge);
            //    }
            //}

            return returnval;
        }

    }; //end of namespace
}; //end of namespace

//
// $Log: MeshGraph3D.cpp,v $
// Revision 1.4  2008/05/29 19:07:39  delisi
// Removed the Write(..) methods, so it is only in the base MeshGraph class. Also, added a line to set the global ID of the geometry object for every element read in.
//
// Revision 1.3  2008/02/08 23:05:52  jfrazier
// More work on 3D components.
//
// Revision 1.2  2008/02/03 05:05:16  jfrazier
// Initial checkin of 3D components.
//
//
