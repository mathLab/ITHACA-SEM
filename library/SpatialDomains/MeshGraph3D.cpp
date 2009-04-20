////////////////////////////////////////////////////////////////////////////////
//
//  File: MeshGraph3D.cpp
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
            /// We know we have it since we made it this far.
            TiXmlHandle docHandle(&doc);
            TiXmlElement* mesh = docHandle.FirstChildElement("NEKTAR").FirstChildElement("GEOMETRY").Element();
            TiXmlElement* field = NULL;

            /// Look for elements in FACE block.
            field = mesh->FirstChildElement("FACE");

            ASSERTL0(field, "Unable to find FACE tag in file.");

            int nextFaceNumber = -1;

            /// All faces are of the form: "<? ID="#"> ... </?>", with
			/// ? being an element type (either Q or T).

            TiXmlElement *element = field->FirstChildElement();

            while (element)
            {
                std::string elementType(element->ValueStr());

                ASSERTL0(elementType == "Q" || elementType == "T",
                    (std::string("Unknown 3D face type: ") + elementType).c_str());

                /// These should be ordered.
                nextFaceNumber++;

                /// Read id attribute.
                int indx;
                int err = element->QueryIntAttribute("ID", &indx);
                ASSERTL0(err == TIXML_SUCCESS, "Unable to read face attribute ID.");
                ASSERTL0(indx == nextFaceNumber, "Face IDs must begin with zero and be sequential.");

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

                ASSERTL0(!elementStr.empty(), "Unable to read face description body.");

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

                        ASSERTL0(!elementDataStrm.fail(), (std::string("Unable to read face data for TRIANGLE: ") + elementStr).c_str());

                        /// Create a TriGeom to hold the new definition.
                        SegGeomSharedPtr edges[TriGeom::kNedges] = 
                        {
                            GetSegGeom(edge1),
                            GetSegGeom(edge2),
                            GetSegGeom(edge3)
                        };

                        StdRegions::EdgeOrientation edgeorient[TriGeom::kNedges] = 
                        {
                            SegGeom::GetEdgeOrientation(*edges[0], *edges[1]),
                            SegGeom::GetEdgeOrientation(*edges[1], *edges[2]), 
                            SegGeom::GetEdgeOrientation(*edges[2], *edges[0])
                        };

                        TriGeomSharedPtr trigeom(MemoryManager<TriGeom>::AllocateSharedPtr(indx, edges, edgeorient));
                        trigeom->SetGlobalID(indx);

                        m_trigeoms.push_back(trigeom);
                    }
                    catch(...)
                    {
                        NEKERROR(ErrorUtil::efatal,
                            (std::string("Unable to read face data for TRIANGLE: ") + elementStr).c_str());
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

                        ASSERTL0(!elementDataStrm.fail(), (std::string("Unable to read face data for QUAD: ") + elementStr).c_str());

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

                        QuadGeomSharedPtr quadgeom(MemoryManager<QuadGeom>::AllocateSharedPtr(indx, edges, edgeorient));
                        quadgeom->SetGlobalID(indx);

                        m_quadgeoms.push_back(quadgeom);

                    }
                    catch(...)
                    {
                        NEKERROR(ErrorUtil::efatal,(std::string("Unable to read face data for QUAD: ") + elementStr).c_str());
                    }
                }

                /// Keep looking
                element = element->NextSiblingElement();
            }
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

                std::istringstream elementDataStrm(elementStr.c_str());

                /// Parse out the element components corresponding to type of element.

                // Tetrahedral
                if (elementType == "A")
                {
                    try
                    {
						/// Create arrays for the tri and quad faces.
						const int kNfaces = TetGeom::kNfaces;
						const int kNtfaces = TetGeom::kNtfaces;
						const int kNqfaces = TetGeom::kNqfaces;
						TriGeomSharedPtr tfaces[kNtfaces];
						//QuadGeomSharedPtr qfaces[kNqfaces];
						int Ntfaces = 0;
						int Nqfaces = 0;

						/// Fill the arrays and make sure there aren't too many faces.
						std::stringstream errorstring;
						errorstring << "Element " << indx << " must have " << kNtfaces << " triangle face(s), and " << kNqfaces << " quadrilateral face(s).";
						for (int i = 0; i < kNfaces; i++)
						{
							int faceID;
							elementDataStrm >> faceID;
							Geometry2DSharedPtr face = GetGeometry2D(faceID);
							if (face == Geometry2DSharedPtr() ||
								(face->GetGeomShapeType() != eTriangle && face->GetGeomShapeType() != eQuadrilateral))
							{
								std::stringstream errorstring;
								errorstring << "Element " << indx << " has invalid face: " << faceID;
								ASSERTL0(false, errorstring.str().c_str());
							}
							else if (face->GetGeomShapeType() == eTriangle)
							{
								ASSERTL0(Ntfaces < kNtfaces, errorstring.str().c_str());
								tfaces[Ntfaces++] = boost::static_pointer_cast<TriGeom>(face);
							}
							else if (face->GetGeomShapeType() == eQuadrilateral)
							{
								ASSERTL0(Nqfaces < kNqfaces, errorstring.str().c_str());
								//qfaces[Nqfaces++] = boost::static_pointer_cast<QuadGeom>(face);
							}
						}

						/// Make sure all of the face indicies could be read, and that there weren't too few.
						ASSERTL0(!elementDataStrm.fail(), (std::string("Unable to read element data for TETRAHEDRON: ") + elementStr).c_str());
						ASSERTL0(Ntfaces == kNtfaces, errorstring.str().c_str());
						ASSERTL0(Nqfaces == kNqfaces, errorstring.str().c_str());

                        StdRegions::FaceOrientation faceorient[kNfaces] = 
                        {
                            //TriGeom::GetFaceOrientation(*faces[0], *faces[1]),
                            //TriGeom::GetFaceOrientation(*faces[1], *faces[2]), 
                            //TriGeom::GetFaceOrientation(*faces[2], *faces[3])
                            //TriGeom::GetFaceOrientation(*faces[3], *faces[0])
                        };

                        //TetGeomSharedPtr tetgeom(MemoryManager<TetGeom>::AllocateSharedPtr(tfaces, qfaces, faceorient));
						TetGeomSharedPtr tetgeom(MemoryManager<TetGeom>::AllocateSharedPtr(tfaces, faceorient));
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
                    try
                    {
						/// Create arrays for the tri and quad faces.
						const int kNfaces = PyrGeom::kNfaces;
						const int kNtfaces = PyrGeom::kNtfaces;
						const int kNqfaces = PyrGeom::kNqfaces;
						TriGeomSharedPtr tfaces[kNtfaces];
						QuadGeomSharedPtr qfaces[kNqfaces];
						int Ntfaces = 0;
						int Nqfaces = 0;

						/// Fill the arrays and make sure there aren't too many faces.
						std::stringstream errorstring;
						errorstring << "Element " << indx << " must have " << kNtfaces << " triangle face(s), and " << kNqfaces << " quadrilateral face(s).";
						for (int i = 0; i < kNfaces; i++)
						{
							int faceID;
							elementDataStrm >> faceID;
							Geometry2DSharedPtr face = GetGeometry2D(faceID);
							if (face == Geometry2DSharedPtr() ||
								(face->GetGeomShapeType() != eTriangle && face->GetGeomShapeType() != eQuadrilateral))
							{
								std::stringstream errorstring;
								errorstring << "Element " << indx << " has invalid face: " << faceID;
								ASSERTL0(false, errorstring.str().c_str());
							}
							else if (face->GetGeomShapeType() == eTriangle)
							{
								ASSERTL0(Ntfaces < kNtfaces, errorstring.str().c_str());
								tfaces[Ntfaces++] = boost::static_pointer_cast<TriGeom>(face);
							}
							else if (face->GetGeomShapeType() == eQuadrilateral)
							{
								ASSERTL0(Nqfaces < kNqfaces, errorstring.str().c_str());
								qfaces[Nqfaces++] = boost::static_pointer_cast<QuadGeom>(face);
							}
						}

						/// Make sure all of the face indicies could be read, and that there weren't too few.
						ASSERTL0(!elementDataStrm.fail(), (std::string("Unable to read element data for PYRAMID: ") + elementStr).c_str());
						ASSERTL0(Ntfaces == kNtfaces, errorstring.str().c_str());
						ASSERTL0(Nqfaces == kNqfaces, errorstring.str().c_str());

                        StdRegions::FaceOrientation faceorient[kNfaces] = 
                        {
                            //TriGeom::GetFaceOrientation(*faces[0], *faces[1]),
                            //TriGeom::GetFaceOrientation(*faces[1], *faces[2]), 
                            //TriGeom::GetFaceOrientation(*faces[2], *faces[3])
                            //TriGeom::GetFaceOrientation(*faces[3], *faces[0])
                        };

                        PyrGeomSharedPtr pyrgeom(MemoryManager<PyrGeom>::AllocateSharedPtr(tfaces, qfaces, faceorient));
						pyrgeom->SetGlobalID(indx);

                        m_pyrgeoms.push_back(pyrgeom);
                    }
                    catch(...)
                    {
                        NEKERROR(ErrorUtil::efatal,
                            (std::string("Unable to read element data for PYRAMID: ") + elementStr).c_str());
                    }
                }
                // Prism
                else if (elementType == "R")
                {
                    try
                    {
						/// Create arrays for the tri and quad faces.
						const int kNfaces = PrismGeom::kNfaces;
						const int kNtfaces = PrismGeom::kNtfaces;
						const int kNqfaces = PrismGeom::kNqfaces;
						TriGeomSharedPtr tfaces[kNtfaces];
						QuadGeomSharedPtr qfaces[kNqfaces];
						int Ntfaces = 0;
						int Nqfaces = 0;

						/// Fill the arrays and make sure there aren't too many faces.
						std::stringstream errorstring;
						errorstring << "Element " << indx << " must have " << kNtfaces << " triangle face(s), and " << kNqfaces << " quadrilateral face(s).";
						for (int i = 0; i < kNfaces; i++)
						{
							int faceID;
							elementDataStrm >> faceID;
							Geometry2DSharedPtr face = GetGeometry2D(faceID);
							if (face == Geometry2DSharedPtr() ||
								(face->GetGeomShapeType() != eTriangle && face->GetGeomShapeType() != eQuadrilateral))
							{
								std::stringstream errorstring;
								errorstring << "Element " << indx << " has invalid face: " << faceID;
								ASSERTL0(false, errorstring.str().c_str());
							}
							else if (face->GetGeomShapeType() == eTriangle)
							{
								ASSERTL0(Ntfaces < kNtfaces, errorstring.str().c_str());
								tfaces[Ntfaces++] = boost::static_pointer_cast<TriGeom>(face);
							}
							else if (face->GetGeomShapeType() == eQuadrilateral)
							{
								ASSERTL0(Nqfaces < kNqfaces, errorstring.str().c_str());
								qfaces[Nqfaces++] = boost::static_pointer_cast<QuadGeom>(face);
							}
						}

						/// Make sure all of the face indicies could be read, and that there weren't too few.
						ASSERTL0(!elementDataStrm.fail(), (std::string("Unable to read element data for PRISM: ") + elementStr).c_str());
						ASSERTL0(Ntfaces == kNtfaces, errorstring.str().c_str());
						ASSERTL0(Nqfaces == kNqfaces, errorstring.str().c_str());

                        StdRegions::FaceOrientation faceorient[kNfaces] = 
                        {
                            //TriGeom::GetFaceOrientation(*faces[0], *faces[1]),
                            //TriGeom::GetFaceOrientation(*faces[1], *faces[2]), 
                            //TriGeom::GetFaceOrientation(*faces[2], *faces[3])
                            //TriGeom::GetFaceOrientation(*faces[3], *faces[0])
                        };

                        PrismGeomSharedPtr prismgeom(MemoryManager<PrismGeom>::AllocateSharedPtr(tfaces, qfaces, faceorient));
						prismgeom->SetGlobalID(indx);

                        m_prismgeoms.push_back(prismgeom);
                    }
                    catch(...)
                    {
                        NEKERROR(ErrorUtil::efatal,
                            (std::string("Unable to read element data for PRISM: ") + elementStr).c_str());
                    }
                }
                // Hexahedral
                else if (elementType == "H")
                {
                    try
                    {
						/// Create arrays for the tri and quad faces.
						const int kNfaces = HexGeom::kNfaces;
						const int kNtfaces = HexGeom::kNtfaces;
						const int kNqfaces = HexGeom::kNqfaces;
						//TriGeomSharedPtr tfaces[kNtfaces];
						QuadGeomSharedPtr qfaces[kNqfaces];
						int Ntfaces = 0;
						int Nqfaces = 0;

						/// Fill the arrays and make sure there aren't too many faces.
						std::stringstream errorstring;
						errorstring << "Element " << indx << " must have " << kNtfaces << " triangle face(s), and " << kNqfaces << " quadrilateral face(s).";
						for (int i = 0; i < kNfaces; i++)
						{
							int faceID;
							elementDataStrm >> faceID;
							Geometry2DSharedPtr face = GetGeometry2D(faceID);
							if (face == Geometry2DSharedPtr() ||
								(face->GetGeomShapeType() != eTriangle && face->GetGeomShapeType() != eQuadrilateral))
							{
								std::stringstream errorstring;
								errorstring << "Element " << indx << " has invalid face: " << faceID;
								ASSERTL0(false, errorstring.str().c_str());
							}
							else if (face->GetGeomShapeType() == eTriangle)
							{
								ASSERTL0(Ntfaces < kNtfaces, errorstring.str().c_str());
								//tfaces[Ntfaces++] = boost::static_pointer_cast<TriGeom>(face);
							}
							else if (face->GetGeomShapeType() == eQuadrilateral)
							{
								ASSERTL0(Nqfaces < kNqfaces, errorstring.str().c_str());
								qfaces[Nqfaces++] = boost::static_pointer_cast<QuadGeom>(face);
							}
						}

						/// Make sure all of the face indicies could be read, and that there weren't too few.
						ASSERTL0(!elementDataStrm.fail(), (std::string("Unable to read element data for HEXAHEDRAL: ") + elementStr).c_str());
						ASSERTL0(Ntfaces == kNtfaces, errorstring.str().c_str());
						ASSERTL0(Nqfaces == kNqfaces, errorstring.str().c_str());

                        StdRegions::FaceOrientation faceorient[kNfaces] = 
                        {
                            //TriGeom::GetFaceOrientation(*faces[0], *faces[1]),
                            //TriGeom::GetFaceOrientation(*faces[1], *faces[2]), 
                            //TriGeom::GetFaceOrientation(*faces[2], *faces[3])
                            //TriGeom::GetFaceOrientation(*faces[3], *faces[0])
                        };

                        //HexGeomSharedPtr hexgeom(MemoryManager<HexGeom>::AllocateSharedPtr(tfaces, qfaces, faceorient));
                        HexGeomSharedPtr hexgeom(MemoryManager<HexGeom>::AllocateSharedPtr(qfaces));
                        hexgeom->SetGlobalID(indx);

                        m_hexgeoms.push_back(hexgeom);
                    }
                    catch(...)
                    {
                        NEKERROR(ErrorUtil::efatal,
                            (std::string("Unable to read element data for HEXAHEDRAL: ") + elementStr).c_str());
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

        Geometry2DSharedPtr MeshGraph3D::GetGeometry2D(int gID)
        {
			for (TriGeomVectorIter iter = m_trigeoms.begin(); iter != m_trigeoms.end(); iter++)
			{
				if ((*iter)->GetGlobalID() == gID)
				{
					return *iter;
				}
			}

			for (QuadGeomVectorIter iter = m_quadgeoms.begin(); iter != m_quadgeoms.end(); iter++)
			{
				if ((*iter)->GetGlobalID() == gID)
				{
					return *iter;
				}
			}

            return Geometry2DSharedPtr();
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

                case 'F':   // Face
                    for (seqIter = seqVector.begin(); seqIter != seqVector.end(); ++seqIter)
                    {
						Geometry2DSharedPtr face = GetGeometry2D(*seqIter);
                        if (face == Geometry2DSharedPtr())
                        {
                            char errStr[16] = "";
                            ::sprintf(errStr, "%d", *seqIter);
                            NEKERROR(ErrorUtil::ewarning, (std::string("Unknown face index: ") + errStr).c_str());
                        }
                        else
                        {
                            m_MeshCompositeVector.back()->push_back(face);
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


        ElementFaceVectorSharedPtr MeshGraph3D::GetElementsFromFace(Geometry2DSharedPtr face)
        {
            // Search tets, prisms, pyramids, and hexahedrons
            // Need to iterate through vectors because there may be multiple
            // occurrences.
            ElementFaceSharedPtr elementFace;

            ElementFaceVectorSharedPtr returnval = MemoryManager<ElementFaceVector>::AllocateSharedPtr();

            CompositeVectorIter compIter;

            TetGeomSharedPtr tetGeomShPtr;
            HexGeomSharedPtr hexGeomShPtr;
            PrismGeomSharedPtr prismGeomShPtr;
            PyrGeomSharedPtr pyrGeomShPtr;

            GeometryVectorIter geomIter;

            for (compIter = m_Domain.begin(); compIter != m_Domain.end(); ++compIter)
            {
                for (geomIter = (*compIter)->begin(); geomIter != (*compIter)->end(); ++geomIter)
                {
                    tetGeomShPtr = boost::dynamic_pointer_cast<TetGeom>(*geomIter);
                    hexGeomShPtr = boost::dynamic_pointer_cast<HexGeom>(*geomIter);
                    prismGeomShPtr = boost::dynamic_pointer_cast<PrismGeom>(*geomIter);
                    pyrGeomShPtr = boost::dynamic_pointer_cast<PyrGeom>(*geomIter);

                    if (tetGeomShPtr || hexGeomShPtr || prismGeomShPtr || pyrGeomShPtr)
                    {
                        int faceNum;                        
                        if (tetGeomShPtr)
                        {
                            if ((faceNum = tetGeomShPtr->WhichFace(face)) > -1)
                            {
                                elementFace = MemoryManager<ElementFace>::AllocateSharedPtr();
                                elementFace->m_Face = tetGeomShPtr;
                                elementFace->m_FaceIndx = faceNum;
                                returnval->push_back(elementFace);
                            }
                        }
                        else if (hexGeomShPtr)
                        {
                            if ((faceNum = hexGeomShPtr->WhichFace(face)) > -1)
                            {
                                elementFace = MemoryManager<ElementFace>::AllocateSharedPtr();
                                elementFace->m_Face = hexGeomShPtr;
                                elementFace->m_FaceIndx = faceNum;
                                returnval->push_back(elementFace);
                            }
                        }
                        else if (prismGeomShPtr)
                        {
                            if ((faceNum = prismGeomShPtr->WhichFace(face)) > -1)
                            {
                                elementFace = MemoryManager<ElementFace>::AllocateSharedPtr();
                                elementFace->m_Face = prismGeomShPtr;
                                elementFace->m_FaceIndx = faceNum;
                                returnval->push_back(elementFace);
                            }
                        }
                        else if (pyrGeomShPtr)
                        {
                            if ((faceNum = pyrGeomShPtr->WhichFace(face)) > -1)
                            {
                                elementFace = MemoryManager<ElementFace>::AllocateSharedPtr();
                                elementFace->m_Face = pyrGeomShPtr;
                                elementFace->m_FaceIndx = faceNum;
                                returnval->push_back(elementFace);
                            }
                        }
                    }
                }         
            }

 
            return returnval;
        }

        LibUtilities::BasisKey MeshGraph3D:: GetFaceBasisKey(Geometry2DSharedPtr face, const int flag)
        {
            ElementFaceVectorSharedPtr elements = GetElementsFromFace(face);
            // Perhaps, a check should be done here to ensure that in case 
            // elements->size!=1, all elements to which the edge belongs have the same type
            // and order of expansion such that no confusion can arise.
            ExpansionShPtr expansion = GetExpansion((*elements)[0]->m_Face);


            ASSERTL0(false,"This method needs reworking");
#if 0 
            int nummodes = (int) expansion->m_NumModesEqn.Evaluate();
            
            switch(expansion->m_ExpansionType)
            {
            case eModified:
                {
                    switch (flag)
                    {
                    case 0:
                        {
                            const LibUtilities::PointsKey pkey(nummodes+1,LibUtilities::eGaussLobattoLegendre);
                            return LibUtilities::BasisKey(LibUtilities::eModified_A,nummodes,pkey);
                        }
                        break;
                    case 1:
                        {
                            const LibUtilities::PointsKey pkey(nummodes,LibUtilities::eGaussRadauMAlpha1Beta0);
                            return LibUtilities::BasisKey(LibUtilities::eModified_B,nummodes,pkey);
                        }
                        break;
                    default:
                        ASSERTL0(false,"invalid value to flag");
                        break;
                    }
                }
                break;
            case eGLL_Lagrange:
                {
                    TriGeomSharedPtr triangle = boost::dynamic_pointer_cast<TriGeom>(face);
                    QuadGeomSharedPtr quadrilateral = boost::dynamic_pointer_cast<QuadGeom>(face);

                    if(quadrilateral)
                    {
                        const LibUtilities::PointsKey pkey(nummodes+1,LibUtilities::eGaussLobattoLegendre);
                        return LibUtilities::BasisKey(LibUtilities::eGLL_Lagrange,nummodes,pkey);
                    }
                    else if(triangle)
                    {
                        switch (flag)
                        {
                        case 0:
                            {
                                const LibUtilities::PointsKey pkey(nummodes+1,LibUtilities::eGaussLobattoLegendre);
                                return LibUtilities::BasisKey(LibUtilities::eOrtho_A,nummodes,pkey);
                            }
                            break;
                        case 1:
                            {
                                const LibUtilities::PointsKey pkey(nummodes,LibUtilities::eGaussRadauMAlpha1Beta0);
                                return LibUtilities::BasisKey(LibUtilities::eOrtho_B,nummodes,pkey);
                            }
                            break;
                        default:
                            ASSERTL0(false,"invalid value to flag");
                            break;
                        }
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a proper Geometry2D failed");
                    }  
                }
                break;
            case eOrthogonal:
                {
                    switch (flag)
                    {
                    case 0:
                        {
                            const LibUtilities::PointsKey pkey(nummodes+1,LibUtilities::eGaussLobattoLegendre);
                            return LibUtilities::BasisKey(LibUtilities::eOrtho_A,nummodes,pkey);
                        }
                        break;
                    case 1:
                        {
                            const LibUtilities::PointsKey pkey(nummodes,LibUtilities::eGaussRadauMAlpha1Beta0);
                            return LibUtilities::BasisKey(LibUtilities::eOrtho_B,nummodes,pkey);
                        }
                        break;
                    default:
                        ASSERTL0(false,"invalid value to flag");
                        break;
                        }
                }
                break;              
            case eGLL_Lagrange_SEM:
                {
                    const LibUtilities::PointsKey pkey(nummodes,LibUtilities::eGaussLobattoLegendre);
                    return LibUtilities::BasisKey(LibUtilities::eGLL_Lagrange,nummodes,pkey);
                }
                break;     
            default:
                ASSERTL0(false,"expansion type unknown");
                return LibUtilities::NullBasisKey; // Keep things happy by returning a value.
                break;
            }            
#endif
        }



    }; //end of namespace
}; //end of namespace

//
// $Log: MeshGraph3D.cpp,v $
// Revision 1.12  2009/01/12 10:26:59  pvos
// Added input tags for nodal expansions
//
// Revision 1.11  2008/09/23 18:19:56  pvos
// Updates for working ProjectContField3D demo
//
// Revision 1.10  2008/09/12 11:26:19  pvos
// Updates for mappings in 3D
//
// Revision 1.9  2008/08/26 02:24:38  ehan
// Added GetElementFromFace()
//
// Revision 1.8  2008/06/30 19:34:26  ehan
// Fixed infinity recursive-loop error.
//
// Revision 1.7  2008/06/12 19:56:05  delisi
// Changed some error handling for reading 3D geometries.
//
// Revision 1.6  2008/06/11 16:10:12  delisi
// Added the 3D reader.
//
// Revision 1.5  2008/06/09 22:36:09  delisi
// Changed Tet representation from 'E' to 'A', since 'E' is used by edges.
//
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
