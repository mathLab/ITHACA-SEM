//////////////////////////////////////////////////////////////////////////////////
//  File:  MeshGraph2D.cpp
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

#include <SpatialDomains/MeshGraph2D.h>
#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/TriGeom.h>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <tinyxml/tinyxml.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        MeshGraph2D::MeshGraph2D()
        {
        }

        MeshGraph2D::~MeshGraph2D()
        {
        }

        MeshGraph2D::MeshGraph2D(const LibUtilities::SessionReaderSharedPtr &pSession)
            : MeshGraph(pSession)
        {
            ReadGeometry(pSession->GetDocument());
            ReadExpansions(pSession->GetDocument());
        }

        void MeshGraph2D::ReadGeometry(const std::string &infilename)
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
        void MeshGraph2D::ReadGeometry(TiXmlDocument &doc)
        {
            // Read mesh first
            MeshGraph::ReadGeometry(doc);
            TiXmlHandle docHandle(&doc);

            TiXmlElement* mesh = NULL;

            /// Look for all geometry related data in GEOMETRY block.
            mesh = docHandle.FirstChildElement("NEKTAR").FirstChildElement("GEOMETRY").Element();

            ASSERTL0(mesh, "Unable to find GEOMETRY tag in file.");
            
            ReadCurves(doc);
            ReadEdges(doc);
            ReadElements(doc);
            ReadComposites(doc);
            ReadDomain(doc);
        }

        void MeshGraph2D::ReadEdges(TiXmlDocument &doc)
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

            /// Since all edge data is one big text block, we need to
            /// accumulate all TEXT data and then parse it.  This
            /// approach effectively skips all comments or other node
            /// types since we only care about the edge list.  We
            /// cannot handle missing edge numbers as we could with
            /// missing element numbers due to the text block format.
            std::string edgeStr;
            int i,indx;
            int nextEdgeNumber = -1;

            // Curved Edges
            map<int, int> edge_curved;
            for(i = 0; i < m_curvedEdges.size(); ++i)
            {
                edge_curved[m_curvedEdges[i]->m_curveID] = i;
            }

            while(edge)
            {
                nextEdgeNumber++;

                int err = edge->QueryIntAttribute("ID",&indx);
                ASSERTL0(err == TIXML_SUCCESS, "Unable to read edge attribute ID.");
//                ASSERTL0(indx == nextEdgeNumber, "Edge IDs must begin with zero and be sequential.");

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

                        // Must check after the read because we may be
                        // at the end and not know it.  If we are at
                        // the end we will add a duplicate of the last
                        // entry if we don't check here.
                        if (!edgeDataStrm.fail())
                        {
                            VertexComponentSharedPtr vertices[2] = {GetVertex(vertex1), GetVertex(vertex2)};

                            SegGeomSharedPtr edge;

                            if(edge_curved.count(indx) == 0)
                            {
                                edge = MemoryManager<SegGeom>::AllocateSharedPtr(indx, m_spaceDimension, vertices);
                                edge->SetGlobalID(indx); // Set global mesh id
                            }
                            else
                            {
                                edge = MemoryManager<SegGeom>::AllocateSharedPtr(indx, m_spaceDimension, vertices, m_curvedEdges[edge_curved.find(indx)->second]);
                                edge->SetGlobalID(indx); //Set global mesh id
                            }

                            m_segGeoms[indx] = edge;
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

        void MeshGraph2D::ReadElements(TiXmlDocument &doc)
        {
            /// We know we have it since we made it this far.
            TiXmlHandle docHandle(&doc);
            TiXmlElement* mesh = docHandle.FirstChildElement("NEKTAR").FirstChildElement("GEOMETRY").Element();
            TiXmlElement* field = NULL;

            /// Look for elements in ELEMENT block.
            field = mesh->FirstChildElement("ELEMENT");

            ASSERTL0(field, "Unable to find ELEMENT tag in file.");

            // Set up curve map for curved elements on an embedded manifold.
            map<int, int> faceCurves;
            map<int,int>::iterator x;
            for (int i = 0; i < m_curvedFaces.size(); ++i)
            {
                faceCurves[m_curvedFaces[i]->m_curveID] = i;
            }

            int nextElementNumber = -1;

            /// All elements are of the form: "<? ID="#"> ... </?>", with
            /// ? being the element type.

            TiXmlElement *element = field->FirstChildElement();

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
//                    ASSERTL0(indx == nextElementNumber, "Element IDs must begin with zero and be sequential.");

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

                            /// Create a TriGeom to hold the new definition.
                            SegGeomSharedPtr edges[TriGeom::kNedges] =
                        {
                            GetSegGeom(edge1),
                            GetSegGeom(edge2),
                            GetSegGeom(edge3)
                        };

                            StdRegions::Orientation edgeorient[TriGeom::kNedges] =
                        {
                            SegGeom::GetEdgeOrientation(*edges[0], *edges[1]),
                            SegGeom::GetEdgeOrientation(*edges[1], *edges[2]),
                            SegGeom::GetEdgeOrientation(*edges[2], *edges[0])
                        };

                            TriGeomSharedPtr trigeom;
                            if ((x = faceCurves.find(indx)) == faceCurves.end())
                            {
                                trigeom = MemoryManager<TriGeom>
                                            ::AllocateSharedPtr(indx,
                                                    edges,
                                                    edgeorient);
                            }
                            else
                            {
                                trigeom = MemoryManager<TriGeom>
                                            ::AllocateSharedPtr(indx,
                                                    edges,
                                                    edgeorient,
                                                    m_curvedFaces[x->second]);
                            }
                            trigeom->SetGlobalID(indx);

                            m_triGeoms[indx] = trigeom;
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

                            /// Create a QuadGeom to hold the new definition.
                            SegGeomSharedPtr edges[QuadGeom::kNedges] =
                        {GetSegGeom(edge1),GetSegGeom(edge2),
                         GetSegGeom(edge3),GetSegGeom(edge4)};

                            StdRegions::Orientation edgeorient[QuadGeom::kNedges] =
                        {
                            SegGeom::GetEdgeOrientation(*edges[0], *edges[1]),
                            SegGeom::GetEdgeOrientation(*edges[1], *edges[2]),
                            SegGeom::GetEdgeOrientation(*edges[2], *edges[3]),
                            SegGeom::GetEdgeOrientation(*edges[3], *edges[0])
                        };

                            QuadGeomSharedPtr quadgeom;
                            if ((x = faceCurves.find(indx)) == faceCurves.end())
                            {
                                quadgeom = MemoryManager<QuadGeom>
                                            ::AllocateSharedPtr(indx,
                                                    edges,
                                                    edgeorient);
                            }
                            else
                            {
                                quadgeom = MemoryManager<QuadGeom>
                                            ::AllocateSharedPtr(indx,
                                                    edges,
                                                    edgeorient,
                                                    m_curvedFaces[x->second]);
                            }
                            quadgeom->SetGlobalID(indx);

                            m_quadGeoms[indx] = quadgeom;

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

        void MeshGraph2D::ReadComposites(TiXmlDocument &doc)
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
//                ASSERTL0(indx == nextCompositeNumber, "Composite IDs must begin with zero and be sequential.");

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

                /// Keep looking
                composite = composite->NextSiblingElement("C");
            }
        }


        SegGeomSharedPtr MeshGraph2D::GetSegGeom(int eID)
        {
            SegGeomSharedPtr returnval;
            SegGeomMap::iterator x = m_segGeoms.find(eID);
            ASSERTL0(x != m_segGeoms.end(), "Segment not found.");
            return x->second;
        };


        // Take the string that is the composite reference and find the
        // pointer to the Geometry object corresponding to it.

        // The only allowable combinations of previous and current items
        // are V (0D); E (1D); and T and Q (2D).  Only elements of the same
        // dimension are allowed to be grouped.
        void MeshGraph2D::ResolveGeomRef(const std::string &prevToken, const std::string &token,
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
                std::vector<unsigned int> seqVector;
                std::vector<unsigned int>::iterator seqIter;

                bool err = ParseUtils::GenerateSeqVector(indxStr.c_str(), seqVector);

                ASSERTL0(err, (std::string("Error reading composite elements: ") + indxStr).c_str());

                prevTokenStream >> prevType;

                // All composites must be of the same dimension.
                bool validSequence = (prevToken.empty() ||         // No previous, then current is just fine.
                                      (type == 'V' && prevType == 'V') ||
                                      (type == 'E' && prevType == 'E') ||
                                      ((type == 'T' || type == 'Q') &&
                                       (prevType == 'T' || prevType == 'Q')));

                ASSERTL0(validSequence, (std::string("Invalid combination of composite items: ")
                                         + type + " and " + prevType + ".").c_str());


                switch(type)
                {
                case 'E':   // Edge
                    for (seqIter = seqVector.begin(); seqIter != seqVector.end(); ++seqIter)
                    {
                        if (m_segGeoms.find(*seqIter) == m_segGeoms.end())
                        {
                            char errStr[16] = "";
                            ::sprintf(errStr, "%d", *seqIter);
                            NEKERROR(ErrorUtil::ewarning, (std::string("Unknown edge index: ") + errStr).c_str());
                        }
                        else
                        {
                            composite->push_back(m_segGeoms[*seqIter]);
                        }
                    }
                    break;

                case 'T':   // Triangle
#if 1
                    {
                        // Set up inverse maps of tris which takes global id
                        // back to local storage in m_triGeoms;
//                        map<int, int> tri_id_map;
//                        for(int i = 0; i < m_triGeoms.size(); ++i)
//                        {
//                            tri_id_map[m_triGeoms[i]->GetGlobalID()] = i;
//                        }

                        for (seqIter = seqVector.begin(); seqIter != seqVector.end(); ++seqIter)
                        {
                            if (m_triGeoms.count(*seqIter) == 0 )
                            {
                                char errStr[16] = "";
                                ::sprintf(errStr, "%d", *seqIter);
                                NEKERROR(ErrorUtil::ewarning, (std::string("Unknown triangle index: ") + errStr).c_str());
                            }
                            else
                            {
                                composite->push_back(m_triGeoms[*seqIter]);
                            }
                        }
                    }
#else
                    for (seqIter = seqVector.begin(); seqIter != seqVector.end(); ++seqIter)
                    {
                        if (*seqIter >= m_triGeoms.size())
                        {
                            char errStr[16] = "";
                            ::sprintf(errStr, "%d", *seqIter);
                            NEKERROR(ErrorUtil::ewarning, (std::string("Unknown triangle index: ") + errStr+std::string(" in Composite section")).c_str());
                        }
                        else
                        {
                            composite->push_back(m_triGeoms[*seqIter]);
                        }
                    }
#endif
                    break;

                case 'Q':   // Quad
#if 1
                    {
                        // Set up inverse maps of tris which takes global id
                        // back to local storage in m_triGeoms;
//                        map<int, int> quad_id_map;
//                        for(int i = 0; i < m_quadGeoms.size(); ++i)
//                        {
//                            quad_id_map[m_quadGeoms[i]->GetGlobalID()] = i;
//                        }

                        for (seqIter = seqVector.begin(); seqIter != seqVector.end(); ++seqIter)
                        {
                            if (m_quadGeoms.count(*seqIter) == 0)
                            {
                                char errStr[16] = "";
                                ::sprintf(errStr, "%d", *seqIter);
                                NEKERROR(ErrorUtil::ewarning, (std::string("Unknown quad index: ") + errStr +std::string(" in Composite section")).c_str());
                            }
                            else
                            {
                                composite->push_back(m_quadGeoms[*seqIter]);
                            }
                        }
                    }
#else
                    for (seqIter = seqVector.begin(); seqIter != seqVector.end(); ++seqIter)
                    {
                        if (*seqIter >= m_quadGeoms.size())
                        {
                            char errStr[16] = "";
                            ::sprintf(errStr, "%d", *seqIter);
                            NEKERROR(ErrorUtil::ewarning, (std::string("Unknown quad index: ") + errStr).c_str());
                        }
                        else
                        {
                            composite->push_back(m_quadGeoms[*seqIter]);
                        }
                    }
#endif
                    break;

                case 'V':   // Vertex
                    for (seqIter = seqVector.begin(); seqIter != seqVector.end(); ++seqIter)
                    {
                        if (*seqIter >= m_vertSet.size())
                        {
                            char errStr[16] = "";
                            ::sprintf(errStr, "%d", *seqIter);
                            NEKERROR(ErrorUtil::ewarning, (std::string("Unknown vertex index: ") + errStr).c_str());
                        }
                        else
                        {
                            composite->push_back(m_vertSet[*seqIter]);
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

        ElementEdgeVectorSharedPtr MeshGraph2D::GetElementsFromEdge(Geometry1DSharedPtr edge)
        {
            SegGeomSharedPtr Sedge;

            if(!(Sedge = boost::dynamic_pointer_cast<SegGeom>(edge)))
            {
                ASSERTL0(false,"Dynamics cast failed");

            }
            return GetElementsFromEdge(Sedge);

        }
        ElementEdgeVectorSharedPtr MeshGraph2D::GetElementsFromEdge(SegGeomSharedPtr edge)
        {
            // Search tris and quads
            // Need to iterate through vectors because there may be multiple
            // occurrences.
            ElementEdgeSharedPtr elementEdge;
            //TriGeomVectorIter triIter;

            ElementEdgeVectorSharedPtr returnval = MemoryManager<ElementEdgeVector>::AllocateSharedPtr();

            CompositeMapIter compIter;
            TriGeomSharedPtr triGeomShPtr;
            QuadGeomSharedPtr quadGeomShPtr;

            GeometryVectorIter geomIter;

            for (compIter = m_domain.begin(); compIter != m_domain.end(); ++compIter)
            {
                for (geomIter = (compIter->second)->begin(); geomIter != (compIter->second)->end(); ++geomIter)
                {
                    triGeomShPtr = boost::dynamic_pointer_cast<TriGeom>(*geomIter);
                    quadGeomShPtr = boost::dynamic_pointer_cast<QuadGeom>(*geomIter);

                    if (triGeomShPtr || quadGeomShPtr)
                    {
                        int edgeNum;
                        if (triGeomShPtr)
                        {
                            if ((edgeNum = triGeomShPtr->WhichEdge(edge)) > -1)
                            {
                                elementEdge = MemoryManager<ElementEdge>::AllocateSharedPtr();
                                elementEdge->m_Element = triGeomShPtr;
                                elementEdge->m_EdgeIndx = edgeNum;
                                returnval->push_back(elementEdge);
                            }
                        }
                        else if (quadGeomShPtr)
                        {
                            if ((edgeNum = quadGeomShPtr->WhichEdge(edge)) > -1)
                            {
                                elementEdge = MemoryManager<ElementEdge>::AllocateSharedPtr();
                                elementEdge->m_Element = quadGeomShPtr;
                                elementEdge->m_EdgeIndx = edgeNum;
                                returnval->push_back(elementEdge);
                            }
                        }
                    }
                }
            }

                //for(triIter = m_triGeoms.begin(); triIter != m_triGeoms.end(); ++triIter)
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

                //for(quadIter = m_quadGeoms.begin(); quadIter != m_quadGeoms.end(); ++quadIter)
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

        LibUtilities::BasisKey MeshGraph2D::GetEdgeBasisKey(SegGeomSharedPtr edge, const std::string variable)
        {
            ElementEdgeVectorSharedPtr elements = GetElementsFromEdge(edge);
            // Perhaps, a check should be done here to ensure that
            // in case elements->size!=1, all elements to which
            // the edge belongs have the same type and order of
            // expansion such that no confusion can arise.
            ExpansionShPtr expansion = GetExpansion((*elements)[0]->m_Element, variable);

            int edge_id = (*elements)[0]->m_EdgeIndx;

            if((*elements)[0]->m_Element->GetShapeType() == LibUtilities::eTriangle)
            {
                edge_id = (edge_id)? 1:0;
            }
            else
            {
                edge_id = edge_id%2;
            }

            int nummodes  = expansion->m_basisKeyVector[edge_id].GetNumModes();
            int numpoints = expansion->m_basisKeyVector[edge_id].GetNumPoints();

            if((*elements)[0]->m_Element->GetShapeType() == LibUtilities::eTriangle)
            {
                // Use edge 0 to define basis of order relevant to edge
                switch(expansion->m_basisKeyVector[edge_id].GetBasisType())
                {
                case LibUtilities::eGLL_Lagrange:
                    {
                        switch(expansion->m_basisKeyVector[edge_id].GetPointsType())
                        {
                        case LibUtilities::eGaussLobattoLegendre:
                            {
                                const  LibUtilities::PointsKey pkey(numpoints+1,LibUtilities::eGaussLobattoLegendre);
                                return LibUtilities::BasisKey(expansion->m_basisKeyVector[0].GetBasisType(),nummodes,pkey);
                            }
                            break;

                        default:
                            ASSERTL0(false,"Unexpected points distribution");

                            // It doesn't matter what we return
                            // here since the ASSERT will stop
                            // execution.  Just return something
                            // to prevent warnings messages.
                            const  LibUtilities::PointsKey pkey(numpoints+1,LibUtilities::eGaussLobattoLegendre);
                            return LibUtilities::BasisKey(expansion->m_basisKeyVector[0].GetBasisType(),nummodes,pkey);
                            break;
                        }
                    }
                    break;
                case LibUtilities::eOrtho_B: // Assume this is called from nodal triangular basis
                    {
                        switch(expansion->m_basisKeyVector[edge_id].GetPointsType())
                        {
                        case LibUtilities::eGaussRadauMAlpha1Beta0:
                            {
                                const  LibUtilities::PointsKey pkey(numpoints+1,LibUtilities::eGaussLobattoLegendre);
                                return LibUtilities::BasisKey(LibUtilities::eGLL_Lagrange,nummodes,pkey);
                            }
                            break;

                        default:
                            ASSERTL0(false,"Unexpected points distribution");

                            // It doesn't matter what we return
                            // here since the ASSERT will stop
                            // execution.  Just return something
                            // to prevent warnings messages.
                            const LibUtilities::PointsKey pkey(numpoints+1,LibUtilities::eGaussLobattoLegendre);
                            return LibUtilities::BasisKey(expansion->m_basisKeyVector[0].GetBasisType(),nummodes,pkey);
                            break;
                        }
                    }
                    break;
                case LibUtilities::eModified_B:
                    {
                        switch(expansion->m_basisKeyVector[edge_id].GetPointsType())
                        {
                        case LibUtilities::eGaussRadauMAlpha1Beta0:
                            {
                                const LibUtilities::PointsKey pkey(numpoints+1,LibUtilities::eGaussLobattoLegendre);
                                return LibUtilities::BasisKey(expansion->m_basisKeyVector[0].GetBasisType(),nummodes,pkey);
                            }
                            break;

                        default:
                            ASSERTL0(false,"Unexpected points distribution");

                            // It doesn't matter what we return
                            // here since the ASSERT will stop
                            // execution.  Just return something
                            // to prevent warnings messages.
                            const LibUtilities::PointsKey pkey(numpoints+1,LibUtilities::eGaussLobattoLegendre);
                            return LibUtilities::BasisKey(expansion->m_basisKeyVector[0].GetBasisType(),nummodes,pkey);
                            break;
                        }
                    }
                    break;
                case LibUtilities::eModified_A:
                    {
                        switch(expansion->m_basisKeyVector[edge_id].GetPointsType())
                        {
                        case LibUtilities::eGaussLobattoLegendre:
                        {
                            const LibUtilities::PointsKey pkey(numpoints,LibUtilities::eGaussLobattoLegendre);
                            return LibUtilities::BasisKey(expansion->m_basisKeyVector[0].GetBasisType(),nummodes,pkey);
                        }
                        break;
                        default:
                            ASSERTL0(false,"Unexpected points distribution");
                            // It doesn't matter what we return here
                            // since the ASSERT will stop execution.
                            // Just return something to prevent
                            // warnings messages.
                            const LibUtilities::PointsKey pkey(numpoints+1,LibUtilities::eGaussLobattoLegendre);
                            return LibUtilities::BasisKey(expansion->m_basisKeyVector[0].GetBasisType(),nummodes,pkey);
                            break;
                        }
                    }
                    break;

                default:
                    ASSERTL0(false,"Unexpected basis distribution");
                    // It doesn't matter what we return here since the
                    // ASSERT will stop execution.  Just return
                    // something to prevent warnings messages.
                    const LibUtilities::PointsKey pkey(numpoints+1,LibUtilities::eGaussLobattoLegendre);
                    return LibUtilities::BasisKey(expansion->m_basisKeyVector[0].GetBasisType(),nummodes,pkey);
                }
            }
            else
            {
                // Quadrilateral
                const LibUtilities::PointsKey pkey(numpoints,expansion->m_basisKeyVector[edge_id].GetPointsType());
                return LibUtilities::BasisKey(expansion->m_basisKeyVector[edge_id].GetBasisType(),nummodes,pkey);
            }
            
            ASSERTL0(false, "Unable to determine edge points type.");
            return LibUtilities::NullBasisKey;
        }
    }; //end of namespace
}; //end of namespace

//
// $Log: MeshGraph2D.cpp,v $
// Revision 1.41  2009/12/17 02:08:04  bnelson
// Fixed visual studio compiler warning.
//
// Revision 1.40  2009/10/22 17:41:47  cbiotto
// Update for variable order expansion
//
// Revision 1.39  2009/07/02 13:32:24  sehunchun
// *** empty log message ***
//
// Revision 1.38  2009/07/02 11:39:44  sehunchun
// Modification for 2D gemoetry embedded in 3D
//
// Revision 1.37  2009/05/01 13:23:21  pvos
// Fixed various bugs
//
// Revision 1.36  2009/04/20 16:13:23  sherwin
// Modified Import and Write functions and redefined how Expansion is used
//
// Revision 1.35  2009/01/21 16:59:03  pvos
// Added additional geometric factors to improve efficiency
//
// Revision 1.34  2009/01/12 10:26:59  pvos
// Added input tags for nodal expansions
//
// Revision 1.33  2008/09/12 11:26:19  pvos
// Updates for mappings in 3D
//
// Revision 1.32  2008/09/09 14:21:44  sherwin
// Updates for first working version of curved edges
//
// Revision 1.31  2008/08/14 22:11:03  sherwin
// Mods for HDG update
//
// Revision 1.30  2008/07/29 22:23:36  sherwin
// various mods for DG advection solver in Multiregions. Added virtual calls to Geometry, Geometry1D, 2D and 3D
//
// Revision 1.29  2008/05/29 19:07:39  delisi
// Removed the Write(..) methods, so it is only in the base MeshGraph class. Also, added a line to set the global ID of the geometry object for every element read in.
//
// Revision 1.28  2008/05/28 21:42:18  jfrazier
// Minor comment spelling change.
//
// Revision 1.27  2008/03/18 14:14:49  pvos
// Update for nodal triangular helmholtz solver
//
// Revision 1.26  2008/03/11 17:02:24  jfrazier
// Modified the GetElementsFromEdge to use the domain.
//
// Revision 1.25  2008/01/21 19:58:14  sherwin
// Updated so that QuadGeom and TriGeom have SegGeoms instead of EdgeComponents
//
// Revision 1.24  2007/12/17 20:27:24  sherwin
// Added normals to GeomFactors
//
// Revision 1.23  2007/12/14 17:39:18  jfrazier
// Fixed composite reader to handle ranges and comma separated lists.
//
// Revision 1.22  2007/12/11 21:51:53  jfrazier
// Updated 2d components so elements could be retrieved from edges.
//
// Revision 1.21  2007/12/04 03:02:26  jfrazier
// Changed to stringstream.
//
// Revision 1.20  2007/09/20 22:25:06  jfrazier
// Added expansion information to meshgraph class.
//
// Revision 1.19  2007/07/26 01:38:33  jfrazier
// Cleanup of some attribute reading code.
//
// Revision 1.18  2007/07/24 16:52:09  jfrazier
// Added domain code.
//
// Revision 1.17  2007/07/05 04:21:10  jfrazier
// Changed id format and propagated from 1d to 2d.
//
// Revision 1.16  2007/06/10 02:27:11  jfrazier
// Another checkin with an incremental completion of the boundary conditions reader
//
// Revision 1.15  2007/06/07 23:55:24  jfrazier
// Intermediate revisions to add parsing for boundary conditions file.
//
// Revision 1.14  2007/05/28 21:48:42  sherwin
// Update for 2D functionality
//
// Revision 1.13  2006/10/17 22:26:01  jfrazier
// Added capability to specify ranges in composite definition.
//
// Revision 1.12  2006/10/17 18:42:54  jfrazier
// Removed "NUMBER" attribute in items.
//
// Revision 1.11  2006/09/26 23:41:53  jfrazier
// Updated to account for highest level NEKTAR tag and changed the geometry tag to GEOMETRY.
//
// Revision 1.10  2006/08/24 18:50:00  jfrazier
// Completed error checking on permissable composite item combinations.
//
// Revision 1.9  2006/08/18 19:37:17  jfrazier
// *** empty log message ***
//
// Revision 1.8  2006/08/17 22:55:00  jfrazier
// Continued adding code to process composites in the mesh2d.
//
// Revision 1.7  2006/08/16 23:34:42  jfrazier
// *** empty log message ***
//
// Revision 1.6  2006/06/02 18:48:40  sherwin
// Modifications to make ProjectLoc2D run bit there are bus errors for order > 3
//
// Revision 1.5  2006/06/01 14:15:30  sherwin
// Added typdef of boost wrappers and made GeoFac a boost shared pointer.
//
// Revision 1.4  2006/05/30 14:00:04  sherwin
// Updates to make MultiRegions and its Demos work
//
// Revision 1.3  2006/05/23 19:56:33  jfrazier
// These build and run, but the expansion pieces are commented out
// because they would not run.
//
// Revision 1.2  2006/05/09 13:37:01  jfrazier
// Removed duplicate definition of shared vertex pointer.
//
// Revision 1.1  2006/05/04 18:59:01  kirby
// *** empty log message ***
//
// Revision 1.11  2006/04/11 23:18:11  jfrazier
// Completed MeshGraph2D for tri's and quads.  Not thoroughly tested.
//
// Revision 1.10  2006/04/09 02:08:35  jfrazier
// Added precompiled header.
//
// Revision 1.9  2006/04/04 23:12:37  jfrazier
// More updates to readers.  Still work to do on MeshGraph2D to store tris and quads.
//
// Revision 1.8  2006/03/25 00:58:29  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.7  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.6  2006/03/12 07:42:03  sherwin
//
// Updated member names and StdRegions call. Still has not been compiled
//
// Revision 1.5  2006/02/26 21:19:43  bnelson
// Fixed a variety of compiler errors caused by updates to the coding standard.
//
// Revision 1.4  2006/02/19 01:37:34  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
