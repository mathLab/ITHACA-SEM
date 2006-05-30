////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/libs/SpatialDomains/MeshGraph2D.cpp,v $
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

#include <SpatialDomains/MeshGraph2D.h>
#include <SpatialDomains/SegGeom.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        MeshGraph2D::MeshGraph2D():
    m_geofac_defined(false)
    {
    }

    MeshGraph2D::~MeshGraph2D()
    {
    }

    void MeshGraph2D::Read(std::string &infilename)
    {
        SetFileName(infilename);
        TiXmlDocument doc(infilename);
        bool loadOkay = doc.LoadFile();

        ASSERTL0(loadOkay, (std::string("Unable to load file: ") + infilename).c_str());

        Read(doc);
    }

    void MeshGraph2D::ReadEdges(TiXmlDocument &doc)
    {
        /// We know we have it since we made it this far.
        TiXmlElement* mesh = doc.FirstChildElement("MESH");
        TiXmlElement* field = NULL;

        /// Look for elements in ELEMENT block.
        field = mesh->FirstChildElement("EDGE");

        ASSERTL0(field, "Unable to find EDGE tag in file.");

        TiXmlAttribute *attr = field->FirstAttribute();
        int numEdges = -1;
        int edgeNumber = 0;
        int nextEdgeNumber = 0;
        int err = 0;

        while (attr)
        {
            std::string attrName(attr->Name());
            if (attrName == "NUMBER")
            {
                err = attr->QueryIntValue(&numEdges);
                if (err)
                {
                    numEdges = -1;
                }
            }
            // Get the next attribute.  Shouldn't be any more, but just in case.
            attr = attr->Next();
        }

        ASSERTL0(numEdges > 0, "Unable to read NUMBER attribute value.");

        /// All elements are of the form: "ID <?> ... </?>", with ? being the element type.

        /// Read the ID field first.
        TiXmlNode *child = field->FirstChild();

        /// Since all edge data is one big text block, we need to accumulate
        /// all TEXT data and then parse it.  This approach effectively skips
        /// all comments or other node types since we only care about the
        /// edge list.  We cannot handle missing edge numbers as we could
        /// with missing element numbers due to the text block format.
        std::string edgeStr;

        while(child)
        {
            if (child->Type() == TiXmlNode::TEXT)
            {
                edgeStr += child->ToText()->ValueStr();
            }

            child = child->NextSibling();
        }

        /// Now parse out the edges, three fields at a time.
        int edgeid, vertex1, vertex2;
        std::istrstream edgeDataStrm(edgeStr.c_str());

        for (int i=0; i<numEdges; ++i)
        {
            try
            {
                edgeDataStrm >> edgeid >> vertex1 >> vertex2;

                ASSERTL0(!edgeDataStrm.fail(), (std::string("Unable to read edge data: ") + edgeStr).c_str());

                VertexComponentSharedPtr vertices[2] = {GetVertex(vertex1), GetVertex(vertex2)};

                EdgeComponentSharedPtr edge(new EdgeComponent(edgeid, m_MeshDimension, vertices));
                m_ecomps.push_back(edge);
            }
            catch(...)
            {
                NEKERROR(ErrorUtil::efatal, (std::string("Unable to read edge data: ") + edgeStr).c_str());
            }
        }
    }

    void MeshGraph2D::ReadElements(TiXmlDocument &doc)
    {
        /// We know we have it since we made it this far.
        TiXmlElement* mesh = doc.FirstChildElement("MESH");
        TiXmlElement* field = NULL;

        /// Look for elements in ELEMENT block.
        field = mesh->FirstChildElement("ELEMENT");

        ASSERTL0(field, "Unable to find ELEMENT tag in file.");

        TiXmlAttribute *attr = field->FirstAttribute();
        int numElements = -1;
        int elementNumber = 0;
        int nextElementNumber = 0;
        int err = 0;

        while (attr)
        {
            std::string attrName(attr->Name());
            if (attrName == "NUMBER")
            {
                err = attr->QueryIntValue(&numElements);
                if (err)
                {
                    numElements = -1;
                }
            }
            // Get the next attribute.  Shouldn't be any more, but just in case.
            attr = attr->Next();
        }

        ASSERTL0(numElements > 0, "Unable to read NUMBER attribute value.");

        /// All elements are of the form: "ID <?> ... </?>", with ? being the element type.

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

                ASSERTL0(elementType == "Q" || elementType == "T",
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

                if (elementType == "T")
                {
                    // Read three edge numbers
                    int edge1, edge2, edge3;
                    std::istrstream elementDataStrm(elementStr.c_str());

                    try
                    {
                        elementDataStrm >> edge1;
                        elementDataStrm >> edge2;
                        elementDataStrm >> edge3;

                        ASSERTL0(!elementDataStrm.fail(), (std::string("Unable to read element data for TRIANGLE: ") + elementStr).c_str());

                        /// Create a TriGeom to hold the new definition.
                        EdgeComponentSharedPtr edges[TriGeom::kNedges] = {GetEdgeComponent(edge1),GetEdgeComponent(edge2),GetEdgeComponent(edge3)};
                        StdRegions::EdgeOrientation edgeorient[TriGeom::kNedges] = {EdgeComponent::GetEdgeOrientation(*edges[0], *edges[1]),
                            EdgeComponent::GetEdgeOrientation(*edges[1], *edges[2]), EdgeComponent::GetEdgeOrientation(*edges[2], *edges[0])};
                        TriGeomSharedPtr trigeom(new TriGeom(edges, edgeorient));

                        m_trigeoms.push_back(trigeom);
                    }
                    catch(...)
                    {
                        NEKERROR(ErrorUtil::efatal,
                            (std::string("Unable to read element data for TRIANGLE: ") + elementStr).c_str());
                    }
                }
                else if (elementType == "Q")
                {
                    // Read four edge numbers
                    int edge1, edge2, edge3, edge4;
                    std::istrstream elementDataStrm(elementStr.c_str());

                    try
                    {
                        elementDataStrm >> edge1;
                        elementDataStrm >> edge2;
                        elementDataStrm >> edge3;
                        elementDataStrm >> edge4;

                        ASSERTL0(!elementDataStrm.fail(), (std::string("Unable to read element data for QUAD: ") + elementStr).c_str());

                        /// Create a QuadGeom to hold the new definition.
                        EdgeComponentSharedPtr edges[QuadGeom::kNedges] = {GetEdgeComponent(edge1),GetEdgeComponent(edge2),GetEdgeComponent(edge3)};
                        StdRegions::EdgeOrientation edgeorient[QuadGeom::kNedges] =
                        {
                            EdgeComponent::GetEdgeOrientation(*edges[0], *edges[1]),
                                EdgeComponent::GetEdgeOrientation(*edges[1], *edges[2]),
                                EdgeComponent::GetEdgeOrientation(*edges[2], *edges[3]),
                                EdgeComponent::GetEdgeOrientation(*edges[3], *edges[0])
                        };
                        QuadGeomSharedPtr quadgeom(new QuadGeom(edges, edgeorient));

                        m_quadgeoms.push_back(quadgeom);

                    }
                    catch(...)
                    {
                        NEKERROR(ErrorUtil::efatal, (std::string("Unable to read element data for QUAD: ") + elementStr).c_str());
                    }
                }
            }

            /// Keep looking
            child = child->NextSibling();
        }
    }

    // \brief Read segments (and general MeshGraph) given TiXmlDocument.
    void MeshGraph2D::Read(TiXmlDocument &doc)
    {
        // Read mesh first
        MeshGraph::Read(doc);

        TiXmlNode* node = NULL;
        TiXmlElement* mesh = NULL;

        /// Look for all geometry related data in MESH block.
        mesh = doc.FirstChildElement("MESH");

        ASSERTL0(mesh, "Unable to find MESH tag in file.");

        ReadEdges(doc);
        ReadElements(doc);

    }

    EdgeComponentSharedPtr MeshGraph2D::GetEdgeComponent(int eID)
    {
        EdgeComponentSharedPtr returnval;

        if (eID >= 0 && eID < int(m_ecomps.size()))
        {
            returnval = m_ecomps[eID];
        }

        return returnval;
    };

    void MeshGraph2D::Write(std::string &outfilename)
    {
    }
    }; //end of namespace
}; //end of namespace

//
// $Log: MeshGraph2D.cpp,v $
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
