////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/libs/SpatialDomains/MeshGraph1D.cpp,v $
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

namespace Nektar
{
    namespace SpatialDomains
    {
        MeshGraph1D::MeshGraph1D():
    m_geofac_defined(false)
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
        TiXmlElement *element = NULL;
        int err;    /// Error value returned by TinyXML.

        /// Look for segments in ELEMENT block.
        element = docHandle.FirstChildElement("MESH").FirstChildElement("ELEMENT").Element();

        ASSERTL0(element, "Unable to find ELEMENT tag in file.");
        TiXmlAttribute *attr = element->FirstAttribute();
        int numSegments = -1;

        while (attr)
        {
            std::string attrName(attr->Name());
            if (attrName == "NUMBER")
            {
                err = attr->QueryIntValue(&numSegments);
                if (err)
                {
                    numSegments = -1;
                }
            }
            // Get the next attribute.  Shouldn't be any more, but just in case.
            attr = attr->Next();
        }

        ASSERTL0(numSegments > 0, "Unable to read NUMBER attribute value.");
        int indx = -1;

        node = element->FirstChild();
        for (int i=0; i<numSegments; ++i)
        {
            /// First read index number
            /// Then read tagged element

            int nodeType = node->Type();
            if (nodeType == TiXmlNode::TEXT)
            {
                // Read index number
                indx = atoi(node->ToText()->Value());
                --i; // Don't count the index
            }
            else if (nodeType == TiXmlNode::COMMENT)
            {
                --numSegments; /// Leave it pointing where it is and skip it.
                continue;
            }
            else if (nodeType == TiXmlNode::ELEMENT)
            {
                /// All segments are embedded in <S></S> tags.
                TiXmlElement* segment = element->FirstChildElement();
                if (std::string(segment->Value())!= "S")
                {
		  NEKERROR(ErrorUtil::ewarning,(std::string("Element type not allowed: ") +
						segment->Value()).c_str());
                }
                else
                {
                    //// Get the entire data block then go through it one piece at a time.
                    std::string segmentData;

                    // Look through all immediate children and
                    // accumulate all text (type TEXT).  This is just
                    // in case it is broken up by other tagged
                    // elements.  Using element->GetText() will not
                    // work if text contains XML comments, for
                    // example.
                    TiXmlNode *child = node->FirstChild();
                    while(child)
                    {
                        if (child->Type() == TiXmlNode::TEXT)
                        {
                            segmentData += child->ToText()->Value();
                            segmentData += " "; // Don't jam together.
                        }
                        child = child->NextSibling();
                    }

                    // Get segment data from the data string.
                    int vertex1, vertex2;
                    std::istrstream vertexDataStrm(segmentData.c_str());

                    try
                    {
                        vertexDataStrm >> vertex1 >> vertex2;

                        ASSERTL0(!vertexDataStrm.fail(), 
				 (std::string("Unable to read segment data: ")
				  + std::string(segmentData) + ".").c_str());

                        SegGeomSharedPtr seg(new SegGeom(indx, 
				      GetVertex(vertex1), GetVertex(vertex2)));
                        m_seggeoms.push_back(seg);
                    }
                    catch(...)
                    {
                        NEKERROR(ErrorUtil::efatal, (std::string("Unable to read segment data: "
                            + segmentData + ".").c_str()));
                    }

                    indx = -1;
                }
            }
            else
            {
                ASSERTL1(false, "Unknown node type.");
            }

            node = node->NextSibling();
        }
    }

    void MeshGraph1D::Write(std::string &outfilename)
    {
    }

    // generate geometric factors based on MeshGraph information. 
    void MeshGraph1D::GenXGeoFac()
    {
        SegGeomVector::const_iterator def; 

        for(def = m_seggeoms.begin(); def != m_seggeoms.end(); ++def)
        {
            (*def)->SetXGeoFac((*def)->GenXGeoFac());
        }

        m_geofac_defined = true;
    }
    }; //end of namespace
}; //end of namespace

//
// $Log: MeshGraph1D.cpp,v $
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
