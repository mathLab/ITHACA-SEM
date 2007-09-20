////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/MeshGraph.cpp,v $
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

#include <boost/foreach.hpp>

#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/ParseUtils.hpp>

// Use the stl version, primarily for string.
#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif

#include <tinyxml/tinyxml.h>
#include <string>
#include <strstream>

namespace Nektar
{
    namespace SpatialDomains
    {

        MeshGraph::MeshGraph():
            m_MeshDimension(3),
            m_SpaceDimension(3)
        {
        }

        // \brief Read will read the meshgraph vertices given a filename.
        void MeshGraph::Read(std::string &infilename)
        {
            SetFileName(infilename);
            TiXmlDocument doc(infilename);
            bool loadOkay = doc.LoadFile();

            std::string errstr = "Unable to load file: ";
            errstr += infilename;
            ASSERTL0(loadOkay, errstr.c_str());

            Read(doc);
        }

        // \brief Read will read the meshgraph vertices given a TiXmlDocument.
        void MeshGraph::Read(TiXmlDocument &doc)
        {
            TiXmlHandle docHandle(&doc);
            TiXmlNode* node = NULL;
            TiXmlElement* mesh = NULL;
            TiXmlElement* master = NULL;    // Master tag within which all data is contained.

            int err;    /// Error value returned by TinyXML.

            master = doc.FirstChildElement("NEKTAR");
            ASSERTL0(master, "Unable to find NEKTAR tag in file");

            // Find the Mesh tag and same the dim and space attributes
            mesh = master->FirstChildElement("GEOMETRY");

            ASSERTL0(mesh, "Unable to find GEOMETRY tag in file.");
            TiXmlAttribute *attr = mesh->FirstAttribute();

            // Initialize the mesh and space dimensions to 3 dimensions.
            // We want to do this each time we read a file, so it should
            // be done here and not just during class initialization.
            m_MeshDimension = 3;
            m_SpaceDimension = 3;

            while (attr)
            {
                std::string attrName(attr->Name());
                if (attrName == "DIM")
                {
                    err = attr->QueryIntValue(&m_MeshDimension);
                    ASSERTL1(err==TIXML_SUCCESS, "Unable to read mesh dimension.");
                }
                else if (attrName == "SPACE")
                {
                    err = attr->QueryIntValue(&m_SpaceDimension);
                    ASSERTL1(err==TIXML_SUCCESS, "Unable to read mesh dimension.");
                }
                else
                {
                    std::string errstr("Unknown attribute: ");
                    errstr += attrName;
                    ASSERTL1(false, errstr.c_str());
                }

                // Get the next attribute.
                attr = attr->Next();
            }

            ASSERTL1(m_MeshDimension<=m_SpaceDimension, "Mesh dimension greater than space dimension");

            // Now read the vertices
            TiXmlElement* element = mesh->FirstChildElement("VERTEX");
            ASSERTL0(element, "Unable to find mesh VERTEX tag in file.");

            TiXmlElement *vertex = element->FirstChildElement("V");

            int indx;
            int nextVertexNumber = -1;

            while (vertex)
            {
                nextVertexNumber++;

                TiXmlAttribute *vertexAttr = vertex->FirstAttribute();
                std::string attrName(vertexAttr->Name());

                ASSERTL0(attrName == "ID", (std::string("Unknown attribute name: ") + attrName).c_str());

                err = vertexAttr->QueryIntValue(&indx);
                ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute ID.");
                ASSERTL0(indx == nextVertexNumber, "Element IDs must begin with zero and be sequential.");

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
                double xval, yval, zval;
                std::istrstream vertexDataStrm(vertexBodyStr.c_str());

                try
                {
                    while(!vertexDataStrm.fail())
                    {
                        vertexDataStrm >> xval >> yval >> zval;

                        // Need to check it here because we may not be good after the read
                        // indicating that there was nothing to read.
                        if (!vertexDataStrm.fail())
                        {
                            VertexComponentSharedPtr vert(MemoryManager<VertexComponent>::AllocateSharedPtr(m_MeshDimension, indx, xval, yval, zval));
                            m_vertset.push_back(vert);
                        }
                    }
                }
                catch(...)
                {
                    ASSERTL0(false, "Unable to read VERTEX data.");
                }
 
                vertex = vertex->NextSiblingElement("V");
            }

       }

        GeometrySharedPtr MeshGraph::GetCompositeItem(int whichComposite, int whichItem)
        {
            GeometrySharedPtr returnval;
            bool error = false;

            if (whichComposite >= 0 && whichComposite < int(m_MeshCompositeVector.size()))
            {
                if (whichItem >= 0 && whichItem < int(m_MeshCompositeVector[whichComposite]->size()))
                {
                    returnval = m_MeshCompositeVector[whichComposite]->at(whichItem);
                }
                else
                {
                    error = true;
                }
            }
            else
            {
                error = true;
            }

            if (error)
            {
                std::ostringstream errStream;
                errStream << "Unable to access composite item [" << whichComposite << "][" << whichItem << "].";
                
                std::string testStr = errStream.str();

                NEKERROR(ErrorUtil::efatal, testStr.c_str());
            }

            return returnval;
        }

        void MeshGraph::ReadDomain(TiXmlDocument &doc)
        {
            TiXmlHandle docHandle(&doc);

            TiXmlElement* mesh = docHandle.FirstChildElement("NEKTAR").FirstChildElement("GEOMETRY").Element();
            TiXmlElement* domain = NULL;

            ASSERTL0(mesh, "Unable to find GEOMETRY tag in file.");

            /// Look for data in DOMAIN block.
            domain = mesh->FirstChildElement("DOMAIN");

            ASSERTL0(domain, "Unable to find DOMAIN tag in file.");

            // find the non comment portion of the body.
            TiXmlNode* elementChild = domain->FirstChild();
            while(elementChild && elementChild->Type() != TiXmlNode::TEXT)
            {
                elementChild = elementChild->NextSibling();
            }

            ASSERTL0(elementChild, "Unable to read DOMAIN body.");
            std::string elementStr = elementChild->ToText()->ValueStr();

            elementStr = elementStr.substr(elementStr.find_first_not_of(" "));

            std::string::size_type indxBeg = elementStr.find_first_of('[') + 1;
            std::string::size_type indxEnd = elementStr.find_last_of(']') - 1;
            std::string indxStr = elementStr.substr(indxBeg, indxEnd - indxBeg + 1);

            ASSERTL0(!indxStr.empty(), "Unable to read domain's composite index (index missing?).");

            // Read the domain composites.
            // Parse the composites into a list.
            GetCompositeList(indxStr, m_Domain);
            ASSERTL0(!m_Domain.empty(), (std::string("Unable to obtain domain's referenced composite: ") + indxStr).c_str());
        }

        void MeshGraph::GetCompositeList(const std::string &compositeStr, CompositeVector &compositeVector) const
        {
            // Parse the composites into a list.
            typedef vector<unsigned int> SeqVector;
            SeqVector seqVector;
            bool parseGood = ParseUtils::GenerateSeqVector(compositeStr.c_str(), seqVector);

            ASSERTL0(parseGood && !seqVector.empty(), (std::string("Unable to read composite index range: ") + compositeStr).c_str());

            SeqVector addedVector;    // Vector of those composites already added to compositeVector;
            for (SeqVector::iterator iter = seqVector.begin(); iter != seqVector.end(); ++iter)
            {
                // Only add a new one if it does not already exist in vector.
                // Can't go back and delete with a vector, so prevent it from
                // being added in the first place.
                if (std::find(addedVector.begin(), addedVector.end(), *iter) == addedVector.end())
                {
                    addedVector.push_back(*iter);
                    Composite composite = GetComposite(*iter);
                    CompositeVector::iterator compIter;
                    if (composite)
                    {
                        compositeVector.push_back(composite);
                    }
                    else
                    {
                        char str[64];
                        ::sprintf(str, "%d", *iter);
                        NEKERROR(ErrorUtil::ewarning, (std::string("Undefined composite: ") + str).c_str());

                    }
                }
            }
        }

        void MeshGraph::Write(std::string &outfilename)
        {
        }

        MeshGraph::~MeshGraph()
        {
        }
    }; //end of namespace
}; //end of namespace

//
// $Log: MeshGraph.cpp,v $
// Revision 1.9  2007/09/03 17:05:01  jfrazier
// Cleanup and addition of composite range in domain specification.
//
// Revision 1.8  2007/07/24 16:52:08  jfrazier
// Added domain code.
//
// Revision 1.7  2007/06/10 02:27:10  jfrazier
// Another checkin with an incremental completion of the boundary conditions reader.
//
// Revision 1.6  2007/03/29 19:25:10  bnelson
// *** empty log message ***
//
// Revision 1.5  2006/10/17 18:42:54  jfrazier
// Removed "NUMBER" attribute in items.
//
// Revision 1.4  2006/09/26 23:41:52  jfrazier
// Updated to account for highest level NEKTAR tag and changed the geometry tag to GEOMETRY.
//
// Revision 1.3  2006/08/24 18:50:00  jfrazier
// Completed error checking on permissable composite item combinations.
//
// Revision 1.2  2006/05/09 13:37:01  jfrazier
// Removed duplicate definition of shared vertex pointer.
//
// Revision 1.1  2006/05/04 18:59:02  kirby
// *** empty log message ***
//
// Revision 1.14  2006/04/09 02:08:35  jfrazier
// Added precompiled header.
//
// Revision 1.13  2006/04/04 23:12:37  jfrazier
// More updates to readers.  Still work to do on MeshGraph2D to store tris and quads.
//
// Revision 1.12  2006/03/25 00:58:29  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.11  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.10  2006/02/26 21:19:43  bnelson
// Fixed a variety of compiler errors caused by updates to the coding standard.
//
// Revision 1.9  2006/02/19 01:37:33  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
