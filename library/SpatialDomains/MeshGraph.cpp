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
#include <cstring>
#include <sstream>

#include <SpatialDomains/MeshGraph1D.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <SpatialDomains/MeshGraph3D.h>

// These are required for the Write(...) and Import(...) functions.
#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/binary_from_base64.hpp>
#include <boost/archive/iterators/transform_width.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/zlib.hpp> 
#include <boost/iostreams/filtering_stream.hpp>

namespace Nektar
{
    namespace SpatialDomains
    {

        MeshGraph::MeshGraph():
    m_MeshDimension(3),
        m_SpaceDimension(3)
    {
    }

    boost::shared_ptr<MeshGraph> MeshGraph::Read(std::string &infilename)
    {
        boost::shared_ptr<MeshGraph> returnval;

        MeshGraph mesh;

        mesh.ReadGeometry(infilename);
        int meshDim = mesh.GetMeshDimension();

        switch(meshDim)
        {
        case 1:
            returnval = MemoryManager<MeshGraph1D>::AllocateSharedPtr();
            break;
            
        case 2:
           returnval = MemoryManager<MeshGraph2D>::AllocateSharedPtr();
           break;

        case 3:
           returnval = MemoryManager<MeshGraph3D>::AllocateSharedPtr();
           break;

        default:
            std::string err = "Invalid mesh dimension: ";
            std::stringstream strstrm;
            strstrm << meshDim;
            err += strstrm.str();
            NEKERROR(ErrorUtil::efatal, err.c_str());
        }

        if (returnval)
        {
            returnval->ReadGeometry(infilename);
            returnval->ReadExpansions(infilename);
        }

        return returnval;
    }


    // \brief Read will read the meshgraph vertices given a filename.
    void MeshGraph::ReadGeometry(std::string &infilename)
    {
        TiXmlDocument doc(infilename);
        bool loadOkay = doc.LoadFile();

        std::string errstr = "Unable to load file: ";
        errstr += infilename;
        ASSERTL0(loadOkay, errstr.c_str());

        ReadGeometry(doc);
    }

    // \brief Read will read the meshgraph vertices given a TiXmlDocument.
    void MeshGraph::ReadGeometry(TiXmlDocument &doc)
    {
        TiXmlHandle docHandle(&doc);
        TiXmlNode* node = NULL;
        TiXmlElement* mesh = NULL;
        TiXmlElement* master = NULL;    // Master tag within which all data is contained.

        int err;    /// Error value returned by TinyXML.

        master = doc.FirstChildElement("NEKTAR");
        ASSERTL0(master, "Unable to find NEKTAR tag in file.");

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
            std::istringstream vertexDataStrm(vertexBodyStr.c_str());

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

    // \brief Read the expansions given the XML file path.
    void MeshGraph::ReadExpansions(std::string &infilename)
    {
        TiXmlDocument doc(infilename);
        bool loadOkay = doc.LoadFile();

        std::string errstr = "Unable to load file: ";
        errstr += infilename;
        ASSERTL0(loadOkay, errstr.c_str());

        ReadExpansions(doc);
    }

    // \brief Read the expansions given the XML document reference.
    void MeshGraph::ReadExpansions(TiXmlDocument &doc)
    {
        TiXmlElement *master = doc.FirstChildElement("NEKTAR");
        ASSERTL0(master, "Unable to find NEKTAR tag in file.");

        // Find the Expansions tag
         TiXmlElement *expansionTypes = master->FirstChildElement("EXPANSIONS");
        ASSERTL0(expansionTypes, "Unable to find EXPANSIONS tag in file.");

        if (expansionTypes)
        {
            /// Expansiontypes will contain composite, nummodes, and expansiontype
            /// (eModified, or eOrthogonal)

            // Need a vector of all elements and their associated expansion information.
            // Default all elements in the domain to linear (2 modes) and Modified type.
            const CompositeVector &domain = this->GetDomain();
            CompositeVector::const_iterator compIter;

            for (compIter = domain.begin(); compIter != domain.end(); ++compIter)
            {
                boost::shared_ptr<GeometryVector> geomVectorShPtr = *compIter;
                GeometryVectorIter geomIter;
                for (geomIter = geomVectorShPtr->begin(); geomIter != geomVectorShPtr->end(); ++geomIter)
                {
                    // Make sure we only have one instance of the GeometrySharedPtr stored in the list.
                    ExpansionVector::iterator elemIter;
                    for (elemIter = m_ExpansionVector.begin(); elemIter != m_ExpansionVector.end(); ++elemIter)
                    {
                        if ((*elemIter)->m_GeomShPtr == *geomIter)
                        {
                            break;
                        }
                    }

                    // Not found in list.
                    if (elemIter == m_ExpansionVector.end())
                    {
                        Equation eqn("2");
                        ExpansionType type = eModified;
                        ExpansionShPtr expansionElementShPtr =
                            MemoryManager<Expansion>::AllocateSharedPtr(*geomIter, eqn, type);
                        m_ExpansionVector.push_back(expansionElementShPtr);
                    }
                }
            }

            // Clear the default linear expansion over the domain.
            TiXmlElement *expansion = expansionTypes->FirstChildElement("E");

            while (expansion)
            {
                /// Mandatory components...optional are to follow later.
                std::string compositeStr = expansion->Attribute("COMPOSITE");
                ASSERTL0(compositeStr.length() > 3, "COMPOSITE must be specified in expansion definition");
                int beg = compositeStr.find_first_of("[");
                int end = compositeStr.find_first_of("]");
                std::string compositeListStr = compositeStr.substr(beg+1,end-beg-1);

                CompositeVector compositeVector;
                GetCompositeList(compositeListStr, compositeVector);

                std::string nummodesStr = expansion->Attribute("NUMMODES");
                ASSERTL0(!nummodesStr.empty(), "NUMMODES must be specified in expansion definition");
                Equation nummodesEqn(nummodesStr);

                std::string typeStr = expansion->Attribute("TYPE");
                ASSERTL0(!typeStr.empty(), "TYPE must be specified in "
                    "expansion definition");

                ExpansionType type;
                const std::string* begStr = kExpansionTypeStr;
                const std::string* endStr = kExpansionTypeStr+eExpansionTypeSize;
                const std::string* expStr = std::find(begStr, endStr, typeStr);

                ASSERTL0(expStr != endStr, "Invalid expansion type.");
                type = (ExpansionType)(expStr - begStr);

                // Now have composite, modes, and type.
                // Cycle through all composites for the geomShPtrs and set the modes and types for the
                // elements contained in the element list.
                CompositeVectorIter compVecIter;
                for (compVecIter = compositeVector.begin(); compVecIter != compositeVector.end(); ++compVecIter)
                {
                    GeometryVectorIter geomVecIter;
                    for (geomVecIter = (*compVecIter)->begin(); geomVecIter != (*compVecIter)->end(); ++geomVecIter)
                    {
                        ExpansionVectorIter expVecIter;
                        for (expVecIter = m_ExpansionVector.begin(); expVecIter != m_ExpansionVector.end(); ++expVecIter)
                        {
                            if (*geomVecIter == (*expVecIter)->m_GeomShPtr)
                            {
                                (*expVecIter)->m_ExpansionType = type;
                                (*expVecIter)->m_NumModesEqn = nummodesEqn;
                                break;
                            }
                        }
                    }
                }

                expansion = expansion->NextSiblingElement("E");
            }
        }
    }


    void MeshGraph::ReadCurves(TiXmlDocument &doc)
    {
        /// We know we have it since we made it this far.
        TiXmlHandle docHandle(&doc);
        TiXmlElement* mesh = docHandle.FirstChildElement("NEKTAR").FirstChildElement("GEOMETRY").Element();
        TiXmlElement* field = NULL;

        int err;

        /// Look for elements in CURVE block.
        field = mesh->FirstChildElement("CURVED");

        if(!field) //return if no curved entities
        {
            return; 
        }

        /// All curves are of the form: "<? ID="#" TYPE="GLL OR other
        /// points type" NUMPOINTS="#"> ... </?>", with ? being an
        /// element type (either E or F).

        TiXmlElement *edgelement = field->FirstChildElement("E");

        int edgeindx, edgeid;
        int nextEdgeNumber = -1;

        while(edgelement)
        {
            /// These should be ordered.
            nextEdgeNumber++;

            std::string edge(edgelement->ValueStr());
            ASSERTL0(edge == "E", (std::string("Unknown 3D curve type:") + edge).c_str());

            /// Read id attribute.
            err = edgelement->QueryIntAttribute("ID", &edgeindx);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read curve attribute ID.");
            ASSERTL0(edgeindx == nextEdgeNumber, "Curve IDs must begin with zero and be sequential.");

            /// Read edge id attribute.
            err = edgelement->QueryIntAttribute("EDGEID", &edgeid);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read curve attribute EDGEID.");

            /// Read text edgelement description.
            std::string elementStr;
            TiXmlNode* elementChild = edgelement->FirstChild();

            while(elementChild)
            {
                // Accumulate all non-comment element data
                if (elementChild->Type() == TiXmlNode::TEXT)
                {
                    elementStr += elementChild->ToText()->ValueStr();
                    elementStr += " ";
                }
                elementChild = elementChild->NextSibling();
            }

            ASSERTL0(!elementStr.empty(), "Unable to read curve description body.");

            /// Parse out the element components corresponding to type of element.
           if (edge == "E")
            {
                int numPts=0;
                // Determine the points type
                std::string typeStr = edgelement->Attribute("TYPE");
                ASSERTL0(!typeStr.empty(), "TYPE must be specified in " "points definition");

                LibUtilities::PointsType type;
                const std::string* begStr = LibUtilities::kPointsTypeStr;
                const std::string* endStr = LibUtilities::kPointsTypeStr + LibUtilities::SIZE_PointsType;
                const std::string* ptsStr = std::find(begStr, endStr, typeStr);

                ASSERTL0(ptsStr != endStr, "Invalid points type.");
                type = (LibUtilities::PointsType)(ptsStr - begStr);

                //Determine the number of points
                err = edgelement->QueryIntAttribute("NUMPOINTS", &numPts);
                ASSERTL0(err == TIXML_SUCCESS, "Unable to read curve attribute NUMPOINTS.");
                CurveSharedPtr curve(MemoryManager<Curve>::AllocateSharedPtr(edgeid, type));
                
                // Read points (x, y, z)
                double xval, yval, zval;
                std::istringstream elementDataStrm(elementStr.c_str());
                try
                {
                    while(!elementDataStrm.fail())
                    {
                        elementDataStrm >> xval >> yval >> zval;

                        // Need to check it here because we may not be
                        // good after the read indicating that there
                        // was nothing to read.
                        if (!elementDataStrm.fail())
                        {
                            VertexComponentSharedPtr vert(MemoryManager<VertexComponent>::AllocateSharedPtr(m_MeshDimension, edgeindx, xval, yval, zval));

                            curve->m_points.push_back(vert);
                        }

                    }
                }
                catch(...)
                {
                    NEKERROR(ErrorUtil::efatal,
                    (std::string("Unable to read curve data for EDGE: ") + elementStr).c_str());

                }

                ASSERTL0(curve->m_points.size() == numPts,"Number of points specificed by attribute NUMPOINTS is different from number of points in list");

                m_curvededges.push_back(curve);

                edgelement = edgelement->NextSiblingElement("E");

            } // end if-loop

        } // end while-loop


        TiXmlElement *facelement = field->FirstChildElement("F");
        int faceindx, faceid;
        int nextFaceNumber = -1;
        
        while(facelement)
        {
            cout << "facelement = " << *facelement << endl;
             
            /// These should be ordered.
            nextFaceNumber++;

            std::string face(facelement->ValueStr());
            ASSERTL0(face == "F", (std::string("Unknown 3D curve type: ") + face).c_str());

            /// Read id attribute.
            err = facelement->QueryIntAttribute("ID", &faceindx);

            ASSERTL0(err == TIXML_SUCCESS, "Unable to read curve attribute ID.");
            ASSERTL0(faceindx == nextFaceNumber, "Face IDs must begin with zero and be sequential.");


            /// Read edge id attribute.
            err = edgelement->QueryIntAttribute("FACEID", &faceid);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read curve attribute FACEID.");

            /// Read text face element description.
            std::string elementStr;
            TiXmlNode* elementChild = facelement->FirstChild();
          
            while(elementChild)
            {
                // Accumulate all non-comment element data
                if (elementChild->Type() == TiXmlNode::TEXT)
                {
                    elementStr += elementChild->ToText()->ValueStr();
                    elementStr += " ";
                }
                elementChild = elementChild->NextSibling();
            }

            ASSERTL0(!elementStr.empty(), "Unable to read curve description body.");
            
            
            /// Parse out the element components corresponding to type of element.
            
            if(face == "F") 
            {
                std::string typeStr = facelement->Attribute("TYPE");
                ASSERTL0(!typeStr.empty(), "TYPE must be specified in " "points definition");
                LibUtilities::PointsType type;
                const std::string* begStr = LibUtilities::kPointsTypeStr;
                const std::string* endStr = LibUtilities::kPointsTypeStr + LibUtilities::SIZE_PointsType;
                const std::string* ptsStr = std::find(begStr, endStr, typeStr);
                
                ASSERTL0(ptsStr != endStr, "Invalid points type.");
                type = (LibUtilities::PointsType)(ptsStr - begStr);
                
                std::string numptsStr = facelement->Attribute("NUMPOINTS");
                ASSERTL0(!numptsStr.empty(), "NUMPOINTS must be specified in points definition");
                int numPts=0;
                std::strstream s;
                s << numptsStr;
                s >> numPts;
                
                
                CurveSharedPtr curve(MemoryManager<Curve>::AllocateSharedPtr(faceid, type));
                
                cout << "numPts = " << numPts << endl;
                ASSERTL0(numPts >= 3, "NUMPOINTS for face must be greater than 2");
                
                if(numPts == 3)
                {
                    ASSERTL0(ptsStr != endStr, "Invalid points type.");
                }
                
                // Read points (x, y, z)
                double xval, yval, zval;
                std::istringstream elementDataStrm(elementStr.c_str());
                try
                {
                    while(!elementDataStrm.fail())
                    {
                        elementDataStrm >> xval >> yval >> zval;
                        
                        // Need to check it here because we may not be good after the read
                        // indicating that there was nothing to read.
                        if (!elementDataStrm.fail())
                        {
                            VertexComponentSharedPtr vert(MemoryManager<VertexComponent>::AllocateSharedPtr(m_MeshDimension, faceindx, xval, yval, zval));
                            curve->m_points.push_back(vert);
                            
                        }
                        
                        cout << "xval = " << xval << "  yval = " << yval <<"  zval = " << zval << endl;
                        cout << endl;
                        
                    }
                }
                catch(...)
                {
                    NEKERROR(ErrorUtil::efatal,
                             (std::string("Unable to read curve data for FACE: ") + elementStr).c_str());
                    
                }
                m_curvedfaces.push_back(curve);
                
                facelement = facelement->NextSiblingElement("F");
                
            } // end if-loop            
        } // end while-loop
        
    } // end of ReadCurves()

        
    void MeshGraph::ReadCurves(std::string &infilename)
    {
        TiXmlDocument doc(infilename);
        bool loadOkay = doc.LoadFile();

        std::string errstr = "Unable to load file: ";
        errstr += infilename;
        ASSERTL0(loadOkay, errstr.c_str());

        ReadCurves(doc);
    }


    LibUtilities::BasisKey MeshGraph::GetBasisKey(ExpansionShPtr in,
        const int flag)
    {
        int order = (int) in->m_NumModesEqn.Evaluate();

        switch(in->m_ExpansionType)
        {
        case eModified:
            switch (flag)
            {
            case 0:
                {
                    const LibUtilities::PointsKey pkey(order+1,LibUtilities::eGaussLobattoLegendre);
                    return LibUtilities::BasisKey(LibUtilities::eModified_A,order,pkey);
                }
                break;
            case 1:
                {
                    const LibUtilities::PointsKey pkey(order,LibUtilities::eGaussRadauMAlpha1Beta0);
                    return LibUtilities::BasisKey(LibUtilities::eModified_B,order,pkey);
                }
                break;
            case 2:
                {
                    const LibUtilities::PointsKey pkey(order,LibUtilities::eGaussRadauMAlpha2Beta0);
                    return LibUtilities::BasisKey(LibUtilities::eModified_C,order,pkey);
                }
                break;
            default:
                ASSERTL0(false,"invalid value to flag");
                break;
            }
            break;
        case eNodal:
            {
                SegGeomSharedPtr segment = boost::dynamic_pointer_cast<SegGeom>(in->m_GeomShPtr);
                TriGeomSharedPtr triangle = boost::dynamic_pointer_cast<TriGeom>(in->m_GeomShPtr);
                QuadGeomSharedPtr quadrilateral = boost::dynamic_pointer_cast<QuadGeom>(in->m_GeomShPtr);

                if(segment || quadrilateral)
                {
                    const LibUtilities::PointsKey pkey(order+1,LibUtilities::eGaussLobattoLegendre);
                    return LibUtilities::BasisKey(LibUtilities::eGLL_Lagrange,order,pkey);
                }
                else if(triangle)
                {
                    switch (flag)
                    {
                    case 0:
                        {
                            const LibUtilities::PointsKey pkey(order+1,LibUtilities::eGaussLobattoLegendre);
                            return LibUtilities::BasisKey(LibUtilities::eOrtho_A,order,pkey);
                        }
                        break;
                    case 1:
                        {
                            const LibUtilities::PointsKey pkey(order,LibUtilities::eGaussRadauMAlpha1Beta0);
                            return LibUtilities::BasisKey(LibUtilities::eOrtho_B,order,pkey);
                        }
                        break;
                    case 2:
                        {
                            const LibUtilities::PointsKey pkey(order,LibUtilities::eGaussRadauMAlpha2Beta0);
                            return LibUtilities::BasisKey(LibUtilities::eOrtho_C,order,pkey);
                        }
                        break;
                    default:
                        ASSERTL0(false,"invalid value to flag");
                        break;
                    }
                }
                else
                {
                    NEKERROR(ErrorUtil::efatal,"Unrecognised Geometry");
                }
            }
            break;
        case eOrthogonal:
            switch (flag)
            {
            case 0:
                {
                    const LibUtilities::PointsKey pkey(order+1,LibUtilities::eGaussLobattoLegendre);
                    return LibUtilities::BasisKey(LibUtilities::eOrtho_A,order,pkey);
                }
                break;
            case 1:
                {
                    const LibUtilities::PointsKey pkey(order,LibUtilities::eGaussRadauMAlpha1Beta0);
                    return LibUtilities::BasisKey(LibUtilities::eOrtho_B,order,pkey);
                }
                break;
            case 2:
                {
                    const LibUtilities::PointsKey pkey(order,LibUtilities::eGaussRadauMAlpha2Beta0);
                    return LibUtilities::BasisKey(LibUtilities::eOrtho_C,order,pkey);
                }
                break;
            default:
                ASSERTL0(false,"invalid value to flag");
                break;
            }
            break;
        default:
            ASSERTL0(false,"expansion type unknonw");
            break;
        }

        // This never gets hit, but keeps the compiler happy.
        // Since the default cases above don't return anything
        // the compiler complains if this is not here.  It is
        // more proper than, say, suppressing the warning.
        return LibUtilities::NullBasisKey;
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

	void MeshGraph::Write(std::string &outfilename, FieldDefinitions &fielddefs, std::vector<double> &fielddata)
    {
		// Calculate serialization size.
		unsigned int size = fielddefs.m_Fields.size() * fielddefs.m_NumModes.size() * fielddefs.m_Elements.size();
		ASSERTL0(fielddefs.m_Elements.size() > 0, "Fielddefs vector must contain at least one element.");
		ASSERTL0(fielddata.size() == size, "Invalid size of fielddata vector.");
		ASSERTL0(fielddata.size() > 0, "Fielddata vector must contain at least one value.");

		// Calculate the attributes.
		GeomShapeType shapeType = fielddefs.m_Elements[0]->GetGeomShapeType();
		std::string typeString = GeomShapeTypeMap[shapeType];
		std::string idString;
		{
			std::stringstream idStringStream;
			bool first = true;
			for (std::vector<Geometry>::size_type i = 0; i < fielddefs.m_Elements.size(); i++)
			{
				ASSERTL0(fielddefs.m_Elements[i]->GetGeomShapeType() == shapeType,
					"All elements must be the same shape.");
				if (!first) idStringStream << ",";
				idStringStream << fielddefs.m_Elements[i]->GetGlobalID();
				first = false;
			}
			idString = idStringStream.str();
		}
		std::string basisString;
		{
			std::stringstream basisStringStream;
			bool first = true;
			for (std::vector<LibUtilities::BasisType>::size_type i = 0; i < fielddefs.m_Basis.size(); i++)
			{
				if (!first) basisStringStream << ",";
				basisStringStream << LibUtilities::BasisTypeMap[fielddefs.m_Basis[i]];
				first = false;
			}
			basisString = basisStringStream.str();
		}
		std::string numModesString;
		{
			std::stringstream numModesStringStream;
			bool first = true;
			for (std::vector<int>::size_type i = 0; i < fielddefs.m_NumModes.size(); i++)
			{
				if (!first) numModesStringStream << ",";
				numModesStringStream << fielddefs.m_NumModes[i];
				first = false;
			}
			numModesString = numModesStringStream.str();
		}
		std::string fieldsString;
		{
			std::stringstream fieldsStringStream;
			bool first = true;
			for (std::vector<int>::size_type i = 0; i < fielddefs.m_Fields.size(); i++)
			{
				if (!first) fieldsStringStream << ",";
				fieldsStringStream << fielddefs.m_Fields[i];
				first = false;
			}
			fieldsString = fieldsStringStream.str();
		}

		std::string compressedDataString;
		{
			// Serialize the fielddata vector to the stringstream.
			std::stringstream archiveStringStream(std::string((char*)&fielddata[0], sizeof(fielddata[0])/sizeof(char)*fielddata.size()));

			// Compress the serialized data.
			std::stringstream compressedData;
			{
				boost::iostreams::filtering_streambuf<boost::iostreams::input> out; 
				out.push(boost::iostreams::zlib_compressor());
				out.push(archiveStringStream);
				boost::iostreams::copy(out, compressedData);
			}

			// If the string length is not divisible by 3, pad it. There is a bug
			// in transform_width that will make it reference past the end and crash.
			switch (compressedData.str().length() % 3)
			{
			case 1:
				compressedData << '\0';
			case 2:
				compressedData << '\0';
				break;
			}
			compressedDataString = compressedData.str();
		}

		// Convert from binary to base64.
		typedef boost::archive::iterators::base64_from_binary<
					boost::archive::iterators::transform_width<
						std::string::const_iterator, 6, 8> > base64_t;
		std::string base64string(base64_t(compressedDataString.begin()), base64_t(compressedDataString.end()));

		// Create the XML output.
		std::ofstream xmlFile(outfilename.c_str(), std::ios::out | std::ios::trunc);
		ASSERTL0(xmlFile, std::string("Unable to save file: ").append(outfilename));
		xmlFile << "<?xml version=\"1.0\" encoding=\"utf-8\" ?>" << std::endl;
		xmlFile << "<nektar>" << std::endl;
		xmlFile << "    <element ";
		xmlFile << "id=\"" << idString << "\" ";
		xmlFile << "type=\"" << typeString << "\" ";
		xmlFile << "basis=\"" << basisString << "\" ";
		xmlFile << "numModes=\"" << numModesString << "\" ";
		xmlFile << "fields=\"" << fieldsString << "\">";
		xmlFile << base64string << "</element>" << std::endl;
		xmlFile << "</nektar>" << std::endl;
		xmlFile.flush();
		xmlFile.close();
    }

	void MeshGraph::Import(std::string &infilename, std::vector<FieldDefinitions> &fielddefs, std::vector<std::vector<double> > &fielddata)
	{
		ASSERTL1(fielddefs.size() == 0, "Expected an empty fielddefs vector.");
		ASSERTL1(fielddata.size() == 0, "Expected an empty fielddata vector.");

		TiXmlDocument doc(infilename);
		bool loadOkay = doc.LoadFile();

		std::stringstream errstr;
		errstr << "Unable to load file: " << infilename << std::endl;
		errstr << "Reason: " << doc.ErrorDesc() << std::endl;
		errstr << "Position: Line " << doc.ErrorRow() << ", Column " << doc.ErrorCol() << std::endl;
		ASSERTL0(loadOkay, errstr.str().c_str());

		TiXmlHandle docHandle(&doc);
		TiXmlElement* master = NULL;    // Master tag within which all data is contained.

		master = doc.FirstChildElement("nektar");
		ASSERTL0(master, "Unable to find nektar tag in file.");

		// Loop through all nektar tags, finding all of the element tags.
		while (master)
		{
			TiXmlElement* element = master->FirstChildElement("element");
			ASSERTL0(element, "Unable to find element tag within nektar tag.");

			while (element)
			{
				// Extract the attributes.
				std::string idString;
				std::string typeString;
				std::string basisString;
				std::string numModesString;
				std::string fieldsString;
				TiXmlAttribute *attr = element->FirstAttribute();
				while (attr)
				{
					std::string attrName(attr->Name());
					if (attrName == "id")
						idString.insert(0, attr->Value());
					else if (attrName == "type")
						typeString.insert(0, attr->Value());
					else if (attrName == "basis")
						basisString.insert(0, attr->Value());
					else if (attrName == "numModes")
						numModesString.insert(0, attr->Value());
					else if (attrName == "fields")
						fieldsString.insert(0, attr->Value());
					else
					{
						std::string errstr("Unknown attribute: ");
						errstr += attrName;
						ASSERTL1(false, errstr.c_str());
					}

					// Get the next attribute.
					attr = attr->Next();
				}

				// Extract the body, which the "data".
				TiXmlNode* elementChild = element->FirstChild();
				ASSERTL0(elementChild, "Unable to extract the data from the element tag.");
				std::string elementStr;
				while(elementChild)
				{
					if (elementChild->Type() == TiXmlNode::TEXT)
					{
						elementStr += elementChild->ToText()->ValueStr();
					}
					elementChild = elementChild->NextSibling();
				}

				// Convert from base64 to binary.
				typedef boost::archive::iterators::transform_width<
					boost::archive::iterators::binary_from_base64<
						std::string::const_iterator>, 8, 6 > binary_t;
				std::stringstream elementCompressedData(std::string(binary_t(elementStr.begin()), binary_t(elementStr.end())));

				// Decompress the binary data.
				std::stringstream elementDecompressedData;
				boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
				in.push(boost::iostreams::zlib_decompressor());
				in.push(elementCompressedData);
				try
				{
					boost::iostreams::copy(in, elementDecompressedData);
				}
				catch (boost::iostreams::zlib_error e)
				{
					if (e.error() == boost::iostreams::zlib::stream_end)
					{
						ASSERTL0(false, "Stream end zlib error");
					}
					else if (e.error() == boost::iostreams::zlib::stream_error)
					{
						ASSERTL0(false, "Stream zlib error");
					}
					else if (e.error() == boost::iostreams::zlib::version_error)
					{
						ASSERTL0(false, "Version zlib error");
					}
					else if (e.error() == boost::iostreams::zlib::data_error)
					{
						ASSERTL0(false, "Data zlib error");
					}
					else if (e.error() == boost::iostreams::zlib::mem_error)
					{
						ASSERTL0(false, "Memory zlib error");
					}
					else if (e.error() == boost::iostreams::zlib::buf_error)
					{
						ASSERTL0(false, "Buffer zlib error");
					}
					else
					{
						ASSERTL0(false, "Unknown zlib error");
					}
				}

				// Deserialize the array.
				double* readFieldData = (double*) elementDecompressedData.str().c_str();
				std::vector<double> elementFieldData(readFieldData, readFieldData + elementDecompressedData.str().length() * sizeof(*elementDecompressedData.str().c_str()) / sizeof(double));
				fielddata.push_back(elementFieldData);

				// Reconstruct the fielddefs.
				std::vector<unsigned int> elementIds;
				bool valid = ParseUtils::GenerateSeqVector(idString.c_str(), elementIds);
				ASSERTL0(valid, "Unable to correctly parse the element ids.");

				GeomShapeType elementType;
				valid = false;
				for (unsigned int i = 0; i < SIZE_GeomShapeType; i++)
				{
					if (GeomShapeTypeMap[i] == typeString)
					{
						elementType = (GeomShapeType) i;
						valid = true;
						break;
					}
				}
				ASSERTL0(valid, "Unable to correctly parse the element type.");

				std::vector<std::string> basisStrings;
				std::vector<LibUtilities::BasisType> basis;
				valid = ParseUtils::GenerateOrderedStringVector(basisString.c_str(), basisStrings);
				ASSERTL0(valid, "Unable to correctly parse the basis types.");
				for (std::vector<std::string>::size_type i = 0; i < basisStrings.size(); i++)
				{
					valid = false;
					for (unsigned int j = 0; j < LibUtilities::SIZE_BasisType; j++)
					{
						if (LibUtilities::BasisTypeMap[j] == basisStrings[i])
						{
							basis.push_back((LibUtilities::BasisType) j);
							valid = true;
							break;
						}
					}
					ASSERTL0(valid, std::string("Unable to correctly parse the basis type: ").append(basisStrings[i]).c_str());
				}

				std::vector<unsigned int> numModes;
				valid = ParseUtils::GenerateOrderedVector(numModesString.c_str(), numModes);
				ASSERTL0(valid, "Unable to correctly parse the number of modes.");

				std::vector<unsigned int> numFields;
				valid = ParseUtils::GenerateOrderedVector(fieldsString.c_str(), numFields);
				ASSERTL0(valid, "Unable to correctly parse the number of fields.");

				// Now the elements vector must be reconstructed from the element IDs.
				GeometryVector elements;
				CompositeVector domain = GetDomain();
				for (std::vector<unsigned int>::size_type elementIdIndex = 0; elementIdIndex < elementIds.size(); elementIdIndex++)
				{
					valid = false;
					for (CompositeVector::size_type compositeIndex = 0; compositeIndex < domain.size(); compositeIndex++)
					{
						for (GeometryVector::size_type geometryIndex = 0; geometryIndex < domain[compositeIndex]->size(); geometryIndex++)
						{
							GeometrySharedPtr geom = (*domain[compositeIndex])[geometryIndex];
							if (geom->GetGeomShapeType() == elementType &&
								geom->GetGlobalID() == elementIds[elementIdIndex])
							{
								valid = true;
								elements.push_back(geom);
							}
						}
						if (valid) break;
					}
					ASSERTL0(valid, "Unable to find element in the domain.");
				}

				fielddefs.push_back(FieldDefinitions(elements, basis, numModes, numFields));

				element = element->NextSiblingElement("element");
			}
			
			master = master->NextSiblingElement("nektar");
		}
	}

    MeshGraph::~MeshGraph()
    {
    }
    }; //end of namespace
}; //end of namespace

//
// $Log: MeshGraph.cpp,v $
// Revision 1.24  2008/08/18 20:54:02  ehan
// Changed name CURVE to CURVED.
//
// Revision 1.23  2008/07/14 21:04:26  ehan
// Added ASSERTL0 to check valid points type and number of points.
//
// Revision 1.22  2008/07/09 23:41:20  ehan
// Added edge component and face component to the curve reader.
//
// Revision 1.21  2008/07/08 18:58:07  ehan
// Added curve reader.
//
// Revision 1.20  2008/06/30 19:34:46  ehan
// Fixed infinity recursive-loop error.
//
// Revision 1.19  2008/06/11 16:10:12  delisi
// Added the 3D reader.
//
// Revision 1.18  2008/05/30 00:33:48  delisi
// Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
//
// Revision 1.17  2008/05/29 21:19:23  delisi
// Added the Write(...) and Import(...) functions which write and read XML files for output.
//
// Revision 1.16  2008/03/18 14:14:49  pvos
// Update for nodal triangular helmholtz solver
//
// Revision 1.15  2007/12/11 18:59:58  jfrazier
// Updated meshgraph so that a generic read could be performed and the proper type read (based on dimension) will be returned.
//
// Revision 1.14  2007/12/04 03:29:56  jfrazier
// Changed to stringstream.
//
// Revision 1.13  2007/11/07 20:31:03  jfrazier
// Added new expansion list to replace the expansion composite list.
//
// Revision 1.12  2007/09/25 04:45:14  jfrazier
// Added default linear expansions for the entire domain.
//
// Revision 1.11  2007/09/20 22:25:05  jfrazier
// Added expansion information to meshgraph class.
//
// Revision 1.10  2007/09/20 02:06:15  jfrazier
// General cleanup.
//
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
