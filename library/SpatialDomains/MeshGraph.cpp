////////////////////////////////////////////////////////////////////////////////
//
//  File: MeshGraph.cpp
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
////////////////////////////////////////////////////////////////////////////////


#include <SpatialDomains/MeshGraph.h>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/Equation.h>
#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdPyrExp.h>
#include <StdRegions/StdPrismExp.h>

// Use the stl version, primarily for string.
#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif

#include <tinyxml.h>
#include <cstring>
#include <sstream>
#include <iomanip>

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

using namespace std;

namespace Nektar
{
    namespace SpatialDomains
    {
        /**
         *
         */
        MeshGraph::MeshGraph():
            m_meshDimension(3),
            m_spaceDimension(3),
            m_domainRange(NullDomainRangeShPtr)
        {
        }


        /**
         *
         */
        MeshGraph::MeshGraph(
                unsigned int meshDimension,
                unsigned int spaceDimension) :
            m_meshDimension(meshDimension),
            m_spaceDimension(spaceDimension),
            m_domainRange(NullDomainRangeShPtr)
        {
        }


        /**
         *
         */
        MeshGraph::MeshGraph(
                       const LibUtilities::SessionReaderSharedPtr &pSession,
                       const DomainRangeShPtr &rng):
            m_session(pSession),
            m_domainRange(rng)
        {
        }



        /**
         *
         */
        MeshGraph::~MeshGraph()
        {
        }


        /**
         *
         */
        boost::shared_ptr<MeshGraph> MeshGraph::Read(
                      const LibUtilities::SessionReaderSharedPtr &pSession,
                      DomainRangeShPtr &rng)
        {
            boost::shared_ptr<MeshGraph> returnval;

            // read the geometry tag to get the dimension

            TiXmlElement* geometry_tag = pSession->GetElement("NEKTAR/GEOMETRY");
            TiXmlAttribute *attr = geometry_tag->FirstAttribute();
            int meshDim = 0;
            while (attr)
            {
                std::string attrName(attr->Name());
                if (attrName == "DIM")
                {
                    int err = attr->QueryIntValue(&meshDim);
                    ASSERTL0(err==TIXML_SUCCESS, "Unable to read mesh dimension.");
                    break;
                }
                else
                {
                    std::string errstr("Unknown attribute: ");
                    errstr += attrName;
                    ASSERTL0(false, errstr.c_str());
                }

                // Get the next attribute.
                attr = attr->Next();
            }

            // instantiate the dimension-specific meshgraph classes

            switch(meshDim)
            {
            case 1:
                returnval = MemoryManager<MeshGraph1D>::AllocateSharedPtr(pSession,rng);
                break;

            case 2:
                returnval = MemoryManager<MeshGraph2D>::AllocateSharedPtr(pSession,rng);
                break;

            case 3:
                returnval = MemoryManager<MeshGraph3D>::AllocateSharedPtr(pSession,rng);
                break;

            default:
                std::string err = "Invalid mesh dimension: ";
                std::stringstream strstrm;
                strstrm << meshDim;
                err += strstrm.str();
                NEKERROR(ErrorUtil::efatal, err.c_str());
            }

            return returnval;
        }



        /*  ====  OUTDATED ROUTINE, PLEASE NOT USE  ==== */
        boost::shared_ptr<MeshGraph> MeshGraph::Read(
                const std::string& infilename,
                bool pReadExpansions)
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
                returnval->ReadGeometryInfo(infilename);
                if (pReadExpansions)
                {
                    returnval->ReadExpansions(infilename);
                }
            }
            return returnval;
        }



        /**
         *
         */
        void MeshGraph::ReadGeometry(const std::string& infilename)
        {
            TiXmlDocument doc(infilename);
            bool loadOkay = doc.LoadFile();

            std::stringstream errstr;
            errstr << "Unable to load file: " << infilename << " (";
            errstr << doc.ErrorDesc() << ", line " << doc.ErrorRow()
                                 << ", column " << doc.ErrorCol() << ")";
            ASSERTL0(loadOkay, errstr.str());

            ReadGeometry(doc);
        }


        /**
         *
         */
        void MeshGraph::ReadGeometry(TiXmlDocument &doc)
        {
            TiXmlHandle docHandle(&doc);
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
            m_meshPartitioned = false;
            m_meshDimension = 3;
            m_spaceDimension = 3;

            while (attr)
            {
                std::string attrName(attr->Name());
                if (attrName == "DIM")
                {
                    err = attr->QueryIntValue(&m_meshDimension);
                    ASSERTL1(err==TIXML_SUCCESS, "Unable to read mesh dimension.");
                }
                else if (attrName == "SPACE")
                {
                    err = attr->QueryIntValue(&m_spaceDimension);
                    ASSERTL1(err==TIXML_SUCCESS, "Unable to read space dimension.");
                }
                else if (attrName == "PARTITION")
                {
                    err = attr->QueryIntValue(&m_partition);
                    ASSERTL1(err==TIXML_SUCCESS, "Unable to read partition.");
                    m_meshPartitioned = true;
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

            ASSERTL1(m_meshDimension<=m_spaceDimension, "Mesh dimension greater than space dimension");

            // Now read the vertices
            TiXmlElement* element = mesh->FirstChildElement("VERTEX");
            ASSERTL0(element, "Unable to find mesh VERTEX tag in file.");

            NekDouble xscale,yscale,zscale;

            // check to see if any scaling parameters are in
            // attributes and determine these values
            LibUtilities::AnalyticExpressionEvaluator expEvaluator;
            //LibUtilities::ExpressionEvaluator expEvaluator;
            const char *xscal =  element->Attribute("XSCALE");
            if(!xscal)
            {
                xscale = 1.0;
            }
            else
            {
                std::string xscalstr = xscal;
                int expr_id = expEvaluator.DefineFunction("",xscalstr);
                xscale = expEvaluator.Evaluate(expr_id);
            }

            const char *yscal =  element->Attribute("YSCALE");
            if(!yscal)
            {
                yscale = 1.0;
            }
            else
            {
                std::string yscalstr = yscal;
                int expr_id = expEvaluator.DefineFunction("",yscalstr);
                yscale = expEvaluator.Evaluate(expr_id);
            }

            const char *zscal = element->Attribute("ZSCALE");
            if(!zscal)
            {
                zscale = 1.0;
            }
            else
            {
                std::string zscalstr = zscal;
                int expr_id = expEvaluator.DefineFunction("",zscalstr);
                zscale = expEvaluator.Evaluate(expr_id);
            }


            NekDouble xmove,ymove,zmove;

            // check to see if any moving parameters are in
            // attributes and determine these values

            //LibUtilities::ExpressionEvaluator expEvaluator;
            const char *xmov =  element->Attribute("XMOVE");
            if(!xmov)
            {
                xmove = 0.0;
            }
            else
            {
                std::string xmovstr = xmov;
                int expr_id = expEvaluator.DefineFunction("",xmovstr);
                xmove = expEvaluator.Evaluate(expr_id);
            }

            const char *ymov =  element->Attribute("YMOVE");
            if(!ymov)
            {
                ymove = 0.0;
            }
            else
            {
                std::string ymovstr = ymov;
                int expr_id = expEvaluator.DefineFunction("",ymovstr);
                ymove = expEvaluator.Evaluate(expr_id);
            }

            const char *zmov = element->Attribute("ZMOVE");
            if(!zmov)
            {
                zmove = 0.0;
            }
            else
            {
                std::string zmovstr = zmov;
                int expr_id = expEvaluator.DefineFunction("",zmovstr);
                zmove = expEvaluator.Evaluate(expr_id);
            }

            string IsCompressed;
            element->QueryStringAttribute("COMPRESSED",&IsCompressed); 

            if(IsCompressed.size()) 
            {
                if(boost::iequals(IsCompressed,
                            LibUtilities::CompressData::GetCompressString()))
                {
                    // Extract the vertex body
                    TiXmlNode* vertexChild = element->FirstChild();
                    ASSERTL0(vertexChild,
                             "Unable to extract the data from the compressed "
                             "vertex tag.");

                    std::string vertexStr;
                    if (vertexChild->Type() == TiXmlNode::TINYXML_TEXT)
                    {
                        vertexStr += vertexChild->ToText()->ValueStr();
                    }

                    std::vector<LibUtilities::MeshVertex> vertData;
                    LibUtilities::CompressData::ZlibDecodeFromBase64Str(
                                                        vertexStr,vertData);

                    int indx;
                    NekDouble xval, yval, zval;
                    for(int i = 0; i < vertData.size(); ++i)
                    {
                        indx = vertData[i].id;
                        xval = vertData[i].x;
                        yval = vertData[i].y;
                        zval = vertData[i].z;

                        xval = xval*xscale + xmove;
                        yval = yval*yscale + ymove;
                        zval = zval*zscale + zmove;

                        PointGeomSharedPtr vert(
                            MemoryManager<PointGeom>::AllocateSharedPtr(
                                m_spaceDimension, indx, xval, yval, zval));

                        vert->SetGlobalID(indx);
                        m_vertSet[indx] = vert;
                    }
                }
                else
                {
                    ASSERTL0(false,"Compressed formats do not match. Expected :"
                             + LibUtilities::CompressData::GetCompressString()
                             + " but got " + std::string(IsCompressed));
                }
            }
            else
            {
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

                    // Now read body of vertex
                    std::string vertexBodyStr;
                    
                    TiXmlNode *vertexBody = vertex->FirstChild();
                    
                    while (vertexBody)
                    {
                        // Accumulate all non-comment body data.
                        if (vertexBody->Type() == TiXmlNode::TINYXML_TEXT)
                        {
                            vertexBodyStr += vertexBody->ToText()->Value();
                            vertexBodyStr += " ";
                        }
                        
                        vertexBody = vertexBody->NextSibling();
                    }
                    
                    ASSERTL0(!vertexBodyStr.empty(), "Vertex definitions must contain vertex data.");
                    
                    // Get vertex data from the data string.
                    NekDouble xval, yval, zval;
                    std::istringstream vertexDataStrm(vertexBodyStr.c_str());
                    
                    try
                    {
                        while(!vertexDataStrm.fail())
                        {
                            vertexDataStrm >> xval >> yval >> zval;
                            
                            xval = xval*xscale + xmove;
                            yval = yval*yscale + ymove;
                            zval = zval*zscale + zmove;
                            
                            // Need to check it here because we may not be
                            // good after the read indicating that there
                            // was nothing to read.
                            if (!vertexDataStrm.fail())
                            {
                                PointGeomSharedPtr vert(MemoryManager<PointGeom>::AllocateSharedPtr(m_spaceDimension, indx, xval, yval, zval));
                                vert->SetGlobalID(indx);
                                m_vertSet[indx] = vert;
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
        }


        /**
         * Read the geometry-related information from the given file. This
         * information is located within the XML tree under
         * <NEKTAR><GEOMETRY><GEOMINFO>.
         * @param   infilename      Filename of XML file.
         */
        void MeshGraph::ReadGeometryInfo(const std::string &infilename)
        {
            TiXmlDocument doc(infilename);
            bool loadOkay = doc.LoadFile();

            std::stringstream errstr;
            errstr << "Unable to load file: " << infilename << std::endl;
            errstr << "Reason: " << doc.ErrorDesc() << std::endl;
            errstr << "Position: Line " << doc.ErrorRow() << ", Column " << doc.ErrorCol() << std::endl;
            ASSERTL0(loadOkay, errstr.str());

            ReadGeometryInfo(doc);
        }


        /**
         * Read the geometry-related information from the given XML document.
         * This information is located within the XML tree under
         * <NEKTAR><GEOMETRY><GEOMINFO>.
         * @param   doc             XML document.
         */
        void MeshGraph::ReadGeometryInfo(TiXmlDocument &doc)
        {
            TiXmlElement *master = doc.FirstChildElement("NEKTAR");
            ASSERTL0(master, "Unable to find NEKTAR tag in file.");

            // Find the Expansions tag
            TiXmlElement *geomTag = master->FirstChildElement("GEOMETRY");
            ASSERTL0(geomTag, "Unable to find GEOMETRY tag in file.");

            // See if we have GEOMINFO. If there is none, it's fine.
            TiXmlElement *geomInfoTag = geomTag->FirstChildElement("GEOMINFO");
            if (!geomInfoTag) return;

            TiXmlElement *infoItem = geomInfoTag->FirstChildElement("I");

            // Multiple nodes will only occur if there is a comment in between
            // definitions.
            while (infoItem)
            {
                std::string geomProperty = infoItem->Attribute("PROPERTY");
                std::string geomValue    = infoItem->Attribute("VALUE");
                GeomInfoMap::iterator x  = m_geomInfo.find(geomProperty);

                ASSERTL0(x == m_geomInfo.end(),
                        "Property " + geomProperty + " already specified.");
                m_geomInfo[geomProperty] = geomValue;
                infoItem = infoItem->NextSiblingElement("I");
            }
        }


        /**
         *
         */
        void MeshGraph::ReadExpansions(const std::string& infilename)
        {
            TiXmlDocument doc(infilename);
            bool loadOkay = doc.LoadFile();

            std::stringstream errstr;
            errstr << "Unable to load file: " << infilename << std::endl;
            errstr << "Reason: " << doc.ErrorDesc() << std::endl;
            errstr << "Position: Line " << doc.ErrorRow() << ", Column " << doc.ErrorCol() << std::endl;
            ASSERTL0(loadOkay, errstr.str());

            ReadExpansions(doc);
        }


        /**
         *
         */
        void MeshGraph::ReadExpansions(TiXmlDocument &doc)
        {
            TiXmlElement *master = doc.FirstChildElement("NEKTAR");
            ASSERTL0(master, "Unable to find NEKTAR tag in file.");

            // Find the Expansions tag
            TiXmlElement *expansionTypes = master->FirstChildElement("EXPANSIONS");
            ASSERTL0(expansionTypes, "Unable to find EXPANSIONS tag in file.");

            if(expansionTypes)
            {
                // Find the Expansion type
                TiXmlElement *expansion = expansionTypes->FirstChildElement();
                std::string   expType   = expansion->Value();

                if(expType == "E")
                {
                    int i;
                    ExpansionMapShPtr expansionMap;

                    /// Expansiontypes will contain composite,
                    /// nummodes, and expansiontype (eModified, or
                    /// eOrthogonal) Or a full list of data of
                    /// basistype, nummodes, pointstype, numpoints;

                    /// Expansiontypes may also contain a list of
                    /// fields that this expansion relates to. If this
                    /// does not exist the variable is only set to
                    /// "DefaultVar".

                    while (expansion)
                    {

                        const char *fStr = expansion->Attribute("FIELDS");
                        std::vector<std::string> fieldStrings;

                        if(fStr) // extract other fields.
                        {
                            std::string fieldStr = fStr;
                            bool  valid = ParseUtils::GenerateOrderedStringVector(fieldStr.c_str(),fieldStrings);
                            ASSERTL0(valid,"Unable to correctly parse the field string in ExpansionTypes.");
                        }

                        // check to see if m_expasionVectorShPtrMap has
                        // already been intiailised and if not intiailse
                        // vector.
                        if(m_expansionMapShPtrMap.count("DefaultVar") == 0) // no previous definitions
                        {
                            expansionMap = SetUpExpansionMap();

                            m_expansionMapShPtrMap["DefaultVar"] = expansionMap;

                            // make sure all fields in this search point
                            // to same expansion vector;
                            for(i = 0; i < fieldStrings.size(); ++i)
                            {
                                m_expansionMapShPtrMap[fieldStrings[i]] = expansionMap;
                            }
                        }
                        else // default variable is defined
                        {

                            if(fieldStrings.size()) // fields are defined
                            {
                                //see if field exists
                                if(m_expansionMapShPtrMap.count(fieldStrings[0]))
                                {
                                    expansionMap = m_expansionMapShPtrMap.find(fieldStrings[0])->second;
                                }
                                else
                                {
                                    expansionMap = SetUpExpansionMap();
                                    // make sure all fields in this search point
                                    // to same expansion vector;
                                    for(i = 0; i < fieldStrings.size(); ++i)
                                    {
                                        if(m_expansionMapShPtrMap.count(fieldStrings[i]) == 0)
                                        {
                                            m_expansionMapShPtrMap[fieldStrings[i]] = expansionMap;
                                        }
                                        else
                                        {
                                            ASSERTL0(false,"Expansion vector for this field is already  setup");
                                        }
                                    }
                                }
                            }
                            else // use default variable list
                            {
                                expansionMap = m_expansionMapShPtrMap.find("DefaultVar")->second;
                            }

                        }

                        /// Mandatory components...optional are to follow later.
                        std::string compositeStr = expansion->Attribute("COMPOSITE");
                        ASSERTL0(compositeStr.length() > 3, "COMPOSITE must be specified in expansion definition");
                        int beg = compositeStr.find_first_of("[");
                        int end = compositeStr.find_first_of("]");
                        std::string compositeListStr = compositeStr.substr(beg+1,end-beg-1);

                        CompositeMap compositeVector;
                        GetCompositeList(compositeListStr, compositeVector);

                        bool          useExpansionType = false;
                        ExpansionType expansion_type;
                        int           num_modes;

                        LibUtilities::BasisKeyVector basiskeyvec;
                        const char * tStr = expansion->Attribute("TYPE");

                        if(tStr) // use type string to define expansion
                        {
                            std::string typeStr = tStr;
                            const std::string* begStr = kExpansionTypeStr;
                            const std::string* endStr = kExpansionTypeStr+eExpansionTypeSize;
                            const std::string* expStr = std::find(begStr, endStr, typeStr);

                            ASSERTL0(expStr != endStr, "Invalid expansion type.");
                            expansion_type = (ExpansionType)(expStr - begStr);


                            /// \todo solvers break the pattern 'instantiate Session -> instantiate MeshGraph'
                            /// and parse command line arguments by themselves; one needs to unify command
                            /// line arguments handling.
                            /// Solvers tend to call MeshGraph::Read statically -> m_session
                            /// is not defined -> no info about command line arguments presented
                            /// ASSERTL0(m_session != 0, "One needs to instantiate SessionReader first");

                            const char *nStr = expansion->Attribute("NUMMODES");
                            ASSERTL0(nStr,"NUMMODES was not defined in EXPANSION section of input");
                            std::string nummodesStr = nStr;

                            // ASSERTL0(m_session,"Session should be defined to evaluate nummodes ");
                            if (m_session)
                            {
                                LibUtilities::Equation nummodesEqn(m_session, nummodesStr);
                                num_modes = (int) nummodesEqn.Evaluate();
                            }
                            else
                            {
                                num_modes = boost::lexical_cast<int>(nummodesStr);
                            }

                            useExpansionType = true;
                        }
                        else // assume expansion is defined individually
                        {
                            // Extract the attributes.
                            const char *bTypeStr = expansion->Attribute("BASISTYPE");
                            ASSERTL0(bTypeStr,"TYPE or BASISTYPE was not defined in EXPANSION section of input");
                            std::string basisTypeStr = bTypeStr;

                            // interpret the basis type string.
                            std::vector<std::string> basisStrings;
                            std::vector<LibUtilities::BasisType> basis;
                            bool valid = ParseUtils::GenerateOrderedStringVector(basisTypeStr.c_str(), basisStrings);
                            ASSERTL0(valid, "Unable to correctly parse the basis types.");
                            for (vector<std::string>::size_type i = 0; i < basisStrings.size(); i++)
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
                            const char *nModesStr = expansion->Attribute("NUMMODES");
                            ASSERTL0(nModesStr,"NUMMODES was not defined in EXPANSION section of input");

                            std::string numModesStr = nModesStr;
                            std::vector<unsigned int> numModes;
                            valid = ParseUtils::GenerateOrderedVector(numModesStr.c_str(), numModes);
                            ASSERTL0(valid, "Unable to correctly parse the number of modes.");
                            ASSERTL0(numModes.size() == basis.size(),"information for num modes does not match the number of basis");

                            const char *pTypeStr =  expansion->Attribute("POINTSTYPE");
                            ASSERTL0(pTypeStr,"POINTSTYPE was not defined in EXPANSION section of input");
                            std::string pointsTypeStr = pTypeStr;
                            // interpret the points type string.
                            std::vector<std::string> pointsStrings;
                            std::vector<LibUtilities::PointsType> points;
                            valid = ParseUtils::GenerateOrderedStringVector(pointsTypeStr.c_str(), pointsStrings);
                            ASSERTL0(valid, "Unable to correctly parse the points types.");
                            for (vector<std::string>::size_type i = 0; i < pointsStrings.size(); i++)
                            {
                                valid = false;
                                for (unsigned int j = 0; j < LibUtilities::SIZE_PointsType; j++)
                                {
                                    if (LibUtilities::kPointsTypeStr[j] == pointsStrings[i])
                                    {
                                        points.push_back((LibUtilities::PointsType) j);
                                        valid = true;
                                        break;
                                    }
                                }
                                ASSERTL0(valid, std::string("Unable to correctly parse the points type: ").append(pointsStrings[i]).c_str());
                            }

                            const char *nPointsStr = expansion->Attribute("NUMPOINTS");
                            ASSERTL0(nPointsStr,"NUMPOINTS was not defined in EXPANSION section of input");
                            std::string numPointsStr = nPointsStr;
                            std::vector<unsigned int> numPoints;
                            valid = ParseUtils::GenerateOrderedVector(numPointsStr.c_str(), numPoints);
                            ASSERTL0(valid, "Unable to correctly parse the number of points.");
                            ASSERTL0(numPoints.size() == numPoints.size(),"information for num points does not match the number of basis");

                            for(int i = 0; i < basis.size(); ++i)
                            {
                                //Generate Basis key  using information
                                const LibUtilities::PointsKey pkey(numPoints[i],points[i]);
                                basiskeyvec.push_back(LibUtilities::BasisKey(basis[i],numModes[i],pkey));
                            }
                        }

                        // Now have composite and basiskeys.  Cycle through
                        // all composites for the geomShPtrs and set the modes
                        // and types for the elements contained in the element
                        // list.
                        CompositeMapIter compVecIter;
                        for (compVecIter = compositeVector.begin(); compVecIter != compositeVector.end(); ++compVecIter)
                        {
                            GeometryVectorIter geomVecIter;
                            for (geomVecIter = (compVecIter->second)->begin(); geomVecIter != (compVecIter->second)->end(); ++geomVecIter)
                            {
                                ExpansionMapIter x = expansionMap->find((*geomVecIter)->GetGlobalID());
                                ASSERTL0(x != expansionMap->end(), "Expansion not found!!");
                                if(useExpansionType)
                                {
                                    (x->second)->m_basisKeyVector = MeshGraph::DefineBasisKeyFromExpansionType(*geomVecIter,expansion_type,num_modes);
                                }
                                else
                                {
                                    ASSERTL0((*geomVecIter)->GetShapeDim() == basiskeyvec.size()," There is an incompatible expansion dimension with geometry dimension");
                                    (x->second)->m_basisKeyVector = basiskeyvec;
                                }
                            }
                        }

                        expansion = expansion->NextSiblingElement("E");
                    }
                }
                else if(expType == "H")
                {
                    int i;
                    ExpansionMapShPtr expansionMap;

                    while (expansion)
                    {

                        const char *fStr = expansion->Attribute("FIELDS");
                        std::vector<std::string> fieldStrings;

                        if(fStr) // extract other fields.
                        {
                            std::string fieldStr = fStr;
                            bool  valid = ParseUtils::GenerateOrderedStringVector(fieldStr.c_str(),fieldStrings);
                            ASSERTL0(valid,"Unable to correctly parse the field string in ExpansionTypes.");
                        }

                        // check to see if m_expasionVectorShPtrMap has
                        // already been intiailised and if not intiailse
                        // vector.
                        if(m_expansionMapShPtrMap.count("DefaultVar") == 0) // no previous definitions
                        {
                            expansionMap = SetUpExpansionMap();

                            m_expansionMapShPtrMap["DefaultVar"] = expansionMap;

                            // make sure all fields in this search point
                            // to same expansion vector;
                            for(i = 0; i < fieldStrings.size(); ++i)
                            {
                                m_expansionMapShPtrMap[fieldStrings[i]] = expansionMap;
                            }
                        }
                        else // default variable is defined
                        {

                            if(fieldStrings.size()) // fields are defined
                            {
                                //see if field exists
                                if(m_expansionMapShPtrMap.count(fieldStrings[0]))
                                {
                                    expansionMap = m_expansionMapShPtrMap.find(fieldStrings[0])->second;
                                }
                                else
                                {
                                    expansionMap = SetUpExpansionMap();
                                    // make sure all fields in this search point
                                    // to same expansion vector;
                                    for(i = 0; i < fieldStrings.size(); ++i)
                                    {
                                        if(m_expansionMapShPtrMap.count(fieldStrings[i]) == 0)
                                        {
                                            m_expansionMapShPtrMap[fieldStrings[i]] = expansionMap;
                                        }
                                        else
                                        {
                                            ASSERTL0(false,"Expansion vector for this field is already  setup");
                                        }
                                    }
                                }
                            }
                            else // use default variable list
                            {
                                expansionMap = m_expansionMapShPtrMap.find("DefaultVar")->second;
                            }

                        }

                        /// Mandatory components...optional are to follow later.
                        std::string compositeStr = expansion->Attribute("COMPOSITE");
                        ASSERTL0(compositeStr.length() > 3, "COMPOSITE must be specified in expansion definition");
                        int beg = compositeStr.find_first_of("[");
                        int end = compositeStr.find_first_of("]");
                        std::string compositeListStr = compositeStr.substr(beg+1,end-beg-1);

                        CompositeMap compositeVector;
                        GetCompositeList(compositeListStr, compositeVector);

                        ExpansionType expansion_type_x = eNoExpansionType;
                        ExpansionType expansion_type_y = eNoExpansionType;
                        ExpansionType expansion_type_z = eNoExpansionType;
                        int           num_modes_x = 0;
                        int           num_modes_y = 0;
                        int           num_modes_z = 0;

                        LibUtilities::BasisKeyVector basiskeyvec;

                        const char * tStr_x = expansion->Attribute("TYPE-X");

                        if(tStr_x) // use type string to define expansion
                        {
                            std::string typeStr = tStr_x;
                            const std::string* begStr = kExpansionTypeStr;
                            const std::string* endStr = kExpansionTypeStr+eExpansionTypeSize;
                            const std::string* expStr = std::find(begStr, endStr, typeStr);

                            ASSERTL0(expStr != endStr, "Invalid expansion type.");
                            expansion_type_x = (ExpansionType)(expStr - begStr);

                            const char *nStr = expansion->Attribute("NUMMODES-X");
                            ASSERTL0(nStr,"NUMMODES-X was not defined in EXPANSION section of input");
                            std::string nummodesStr = nStr;

                            // ASSERTL0(m_session,"Session should be defined to evaluate nummodes ");

                            if (m_session)
                            {
                                LibUtilities::Equation nummodesEqn(m_session, nummodesStr);
                                num_modes_x = (int) nummodesEqn.Evaluate();
                            }
                            else
                            {
                                num_modes_x = boost::lexical_cast<int>(nummodesStr);
                            }

                        }

                        const char * tStr_y = expansion->Attribute("TYPE-Y");

                        if(tStr_y) // use type string to define expansion
                        {
                            std::string typeStr = tStr_y;
                            const std::string* begStr = kExpansionTypeStr;
                            const std::string* endStr = kExpansionTypeStr+eExpansionTypeSize;
                            const std::string* expStr = std::find(begStr, endStr, typeStr);

                            ASSERTL0(expStr != endStr, "Invalid expansion type.");
                            expansion_type_y = (ExpansionType)(expStr - begStr);

                            const char *nStr = expansion->Attribute("NUMMODES-Y");
                            ASSERTL0(nStr,"NUMMODES-Y was not defined in EXPANSION section of input");
                            std::string nummodesStr = nStr;

                            // ASSERTL0(m_session,"Session should be defined to evaluate nummodes ");
                            if (m_session)
                            {
                                LibUtilities::Equation nummodesEqn(m_session, nummodesStr);
                                num_modes_y = (int) nummodesEqn.Evaluate();
                            }
                            else
                            {
                                num_modes_y = boost::lexical_cast<int>(nummodesStr);
                            }

                        }

                        const char * tStr_z = expansion->Attribute("TYPE-Z");

                        if(tStr_z) // use type string to define expansion
                        {
                            std::string typeStr = tStr_z;
                            const std::string* begStr = kExpansionTypeStr;
                            const std::string* endStr = kExpansionTypeStr+eExpansionTypeSize;
                            const std::string* expStr = std::find(begStr, endStr, typeStr);

                            ASSERTL0(expStr != endStr, "Invalid expansion type.");
                            expansion_type_z = (ExpansionType)(expStr - begStr);

                            const char *nStr = expansion->Attribute("NUMMODES-Z");
                            ASSERTL0(nStr,"NUMMODES-Z was not defined in EXPANSION section of input");
                            std::string nummodesStr = nStr;

                            // ASSERTL0(m_session,"Session should be defined to evaluate nummodes ");
                            if (m_session)
                            {
                                LibUtilities::Equation nummodesEqn(m_session, nummodesStr);
                                num_modes_z = (int) nummodesEqn.Evaluate();
                            }
                            else
                            {
                                num_modes_z = boost::lexical_cast<int>(nummodesStr);
                            }

                        }

                        CompositeMapIter compVecIter;
                        for (compVecIter = compositeVector.begin(); compVecIter != compositeVector.end(); ++compVecIter)
                        {
                            GeometryVectorIter geomVecIter;
                            for (geomVecIter = (compVecIter->second)->begin(); geomVecIter != (compVecIter->second)->end(); ++geomVecIter)
                            {
                                ExpansionMapIter expVecIter;
                                for (expVecIter = expansionMap->begin(); expVecIter != expansionMap->end(); ++expVecIter)
                                {

                                    (expVecIter->second)->m_basisKeyVector = DefineBasisKeyFromExpansionTypeHomo(*geomVecIter,
                                            expansion_type_x,
                                            expansion_type_y,
                                            expansion_type_z,
                                            num_modes_x,
                                            num_modes_y,
                                            num_modes_z);
                                }
                            }
                        }

                        expansion = expansion->NextSiblingElement("H");
                    }
                }
                else if(expType == "ELEMENTS")  // Reading a file with the expansion definition
                {
                    std::vector<LibUtilities::FieldDefinitionsSharedPtr> fielddefs;
                    LibUtilities::FieldIO f(m_session->GetComm());
                    f.ImportFieldDefs(doc, fielddefs, true);
                    cout << "    Number of elements: " << fielddefs.size() << endl;
                    SetExpansions(fielddefs);
                }
                else
                {
                    ASSERTL0(false,"Expansion type not defined");
                }
            }
        }


        /**
         *
         */
        void MeshGraph::ReadDomain(TiXmlDocument &doc)
        {
            TiXmlHandle docHandle(&doc);

            TiXmlElement* mesh = docHandle.FirstChildElement("NEKTAR").FirstChildElement("GEOMETRY").Element();
            TiXmlElement* domain = NULL;

            ASSERTL0(mesh, "Unable to find GEOMETRY tag in file.");

            /// Look for data in DOMAIN block.
            domain = mesh->FirstChildElement("DOMAIN");

            ASSERTL0(domain, "Unable to find DOMAIN tag in file.");

            /// Elements are of the form: "<D ID = "N"> ... </D>".
            /// Read the ID field first.
            TiXmlElement *multidomains = domain->FirstChildElement("D");

            if(multidomains)
            {
                int nextDomainNumber = 0;
                while (multidomains)
                {
                    int indx;
                    int err = multidomains->QueryIntAttribute("ID", &indx);
                    ASSERTL0(err == TIXML_SUCCESS,
                             "Unable to read attribute ID in Domain.");


                    TiXmlNode* elementChild = multidomains->FirstChild();
                    while(elementChild && elementChild->Type() != TiXmlNode::TINYXML_TEXT)
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
                    CompositeMap unrollDomain;
                    GetCompositeList(indxStr, unrollDomain);
                    m_domain.push_back(unrollDomain);

                    ASSERTL0(!m_domain[nextDomainNumber++].empty(), (std::string("Unable to obtain domain's referenced composite: ") + indxStr).c_str());

                    /// Keep looking
                    multidomains = multidomains->NextSiblingElement("D");
                }

            }
            else // previous definition of just one composite
            {

                // find the non comment portion of the body.
                TiXmlNode* elementChild = domain->FirstChild();
                while(elementChild && elementChild->Type() != TiXmlNode::TINYXML_TEXT)
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
                CompositeMap fullDomain;
                GetCompositeList(indxStr, fullDomain);
                m_domain.push_back(fullDomain);

                ASSERTL0(!m_domain[0].empty(), (std::string("Unable to obtain domain's referenced composite: ") + indxStr).c_str());
            }
        }


        /**
         *
         */
        void MeshGraph::ReadCurves(TiXmlDocument &doc)
        {
            /// We know we have it since we made it this far.
            TiXmlHandle docHandle(&doc);
            TiXmlElement* mesh = docHandle.FirstChildElement("NEKTAR").FirstChildElement("GEOMETRY").Element();
            TiXmlElement* field = NULL;

            // check to see if any scaling parameters are in
            // attributes and determine these values
            TiXmlElement* element = mesh->FirstChildElement("VERTEX");
            ASSERTL0(element, "Unable to find mesh VERTEX tag in file.");

            NekDouble xscale,yscale,zscale;

            LibUtilities::AnalyticExpressionEvaluator expEvaluator;
            const char *xscal =  element->Attribute("XSCALE");
            if(!xscal)
            {
                xscale = 1.0;
            }
            else
            {
                std::string xscalstr = xscal;
                int expr_id = expEvaluator.DefineFunction("",xscalstr);
                xscale = expEvaluator.Evaluate(expr_id);
            }

            const char *yscal =  element->Attribute("YSCALE");
            if(!yscal)
            {
                yscale = 1.0;
            }
            else
            {
                std::string yscalstr = yscal;
                int expr_id = expEvaluator.DefineFunction("",yscalstr);
                yscale = expEvaluator.Evaluate(expr_id);
            }

            const char *zscal = element->Attribute("ZSCALE");
            if(!zscal)
            {
                zscale = 1.0;
            }
            else
            {
                std::string zscalstr = zscal;
                int expr_id = expEvaluator.DefineFunction("",zscalstr);
                zscale = expEvaluator.Evaluate(expr_id);
            }

            NekDouble xmove,ymove,zmove;

            // check to see if any moving parameters are in
            // attributes and determine these values

            //LibUtilities::ExpressionEvaluator expEvaluator;
            const char *xmov =  element->Attribute("XMOVE");
            if(!xmov)
            {
                xmove = 0.0;
            }
            else
            {
                std::string xmovstr = xmov;
                int expr_id = expEvaluator.DefineFunction("",xmovstr);
                xmove = expEvaluator.Evaluate(expr_id);
            }

            const char *ymov =  element->Attribute("YMOVE");
            if(!ymov)
            {
                ymove = 0.0;
            }
            else
            {
                std::string ymovstr = ymov;
                int expr_id = expEvaluator.DefineFunction("",ymovstr);
                ymove = expEvaluator.Evaluate(expr_id);
            }

            const char *zmov = element->Attribute("ZMOVE");
            if(!zmov)
            {
                zmove = 0.0;
            }
            else
            {
                std::string zmovstr = zmov;
                int expr_id = expEvaluator.DefineFunction("",zmovstr);
                zmove = expEvaluator.Evaluate(expr_id);
            }

            int err;

            /// Look for elements in CURVE block.
            field = mesh->FirstChildElement("CURVED");

            if(!field) //return if no curved entities
            {
                return;
            }

            string IsCompressed;
            field->QueryStringAttribute("COMPRESSED",&IsCompressed); 
            
            if(IsCompressed.size()) 
            {
                ASSERTL0(boost::iequals(IsCompressed,
                            LibUtilities::CompressData::GetCompressString()),
                         "Compressed formats do not match. Expected :"
                         + LibUtilities::CompressData::GetCompressString()
                         + " but got "
                         + boost::lexical_cast<std::string>(IsCompressed));

                std::vector<LibUtilities::MeshCurvedInfo> edginfo;
                std::vector<LibUtilities::MeshCurvedInfo> facinfo;
                LibUtilities::MeshCurvedPts cpts;

                // read edge, face info and curved poitns.
                TiXmlElement *x = field->FirstChildElement();
                while(x)
                {
                    const char *entitytype = x->Value();
                    // read in edge or face info
                    if(boost::iequals(entitytype,"E"))
                    {
                        // read in data
                        std::string elmtStr;
                        TiXmlNode* child = x->FirstChild();

                        if (child->Type() == TiXmlNode::TINYXML_TEXT)
                        {
                            elmtStr += child->ToText()->ValueStr();
                        }

                        LibUtilities::CompressData::ZlibDecodeFromBase64Str(
                                                        elmtStr,edginfo);
                    }
                    else if(boost::iequals(entitytype,"F"))
                    {
                        // read in data
                        std::string elmtStr;
                        TiXmlNode* child = x->FirstChild();

                        if (child->Type() == TiXmlNode::TINYXML_TEXT)
                        {
                            elmtStr += child->ToText()->ValueStr();
                        }

                        LibUtilities::CompressData::ZlibDecodeFromBase64Str(
                                                        elmtStr,facinfo);
                    }
                    else if(boost::iequals(entitytype,"DATAPOINTS"))
                    {
                        NekInt id;
                        ASSERTL0(x->Attribute("ID", &id),
                                 "Failed to get ID from PTS section");
                        cpts.id = id;

                        // read in data
                        std::string elmtStr;

                        TiXmlElement* DataIdx =
                            x->FirstChildElement("INDEX");
                        ASSERTL0(DataIdx,
                                 "Cannot read data index tag in compressed "
                                 "curved section");

                        TiXmlNode* child = DataIdx->FirstChild();
                        if (child->Type() == TiXmlNode::TINYXML_TEXT)
                        {
                            elmtStr = child->ToText()->ValueStr();
                        }

                        LibUtilities::CompressData::ZlibDecodeFromBase64Str(
                                                    elmtStr,cpts.index);

                        TiXmlElement* DataPts =
                            x->FirstChildElement("POINTS");
                        ASSERTL0(DataPts,
                                 "Cannot read data pts tag in compressed "
                                 "curved section");

                        child = DataPts->FirstChild();
                        if (child->Type() == TiXmlNode::TINYXML_TEXT)
                        {
                            elmtStr = child->ToText()->ValueStr();
                        }

                        LibUtilities::CompressData::ZlibDecodeFromBase64Str(
                                                    elmtStr,cpts.pts);
                    }
                    else
                    {
                        ASSERTL0(false,"Unknown tag in curved section");
                    }
                    x = x->NextSiblingElement();
                }

                // rescale (x,y,z) points;
                for(int i = 0; i > cpts.pts.size(); ++i)
                {
                    cpts.pts[i].x = xscale*cpts.pts[i].x + xmove;
                    cpts.pts[i].y = yscale*cpts.pts[i].y + ymove;
                    cpts.pts[i].z = zscale*cpts.pts[i].z + zmove;
                }

                for(int i = 0; i < edginfo.size(); ++i)
                {
                    int edgeid = edginfo[i].entityid;
                    LibUtilities::PointsType ptype;

                    CurveSharedPtr curve(
                            MemoryManager<Curve>::AllocateSharedPtr(
                                edgeid, ptype = (LibUtilities::PointsType)
                                                    edginfo[i].ptype));

                    // load points
                    int offset = edginfo[i].ptoffset;
                    for(int j = 0; j < edginfo[i].npoints; ++j)
                    {
                        int idx = cpts.index[offset+j];

                        PointGeomSharedPtr vert(
                            MemoryManager<PointGeom>::AllocateSharedPtr(
                                m_meshDimension, edginfo[i].id,
                                cpts.pts[idx].x, cpts.pts[idx].y,
                                cpts.pts[idx].z));
                        curve->m_points.push_back(vert);
                    }

                    m_curvedEdges[edgeid] = curve;
                }

                for(int i = 0; i < facinfo.size(); ++i)
                {
                    int faceid = facinfo[i].entityid;
                    LibUtilities::PointsType ptype;

                    CurveSharedPtr curve(
                            MemoryManager<Curve>::AllocateSharedPtr(
                                faceid, ptype = (LibUtilities::PointsType)
                                                    facinfo[i].ptype));

                    int offset = facinfo[i].ptoffset;
                    for(int j = 0; j < facinfo[i].npoints; ++j)
                    {
                        int idx = cpts.index[offset+j];

                        PointGeomSharedPtr vert(MemoryManager<PointGeom>::
                                   AllocateSharedPtr(m_meshDimension,
                                                     facinfo[i].id,
                                                     cpts.pts[idx].x,
                                                     cpts.pts[idx].y,
                                                     cpts.pts[idx].z));
                        curve->m_points.push_back(vert);
                    }

                    m_curvedFaces[faceid] = curve;
                }
            }
            else
            {
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
                    
                    /// Read edge id attribute.
                    err = edgelement->QueryIntAttribute("EDGEID", &edgeid);
                    ASSERTL0(err == TIXML_SUCCESS, "Unable to read curve attribute EDGEID.");

                    /// Read text edgelement description.
                    std::string elementStr;
                    TiXmlNode* elementChild = edgelement->FirstChild();
                    
                    while(elementChild)
                    {
                        // Accumulate all non-comment element data
                        if (elementChild->Type() == TiXmlNode::TINYXML_TEXT)
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
                        NekDouble xval, yval, zval;
                        std::istringstream elementDataStrm(elementStr.c_str());
                        try
                        {
                            while(!elementDataStrm.fail())
                            {
                                elementDataStrm >> xval >> yval >> zval;
                                
                                xval = xval*xscale + xmove;
                                yval = yval*yscale + ymove;
                                zval = zval*zscale + zmove;
                                
                                // Need to check it here because we may not be
                                // good after the read indicating that there
                                // was nothing to read.
                                if (!elementDataStrm.fail())
                                {
                                    PointGeomSharedPtr vert(MemoryManager<PointGeom>::AllocateSharedPtr(m_meshDimension, edgeindx, xval, yval, zval));
                                    
                                    curve->m_points.push_back(vert);
                                }
                                
                            }
                        }
                        catch(...)
                        {
                            NEKERROR(ErrorUtil::efatal,
                                     (std::string("Unable to read curve data for EDGE: ") + elementStr).c_str());
                            
                        }
                        
                        ASSERTL0(curve->m_points.size() == numPts,
                                 "Number of points specificed by attribute "
                                 "NUMPOINTS is different from number of points "
                                 "in list (edgeid = " +
                                 boost::lexical_cast<string>(edgeid));

                        m_curvedEdges[edgeid] = curve;
                        
                        edgelement = edgelement->NextSiblingElement("E");

                    } // end if-loop

                } // end while-loop

                TiXmlElement *facelement = field->FirstChildElement("F");
                int faceindx, faceid;
                
                while(facelement)
                {
                    std::string face(facelement->ValueStr());
                    ASSERTL0(face == "F", (std::string("Unknown 3D curve type: ") + face).c_str());
                    
                    /// Read id attribute.
                    err = facelement->QueryIntAttribute("ID", &faceindx);
                    ASSERTL0(err == TIXML_SUCCESS, "Unable to read curve attribute ID.");
                    
                    /// Read face id attribute.
                    err = facelement->QueryIntAttribute("FACEID", &faceid);
                    ASSERTL0(err == TIXML_SUCCESS, "Unable to read curve attribute FACEID.");
                    
                    /// Read text face element description.
                    std::string elementStr;
                    TiXmlNode* elementChild = facelement->FirstChild();
                    
                    while(elementChild)
                    {
                        // Accumulate all non-comment element data
                        if (elementChild->Type() == TiXmlNode::TINYXML_TEXT)
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
                        std::stringstream s;
                        s << numptsStr;
                        s >> numPts;
                        
                        CurveSharedPtr curve(MemoryManager<Curve>::AllocateSharedPtr(faceid, type));
                        
                        ASSERTL0(numPts >= 3, "NUMPOINTS for face must be greater than 2");
                        
                        if(numPts == 3)
                        {
                            ASSERTL0(ptsStr != endStr, "Invalid points type.");
                        }
                        
                        // Read points (x, y, z)
                        NekDouble xval, yval, zval;
                        std::istringstream elementDataStrm(elementStr.c_str());
                        try
                        {
                            while(!elementDataStrm.fail())
                            {
                                elementDataStrm >> xval >> yval >> zval;

                                // Need to check it here because we
                                // may not be good after the read
                                // indicating that there was nothing
                                // to read.
                                if (!elementDataStrm.fail())
                                {
                                    PointGeomSharedPtr vert(MemoryManager<PointGeom>::AllocateSharedPtr(m_meshDimension, faceindx, xval, yval, zval));
                                    curve->m_points.push_back(vert);
                                }
                            }
                        }
                        catch(...)
                        {
                            NEKERROR(ErrorUtil::efatal,
                                     (std::string("Unable to read curve data for FACE: ")
                                      + elementStr).c_str());
                        }
                        m_curvedFaces[faceid] = curve;
                        
                        facelement = facelement->NextSiblingElement("F");
                        
                    } // end if-loop
                } // end while-loop
            } // end of compressed else
        } // end of ReadCurves()

        /**
         *
         */
        void MeshGraph::ReadCurves(std::string &infilename)
        {
            TiXmlDocument doc(infilename);
            bool loadOkay = doc.LoadFile();

            std::stringstream errstr;
            errstr << "Unable to load file: " << infilename << std::endl;
            errstr << "Reason: " << doc.ErrorDesc() << std::endl;
            errstr << "Position: Line " << doc.ErrorRow() << ", Column " << doc.ErrorCol() << std::endl;
            ASSERTL0(loadOkay, errstr.str());

            ReadCurves(doc);
        }

        /**
         * @brief Populate a TinyXML document with a GEOMETRY tag inside the
         * NEKTAR tag.
         *
         * This routine will create a GEOMETRY XML tag which represents the
         * MeshGraph object. If a NEKTAR tag does not already exist, it will
         * create one. If a GEOMETRY block already exists inside the NEKTAR tag,
         * it will overwrite it.
         */
        void MeshGraph::WriteGeometry(TiXmlDocument &doc)
        {
            TiXmlElement *root = doc.FirstChildElement("NEKTAR");
            TiXmlElement *geomTag;

            // Try to find existing NEKTAR tag.
            if (!root)
            {
                root = new TiXmlElement("NEKTAR");
                doc.LinkEndChild(root);

                geomTag = new TiXmlElement("GEOMETRY");
                root->LinkEndChild(geomTag);
            }
            else
            {
                // Try to find existing GEOMETRY tag.
                geomTag = root->FirstChildElement("GEOMETRY");

                if (!geomTag)
                {
                    geomTag = new TiXmlElement("GEOMETRY");
                    root->LinkEndChild(geomTag);
                }
            }

            // Update attributes with dimensions.
            geomTag->SetAttribute("DIM",   m_meshDimension);
            geomTag->SetAttribute("SPACE", m_spaceDimension);

            // Clear existing elements.
            geomTag->Clear();

            // Construct <VERTEX> block
            TiXmlElement *vertTag = new TiXmlElement("VERTEX");
            PointGeomMap::iterator pIt;

            for (pIt = m_vertSet.begin(); pIt != m_vertSet.end(); ++pIt)
            {
                stringstream s;
                s << scientific << setprecision(8)
                  << (*pIt->second)(0) << " " << (*pIt->second)(1) << " "
                  << (*pIt->second)(2);
                TiXmlElement * v = new TiXmlElement("V");
                v->SetAttribute("ID", pIt->second->GetVid());
                v->LinkEndChild(new TiXmlText(s.str()));
                vertTag->LinkEndChild(v);
            }

            geomTag->LinkEndChild(vertTag);

            // Construct <EDGE> or <ELEMENT> block
            TiXmlElement *edgeTag = new TiXmlElement(
                m_meshDimension == 1 ? "ELEMENT" : "EDGE");
            SegGeomMap::iterator sIt;
            string tag = m_meshDimension == 1 ? "S" : "E";

            for (sIt = m_segGeoms.begin(); sIt != m_segGeoms.end(); ++sIt)
            {
                stringstream s;
                SegGeomSharedPtr seg = sIt->second;
                s << seg->GetVid(0) << " " << seg->GetVid(1);
                TiXmlElement *e = new TiXmlElement(tag);
                e->SetAttribute("ID", sIt->first);
                e->LinkEndChild(new TiXmlText(s.str()));
                edgeTag->LinkEndChild(e);
            }

            geomTag->LinkEndChild(edgeTag);

            // Construct <FACE> or <ELEMENT> block
            if (m_meshDimension > 1)
            {
                TiXmlElement *faceTag = new TiXmlElement(
                    m_meshDimension == 2 ? "ELEMENT" : "FACE");

                TriGeomMap::iterator tIt;
                tag = "T";

                for (tIt = m_triGeoms.begin(); tIt != m_triGeoms.end(); ++tIt)
                {
                    stringstream s;
                    TriGeomSharedPtr tri = tIt->second;
                    s << tri->GetEid(0) << " " << tri->GetEid(1) << " "
                      << tri->GetEid(2);
                    TiXmlElement *t = new TiXmlElement(tag);
                    t->SetAttribute("ID", tIt->first);
                    t->LinkEndChild(new TiXmlText(s.str()));
                    faceTag->LinkEndChild(t);
                }

                QuadGeomMap::iterator qIt;
                tag = "Q";

                for (qIt = m_quadGeoms.begin(); qIt != m_quadGeoms.end(); ++qIt)
                {
                    stringstream s;
                    QuadGeomSharedPtr quad = qIt->second;
                    s << quad->GetEid(0) << " " << quad->GetEid(1) << " "
                      << quad->GetEid(2) << " " << quad->GetEid(3);
                    TiXmlElement *q = new TiXmlElement(tag);
                    q->SetAttribute("ID", qIt->first);
                    q->LinkEndChild(new TiXmlText(s.str()));
                    faceTag->LinkEndChild(q);
                }

                geomTag->LinkEndChild(faceTag);
            }

            if (m_meshDimension > 2)
            {
                TiXmlElement *elmtTag = new TiXmlElement("ELEMENT");

                HexGeomMap::iterator hIt;
                tag = "H";

                for (hIt = m_hexGeoms.begin(); hIt != m_hexGeoms.end(); ++hIt)
                {
                    stringstream s;
                    HexGeomSharedPtr hex = hIt->second;
                    s << hex->GetFid(0) << " " << hex->GetFid(1) << " "
                      << hex->GetFid(2) << " " << hex->GetFid(3) << " "
                      << hex->GetFid(4) << " " << hex->GetFid(5) << " ";
                    TiXmlElement *h = new TiXmlElement(tag);
                    h->SetAttribute("ID", hIt->first);
                    h->LinkEndChild(new TiXmlText(s.str()));
                    elmtTag->LinkEndChild(h);
                }

                PrismGeomMap::iterator rIt;
                tag = "R";

                for (rIt = m_prismGeoms.begin(); rIt != m_prismGeoms.end(); ++rIt)
                {
                    stringstream s;
                    PrismGeomSharedPtr prism = rIt->second;
                    s << prism->GetFid(0) << " " << prism->GetFid(1) << " "
                      << prism->GetFid(2) << " " << prism->GetFid(3) << " "
                      << prism->GetFid(4) << " ";
                    TiXmlElement *p = new TiXmlElement(tag);
                    p->SetAttribute("ID", rIt->first);
                    p->LinkEndChild(new TiXmlText(s.str()));
                    elmtTag->LinkEndChild(p);
                }

                PyrGeomMap::iterator pIt;
                tag = "P";

                for (pIt = m_pyrGeoms.begin(); pIt != m_pyrGeoms.end(); ++pIt)
                {
                    stringstream s;
                    PyrGeomSharedPtr pyr = pIt->second;
                    s << pyr->GetFid(0) << " " << pyr->GetFid(1) << " "
                      << pyr->GetFid(2) << " " << pyr->GetFid(3) << " "
                      << pyr->GetFid(4) << " ";
                    TiXmlElement *p = new TiXmlElement(tag);
                    p->SetAttribute("ID", pIt->first);
                    p->LinkEndChild(new TiXmlText(s.str()));
                    elmtTag->LinkEndChild(p);
                }

                TetGeomMap::iterator tIt;
                tag = "A";

                for (tIt = m_tetGeoms.begin(); tIt != m_tetGeoms.end(); ++tIt)
                {
                    stringstream s;
                    TetGeomSharedPtr tet = tIt->second;
                    s << tet->GetFid(0) << " " << tet->GetFid(1) << " "
                      << tet->GetFid(2) << " " << tet->GetFid(3) << " ";
                    TiXmlElement *t = new TiXmlElement(tag);
                    t->SetAttribute("ID", tIt->first);
                    t->LinkEndChild(new TiXmlText(s.str()));
                    elmtTag->LinkEndChild(t);
                }

                geomTag->LinkEndChild(elmtTag);
            }

            // Construct <CURVED> block
            TiXmlElement *curveTag = new TiXmlElement("CURVED");
            CurveMap::iterator curveIt;
            int curveId = 0;

            for (curveIt  = m_curvedEdges.begin();
                 curveIt != m_curvedEdges.end(); ++curveIt)
            {
                CurveSharedPtr curve = curveIt->second;
                TiXmlElement *c = new TiXmlElement("E");
                stringstream s;
                s.precision(8);

                for (int j = 0; j < curve->m_points.size(); ++j)
                {
                    SpatialDomains::PointGeomSharedPtr p = curve->m_points[j];
                    s << scientific << (*p)(0) << " " << (*p)(1) << " " << (*p)(2) << "   ";
                }

                c->SetAttribute("ID", curveId++);
                c->SetAttribute("EDGEID", curve->m_curveID);
                c->SetAttribute("NUMPOINTS", curve->m_points.size());
                c->SetAttribute("TYPE", LibUtilities::kPointsTypeStr[curve->m_ptype]);
                c->LinkEndChild(new TiXmlText(s.str()));
                curveTag->LinkEndChild(c);
            }

            for (curveIt  = m_curvedFaces.begin();
                 curveIt != m_curvedFaces.end(); ++curveIt)
            {
                CurveSharedPtr curve = curveIt->second;
                TiXmlElement *c = new TiXmlElement("F");
                stringstream s;
                s.precision(8);

                for (int j = 0; j < curve->m_points.size(); ++j)
                {
                    SpatialDomains::PointGeomSharedPtr p = curve->m_points[j];
                    s << scientific << (*p)(0) << " " << (*p)(1) << " " << (*p)(2) << "   ";
                }

                c->SetAttribute("ID", curveId++);
                c->SetAttribute("FACEID", curve->m_curveID);
                c->SetAttribute("NUMPOINTS", curve->m_points.size());
                c->SetAttribute("TYPE", LibUtilities::kPointsTypeStr[curve->m_ptype]);
                c->LinkEndChild(new TiXmlText(s.str()));
                curveTag->LinkEndChild(c);
            }

            geomTag->LinkEndChild(curveTag);

            // Construct <COMPOSITE> blocks
            TiXmlElement *compTag = new TiXmlElement("COMPOSITE");
            CompositeMap::iterator cIt;

            // Create a map that gets around the issue of mapping faces -> F and
            // edges -> E inside the tag.
            map<LibUtilities::ShapeType, pair<string, string> > compMap;
            compMap[LibUtilities::eSegment]       = make_pair("S", "E");
            compMap[LibUtilities::eQuadrilateral] = make_pair("Q", "F");
            compMap[LibUtilities::eTriangle]      = make_pair("T", "F");
            compMap[LibUtilities::eTetrahedron]   = make_pair("A", "A");
            compMap[LibUtilities::ePyramid]       = make_pair("P", "P");
            compMap[LibUtilities::ePrism]         = make_pair("R", "R");
            compMap[LibUtilities::eHexahedron]    = make_pair("H", "H");

            std::vector<unsigned int> idxList;

            for (cIt = m_meshComposites.begin(); cIt != m_meshComposites.end(); ++cIt)
            {
                stringstream s;
                TiXmlElement *c = new TiXmlElement("C");
                GeometrySharedPtr firstGeom = cIt->second->at(0);
                int shapeDim = firstGeom->GetShapeDim();
                string tag = (shapeDim < m_meshDimension) ?
                    compMap[firstGeom->GetShapeType()].second :
                    compMap[firstGeom->GetShapeType()].first;

                idxList.clear();
                s << " " << tag << "[";

                for (int i = 0; i < cIt->second->size(); ++i)
                {
                    idxList.push_back((*cIt->second)[i]->GetGlobalID());
                }

                s << ParseUtils::GenerateSeqString(idxList) << "] ";

                c->SetAttribute("ID", cIt->first);
                c->LinkEndChild(new TiXmlText(s.str()));
                compTag->LinkEndChild(c);
            }

            geomTag->LinkEndChild(compTag);

            // Construct <DOMAIN> block
            TiXmlElement *domTag = new TiXmlElement("DOMAIN");
            stringstream domString;

            // @todo Fix this to accomodate multi domain output
            idxList.clear();
            for (cIt = m_domain[0].begin(); cIt != m_domain[0].end(); ++cIt)
            {
                idxList.push_back(cIt->first);
            }

            domString << " C[" << ParseUtils::GenerateSeqString(idxList) << "] ";
            domTag->LinkEndChild(new TiXmlText(domString.str()));
            geomTag->LinkEndChild(domTag);
        }

        /**
         * @brief Write out an XML file containing the GEOMETRY block
         * representing this MeshGraph instance inside a NEKTAR tag.
         */
        void MeshGraph::WriteGeometry(std::string &outfilename)
        {
            // Create empty TinyXML document.
            TiXmlDocument doc;
            TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "utf-8", "");
            doc.LinkEndChild(decl);

            // Write out geometry information.
            WriteGeometry(doc);

            // Save file.
            doc.SaveFile(outfilename);
        }

        void MeshGraph::SetDomainRange (
            NekDouble xmin, NekDouble xmax, NekDouble ymin,
            NekDouble ymax, NekDouble zmin, NekDouble zmax)
        {
            m_domainRange->m_checkShape = false;

            if(m_domainRange == NullDomainRangeShPtr)
            {
                m_domainRange = MemoryManager<DomainRange>::AllocateSharedPtr();
                m_domainRange->m_doXrange = true;
            }

            m_domainRange->m_xmin = xmin;
            m_domainRange->m_xmax = xmax;

            if(ymin == NekConstants::kNekUnsetDouble)
            {
                m_domainRange->m_doYrange = false;
            }
            else
            {
                m_domainRange->m_doYrange = true;
                m_domainRange->m_ymin = ymin;
                m_domainRange->m_ymax = ymax;
            }

            if(zmin == NekConstants::kNekUnsetDouble)
            {
                m_domainRange->m_doZrange = false;
            }
            else
            {
                m_domainRange->m_doZrange = true;
                m_domainRange->m_zmin = zmin;
                m_domainRange->m_zmax = zmax;
            }
        }

        bool MeshGraph::CheckRange(Geometry2D &geom)
        {
            bool returnval = true;

            if(m_domainRange  != NullDomainRangeShPtr)
            {
                int nverts  = geom.GetNumVerts();
                int coordim = geom.GetCoordim();

                // exclude elements outside x range if all vertices not in region
                if(m_domainRange->m_doXrange)
                {
                    int ncnt_low = 0;
                    int ncnt_up = 0;
                    for(int i = 0; i < nverts; ++i)
                    {
                        NekDouble xval = (*geom.GetVertex(i))[0];
                        if(xval < m_domainRange->m_xmin)
                        {
                            ncnt_low++;
                        }

                        if(xval > m_domainRange->m_xmax)
                        {
                            ncnt_up++;
                        }
                    }

                    // check for all verts to be less or greater than
                    // range so that if element spans thin range then
                    // it is still included
                    if((ncnt_up == nverts)||(ncnt_low == nverts))
                    {
                        returnval = false;
                    }
                }

                // exclude elements outside y range if all vertices not in region
                if(m_domainRange->m_doYrange)
                {
                    int ncnt_low = 0;
                    int ncnt_up  = 0;
                    for(int i = 0; i < nverts; ++i)
                    {
                        NekDouble yval = (*geom.GetVertex(i))[1];
                        if(yval < m_domainRange->m_ymin)
                        {
                            ncnt_low++;
                        }

                        if(yval > m_domainRange->m_ymax)
                        {
                            ncnt_up++;
                        }
                    }

                    // check for all verts to be less or greater than
                    // range so that if element spans thin range then
                    // it is still included
                    if((ncnt_up == nverts)||(ncnt_low == nverts))
                    {
                        returnval = false;
                    }
                }

                if(coordim > 2)
                {
                    // exclude elements outside z range if all vertices not in region
                    if(m_domainRange->m_doZrange)
                    {
                        int ncnt_low = 0;
                        int ncnt_up  = 0;

                        for(int i = 0; i < nverts; ++i)
                        {
                            NekDouble zval = (*geom.GetVertex(i))[2];

                            if(zval < m_domainRange->m_zmin)
                            {
                                ncnt_low++;
                            }

                            if(zval > m_domainRange->m_zmax)
                            {
                                ncnt_up++;
                            }
                        }

                        // check for all verts to be less or greater than
                        // range so that if element spans thin range then
                        // it is still included
                        if((ncnt_up == nverts)||(ncnt_low == nverts))
                        {
                            returnval = false;
                        }
                    }
                }
            }
            return returnval;
        }


        /* Domain checker for 3D geometries */
        bool MeshGraph::CheckRange(Geometry3D &geom)
        {
            bool returnval = true;

            if(m_domainRange  != NullDomainRangeShPtr)
            {
                int nverts  = geom.GetNumVerts();

                if(m_domainRange->m_doXrange)
                {
                    int ncnt_low = 0;
                    int ncnt_up = 0;

                    for(int i = 0; i < nverts; ++i)
                    {
                        NekDouble xval = (*geom.GetVertex(i))[0];
                        if(xval < m_domainRange->m_xmin)
                        {
                            ncnt_low++;
                        }

                        if(xval > m_domainRange->m_xmax)
                        {
                            ncnt_up++;
                        }
                    }

                    // check for all verts to be less or greater than
                    // range so that if element spans thin range then
                    // it is still included
                    if((ncnt_up == nverts)||(ncnt_low == nverts))
                    {
                        returnval = false;
                    }
                }

                if(m_domainRange->m_doYrange)
                {
                    int ncnt_low = 0;
                    int ncnt_up = 0;
                    for(int i = 0; i < nverts; ++i)
                    {
                        NekDouble yval = (*geom.GetVertex(i))[1];
                        if(yval < m_domainRange->m_ymin)
                        {
                            ncnt_low++;
                        }

                        if(yval > m_domainRange->m_ymax)
                        {
                            ncnt_up++;
                        }
                    }

                    // check for all verts to be less or greater than
                    // range so that if element spans thin range then
                    // it is still included
                    if((ncnt_up == nverts)||(ncnt_low == nverts))
                    {
                        returnval = false;
                    }
                }

                if(m_domainRange->m_doZrange)
                {
                    int ncnt_low = 0;
                    int ncnt_up  = 0;
                    for(int i = 0; i < nverts; ++i)
                    {
                        NekDouble zval = (*geom.GetVertex(i))[2];

                        if(zval < m_domainRange->m_zmin)
                        {
                            ncnt_low++;
                        }

                        if(zval > m_domainRange->m_zmax)
                        {
                            ncnt_up++;
                        }
                    }

                    // check for all verts to be less or greater than
                    // range so that if element spans thin range then
                    // it is still included
                    if((ncnt_up == nverts)||(ncnt_low == nverts))
                    {
                        returnval = false;
                    }
                }

                if(m_domainRange->m_checkShape)
                {
                    if(geom.GetShapeType() != m_domainRange->m_shapeType)
                    {
                        returnval = false;
                    }
                }

            }

            return returnval;
        }

        /**
         *
         */
        GeometrySharedPtr MeshGraph::GetCompositeItem(int whichComposite, int whichItem)
        {
            GeometrySharedPtr returnval;
            bool error = false;

            if (whichComposite >= 0 && whichComposite < int(m_meshComposites.size()))
            {
                if (whichItem >= 0 && whichItem < int(m_meshComposites[whichComposite]->size()))
                {
                    returnval = m_meshComposites[whichComposite]->at(whichItem);
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


        /**
         *
         */
        void MeshGraph::GetCompositeList(const std::string &compositeStr, CompositeMap &compositeVector) const
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

                    // If the composite listed is not found and we are working
                    // on a partitioned mesh, silently ignore it.
                    if (m_meshComposites.find(*iter) == m_meshComposites.end()
                            && m_meshPartitioned)
                    {
                        continue;
                    }

                    addedVector.push_back(*iter);
                    Composite composite = GetComposite(*iter);
                    CompositeMap::iterator compIter;
                    if (composite)
                    {
                        compositeVector[*iter] = composite;
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


        /**
         *
         */
        const ExpansionMap &MeshGraph::GetExpansions(const std::string variable)
        {
            ExpansionMapShPtr returnval;

            if(m_expansionMapShPtrMap.count(variable))
            {
                returnval = m_expansionMapShPtrMap.find(variable)->second;
            }
            else
            {
                if(m_expansionMapShPtrMap.count("DefaultVar") == 0)
                {
                    NEKERROR(ErrorUtil::efatal, (std::string("Unable to find expansion vector definition for field: ")+variable).c_str());
                }
                returnval = m_expansionMapShPtrMap.find("DefaultVar")->second;
                m_expansionMapShPtrMap[variable] = returnval;

                NEKERROR(ErrorUtil::ewarning, (std::string("Using Default variable expansion definition for field: ")+variable).c_str());
            }

            return *returnval;
        }


        /**
         *
         */
        ExpansionShPtr MeshGraph::GetExpansion(GeometrySharedPtr geom, const std::string variable)
        {
            ExpansionMapIter iter;
            ExpansionShPtr returnval;

            ExpansionMapShPtr expansionMap = m_expansionMapShPtrMap.find(variable)->second;

            for (iter = expansionMap->begin(); iter!=expansionMap->end(); ++iter)
            {
                if ((iter->second)->m_geomShPtr == geom)
                {
                    returnval = iter->second;
                    break;
                }
            }
            return returnval;
        }


        /**
         *
         */
        void MeshGraph::SetExpansions(
                std::vector<LibUtilities::FieldDefinitionsSharedPtr> &fielddef)
        {
            int i, j, k, cnt, id;
            GeometrySharedPtr geom;

            ExpansionMapShPtr expansionMap;

            // Loop over fields and determine unique fields string and
            // declare whole expansion list
            for(i = 0; i < fielddef.size(); ++i)
            {
                for(j = 0; j < fielddef[i]->m_fields.size(); ++j)
                {
                    std::string field = fielddef[i]->m_fields[j];
                    if(m_expansionMapShPtrMap.count(field) == 0)
                    {
                        expansionMap = MemoryManager<ExpansionMap>::AllocateSharedPtr();
                        m_expansionMapShPtrMap[field] = expansionMap;

                        // check to see if DefaultVar also not set and
                        // if so assign it to this expansion
                        if(m_expansionMapShPtrMap.count("DefaultVar") == 0)
                        {
                            m_expansionMapShPtrMap["DefaultVar"] = expansionMap;
                        }

                        // loop over all elements and set expansion
                        for(k = 0; k < fielddef.size(); ++k)
                        {
                            for(int h = 0; h < fielddef[k]->m_fields.size(); ++h)
                            {
                                if(fielddef[k]->m_fields[h] == field)
                                {
                                    expansionMap = m_expansionMapShPtrMap.find(field)->second;
                                    LibUtilities::BasisKeyVector def;

                                    for(int g = 0; g < fielddef[k]->m_elementIDs.size(); ++g)
                                    {
                                        ExpansionShPtr tmpexp =
                                                MemoryManager<Expansion>::AllocateSharedPtr(geom, def);
                                        (*expansionMap)[fielddef[k]->m_elementIDs[g]] = tmpexp;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // loop over all elements find the geometry shared ptr and
            // set up basiskey vector
            for(i = 0; i < fielddef.size(); ++i)
            {
                cnt = 0;
                std::vector<std::string>  fields = fielddef[i]->m_fields;
                std::vector<unsigned int> nmodes = fielddef[i]->m_numModes;
                std::vector<LibUtilities::BasisType> basis = fielddef[i]->m_basis;
                bool pointDef = fielddef[i]->m_pointsDef;
                bool numPointDef = fielddef[i]->m_numPointsDef;

                // Check points and numpoints
                std::vector<unsigned int> npoints = fielddef[i]->m_numPoints;
                std::vector<LibUtilities::PointsType> points = fielddef[i]->m_points;

                bool UniOrder =  fielddef[i]->m_uniOrder;

                for (j = 0; j < fielddef[i]->m_elementIDs.size(); ++j)
                {

                    LibUtilities::BasisKeyVector bkeyvec;
                    id = fielddef[i]->m_elementIDs[j];

                    switch (fielddef[i]->m_shapeType)
                    {
                    case LibUtilities::eSegment:
                        {
                            if(m_segGeoms.count(fielddef[i]->m_elementIDs[j]) == 0)
                            {
                                // skip element likely from parallel read
                                continue;
                            }
                            geom = m_segGeoms[fielddef[i]->m_elementIDs[j]];

                            LibUtilities::PointsKey pkey(nmodes[cnt]+1, LibUtilities::eGaussLobattoLegendre);

                            if(numPointDef&&pointDef)
                            {
                                const LibUtilities::PointsKey pkey1(npoints[cnt], points[0]);
                                pkey = pkey1;
                            }
                            else if(!numPointDef&&pointDef)
                            {
                                const LibUtilities::PointsKey pkey1(nmodes[cnt]+1, points[0]);
                                pkey = pkey1;
                            }
                            else if(numPointDef&&!pointDef)
                            {
                                const LibUtilities::PointsKey pkey1(npoints[cnt], LibUtilities::eGaussLobattoLegendre);
                                pkey = pkey1;
                            }

                            LibUtilities::BasisKey bkey(basis[0], nmodes[cnt], pkey);

                            if(!UniOrder)
                            {
                                cnt++;
                                cnt += fielddef[i]->m_numHomogeneousDir;
                            }
                            bkeyvec.push_back(bkey);
                        }
                        break;
                    case LibUtilities::eTriangle:
                        {
                            if(m_triGeoms.count(fielddef[i]->m_elementIDs[j]) == 0)
                            {
                                // skip element likely from parallel read
                                continue;
                            }
                            geom = m_triGeoms[fielddef[i]->m_elementIDs[j]];

                            LibUtilities::PointsKey pkey(nmodes[cnt]+1, LibUtilities::eGaussLobattoLegendre);
                            if(numPointDef&&pointDef)
                            {
                                const LibUtilities::PointsKey pkey2(npoints[cnt], points[0]);
                                pkey = pkey2;
                            }
                            else if(!numPointDef&&pointDef)
                            {
                                const LibUtilities::PointsKey pkey2(nmodes[cnt]+1, points[0]);
                                pkey = pkey2;
                            }
                            else if(numPointDef&&!pointDef)
                            {
                                const LibUtilities::PointsKey pkey2(npoints[cnt], LibUtilities::eGaussLobattoLegendre);
                                pkey = pkey2;
                            }
                            LibUtilities::BasisKey bkey(basis[0], nmodes[cnt], pkey);

                            bkeyvec.push_back(bkey);

                            LibUtilities::PointsKey pkey1(nmodes[cnt+1], LibUtilities::eGaussRadauMAlpha1Beta0);
                            if(numPointDef&&pointDef)
                            {
                                const LibUtilities::PointsKey pkey2(npoints[cnt+1], points[1]);
                                pkey1 = pkey2;
                            }
                            else if(!numPointDef&&pointDef)
                            {
                                const LibUtilities::PointsKey pkey2(nmodes[cnt+1], points[1]);
                                pkey1 = pkey2;
                            }
                            else if(numPointDef&&!pointDef)
                            {
                                const LibUtilities::PointsKey pkey2(npoints[cnt+1], LibUtilities::eGaussRadauMAlpha1Beta0);
                                pkey1 = pkey2;
                            }
                            LibUtilities::BasisKey bkey1(basis[1], nmodes[cnt+1], pkey1);
                            bkeyvec.push_back(bkey1);

                            if(!UniOrder)
                            {
                                cnt += 2;
                                cnt += fielddef[i]->m_numHomogeneousDir;
                            }
                        }
                        break;
                    case LibUtilities::eQuadrilateral:
                        {
                            if(m_quadGeoms.count(fielddef[i]->m_elementIDs[j]) == 0)
                            {
                                // skip element likely from parallel read
                                continue;
                            }

                            geom = m_quadGeoms[fielddef[i]->m_elementIDs[j]];

                            for(int b = 0; b < 2; ++b)
                            {
                                LibUtilities::PointsKey pkey(nmodes[cnt+b]+1, LibUtilities::eGaussLobattoLegendre);

                                if(numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(npoints[cnt+b],points[b]);
                                    pkey = pkey2;
                                }
                                else if(!numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(nmodes[cnt+b]+1,points[b]);
                                    pkey = pkey2;
                                }
                                else if(numPointDef&&!pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(npoints[cnt+b],LibUtilities::eGaussLobattoLegendre);
                                    pkey = pkey2;
                                }
                                LibUtilities::BasisKey bkey(basis[b],nmodes[cnt+b],pkey);
                                bkeyvec.push_back(bkey);
                            }

                            if(!UniOrder)
                            {
                                cnt += 2;
                                cnt += fielddef[i]->m_numHomogeneousDir;
                            }
                        }
                        break;

                        case LibUtilities::eTetrahedron:
                        {
                            k = fielddef[i]->m_elementIDs[j];

                            // allow for possibility that fielddef is
                            // larger than m_graph which can happen in
                            // parallel runs
                            if(m_tetGeoms.count(k) == 0)
                            {
                                continue;
                            }
                            geom = m_tetGeoms[k];

                            {
                                LibUtilities::PointsKey pkey(nmodes[cnt]+1, LibUtilities::eGaussLobattoLegendre);

                                if(numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(npoints[cnt],points[0]);
                                    pkey = pkey2;
                                }
                                else if(!numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(nmodes[cnt]+1,points[0]);
                                    pkey = pkey2;
                                }
                                else if(numPointDef&&!pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(npoints[cnt],LibUtilities::eGaussLobattoLegendre);
                                    pkey = pkey2;
                                }

                                LibUtilities::BasisKey bkey(basis[0],nmodes[cnt],pkey);

                                bkeyvec.push_back(bkey);
                            }
                            {
                                LibUtilities::PointsKey pkey(nmodes[cnt+1], LibUtilities::eGaussRadauMAlpha1Beta0);

                                if(numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(npoints[cnt+1],points[1]);
                                    pkey = pkey2;
                                }
                                else if(!numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(nmodes[cnt+1]+1,points[1]);
                                    pkey = pkey2;
                                }
                                else if(numPointDef&&!pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(npoints[cnt+1],LibUtilities::eGaussRadauMAlpha1Beta0);
                                    pkey = pkey2;
                                }

                                LibUtilities::BasisKey bkey(basis[1],nmodes[cnt+1],pkey);

                                bkeyvec.push_back(bkey);
                            }

                            {
                                LibUtilities::PointsKey pkey(nmodes[cnt+2], LibUtilities::eGaussRadauMAlpha2Beta0);

                                if(numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(npoints[cnt+2],points[2]);
                                    pkey = pkey2;
                                }
                                else if(!numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(nmodes[cnt+2]+1,points[2]);
                                    pkey = pkey2;
                                }
                                else if(numPointDef&&!pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(npoints[cnt+2],LibUtilities::eGaussRadauMAlpha1Beta0);
                                    pkey = pkey2;
                                }

                                LibUtilities::BasisKey bkey(basis[2],nmodes[cnt+2],pkey);

                                bkeyvec.push_back(bkey);
                            }

                            if(!UniOrder)
                            {
                                cnt += 3;
                            }
                        }
                        break;
                        case LibUtilities::ePrism:
                        {
                            k = fielddef[i]->m_elementIDs[j];
                            if(m_prismGeoms.count(k) == 0)
                            {
                                continue;
                            }
                            geom = m_prismGeoms[k];

                            for(int b = 0; b < 2; ++b)
                            {
                                LibUtilities::PointsKey pkey(nmodes[cnt+b]+1,LibUtilities::eGaussLobattoLegendre);

                                if(numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(npoints[cnt+b],points[b]);
                                    pkey = pkey2;
                                }
                                else if(!numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(nmodes[cnt+b]+1,points[b]);
                                    pkey = pkey2;
                                }
                                else if(numPointDef&&!pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(npoints[cnt+b],LibUtilities::eGaussLobattoLegendre);
                                    pkey = pkey2;
                                }

                                LibUtilities::BasisKey bkey(basis[b],nmodes[cnt+b],pkey);
                                bkeyvec.push_back(bkey);
                            }

                            {
                                LibUtilities::PointsKey pkey(nmodes[cnt+2],LibUtilities::eGaussRadauMAlpha1Beta0);

                                if(numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(npoints[cnt+2],points[2]);
                                    pkey = pkey2;
                                }
                                else if(!numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(nmodes[cnt+2]+1,points[2]);
                                    pkey = pkey2;
                                }
                                else if(numPointDef&&!pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(npoints[cnt+2],LibUtilities::eGaussLobattoLegendre);
                                    pkey = pkey2;
                                }

                                LibUtilities::BasisKey bkey(basis[2],nmodes[cnt+2],pkey);
                                bkeyvec.push_back(bkey);
                            }

                            if(!UniOrder)
                            {
                                cnt += 3;
                            }
                        }
                        break;
                        case LibUtilities::ePyramid:
                        {
                            k = fielddef[i]->m_elementIDs[j];
                            ASSERTL0(m_pyrGeoms.find(k) != m_pyrGeoms.end(),
                                    "Failed to find geometry with same global id");
                            geom = m_pyrGeoms[k];


                            for(int b = 0; b < 2; ++b)
                            {
                                LibUtilities::PointsKey pkey(nmodes[cnt+b]+1,LibUtilities::eGaussLobattoLegendre);

                                if(numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(npoints[cnt+b],points[b]);
                                    pkey = pkey2;
                                }
                                else if(!numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(nmodes[cnt+b]+1,points[b]);
                                    pkey = pkey2;
                                }
                                else if(numPointDef&&!pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(npoints[cnt+b],LibUtilities::eGaussLobattoLegendre);
                                    pkey = pkey2;
                                }

                                LibUtilities::BasisKey bkey(basis[b],nmodes[cnt+b],pkey);
                                bkeyvec.push_back(bkey);
                            }

                            {
                                LibUtilities::PointsKey pkey(nmodes[cnt+2],LibUtilities::eGaussRadauMAlpha2Beta0);

                                if(numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(npoints[cnt+2],points[2]);
                                    pkey = pkey2;
                                }
                                else if(!numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(nmodes[cnt+2]+1,points[2]);
                                    pkey = pkey2;
                                }
                                else if(numPointDef&&!pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(npoints[cnt+2],LibUtilities::eGaussLobattoLegendre);
                                    pkey = pkey2;
                                }

                                LibUtilities::BasisKey bkey(basis[2],nmodes[cnt+2],pkey);
                                bkeyvec.push_back(bkey);
                            }

                            if(!UniOrder)
                            {
                                cnt += 3;
                            }
                        }
                        break;
                        case LibUtilities::eHexahedron:
                        {
                            k = fielddef[i]->m_elementIDs[j];
                            if(m_hexGeoms.count(k) == 0)
                            {
                                continue;
                            }

                            geom = m_hexGeoms[k];

                            for(int b = 0; b < 3; ++b)
                            {
                                LibUtilities::PointsKey pkey(nmodes[cnt+b],LibUtilities::eGaussLobattoLegendre);

                                if(numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(npoints[cnt+b],points[b]);
                                    pkey = pkey2;
                                }
                                else if(!numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(nmodes[cnt+b]+1,points[b]);
                                    pkey = pkey2;
                                }
                                else if(numPointDef&&!pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(npoints[cnt+b],LibUtilities::eGaussLobattoLegendre);
                                    pkey = pkey2;
                                }

                                LibUtilities::BasisKey bkey(basis[b],nmodes[cnt+b],pkey);
                                bkeyvec.push_back(bkey);
                            }

                            if(!UniOrder)
                            {
                                cnt += 3;
                            }
                        }
                        break;
                        default:
                            ASSERTL0(false,"Need to set up for pyramid and prism 3D Expansions");
                            break;
                        }

                        for(k = 0; k < fields.size(); ++k)
                        {
                            expansionMap = m_expansionMapShPtrMap.find(fields[k])->second;
                            if((*expansionMap).find(id) != (*expansionMap).end())
                            {
                                (*expansionMap)[id]->m_geomShPtr = geom;
                                (*expansionMap)[id]->m_basisKeyVector = bkeyvec;
                            }
                        }
                }
            }
        }


        /**
         *
         */
        void MeshGraph::SetExpansions(
                std::vector<LibUtilities::FieldDefinitionsSharedPtr> &fielddef,
                std::vector< std::vector<LibUtilities::PointsType> > &pointstype)
        {
            int i,j,k,g,h,cnt,id;
            GeometrySharedPtr geom;

            ExpansionMapShPtr expansionMap;

            // Loop over fields and determine unique fields string and
            // declare whole expansion list
            for(i = 0; i < fielddef.size(); ++i)
            {
                for(j = 0; j < fielddef[i]->m_fields.size(); ++j)
                {
                    std::string field = fielddef[i]->m_fields[j];
                    if(m_expansionMapShPtrMap.count(field) == 0)
                    {
                        expansionMap = MemoryManager<ExpansionMap>::AllocateSharedPtr();
                        m_expansionMapShPtrMap[field] = expansionMap;

                        // check to see if DefaultVar also not set and
                        // if so assign it to this expansion
                        if(m_expansionMapShPtrMap.count("DefaultVar") == 0)
                        {
                            m_expansionMapShPtrMap["DefaultVar"] = expansionMap;
                        }

                        // loop over all elements and set expansion
                        for(k = 0; k < fielddef.size(); ++k)
                        {
                            for(h = 0; h < fielddef[k]->m_fields.size(); ++h)
                            {
                                if(fielddef[k]->m_fields[h] == field)
                                {
                                    expansionMap = m_expansionMapShPtrMap.find(field)->second;
                                    LibUtilities::BasisKeyVector def;

                                    for(g = 0; g < fielddef[k]->m_elementIDs.size(); ++g)
                                    {
                                        ExpansionShPtr tmpexp =
                                                MemoryManager<Expansion>::AllocateSharedPtr(geom, def);
                                        (*expansionMap)[fielddef[k]->m_elementIDs[g]] = tmpexp;
                                    }
                                }
                            }
                        }
                    }
                }
            }


            // loop over all elements find the geometry shared ptr and
            // set up basiskey vector
            for(i = 0; i < fielddef.size(); ++i)
            {
                cnt = 0;
                std::vector<std::string>  fields = fielddef[i]->m_fields;
                std::vector<unsigned int> nmodes = fielddef[i]->m_numModes;
                std::vector<LibUtilities::BasisType> basis = fielddef[i]->m_basis;
                bool UniOrder =  fielddef[i]->m_uniOrder;

                for(j = 0; j < fielddef[i]->m_elementIDs.size(); ++j)
                {
                    LibUtilities::BasisKeyVector bkeyvec;
                    id = fielddef[i]->m_elementIDs[j];

                    switch(fielddef[i]->m_shapeType)
                    {
                    case LibUtilities::eSegment:
                    {
                        k = fielddef[i]->m_elementIDs[j];
                        ASSERTL0(m_segGeoms.find(k) != m_segGeoms.end(),
                                "Failed to find geometry with same global id.");
                        geom = m_segGeoms[k];

                        const LibUtilities::PointsKey pkey(nmodes[cnt], pointstype[i][0]);
                        LibUtilities::BasisKey bkey(basis[0], nmodes[cnt], pkey);
                        if(!UniOrder)
                        {
                            cnt++;
                        }
                        bkeyvec.push_back(bkey);
                    }
                    break;
                    case  LibUtilities::eTriangle:
                    {
                        k = fielddef[i]->m_elementIDs[j];
                        ASSERTL0(m_triGeoms.find(k) != m_triGeoms.end(),
                                "Failed to find geometry with same global id.");
                        geom = m_triGeoms[k];
                        for(int b = 0; b < 2; ++b)
                        {
                            const LibUtilities::PointsKey pkey(nmodes[cnt+b],pointstype[i][b]);
                            LibUtilities::BasisKey bkey(basis[b],nmodes[cnt+b],pkey);
                            bkeyvec.push_back(bkey);
                        }

                        if(!UniOrder)
                        {
                            cnt += 2;
                        }
                    }
                    break;
                    case  LibUtilities::eQuadrilateral:
                    {
                        k = fielddef[i]->m_elementIDs[j];
                        ASSERTL0(m_quadGeoms.find(k) != m_quadGeoms.end(),
                                "Failed to find geometry with same global id");
                        geom = m_quadGeoms[k];

                        for(int b = 0; b < 2; ++b)
                        {
                            const LibUtilities::PointsKey pkey(nmodes[cnt+b],pointstype[i][b]);
                            LibUtilities::BasisKey bkey(basis[b],nmodes[cnt+b],pkey);
                            bkeyvec.push_back(bkey);
                        }

                        if(!UniOrder)
                        {
                            cnt += 2;
                        }
                    }
                    break;
                    case  LibUtilities::eTetrahedron:
                    {
                        k = fielddef[i]->m_elementIDs[j];
                        ASSERTL0(m_tetGeoms.find(k) != m_tetGeoms.end(),
                                "Failed to find geometry with same global id");
                        geom = m_tetGeoms[k];

                        for(int b = 0; b < 3; ++b)
                        {
                            const LibUtilities::PointsKey pkey(nmodes[cnt+b],pointstype[i][b]);
                            LibUtilities::BasisKey bkey(basis[b],nmodes[cnt+b],pkey);
                            bkeyvec.push_back(bkey);
                        }

                        if(!UniOrder)
                        {
                            cnt += 2;
                        }
                    }
                    break;
                    case  LibUtilities::ePyramid:
                    {
                        k = fielddef[i]->m_elementIDs[j];
                        ASSERTL0(m_pyrGeoms.find(k) != m_pyrGeoms.end(),
                                "Failed to find geometry with same global id");
                        geom = m_pyrGeoms[k];

                        for(int b = 0; b < 3; ++b)
                        {
                            const LibUtilities::PointsKey pkey(nmodes[cnt+b],pointstype[i][b]);
                            LibUtilities::BasisKey bkey(basis[b],nmodes[cnt+b],pkey);
                            bkeyvec.push_back(bkey);
                        }

                        if(!UniOrder)
                        {
                            cnt += 2;
                        }
                    }
                    break;
                    case  LibUtilities::ePrism:
                    {
                        k = fielddef[i]->m_elementIDs[j];
                        ASSERTL0(m_prismGeoms.find(k) != m_prismGeoms.end(),
                                "Failed to find geometry with same global id");
                        geom = m_prismGeoms[k];

                        for(int b = 0; b < 3; ++b)
                        {
                            const LibUtilities::PointsKey pkey(nmodes[cnt+b],pointstype[i][b]);
                            LibUtilities::BasisKey bkey(basis[b],nmodes[cnt+b],pkey);
                            bkeyvec.push_back(bkey);
                        }

                        if(!UniOrder)
                        {
                            cnt += 2;
                        }
                    }
                    break;
                    case  LibUtilities::eHexahedron:
                    {
                        k = fielddef[i]->m_elementIDs[j];
                        ASSERTL0(m_hexGeoms.find(k) != m_hexGeoms.end(),
                                "Failed to find geometry with same global id");
                        geom = m_hexGeoms[k];

                        for(int b = 0; b < 3; ++b)
                        {
                            const LibUtilities::PointsKey pkey(nmodes[cnt+b],pointstype[i][b]);
                            LibUtilities::BasisKey bkey(basis[b],nmodes[cnt+b],pkey);
                            bkeyvec.push_back(bkey);
                        }

                        if(!UniOrder)
                        {
                            cnt += 2;
                        }
                    }
                    break;
                    default:
                        ASSERTL0(false,"Need to set up for pyramid and prism 3D Expansions");
                        break;
                    }

                    for(k = 0; k < fields.size(); ++k)
                    {
                        expansionMap = m_expansionMapShPtrMap.find(fields[k])->second;
                        if((*expansionMap).find(id) != (*expansionMap).end())
                        {
                            (*expansionMap)[id]->m_geomShPtr = geom;
                            (*expansionMap)[id]->m_basisKeyVector = bkeyvec;
                        }
                    }
                }
            }
        }

        /**
         * \brief Reset all points keys to have equispaced points with
         * optional arguemt of \a npoints which redefines how many
         * points are to be used.
         */
        void MeshGraph::SetExpansionsToEvenlySpacedPoints(int npoints)
        {
            ExpansionMapShPtrMapIter   it;

            // iterate over all defined expansions
            for(it = m_expansionMapShPtrMap.begin(); it != m_expansionMapShPtrMap.end(); ++it)
            {
                ExpansionMapIter expIt;

                for(expIt = it->second->begin(); expIt != it->second->end(); ++expIt)
                {
                    for(int i = 0; i < expIt->second->m_basisKeyVector.size(); ++i)
                    {
                        LibUtilities::BasisKey  bkeyold = expIt->second->m_basisKeyVector[i];

                        int npts;

                        if(npoints) // use input
                        {
                            npts = npoints;
                        }
                        else
                        {
                            npts = bkeyold.GetNumModes();
                        }


                        const LibUtilities::PointsKey pkey(npts,LibUtilities::ePolyEvenlySpaced);
                        LibUtilities::BasisKey bkeynew(bkeyold.GetBasisType(),bkeyold.GetNumModes(), pkey);
                        expIt->second->m_basisKeyVector[i] = bkeynew;

                    }
                }
            }
        }

        /** 
         * \brief Reset all points keys to have expansion order of \a
         *  nmodes.  we keep the point distribution the same and make
         *  the number of points the same difference from the number
         *  of modes as the original expansion definition
         */
        void MeshGraph::SetExpansionsToPolyOrder(int nmodes)
        {
            ExpansionMapShPtrMapIter   it;

            // iterate over all defined expansions
            for(it = m_expansionMapShPtrMap.begin(); it != m_expansionMapShPtrMap.end(); ++it)
            {
                ExpansionMapIter expIt;
                
                for(expIt = it->second->begin(); expIt != it->second->end(); ++expIt)
                {
                    for(int i = 0; i < expIt->second->m_basisKeyVector.size(); ++i)
                    {
                        LibUtilities::BasisKey  bkeyold = expIt->second->m_basisKeyVector[i]; 
                        
                        int npts = nmodes + (bkeyold.GetNumPoints() - bkeyold.GetNumModes());
                        
                        const LibUtilities::PointsKey pkey(npts,bkeyold.GetPointsType());
                        LibUtilities::BasisKey bkeynew(bkeyold.GetBasisType(),nmodes, pkey);
                        expIt->second->m_basisKeyVector[i] = bkeynew; 
                        
                    }
                }
            }
        }


        /**
         * \brief Reset all points keys to have expansion order of \a
         *  nmodes.  we keep the point distribution the same and make
         *  the number of points the same difference from the number
         *  of modes as the original expansion definition
         */
        void MeshGraph::SetExpansionsToPointOrder(int npts)
        {
            ExpansionMapShPtrMapIter   it;

            // iterate over all defined expansions
            for (it = m_expansionMapShPtrMap.begin();
                 it != m_expansionMapShPtrMap.end();
                 ++it)
            {
                ExpansionMapIter expIt;

                for (expIt = it->second->begin();
                     expIt != it->second->end();
                     ++expIt)
                {
                    for(int i = 0;
                        i < expIt->second->m_basisKeyVector.size();
                        ++i)
                    {
                        LibUtilities::BasisKey  bkeyold =
                            expIt->second->m_basisKeyVector[i];

                        const LibUtilities::PointsKey pkey(
                            npts, bkeyold.GetPointsType());

                        LibUtilities::BasisKey bkeynew(bkeyold.GetBasisType(),
                                                       bkeyold.GetNumModes(),
                                                       pkey);
                        expIt->second->m_basisKeyVector[i] = bkeynew;
                    }
                }
            }
        }


        /**
         * For each element of shape given by \a shape in field \a
         * var, replace the current BasisKeyVector describing the
         * expansion in each dimension, with the one provided by \a
         * keys.
         *
         * @TODO: Allow selection of elements through a CompositeVector,
         * as well as by type.
         *
         * @param   shape     The shape of elements to be changed.
         * @param   keys      The new basis vector to apply to those elements.
         */
        void MeshGraph::SetBasisKey(LibUtilities::ShapeType   shape,
            LibUtilities::BasisKeyVector    &keys,
            std::string                     var)
        {
            ExpansionMapIter elemIter;

            ExpansionMapShPtr expansionMap = m_expansionMapShPtrMap.find(var)->second;

            for (elemIter = expansionMap->begin(); elemIter != expansionMap->end(); ++elemIter)
            {
                if ((elemIter->second)->m_geomShPtr->GetShapeType() == shape)
                {
                    (elemIter->second)->m_basisKeyVector = keys;
                }
            }
        }


        /**
         *
         */
        LibUtilities::BasisKeyVector MeshGraph::DefineBasisKeyFromExpansionType(
            GeometrySharedPtr   in,
            ExpansionType       type,
            const int           nummodes)
        {
            LibUtilities::BasisKeyVector returnval;

            LibUtilities::ShapeType shape= in->GetShapeType();

            int quadoffset = 1;
            switch(type)
            {
            case eModified:
            case eModifiedGLLRadau10:
                quadoffset = 1;
                break;
            case eModifiedQuadPlus1:
                quadoffset = 2;
                break;
            case eModifiedQuadPlus2:
                quadoffset = 3;
                break;
            default:
                break;
            }

            switch(type)
            {
            case eModified:
            case eModifiedQuadPlus1:
            case eModifiedQuadPlus2:
            case eModifiedGLLRadau10:
                {
                    switch (shape)
                    {
                    case LibUtilities::eSegment:
                        {
                            const LibUtilities::PointsKey pkey(nummodes+quadoffset, LibUtilities::eGaussLobattoLegendre);
                            LibUtilities::BasisKey bkey(LibUtilities::eModified_A, nummodes, pkey);
                            returnval.push_back(bkey);
                        }
                        break;
                    case LibUtilities::eQuadrilateral:
                        {
                            const LibUtilities::PointsKey pkey(nummodes+quadoffset, LibUtilities::eGaussLobattoLegendre);
                            LibUtilities::BasisKey bkey(LibUtilities::eModified_A, nummodes, pkey);
                            returnval.push_back(bkey);
                            returnval.push_back(bkey);
                        }
                        break;
                    case LibUtilities::eHexahedron:
                        {
                            const LibUtilities::PointsKey pkey(nummodes+quadoffset, LibUtilities::eGaussLobattoLegendre);
                            LibUtilities::BasisKey bkey(LibUtilities::eModified_A, nummodes, pkey);
                            returnval.push_back(bkey);
                            returnval.push_back(bkey);
                            returnval.push_back(bkey);
                        }
                        break;
                    case LibUtilities::eTriangle:
                        {
                            const LibUtilities::PointsKey pkey(nummodes+quadoffset, LibUtilities::eGaussLobattoLegendre);
                            LibUtilities::BasisKey bkey(LibUtilities::eModified_A, nummodes, pkey);
                            returnval.push_back(bkey);

                            const LibUtilities::PointsKey pkey1(nummodes+quadoffset-1, LibUtilities::eGaussRadauMAlpha1Beta0);
                            LibUtilities::BasisKey bkey1(LibUtilities::eModified_B, nummodes, pkey1);

                            returnval.push_back(bkey1);
                        }
                        break;
                    case LibUtilities::eTetrahedron:
                        {
                            const LibUtilities::PointsKey pkey(nummodes+quadoffset, LibUtilities::eGaussLobattoLegendre);
                            LibUtilities::BasisKey bkey(LibUtilities::eModified_A, nummodes, pkey);
                            returnval.push_back(bkey);

                            const LibUtilities::PointsKey pkey1(nummodes+quadoffset-1, LibUtilities::eGaussRadauMAlpha1Beta0);
                            LibUtilities::BasisKey bkey1(LibUtilities::eModified_B, nummodes, pkey1);
                            returnval.push_back(bkey1);

                            if(type == eModifiedGLLRadau10)
                            {
                                const LibUtilities::PointsKey pkey2(nummodes+quadoffset-1, LibUtilities::eGaussRadauMAlpha1Beta0);
                                LibUtilities::BasisKey bkey2(LibUtilities::eModified_C, nummodes, pkey2); 
                                returnval.push_back(bkey2);
                            }
                            else
                            {
                                const LibUtilities::PointsKey pkey2(nummodes+quadoffset-1, LibUtilities::eGaussRadauMAlpha2Beta0);
                                LibUtilities::BasisKey bkey2(LibUtilities::eModified_C, nummodes, pkey2);
                                returnval.push_back(bkey2);
                            }
                        }
                        break;
                    case LibUtilities::ePyramid:
                        {
                            const LibUtilities::PointsKey pkey(nummodes+quadoffset, LibUtilities::eGaussLobattoLegendre);
                            LibUtilities::BasisKey bkey(LibUtilities::eModified_A, nummodes, pkey);
                            returnval.push_back(bkey);
                            returnval.push_back(bkey);

                            const LibUtilities::PointsKey pkey1(nummodes+quadoffset-1, LibUtilities::eGaussRadauMAlpha2Beta0);
                            LibUtilities::BasisKey bkey1(LibUtilities::eModified_C, nummodes, pkey1);
                            returnval.push_back(bkey1);
                        }
                        break;
                    case LibUtilities::ePrism:
                        {
                            const LibUtilities::PointsKey pkey(nummodes+quadoffset, LibUtilities::eGaussLobattoLegendre);
                            LibUtilities::BasisKey bkey(LibUtilities::eModified_A, nummodes, pkey);
                            returnval.push_back(bkey);
                            returnval.push_back(bkey);

                            const LibUtilities::PointsKey pkey1(nummodes+quadoffset-1, LibUtilities::eGaussRadauMAlpha1Beta0);
                            LibUtilities::BasisKey bkey1(LibUtilities::eModified_B, nummodes, pkey1);
                            returnval.push_back(bkey1);

                        }
                        break;
                    default:
                        {
                            ASSERTL0(false,"Expansion not defined in switch for this shape");
                        }
                        break;
                    }
                }
                break;

            case eGLL_Lagrange:
                {
                    switch(shape)
                    {
                    case LibUtilities::eSegment:
                        {
                            const LibUtilities::PointsKey pkey(nummodes+1, LibUtilities::eGaussLobattoLegendre);
                            LibUtilities::BasisKey bkey(LibUtilities::eGLL_Lagrange, nummodes, pkey);
                        returnval.push_back(bkey);
                        }
                        break;
                    case LibUtilities::eQuadrilateral:
                        {
                            const LibUtilities::PointsKey pkey(nummodes+1, LibUtilities::eGaussLobattoLegendre);
                            LibUtilities::BasisKey bkey(LibUtilities::eGLL_Lagrange, nummodes, pkey);
                            returnval.push_back(bkey);
                            returnval.push_back(bkey);
                        }
                        break;
                    case LibUtilities::eTriangle: // define with corrects points key
                        // and change to Ortho on construction
                        {
                            const LibUtilities::PointsKey pkey(nummodes+1, LibUtilities::eGaussLobattoLegendre);
                            LibUtilities::BasisKey bkey(LibUtilities::eGLL_Lagrange, nummodes, pkey);
                            returnval.push_back(bkey);

                            const LibUtilities::PointsKey pkey1(nummodes, LibUtilities::eGaussRadauMAlpha1Beta0);
                            LibUtilities::BasisKey bkey1(LibUtilities::eOrtho_B, nummodes, pkey1);
                            returnval.push_back(bkey1);
                        }
                        break;
                    case LibUtilities::eHexahedron:
                        {
                            const LibUtilities::PointsKey pkey(nummodes+1,LibUtilities::eGaussLobattoLegendre);
                            LibUtilities::BasisKey bkey(LibUtilities::eGLL_Lagrange, nummodes, pkey);

                            returnval.push_back(bkey);
                            returnval.push_back(bkey);
                            returnval.push_back(bkey);
                        }
                        break;
                    default:
                        {
                            ASSERTL0(false, "Expansion not defined in switch  for this shape");
                        }
                        break;
                    }
                }
                break;

            case eGauss_Lagrange:
                {
                    switch (shape)
                    {
                    case LibUtilities::eSegment:
                        {
                            const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eGaussGaussLegendre);
                            LibUtilities::BasisKey bkey(LibUtilities::eGauss_Lagrange, nummodes, pkey);

                            returnval.push_back(bkey);
                        }
                        break;
                    case LibUtilities::eQuadrilateral:
                        {
                            const LibUtilities::PointsKey pkey(nummodes,LibUtilities::eGaussGaussLegendre);
                            LibUtilities::BasisKey bkey(LibUtilities::eGauss_Lagrange, nummodes, pkey);

                            returnval.push_back(bkey);
                            returnval.push_back(bkey);
                        }
                        break;
                    case LibUtilities::eHexahedron:
                        {
                            const LibUtilities::PointsKey pkey(nummodes,LibUtilities::eGaussGaussLegendre);
                            LibUtilities::BasisKey bkey(LibUtilities::eGauss_Lagrange, nummodes, pkey);

                            returnval.push_back(bkey);
                            returnval.push_back(bkey);
                            returnval.push_back(bkey);
                        }
                        break;
                    default:
                        {
                            ASSERTL0(false, "Expansion not defined in switch  for this shape");
                        }
                        break;
                    }
                }
                break;

            case eOrthogonal:
                {
                    switch (shape)
                    {
                    case LibUtilities::eSegment:
                        {
                            const LibUtilities::PointsKey pkey(nummodes+1, LibUtilities::eGaussLobattoLegendre);
                            LibUtilities::BasisKey bkey(LibUtilities::eOrtho_A, nummodes, pkey);

                            returnval.push_back(bkey);
                        }
                        break;
                    case LibUtilities::eTriangle:
                        {
                            const LibUtilities::PointsKey pkey(nummodes+1, LibUtilities::eGaussLobattoLegendre);
                            LibUtilities::BasisKey bkey(LibUtilities::eOrtho_A, nummodes, pkey);

                            returnval.push_back(bkey);

                            const LibUtilities::PointsKey pkey1(nummodes, LibUtilities::eGaussRadauMAlpha1Beta0);
                            LibUtilities::BasisKey bkey1(LibUtilities::eOrtho_B, nummodes, pkey1);

                            returnval.push_back(bkey1);
                        }
                        break;
                    case LibUtilities::eQuadrilateral:
                        {
                            const LibUtilities::PointsKey pkey(nummodes+1, LibUtilities::eGaussLobattoLegendre);
                            LibUtilities::BasisKey bkey(LibUtilities::eOrtho_A, nummodes, pkey);

                            returnval.push_back(bkey);
                            returnval.push_back(bkey);
                        }
                        break;
                    case LibUtilities::eTetrahedron:
                        {
                            const LibUtilities::PointsKey pkey(nummodes+1, LibUtilities::eGaussLobattoLegendre);
                            LibUtilities::BasisKey bkey(LibUtilities::eOrtho_A, nummodes, pkey);

                            returnval.push_back(bkey);

                            const LibUtilities::PointsKey pkey1(nummodes, LibUtilities::eGaussRadauMAlpha1Beta0);
                            LibUtilities::BasisKey bkey1(LibUtilities::eOrtho_B, nummodes, pkey1);

                            returnval.push_back(bkey1);

                            const LibUtilities::PointsKey pkey2(nummodes, LibUtilities::eGaussRadauMAlpha2Beta0);
                            LibUtilities::BasisKey bkey2(LibUtilities::eOrtho_C, nummodes, pkey2);
                        }
                        break;
                    default:
                        {
                    ASSERTL0(false,"Expansion not defined in switch  for this shape");
                        }
                        break;
                    }
                }
                break;

            case eGLL_Lagrange_SEM:
                {
                    switch (shape)
                    {
                    case LibUtilities::eSegment:
                        {
                            const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eGLL_Lagrange, nummodes, pkey);

                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eGLL_Lagrange, nummodes, pkey);

                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eHexahedron:
                {
                    const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eGLL_Lagrange, nummodes, pkey);

                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                default:
                {
                    ASSERTL0(false,"Expansion not defined in switch  for this shape");
                }
                break;
                }
            }
            break;


            case eFourier:
            {
                switch (shape)
                {
                case LibUtilities::eSegment:
                {
                    const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eFourierEvenlySpaced);
                    LibUtilities::BasisKey bkey(LibUtilities::eFourier, nummodes, pkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eFourierEvenlySpaced);
                    LibUtilities::BasisKey bkey(LibUtilities::eFourier, nummodes, pkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eHexahedron:
                {
                    const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eFourierEvenlySpaced);
                    LibUtilities::BasisKey bkey(LibUtilities::eFourier, nummodes, pkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                default:
                {
                    ASSERTL0(false,"Expansion not defined in switch  for this shape");
                }
                break;
                }
            }
            break;


            case eFourierSingleMode:
            {
                switch (shape)
                {
                    case LibUtilities::eSegment:
                    {
                        const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eFourierSingleModeSpaced);
                        LibUtilities::BasisKey bkey(LibUtilities::eFourierSingleMode, nummodes, pkey);
                        returnval.push_back(bkey);
                    }
                    break;
                    case LibUtilities::eQuadrilateral:
                    {
                        const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eFourierSingleModeSpaced);
                        LibUtilities::BasisKey bkey(LibUtilities::eFourierSingleMode, nummodes, pkey);
                        returnval.push_back(bkey);
                        returnval.push_back(bkey);
                    }
                    break;
                    case LibUtilities::eHexahedron:
                    {
                        const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eFourierSingleModeSpaced);
                        LibUtilities::BasisKey bkey(LibUtilities::eFourierSingleMode, nummodes, pkey);
                        returnval.push_back(bkey);
                        returnval.push_back(bkey);
                        returnval.push_back(bkey);
                    }
                    break;
                    default:
                    {
                        ASSERTL0(false,"Expansion not defined in switch  for this shape");
                    }
                    break;
                }
            }
            break;

            case eFourierHalfModeRe:
            {
                switch (shape)
                {
                    case LibUtilities::eSegment:
                    {
                        const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eFourierSingleModeSpaced);
                        LibUtilities::BasisKey bkey(LibUtilities::eFourierHalfModeRe, nummodes, pkey);
                        returnval.push_back(bkey);
                    }
                    break;
                    case LibUtilities::eQuadrilateral:
                    {
                        const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eFourierSingleModeSpaced);
                        LibUtilities::BasisKey bkey(LibUtilities::eFourierHalfModeRe, nummodes, pkey);
                        returnval.push_back(bkey);
                        returnval.push_back(bkey);
                    }
                    break;
                    case LibUtilities::eHexahedron:
                    {
                        const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eFourierSingleModeSpaced);
                        LibUtilities::BasisKey bkey(LibUtilities::eFourierHalfModeRe, nummodes, pkey);
                        returnval.push_back(bkey);
                        returnval.push_back(bkey);
                        returnval.push_back(bkey);
                    }
                        break;
                    default:
                    {
                        ASSERTL0(false,"Expansion not defined in switch  for this shape");
                    }
                    break;
                }
            }
            break;

            case eFourierHalfModeIm:
            {
                switch (shape)
                {
                    case LibUtilities::eSegment:
                    {
                        const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eFourierSingleModeSpaced);
                        LibUtilities::BasisKey bkey(LibUtilities::eFourierHalfModeIm, nummodes, pkey);
                        returnval.push_back(bkey);
                    }
                    break;
                    case LibUtilities::eQuadrilateral:
                    {
                        const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eFourierSingleModeSpaced);
                        LibUtilities::BasisKey bkey(LibUtilities::eFourierHalfModeIm, nummodes, pkey);
                        returnval.push_back(bkey);
                        returnval.push_back(bkey);
                    }
                    break;
                    case LibUtilities::eHexahedron:
                    {
                        const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eFourierSingleModeSpaced);
                        LibUtilities::BasisKey bkey(LibUtilities::eFourierHalfModeIm, nummodes, pkey);
                        returnval.push_back(bkey);
                        returnval.push_back(bkey);
                        returnval.push_back(bkey);
                    }
                    break;
                    default:
                    {
                        ASSERTL0(false,"Expansion not defined in switch  for this shape");
                    }
                    break;
                }
            }
            break;

            case eChebyshev:
            {
                switch (shape)
                {
                case LibUtilities::eSegment:
                {
                    const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eGaussGaussChebyshev);
                    LibUtilities::BasisKey bkey(LibUtilities::eChebyshev, nummodes, pkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eGaussGaussChebyshev);
                    LibUtilities::BasisKey bkey(LibUtilities::eChebyshev, nummodes, pkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eHexahedron:
                {
                    const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eGaussGaussChebyshev);
                    LibUtilities::BasisKey bkey(LibUtilities::eChebyshev, nummodes, pkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                default:
                {
                    ASSERTL0(false,"Expansion not defined in switch  for this shape");
                }
                break;
                }
            }
            break;

            case eFourierChebyshev:
            {
                switch (shape)
                {
                case LibUtilities::eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eFourierEvenlySpaced);
                    LibUtilities::BasisKey bkey(LibUtilities::eFourier, nummodes, pkey);
                    returnval.push_back(bkey);

                    const LibUtilities::PointsKey pkey1(nummodes, LibUtilities::eGaussGaussChebyshev);
                    LibUtilities::BasisKey bkey1(LibUtilities::eChebyshev, nummodes, pkey1);
                    returnval.push_back(bkey1);
                }
                break;
                default:
                {
                    ASSERTL0(false,"Expansion not defined in switch  for this shape");
                }
                break;
                }
            }
            break;

            case eChebyshevFourier:
            {
                switch (shape)
                {
                case LibUtilities::eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey1(nummodes, LibUtilities::eGaussGaussChebyshev);
                    LibUtilities::BasisKey bkey1(LibUtilities::eChebyshev, nummodes, pkey1);
                    returnval.push_back(bkey1);

                    const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eFourierEvenlySpaced);
                    LibUtilities::BasisKey bkey(LibUtilities::eFourier, nummodes, pkey);
                    returnval.push_back(bkey);
                }
                break;
                default:
                {
                    ASSERTL0(false,"Expansion not defined in switch  for this shape");
                }
                break;
                }
            }
            break;

            case eFourierModified:
            {
                switch (shape)
                {
                case LibUtilities::eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey(nummodes, LibUtilities::eFourierEvenlySpaced);
                    LibUtilities::BasisKey bkey(LibUtilities::eFourier, nummodes, pkey);
                    returnval.push_back(bkey);

                    const LibUtilities::PointsKey pkey1(nummodes+1, LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey1(LibUtilities::eModified_A, nummodes, pkey1);
                    returnval.push_back(bkey1);
                }
                break;
                default:
                {
                    ASSERTL0(false,"Expansion not defined in switch  for this shape");
                }
                break;
                }
            }
            break;

            default:
            {
                ASSERTL0(false,"Expansion type not defined");
            }
            break;

            }

            return returnval;
        }

        /**
         *
         */
        LibUtilities::BasisKeyVector
                            MeshGraph::DefineBasisKeyFromExpansionTypeHomo(
                GeometrySharedPtr in,
                ExpansionType type_x,
                ExpansionType type_y,
                ExpansionType type_z,
                const int nummodes_x,
                const int nummodes_y,
                const int nummodes_z)
        {
            LibUtilities::BasisKeyVector returnval;

            LibUtilities::ShapeType shape = in->GetShapeType();

            switch (shape)
            {
                case LibUtilities::eSegment:
                {
                    ASSERTL0(false,"Homogeneous expansion not defined for this shape");
                }
                break;

                case LibUtilities::eQuadrilateral:
                {
                    ASSERTL0(false,"Homogeneous expansion not defined for this shape");
                }
                break;

                case LibUtilities::eHexahedron:
                {
                    switch(type_x)
                    {
                        case eFourier:
                        {
                            const LibUtilities::PointsKey pkey1(nummodes_x,LibUtilities::eFourierEvenlySpaced);
                            LibUtilities::BasisKey bkey1(LibUtilities::eFourier,nummodes_x,pkey1);
                            returnval.push_back(bkey1);
                        }
                            break;

						case eFourierSingleMode:
                        {
                            const LibUtilities::PointsKey pkey1(nummodes_x,LibUtilities::eFourierSingleModeSpaced);
                            LibUtilities::BasisKey bkey1(LibUtilities::eFourierSingleMode,nummodes_x,pkey1);
                            returnval.push_back(bkey1);
                        }
                            break;

						case eFourierHalfModeRe:
                        {
                            const LibUtilities::PointsKey pkey1(nummodes_x,LibUtilities::eFourierSingleModeSpaced);
                            LibUtilities::BasisKey bkey1(LibUtilities::eFourierHalfModeRe,nummodes_x,pkey1);
                            returnval.push_back(bkey1);
                        }
                            break;

						case eFourierHalfModeIm:
                        {
                            const LibUtilities::PointsKey pkey1(nummodes_x,LibUtilities::eFourierSingleModeSpaced);
                            LibUtilities::BasisKey bkey1(LibUtilities::eFourierHalfModeIm,nummodes_x,pkey1);
                            returnval.push_back(bkey1);
                        }
                            break;


                        case eChebyshev:
                        {
                            const LibUtilities::PointsKey pkey1(nummodes_x,LibUtilities::eGaussGaussChebyshev);
                            LibUtilities::BasisKey bkey1(LibUtilities::eChebyshev,nummodes_x,pkey1);
                            returnval.push_back(bkey1);
                        }
                            break;



                        default:
                        {
                            ASSERTL0(false,"Homogeneous expansion can be of Fourier or Chebyshev type only");
                        }
                            break;
                    }


                    switch(type_y)
                    {
                        case eFourier:
                        {
                            const LibUtilities::PointsKey pkey2(nummodes_y,LibUtilities::eFourierEvenlySpaced);
                            LibUtilities::BasisKey bkey2(LibUtilities::eFourier,nummodes_y,pkey2);
                            returnval.push_back(bkey2);
                        }
                            break;


						case eFourierSingleMode:
                        {
                            const LibUtilities::PointsKey pkey2(nummodes_y,LibUtilities::eFourierSingleModeSpaced);
                            LibUtilities::BasisKey bkey2(LibUtilities::eFourierSingleMode,nummodes_y,pkey2);
                            returnval.push_back(bkey2);
                        }
                            break;

						case eFourierHalfModeRe:
                        {
                            const LibUtilities::PointsKey pkey2(nummodes_y,LibUtilities::eFourierSingleModeSpaced);
                            LibUtilities::BasisKey bkey2(LibUtilities::eFourierHalfModeRe,nummodes_y,pkey2);
                            returnval.push_back(bkey2);
                        }
                            break;

						case eFourierHalfModeIm:
                        {
                            const LibUtilities::PointsKey pkey2(nummodes_y,LibUtilities::eFourierSingleModeSpaced);
                            LibUtilities::BasisKey bkey2(LibUtilities::eFourierHalfModeIm,nummodes_y,pkey2);
                            returnval.push_back(bkey2);
                        }
                            break;

                        case eChebyshev:
                        {
                            const LibUtilities::PointsKey pkey2(nummodes_y,LibUtilities::eGaussGaussChebyshev);
                            LibUtilities::BasisKey bkey2(LibUtilities::eChebyshev,nummodes_y,pkey2);
                            returnval.push_back(bkey2);
                        }
                            break;

                        default:
                        {
                            ASSERTL0(false,"Homogeneous expansion can be of Fourier or Chebyshev type only");
                        }
                            break;
                    }

                    switch(type_z)
                    {
                        case eFourier:
                        {
                            const LibUtilities::PointsKey pkey3(nummodes_z,LibUtilities::eFourierEvenlySpaced);
                            LibUtilities::BasisKey bkey3(LibUtilities::eFourier,nummodes_z,pkey3);
                            returnval.push_back(bkey3);
                        }
                            break;

						case eFourierSingleMode:
                        {
                            const LibUtilities::PointsKey pkey3(nummodes_z,LibUtilities::eFourierSingleModeSpaced);
                            LibUtilities::BasisKey bkey3(LibUtilities::eFourierSingleMode,nummodes_z,pkey3);
                            returnval.push_back(bkey3);
                        }
                            break;

						case eFourierHalfModeRe:
                        {
                            const LibUtilities::PointsKey pkey3(nummodes_z,LibUtilities::eFourierSingleModeSpaced);
                            LibUtilities::BasisKey bkey3(LibUtilities::eFourierHalfModeRe,nummodes_z,pkey3);
                            returnval.push_back(bkey3);
                        }
                            break;

						case eFourierHalfModeIm:
                        {
                            const LibUtilities::PointsKey pkey3(nummodes_z,LibUtilities::eFourierSingleModeSpaced);
                            LibUtilities::BasisKey bkey3(LibUtilities::eFourierHalfModeIm,nummodes_z,pkey3);
                            returnval.push_back(bkey3);
                        }
                            break;

                        case eChebyshev:
                        {
                            const LibUtilities::PointsKey pkey3(nummodes_z,LibUtilities::eGaussGaussChebyshev);
                            LibUtilities::BasisKey bkey3(LibUtilities::eChebyshev,nummodes_z,pkey3);
                            returnval.push_back(bkey3);
                        }
                            break;

                        default:
                        {
                            ASSERTL0(false,"Homogeneous expansion can be of Fourier or Chebyshev type only");
                        }
                            break;
                    }
                }
                break;

                case LibUtilities::eTriangle:
                {
                    ASSERTL0(false,"Homogeneous expansion not defined for this shape");
                }
                break;

                case LibUtilities::eTetrahedron:
                {
                    ASSERTL0(false,"Homogeneous expansion not defined for this shape");
                }
                break;

                default:
                    ASSERTL0(false,"Expansion not defined in switch  for this shape");
                break;
            }

            return returnval;
        }


        /**
         *
         */
        PointGeomSharedPtr MeshGraph::AddVertex(NekDouble x, NekDouble y, NekDouble z)
        {
            unsigned int nextId = m_vertSet.rbegin()->first + 1;
            PointGeomSharedPtr vert(MemoryManager<PointGeom>::AllocateSharedPtr(m_spaceDimension, nextId, x, y, z));
            m_vertSet[nextId] = vert;
            return vert;
        }

        /**
         *
         */
        SegGeomSharedPtr MeshGraph::AddEdge(PointGeomSharedPtr v0, PointGeomSharedPtr v1,
                CurveSharedPtr curveDefinition)
        {
            PointGeomSharedPtr vertices[] = {v0, v1};
            SegGeomSharedPtr edge;
            int edgeId = m_segGeoms.rbegin()->first + 1;

            if( curveDefinition )
            {
                edge = MemoryManager<SegGeom>::AllocateSharedPtr(edgeId, m_spaceDimension, vertices, curveDefinition);
            }
            else
            {
                edge = MemoryManager<SegGeom>::AllocateSharedPtr(edgeId, m_spaceDimension, vertices);
            }
            m_segGeoms[edgeId] = edge;
            return edge;
        }


        /**
         *
         */
        TriGeomSharedPtr MeshGraph::AddTriangle(SegGeomSharedPtr edges[], StdRegions::Orientation orient[])
        {
            int indx = m_triGeoms.rbegin()->first + 1;
            TriGeomSharedPtr trigeom(MemoryManager<TriGeom>::AllocateSharedPtr(indx, edges, orient));
            trigeom->SetGlobalID(indx);

            m_triGeoms[indx] = trigeom;

            return trigeom;
        }


        /**
         *
         */
        QuadGeomSharedPtr MeshGraph::AddQuadrilateral(SegGeomSharedPtr edges[], StdRegions::Orientation orient[])
        {
            int indx = m_quadGeoms.rbegin()->first + 1;
            QuadGeomSharedPtr quadgeom(MemoryManager<QuadGeom>::AllocateSharedPtr(indx, edges, orient));
            quadgeom->SetGlobalID(indx);

            m_quadGeoms[indx] = quadgeom;
            return quadgeom;
        }


        /**
         *
         */
        PrismGeomSharedPtr MeshGraph::AddPrism(TriGeomSharedPtr tfaces[PrismGeom::kNtfaces],
                QuadGeomSharedPtr qfaces[PrismGeom::kNqfaces])
        {
            // Setting the orientation is disabled in the reader.  Why?
            Geometry2DSharedPtr faces[] = { qfaces[0], tfaces[0], qfaces[1], tfaces[1], qfaces[2] };
            unsigned int index = m_prismGeoms.rbegin()->first + 1;
            PrismGeomSharedPtr prismgeom(MemoryManager<PrismGeom>::AllocateSharedPtr(faces));
            prismgeom->SetGlobalID(index);

            m_prismGeoms[index] = prismgeom;
            return prismgeom;
        }


        /**
         *
         */
        TetGeomSharedPtr MeshGraph::AddTetrahedron(TriGeomSharedPtr tfaces[TetGeom::kNtfaces])
        {
            unsigned int index = m_tetGeoms.rbegin()->first + 1;
            TetGeomSharedPtr tetgeom(MemoryManager<TetGeom>::AllocateSharedPtr(tfaces));
            tetgeom->SetGlobalID(index);

            m_tetGeoms[index] = tetgeom;
            return tetgeom;
        }


        /**
         *
         */
        PyrGeomSharedPtr MeshGraph::AddPyramid(TriGeomSharedPtr tfaces[PyrGeom::kNtfaces],
                QuadGeomSharedPtr qfaces[PyrGeom::kNqfaces])
        {
            Geometry2DSharedPtr faces[] = { qfaces[0], tfaces[0], tfaces[1], tfaces[2], tfaces[3] };
            unsigned int index = m_pyrGeoms.rbegin()->first + 1;

            PyrGeomSharedPtr pyrgeom(MemoryManager<PyrGeom>::AllocateSharedPtr(faces));
            pyrgeom->SetGlobalID(index);

            m_pyrGeoms[index] = pyrgeom;
            return pyrgeom;
        }


        /**
         *
         */
        HexGeomSharedPtr MeshGraph::AddHexahedron(QuadGeomSharedPtr qfaces[HexGeom::kNqfaces])
        {
            unsigned int index = m_hexGeoms.rbegin()->first + 1;
            HexGeomSharedPtr hexgeom(MemoryManager<HexGeom>::AllocateSharedPtr(qfaces));
            hexgeom->SetGlobalID(index);
            m_hexGeoms[index] = hexgeom;
            return hexgeom;
        }


        /**
         * Generate a single vector of Expansion structs mapping global element
         * ID to a corresponding Geometry shared pointer and basis key.
         *
         * Expansion map ensures elements which appear in multiple composites
         * within the domain are only listed once.
         */
        ExpansionMapShPtr MeshGraph::SetUpExpansionMap(void)
        {
            ExpansionMapShPtr returnval;
            returnval = MemoryManager<ExpansionMap>::AllocateSharedPtr();

            for(int d = 0; d < m_domain.size(); ++d)
            {
                CompositeMap::const_iterator compIter;

                for (compIter = m_domain[d].begin(); compIter != m_domain[d].end(); ++compIter)
                {
                    GeometryVector::const_iterator x;
                    for (x = compIter->second->begin(); x != compIter->second->end(); ++x)
                    {
                        LibUtilities::BasisKeyVector def;
                        ExpansionShPtr expansionElementShPtr =
                            MemoryManager<Expansion>::AllocateSharedPtr(*x, def);
                        int id = (*x)->GetGlobalID();
                        (*returnval)[id] = expansionElementShPtr;
                    }
                }
            }

            return returnval;
        }
    }; //end of namespace
}; //end of namespace
