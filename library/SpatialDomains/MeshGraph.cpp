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
#include <cmath>

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
            m_meshDimension(3),
            m_spaceDimension(3)
        {
        }

        MeshGraph::MeshGraph(unsigned int meshDimension, unsigned int spaceDimension) :
            m_meshDimension(meshDimension),
            m_spaceDimension(spaceDimension)
        {
        }

        boost::shared_ptr<MeshGraph> MeshGraph::Read(const std::string& infilename)
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
                returnval->ReadExpansions(infilename);
            }
            return returnval;
        }

        void MeshGraph::SetExpansions(std::vector<SpatialDomains::FieldDefinitionsSharedPtr> &fielddef)
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
                        
                        // check to see if DefaultVar also not set and if so assign it to this expansion
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
                bool pointDef = fielddef[i]->m_pointsDef;
                bool numPointDef = fielddef[i]->m_numPointsDef;

                // Check points and numpoints
                std::vector<unsigned int> npoints = fielddef[i]->m_numPoints;
                std::vector<LibUtilities::PointsType> points = fielddef[i]->m_points;

                bool UniOrder =  fielddef[i]->m_uniOrder;

                int check = 0;
                for(j=0; j< basis.size(); ++j)
                {
                    if(LibUtilities::BasisTypeMap[basis[j]]=="Modified_A" ||
                        LibUtilities::BasisTypeMap[basis[j]]=="Modified_B")
                        check++;
                }

                if(check==basis.size())
                {
                    for(j = 0; j < fielddef[i]->m_elementIDs.size(); ++j)
                    {

                        LibUtilities::BasisKeyVector bkeyvec;
                        id = fielddef[i]->m_elementIDs[j];
                        
                        switch(fielddef[i]->m_shapeType)
                        {
                        case eSegment:
                            {
                                for(k = 0; k < m_segGeoms.size();++k)
                                {
                                    if(m_segGeoms[k]->GetGlobalID() == fielddef[i]->m_elementIDs[j])
                                    {
                                        geom = m_segGeoms[k];
                                        break;
                                    }
                                }
                                ASSERTL0(k != m_segGeoms.size(),"Failed to find geometry with same global id");

                                LibUtilities::PointsKey pkey(nmodes[cnt]+1,LibUtilities::eGaussLobattoLegendre);

                                if(numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey1(npoints[cnt],points[0]);
                                    pkey = pkey1;
                                }
                                else if(!numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey1(nmodes[cnt]+1,points[0]);
                                    pkey = pkey1;
                                }
                                else if(numPointDef&&!pointDef)
                                {
                                    const LibUtilities::PointsKey pkey1(npoints[cnt],LibUtilities::eGaussLobattoLegendre);
                                    pkey = pkey1;
                                }

                                LibUtilities::BasisKey bkey(basis[0],nmodes[cnt],pkey);

                                if(!UniOrder)
                                {
                                    cnt++;
                                }
                                bkeyvec.push_back(bkey);
                            }
                            break;
                        case eTriangle:
                            {
                                for(k = 0; k < m_triGeoms.size();++k)
                                {
                                    if(m_triGeoms[k]->GetGlobalID() == fielddef[i]->m_elementIDs[j])
                                    {
                                        geom = m_triGeoms[k];
                                        break;
                                    }
                                }
                                ASSERTL0(k != m_triGeoms.size(),"Failed to find geometry with same global id");

                                LibUtilities::PointsKey pkey(nmodes[cnt]+1,LibUtilities::eGaussLobattoLegendre);
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

                                LibUtilities::PointsKey pkey1(nmodes[cnt+1],LibUtilities::eGaussRadauMAlpha1Beta0);
                                if(numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(npoints[cnt+1],points[1]);
                                    pkey1 = pkey2;
                                }
                                else if(!numPointDef&&pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(nmodes[cnt+1],points[1]);
                                    pkey1 = pkey2;
                                }
                                else if(numPointDef&&!pointDef)
                                {
                                    const LibUtilities::PointsKey pkey2(npoints[cnt+1],LibUtilities::eGaussRadauMAlpha1Beta0);
                                    pkey1 = pkey2;
                                }
                                LibUtilities::BasisKey bkey1(basis[1],nmodes[cnt+1],pkey1);
                                bkeyvec.push_back(bkey1);

                                if(!UniOrder)
                                {
                                    cnt += 2;
                                }
                            }
                            break;
                        case eQuadrilateral:
                            {
                                for(k = 0; k < m_quadGeoms.size();++k)
                                {
                                    if(m_quadGeoms[k]->GetGlobalID() == fielddef[i]->m_elementIDs[j])
                                    {
                                        geom = m_quadGeoms[k];
                                        break;
                                    }
                                }
                                ASSERTL0(k != m_quadGeoms.size(),"Failed to find geometry with same global id");
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

                                if(!UniOrder)
                                {
                                    cnt += 2;
                                }
                            }
                            break;
                        default:
                            ASSERTL0(false,"Need to set up for 3D Expansions");
                            break;
                        }

                        for(k = 0; k < fields.size(); ++k)
                        {
                            expansionMap = m_expansionMapShPtrMap.find(fields[k])->second;
                            (*expansionMap)[id]->m_geomShPtr = geom;
                            (*expansionMap)[id]->m_basisKeyVector = bkeyvec;
                        }
                    }
                }
                else
                {
                    ASSERTL0(false,"Need to set up for non Modified basis");
                }
            }
        }


        void MeshGraph::SetExpansions(std::vector<SpatialDomains::FieldDefinitionsSharedPtr> &fielddef, std::vector< std::vector<LibUtilities::PointsType> > &pointstype)
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
                        
                        // check to see if DefaultVar also not set and if so assign it to this expansion
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
                    case eSegment:
                        {
                            k = fielddef[i]->m_elementIDs[j];
                            ASSERTL0(m_segGeoms.find(k) != m_segGeoms.end(),
                                     "Failed to find geometry with same global id.");
                            geom = m_segGeoms[k];
//                            for(k = 0; k < m_segGeoms.size();++k)
//                            {
//                                if(m_segGeoms[k]->GetGlobalID() == fielddef[i]->m_elementIDs[j])
//                                {
//                                    geom = m_segGeoms[k];
//                                    break;
//                                }
//                            }
//                            ASSERTL0(k != m_segGeoms.size(),"Failed to find geometry with same global id");
                            const LibUtilities::PointsKey pkey(nmodes[cnt],pointstype[i][0]);
                            LibUtilities::BasisKey bkey(basis[0],nmodes[cnt],pkey);
                            if(!UniOrder)
                            {
                                cnt++;
                            }
                            bkeyvec.push_back(bkey);
                        }
                        break;
                    case eTriangle:
                        {
                            k = fielddef[i]->m_elementIDs[j];
                            ASSERTL0(m_triGeoms.find(k) != m_triGeoms.end(),
                                     "Failed to find geometry with same global id.");
                            geom = m_triGeoms[k];
//                            for(k = 0; k < m_triGeoms.size();++k)
//                            {
//                                if(m_triGeoms[k]->GetGlobalID() == fielddef[i]->m_elementIDs[j])
//                                {
//                                    geom = m_triGeoms[k];
//                                    break;
//                                }
//                            }
//                            ASSERTL0(k != m_triGeoms.size(),"Failed to find geometry with same global id");
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
                    case eQuadrilateral:
                        {
                            k = fielddef[i]->m_elementIDs[j];
                            ASSERTL0(m_quadGeoms.find(k) != m_quadGeoms.end(),
                                    "Failed to find geometry with same global id");
                            geom = m_quadGeoms[k];
//                            for(k = 0; k < m_quadGeoms.size();++k)
//                            {
//                                if(m_quadGeoms[k]->GetGlobalID() == fielddef[i]->m_elementIDs[j])
//                                {
//                                    geom = m_quadGeoms[k];
//                                    break;
//                                }
//                            }
//                            ASSERTL0(k != m_quadGeoms.size(),"Failed to find geometry with same global id");

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
                    case eTetrahedron:
                        {
                            k = fielddef[i]->m_elementIDs[j];
                            ASSERTL0(m_tetGeoms.find(k) != m_tetGeoms.end(),
                                    "Failed to find geometry with same global id");
                            geom = m_tetGeoms[k];

//                            for(k = 0; k < m_tetGeoms.size();++k)
//                            {
//                                if(m_tetGeoms[k]->GetGlobalID() == fielddef[i]->m_elementIDs[j])
//                                {
//                                    geom = m_tetGeoms[k];
//                                    break;
//                                }
//                            }
//                            ASSERTL0(k != m_tetGeoms.size(),"Failed to find geometry with same global id");
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
                    case eHexahedron:
                        {
                            k = fielddef[i]->m_elementIDs[j];
                            ASSERTL0(m_hexGeoms.find(k) != m_hexGeoms.end(),
                                    "Failed to find geometry with same global id");
                            geom = m_hexGeoms[k];
//                            for(k = 0; k < m_quadGeoms.size();++k)
//                            {
//                                if(m_hexGeoms[k]->GetGlobalID() == fielddef[i]->m_elementIDs[j])
//                                {
//                                    geom = m_hexGeoms[k];
//                                    break;
//                                }
//                            }
//                            ASSERTL0(k != m_hexGeoms.size(),"Failed to find geometry with same global id");

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
                            (*expansionMap)[id]->m_geomShPtr = geom;
                            (*expansionMap)[id]->m_basisKeyVector = bkeyvec;
                        }
                }
            }

        }


        // \brief Read will read the meshgraph vertices given a filename.
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
//                ASSERTL0(indx == nextVertexNumber, "Vertex IDs must begin with zero and be sequential.");

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
                NekDouble xval, yval, zval;
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
                            VertexComponentSharedPtr vert(MemoryManager<VertexComponent>::AllocateSharedPtr(m_spaceDimension, indx, xval, yval, zval));
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


        void MeshGraph::WriteGeometry(const std::string& fileName)
        {
            TiXmlDocument doc(fileName);
            WriteGeometry(doc);
        }

        void MeshGraph::WriteGeometry(TiXmlDocument& doc)
        {

        }

        VertexComponentSharedPtr MeshGraph::AddVertex(NekDouble x, NekDouble y, NekDouble z)
        {
            unsigned int nextId = m_vertSet.rbegin()->first + 1;
            VertexComponentSharedPtr vert(MemoryManager<VertexComponent>::AllocateSharedPtr(m_spaceDimension, nextId, x, y, z));
            m_vertSet[nextId] = vert;
            return vert;
        }

        SegGeomSharedPtr MeshGraph::AddEdge(VertexComponentSharedPtr v0, VertexComponentSharedPtr v1,
            CurveSharedPtr curveDefinition)
        {
            VertexComponentSharedPtr vertices[] = {v0, v1};
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

        TriGeomSharedPtr MeshGraph::AddTriangle(SegGeomSharedPtr edges[], StdRegions::EdgeOrientation orient[])
        {
            int indx = m_triGeoms.rbegin()->first + 1;
            TriGeomSharedPtr trigeom(MemoryManager<TriGeom>::AllocateSharedPtr(indx, edges, orient));
            trigeom->SetGlobalID(indx);

            m_triGeoms[indx] = trigeom;

            return trigeom;
        }

        QuadGeomSharedPtr MeshGraph::AddQuadrilateral(SegGeomSharedPtr edges[], StdRegions::EdgeOrientation orient[])
        {
            int indx = m_quadGeoms.rbegin()->first + 1;
            QuadGeomSharedPtr quadgeom(MemoryManager<QuadGeom>::AllocateSharedPtr(indx, edges, orient));
            quadgeom->SetGlobalID(indx);

            m_quadGeoms[indx] = quadgeom;
            return quadgeom;
        }

        PrismGeomSharedPtr MeshGraph::AddPrism(TriGeomSharedPtr tfaces[PrismGeom::kNtfaces],
            QuadGeomSharedPtr qfaces[PrismGeom::kNqfaces])
        {
            // Setting the orientation is disabled in the reader.  Why?
            StdRegions::FaceOrientation faceorient[PrismGeom::kNtfaces + PrismGeom::kNqfaces];

            unsigned int index = m_prismGeoms.rbegin()->first + 1;
            PrismGeomSharedPtr prismgeom(MemoryManager<PrismGeom>::AllocateSharedPtr(tfaces, qfaces, faceorient));
            prismgeom->SetGlobalID(index);

            m_prismGeoms[index] = prismgeom;
            return prismgeom;
        }

        TetGeomSharedPtr MeshGraph::AddTetrahedron(TriGeomSharedPtr tfaces[TetGeom::kNtfaces])
        {
            unsigned int index = m_tetGeoms.rbegin()->first + 1;
            TetGeomSharedPtr tetgeom(MemoryManager<TetGeom>::AllocateSharedPtr(tfaces));
            tetgeom->SetGlobalID(index);

            m_tetGeoms[index] = tetgeom;
            return tetgeom;
        }

        PyrGeomSharedPtr MeshGraph::AddPyramid(TriGeomSharedPtr tfaces[PyrGeom::kNtfaces],
            QuadGeomSharedPtr qfaces[PyrGeom::kNqfaces])
        {
            unsigned int index = m_pyrGeoms.rbegin()->first + 1;

            StdRegions::FaceOrientation faceorient[PyrGeom::kNtfaces + PyrGeom::kNqfaces];
            PyrGeomSharedPtr pyrgeom(MemoryManager<PyrGeom>::AllocateSharedPtr(tfaces, qfaces, faceorient));
            pyrgeom->SetGlobalID(index);

            m_pyrGeoms[index] = pyrgeom;
            return pyrgeom;
        }

        HexGeomSharedPtr MeshGraph::AddHexahedron(QuadGeomSharedPtr qfaces[HexGeom::kNqfaces])
        {
            unsigned int index = m_hexGeoms.rbegin()->first + 1;
            HexGeomSharedPtr hexgeom(MemoryManager<HexGeom>::AllocateSharedPtr(qfaces));
            hexgeom->SetGlobalID(index);
            m_hexGeoms[index] = hexgeom;
            return hexgeom;
        }

        // \brief Read the expansions given the XML file path.
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

        LibUtilities::BasisKeyVector DefineBasisKeyFromExpansionType(GeometrySharedPtr in,
                                     ExpansionType type,
                                     const int order)
        {
            LibUtilities::BasisKeyVector returnval;
            
            GeomShapeType shape= in->GetGeomShapeType();
            
            switch(type)
            {
            case eModified:
                switch (shape)
                {
                case eSegment:
                    {
                        const LibUtilities::PointsKey pkey(order+1,LibUtilities::eGaussLobattoLegendre);
                        LibUtilities::BasisKey bkey(LibUtilities::eModified_A,order,pkey);
                        returnval.push_back(bkey);
                    }
                    break;
                case eQuadrilateral:
                    {
                        const LibUtilities::PointsKey pkey(order+1,LibUtilities::eGaussLobattoLegendre);
                        LibUtilities::BasisKey bkey(LibUtilities::eModified_A,order,pkey);
                        returnval.push_back(bkey);
                        returnval.push_back(bkey);
                    }
                    break;
                case eHexahedron:
                    {
                        const LibUtilities::PointsKey pkey(order+1,LibUtilities::eGaussLobattoLegendre);
                        LibUtilities::BasisKey bkey(LibUtilities::eModified_A,order,pkey);
                        returnval.push_back(bkey);
                        returnval.push_back(bkey);
                        returnval.push_back(bkey);
                    }
                    break;
                case eTriangle:
                    {
                        const LibUtilities::PointsKey pkey(order+1,LibUtilities::eGaussLobattoLegendre);
                        LibUtilities::BasisKey bkey(LibUtilities::eModified_A,order,pkey);
                        returnval.push_back(bkey);
                        
                        const LibUtilities::PointsKey pkey1(order,LibUtilities::eGaussRadauMAlpha1Beta0);
                        LibUtilities::BasisKey bkey1(LibUtilities::eModified_B,order,pkey1);
                        
                        returnval.push_back(bkey1);
                    }
                    break;
                case eTetrahedron:
                    {
                        const LibUtilities::PointsKey pkey(order+1,LibUtilities::eGaussLobattoLegendre);
                        LibUtilities::BasisKey bkey(LibUtilities::eModified_A,order,pkey);
                        returnval.push_back(bkey);
                        
                        const LibUtilities::PointsKey pkey1(order,LibUtilities::eGaussRadauMAlpha1Beta0);
                        LibUtilities::BasisKey bkey1(LibUtilities::eModified_B,order,pkey1);
                        returnval.push_back(bkey1);
                        
                        const LibUtilities::PointsKey pkey2(order,LibUtilities::eGaussRadauMAlpha2Beta0);
                        LibUtilities::BasisKey bkey2(LibUtilities::eModified_C,order,pkey2);
                        returnval.push_back(bkey2);
                    }
                    break;
                default:
                    ASSERTL0(false,"Expansion not defined in switch  for this shape");
                    break;
                }
                break;
                
            case eGLL_Lagrange:
            {
                
                switch(shape)
                {
                case eSegment:
                    {
                        const LibUtilities::PointsKey pkey(order+1,LibUtilities::eGaussLobattoLegendre);
                        LibUtilities::BasisKey bkey(LibUtilities::eGLL_Lagrange,order,pkey);
                        returnval.push_back(bkey);
            }
                    break;
                case eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey(order+1,LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eGLL_Lagrange,order,pkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                case eTriangle: // define with corrects points key
                    // and change to Ortho on
                    // construction
                {

                    const LibUtilities::PointsKey pkey(order+1,LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eGLL_Lagrange,order,pkey);
                    
                    returnval.push_back(bkey);
                    
                    const LibUtilities::PointsKey pkey1(order,LibUtilities::eGaussRadauMAlpha1Beta0);
                    LibUtilities::BasisKey bkey1(LibUtilities::eOrtho_B,order,pkey1);
                    
                    returnval.push_back(bkey1);
                }
                break;
                default:
                    ASSERTL0(false,"Expansion not defined in switch  for this shape");
                    break;
                }
            }
            break;
            case eOrthogonal:
                switch (shape)
                {
                case eSegment:
                    {
                        const LibUtilities::PointsKey pkey(order+1,LibUtilities::eGaussLobattoLegendre);
                        LibUtilities::BasisKey bkey(LibUtilities::eOrtho_A,order,pkey);
                        
                        returnval.push_back(bkey);
                    }
                    break;
                case eTriangle:
                {
                    const LibUtilities::PointsKey pkey(order+1,LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eOrtho_A,order,pkey);
                    
                    returnval.push_back(bkey);
                    
                    const LibUtilities::PointsKey pkey1(order,LibUtilities::eGaussRadauMAlpha1Beta0);
                    LibUtilities::BasisKey bkey1(LibUtilities::eOrtho_B,order,pkey1);
                    
                    returnval.push_back(bkey1);
                }
                break;
                case eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey(order+1,LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eOrtho_A,order,pkey);
                    
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                case eTetrahedron:
                {
                    const LibUtilities::PointsKey pkey(order+1,LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eOrtho_A,order,pkey);
                    
                    returnval.push_back(bkey);
                    
                    const LibUtilities::PointsKey pkey1(order,LibUtilities::eGaussRadauMAlpha1Beta0);
                    LibUtilities::BasisKey bkey1(LibUtilities::eOrtho_B,order,pkey1);
                    
                    returnval.push_back(bkey1);
                    
                    const LibUtilities::PointsKey pkey2(order,LibUtilities::eGaussRadauMAlpha2Beta0);
                    LibUtilities::BasisKey bkey2(LibUtilities::eOrtho_C,order,pkey2);
                }
                break;
                default:
                    ASSERTL0(false,"Expansion not defined in switch  for this shape");
                    break;
                }
                break;
            case eGLL_Lagrange_SEM:
            {
                switch (shape)
                {
                case eSegment:
                {
                    const LibUtilities::PointsKey pkey(order,LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eGLL_Lagrange,order,pkey);
                    
                    returnval.push_back(bkey);
                }
                break;
                case eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey(order,LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eGLL_Lagrange,order,pkey);
                    
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                case eHexahedron:
                {
                    const LibUtilities::PointsKey pkey(order,LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eGLL_Lagrange,order,pkey);
                    
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                default:
                    ASSERTL0(false,"Expansion not defined in switch  for this shape");
                    break;
                }
            }
            break;
            default:
                break;
            }
            
            return returnval;
        }

        ExpansionMapShPtr MeshGraph::SetUpExpansionMap(void)
        {
            ExpansionMapShPtr returnval; 
            returnval = MemoryManager<ExpansionMap>::AllocateSharedPtr();
            
            // Need a vector of all elements and their associated
            // expansion information.
            const CompositeMap &domain = this->GetDomain();
            CompositeMap::const_iterator compIter;
            
            for (compIter = domain.begin(); compIter != domain.end(); ++compIter)
            {
                boost::shared_ptr<GeometryVector> geomVectorShPtr = compIter->second;
                GeometryVectorIter geomIter;
                for (geomIter = geomVectorShPtr->begin(); geomIter != geomVectorShPtr->end(); ++geomIter)
                {
                    // Make sure we only have one instance of the
                    // GeometrySharedPtr stored in the list.
                    ExpansionMap::iterator elemIter;
                    for (elemIter = returnval->begin(); elemIter != returnval->end(); ++elemIter)
                    {
                        if (elemIter->second->m_geomShPtr == *geomIter)
                        {
                            break;
                        }
                    }
                    
                    // Not found in list.
                    if (elemIter == returnval->end())
                    {
                        LibUtilities::BasisKeyVector def;
                        ExpansionShPtr expansionElementShPtr =
                            MemoryManager<Expansion>::AllocateSharedPtr(*geomIter, def);
                        int id = (*geomIter)->GetGlobalID();
                        (*returnval)[id] = expansionElementShPtr;
                    }
                }
            }
            
            return returnval;
        }
        
        // \brief Read the expansions given the XML document reference.
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
                        int           expansion_order;

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

                            const char *nStr = expansion->Attribute("NUMMODES");
                            ASSERTL0(nStr,"NUMMODES was not defined in EXPANSION section of input");
                            std::string nummodesStr = nStr;

                            Equation nummodesEqn(nummodesStr);

                            expansion_order = (int) nummodesEqn.Evaluate();

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
                                ExpansionMapIter expVecIter;
                                for (expVecIter = expansionMap->begin(); expVecIter != expansionMap->end(); ++expVecIter)
                                {
                                    if (*geomVecIter == (expVecIter->second)->m_geomShPtr)
                                    {
                                        if(useExpansionType)
                                        {
                                            (expVecIter->second)->m_basisKeyVector = DefineBasisKeyFromExpansionType(*geomVecIter,expansion_type,expansion_order);
                                        }
                                        else
                                        {
                                            ASSERTL0((*geomVecIter)->GetShapeDim() == basiskeyvec.size()," There is an incompatible expansion dimension with geometry dimension");
                                            (expVecIter->second)->m_basisKeyVector = basiskeyvec;
                                        }
                                        break;
                                    }
                                }
                            }
                        }

                        expansion = expansion->NextSiblingElement("E");
                    }
                }
                else if(expType == "ELEMENTS")  // Reading a file with the expansion definition
                {
                    std::vector<FieldDefinitionsSharedPtr> fielddefs;
                    ImportFieldDefs(doc, fielddefs, true);
                    cout << "    Number of elements: " << fielddefs.size() << endl;
                    SetExpansions(fielddefs);
                }
                else
                {
                    ASSERTL0(false,"Expansion type not defined");
                }
            }
        }

#ifdef NOTNEEDED // SJS: Think this needs depracating
        /**
         * For each element of shape given by \a shape in field \a
         * var, replace the current BasisKeyVector describing the
         * expansion in each dimension, with the one provided by \a
         * keys.
         *
         * @TODO: Allow selection of elements through a CompositeVector, as well as
         * by type.
         *
         * @param   shape           The shape of elements to be changed.
         * @param   keys            The new basis vector to apply to those elements.
         */
        void MeshGraph::SetBasisKey(std::string var,
                                    SpatialDomains::GeomShapeType shape,
                                    LibUtilities::BasisKeyVector &keys)
        {
            ExpansionMapIter elemIter;
            for (elemIter = m_expansions.begin(); elemIter != m_expansions.end(); ++elemIter)
            {
                if ((elemIter->second)->m_geomShPtr->GetGeomShapeType() == shape)
                {
                    (elemIter->second)->m_basisKeyVector = keys;
                }
            }
        }
#endif

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
                    NekDouble xval, yval, zval;
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
                                VertexComponentSharedPtr vert(MemoryManager<VertexComponent>::AllocateSharedPtr(m_meshDimension, edgeindx, xval, yval, zval));

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

                    m_curvedEdges.push_back(curve);

                    edgelement = edgelement->NextSiblingElement("E");

                } // end if-loop

            } // end while-loop


            TiXmlElement *facelement = field->FirstChildElement("F");
            int faceindx, faceid;
            int nextFaceNumber = -1;

            while(facelement)
            {
                /// These should be ordered.
                nextFaceNumber++;

                std::string face(facelement->ValueStr());
                ASSERTL0(face == "F", (std::string("Unknown 3D curve type: ") + face).c_str());

                /// Read id attribute.
                err = facelement->QueryIntAttribute("ID", &faceindx);

                ASSERTL0(err == TIXML_SUCCESS, "Unable to read curve attribute ID.");
                ASSERTL0(faceindx == nextFaceNumber, "Face IDs must begin with zero and be sequential.");

                /// Read face id attribute.
                err = facelement->QueryIntAttribute("FACEID", &faceid);
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

                            // Need to check it here because we may not be good after the read
                            // indicating that there was nothing to read.
                            if (!elementDataStrm.fail())
                            {
                                VertexComponentSharedPtr vert(MemoryManager<VertexComponent>::AllocateSharedPtr(m_meshDimension, faceindx, xval, yval, zval));
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
                    m_curvedFaces.push_back(curve);

                    facelement = facelement->NextSiblingElement("F");

                } // end if-loop
            } // end while-loop
        } // end of ReadCurves()


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
          GetCompositeList(indxStr, m_domain);
          ASSERTL0(!m_domain.empty(), (std::string("Unable to obtain domain's referenced composite: ") + indxStr).c_str());
        }

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

        int MeshGraph::CheckFieldDefinition(const FieldDefinitionsSharedPtr &fielddefs)
          {
          int i;
          ASSERTL0(fielddefs->m_elementIDs.size() > 0, "Fielddefs vector must contain at least one element of data .");

          unsigned int numbasis = 0;

          // Determine nummodes vector lists are correct length
          switch(fielddefs->m_shapeType)
          {
          case eSegment:
              numbasis = 1;
              if(fielddefs->m_numHomogeneousDir)
              {
                  numbasis += fielddefs->m_numHomogeneousDir;
              }

              break;
          case eTriangle:  case eQuadrilateral:
              if(fielddefs->m_numHomogeneousDir)
              {
                  numbasis = 3;
              }
              else
              {
                  numbasis = 2;
              }
              break;
          case eTetrahedron:
          case ePyramid:
          case ePrism:
          case eHexahedron:
              numbasis = 3;
              break;
          }

          unsigned int datasize = 0;

          ASSERTL0(fielddefs->m_basis.size() == numbasis, "Length of basis vector is incorrect");

          if(fielddefs->m_uniOrder == true)
          {
              unsigned int cnt = 0;
              // calculate datasize
              switch(fielddefs->m_shapeType)
              {
              case eSegment:
                  if(fielddefs->m_numHomogeneousDir == 1)
                  {
                      datasize += fielddefs->m_numModes[cnt++]*fielddefs->m_numModes[cnt++];
                  }
                  else if(fielddefs->m_numHomogeneousDir == 2)
                  {
                      datasize += fielddefs->m_numModes[cnt++]*fielddefs->m_numModes[cnt++]*fielddefs->m_numModes[cnt++];
                  }
                  else
                  {
                      datasize += fielddefs->m_numModes[cnt++];
                  }
                  break;
              case eTriangle:
                  {
                      int l = fielddefs->m_numModes[cnt++];
                      int m = fielddefs->m_numModes[cnt++];

                      if(fielddefs->m_numHomogeneousDir == 1)
                      {
                          datasize += StdRegions::StdTriData::getNumberOfCoefficients(l,m)*fielddefs->m_numModes[cnt++];
                      }
                      else
                      {
                          datasize += StdRegions::StdTriData::getNumberOfCoefficients(l,m);
                      }
                  }
                  break;
              case eQuadrilateral:
                  {
                      if(fielddefs->m_numHomogeneousDir == 1)
                      {
                          datasize += fielddefs->m_numModes[cnt++]*
                              fielddefs->m_numModes[cnt++]*
                              fielddefs->m_numModes[cnt++];
                      }
                      else
                      {
                          datasize += fielddefs->m_numModes[cnt++]*
                              fielddefs->m_numModes[cnt++];
                      }
                  }
                  break;
              case eTetrahedron:
                  {
                      int l = fielddefs->m_numModes[cnt++];
                      int m = fielddefs->m_numModes[cnt++];
                      int n = fielddefs->m_numModes[cnt++];
                      datasize += StdRegions::StdTetData::getNumberOfCoefficients(l,m,n);
                  }
                  break;
              case ePyramid:
                  {
                      int l = fielddefs->m_numModes[cnt++];
                      int m = fielddefs->m_numModes[cnt++];
                      int n = fielddefs->m_numModes[cnt++];
                      datasize += StdRegions::StdPyrData::getNumberOfCoefficients(l,m,n);
                  }
                  break;
              case  ePrism:
                  {
                      int l = fielddefs->m_numModes[cnt++];
                      int m = fielddefs->m_numModes[cnt++];
                      int n = fielddefs->m_numModes[cnt++];
                      datasize += StdRegions::StdPrismData::getNumberOfCoefficients(l,m,n);
                  }
                  break;
              case eHexahedron:
                  datasize += fielddefs->m_numModes[cnt++]*
                      fielddefs->m_numModes[cnt++]*
                      fielddefs->m_numModes[cnt++];
                  break;
              }

              datasize *= fielddefs->m_elementIDs.size();
          }
          else
          {
              unsigned int cnt = 0;
              // calculate data length
              for(i = 0; i < fielddefs->m_elementIDs.size(); ++i)
              {
                  switch(fielddefs->m_shapeType)
                  {
                  case eSegment:
                      datasize += fielddefs->m_numModes[cnt++];
                      break;
                  case eTriangle:
                      {
                          int l = fielddefs->m_numModes[cnt++];
                          int m = fielddefs->m_numModes[cnt++];
                          datasize += StdRegions::StdTriData::getNumberOfCoefficients(l,m);
                      }
                      break;
                  case eQuadrilateral:
                      datasize += fielddefs->m_numModes[cnt++]*
                          fielddefs->m_numModes[cnt++];
                      break;
                  case eTetrahedron:
                      {
                          int l = fielddefs->m_numModes[cnt++];
                          int m = fielddefs->m_numModes[cnt++];
                          int n = fielddefs->m_numModes[cnt++];
                          datasize += StdRegions::StdTetData::getNumberOfCoefficients(l,m,n);
                      }
                      break;
                  case ePyramid:
                      {
                          int l = fielddefs->m_numModes[cnt++];
                          int m = fielddefs->m_numModes[cnt++];
                          int n = fielddefs->m_numModes[cnt++];
                          datasize += StdRegions::StdPyrData::getNumberOfCoefficients(l,m,n);
                      }
                      break;
                  case  ePrism:
                      {
                          int l = fielddefs->m_numModes[cnt++];
                          int m = fielddefs->m_numModes[cnt++];
                          int n = fielddefs->m_numModes[cnt++];
                          datasize += StdRegions::StdPrismData::getNumberOfCoefficients(l,m,n);
                      }
                      break;
                  case eHexahedron:
                      datasize += fielddefs->m_numModes[cnt++]*
                          fielddefs->m_numModes[cnt++]*
                          fielddefs->m_numModes[cnt++];
                      break;
                  }
              }
          }

          return datasize;
        }

        void MeshGraph::Write(const std::string &outFile, std::vector<
    		    FieldDefinitionsSharedPtr> &fielddefs,
    		    std::vector<std::vector<NekDouble> > &fielddata)
        {
    	    ASSERTL1(fielddefs.size() == fielddata.size(),
    			    "Length of fielddefs and fielddata incompatible");

    	    TiXmlDocument doc;
    	    TiXmlDeclaration * decl = new TiXmlDeclaration("1.0", "utf-8", "");
    	    doc.LinkEndChild(decl);

    	    cout << "Writing outfile: " << outFile << endl;

    	    TiXmlElement * root = new TiXmlElement("NEKTAR");
    	    doc.LinkEndChild(root);

    	    for (int f = 0; f < fielddefs.size(); ++f)
    	    {

    		    ASSERTL1(fielddata[f].size() > 0,
    				    "Fielddata vector must contain at least one value.");

    		    int datasize = CheckFieldDefinition(fielddefs[f]);
    		    ASSERTL1(fielddata[f].size() == fielddefs[f]->m_fields.size()
    				    * datasize, "Invalid size of fielddata vector.");

    		    //---------------------------------------------
    		    // Write ELEMENTS
    		    TiXmlElement * elemTag = new TiXmlElement("ELEMENTS");
    		    root->LinkEndChild(elemTag);

    		    // Write FIELDS
    		    std::string fieldsString;
    		    {
                        std::stringstream fieldsStringStream;
                        bool first = true;
                        for (std::vector<int>::size_type i = 0; i
                                 < fielddefs[f]->m_fields.size(); i++)
                        {
    				    if (!first)
                                        fieldsStringStream << ",";
    				    fieldsStringStream << fielddefs[f]->m_fields[i];
    				    first = false;
    			    }
    			    fieldsString = fieldsStringStream.str();
    		    }
    		    elemTag->SetAttribute("FIELDS", fieldsString);

    		    // Write SHAPE
    		    std::string shapeString;
    		    {
                        std::stringstream shapeStringStream;
                        shapeStringStream << GeomShapeTypeMap[fielddefs[f]->m_shapeType];
                        if(fielddefs[f]->m_numHomogeneousDir == 1)
                        {
                            shapeStringStream << "-HomogenousExp1D";
                        }
                        else if (fielddefs[f]->m_numHomogeneousDir == 2)
                        {
                            shapeStringStream << "-HomogenousExp2D";
                        }

                        shapeString = shapeStringStream.str();
    		    }
    		    elemTag->SetAttribute("SHAPE", shapeString);

    		    // Write BASIS
    		    std::string basisString;
    		    {
    			    std::stringstream basisStringStream;
    			    bool first = true;
    			    for (std::vector<LibUtilities::BasisType>::size_type i = 0; i < fielddefs[f]->m_basis.size(); i++)
    			    {
    				    if (!first)
    					    basisStringStream << ",";
    				    basisStringStream
    				    << LibUtilities::BasisTypeMap[fielddefs[f]->m_basis[i]];
    				    first = false;
    			    }
    			    basisString = basisStringStream.str();
    		    }
    		    elemTag->SetAttribute("BASIS", basisString);

                    // Write homogeneuous length details
                    if(fielddefs[f]->m_numHomogeneousDir)
                    {
                        std::string homoLenString;
                        {
    			    std::stringstream homoLenStringStream;
    			    bool first = true;
    			    for (int i = 0; i < fielddefs[f]->m_numHomogeneousDir; ++i)
    			    {
                                if (!first)
                                    homoLenStringStream << ",";
                                homoLenStringStream
    				    << fielddefs[f]->m_homogeneousLengths[i];
                                first = false;
    			    }
    			    homoLenString = homoLenStringStream.str();
                        }
                        elemTag->SetAttribute("HOMOGENEOUSLENGTHS", homoLenString);
                    }

    		    // Write NUMMODESPERDIR
    		    std::string numModesString;
    		    {
                        std::stringstream numModesStringStream;

                        if (fielddefs[f]->m_uniOrder)
                        {
                            numModesStringStream << "UNIORDER:";
                            // Just dump single definition
                            bool first = true;
                            for (std::vector<int>::size_type i = 0; i
                                     < fielddefs[f]->m_basis.size(); i++)
                            {
                                if (!first)
                                    numModesStringStream << ",";
                                numModesStringStream << fielddefs[f]->m_numModes[i];
                                first = false;
                            }
                        }
                        else
                        {
                            numModesStringStream << "MIXORDER:";
                            bool first = true;
                            for (std::vector<int>::size_type i = 0; i
                                     < fielddefs[f]->m_numModes.size(); i++)
                            {
                                if (!first)
                                    numModesStringStream << ",";
                                numModesStringStream << fielddefs[f]->m_numModes[i];
                                first = false;
                            }
                        }

                        numModesString = numModesStringStream.str();
    		    }
    		    elemTag->SetAttribute("NUMMODESPERDIR", numModesString);

    		    //Write ID
    		    // Should ideally look at ways of compressing this stream
    		    // if just sequential;
    		    std::string idString;
    		    {
                        std::stringstream idStringStream;
                        bool first = true;
                        for (std::vector<Geometry>::size_type i = 0; i
                                 < fielddefs[f]->m_elementIDs.size(); i++)
                        {
                            if (!first)
                                idStringStream << ",";
                            idStringStream << fielddefs[f]->m_elementIDs[i];
                            first = false;
                        }
                        idString = idStringStream.str();
    		    }
    		    elemTag->SetAttribute("ID", idString);

    		    // Write binary data
    		    std::string compressedDataString;
    		    {
                        // Serialize the fielddata vector to the stringstream.
                        std::stringstream archiveStringStream(std::string((char*) &fielddata[f][0], sizeof(fielddata[f][0]) / sizeof(char) * fielddata[f].size()));

                        // Compress the serialized data.
                        std::stringstream compressedData;
                        {
                            boost::iostreams::filtering_streambuf<boost::iostreams::input>   out;
                            out.push(boost::iostreams::zlib_compressor());
                            out.push(archiveStringStream);
                            boost::iostreams::copy(out, compressedData);
                        }

                        // If the string length is not divisible by 3,
                        // pad it. There is a bug in transform_width
                        // that will make it reference past the end
                        // and crash.
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
    		    std::string base64string(base64_t(compressedDataString.begin()),
                                             base64_t(compressedDataString.end()));
    		    elemTag->LinkEndChild(new TiXmlText(base64string));

    	    }
    	    doc.SaveFile(outFile);
        }

        // Allow global id string
        // Use singe shape per block, single expansion per block
        // Change numModes to numModesPerDir, add initial identifier
        // to be "UniOrder" or "MixOrder"

        void MeshGraph::Import(const std::string& infilename, std::vector<FieldDefinitionsSharedPtr> &fielddefs, std::vector<std::vector<NekDouble> > &fielddata)
        {
          TiXmlDocument doc(infilename);
          bool loadOkay = doc.LoadFile();

          std::stringstream errstr;
          errstr << "Unable to load file: " << infilename << std::endl;
          errstr << "Reason: " << doc.ErrorDesc() << std::endl;
          errstr << "Position: Line " << doc.ErrorRow() << ", Column " << doc.ErrorCol() << std::endl;
          ASSERTL0(loadOkay, errstr.str());

          ImportFieldDefs(doc, fielddefs, false);
          ImportFieldData(doc, fielddefs, fielddata);
        }

        /**
         * The bool decides if the FieldDefs are in <EXPANSIONS> or in <NEKTAR>.
         */
        void MeshGraph::ImportFieldDefs(TiXmlDocument &doc, std::vector<FieldDefinitionsSharedPtr> &fielddefs, bool expChild)
          {
              ASSERTL1(fielddefs.size() == 0, "Expected an empty fielddefs vector.");

              TiXmlHandle docHandle(&doc);
              TiXmlElement* master = NULL;    // Master tag within which all data is contained.

              master = doc.FirstChildElement("NEKTAR");
              ASSERTL0(master, "Unable to find NEKTAR tag in file.");
              std::string strLoop = "NEKTAR";
              TiXmlElement* loopXml = master;

              TiXmlElement *expansionTypes;
              if(expChild)
              {
                  expansionTypes = master->FirstChildElement("EXPANSIONS");
                  ASSERTL0(expansionTypes, "Unable to find EXPANSIONS tag in file.");
                  loopXml = expansionTypes;
                  strLoop = "EXPANSIONS";
              }

              // Loop through all nektar tags, finding all of the element tags.
              while (loopXml)
              {
                  TiXmlElement* element = loopXml->FirstChildElement("ELEMENTS");
                  ASSERTL0(element, "Unable to find ELEMENTS tag within nektar tag.");

                  while (element)
                  {
                      // Extract the attributes.
                      std::string idString;
                      std::string shapeString;
                      std::string basisString;
                      std::string homoLengthsString;
                      std::string numModesString;
                      std::string numPointsString;
                      std::string fieldsString;
                      std::string pointsString;
                      bool pointDef = false;
                      bool numPointDef = false;
                      TiXmlAttribute *attr = element->FirstAttribute();
                      while (attr)
                      {
                          std::string attrName(attr->Name());
                          if (attrName == "FIELDS")
                          {
                              fieldsString.insert(0, attr->Value());
                          }
                          else if (attrName == "SHAPE")
                          {
                              shapeString.insert(0, attr->Value());
                          }
                          else if (attrName == "BASIS")
                          {
                              basisString.insert(0, attr->Value());
                          }
                          else if (attrName == "HOMOGENEOUSLENGTHS")
                          {
                              homoLengthsString.insert(0,attr->Value());
                          }
                          else if (attrName == "NUMMODESPERDIR")
                          {
                              numModesString.insert(0, attr->Value());
                          }
                          else if (attrName == "ID")
                          {
                              idString.insert(0, attr->Value());
                          }
                          else if (attrName == "POINTSTYPE")
                          {
                              pointsString.insert(0, attr->Value());
                              pointDef = true;
                          }
                          else if (attrName == "NUMPOINTSPERDIR")
                          {
                              numPointsString.insert(0, attr->Value());
                              numPointDef = true;
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

                      // Check to see if homogeneous expansion and if so
                      // strip down the shapeString definition
                      int numHomoDir = 0;
                      size_t loc;
                      //---> This finds the first location of  'n'!
                      if((loc = shapeString.find_first_of("-"))!=string::npos)
                      {
                          if(shapeString.find("Exp1D")!=string::npos)
                          {
                              numHomoDir = 1;
                          }
                          else // HomogeneousExp1D
                          {
                              numHomoDir = 2;
                          }

                          shapeString.erase(loc,shapeString.length());
                      }

                      // Reconstruct the fielddefs.
                      std::vector<unsigned int> elementIds;
                      {
                          bool valid = ParseUtils::GenerateSeqVector(idString.c_str(), elementIds);
                          ASSERTL0(valid, "Unable to correctly parse the element ids.");
                      }

                      // Get the geometrical shape
                      SpatialDomains::GeomShapeType shape;
                      bool valid = false;
                      for (unsigned int j = 0; j < SpatialDomains::SIZE_GeomShapeType; j++)
                      {
                          if (SpatialDomains::GeomShapeTypeMap[j] == shapeString)
                          {
                              shape = (SpatialDomains::GeomShapeType) j;
                              valid = true;
                              break;
                          }
                      }

                      ASSERTL0(valid, std::string("Unable to correctly parse the shape type: ").append(shapeString).c_str());

                      // Get the basis
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

                      // Get homoLengths
                      std::vector<NekDouble> homoLengths;
                      if(numHomoDir)
                      {
                          valid = ParseUtils::GenerateOrderedVector(homoLengthsString.c_str(), homoLengths);
                          ASSERTL0(valid, "Unable to correctly parse the number of homogeneous lengths.");
                      }


                      // Get points type
                      std::vector<LibUtilities::PointsType> points;

                      if(pointDef)
                      {
                          std::vector<std::string> pointsStrings;
                          valid = ParseUtils::GenerateOrderedStringVector(pointsString.c_str(), pointsStrings);
                          ASSERTL0(valid, "Unable to correctly parse the points types.");
                          for (std::vector<std::string>::size_type i = 0; i < pointsStrings.size(); i++)
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
                      }

                      // Get numModes
                      std::vector<unsigned int> numModes;
                      bool UniOrder = false;

                      if(strstr(numModesString.c_str(),"UNIORDER:"))
                      {
                          UniOrder  = true;
                      }

                      valid = ParseUtils::GenerateOrderedVector(numModesString.c_str()+9, numModes);
                      ASSERTL0(valid, "Unable to correctly parse the number of modes.");

                      // Get numPoints
                      std::vector<unsigned int> numPoints;
                      if(numPointDef)
                      {
                          valid = ParseUtils::GenerateOrderedVector(numPointsString.c_str(), numPoints);
                          ASSERTL0(valid, "Unable to correctly parse the number of points.");
                      }

                      // Get fields names
                      std::vector<std::string> Fields;
                      valid = ParseUtils::GenerateOrderedStringVector(fieldsString.c_str(), Fields);
                      ASSERTL0(valid, "Unable to correctly parse the number of fields.");

                      SpatialDomains::FieldDefinitionsSharedPtr fielddef  = MemoryManager<SpatialDomains::FieldDefinitions>::AllocateSharedPtr(shape, elementIds, basis, UniOrder, numModes, Fields, numHomoDir, homoLengths, points, pointDef, numPoints, numPointDef);
                      int datasize = CheckFieldDefinition(fielddef);

                      fielddefs.push_back(fielddef);

                      element = element->NextSiblingElement("ELEMENTS");
                  }
                  loopXml = loopXml->NextSiblingElement(strLoop);
              }
          }

        void MeshGraph::ImportFieldData(TiXmlDocument &doc, const std::vector<FieldDefinitionsSharedPtr> &fielddefs, std::vector<std::vector<NekDouble> > &fielddata)
        {
            int cntdumps = 0;
            ASSERTL1(fielddata.size() == 0, "Expected an empty fielddata vector.");

            TiXmlHandle docHandle(&doc);
            TiXmlElement* master = NULL;    // Master tag within which all data is contained.

            master = doc.FirstChildElement("NEKTAR");
            ASSERTL0(master, "Unable to find NEKTAR tag in file.");

            // Loop through all nektar tags, finding all of the element tags.
            while (master)
            {
                TiXmlElement* element = master->FirstChildElement("ELEMENTS");
                ASSERTL0(element, "Unable to find ELEMENTS tag within nektar tag.");
                while (element)
                {
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
                    std::string vData = elementDecompressedData.str();
                    NekDouble* readFieldData = (NekDouble*) vData.c_str();
                    std::vector<NekDouble> elementFieldData(readFieldData, readFieldData + elementDecompressedData.str().length() * sizeof(*elementDecompressedData.str().c_str()) / sizeof(NekDouble));
                    fielddata.push_back(elementFieldData);

                    int datasize = CheckFieldDefinition(fielddefs[cntdumps]);
                    ASSERTL0(fielddata[cntdumps].size() == datasize*fielddefs[cntdumps]->m_fields.size(),"Input data is not the same length as header information");

                    cntdumps++;

                    element = element->NextSiblingElement("ELEMENTS");
                }
                master = master->NextSiblingElement("NEKTAR");
            }
        }

        MeshGraph::~MeshGraph()
          {
          }
  }; //end of namespace
}; //end of namespace

//
// $Log: MeshGraph.cpp,v $
// Revision 1.41  2009/12/16 21:10:52  bnelson
// Updated file read methods to take const std::string& instead of std::string&
//
// Revision 1.40  2009/12/15 18:09:02  cantwell
// Split GeomFactors into 1D, 2D and 3D
// Added generation of tangential basis into GeomFactors
// Updated ADR2DManifold solver to use GeomFactors for tangents
// Added <GEOMINFO> XML session section support in MeshGraph
// Fixed const-correctness in VmathArray
// Cleaned up LocalRegions code to generate GeomFactors
// Removed GenSegExp
// Temporary fix to SubStructuredGraph
// Documentation for GlobalLinSys and GlobalMatrix classes
//
// Revision 1.39  2009/11/18 22:31:46  bnelson
// Changed Write parameter list to accept a const string& as a first parameter.
//
// Revision 1.38  2009/10/22 17:34:20  cbiotto
// Fixing bug for variable order expansion
//
// Revision 1.37  2009/09/24 10:58:12  cbiotto
// Updates for variable order expansions
//
// Revision 1.36  2009/09/06 21:56:42  sherwin
// Updates for Navier-Stokes solver including restart details
//
// Revision 1.35  2009/08/26 11:42:47  rcantao
// - Adequation to new namespace scheme of Spirit version 2. Tested against Boost
// 1.35 and 1.39.
// - Cleared some deprecated headers warnings.
//
// Revision 1.34  2009/08/19 14:13:34  claes
// Removed Gauss-Kronrod parts
//
// Revision 1.33  2009/07/02 13:26:12  sehunchun
// Changes for 2D manifold imbedded in 3D
//
// Revision 1.32  2009/06/18 11:47:24  claes
// changes supporting the static use of Kronrod points
//
// Revision 1.31  2009/06/15 01:59:21  claes
// Gauss-Kronrod updates
//
// Revision 1.30  2009/04/29 10:55:17  pvos
// made demos working with nodal expansions
//
// Revision 1.29  2009/04/27 21:33:35  sherwin
// Added SolverInfoExists and others
//
// Revision 1.28  2009/04/20 16:13:23  sherwin
// Modified Import and Write functions and redefined how Expansion is used
//
// Revision 1.27  2009/01/12 10:26:59  pvos
// Added input tags for nodal expansions
//
// Revision 1.26  2008/10/04 19:32:46  sherwin
// Added SharedPtr Typedef and replaced MeshDimension with SpaceDimension
//
// Revision 1.25  2008/09/09 14:20:30  sherwin
// Updated to handle curved edges (first working version)
//
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
