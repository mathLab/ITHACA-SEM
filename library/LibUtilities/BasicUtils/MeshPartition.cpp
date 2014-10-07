////////////////////////////////////////////////////////////////////////////////
//
//  File: MeshPartition.cpp
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

#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif

#include <LibUtilities/BasicUtils/MeshPartition.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>

#include <tinyxml.h>

#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/BasicUtils/FileSystem.h>

#include <boost/algorithm/string.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/graph/detail/edge.hpp>
#include <boost/format.hpp>

namespace Nektar
{
    namespace LibUtilities
    {
        MeshPartitionFactory& GetMeshPartitionFactory()
        {
            typedef Loki::SingletonHolder<MeshPartitionFactory,
                Loki::CreateUsingNew,
                Loki::NoDestroy > Type;
            return Type::Instance();
        }

        MeshPartition::MeshPartition(const LibUtilities::SessionReaderSharedPtr& pSession) :
                m_numFields(0),
                m_fieldNameToId(),
                m_comm(pSession->GetComm()),
                m_weightingRequired(false)
        {
            ReadConditions(pSession);
            ReadGeometry(pSession);
            ReadExpansions(pSession);
        }

        MeshPartition::~MeshPartition()
        {

        }

        void MeshPartition::PartitionMesh(int nParts, bool shared)
        {
            ASSERTL0(m_meshElements.size() >= nParts,
                     "Too few elements for this many processes.");
            m_shared = shared;

            if (m_weightingRequired)  WeightElements();
            CreateGraph(m_mesh);
            PartitionGraph(m_mesh, nParts, m_localPartition);
        }

        void MeshPartition::WriteLocalPartition(LibUtilities::SessionReaderSharedPtr& pSession)
        {
            TiXmlDocument vNew;
            TiXmlDeclaration * decl = new TiXmlDeclaration("1.0", "utf-8", "");
            vNew.LinkEndChild(decl);

            TiXmlElement* vElmtNektar;
            vElmtNektar = new TiXmlElement("NEKTAR");

            int rank = m_comm->GetRowComm()->GetRank();
            OutputPartition(pSession, m_localPartition[rank], vElmtNektar);

            vNew.LinkEndChild(vElmtNektar);

            std::string  dirname = pSession->GetSessionName() + "_xml"; 
            fs::path    pdirname(dirname);
            
            boost::format pad("P%1$07d.xml");
            pad % rank;
            fs::path    pFilename(pad.str());
            
            if(!fs::is_directory(dirname))
            {
                fs::create_directory(dirname);
            }
            
            fs::path fullpath = pdirname / pFilename; 
            vNew.SaveFile(PortablePath(fullpath));
        }

        void MeshPartition::WriteAllPartitions(LibUtilities::SessionReaderSharedPtr& pSession)
        {
            for (int i = 0; i < m_localPartition.size(); ++i)
            {
                TiXmlDocument vNew;
                TiXmlDeclaration * decl = new TiXmlDeclaration("1.0", "utf-8", "");
                vNew.LinkEndChild(decl);

                TiXmlElement* vElmtNektar;
                vElmtNektar = new TiXmlElement("NEKTAR");

                OutputPartition(pSession, m_localPartition[i], vElmtNektar);

                vNew.LinkEndChild(vElmtNektar);

                std::string  dirname = pSession->GetSessionName() + "_xml"; 
                fs::path    pdirname(dirname);
                
                boost::format pad("P%1$07d.xml");
                pad % i;
                fs::path    pFilename(pad.str());
                
                fs::path fullpath = pdirname / pFilename; 
                
                if(!fs::is_directory(dirname))
                {
                    fs::create_directory(dirname);
                }

                vNew.SaveFile(PortablePath(fullpath));
            }
        }

        void MeshPartition::GetCompositeOrdering(CompositeOrdering &composites)
        {
            std::map<int, MeshEntity>::iterator it;
            for (it  = m_meshComposites.begin();
                 it != m_meshComposites.end(); ++it)
            {
                composites[it->first] = it->second.list;
            }
        }

        void MeshPartition::GetBndRegionOrdering(BndRegionOrdering &bndRegs)
        {
            bndRegs = m_bndRegOrder;
        }


        void MeshPartition::ReadExpansions(const LibUtilities::SessionReaderSharedPtr& pSession)
        {
            // Find the Expansions tag
            TiXmlElement *expansionTypes = pSession->GetElement("Nektar/Expansions");

            // Find the Expansion type
            TiXmlElement *expansion = expansionTypes->FirstChildElement();
            std::string   expType   = expansion->Value();


            if(expType != "E")
            {
                ASSERTL0(false,"Expansion type not defined or not supported at the moment");
            }

            /// Expansiontypes will contain plenty of data,
            /// where relevant at this stage are composite
            /// ID(s) that this expansion type describes,
            /// nummodes and a list of fields that this
            /// expansion relates to. If this does not exist
            /// the variable is only set to "DefaultVar".

            while (expansion)
            {
                std::vector<unsigned int> composite;
                std::vector<unsigned int> nummodes;
                std::vector<std::string>  fieldName;

                const char *nModesStr = expansion->Attribute("NUMMODES");
                ASSERTL0(nModesStr,"NUMMODES was not defined in EXPANSION section of input");
                std::string numModesStr = nModesStr;
                bool valid = ParseUtils::GenerateOrderedVector(numModesStr.c_str(), nummodes);
                ASSERTL0(valid, "Unable to correctly parse the number of modes.");

                if (nummodes.size() == 1)
                {
                    for (int i = 1; i < m_dim; i++)
                    {
                        nummodes.push_back( nummodes[0] );
                    }
                }
                ASSERTL0(nummodes.size() == m_dim,"Number of modes should match mesh dimension");


                const char *fStr = expansion->Attribute("FIELDS");
                if(fStr)
                {
                    std::string fieldStr = fStr;
                    bool  valid = ParseUtils::GenerateOrderedStringVector(fieldStr.c_str(),fieldName);
                    ASSERTL0(valid,"Unable to correctly parse the field string in ExpansionTypes.");

                    for (int i = 0; i < fieldName.size(); ++i)
                    {
                        if (m_fieldNameToId.count(fieldName[i]) == 0)
                        {
                            int k = m_fieldNameToId.size();
                            m_fieldNameToId[ fieldName[i] ] = k;
                            m_numFields++;
                        }
                    }
                }
                else
                {
                    fieldName.push_back("DefaultVar");
                    int k = m_fieldNameToId.size();

                    if (m_fieldNameToId.count("DefaultVar") == 0)
                    {
                        ASSERTL0(k == 0,
                                 "Omitting field variables and explicitly listing " \
                                 "them in different ExpansionTypes is wrong practise");

                        m_fieldNameToId[ "DefaultVar" ] = k;
                        m_numFields++;
                    }
                }

                std::string compositeStr = expansion->Attribute("COMPOSITE");
                ASSERTL0(compositeStr.length() > 3, "COMPOSITE must be specified in expansion definition");
                int beg = compositeStr.find_first_of("[");
                int end = compositeStr.find_first_of("]");
                std::string compositeListStr = compositeStr.substr(beg+1,end-beg-1);
                bool parseGood = ParseUtils::GenerateSeqVector(compositeListStr.c_str(), composite);
                ASSERTL0(parseGood && !composite.empty(),
                    (std::string("Unable to read composite index range: ") + compositeListStr).c_str());


                // construct mapping (field name, CompositeID) -> nummodes
                for (int i = 0; i < composite.size(); ++i)
                {
                    for (int j = 0; j < fieldName.size(); j++)
                    {
                        m_expansions[composite[i]][fieldName[j]] = nummodes;
                    }
                }

                expansion = expansion->NextSiblingElement("E");
            }
        }




        void MeshPartition::ReadGeometry(const LibUtilities::SessionReaderSharedPtr& pSession)
        {
            TiXmlElement* x;
            TiXmlElement *vGeometry, *vSubElement;
            int i;

            vGeometry = pSession->GetElement("Nektar/Geometry");
            m_dim = atoi(vGeometry->Attribute("DIM"));

            // Read mesh vertices
            vSubElement = pSession->GetElement("Nektar/Geometry/Vertex");

            // Retrieve any VERTEX attributes specifying mesh transforms
            std::string attr[] = {"XSCALE", "YSCALE", "ZSCALE",
                                  "XMOVE",  "YMOVE",  "ZMOVE" };
            for (i = 0; i < 6; ++i)
            {
                const char *val =  vSubElement->Attribute(attr[i].c_str());
                if (val)
                {
                    m_vertexAttributes[attr[i]] = std::string(val);
                }
            }

            x = vSubElement->FirstChildElement();
            i = 0;
            if (x->FirstAttribute())
            {
                i = x->FirstAttribute()->IntValue();
            }
            while(x)
            {
                TiXmlAttribute* y = x->FirstAttribute();
                ASSERTL0(y, "Failed to get attribute.");
                MeshVertex v;
                v.id = y->IntValue();
                ASSERTL0(v.id == i++, "Vertex IDs not sequential.");
                std::vector<std::string> vCoords;
                std::string vCoordStr = x->FirstChild()->ToText()->Value();
                boost::split(vCoords, vCoordStr, boost::is_any_of("\t "));
                v.x = atof(vCoords[0].c_str());
                v.y = atof(vCoords[1].c_str());
                v.z = atof(vCoords[2].c_str());
                m_meshVertices[v.id] = v;
                x = x->NextSiblingElement();
            }

            // Read mesh edges
            if (m_dim >= 2)
            {
                vSubElement = pSession->GetElement("Nektar/Geometry/Edge");
                ASSERTL0(vSubElement, "Cannot read edges");
                x = vSubElement->FirstChildElement();
                i = 0;
                while(x)
                {
                    TiXmlAttribute* y = x->FirstAttribute();
                    ASSERTL0(y, "Failed to get attribute.");
                    MeshEntity e;
                    e.id = y->IntValue();
                    e.type = 'E';
                    ASSERTL0(e.id == i++, "Edge IDs not sequential.");
                    std::vector<std::string> vVertices;
                    std::string vVerticesString = x->FirstChild()->ToText()->Value();
                    boost::split(vVertices, vVerticesString, boost::is_any_of("\t "));
                    e.list.push_back(atoi(vVertices[0].c_str()));
                    e.list.push_back(atoi(vVertices[1].c_str()));
                    m_meshEdges[e.id] = e;
                    x = x->NextSiblingElement();
                }
            }

            // Read mesh faces
            if (m_dim == 3)
            {
                vSubElement = pSession->GetElement("Nektar/Geometry/Face");
                ASSERTL0(vSubElement, "Cannot read faces.");
                x = vSubElement->FirstChildElement();
                i = 0;
                while(x)
                {
                    TiXmlAttribute* y = x->FirstAttribute();
                    ASSERTL0(y, "Failed to get attribute.");
                    MeshEntity f;
                    f.id = y->IntValue();
                    f.type = x->Value()[0];
                    ASSERTL0(f.id == i++, "Face IDs not sequential.");
                    std::vector<std::string> vEdges;
                    std::string vEdgeStr = x->FirstChild()->ToText()->Value();
                    boost::split(vEdges, vEdgeStr, boost::is_any_of("\t "));
                    for (int i = 0; i < vEdges.size(); ++i)
                    {
                        f.list.push_back(atoi(vEdges[i].c_str()));
                    }
                    m_meshFaces[f.id] = f;
                    x = x->NextSiblingElement();
                }
            }

            // Read mesh elements
            vSubElement = pSession->GetElement("Nektar/Geometry/Element");
            ASSERTL0(vSubElement, "Cannot read elements.");
            x = vSubElement->FirstChildElement();
            i = 0;
            while(x)
            {
                TiXmlAttribute* y = x->FirstAttribute();
                ASSERTL0(y, "Failed to get attribute.");
                MeshEntity e;
                e.id = y->IntValue();
                ASSERTL0(e.id == i++, "Element IDs not sequential.");
                std::vector<std::string> vItems;
                std::string vItemStr = x->FirstChild()->ToText()->Value();
                boost::split(vItems, vItemStr, boost::is_any_of("\t "));
                for (int i = 0; i < vItems.size(); ++i)
                {
                    e.list.push_back(atoi(vItems[i].c_str()));
                }
                e.type = x->Value()[0];
                m_meshElements[e.id] = e;
                x = x->NextSiblingElement();
            }

            // Read mesh curves
            if (pSession->DefinesElement("Nektar/Geometry/Curved"))
            {
                vSubElement = pSession->GetElement("Nektar/Geometry/Curved");
                x = vSubElement->FirstChildElement();
                i = 0;
                while(x)
                {
                    MeshCurved c;
                    ASSERTL0(x->Attribute("ID", &c.id),
                             "Failed to get attribute ID");
                    c.type = std::string(x->Attribute("TYPE"));
                    ASSERTL0(!c.type.empty(),
                             "Failed to get attribute TYPE");
                    ASSERTL0(x->Attribute("NUMPOINTS", &c.npoints),
                             "Failed to get attribute NUMPOINTS");
                    c.data = x->FirstChild()->ToText()->Value();
                    c.entitytype = x->Value()[0];
                    if (c.entitytype == "E")
                    {
                        ASSERTL0(x->Attribute("EDGEID", &c.entityid),
                             "Failed to get attribute EDGEID");
                    }
                    else if (c.entitytype == "F")
                    {
                        ASSERTL0(x->Attribute("FACEID", &c.entityid),
                             "Failed to get attribute FACEID");
                    }
                    else
                    {
                        ASSERTL0(false, "Unknown curve type.");
                    }
                    m_meshCurved[std::make_pair(c.entitytype, c.id)] = c;
                    x = x->NextSiblingElement();
                }
            }

            // Read composites
            vSubElement = pSession->GetElement("Nektar/Geometry/Composite");
            ASSERTL0(vSubElement, "Cannot read composites.");
            x = vSubElement->FirstChildElement();
            i = 0;
            while(x)
            {
                TiXmlAttribute* y = x->FirstAttribute();
                ASSERTL0(y, "Failed to get attribute.");
                MeshEntity c;
                c.id = y->IntValue();
                std::string vSeqStr = x->FirstChild()->ToText()->Value();
                c.type = vSeqStr[0];
                std::string::size_type indxBeg = vSeqStr.find_first_of('[') + 1;
                std::string::size_type indxEnd = vSeqStr.find_last_of(']') - 1;
                vSeqStr = vSeqStr.substr(indxBeg, indxEnd - indxBeg + 1);

                std::vector<unsigned int> vSeq;
                ParseUtils::GenerateSeqVector(vSeqStr.c_str(), vSeq);

                for (int i = 0; i < vSeq.size(); ++i)
                {
                    c.list.push_back(vSeq[i]);
                }
                m_meshComposites[c.id] = c;
                x = x->NextSiblingElement();
            }

            // Read Domain
            vSubElement = pSession->GetElement("Nektar/Geometry/Domain");
            ASSERTL0(vSubElement, "Cannot read domain");
            std::string vSeqStr = vSubElement->FirstChild()->ToText()->Value();
            std::string::size_type indxBeg = vSeqStr.find_first_of('[') + 1;
            std::string::size_type indxEnd = vSeqStr.find_last_of(']') - 1;
            vSeqStr = vSeqStr.substr(indxBeg, indxEnd - indxBeg + 1);
            ParseUtils::GenerateSeqVector(vSeqStr.c_str(), m_domain);
        }

        void MeshPartition::PrintPartInfo(std::ostream &out)
        {
            int nElmt = boost::num_vertices(m_mesh);
            int nPart = m_localPartition.size();

            out << "# Partition information:" << std::endl;
            out << "# No. elements  : " << nElmt << std::endl;
            out << "# No. partitions: " << nPart << std::endl;
            out << "# ID  nElmt  nLocDof  nBndDof" << std::endl;

            BoostVertexIterator vertit, vertit_end;
            std::vector<int> partElmtCount(nPart, 0);
            std::vector<int> partLocCount (nPart, 0);
            std::vector<int> partBndCount (nPart, 0);

            std::map<int, int> elmtSizes;
            std::map<int, int> elmtBndSizes;
            
            for (unsigned int i = 0; i < m_domain.size(); ++i)
            {
                int cId = m_domain[i];
                NummodesPerField npf = m_expansions[cId];

                for (NummodesPerField::iterator it = npf.begin(); it != npf.end(); ++it)
                {
                    ASSERTL0(it->second.size() == m_dim,
                        " Number of directional" \
                        " modes in expansion spec for composite id = " + 
                        boost::lexical_cast<std::string>(cId) +
                        " and field " +
                        boost::lexical_cast<std::string>(it->first) +
                        " does not correspond to mesh dimension");

                    int na = it->second[0];
                    int nb = it->second[1];
                    int nc = 0;
                    if (m_dim == 3)
                    {
                        nc = it->second[2];
                    }

                    int weight = 0;
                    int bndWeight = 0;
                    switch (m_meshComposites[cId].type)
                    {
                        case 'A':
                            weight    = StdTetData::getNumberOfCoefficients(na, nb, nc);
                            bndWeight = StdTetData::getNumberOfBndCoefficients(na, nb, nc);
                            break;
                        case 'R':
                            weight    = StdPrismData::getNumberOfCoefficients(na, nb, nc);
                            bndWeight = StdPrismData::getNumberOfBndCoefficients(na, nb, nc);
                            break;
                        case 'H':
                            weight    = StdHexData::getNumberOfCoefficients(na, nb, nc);
                            bndWeight = StdHexData::getNumberOfBndCoefficients(na, nb, nc);
                            break;
                        case 'P':
                            weight    = StdPyrData::getNumberOfCoefficients(na, nb, nc);
                            bndWeight = StdPyrData::getNumberOfBndCoefficients(na, nb, nc);
                            break;
                        case 'Q':
                            weight    = StdQuadData::getNumberOfCoefficients(na, nb);
                            bndWeight = StdQuadData::getNumberOfBndCoefficients(na, nb);
                            break;
                        case 'T':
                            weight    = StdTriData::getNumberOfCoefficients(na, nb);
                            bndWeight = StdTriData::getNumberOfBndCoefficients(na, nb);
                            break;
                        case 'S':
                            weight    = StdSegData::getNumberOfCoefficients(na);
                            bndWeight = StdSegData::getNumberOfBndCoefficients(na);
                            break;
                        case 'V':
                            weight    = 1;
                            bndWeight = 1;
                            break;
                        default:
                            break;
                    }

                    for (unsigned int j = 0; j < m_meshComposites[cId].list.size(); ++j)
                    {
                        int elid = m_meshComposites[cId].list[j]; 
                        elmtSizes[elid] = weight;
                        elmtBndSizes[elid] = bndWeight;
                    }
                }
            }

            for (boost::tie(vertit, vertit_end) = boost::vertices(m_mesh);
                 vertit != vertit_end; ++vertit)
            {
                int partId = m_mesh[*vertit].partition;
                partElmtCount[partId]++;
                partLocCount [partId] += elmtSizes[m_mesh[*vertit].id];
                partBndCount [partId] += elmtBndSizes[m_mesh[*vertit].id];
            }

            for (int i = 0; i < nPart; ++i)
            {
                out << i << " " << partElmtCount[i] << " " << partLocCount[i] << " " << partBndCount[i] << std::endl;
            }
        }

        void MeshPartition::ReadConditions(const SessionReaderSharedPtr& pSession)
        {
            if (!pSession->DefinesElement("Nektar/Conditions/SolverInfo"))
            {
                // No SolverInfo = no change of default action to weight
                // mesh graph.
                return;
            }

            TiXmlElement* solverInfoElement = 
                    pSession->GetElement("Nektar/Conditions/SolverInfo");

            TiXmlElement* solverInfo = 
                    solverInfoElement->FirstChildElement("I");
            ASSERTL0(solverInfo, "Cannot read SolverInfo tags");

            while (solverInfo)
            {
                // read the property name
                ASSERTL0(solverInfo->Attribute("PROPERTY"),
                         "Missing PROPERTY attribute in solver info "
                         "section. ");
                std::string solverProperty = 
                    solverInfo->Attribute("PROPERTY");
                ASSERTL0(!solverProperty.empty(),
                         "Solver info properties must have a non-empty "
                         "name. ");
                // make sure that solver property is capitalised
                std::string solverPropertyUpper =
                    boost::to_upper_copy(solverProperty);


                // read the value
                ASSERTL0(solverInfo->Attribute("VALUE"),
                        "Missing VALUE attribute in solver info section. ");
                std::string solverValue    = solverInfo->Attribute("VALUE");
                ASSERTL0(!solverValue.empty(),
                         "Solver info properties must have a non-empty value");
                // make sure that property value is capitalised
                std::string propertyValueUpper =
                    boost::to_upper_copy(solverValue);

                if (solverPropertyUpper == "WEIGHTPARTITIONS") 
                {
                    if (propertyValueUpper != "UNIFORM")
                    {
                        m_weightingRequired = true;
                    }
                    return;
                }
                solverInfo = solverInfo->NextSiblingElement("I");
            }
        }


        /*
         * Calculate element weights based on
         *   - element type (Q,T,H,P,R,A)
         *   - nummodes in expansion which this element belongs to via composite.
         *
         * For each element we prepare two vertex weightings, one associated
         * with the number of matrix elements associated with it (to balance
         * matrix multiplication work) and another associated
         * with all work which scales linearly with the number of its 
         * coefficients: communication, vector updates etc.
         *
         * \todo Refactor this code to explicitly represent performance model
         * and flexibly generate graph vertex weights depending on perf data.
         */
        void MeshPartition::WeightElements()
        {
            std::vector<unsigned int> weight(m_numFields, 1);
            for (int i = 0; i < m_meshElements.size(); ++i)
            {
                m_vertWeights.push_back( weight );
            }

            for (unsigned int i = 0; i < m_domain.size(); ++i)
            {
                int cId = m_domain[i];
                NummodesPerField npf = m_expansions[cId];

                for (NummodesPerField::iterator it = npf.begin(); it != npf.end(); ++it)
                {
                    ASSERTL0(it->second.size() == m_dim,
                        " Number of directional" \
                        " modes in expansion spec for composite id = " + 
                        boost::lexical_cast<std::string>(cId) +
                        " and field " +
                        boost::lexical_cast<std::string>(it->first) +
                        " does not correspond to mesh dimension");

                    int na = it->second[0];
                    int nb = 0;
                    int nc = 0;
                    if (m_dim >= 2)
                    {
                        nb = it->second[1];
                    }
                    if (m_dim == 3)
                    {
                        nc = it->second[2];
                    }

                    int weight = 0;
                    int bndWeight = 0;
                    switch (m_meshComposites[cId].type)
                    {
                        case 'A':
                            weight    = StdTetData::getNumberOfCoefficients(na, nb, nc);
                            bndWeight = StdTetData::getNumberOfBndCoefficients(na, nb, nc);
                            break;
                        case 'R':
                            weight    = StdPrismData::getNumberOfCoefficients(na, nb, nc);
                            bndWeight = StdPrismData::getNumberOfBndCoefficients(na, nb, nc);
                            break;
                        case 'H':
                            weight    = StdHexData::getNumberOfCoefficients(na, nb, nc);
                            bndWeight = StdHexData::getNumberOfBndCoefficients(na, nb, nc);
                            break;
                        case 'P':
                            weight    = StdPyrData::getNumberOfCoefficients(na, nb, nc);
                            bndWeight = StdPyrData::getNumberOfBndCoefficients(na, nb, nc);
                            break;
                        case 'Q':
                            weight    = StdQuadData::getNumberOfCoefficients(na, nb);
                            bndWeight = StdQuadData::getNumberOfBndCoefficients(na, nb);
                            break;
                        case 'T':
                            weight    = StdTriData::getNumberOfCoefficients(na, nb);
                            bndWeight = StdTriData::getNumberOfBndCoefficients(na, nb);
                            break;
                        case 'S':
                            weight    = StdSegData::getNumberOfCoefficients(na);
                            bndWeight = StdSegData::getNumberOfBndCoefficients(na);
                            break;
                        case 'V':
                            weight    = 1;
                            bndWeight = 1;
                            break;
                        default:
                            break;
                    }

                    for (unsigned int j = 0; j < m_meshComposites[cId].list.size(); ++j)
                    {
                        int elmtId = m_meshComposites[cId].list[j];
                        m_vertWeights[elmtId][ m_fieldNameToId[it->first]] = bndWeight;
                        //m_vertWeights[elmtId][ m_fieldNameToId[ it->first ] * 2 + 1 ] = weight*weight;
                    }
                }
            } // for i
        }

        void MeshPartition::CreateGraph(BoostSubGraph& pGraph)
        {
            // Maps edge/face to first mesh element id.
            // On locating second mesh element id, graph edge is created instead.
            std::map<int, int> vGraphEdges;

            for (unsigned int i = 0; i < m_meshElements.size(); ++i)
            {
                int p = m_meshElements[i].id;
                BoostVertex v = boost::add_vertex(pGraph);
                pGraph[v].id = p;
                pGraph[v].partition = 0;
                if (m_weightingRequired)
                {
                    pGraph[v].weight = m_vertWeights[i];
                }

                // Process element entries and add graph edges
                for (unsigned j = 0; j < m_meshElements[i].list.size(); ++j)
                {
                    int eId = m_meshElements[i].list[j];

                    // Look to see if we've examined this edge/face before
                    // If so, we've got both graph vertices so add edge
                    if (vGraphEdges.find(eId) != vGraphEdges.end())
                    {
                        BoostEdge e = boost::add_edge( p, vGraphEdges[eId], pGraph).first;
                        pGraph[e].id = eId;
                    }
                    else
                    {
                        vGraphEdges[eId] = p;
                    }
                }
            }
        }

        void MeshPartition::PartitionGraph(BoostSubGraph& pGraph,
                                           int nParts,
                                           std::vector<BoostSubGraph>& pLocalPartition)
        {
            int i;
            int nGraphVerts = boost::num_vertices(pGraph);
            int nGraphEdges = boost::num_edges(pGraph);

            // Convert boost graph into CSR format
            BoostVertexIterator    vertit, vertit_end;
            Array<OneD, int> part(nGraphVerts,0);

            if (m_comm->GetRowComm()->TreatAsRankZero())
            {
                int acnt = 0;
                int vcnt = 0;
                int nWeight = nGraphVerts;
                BoostAdjacencyIterator adjvertit, adjvertit_end;
                Array<OneD, int> xadj(nGraphVerts+1,0);
                Array<OneD, int> adjncy(2*nGraphEdges);
                Array<OneD, int> vwgt(nWeight, 1);
                Array<OneD, int> vsize(nGraphVerts, 1);
                for ( boost::tie(vertit, vertit_end) = boost::vertices(pGraph);
                      vertit != vertit_end;
                      ++vertit)
                {
                    for ( boost::tie(adjvertit, adjvertit_end) = boost::adjacent_vertices(*vertit,pGraph);
                          adjvertit != adjvertit_end;
                          ++adjvertit)
                    {
                        adjncy[acnt++] = *adjvertit;

                    }
                    xadj[++vcnt] = acnt;

                    if (m_weightingRequired)
                    {
                        vwgt[pGraph[*vertit].id ] = pGraph[*vertit].weight[0];
                    }
                    else
                    {
                        vwgt[pGraph[*vertit].id] = 1;
                    }
                }

                // Call Metis and partition graph
                int vol = 0;

                try
                {
                    //////////////////////////////////////////////////////
                    // On a cartesian communicator do mesh partiotion just on the first column
                    // so there is no doubt the partitions are all the same in all the columns
                    if(m_comm->GetColumnComm()->GetRank() == 0)
                    {
                        // Attempt partitioning using METIS.
                        int ncon = 1;
                        PartitionGraphImpl(nGraphVerts, ncon, xadj, adjncy, vwgt, vsize, nParts, vol, part);

                        // Check METIS produced a valid partition and fix if not.
                        CheckPartitions(nParts, part);
                        if (!m_shared)
                        {
                            // distribute among columns
                            for (i = 1; i < m_comm->GetColumnComm()->GetSize(); ++i)
                            {
                                m_comm->GetColumnComm()->Send(i, part);
                            }
                        }
                    }
                    else 
                    {
                        m_comm->GetColumnComm()->Recv(0, part);
                    }
                    if (!m_shared)
                    {
                        m_comm->GetColumnComm()->Block();

                        //////////////////////////////////
                        // distribute among rows
                        for (i = 1; i < m_comm->GetRowComm()->GetSize(); ++i)
                        {
                            m_comm->GetRowComm()->Send(i, part);
                        }
                    }
                }
                catch (...)
                {
                    NEKERROR(ErrorUtil::efatal,
                             "Error in calling metis to partition graph.");
                }
            }
            else
            {
                m_comm->GetRowComm()->Recv(0, part);
            }

            // Create boost subgraph for this process's partitions
            int nCols = nParts;
            pLocalPartition.resize(nCols);
            for (i = 0; i < nCols; ++i)
            {
                pLocalPartition[i] = pGraph.create_subgraph();
            }

            // Populate subgraph
            i = 0;
            for ( boost::tie(vertit, vertit_end) = boost::vertices(pGraph);
                  vertit != vertit_end;
                  ++vertit, ++i)
            {
                pGraph[*vertit].partition = part[i];
                pGraph[*vertit].partid = boost::num_vertices(pLocalPartition[part[i]]);
                boost::add_vertex(i, pLocalPartition[part[i]]);
            }
        }


        void MeshPartition::CheckPartitions(int nParts, Array<OneD, int> &pPart)
        {
            unsigned int       i     = 0;
            unsigned int       cnt   = 0;
            const unsigned int npart = m_comm->GetRowComm()->GetSize();
            bool               valid = true;

            // Check that every process has at least one element assigned
            for (i = 0; i < nParts; ++i)
            {
                cnt = std::count(pPart.begin(), pPart.end(), i);
                if (cnt == 0)
                {
                    valid = false;
                }
            }

            // If METIS produced an invalid partition, repartition naively.
            // Elements are assigned to processes in a round-robin fashion.
            // It is assumed that METIS failure only occurs when the number of
            // elements is approx. the number of processes, so this approach
            // should not be too inefficient communication-wise.
            if (!valid)
            {
                for (i = 0; i < pPart.num_elements(); ++i)
                {
                    pPart[i] = i % nParts;
                }
            }
        }


        void MeshPartition::OutputPartition(
                LibUtilities::SessionReaderSharedPtr& pSession,
                BoostSubGraph& pGraph,
                TiXmlElement* pNektar)
        {
            // Write Geometry data
            std::string vDim   = pSession->GetElement("Nektar/Geometry")->Attribute("DIM");
            std::string vSpace = pSession->GetElement("Nektar/Geometry")->Attribute("SPACE");
            std::string vPart  = boost::lexical_cast<std::string>(pGraph[*boost::vertices(pGraph).first].partition);
            TiXmlElement* vElmtGeometry = new TiXmlElement("GEOMETRY");
            vElmtGeometry->SetAttribute("DIM", vDim);
            vElmtGeometry->SetAttribute("SPACE", vSpace);
            vElmtGeometry->SetAttribute("PARTITION", vPart);

            TiXmlElement *vVertex  = new TiXmlElement("VERTEX");
            TiXmlElement *vEdge    = new TiXmlElement("EDGE");
            TiXmlElement *vFace    = new TiXmlElement("FACE");
            TiXmlElement *vElement = new TiXmlElement("ELEMENT");
            TiXmlElement *vCurved  = new TiXmlElement("CURVED");
            TiXmlElement *vComposite = new TiXmlElement("COMPOSITE");
            TiXmlElement *vDomain  = new TiXmlElement("DOMAIN");

            TiXmlElement *x;
            TiXmlText    *y;

            BoostVertexIterator    vertit, vertit_end;
            int id;

            std::map<int, MeshEntity> vComposites;
            std::map<int, MeshEntity> vElements;
            std::map<int, MeshEntity> vEdges;
            std::map<int, MeshEntity> vFaces;
            std::map<int, MeshVertex> vVertices;
            std::map<int, MeshEntity>::iterator vIt;
            std::map<int, MeshVertex>::iterator vVertIt;
            std::map<std::string, std::string>::iterator vAttrIt;

            // Populate lists of elements, edges and vertices required.
            for ( boost::tie(vertit, vertit_end) = boost::vertices(pGraph);
                  vertit != vertit_end;
                  ++vertit)
            {
                id = pGraph[*vertit].id;
                vElements[id] = m_meshElements[pGraph[*vertit].id];
            }

            std::map<int, MeshEntity> * vNext = &vElements;
            switch (m_dim)
            {
                case 3:
                {
                    // Compile list of faces
                    for (vIt = vNext->begin(); vIt != vNext->end(); vIt++)
                    {
                        for (unsigned int j = 0; j < vIt->second.list.size(); ++j)
                        {
                            id = vIt->second.list[j];
                            vFaces[id] = m_meshFaces[id];
                        }
                    }
                    vNext = &vFaces;
                }
                case 2:
                {
                    // Compile list of edges
                    for (vIt = vNext->begin(); vIt != vNext->end(); vIt++)
                    {
                        for (unsigned int j = 0; j < vIt->second.list.size(); ++j)
                        {
                            id = vIt->second.list[j];
                            vEdges[id] = m_meshEdges[id];
                        }
                    }
                    vNext = &vEdges;
                }
                case 1:
                {
                    // Compile list of vertices
                    for (vIt = vNext->begin(); vIt != vNext->end(); vIt++)
                    {
                        for (unsigned int j = 0; j < vIt->second.list.size(); ++j)
                        {
                            id = vIt->second.list[j];
                            vVertices[id] = m_meshVertices[id];
                        }
                    }
                }
            }

            // Generate XML data for these mesh entities
            for (vVertIt = vVertices.begin(); vVertIt != vVertices.end(); vVertIt++)
            {
                x = new TiXmlElement("V");
                x->SetAttribute("ID", vVertIt->first);
                std::stringstream vCoords;
                vCoords.precision(12);
                vCoords << std::setw(15) << vVertIt->second.x
                        << std::setw(15) << vVertIt->second.y
                        << std::setw(15) << vVertIt->second.z << " ";
                y = new TiXmlText(vCoords.str());
                x->LinkEndChild(y);
                vVertex->LinkEndChild(x);
            }

            // Apply transformation attributes to VERTEX section
            for (vAttrIt  = m_vertexAttributes.begin();
                 vAttrIt != m_vertexAttributes.end();
                 ++ vAttrIt)
            {
                vVertex->SetAttribute(vAttrIt->first, vAttrIt->second);
            }

            if (m_dim >= 2)
            {
                for (vIt = vEdges.begin(); vIt != vEdges.end(); vIt++)
                {
                    x = new TiXmlElement("E");
                    x->SetAttribute("ID", vIt->first);
                    std::stringstream vVertices;
                    vVertices << std::setw(10) << vIt->second.list[0]
                            << std::setw(10) << vIt->second.list[1] << " ";
                    y = new TiXmlText(vVertices.str());
                    x->LinkEndChild(y);
                    vEdge->LinkEndChild(x);
                }
            }

            if (m_dim >= 3)
            {
                for (vIt = vFaces.begin(); vIt != vFaces.end(); vIt++)
                {
                    std::string vType("F");
                    vType[0] = vIt->second.type;
                    x = new TiXmlElement(vType);
                    x->SetAttribute("ID", vIt->first);
                    std::stringstream vListStr;
                    for (unsigned int i = 0; i < vIt->second.list.size(); ++i)
                    {
                        vListStr << std::setw(10) << vIt->second.list[i];
                    }
                    vListStr << " ";
                    y = new TiXmlText(vListStr.str());
                    x->LinkEndChild(y);
                    vFace->LinkEndChild(x);
                }
            }

            for (vIt = vElements.begin(); vIt != vElements.end(); vIt++)
            {
                std::string vType("T");
                vType[0] = vIt->second.type;
                x = new TiXmlElement(vType.c_str());
                x->SetAttribute("ID", vIt->first);
                std::stringstream vEdges;
                for (unsigned i = 0; i < vIt->second.list.size(); ++i)
                {
                    vEdges << std::setw(10) << vIt->second.list[i];
                }
                vEdges << " ";
                y = new TiXmlText(vEdges.str());
                x->LinkEndChild(y);
                vElement->LinkEndChild(x);
            }

            if (m_dim >= 2)
            {
                std::map<MeshCurvedKey, MeshCurved>::const_iterator vItCurve;
                for (vItCurve  = m_meshCurved.begin(); 
                     vItCurve != m_meshCurved.end(); 
                     ++vItCurve)
                {
                    MeshCurved c = vItCurve->second;
                    
                    if (vEdges.find(c.entityid) != vEdges.end() || 
                        vFaces.find(c.entityid) != vFaces.end())
                    {
                        x = new TiXmlElement(c.entitytype);
                        x->SetAttribute("ID", c.id);
                        if (c.entitytype == "E")
                        {
                            x->SetAttribute("EDGEID", c.entityid);
                        }
                        else
                        {
                            x->SetAttribute("FACEID", c.entityid);
                        }
                        x->SetAttribute("TYPE", c.type);
                        x->SetAttribute("NUMPOINTS", c.npoints);
                        y = new TiXmlText(c.data);
                        x->LinkEndChild(y);
                        vCurved->LinkEndChild(x);
                    }
                }
            }

            // Generate composites section comprising only those mesh entities
            // which belong to this partition.
            for (vIt = m_meshComposites.begin(); vIt != m_meshComposites.end(); ++vIt)
            {
                bool comma = false; // Set to true after first entity output
                bool range = false; // True when entity IDs form a range
                int last_idx = -2;  // Last entity ID output
                std::string vCompositeStr = "";
                for (unsigned int j = 0; j < vIt->second.list.size(); ++j)
                {
                    // Based on entity type, check if in this partition
                    switch (vIt->second.type)
                    {
                    case 'V':
                        if (vVertices.find(vIt->second.list[j]) == vVertices.end())
                        {
                            continue;
                        }
                        break;
                    case 'E':
                        if (vEdges.find(vIt->second.list[j]) == vEdges.end())
                        {
                            continue;
                        }
                        break;
                    case 'F':
                        if (vFaces.find(vIt->second.list[j]) == vFaces.end())
                        {
                            continue;
                        }
                        break;
                    default:
                        if (vElements.find(vIt->second.list[j]) == vElements.end())
                        {
                            continue;
                        }
                        break;
                    }

                    // Condense consecutive entity IDs into ranges
                    // last_idx initially -2 to avoid error for ID=0
                    if (last_idx + 1 == vIt->second.list[j])
                    {
                        last_idx++;
                        range = true;
                        continue;
                    }
                    // This entity is not in range, so close previous range with
                    // last_idx
                    if (range)
                    {
                        vCompositeStr += "-" + boost::lexical_cast<std::string>(last_idx);
                        range = false;
                    }
                    // Output ID, which is either standalone, or will start a
                    // range.
                    vCompositeStr += comma ? "," : "";
                    vCompositeStr += boost::lexical_cast<std::string>(vIt->second.list[j]);
                    last_idx = vIt->second.list[j];
                    comma = true;
                }
                // If last entity is part of a range, it must be output now
                if (range)
                {
                    vCompositeStr += "-" + boost::lexical_cast<std::string>(last_idx);
                }

                if (vCompositeStr.length() > 0)
                {
                    vComposites[vIt->first] = vIt->second;
                    x = new TiXmlElement("C");
                    x->SetAttribute("ID", vIt->first);
                    vCompositeStr = "X[" + vCompositeStr + "]";
                    vCompositeStr[0] = vIt->second.type;
                    y = new TiXmlText(vCompositeStr.c_str());
                    x->LinkEndChild(y);
                    vComposite->LinkEndChild(x);
                }
            }

            std::string vDomainListStr;
            bool comma = false;
            for (unsigned int i = 0; i < m_domain.size(); ++i)
            {
                if (vComposites.find(m_domain[i]) != vComposites.end())
                {
                    vDomainListStr += comma ? "," : "";
                    comma = true;
                    vDomainListStr += boost::lexical_cast<std::string>(m_domain[i]);
                }
            }
            vDomainListStr = "C[" + vDomainListStr + "]";
            TiXmlText* vDomainList = new TiXmlText(vDomainListStr);
            vDomain->LinkEndChild(vDomainList);

            vElmtGeometry->LinkEndChild(vVertex);
            if (m_dim >= 2)
            {
                vElmtGeometry->LinkEndChild(vEdge);
            }
            if (m_dim >= 3)
            {
                vElmtGeometry->LinkEndChild(vFace);
            }
            vElmtGeometry->LinkEndChild(vElement);
            if (m_dim >= 2)
            {
                vElmtGeometry->LinkEndChild(vCurved);
            }
            vElmtGeometry->LinkEndChild(vComposite);
            vElmtGeometry->LinkEndChild(vDomain);

            pNektar->LinkEndChild(vElmtGeometry);

            if (pSession->DefinesElement("Nektar/Conditions"))
            {
                std::set<int> vBndRegionIdList;
                TiXmlElement* vConditions    = new TiXmlElement(*pSession->GetElement("Nektar/Conditions"));
                TiXmlElement* vBndRegions    = vConditions->FirstChildElement("BOUNDARYREGIONS");
                TiXmlElement* vBndConditions = vConditions->FirstChildElement("BOUNDARYCONDITIONS");
                TiXmlElement* vItem;

                if (vBndRegions)
                {
                    TiXmlElement* vNewBndRegions = new TiXmlElement("BOUNDARYREGIONS");
                    vItem = vBndRegions->FirstChildElement();
                    while (vItem)
                    {
                        std::string vSeqStr = vItem->FirstChild()->ToText()->Value();
                        std::string::size_type indxBeg = vSeqStr.find_first_of('[') + 1;
                        std::string::size_type indxEnd = vSeqStr.find_last_of(']') - 1;
                        vSeqStr = vSeqStr.substr(indxBeg, indxEnd - indxBeg + 1);
                        std::vector<unsigned int> vSeq;
                        ParseUtils::GenerateSeqVector(vSeqStr.c_str(), vSeq);
                        std::string vListStr;
                        bool comma = false;
                        for (unsigned int i = 0; i < vSeq.size(); ++i)
                        {
                            if (vComposites.find(vSeq[i]) != vComposites.end())
                            {
                                vListStr += comma ? "," : "";
                                comma = true;
                                vListStr += boost::lexical_cast<std::string>(vSeq[i]);
                            }
                        }
                        int p = atoi(vItem->Attribute("ID"));
                        if (vListStr.length() == 0)
                        {
                            vBndRegions->RemoveChild(vItem);
                        }
                        else
                        {
                            vListStr = "C[" + vListStr + "]";
                            TiXmlText* vList = new TiXmlText(vListStr);
                            TiXmlElement* vNewElement = new TiXmlElement("B");
                            vNewElement->SetAttribute("ID", p);
                            vNewElement->LinkEndChild(vList);
                            vNewBndRegions->LinkEndChild(vNewElement);
                            vBndRegionIdList.insert(p);
                        }

                        // Store original order of boundary region.
                        m_bndRegOrder[p] = vSeq;
                        
                        vItem = vItem->NextSiblingElement();
                    }
                    vConditions->ReplaceChild(vBndRegions, *vNewBndRegions);
                }

                if (vBndConditions)
                {
                    vItem = vBndConditions->FirstChildElement();
                    while (vItem)
                    {
                        std::set<int>::iterator x;
                        if ((x = vBndRegionIdList.find(atoi(vItem->Attribute("REF")))) != vBndRegionIdList.end())
                        {
                            vItem->SetAttribute("REF", *x);
                        }
                        else
                        {
                            vBndConditions->RemoveChild(vItem);
                        }
                        vItem = vItem->NextSiblingElement();
                    }
                }
                pNektar->LinkEndChild(vConditions);
            }

            // Distribute other sections of the XML to each process as is.
            TiXmlElement* vSrc = pSession->GetElement("Nektar")
                                                    ->FirstChildElement();
            while (vSrc)
            {
                std::string vName = boost::to_upper_copy(vSrc->ValueStr());
                if (vName != "GEOMETRY" && vName != "CONDITIONS")
                {
                    pNektar->LinkEndChild(new TiXmlElement(*vSrc));
                }
                vSrc = vSrc->NextSiblingElement();
            }
        }

        void MeshPartition::GetElementIDs(const int procid, std::vector<unsigned int> &elmtid)
        {
            BoostVertexIterator    vertit, vertit_end;

            ASSERTL0(procid < m_localPartition.size(),"procid is less than the number of partitions");
            
            // Populate lists of elements, edges and vertices required.
            for ( boost::tie(vertit, vertit_end) = boost::vertices(m_localPartition[procid]);
                  vertit != vertit_end;
                  ++vertit)
            {
                elmtid.push_back(m_meshElements[m_localPartition[procid][*vertit].id].id);
            }
        }
    }
}
