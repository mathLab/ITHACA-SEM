///////////////////////////////////////////////////////////////////////////////
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

#include <SpatialDomains/MeshPartition.h>

#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

#include <tinyxml.h>

#include <LibUtilities/BasicUtils/FileSystem.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <SpatialDomains/MeshGraph.h>

#include <LibUtilities/Foundations/Foundations.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/detail/edge.hpp>

namespace Nektar
{
namespace SpatialDomains
{

SPATIAL_DOMAINS_EXPORT MeshPartitionFactory &GetMeshPartitionFactory()
{
    static MeshPartitionFactory instance;
    return instance;
}

MeshPartition::MeshPartition(const LibUtilities::SessionReaderSharedPtr session,
                             const MeshGraphSharedPtr m)
    :  m_dim(m->GetMeshDimension()), m_numFields(0), m_session(session),
       m_meshgraph(m), m_fieldNameToId(), m_comm(session->GetComm()),
       m_weightingRequired(false), m_weightBnd(false), m_weightDofs(false)
{
    // leave the meshpartition method of reading expansions and conditions
    ReadConditions();
    ReadExpansions();
}

MeshPartition::~MeshPartition()
{
}

void MeshPartition::PartitionMesh(int nParts, bool shared, bool overlapping)
{
    ASSERTL0(m_meshgraph->GetNumElements() >= nParts,
             "Too few elements for this many processes.");
    m_shared = shared;

    TransferElements();

    if (m_weightingRequired)
    {
        WeightElements();
    }
    CreateGraph();

    PartitionGraph(nParts, overlapping);
}

void MeshPartition::WriteLocalPartition()
{
    int rank              = m_comm->GetRowComm()->GetRank();
    BoostSubGraph &pGraph = m_localPartition[rank];
    unsigned int part     = pGraph[*boost::vertices(pGraph).first].partition;

    vector<unsigned int> parts;
    parts.push_back(part);
    vector<set<unsigned int>> elementLists;
    vector<unsigned int> ell;
    set<unsigned int> elementList;

    GetElementIDs(rank, ell);

    elementList.insert(ell.begin(), ell.end());
    elementLists.push_back(elementList);

    m_meshgraph->WriteGeometry(m_session->GetSessionName(),
                               elementLists, parts);
}

void MeshPartition::WriteAllPartitions()
{
    vector<unsigned int> parts;
    vector<set<unsigned int>> elementLists;

    for (int i = 0; i < m_localPartition.size(); ++i)
    {
        BoostSubGraph &pGraph = m_localPartition[i];
        unsigned int part = pGraph[*boost::vertices(pGraph).first].partition;

        parts.push_back(part);

        vector<unsigned int> ell;

        GetElementIDs(i, ell);

        set<unsigned int> elementList;

        elementList.insert(ell.begin(), ell.end());

        elementLists.push_back(elementList);
    }

    m_meshgraph->WriteGeometry(m_session->GetSessionName(),
                               elementLists, parts);
}

void MeshPartition::TransferElements()
{
    // elements are needed so process them from MeshGraph
    switch (m_dim)
    {
        case 1:
        {
            SegGeomMap &segGeoms = m_meshgraph->GetAllSegGeoms();
            for (auto &i : segGeoms)
            {
                MeshEntity e;
                e.id = i.first;
                e.list.push_back(i.second->GetVertex(0)->GetGlobalID());
                e.list.push_back(i.second->GetVertex(1)->GetGlobalID());
                m_elements[e.id] = e;
            }
        }
        break;
        case 2:
        {
            TriGeomMap &triGeoms = m_meshgraph->GetAllTriGeoms();
            for (auto &i : triGeoms)
            {
                MeshEntity e;
                e.id = i.first;
                e.list.push_back(i.second->GetEdge(0)->GetGlobalID());
                e.list.push_back(i.second->GetEdge(1)->GetGlobalID());
                e.list.push_back(i.second->GetEdge(2)->GetGlobalID());
                m_elements[e.id] = e;
            }
            QuadGeomMap &quadGeoms = m_meshgraph->GetAllQuadGeoms();
            for (auto &i : quadGeoms)
            {
                MeshEntity e;
                e.id = i.first;
                e.list.push_back(i.second->GetEdge(0)->GetGlobalID());
                e.list.push_back(i.second->GetEdge(1)->GetGlobalID());
                e.list.push_back(i.second->GetEdge(2)->GetGlobalID());
                e.list.push_back(i.second->GetEdge(3)->GetGlobalID());
                m_elements[e.id] = e;
            }
        }
        break;
        case 3:
        {
            TetGeomMap &tetGeoms = m_meshgraph->GetAllTetGeoms();
            for (auto &i : tetGeoms)
            {
                MeshEntity e;
                e.id = i.first;
                e.list.push_back(i.second->GetFace(0)->GetGlobalID());
                e.list.push_back(i.second->GetFace(1)->GetGlobalID());
                e.list.push_back(i.second->GetFace(2)->GetGlobalID());
                e.list.push_back(i.second->GetFace(3)->GetGlobalID());
                m_elements[e.id] = e;
            }
            PyrGeomMap &pyrGeoms = m_meshgraph->GetAllPyrGeoms();
            for (auto &i : pyrGeoms)
            {
                MeshEntity e;
                e.id = i.first;
                e.list.push_back(i.second->GetFace(0)->GetGlobalID());
                e.list.push_back(i.second->GetFace(1)->GetGlobalID());
                e.list.push_back(i.second->GetFace(2)->GetGlobalID());
                e.list.push_back(i.second->GetFace(3)->GetGlobalID());
                e.list.push_back(i.second->GetFace(4)->GetGlobalID());
                m_elements[e.id] = e;
            }
            PrismGeomMap &prismGeoms = m_meshgraph->GetAllPrismGeoms();
            for (auto &i : prismGeoms)
            {
                MeshEntity e;
                e.id = i.first;
                e.list.push_back(i.second->GetFace(0)->GetGlobalID());
                e.list.push_back(i.second->GetFace(1)->GetGlobalID());
                e.list.push_back(i.second->GetFace(2)->GetGlobalID());
                e.list.push_back(i.second->GetFace(3)->GetGlobalID());
                e.list.push_back(i.second->GetFace(4)->GetGlobalID());
                m_elements[e.id] = e;
            }
            HexGeomMap &hexGeoms = m_meshgraph->GetAllHexGeoms();
            for (auto &i : hexGeoms)
            {
                MeshEntity e;
                e.id = i.first;
                e.list.push_back(i.second->GetFace(0)->GetGlobalID());
                e.list.push_back(i.second->GetFace(1)->GetGlobalID());
                e.list.push_back(i.second->GetFace(2)->GetGlobalID());
                e.list.push_back(i.second->GetFace(3)->GetGlobalID());
                e.list.push_back(i.second->GetFace(4)->GetGlobalID());
                e.list.push_back(i.second->GetFace(5)->GetGlobalID());
                m_elements[e.id] = e;
            }
        }
        break;
    }
}

CompositeOrdering MeshPartition::GetCompositeOrdering()
{
    map<int, CompositeSharedPtr> &cmap = m_meshgraph->GetComposites();

    CompositeOrdering ret;

    for (auto &i : cmap)
    {
        vector<unsigned int> list;
        for (int j = 0; j < i.second->m_geomVec.size(); j++)
        {
            list.push_back(i.second->m_geomVec[j]->GetGlobalID());
        }
        ret[i.first] = list;
    }

    return ret;
}

BndRegionOrdering MeshPartition::GetBndRegionOrdering()
{
    return m_meshgraph->GetBndRegionOrdering();
}

void MeshPartition::ReadExpansions()
{
    // Find the Expansions tag
    TiXmlElement *expansionTypes = m_session->GetElement("Nektar/Expansions");

    // Find the Expansion type
    TiXmlElement *expansion = expansionTypes->FirstChildElement();
    std::string expType     = expansion->Value();

    /// Expansiontypes will contain plenty of data,
    /// where relevant at this stage are composite
    /// ID(s) that this expansion type describes,
    /// nummodes and a list of fields that this
    /// expansion relates to. If this does not exist
    /// the variable is only set to "DefaultVar".

    if (expType == "E")
    {
        while (expansion)
        {
            std::vector<unsigned int> composite;
            std::vector<unsigned int> nummodes;
            std::vector<std::string> fieldName;

            const char *nModesStr = expansion->Attribute("NUMMODES");
            ASSERTL0(nModesStr,
                     "NUMMODES was not defined in EXPANSION section of input");
            std::string numModesStr = nModesStr;
            bool valid = ParseUtils::GenerateVector(numModesStr.c_str(),
                                                           nummodes);
            ASSERTL0(valid, "Unable to correctly parse the number of modes.");

            if (nummodes.size() == 1)
            {
                for (int i = 1; i < m_dim; i++)
                {
                    nummodes.push_back(nummodes[0]);
                }
            }
            ASSERTL0(nummodes.size() == m_dim,
                     "Number of modes should match mesh dimension");

            const char *fStr = expansion->Attribute("FIELDS");
            if (fStr)
            {
                std::string fieldStr = fStr;
                bool valid           = ParseUtils::GenerateVector(
                    fieldStr.c_str(), fieldName);
                ASSERTL0(valid, "Unable to correctly parse the field string in "
                                "ExpansionTypes.");

                for (int i = 0; i < fieldName.size(); ++i)
                {
                    if (m_fieldNameToId.count(fieldName[i]) == 0)
                    {
                        int k                         = m_fieldNameToId.size();
                        m_fieldNameToId[fieldName[i]] = k;
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
                    ASSERTL0(
                        k == 0,
                        "Omitting field variables and explicitly listing "
                        "them in different ExpansionTypes is wrong practise");

                    m_fieldNameToId["DefaultVar"] = k;
                    m_numFields++;
                }
            }

            std::string compositeStr = expansion->Attribute("COMPOSITE");
            ASSERTL0(compositeStr.length() > 3,
                     "COMPOSITE must be specified in expansion definition");
            int beg = compositeStr.find_first_of("[");
            int end = compositeStr.find_first_of("]");
            std::string compositeListStr =
                compositeStr.substr(beg + 1, end - beg - 1);
            bool parseGood = ParseUtils::GenerateSeqVector(
                compositeListStr.c_str(), composite);
            ASSERTL0(parseGood && !composite.empty(),
                     (std::string("Unable to read composite index range: ") +
                      compositeListStr)
                         .c_str());

            // construct mapping (elmt id, field name) -> nummodes
            map<int, CompositeSharedPtr> &compMap =
                m_meshgraph->GetComposites();
            for (auto &i : compMap)
            {
                for (int j = 0; j < fieldName.size(); j++)
                {
                    for (int k = 0; k < i.second->m_geomVec.size(); k++)
                    {
                        int elid = i.second->m_geomVec[k]->GetGlobalID();
                        m_expansions[elid][fieldName[j]] = nummodes;
                        m_shape[elid] = i.second->m_geomVec[k]->GetShapeType();
                    }
                }
            }

            expansion = expansion->NextSiblingElement("E");
        }
    }
    else if (expType == "F")
    {
        ASSERTL0(expansion->Attribute("FILE"),
                 "Attribute FILE expected for type F expansion");
        std::string filenameStr = expansion->Attribute("FILE");
        ASSERTL0(!filenameStr.empty(),
                 "A filename must be specified for the FILE "
                 "attribute of expansion");

        // Create fieldIO object to load file
        //    need a serial communicator to avoid problems with
        //    shared file system
        LibUtilities::CommSharedPtr comm =
            LibUtilities::GetCommFactory().CreateInstance("Serial", 0, 0);
        std::string iofmt  = LibUtilities::FieldIO::GetFileType(filenameStr, comm);
        LibUtilities::FieldIOSharedPtr f = LibUtilities::GetFieldIOFactory().CreateInstance(
            iofmt, comm, m_session->GetSharedFilesystem());

        // Load field definitions from file
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> fielddefs;
        f->Import(filenameStr, fielddefs);

        // Parse field definitions
        for (int i = 0; i < fielddefs.size(); ++i)
        {
            // Name of fields
            for (int j = 0; j < fielddefs[i]->m_fields.size(); ++j)
            {
                std::string fieldName = fielddefs[i]->m_fields[j];
                if (m_fieldNameToId.count(fieldName) == 0)
                {
                    int k                      = m_fieldNameToId.size();
                    m_fieldNameToId[fieldName] = k;
                    m_numFields++;
                }
            }
            // Number of modes and shape for each element
            int numHomoDir = fielddefs[i]->m_numHomogeneousDir;
            int cnt        = 0;
            for (int j = 0; j < fielddefs[i]->m_elementIDs.size(); ++j)
            {
                int elid = fielddefs[i]->m_elementIDs[j];
                std::vector<unsigned int> nummodes;
                for (int k = 0; k < m_dim; k++)
                {
                    nummodes.push_back(fielddefs[i]->m_numModes[cnt++]);
                }
                if (fielddefs[i]->m_uniOrder)
                {
                    cnt = 0;
                }
                else
                {
                    cnt += numHomoDir;
                }
                for (int k = 0; k < fielddefs[i]->m_fields.size(); k++)
                {
                    std::string fieldName         = fielddefs[i]->m_fields[k];
                    m_expansions[elid][fieldName] = nummodes;
                }
                switch (fielddefs[i]->m_shapeType)
                {
                    case LibUtilities::eSegment:
                    {
                        m_shape[elid] = 'S';
                        break;
                    }
                    case LibUtilities::eTriangle:
                    {
                        m_shape[elid] = 'T';
                        break;
                    }
                    case LibUtilities::eQuadrilateral:
                    {
                        m_shape[elid] = 'Q';
                        break;
                    }
                    case LibUtilities::eTetrahedron:
                    {
                        m_shape[elid] = 'A';
                        break;
                    }
                    case LibUtilities::ePyramid:
                    {
                        m_shape[elid] = 'R';
                        break;
                    }
                    case LibUtilities::ePrism:
                    {
                        m_shape[elid] = 'P';
                        break;
                    }
                    case LibUtilities::eHexahedron:
                    {
                        m_shape[elid] = 'H';
                        break;
                    }
                    default:
                        ASSERTL0(false, "Shape not recognized.");
                        break;
                }
            }
        }
    }
    else
    {
        ASSERTL0(false,
                 "Expansion type not defined or not supported at the moment");
    }
}

void MeshPartition::PrintPartInfo(std::ostream &out)
{
    int nElmt = boost::num_vertices(m_graph);
    int nPart = m_localPartition.size();

    out << "# Partition information:" << std::endl;
    out << "# No. elements  : " << nElmt << std::endl;
    out << "# No. partitions: " << nPart << std::endl;
    out << "# ID  nElmt  nLocDof  nBndDof" << std::endl;

    BoostVertexIterator vertit, vertit_end;
    std::vector<int> partElmtCount(nPart, 0);
    std::vector<int> partLocCount(nPart, 0);
    std::vector<int> partBndCount(nPart, 0);

    std::map<int, int> elmtSizes;
    std::map<int, int> elmtBndSizes;

    for (std::map<int, NummodesPerField>::iterator expIt = m_expansions.begin();
         expIt != m_expansions.end(); ++expIt)
    {
        int elid             = expIt->first;
        NummodesPerField npf = expIt->second;

        for (NummodesPerField::iterator it = npf.begin(); it != npf.end(); ++it)
        {
            ASSERTL0(it->second.size() == m_dim,
                     " Number of directional"
                     " modes in expansion spec for element id = " +
                         boost::lexical_cast<std::string>(elid) +
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

            elmtSizes[elid] =
                CalculateElementWeight(m_shape[elid], false, na, nb, nc);
            elmtBndSizes[elid] =
                CalculateElementWeight(m_shape[elid], true, na, nb, nc);
        }
    }

    for (boost::tie(vertit, vertit_end) = boost::vertices(m_graph);
         vertit != vertit_end; ++vertit)
    {
        int partId = m_graph[*vertit].partition;
        partElmtCount[partId]++;
        partLocCount[partId] += elmtSizes[m_graph[*vertit].id];
        partBndCount[partId] += elmtBndSizes[m_graph[*vertit].id];
    }

    for (int i = 0; i < nPart; ++i)
    {
        out << i << " " << partElmtCount[i] << " " << partLocCount[i] << " "
            << partBndCount[i] << std::endl;
    }
}

void MeshPartition::ReadConditions()
{
    if (!m_session->DefinesElement("Nektar/Conditions/SolverInfo"))
    {
        // No SolverInfo = no change of default action to weight
        // mesh graph.
        return;
    }

    TiXmlElement *solverInfoElement =
        m_session->GetElement("Nektar/Conditions/SolverInfo");

    TiXmlElement *solverInfo = solverInfoElement->FirstChildElement("I");
    ASSERTL0(solverInfo, "Cannot read SolverInfo tags");

    while (solverInfo)
    {
        // read the property name
        ASSERTL0(solverInfo->Attribute("PROPERTY"),
                 "Missing PROPERTY attribute in solver info "
                 "section. ");
        std::string solverProperty = solverInfo->Attribute("PROPERTY");
        ASSERTL0(!solverProperty.empty(),
                 "Solver info properties must have a non-empty "
                 "name. ");
        // make sure that solver property is capitalised
        std::string solverPropertyUpper = boost::to_upper_copy(solverProperty);

        // read the value
        ASSERTL0(solverInfo->Attribute("VALUE"),
                 "Missing VALUE attribute in solver info section. ");
        std::string solverValue = solverInfo->Attribute("VALUE");
        ASSERTL0(!solverValue.empty(),
                 "Solver info properties must have a non-empty value");
        // make sure that property value is capitalised
        std::string propertyValueUpper = boost::to_upper_copy(solverValue);

        if (solverPropertyUpper == "WEIGHTPARTITIONS")
        {
            if (propertyValueUpper == "DOF")
            {
                m_weightingRequired = true;
                m_weightDofs        = true;
            }
            else if (propertyValueUpper == "BOUNDARY")
            {
                m_weightingRequired = true;
                m_weightBnd         = true;
            }
            else if (propertyValueUpper == "BOTH")
            {
                m_weightingRequired = true;
                m_weightDofs        = true;
                m_weightBnd         = true;
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
    std::map<int, MeshEntity>::iterator eIt;
    for (eIt = m_elements.begin(); eIt != m_elements.end(); ++eIt)
    {
        m_vertWeights[eIt->first]    = weight;
        m_vertBndWeights[eIt->first] = weight;
        m_edgeWeights[eIt->first]    = weight;
    }

    for (std::map<int, NummodesPerField>::iterator expIt = m_expansions.begin();
         expIt != m_expansions.end(); ++expIt)
    {
        int elid             = expIt->first;
        NummodesPerField npf = expIt->second;

        for (NummodesPerField::iterator it = npf.begin(); it != npf.end(); ++it)
        {
            ASSERTL0(it->second.size() == m_dim,
                     " Number of directional"
                     " modes in expansion spec for element id = " +
                         boost::lexical_cast<std::string>(elid) +
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

            m_vertWeights[elid][m_fieldNameToId[it->first]] =
                CalculateElementWeight(m_shape[elid], false, na, nb, nc);
            m_vertBndWeights[elid][m_fieldNameToId[it->first]] =
                CalculateElementWeight(m_shape[elid], true, na, nb, nc);
            m_edgeWeights[elid][m_fieldNameToId[it->first]] =
                CalculateEdgeWeight(m_shape[elid], na, nb, nc);
        }
    } // for i
}

void MeshPartition::CreateGraph()
{
    // Maps edge/face to first mesh element id.
    // On locating second mesh element id, graph edge is created instead.
    std::map<int, int> vGraphEdges;
    map<int, MeshEntity>::iterator eIt;
    int vcnt = 0;

    for (eIt = m_elements.begin(); eIt != m_elements.end(); ++eIt, ++vcnt)
    {
        BoostVertex v       = boost::add_vertex(m_graph);
        m_graph[v].id        = eIt->first;
        m_graph[v].partition = 0;

        if (m_weightingRequired)
        {
            m_graph[v].weight     = m_vertWeights[eIt->first];
            m_graph[v].bndWeight  = m_vertBndWeights[eIt->first];
            m_graph[v].edgeWeight = m_edgeWeights[eIt->first];
        }

        // Process element entries and add graph edges
        for (unsigned j = 0; j < eIt->second.list.size(); ++j)
        {
            int eId = eIt->second.list[j];

            // Look to see if we've examined this edge/face before
            // If so, we've got both graph vertices so add edge
            if (vGraphEdges.find(eId) != vGraphEdges.end())
            {
                BoostEdge e =
                    boost::add_edge(vcnt, vGraphEdges[eId], m_graph).first;
                m_graph[e].id = vcnt;
            }
            else
            {
                vGraphEdges[eId] = vcnt;
            }
        }
    }
}

/**
 * @brief Partition the graph.
 *
 * This routine partitions the graph @p pGraph into @p nParts, producing
 * subgraphs that are populated in @p pLocalPartition. If the @p
 * overlapping option is set (which is used for post-processing
 * purposes), the resulting partitions are extended to cover
 * neighbouring elements by additional vertex on the dual graph, which
 * produces overlapping partitions (i.e. the intersection of two
 * connected partitions is non-empty).
 *
 * @param nParts           Number of partitions.
 * @param pLocalPartition  Vector of sub-graphs representing each
 * @param overlapping      True if resulting partitions should overlap.
 */
void MeshPartition::PartitionGraph(int nParts, bool overlapping)
{
    int i;
    int nGraphVerts = boost::num_vertices(m_graph);
    int nGraphEdges = boost::num_edges(m_graph);

    int ncon = 1;
    if (m_weightDofs && m_weightBnd)
    {
        ncon = 2;
    }
    // Convert boost graph into CSR format
    BoostVertexIterator vertit, vertit_end;
    BoostAdjacencyIterator adjvertit, adjvertit_end;
    Array<OneD, int> part(nGraphVerts, 0);

    if (m_comm->GetRowComm()->TreatAsRankZero())
    {
        int acnt    = 0;
        int vcnt    = 0;
        int nWeight = ncon * nGraphVerts;

        Array<OneD, int> xadj(nGraphVerts + 1, 0);
        Array<OneD, int> adjncy(2 * nGraphEdges);
        Array<OneD, int> adjwgt(2 * nGraphEdges, 1);
        Array<OneD, int> vwgt(nWeight, 1);
        Array<OneD, int> vsize(nGraphVerts, 1);

        for (boost::tie(vertit, vertit_end) = boost::vertices(m_graph);
             vertit != vertit_end; ++vertit)
        {
            for (boost::tie(adjvertit, adjvertit_end) =
                     boost::adjacent_vertices(*vertit, m_graph);
                 adjvertit != adjvertit_end; ++adjvertit)
            {
                adjncy[acnt++] = *adjvertit;
                if (m_weightingRequired)
                {
                    adjwgt[acnt - 1] = m_graph[*vertit].edgeWeight[0];
                }
            }

            xadj[++vcnt] = acnt;

            if (m_weightingRequired)
            {
                int ccnt = 0;
                if (m_weightDofs)
                {
                    vwgt[ncon * (vcnt - 1) + ccnt] = m_graph[*vertit].weight[0];
                    ccnt++;
                }
                if (m_weightBnd)
                {
                    vwgt[ncon * (vcnt - 1) + ccnt] =
                        m_graph[*vertit].bndWeight[0];
                }
            }
        }

        // Call partitioner to partition graph
        int vol = 0;

        try
        {
            //////////////////////////////////////////////////////
            // On a cartesian communicator do mesh partiotion just on the first
            // column
            // so there is no doubt the partitions are all the same in all the
            // columns
            if (m_comm->GetColumnComm()->GetRank() == 0)
            {
                // Attempt partitioning.
                PartitionGraphImpl(nGraphVerts, ncon, xadj, adjncy, vwgt, vsize,
                                   adjwgt, nParts, vol, part);

                // Check the partitioner produced a valid partition and fix if
                // not.
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
                     "Error in calling graph partitioner.");
        }
    }
    else
    {
        m_comm->GetRowComm()->Recv(0, part);
    }

    // Create boost subgraph for this process's partitions;
    m_localPartition.resize(nParts);
    for (i = 0; i < nParts; ++i)
    {
        m_localPartition[i] = m_graph.create_subgraph();
    }

    // Populate subgraph
    i = 0;
    for (boost::tie(vertit, vertit_end) = boost::vertices(m_graph);
         vertit != vertit_end; ++vertit, ++i)
    {
        m_graph[*vertit].partition = part[i];
        boost::add_vertex(i, m_localPartition[part[i]]);
    }

    // If the overlapping option is set (for post-processing purposes),
    // add vertices that correspond to the neighbouring elements.
    if (overlapping)
    {
        for (boost::tie(vertit, vertit_end) = boost::vertices(m_graph);
             vertit != vertit_end; ++vertit)
        {
            for (boost::tie(adjvertit, adjvertit_end) =
                     boost::adjacent_vertices(*vertit, m_graph);
                 adjvertit != adjvertit_end; ++adjvertit)
            {
                if (part[*adjvertit] != part[*vertit])
                {
                    boost::add_vertex(*adjvertit,
                                      m_localPartition[part[*vertit]]);
                }
            }
        }
    }
}

void MeshPartition::CheckPartitions(int nParts, Array<OneD, int> &pPart)
{
    unsigned int i   = 0;
    unsigned int cnt = 0;
    bool valid       = true;

    // Check that every process has at least one element assigned
    for (i = 0; i < nParts; ++i)
    {
        cnt = std::count(pPart.begin(), pPart.end(), i);
        if (cnt == 0)
        {
            valid = false;
        }
    }

    // If the graph partitioner produced an invalid partition, repartition
    // naively.  Elements are assigned to processes in a round-robin fashion.
    // It is assumed that graph partitioner failure only occurs when the number
    // of elements is approx. the number of processes, so this approach should
    // not be too inefficient communication-wise.
    if (!valid)
    {
        for (i = 0; i < pPart.num_elements(); ++i)
        {
            pPart[i] = i % nParts;
        }
    }
}

void MeshPartition::GetElementIDs(const int procid,
                                  std::vector<unsigned int> &elmtid)
{
    BoostVertexIterator vertit, vertit_end;

    ASSERTL0(procid < m_localPartition.size(),
             "procid is less than the number of partitions");

    // Populate lists of elements, edges and vertices required.
    for (boost::tie(vertit, vertit_end) =
             boost::vertices(m_localPartition[procid]);
         vertit != vertit_end; ++vertit)
    {
        elmtid.push_back(m_elements[m_localPartition[procid][*vertit].id].id);
    }
}

int MeshPartition::CalculateElementWeight(char elmtType, bool bndWeight, int na,
                                          int nb, int nc)
{
    int weight = 0;

    switch (elmtType)
    {
        case 'A':
            weight = bndWeight
                         ? LibUtilities::StdTetData::getNumberOfBndCoefficients(
                               na, nb, nc)
                         : LibUtilities::StdTetData::getNumberOfCoefficients(
                               na, nb, nc);
            break;
        case 'R':
            weight =
                bndWeight
                    ? LibUtilities::StdPrismData::getNumberOfBndCoefficients(
                          na, nb, nc)
                    : LibUtilities::StdPrismData::getNumberOfCoefficients(
                          na, nb, nc);
            break;
        case 'H':
            weight = bndWeight
                         ? LibUtilities::StdHexData::getNumberOfBndCoefficients(
                               na, nb, nc)
                         : LibUtilities::StdHexData::getNumberOfCoefficients(
                               na, nb, nc);
            break;
        case 'P':
            weight = bndWeight
                         ? LibUtilities::StdPyrData::getNumberOfBndCoefficients(
                               na, nb, nc)
                         : LibUtilities::StdPyrData::getNumberOfCoefficients(
                               na, nb, nc);
            break;
        case 'Q':
            weight =
                bndWeight
                    ? LibUtilities::StdQuadData::getNumberOfBndCoefficients(na,
                                                                            nb)
                    : LibUtilities::StdQuadData::getNumberOfCoefficients(na,
                                                                         nb);
            break;
        case 'T':
            weight =
                bndWeight
                    ? LibUtilities::StdTriData::getNumberOfBndCoefficients(na,
                                                                           nb)
                    : LibUtilities::StdTriData::getNumberOfCoefficients(na, nb);
            break;
        case 'S':
            weight =
                bndWeight
                    ? LibUtilities::StdSegData::getNumberOfBndCoefficients(na)
                    : LibUtilities::StdSegData::getNumberOfCoefficients(na);
            break;
        case 'V':
            weight = 1;
            break;
        default:
            break;
    }

    return weight;
}

/**
 *     Calculate the number of modes needed for communication when
 *        in partition boundary, to be used as weighting for edges.
 *     Since we do not know exactly which face this refers to, assume
 *        the max order and quad face (for prisms) as arbitrary choices
 */
int MeshPartition::CalculateEdgeWeight(char elmtType, int na, int nb, int nc)
{
    int weight = 0;
    int n      = std::max(na, std::max(nb, nc));
    switch (elmtType)
    {
        case 'A':
            weight = LibUtilities::StdTriData::getNumberOfCoefficients(n, n);
            break;
        case 'R':
            weight = LibUtilities::StdQuadData::getNumberOfCoefficients(n, n);
            break;
        case 'H':
            weight = LibUtilities::StdQuadData::getNumberOfCoefficients(n, n);
            break;
        case 'P':
            weight = LibUtilities::StdQuadData::getNumberOfCoefficients(n, n);
            break;
        case 'Q':
        case 'T':
            weight = n;
            break;
        case 'S':
            weight = 1;
            break;
        default:
            break;
    }

    return weight;
}
}
}
