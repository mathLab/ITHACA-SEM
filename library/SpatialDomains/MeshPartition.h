////////////////////////////////////////////////////////////////////////////////
//
//  File: MeshPartition.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
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
#ifndef NEKTAR_SPATIALDOMAINS_MESHPARTITION_H
#define NEKTAR_SPATIALDOMAINS_MESHPARTITION_H

#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
#include <SpatialDomains/MeshEntities.hpp>
#include <boost/graph/adjacency_list.hpp>

class TiXmlElement;

namespace Nektar
{
namespace SpatialDomains
{
class MeshPartition;

typedef std::map<int, std::pair<LibUtilities::ShapeType, std::vector<int>>>
CompositeDescriptor;

/// Datatype of the NekFactory used to instantiate classes derived from
/// the EquationSystem class.
typedef LibUtilities::NekFactory<std::string, MeshPartition,
                                 const LibUtilities::SessionReaderSharedPtr,
                                 int, std::map<int, MeshEntity>,
                                 CompositeDescriptor>
    MeshPartitionFactory;

SPATIAL_DOMAINS_EXPORT MeshPartitionFactory &GetMeshPartitionFactory();

class MeshPartition
{

public:
    MeshPartition(const LibUtilities::SessionReaderSharedPtr session,
                  int                                        meshDim,
                  std::map<int, MeshEntity>                  element,
                  CompositeDescriptor                        compMap);
    virtual ~MeshPartition();

    SPATIAL_DOMAINS_EXPORT void PartitionMesh(
        int  nParts,
        bool shared      = false,
        bool overlapping = false,
        int  nLocal      = 0);

    SPATIAL_DOMAINS_EXPORT void PrintPartInfo(std::ostream &out);

    SPATIAL_DOMAINS_EXPORT void GetElementIDs(
        const int                  procid,
        std::vector<unsigned int> &tmp);

protected:
    typedef std::vector<unsigned int> MultiWeight;

    // Element in a mesh
    struct GraphVertexProperties
    {
        int id;             ///< Universal ID of the vertex
        int partition;      ///< Index of the partition to which it belongs
        MultiWeight weight; ///< Weightings to this graph vertex
        MultiWeight bndWeight;
        MultiWeight edgeWeight;
    };

    // Face/Edge/Vertex between two adjacent elements
    struct GraphEdgeProperties
    {
        int id;
        std::vector<MeshVertex> vertices;
        std::vector<MeshEdge> edges;
    };

    // Basic graph definition
    typedef boost::adjacency_list<
        boost::setS, boost::vecS, boost::undirectedS, GraphVertexProperties,
        boost::property<boost::edge_index_t, unsigned int, GraphEdgeProperties>>
        BoostGraph;

    typedef boost::graph_traits<BoostGraph>::vertex_descriptor BoostVertex;
    typedef boost::graph_traits<BoostGraph>::edge_descriptor BoostEdge;
    typedef boost::graph_traits<BoostGraph>::edge_iterator BoostEdgeIterator;
    typedef boost::graph_traits<BoostGraph>::vertex_iterator
        BoostVertexIterator;
    typedef boost::graph_traits<BoostGraph>::adjacency_iterator
        BoostAdjacencyIterator;

    typedef std::vector<unsigned int> NumModes;
    typedef std::map<std::string, NumModes> NummodesPerField;

    LibUtilities::SessionReaderSharedPtr m_session;

    int m_dim;
    int m_numFields;

    std::map<int, MeshEntity> m_elements;
    std::map<int, MeshEntity> m_ghostElmts;
    CompositeDescriptor m_compMap;

    // hierarchial mapping: elmt id -> field name -> integer list
    // of directional nummodes described by expansion type clause.
    std::map<int, NummodesPerField> m_expansions;

    // map of each elements shape
    std::map<int, LibUtilities::ShapeType> m_shape;

    std::map<std::string, int> m_fieldNameToId;
    std::map<int, MultiWeight> m_vertWeights;
    std::map<int, MultiWeight> m_vertBndWeights;
    std::map<int, MultiWeight> m_edgeWeights;

    BoostGraph m_graph;
    std::vector<std::vector<unsigned int>> m_localPartition;

    LibUtilities::CommSharedPtr m_comm;

    bool m_weightingRequired;
    bool m_weightBnd;
    bool m_weightDofs;
    bool m_shared;
    bool m_parallel;

    void ReadExpansions();
    void ReadConditions();
    void WeightElements();
    void CreateGraph();
    void PartitionGraph(int nParts, bool overlapping = false);

    virtual void PartitionGraphImpl(int &nVerts, int &nVertConds,
                                    Nektar::Array<Nektar::OneD, int> &xadj,
                                    Nektar::Array<Nektar::OneD, int> &adjcy,
                                    Nektar::Array<Nektar::OneD, int> &vertWgt,
                                    Nektar::Array<Nektar::OneD, int> &vertSize,
                                    Nektar::Array<Nektar::OneD, int> &edgeWgt,
                                    int &nparts, int &volume,
                                    Nektar::Array<Nektar::OneD, int> &part) = 0;

    void CheckPartitions(int nParts, Array<OneD, int> &pPart);
    int CalculateElementWeight(LibUtilities::ShapeType elmtType, bool bndWeight,
                               int na, int nb, int nc);
    int CalculateEdgeWeight(LibUtilities::ShapeType elmtType,
                            int na, int nb, int nc);
};

typedef std::shared_ptr<MeshPartition> MeshPartitionSharedPtr;
}
}

#endif
