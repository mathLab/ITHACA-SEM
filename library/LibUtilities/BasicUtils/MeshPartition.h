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
#ifndef NEKTAR_LIBUTILITIES_BASICUTILS_MESHPARTITION_H
#define NEKTAR_LIBUTILITIES_BASICUTILS_MESHPARTITION_H

#include <boost/graph/subgraph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/MeshEntities.hpp>

class TiXmlElement;

namespace Nektar
{
    namespace LibUtilities
    {
        class MeshPartition;
        class SessionReader;
        typedef boost::shared_ptr<SessionReader> SessionReaderSharedPtr;
        typedef std::map<int, std::vector<unsigned int> > CompositeOrdering;
        typedef std::map<int, std::vector<unsigned int> > BndRegionOrdering;

        /// Datatype of the NekFactory used to instantiate classes derived from
        /// the EquationSystem class.
        typedef LibUtilities::NekFactory< std::string, MeshPartition, const SessionReaderSharedPtr& > MeshPartitionFactory;

        LIB_UTILITIES_EXPORT MeshPartitionFactory& GetMeshPartitionFactory();

        class MeshPartition
        {

        public:
            LIB_UTILITIES_EXPORT MeshPartition(const SessionReaderSharedPtr& pSession);
            LIB_UTILITIES_EXPORT virtual ~MeshPartition();

            LIB_UTILITIES_EXPORT void PartitionMesh(int nParts, bool shared = false);
            LIB_UTILITIES_EXPORT void WriteLocalPartition(
                    SessionReaderSharedPtr& pSession);
            LIB_UTILITIES_EXPORT void WriteAllPartitions(
                    SessionReaderSharedPtr& pSession);

            LIB_UTILITIES_EXPORT void PrintPartInfo(std::ostream &out);
            LIB_UTILITIES_EXPORT void GetCompositeOrdering(
                    CompositeOrdering &composites);
            LIB_UTILITIES_EXPORT void GetBndRegionOrdering(
                    BndRegionOrdering &composites);

            LIB_UTILITIES_EXPORT void GetElementIDs(const int procid,
                                                    std::vector<unsigned int> &tmp);

        private:
            struct MeshEntity
            {
                int id;
                char type;
                std::vector<unsigned int> list;
            };
            
            struct MeshFace
            {
                int id;
                std::vector<int> edgeList;
            };
            
            struct MeshElement
            {
                int id;
                char type;
                std::vector<int> list;
            };
            
            struct MeshCurved
            {
                int id;
                std::string entitytype;
                int entityid;
                std::string type;
                int npoints;
                std::string data;
                int ptid;
                int ptoffset;
            };

            struct MeshComposite
            {
                int id;
                char type;
                std::vector<int> list;
            };

            bool m_isCompressed; // Idenfity if input is compressed and if so set output to be compressed
            typedef std::pair<std::string, int> MeshCurvedKey;
            typedef std::vector<unsigned int>   MultiWeight;

            // Element in a mesh
            struct GraphVertexProperties
            {
                int id;             ///< Universal ID of the vertex
                int partition;      ///< Index of the partition to which it belongs
                MultiWeight weight; ///< Weightings to this graph vertex
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
                        boost::setS,
                        boost::vecS,
                        boost::undirectedS,
                        GraphVertexProperties,
                        boost::property<
                            boost::edge_index_t,
                            unsigned int,
                            GraphEdgeProperties
                        >
                    > BoostGraph;

            // Use induced subgraphs to manage the partitions
            typedef boost::subgraph< BoostGraph > BoostSubGraph;

            typedef boost::graph_traits<
                        BoostGraph
                    >::vertex_descriptor BoostVertex;
            typedef boost::graph_traits<
                        BoostGraph
                    >::edge_descriptor BoostEdge;
            typedef boost::graph_traits<
                        BoostGraph
                    >::edge_iterator BoostEdgeIterator;
            typedef boost::graph_traits<
                        BoostGraph
                    >::vertex_iterator BoostVertexIterator;
            typedef boost::graph_traits<
                        BoostGraph
                    >::adjacency_iterator BoostAdjacencyIterator;

            typedef std::vector<unsigned int>       NumModes;
            typedef std::map<std::string, NumModes> NummodesPerField;

            int                                 m_dim;
            int                                 m_numFields;

            std::map<int, MeshVertex>           m_meshVertices;
            std::map<int, MeshEntity>           m_meshEdges;
            std::map<int, MeshEntity>           m_meshFaces;
            std::map<int, MeshEntity>           m_meshElements;
            std::map<MeshCurvedKey, MeshCurved> m_meshCurved;
            std::map<int, MeshCurvedPts>        m_meshCurvedPts;
            std::map<int, MeshEntity>           m_meshComposites;
            std::vector<unsigned int>           m_domain;
            std::map<std::string, std::string>  m_vertexAttributes;

            // hierarchial mapping: composite id -> field name -> integer list
            // of directional nummodes described by expansion type clause.
            std::map<int, NummodesPerField>     m_expansions;

            std::map<std::string, int>          m_fieldNameToId;
            std::map<int, MultiWeight>          m_vertWeights;

            BndRegionOrdering                   m_bndRegOrder;

            BoostSubGraph                       m_mesh;
            std::vector<BoostSubGraph>          m_localPartition;

            CommSharedPtr                       m_comm;

            bool                                m_weightingRequired;
            bool                                m_shared;

            void ReadExpansions(const SessionReaderSharedPtr& pSession);
            void ReadGeometry(const SessionReaderSharedPtr& pSession);
            void ReadConditions(const SessionReaderSharedPtr& pSession);
            void WeightElements();
            void CreateGraph(BoostSubGraph& pGraph);
            void PartitionGraph(BoostSubGraph& pGraph,
                                int nParts,
                                std::vector<BoostSubGraph>& pLocalPartition);

            virtual void PartitionGraphImpl(
                    int&                              nVerts,
                    int&                              nVertConds,
                    Nektar::Array<Nektar::OneD, int>& xadj,
                    Nektar::Array<Nektar::OneD, int>& adjcy,
                    Nektar::Array<Nektar::OneD, int>& vertWgt,
                    Nektar::Array<Nektar::OneD, int>& vertSize,
                    int&                              nparts,
                    int&                              volume,
                    Nektar::Array<Nektar::OneD, int>& part) = 0;

            void OutputPartition(SessionReaderSharedPtr& pSession, BoostSubGraph& pGraph, TiXmlElement* pGeometry);
            void CheckPartitions(int nParts, Array<OneD, int> &pPart);
            int CalculateElementWeight(char elmtType, bool bndWeight, int na, int nb, int nc);
        };

        typedef boost::shared_ptr<MeshPartition> MeshPartitionSharedPtr;
    }
}

#endif
