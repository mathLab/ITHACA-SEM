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

class TiXmlElement;

namespace Nektar
{
    namespace LibUtilities
    {
        class SessionReader;

        class MeshPartition
        {
        typedef boost::shared_ptr<SessionReader> SessionReaderSharedPtr;

        public:
            LIB_UTILITIES_EXPORT MeshPartition(const SessionReaderSharedPtr& pSession);
            LIB_UTILITIES_EXPORT ~MeshPartition();

            LIB_UTILITIES_EXPORT void PartitionMesh();
            LIB_UTILITIES_EXPORT void WriteLocalPartition(
                    SessionReaderSharedPtr& pSession);

        private:
            struct MeshEntity
            {
                int id;
                char type;
                std::vector<int> list;
            };

            struct MeshVertex
            {
                int id;
                NekDouble x;
                NekDouble y;
                NekDouble z;
            };

            struct MeshEdge
            {
                int id;
                int v0;
                int v1;
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
                int edgeid;
                std::string type;
                int npoints;
                std::string data;
            };

            struct MeshComposite
            {
                int id;
                char type;
                std::vector<int> list;
            };

            // Element in a mesh
            struct GraphVertexProperties
            {
                int id;         ///< Universal ID of the vertex
                int partition;  ///< Index of the partition to which it belongs
                int partid;     ///< Global ID of the vertex in the partition
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

            int                        m_dim;

            std::map<int, MeshVertex>  m_meshVertices;
            std::map<int, MeshEntity>  m_meshEdges;
            std::map<int, MeshEntity>  m_meshFaces;
            std::map<int, MeshEntity>  m_meshElements;
            std::map<int, MeshCurved>  m_meshCurved;
            std::map<int, MeshEntity>  m_meshComposites;
            std::vector<unsigned int>  m_domain;

            BoostSubGraph              m_mesh;
            BoostSubGraph              m_localPartition;

            CommSharedPtr              m_comm;

            void ReadMesh(const SessionReaderSharedPtr& pSession);
            void CreateGraph(BoostSubGraph& pGraph);
            void PartitionGraph(BoostSubGraph& pGraph,
                                BoostSubGraph& pLocalPartition);
            void OutputPartition(SessionReaderSharedPtr& pSession, BoostSubGraph& pGraph, TiXmlElement* pGeometry);
            void CheckPartitions(Array<OneD, int> &pPart);
        };

        typedef boost::shared_ptr<MeshPartition> MeshPartitionSharedPtr;
    }
}

#endif
