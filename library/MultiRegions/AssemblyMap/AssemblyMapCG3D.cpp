///////////////////////////////////////////////////////////////////////////////
//
// File AssemblyMapCG3D.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: C0-continuous assembly mappings specific to 3D
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/AssemblyMap/AssemblyMapCG3D.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class AssemblyMapCG3D
         * Mappings are created for three possible global solution types:
         *  - Direct full matrix
         *  - Direct static condensation
         *  - Direct multi-level static condensation
         * In the latter case, mappings are created recursively for the
         * different levels of static condensation.
         *
         * These mappings are used by GlobalLinSys to generate the global
         * system.
         */

        /**
         *
         */
        AssemblyMapCG3D::AssemblyMapCG3D(
                const LibUtilities::SessionReaderSharedPtr &pSession):
            AssemblyMapCG(pSession)
        {
        }



        /**
         *
         */
        AssemblyMapCG3D::AssemblyMapCG3D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const int numLocalCoeffs,
                const ExpList &locExp,
                const Array<OneD, const ExpListSharedPtr> &bndCondExp,
                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
                                                                &bndConditions,
                const map<int,int>& periodicVerticesId,
                const map<int,int>& periodicEdgesId,
                const map<int,pair<int, StdRegions::Orientation> >& periodicFacesId):
            AssemblyMapCG(pSession)
        {
            SetUp3DExpansionC0ContMap(numLocalCoeffs,
                                      locExp,
                                      bndCondExp,
                                      bndConditions,
                                      periodicVerticesId,
                                      periodicEdgesId,
                                      periodicFacesId);

            CalculateBndSystemBandWidth();
            CalculateFullSystemBandWidth();
        }


        /**
         *
         */
        AssemblyMapCG3D::AssemblyMapCG3D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const int numLocalCoeffs,
                const ExpList &locExp):
            AssemblyMapCG(pSession)
        {
            SetUp3DExpansionC0ContMap(numLocalCoeffs, locExp);
            CalculateBndSystemBandWidth();
            CalculateFullSystemBandWidth();
        }


        /**
         *
         */
        AssemblyMapCG3D::~AssemblyMapCG3D()
        {
        }




        /**
         * Construction of the local->global map is achieved in several stages.
         * A mesh vertex, mesh edge and mesh face renumbering is constructed
         * in #vertReorderedGraphVertId, #edgeReorderedGraphVertId and
         * #faceReorderedGraphVertId
         *
         * The only unique identifiers of the vertices, edges and faces of the
         * mesh are the vertex id and the mesh id (stored in their corresponding
         * Geometry object).  However, setting up a global numbering based on
         * these id's will not lead to a suitable or optimal numbering. Mainly
         * because:
         *  - we want the Dirichlet DOF's to be listed first
         *  - we want an optimal global numbering of the remaining DOF's
         *    (strategy still need to be defined but can for example be:
         *    minimum bandwith or minimum fill-in of the resulting global
         *    system matrix)
         *
         * That's why the vertices, edges and faces will be rearranged. This is
         * done in the following way: The vertices, edges and faces of the mesh
         * are considered as vertices of a graph (in a computer science way)
         * (equivalently, they can also be considered as boundary degrees of
         * freedom, whereby all boundary modes of a single edge are considered
         * as a single DOF). We then will use algorithms to reorder these
         * graph-vertices (or boundary DOF's).
         *
         * We will use a boost graph object to store this graph the first
         * template parameter (=OutEdgeList) is chosen to be of type std::set
         * as in the set up of the adjacency, a similar edge might be created
         * multiple times.  And to prevent the definition of parallel edges,
         * we use std::set (=boost::setS) rather than std::vector
         * (=boost::vecS).
         *
         * Two different containers are used to store the graph vertex id's of
         * the different mesh vertices and edges. They are implemented as a STL
         * map such that the graph vertex id can later be retrieved by the
         * unique mesh vertex or edge id's which serve as the key of the map.
         *
         * Therefore, the algorithm proceeds as follows:
         */
        void AssemblyMapCG3D::SetUp3DExpansionC0ContMap(
            const int numLocalCoeffs,
            const ExpList &locExp,
            const Array<OneD, const ExpListSharedPtr> &bndCondExp,
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions,
            const map<int,int>& periodicVerticesId,
            const map<int,int>& periodicEdgesId,
            const map<int,pair<int, StdRegions::Orientation> >& periodicFacesId)
        {
            int i,j,k,l;
            int cnt = 0;
            int intDofCnt;
            int meshVertId;
            int meshVertId2;
            int meshEdgeId;
            int meshEdgeId2;
            int meshFaceId;
            int meshFaceId2;
            int globalId;
            int nEdgeInteriorCoeffs;
            int nFaceInteriorCoeffs;
            int firstNonDirGraphVertId;
            int nLocBndCondDofs = 0;
            int nLocDirBndCondDofs = 0;
            int graphVertId = 0;
            StdRegions::StdExpansion3DSharedPtr locExpansion;
            StdRegions::StdExpansion2DSharedPtr bndCondFaceExp;
            LibUtilities::BasisType             bType;
            StdRegions::Orientation             edgeOrient;
            StdRegions::Orientation             faceOrient;
            Array<OneD, unsigned int>           edgeInteriorMap;
            Array<OneD, int>                    edgeInteriorSign;
            Array<OneD, unsigned int>           faceInteriorMap;
            Array<OneD, int>                    faceInteriorSign;

            const StdRegions::StdExpansionVector &locExpVector = *(locExp.GetExp());

            m_signChange = false;
            //m_systemSingular = false;

            map<int,int> vertReorderedGraphVertId;
            map<int,int> edgeReorderedGraphVertId;
            map<int,int> faceReorderedGraphVertId;
            map<int,int>::iterator mapIt;
            map<int,int>::const_iterator mapConstIt;
            map<int,pair<int, StdRegions::Orientation> >::const_iterator mapFaceIt;

            bool systemSingular = true;

            /**
             * STEP 1: Order the Dirichlet vertices and edges first
             */
            for(i = 0; i < bndCondExp.num_elements(); i++)
            {
                // Check to see if any value on face has Dirichlet value.
                cnt = 0;
                for(j = 0; j < bndCondExp[i]->GetNumElmts(); j++)
                {
                    bndCondFaceExp = boost::dynamic_pointer_cast<
                        StdRegions::StdExpansion2D>(
                            bndCondExp[i]->GetExp(j));
                    if (bndConditions[i]->GetBoundaryConditionType() == 
                        SpatialDomains::eDirichlet)
                    {
                        for(k = 0; k < bndCondFaceExp->GetNverts(); k++)
                        {
                            meshVertId = (bndCondFaceExp->GetGeom2D())->GetVid(k);
                            if(vertReorderedGraphVertId.count(meshVertId) == 0)
                            {
                                vertReorderedGraphVertId[meshVertId] = graphVertId++;
                            }
                        }

                        for(k = 0; k < bndCondFaceExp->GetNedges(); k++)
                        {
                            meshEdgeId = (bndCondFaceExp->GetGeom2D())->GetEid(k);
                            if(edgeReorderedGraphVertId.count(meshEdgeId) == 0)
                            {
                                edgeReorderedGraphVertId[meshEdgeId] = graphVertId++;
                            }
                        }
                        meshFaceId = (bndCondFaceExp->GetGeom2D())->GetFid();
                        faceReorderedGraphVertId[meshFaceId] = graphVertId++;
                        nLocDirBndCondDofs += bndCondFaceExp->GetNcoeffs();
                    }
                    if (bndConditions[i]->GetBoundaryConditionType() !=
                        SpatialDomains::eNeumann)
                    {
                        systemSingular = false;
                    }
                    nLocBndCondDofs += bndCondFaceExp->GetNcoeffs();
                }
            }

            /**
             * STEP 1.5: Exchange Dirichlet mesh vertices between processes and
             * check for singular problems.
             */
            
            // Collate information on Dirichlet vertices from all processes
            int n = m_comm->GetSize();
            int p = m_comm->GetRank();

            Array<OneD, int> vertcounts (n, 0);
            Array<OneD, int> vertoffsets(n, 0);
            Array<OneD, int> edgecounts (n, 0);
            Array<OneD, int> edgeoffsets(n, 0);
            vertcounts[p] = vertReorderedGraphVertId.size();
            edgecounts[p] = edgeReorderedGraphVertId.size();
            m_comm->AllReduce(vertcounts, LibUtilities::ReduceSum);
            m_comm->AllReduce(edgecounts, LibUtilities::ReduceSum);

            for (i = 1; i < n; ++i)
            {
                vertoffsets[i] = vertoffsets[i-1] + vertcounts[i-1];
                edgeoffsets[i] = edgeoffsets[i-1] + edgecounts[i-1];
            }

            int nTotVerts = Vmath::Vsum(n,vertcounts,1);
            int nTotEdges = Vmath::Vsum(n,edgecounts,1);
            
            Array<OneD, int> vertlist(nTotVerts, 0);
            Array<OneD, int> edgelist(nTotEdges, 0);
            std::map<int, int>::iterator it;
            for (it  = vertReorderedGraphVertId.begin(), i = 0;
                 it != vertReorderedGraphVertId.end();
                 ++it, ++i)
            {
                vertlist[vertoffsets[p] + i] = it->first;
            }
            
            for (it  = edgeReorderedGraphVertId.begin(), i = 0;
                 it != edgeReorderedGraphVertId.end();
                 ++it, ++i)
            {
                edgelist[edgeoffsets[p] + i] = it->first;
            }
            m_comm->AllReduce(vertlist, LibUtilities::ReduceSum);
            m_comm->AllReduce(edgelist, LibUtilities::ReduceSum);
            
            int nExtraDirichlet = 0;
            map<int, int> extraDirVertIds, extraDirEdgeIds;
            
            // Ensure Dirchlet vertices are consistently recorded between
            // processes (e.g. Dirichlet region meets Neumann region across a
            // partition boundary requires vertex on partition to be Dirichlet).
            for (i = 0; i < n; ++i)
            {
                if (i == p)
                {
                    continue;
                }
                
                for(j = 0; j < locExpVector.size(); j++)
                {
                    locExpansion = boost::dynamic_pointer_cast<
                    StdRegions::StdExpansion3D>(
                            locExpVector[locExp.GetOffset_Elmt_Id(j)]);
                    
                    for(k = 0; k < locExpansion->GetNverts(); k++)
                    {
                        meshVertId = locExpansion->GetGeom3D()->GetVid(k);
                        if(vertReorderedGraphVertId.count(meshVertId) == 0)
                        {
                            for (l = 0; l < vertcounts[i]; ++l)
                            {
                                if (vertlist[vertoffsets[i]+l] == meshVertId)
                                {
                                    extraDirVertIds[meshVertId] = i;
                                    vertReorderedGraphVertId[meshVertId] = graphVertId++;
                                    nExtraDirichlet++;
                                }
                            }
                        }
                    }
                    
                    for(k = 0; k < locExpansion->GetNedges(); k++)
                    {
                        meshEdgeId = locExpansion->GetGeom3D()->GetEid(k);
                        if(edgeReorderedGraphVertId.count(meshEdgeId) == 0)
                        {
                            for (l = 0; l < edgecounts[i]; ++l)
                            {
                                if (edgelist[edgeoffsets[i]+l] == meshEdgeId)
                                {
                                    extraDirEdgeIds[meshEdgeId] = i;
                                    edgeReorderedGraphVertId[meshEdgeId] = graphVertId++;
                                    nExtraDirichlet += locExpansion->GetEdgeNcoeffs(k) - 2;
                                }
                            }
                        }
                    }
                }
            }

            for (i = 0; i < n; ++i)
            {
                vertcounts [i] = 0;
                vertoffsets[i] = 0;
                edgecounts [i] = 0;
                edgeoffsets[i] = 0;
            }
            
            vertcounts[p] = extraDirVertIds.size();
            edgecounts[p] = extraDirEdgeIds.size();
            m_comm->AllReduce(vertcounts, LibUtilities::ReduceSum);
            m_comm->AllReduce(edgecounts, LibUtilities::ReduceSum);
            nTotVerts = Vmath::Vsum(n, vertcounts, 1);
            nTotEdges = Vmath::Vsum(n, edgecounts, 1);
            
            vertoffsets[0] = edgeoffsets[0] = 0;
            
            for (i = 1; i < n; ++i)
            {
                vertoffsets[i] = vertoffsets[i-1] + vertcounts[i-1];
                edgeoffsets[i] = edgeoffsets[i-1] + edgecounts[i-1];
            }

            Array<OneD, int> vertids  (nTotVerts, 0);
            Array<OneD, int> edgeids  (nTotEdges, 0);
            Array<OneD, int> vertprocs(nTotVerts, 0);
            Array<OneD, int> edgeprocs(nTotEdges, 0);
            
            for (it  = extraDirVertIds.begin(), i = 0; 
                 it != extraDirVertIds.end(); ++it, ++i)
            {
                vertids  [vertoffsets[p]+i] = it->first;
                vertprocs[vertoffsets[p]+i] = it->second;
            }

            for (it  = extraDirEdgeIds.begin(), i = 0; 
                 it != extraDirEdgeIds.end(); ++it, ++i)
            {
                edgeids  [edgeoffsets[p]+i] = it->first;
                edgeprocs[edgeoffsets[p]+i] = it->second;
            }

            m_comm->AllReduce(vertids,   LibUtilities::ReduceSum);
            m_comm->AllReduce(vertprocs, LibUtilities::ReduceSum);
            m_comm->AllReduce(edgeids,   LibUtilities::ReduceSum);
            m_comm->AllReduce(edgeprocs, LibUtilities::ReduceSum);
            
            set<int> extraDirVerts;
            set<int> extraDirEdges;

            for (i = 0; i < nTotVerts; ++i)
            {
                if (m_comm->GetRank() != vertprocs[i])
                {
                    continue;
                }
                
                extraDirVerts.insert(vertids[i]);
            }
            
            for (i = 0; i < nTotEdges; ++i)
            {
                if (m_comm->GetRank() != edgeprocs[i])
                {
                    continue;
                }
                
                extraDirEdges.insert(edgeids[i]);
            }

            // Check between processes if the whole system is singular
            int s = (systemSingular ? 1 : 0);
            m_comm->AllReduce(s, LibUtilities::ReduceMin);
            systemSingular = (s == 1 ? true : false);

            // Count the number of boundary regions on each process
            Array<OneD, int> bccounts(n, 0);
            bccounts[p] = bndCondExp.num_elements();
            m_comm->AllReduce(bccounts, LibUtilities::ReduceSum);

            // Find the process rank with the maximum number of boundary regions
            int maxBCIdx = Vmath::Imax(n, bccounts, 1);

            // If the system is singular, the process with the maximum number of
            // BCs will set a Dirichlet vertex to make system non-singular.
            // Note: we find the process with maximum boundary regions to ensure
            // we do not try to set a Dirichlet vertex on a partition with no
            // intersection with the boundary.
            if(systemSingular == true && maxBCIdx == p)
            {
                if(m_session->DefinesParameter("SingularElement"))
                {
                    int s_eid;
                    m_session->LoadParameter("SingularElement", s_eid);

                    ASSERTL1(s_eid < locExpVector.size(),
                             "SingularElement Parameter is too large");
                    
                    meshVertId = locExpVector[s_eid]->GetGeom2D()->GetVid(0);
                }
                else
                {
                    //last region i and j=0 edge
                    bndCondFaceExp = boost::dynamic_pointer_cast<
                        StdRegions::StdExpansion2D>(
                            bndCondExp[bndCondExp.num_elements()-1]->GetExp(0));
                    
                    //first vertex 0 of the edge
                    meshVertId = bndCondFaceExp->GetGeom2D()->GetVid(0);
                }

                if(vertReorderedGraphVertId.count(meshVertId) == 0)
                {
                    vertReorderedGraphVertId[meshVertId] = graphVertId++;
                }
            }
            

            m_numLocalDirBndCoeffs = nLocDirBndCondDofs + nExtraDirichlet;
            firstNonDirGraphVertId = graphVertId;

            /**
             * STEP 2: Now order all other vertices and edges in the graph and
             * create a temporary numbering of domain-interior vertices and
             * edges.
             */
            typedef boost::adjacency_list<
                boost::setS, boost::vecS, boost::undirectedS> BoostGraph;
            typedef boost::graph_traits<
                BoostGraph>::vertex_descriptor BoostVertex;
            BoostGraph boostGraphObj;

            map<int, int>    vertTempGraphVertId;
            map<int, int>    edgeTempGraphVertId;
            map<int, int>    faceTempGraphVertId;
            map<int, int>    vwgts_map;
            Array<OneD, int> localVerts;
            Array<OneD, int> localEdges;
            Array<OneD, int> localFaces;

            int tempGraphVertId = 0;
            int localVertOffset = 0;
            int localEdgeOffset = 0;
            int localFaceOffset = 0;
            int nTotalVerts     = 0;
            int nTotalEdges     = 0;
            int nTotalFaces     = 0;
            int nVerts;
            int nEdges;
            int nFaces;
            int vertCnt;
            int edgeCnt;
            int faceCnt;

            m_numLocalBndCoeffs = 0;

            /// - Periodic vertices
            for(mapConstIt  = periodicVerticesId.begin(); 
                mapConstIt != periodicVerticesId.end(); mapConstIt++)
            {
                meshVertId  = mapConstIt->first;
                meshVertId2 = mapConstIt->second;

                if(vertReorderedGraphVertId.count(meshVertId) == 0)
                {

                    if(vertReorderedGraphVertId.count(meshVertId2) == 0)
                    {

                        if(vertTempGraphVertId.count(meshVertId) == 0)
                        {
                            vertTempGraphVertId[meshVertId]  = tempGraphVertId;
                            if(vertTempGraphVertId.count(meshVertId2) == 0)
                            {
                                vertTempGraphVertId[meshVertId2] = tempGraphVertId++;
                            }
                            else
                            {
                                ASSERTL0(false,"Unexplained Periodicity connectivity");
                            }
                        }
                        else
                        {
                            if(vertTempGraphVertId.count(meshVertId2) == 0)
                            {
                                ASSERTL0(false,"Unexplained Periodicity connectivity");
                            }
                            else // Doubly periodic region
                            {
                                int id1 = vertTempGraphVertId[meshVertId];
                                int id2 = vertTempGraphVertId[meshVertId2];
                                int id;

                                if(id1 != id2)
                                {
                                    // Reset any values set to
                                    // id2 to id1. In addition
                                    // if local id is greater
                                    // than id2 decrement list
                                    for(mapIt = vertTempGraphVertId.begin();
                                        mapIt != vertTempGraphVertId.end(); mapIt++)
                                    {
                                        id = mapIt->second;
                                        if(id == id2)
                                        {
                                            vertTempGraphVertId[mapIt->first] = id1;
                                        }
                                        else if (id > id2)
                                        {
                                            vertTempGraphVertId[mapIt->first] = id-1;
                                        }
                                    }
                                    tempGraphVertId--;
                                }
                            }
                        }

                    }
                    else
                    {
                        vertReorderedGraphVertId[meshVertId] = vertReorderedGraphVertId[meshVertId2];
                    }
                }
                else
                {
                    if(vertReorderedGraphVertId.count(meshVertId2) == 0)
                    {
                        vertReorderedGraphVertId[meshVertId2] = vertReorderedGraphVertId[meshVertId];
                    }
                }
            }

            /// - Periodic edges
            for(mapConstIt  = periodicEdgesId.begin(); 
                mapConstIt != periodicEdgesId.end(); 
                mapConstIt++)
            {
                meshEdgeId  = mapConstIt->first;
                meshEdgeId2 = mapConstIt->second;

                if(edgeReorderedGraphVertId.count(meshEdgeId) == 0)
                {
                    if(edgeReorderedGraphVertId.count(meshEdgeId2) == 0)
                    {
                        if(edgeTempGraphVertId.count(meshEdgeId) == 0)
                        {
                            edgeTempGraphVertId[meshEdgeId]  = tempGraphVertId;
                            if(edgeTempGraphVertId.count(meshEdgeId2) == 0)
                            {
                                edgeTempGraphVertId[meshEdgeId2] = tempGraphVertId++;
                            }
                            else
                            {
                                ASSERTL0(false,"Unexplained Periodicity connectivity");
                            }
                        }
                        else
                        {
                            if(edgeTempGraphVertId.count(meshEdgeId2) == 0)
                            {
                                ASSERTL0(false,"Unexplained Periodicity connectivity");
                            }
                            else // Doubly periodic region
                            {
                                int id1 = edgeTempGraphVertId[meshEdgeId];
                                int id2 = edgeTempGraphVertId[meshEdgeId2];
                                int id;

                                if(id1 != id2)
                                {
                                    // Reset any values set to
                                    // id2 to id1. In addition
                                    // if local id is greater
                                    // than id2 decrement list
                                    for(mapIt = edgeTempGraphVertId.begin();
                                        mapIt != edgeTempGraphVertId.end(); mapIt++)
                                    {
                                        id = mapIt->second;
                                        if(id == id2)
                                        {
                                            edgeTempGraphVertId[mapIt->first] = id1;
                                        }
                                        else if (id > id2)
                                        {
                                            edgeTempGraphVertId[mapIt->first] = id-1;
                                        }
                                    }
                                    tempGraphVertId--;
                                }
                            }
                        }
                    }
                    else
                    {
                        edgeReorderedGraphVertId[meshEdgeId] = edgeReorderedGraphVertId[meshEdgeId2];
                    }
                }
                else
                {
                    if(edgeReorderedGraphVertId.count(meshEdgeId2) == 0)
                    {
                        edgeReorderedGraphVertId[meshEdgeId2] = edgeReorderedGraphVertId[meshEdgeId];
                    }
                    /*
                    else
                    {
                        ASSERTL0(edgeReorderedGraphVertId[meshEdgeId2] == edgeReorderedGraphVertId[meshEdgeId],
                                 "These values should be equal");
                    }
                    */
                }
            }

            /// - Periodic faces
            for(mapFaceIt  = periodicFacesId.begin(); 
                mapFaceIt != periodicFacesId.end(); 
                mapFaceIt++)
            {
                meshFaceId  = mapFaceIt->first;
                meshFaceId2 = mapFaceIt->second.first;

                if(meshFaceId < meshFaceId2)
                {
                    ASSERTL0(faceReorderedGraphVertId.count(meshFaceId) == 0,
                             "This periodic boundary face has been specified before");
                    ASSERTL0(faceReorderedGraphVertId.count(meshFaceId2) == 0,
                             "This periodic boundary face has been specified before");

                    faceTempGraphVertId[meshFaceId]  = tempGraphVertId;
                    faceTempGraphVertId[meshFaceId2] = tempGraphVertId++;
                }
            }


            /// - All other vertices and edges
            for(i = 0; i < locExpVector.size(); ++i)
            {
                if((locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion3D>(
                        locExpVector[locExp.GetOffset_Elmt_Id(i)])))
                {
                    nTotalVerts += locExpansion->GetNverts();
                    nTotalEdges += locExpansion->GetNedges();
                    nTotalFaces += locExpansion->GetNfaces();
                }
            }

            // Store the temporary graph vertex
            // id's of all element edges and
            // vertices in these 3 arrays below
            localVerts = Array<OneD, int>(nTotalVerts,-1);
            localEdges = Array<OneD, int>(nTotalEdges,-1);
            localFaces = Array<OneD, int>(nTotalFaces,-1);

            for(i = 0; i < locExpVector.size(); ++i)
            {
                if((locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion3D>(
                        locExpVector[locExp.GetOffset_Elmt_Id(i)])))
                {
                    vertCnt = 0;
                    nVerts = locExpansion->GetNverts();
                    for(j = 0; j < nVerts; ++j)
                    {
                        meshVertId = (locExpansion->GetGeom3D())->GetVid(j);
                        if(vertReorderedGraphVertId.count(meshVertId) == 0)
                        {
                            if(vertTempGraphVertId.count(meshVertId) == 0)
                            {
                                boost::add_vertex(boostGraphObj);
                                vertTempGraphVertId[meshVertId] = tempGraphVertId++;
                                m_numNonDirVertexModes+=1;
                            }
                            localVerts[localVertOffset+vertCnt++] = vertTempGraphVertId[meshVertId];
                            vwgts_map[ vertTempGraphVertId[meshVertId] ] = 1;
                        }
                    }
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a local 3D expansion failed");
                }
                localVertOffset+=nVerts;
            }



            for(i = 0; i < locExpVector.size(); ++i)
            {
                if((locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion3D>(
                        locExpVector[locExp.GetOffset_Elmt_Id(i)])))
                {
                    edgeCnt = 0;
                    nEdges = locExpansion->GetNedges();

                    for(j = 0; j < nEdges; ++j)
                    {
                        nEdgeInteriorCoeffs = locExpansion->GetEdgeNcoeffs(j) - 2;
                        meshEdgeId = (locExpansion->GetGeom3D())->GetEid(j);
                        if(edgeReorderedGraphVertId.count(meshEdgeId) == 0)
                        {
                            if(edgeTempGraphVertId.count(meshEdgeId) == 0)
                            {
                                boost::add_vertex(boostGraphObj);
                                edgeTempGraphVertId[meshEdgeId] = tempGraphVertId++;
                                m_numNonDirEdgeModes+=nEdgeInteriorCoeffs;
                            }
                            localEdges[localEdgeOffset+edgeCnt++] = edgeTempGraphVertId[meshEdgeId];
                            vwgts_map[ edgeTempGraphVertId[meshEdgeId] ] = nEdgeInteriorCoeffs;
                        }
                    }
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a local 3D expansion failed");
                }
                localEdgeOffset+=nEdges;
            }

            for(i = 0; i < locExpVector.size(); ++i)
            {
                if((locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion3D>(
                        locExpVector[locExp.GetOffset_Elmt_Id(i)])))
                {
                    nFaces = locExpansion->GetNfaces();
                    faceCnt = 0;
                    for(j = 0; j < nFaces; ++j)
                    {
                        nFaceInteriorCoeffs = locExpansion->GetFaceIntNcoeffs(j);
                        meshFaceId = (locExpansion->GetGeom3D())->GetFid(j);
                        if(faceReorderedGraphVertId.count(meshFaceId) == 0)
                        {
                            if(faceTempGraphVertId.count(meshFaceId) == 0)
                            {
                                boost::add_vertex(boostGraphObj);
                                faceTempGraphVertId[meshFaceId] = tempGraphVertId++;
                                m_numNonDirFaceModes+=nFaceInteriorCoeffs;
                            }
                            localFaces[localFaceOffset+faceCnt++] = faceTempGraphVertId[meshFaceId];
                            vwgts_map[ faceTempGraphVertId[meshFaceId] ] = nFaceInteriorCoeffs;
                        }
                    }
                    m_numLocalBndCoeffs += locExpansion->NumBndryCoeffs();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a local 3D expansion failed");
                }
                localFaceOffset+=nFaces;
            }

            localVertOffset=0;
            localEdgeOffset=0;
            localFaceOffset=0;
            for(i = 0; i < locExpVector.size(); ++i)
            {
                if((locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion3D>(
                        locExpVector[locExp.GetOffset_Elmt_Id(i)])))
                {
                    nVerts = locExpansion->GetNverts();
                    nEdges = locExpansion->GetNedges();
                    nFaces = locExpansion->GetNfaces();
                    // Now loop over all local faces, edges and vertices
                    // of this element and define that all other
                    // faces, edges and verices of this element are
                    // adjacent to them.

                    // Vertices
                    for(j = 0; j < nVerts; j++)
                    {
                        if(localVerts[j+localVertOffset]==-1)
                        {
                            break;
                        }
                        // associate to other vertices
                        for(k = 0; k < nVerts; k++)
                        {
                            if(localVerts[k+localVertOffset]==-1)
                            {
                                break;
                            }
                            if(k!=j)
                            {
                                boost::add_edge( (size_t) localVerts[j+localVertOffset], 
                                                 (size_t) localVerts[k+localVertOffset],boostGraphObj);
                            }
                        }
                        // associate to other edges
                        for(k = 0; k < nEdges; k++)
                        {
                            if(localEdges[k+localEdgeOffset]==-1)
                            {
                                break;
                            }
                            boost::add_edge( (size_t) localVerts[j+localVertOffset], 
                                             (size_t) localEdges[k+localEdgeOffset],boostGraphObj);
                        }
                        // associate to other faces
                        for(k = 0; k < nFaces; k++)
                        {
                            if(localFaces[k+localFaceOffset]==-1)
                            {
                                break;
                            }
                            boost::add_edge( (size_t) localVerts[j+localVertOffset], 
                                             (size_t) localFaces[k+localFaceOffset],boostGraphObj);
                        }
                    }

                    // Edges
                    for(j = 0; j < nEdges; j++)
                    {
                        if(localEdges[j+localEdgeOffset]==-1)
                        {
                            break;
                        }
                        // Associate to other edges
                        for(k = 0; k < nEdges; k++)
                        {
                            if(localEdges[k+localEdgeOffset]==-1)
                            {
                                break;
                            }
                            if(k!=j)
                            {
                                boost::add_edge( (size_t) localEdges[j+localEdgeOffset], 
                                                 (size_t) localEdges[k+localEdgeOffset],boostGraphObj);
                            }
                        }
                        // Associate to vertices
                        for(k = 0; k < nVerts; k++)
                        {
                            if(localVerts[k+localVertOffset]==-1)
                            {
                                break;
                            }
                            boost::add_edge( (size_t) localEdges[j+localEdgeOffset], 
                                             (size_t) localVerts[k+localVertOffset],boostGraphObj);
                        }
                        // Associate to faces
                        for(k = 0; k < nFaces; k++)
                        {
                            if(localFaces[k+localFaceOffset]==-1)
                            {
                                break;
                            }
                            boost::add_edge( (size_t) localEdges[j+localEdgeOffset], 
                                             (size_t) localFaces[k+localFaceOffset],boostGraphObj);
                        }
                    }

                    // Faces
                    for(j = 0; j < nFaces; j++)
                    {
                        if(localFaces[j+localFaceOffset]==-1)
                        {
                            break;
                        }
                        // Associate to other faces
                        for(k = 0; k < nFaces; k++)
                        {
                            if(localFaces[k+localFaceOffset]==-1)
                            {
                                break;
                            }
                            if(k!=j)
                            {
                                boost::add_edge( (size_t) localFaces[j+localFaceOffset], 
                                                 (size_t) localFaces[k+localFaceOffset],boostGraphObj);
                            }
                        }
                        // Associate to vertices
                        for(k = 0; k < nVerts; k++)
                        {
                            if(localVerts[k+localVertOffset]==-1)
                            {
                                break;
                            }
                            boost::add_edge( (size_t) localFaces[j+localFaceOffset], 
                                             (size_t) localVerts[k+localVertOffset],boostGraphObj);
                        }
                        // Associate to edges
                        for(k = 0; k < nEdges; k++)
                        {
                            if(localEdges[k+localEdgeOffset]==-1)
                            {
                                break;
                            }
                            boost::add_edge( (size_t) localFaces[j+localFaceOffset], 
                                             (size_t) localEdges[k+localEdgeOffset],boostGraphObj);
                        }
                    }
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a local 3D expansion failed");
                }
                localVertOffset+=nVerts;
                localEdgeOffset+=nEdges;
                localFaceOffset+=nFaces;
            }

            // Container to store vertices of the graph which correspond to
            // degrees of freedom along the boundary.
            set<int> partVerts;
            
            if (m_solnType == eIterativeMultiLevelStaticCond)
            {
                vector<long> procVerts,  procEdges,  procFaces;
                set   <int>  foundVerts, foundEdges, foundFaces;
                
                // Loop over element and construct the procVerts and procEdges
                // vectors, which store the geometry IDs of mesh vertices and
                // edges respectively which are local to this process.
                for(i = cnt = 0; i < locExpVector.size(); ++i)
                {
                    int elmtid = locExp.GetOffset_Elmt_Id(i);
                    if((locExpansion = boost::dynamic_pointer_cast<
                            StdRegions::StdExpansion3D>(locExpVector[elmtid])))
                    {
                        for (j = 0; j < locExpansion->GetNverts(); ++j)
                        {
                            int vid = locExpansion->GetGeom3D()->GetVid(j)+1;
                            
                            if (foundVerts.count(vid) == 0)
                            {
                                procVerts.push_back(vid);
                                foundVerts.insert(vid);
                            }
                        }

                        for (j = 0; j < locExpansion->GetNedges(); ++j)
                        {
                            int eid = locExpansion->GetGeom3D()->GetEid(j)+1;

                            if (foundEdges.count(eid) == 0)
                            {
                                procEdges.push_back(eid);
                                foundEdges.insert(eid);
                            }
                        }

                        for (j = 0; j < locExpansion->GetNfaces(); ++j)
                        {
                            int fid = locExpansion->GetGeom3D()->GetFid(j)+1;
                            
                            if (foundFaces.count(fid) == 0)
                            {
                                procFaces.push_back(fid);
                                foundFaces.insert(fid);
                            }
                        }
                    }
                    else
                    {
                        ASSERTL0(false,
                                 "dynamic cast to a local 3D expansion failed");
                    }
                }

                int unique_verts = foundVerts.size();
                int unique_edges = foundEdges.size();
                int unique_faces = foundFaces.size();

                // Now construct temporary GS objects. These will be used to
                // populate the arrays tmp3 and tmp4 with the multiplicity of
                // the vertices and edges respectively to identify those
                // vertices and edges which are located on partition boundary.
                Array<OneD, long> vertArray(unique_verts, &procVerts[0]);
                Array<OneD, long> edgeArray(unique_edges, &procEdges[0]);
                Array<OneD, long> faceArray(unique_faces, &procFaces[0]);
                Gs::gs_data *tmp1 = Gs::Init(vertArray, m_comm);
                Gs::gs_data *tmp2 = Gs::Init(edgeArray, m_comm);
                Gs::gs_data *tmp3 = Gs::Init(faceArray, m_comm);
                Array<OneD, NekDouble> tmp4(unique_verts, 1.0);
                Array<OneD, NekDouble> tmp5(unique_edges, 1.0);
                Array<OneD, NekDouble> tmp6(unique_faces, 1.0);
                Gs::Gather(tmp4, Gs::gs_add, tmp1);
                Gs::Gather(tmp5, Gs::gs_add, tmp2);
                Gs::Gather(tmp6, Gs::gs_add, tmp3);

                // Finally, fill the partVerts set with all non-Dirichlet
                // vertices which lie on a partition boundary.
                for (i = 0; i < unique_verts; ++i)
                {
                    if (tmp4[i] > 1.0)
                    {
                        if (vertReorderedGraphVertId.count(procVerts[i]-1) == 0)
                        {
                            partVerts.insert(vertTempGraphVertId[procVerts[i]-1]);
                        }
                    }
                }
            
                for (i = 0; i < unique_edges; ++i)
                {
                    if (tmp5[i] > 1.0)
                    {
                        if (edgeReorderedGraphVertId.count(procEdges[i]-1) == 0)
                        {
                            partVerts.insert(edgeTempGraphVertId[procEdges[i]-1]);
                        }
                    }
                }

                for (i = 0; i < unique_faces; ++i)
                {
                    if (tmp6[i] > 1.0)
                    {
                        if (faceReorderedGraphVertId.count(procFaces[i]-1) == 0)
                        {
                            partVerts.insert(faceTempGraphVertId[procFaces[i]-1]);
                        }
                    }
                }
            }
            
            /**
             * STEP 3: Reorder graph for optimisation.
             */
            BottomUpSubStructuredGraphSharedPtr bottomUpGraph;
            int nGraphVerts = tempGraphVertId;
            Array<OneD, int> perm(nGraphVerts);
            Array<OneD, int> iperm(nGraphVerts);
            Array<OneD, int> vwgts(nGraphVerts);
            ASSERTL1(vwgts_map.size()==nGraphVerts,"Non matching dimensions");
            for(i = 0; i < nGraphVerts; ++i)
            {
                vwgts[i] = vwgts_map[i];
            }

            if(nGraphVerts)
            {
                switch(m_solnType)
                {
                case eDirectFullMatrix:
                case eIterativeFull:
                case eIterativeStaticCond:
                case eXxtFullMatrix:
                case eXxtStaticCond:
                    {
                        NoReordering(boostGraphObj,perm,iperm);
                    }
                    break;
                case eDirectStaticCond:
                    {
                        CuthillMckeeReordering(boostGraphObj,perm,iperm);
                    }
                    break;
                case eDirectMultiLevelStaticCond:
                case eIterativeMultiLevelStaticCond:
                case eXxtMultiLevelStaticCond:
                    {
                        MultiLevelBisectionReordering(boostGraphObj,perm,iperm,bottomUpGraph,partVerts);
                    }
                    break;
                default:
                    {
                        ASSERTL0(false,"Unrecognised solution type: " + std::string(MultiRegions::GlobalSysSolnTypeMap[m_solnType]));
                    }
                }
            }

            // For parallel multi-level static condensation determine the lowest
            // static condensation level amongst processors.
            if (m_solnType == eIterativeMultiLevelStaticCond)
            {
                m_lowestStaticCondLevel = bottomUpGraph->GetNlevels()-1;
                m_comm->AllReduce(m_lowestStaticCondLevel, 
                                  LibUtilities::ReduceMax);
            }
            else
            {
                m_lowestStaticCondLevel = 0;
            }

            /**
             * STEP 4: Fill the #vertReorderedGraphVertId and
             * #edgeReorderedGraphVertId with the optimal ordering from boost.
             */
            for(mapIt = vertTempGraphVertId.begin(); mapIt != vertTempGraphVertId.end(); mapIt++)
            {
                vertReorderedGraphVertId[mapIt->first] = iperm[mapIt->second] + graphVertId;
            }
            for(mapIt = edgeTempGraphVertId.begin(); mapIt != edgeTempGraphVertId.end(); mapIt++)
            {
                edgeReorderedGraphVertId[mapIt->first] = iperm[mapIt->second] + graphVertId;
            }
            for(mapIt = faceTempGraphVertId.begin(); mapIt != faceTempGraphVertId.end(); mapIt++)
            {
                faceReorderedGraphVertId[mapIt->first] = iperm[mapIt->second] + graphVertId;
            }


            /**
             * STEP 5: Set up an array which contains the offset information of
             * the different graph vertices.
             *
             * This basically means to identify to how many global degrees of
             * freedom the individual graph vertices correspond. Obviously,
             * the graph vertices corresponding to the mesh-vertices account
             * for a single global DOF. However, the graph vertices
             * corresponding to the element edges correspond to N-2 global DOF
             * where N is equal to the number of boundary modes on this edge.
             */
            Array<OneD, int> graphVertOffset(vertReorderedGraphVertId.size()+
                                             edgeReorderedGraphVertId.size()+
                                             faceReorderedGraphVertId.size()+1);
            graphVertOffset[0] = 0;

            for(i = 0; i < locExpVector.size(); ++i)
            {
                locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion3D>(locExpVector[locExp.GetOffset_Elmt_Id(i)]);

                for(j = 0; j < locExpansion->GetNverts(); ++j)
                {
                    meshVertId = (locExpansion->GetGeom3D())->GetVid(j);
                    graphVertOffset[vertReorderedGraphVertId[meshVertId]+1] = 1;
                }

                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    nEdgeInteriorCoeffs = locExpansion->GetEdgeNcoeffs(j) - 2;
                    meshEdgeId = (locExpansion->GetGeom3D())->GetEid(j);
                    graphVertOffset[edgeReorderedGraphVertId[meshEdgeId]+1] = nEdgeInteriorCoeffs;

                    bType = locExpansion->GetEdgeBasisType(j);
                    // need a sign vector for modal expansions if nEdgeCoeffs >=4
                    if( (nEdgeInteriorCoeffs+2 >= 4)&&
                        ( (bType == LibUtilities::eModified_A)||
                          (bType == LibUtilities::eModified_B) ) )
                    {
                        m_signChange = true;
                    }
                }

                for(j = 0; j < locExpansion->GetNfaces(); ++j)
                {
                    nFaceInteriorCoeffs = locExpansion->GetFaceIntNcoeffs(j);
                    meshFaceId = (locExpansion->GetGeom3D())->GetFid(j);
                    graphVertOffset[faceReorderedGraphVertId[meshFaceId]+1] = nFaceInteriorCoeffs;
                }
            }

            for(i = 1; i < graphVertOffset.num_elements(); i++)
            {
                graphVertOffset[i] += graphVertOffset[i-1];
            }
            
            // Allocate the proper amount of space for the class-data
            m_numLocalCoeffs                 = numLocalCoeffs;
            m_numGlobalDirBndCoeffs          = graphVertOffset[firstNonDirGraphVertId];
            m_localToGlobalMap               = Array<OneD, int>(m_numLocalCoeffs,-1);
            m_localToGlobalBndMap            = Array<OneD, int>(m_numLocalBndCoeffs,-1);
            m_bndCondCoeffsToGlobalCoeffsMap = Array<OneD,int>(nLocBndCondDofs,-1);
            // If required, set up the sign-vector
            if(m_signChange)
            {
                m_localToGlobalSign = Array<OneD, NekDouble>(m_numLocalCoeffs,1.0);
                m_localToGlobalBndSign = Array<OneD, NekDouble>(m_numLocalBndCoeffs,1.0);
                m_bndCondCoeffsToGlobalCoeffsSign = Array<OneD,NekDouble>(nLocBndCondDofs,1.0);
            }

            m_staticCondLevel = 0;
            m_numPatches =  locExpVector.size();
            m_numLocalBndCoeffsPerPatch = Array<OneD, unsigned int>(m_numPatches);
            m_numLocalIntCoeffsPerPatch = Array<OneD, unsigned int>(m_numPatches);
            for(i = 0; i < m_numPatches; ++i)
            {
                m_numLocalBndCoeffsPerPatch[i] = (unsigned int) 
                    locExpVector[locExp.GetOffset_Elmt_Id(i)]->NumBndryCoeffs();
                m_numLocalIntCoeffsPerPatch[i] = (unsigned int) 
                    locExpVector[locExp.GetOffset_Elmt_Id(i)]->GetNcoeffs() - 
                    locExpVector[locExp.GetOffset_Elmt_Id(i)]->NumBndryCoeffs();
            }

            /**
             * STEP 6: Now, all ingredients are ready to set up the actual
             * local to global mapping.
             *
             * The remainder of the map consists of the element-interior
             * degrees of freedom. This leads to the block-diagonal submatrix
             * as each element-interior mode is globally orthogonal to modes
             * in all other elements.
             */
            cnt = 0;
            // Loop over all the elements in the domain
            for(i = 0; i < locExpVector.size(); ++i)
            {
                locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion3D>(locExpVector[i]);
                cnt = locExp.GetCoeff_Offset(i);
                for(j = 0; j < locExpansion->GetNverts(); ++j)
                {
                    meshVertId          = (locExpansion->GetGeom3D())->GetVid(j);

                    // Set the global DOF for vertex j of element i
                    m_localToGlobalMap[cnt+locExpansion->GetVertexMap(j)] =
                        graphVertOffset[vertReorderedGraphVertId[meshVertId]];
                }

                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    nEdgeInteriorCoeffs = locExpansion->GetEdgeNcoeffs(j)-2;
                    edgeOrient          = (locExpansion->GetGeom3D())->GetEorient(j);
                    meshEdgeId          = (locExpansion->GetGeom3D())->GetEid(j);

                    locExpansion->GetEdgeInteriorMap(j,edgeOrient,edgeInteriorMap,edgeInteriorSign);

                    // Set the global DOF's for the interior modes of edge j
                    for(k = 0; k < nEdgeInteriorCoeffs; ++k)
                    {
                        m_localToGlobalMap[cnt+edgeInteriorMap[k]] =
                            graphVertOffset[edgeReorderedGraphVertId[meshEdgeId]]+k;
                    }

                    // Fill the sign vector if required
                    if(m_signChange)
                    {
                        for(k = 0; k < nEdgeInteriorCoeffs; ++k)
                        {
                            m_localToGlobalSign[cnt+edgeInteriorMap[k]] = (NekDouble) edgeInteriorSign[k];
                        }
                    }
                }

                for(j = 0; j < locExpansion->GetNfaces(); ++j)
                {
                    map<int, pair<int, StdRegions::Orientation> >::const_iterator it;
                    
                    nFaceInteriorCoeffs = locExpansion->GetFaceIntNcoeffs(j);
                    faceOrient          = (locExpansion->GetGeom3D())->GetFaceOrient(j);
                    meshFaceId          = (locExpansion->GetGeom3D())->GetFid(j);
                    
                    /*
                    it = periodicFaces.find(meshFaceId);
                    if (it == periodicFaces.begin())
                    {
                        
                    }
                    */
                    
                    locExpansion->GetFaceInteriorMap(j,faceOrient,faceInteriorMap,faceInteriorSign);

                    // Set the global DOF's for the interior modes of face j
                    for(k = 0; k < nFaceInteriorCoeffs; ++k)
                    {
                        m_localToGlobalMap[cnt+faceInteriorMap[k]] =
                            graphVertOffset[faceReorderedGraphVertId[meshFaceId]]+k;
                    }

                    if(m_signChange)
                    {
                        for(k = 0; k < nFaceInteriorCoeffs; ++k)
                        {
                            m_localToGlobalSign[cnt+faceInteriorMap[k]] = (NekDouble) faceInteriorSign[k];
                        }
                    }
                }
            }

            // Set up the mapping for the boundary conditions
            cnt = 0;
            int offset = 0;
            for(i = 0; i < bndCondExp.num_elements(); i++)
            {
                set<int> foundExtraVerts, foundExtraEdges;
                for(j = 0; j < bndCondExp[i]->GetNumElmts(); j++)
                {
                    bndCondFaceExp  = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(bndCondExp[i]->GetExp(j));
                    cnt = offset + bndCondExp[i]->GetCoeff_Offset(j);
                    for(k = 0; k < bndCondFaceExp->GetNverts(); k++)
                    {
                        meshVertId = (bndCondFaceExp->GetGeom2D())->GetVid(k);
                        m_bndCondCoeffsToGlobalCoeffsMap[cnt+bndCondFaceExp->GetVertexMap(k)] = graphVertOffset[vertReorderedGraphVertId[meshVertId]];

                        if (bndConditions[i]->GetBoundaryConditionType() != 
                                SpatialDomains::eDirichlet)
                        {
                            continue;
                        }
                        
                        set<int>::iterator iter = extraDirVerts.find(meshVertId);
                        if (iter != extraDirVerts.end() && 
                            foundExtraVerts.count(meshVertId) == 0)
                        {
                            int loc = bndCondExp[i]->GetCoeff_Offset(j) + 
                                bndCondFaceExp->GetVertexMap(k);
                            int gid = graphVertOffset[
                                vertReorderedGraphVertId[meshVertId]];
                            m_extraDirDofs[i].push_back(make_pair(loc,gid));
                            foundExtraVerts.insert(meshVertId);
                        }
                    }

                    for(k = 0; k < bndCondFaceExp->GetNedges(); k++)
                    {
                        nEdgeInteriorCoeffs = bndCondFaceExp->GetEdgeNcoeffs(k)-2;
                        edgeOrient          = (bndCondFaceExp->GetGeom2D())->GetEorient(k);
                        meshEdgeId          = (bndCondFaceExp->GetGeom2D())->GetEid(k);

                        bndCondFaceExp->GetEdgeInteriorMap(
                            k,edgeOrient,edgeInteriorMap,edgeInteriorSign);

                        for(l = 0; l < nEdgeInteriorCoeffs; ++l)
                        {
                            m_bndCondCoeffsToGlobalCoeffsMap[cnt+edgeInteriorMap[l]] =
                                graphVertOffset[edgeReorderedGraphVertId[meshEdgeId]]+l;
                        }

                        // Fill the sign vector if required
                        if(m_signChange)
                        {
                            for(l = 0; l < nEdgeInteriorCoeffs; ++l)
                            {
                                m_bndCondCoeffsToGlobalCoeffsSign[cnt+edgeInteriorMap[l]] = (NekDouble) edgeInteriorSign[l];
                            }
                        }
                        
                        if (bndConditions[i]->GetBoundaryConditionType() != 
                                SpatialDomains::eDirichlet)
                        {
                            continue;
                        }
                        
                        set<int>::iterator iter = extraDirEdges.find(meshEdgeId);
                        if (iter != extraDirEdges.end()            && 
                            foundExtraEdges.count(meshEdgeId) == 0 &&
                            nEdgeInteriorCoeffs > 0)
                        {
                            for(l = 0; l < nEdgeInteriorCoeffs; ++l)
                            {
                                int loc = bndCondExp[i]->GetCoeff_Offset(j) + 
                                    edgeInteriorMap[l];
                                int gid = graphVertOffset[
                                    edgeReorderedGraphVertId[meshEdgeId]]+l;
                                m_extraDirDofs[i].push_back(make_pair(loc,gid));
                            }
                            foundExtraEdges.insert(meshEdgeId);
                        }
                    }

                    meshFaceId = (bndCondFaceExp->GetGeom2D())->GetFid();
                    intDofCnt = 0;
                    for(k = 0; k < bndCondFaceExp->GetNcoeffs(); k++)
                    {
                        if(m_bndCondCoeffsToGlobalCoeffsMap[cnt+k] == -1)
                        {
                            m_bndCondCoeffsToGlobalCoeffsMap[cnt+k] =
                                graphVertOffset[faceReorderedGraphVertId[meshFaceId]]+intDofCnt;
                            intDofCnt++;
                        }
                    }
                }
                offset += bndCondExp[i]->GetNcoeffs();
            }

            globalId = Vmath::Vmax(m_numLocalCoeffs,&m_localToGlobalMap[0],1)+1;
            m_numGlobalBndCoeffs = globalId;


            /**
             * STEP 7: The boundary condition mapping is generated from the
             * same vertex renumbering.
             */
            cnt=0;
            for(i = 0; i < m_numLocalCoeffs; ++i)
            {
                if(m_localToGlobalMap[i] == -1)
                {
                    m_localToGlobalMap[i] = globalId++;
                }
                else
                {
                    if(m_signChange)
                    {
                        m_localToGlobalBndSign[cnt]=m_localToGlobalSign[i];
                    }
                    m_localToGlobalBndMap[cnt++]=m_localToGlobalMap[i];
                }
            }

            m_numGlobalCoeffs = globalId;

            SetUpUniversalC0ContMap(locExp);

            // Set up the local to global map for the next level when using
            // multi-level static condensation
            if ((m_solnType == eDirectMultiLevelStaticCond ||
                 m_solnType == eIterativeMultiLevelStaticCond) && nGraphVerts)
            {
                if (m_staticCondLevel < (bottomUpGraph->GetNlevels()-1))
                {
                    Array<OneD, int> vwgts_perm(nGraphVerts);
                    for(i = 0; i < nGraphVerts; ++i)
                    {
                        vwgts_perm[i] = vwgts[perm[i]];
                    }

                    bottomUpGraph->ExpandGraphWithVertexWeights(vwgts_perm);
                    m_nextLevelLocalToGlobalMap = MemoryManager<
                        AssemblyMap>::AllocateSharedPtr(this,bottomUpGraph);
                }
            }

            m_hash = boost::hash_range(m_localToGlobalMap.begin(), 
                                       m_localToGlobalMap.end());

            // Add up hash values if parallel
            int hash = m_hash;
            m_comm->GetRowComm()->AllReduce(hash, 
                              LibUtilities::ReduceSum);
            m_hash = hash;
        }
    } // namespace
} // namespace
