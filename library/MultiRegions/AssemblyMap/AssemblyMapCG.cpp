///////////////////////////////////////////////////////////////////////////////
//
// File AssemblyMapCG.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: C0-continuous Local to Global mapping routines, base class
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>
#include <MultiRegions/ExpList.h>
#include <LocalRegions/Expansion.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/BasicUtils/HashUtils.hpp>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class AssemblyMapCG
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
        AssemblyMapCG::AssemblyMapCG(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const std::string variable):
            AssemblyMap(pSession,variable)
        {
            pSession->LoadParameter(
                "MaxStaticCondLevel", m_maxStaticCondLevel, 100);
        }

        int AssemblyMapCG::CreateGraph(
            const ExpList                       &locExp,
            const BndCondExp                    &bndCondExp,
            const Array<OneD, const BndCond>    &bndConditions,
            const bool                           checkIfSystemSingular,
            const PeriodicMap                   &periodicVerts,
            const PeriodicMap                   &periodicEdges,
            const PeriodicMap                   &periodicFaces,
            DofGraph                            &graph,
            BottomUpSubStructuredGraphSharedPtr &bottomUpGraph,
            set<int>                            &extraDirVerts,
            set<int>                            &extraDirEdges,
            int                                 &firstNonDirGraphVertId,
            int                                 &nExtraDirichlet,
            int                                  mdswitch)
        {
            int graphVertId = 0;
            int vMaxVertId = -1;
            int i, j, k, l, cnt;
            int meshVertId, meshEdgeId, meshFaceId;
            int meshVertId2, meshEdgeId2;

            LocalRegions::ExpansionSharedPtr exp, bndExp;
            const LocalRegions::ExpansionVector &locExpVector =
                *(locExp.GetExp());
            LibUtilities::CommSharedPtr vComm = m_comm->GetRowComm();

            m_numLocalBndCondCoeffs = 0;
            m_systemSingular = checkIfSystemSingular;

            for(i = 0; i < bndCondExp.size(); i++)
            {

                m_numLocalBndCondCoeffs += bndCondExp[i]->GetNcoeffs();

                if (bndConditions[0][i]->GetBoundaryConditionType() ==
                       SpatialDomains::ePeriodic)
                {
                    continue;
                }

                // Check to see if any value on boundary has Dirichlet
                // value.  note this is a vector to manage coupled
                // solver but for scalar will just be a vector of size 11
                cnt = 0;
                for(k = 0; k < bndConditions.size(); ++k)
                {
                    if (bndConditions[k][i]->GetBoundaryConditionType() ==
                            SpatialDomains::eDirichlet)
                    {
                        cnt++;
                    }
                    if (bndConditions[k][i]->GetBoundaryConditionType() !=
                            SpatialDomains::eNeumann)
                    {
                        m_systemSingular = false;
                    }
                }

                // Find the maximum boundary vertex ID on this process. This is
                // used later to pin a vertex if the system is singular.
                for (j = 0; j < bndCondExp[i]->GetNumElmts(); ++j)
                {
                    bndExp = bndCondExp[i]->GetExp(j)->as<LocalRegions::Expansion>();
                    for (k = 0; k < bndExp->GetNverts(); ++k)
                    {
                        if (vMaxVertId < bndExp->GetGeom()->GetVid(k))
                        {
                            vMaxVertId = bndExp->GetGeom()->GetVid(k);
                        }
                    }
                }

                // If all boundaries are Dirichlet fill in graph
                if(cnt == bndConditions.size())
                {
                    for(j = 0; j < bndCondExp[i]->GetNumElmts(); j++)
                    {
                        bndExp = bndCondExp[i]->GetExp(j);

                        for (k = 0; k < bndExp->GetNverts(); k++)
                        {
                            meshVertId = bndExp->GetGeom()->GetVid(k);
                            if (graph[0].count(meshVertId) == 0)
                            {
                                graph[0][meshVertId] = graphVertId++;
                            }
                        }

                        for (k = 0; k < bndExp->GetNedges(); k++)
                        {
                            meshEdgeId = bndExp->GetGeom()->GetEid(k);
                            if (graph[1].count(meshEdgeId) == 0)
                            {
                                graph[1][meshEdgeId] = graphVertId++;
                            }
                        }

                        // Possibility of a face in 3D or edge in 2D
                        meshFaceId = bndExp->GetGeom()->GetGlobalID();
                        const int bndDim = bndExp->GetNumBases();
                        if (graph[bndDim].count(meshFaceId) == 0)
                        {
                            graph[bndDim][meshFaceId] = graphVertId++;
                        }
                        m_numLocalDirBndCoeffs += bndExp->GetNcoeffs();
                    }
                }
            }

            // Number of dirichlet edges and faces (not considering periodic
            // BCs)
            m_numDirEdges = graph[1].size();
            m_numDirFaces = graph[2].size();

            /*
             * The purpose of this routine is to deal with those degrees of
             * freedom that are Dirichlet, but do not have a local Dirichlet
             * boundary condition expansion set.
             *
             * For example, in 2D, consider a triangulation of a square into two
             * triangles. Now imagine one edge of the square is Dirichlet and
             * the problem is run on two processors. On one processor, one
             * triangle vertex is Dirichlet, but doesn't know this since the
             * Dirichlet composite lives on the other processor.
             *
             * When the global linear system is solved therefore, there is an
             * inconsistency that at best leads to an inaccurate answer or a
             * divergence of the system.
             *
             * This routine identifies such cases for 2D, and also for 3D where
             * e.g. edges may have the same problem (consider an extrusion of
             * the case above, for example).
             */

            // Collate information on Dirichlet vertices from all processes
            int n = vComm->GetSize();
            int p = vComm->GetRank();

            if(vComm->IsSerial())
            {
                // for FieldConvert Comm this is true and it resets
                // parallel processing back to serial case
                n = 1;
                p = 0;
            }
            // At this point, graph only contains information from Dirichlet
            // boundaries. Therefore make a global list of the vert and edge
            // information on all processors.
            Array<OneD, int> vertcounts (n, 0);
            Array<OneD, int> vertoffsets(n, 0);
            Array<OneD, int> edgecounts (n, 0);
            Array<OneD, int> edgeoffsets(n, 0);
            vertcounts[p] = graph[0].size();
            edgecounts[p] = graph[1].size();
            vComm->AllReduce(vertcounts, LibUtilities::ReduceSum);
            vComm->AllReduce(edgecounts, LibUtilities::ReduceSum);

            for (i = 1; i < n; ++i)
            {
                vertoffsets[i] = vertoffsets[i-1] + vertcounts[i-1];
                edgeoffsets[i] = edgeoffsets[i-1] + edgecounts[i-1];
            }

            int nTotVerts = Vmath::Vsum(n,vertcounts,1);
            int nTotEdges = Vmath::Vsum(n,edgecounts,1);

            Array<OneD, int> vertlist(nTotVerts, 0);
            Array<OneD, int> edgelist(nTotEdges, 0);

            // construct list of global ids of global vertices
            i = 0;
            for (auto &it : graph[0])
            {
                vertlist[vertoffsets[p] + i++] = it.first;
            }

            // construct list of global ids of global edges
            i = 0;
            for (auto &it : graph[1])
            {
                edgelist[edgeoffsets[p] + i++] = it.first;
            }
            vComm->AllReduce(vertlist, LibUtilities::ReduceSum);
            vComm->AllReduce(edgelist, LibUtilities::ReduceSum);

            // Now we have a list of all Dirichlet vertices and edges on all
            // processors.
            nExtraDirichlet = 0;
            map<int, int> extraDirVertIds, extraDirEdgeIds;

            // Ensure Dirchlet vertices are consistently recorded between
            // processes (e.g. Dirichlet region meets Neumann region across a
            // partition boundary requires vertex on partition to be Dirichlet).
            //
            // To do this we look over all elements and vertices in local
            // partition and see if they match the values stored in the vertlist
            // from other processors and if so record the meshVertId/meshEdgeId
            // and the processor it comes from.
            for (i = 0; i < n; ++i)
            {
                if (i == p)
                {
                    continue;
                }

                for(j = 0; j < locExpVector.size(); j++)
                {
                    exp = locExpVector[j];

                    for(k = 0; k < exp->GetNverts(); k++)
                    {
                        meshVertId = exp->GetGeom()->GetVid(k);
                        if(graph[0].count(meshVertId) == 0)
                        {
                            for (l = 0; l < vertcounts[i]; ++l)
                            {
                                if (vertlist[vertoffsets[i]+l] == meshVertId)
                                {
                                    extraDirVertIds[meshVertId] = i;
                                    graph[0][meshVertId] = graphVertId++;
                                    nExtraDirichlet++;
                                }
                            }
                        }
                    }

                    for(k = 0; k < exp->GetNedges(); k++)
                    {
                        meshEdgeId = exp->GetGeom()->GetEid(k);
                        if(graph[1].count(meshEdgeId) == 0)
                        {
                            for (l = 0; l < edgecounts[i]; ++l)
                            {
                                if (edgelist[edgeoffsets[i]+l] == meshEdgeId)
                                {
                                    extraDirEdgeIds[meshEdgeId] = i;
                                    graph[1][meshEdgeId] = graphVertId++;
                                    nExtraDirichlet += exp->GetEdgeNcoeffs(k)-2;
                                }
                            }
                        }
                    }
                }
            }

            // Low Energy preconditioner needs to know how many extra Dirichlet
            // edges are on this process so store map in array.
            m_extraDirEdges = Array<OneD, int>(extraDirEdgeIds.size(), -1);
            i = 0;
            for (auto &it : extraDirEdgeIds)
            {
                meshEdgeId = it.first;
                m_extraDirEdges[i++] = meshEdgeId;
            }

            // Now we have a list of all vertices and edges that are Dirichlet
            // and not defined on the local partition as well as which processor
            // they are stored on.
            //
            // Make a full list of all such entities on all processors and which
            // processor they belong to.
            for (i = 0; i < n; ++i)
            {
                vertcounts [i] = 0;
                vertoffsets[i] = 0;
                edgecounts [i] = 0;
                edgeoffsets[i] = 0;
            }

            vertcounts[p] = extraDirVertIds.size();
            edgecounts[p] = extraDirEdgeIds.size();
            vComm->AllReduce(vertcounts, LibUtilities::ReduceSum);
            vComm->AllReduce(edgecounts, LibUtilities::ReduceSum);
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

            i = 0;
            for (auto &it : extraDirVertIds)
            {
                vertids  [vertoffsets[p]+i] = it.first;
                vertprocs[vertoffsets[p]+i] = it.second;
                ++i;
            }

            i = 0;
            for (auto &it : extraDirEdgeIds)
            {
                edgeids  [edgeoffsets[p]+i] = it.first;
                edgeprocs[edgeoffsets[p]+i] = it.second;
                ++i;
            }

            vComm->AllReduce(vertids,   LibUtilities::ReduceSum);
            vComm->AllReduce(vertprocs, LibUtilities::ReduceSum);
            vComm->AllReduce(edgeids,   LibUtilities::ReduceSum);
            vComm->AllReduce(edgeprocs, LibUtilities::ReduceSum);

            // Set up list of vertices that need to be shared to other
            // partitions
            for (i = 0; i < nTotVerts; ++i)
            {
                if (p == vertprocs[i]) // rank = vertproc[i]
                {
                    extraDirVerts.insert(vertids[i]);
                }
            }

            // Set up list of edges that need to be shared to other partitions
            for (i = 0; i < nTotEdges; ++i)
            {
                if (p == edgeprocs[i]) // rank = vertproc[i]
                {
                    extraDirEdges.insert(edgeids[i]);
                }
            }

            // Check between processes if the whole system is singular
            int s = m_systemSingular ? 1 : 0;
            vComm->AllReduce(s, LibUtilities::ReduceMin);
            m_systemSingular = s == 1 ? true : false;

            // Find the minimum boundary vertex ID on each process
            Array<OneD, int> bcminvertid(n, 0);
            bcminvertid[p] = vMaxVertId;
            vComm->AllReduce(bcminvertid, LibUtilities::ReduceMax);

            // Find the process rank with the minimum boundary vertex ID
            int maxIdx = Vmath::Imax(n, bcminvertid, 1);

            // If the system is singular, the process with the maximum
            // number of BCs will set a Dirichlet vertex to make
            // system non-singular.  Note: we find the process with
            // maximum boundary regions to ensure we do not try to set
            // a Dirichlet vertex on a partition with no intersection
            // with the boundary.
            meshVertId = 0;

            if (m_systemSingular && checkIfSystemSingular && maxIdx == p)
            {
                if (m_session->DefinesParameter("SingularVertex"))
                {
                    m_session->LoadParameter("SingularVertex", meshVertId);
                }
                else if (vMaxVertId == -1)
                {
                    // All boundaries are periodic.
                    meshVertId = locExpVector[0]->GetGeom()->GetVid(0);
                }
                else
                {
                    // Set pinned vertex to that with minimum vertex ID to
                    // ensure consistency in parallel.
                    meshVertId = bcminvertid[p];
                }

                if (graph[0].count(meshVertId) == 0)
                {
                    graph[0][meshVertId] = graphVertId++;
                }
            }

            vComm->AllReduce(meshVertId, LibUtilities::ReduceSum);

            // When running in parallel, we need to ensure that the singular
            // mesh vertex is communicated to any periodic vertices, otherwise
            // the system may diverge.
            if(m_systemSingular && checkIfSystemSingular)
            {
                // Firstly, we check that no other processors have this
                // vertex. If they do, then we mark the vertex as also being
                // Dirichlet.
                if (maxIdx != p)
                {
                    for (i = 0; i < locExpVector.size(); ++i)
                    {
                        for (j = 0; j < locExpVector[i]->GetNverts(); ++j)
                        {
                            if (locExpVector[i]->GetGeom()->GetVid(j) !=
                                    meshVertId)
                            {
                                continue;
                            }

                            if (graph[0].count(meshVertId) == 0)
                            {
                                graph[0][meshVertId] =
                                    graphVertId++;
                            }
                        }
                    }
                }

                // In the case that meshVertId is periodic with other vertices,
                // this process and all other processes need to make sure that
                // the periodic vertices are also marked as Dirichlet.
                int gId;

                // At least one process (maxBCidx) will have already associated
                // a graphVertId with meshVertId. Others won't even have any of
                // the vertices. The logic below is designed to handle both
                // cases.
                if (graph[0].count(meshVertId) == 0)
                {
                    gId = -1;
                }
                else
                {
                    gId = graph[0][meshVertId];
                }

                for (auto &pIt : periodicVerts)
                {
                    // Either the vertex is local to this processor (in which
                    // case it will be in the pIt.first position) or else
                    // meshVertId might be contained within another processor's
                    // vertex list. The if statement below covers both cases. If
                    // we find it, set as Dirichlet with the vertex id gId.
                    if (pIt.first == meshVertId)
                    {
                        gId = gId < 0 ? graphVertId++ : gId;
                        graph[0][meshVertId] = gId;

                        for (i = 0; i < pIt.second.size(); ++i)
                        {
                            if (pIt.second[i].isLocal)
                            {
                                graph[0][pIt.second[i].id] = graph[0][meshVertId];
                            }
                        }
                    }
                    else
                    {
                        bool found = false;
                        for (i = 0; i < pIt.second.size(); ++i)
                        {
                            if (pIt.second[i].id == meshVertId)
                            {
                                found = true;
                                break;
                            }
                        }

                        if (found)
                        {
                            gId = gId < 0 ? graphVertId++ : gId;
                            graph[0][pIt.first] = gId;

                            for (i = 0; i < pIt.second.size(); ++i)
                            {
                                if (pIt.second[i].isLocal)
                                {
                                    graph[0][pIt.second[i].id] = graph[0][pIt.first];
                                }
                            }
                        }
                    }
                }
            }

            // Add extra dirichlet boundary conditions to count.
            m_numLocalDirBndCoeffs += nExtraDirichlet;
            firstNonDirGraphVertId  = graphVertId;

            typedef boost::adjacency_list<
                boost::setS, boost::vecS, boost::undirectedS> BoostGraph;
            BoostGraph boostGraphObj;

            vector<map<int,int> > tempGraph(3);
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

            m_numNonDirVertexModes = 0;
            m_numNonDirEdges       = 0;
            m_numNonDirFaces       = 0;
            m_numNonDirFaceModes   = 0;
            m_numNonDirFaceModes   = 0;
            m_numLocalBndCoeffs    = 0;

            map<int,int> EdgeSize;
            map<int,int> FaceSize;

            /// -  Count verts, edges, face and add up edges and face sizes
            for(i = 0; i < locExpVector.size(); ++i)
            {
                exp = locExpVector[i];
                nTotalVerts += exp->GetNverts();
                nTotalEdges += exp->GetNedges();
                nTotalFaces += exp->GetNfaces();

                nEdges = exp->GetNedges();
                for(j = 0; j < nEdges; ++j)
                {
                    meshEdgeId = exp->GetGeom()->GetEid(j);
                    if (EdgeSize.count(meshEdgeId) > 0)
                    {
                        EdgeSize[meshEdgeId] =
                                min(EdgeSize[meshEdgeId],
                                    exp->GetEdgeNcoeffs(j) - 2);
                    }
                    else
                    {
                        EdgeSize[meshEdgeId] = exp->GetEdgeNcoeffs(j) - 2;
                    }
                }

                nFaces = exp->GetNfaces();
                faceCnt = 0;
                for(j = 0; j < nFaces; ++j)
                {
                    meshFaceId = exp->GetGeom()->GetFid(j);
                    if (FaceSize.count(meshFaceId) > 0)
                    {
                        FaceSize[meshFaceId] =
                                min(FaceSize[meshFaceId],
                                    exp->GetFaceIntNcoeffs(j));
                    }
                    else
                    {
                        FaceSize[meshFaceId] = exp->GetFaceIntNcoeffs(j);
                    }
                    FaceSize[meshFaceId] = exp->GetFaceIntNcoeffs(j);
                }
            }

            /// - Periodic vertices
            for (auto &pIt : periodicVerts)
            {
                meshVertId = pIt.first;

                // This periodic vertex is joined to a Dirichlet condition.
                if (graph[0].count(pIt.first) != 0)
                {
                    for (i = 0; i < pIt.second.size(); ++i)
                    {
                        meshVertId2 = pIt.second[i].id;
                        if (graph[0].count(meshVertId2) == 0 &&
                            pIt.second[i].isLocal)
                        {
                            graph[0][meshVertId2] =
                                graph[0][meshVertId];
                        }
                    }
                    continue;
                }

                // One of the attached vertices is Dirichlet.
                bool isDirichlet = false;
                for (i = 0; i < pIt.second.size(); ++i)
                {
                    if (!pIt.second[i].isLocal)
                    {
                        continue;
                    }

                    meshVertId2 = pIt.second[i].id;
                    if (graph[0].count(meshVertId2) > 0)
                    {
                        isDirichlet = true;
                        break;
                    }
                }

                if (isDirichlet)
                {
                    graph[0][meshVertId] =
                        graph[0][pIt.second[i].id];

                    for (j = 0; j < pIt.second.size(); ++j)
                    {
                        meshVertId2 = pIt.second[i].id;
                        if (j == i || !pIt.second[j].isLocal ||
                            graph[0].count(meshVertId2) > 0)
                        {
                            continue;
                        }

                        graph[0][meshVertId2] =
                            graph[0][pIt.second[i].id];
                    }

                    continue;
                }

                // Otherwise, see if a vertex ID has already been set.
                for (i = 0; i < pIt.second.size(); ++i)
                {
                    if (!pIt.second[i].isLocal)
                    {
                        continue;
                    }

                    if (tempGraph[0].count(pIt.second[i].id) > 0)
                    {
                        break;
                    }
                }

                if (i == pIt.second.size())
                {
                    boost::add_vertex(boostGraphObj);
                    tempGraph[0][meshVertId] = tempGraphVertId++;
                    m_numNonDirVertexModes++;
                }
                else
                {
                    tempGraph[0][meshVertId] = tempGraph[0][pIt.second[i].id];
                }
            }

            // Store the temporary graph vertex id's of all element edges and
            // vertices in these 3 arrays below
            localVerts = Array<OneD, int>(nTotalVerts,-1);
            localEdges = Array<OneD, int>(nTotalEdges,-1);
            localFaces = Array<OneD, int>(nTotalFaces,-1);

            // Set up vertex numbering
            for(i = 0; i < locExpVector.size(); ++i)
            {
                exp = locExpVector[i];
                vertCnt = 0;
                nVerts = exp->GetNverts();
                for(j = 0; j < nVerts; ++j)
                {
                    meshVertId = exp->GetGeom()->GetVid(j);
                    if(graph[0].count(meshVertId) == 0)
                    {
                        if(tempGraph[0].count(meshVertId) == 0)
                        {
                            boost::add_vertex(boostGraphObj);
                            tempGraph[0][meshVertId] = tempGraphVertId++;
                            m_numNonDirVertexModes+=1;
                        }
                        localVerts[localVertOffset+vertCnt++] = tempGraph[0][meshVertId];
                        vwgts_map[ tempGraph[0][meshVertId] ] = 1;
                    }
                }

                localVertOffset+=nVerts;
            }

            /// - Periodic edges
            for (auto &pIt : periodicEdges)
            {
                meshEdgeId = pIt.first;

                // This periodic edge is joined to a Dirichlet condition.
                if (graph[1].count(pIt.first) != 0)
                {
                    for (i = 0; i < pIt.second.size(); ++i)
                    {
                        meshEdgeId2 = pIt.second[i].id;
                        if (graph[1].count(meshEdgeId2) == 0 &&
                            pIt.second[i].isLocal)
                        {
                            graph[1][meshEdgeId2] =
                                graph[1][meshEdgeId];
                        }
                    }
                    continue;
                }

                // One of the attached edges is Dirichlet.
                bool isDirichlet = false;
                for (i = 0; i < pIt.second.size(); ++i)
                {
                    if (!pIt.second[i].isLocal)
                    {
                        continue;
                    }

                    meshEdgeId2 = pIt.second[i].id;
                    if (graph[1].count(meshEdgeId2) > 0)
                    {
                        isDirichlet = true;
                        break;
                    }
                }

                if (isDirichlet)
                {
                    graph[1][meshEdgeId] =
                        graph[1][pIt.second[i].id];

                    for (j = 0; j < pIt.second.size(); ++j)
                    {
                        meshEdgeId2 = pIt.second[i].id;
                        if (j == i || !pIt.second[j].isLocal ||
                            graph[1].count(meshEdgeId2) > 0)
                        {
                            continue;
                        }

                        graph[1][meshEdgeId2] =
                            graph[1][pIt.second[i].id];
                    }

                    continue;
                }

                // Otherwise, see if a edge ID has already been set.
                for (i = 0; i < pIt.second.size(); ++i)
                {
                    if (!pIt.second[i].isLocal)
                    {
                        continue;
                    }

                    if (tempGraph[1].count(pIt.second[i].id) > 0)
                    {
                        break;
                    }
                }

                if (i == pIt.second.size())
                {
                    boost::add_vertex(boostGraphObj);
                    tempGraph[1][meshEdgeId] = tempGraphVertId++;
                    m_numNonDirEdgeModes += EdgeSize[meshEdgeId];
                    m_numNonDirEdges++;
                }
                else
                {
                    tempGraph[1][meshEdgeId] = tempGraph[1][pIt.second[i].id];
                }
            }

            int nEdgeIntCoeffs, nFaceIntCoeffs;

            // Set up edge numbering
            for(i = 0; i < locExpVector.size(); ++i)
            {
                exp = locExpVector[i];
                edgeCnt = 0;
                nEdges = exp->GetNedges();

                for(j = 0; j < nEdges; ++j)
                {
                    meshEdgeId = exp->GetGeom()->GetEid(j);
                    nEdgeIntCoeffs = EdgeSize[meshEdgeId];
                    if(graph[1].count(meshEdgeId) == 0)
                    {
                        if(tempGraph[1].count(meshEdgeId) == 0)
                        {
                            boost::add_vertex(boostGraphObj);
                            tempGraph[1][meshEdgeId] = tempGraphVertId++;
                            m_numNonDirEdgeModes+=nEdgeIntCoeffs;

                            m_numNonDirEdges++;
                        }
                        localEdges[localEdgeOffset+edgeCnt++] = tempGraph[1][meshEdgeId];
                        vwgts_map[ tempGraph[1][meshEdgeId] ] = nEdgeIntCoeffs;
                    }
                }

                localEdgeOffset+=nEdges;
            }

            /// - Periodic faces
            for (auto &pIt : periodicFaces)
            {
                if (!pIt.second[0].isLocal)
                {
                    // The face mapped to is on another process.
                    meshFaceId = pIt.first;
                    ASSERTL0(graph[2].count(meshFaceId) == 0,
                             "This periodic boundary edge has been specified before");
                    boost::add_vertex(boostGraphObj);
                    tempGraph[2][meshFaceId]  = tempGraphVertId++;
                    nFaceIntCoeffs  = FaceSize[meshFaceId];
                    m_numNonDirFaceModes+=nFaceIntCoeffs;
                    m_numNonDirFaces++;
                }
                else if (pIt.first < pIt.second[0].id)
                {
                    ASSERTL0(graph[2].count(pIt.first) == 0,
                             "This periodic boundary face has been specified before");
                    ASSERTL0(graph[2].count(pIt.second[0].id) == 0,
                             "This periodic boundary face has been specified before");

                    boost::add_vertex(boostGraphObj);
                    tempGraph[2][pIt.first]        = tempGraphVertId;
                    tempGraph[2][pIt.second[0].id] = tempGraphVertId++;
                    nFaceIntCoeffs  = FaceSize[pIt.first];
                    m_numNonDirFaceModes+=nFaceIntCoeffs;
                    m_numNonDirFaces++;
                }
            }

            // setup face numbering
            for(i = 0; i < locExpVector.size(); ++i)
            {
                exp = locExpVector[i];
                nFaces = exp->GetNfaces();
                faceCnt = 0;
                for(j = 0; j < nFaces; ++j)
                {
                    nFaceIntCoeffs = exp->GetFaceIntNcoeffs(j);
                    meshFaceId = exp->GetGeom()->GetFid(j);
                    if(graph[2].count(meshFaceId) == 0)
                    {
                        if(tempGraph[2].count(meshFaceId) == 0)
                        {
                            boost::add_vertex(boostGraphObj);
                            tempGraph[2][meshFaceId] = tempGraphVertId++;
                            m_numNonDirFaceModes+=nFaceIntCoeffs;

                            m_numNonDirFaces++;
                        }
                        localFaces[localFaceOffset+faceCnt++] = tempGraph[2][meshFaceId];
                        vwgts_map[ tempGraph[2][meshFaceId] ] = nFaceIntCoeffs;
                    }
                }
                m_numLocalBndCoeffs += exp->NumBndryCoeffs();

                localFaceOffset+=nFaces;
            }

            localVertOffset=0;
            localEdgeOffset=0;
            localFaceOffset=0;
            for(i = 0; i < locExpVector.size(); ++i)
            {
                exp = locExpVector[i];
                nVerts = exp->GetNverts();
                nEdges = exp->GetNedges();
                nFaces = exp->GetNfaces();

                // Now loop over all local faces, edges and vertices of this
                // element and define that all other faces, edges and verices of
                // this element are adjacent to them.

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

                localVertOffset+=nVerts;
                localEdgeOffset+=nEdges;
                localFaceOffset+=nFaces;
            }

            // Container to store vertices of the graph which correspond to
            // degrees of freedom along the boundary and periodic BCs.
            set<int> partVerts;

            if (m_solnType == eIterativeMultiLevelStaticCond ||
                m_solnType == eXxtMultiLevelStaticCond)
            {
                vector<long> procVerts,  procEdges,  procFaces;
                set   <int>  foundVerts, foundEdges, foundFaces;

                // Loop over element and construct the procVerts and procEdges
                // vectors, which store the geometry IDs of mesh vertices and
                // edges respectively which are local to this process.
                for(i = cnt = 0; i < locExpVector.size(); ++i)
                {
                    int elmtid = i;
                    exp = locExpVector[elmtid];
                    for (j = 0; j < exp->GetNverts(); ++j)
                    {
                        int vid = exp->GetGeom()->GetVid(j)+1;
                        if (foundVerts.count(vid) == 0)
                        {
                            procVerts.push_back(vid);
                            foundVerts.insert(vid);
                        }
                    }

                    for (j = 0; j < exp->GetNedges(); ++j)
                    {
                        int eid = exp->GetGeom()->GetEid(j)+1;

                        if (foundEdges.count(eid) == 0)
                        {
                            procEdges.push_back(eid);
                            foundEdges.insert(eid);
                        }
                    }

                    for (j = 0; j < exp->GetNfaces(); ++j)
                    {
                        int fid = exp->GetGeom()->GetFid(j)+1;

                        if (foundFaces.count(fid) == 0)
                        {
                            procFaces.push_back(fid);
                            foundFaces.insert(fid);
                        }
                    }
                }

                int unique_verts = foundVerts.size();
                int unique_edges = foundEdges.size();
                int unique_faces = foundFaces.size();

                bool verbose = m_session->DefinesCmdLineArgument("verbose");

                // Now construct temporary GS objects. These will be used to
                // populate the arrays tmp3 and tmp4 with the multiplicity of
                // the vertices and edges respectively to identify those
                // vertices and edges which are located on partition boundary.
                Array<OneD, long> vertArray(unique_verts, &procVerts[0]);
                Gs::gs_data *tmp1 = Gs::Init(vertArray, vComm, verbose);
                Array<OneD, NekDouble> tmp4(unique_verts, 1.0);
                Array<OneD, NekDouble> tmp5(unique_edges, 1.0);
                Array<OneD, NekDouble> tmp6(unique_faces, 1.0);
                Gs::Gather(tmp4, Gs::gs_add, tmp1);
                Gs::Finalise(tmp1);

                if (unique_edges > 0)
                {
                    Array<OneD, long> edgeArray(unique_edges, &procEdges[0]);
                    Gs::gs_data *tmp2 = Gs::Init(edgeArray, vComm, verbose);
                    Gs::Gather(tmp5, Gs::gs_add, tmp2);
                    Gs::Finalise(tmp2);
                }

                if (unique_faces > 0)
                {
                    Array<OneD, long> faceArray(unique_faces, &procFaces[0]);
                    Gs::gs_data *tmp3 = Gs::Init(faceArray, vComm, verbose);
                    Gs::Gather(tmp6, Gs::gs_add, tmp3);
                    Gs::Finalise(tmp3);
                }

                // Finally, fill the partVerts set with all non-Dirichlet
                // vertices which lie on a partition boundary.
                for (i = 0; i < unique_verts; ++i)
                {
                    if (tmp4[i] > 1.0)
                    {
                        if (graph[0].count(procVerts[i]-1) == 0)
                        {
                            partVerts.insert(tempGraph[0][procVerts[i]-1]);
                        }
                    }
                }

                for (i = 0; i < unique_edges; ++i)
                {
                    if (tmp5[i] > 1.0)
                    {
                        if (graph[1].count(procEdges[i]-1) == 0)
                        {
                            partVerts.insert(tempGraph[1][procEdges[i]-1]);
                        }
                    }
                }

                for (i = 0; i < unique_faces; ++i)
                {
                    if (tmp6[i] > 1.0)
                    {
                        if (graph[2].count(procFaces[i]-1) == 0)
                        {
                            partVerts.insert(tempGraph[2][procFaces[i]-1]);
                        }
                    }
                }

                // Now fill with all vertices on periodic BCs
                for (auto &pIt : periodicVerts)
                {
                    if (graph[0].count(pIt.first) == 0)
                    {
                        partVerts.insert(tempGraph[0][pIt.first]);
                    }
                }
                for (auto &pIt : periodicEdges)
                {
                    if (graph[1].count(pIt.first) == 0)
                    {
                        partVerts.insert(tempGraph[1][pIt.first]);
                    }
                }
                for (auto &pIt : periodicFaces)
                {
                    if (graph[2].count(pIt.first) == 0)
                    {
                        partVerts.insert(tempGraph[2][pIt.first]);
                    }
                }
            }

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
                    case ePETScStaticCond:
                    case ePETScFullMatrix:
                    case eXxtFullMatrix:
                    case eXxtStaticCond:
                    {
                        NoReordering(boostGraphObj,perm,iperm);
                        break;
                    }

                    case eDirectStaticCond:
                    {
                        CuthillMckeeReordering(boostGraphObj,perm,iperm);
                        break;
                    }

                    case ePETScMultiLevelStaticCond:
                    case eDirectMultiLevelStaticCond:
                    case eIterativeMultiLevelStaticCond:
                    case eXxtMultiLevelStaticCond:
                    {
                        MultiLevelBisectionReordering(
                            boostGraphObj, perm, iperm, bottomUpGraph,
                            partVerts, mdswitch);
                        break;
                    }
                    default:
                    {
                        ASSERTL0(false,
                                 "Unrecognised solution type: " + std::string(
                                     GlobalSysSolnTypeMap[m_solnType]));
                    }
                }
            }

            // For parallel multi-level static condensation determine the lowest
            // static condensation level amongst processors.
            if ((m_solnType == eDirectMultiLevelStaticCond ||
                 m_solnType == ePETScMultiLevelStaticCond ||
                 m_solnType == eIterativeMultiLevelStaticCond ||
                 m_solnType == eXxtMultiLevelStaticCond) && bottomUpGraph)
            {
                m_lowestStaticCondLevel = bottomUpGraph->GetNlevels()-1;
                vComm->AllReduce(m_lowestStaticCondLevel,
                                 LibUtilities::ReduceMax);
            }
            else
            {
                m_lowestStaticCondLevel = 0;
            }

            /**
             * STEP 4: Fill the #graph[0] and
             * #graph[1] with the optimal ordering from boost.
             */
            for(auto &mapIt : tempGraph[0])
            {
                graph[0][mapIt.first] = iperm[mapIt.second] + graphVertId;
            }
            for(auto &mapIt : tempGraph[1])
            {
                graph[1][mapIt.first] = iperm[mapIt.second] + graphVertId;
            }
            for(auto &mapIt : tempGraph[2])
            {
                graph[2][mapIt.first] = iperm[mapIt.second] + graphVertId;
            }

            return nGraphVerts;
        }

        /**
         *
         */
        AssemblyMapCG::AssemblyMapCG(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const int                                   numLocalCoeffs,
            const ExpList                              &locExp,
            const BndCondExp                           &bndCondExp,
            const BndCond                              &bndConditions,
            const bool                                  checkIfSystemSingular,
            const std::string                           variable,
            const PeriodicMap                          &periodicVerts,
            const PeriodicMap                          &periodicEdges,
            const PeriodicMap                          &periodicFaces)
            : AssemblyMap(pSession, variable)
        {
            int i, j, k;
            int p, q, numModes0, numModes1;
            int cnt = 0;
            int meshVertId, meshEdgeId, meshEdgeId2, meshFaceId, meshFaceId2;
            int globalId;
            int nEdgeInteriorCoeffs;
            int firstNonDirGraphVertId;
            LibUtilities::CommSharedPtr         vComm = m_comm->GetRowComm();
            LocalRegions::ExpansionSharedPtr    exp, bndExp;
            LibUtilities::BasisType             bType;
            StdRegions::Orientation             edgeOrient;
            StdRegions::Orientation             faceOrient;
            Array<OneD, unsigned int>           edgeInteriorMap;
            Array<OneD, int>                    edgeInteriorSign;
            Array<OneD, unsigned int>           faceInteriorMap;
            Array<OneD, int>                    faceInteriorSign;

            const LocalRegions::ExpansionVector &locExpVector = *(locExp.GetExp());

            bool verbose = m_session->DefinesCmdLineArgument("verbose");

            m_signChange = false;

            // Stores vertex, edge and face reordered vertices.
            DofGraph graph(3);
            DofGraph dofs(3);
            vector<map<int, int> > faceModes(2);
            map<int, LibUtilities::ShapeType> faceType;

            set<int> extraDirVerts, extraDirEdges;
            BottomUpSubStructuredGraphSharedPtr bottomUpGraph;

            // Construct list of number of degrees of freedom for each vertex,
            // edge and face.
            for (i = 0; i < locExpVector.size(); ++i)
            {
                exp = locExpVector[i];

                for(j = 0; j < exp->GetNverts(); ++j)
                {
                    dofs[0][exp->GetGeom()->GetVid(j)] = 1;
                }

                for(j = 0; j < exp->GetNedges(); ++j)
                {
                    if (dofs[1].count(exp->GetGeom()->GetEid(j)) > 0)
                    {
                        if (dofs[1][exp->GetGeom()->GetEid(j)] !=
                                    exp->GetEdgeNcoeffs(j)-2)
                        {
                            ASSERTL0( (exp->GetEdgeBasisType(j) == LibUtilities::eModified_A) ||
                                      (exp->GetEdgeBasisType(j) == LibUtilities::eModified_B) ||
                                      (exp->GetEdgeBasisType(j) == LibUtilities::eModified_C) ||
                                      (exp->GetEdgeBasisType(j) == LibUtilities::eModifiedPyr_C),
                                      "CG with variable order only available with modal expansion");
                        }
                        dofs[1][exp->GetGeom()->GetEid(j)] =
                                min(dofs[1][exp->GetGeom()->GetEid(j)],
                                    exp->GetEdgeNcoeffs(j)-2);
                    }
                    else
                    {
                        dofs[1][exp->GetGeom()->GetEid(j)] =
                            exp->GetEdgeNcoeffs(j) - 2;
                    }
                }

                for(j = 0; j < exp->GetNfaces(); ++j)
                {
                    faceOrient  = exp->GetGeom()->GetForient(j);
                    meshFaceId  = exp->GetGeom()->GetFid(j);
                    exp->GetFaceNumModes(j, faceOrient, numModes0, numModes1);

                    if (faceModes[0].count(meshFaceId) > 0)
                    {
                        faceModes[0][meshFaceId] =
                                min(faceModes[0][meshFaceId], numModes0);

                        faceModes[1][meshFaceId] =
                                min(faceModes[1][meshFaceId], numModes1);
                    }
                    else
                    {
                        faceModes[0][meshFaceId] = numModes0;
                        faceModes[1][meshFaceId] = numModes1;

                        // Get shape of this face
                        SpatialDomains::Geometry3DSharedPtr geom;
                        geom = std::dynamic_pointer_cast<SpatialDomains::
                                Geometry3D> (exp->GetGeom());
                        faceType[meshFaceId] =
                                geom->GetFace(j)->GetShapeType();
                    }
                }
            }

            // Add non-local periodic dofs to the map
            for (auto &pIt : periodicEdges)
            {
                for (i = 0; i < pIt.second.size(); ++i)
                {
                    meshEdgeId2 = pIt.second[i].id;
                    if (dofs[1].count(meshEdgeId2) == 0)
                    {
                        dofs[1][meshEdgeId2] = 1e6;
                    }
                }
            }
            for (auto &pIt : periodicFaces)
            {
                for (i = 0; i < pIt.second.size(); ++i)
                {
                    meshFaceId2 = pIt.second[i].id;
                    if (faceModes[0].count(meshFaceId2) == 0)
                    {
                        faceModes[0][meshFaceId2] = 1e6;
                        faceModes[1][meshFaceId2] = 1e6;
                    }
                }
            }

            // Now use information from all partitions to determine the correct
            // size

            // edges
            Array<OneD, long> edgeId (dofs[1].size());
            Array<OneD, NekDouble> edgeDof (dofs[1].size());
            i = 0;
            for(auto &dofIt : dofs[1])
            {
                edgeId [i  ] = dofIt.first + 1;
                edgeDof[i++] = (NekDouble) dofIt.second;
            }
            Gs::gs_data *tmp = Gs::Init(edgeId, vComm, verbose);
            Gs::Gather(edgeDof, Gs::gs_min, tmp);
            Gs::Finalise(tmp);
            for (i = 0; i < dofs[1].size(); i++)
            {
                dofs[1][edgeId[i]-1] = (int) (edgeDof[i]+0.5);
            }
            // Periodic edges
            for (auto &pIt : periodicEdges)
            {
                meshEdgeId = pIt.first;
                for (i = 0; i < pIt.second.size(); ++i)
                {
                    meshEdgeId2 = pIt.second[i].id;
                    if (dofs[1][meshEdgeId2] < dofs[1][meshEdgeId])
                    {
                        dofs[1][meshEdgeId] = dofs[1][meshEdgeId2];
                    }
                }
            }
            // faces
            Array<OneD, long> faceId (faceModes[0].size());
            Array<OneD, NekDouble> faceP (faceModes[0].size());
            Array<OneD, NekDouble> faceQ (faceModes[0].size());

            i = 0;
            for(auto dofIt = faceModes[0].begin(), dofIt2 = faceModes[1].begin();
                dofIt != faceModes[0].end(); dofIt++, dofIt2++, i++)
            {
                faceId[i] = dofIt->first+1;
                faceP[i] = (NekDouble) dofIt->second;
                faceQ[i] = (NekDouble) dofIt2->second;
            }
            Gs::gs_data *tmp2 = Gs::Init(faceId, vComm, verbose);
            Gs::Gather(faceP, Gs::gs_min, tmp2);
            Gs::Gather(faceQ, Gs::gs_min, tmp2);
            Gs::Finalise(tmp2);
            for (i=0; i < faceModes[0].size(); i++)
            {
                faceModes[0][faceId[i]-1] = (int) (faceP[i]+0.5);
                faceModes[1][faceId[i]-1] = (int) (faceQ[i]+0.5);
            }
            // Periodic faces
            for (auto &pIt : periodicFaces)
            {
                meshFaceId = pIt.first;
                for (i = 0; i < pIt.second.size(); ++i)
                {
                    meshFaceId2 = pIt.second[i].id;
                    if (faceModes[0][meshFaceId2] < faceModes[0][meshFaceId])
                    {
                        faceModes[0][meshFaceId] = faceModes[0][meshFaceId2];
                    }
                    if (faceModes[1][meshFaceId2] < faceModes[1][meshFaceId])
                    {
                        faceModes[1][meshFaceId] = faceModes[1][meshFaceId2];
                    }
                }
            }
            // Calculate number of dof in each face
            int P, Q;
            for (i=0; i < faceModes[0].size(); i++)
            {
                P = faceModes[0][faceId[i]-1];
                Q = faceModes[1][faceId[i]-1];
                if (faceType[faceId[i]-1] == LibUtilities::eQuadrilateral)
                {
                    // Quad face
                    dofs[2][faceId[i]-1] =
                      LibUtilities::StdQuadData::getNumberOfCoefficients(P,Q) -
                      LibUtilities::StdQuadData::getNumberOfBndCoefficients(P,Q);
                }
                else
                {
                    // Tri face
                    dofs[2][faceId[i]-1] =
                      LibUtilities::StdTriData::getNumberOfCoefficients(P,Q) -
                      LibUtilities::StdTriData::getNumberOfBndCoefficients(P,Q);
                }
            }

            Array<OneD, const BndCond> bndCondVec(1, bndConditions);

            // Note that nExtraDirichlet is not used in the logic below; it just
            // needs to be set so that the coupled solver in
            // IncNavierStokesSolver can work.
            int nExtraDirichlet;
            int mdswitch;
            m_session->LoadParameter(
                "MDSwitch", mdswitch, 10);

            int nGraphVerts =
                CreateGraph(locExp, bndCondExp, bndCondVec,
                            checkIfSystemSingular, periodicVerts, periodicEdges,
                            periodicFaces, graph, bottomUpGraph, extraDirVerts,
                            extraDirEdges, firstNonDirGraphVertId,
                            nExtraDirichlet, mdswitch);

            /*
             * Set up an array which contains the offset information of the
             * different graph vertices.
             *
             * This basically means to identify to how many global degrees of
             * freedom the individual graph vertices correspond. Obviously,
             * the graph vertices corresponding to the mesh-vertices account
             * for a single global DOF. However, the graph vertices
             * corresponding to the element edges correspond to N-2 global DOF
             * where N is equal to the number of boundary modes on this edge.
             */
            Array<OneD, int> graphVertOffset(
                graph[0].size() + graph[1].size() + graph[2].size() + 1);

            graphVertOffset[0] = 0;

            for(i = 0; i < locExpVector.size(); ++i)
            {
                exp = locExpVector[i];
                
                for(j = 0; j < exp->GetNverts(); ++j)
                {
                    meshVertId = exp->GetGeom()->GetVid(j);
                    graphVertOffset[graph[0][meshVertId]+1] = 1;
                }

                for(j = 0; j < exp->GetNedges(); ++j)
                {
                    nEdgeInteriorCoeffs = exp->GetEdgeNcoeffs(j) - 2;
                    meshEdgeId = exp->GetGeom()->GetEid(j);
                    graphVertOffset[graph[1][meshEdgeId]+1]
                        = dofs[1][meshEdgeId];

                    bType = exp->GetEdgeBasisType(j);

                    // need a sign vector for modal expansions if nEdgeCoeffs
                    // >=3 (not 4 because of variable order case)
                    if(nEdgeInteriorCoeffs &&
                       (bType == LibUtilities::eModified_A ||
                        bType == LibUtilities::eModified_B))
                    {
                        m_signChange = true;
                    }
                }

                for(j = 0; j < exp->GetNfaces(); ++j)
                {
                    meshFaceId = exp->GetGeom()->GetFid(j);
                    graphVertOffset[graph[2][meshFaceId]+1] =
                        dofs[2][meshFaceId];
                }
            }

            for(i = 1; i < graphVertOffset.size(); i++)
            {
                graphVertOffset[i] += graphVertOffset[i-1];
            }

            // Allocate the proper amount of space for the class-data
            m_numLocalCoeffs        = numLocalCoeffs;
            m_numGlobalDirBndCoeffs = graphVertOffset[firstNonDirGraphVertId];
            m_localToGlobalMap      = Array<OneD, int>(m_numLocalCoeffs,-1);
            m_localToGlobalBndMap   = Array<OneD, int>(m_numLocalBndCoeffs,-1);
            m_localToLocalBndMap    = Array<OneD, int>(m_numLocalBndCoeffs,-1);
            m_localToLocalIntMap    = Array<OneD, int>(m_numLocalCoeffs-
                                                       m_numLocalBndCoeffs,-1);
            m_bndCondCoeffsToLocalCoeffsMap  = Array<OneD, int>
                (m_numLocalBndCondCoeffs,-1);
            
            // If required, set up the sign-vector
            if(m_signChange)
            {
                m_localToGlobalSign              = Array<OneD, NekDouble>
                    (m_numLocalCoeffs,1.0);
                m_localToGlobalBndSign           = Array<OneD, NekDouble>
                    (m_numLocalBndCoeffs,1.0);
                m_bndCondCoeffsToLocalCoeffsSign = Array<OneD, NekDouble>
                    (m_numLocalBndCondCoeffs,1.0);
            }

            m_staticCondLevel = 0;
            m_numPatches =  locExpVector.size();
            m_numLocalBndCoeffsPerPatch = Array<OneD, unsigned int>(m_numPatches);
            m_numLocalIntCoeffsPerPatch = Array<OneD, unsigned int>(m_numPatches);
            for(i = 0; i < m_numPatches; ++i)
            {
                m_numLocalBndCoeffsPerPatch[i] = (unsigned int)
                    locExpVector[i]->NumBndryCoeffs();
                m_numLocalIntCoeffsPerPatch[i] = (unsigned int)
                    locExpVector[i]->GetNcoeffs() -
                    locExpVector[i]->NumBndryCoeffs();
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
            int cntbdry = 0; 
            int cntint  = 0; 
            for(i = 0; i < locExpVector.size(); ++i)
            {
                exp = locExpVector[i];
                cnt = locExp.GetCoeff_Offset(i);

                int nbdry = exp->NumBndryCoeffs();
                int nint  = exp->GetNcoeffs() - nbdry;

                Array<OneD,unsigned int> bmap(nbdry);
                Array<OneD,unsigned int> imap(nint);

                exp->GetBoundaryMap(bmap);
                exp->GetInteriorMap(imap);

                for(j = 0; j < nbdry; ++j)
                {
                    m_localToLocalBndMap[cntbdry++] = cnt + bmap[j];
                }

                for(j = 0; j < nint; ++j)
                {
                    m_localToLocalIntMap[cntint++] = cnt + imap[j];
                }

                for(j = 0; j < exp->GetNverts(); ++j)
                {
                    meshVertId = exp->GetGeom()->GetVid(j);

                    // Set the global DOF for vertex j of element i
                    m_localToGlobalMap[cnt+exp->GetVertexMap(j)] =
                        graphVertOffset[graph[0][meshVertId]];
                }

                for(j = 0; j < exp->GetNedges(); ++j)
                {
                    nEdgeInteriorCoeffs = exp->GetEdgeNcoeffs(j)-2;
                    edgeOrient          = exp->GetGeom()->GetEorient(j);
                    meshEdgeId          = exp->GetGeom()->GetEid(j);

                    auto pIt = periodicEdges.find(meshEdgeId);

                    // See if this edge is periodic. If it is, then we map all
                    // edges to the one with lowest ID, and align all
                    // coefficients to this edge orientation.
                    if (pIt != periodicEdges.end())
                    {
                        pair<int, StdRegions::Orientation> idOrient =
                            DeterminePeriodicEdgeOrientId(
                                meshEdgeId, edgeOrient, pIt->second);
                        edgeOrient = idOrient.second;
                    }

                    exp->GetEdgeInteriorMap(j,edgeOrient,edgeInteriorMap,edgeInteriorSign);

                    // Set the global DOF's for the interior modes of edge j
                    for(k = 0; k < dofs[1][meshEdgeId]; ++k)
                    {
                        m_localToGlobalMap[cnt+edgeInteriorMap[k]] =
                            graphVertOffset[graph[1][meshEdgeId]]+k;
                    }
                    for(k = dofs[1][meshEdgeId]; k < nEdgeInteriorCoeffs; ++k)
                    {
                        m_localToGlobalMap[cnt+edgeInteriorMap[k]] = 0;
                    }

                    // Fill the sign vector if required
                    if(m_signChange)
                    {
                        for(k = 0; k < dofs[1][meshEdgeId]; ++k)
                        {
                            m_localToGlobalSign[cnt+edgeInteriorMap[k]] =
                                (NekDouble) edgeInteriorSign[k];
                        }
                        for(k = dofs[1][meshEdgeId]; k < nEdgeInteriorCoeffs; ++k)
                        {
                            m_localToGlobalSign[cnt+edgeInteriorMap[k]] = 0.0;
                        }
                    }
                }

                for(j = 0; j < exp->GetNfaces(); ++j)
                {
                    faceOrient          = exp->GetGeom()->GetForient(j);
                    meshFaceId          = exp->GetGeom()->GetFid(j);

                    auto pIt = periodicFaces.find(meshFaceId);

                    if (pIt != periodicFaces.end() &&
                        meshFaceId == min(meshFaceId, pIt->second[0].id))
                    {
                        faceOrient = DeterminePeriodicFaceOrient(faceOrient,pIt->second[0].orient);
                    }

                    exp->GetFaceInteriorMap(j,faceOrient,faceInteriorMap,faceInteriorSign);

                    // Set the global DOF's for the interior modes of face j
                    exp->GetFaceNumModes(j, faceOrient, numModes0, numModes1);
                    switch(faceType[meshFaceId])
                    {
                    case LibUtilities::eQuadrilateral:
                    {
                        int kLoc=0;
                        k = 0;
                        for( q = 2; q < numModes1; q++)
                        {
                            for( p = 2; p < numModes0; p++)
                            {
                                if( (p < faceModes[0][meshFaceId]) &&
                                    (q < faceModes[1][meshFaceId]))
                                {
                                    m_localToGlobalMap[cnt+faceInteriorMap[kLoc]] =
                                        graphVertOffset[graph[2][meshFaceId]]+k;
                                    if(m_signChange)
                                    {
                                        m_localToGlobalSign[cnt+faceInteriorMap[kLoc]] =
                                            (NekDouble) faceInteriorSign[kLoc];
                                    }
                                    k++;
                                }
                                else
                                {
                                    m_localToGlobalMap[cnt+faceInteriorMap[kLoc]] =  0;
                                    if(m_signChange)
                                    {
                                        m_localToGlobalSign[cnt+faceInteriorMap[kLoc]] = 0.0;
                                    }
                                }
                                kLoc++;
                            }
                        }
                    }
                    break;
                    case LibUtilities::eTriangle:
                    {
                        int kLoc=0;
                        k = 0;
                        for( p = 2; p < numModes0; p++)
                        {
                            for( q = 1; q < numModes1-p; q++)
                            {
                                if( (p   < faceModes[0][meshFaceId]) &&
                                    (p+q < faceModes[1][meshFaceId]))
                                {
                                    m_localToGlobalMap[cnt+faceInteriorMap[kLoc]] =
                                        graphVertOffset[graph[2][meshFaceId]]+k;
                                    if(m_signChange)
                                    {
                                        m_localToGlobalSign[cnt+faceInteriorMap[kLoc]] =
                                            (NekDouble) faceInteriorSign[kLoc];
                                    }
                                    k++;
                                }
                                else
                                {
                                    m_localToGlobalMap[cnt+faceInteriorMap[kLoc]] = 0;
                                    if(m_signChange)
                                    {
                                        m_localToGlobalSign[cnt+faceInteriorMap[kLoc]] = 0.0;
                                    }
                                }
                                kLoc++;
                            }
                        }
                    }
                    break;
                    default:
                        ASSERTL0(false,"Shape not recognised");
                        break;
                    }
                }
            }

            // Set up the mapping for the boundary conditions
            // Set up boundary mapping
            map<int, pair<int,int> > traceToElmtTraceMap;
            int id;
            
            for(cnt = i = 0; i < locExpVector.size(); ++i)
            {
                exp = locExpVector[i];

                for(j = 0; j < exp->GetNtrace(); ++j)
                {
                    id = exp->GetGeom()->GetTid(j);
                    
                    traceToElmtTraceMap[id] = pair<int,int>(i,j);
                }
            }

            Array<OneD, unsigned int> maparray;
            Array<OneD, int>          signarray;
            map<int,pair<int,NekDouble>> GloDirBndCoeffToLocalCoeff; 
            set<int> CoeffOnDirTrace;
            
            cnt = 0;
            int offset = 0;
            for(i = 0; i < bndCondExp.size(); i++)
            {
                set<int> foundExtraVerts, foundExtraEdges;
                for(j = 0; j < bndCondExp[i]->GetNumElmts(); j++)
                {
                    bndExp  = bndCondExp[i]->GetExp(j);
                    cnt = offset + bndCondExp[i]->GetCoeff_Offset(j);

                    int id = bndExp->GetGeom()->GetGlobalID();

                    ASSERTL1(traceToElmtTraceMap.count(id) > 0,
                             "Failed to find trace id");

                    int eid = traceToElmtTraceMap[id].first;
                    int tid = traceToElmtTraceMap[id].second;

                    exp = locExpVector[eid]; 
                    int dim = exp->GetShapeDimension();
                    
                    if(dim == 1)
                    {
                        m_bndCondCoeffsToLocalCoeffsMap [cnt] =
                            locExp.GetCoeff_Offset(eid) + exp->GetVertexMap(tid);
                    }
                    else
                    {
                        if(dim == 2)
                        {
                            exp->GetEdgeToElementMap(tid,exp->GetEorient(tid),
                                                     maparray,signarray,
                                                     bndExp->GetBasisNumModes(0));
                        }
                        else if (dim == 3)
                        {
                            exp->GetFaceToElementMap(tid,exp->GetForient(tid),
                                                     maparray,signarray,
                                                     bndExp->GetBasisNumModes(0),
                                                     bndExp->GetBasisNumModes(1));
                        }
                        
                        for(k = 0; k < bndExp->GetNcoeffs(); k++)
                        {
                            m_bndCondCoeffsToLocalCoeffsMap [cnt+k] =
                                locExp.GetCoeff_Offset(eid) + maparray[k];
                            if(m_signChange) 
                            {
                                m_bndCondCoeffsToLocalCoeffsSign[cnt+k] = signarray[k];
                            }
                        }

                    }

                    // we now need some information to work out how to
                    // handle vertices and edges that are only just
                    // touching a dirichlet boundary (and not the
                    // whole edge/face)
                        
                    for(k = 0; k < bndExp->GetNcoeffs(); k++)
                    {
                        int locid = m_bndCondCoeffsToLocalCoeffsMap [cnt+k];
                        int gloid = m_localToGlobalMap[locid];
                        NekDouble sign = 1.0;

                        if(m_signChange)
                        {
                            sign = m_bndCondCoeffsToLocalCoeffsSign[cnt+k];
                        }

                        if (bndConditions[i]->GetBoundaryConditionType() ==
                            SpatialDomains::eDirichlet)
                        {
                            CoeffOnDirTrace.insert(locid);
                        
                            // store the local id and sign from global id
                            // back to local space;
                            GloDirBndCoeffToLocalCoeff[gloid] =
                                pair<int,NekDouble> (locid,sign);
                        }
                    }
                }
                offset += bndCondExp[i]->GetNcoeffs();
            }
            
            globalId = Vmath::Vmax(m_numLocalCoeffs,&m_localToGlobalMap[0],1)+1;
            m_numGlobalBndCoeffs = globalId;

            // Set up a mapping list of Dirichlet Local Dofs that
            // arise due to one vertex or edge just touching a
            // Dirichlet boundary and need the value from another
            // local coeff that has been filled by the boundary
            // coeffs.

            Array<OneD, NekDouble> gloParaDirBnd(m_numLocalCoeffs,-1.0);

            Array<OneD, unsigned int> bndmap; 
            cnt = 0; 
            for(i = 0; i < locExpVector.size(); ++i)
            {
                int gloid;
                
                exp = locExpVector[i];

                exp->GetBoundaryMap(bndmap);
                
                for(j = 0; j < bndmap.size(); ++j)
                {
                    k = cnt + bndmap[j];
                    
                    if(CoeffOnDirTrace.count(k) == 0) 
                    {
                        gloid = m_localToGlobalMap[k];
                
                        if(gloid < m_numGlobalDirBndCoeffs) // point on Dir BC
                        {
                            if(GloDirBndCoeffToLocalCoeff.count(gloid))
                            {
                                int locid = GloDirBndCoeffToLocalCoeff[gloid].
                                    first;
                                NekDouble sign = 1.0;
                                
                                if(m_signChange)
                                {
                                    sign = m_localToGlobalSign[locid]*
                                    m_localToGlobalSign[k];
                                }

                                ExtraDirDof  DirDofs(k,locid,sign);
                                // could make same `structure as extraDirDof
                                m_copyLocalDirDofs.insert(DirDofs);
                            } 
                            else // else could be on another parallel partition.
                            {
                                gloParaDirBnd[gloid] = gloid;
                            }
                        }
                    }
                }
                
                cnt += exp->GetNcoeffs();
            }
                        
            /*
             * The boundary condition mapping is generated from the same vertex
             * renumbering.
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

            SetUpUniversalC0ContMap(locExp, periodicVerts, periodicEdges,
                                    periodicFaces);

            // Now that universal map is setup reset gloParaDirBnd to
            // 0 if no point communicated or universal value of not
            // equal to -1.0
            for(i = 0; i < m_numGlobalBndCoeffs; ++i)
            {
                int gloid = gloParaDirBnd[i];
                if(gloid  == -1)
                {
                    gloParaDirBnd[i] = 0.0;
                }
                else
                {
                    gloParaDirBnd[i] = m_globalToUniversalMap[gloid]; 
                }
            }
            
            // Use parallel boundary communication to set parallel
            // dirichlet values on all processors Needs to be after
            // SetupUuniversialC0ContMap
            Gs::Gather(gloParaDirBnd,Gs::gs_max,m_bndGsh);
 
            // copy global ids back to local values in partition to initialise gs communicator. 
            Array<OneD, long> paraDirBnd(m_numLocalCoeffs);
            for(i = 0; i < numLocalCoeffs; ++i)
            {
                paraDirBnd[i] = 0.0; 

                int id = m_localToGlobalMap[i];

                if(id >= m_numGlobalDirBndCoeffs)
                {
                    continue;
                }

                paraDirBnd[i] = gloParaDirBnd[id];

                if(gloParaDirBnd[id] > 0.0)
                {
                    // gather any sign changes due to edge modes
                    if(m_signChange)
                    {
                        if(m_localToGlobalSign[i] < 0) 
                        {
                            m_parallelDirBndSign.insert(i);
                        }
                    }
                }
            }

            m_dirBndGsh = Gs::Init(paraDirBnd,vComm,verbose);
            
            // Set up the local to global map for the next level when using
            // multi-level static condensation
            if ((m_solnType == eDirectMultiLevelStaticCond    ||
                 m_solnType == eIterativeMultiLevelStaticCond ||
                 m_solnType == eXxtMultiLevelStaticCond       ||
                 m_solnType == ePETScMultiLevelStaticCond) && nGraphVerts)
            {
                if (m_staticCondLevel < (bottomUpGraph->GetNlevels()-1))
                {
                    Array<OneD, int> vwgts_perm(
                         graph[0].size() + graph[1].size() + graph[2].size()
                         - firstNonDirGraphVertId);
                    
                    for (i = 0; i < locExpVector.size(); ++i)
                    {
                        exp = locExpVector[i];
                        
                        for (j = 0; j < exp->GetNverts(); ++j)
                        {
                            meshVertId = exp->GetGeom()->GetVid(j);

                            if (graph[0][meshVertId] >= firstNonDirGraphVertId)
                            {
                                vwgts_perm[graph[0][meshVertId] -
                                           firstNonDirGraphVertId] =
                                    dofs[0][meshVertId];
                            }
                        }

                        for (j = 0; j < exp->GetNedges(); ++j)
                        {
                            meshEdgeId = exp->GetGeom()->GetEid(j);

                            if (graph[1][meshEdgeId] >= firstNonDirGraphVertId)
                            {
                                vwgts_perm[graph[1][meshEdgeId] -
                                           firstNonDirGraphVertId] =
                                    dofs[1][meshEdgeId];
                            }
                        }

                        for (j = 0; j < exp->GetNfaces(); ++j)
                        {
                            meshFaceId = exp->GetGeom()->GetFid(j);
                            
                            if (graph[2][meshFaceId] >= firstNonDirGraphVertId)
                            {
                                vwgts_perm[graph[2][meshFaceId] -
                                           firstNonDirGraphVertId] =
                                    dofs[2][meshFaceId];
                            }
                        }
                    }
                    
                    bottomUpGraph->ExpandGraphWithVertexWeights(vwgts_perm);
                    m_nextLevelLocalToGlobalMap = MemoryManager<AssemblyMap>::
                        AllocateSharedPtr(this, bottomUpGraph);
                }
            }
            
            m_hash = hash_range(m_localToGlobalMap.begin(),
                                m_localToGlobalMap.end());
            
            // Add up hash values if parallel
            int hash = m_hash;
            vComm->AllReduce(hash, LibUtilities::ReduceSum);
            m_hash = hash;
            
            CalculateBndSystemBandWidth();
            CalculateFullSystemBandWidth();
        }
        
        /**
         *
         */
        AssemblyMapCG::~AssemblyMapCG()
        {
            Gs::Finalise(m_gsh);
            Gs::Finalise(m_bndGsh);
        }

        /**
         * @brief Determine orientation of an edge to its periodic equivalents,
         * as well as the ID of the representative edge.
         *
         * Since an edge may be periodic with more than one other edge (e.g. a
         * periodic cube has sets of four periodic edges in each coordinate
         * direction), we have to define a 'representative' edge. In this
         * assembly map we define it to be the one with the minimum ID. This
         * routine is set up to calculate the orientation of a given edge with
         * ID @p meshEdgeId with respect to the edge ID.
         *
         * @param meshEdgeId     ID of a periodic edge.
         * @param edgeOrient     Edge orientation of meshEdgeId with respect to
         *                       its parent element.
         * @param periodicEdges  The map of all periodic edges.
         *
         * @return Pair containing the ID of the periodic edge and the
         *         orientation of @p meshEdgeID with respect to this edge.
         */
        pair<int, StdRegions::Orientation> DeterminePeriodicEdgeOrientId(
            int                           meshEdgeId,
            StdRegions::Orientation       edgeOrient,
            const vector<PeriodicEntity> &periodicEdges)
        {
            int minId  = periodicEdges[0].id;
            int minIdK = 0;
            int k;

            for (k = 1; k < periodicEdges.size(); ++k)
            {
                if (periodicEdges[k].id < minId)
                {
                    minId  = min(minId, periodicEdges[k].id);
                    minIdK = k;
                }
            }

            minId = min(minId, meshEdgeId);

            if (meshEdgeId != minId)
            {
                if (periodicEdges[minIdK].orient == StdRegions::eBackwards)
                {
                    // Swap edge orientation
                    edgeOrient = (edgeOrient == StdRegions::eForwards) ?
                        StdRegions::eBackwards : StdRegions::eForwards;
                }
            }

            return make_pair(minId, edgeOrient);
        }

        /**
         * @brief Determine relative orientation between two faces.
         *
         * Given the orientation of a local element to its local face, defined
         * as @p faceOrient, and @p perFaceOrient which states the alignment of
         * one periodic face to the other global face, this routine determines
         * the orientation that takes this local element face to the
         * global/unique face.
         *
         * @param faceOrient     Orientation of the face with respect to its
         *                       parent element.
         * @param perFaceOrient  Orientation of the representative/global face.
         *
         * @return Orientation between the two faces.
         */
        StdRegions::Orientation DeterminePeriodicFaceOrient(
            StdRegions::Orientation faceOrient,
            StdRegions::Orientation perFaceOrient)
        {
            StdRegions::Orientation  returnval = faceOrient;

            if(perFaceOrient != StdRegions::eDir1FwdDir1_Dir2FwdDir2)
            {
                int tmp1 = (int)faceOrient    - 5;
                int tmp2 = (int)perFaceOrient - 5;

                int flipDir1Map [8] = {2,3,0,1,6,7,4,5};
                int flipDir2Map [8] = {1,0,3,2,5,4,7,6};
                int transposeMap[8] = {4,5,6,7,0,2,1,3};

                // Transpose orientation
                if (tmp2 > 3)
                {
                    tmp1 = transposeMap[tmp1];
                }

                // Reverse orientation in direction 1.
                if (tmp2 == 2 || tmp2 == 3 || tmp2 == 6 || tmp2 == 7)
                {
                    tmp1 = flipDir1Map[tmp1];
                }

                // Reverse orientation in direction 2
                if (tmp2 % 2 == 1)
                {
                    tmp1 = flipDir2Map[tmp1];
                }

                returnval = (StdRegions::Orientation)(tmp1+5);
            }
            return returnval;
        }


        /**
         * Sets up the global to universal mapping of degrees of freedom across
         * processors.
         */
        void AssemblyMapCG::SetUpUniversalC0ContMap(
            const ExpList     &locExp,
            const PeriodicMap &perVerts,
            const PeriodicMap &perEdges,
            const PeriodicMap &perFaces)
        {
            LocalRegions::ExpansionSharedPtr exp;
            int nVert = 0;
            int nEdge = 0;
            int nFace = 0;
            int maxEdgeDof = 0;
            int maxFaceDof = 0;
            int maxIntDof = 0;
            int dof = 0;
            int cnt;
            int i,j,k,l;
            int meshVertId;
            int meshEdgeId;
            int meshFaceId;
            int elementId;
            int vGlobalId;
            int maxBndGlobalId = 0;
            StdRegions::Orientation     edgeOrient;
            StdRegions::Orientation     faceOrient;
            Array<OneD, unsigned int>   edgeInteriorMap;
            Array<OneD, int>            edgeInteriorSign;
            Array<OneD, unsigned int>   faceInteriorMap;
            Array<OneD, int>            faceInteriorSign;
            Array<OneD, unsigned int>   interiorMap;

            const LocalRegions::ExpansionVector &locExpVector = *(locExp.GetExp());
            LibUtilities::CommSharedPtr vCommRow = m_comm->GetRowComm();
            const bool verbose = locExp.GetSession()->DefinesCmdLineArgument("verbose");

            m_globalToUniversalMap = Nektar::Array<OneD, int>(m_numGlobalCoeffs, -1);
            m_globalToUniversalMapUnique = Nektar::Array<OneD, int>(m_numGlobalCoeffs, -1);
            m_globalToUniversalBndMap = Nektar::Array<OneD, int>(m_numGlobalBndCoeffs, -1);
            m_globalToUniversalBndMapUnique = Nektar::Array<OneD, int>(m_numGlobalBndCoeffs, -1);

            // Loop over all the elements in the domain to gather mesh data
            for(i = 0; i < locExpVector.size(); ++i)
            {
                exp = locExpVector[i];
                nVert += exp->GetNverts();
                nEdge += exp->GetNedges();
                nFace += exp->GetNfaces();
                // Loop over all edges (and vertices) of element i
                for(j = 0; j < exp->GetNedges(); ++j)
                {
                    dof = exp->GetEdgeNcoeffs(j)-2;
                    maxEdgeDof = (dof > maxEdgeDof ? dof : maxEdgeDof);
                }
                for(j = 0; j < exp->GetNfaces(); ++j)
                {
                    dof = exp->GetFaceIntNcoeffs(j);
                    maxFaceDof = (dof > maxFaceDof ? dof : maxFaceDof);
                }
                exp->GetInteriorMap(interiorMap);
                dof = interiorMap.size();
                maxIntDof = (dof > maxIntDof ? dof : maxIntDof);
            }

            // Tell other processes about how many dof we have
            vCommRow->AllReduce(nVert, LibUtilities::ReduceSum);
            vCommRow->AllReduce(nEdge, LibUtilities::ReduceSum);
            vCommRow->AllReduce(nFace, LibUtilities::ReduceSum);
            vCommRow->AllReduce(maxEdgeDof, LibUtilities::ReduceMax);
            vCommRow->AllReduce(maxFaceDof, LibUtilities::ReduceMax);
            vCommRow->AllReduce(maxIntDof,  LibUtilities::ReduceMax);

            // Assemble global to universal mapping for this process
            for(i = 0; i < locExpVector.size(); ++i)
            {
                exp = locExpVector[i];
                cnt = locExp.GetCoeff_Offset(i);

                // Loop over all vertices of element i
                for(j = 0; j < exp->GetNverts(); ++j)
                {
                    meshVertId = exp->GetGeom()->GetVid(j);
                    vGlobalId  = m_localToGlobalMap[cnt+exp->GetVertexMap(j)];

                    auto pIt = perVerts.find(meshVertId);
                    if (pIt != perVerts.end())
                    {
                        for (k = 0; k < pIt->second.size(); ++k)
                        {
                            meshVertId = min(meshVertId, pIt->second[k].id);
                        }
                    }

                    m_globalToUniversalMap[vGlobalId] = meshVertId + 1;
                    m_globalToUniversalBndMap[vGlobalId]=m_globalToUniversalMap[vGlobalId];
                    maxBndGlobalId = (vGlobalId > maxBndGlobalId ? vGlobalId : maxBndGlobalId);
                }

                // Loop over all edges of element i
                for(j = 0; j < exp->GetNedges(); ++j)
                {
                    meshEdgeId = exp->GetGeom()->GetEid(j);
                    auto pIt = perEdges.find(meshEdgeId);
                    edgeOrient = exp->GetGeom()->GetEorient(j);

                    if (pIt != perEdges.end())
                    {
                        pair<int, StdRegions::Orientation> idOrient =
                            DeterminePeriodicEdgeOrientId(
                                meshEdgeId, edgeOrient, pIt->second);
                        meshEdgeId = idOrient.first;
                        edgeOrient = idOrient.second;
                    }

                    exp->GetEdgeInteriorMap(j,edgeOrient,edgeInteriorMap,edgeInteriorSign);
                    dof = exp->GetEdgeNcoeffs(j)-2;

                    // Set the global DOF's for the interior modes of edge j
                    //    for varP, ignore modes with sign == 0
                    for(k = 0, l = 0; k < dof; ++k)
                    {
                        if (m_signChange)
                        {
                            if (m_localToGlobalSign[cnt+edgeInteriorMap[k]]==0)
                            {
                                continue;
                            }
                        }
                        vGlobalId = m_localToGlobalMap[cnt+edgeInteriorMap[k]];
                        m_globalToUniversalMap[vGlobalId]
                           = nVert + meshEdgeId * maxEdgeDof + l + 1;
                        m_globalToUniversalBndMap[vGlobalId]=m_globalToUniversalMap[vGlobalId];
                        maxBndGlobalId = (vGlobalId > maxBndGlobalId ? vGlobalId : maxBndGlobalId);
                        l++;
                    }
                }

                // Loop over all faces of element i
                for(j = 0; j < exp->GetNfaces(); ++j)
                {
                    faceOrient = exp->GetGeom()->GetForient(j);

                    meshFaceId = exp->GetGeom()->GetFid(j);

                    auto pIt = perFaces.find(meshFaceId);
                    if (pIt != perFaces.end())
                    {
                        if(meshFaceId == min(meshFaceId, pIt->second[0].id))
                        {
                            faceOrient = DeterminePeriodicFaceOrient(faceOrient,pIt->second[0].orient);
                        }
                        meshFaceId = min(meshFaceId, pIt->second[0].id);
                    }


                    exp->GetFaceInteriorMap(j,faceOrient,faceInteriorMap,faceInteriorSign);
                    dof = exp->GetFaceIntNcoeffs(j);

                    for(k = 0, l = 0; k < dof; ++k)
                    {
                        if (m_signChange)
                        {
                            if (m_localToGlobalSign[cnt+faceInteriorMap[k]]==0)
                            {
                                continue;
                            }
                        }
                        vGlobalId = m_localToGlobalMap[cnt+faceInteriorMap[k]];
                        m_globalToUniversalMap[vGlobalId]
                           = nVert + nEdge*maxEdgeDof + meshFaceId * maxFaceDof
                                   + l + 1;
                        m_globalToUniversalBndMap[vGlobalId]=m_globalToUniversalMap[vGlobalId];

                        maxBndGlobalId = (vGlobalId > maxBndGlobalId ? vGlobalId : maxBndGlobalId);
                        l++;
                    }
                }

                // Add interior DOFs to complete universal numbering
                exp->GetInteriorMap(interiorMap);
                dof = interiorMap.size();
                elementId = (exp->GetGeom())->GetGlobalID();
                for (k = 0; k < dof; ++k)
                {
                    vGlobalId = m_localToGlobalMap[cnt+interiorMap[k]];
                    m_globalToUniversalMap[vGlobalId]
                           = nVert + nEdge*maxEdgeDof + nFace*maxFaceDof + elementId*maxIntDof + k + 1;
                }
            }

            // Set up the GSLib universal assemble mapping
            // Internal DOF do not participate in any data
            // exchange, so we keep these set to the special GSLib id=0 so
            // they are ignored.
            Nektar::Array<OneD, long> tmp(m_numGlobalCoeffs);
            Vmath::Zero(m_numGlobalCoeffs, tmp, 1);
            Nektar::Array<OneD, long> tmp2(m_numGlobalBndCoeffs, tmp);
            for (unsigned int i = 0; i < m_numGlobalBndCoeffs; ++i)
            {
                tmp[i] = m_globalToUniversalMap[i];
            }

            m_gsh = Gs::Init(tmp, vCommRow, verbose);
            m_bndGsh = Gs::Init(tmp2, vCommRow, verbose);
            Gs::Unique(tmp, vCommRow);
            for (unsigned int i = 0; i < m_numGlobalCoeffs; ++i)
            {
                m_globalToUniversalMapUnique[i] = (tmp[i] >= 0 ? 1 : 0);
            }
            for (unsigned int i = 0; i < m_numGlobalBndCoeffs; ++i)
            {
                m_globalToUniversalBndMapUnique[i] = (tmp2[i] >= 0 ? 1 : 0);
            }
        }

        /**
         * @brief Construct an AssemblyMapCG object which corresponds to the
         * linear space of the current object.
         *
         * This function is used to create a linear-space assembly map, which is
         * then used in the linear space preconditioner in the conjugate
         * gradient solve.
         */
        AssemblyMapSharedPtr AssemblyMapCG::v_LinearSpaceMap(
            const ExpList &locexp, GlobalSysSolnType solnType)
        {
            AssemblyMapCGSharedPtr returnval;

            int i, j;
            int nverts = 0;
            const std::shared_ptr<LocalRegions::ExpansionVector> exp
                = locexp.GetExp();
            int nelmts = exp->size();
            const bool verbose = locexp.GetSession()->DefinesCmdLineArgument("verbose");

            // Get Default Map and turn off any searched values.
            returnval = MemoryManager<AssemblyMapCG>
                ::AllocateSharedPtr(m_session);
            returnval->m_solnType           = solnType;
            returnval->m_preconType         = eNull;
            returnval->m_maxStaticCondLevel = 0;
            returnval->m_signChange         = false;
            returnval->m_comm               = m_comm;

            // Count the number of vertices
            for (i = 0; i < nelmts; ++i)
            {
                nverts += (*exp)[i]->GetNverts();
            }

            returnval->m_numLocalCoeffs   = nverts;
            returnval->m_localToGlobalMap = Array<OneD, int>(nverts, -1);

            // Store original global ids in this map
            returnval->m_localToGlobalBndMap = Array<OneD, int>(nverts, -1);

            int cnt  = 0;
            int cnt1 = 0;
            Array<OneD, int> GlobCoeffs(m_numGlobalCoeffs, -1);

            // Set up local to global map;
            for (i = 0; i < nelmts; ++i)
            {
                for (j = 0; j < (*exp)[i]->GetNverts(); ++j)
                {
                    returnval->m_localToGlobalMap[cnt] =
                        returnval->m_localToGlobalBndMap[cnt] =
                        m_localToGlobalMap[cnt1 + (*exp)[i]->GetVertexMap(j,true)];
                    GlobCoeffs[returnval->m_localToGlobalMap[cnt]] = 1;

                    // Set up numLocalDirBndCoeffs
                    if ((returnval->m_localToGlobalMap[cnt]) <
                            m_numGlobalDirBndCoeffs)
                    {
                            returnval->m_numLocalDirBndCoeffs++;
                    }
                    cnt++;
                }
                cnt1 += (*exp)[i]->GetNcoeffs();
            }

            cnt = 0;
            // Reset global numbering and count number of dofs
            for (i = 0; i < m_numGlobalCoeffs; ++i)
            {
                if (GlobCoeffs[i] != -1)
                {
                    GlobCoeffs[i] = cnt++;
                }
            }

            // Set up number of globalCoeffs;
            returnval->m_numGlobalCoeffs = cnt;

            // Set up number of global Dirichlet boundary coefficients
            for (i = 0; i < m_numGlobalDirBndCoeffs; ++i)
            {
                if (GlobCoeffs[i] != -1)
                {
                    returnval->m_numGlobalDirBndCoeffs++;
                }
            }

            // Set up global to universal map
            if (m_globalToUniversalMap.size())
            {
                LibUtilities::CommSharedPtr vCommRow
                    = m_session->GetComm()->GetRowComm();
                int nglocoeffs = returnval->m_numGlobalCoeffs;
                returnval->m_globalToUniversalMap
                    = Array<OneD, int> (nglocoeffs);
                returnval->m_globalToUniversalMapUnique
                    = Array<OneD, int> (nglocoeffs);

                // Reset local to global map and setup universal map
                for (i = 0; i < nverts; ++i)
                {
                    cnt = returnval->m_localToGlobalMap[i];
                    returnval->m_localToGlobalMap[i] = GlobCoeffs[cnt];

                    returnval->m_globalToUniversalMap[GlobCoeffs[cnt]] =
                        m_globalToUniversalMap[cnt];
                }

                Nektar::Array<OneD, long> tmp(nglocoeffs);
                Vmath::Zero(nglocoeffs, tmp, 1);
                for (unsigned int i = 0; i < nglocoeffs; ++i)
                {
                    tmp[i] = returnval->m_globalToUniversalMap[i];
                }
                returnval->m_gsh = Gs::Init(tmp, vCommRow, verbose);
                Gs::Unique(tmp, vCommRow);
                for (unsigned int i = 0; i < nglocoeffs; ++i)
                {
                    returnval->m_globalToUniversalMapUnique[i]
                        = (tmp[i] >= 0 ? 1 : 0);
                }
            }
            else // not sure this option is ever needed.
            {
                for (i = 0; i < nverts; ++i)
                {
                    cnt = returnval->m_localToGlobalMap[i];
                    returnval->m_localToGlobalMap[i] = GlobCoeffs[cnt];
                }
            }

            return returnval;
        }

        /**
         * The bandwidth calculated here corresponds to what is referred to as
         * half-bandwidth.  If the elements of the matrix are designated as
         * a_ij, it corresponds to the maximum value of |i-j| for non-zero
         * a_ij.  As a result, the value also corresponds to the number of
         * sub- or super-diagonals.
         *
         * The bandwith can be calculated elementally as it corresponds to the
         * maximal elemental bandwith (i.e. the maximal difference in global
         * DOF index for every element).
         *
         * We caluclate here the bandwith of the full global system.
         */
        void AssemblyMapCG::CalculateFullSystemBandWidth()
        {
            int i,j;
            int cnt = 0;
            int locSize;
            int maxId;
            int minId;
            int bwidth = -1;
            for(i = 0; i < m_numPatches; ++i)
            {
                locSize = m_numLocalBndCoeffsPerPatch[i]+m_numLocalIntCoeffsPerPatch[i];
                maxId = -1;
                minId = m_numLocalCoeffs+1;
                for(j = 0; j < locSize; j++)
                {
                    if(m_localToGlobalMap[cnt+j] >= m_numGlobalDirBndCoeffs)
                    {
                        if(m_localToGlobalMap[cnt+j] > maxId)
                        {
                            maxId = m_localToGlobalMap[cnt+j];
                        }

                        if(m_localToGlobalMap[cnt+j] < minId)
                        {
                            minId = m_localToGlobalMap[cnt+j];
                        }
                    }
                }
                bwidth = (bwidth>(maxId-minId))?bwidth:(maxId-minId);

                cnt+=locSize;
            }

            m_fullSystemBandWidth = bwidth;
        }


        int AssemblyMapCG::v_GetLocalToGlobalMap(const int i) const
        {
            return m_localToGlobalMap[i];
        }

        int AssemblyMapCG::v_GetGlobalToUniversalMap(const int i) const
        {
            return m_globalToUniversalMap[i];
        }

        int AssemblyMapCG::v_GetGlobalToUniversalMapUnique(const int i) const
        {
            return m_globalToUniversalMapUnique[i];
        }

        const Array<OneD,const int>&
                    AssemblyMapCG::v_GetLocalToGlobalMap(void)
        {
            return m_localToGlobalMap;
        }

        const Array<OneD,const int>&
                    AssemblyMapCG::v_GetGlobalToUniversalMap(void)
        {
            return m_globalToUniversalMap;
        }

        const Array<OneD,const int>&
                    AssemblyMapCG::v_GetGlobalToUniversalMapUnique(void)
        {
            return m_globalToUniversalMapUnique;
        }

        NekDouble AssemblyMapCG::v_GetLocalToGlobalSign(
                    const int i) const
        {
            if(m_signChange)
            {
                return m_localToGlobalSign[i];
            }
            else
            {
                return 1.0;
            }
        }

        const Array<OneD, NekDouble>& AssemblyMapCG::v_GetLocalToGlobalSign() const
        {
            return m_localToGlobalSign;
        }

        void AssemblyMapCG::v_LocalToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                    Array<OneD,       NekDouble>& global,
                    bool useComm) const
        {
            Array<OneD, const NekDouble> local;
            if(global.data() == loc.data())
            {
                local = Array<OneD, NekDouble>(m_numLocalCoeffs,loc.data());
            }
            else
            {
                local = loc; // create reference
            }


            if(m_signChange)
            {
                Vmath::Scatr(m_numLocalCoeffs, m_localToGlobalSign.get(), local.get(), m_localToGlobalMap.get(), global.get());
            }
            else
            {
                Vmath::Scatr(m_numLocalCoeffs, local.get(), m_localToGlobalMap.get(), global.get());
            }

            // ensure all values are unique by calling a max
            if(useComm)
            {
                Gs::Gather(global, Gs::gs_max, m_gsh);
            }
        }


        void AssemblyMapCG::v_LocalToGlobal(
                    const NekVector<NekDouble>& loc,
                    NekVector<      NekDouble>& global,
                    bool useComm) const
        {
            LocalToGlobal(loc.GetPtr(),global.GetPtr(),useComm);
        }

        void AssemblyMapCG::v_GlobalToLocal(
                    const Array<OneD, const NekDouble>& global,
                          Array<OneD,       NekDouble>& loc) const
        {
            Array<OneD, const NekDouble> glo;
            if(global.data() == loc.data())
            {
                glo = Array<OneD, NekDouble>(m_numGlobalCoeffs,global.data());
            }
            else
            {
                glo = global; // create reference
            }


            if(m_signChange)
            {
                Vmath::Gathr(m_numLocalCoeffs, m_localToGlobalSign.get(), glo.get(), m_localToGlobalMap.get(), loc.get());
            }
            else
            {
                Vmath::Gathr(m_numLocalCoeffs, glo.get(), m_localToGlobalMap.get(), loc.get());
            }
        }

        void AssemblyMapCG::v_GlobalToLocal(
                    const NekVector<NekDouble>& global,
                          NekVector<      NekDouble>& loc) const
        {
            GlobalToLocal(global.GetPtr(),loc.GetPtr());
        }

        void AssemblyMapCG::v_Assemble(
                    const Array<OneD, const NekDouble> &loc,
                          Array<OneD,       NekDouble> &global) const
        {
            Array<OneD, const NekDouble> local;
            if(global.data() == loc.data())
            {
                local = Array<OneD, NekDouble>(m_numLocalCoeffs,loc.data());
            }
            else
            {
                local = loc; // create reference
            }

            Vmath::Zero(m_numGlobalCoeffs, global.get(), 1);

            if(m_signChange)
            {
                Vmath::Assmb(m_numLocalCoeffs, m_localToGlobalSign.get(), local.get(), m_localToGlobalMap.get(), global.get());
            }
            else
            {
                Vmath::Assmb(m_numLocalCoeffs, local.get(), m_localToGlobalMap.get(), global.get());
            }
            UniversalAssemble(global);
        }

        void AssemblyMapCG::v_Assemble(
                    const NekVector<NekDouble>& loc,
                          NekVector<      NekDouble>& global) const
        {
            Assemble(loc.GetPtr(),global.GetPtr());
        }

        void AssemblyMapCG::v_UniversalAssemble(
                      Array<OneD,     NekDouble>& pGlobal) const
        {
            Gs::Gather(pGlobal, Gs::gs_add, m_gsh);
        }

        void AssemblyMapCG::v_UniversalAssemble(
                      NekVector<      NekDouble>& pGlobal) const
        {
            UniversalAssemble(pGlobal.GetPtr());
        }

        void AssemblyMapCG::v_UniversalAssemble(
                      Array<OneD,     NekDouble>& pGlobal,
                      int                         offset) const
        {
            Array<OneD, NekDouble> tmp(offset);
            Vmath::Vcopy(offset, pGlobal, 1, tmp, 1);
            UniversalAssemble(pGlobal);
            Vmath::Vcopy(offset, tmp, 1, pGlobal, 1);
        }

        int AssemblyMapCG::v_GetFullSystemBandWidth() const
        {
            return m_fullSystemBandWidth;
        }

        int AssemblyMapCG::v_GetNumNonDirVertexModes() const
        {
            return m_numNonDirVertexModes;
        }

        int AssemblyMapCG::v_GetNumNonDirEdgeModes() const
        {
            return m_numNonDirEdgeModes;
        }

        int AssemblyMapCG::v_GetNumNonDirFaceModes() const
        {
            return m_numNonDirFaceModes;
        }

        int AssemblyMapCG::v_GetNumDirEdges() const
        {
            return m_numDirEdges;
        }

        int AssemblyMapCG::v_GetNumDirFaces() const
        {
            return m_numDirFaces;
        }

        int AssemblyMapCG::v_GetNumNonDirEdges() const
        {
            return m_numNonDirEdges;
        }

        int AssemblyMapCG::v_GetNumNonDirFaces() const
        {
            return m_numNonDirFaces;
        }

        const Array<OneD, const int>& AssemblyMapCG::v_GetExtraDirEdges()
        {
            return m_extraDirEdges;
        }
    } // namespace
} // namespace
