///////////////////////////////////////////////////////////////////////////////
//
// File AssemblyCommDG.cpp
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
// Description: Parallel communication methods for DG with MPI
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/Expansion2D.h>
#include <MultiRegions/AssemblyMap/AssemblyCommDG.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

#include <utility>

namespace Nektar
{
namespace MultiRegions
{

AllToAll::AllToAll(const LibUtilities::CommSharedPtr &comm, const int &maxQuad,
                   const int &nRanks,
                   const std::map<int, std::vector<int>> &rankSharedEdges,
                   const std::map<int, std::vector<int>> &edgeToTrace)
    : m_comm(comm), m_maxQuad(maxQuad), m_nRanks(nRanks)
{
    for (size_t i = 0; i < nRanks; ++i)
    {
        if (rankSharedEdges.find(i) != rankSharedEdges.end())
        {
            m_maxCount = (rankSharedEdges.at(i).size() > m_maxCount)
                             ? static_cast<int>(rankSharedEdges.at(i).size())
                             : m_maxCount;
        }
    }

    comm->AllReduce(m_maxCount, LibUtilities::ReduceMax);

    // Creates the edge index vector where value -1 indicates
    // padding of value 0 to be inserted instead of value from Fwd
    for (size_t i = 0; i < nRanks; ++i)
    {
        if (rankSharedEdges.find(i) != rankSharedEdges.end())
        {
            for (size_t j = 0; j < rankSharedEdges.at(i).size(); ++j)
            {
                std::vector<int> edgeIndex =
                    edgeToTrace.at(rankSharedEdges.at(i)[j]);
                if (edgeIndex.size() < maxQuad)
                {
                    std::vector<int> diff(maxQuad - edgeIndex.size(), -1);
                    edgeIndex.insert(edgeIndex.end(), diff.begin(), diff.end());
                }

                m_allEdgeIndex.insert(m_allEdgeIndex.end(), edgeIndex.begin(),
                                      edgeIndex.end());
            }

            if (rankSharedEdges.at(i).size() < m_maxCount)
            {
                std::vector<int> edgeIndex(
                    maxQuad * (m_maxCount - rankSharedEdges.at(i).size()), -1);
                m_allEdgeIndex.insert(m_allEdgeIndex.end(), edgeIndex.begin(),
                                      edgeIndex.end());
            }
        }
        else
        {
            std::vector<int> edgeIndex(maxQuad * m_maxCount, -1);
            m_allEdgeIndex.insert(m_allEdgeIndex.end(), edgeIndex.begin(),
                                  edgeIndex.end());
        }
    }
}

AllToAllV::AllToAllV(const LibUtilities::CommSharedPtr &comm,
                     const std::map<int, std::vector<int>> &rankSharedEdges,
                     const std::map<int, std::vector<int>> &edgeToTrace,
                     const int &nRanks)
    : m_comm(comm)
{
    m_allVSendCount = Array<OneD, int>(nRanks, 0);
    for (size_t i = 0; i < nRanks; ++i)
    {
        if (rankSharedEdges.find(i) != rankSharedEdges.end())
        {
            for (size_t j = 0; j < rankSharedEdges.at(i).size(); ++j)
            {
                std::vector<int> edgeIndex =
                    edgeToTrace.at(rankSharedEdges.at(i)[j]);
                m_allVEdgeIndex.insert(m_allVEdgeIndex.end(), edgeIndex.begin(),
                                       edgeIndex.end());
                m_allVSendCount[i] += edgeIndex.size();
            }
        }
        else
        {
            m_allVSendCount[i] = 0;
        }
    }

    m_allVSendDisp = Array<OneD, int>(nRanks, 0);
    for (size_t i = 1; i < nRanks; ++i)
    {
        m_allVSendDisp[i] = m_allVSendDisp[i - 1] + m_allVSendCount[i - 1];
    }
}

NeighborAllToAllV::NeighborAllToAllV(
    const LibUtilities::CommSharedPtr &comm,
    const std::map<int, std::vector<int>> &rankSharedEdges,
    const std::map<int, std::vector<int>> &edgeToTrace)
    : m_comm(comm)
{
    int nNeighbours = rankSharedEdges.size();
    Array<OneD, int> destinations(nNeighbours, 0);
    Array<OneD, int> weights(nNeighbours, 0);
    int cnt = 0;
    for (auto &rankEdgeVec : rankSharedEdges)
    {
        destinations[cnt] = rankEdgeVec.first;
        weights[cnt]      = rankEdgeVec.second.size();
        ++cnt;
    }

    comm->DistGraphCreateAdjacent(destinations, weights, 1);

    // Setting up indices
    m_sendCount = Array<OneD, int>(nNeighbours, 0);
    cnt         = 0;

    for (const auto &rankEdgeSet : rankSharedEdges)
    {
        for (size_t i : rankEdgeSet.second)
        {
            std::vector<int> edgeIndex = edgeToTrace.at(i);
            m_edgeTraceIndex.insert(m_edgeTraceIndex.end(), edgeIndex.begin(),
                                    edgeIndex.end());
            m_sendCount[cnt] += edgeIndex.size();
        }

        ++cnt;
    }

    m_sendDisp = Array<OneD, int>(nNeighbours, 0);
    for (size_t i = 1; i < nNeighbours; ++i)
    {
        m_sendDisp[i] = m_sendDisp[i - 1] + m_sendCount[i - 1];
    }
}

Pairwise::Pairwise(const LibUtilities::CommSharedPtr &comm,
                   const std::map<int, std::vector<int>> &rankSharedEdges,
                   const std::map<int, std::vector<int>> &edgeToTrace)
    : m_comm(comm)
{
    int cnt         = 0;
    int nNeighbours = rankSharedEdges.size();
    Array<OneD, int> sendCount(nNeighbours, -1);

    for (const auto &rankEdgeSet : rankSharedEdges)
    {
        std::vector<int> edgeTraceIndex;
        for (size_t i : rankEdgeSet.second)
        {
            std::vector<int> edgeIndex = edgeToTrace.at(i);
            edgeTraceIndex.insert(edgeTraceIndex.end(), edgeIndex.begin(),
                                  edgeIndex.end());
        }

        m_vecPairPartitionTrace.emplace_back(
            std::make_pair(rankEdgeSet.first, edgeTraceIndex));

        sendCount[cnt++] = edgeTraceIndex.size();
    }

    m_sendDisp = Array<OneD, int>(nNeighbours, 0);

    for (size_t i = 1; i < nNeighbours; ++i)
    {
        m_sendDisp[i] = m_sendDisp[i - 1] + sendCount[i - 1];
    }

    size_t totSends = std::accumulate(sendCount.begin(), sendCount.end(), 0);

    m_recvBuff = Array<OneD, NekDouble>(totSends, -1);
    m_sendBuff = Array<OneD, NekDouble>(totSends, -1);

    m_recvRequest = m_comm->CreateRequest(m_vecPairPartitionTrace.size());
    m_sendRequest = m_comm->CreateRequest(m_vecPairPartitionTrace.size());

    // Construct persistent requests
    for (size_t i = 0; i < m_vecPairPartitionTrace.size(); ++i)
    {
        size_t len = m_vecPairPartitionTrace[i].second.size();

        // Initialise receive requests
        m_comm->RecvInit(m_vecPairPartitionTrace[i].first,
                         m_recvBuff[m_sendDisp[i]], len, m_recvRequest, i);

        // Initialise send requests
        m_comm->SendInit(m_vecPairPartitionTrace[i].first,
                         m_sendBuff[m_sendDisp[i]], len, m_sendRequest, i);
    }
}

void AllToAll::PerformExchange(const Array<OneD, NekDouble> &testFwd,
                               Array<OneD, NekDouble> &testBwd)
{
    int size = m_maxQuad * m_maxCount * m_nRanks;
    Array<OneD, NekDouble> sendBuff(size, -1);
    Array<OneD, NekDouble> recvBuff(size, -1);

    for (size_t j = 0; j < size; ++j)
    {
        if (m_allEdgeIndex[j] == -1)
        {
            sendBuff[j] = 0;
        }
        else
        {
            sendBuff[j] = testFwd[m_allEdgeIndex[j]];
        }
    }

    m_comm->AlltoAll(sendBuff, recvBuff);

    for (size_t j = 0; j < size; ++j)
    {
        if (m_allEdgeIndex[j] != -1)
        {
            testBwd[m_allEdgeIndex[j]] = recvBuff[j];
        }
    }
}

void AllToAllV::PerformExchange(const Array<OneD, NekDouble> &testFwd,
                                Array<OneD, NekDouble> &testBwd)
{
    Array<OneD, NekDouble> sendBuff(m_allVEdgeIndex.size(), -1);
    Array<OneD, NekDouble> recvBuff(m_allVEdgeIndex.size(), -1);

    for (size_t i = 0; i < m_allVEdgeIndex.size(); ++i)
    {
        sendBuff[i] = testFwd[m_allVEdgeIndex[i]];
    }

    m_comm->AlltoAllv(sendBuff, m_allVSendCount, m_allVSendDisp, recvBuff,
                      m_allVSendCount, m_allVSendDisp);

    for (size_t i = 0; i < m_allVEdgeIndex.size(); ++i)
    {
        testBwd[m_allVEdgeIndex[i]] = recvBuff[i];
    }
}

void NeighborAllToAllV::PerformExchange(const Array<OneD, NekDouble> &testFwd,
                                        Array<OneD, NekDouble> &testBwd)
{
    Array<OneD, NekDouble> sendBuff(m_edgeTraceIndex.size(), -1);
    Array<OneD, NekDouble> recvBuff(m_edgeTraceIndex.size(), -1);
    for (size_t i = 0; i < m_edgeTraceIndex.size(); ++i)
    {
        sendBuff[i] = testFwd[m_edgeTraceIndex[i]];
    }

    m_comm->NeighborAlltoAllv(sendBuff, m_sendCount, m_sendDisp, recvBuff,
                              m_sendCount, m_sendDisp);

    for (size_t i = 0; i < m_edgeTraceIndex.size(); ++i)
    {
        testBwd[m_edgeTraceIndex[i]] = recvBuff[i];
    }
}

void Pairwise::PerformExchange(const Array<OneD, NekDouble> &testFwd,
                               Array<OneD, NekDouble> &testBwd)
{
    // Perform receive posts
    m_comm->StartAll(m_recvRequest);

    // Fill send buffer from Fwd trace
    for (size_t i = 0; i < m_vecPairPartitionTrace.size(); ++i)
    {
        size_t len = m_vecPairPartitionTrace[i].second.size();
        for (size_t j = 0; j < len; ++j)
        {
            m_sendBuff[m_sendDisp[i] + j] =
                testFwd[m_vecPairPartitionTrace[i].second[j]];
        }
    }

    // Perform send posts
    m_comm->StartAll(m_sendRequest);

    // Wait for all send/recvs to complete
    m_comm->WaitAll(m_sendRequest);
    m_comm->WaitAll(m_recvRequest);

    // Fill Bwd trace from recv buffer
    for (size_t i = 0; i < m_vecPairPartitionTrace.size(); ++i)
    {
        size_t len = m_vecPairPartitionTrace[i].second.size();
        for (size_t j = 0; j < len; ++j)
        {
            testBwd[m_vecPairPartitionTrace[i].second[j]] =
                m_recvBuff[m_sendDisp[i] + j];
        }
    }
}

AssemblyCommDG::AssemblyCommDG(
    const ExpList &locExp, const ExpListSharedPtr &trace,
    const Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr>>
        &elmtToTrace,
    const Array<OneD, const ExpListSharedPtr> &bndCondExp,
    const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndCond,
    const PeriodicMap &perMap)
{
    auto comm = locExp.GetSession()->GetComm();

    // If serial then skip initialising graph structure and the MPI timing
    if (comm->IsSerial())
    {
        m_exchange =
            ExchangeMethodSharedPtr(MemoryManager<Serial>::AllocateSharedPtr());
    }
    else
    {
        // Initialise graph structure and link processes across partition
        // boundaries
        AssemblyCommDG::InitialiseStructure(locExp, trace, elmtToTrace,
                                            bndCondExp, bndCond, perMap, comm);

        // Timing MPI comm methods, warm up with 10 iterations then time over 50
        std::vector<ExchangeMethodSharedPtr> MPIFuncs;
        std::vector<std::string> MPIFuncsNames;

        MPIFuncs.emplace_back(
            ExchangeMethodSharedPtr(MemoryManager<AllToAll>::AllocateSharedPtr(
                comm, m_maxQuad, m_nRanks, m_rankSharedEdges, m_edgeToTrace)));
        MPIFuncsNames.emplace_back("AllToAll");

        MPIFuncs.emplace_back(
            ExchangeMethodSharedPtr(MemoryManager<AllToAllV>::AllocateSharedPtr(
                comm, m_rankSharedEdges, m_edgeToTrace, m_nRanks)));
        MPIFuncsNames.emplace_back("AllToAllV");

        MPIFuncs.emplace_back(
            ExchangeMethodSharedPtr(MemoryManager<Pairwise>::AllocateSharedPtr(
                comm, m_rankSharedEdges, m_edgeToTrace)));
        MPIFuncsNames.emplace_back("PairwiseSendRecv");

        // Disable neighbor MPI method on unsupported MPI version (below 3.0)
        if (std::get<0>(comm->GetVersion()) >= 3)
        {
            MPIFuncs.emplace_back(ExchangeMethodSharedPtr(
                MemoryManager<NeighborAllToAllV>::AllocateSharedPtr(
                    comm, m_rankSharedEdges, m_edgeToTrace)));
            MPIFuncsNames.emplace_back("NeighborAllToAllV");
        }

        int numPoints = trace->GetNpoints();
        int warmup = 10, iter = 50;
        NekDouble min, max;
        std::vector<NekDouble> avg(MPIFuncs.size(), -1);
        bool verbose = locExp.GetSession()->DefinesCmdLineArgument("verbose");

        if (verbose && comm->GetRank() == 0)
        {
            std::cout << "MPI setup for trace exchange: " << std::endl;
        }

        for (size_t i = 0; i < MPIFuncs.size(); ++i)
        {
            Timing(comm, warmup, numPoints, MPIFuncs[i]);
            std::tie(avg[i], min, max) =
                Timing(comm, iter, numPoints, MPIFuncs[i]);
            if (verbose && comm->GetRank() == 0)
            {
                std::cout << "  " << MPIFuncsNames[i]
                          << " times (avg, min, max): " << avg[i] << " " << min
                          << " " << max << std::endl;
            }
        }

        // Gets the fastest MPI method
        int fastestMPI = std::distance(
            avg.begin(), std::min_element(avg.begin(), avg.end()));

        if (verbose && comm->GetRank() == 0)
        {
            std::cout << "  Chosen fastest method: "
                      << MPIFuncsNames[fastestMPI] << std::endl;
        }

        m_exchange = MPIFuncs[fastestMPI];
    }
}

/**
 * This function sets up the initial structure to allow for the exchange methods
 * to be created. This structure is contained within the member variable
 * #m_rankSharedEdges which is a map of rank to vector of the shared edges with
 * that rank. This is filled by:
 *
 * - Create an edge to trace mapping, and realign periodic edges within this
 *   mapping so that they have the same data layout for ranks sharing periodic
 *   boundaries.
 * - Create a list of all local edge IDs and calculate the maximum number of
 *   quadrature points used locally, then perform an AllReduce to find the
 *   maximum number of quadrature points across all ranks (for the AllToAll
 *   method).
 * - Create a list of all boundary edge IDs except for those which are periodic
 * - Using the boundary ID list, and all local ID list we can construct a unique
 *   list of IDs which are on a partition boundary (e.g. if doesn't occur in the
 *   local list twice, and doesn't occur in the boundary list it is on a
 *   partition boundary). We also check, if it is a periodic edge, whether the
 *   other side is local, if not we add the minimum of the two periodic IDs to
 *   the unique list as we must have a consistent numbering scheme across ranks.
 * - We send the unique list to all other ranks/partitions. Each ranks unique
 *   list is then compared with the local unique edge ID list, if a match is
 *   found then the member variable #m_rankSharedEdges is filled with the
 *   matching rank and unique edge ID.
 */
void AssemblyCommDG::InitialiseStructure(
    const ExpList &locExp, const ExpListSharedPtr &trace,
    const Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr>>
        &elmtToTrace,
    const Array<OneD, const ExpListSharedPtr> &bndCondExp,
    const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndCond,
    const PeriodicMap &perMap, const LibUtilities::CommSharedPtr &comm)
{
    Array<OneD, int> tmp;
    int quad = 0, nDim = 0, eid = 0, offset = 0;
    const LocalRegions::ExpansionVector &locExpVector = *(locExp.GetExp());

    // Assume that each element of the expansion is of the same
    // dimension.
    nDim = locExpVector[0]->GetShapeDimension();

    // This sets up the edge to trace mapping and realigns periodic edges
    if (nDim == 1)
    {
        for (size_t i = 0; i < trace->GetExpSize(); ++i)
        {
            eid    = trace->GetExp(i)->GetGeom()->GetGlobalID();
            offset = trace->GetPhys_Offset(i);

            // Check to see if this vert is periodic. If it is, then we
            // need use the unique eid of the two points
            auto it = perMap.find(eid);
            if (perMap.count(eid) > 0)
            {
                PeriodicEntity ent = it->second[0];
                if (!ent.isLocal) // Not sure if true in 1D
                {
                    eid = std::min(eid, ent.id);
                }
            }

            m_edgeToTrace[eid].emplace_back(offset);
        }
    }
    else
    {
        for (size_t i = 0; i < trace->GetExpSize(); ++i)
        {
            eid    = trace->GetExp(i)->GetGeom()->GetGlobalID();
            offset = trace->GetPhys_Offset(i);
            quad   = trace->GetExp(i)->GetTotPoints();

            // Check to see if this edge is periodic. If it is, then we
            // need to reverse the trace order of one edge only in the
            // edge to trace map so that the data are reversed w.r.t each
            // other. We do this by using the minimum of the two IDs.
            auto it      = perMap.find(eid);
            bool realign = false;
            if (perMap.count(eid) > 0)
            {
                PeriodicEntity ent = it->second[0];
                if (!ent.isLocal)
                {
                    realign = eid == std::min(eid, ent.id);
                    eid     = std::min(eid, ent.id);
                }
            }

            for (size_t j = 0; j < quad; ++j)
            {
                m_edgeToTrace[eid].emplace_back(offset + j);
            }

            if (realign)
            {
                // Realign some periodic edges in m_edgeToTrace
                Array<OneD, int> tmpArray(m_edgeToTrace[eid].size());
                for (size_t j = 0; j < m_edgeToTrace[eid].size(); ++j)
                {
                    tmpArray[j] = m_edgeToTrace[eid][j];
                }

                StdRegions::Orientation orient = it->second[0].orient;

                if (nDim == 2)
                {
                    AssemblyMapDG::RealignTraceElement(tmpArray, orient, quad);
                }
                else
                {
                    // Orient is going from face 2 -> face 1 but we want face 1
                    // -> face 2; in all cases except below these are
                    // equivalent. However below is not equivalent so we use the
                    // reverse of the mapping.
                    if (orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1)
                    {
                        orient = StdRegions::eDir1BwdDir2_Dir2FwdDir1;
                    }
                    else if (orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1)
                    {
                        orient = StdRegions::eDir1FwdDir2_Dir2BwdDir1;
                    }

                    AssemblyMapDG::RealignTraceElement(
                        tmpArray, orient, trace->GetExp(i)->GetNumPoints(0),
                        trace->GetExp(i)->GetNumPoints(1));
                }

                for (size_t j = 0; j < m_edgeToTrace[eid].size(); ++j)
                {
                    m_edgeToTrace[eid][j] = tmpArray[j];
                }
            }
        }
    }

    // This creates a list of all geometry of problem dimension - 1
    // and populates the maxQuad member variable
    std::vector<int> localEdgeIds;
    for (eid = 0; eid < locExpVector.size(); ++eid)
    {
        LocalRegions::ExpansionSharedPtr locExpansion = locExpVector[eid];
        nDim = locExpansion->GetShapeDimension();

        if (nDim == 1)
        {
            int nVerts = locExpansion->GetNverts();
            for (size_t j = 0; j < nVerts; ++j)
            {
                LocalRegions::PointExpSharedPtr locPointExp =
                    elmtToTrace[eid][j]->as<LocalRegions::PointExp>();
                int id = locPointExp->GetGeom()->GetGlobalID();
                localEdgeIds.emplace_back(id);
            }

            m_maxQuad = (1 > m_maxQuad ? 1 : m_maxQuad);
        }
        else
        {
            if (nDim == 2)
            {
                for (size_t j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    LocalRegions::SegExpSharedPtr locSegExp =
                        elmtToTrace[eid][j]->as<LocalRegions::SegExp>();
                    int id = locSegExp->GetGeom()->GetGlobalID();
                    localEdgeIds.emplace_back(id);
                }
            }
            else if (nDim == 3)
            {
                for (size_t j = 0; j < locExpansion->GetNfaces(); ++j)
                {
                    LocalRegions::Expansion2DSharedPtr locFaceExp =
                        elmtToTrace[eid][j]->as<LocalRegions::Expansion2D>();
                    int id = locFaceExp->GetGeom()->GetGlobalID();
                    localEdgeIds.emplace_back(id);
                }
            }

            quad = locExpansion->GetTotPoints();
            if (quad > m_maxQuad)
            {
                m_maxQuad = quad;
            }
        }
    }

    // Find max quadrature points across all processes
    comm->AllReduce(m_maxQuad, LibUtilities::ReduceMax);

    // Create list of boundary edge IDs
    std::set<int> bndIdList;
    for (size_t i = 0; i < bndCond.size(); ++i)
    {
        for (size_t j = 0; j < bndCondExp[i]->GetExpSize(); ++j)
        {
            eid = bndCondExp[i]->GetExp(j)->GetGeom()->GetGlobalID();
            if (perMap.find(eid) ==
                perMap.end()) // Don't add if periodic boundary
            {
                bndIdList.insert(eid);
            }
        }
    }

    // Get unique edges to send
    std::vector<int> uniqueEdgeIds;
    std::vector<bool> duplicated(localEdgeIds.size(), false);
    for (size_t i = 0; i < localEdgeIds.size(); ++i)
    {
        eid = localEdgeIds[i];
        for (size_t j = i + 1; j < localEdgeIds.size(); ++j)
        {
            if (eid == localEdgeIds[j])
            {
                duplicated[i] = duplicated[j] = true;
            }
        }

        if (!duplicated[i]) // Not duplicated in local partition
        {
            if (bndIdList.find(eid) == bndIdList.end()) // Not a boundary edge
            {
                // Check if periodic and if not local set eid to other side
                auto it = perMap.find(eid);
                if (it != perMap.end())
                {
                    if (!it->second[0].isLocal)
                    {
                        uniqueEdgeIds.emplace_back(
                            std::min(eid, it->second[0].id));
                    }
                }
                else
                {
                    uniqueEdgeIds.emplace_back(eid);
                }
            }
        }
    }

    // Send uniqueEdgeIds size so all partitions can prepare buffers
    m_nRanks = comm->GetSize();
    Array<OneD, int> rankNumEdges(m_nRanks);
    Array<OneD, int> localEdgeSize(1, uniqueEdgeIds.size());
    comm->AllGather(localEdgeSize, rankNumEdges);

    Array<OneD, int> rankLocalEdgeDisp(m_nRanks, 0);
    for (size_t i = 1; i < m_nRanks; ++i)
    {
        rankLocalEdgeDisp[i] = rankLocalEdgeDisp[i - 1] + rankNumEdges[i - 1];
    }

    Array<OneD, int> localEdgeIdsArray(uniqueEdgeIds.size());
    for (size_t i = 0; i < uniqueEdgeIds.size(); ++i)
    {
        localEdgeIdsArray[i] = uniqueEdgeIds[i];
    }

    // Sort localEdgeIdsArray before sending (this is important!)
    std::sort(localEdgeIdsArray.begin(), localEdgeIdsArray.end());

    Array<OneD, int> rankLocalEdgeIds(
        std::accumulate(rankNumEdges.begin(), rankNumEdges.end(), 0), 0);

    // Send all unique edge IDs to all partitions
    comm->AllGatherv(localEdgeIdsArray, rankLocalEdgeIds, rankNumEdges,
                     rankLocalEdgeDisp);

    // Find what edge Ids match with other ranks
    size_t myRank = comm->GetRank();
    Array<OneD, int> perTraceSend(m_nRanks, 0);
    for (size_t i = 0; i < m_nRanks; ++i)
    {
        if (i == myRank)
        {
            continue;
        }

        for (size_t j = 0; j < rankNumEdges[i]; ++j)
        {
            int edgeId = rankLocalEdgeIds[rankLocalEdgeDisp[i] + j];
            if (std::find(uniqueEdgeIds.begin(), uniqueEdgeIds.end(), edgeId) !=
                uniqueEdgeIds.end())
            {
                m_rankSharedEdges[i].emplace_back(edgeId);
            }
        }
    }
}

/**
 * Timing of the exchange method @p f, performing the exchange @p count times
 * for array of length @p num.
 *
 * @param comm   Communicator
 * @param count  Number of timing iterations to run
 * @param num    Number of quadrature points to communicate
 * @param f      #ExchangeMethod to time
 *
 * @return tuple of loop times {avg, min, max}
 */
std::tuple<NekDouble, NekDouble, NekDouble> AssemblyCommDG::Timing(
    const LibUtilities::CommSharedPtr &comm, const int &count, const int &num,
    const ExchangeMethodSharedPtr& f)
{
    Array<OneD, NekDouble> testFwd(num, 1);
    Array<OneD, NekDouble> testBwd(num, -2);

    LibUtilities::Timer t;
    t.Start();
    for (size_t i = 0; i < count; ++i)
    {
        f->PerformExchange(testFwd, testBwd);
    }
    t.Stop();

    // These can just be 'reduce' but need to setup the wrapper in comm.h
    Array<OneD, NekDouble> minTime(1, t.TimePerTest(count));
    comm->AllReduce(minTime, LibUtilities::ReduceMin);

    Array<OneD, NekDouble> maxTime(1, t.TimePerTest(count));
    comm->AllReduce(maxTime, LibUtilities::ReduceMax);

    Array<OneD, NekDouble> sumTime(1, t.TimePerTest(count));
    comm->AllReduce(sumTime, LibUtilities::ReduceSum);

    NekDouble avgTime = sumTime[0] / comm->GetSize();
    return std::make_tuple(avgTime, minTime[0], maxTime[0]);
}

} // namespace MultiRegions
} // namespace Nektar
