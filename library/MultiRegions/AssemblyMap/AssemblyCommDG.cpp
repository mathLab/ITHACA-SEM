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
// Description: Local to Global DG mapping routines, header file
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/AssemblyMap/AssemblyCommDG.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <LocalRegions/Expansion2D.h>

namespace Nektar
{
namespace MultiRegions
{
AllToAll::AllToAll(
    const LibUtilities::CommSharedPtr &comm,
    const int &maxQuad,
    const int &nRanks,
    const std::map<int, std::vector<int>> &rankSharedEdges,
    const std::map<int, std::vector<int>> &edgeToTrace)
    : m_comm(comm),
    m_maxQuad(maxQuad),
    m_nRanks(nRanks)
{
    for (size_t i = 0; i < nRanks; ++i)
    {
        if (rankSharedEdges.find(i) != rankSharedEdges.end())
        {
            m_maxCount = (rankSharedEdges.at(i).size() > m_maxCount)
                ? rankSharedEdges.at(i).size() : m_maxCount;
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
                std::vector<int> edgeIndex
                = edgeToTrace.at(rankSharedEdges.at(i)[j]);
                if (edgeIndex.size() < maxQuad)
                {
                    std::vector<int> diff(maxQuad - edgeIndex.size(), -1);
                    edgeIndex.insert(edgeIndex.end(), diff.begin(), diff.end());
                }

                m_allEdgeIndex.insert(m_allEdgeIndex.end(),
                    edgeIndex.begin(), edgeIndex.end());
            }

            if (rankSharedEdges.at(i).size() < m_maxCount)
            {
                std::vector<int> edgeIndex(
                    maxQuad * (m_maxCount - rankSharedEdges.at(i).size()), -1);
                m_allEdgeIndex.insert(
                    m_allEdgeIndex.end(), edgeIndex.begin(), edgeIndex.end());
            }
        }
        else
        {
            std::vector<int> edgeIndex(maxQuad * m_maxCount, -1);
            m_allEdgeIndex.insert(
                m_allEdgeIndex.end(), edgeIndex.begin(), edgeIndex.end());
        }
    }
}

AllToAllV::AllToAllV(
    const LibUtilities::CommSharedPtr &comm,
    const std::map<int, std::vector<int>> &rankSharedEdges,
    const std::map<int, std::vector<int>> &edgeToTrace,
    const int &nRanks) : m_comm(comm)
{
    m_allVSendCount = Nektar::Array<OneD, int>(nRanks, 0);
    for (size_t i = 0; i < nRanks; ++i)
    {
        if (rankSharedEdges.find(i) != rankSharedEdges.end())
        {
            for (size_t j = 0; j < rankSharedEdges.at(i).size(); ++j)
            {
                std::vector<int> edgeIndex = edgeToTrace.at(
                    rankSharedEdges.at(i)[j]);
                m_allVEdgeIndex.insert(
                    m_allVEdgeIndex.end(), edgeIndex.begin(), edgeIndex.end());
                m_allVSendCount[i] += edgeIndex.size();
            }
        }
        else
        {
            m_allVSendCount[i] = 0;
        }
    }

    m_allVSendDisp = Nektar::Array<OneD, int>(nRanks, 0);
    for (size_t i = 1; i < nRanks; ++i)
    {
        m_allVSendDisp[i] = m_allVSendDisp[i-1] + m_allVSendCount[i-1];
    }
}

NeighborAllToAllV::NeighborAllToAllV(
    const LibUtilities::CommSharedPtr &comm,
    const std::map<int, std::vector<int>> &rankSharedEdges,
    const std::map<int, std::vector<int>> &edgeToTrace) : m_comm(comm)
{
    int nNeighbours = rankSharedEdges.size();
    Array<OneD,int> destinations(nNeighbours, 0);
    Array<OneD,int> weights(nNeighbours, 0);
    int cnt = 0;
    for (auto &rankEdgeVec : rankSharedEdges)
    {
        destinations[cnt] = rankEdgeVec.first;
        weights[cnt] = rankEdgeVec.second.size();
        ++cnt;
    }

    comm->DistGraphCreateAdjacent(destinations, weights, 1);

    //Setting up indices
    m_sendCount = Array<OneD, int>(nNeighbours, 0);
    cnt = 0;
    for (const auto &rankEdgeSet : rankSharedEdges)
    {
        for (size_t i : rankEdgeSet.second)
        {
            std::vector<int> edgeIndex = edgeToTrace.at(i);
            m_edgeTraceIndex.insert(
                m_edgeTraceIndex.end(), edgeIndex.begin(), edgeIndex.end());
            m_sendCount[cnt] += edgeIndex.size();
        }

        ++cnt;
    }

    m_sendDisp = Nektar::Array<OneD, int>(nNeighbours, 0);
    for (size_t i = 1; i < nNeighbours; ++i)
    {
        m_sendDisp[i] = m_sendDisp[i-1] + m_sendCount[i-1];
    }
}

Pairwise::Pairwise(
    const LibUtilities::CommSharedPtr &comm,
    const std::map<int, std::vector<int>> &rankSharedEdges,
    const std::map<int, std::vector<int>> &edgeToTrace) : m_comm(comm)
{
    int cnt = 0;
    int nNeighbours = rankSharedEdges.size();
    Array<OneD,int> sendCount(nNeighbours, -1);
    for (const auto& rankEdgeSet : rankSharedEdges)
    {
        std::vector<int> edgeTraceIndex;
        for (size_t i : rankEdgeSet.second)
        {
            std::vector<int> edgeIndex = edgeToTrace.at(i);
            edgeTraceIndex.insert(
                edgeTraceIndex.end(), edgeIndex.begin(), edgeIndex.end());
        }

        m_vecPairPartitionTrace.emplace_back(
            std::make_pair(rankEdgeSet.first, edgeTraceIndex));

        sendCount[cnt++] = edgeTraceIndex.size();
    }

    m_totSends = std::accumulate(sendCount.begin(), sendCount.end(), 0);

    m_sendDisp = Nektar::Array<OneD, int>(nNeighbours, 0);
    for (size_t i = 1; i < nNeighbours; ++i)
    {
        m_sendDisp[i] = m_sendDisp[i-1] + sendCount[i-1];
    }

    m_requests = m_comm->CreateRequest(m_vecPairPartitionTrace.size() * 2);
}

void AllToAll::PerformExchange(
    const Array<OneD, NekDouble> &testFwd,
    Array<OneD, NekDouble> &testBwd)
{
    int size = m_maxQuad * m_maxCount * m_nRanks;
    Array<OneD, double> sendBuff(size, -1);
    Array<OneD, double> recvBuff(size, -1);

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

void AllToAllV::PerformExchange(
    const Array<OneD, NekDouble> &testFwd,
    Array<OneD, NekDouble> &testBwd)
{
    Array<OneD, double> sendBuff(m_allVEdgeIndex.size(), -1);
    Array<OneD, double> recvBuff(m_allVEdgeIndex.size(), -1);

    for (size_t i = 0; i < m_allVEdgeIndex.size(); ++i)
    {
        sendBuff[i] = testFwd[m_allVEdgeIndex[i]];
    }

    m_comm->AlltoAllv(sendBuff, m_allVSendCount, m_allVSendDisp,
                      recvBuff, m_allVSendCount, m_allVSendDisp);

    for (size_t i = 0; i < m_allVEdgeIndex.size(); ++i)
    {
        testBwd[m_allVEdgeIndex[i]] = recvBuff[i];
    }
}

void NeighborAllToAllV::PerformExchange(
    const Array<OneD, NekDouble> &testFwd,
    Array<OneD, NekDouble> &testBwd)
{
    Array<OneD, double> sendBuff(m_edgeTraceIndex.size(), -1);
    Array<OneD, double> recvBuff(m_edgeTraceIndex.size(), -1);
    for (size_t i = 0; i < m_edgeTraceIndex.size(); ++i)
    {
        sendBuff[i] = testFwd[m_edgeTraceIndex[i]];
    }

    m_comm->NeighborAlltoAllv(
        sendBuff, m_sendCount, m_sendDisp, recvBuff, m_sendCount, m_sendDisp);

    for (size_t i = 0; i < m_edgeTraceIndex.size(); ++i)
    {
        testBwd[m_edgeTraceIndex[i]] = recvBuff[i];
    }
}

void Pairwise::PerformExchange(
    const Array<OneD, NekDouble> &testFwd,
    Array<OneD, NekDouble> &testBwd)
{
    Array<OneD, NekDouble> recvBuff(m_totSends, -1);
    // Perform receive operations
    for (size_t i = 0; i < m_vecPairPartitionTrace.size(); ++i)
    {
        size_t len = m_vecPairPartitionTrace[i].second.size();
        m_comm->Irecv(m_vecPairPartitionTrace[i].first,
                      recvBuff[m_sendDisp[i]], len, m_requests, 2*i);
    }

    // Can't reuse buffer until previous send has completed so we store all
    Array<OneD, NekDouble> sendBuff(m_totSends, -1);
    for (size_t i = 0; i < m_vecPairPartitionTrace.size(); ++i)
    {
        size_t len = m_vecPairPartitionTrace[i].second.size();
        for (size_t j = 0; j < len; ++j)
        {
            sendBuff[m_sendDisp[i] + j] =
                testFwd[m_vecPairPartitionTrace[i].second[j]];
        }

        // Perform send operations
        m_comm->Irsend(m_vecPairPartitionTrace[i].first,
                       sendBuff[m_sendDisp[i]], len, m_requests, 2*i + 1);
    }

    //Wait for all send/recvs to complete
    m_comm->WaitAll(m_requests);

    for (size_t i = 0; i < m_vecPairPartitionTrace.size(); ++i)
    {
        size_t len = m_vecPairPartitionTrace[i].second.size();
        for (size_t j = 0; j < len; ++j)
        {
            testBwd[m_vecPairPartitionTrace[i].second[j]] =
                recvBuff[m_sendDisp[i] + j];
        }
    }
}

AssemblyCommDG::AssemblyCommDG(
    const ExpList         &locExp,
    const ExpListSharedPtr &trace,
    const Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr>> elmtToTrace,
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &bndCondExp,
    const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndCond,
    const PeriodicMap &perMap)
{
    auto comm = locExp.GetSession()->GetComm();

    // Initialise graph structure and link processes across partition boundaries
    AssemblyCommDG::InitialiseStructure(
        locExp, trace, elmtToTrace, bndCondExp, bndCond, perMap, comm);

    //Timing MPI comm methods, warm up with 10 iterations then time over 50
    std::map<int, ExchangeMethodSharedPtr> MPIFuncMap;

    MPIFuncMap[0] = ExchangeMethodSharedPtr(
        MemoryManager<AllToAll>::AllocateSharedPtr(
            comm, m_maxQuad, m_nRanks, m_rankSharedEdges, m_edgeToTrace));

    MPIFuncMap[1] = ExchangeMethodSharedPtr(
        MemoryManager<AllToAllV>::AllocateSharedPtr(
            comm, m_rankSharedEdges, m_edgeToTrace, m_nRanks));

    MPIFuncMap[2] = ExchangeMethodSharedPtr(
        MemoryManager<NeighborAllToAllV>::AllocateSharedPtr(
            comm, m_rankSharedEdges, m_edgeToTrace));

    MPIFuncMap[3] = ExchangeMethodSharedPtr(
        MemoryManager<Pairwise>::AllocateSharedPtr(
            comm, m_rankSharedEdges, m_edgeToTrace));

    const char* MPITypeMap[] =
    {
        "AllToAll",
        "AllToAllV",
        "NeighborAllToAllV",
        "PairwiseSendRecv"
    };

    int numPoints = trace->GetNpoints();
    int warmup = 10, iter = 50;
    NekDouble min, max;
    std::vector<NekDouble> avg(4, -1);

    if (comm->GetRank() == 0)
    {
        std::cout << "MPI setup: " << std::endl;
    }

    for (auto &func : MPIFuncMap)
    {
        Timing(comm, warmup, numPoints, func.second);
        std::tie(avg[func.first], min, max) =
            Timing(comm, iter, numPoints, func.second);
        if (comm->GetRank() == 0)
        {
            std::cout << "  " << MPITypeMap[func.first]
                      << " times (avg, min, max): "
                      << avg[func.first] << " " << min
                      << " " << max << std::endl;
        }
    }

    // Gets the fastest MPI method greater than 0
    int fastestMPI = std::distance(
        avg.begin(), std::min_element(avg.begin(), avg.end(),
            [](NekDouble a, NekDouble b)
            { return (a < 0) ? false : (b < 0) ? true : (a < b); }));

    if (comm->GetRank() == 0)
    {
        std::cout << "  Chosen fastest method: " << MPITypeMap[fastestMPI]
                  << std::endl;
    }

    m_exchange = MPIFuncMap[fastestMPI];
}

void AssemblyCommDG::InitialiseStructure(
    const Nektar::MultiRegions::ExpList &locExp,
    const Nektar::MultiRegions::ExpListSharedPtr &trace,
    const Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr>> elmtToTrace,
    const Nektar::Array<Nektar::OneD, const Nektar::MultiRegions::ExpListSharedPtr> &bndCondExp,
    const Nektar::Array<Nektar::OneD, const Nektar::SpatialDomains::BoundaryConditionShPtr> &bndCond,
    const Nektar::MultiRegions::PeriodicMap &perMap,
    const Nektar::LibUtilities::CommSharedPtr &comm)
{
    Array<OneD, int> tmp;
    int quad = 0, nDim = 0, eid = 0, offset = 0;
    const LocalRegions::ExpansionVector &locExpVector = *(locExp.GetExp());

    // Assume that each element of the expansion is of the same
    // dimension.
    nDim = locExpVector[0]->GetShapeDimension();

    // @todo combine the two GetExpSize() loops into one.
    // This sets up the edge to trace mapping and realigns periodic edges
    if (nDim == 1)
    {
        for (size_t i = 0; i < trace->GetExpSize(); ++i)
        {
            eid = trace->GetExp(i)->GetGeom()->GetGlobalID();
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
            auto it = perMap.find(eid);
            bool realign = false;
            if (perMap.count(eid) > 0)
            {
                PeriodicEntity ent = it->second[0];
                if (!ent.isLocal)
                {
                    realign = eid == std::min(eid, ent.id);
                    eid = std::min(eid, ent.id);
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

                if (nDim == 2)
                {
                    AssemblyMapDG::RealignTraceElement(
                        tmpArray,
                        it->second[0].orient, quad);
                }
                else
                {
                    AssemblyMapDG::RealignTraceElement(
                        tmpArray,
                        it->second[0].orient,
                        trace->GetExp(i)->GetNumPoints(0),
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

    //Create list of boundary edge IDs
    std::set<int> bndIdList;
    for (size_t i = 0; i < bndCond.num_elements(); ++i)
    {
        for (size_t j = 0; j < bndCondExp[i]->GetExpSize(); ++j)
        {
            eid = bndCondExp[i]->GetExp(j)->GetGeom()->GetGlobalID();
            if (perMap.find(eid) == perMap.end())  // Don't add if periodic boundary
            {
                bndIdList.insert(eid);
            }
        }
    }

    //Get unique edges to send
    std::vector<int> uniqueEdgeIds, uniqueEdgeIdsLocal;
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

        if (!duplicated[i])  // Not duplicated in local partition
        {
            if (bndIdList.find(eid) == bndIdList.end())  // Not a boundary edge
            {
                // Check if periodic and if not local set eid to other side
                auto it = perMap.find(eid);
                if (it != perMap.end())
                {
                    if (!it->second[0].isLocal)
                    {
                        uniqueEdgeIds.emplace_back(it->second[0].id);
                    }
                }
                else
                {
                    uniqueEdgeIds.emplace_back(eid);
                }

                // Separate local edge ID list where periodic edge IDs
                // aren't swapped
                uniqueEdgeIdsLocal.emplace_back(eid);
            }
        }
    }

    // Send uniqueEdgeIds size so all partitions can prepare buffers
    m_nRanks =  comm->GetSize();
    Array<OneD, int> rankNumEdges(m_nRanks);
    Array<OneD, int> localEdgeSize(1, uniqueEdgeIds.size());
    comm->AllGather(localEdgeSize, rankNumEdges);

    Array<OneD, int> rankLocalEdgeDisp(m_nRanks, 0);
    for (size_t i = 1; i < m_nRanks; ++i)
    {
        rankLocalEdgeDisp[i] = rankLocalEdgeDisp[i-1] + rankNumEdges[i-1];
    }

    Array<OneD, int> localEdgeIdsArray(uniqueEdgeIds.size());
    for (size_t i = 0; i < uniqueEdgeIds.size(); ++i)
    {
        localEdgeIdsArray[i] = uniqueEdgeIds[i];
    }

    //Sort localEdgeIdsArray before sending (this is important!)
    std::sort(localEdgeIdsArray.begin(), localEdgeIdsArray.end());

    // If only one rank set array size to 1 to avoid undefined behaviour
    int size = (m_nRanks != 1) ?
        std::accumulate(rankNumEdges.begin(), rankNumEdges.end(), 0) : 1;
    Array<OneD, int> rankLocalEdgeIds(size, 0);

    //Send all unique edge IDs to all partitions
    comm->AllGatherv(localEdgeIdsArray, rankLocalEdgeIds,
                     rankNumEdges, rankLocalEdgeDisp);

    //Create an array of other rank IDs (all except mine)
    int myRank = comm->GetRank();
    int cnt = 0;
    Array<OneD, int> otherRanks(m_nRanks - 1);
    for(int i = 0; i < m_nRanks; ++i)
    {
        if (i != myRank)
        {
            otherRanks[cnt] = i;
            ++cnt;
        }
    }

    // Find what edge Ids match with other ranks
    Array<OneD, int> perTraceSend(m_nRanks, 0);
    for (auto &rank : otherRanks)
    {
        std::vector<int> periodicEdgeList;
        for (size_t j = 0; j < rankNumEdges[rank]; ++j)
        {
            int edgeId = rankLocalEdgeIds[rankLocalEdgeDisp[rank] + j];
            if (std::find(
                uniqueEdgeIdsLocal.begin(), uniqueEdgeIdsLocal.end(), edgeId)
                != uniqueEdgeIdsLocal.end())
            {
                // If periodic then create separate list of minimum
                // of the two ids to be appended on to the end
                auto it = perMap.find(edgeId);
                if (it != perMap.end())
                {
                    int locVal = std::min(edgeId, it->second[0].id);
                    periodicEdgeList.emplace_back(locVal);
                }
                else
                {
                    m_rankSharedEdges[rank].emplace_back(edgeId);
                }
            }
        }

        // Sort periodic edges, keep IDs as is becausethe edgeToTrace
        // is constructed using the min IDs also
        std::sort(periodicEdgeList.begin(), periodicEdgeList.end());
        if(!periodicEdgeList.empty())
        {
            m_rankSharedEdges[rank].insert(m_rankSharedEdges[rank].end(),
                                           periodicEdgeList.begin(),
                                           periodicEdgeList.end());
        }

        // List of number of quad points in periodic conditions for each rank
        for (auto edgeId : periodicEdgeList)
        {
            perTraceSend[rank] += m_edgeToTrace[edgeId].size();
        }
    }

    // Check that periodic trace edges being communicated are of same order
    Array<OneD, int> perTraceRecv(m_nRanks, 0);
    comm->AlltoAll(perTraceSend, perTraceRecv);

    ASSERTL0(perTraceSend == perTraceRecv,
        "Periodic boundary conditions require the same basis order.")
}

std::tuple<NekDouble, NekDouble, NekDouble> AssemblyCommDG::Timing(
    const LibUtilities::CommSharedPtr &comm,
    const int &count,
    const int &num,
    ExchangeMethodSharedPtr f)
{
    Array<OneD, double> testFwd(num, 1);
    Array<OneD, double> testBwd(num, -2);

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

    return {sumTime[0]/comm->GetSize(), minTime[0], maxTime[0]};
}

} // namespace
} // namespace
