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
#include <LibUtilities/Communication/CommMpi.h>  // TODO: Implement using wrappers so can remove this include
#include <MultiRegions/ExpList.h>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/PointExp.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>

namespace Nektar
{
namespace MultiRegions
{
AssemblyCommDG::AssemblyCommDG()
{
}

AssemblyCommDG::~AssemblyCommDG()
{
}

AssemblyCommDG::AssemblyCommDG(
    const ExpList         &locExp,
    const ExpListSharedPtr &trace,
    const Array<OneD, const MultiRegions::ExpListSharedPtr>         &bndCondExp,
    const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>  &bndCond,
    const PeriodicMap     &perMap)
{
    // Initialise graph structure and link processes across partition boundaries
    AssemblyCommDG::InitialiseStructure(locExp, trace, bndCondExp, bndCond,
                                        perMap);

    // Setting up MPI comm methods
    auto allToAll = AllToAll();
    auto allToAllV = AllToAllV();
    auto neighborAllToAllV = NeighborAllToAllV();
    auto pairwise = Pairwise();

    //Timing MPI comm methods, warm up with 10 iterations then time over 50
    using std::placeholders::_1;
    using std::placeholders::_2;
    std::map<int, func_t> MPIFuncMap;
    MPIFuncMap[0] = std::bind(&AllToAll::PerformExchange, allToAll, _1, _2);
    MPIFuncMap[1] = std::bind(&AllToAllV::PerformExchange, allToAllV, _1, _2);
    MPIFuncMap[2] = std::bind(&NeighborAllToAllV::PerformExchange, neighborAllToAllV, _1, _2);
    MPIFuncMap[3] = std::bind(&Pairwise::PerformExchange, pairwise, _1, _2);

    int numPoints = trace->GetNpoints();
    int warmup = 10, iter = 50;
    NekDouble min, max;
    std::vector<NekDouble> avg(4);

    if (m_comm->GetRank() == 0)
    {
        std::cout << "MPI setup: " << std::endl;
    }

    for (auto const &func : MPIFuncMap)
    {
        MPITiming(m_comm, warmup, numPoints, func.second);
        std::tie(avg[func.first], min, max) = MPITiming(m_comm, iter, numPoints, func.second);
        if (m_comm->GetRank() == 0)
        {
            std::cout << "  " << MPITypeMap[func.first] <<" times (avg, min, max): "
                      << avg[func.first] << " " << min << " " << max << std::endl;
        }
    }

    int fastestMPI = std::distance(avg.begin(), std::min_element(avg.begin(), avg.end()));

    if (m_comm->GetRank() == 0)
    {
        std::cout << "  Chosen fastest method: " << MPITypeMap[fastestMPI] << std::endl;
    }

    m_MPITraceAssemble = MPIFuncMap[fastestMPI];

}

void AssemblyCommDG::InitialiseStructure(
    const ExpList &locExp,
    const ExpListSharedPtr &trace,
    const Nektar::Array<Nektar::OneD, const Nektar::MultiRegions::ExpListSharedPtr> &bndCondExp,
    const Nektar::Array< Nektar::OneD, const Nektar::SpatialDomains::BoundaryConditionShPtr> &bndCond,
    const Nektar::MultiRegions::PeriodicMap &perMap)
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
                    m_elmtToTrace[eid][j]->as<LocalRegions::PointExp>();
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
                        m_elmtToTrace[eid][j]->as<LocalRegions::SegExp>();
                    int id = locSegExp->GetGeom()->GetGlobalID();
                    localEdgeIds.emplace_back(id);
                }
            }
            else if (nDim == 3)
            {
                for (size_t j = 0; j < locExpansion->GetNfaces(); ++j)
                {
                    LocalRegions::Expansion2DSharedPtr locFaceExp =
                        m_elmtToTrace[eid][j]->as<LocalRegions::Expansion2D>();
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
    m_comm->AllReduce(m_maxQuad, LibUtilities::ReduceMax);

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

                uniqueEdgeIdsLocal.emplace_back(eid); // Separate local edge ID list where periodic edge IDs aren't swapped
            }
        }
    }

    // Send uniqueEdgeIds size so all partitions can prepare buffers
    m_nRanks =  m_comm->GetSize();
    Array<OneD, int> rankNumEdges(m_nRanks);
    Array<OneD, int> localEdgeSize(1, uniqueEdgeIds.size());
    m_comm->AllGather(localEdgeSize, rankNumEdges);

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

    Array<OneD, int> rankLocalEdgeIds(std::accumulate(
        rankNumEdges.begin(), rankNumEdges.end(), 0), 0);

    //Send all unique edge IDs to all partitions
    m_comm->AllGatherv(localEdgeIdsArray, rankLocalEdgeIds, rankNumEdges, rankLocalEdgeDisp);

    //Create an array of other rank IDs (all except mine)
    int myRank = m_comm->GetRank();
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
            if (std::find(uniqueEdgeIdsLocal.begin(), uniqueEdgeIdsLocal.end(), edgeId) != uniqueEdgeIdsLocal.end())
            {
                // If periodic then create separate list of minimum of the two ids to be appended on to the end
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

        // Sort periodic edges, keep IDs as is, the m_edgeToTrace is constructed using the min IDs also
        std::sort(periodicEdgeList.begin(), periodicEdgeList.end());
        if(!periodicEdgeList.empty())
        {
            m_rankSharedEdges[rank].insert(m_rankSharedEdges[rank].end(),periodicEdgeList.begin(), periodicEdgeList.end());
        }

        // List of number of quad points in periodic conditions for each rank
        for (auto edgeId : periodicEdgeList)
        {
            perTraceSend[rank] += m_edgeToTrace[edgeId].size();
        }
    }

    // Check that periodic trace edges being communicated are of same order
    Array<OneD, int> perTraceRecv(m_nRanks);
    m_comm->AlltoAll(perTraceSend, perTraceRecv);

    ASSERTL0(perTraceSend == perTraceRecv, "Periodic boundary conditions require the same basis order.")
}

void AssemblyCommDG::PerformExchange(
    const Array<OneD, NekDouble> &inArray,
    Array<OneD, NekDouble> &outArray)
{
    m_MPITraceAssemble(inArray, outArray);
}

AllToAll::AllToAll()
{
    // Get maxCount which is the largest shared partition edge
    for (size_t i = 0; i < m_nRanks; ++i)
    {
        if (m_rankSharedEdges.find(i) != m_rankSharedEdges.end())
        {
            m_maxCount = (m_rankSharedEdges[i].size() > m_maxCount) ? m_rankSharedEdges[i].size() : m_maxCount;
        }
    }

    m_comm->AllReduce(m_maxCount, LibUtilities::ReduceMax);

    // Creates the edge index vector where value -1 indicates
    // padding of value 0 to be inserted instead of value from Fwd
    for (size_t i = 0; i < m_nRanks; ++i)
    {
        if (m_rankSharedEdges.find(i) != m_rankSharedEdges.end())
        {
            for (size_t j = 0; j < m_rankSharedEdges[i].size(); ++j)
            {
                std::vector<int> edgeIndex = m_edgeToTrace[m_rankSharedEdges[i][j]];
                if (edgeIndex.size() < m_maxQuad)
                {
                    std::vector<int> diff(m_maxQuad - edgeIndex.size(), -1);
                    edgeIndex.insert(edgeIndex.end(), diff.begin(), diff.end());
                }

                m_allEdgeIndex.insert(m_allEdgeIndex.end(), edgeIndex.begin(), edgeIndex.end());
            }

            if (m_rankSharedEdges[i].size() < m_maxCount)
            {
                std::vector<int> edgeIndex(m_maxQuad * (m_maxCount - m_rankSharedEdges[i].size()), -1);
                m_allEdgeIndex.insert(m_allEdgeIndex.end(), edgeIndex.begin(), edgeIndex.end());
            }
        }
        else
        {
            std::vector<int> edgeIndex(m_maxQuad * m_maxCount, -1);
            m_allEdgeIndex.insert(m_allEdgeIndex.end(), edgeIndex.begin(), edgeIndex.end());
        }
    }
}

AllToAllV::AllToAllV()
{
    m_allVSendCount = Nektar::Array<OneD, int>(m_nRanks, 0);
    for (size_t i = 0; i < m_nRanks; ++i)
    {
        if (m_rankSharedEdges.find(i) != m_rankSharedEdges.end())
        {
            for (size_t j = 0; j < m_rankSharedEdges[i].size(); ++j)
            {
                std::vector<int> edgeIndex = m_edgeToTrace[m_rankSharedEdges[i][j]];
                m_allVEdgeIndex.insert(m_allVEdgeIndex.end(), edgeIndex.begin(), edgeIndex.end());
                m_allVSendCount[i] += edgeIndex.size();
            }
        }
        else
        {
            m_allVSendCount[i] = 0;
        }
    }

    m_allVSendDisp = Nektar::Array<OneD, int>(m_nRanks, 0);
    for (size_t i = 1; i < m_nRanks; ++i)
    {
        m_allVSendDisp[i] = m_allVSendDisp[i-1] + m_allVSendCount[i-1];
    }
}

NeighborAllToAllV::NeighborAllToAllV()
{
    int nNeighbours = m_rankSharedEdges.size();
    Array<OneD,int> destinations(nNeighbours, 0);
    Array<OneD,int> weights(nNeighbours, 0);
    int cnt = 0;
    for (auto &rankEdgeVec : m_rankSharedEdges)
    {
        destinations[cnt] = rankEdgeVec.first;
        weights[cnt] = rankEdgeVec.second.size();
        ++cnt;
    }

    int retval = MPI_Dist_graph_create_adjacent(MPI_COMM_WORLD,
                                                nNeighbours, destinations.get(), weights.get(),  // Sources
                                                nNeighbours, destinations.get(), weights.get(),  // Destinations
                                                MPI_INFO_NULL, 1, &m_commGraph);

    ASSERTL0(retval == MPI_SUCCESS, "MPI error creating the distributed graph.")

    //Setting up indices
    m_sendCount = Array<OneD, int>(nNeighbours, 0);
    cnt = 0;
    for (const auto &rankEdgeSet : m_rankSharedEdges)
    {
        for (size_t i : rankEdgeSet.second)
        {
            std::vector<int> edgeIndex = m_edgeToTrace[i];
            m_edgeTraceIndex.insert(m_edgeTraceIndex.end(), edgeIndex.begin(), edgeIndex.end());
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

Pairwise::Pairwise()
{
    int cnt = 0;
    Array<OneD,int> sendCount;
    for (const auto& rankEdgeSet : m_rankSharedEdges)
    {
        std::vector<int> edgeTraceIndex;
        for (size_t i : rankEdgeSet.second)
        {
            std::vector<int> edgeIndex = m_edgeToTrace[i];
            edgeTraceIndex.insert(edgeTraceIndex.end(), edgeIndex.begin(), edgeIndex.end());
            sendCount[cnt] += edgeIndex.size();
            m_totSends += m_edgeToTrace[i].size();
        }

        m_vecPairPartitionTrace.emplace_back(std::make_pair(rankEdgeSet.first, edgeTraceIndex));
        ++cnt;
    }

    int nNeighbours = m_rankSharedEdges.size();
    m_sendDisp = Nektar::Array<OneD, int>(nNeighbours, 0);
    for (size_t i = 1; i < nNeighbours; ++i)
    {
        m_sendDisp[i] = m_sendDisp[i-1] + sendCount[i-1];
    }
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


    MPI_Neighbor_alltoallv(sendBuff.get(), m_sendCount.get(), m_sendDisp.get(), MPI_DOUBLE,
                           recvBuff.get(), m_sendCount.get(), m_sendDisp.get(), MPI_DOUBLE,
                           m_commGraph);

    for (size_t i = 0; i < m_edgeTraceIndex.size(); ++i)
    {
        testBwd[m_edgeTraceIndex[i]] = recvBuff[i];
    }
}

void Pairwise::PerformExchange(
    const Array<OneD, NekDouble> &testFwd,
    Array<OneD, NekDouble> &testBwd)
{
    Array<OneD, MPI_Request> request(m_vecPairPartitionTrace.size() * 2);
    Array<OneD, MPI_Status> status(m_vecPairPartitionTrace.size() * 2);
    size_t count = 0, count2 = 0;

    Array<OneD, NekDouble> recvBuff(m_totSends, -1);
    for (auto &pairPartitionTrace : m_vecPairPartitionTrace)
    {
        size_t len = pairPartitionTrace.second.size();

        MPI_Irecv(static_cast<void *>(recvBuff.get() +
                                      m_sendDisp[count++]),
                  len,
                  MPI_DOUBLE,
                  pairPartitionTrace.first,  // rank of source
                  0,                         // message tag
                  MPI_COMM_WORLD,
                  &request[count2++]);
    }

    for (auto &pairPartitionTrace : m_vecPairPartitionTrace)
    {
        size_t len = pairPartitionTrace.second.size();
        Array <OneD, NekDouble> sendBuff(len, -1);

        for (size_t j = 0; j < len; ++j)
        {
            sendBuff[j] = testFwd[pairPartitionTrace.second[j]];
        }

        MPI_Isend(sendBuff.get(),
                  len,
                  MPI_DOUBLE,
                  pairPartitionTrace.first,  // rank of destination
                  0,                         // message tag
                  MPI_COMM_WORLD,
                  &request[count2++]);
    }

    MPI_Waitall(m_vecPairPartitionTrace.size() * 2, request.get(), status.get());

    count = 0;
    for (auto &pairPartitionTrace : m_vecPairPartitionTrace)
    {
        size_t len = pairPartitionTrace.second.size();
        for (size_t j = 0; j < len; ++j)
        {
            testBwd[pairPartitionTrace.second[j]] = recvBuff[m_sendDisp[count] + j];
        }
        count++;
    }
}
} // namespace
} // namespace
