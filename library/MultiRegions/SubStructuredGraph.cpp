///////////////////////////////////////////////////////////////////////////////
//
// File SubStructuredGraph.cpp
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
// Description: a collection of classes that facilitates the implementation
//              of the multi-level static condensation routines
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/SubStructuredGraph.h>
#include <LibUtilities/BasicUtils/VmathArray.hpp>

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/algorithm/string/replace.hpp>

using std::max;
using std::cout;
using std::endl;

#ifdef NEKTAR_USE_SCOTCH
#ifdef NEKTAR_USE_MPI
#include <ptscotch.h>
#else
#include <scotch.h>
#endif

#define SCOTCH_CALL(scotchFunc, args)                                   \
    {                                                                   \
        ASSERTL0(scotchFunc args == 0,                                  \
                 std::string("Error in Scotch calling function ")       \
                 + std::string(#scotchFunc));                           \
    }
#endif

namespace Nektar
{
    namespace MultiRegions
    {
        PatchMap::PatchMap(void)
        {

        }

        PatchMap::PatchMap(const int nvals)
        {
            m_patchId  = Array<OneD, int>         (nvals);
            m_dofId    = Array<OneD, int>         (nvals);
            m_bndPatch = Array<OneD, unsigned int>(nvals);
            m_sign     = Array<OneD, NekDouble>   (nvals);
        }

        PatchMap::~PatchMap(void)
        {

        }

        void PatchMap::SetPatchMap(
            const int          n,
            const int          patchId,
            const int          dofId,
            const unsigned int bndPatch,
            const NekDouble    sign)
        {
            m_patchId [n] = patchId;
            m_dofId   [n] = dofId;
            m_bndPatch[n] = bndPatch;
            m_sign    [n] = sign;
        }
        
        void PatchMap::SetNewLevelMap(Array<OneD, const unsigned int> numLocalBndCondPerPatch,
                                      Array<OneD, const unsigned int> numLocalIntCondPerPatch)
        {
            int npatch = numLocalBndCondPerPatch.size();
            
            Array<OneD, int> bndoffset(npatch+1);
            Array<OneD, int> intoffset(npatch+1);

            bndoffset[0] = intoffset[0] = 0;
            for(int i = 1; i <= npatch; ++i)
            {
                bndoffset[i] = bndoffset[i-1] + numLocalBndCondPerPatch[i-1];
                intoffset[i] = intoffset[i-1] + numLocalIntCondPerPatch[i-1];
            }
            
            m_newLevelMap = Array<OneD, int>(m_dofId.size());

            for(int i = 0; i < m_dofId.size(); ++i)
            {
                if(m_bndPatch[i] == true)
                {
                    m_newLevelMap[i] = bndoffset[m_patchId[i]] + m_dofId[i];
                }
                else
                {
                    m_newLevelMap[i] = bndoffset[npatch] + intoffset[m_patchId[i]] + m_dofId[i];
                }
            }
        }
        
        MultiLevelBisectedGraph::MultiLevelBisectedGraph(
            MultiLevelBisectedGraphSharedPtr oldLevel,
            const int                        nPartition)
        {
            m_daughterGraphs.push_back(oldLevel);
            m_BndDofs = MemoryManager<SubGraph>::AllocateSharedPtr(nPartition);
        }

        MultiLevelBisectedGraph::MultiLevelBisectedGraph(const int nBndDofs):
            m_BndDofs(MemoryManager<SubGraph>::AllocateSharedPtr(nBndDofs))
        {
        }

        MultiLevelBisectedGraph::~MultiLevelBisectedGraph(void)
        {
        }

        int MultiLevelBisectedGraph::GetTotDofs() const
        {
            int returnval = 0;

            for (auto &g : m_daughterGraphs)
            {
                returnval += g->GetTotDofs();
            }

            returnval += m_BndDofs->GetNverts();
            return returnval;
        }

        void MultiLevelBisectedGraph::SetGlobalNumberingOffset()
        {
            static int level = 0;
            static int offset = 0;
            level++;

            for (auto &g : m_daughterGraphs)
            {
                g->SetGlobalNumberingOffset();
            }

            m_BndDofs->SetIdOffset(offset);
            offset += m_BndDofs->GetNverts();

            level--;
            if(level == 0)
            {
                offset = 0;
            }
        }

        void MultiLevelBisectedGraph::DumpNBndDofs(void) const
        {
            static int level = 0;
            level++;
            cout << "LEVEL " << level << " " << m_BndDofs->GetNverts() << endl;

            for (auto &g : m_daughterGraphs)
            {
                g->DumpNBndDofs();
            }

            level--;
        }

        void MultiLevelBisectedGraph::CollectLeaves(
            std::vector<SubGraphSharedPtr>& leaves) const
        {
            int cnt = 0;

            for (auto &g : m_daughterGraphs)
            {
                g->CollectLeaves(leaves);
                cnt++;
            }

            if (cnt == 0)
            {
                SubGraphSharedPtr leave = m_BndDofs;
                leaves.push_back(leave);
            }
        }

        int MultiLevelBisectedGraph::CutEmptyLeaves()
        {
            int        returnval;
            static int level = 0;
            static int nLeaves = 0;
            level++;

            if (level == 1 && m_daughterGraphs.size() == 0)
            {
                level   = 0;
                nLeaves = 0;
                return 0;
            }

            for (auto it = m_daughterGraphs.begin(); it != m_daughterGraphs.end();)
            {
                auto g = *it;
                if (g->GetNdaughterGraphs() == 0 &&
                    g->GetBndDofsGraph()->GetNverts() == 0)
                {
                    it = m_daughterGraphs.erase(it);
                    nLeaves++;
                }
                else
                {
                    g->CutEmptyLeaves();
                    ++it;
                }
            }

            returnval = nLeaves;

            level--;
            if(level == 0)
            {
                nLeaves = 0;
            }

            return returnval;
        }

        int MultiLevelBisectedGraph::CutLeaves()
        {
            int returnval;
            static int level = 0;
            static int nLeaves = 0;
            level++;

            if (level == 1 && m_daughterGraphs.size() == 0)
            {
                level = 0;
                nLeaves = 0;
                return 0;
            }

            for (auto it = m_daughterGraphs.begin(); it != m_daughterGraphs.end();)
            {
                auto g = *it;
                if (g->GetNdaughterGraphs() == 0)
                {
                    it = m_daughterGraphs.erase(it);
                    nLeaves++;
                }
                else
                {
                    g->CutLeaves();
                    ++it;
                }
            }

            returnval = nLeaves;

            level--;
            if(level == 0)
            {
                nLeaves = 0;
            }

            return returnval;
        }


        BottomUpSubStructuredGraph::BottomUpSubStructuredGraph(
            const int nVerts) : m_IntBlocks(), m_daughterGraph()
        {
            // This constructor should only be used in the very special case
            // when there is a graph consisting of nVerts vertices but without
            // any connectivity whatsowever (i.e. no edges)

            SubGraphSharedPtr subgraph;
            for (int i = 0 ; i < nVerts; i++)
            {
                subgraph = MemoryManager<SubGraph>::AllocateSharedPtr(1,i);
                m_IntBlocks.push_back(subgraph);
            }
        }

        BottomUpSubStructuredGraph::BottomUpSubStructuredGraph(
            MultiLevelBisectedGraphSharedPtr graph,
            int nPartition,
            bool globaloffset) :
            m_IntBlocks    (),
            m_daughterGraph()
        {
            int ncuts;

            if (nPartition > 0)
            {
                graph = MemoryManager<MultiLevelBisectedGraph>::
                    AllocateSharedPtr(graph, nPartition);
            }

            if (globaloffset)
            {
                graph->SetGlobalNumberingOffset();
            }

            graph->CutEmptyLeaves();
            graph->CollectLeaves(m_IntBlocks);
            ncuts = graph->CutLeaves();

            if(ncuts)
            {
                m_daughterGraph = MemoryManager<BottomUpSubStructuredGraph>::
                    AllocateSharedPtr(graph);
            }
        }

        BottomUpSubStructuredGraph::~BottomUpSubStructuredGraph(void)
        {
        }

        int BottomUpSubStructuredGraph::GetTotDofs() const
        {
            static int nIntDofs = 0;
            static int level = 0;
            level++;

            int returnval;

            if( m_daughterGraph.get())
            {
                m_daughterGraph->GetTotDofs();
            }

            for(int i = 0; i < m_IntBlocks.size(); i++)
            {
                nIntDofs += m_IntBlocks[i]->GetNverts();
            }

            returnval = nIntDofs;

            level--;
            if(level == 0)
            {
                nIntDofs = 0;
            }

            return returnval;
        }

        void BottomUpSubStructuredGraph::Dump() const
        {
            static int level = 0;
            level++;

            cout << "LEVEL " << level << endl;
            cout << "interior blocks" << endl;
            for (int i = 0; i < m_IntBlocks.size(); i++)
            {
                cout << "  " << i
                     << "/"  << m_IntBlocks.size()-1
                     << ": " << m_IntBlocks[i]->GetNverts()
                     << ", " << m_IntBlocks[i]->GetIdOffset() << endl;
            }
            cout << endl;

            if (m_daughterGraph.get())
            {
                m_daughterGraph->Dump();
            }

            level--;
        }

        void BottomUpSubStructuredGraph::UpdateBottomUpReordering(
            Array<OneD, int> &perm,
            Array<OneD, int> &iperm) const
        {
            int nDofs = GetTotDofs();

            // Step 1: make a permutation array that goes from the current
            // reordering in the bottom-up graph to an ordering where the
            // interior dofs of the first (=bottom) level are ordered last,
            // followed interior dofs of the second level, ...
            Array<OneD, int> iperm1(nDofs);
            SetBottomUpReordering(iperm1);

            // Now, based upon the input permutation array, update the
            // permutation arrays between the original ordering of the dofs
            // (defined somewhere outside) and the final reordering as given by
            // iperm1
            for (int i=0; i < nDofs; i++)
            {
                iperm[i]       = iperm1[iperm[i]];
                perm[iperm[i]] = i;
            }
        }

        void BottomUpSubStructuredGraph::SetBottomUpReordering(
            Array<OneD, int>& iperm) const
        {
            static int offset = 0;
            static int level = 0;
            level++;

            if( m_daughterGraph.get() )
            {
                m_daughterGraph->SetBottomUpReordering(iperm);
            }

            for(int i = 0; i < m_IntBlocks.size(); i++)
            {
                int GlobIdOffset = m_IntBlocks[i]->GetIdOffset();
                m_IntBlocks[i]->SetIdOffset(offset);
                for(int j = 0; j < m_IntBlocks[i]->GetNverts(); j++)
                {
                    iperm[GlobIdOffset+j] = offset;
                    offset++;
                }
            }

            level--;
            if(level == 0)
            {
                offset = 0;
            }
        }

        void BottomUpSubStructuredGraph::ExpandGraphWithVertexWeights(
            const Array<OneD, const int>& wgts)
        {
            static int offset = 0;
            static int level = 0;
            level++;

            if( m_daughterGraph.get())
            {
                m_daughterGraph->ExpandGraphWithVertexWeights(wgts);
            }

            for(int i = 0; i < m_IntBlocks.size(); i++)
            {
                int OrigGlobIdOffset = m_IntBlocks[i]->GetIdOffset();
                int newNverts = 0;

                m_IntBlocks[i]->SetIdOffset(offset);

                for(int j = 0; j < m_IntBlocks[i]->GetNverts(); j++)
                {
                    newNverts += wgts[OrigGlobIdOffset+j];
                    offset    += wgts[OrigGlobIdOffset+j];
                }

                m_IntBlocks[i]->SetNverts(newNverts);
            }

            // erase the blocks that do not have any vertices
            m_IntBlocks.erase(
                remove_if(m_IntBlocks.begin(),
                          m_IntBlocks.end(),
                          SubGraphWithoutVerts),
                m_IntBlocks.end());

            // remove the current level if there are no interior blocks
            if(m_IntBlocks.size() == 0 && m_daughterGraph.get())
            {
                m_IntBlocks     = m_daughterGraph->GetInteriorBlocks();
                m_daughterGraph = m_daughterGraph->GetDaughterGraph();
            }

            level--;
            if(level == 0)
            {
                offset = 0;
            }
        }

        void BottomUpSubStructuredGraph::MaskPatches(
            const int               leveltomask,
            Array<OneD, NekDouble> &maskarray) const
        {
            static int level = 0;
            level++;

            if(level < leveltomask)
            {
                m_daughterGraph->MaskPatches(leveltomask,maskarray);
            }
            else if(level == leveltomask)
            {
                int GlobIdOffset;
                int nVerts;
                for(int i = 0; i < m_IntBlocks.size(); i++)
                {
                    GlobIdOffset = m_IntBlocks[i]->GetIdOffset();
                    nVerts       = m_IntBlocks[i]->GetNverts();
                    for(int j = 0; j < nVerts; j++)
                    {
                        maskarray[GlobIdOffset+j] = (NekDouble) i;
                    }
                }
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                         "If statement should not arrive here");
            }

            level--;
        }

        int BottomUpSubStructuredGraph::GetNpatchesWithInterior(
            const int whichlevel) const
        {
            int returnval = -1;
            static int level = 0;
            level++;

            if(level < whichlevel)
            {
                returnval = m_daughterGraph->GetNpatchesWithInterior(
                    whichlevel);
            }
            else if(level == whichlevel)
            {
                returnval = m_IntBlocks.size();
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                         "If statement should not arrive here");
            }

            level--;

            return returnval;
        }

        void BottomUpSubStructuredGraph::GetNintDofsPerPatch(
            const int                  whichlevel,
            Array<OneD, unsigned int> &outarray) const
        {
            static int level = 0;
            level++;

            if(level < whichlevel)
            {
                m_daughterGraph->GetNintDofsPerPatch(whichlevel,outarray);
            }
            else if(level == whichlevel)
            {
                ASSERTL1(outarray.size() >= m_IntBlocks.size(),
                         "Array dimension not sufficient");

                for(int i = 0; i < m_IntBlocks.size(); i++)
                {
                    outarray[i] = (unsigned int) m_IntBlocks[i]->GetNverts();
                }
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                         "If statement should not arrive here");
            }

            level--;
        }

        int BottomUpSubStructuredGraph::GetInteriorOffset(
            const int whichlevel, const int patch) const
        {
            int retval = -1;
            static int level = 0;
            level++;

            if(level < whichlevel)
            {
                retval = m_daughterGraph->GetInteriorOffset(whichlevel,patch);
            }
            else if(level == whichlevel)
            {
                retval = m_IntBlocks[patch]->GetIdOffset();
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                         "If statement should not arrive here");
            }

            level--;

            return retval;
        }

        std::vector<SubGraphSharedPtr> BottomUpSubStructuredGraph::
            GetInteriorBlocks(const int whichlevel) const
        {
            std::vector<SubGraphSharedPtr> returnval;
            static int level = 0;
            level++;

            if(level < whichlevel)
            {
                returnval = m_daughterGraph->GetInteriorBlocks(whichlevel);
            }
            else if(level == whichlevel)
            {
                returnval = m_IntBlocks;
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                         "If statement should not arrive here");
            }

            level--;
            return returnval;
        }

        int BottomUpSubStructuredGraph::GetNumGlobalDofs(
            const int whichlevel) const
        {
            int returnval = -1;
            static int level = 0;
            level++;

            if(level < whichlevel)
            {
                returnval = m_daughterGraph->GetNumGlobalDofs(whichlevel);
            }
            else if(level == whichlevel)
            {
                returnval = m_IntBlocks[m_IntBlocks.size()-1]->GetIdOffset()+
                    m_IntBlocks[m_IntBlocks.size()-1]->GetNverts();
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                         "If statement should not arrive here");
            }

            level--;

            return returnval;
        }

        int BottomUpSubStructuredGraph::GetNlevels() const
        {
            int returnval = 0;
            static int level = 0;
            level++;

            if( m_daughterGraph.get())
            {
                returnval = m_daughterGraph->GetNlevels();
            }
            returnval = max(returnval,level);

            level--;

            return returnval;

        }

        bool SubGraphWithoutVerts(const SubGraphSharedPtr g)
        {
            return g->GetNverts() == 0;
        }

        namespace
        {
            // the first template parameter (=OutEdgeList) is chosen to be of
            // type std::set as in the set up of the adjacency, a similar edge
            // might be created multiple times.  And to prevent the definition
            // of parallell edges, we use std::set (=boost::setS) rather than
            // std::vector (=boost::vecS)
            typedef boost::adjacency_list<
                boost::setS, boost::vecS, boost::undirectedS> BoostGraph;
            typedef boost::graph_traits<BoostGraph>::vertex_descriptor
                BoostVertex;
        }

        void CuthillMckeeReordering(const BoostGraph& graph,
                                    Array<OneD, int>& perm,
                                    Array<OneD, int>& iperm)
        {
            int nGraphVerts = boost::num_vertices(graph);

            ASSERTL1(perm. size() >= nGraphVerts &&
                     iperm.size() >= nGraphVerts,
                     "Non-matching dimensions");

            // Call boost::cuthill_mckee_ordering to reorder the graph-vertices
            // using the reverse Cuthill-Mckee algorithm
            std::vector<BoostVertex> reorderedVerts(nGraphVerts);
            boost::cuthill_mckee_ordering(graph, reorderedVerts.rbegin());

            //copy the reordering to the Arrays perm and iperm
            for(int i = 0; i < nGraphVerts; i++)
            {
                perm[i] = reorderedVerts[i];
                iperm[ reorderedVerts[i] ] = i;
            }
        }

        void MultiLevelBisectionReordering(
            const BoostGraph                    &graph,
            Array<OneD, int>                    &perm,
            Array<OneD, int>                    &iperm,
            BottomUpSubStructuredGraphSharedPtr &substructgraph,
            std::set<int>                        partVerts,
            int                                  mdswitch)
        {
#ifndef NEKTAR_USE_SCOTCH
            boost::ignore_unused(graph, perm, iperm, substructgraph, partVerts,
                                 mdswitch);
            ASSERTL0(false, "Multi-level static condensation requires Nektar++"
                            " to be built with SCOTCH.");
#else
            int nGraphVerts = boost::num_vertices(graph);
            int nGraphEdges = boost::num_edges   (graph);

            ASSERTL1(perm. size() >= nGraphVerts &&
                     iperm.size() >= nGraphVerts,
                     "Non-matching dimensions");

            // We will now use Scotch to reorder the graph.  For the purpose of
            // multi-level static condensation, we will use a Scotch routine
            // that partitions the graph recursively using a multi-level nested
            // bisection algorithm.  The name of this routine is
            // SCOTCH_graphOrder and it was originally designed to reorder the
            // DOFs in a matrix in order to minimise the fill-in when applying a
            // factorisation technique (such as Cholesky).  However, this
            // reordering of DOFs also seems to be perfectly suited in the
            // context of multilevel substructuring. Therefore, we will use this
            // Scotch routine instead of the more well-known graph-partitioning
            // routines.
            if(nGraphEdges)
            {
                //
                // Step 1: Convert boost graph to a graph in adjncy-list format
                // as required by Scotch.
                //
                int acnt = 0, vcnt = 0, i, cnt;
                int nPartition    = partVerts.size();
                int nNonPartition = nGraphVerts - partVerts.size();

                Array<OneD, int> xadj(nNonPartition+1,0);
                Array<OneD, int> adjncy(2*nGraphEdges);
                Array<OneD, int> initial_perm(nGraphVerts);
                Array<OneD, int> iinitial_perm(nGraphVerts);
                Array<OneD, int> perm_tmp (nNonPartition);
                Array<OneD, int> iperm_tmp(nNonPartition);

                // Perform an initial reordering of the vertices, so that
                // partition nodes are at the end. This allows Scotch to
                // partition the interior nodes from values starting at zero.
                for (i = cnt = 0; i < nGraphVerts; ++i)
                {
                    if (partVerts.count(i) == 0)
                    {
                        initial_perm [i]     = cnt;
                        iinitial_perm[cnt++] = i;
                    }
                }

                for (i = 0; i < nGraphVerts; ++i)
                {
                    if (partVerts.count(i) > 0)
                    {
                        initial_perm [i]     = cnt;
                        iinitial_perm[cnt++] = i;
                    }
                }

                // Apply this reordering to the graph.
                boost::property_map<BoostGraph, boost::vertex_index_t>::type
                    index = get(boost::vertex_index, graph);

                // Now construct the adjaceny list using
                // boost::adjacent_vertices.
                auto verts = boost::vertices(graph);
                for (auto vertit = verts.first; vertit != verts.second; ++vertit)
                {
                    if (partVerts.count(index[*vertit]) > 0)
                    {
                        continue;
                    }

                    auto adjverts = boost::adjacent_vertices(*vertit,graph);
                    for (auto adjvertit = adjverts.first;
                         adjvertit != adjverts.second; ++adjvertit)
                    {
                        if (partVerts.count(index[*adjvertit]) > 0)
                        {
                            continue;
                        }
                        adjncy[acnt++] = initial_perm[*adjvertit];
                    }
                    xadj[++vcnt] = acnt;
                }

                //
                // Step 2: pass the graph to Scotch and perform the nested
                // dissection to obtain a separator tree.
                //

                // Pass the adjaceny graph into Scotch.
                SCOTCH_Graph scGraph;
                SCOTCH_CALL(SCOTCH_graphBuild,
                            (&scGraph, 0, nNonPartition, &xadj[0], &xadj[1],
                             NULL, NULL, xadj[nNonPartition], &adjncy[0], NULL));

                // This horrible looking string defines the Scotch graph
                // reordering strategy, which essentially does a nested
                // dissection + compression. We take this almost directly from
                // the SCOTCH_stratGraphOrderBuild function (defined in
                // library_graph_order.c), but by specifying the string
                // manually, we can replace the subdivision strategy to allow us
                // to control the number of vertices used to determine whether
                // to perform another dissection using the mdswitch
                // parameter. The below is essentially equivalent to calling
                // SCOTCH_stratGraphOrderBuild with the flags
                // SCOTCH_STRATLEAFSIMPLE and SCOTCH_STRATSEPASIMPLE to make
                // sure leaf nodes do not have any reordering applied to them.
                std::string strat_str =
                    "c{rat=0.7,cpr=n{sep=/(<TSTS>)?m{rat=0.7,vert=100,low="
                    "h{pass=10},asc=b{width=3,bnd=f{bal=<BBAL>},"
                    "org=(|h{pass=10})f{bal=<BBAL>}}}<SEPA>;,"
                    "ole=<OLEA>,ose=<OSEP>},unc=n{sep=/(<TSTS>)?m{rat=0.7,"
                    "vert=100,low=h{pass=10},asc=b{width=3,bnd=f{bal=<BBAL>},"
                    "org=(|h{pass=10})f{bal=<BBAL>}}}<SEPA>;"
                    ",ole=<OLEA>,ose=<OSEP>}}";

                // Replace flags in the string with appropriate values.
                boost::replace_all(
                    strat_str, "<SEPA>", "|m{rat=0.7,vert=100,low=h{pass=10},"
                    "asc=b{width=3,bnd=f{bal=<BBAL>},"
                    "org=(|h{pass=10})f{bal=<BBAL>}}}");
                boost::replace_all(strat_str, "<OSEP>", "s");
                boost::replace_all(strat_str, "<OLEA>", "s");
                boost::replace_all(strat_str, "<BBAL>", "0.1");
                boost::replace_all(
                    strat_str, "<TSTS>",
                    "vert>"+std::to_string(mdswitch));

                // Set up the re-ordering strategy.
                SCOTCH_Strat strat;
                SCOTCH_CALL(SCOTCH_stratInit, (&strat));
                SCOTCH_CALL(SCOTCH_stratGraphOrder, (&strat, strat_str.c_str()));

                // As output, Scotch will give us the total number of 'blocks'
                // (i.e. the separators and all of the leaves), the separator
                // tree as a mapping of block to parent block, and the range of
                // indices that is contained within each block. Reordering of
                // the vertices goes from largest index (at the top level) to
                // smallest (at the bottom level). The precise ordering is given
                // in the Scotch user guide.
                //
                // Note that we pass in iperm into the 'permtab' field of
                // graphOrder and 'perm' into the 'peritab' field; this is
                // because our definition of the permutation is the opposite of
                // what's defined in Scotch (this is leftover from a previous
                // implementation that used Metis).
                Array<OneD, int> treetab(nNonPartition);
                Array<OneD, int> rangtab(nNonPartition + 1);
                int cblknbr = 0;
                SCOTCH_CALL(SCOTCH_graphOrder,
                            (&scGraph, &strat, &iperm_tmp[0], &perm_tmp[0],
                             &cblknbr, &rangtab[0], &treetab[0]));

                // We're now done with Scotch: clean up the created structures.
                SCOTCH_graphExit(&scGraph);
                SCOTCH_stratExit(&strat);

                //
                // Step 3: create a MultiLevelBisectedGraph by reading the
                // separator tree we obtained from Scotch.
                //

                // Setup root block, which lies at the end of the blocks
                // described in treetab[].
                std::vector<MultiLevelBisectedGraphSharedPtr> graphs(cblknbr);

                // The strategy now is to traverse backwards over the blocks
                // described in treetab to set up the levels of the top-down
                // graph. rangtab allows us to calculate how many degrees of
                // freedom lie in the separator.
                for (i = cblknbr-1; i >= 0; --i)
                {
                    // Set up this block.
                    graphs[i] = MemoryManager<MultiLevelBisectedGraph>
                        ::AllocateSharedPtr(rangtab[i+1] - rangtab[i]);

                    // If we're a root block (treetab[i] == -1) we don't need to
                    // do anything, just move onto the next block.
                    if (treetab[i] == -1)
                    {
                        continue;
                    }

                    // Now use treetab[i] to figure out the parent block.  We
                    // have to be a bit careful in setting left/right daughters
                    // here. The left daughter's degrees of freedom are ordered
                    // _first_ in the iperm/perm arrays returned from Scotch,
                    // but if there is both a left and right daughter, we'll
                    // come across the right daughter first because the
                    // separators are being traversed backwards. We'll therefore
                    // insert this at the beginning of the daughter graphs
                    // vector.
                    MultiLevelBisectedGraphSharedPtr tmp = graphs[treetab[i]];
                    std::vector<MultiLevelBisectedGraphSharedPtr> &daughters =
                        tmp->GetDaughterGraphs();
                    daughters.insert(daughters.begin(), graphs[i]);
                }

                // Change permutations from Scotch to account for initial offset
                // in case we had partition vertices.
                for (i = 0; i < nGraphVerts; ++i)
                {
                    if (partVerts.count(i) == 0)
                    {
                        iperm[i] = iperm_tmp[initial_perm[i]];
                        perm[iperm[i]] = i;
                    }
                }

                auto it = partVerts.begin(), it2 = partVerts.end();
                for (i = nNonPartition; it != it2; ++it, ++i)
                {
                    perm [i]   = *it;
                    iperm[*it] = i;
                }

                for (i = 0; i < nGraphVerts; ++i)
                {
                    ASSERTL1(perm[iperm[i]] == i, "Perm error "
                             + std::to_string(i));
                }

                // If we were passed a graph with disconnected regions, we need
                // to create a bisected graph with the appropriate roots.
                std::vector<int> rootBlocks;
                for (i = 0; i < cblknbr; ++i)
                {
                    if (treetab[i] == -1)
                    {
                        rootBlocks.push_back(i);
                    }
                }

                MultiLevelBisectedGraphSharedPtr root;
                if (rootBlocks.size() == 1)
                {
                    root = graphs[rootBlocks[0]];
                }
                else
                {
                    root = MemoryManager<MultiLevelBisectedGraph>
                        ::AllocateSharedPtr(0);

                    for (int i = 0; i < rootBlocks.size(); ++i)
                    {
                        root->GetDaughterGraphs().push_back(graphs[rootBlocks[i]]);
                    }
                }

                // Check that our degree of freedom count in the constructed
                // graph is the same as the number of degrees of freedom that we
                // were given in the function input.
                ASSERTL0(root->GetTotDofs() == nNonPartition,
                         "Error in constructing Scotch graph for multi-level"
                         " static condensation.");

                //
                // Step 4: Set up the bottom-up graph from the top-down graph,
                // and reorder the permutation from Scotch.
                //
                substructgraph = MemoryManager<BottomUpSubStructuredGraph>::
                    AllocateSharedPtr(root, nPartition, true);

                // Important: we cannot simply use the ordering given by Scotch
                // as it does not order the different blocks as we would like
                // it. Therefore, we use following command to re-order them
                // again in the context of the bottom-up substructuring. As a
                // result, we will now obtain an ordering where the interior
                // degrees of freedom of the first (=bottom) level will be
                // ordered last (block by block ofcoarse). The interior degrees
                // of freedom of the second level will be ordered second to
                // last, etc ... As a result, the boundary degrees of freedom of
                // the last level (i.e. the dofs that will have to solved
                // non-recursively) will be ordered first (after the Dirichlet
                // Dofs that is).  (this way, we actually follow the same idea
                // and convention in the standard (non-multi-level) static
                // condensation approach).
                substructgraph->UpdateBottomUpReordering(perm,iperm);
            }
            else
            {
                // This is the very special case of a graph without any
                // connectivity i.e. a collection of vertices without any edges
                for(int i = 0; i < nGraphVerts; i++)
                {
                    perm[i] = i;
                    iperm[i] = i;
                }
                substructgraph = MemoryManager<BottomUpSubStructuredGraph>::
                    AllocateSharedPtr(nGraphVerts);
            }
#endif
        }

        void NoReordering(const BoostGraph& graph,
                          Array<OneD, int>& perm,
                          Array<OneD, int>& iperm)
        {
            int nGraphVerts = boost::num_vertices(graph);

            ASSERTL1(perm. size() >= nGraphVerts &&
                     iperm.size() >= nGraphVerts,
                     "Non-matching dimensions");

            for (int i = 0; i < nGraphVerts; i++)
            {
                perm [i] = i;
                iperm[i] = i;
            }
        }
    }
}
