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
// Description: a collection of classes that facilitates the implementation
//              of the multi-level static condensation routines
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/SubStructuredGraph.h>
#include <LibUtilities/BasicUtils/Metis.hpp>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

namespace Nektar
{
    namespace MultiRegions
    {
       
        PatchMap::PatchMap(void)
        {
        };

        PatchMap::PatchMap(const int nvals)
        {
            m_patchId  = Array<OneD, int> (nvals);
            m_dofId    = Array<OneD, int> (nvals);
#if 0 
            m_bndPatch = Array<OneD, bool>(nvals);
#else
            m_bndPatch = Array<OneD, unsigned int>(nvals);
#endif
            m_sign     = Array<OneD, NekDouble>(nvals);
        }
            
        PatchMap::~PatchMap(void)
        {
        }
        

#if 0
        void PatchMap::SetPatchMap(const int n, const int patchId, const int dofId,const bool bndPatch,const NekDouble sign)
#else
        void PatchMap::SetPatchMap(const int n, const int patchId, const int dofId,const unsigned int bndPatch,const NekDouble sign)
#endif
        {
            m_patchId[n]  = patchId;
            m_dofId[n]    = dofId;
            m_bndPatch[n] = bndPatch;
            m_sign[n]     = sign;
        }
        
        MultiLevelBisectedGraph::MultiLevelBisectedGraph(const Array<OneD, const int> sepTree):
            m_BndDofs(),
            m_leftDaughterGraph(),
            m_rightDaughterGraph()
        {
            static int offset = -5;
            offset += 5;

            int recurLevel    = sepTree[offset+0];
            int branch        = sepTree[offset+1];
            int nLeftIntDofs  = sepTree[offset+2];
            int nRightIntDofs = sepTree[offset+3];
            int nBndDofs      = sepTree[offset+4];

            bool daughtersConstructed[2] = {false,false};
            
            if( (offset + 5) < sepTree.num_elements() )
            {
                while( sepTree[offset+5] > recurLevel ) 
                {         
                    switch(sepTree[offset+6])
                    {
                    case 1:
                        {
                            m_leftDaughterGraph = MemoryManager<MultiLevelBisectedGraph>::AllocateSharedPtr(sepTree);
                            daughtersConstructed[0] = true;
                        }
                        break;
                    case 2:
                        {
                            m_rightDaughterGraph = MemoryManager<MultiLevelBisectedGraph>::AllocateSharedPtr(sepTree);
                            daughtersConstructed[1] = true;
                        }
                        break;
                    default:
                        {
                            NEKERROR(ErrorUtil::efatal,"Invalid branch id");
                        }
                    }
                    if( (offset + 5) >= sepTree.num_elements() )
                    {
                        break;
                    }
                }
            }

            m_BndDofs = MemoryManager<SubGraph>::AllocateSharedPtr(nBndDofs);
                
            if( !daughtersConstructed[0] && nLeftIntDofs )
            {
                m_leftDaughterGraph = MemoryManager<MultiLevelBisectedGraph>::AllocateSharedPtr(nLeftIntDofs);
            }
            
            if( !daughtersConstructed[1] && nRightIntDofs )
            {
                m_rightDaughterGraph = MemoryManager<MultiLevelBisectedGraph>::AllocateSharedPtr(nRightIntDofs);
            }

            if( recurLevel == 1 )
            {
                offset = -5;
            }
        }

        MultiLevelBisectedGraph::MultiLevelBisectedGraph(const int nBndDofs):
            m_BndDofs(MemoryManager<SubGraph>::AllocateSharedPtr(nBndDofs)),
            m_leftDaughterGraph(),
            m_rightDaughterGraph()
        {
        }
        
        MultiLevelBisectedGraph::~MultiLevelBisectedGraph(void)
        {
        }

        int MultiLevelBisectedGraph::GetTotDofs() const
        {
            static int nBndDofs = 0;
            static int level = 0;
            level++;

            int returnval;

            if(m_leftDaughterGraph.get())
            {
                m_leftDaughterGraph->GetTotDofs();
            }
            if(m_rightDaughterGraph.get())
            {
                m_rightDaughterGraph->GetTotDofs();
            }

            nBndDofs += m_BndDofs->GetNverts();
            returnval = nBndDofs;
 
            level--;
            if(level == 0)
            {
                nBndDofs = 0;
            }

            return returnval;
        }

        void MultiLevelBisectedGraph::SetGlobalNumberingOffset()
        {
            static int level = 0;
            static int offset = 0;
            level++;

            if(m_leftDaughterGraph.get())
            {
                m_leftDaughterGraph->SetGlobalNumberingOffset();
            }
            if(m_rightDaughterGraph.get())
            {
                m_rightDaughterGraph->SetGlobalNumberingOffset();
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
            if(m_leftDaughterGraph.get())
            {
                m_leftDaughterGraph->DumpNBndDofs();
            }
            if(m_rightDaughterGraph.get())
            {
                m_rightDaughterGraph->DumpNBndDofs();
            }
            cout << "  "  << m_BndDofs->GetNverts() << "   " << m_BndDofs->GetIdOffset() << endl;
        }

        void MultiLevelBisectedGraph::CollectLeaves(vector<SubGraphSharedPtr>& leaves) const
        {
            int cnt = 0;
            if(m_leftDaughterGraph.get())
            {
                m_leftDaughterGraph->CollectLeaves(leaves);
                cnt++;
            }
            if(m_rightDaughterGraph.get())
            {
                m_rightDaughterGraph->CollectLeaves(leaves);
                cnt++;
            }

            if(cnt == 0)
            {
                SubGraphSharedPtr leave = m_BndDofs;
                leaves.push_back(leave);
            }
        }

        inline int MultiLevelBisectedGraph::GetNdaughterGraphs() const
        {
            int cnt = 0;
            if(m_leftDaughterGraph.get())
            {
                cnt++;
            }
            if(m_rightDaughterGraph.get())
            {
                cnt++;
            }
            return cnt;
        } 

        int MultiLevelBisectedGraph::CutEmptyLeaves()
        {
            int returnval;
            static int level = 0;
            static int nLeaves = 0;
            level++;

            if( (level == 1) && 
                (!m_leftDaughterGraph.get()) && 
                (!m_rightDaughterGraph.get()) )
            {
                level = 0;
                nLeaves = 0;
                return 0;
            }

            if(m_leftDaughterGraph.get())
            {
                if( (m_leftDaughterGraph->GetNdaughterGraphs() == 0) && 
                    (m_leftDaughterGraph->GetBndDofsGraph()->GetNverts() == 0) )
                {
                    m_leftDaughterGraph = MultiLevelBisectedGraphSharedPtr();
                    nLeaves++;
                }
                else
                {
                    m_leftDaughterGraph->CutEmptyLeaves();
                }
            }
            if(m_rightDaughterGraph.get())
            {
                if( (m_rightDaughterGraph->GetNdaughterGraphs() == 0) &&
                    (m_rightDaughterGraph->GetBndDofsGraph()->GetNverts() == 0) )
                {
                    m_rightDaughterGraph = MultiLevelBisectedGraphSharedPtr();
                    nLeaves++;
                }
                else
                {
                    m_rightDaughterGraph->CutEmptyLeaves();
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

            if( (level == 1) && 
                (!m_leftDaughterGraph.get()) && 
                (!m_rightDaughterGraph.get()) )
            {
                level = 0;
                nLeaves = 0;
                return 0;
            }

            if(m_leftDaughterGraph.get())
            {
                if(m_leftDaughterGraph->GetNdaughterGraphs() == 0)
                {
                    m_leftDaughterGraph = MultiLevelBisectedGraphSharedPtr();
                    nLeaves++;
                }
                else
                {
                    m_leftDaughterGraph->CutLeaves();
                }
            }
            if(m_rightDaughterGraph.get())
            {
                if(m_rightDaughterGraph->GetNdaughterGraphs() == 0)
                {
                    m_rightDaughterGraph = MultiLevelBisectedGraphSharedPtr();
                    nLeaves++;
                }
                else
                {
                    m_rightDaughterGraph->CutLeaves();
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


        BottomUpSubStructuredGraph::BottomUpSubStructuredGraph(const int nVerts):
            m_IntBlocks(),
            m_daughterGraph()
        {
            // This constructor should only be used in the very special case
            // when there is a graph consisting of nVerts vertices but without
            // any connectivity whatsowever (i.e. no edges)

            SubGraphSharedPtr subgraph;
            const int one = 1;
            for(int i = 0 ; i < nVerts ; i++)
            {
                subgraph = MemoryManager<SubGraph>::AllocateSharedPtr(one,i);
                m_IntBlocks.push_back(subgraph);
            } 
        }

        BottomUpSubStructuredGraph::BottomUpSubStructuredGraph(const Array<OneD, const int> septree):
            m_IntBlocks(),
            m_daughterGraph()
        {
            // First, create a top-down graph structure based upon the
            // separator tree. This is easier as separation tree is also 
            // structured following a top-down approach
            MultiLevelBisectedGraphSharedPtr topDownGraph = MemoryManager<MultiLevelBisectedGraph>::
                AllocateSharedPtr(septree);
            // set the global numbering of the top-down graph
            topDownGraph->SetGlobalNumberingOffset();

            topDownGraph->CutEmptyLeaves();
            // Secondly, recursively construct the subgraphs of the
            // bottom up point of view 1. Collect all the leaves of
            // the topdown graph this will be the first level of the
            // bottom up graph
            topDownGraph->CollectLeaves(m_IntBlocks);
            // 2. Reduce the topdown graph by cutting the leaves (this
            //    will allow a recursive approach)
            int ncuts = topDownGraph->CutLeaves();
            // 3. If there were leaves to cut, proceed recursively
            if(ncuts)
            {
                m_daughterGraph = MemoryManager<BottomUpSubStructuredGraph>::AllocateSharedPtr(topDownGraph);
            }
        }

        BottomUpSubStructuredGraph::BottomUpSubStructuredGraph(const MultiLevelBisectedGraphSharedPtr& graph):
            m_IntBlocks(),
            m_daughterGraph()
        {
            int ncuts;
            graph->CutEmptyLeaves();
            graph->CollectLeaves(m_IntBlocks);
            ncuts = graph->CutLeaves();

            if(ncuts)
            {
                m_daughterGraph = MemoryManager<BottomUpSubStructuredGraph>::AllocateSharedPtr(graph);
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
            for(int i = 0; i < m_IntBlocks.size(); i++)
            {
                cout << "  " << i << "/"<< m_IntBlocks.size()-1<< ": " << m_IntBlocks[i]->GetNverts() << ", " << m_IntBlocks[i]->GetIdOffset() << endl;
            }
            cout << endl;
            
            if( m_daughterGraph.get())
            {
                m_daughterGraph->Dump();
            }
 
            level--;
        }

        void BottomUpSubStructuredGraph::UpdateBottomUpReordering(Array<OneD,       int>& perm, 
                                                                  Array<OneD,       int>& iperm) const
        {
            int nDofs = GetTotDofs();            
            
            // Step 1: make a permutation array
            // that goes from the current reordering 
            // in the bottom-up graph to an ordering
            // where the interior dofs of the first (=bottom)
            // level are ordered last, followed interior dofs
            // of the second level, ...
            Array<OneD, int> iperm1(nDofs);
            SetBottomUpReordering(iperm1);
            
            // now, based upon the input permutation array, update
            // the permutation arrays between the original ordering of the dofs 
            // (defined somewhere outside) and the final reordering as
            // given by iperm1
            for(int i=0; i < nDofs; i++)
            {
                iperm[i]       = iperm1[ iperm[i] ];
                perm[iperm[i]] = i;
            }
        }
        
        void BottomUpSubStructuredGraph::SetBottomUpReordering(Array<OneD, int>& iperm) const
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

        void BottomUpSubStructuredGraph::ExpandGraphWithVertexWeights(const Array<OneD, const int>& wgts)
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
                int OrigNverts       = m_IntBlocks[i]->GetNverts();
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
            m_IntBlocks.erase(remove_if(m_IntBlocks.begin(), m_IntBlocks.end(), SubGraphWithoutVerts), m_IntBlocks.end());
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

        void BottomUpSubStructuredGraph::MaskPatches(const int leveltomask, Array<OneD, NekDouble>& maskarray) const
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
                NEKERROR(ErrorUtil::efatal,"If statement should not arrive here");
            }
 
            level--;
        }

        int BottomUpSubStructuredGraph::GetNpatchesWithInterior(const int whichlevel) const
        {
            int returnval;
            static int level = 0;
            level++;

            if(level < whichlevel)
            {
                returnval = m_daughterGraph->GetNpatchesWithInterior(whichlevel);
            }
            else if(level == whichlevel)
            {
                returnval = m_IntBlocks.size();
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,"If statement should not arrive here");
            }
 
            level--;

            return returnval;
        }

        void BottomUpSubStructuredGraph::GetNintDofsPerPatch(const int whichlevel, Array<OneD, unsigned int>& outarray) const
        {
            static int level = 0;
            level++;

            if(level < whichlevel)
            {
                m_daughterGraph->GetNintDofsPerPatch(whichlevel,outarray);
            }
            else if(level == whichlevel)
            {
                ASSERTL1(outarray.num_elements() >= m_IntBlocks.size(),"Array dimension not sufficient");
                for(int i = 0; i < m_IntBlocks.size(); i++)
                {
                    outarray[i] = (unsigned int) m_IntBlocks[i]->GetNverts();
                }
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,"If statement should not arrive here");
            }
 
            level--;
        }
            
        int BottomUpSubStructuredGraph::GetInteriorOffset(const int whichlevel, const int patch) const
        {
            int returnval;
            static int level = 0;
            level++;

            if(level < whichlevel)
            {
                returnval = m_daughterGraph->GetInteriorOffset(whichlevel,patch);
            }
            else if(level == whichlevel)
            {
                returnval = m_IntBlocks[patch]->GetIdOffset();
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,"If statement should not arrive here");
            }
 
            level--;

            return returnval;
        }

        vector<SubGraphSharedPtr> BottomUpSubStructuredGraph::GetInteriorBlocks(const int whichlevel) const
        {
            vector<SubGraphSharedPtr> returnval;
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
                NEKERROR(ErrorUtil::efatal,"If statement should not arrive here");
            }
 
            level--;
            return returnval; 
        }

        int BottomUpSubStructuredGraph::GetNumGlobalDofs(const int whichlevel) const
        {
            int returnval;
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
                NEKERROR(ErrorUtil::efatal,"If statement should not arrive here");
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
            bool returnvalue = false;
            if( g->GetNverts() == 0 )
            {
                returnvalue = true;
            }
            return returnvalue;
        } 

        namespace
        {
                // the first template parameter (=OutEdgeList) is
                // chosen to be of type std::set as in the set up of
                // the adjacency, a similar edge might be created
                // multiple times.  And to prevent the definition of
                // parallell edges, we use std::set (=boost::setS)
                // rather than std::vector (=boost::vecS)
                typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> BoostGraph;
                typedef boost::graph_traits<BoostGraph>::vertex_descriptor                  BoostVertex;
                typedef boost::graph_traits<BoostGraph>::vertex_iterator                    BoostVertexIterator;
                typedef boost::graph_traits<BoostGraph>::adjacency_iterator                 BoostAdjacencyIterator;
        }

        void CuthillMckeeReordering(const BoostGraph& graph,
                                    Array<OneD, int>& perm,
                                    Array<OneD, int>& iperm)
        {
            int nGraphVerts = boost::num_vertices(graph);

            ASSERTL1(perm.num_elements()>=nGraphVerts,"Non-matching dimensions");
            ASSERTL1(iperm.num_elements()>=nGraphVerts,"Non-matching dimensions");

            // Call boost::cuthill_mckee_ordering to reorder the
            // graph-vertices using the reverse Cuthill-Mckee
            // algorithm
            vector<BoostVertex> reorderedVerts(nGraphVerts);
            boost::cuthill_mckee_ordering(graph, reorderedVerts.rbegin());

            //copy the reordering to the Arrays perm and iperm
            for(int i = 0; i < nGraphVerts; i++)
            {
                perm[i] = reorderedVerts[i];
                iperm[ reorderedVerts[i] ] = i;
            }
        }

        void MultiLevelBisectionReordering(const BoostGraph& graph,
                                           const Array<OneD, const int>& vwgts,
                                           Array<OneD, int>& perm,
                                           Array<OneD, int>& iperm,
                                           BottomUpSubStructuredGraphSharedPtr& substructgraph,
                                           const int mdswitch)
        {
            int nGraphVerts = boost::num_vertices(graph);
            int nGraphEdges = boost::num_edges(graph);

            ASSERTL1(perm.num_elements()>=nGraphVerts,"Non-matching dimensions");
            ASSERTL1(iperm.num_elements()>=nGraphVerts,"Non-matching dimensions");

            // We will now use METIS to reorder the graph.  For the
            // purpose of multi-level static condensation, we will use
            // a METIS routine that partitions the graph recursively
            // using a multi-level bisection algorithm.  The name of
            // this routine is METIS_NodeND and it was originally
            // designed to reorder the DOFs in a matrix in order to
            // minimise the fill-in when applying a factorisation
            // technique (such as Cholesky).  However, this reordering
            // of DOFs also seems to be perfectly suited in the
            // context of multilevel substructering. Therefore, we
            // will use this metis routine instead of the more
            // well-known graph-partitioning routines. However, as the
            // standard metis implementation of METIS_NodeND only
            // gives the resulting re-ordering as an output, we we
            // will use an modified version of this routine that also
            // returns information about the structure of the
            // multi-level bisected partitioning.

            // This modified implementation has been written by W. GAO
            // and collaborators and it additionally returns the
            // separator tree compared to the standard implementation.

            // The name of this modified routine AS_METIS_NodeND (where AS stands for 
            // automated substructering)
            // More information can be found in the paper:
            //  W. Gao, S. Li Xiaoye, C. Yang and Z. Bai
            //  'An implementation and evaluation of the AMLS method 
            //   for sparse eigenvalue problems'
            //  ACM Trans. Math. Softw. 34, 4, Article 20 (July 2008)
            
            if(nGraphEdges)
            {
                // Step 1: Convert boost graph to a graph in adjncy-list format
                // as required by METIS
                int acnt = 0;
                int vcnt = 0;
                BoostVertexIterator    vertit, vertit_end;
                BoostAdjacencyIterator adjvertit, adjvertit_end;
                Array<OneD, int> xadj(nGraphVerts+1,0);
                Array<OneD, int> adjncy(2*nGraphEdges);
                
                for ( boost::tie(vertit, vertit_end) = boost::vertices(graph); 
                      vertit != vertit_end; 
                      ++vertit) 
                {
                    for ( boost::tie(adjvertit, adjvertit_end) = boost::adjacent_vertices(*vertit,graph); 
                          adjvertit != adjvertit_end;
                          ++adjvertit) 
                    {
                        adjncy[acnt++] = *adjvertit;
                        
                    }
                    xadj[++vcnt] = acnt;
                }
                
                // Step 2: use metis to reorder the dofs
                // We do not know on forehand the size of the
                // separator tree that METIS will return, so we just 
                // assume a really big value and try with that
                int sizeSeparatorTree = nGraphVerts*10;
                Array<OneD,int> septreeTmp(sizeSeparatorTree,-1);  
                
                // The separatortree returned by metis has the following structure:
                // It is a one dimensional array and information per level
                // is contained per 5 elements:
                // m_septree[i*5 + 0]: the level of recursion (top-level = 1)
                // m_septree[i*5 + 1]: is this substructure a left or right branch? (1 or 2)
                //                                1: left branch
                //                                2: right branch
                // m_septree[i*5 + 2]: the number of 'interior' DOFs in left branch
                // m_septree[i*5 + 3]: the number of 'interior' DOFs in right branch
                // m_septree[i*5 + 4]: the number of 'boundary' DOFs     
                
                // now try to call Call METIS 
                try
                {
                    Metis::as_onmetis(nGraphVerts,xadj,adjncy,perm,iperm,septreeTmp,mdswitch);
                }
                catch(...)
                {
                    NEKERROR(ErrorUtil::efatal,
                             "Error in calling metis (the size of the separator tree might not be sufficient)");
                }
                
                // Post-process the separatortree
                int trueSizeSepTree=0;
                for(int i = 0 ; septreeTmp[i] != -1; i++)
                {
                    trueSizeSepTree++;
                }            
                Array<OneD,int> septree(trueSizeSepTree);
                Vmath::Vcopy(trueSizeSepTree,septreeTmp,1,septree,1);
                
                // Based upon the separator tree, where are going to
                // set up an object of the class
                // BottomUpSubStructuredGraph. The constructor will
                // read the separatortree and will interprete the
                // information from a bottom-up point of view.
                substructgraph = MemoryManager<BottomUpSubStructuredGraph>::AllocateSharedPtr(septree);
                // Important, we cannot simply use the ordering given
                // by metis as it does not order the different blocks
                // as we would like it. Therefore, we use following
                // command to re-order them again in the context of
                // the bottom-up substructering. As a result, we will
                // now obtain an ordering where the interior degrees
                // of freedom of the first (=bottom) level will be
                // ordered last (block by block ofcoarse). The
                // interior degrees of freedom of the second level
                // will be ordered second to last, etc ... As a
                // result, the boundary degrees of freedom of the last
                // level (i.e. the dofs that will have to solved
                // non-recursively) will be ordered first (after the
                // Dirichlet Dofs that is).  (this way, we actually
                // follow the same idea and convention in the standard
                // (non-multi-level) static condensation approach)
                substructgraph->UpdateBottomUpReordering(perm,iperm);
            }
            else
            {
                // This is the very special case of a graph without any connectivity
                // i.e. a collection of vertices without any edges
                for(int i = 0; i < nGraphVerts; i++)
                {
                    perm[i] = i;
                    iperm[i] = i;
                }
                substructgraph = MemoryManager<BottomUpSubStructuredGraph>::AllocateSharedPtr(nGraphVerts);
            }

#if 0 // Removed this operation to higher routine. 
            Array<OneD,int> vwgts_perm(nGraphVerts);
            for(int i = 0; i < nGraphVerts; i++)
            {
                vwgts_perm[i] = vwgts[perm[i]];
            }
            substructgraph->ExpandGraphWithVertexWeights(vwgts_perm);
            cout << "Expanded graph" <<endl;
            substructgraph->Dump();
#endif
        }

        void NoReordering(const BoostGraph& graph,
                          Array<OneD, int>& perm,
                          Array<OneD, int>& iperm)
        {
            int nGraphVerts = boost::num_vertices(graph);

            ASSERTL1(perm.num_elements()>=nGraphVerts,"Non-matching dimensions");
            ASSERTL1(iperm.num_elements()>=nGraphVerts,"Non-matching dimensions");

            for(int i = 0; i < nGraphVerts; i++)
            {
                perm[i] = i;
                iperm[i] = i;
            }

        }

    }
}

/**
 * $Log: SubStructuredGraph.cpp,v $
 * Revision 1.7  2010/01/03 12:18:26  cantwell
 * Removed unnecessary header file inclusions.
 *
 * Revision 1.6  2009/12/21 17:50:09  cantwell
 * Restored CutEmptyLeaves calls.
 *
 * Revision 1.5  2009/12/15 18:09:03  cantwell
 * Split GeomFactors into 1D, 2D and 3D
 * Added generation of tangential basis into GeomFactors
 * Updated ADR2DManifold solver to use GeomFactors for tangents
 * Added <GEOMINFO> XML session section support in MeshGraph
 * Fixed const-correctness in VmathArray
 * Cleaned up LocalRegions code to generate GeomFactors
 * Removed GenSegExp
 * Temporary fix to SubStructuredGraph
 * Documentation for GlobalLinSys and GlobalMatrix classes
 *
 * Revision 1.4  2009/11/19 11:41:07  pvos
 * Fixed various bugs
 *
 * Revision 1.3  2009/11/09 15:57:11  pvos
 * multi-level recursion bug fixes
 *
 * Revision 1.2  2009/11/02 11:19:44  pvos
 * Fixed a bug for reordering a graph without edges
 *
 * Revision 1.1  2009/10/30 14:02:55  pvos
 * Multi-level static condensation updates
 *
 **/

