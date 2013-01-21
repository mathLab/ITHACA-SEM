///////////////////////////////////////////////////////////////////////////////
//
// File AssemblyMapCG2D.cpp
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
// Description: C0-continuous assembly mappings specific to 2D
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/AssemblyMap/AssemblyMapCG2D.h>
#include <MultiRegions/ExpList.h>
#include <LocalRegions/SegExp.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

#include <iomanip>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class AssemblyMapCG2D
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
        AssemblyMapCG2D::AssemblyMapCG2D(
                const LibUtilities::SessionReaderSharedPtr &pSession):
            AssemblyMapCG(pSession)
        {
        }




        /**
         *
         */
        AssemblyMapCG2D::AssemblyMapCG2D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const int numLocalCoeffs,
                const ExpList &locExp,
                const Array<OneD, const ExpListSharedPtr> &bndCondExp,
                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
                                                            &bndConditions,
                const vector<map<int,int> >& periodicVerticesId,
                const map<int,int>& periodicEdgesId,
                const bool checkIfSystemSingular) :
            AssemblyMapCG(pSession)
        {
            SetUp2DExpansionC0ContMap(numLocalCoeffs,
                                      locExp,
                                      bndCondExp,
                                      bndConditions,
                                      periodicVerticesId,
                                      periodicEdgesId,
                                      checkIfSystemSingular);

            CalculateBndSystemBandWidth();
            CalculateFullSystemBandWidth();
        }


        /**
         *
         */
        AssemblyMapCG2D::AssemblyMapCG2D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const int numLocalCoeffs,
                const ExpList &locExp):
            AssemblyMapCG(pSession)
        {
            SetUp2DExpansionC0ContMap(numLocalCoeffs, locExp);
            CalculateBndSystemBandWidth();
            CalculateFullSystemBandWidth();
        }


        /**
         *
         */
        AssemblyMapCG2D::~AssemblyMapCG2D()
        {
        }


        /**
         * Construction of the local->global map is achieved in
         * several stages.  A mesh vertex and mesh edge renumbering is
         * constructed in #vertReorderedGraphVertId and
         * #edgeReorderedGraphVertId through a call ot
         * #SetUp2DGraphC0ContMap. 
         *
         * The local numbering is then deduced in the following steps:
         */
        void AssemblyMapCG2D::SetUp2DExpansionC0ContMap(
                const int numLocalCoeffs,
                const ExpList &locExp,
                const Array<OneD, const ExpListSharedPtr> &bndCondExp,
                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
                                                                &bndConditions,
                const vector<map<int,int> >& periodicVerticesId,
                const map<int,int>& periodicEdgesId,
                const bool checkIfSystemSingular)
        {
            int i,j,k;
            int cnt = 0,offset=0;
            int meshVertId;
            int meshEdgeId;
            int bndEdgeCnt;
            int globalId, nGraphVerts;
            int nEdgeCoeffs;
            int nEdgeInteriorCoeffs;
            int firstNonDirGraphVertId;
            int nLocBndCondDofs = 0;
            int nLocDirBndCondDofs = 0;
            int nExtraDirichlet = 0;
            StdRegions::StdExpansion2DSharedPtr locExpansion;
            LocalRegions::SegExpSharedPtr       bndSegExp;
            LibUtilities::BasisType             bType;
            StdRegions::Orientation         edgeOrient;
            Array<OneD, unsigned int>           edgeInteriorMap;
            Array<OneD, int>                    edgeInteriorSign;
            const StdRegions::StdExpansionVector &locExpVector = *(locExp.GetExp());
            m_signChange = false;
            m_systemSingular = false;
            Array<OneD, map<int,int> > ReorderedGraphVertId(2);
            Array<OneD, map<int,int> > Dofs(2);
            BottomUpSubStructuredGraphSharedPtr bottomUpGraph;
            set<int> extraDirVerts;
            
            for(i = 0; i < locExpVector.size(); ++i)
            {
                for(j = 0; j < locExpVector[i]->GetNverts(); ++j)
                {
                    // Vert Dofs
                    Dofs[0][(locExpVector[i]->GetGeom2D())->GetVid(j)] = 1;
                    // Edge Dofs
                    Dofs[1][(locExpVector[i]->GetGeom2D())->GetEid(j)] = 
                        locExpVector[i]->GetEdgeNcoeffs(j)-2;
                }
            }

            /**
             * STEP 1: Wrap boundary conditions in an array (since
             * routine is set up for multiple fields) and call the
             * graph re-ordering subroutine to obtain the reordered
             * values
             */
            Array<OneD, Array<OneD, const SpatialDomains::BoundaryConditionShPtr> > bndConditionsVec(1,bndConditions);
            nGraphVerts = SetUp2DGraphC0ContMap(locExp,
                                                bndCondExp,bndConditionsVec,
                                                periodicVerticesId,periodicEdgesId,
                                                Dofs,
                                                ReorderedGraphVertId,
                                                firstNonDirGraphVertId,
                                                nExtraDirichlet,
                                                bottomUpGraph,
                                                extraDirVerts,
                                                checkIfSystemSingular);

            /**
             * STEP 2: Count out the number of Dirichlet vertices and
             * edges first
             */
            for(i = 0; i < bndCondExp.num_elements(); i++)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); j++)
                {
                    bndSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[i]->GetExp(j));

                    if(bndConditions[i]->GetBoundaryConditionType()==SpatialDomains::eDirichlet)
                    {
                        nLocDirBndCondDofs += bndSegExp->GetNcoeffs();
                    }
                    if( bndCondExp[i]->GetExp(j)==bndCondExp[bndCondExp.num_elements()-1]->GetExp(bndCondExp[bndCondExp.num_elements()-1]->GetExpSize()-1)
                            && i==(bndCondExp.num_elements()-1)
                            && nLocDirBndCondDofs ==0
                            && checkIfSystemSingular==true)
                    {
                        nLocDirBndCondDofs =1;
                        m_systemSingular=true;
                    }

                    nLocBndCondDofs += bndSegExp->GetNcoeffs();
                }
            }
            m_numLocalDirBndCoeffs = nLocDirBndCondDofs + nExtraDirichlet;

            /**
             * STEP 3: Set up an array which contains the offset information of
             * the different graph vertices.
             *
             * This basically means to identify to how many global degrees of
             * freedom the individual graph vertices correspond. Obviously,
             * the graph vertices corresponding to the mesh-vertices account
             * for a single global DOF. However, the graph vertices
             * corresponding to the element edges correspond to N-2 global DOF
             * where N is equal to the number of boundary modes on this edge.
             */
            Array<OneD, int> graphVertOffset(ReorderedGraphVertId[0].size()+
                                             ReorderedGraphVertId[1].size()+1);
            graphVertOffset[0] = 0;
            m_numLocalBndCoeffs = 0;
            for(i = 0; i < locExpVector.size(); ++i)
            {
                locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(locExpVector[locExp.GetOffset_Elmt_Id(i)]);
                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    nEdgeCoeffs = locExpansion->GetEdgeNcoeffs(j);
                    meshEdgeId = (locExpansion->GetGeom2D())->GetEid(j);
                    meshVertId = (locExpansion->GetGeom2D())->GetVid(j);
                    graphVertOffset[ReorderedGraphVertId[0][meshVertId]+1] = 1;
                    graphVertOffset[ReorderedGraphVertId[1][meshEdgeId]+1] = nEdgeCoeffs-2;

                    bType = locExpansion->GetEdgeBasisType(j);
                    // need a sign vector for modal expansions if nEdgeCoeffs >=4
                    if( (nEdgeCoeffs >= 4)&&
                        ( (bType == LibUtilities::eModified_A)||
                          (bType == LibUtilities::eModified_B) ) )
                    {
                        m_signChange = true;
                    }
                }
                m_numLocalBndCoeffs += locExpVector[i]->NumBndryCoeffs();
            }

            for(i = 1; i < graphVertOffset.num_elements(); i++)
            {
                graphVertOffset[i] += graphVertOffset[i-1];
            }

            // Allocate the proper amount of space for the class-data and fill
            // information that is already known
            m_numLocalCoeffs                 = numLocalCoeffs;
            m_numGlobalDirBndCoeffs          = graphVertOffset[firstNonDirGraphVertId];
            m_localToGlobalMap               = Array<OneD, int>(m_numLocalCoeffs,-1);
            m_localToGlobalBndMap            = Array<OneD, int>(m_numLocalBndCoeffs,-1);
            m_bndCondCoeffsToGlobalCoeffsMap = Array<OneD, int>(nLocBndCondDofs,-1);

            // If required, set up the sign-vector
            if(m_signChange)
            {
                m_localToGlobalSign = Array<OneD, NekDouble>(m_numLocalCoeffs,1.0);
                m_localToGlobalBndSign = Array<OneD, NekDouble>(m_numLocalBndCoeffs,1.0);
                m_bndCondCoeffsToGlobalCoeffsSign = Array<OneD, NekDouble>(nLocBndCondDofs,1.0);
            }
            else
            {
                m_localToGlobalSign = NullNekDouble1DArray;
                m_localToGlobalBndSign = NullNekDouble1DArray;
                m_bndCondCoeffsToGlobalCoeffsSign = NullNekDouble1DArray;
            }

            // Set up information for multi-level static condensation.
            m_staticCondLevel = 0;
            m_numPatches =  locExpVector.size();
            m_numLocalBndCoeffsPerPatch = Array<OneD, unsigned int>(m_numPatches);
            m_numLocalIntCoeffsPerPatch = Array<OneD, unsigned int>(m_numPatches);
            for(i = 0; i < m_numPatches; ++i)
            {
                int elmtid = locExp.GetOffset_Elmt_Id(i);
                locExpansion = boost::dynamic_pointer_cast<
                    StdRegions::StdExpansion2D>(locExpVector[elmtid]);
                m_numLocalBndCoeffsPerPatch[i] = (unsigned int) 
                    locExpVector[elmtid]->NumBndryCoeffs();
                m_numLocalIntCoeffsPerPatch[i] = (unsigned int) 
                    locExpVector[elmtid]->GetNcoeffs() - 
                    locExpVector[elmtid]->NumBndryCoeffs();
           }

            /**
             * STEP 4: Now, all ingredients are ready to set up the actual
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
                locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(locExpVector[i]);

                cnt = locExp.GetCoeff_Offset(i);

                // Loop over all edges (and vertices) of element i
                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    map<int,int>::const_iterator it;
                    
                    nEdgeInteriorCoeffs = locExpansion->GetEdgeNcoeffs(j)-2;
                    edgeOrient          = (locExpansion->GetGeom2D())->GetEorient(j);
                    meshEdgeId          = (locExpansion->GetGeom2D())->GetEid(j);
                    meshVertId          = (locExpansion->GetGeom2D())->GetVid(j);

                    /*
                     * Where the edge orientations of periodic edges matches,
                     * we reverse the vertices of each edge in
                     * DisContField2D::GetPeriodicEdges. We must therefore
                     * reverse the orientation of precisely one of the two
                     * edges so that the sign array is correctly populated.
                     */ 
                    it = periodicEdgesId.find(meshEdgeId);
                    if (it != periodicEdgesId.end() && it->second < 0)
                    {
                        if (edgeOrient == StdRegions::eForwards)
                        {
                            edgeOrient = StdRegions::eBackwards;
                        }
                        else
                        {
                            edgeOrient = StdRegions::eForwards;
                        }
                    }

                    locExpansion->GetEdgeInteriorMap(j,edgeOrient,edgeInteriorMap,edgeInteriorSign);
                    // Set the global DOF for vertex j of element i
                    m_localToGlobalMap[cnt+locExpansion->GetVertexMap(j)] =
                        graphVertOffset[ReorderedGraphVertId[0][meshVertId]];

                    // Set the global DOF's for the interior modes of edge j
                    for(k = 0; k < nEdgeInteriorCoeffs; ++k)
                    {
                        m_localToGlobalMap[cnt+edgeInteriorMap[k]] =
                            graphVertOffset[ReorderedGraphVertId[1][meshEdgeId]]+k;
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
                
                cnt += locExpVector[locExp.GetOffset_Elmt_Id(i)]->GetNcoeffs();
            }


            // Set up the mapping for the boundary conditions
            offset = cnt = 0;
            for(i = 0; i < bndCondExp.num_elements(); i++)
            {
                set<int> foundExtraVerts;
                for(j = 0; j < bndCondExp[i]->GetExpSize(); j++)
                {
                    bndSegExp  = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[i]->GetExp(j));

                    cnt = offset + bndCondExp[i]->GetCoeff_Offset(j);
                    for(k = 0; k < 2; k++)
                    {
                        meshVertId = (bndSegExp->GetGeom1D())->GetVid(k);
                        m_bndCondCoeffsToGlobalCoeffsMap[cnt+bndSegExp->GetVertexMap(k)] = graphVertOffset[ReorderedGraphVertId[0][meshVertId]];

                        set<int>::iterator iter = extraDirVerts.find(meshVertId);
                        if (iter != extraDirVerts.end() && 
                            foundExtraVerts.count(meshVertId) == 0)
                        {
                            int loc = bndCondExp[i]->GetCoeff_Offset(j) + 
                                bndSegExp->GetVertexMap(k);
                            int gid = graphVertOffset[
                                ReorderedGraphVertId[0][meshVertId]];
                            m_extraDirDofs[i].push_back(make_pair(loc,gid));
                            foundExtraVerts.insert(meshVertId);
                        }
                    }

                    meshEdgeId = (bndSegExp->GetGeom1D())->GetEid();
                    bndEdgeCnt = 0;
                    for(k = 0; k < bndSegExp->GetNcoeffs(); k++)
                    {
                        if(m_bndCondCoeffsToGlobalCoeffsMap[cnt+k] == -1)
                        {
                            m_bndCondCoeffsToGlobalCoeffsMap[cnt+k] =
                                graphVertOffset[ReorderedGraphVertId[1][meshEdgeId]]+bndEdgeCnt;
                            bndEdgeCnt++;
                        }
                    }
                }
                offset += bndCondExp[i]->GetNcoeffs();
            }

            globalId = Vmath::Vmax(m_numLocalCoeffs,&m_localToGlobalMap[0],1)+1;
            m_numGlobalBndCoeffs = globalId;

            /**
             * STEP 5: The boundary condition mapping is generated from the
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
                 m_solnType == eIterativeMultiLevelStaticCond ||
                 m_solnType == eXxtMultiLevelStaticCond) && nGraphVerts)
            {
                if (m_staticCondLevel < (bottomUpGraph->GetNlevels()-1) &&
                    m_staticCondLevel < m_maxStaticCondLevel)
                {
                    Array<OneD, int> vwgts_perm(
                        Dofs[0].size()+Dofs[1].size()-firstNonDirGraphVertId);
                    for(i = 0; i < locExpVector.size(); ++i)
                    {
                        int eid = locExp.GetOffset_Elmt_Id(i);
                        locExpansion = boost::dynamic_pointer_cast<
                            StdRegions::StdExpansion2D>(locExpVector[eid]);
                        for(j = 0; j < locExpansion->GetNverts(); ++j)
                        {
                            meshEdgeId = (locExpansion->GetGeom2D())->GetEid(j);
                            meshVertId = (locExpansion->GetGeom2D())->GetVid(j);
                            
                            if(ReorderedGraphVertId[0][meshVertId] >= 
                                   firstNonDirGraphVertId)
                            {
                                vwgts_perm[ReorderedGraphVertId[0][meshVertId]-
                                           firstNonDirGraphVertId] = 
                                    Dofs[0][meshVertId];
                            }
                            
                            if(ReorderedGraphVertId[1][meshEdgeId] >= 
                                   firstNonDirGraphVertId)
                            {
                                vwgts_perm[ReorderedGraphVertId[1][meshEdgeId]-
                                           firstNonDirGraphVertId] = 
                                    Dofs[1][meshEdgeId];
                            }
                        }
                    }

                    bottomUpGraph->ExpandGraphWithVertexWeights(vwgts_perm);
                    
                    m_nextLevelLocalToGlobalMap = MemoryManager<AssemblyMap>::
                        AllocateSharedPtr(this,bottomUpGraph);
                }
            }
            
            m_hash = boost::hash_range(
                m_localToGlobalMap.begin(), m_localToGlobalMap.end());
        }

        /**
         * The only unique identifiers of the vertices and edges of the mesh are
         * the vertex id and the mesh id (stored in their corresponding Geometry
         * object).  However, setting up a global numbering based on these id's
         * will not lead to a suitable or optimal numbering. Mainly because:
         *  - we want the Dirichlet DOF's to be listed first
         *  - we want an optimal global numbering of the remaining DOF's
         *    (strategy still need to be defined but can for example be: minimum
         *    bandwith or minimum fill-in of the resulting global system matrix)
         *
         * The vertices and egdes therefore need to be rearranged which is
         * perofrmed in in the following way: The vertices and edges of the mesh
         * are considered as vertices of a graph (in a computer science
         * terminology, equivalently, they can also be considered as boundary
         * degrees of freedom, whereby all boundary modes of a single edge are
         * considered as a single DOF). We then will use different algorithms to
         * reorder the graph-vertices.
         *
         * In the following we use a boost graph object to store this graph the
         * first template parameter (=OutEdgeList) is chosen to be of type
         * std::set. Similarly we also use a std::set to hold the adjacency
         * information. A similar edge might exist multiple times and so to
         * prevent the definition of parallel edges, we use std::set
         * (=boost::setS) rather than std::vector (=boost::vecS).
         *
         * Two different containers are used to store the graph vertex id's of
         * the different mesh vertices and edges. They are implemented as a STL
         * map such that the graph vertex id can later be retrieved by the
         * unique mesh vertex or edge id's which serve as the key of the map.
         *
         * Therefore, the algorithm proceeds as follows:
         */
        int AssemblyMapCG2D::SetUp2DGraphC0ContMap(
                const ExpList  &locExp,
                const Array<OneD, const ExpListSharedPtr> &bndCondExp,
                const Array<OneD, Array<OneD, const SpatialDomains::BoundaryConditionShPtr> >  &bndConditions,
                const vector<map<int,int> >& periodicVerticesId,
                const map<int,int>& periodicEdgesId,
                Array<OneD, map<int,int> > &Dofs,
                Array<OneD, map<int,int> > &ReorderedGraphVertId,
                int          &firstNonDirGraphVertId,
                int          &nExtraDirichlet,
                BottomUpSubStructuredGraphSharedPtr &bottomUpGraph, 
                set<int> &extraDirVerts,
                const bool checkIfSystemSingular,
                int mdswitch, 
                bool doInteriorMap)
        {
            int i,j,k,l,m;
            int cnt = 0;
            int meshVertId, meshVertId2;
            int meshEdgeId, meshEdgeId2;
            int graphVertId = 0;
            StdRegions::StdExpansion2DSharedPtr  locExpansion;
            LocalRegions::SegExpSharedPtr        bndSegExp;
            MultiRegions::ExpList0DSharedPtr     bndVertExp;
            const StdRegions::StdExpansionVector &locExpVector = *(locExp.GetExp());
            map<int,int>::iterator mapIt;
            map<int,int>::const_iterator mapConstIt;
            bool systemSingular = true;
            LibUtilities::CommSharedPtr vCommRow = m_comm->GetRowComm();

            /**
             * STEP 1: Order the Dirichlet vertices and edges first
             */
            for(i = 0; i < bndCondExp.num_elements(); i++)
            {
                // Check to see if any value on edge has Dirichlet value.
                cnt = 0;
                for(k = 0; k < bndConditions.num_elements(); ++k)
                {
                    if(bndConditions[k][i]->GetBoundaryConditionType()==SpatialDomains::eDirichlet)
                    {
                        cnt++;
                    }
                    if(bndConditions[k][i]->GetBoundaryConditionType()!=SpatialDomains::eNeumann)
                    {
                        systemSingular = false;
                    }
                }
                
                // If all boundaries are Dirichlet take out of mask 
                if(cnt == bndConditions.num_elements())
                {
                    for(j = 0; j < bndCondExp[i]->GetNumElmts(); j++)
                    {
                        bndSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[i]->GetExp(j));
                        meshEdgeId = (bndSegExp->GetGeom1D())->GetEid();
                        ReorderedGraphVertId[1][meshEdgeId] = graphVertId++;
                        for(k = 0; k < 2; k++)
                        {
                            meshVertId = (bndSegExp->GetGeom1D())->GetVid(k);
                            if(ReorderedGraphVertId[0].count(meshVertId) == 0)
                            {
                                ReorderedGraphVertId[0][meshVertId] = graphVertId++;
                            }
                        }
                    }
                }
            }


            /**
             * STEP 1.5: Exchange Dirichlet mesh vertices between processes and
             * check for singular problems.
             */
            // Collate information on Dirichlet vertices from all processes
            int n = m_comm->GetSize();
            int p  = m_comm->GetRank();
            Array<OneD, int> counts (n, 0);
            Array<OneD, int> offsets(n, 0);
            counts[p] = ReorderedGraphVertId[0].size();
            vCommRow->AllReduce(counts, LibUtilities::ReduceSum);
            
            for (i = 1; i < n; ++i)
            {
                offsets[i] = offsets[i-1] + counts[i-1];
            }

            int nTot = Vmath::Vsum(n,counts,1);
            Array<OneD, int> vertexlist(nTot, 0);
            std::map<int, int>::iterator it;
            for (it = ReorderedGraphVertId[0].begin(), i = 0;
                 it != ReorderedGraphVertId[0].end();
                 ++it, ++i)
            {
                vertexlist[offsets[p] + i] = it->first;
            }
            vCommRow->AllReduce(vertexlist, LibUtilities::ReduceSum);

            map<int, int> extraDirVertIds;

            // Ensure Dirchlet vertices are consistently recorded between
            // processes (e.g. Dirichlet region meets Neumann region across a
            // partition boundary requires vertex on partition to be Dirichlet).
            for (i = 0; i < n; ++i)
            {
                if (i == p)
                {
                    continue;
                }

                for(j = 0; j < bndCondExp.num_elements(); j++)
                {
                    for(k = 0; k < bndCondExp[j]->GetNumElmts(); k++)
                    {
                        bndSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[j]->GetExp(k));
                        for(l = 0; l < 2; l++)
                        {
                            meshVertId = (bndSegExp->GetGeom1D())->GetVid(l);
                            if(ReorderedGraphVertId[0].count(meshVertId) == 0)
                            {
                                for (m = 0; m < counts[i]; ++m)
                                {
                                    if (vertexlist[offsets[i]+m] == meshVertId)
                                    {
                                        extraDirVertIds[meshVertId] = i;
                                        ReorderedGraphVertId[0][meshVertId] = graphVertId++;
                                        nExtraDirichlet++;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            for (i = 0; i < n; ++i)
            {
                counts [i] = 0;
                offsets[i] = 0;
            }

            counts[p] = extraDirVertIds.size();
            vCommRow->AllReduce(counts, LibUtilities::ReduceSum);
            nTot = Vmath::Vsum(n, counts, 1);
            
            offsets[0] = 0;
            
            for (i = 1; i < n; ++i)
            {
                offsets[i] = offsets[i-1] + counts[i-1];
            }

            Array<OneD, int> vertids  (nTot, 0);
            Array<OneD, int> vertprocs(nTot, 0);
            
            for (it  = extraDirVertIds.begin(), i = 0; 
                 it != extraDirVertIds.end(); ++it, ++i)
            {
                vertids  [offsets[p]+i] = it->first;
                vertprocs[offsets[p]+i] = it->second;
            }

            vCommRow->AllReduce(vertids,   LibUtilities::ReduceSum);
            vCommRow->AllReduce(vertprocs, LibUtilities::ReduceSum);
            
            for (i = 0; i < nTot; ++i)
            {
                if (m_comm->GetRank() != vertprocs[i])
                {
                    continue;
                }
                
                extraDirVerts.insert(vertids[i]);
            }

            // Check between processes if the whole system is singular
            int s = (systemSingular ? 1 : 0);
            vCommRow->AllReduce(s, LibUtilities::ReduceMin);
            systemSingular = (s == 1 ? true : false);

            // Count the number of boundary regions on each process
            Array<OneD, int> bccounts(n, 0);
            bccounts[p] = bndCondExp.num_elements();
            vCommRow->AllReduce(bccounts, LibUtilities::ReduceSum);

            // Find the process rank with the maximum number of boundary regions
            int maxBCIdx = Vmath::Imax(n, bccounts, 1);

            // If the system is singular, the process with the maximum number of
            // BCs will set a Dirichlet vertex to make system non-singular.
            // Note: we find the process with maximum boundary regions to ensure
            // we do not try to set a Dirichlet vertex on a partition with no
            // intersection with the boundary.
            if(systemSingular == true && checkIfSystemSingular && maxBCIdx == p)
            {
                if(m_session->DefinesParameter("SingularElement"))
                {
                    int s_eid;
                    m_session->LoadParameter("SingularElement", s_eid);

                    ASSERTL1(s_eid < locExpVector.size(),"SingularElement Parameter is too large");
                    
                    meshVertId = locExpVector[s_eid]->GetGeom2D()->GetVid(0);
                }
                else
                {
                    //last region i and j=0 edge
                    bndSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[bndCondExp.num_elements()-1]->GetExp(0));
                    
                    //first vertex 0 of the edge
                    meshVertId = (bndSegExp->GetGeom1D())->GetVid(0);
                }

                if(ReorderedGraphVertId[0].count(meshVertId) == 0)
                {
                    ReorderedGraphVertId[0][meshVertId] = graphVertId++;
                }
            }

            firstNonDirGraphVertId = graphVertId;

            /**
             * STEP 2: Now order all other vertices and edges in the graph and
             * create a temporary numbering of domain-interior vertices and
             * edges.
             */
            typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> BoostGraph;
            typedef boost::graph_traits<BoostGraph>::vertex_descriptor BoostVertex;
            BoostGraph boostGraphObj;

            int tempGraphVertId = 0;
            int nVerts;
            int vertCnt;
            int edgeCnt;
            int localOffset=0;
            int nTotalVerts=0;
            map<int, int>    vertTempGraphVertId;
            map<int, int>    edgeTempGraphVertId;
            map<int, int>    intTempGraphVertId;
            Array<OneD, int> localVerts;
            Array<OneD, int> localEdges;
            Array<OneD, int> localinterior;

            m_numLocalBndCoeffs = 0;

            /// - Periodic vertices
            for(k = 0; k < periodicVerticesId.size(); ++k)
            {
                for(mapConstIt = periodicVerticesId[k].begin(); mapConstIt != periodicVerticesId[k].end(); mapConstIt++)
                {
                    meshVertId  = mapConstIt->first;
                    meshVertId2 = mapConstIt->second;

                    if(ReorderedGraphVertId[0].count(meshVertId) == 0)
                    {

                        if(ReorderedGraphVertId[0].count(meshVertId2) == 0)
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
                            ReorderedGraphVertId[0][meshVertId] = ReorderedGraphVertId[0][meshVertId2];
                        }
                    }
                    else
                    {
                        if(ReorderedGraphVertId[0].count(meshVertId2) == 0)
                        {
                            ReorderedGraphVertId[0][meshVertId2] = ReorderedGraphVertId[0][meshVertId];
                        }
                    }
                }
            }

            /// - Periodic edges
            for(mapConstIt = periodicEdgesId.begin(); mapConstIt != periodicEdgesId.end(); mapConstIt++)
            {
                meshEdgeId  = mapConstIt->first;
                meshEdgeId2 = abs(mapConstIt->second);

                if(meshEdgeId < meshEdgeId2)
                {
                    ASSERTL0(ReorderedGraphVertId[1].count(meshEdgeId) == 0,
                             "This periodic boundary edge has been specified before");
                    ASSERTL0(ReorderedGraphVertId[1].count(meshEdgeId2) == 0,
                             "This periodic boundary edge has been specified before");

                    edgeTempGraphVertId[meshEdgeId]  = tempGraphVertId;
                    edgeTempGraphVertId[meshEdgeId2] = tempGraphVertId++;
                }
            }

            /// - All other vertices and edges
            int nEdgeCoeffs;
            int elmtid;
            for(i = 0; i < locExpVector.size(); ++i)
            {
                elmtid = locExp.GetOffset_Elmt_Id(i);
                if((locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(
                        locExpVector[elmtid])))
                {
                    m_numLocalBndCoeffs += locExpansion->NumBndryCoeffs();

                    nTotalVerts += locExpansion->GetNverts();

                }
            }

            // Store the temporary graph vertex
            // id's of all element edges and
            // vertices in these 2 arrays below
            localVerts = Array<OneD, int>(nTotalVerts,-1);
            localEdges = Array<OneD, int>(nTotalVerts,-1);

            for(i = 0; i < locExpVector.size(); ++i)
            {
                elmtid = locExp.GetOffset_Elmt_Id(i);
                if((locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(
                        locExpVector[elmtid])))
                {
                    vertCnt = 0;
                    nVerts = locExpansion->GetNverts();
                    for(j = 0; j < nVerts; ++j)
                    {
                        meshVertId = (locExpansion->GetGeom2D())->GetVid(j);
                        if(ReorderedGraphVertId[0].count(meshVertId) == 0)
                        {
                            // non-periodic & non-Dirichlet vertex 
                            if(vertTempGraphVertId.count(meshVertId) == 0)
                            {
                                boost::add_vertex(boostGraphObj);
                                vertTempGraphVertId[meshVertId] = tempGraphVertId++;
                            }
                            localVerts[localOffset + vertCnt++] = vertTempGraphVertId[meshVertId];
                        }
                    }
                }
                localOffset+=nVerts;
            }

            m_numNonDirVertexModes=tempGraphVertId;

            localOffset=0;
            for(i = 0; i < locExpVector.size(); ++i)
            {
                elmtid = locExp.GetOffset_Elmt_Id(i);
                if((locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(
                        locExpVector[elmtid])))
                {
                    edgeCnt = 0;
                    nVerts = locExpansion->GetNverts();
                    for(j = 0; j < nVerts; ++j)
                    {
                        meshEdgeId = (locExpansion->GetGeom2D())->GetEid(j);
                        if(ReorderedGraphVertId[1].count(meshEdgeId) == 0)
                        {
                            // non-periodic & non-Dirichlet edge
                            if(edgeTempGraphVertId.count(meshEdgeId) == 0)
                            {
                                boost::add_vertex(boostGraphObj);
                                edgeTempGraphVertId[meshEdgeId] = tempGraphVertId++;
                            }
                            localEdges[localOffset + edgeCnt++] = edgeTempGraphVertId[meshEdgeId];
                        }
                    }
                }
                localOffset+=nVerts;
            }

            if(doInteriorMap)
            {
                for(i = 0; i < locExpVector.size(); ++i)
                {
                    elmtid = locExp.GetOffset_Elmt_Id(i);
                    if((locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(
                            locExpVector[elmtid])))
                    {

                        boost::add_vertex(boostGraphObj);
                        intTempGraphVertId[elmtid] = tempGraphVertId++;
                    }
                }
            }

            localOffset=0;
            for(i = 0; i < locExpVector.size(); ++i)
            {
                elmtid = locExp.GetOffset_Elmt_Id(i);
                if((locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(
                        locExpVector[elmtid])))
                {
                    nVerts = locExpansion->GetNverts();
                    // Now loop over all local edges and vertices of this
                    // element and define that all other edges and vertices of
                    // this element are adjacent to them.
                    for(j = 0; j < nVerts; j++)
                    {
                        if(localVerts[j+localOffset]==-1)
                        {
                            break;
                        }
                        
                        for(k = 0; k < nVerts; k++)
                        {
                            if(localVerts[k+localOffset]==-1)
                            {
                                break;
                            }
                            if(k!=j)
                            {
                                boost::add_edge( (size_t) localVerts[j+localOffset], 
                                                 (size_t) localVerts[k+localOffset],boostGraphObj);
                            }
                        }

                        for(k = 0; k < nVerts; k++)
                        {
                            if(localEdges[k+localOffset]==-1)
                            {
                                break;
                            }
                            boost::add_edge( (size_t) localVerts[j+localOffset], 
                                             (size_t) localEdges[k+localOffset],boostGraphObj);
                        }

                        if(doInteriorMap)
                        {
                            boost::add_edge( (size_t)  localVerts[j+localOffset], 
                                             (size_t) intTempGraphVertId[elmtid],boostGraphObj);
                        }
                    }

                    for(j = 0; j < nVerts; j++)
                    {
                        if(localEdges[j+localOffset]==-1)
                        {
                            break;
                        }
                        for(k = 0; k < nVerts; k++)
                        {
                            if(localEdges[k+localOffset]==-1)
                            {
                                break;
                            }
                            if(k!=j)
                            {
                                boost::add_edge( (size_t) localEdges[j+localOffset], 
                                                 (size_t) localEdges[k+localOffset],boostGraphObj);
                            }
                        }
                        for(k = 0; k < nVerts; k++)
                        {
                            if(localVerts[k+localOffset]==-1)
                            {
                                break;
                            }
                            boost::add_edge( (size_t) localEdges[j+localOffset], 
                                             (size_t) localVerts[k+localOffset],boostGraphObj);
                        }

                        if(doInteriorMap)
                        {
                            boost::add_edge( (size_t) localEdges[j+localOffset],  
                                             (size_t) intTempGraphVertId[elmtid], boostGraphObj);
                        }
                    }
                    
                    if(doInteriorMap)
                    {
                        for(j = 0; j < nVerts; j++)
                        {
                            if(localVerts[j+localOffset]==-1)
                            {
                                break;
                            }
                            boost::add_edge( (size_t) intTempGraphVertId[elmtid], 
                                             (size_t) localVerts[j+localOffset], boostGraphObj);
                        }
                        
                        for(j = 0; j < nVerts; j++)
                        {
                            if(localEdges[j+localOffset]==-1)
                            {
                                break;
                            }

                            boost::add_edge( (size_t) intTempGraphVertId[elmtid], 
                                             (size_t) localEdges[j+localOffset], boostGraphObj);
                        }
                    }


                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a local 2D expansion failed");
                }
                localOffset+=nVerts;
            }

            // Container to store vertices of the graph which correspond to
            // degrees of freedom along the boundary.
            set<int> partVerts;
            
            if (m_solnType == eIterativeMultiLevelStaticCond)
            {
                vector<long> procVerts,  procEdges;
                set   <int>  foundVerts, foundEdges;
                
                // Loop over element and construct the procVerts and procEdges
                // vectors, which store the geometry IDs of mesh vertices and
                // edges respectively which are local to this process.
                for(i = cnt = 0; i < locExpVector.size(); ++i)
                {
                    elmtid = locExp.GetOffset_Elmt_Id(i);
                    if((locExpansion = boost::dynamic_pointer_cast<
                            StdRegions::StdExpansion2D>(locExpVector[elmtid])))
                    {
                        for (j = 0; j < locExpansion->GetNverts(); ++j, ++cnt)
                        {
                            int vid = locExpansion->GetGeom2D()->GetVid(j)+1;
                            int eid = locExpansion->GetGeom2D()->GetEid(j)+1;
                        
                            if (foundVerts.count(vid) == 0)
                            {
                                procVerts.push_back(vid);
                                foundVerts.insert(vid);
                            }
                        
                            if (foundEdges.count(eid) == 0)
                            {
                                procEdges.push_back(eid);
                                foundEdges.insert(eid);
                            }
                        }
                    }
                    else
                    {
                        ASSERTL0(false,
                                 "dynamic cast to a local 2D expansion failed");
                    }
                }

                int unique_verts = foundVerts.size();
                int unique_edges = foundEdges.size();

                // Now construct temporary GS objects. These will be used to
                // populate the arrays tmp3 and tmp4 with the multiplicity of
                // the vertices and edges respectively to identify those
                // vertices and edges which are located on partition boundary.
                Array<OneD, long> vertArray(unique_verts, &procVerts[0]);
                Array<OneD, long> edgeArray(unique_edges, &procEdges[0]);
                Gs::gs_data *tmp1 = Gs::Init(vertArray, m_comm);
                Gs::gs_data *tmp2 = Gs::Init(edgeArray, m_comm);
                Array<OneD, NekDouble> tmp3(unique_verts, 1.0);
                Array<OneD, NekDouble> tmp4(unique_edges, 1.0);
                Gs::Gather(tmp3, Gs::gs_add, tmp1);
                Gs::Gather(tmp4, Gs::gs_add, tmp2);

                // Finally, fill the partVerts set with all non-Dirichlet
                // vertices which lie on a partition boundary.
                for (i = 0; i < unique_verts; ++i)
                {
                    if (tmp3[i] > 1.0)
                    {
                        if (ReorderedGraphVertId[0].count(procVerts[i]-1) == 0)
                        {
                            partVerts.insert(vertTempGraphVertId[procVerts[i]-1]);
                        }
                    }
                }
            
                for (i = 0; i < unique_edges; ++i)
                {
                    if (tmp4[i] > 1.0)
                    {
                        if (ReorderedGraphVertId[1].count(procEdges[i]-1) == 0)
                        {
                            partVerts.insert(edgeTempGraphVertId[procEdges[i]-1]);
                        }
                    }
                }
            }

            /**
             * STEP 3: Reorder graph for optimisation.
             */
            int nGraphVerts = tempGraphVertId;
            Array<OneD, int> perm(nGraphVerts);
            Array<OneD, int> iperm(nGraphVerts);

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
                        MultiLevelBisectionReordering(boostGraphObj,perm,iperm,bottomUpGraph,mdswitch,partVerts);
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
                vCommRow->AllReduce(m_lowestStaticCondLevel, 
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
                ReorderedGraphVertId[0][mapIt->first] = iperm[mapIt->second] + graphVertId;
            }
            for(mapIt = edgeTempGraphVertId.begin(); mapIt != edgeTempGraphVertId.end(); mapIt++)
            {
                ReorderedGraphVertId[1][mapIt->first] = iperm[mapIt->second] + graphVertId;
            }

            if(doInteriorMap)
            {
                for(mapIt = intTempGraphVertId.begin(); mapIt != intTempGraphVertId.end(); mapIt++)
                {
                    ReorderedGraphVertId[2][mapIt->first] = iperm[mapIt->second] + graphVertId;
                }
            }
            
            return nGraphVerts;
        }
    }
}
