///////////////////////////////////////////////////////////////////////////////
//
// File LocToGlobalC0ContMap.cpp
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
// Description: C0-continuous Local to Global mapping routines
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/LocalToGlobalC0ContMap.h>

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
         * @class LocalToGlobalC0ContMap
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
        LocalToGlobalC0ContMap::LocalToGlobalC0ContMap(
                const LibUtilities::CommSharedPtr &pComm):
            LocalToGlobalBaseMap(pComm)
        {
        }


        /**
         *
         */
        LocalToGlobalC0ContMap::LocalToGlobalC0ContMap(
                const LibUtilities::CommSharedPtr &pComm,
                const int numLocalCoeffs,
                const ExpList &locExp,
                const GlobalSysSolnType solnType):
            LocalToGlobalBaseMap(pComm)
        {
            switch(locExp.GetExp(0)->GetShapeDimension())
            {
            case 1:
                {
                    SetUp1DExpansionC0ContMap(numLocalCoeffs, locExp, solnType);
                }
                break;
            case 2:
                {
                    SetUp2DExpansionC0ContMap(numLocalCoeffs, locExp, solnType);
                }
                break;
            case 3:
                {
                    SetUp3DExpansionC0ContMap(numLocalCoeffs, locExp,  solnType);
                }
                break;
            default:
                {
                    ASSERTL0(false,"Local To Global map not defined for this dimension");
                }
            }

            Nektar::Array<OneD, long> tmp(m_globalToUniversalMap.num_elements());
            for (unsigned int i = 0; i < m_globalToUniversalMap.num_elements(); ++i)
            {
                tmp[i] = m_globalToUniversalMap[i];
            }
            m_gsh = Gs::Init(tmp, pComm);
            Gs::Unique(tmp, pComm);
            for (unsigned int i = 0; i < m_globalToUniversalMap.num_elements(); ++i)
            {
                m_globalToUniversalMapUnique[i] = (tmp[i] >= 0 ? 1 : 0);
            }

            CalculateBndSystemBandWidth();
            CalculateFullSystemBandWidth();
        }


        /**
         *
         */
        LocalToGlobalC0ContMap::LocalToGlobalC0ContMap(
                const LibUtilities::CommSharedPtr &pComm,
                const int numLocalCoeffs,
                const ExpList &locExp,
                const GlobalSysSolnType solnType,
                const Array<OneD, const ExpListSharedPtr>
                                                            &bndCondExp,
                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
                                                            &bndConditions,
                const map<int,int>& periodicVerticesId):
            LocalToGlobalBaseMap(pComm)
        {
            SetUp1DExpansionC0ContMap(numLocalCoeffs,
                                      locExp,
                                      solnType,
                                      bndCondExp,
                                      bndConditions,
                                      periodicVerticesId);

            Nektar::Array<OneD, long> tmp(m_globalToUniversalMap.num_elements());
            for (unsigned int i = 0; i < m_globalToUniversalMap.num_elements(); ++i)
            {
                tmp[i] = m_globalToUniversalMap[i];
            }
            m_gsh = Gs::Init(tmp, pComm);
            Gs::Unique(tmp, pComm);
            for (unsigned int i = 0; i < m_globalToUniversalMap.num_elements(); ++i)
            {
                m_globalToUniversalMapUnique[i] = (tmp[i] >= 0 ? 1 : 0);
            }

            CalculateBndSystemBandWidth();
            CalculateFullSystemBandWidth();
        }


        /**
         *
         */
        LocalToGlobalC0ContMap::LocalToGlobalC0ContMap(
                const LibUtilities::CommSharedPtr &pComm,
                const int numLocalCoeffs,
                const ExpList &locExp,
                const GlobalSysSolnType solnType,
                const Array<OneD, const ExpListSharedPtr> &bndCondExp,
                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
                                                            &bndConditions,
                const vector<map<int,int> >& periodicVerticesId,
                const map<int,int>& periodicEdgesId,
                const bool checkIfSystemSingular) :
            LocalToGlobalBaseMap(pComm)
        {
            SetUp2DExpansionC0ContMap(numLocalCoeffs,
                                      locExp,
                                      solnType,
                                      bndCondExp,
                                      bndConditions,
                                      periodicVerticesId,
                                      periodicEdgesId,
                                      checkIfSystemSingular);


            Nektar::Array<OneD, long> tmp(m_globalToUniversalMap.num_elements());
            for (unsigned int i = 0; i < m_globalToUniversalMap.num_elements(); ++i)
            {
                tmp[i] = m_globalToUniversalMap[i];
            }
            m_gsh = Gs::Init(tmp, pComm);
            Gs::Unique(tmp, pComm);
            for (unsigned int i = 0; i < m_globalToUniversalMap.num_elements(); ++i)
            {
                m_globalToUniversalMapUnique[i] = (tmp[i] >= 0 ? 1 : 0);
            }


            CalculateBndSystemBandWidth();
            CalculateFullSystemBandWidth();
        }


        /**
         *
         */
        LocalToGlobalC0ContMap::LocalToGlobalC0ContMap(
                const LibUtilities::CommSharedPtr &pComm,
                const int numLocalCoeffs,
                const ExpList &locExp,
                const GlobalSysSolnType solnType,
                const Array<OneD, const ExpListSharedPtr> &bndCondExp,
                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
                                                                &bndConditions,
                const map<int,int>& periodicVerticesId,
                const map<int,int>& periodicEdgesId,
                const map<int,int>& periodicFacesId):
            LocalToGlobalBaseMap(pComm)
        {
            SetUp3DExpansionC0ContMap(numLocalCoeffs,
                                      locExp,
                                      solnType,
                                      bndCondExp,
                                      bndConditions,
                                      periodicVerticesId,
                                      periodicEdgesId,
                                      periodicFacesId);

            Nektar::Array<OneD, long> tmp(m_globalToUniversalMap.num_elements());
            for (unsigned int i = 0; i < m_globalToUniversalMap.num_elements(); ++i)
            {
                tmp[i] = m_globalToUniversalMap[i];
            }
            m_gsh = Gs::Init(tmp, pComm);
            Gs::Unique(tmp, pComm);
            for (unsigned int i = 0; i < m_globalToUniversalMap.num_elements(); ++i)
            {
                m_globalToUniversalMapUnique[i] = (tmp[i] >= 0 ? 1 : 0);
            }

            CalculateBndSystemBandWidth();
            CalculateFullSystemBandWidth();
        }


        /**
         *
         */
        LocalToGlobalC0ContMap::~LocalToGlobalC0ContMap()
        {
        }


        /**
         * Construction of the local->global map is achieved in several stages.
         * A mesh vertex renumbering is constructed as follows
         *  - Domain Dirichlet boundaries are numbered first.
         *  - Domain periodic boundaries are numbered next and given
         *    consistent global indices.
         *  - Remaining vertices are ordered last.
         * This mesh vertex renumbering is then used to generate the first
         * part of the local to global mapping of the degrees of freedom. The
         * remainder of the map consists of the element-interior degrees of
         * freedom. This leads to the block-diagonal submatrix as each
         * element-interior mode is globally orthogonal to modes in all other
         * elements.
         *
         * The boundary condition mapping is generated from the same vertex
         * renumbering.
         */
        void LocalToGlobalC0ContMap::SetUp1DExpansionC0ContMap(
                const int numLocalCoeffs,
                const ExpList &locExp,
                const GlobalSysSolnType solnType,
                const Array<OneD, const MultiRegions::ExpListSharedPtr>
                                                                &bndCondExp,
                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
                                                                &bndConditions,
                const map<int,int>& periodicVerticesId)
        {
            int i,j;
            int cnt = 0;
            int meshVertId;
            int meshVertId2;
            int graphVertId = 0;
            int globalId;

            const StdRegions::StdExpansionVector &locExpVector = *(locExp.GetExp());

            LocalRegions::SegExpSharedPtr locSegExp;

            int nbnd = bndCondExp.num_elements();

            m_staticCondLevel                = 0;
            m_solnType                       = solnType;
            m_signChange                     = false;
            m_numPatches                     = locExpVector.size();
            m_numLocalCoeffs                 = numLocalCoeffs;
            m_numLocalBndCoeffs              = 2*locExpVector.size();
            m_localToGlobalMap               = Array<OneD, int>(m_numLocalCoeffs,-1);
            m_localToGlobalBndMap            = Array<OneD, int>(m_numLocalBndCoeffs,-1);
            m_bndCondCoeffsToGlobalCoeffsMap = Array<OneD, int>(nbnd);
            m_numLocalBndCoeffsPerPatch      = Array<OneD, unsigned int>(m_numPatches);
            m_numLocalIntCoeffsPerPatch      = Array<OneD, unsigned int>(m_numPatches);
            for(i = 0; i < m_numPatches; ++i)
            {
                m_numLocalBndCoeffsPerPatch[i] = (unsigned int) locExpVector[locExp.GetOffset_Elmt_Id(i)]->NumBndryCoeffs();
                m_numLocalIntCoeffsPerPatch[i] = (unsigned int) locExpVector[locExp.GetOffset_Elmt_Id(i)]->GetNcoeffs() - locExpVector[locExp.GetOffset_Elmt_Id(i)]->NumBndryCoeffs();
            }

            // The only unique identifiers of the vertices of the mesh
            // are the vertex id (stored in their corresponding
            // Geometry object).  However, setting up a global
            // numbering based on these id's will not lead to a
            // suitable or optimal numbering. Mainly because: - we
            // want the Dirichlet DOF's to be listed first - we want
            // an optimal global numbering of the remaining DOF's
            // (strategy still need to be defined but can for example
            // be: minimum bandwith or minimum fill-in of the
            // resulting global system matrix)
            //
            // That's why the vertices be rearranged. Currently, this
            // is done in the following way: The vertices of the mesh
            // are considered as vertices of a graph. (We then will
            // use algorithms to reorder these graph-vertices -
            // although not implemented yet in 1D).
            //
            // A container is used to store the graph vertex id's of
            // the different mesh vertices. This is implemented as a
            // STL map such that the graph vertex id can later be
            // retrieved by the unique mesh vertex id's which serve as
            // the key of the map.
            map<int, int> vertReorderedGraphVertId;
            map<int,int>::const_iterator mapConstIt;

            // STEP 1: Order the Dirichlet vertices first
            m_numGlobalDirBndCoeffs = 0;
            for(i = 0; i < nbnd; i++)
            {
                if(bndConditions[i]->GetBoundaryConditionType()==SpatialDomains::eDirichlet)
                {
                    meshVertId = ((bndCondExp[i])->GetVertex())->GetVid();
                    vertReorderedGraphVertId[meshVertId] = graphVertId++;
                    m_numGlobalDirBndCoeffs++;
                    m_numLocalDirBndCoeffs++;
                }
            }

            // STEP 2: Order the periodic vertices next
            // This allows to give corresponding DOF's the same
            // global ID
            for(mapConstIt = periodicVerticesId.begin(); mapConstIt != periodicVerticesId.end(); mapConstIt++)
            {
                meshVertId  = mapConstIt->first;
                meshVertId2 = mapConstIt->second;

                ASSERTL0(vertReorderedGraphVertId.count(meshVertId) == 0,
                         "This periodic boundary vertex has been specified before");
                ASSERTL0(vertReorderedGraphVertId.count(meshVertId) == 0,
                         "This periodic boundary vertex has been specified before");

                vertReorderedGraphVertId[meshVertId]  = graphVertId;
                vertReorderedGraphVertId[meshVertId2] = graphVertId++;
            }


            // STEP 3: List the other vertices
            // STEP 4: Set up simple map based on vertex and edge id's
            for(i = 0; i < locExpVector.size(); ++i)
            {
                cnt = locExp.GetCoeff_Offset(i);
                if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(locExpVector[i]))
                {
                    for(j = 0; j < locSegExp->GetNverts(); ++j)
                    {
                        meshVertId = (locSegExp->GetGeom1D())->GetVid(j);
                        if(vertReorderedGraphVertId.count(meshVertId) == 0)
                        {
                            vertReorderedGraphVertId[meshVertId] = graphVertId++;
                        }

                        m_localToGlobalMap[cnt+locSegExp->GetVertexMap(j)] =
                            vertReorderedGraphVertId[meshVertId];
                    }
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a segment expansion failed");
                }
            }

            // Set up the mapping for the boundary conditions
            for(i = 0; i < nbnd; i++)
            {
                meshVertId = ((bndCondExp[i])->GetVertex())->GetVid();
                m_bndCondCoeffsToGlobalCoeffsMap[i] = vertReorderedGraphVertId[meshVertId];
            }

            // Setup interior mapping and the boundary map
            globalId = graphVertId;
            m_numGlobalBndCoeffs = globalId;

            cnt=0;
            for(i = 0; i < m_numLocalCoeffs; ++i)
            {
                if(m_localToGlobalMap[i] == -1)
                {
                    m_localToGlobalMap[i] = globalId++;
                }
                else
                {
                    m_localToGlobalBndMap[cnt++]=m_localToGlobalMap[i];
                }
            }
            m_numGlobalCoeffs = globalId;

            SetUpUniversalC0ContMap(locExp);
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
        void LocalToGlobalC0ContMap::SetUp2DExpansionC0ContMap(
                const int numLocalCoeffs,
                const ExpList &locExp,
                const GlobalSysSolnType solnType,
                const Array<OneD, const ExpListSharedPtr> &bndCondExp,
                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
                                                                &bndConditions,
                const vector<map<int,int> >& periodicVerticesId,
                const map<int,int>& periodicEdgesId,
                const bool checkIfSystemSingular)
        {
            int i,j,k;
            int cnt = 0,offset=0;
            int meshVertId, meshVertId2;
            int meshEdgeId, meshEdgeId2;
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
            StdRegions::EdgeOrientation         edgeOrient;
            Array<OneD, unsigned int>           edgeInteriorMap;
            Array<OneD, int>                    edgeInteriorSign;
            const StdRegions::StdExpansionVector &locExpVector = *(locExp.GetExp());
            m_signChange = false;
            m_systemSingular = false;
            Array<OneD, map<int,int> > ReorderedGraphVertId(2);
            Array<OneD, map<int,int> > Dofs(2);
            BottomUpSubStructuredGraphSharedPtr bottomUpGraph;

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
            nGraphVerts = SetUp2DGraphC0ContMap(locExp,solnType,
                                                bndCondExp,bndConditionsVec,
                                                periodicVerticesId,periodicEdgesId,
                                                Dofs,
                                                ReorderedGraphVertId,
                                                firstNonDirGraphVertId,
                                                nExtraDirichlet,
                                                bottomUpGraph,
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
                    if( bndCondExp[i]->GetExp(j)==bndCondExp[(bndCondExp.num_elements()-1)]->GetExp(bndCondExp[i]->GetExpSize()-1)
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
            m_solnType = solnType;
            m_staticCondLevel = 0;
            m_numPatches =  locExpVector.size();
            m_numLocalBndCoeffsPerPatch = Array<OneD, unsigned int>(m_numPatches);
            m_numLocalIntCoeffsPerPatch = Array<OneD, unsigned int>(m_numPatches);
            for(i = 0; i < m_numPatches; ++i)
            {
                m_numLocalBndCoeffsPerPatch[i] = (unsigned int) locExpVector[locExp.GetOffset_Elmt_Id(i)]->NumBndryCoeffs();
                m_numLocalIntCoeffsPerPatch[i] = (unsigned int) locExpVector[locExp.GetOffset_Elmt_Id(i)]->GetNcoeffs() - m_numLocalBndCoeffsPerPatch[i];
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
                    nEdgeInteriorCoeffs = locExpansion->GetEdgeNcoeffs(j)-2;
                    edgeOrient          = (locExpansion->GetGeom2D())->GetEorient(j);
                    meshEdgeId          = (locExpansion->GetGeom2D())->GetEid(j);
                    meshVertId          = (locExpansion->GetGeom2D())->GetVid(j);

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
                for(j = 0; j < bndCondExp[i]->GetExpSize(); j++)
                {
                    bndSegExp  = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[i]->GetExp(j));

                    cnt = offset + bndCondExp[i]->GetCoeff_Offset(j);
                    for(k = 0; k < 2; k++)
                    {
                        meshVertId = (bndSegExp->GetGeom1D())->GetVid(k);
                        m_bndCondCoeffsToGlobalCoeffsMap[cnt+bndSegExp->GetVertexMap(k)] = graphVertOffset[ReorderedGraphVertId[0][meshVertId]];
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

            // Set up the local to global map for the next level when using
            // multi-level static condensation
            if( (solnType == eDirectMultiLevelStaticCond || solnType == eIterativeMultiLevelStaticCond) && nGraphVerts )
            {
                if(m_staticCondLevel < (bottomUpGraph->GetNlevels()-1))
                {
                    Array<OneD, int> vwgts_perm(Dofs[0].size()+Dofs[1].size()+-firstNonDirGraphVertId);
                    for(i = 0; i < Dofs[0].size(); ++i)
                    {
                        if(ReorderedGraphVertId[0][i] >= firstNonDirGraphVertId)
                        {
                            vwgts_perm[ReorderedGraphVertId[0][i]-firstNonDirGraphVertId] = Dofs[0][i];
                        }
                    }
                        
                    for(i = 0; i < Dofs[1].size(); ++i)
                    {
                        if(ReorderedGraphVertId[1][i] >= firstNonDirGraphVertId)
                        {
                            vwgts_perm[ReorderedGraphVertId[1][i]-firstNonDirGraphVertId] = Dofs[1][i];
                        }
                    }

                    bottomUpGraph->ExpandGraphWithVertexWeights(vwgts_perm);
                    
                    m_nextLevelLocalToGlobalMap = MemoryManager<LocalToGlobalBaseMap>::
                        AllocateSharedPtr(this,bottomUpGraph);
                }
            }

            SetUpUniversalC0ContMap(locExp);
        }


        /**
         * Sets up the global to universal mapping of degrees of freedom across
         * processors.
         */
        void LocalToGlobalC0ContMap::SetUpUniversalC0ContMap(
                const ExpList &locExp)
        {
            StdRegions::StdExpansionSharedPtr locExpansion;
            int nDim = 0;
            int nVert = 0;
            int nEdge = 0;
            int nFace = 0;
            int maxEdgeDof = 0;
            int maxFaceDof = 0;
            int dof = 0;
            int cnt;
            int i,j,k;
            int meshVertId;
            int meshEdgeId;
            int meshFaceId;
            int vGlobalId;
            int maxBndGlobalId = 0;
            StdRegions::EdgeOrientation         edgeOrient;
            StdRegions::FaceOrientation         faceOrient;
            Array<OneD, unsigned int>           edgeInteriorMap;
            Array<OneD, int>                    edgeInteriorSign;
            Array<OneD, unsigned int>           faceInteriorMap;
            Array<OneD, int>                    faceInteriorSign;

            const StdRegions::StdExpansionVector &locExpVector = *(locExp.GetExp());

            m_globalToUniversalMap = Nektar::Array<OneD, int>(m_numGlobalCoeffs, -1);
            m_globalToUniversalMapUnique = Nektar::Array<OneD, int>(m_numGlobalCoeffs, -1);

            // Loop over all the elements in the domain to gather mesh data
            for(i = 0; i < locExpVector.size(); ++i)
            {
                locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion>(locExpVector[i]);
                nVert += locExpansion->GetNverts();
                nEdge += locExpansion->GetNedges();
                nFace += locExpansion->GetNfaces();
                // Loop over all edges (and vertices) of element i
                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    dof = locExpansion->GetEdgeNcoeffs(j)-2;
                    maxEdgeDof = (dof > maxEdgeDof ? dof : maxEdgeDof);
                }
                for(j = 0; j < locExpansion->GetNfaces(); ++j)
                {
                    dof = locExpansion->GetFaceIntNcoeffs(j);
                    maxFaceDof = (dof > maxFaceDof ? dof : maxFaceDof);
                }
            }

            // Tell other processes about how many dof we have
            m_comm->AllReduce(nVert, LibUtilities::ReduceSum);
            m_comm->AllReduce(nEdge, LibUtilities::ReduceSum);
            m_comm->AllReduce(nFace, LibUtilities::ReduceSum);
            m_comm->AllReduce(maxEdgeDof, LibUtilities::ReduceMax);
            m_comm->AllReduce(maxFaceDof, LibUtilities::ReduceMax);

            // Assemble global to universal mapping for this process
            for(i = 0; i < locExpVector.size(); ++i)
            {
                locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion>(locExpVector[i]);
                nDim = locExpansion->GetShapeDimension();
                cnt = locExp.GetCoeff_Offset(i);

                // Loop over all vertices of element i
                for(j = 0; j < locExpansion->GetNverts(); ++j)
                {
                    meshVertId   = (locExpansion->GetGeom())->GetVid(j);
                    vGlobalId    = m_localToGlobalMap[cnt+locExpansion->GetVertexMap(j)];
                    m_globalToUniversalMap[vGlobalId] = meshVertId + 1;
                    maxBndGlobalId = (vGlobalId > maxBndGlobalId ? vGlobalId : maxBndGlobalId);
                }

                // Loop over all edges of element i
                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    locExpansion->GetEdgeInteriorMap(j,edgeOrient,edgeInteriorMap,edgeInteriorSign);
                    dof = locExpansion->GetEdgeNcoeffs(j)-2;
                    if (nDim == 2)
                    {
                        meshEdgeId   = (locExpansion->GetGeom2D())->GetEid(j);
                    }
                    else
                    {
                        meshEdgeId   = (locExpansion->GetGeom3D())->GetEid(j);
                    }

                    // Set the global DOF's for the interior modes of edge j
                    for(k = 0; k < dof; ++k)
                    {
                        vGlobalId = m_localToGlobalMap[cnt+edgeInteriorMap[k]];
                        m_globalToUniversalMap[vGlobalId]
                           = nVert + meshEdgeId * maxEdgeDof + k + 1;
                        maxBndGlobalId = (vGlobalId > maxBndGlobalId ? vGlobalId : maxBndGlobalId);
                    }
                }

                // Loop over all faces of element i
                for(j = 0; j < locExpansion->GetNfaces(); ++j)
                {
                    faceOrient          = (locExpansion->GetGeom3D())->GetFaceorient(j);

                    locExpansion->GetFaceInteriorMap(j,faceOrient,faceInteriorMap,faceInteriorSign);
                    dof = locExpansion->GetFaceIntNcoeffs(j);
                    meshFaceId = (locExpansion->GetGeom3D())->GetFid(j);

                    for(k = 0; k < dof; ++k)
                    {
                        vGlobalId = m_localToGlobalMap[cnt+faceInteriorMap[k]];
                        m_globalToUniversalMap[vGlobalId]
                           = nVert + nEdge*maxEdgeDof + meshFaceId * maxFaceDof
                                   + k + 1;
                        maxBndGlobalId = (vGlobalId > maxBndGlobalId ? vGlobalId : maxBndGlobalId);
                    }
                }

            }

            // Finally, internal DOF do not participate in any data
            // exchange, so we set these to the special GSLib id=0 so
            // they are ignored.
            for (k = maxBndGlobalId + 1; k < m_globalToUniversalMap.num_elements(); ++k)
            {
                m_globalToUniversalMap[k] = 0;
            }
        }


        /**
         * The only unique identifiers of the vertices and edges of the mesh
         * are the vertex id and the mesh id (stored in their corresponding
         * Geometry object).  However, setting up a global numbering based on
         * these id's will not lead to a suitable or optimal numbering. Mainly
         * because:
         *  - we want the Dirichlet DOF's to be listed first
         *  - we want an optimal global numbering of the remaining DOF's
         *    (strategy still need to be defined but can for example be:
         *    minimum bandwith or minimum fill-in of the resulting global
         *    system matrix)
         *
         * The vertices and egdes therefore need to be rearranged
         * which is perofrmed in in the following way: The vertices
         * and edges of the mesh are considered as vertices of a graph
         * (in a computer science terminology, equivalently, they can
         * also be considered as boundary degrees of freedom, whereby
         * all boundary modes of a single edge are considered as a
         * single DOF). We then will use different algorithms to
         * reorder the graph-vertices.
         *
         * In the following we use a boost graph object to store this
         * graph the first template parameter (=OutEdgeList) is chosen
         * to be of type std::set. Similarly we also use a std::set to
         * hold the adjacency information. A similar edge might exist
         * multiple times and so to prevent the definition of parallel
         * edges, we use std::set (=boost::setS) rather than
         * std::vector (=boost::vecS).
         *
         * Two different containers are used to store the graph vertex id's of
         * the different mesh vertices and edges. They are implemented as a STL
         * map such that the graph vertex id can later be retrieved by the
         * unique mesh vertex or edge id's which serve as the key of the map.
         *
         * Therefore, the algorithm proceeds as follows:
         */

        int LocalToGlobalC0ContMap::SetUp2DGraphC0ContMap(
                const ExpList  &locExp,
                const GlobalSysSolnType solnType,
                const Array<OneD, const ExpListSharedPtr> &bndCondExp,
                const Array<OneD, Array<OneD, const SpatialDomains::BoundaryConditionShPtr> >  &bndConditions,
                const vector<map<int,int> >& periodicVerticesId,
                const map<int,int>& periodicEdgesId,
                Array<OneD, map<int,int> > &Dofs,
                Array<OneD, map<int,int> > &ReorderedGraphVertId,
                int          &firstNonDirGraphVertId,
                int          &nExtraDirichlet,
                BottomUpSubStructuredGraphSharedPtr &bottomUpGraph, 
                const bool checkIfSystemSingular,
                int mdswitch, 
                bool doInteriorMap)
        {
            int i,j,k,l,m;
            int cnt = 0;
            int meshVertId, meshVertId2;
            int meshEdgeId, meshEdgeId2;
            int graphVertId = 0;
            StdRegions::StdExpansion2DSharedPtr locExpansion;
            LocalRegions::SegExpSharedPtr       bndSegExp;
            MultiRegions::ExpList0DSharedPtr     bndVertExp;
            const StdRegions::StdExpansionVector &locExpVector = *(locExp.GetExp());
            map<int,int>::iterator mapIt;
            map<int,int>::const_iterator mapConstIt;
            bool systemSingular = true;
            
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
                    for(j = 0; j < bndCondExp[i]->GetExpSize(); j++)
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
             * STEP 1.5: Exchange Dirichlet mesh vertices between processes.
             */
            int n = m_comm->GetSize();
            int p  = m_comm->GetRank();
            Array<OneD, int> counts (n, 0);
            Array<OneD, int> offsets(n, 0);
            counts[p] = ReorderedGraphVertId[0].size();
            m_comm->AllReduce(counts, LibUtilities::ReduceSum);
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
            m_comm->AllReduce(vertexlist, LibUtilities::ReduceSum);

            for (i = 0; i < n; ++i)
            {
                if (i == p)
                {
                    continue;
                }

                for(j = 0; j < bndCondExp.num_elements(); j++)
                {
                    for(k = 0; k < bndCondExp[j]->GetExpSize(); k++)
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
                                        ReorderedGraphVertId[0][meshVertId] = graphVertId++;
                                        nExtraDirichlet++;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            int s = (systemSingular ? 1 : 0);
            m_comm->AllReduce(s, LibUtilities::ReduceMin);
            systemSingular = (s == 1 ? true : false);
            if(systemSingular == true && checkIfSystemSingular)
            {
                //last region i and j=0 edge
                bndSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[bndCondExp.num_elements()-1]->GetExp(0));

                //first vertex 0 of the edge
                meshVertId = (bndSegExp->GetGeom1D())->GetVid(0);
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
            map<int, int>    vertTempGraphVertId;
            map<int, int>    edgeTempGraphVertId;
            map<int, int>    intTempGraphVertId;
            map<int, int>    vwgts_map;
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
                meshEdgeId2 = mapConstIt->second;

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
                if(locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(locExpVector[elmtid]))
                {
                    m_numLocalBndCoeffs += locExpansion->NumBndryCoeffs();

                    nVerts = locExpansion->GetNverts();
                    // For element i, store the temporary graph
                    // vertex id's of all element edges and
                    // verices in these 2 arrays below
                    localVerts = Array<OneD, int>(nVerts,-1);
                    localEdges = Array<OneD, int>(nVerts,-1);
                    vertCnt = 0;
                    edgeCnt = 0;

                    for(j = 0; j < nVerts; ++j)
                    {
                        nEdgeCoeffs = locExpansion->GetEdgeNcoeffs(j);

                        meshVertId = (locExpansion->GetGeom2D())->GetVid(j);
                        if(ReorderedGraphVertId[0].count(meshVertId) == 0)
                        {
                            // non-periodic & non-Dirichlet vertex 
                            if(vertTempGraphVertId.count(meshVertId) == 0)
                            {
                                boost::add_vertex(boostGraphObj);
                                vertTempGraphVertId[meshVertId] = tempGraphVertId++;
                            }
                            localVerts[vertCnt++] = vertTempGraphVertId[meshVertId];
                            vwgts_map[ vertTempGraphVertId[meshVertId] ] = Dofs[0][meshVertId];
                        }

                        meshEdgeId = (locExpansion->GetGeom2D())->GetEid(j);
                        if(ReorderedGraphVertId[1].count(meshEdgeId) == 0)
                        {
                            // non-periodic & non-Dirichlet edge
                            if(edgeTempGraphVertId.count(meshEdgeId) == 0)
                            {
                                boost::add_vertex(boostGraphObj);
                                edgeTempGraphVertId[meshEdgeId] = tempGraphVertId++;
                            }
                            localEdges[edgeCnt++] = edgeTempGraphVertId[meshEdgeId];
                            vwgts_map[ edgeTempGraphVertId[meshEdgeId] ] = Dofs[1][meshEdgeId];
                        }
                    }

                    if(doInteriorMap)
                    {
                        boost::add_vertex(boostGraphObj);
                        intTempGraphVertId[elmtid] = tempGraphVertId++;
                        vwgts_map[ intTempGraphVertId[elmtid] ] = Dofs[2][elmtid];
                    }

                    // Now loop over all local edges and vertices
                    // of this element and define that all other
                    // edges and vertices of this element are
                    // adjacent to them.
                    for(j = 0; j < nVerts; j++)
                    {
                        if(localVerts[j]==-1)
                        {
                            break;
                        }
                        
                        for(k = 0; k < nVerts; k++)
                        {
                            if(localVerts[k]==-1)
                            {
                                break;
                            }
                            if(k!=j)
                            {
                                boost::add_edge( (size_t) localVerts[j], (size_t) localVerts[k],boostGraphObj);
                            }
                        }

                        for(k = 0; k < nVerts; k++)
                        {
                            if(localEdges[k]==-1)
                            {
                                break;
                            }
                            boost::add_edge( (size_t) localVerts[j], (size_t) localEdges[k],boostGraphObj);
                        }

                        if(doInteriorMap)
                        {
                            boost::add_edge( (size_t)  localVerts[j], (size_t) intTempGraphVertId[elmtid],boostGraphObj);
                        }
                    }
                    for(j = 0; j < nVerts; j++)
                    {
                        if(localEdges[j]==-1)
                        {
                            break;
                        }
                        for(k = 0; k < nVerts; k++)
                        {
                            if(localEdges[k]==-1)
                            {
                                break;
                            }
                            if(k!=j)
                            {
                                boost::add_edge( (size_t) localEdges[j], (size_t) localEdges[k],boostGraphObj);
                            }
                        }
                        for(k = 0; k < nVerts; k++)
                        {
                            if(localVerts[k]==-1)
                            {
                                break;
                            }
                            boost::add_edge( (size_t) localEdges[j], (size_t) localVerts[k],boostGraphObj);
                        }

                        if(doInteriorMap)
                        {
                            boost::add_edge( (size_t) localEdges[j],  (size_t) intTempGraphVertId[elmtid], boostGraphObj);
                        }
                    }
                    
                    if(doInteriorMap)
                    {
                        for(j = 0; j < nVerts; j++)
                        {
                            if(localVerts[j]==-1)
                            {
                                break;
                            }
                            boost::add_edge( (size_t) intTempGraphVertId[elmtid], (size_t) localVerts[j], boostGraphObj);
                        }
                        
                        for(j = 0; j < nVerts; j++)
                        {
                            if(localEdges[j]==-1)
                            {
                                break;
                            }

                            boost::add_edge( (size_t) intTempGraphVertId[elmtid], (size_t) localEdges[j], boostGraphObj);
                        }
                    }


                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a local 2D expansion failed");
                }
            }


            /**
             * STEP 3: Reorder graph for optimisation.
             */
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
                switch(solnType)
                {
                case eDirectFullMatrix:
                case eIterativeFull:
                    {
                        NoReordering(boostGraphObj,perm,iperm);
                    }
                    break;
                case eDirectStaticCond:
                case eIterativeStaticCond:
                    {
                        CuthillMckeeReordering(boostGraphObj,perm,iperm);
                    }
                    break;
                case eDirectMultiLevelStaticCond:
                case eIterativeMultiLevelStaticCond:
                    {
                        MultiLevelBisectionReordering(boostGraphObj,vwgts,perm,iperm,bottomUpGraph, mdswitch);
                    }
                    break;
                default:
                    {
                        ASSERTL0(false,"Unrecognised solution type");
                    }
                }
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
        void LocalToGlobalC0ContMap::SetUp3DExpansionC0ContMap(const int numLocalCoeffs,
                                                               const ExpList &locExp,
                                                               const GlobalSysSolnType solnType,
                                                               const Array<OneD, const ExpListSharedPtr> &bndCondExp,
                                                               const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions,
                                                               const map<int,int>& periodicVerticesId,
                                                               const map<int,int>& periodicEdgesId,
                                                               const map<int,int>& periodicFacesId)
        {
            int i,j,k,l;
            int cnt = 0,cnt1=0;
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
            StdRegions::EdgeOrientation         edgeOrient;
            StdRegions::FaceOrientation         faceOrient;
            Array<OneD, unsigned int>           edgeInteriorMap;
            Array<OneD, int>                    edgeInteriorSign;
            Array<OneD, unsigned int>           faceInteriorMap;
            Array<OneD, int>                    faceInteriorSign;

            const StdRegions::StdExpansionVector &locExpVector = *(locExp.GetExp());

            m_signChange = false;

            map<int,int> vertReorderedGraphVertId;
            map<int,int> edgeReorderedGraphVertId;
            map<int,int> faceReorderedGraphVertId;
            map<int,int>::iterator mapIt;
            map<int,int>::const_iterator mapConstIt;

            /**
             * STEP 1: Order the Dirichlet vertices and edges first
             */
            for(i = 0; i < bndCondExp.num_elements(); i++)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); j++)
                {
                    bndCondFaceExp = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(bndCondExp[i]->GetExp(j));
                    if(bndConditions[i]->GetBoundaryConditionType()==SpatialDomains::eDirichlet)
                    {
                        meshFaceId = (bndCondFaceExp->GetGeom2D())->GetFid();
                        faceReorderedGraphVertId[meshFaceId] = graphVertId++;
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
                        nLocDirBndCondDofs += bndCondFaceExp->GetNcoeffs();
                    }
                    nLocBndCondDofs += bndCondFaceExp->GetNcoeffs();
                }
            }
            m_numLocalDirBndCoeffs = nLocDirBndCondDofs;
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
            int nEdges;
            int nFaces;
            int vertCnt;
            int edgeCnt;
            int faceCnt;
            map<int, int>          vertTempGraphVertId;
            map<int, int>          edgeTempGraphVertId;
            map<int, int>          faceTempGraphVertId;
            map<int, int>          vwgts_map;
            Array<OneD, int>       localVerts;
            Array<OneD, int>       localEdges;
            Array<OneD, int>       localFaces;

            m_numLocalBndCoeffs = 0;

            /// - Periodic vertices
            for(mapConstIt = periodicVerticesId.begin(); mapConstIt != periodicVerticesId.end(); mapConstIt++)
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
            for(mapConstIt = periodicEdgesId.begin(); mapConstIt != periodicEdgesId.end(); mapConstIt++)
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
                    else
                    {
                        ASSERTL0(edgeReorderedGraphVertId[meshEdgeId2] == edgeReorderedGraphVertId[meshEdgeId],
                                 "These values should be equal");
                    }
                }
            }

            /// - Periodic faces
            for(mapConstIt = periodicFacesId.begin(); mapConstIt != periodicFacesId.end(); mapConstIt++)
            {
                meshFaceId  = mapConstIt->first;
                meshFaceId2 = mapConstIt->second;

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
                if(locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion3D>(locExpVector[locExp.GetOffset_Elmt_Id(i)]))
                {
                	nVerts = locExpansion->GetNverts();
                	nEdges = locExpansion->GetNedges();
                	nFaces = locExpansion->GetNfaces();
                    // For element i, store the temporary graph
                    // vertex id's of all element edges and
                    // verices in these 2 arrays below
                    localVerts = Array<OneD, int>(nVerts,-1);
                    localEdges = Array<OneD, int>(nEdges,-1);
                    localFaces = Array<OneD, int>(nFaces,-1);
                    vertCnt = 0;
                    edgeCnt = 0;
                    faceCnt = 0;

                    for(j = 0; j < nVerts; ++j)
                    {
                        meshVertId = (locExpansion->GetGeom3D())->GetVid(j);
                        if(vertReorderedGraphVertId.count(meshVertId) == 0)
                        {
                        	if(vertTempGraphVertId.count(meshVertId) == 0)
                        	{
                        		boost::add_vertex(boostGraphObj);
                        		vertTempGraphVertId[meshVertId] = tempGraphVertId++;
                        	}
                        	localVerts[vertCnt++] = vertTempGraphVertId[meshVertId];
                        	vwgts_map[ vertTempGraphVertId[meshVertId] ] = 1;
                        }
                    }
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
                        	}
                        	localEdges[edgeCnt++] = edgeTempGraphVertId[meshEdgeId];
                        	vwgts_map[ edgeTempGraphVertId[meshEdgeId] ] = nEdgeInteriorCoeffs;
                        }
                    }
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
                        	}
                        	localFaces[faceCnt++] = faceTempGraphVertId[meshFaceId];
                        	vwgts_map[ faceTempGraphVertId[meshFaceId] ] = nFaceInteriorCoeffs;
                        }
                    }
                    m_numLocalBndCoeffs += locExpansion->NumBndryCoeffs();

                    // Now loop over all local faces, edges and vertices
                    // of this element and define that all other
                    // faces, edges and verices of this element are
                    // adjacent to them.

                    // Vertices
                    for(j = 0; j < nVerts; j++)
                    {
                        if(localVerts[j]==-1)
                        {
                            break;
                        }
                        // associate to other vertices
                        for(k = 0; k < nVerts; k++)
                        {
                            if(localVerts[k]==-1)
                            {
                                break;
                            }
                            if(k!=j)
                            {
                                boost::add_edge( (size_t) localVerts[j], (size_t) localVerts[k],boostGraphObj);
                            }
                        }
                        // associate to other edges
                        for(k = 0; k < nEdges; k++)
                        {
                            if(localEdges[k]==-1)
                            {
                                break;
                            }
                            boost::add_edge( (size_t) localVerts[j], (size_t) localEdges[k],boostGraphObj);
                        }
                        // associate to other faces
                        for(k = 0; k < nFaces; k++)
                        {
                            if(localFaces[k]==-1)
                            {
                                break;
                            }
                            boost::add_edge( (size_t) localVerts[j], (size_t) localFaces[k],boostGraphObj);
                        }
                    }

                    // Edges
                    for(j = 0; j < nEdges; j++)
                    {
                        if(localEdges[j]==-1)
                        {
                            break;
                        }
                        // Associate to other edges
                        for(k = 0; k < nEdges; k++)
                        {
                            if(localEdges[k]==-1)
                            {
                                break;
                            }
                            if(k!=j)
                            {
                                boost::add_edge( (size_t) localEdges[j], (size_t) localEdges[k],boostGraphObj);
                            }
                        }
                        // Associate to vertices
                        for(k = 0; k < nVerts; k++)
                        {
                            if(localVerts[k]==-1)
                            {
                                break;
                            }
                            boost::add_edge( (size_t) localEdges[j], (size_t) localVerts[k],boostGraphObj);
                        }
                        // Associate to faces
                        for(k = 0; k < nFaces; k++)
                        {
                            if(localFaces[k]==-1)
                            {
                                break;
                            }
                            boost::add_edge( (size_t) localEdges[j], (size_t) localFaces[k],boostGraphObj);
                        }
                    }

                    // Faces
                    for(j = 0; j < nFaces; j++)
                    {
                        if(localFaces[j]==-1)
                        {
                            break;
                        }
                        // Associate to other faces
                        for(k = 0; k < nFaces; k++)
                        {
                            if(localFaces[k]==-1)
                            {
                                break;
                            }
                            if(k!=j)
                            {
                                boost::add_edge( (size_t) localFaces[j], (size_t) localFaces[k],boostGraphObj);
                            }
                        }
                        // Associate to vertices
                        for(k = 0; k < nVerts; k++)
                        {
                            if(localVerts[k]==-1)
                            {
                                break;
                            }
                            boost::add_edge( (size_t) localFaces[j], (size_t) localVerts[k],boostGraphObj);
                        }
                        // Associate to edges
                        for(k = 0; k < nEdges; k++)
                        {
                            if(localEdges[k]==-1)
                            {
                                break;
                            }
                            boost::add_edge( (size_t) localFaces[j], (size_t) localEdges[k],boostGraphObj);
                        }
                    }
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a local 3D expansion failed");
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
                switch(solnType)
                {
                case eDirectFullMatrix:
                case eIterativeFull:
                    {
                        NoReordering(boostGraphObj,perm,iperm);
                    }
                    break;
                case eDirectStaticCond:
                case eIterativeStaticCond:
                    {
                        CuthillMckeeReordering(boostGraphObj,perm,iperm);
                    }
                    break;
                case eDirectMultiLevelStaticCond:
                case eIterativeMultiLevelStaticCond:
                    {
                        MultiLevelBisectionReordering(boostGraphObj,vwgts,perm,iperm,bottomUpGraph);
                    }
                    break;
                default:
                    {
                        ASSERTL0(false,"Unrecognised solution type");
                    }
                }
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

            m_solnType = solnType;
            m_staticCondLevel = 0;
            m_numPatches =  locExpVector.size();
            m_numLocalBndCoeffsPerPatch = Array<OneD, unsigned int>(m_numPatches);
            m_numLocalIntCoeffsPerPatch = Array<OneD, unsigned int>(m_numPatches);
            for(i = 0; i < m_numPatches; ++i)
            {
                m_numLocalBndCoeffsPerPatch[i] = (unsigned int) locExpVector[locExp.GetOffset_Elmt_Id(i)]->NumBndryCoeffs();
                m_numLocalIntCoeffsPerPatch[i] = (unsigned int) locExpVector[locExp.GetOffset_Elmt_Id(i)]->GetNcoeffs() - locExpVector[locExp.GetOffset_Elmt_Id(i)]->NumBndryCoeffs();
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
                    nFaceInteriorCoeffs = locExpansion->GetFaceIntNcoeffs(j);
                    faceOrient          = (locExpansion->GetGeom3D())->GetFaceorient(j);
                    meshFaceId          = (locExpansion->GetGeom3D())->GetFid(j);

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
                for(j = 0; j < bndCondExp[i]->GetExpSize(); j++)
                {
                    bndCondFaceExp  = boost::dynamic_pointer_cast<StdRegions::StdExpansion2D>(bndCondExp[i]->GetExp(j));
                    cnt = offset + bndCondExp[i]->GetCoeff_Offset(j);
                    for(k = 0; k < bndCondFaceExp->GetNverts(); k++)
                    {
                        meshVertId = (bndCondFaceExp->GetGeom2D())->GetVid(k);
                        m_bndCondCoeffsToGlobalCoeffsMap[cnt+bndCondFaceExp->GetVertexMap(k)] = graphVertOffset[vertReorderedGraphVertId[meshVertId]];
                    }

                    for(k = 0; k < bndCondFaceExp->GetNedges(); k++)
                    {
                        nEdgeInteriorCoeffs = bndCondFaceExp->GetEdgeNcoeffs(k)-2;
                        edgeOrient          = (bndCondFaceExp->GetGeom2D())->GetEorient(k);
                        meshEdgeId          = (bndCondFaceExp->GetGeom2D())->GetEid(k);

                        bndCondFaceExp->GetEdgeInteriorMap(k,edgeOrient,edgeInteriorMap,edgeInteriorSign);

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
            if( (solnType == eDirectMultiLevelStaticCond || solnType == eIterativeMultiLevelStaticCond) && nGraphVerts )
            {
                if(m_staticCondLevel < (bottomUpGraph->GetNlevels()-1))
                {
                    Array<OneD, int> vwgts_perm(nGraphVerts);
                    for(i = 0; i < nGraphVerts; ++i)
                    {
                        vwgts_perm[i] = vwgts[perm[i]];
                    }
                    
                    bottomUpGraph->ExpandGraphWithVertexWeights(vwgts_perm);
                    m_nextLevelLocalToGlobalMap = MemoryManager<LocalToGlobalBaseMap>::AllocateSharedPtr(this,bottomUpGraph);
                }
            }
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
        void LocalToGlobalC0ContMap::CalculateFullSystemBandWidth()
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


        int LocalToGlobalC0ContMap::v_GetLocalToGlobalMap(const int i) const
        {
            return m_localToGlobalMap[i];
        }

        int LocalToGlobalC0ContMap::v_GetGlobalToUniversalMap(const int i) const
        {
            return m_globalToUniversalMap[i];
        }

        int LocalToGlobalC0ContMap::v_GetGlobalToUniversalMapUnique(const int i) const
        {
            return m_globalToUniversalMapUnique[i];
        }

        const Array<OneD,const int>&
                    LocalToGlobalC0ContMap::v_GetLocalToGlobalMap(void)
        {
            return m_localToGlobalMap;
        }

        void LocalToGlobalC0ContMap::v_SetLocalToGlobalMap(const Array<OneD, int>& inarray)
        {
            m_localToGlobalMap = inarray;
        }

        const Array<OneD,const int>&
                    LocalToGlobalC0ContMap::v_GetGlobalToUniversalMap(void)
        {
            return m_globalToUniversalMap;
        }

        const Array<OneD,const int>&
                    LocalToGlobalC0ContMap::v_GetGlobalToUniversalMapUnique(void)
        {
            return m_globalToUniversalMapUnique;
        }

        NekDouble LocalToGlobalC0ContMap::v_GetLocalToGlobalSign(
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

        const Array<OneD, NekDouble>& LocalToGlobalC0ContMap::v_GetLocalToGlobalSign() const
        {
            if (m_signChange)
            {
                return m_localToGlobalSign;
            }
            else
            {
                return NullNekDouble1DArray;
            }
        }

        void LocalToGlobalC0ContMap::v_SetLocalToGlobalSign(const Array<OneD, NekDouble>& inarray)
        {
            m_localToGlobalSign = inarray;
        }

        const void LocalToGlobalC0ContMap::v_LocalToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                          Array<OneD,       NekDouble>& global) const
        {
            if(m_signChange)
            {
                Vmath::Scatr(m_numLocalCoeffs, m_localToGlobalSign.get(), loc.get(), m_localToGlobalMap.get(), global.get());
            }
            else
            {
                Vmath::Scatr(m_numLocalCoeffs, loc.get(), m_localToGlobalMap.get(), global.get());
            }
        }

        const void LocalToGlobalC0ContMap::v_LocalToGlobal(
                    const NekVector<const NekDouble>& loc,
                          NekVector<      NekDouble>& global) const
        {
            LocalToGlobal(loc.GetPtr(),global.GetPtr());
        }

        const void LocalToGlobalC0ContMap::v_GlobalToLocal(
                    const Array<OneD, const NekDouble>& global,
                          Array<OneD,       NekDouble>& loc) const
        {
            if(m_signChange)
            {
                Vmath::Gathr(m_numLocalCoeffs, m_localToGlobalSign.get(), global.get(), m_localToGlobalMap.get(), loc.get());
            }
            else
            {
                Vmath::Gathr(m_numLocalCoeffs, global.get(), m_localToGlobalMap.get(), loc.get());
            }
        }

        const void LocalToGlobalC0ContMap::v_GlobalToLocal(
                    const NekVector<const NekDouble>& global,
                          NekVector<      NekDouble>& loc) const
        {
            GlobalToLocal(global.GetPtr(),loc.GetPtr());
        }

        const void LocalToGlobalC0ContMap::v_Assemble(
                    const Array<OneD, const NekDouble> &loc,
                          Array<OneD,       NekDouble> &global) const
        {
            ASSERTL1(loc.get() != global.get(),"Local and Global Arrays cannot be the same");

            Vmath::Zero(m_numGlobalCoeffs, global.get(), 1);

            if(m_signChange)
            {
                Vmath::Assmb(m_numLocalCoeffs, m_localToGlobalSign.get(), loc.get(), m_localToGlobalMap.get(), global.get());
            }
            else
            {
                Vmath::Assmb(m_numLocalCoeffs, loc.get(), m_localToGlobalMap.get(), global.get());
            }
            UniversalAssemble(global);
        }

        const void LocalToGlobalC0ContMap::v_Assemble(
                    const NekVector<const NekDouble>& loc,
                          NekVector<      NekDouble>& global) const
        {
            Assemble(loc.GetPtr(),global.GetPtr());
        }

        const void LocalToGlobalC0ContMap::v_UniversalAssemble(
                      Array<OneD,     NekDouble>& pGlobal) const
        {
            Gs::Gather(pGlobal, Gs::gs_add, m_gsh);
        }

        const void LocalToGlobalC0ContMap::v_UniversalAssemble(
                      NekVector<      NekDouble>& pGlobal) const
        {
            UniversalAssemble(pGlobal.GetPtr());
        }

        const int LocalToGlobalC0ContMap::v_GetFullSystemBandWidth() const
        {
            return m_fullSystemBandWidth;
        }

    }
}


/**
 * $Log: LocalToGlobalC0ContMap.cpp,v $
 * Revision 1.14  2010/03/01 17:57:28  cantwell
 * Fixed 3D global matrix operations.
 * Fixed ProjectCont{1,2,3}D demos.
 * Fixed incorrectly placed ASSERT in boundary conditions.
 * Updated TimingGeneralMatrixOp3D to use contfield3d rather than explist3d.
 *
 * Revision 1.13  2009/11/02 11:19:44  pvos
 * Fixed a bug for reordering a graph without edges
 *
 * Revision 1.12  2009/10/30 14:02:55  pvos
 * Multi-level static condensation updates
 *
 * Revision 1.11  2009/05/10 23:17:12  sherwin
 * Updated mainly to handle doubly periodic meshes which required modification to vertex handling from a numbering perspective
 *
 * Revision 1.10  2009/04/27 21:34:58  sherwin
 * Modified WriteToField Method
 *
 * Revision 1.9  2009/04/20 16:14:06  sherwin
 * Updates for optimising bandwidth of DG solver and allowing write import on explist
 *
 * Revision 1.8  2009/04/03 15:06:35  sherwin
 * Added LinearAdvection, Generalised GenGlobalLinSys and added in bandwidth routine for DG solver
 *
 * Revision 1.7  2009/04/02 13:06:42  sherwin
 * Modified to take symmetric banded system for HDH solver
 *
 * Revision 1.6  2009/02/27 15:47:13  sherwin
 * Corrected bug in initialising just Neumann boundary conditions
 *
 * Revision 1.5  2008/12/19 15:12:39  pvos
 * Updates for precomputed dirichlet forcing functionality
 *
 * Revision 1.4  2008/11/05 16:15:24  pvos
 * Added bandwith calculation routines
 *
 * Revision 1.3  2008/09/23 18:21:00  pvos
 * Updates for working ProjectContField3D demo
 *
 * Revision 1.2  2008/09/17 13:46:40  pvos
 * Added LocalToGlobalC0ContMap for 3D expansions
 *
 * Revision 1.1  2008/09/16 13:36:06  pvos
 * Restructured the LocalToGlobalMap classes
 *
 **/
