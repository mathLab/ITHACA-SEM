///////////////////////////////////////////////////////////////////////////////
//
// File CoupledLcoalToGlobalC0ContMap.cpp
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
// Description: Wrapper class around the library
// LocalToGlobalC0ContMap class for use in the Couplied Linearised NS
// solver.
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/CoupledLocalToGlobalC0ContMap.h>
#include <SpatialDomains/MeshGraph.h>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/Expansion2D.h>
#include <MultiRegions/GlobalLinSysDirectStaticCond.h>

using namespace std;

namespace Nektar
{
    /**
     * This is an vector extension of
     * MultiRegions::AssemblyMapCG::SetUp2DExpansionC0ContMap related to the
     * Linearised Navier Stokes problem
     */
    CoupledLocalToGlobalC0ContMap::CoupledLocalToGlobalC0ContMap(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &graph,
        const SpatialDomains::BoundaryConditionsSharedPtr &boundaryConditions,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const MultiRegions::ExpListSharedPtr &pressure,
        const int nz_loc,
        const bool CheckforSingularSys):
        AssemblyMapCG(pSession)
    {
        int i,j,k,n;
        int cnt = 0,offset=0;
        int meshVertId;
        int meshEdgeId;
        int bndEdgeCnt;
        int globalId;
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
        int nvel = fields.size();

        const LocalRegions::ExpansionVector &locExpVector = *(fields[0]->GetExp());
        int id, diff;
        int nel = fields[0]->GetNumElmts();

        MultiRegions::PeriodicMap periodicVerts;
        MultiRegions::PeriodicMap periodicEdges;
        MultiRegions::PeriodicMap periodicFaces;
        vector<map<int,int> > ReorderedGraphVertId(3);
        MultiRegions::BottomUpSubStructuredGraphSharedPtr bottomUpGraph;
        int staticCondLevel = 0;

        if(CheckforSingularSys) //all singularity checking by setting flag to true
        {
            m_systemSingular = true;
        }
        else  // Turn off singular checking by setting flag to false
        {
            m_systemSingular = false;
        }

        /**
         * STEP 1: Wrap boundary conditions vector in an array
         * (since routine is set up for multiple fields) and call
         * the graph re-odering subroutine to obtain the reordered
         * values
         */

        // Obtain any periodic information and allocate default mapping array
        fields[0]->GetPeriodicEntities(periodicVerts,periodicEdges,periodicFaces);


        const Array<OneD, const MultiRegions::ExpListSharedPtr> bndCondExp = fields[0]->GetBndCondExpansions();

        Array<OneD, Array<OneD, const SpatialDomains::BoundaryConditionShPtr> > bndConditionsVec(nvel);
        for(i = 0; i < nvel; ++i)
        {
            bndConditionsVec[i] = fields[i]->GetBndConditions();
        }

        map<int,int> IsDirVertDof;
        map<int,int> IsDirEdgeDof;

        SpatialDomains::Geometry1DSharedPtr g;
        for(j = 0; j < bndCondExp.size(); ++j)
        {
            map<int,int> BndExpVids;
            // collect unique list of vertex ids for this expansion
            for(k = 0; k < bndCondExp[j]->GetNumElmts(); ++k)
            {
                g = bndCondExp[j]->GetExp(k)->as<LocalRegions::Expansion1D>()
                                            ->GetGeom1D();
                BndExpVids[g->GetVid(0)] = g->GetVid(0);
                BndExpVids[g->GetVid(1)] = g->GetVid(1);
            }

            for(i = 0; i < nvel; ++i)
            {
                if(bndConditionsVec[i][j]->GetBoundaryConditionType()==SpatialDomains::eDirichlet)
                {
                    // set number of Dirichlet conditions along edge
                    for(k = 0; k < bndCondExp[j]->GetNumElmts(); ++k)
                    {
                        IsDirEdgeDof[bndCondExp[j]->GetExp(k)
                                     ->as<LocalRegions::Expansion1D>()
                                     ->GetGeom1D()->GetGlobalID()] += 1;
                    }


                    // Set number of Dirichlet conditions at vertices
                    // with a clamp on its maximum value being nvel to
                    // handle corners between expansions
                    for(auto &mapIt : BndExpVids)
                    {
                        id = IsDirVertDof[mapIt.second]+1;
                        IsDirVertDof[mapIt.second] = (id > nvel)?nvel:id;
                    }
                }
                else
                {
                    // Check to see that edge normals have non-zero
                    // component in this direction since otherwise
                    // also can be singular.
                    /// @TODO: Fix this so that we can extract normals from edges
                    for(k = 0; k < bndCondExp[j]->GetNumElmts(); ++k)
                    {
                        Array<OneD,Array<OneD,NekDouble> > locnorm;
                        LocalRegions::Expansion1DSharedPtr loc_exp
                            = bndCondExp[j]->GetExp(k)
                                ->as<LocalRegions::Expansion1D>();
                        locnorm = loc_exp->GetLeftAdjacentElementExp()->GetEdgeNormal(loc_exp->GetLeftAdjacentElementEdge());
                        //locnorm = bndCondExp[j]->GetExp(k)->Get GetMetricInfo()->GetNormal();

                        int ndir = locnorm.size();
                        if(i < ndir) // account for Fourier version where n can be larger then ndir
                        {
                            for(int l = 0; l < locnorm[0].size(); ++l)
                            {
                                if(fabs(locnorm[i][l]) > NekConstants::kNekZeroTol)
                                {
                                    m_systemSingular = false;
                                    break;
                                }
                            }
                        }
                        if(m_systemSingular == false)
                        {
                            break;
                        }
                    }
                }
            }
        }

        Array<OneD, map<int,int> >Dofs(2);

        Array<OneD, int> AddMeanPressureToEdgeId(nel,-1);
        int edgeId,vertId;


        // special case of singular problem - need to fix one pressure
        // dof to a dirichlet edge. Since we attached pressure dof to
        // last velocity component of edge need to make sure this
        // component is Dirichlet
        if(m_systemSingular)
        {
            id  = -1;
            for(i = 0; i < bndConditionsVec[0].size(); ++i)
            {
                if(bndConditionsVec[nvel-1][i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    id = bndCondExp[i]->GetExp(0)
                                ->as<LocalRegions::Expansion1D>()->GetGeom1D()
                                ->GetGlobalID();
                    break;
                }
            }

            ASSERTL0(id != -1," Did not find an edge to attach singular pressure degree of freedom");

            // determine element with this edge id. There may be a
            // more direct way of getting element from spatialDomains
            for(i = 0; i < nel; ++i)
            {
                for(j = 0; j < locExpVector[i]->GetNverts(); ++j)
                {
                    edgeId = (locExpVector[i]->as<LocalRegions::Expansion2D>()
                                             ->GetGeom2D())->GetEid(j);

                    if(edgeId == id)
                    {
                        AddMeanPressureToEdgeId[i] = id;
                        break;
                    }
                }

                if(AddMeanPressureToEdgeId[i] != -1)
                {
                    break;
                }
            }
        }


        for(i = 0; i < nel; ++i)
        {
            for(j = 0; j < locExpVector[i]->GetNverts(); ++j)
            {
                vertId = (locExpVector[i]->as<LocalRegions::Expansion2D>()
                                           ->GetGeom2D())->GetVid(j);
                if(Dofs[0].count(vertId) == 0)
                {
                    Dofs[0][vertId] = nvel*nz_loc;

                    // Adjust for a Dirichlet boundary condition to give number to be solved
                    if(IsDirVertDof.count(vertId) != 0)
                    {
                        Dofs[0][vertId] -= IsDirVertDof[vertId]*nz_loc;
                    }
                }

                edgeId = (locExpVector[i]->as<LocalRegions::Expansion2D>()
                                           ->GetGeom2D())->GetEid(j);
                if(Dofs[1].count(edgeId) == 0)
                {
                    Dofs[1][edgeId] = nvel*(locExpVector[i]->GetEdgeNcoeffs(j)-2)*nz_loc;
                }

                // Adjust for Dirichlet boundary conditions to give number to be solved
                if(IsDirEdgeDof.count(edgeId) != 0)
                {
                    Dofs[1][edgeId] -= IsDirEdgeDof[edgeId]*nz_loc*(locExpVector[i]->GetEdgeNcoeffs(j)-2);
                }
            }
        }

        set<int> extraDirVerts, extraDirEdges;

        CreateGraph(*fields[0], bndCondExp, bndConditionsVec, false,
                    periodicVerts, periodicEdges, periodicFaces,
                    ReorderedGraphVertId, bottomUpGraph, extraDirVerts,
                    extraDirEdges, firstNonDirGraphVertId, nExtraDirichlet, 4);
        /*
        SetUp2DGraphC0ContMap(*fields[0],
                              bndCondExp,
                              bndConditionsVec,
                              periodicVerts,          periodicEdges,
                              Dofs,                   ReorderedGraphVertId,
                              firstNonDirGraphVertId, nExtraDirichlet,
                              bottomUpGraph, extraDir,  false,  4);
        */

        /**
         * STEP 2a: Set the mean pressure modes to edges depending on
         * type of direct solver technique;
         */

        // determine which edge to add mean pressure dof based on
        // ensuring that at least one pressure dof from an internal
        // patch is associated with its boundary system
        if(m_session->MatchSolverInfoAsEnum("GlobalSysSoln", MultiRegions::eDirectMultiLevelStaticCond))
        {


			FindEdgeIdToAddMeanPressure(ReorderedGraphVertId,
				 nel, locExpVector,
				 edgeId, vertId, firstNonDirGraphVertId, IsDirEdgeDof,
				 bottomUpGraph,
				 AddMeanPressureToEdgeId);
		}

        // Set unset elmts to non-Dirichlet edges.
        // special case of singular problem - need to fix one
        // pressure dof to a dirichlet edge
        for(i = 0; i < nel; ++i)
        {
            for(j = 0; j < locExpVector[i]->GetNverts(); ++j)
            {
                edgeId = (locExpVector[i]->as<LocalRegions::Expansion2D>()
                                           ->GetGeom2D())->GetEid(j);

                if(IsDirEdgeDof.count(edgeId) == 0) // interior edge
                {
                    // setup AddMeanPressureToEdgeId to decide where to
                    // put pressure
                    if(AddMeanPressureToEdgeId[i] == -1)
                    {
                        AddMeanPressureToEdgeId[i] = edgeId;
                    }
                }
            }
            ASSERTL0((AddMeanPressureToEdgeId[i] != -1),"Did not determine "
                     "an edge to attach mean pressure dof");
            // Add the mean pressure degree of freedom to this edge
            Dofs[1][AddMeanPressureToEdgeId[i]] += nz_loc;
        }

        map<int,int> pressureEdgeOffset;

        /**
         * STEP 2: Count out the number of Dirichlet vertices and edges first
         */
        for(i = 0; i < bndCondExp.size(); i++)
        {
            for(j = 0; j < bndCondExp[i]->GetNumElmts(); j++)
            {
                bndSegExp = bndCondExp[i]->GetExp(j)
                                         ->as<LocalRegions::SegExp>();
                for(k = 0; k < nvel; ++k)
                {
                    if(bndConditionsVec[k][i]->GetBoundaryConditionType()==SpatialDomains::eDirichlet)
                    {
                        nLocDirBndCondDofs += bndSegExp->GetNcoeffs()*nz_loc;
                    }

                    if (bndConditionsVec[k][i]->GetBoundaryConditionType()!=SpatialDomains::ePeriodic)
                    {
                        nLocBndCondDofs += bndSegExp->GetNcoeffs()*nz_loc;
                    }
                }
            }
        }

        if(m_systemSingular)
        {
            m_numLocalDirBndCoeffs = nLocDirBndCondDofs+nExtraDirichlet+nz_loc;
        }
        else
        {
            m_numLocalDirBndCoeffs = nLocDirBndCondDofs+nExtraDirichlet;
        }

        /**
         * STEP 3: Set up an array which contains the offset information of
         * the different graph vertices.
         *
         * This basically means to identify how many global degrees of
         * freedom the individual graph vertices correspond. Obviously,
         * the graph vertices corresponding to the mesh-vertices account
         * for a single global DOF. However, the graph vertices
         * corresponding to the element edges correspond to 2*(N-2) global DOF
         * where N is equal to the number of boundary modes on this edge.
         */
        Array<OneD, int> graphVertOffset(nvel*nz_loc*(ReorderedGraphVertId[0].size() + ReorderedGraphVertId[1].size()),0);
        graphVertOffset[0] = 0;

        m_signChange = false;

        for(i = 0; i < nel; ++i)
        {
            locExpansion = locExpVector[i]->as<StdRegions::StdExpansion2D>();

            for(j = 0; j < locExpansion->GetNedges(); ++j)
            {
                nEdgeCoeffs = locExpansion->GetEdgeNcoeffs(j);
                meshEdgeId = (locExpansion->as<LocalRegions::Expansion2D>()
                                          ->GetGeom2D())->GetEid(j);
                meshVertId = (locExpansion->as<LocalRegions::Expansion2D>()
                                          ->GetGeom2D())->GetVid(j);

                for(k = 0; k < nvel*nz_loc; ++k)
                {
                    graphVertOffset[ReorderedGraphVertId[0][meshVertId]*nvel*nz_loc+k] = 1;
                    graphVertOffset[ReorderedGraphVertId[1][meshEdgeId]*nvel*nz_loc+k] = (nEdgeCoeffs-2);
                }

                bType = locExpansion->GetEdgeBasisType(j);
                // need a sign vector for modal expansions if nEdgeCoeffs >=4
                if( (nEdgeCoeffs >= 4)&&
                    ( (bType == LibUtilities::eModified_A)||
                      (bType == LibUtilities::eModified_B) ) )
                {
                    m_signChange = true;
                }
            }
        }

        // Add mean pressure modes;
        for(i = 0; i < nel; ++i)
        {
            graphVertOffset[(ReorderedGraphVertId[1][AddMeanPressureToEdgeId[i]]+1)*nvel*nz_loc-1] += nz_loc;
            //graphVertOffset[(ReorderedGraphVertId[1][AddMeanPressureToEdgeId[i]])*nvel*nz_loc] += nz_loc;
        }

        // Negate the vertices and edges with only a partial
        // Dirichlet conditon. Essentially we check to see if an edge
        // has a mixed Dirichlet with Neumann/Robin Condition and if
        // so negate the offset associated with this vertex.

        map<int,int> DirVertChk;

        for(i = 0; i < bndConditionsVec[0].size(); ++i)
        {
            cnt = 0;
            for(j = 0; j < nvel; ++j)
            {
                if(bndConditionsVec[j][i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    cnt ++;
                }
            }

            // Case where partial Dirichlet boundary condition
            if((cnt > 0)&&(cnt < nvel))
            {
                for(j  = 0; j < nvel; ++j)
                {
                    if(bndConditionsVec[j][i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        //negate graph offsets which should be
                        //Dirichlet conditions
                        for(k = 0; k < bndCondExp[i]->GetNumElmts(); ++k)
                        {
                            // vertices with mix condition;
                            id = bndCondExp[i]->GetExp(k)
                                              ->as<LocalRegions::Expansion1D>()
                                              ->GetGeom1D()->GetVid(0);
                            if(DirVertChk.count(id*nvel+j) == 0)
                            {
                                DirVertChk[id*nvel+j] = 1;
                                for(n = 0; n < nz_loc; ++n)
                                {
                                    graphVertOffset[ReorderedGraphVertId[0][id]*nvel*nz_loc+j*nz_loc + n] *= -1;
                                }
                            }

                            id = bndCondExp[i]->GetExp(k)
                                              ->as<LocalRegions::Expansion1D>()
                                              ->GetGeom1D()->GetVid(1);
                            if(DirVertChk.count(id*nvel+j) == 0)
                            {
                                DirVertChk[id*nvel+j] = 1;
                                for(n = 0; n < nz_loc; ++n)
                                {
                                    graphVertOffset[ReorderedGraphVertId[0][id]*nvel*nz_loc+j*nz_loc+n] *= -1;
                                }
                            }

                            // edges with mixed id;
                            id = bndCondExp[i]->GetExp(k)
                                              ->as<LocalRegions::Expansion1D>()
                                              ->GetGeom1D()->GetGlobalID();
                            for(n = 0; n < nz_loc; ++n)
                            {
                                graphVertOffset[ReorderedGraphVertId[1][id]*nvel*nz_loc+j*nz_loc +n] *= -1;
                            }
                        }
                    }
                }
            }
        }


        cnt = 0;
        // assemble accumulative list of full Dirichlet values.
        for(i = 0; i < firstNonDirGraphVertId*nvel*nz_loc; ++i)
        {
            diff = abs(graphVertOffset[i]);
            graphVertOffset[i] = cnt;
            cnt += diff;
        }

        // set Dirichlet values with negative values to Dirichlet value
        for(i = firstNonDirGraphVertId*nvel*nz_loc; i <  graphVertOffset.size(); ++i)
        {
            if(graphVertOffset[i] < 0)
            {
                diff = -graphVertOffset[i];
                graphVertOffset[i] = -cnt;
                cnt += diff;
            }
        }

        // Accumulate all interior degrees of freedom with positive values
        m_numGlobalDirBndCoeffs = cnt;

        // offset values
        for(i = firstNonDirGraphVertId*nvel*nz_loc; i < graphVertOffset.size(); ++i)
        {
            if(graphVertOffset[i] >= 0)
            {
                diff = graphVertOffset[i];
                graphVertOffset[i] = cnt;
                cnt += diff;
            }
        }

        // Finally set negative entries (corresponding to Dirichlet
        // values ) to be positive
        for(i = firstNonDirGraphVertId*nvel*nz_loc; i < graphVertOffset.size(); ++i)
        {
            if(graphVertOffset[i] < 0)
            {
                graphVertOffset[i] = -graphVertOffset[i];
            }
        }


        // Allocate the proper amount of space for the class-data and fill
        // information that is already known
        cnt = 0;
        m_numLocalBndCoeffs = 0;
        m_numLocalCoeffs = 0;

        for(i = 0; i < nel; ++i)
        {
            m_numLocalBndCoeffs += nz_loc*(nvel*locExpVector[i]->NumBndryCoeffs() + 1);
            // add these coeffs up separately since
            // pressure->GetNcoeffs can include the coefficient in
            // multiple planes.
            m_numLocalCoeffs += (pressure->GetExp(i)->GetNcoeffs()-1)*nz_loc;
        }

        m_numLocalCoeffs += m_numLocalBndCoeffs;


        m_localToGlobalMap    = Array<OneD, int>(m_numLocalCoeffs,-1);
        m_localToGlobalBndMap = Array<OneD, int>(m_numLocalBndCoeffs,-1);
        m_bndCondIDToGlobalTraceID = Array<OneD, int>(nLocBndCondDofs,-1);


        // Set default sign array.
        m_localToGlobalSign    = Array<OneD, NekDouble>(m_numLocalCoeffs,1.0);
        m_localToGlobalBndSign = Array<OneD, NekDouble>(m_numLocalBndCoeffs,1.0);

        m_staticCondLevel = staticCondLevel;
        m_numPatches = nel;
        
        m_numLocalBndCoeffsPerPatch = Array<OneD, unsigned int>(nel);
        m_numLocalIntCoeffsPerPatch = Array<OneD, unsigned int>(nel);

        for(i = 0; i < nel; ++i)
        {
            m_numLocalBndCoeffsPerPatch[i] = (unsigned int) nz_loc*(nvel*locExpVector[i]->NumBndryCoeffs() + 1);
            m_numLocalIntCoeffsPerPatch[i] = (unsigned int) nz_loc*(pressure->GetExp(i)->GetNcoeffs()-1);
        }

        // Set up local to local bnd and local int maps
        m_localToLocalBndMap    = Array<OneD, int>(m_numLocalBndCoeffs,-1);
        m_localToLocalIntMap    = Array<OneD, int>(m_numLocalCoeffs-
                                                   m_numLocalBndCoeffs,-1);

        int bndcnt = 0;
        int intcnt = 0;
        cnt = 0; 
        for(i = 0; i < nel; ++i)
        {
            for(j = 0; j < nz_loc*(nvel*locExpVector[i]->NumBndryCoeffs()); ++j)
            {
                m_localToLocalBndMap[bndcnt++]  = cnt++;
            }

            for(n = 0; n < nz_loc; ++n)
            {
                m_localToLocalBndMap[bndcnt++]  = cnt++;
                for(j = 1; j < pressure->GetExp(i)->GetNcoeffs(); ++j)
                {
                    m_localToLocalIntMap[intcnt++]  = cnt++;
                }
            }

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
        int nv,velnbndry;
        Array<OneD, unsigned int> bmap;


        // Loop over all the elements in the domain in shuffled
        // ordering (element type consistency)
        for(i = 0; i < nel; ++i)
        {
            locExpansion = locExpVector[i]->as<StdRegions::StdExpansion2D>();

            velnbndry = locExpansion->NumBndryCoeffs();

            // Require an inverse ordering of the bmap system to store
            // local numbering system. Therefore get hold of elemental
            // bmap and set up an inverse map
            map<int,int> inv_bmap;
            locExpansion->GetBoundaryMap(bmap);
            for(j = 0; j < bmap.size(); ++j)
            {
                inv_bmap[bmap[j]] = j;
            }

            // Loop over all edges (and vertices) of element i
            for(j = 0; j < locExpansion->GetNedges(); ++j)
            {
                nEdgeInteriorCoeffs = locExpansion->GetEdgeNcoeffs(j)-2;
                edgeOrient = (locExpansion->as<LocalRegions::Expansion2D>()
                                            ->GetGeom2D())->GetEorient(j);
                meshEdgeId = (locExpansion->as<LocalRegions::Expansion2D>()
                                            ->GetGeom2D())->GetEid(j);
                meshVertId = (locExpansion->as<LocalRegions::Expansion2D>()
                                            ->GetGeom2D())->GetVid(j);

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

                locExpansion->GetEdgeInteriorMap(j,edgeOrient,
                                          edgeInteriorMap,edgeInteriorSign);

                // Set the global DOF for vertex j of element i
                for(nv = 0; nv < nvel*nz_loc; ++nv)
                {
                    m_localToGlobalMap[cnt+nv*velnbndry+
                               inv_bmap[locExpansion->GetVertexMap(j)]]
                        = graphVertOffset[ReorderedGraphVertId[0]
                                          [meshVertId]*nvel*nz_loc+ nv];
                    
                    // Set the global DOF's for the interior modes of edge j
                    for(k = 0; k < nEdgeInteriorCoeffs; ++k)
                    {
                        m_localToGlobalMap[cnt+nv*velnbndry+
                                        inv_bmap[edgeInteriorMap[k]]] =
                          graphVertOffset[ReorderedGraphVertId[1]
                                          [meshEdgeId]*nvel*nz_loc+nv]+k;
                    }
                }

                // Fill the sign vector if required
                if(m_signChange)
                {
                    for(nv = 0; nv < nvel*nz_loc; ++nv)
                    {
                        for(k = 0; k < nEdgeInteriorCoeffs; ++k)
                        {
                            m_localToGlobalSign[cnt+nv*velnbndry +
                                                inv_bmap[edgeInteriorMap[k]]]
                                = (NekDouble) edgeInteriorSign[k];
                        }
                    }
                }
            }

            // use difference between two edges of the
            // AddMeanPressureEdgeId to det nEdgeInteriorCoeffs.
            nEdgeInteriorCoeffs = graphVertOffset[(ReorderedGraphVertId[1][AddMeanPressureToEdgeId[i]])*nvel*nz_loc+1] - graphVertOffset[(ReorderedGraphVertId[1][AddMeanPressureToEdgeId[i]])*nvel*nz_loc];

            int psize = pressure->GetExp(i)->GetNcoeffs();
            for(n = 0; n < nz_loc; ++n)
            {
                m_localToGlobalMap[cnt + nz_loc*nvel*velnbndry + n*psize] =
                    graphVertOffset[(ReorderedGraphVertId[1]
                             [AddMeanPressureToEdgeId[i]]+1)*nvel*nz_loc-1]+
                    nEdgeInteriorCoeffs +
                    pressureEdgeOffset[AddMeanPressureToEdgeId[i]];

                pressureEdgeOffset[AddMeanPressureToEdgeId[i]] += 1;
            }

            cnt += (velnbndry*nvel+ psize)*nz_loc;
        }

        // Set up the mapping for the boundary conditions
        offset = cnt = 0;
        for(nv = 0; nv < nvel; ++nv)
        {
            for(i = 0; i < bndCondExp.size(); i++)
            {
                if (bndConditionsVec[nv][i]->GetBoundaryConditionType()==
                    SpatialDomains::ePeriodic)
                {
                    continue;
                }

                for(n = 0; n < nz_loc; ++n)
                {
                    int ncoeffcnt = 0;
                    for(j = 0; j < bndCondExp[i]->GetNumElmts(); j++)
                    {
                        bndSegExp  = bndCondExp[i]->GetExp(j)
                                                  ->as<LocalRegions::SegExp>();

                        cnt = offset + bndCondExp[i]->GetCoeff_Offset(j);
                        for(k = 0; k < 2; k++)
                        {
                            meshVertId = (bndSegExp->GetGeom1D())->GetVid(k);
                            m_bndCondIDToGlobalTraceID[cnt+
                                         bndSegExp->GetVertexMap(k)] =
                               graphVertOffset[ReorderedGraphVertId[0]
                                      [meshVertId]*nvel*nz_loc+nv*nz_loc+n];
                        }

                        meshEdgeId = (bndSegExp->GetGeom1D())->GetGlobalID();
                        bndEdgeCnt = 0;
                        nEdgeCoeffs = bndSegExp->GetNcoeffs();
                        for(k = 0; k < nEdgeCoeffs; k++)
                        {
                            if(m_bndCondIDToGlobalTraceID[cnt+k] == -1)
                            {
                                m_bndCondIDToGlobalTraceID[cnt+k] =
                                    graphVertOffset[ReorderedGraphVertId[1]
                                     [meshEdgeId]*nvel*nz_loc+nv*nz_loc+n]+bndEdgeCnt;
                                bndEdgeCnt++;
                            }
                        }
                        ncoeffcnt += nEdgeCoeffs;
                    }
                    // Note: Can not use bndCondExp[i]->GetNcoeffs()
                    // due to homogeneous extension not returning just
                    // the value per plane
                    offset += ncoeffcnt;
                }
            }
        }

        globalId = Vmath::Vmax(m_numLocalCoeffs,&m_localToGlobalMap[0],1)+1;
        m_numGlobalBndCoeffs = globalId;

        /**
         * STEP 5: The boundary condition mapping is generated from the
         * same vertex renumbering and fill in a unique interior map.
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
        if( m_session->MatchSolverInfoAsEnum("GlobalSysSoln", MultiRegions::eDirectMultiLevelStaticCond) )
        {
            if (m_staticCondLevel < (bottomUpGraph->GetNlevels()-1))
            {
                Array<OneD, int> vwgts_perm(
                    Dofs[0].size()+Dofs[1].size()-firstNonDirGraphVertId);
                for(i = 0; i < locExpVector.size(); ++i)
                {
                    locExpansion = locExpVector[i]
                                            ->as<StdRegions::StdExpansion2D>();
                    for(j = 0; j < locExpansion->GetNverts(); ++j)
                    {
                        meshEdgeId = (locExpansion
                                        ->as<LocalRegions::Expansion2D>()
                                        ->GetGeom2D())->GetEid(j);
                        meshVertId = (locExpansion
                                        ->as<LocalRegions::Expansion2D>()
                                        ->GetGeom2D())->GetVid(j);

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
    }



void CoupledLocalToGlobalC0ContMap::FindEdgeIdToAddMeanPressure(vector<map<int,int> > &ReorderedGraphVertId,
										 int &nel, const LocalRegions::ExpansionVector &locExpVector,
										 int &edgeId, int &vertId, int &firstNonDirGraphVertId, map<int,int> &IsDirEdgeDof,
										 MultiRegions::BottomUpSubStructuredGraphSharedPtr &bottomUpGraph,
										 Array<OneD, int> &AddMeanPressureToEdgeId)
{

	int i,j,k;

	// Make list of homogeneous graph edges to elmt mappings
	Array<TwoD, int> EdgeIdToElmts(ReorderedGraphVertId[1].size(),2,-1);
	map<int,int> HomGraphEdgeIdToEdgeId;

	for(i = 0; i < nel; ++i)
	{
		for(j = 0; j < locExpVector[i]->GetNverts(); ++j)
		{
			edgeId = (locExpVector[i]->as<LocalRegions::Expansion2D>()
			                         ->GetGeom2D())->GetEid(j);

			// note second condition stops us using mixed boundary condition
			if((ReorderedGraphVertId[1][edgeId] >= firstNonDirGraphVertId)
			   && (IsDirEdgeDof.count(edgeId) == 0))
			{
				HomGraphEdgeIdToEdgeId[ReorderedGraphVertId[1][edgeId]-firstNonDirGraphVertId] = edgeId;

				if(EdgeIdToElmts[edgeId][0] == -1)
				{
					EdgeIdToElmts[edgeId][0] = i;
				}
				else
				{
					EdgeIdToElmts[edgeId][1] = i;
				}
			}
		}
	}

	// Start at second to last level and find edge on boundary
	// to attach element
	int nlevels = bottomUpGraph->GetNlevels();

	// determine a default edge to attach pressure modes to
	// which is part of the inner solve;
	int defedge = -1;

	vector<MultiRegions::SubGraphSharedPtr> bndgraphs = bottomUpGraph->GetInteriorBlocks(nlevels);
	for(i = 0; i < bndgraphs.size(); ++i)
	{
		int GlobIdOffset = bndgraphs[i]->GetIdOffset();

		for(j = 0; j < bndgraphs[i]->GetNverts(); ++j)
		{
			// find edge in graph vert list
			if(HomGraphEdgeIdToEdgeId.count(GlobIdOffset+j) != 0)
			{
				edgeId = HomGraphEdgeIdToEdgeId[GlobIdOffset+j];

				if(defedge == -1)
				{
					defedge = edgeId;
					break;
				}
			}
		}
		if(defedge != -1)
		{
			break;
		}
	}

	for(int n = 1; n < nlevels; ++n)
	{
		// produce a map with a key that is the element id
		// that contains which next level patch it belongs to
		vector<MultiRegions::SubGraphSharedPtr> bndgraphs = bottomUpGraph->GetInteriorBlocks(n+1);

		// Fill next level graph  of adjacent elements and their level
		map<int,int> ElmtInBndry;

		for(i = 0; i < bndgraphs.size(); ++i)
		{
			int GlobIdOffset = bndgraphs[i]->GetIdOffset();

			for(j = 0; j < bndgraphs[i]->GetNverts(); ++j)
			{
				// find edge in graph vert list
				if(HomGraphEdgeIdToEdgeId.count(GlobIdOffset+j) != 0)
				{
					edgeId = HomGraphEdgeIdToEdgeId[GlobIdOffset+j];

					if(EdgeIdToElmts[edgeId][0] != -1)
					{
						ElmtInBndry[EdgeIdToElmts[edgeId][0]] = i;
					}
					if(EdgeIdToElmts[edgeId][1] != -1)
					{
						ElmtInBndry[EdgeIdToElmts[edgeId][1]] = i;
					}
				}
			}
		}

		// Now search interior patches in this level for edges
		// that share the same element as a boundary edge and
		// assign this elmt that boundary edge
		vector<MultiRegions::SubGraphSharedPtr> intgraphs = bottomUpGraph->GetInteriorBlocks(n);
		for(i = 0; i < intgraphs.size(); ++i)
		{
			int GlobIdOffset = intgraphs[i]->GetIdOffset();
			bool SetEdge = false;
			int elmtid = 0;
			for(j = 0; j < intgraphs[i]->GetNverts(); ++j)
			{
				// Check to see if graph vert is an edge
				if(HomGraphEdgeIdToEdgeId.count(GlobIdOffset+j) != 0)
				{
					edgeId = HomGraphEdgeIdToEdgeId[GlobIdOffset+j];

					for(k = 0; k < 2; ++k)
					{
						// relevant edge id
						elmtid = EdgeIdToElmts[edgeId][k];

						if(elmtid != -1)
						{
							auto mapIt = ElmtInBndry.find(elmtid);

							if(mapIt != ElmtInBndry.end())
							{
								// now find a edge in the next level boundary graph
								int GlobIdOffset1 = bndgraphs[mapIt->second]->GetIdOffset();
								for(int l = 0; l < bndgraphs[mapIt->second]->GetNverts(); ++l)
								{
									// find edge in graph vert list
									if(HomGraphEdgeIdToEdgeId.count(GlobIdOffset1+l) != 0)
									{
										//June 2012: commenting this condition apparently
										//solved the bug caused by the edge reordering procedure

										//if(AddMeanPressureToEdgeId[elmtid] == -1)
										//{

										//AddMeanPressureToEdgeId[elmtid] = HomGraphEdgeIdToEdgeId[GlobIdOffset1+l];
										AddMeanPressureToEdgeId[elmtid] = defedge;

										//}
										SetEdge = true;
										break;
									}
								}
							}
						}
					}
				}
			}


			// if we have failed to find matching edge in next
			// level patch boundary then set last found elmt
			// associated to this interior patch to the
			// default edget value
			if(SetEdge == false)
			{
				if(elmtid == -1) // find an elmtid in patch
				{
					for(j = 0; j < intgraphs[i]->GetNverts(); ++j)
					{
						if(HomGraphEdgeIdToEdgeId.count(GlobIdOffset+j) != 0)
						{
							edgeId = HomGraphEdgeIdToEdgeId[GlobIdOffset+j];
							for(k = 0; k < 2; ++k)
							{
								// relevant edge id
								elmtid = EdgeIdToElmts[edgeId][k];
								if(elmtid != -1)
								{
									break;
								}
							}
						}
						if(elmtid != -1)
						{
							break;
						}
					}
				}
				if(AddMeanPressureToEdgeId[elmtid] == -1)
				{
					AddMeanPressureToEdgeId[elmtid] = defedge;
				}
			}
		}
	}

}




}

