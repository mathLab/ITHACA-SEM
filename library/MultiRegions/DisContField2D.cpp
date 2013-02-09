//////////////////////////////////////////////////////////////////////////////
//
// File DisContField2D.cpp
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
// Description: Field definition for 2D domain with boundary conditions using
// LDG flux.
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/DisContField2D.h>
#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion.h>     // for Expansion
#include <SpatialDomains/MeshGraph2D.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>


namespace Nektar
{
    namespace MultiRegions
    {
        DisContField2D::DisContField2D(void):
            ExpList2D          (),
            m_bndCondExpansions(),
            m_bndConditions    (),
            m_trace            (NullExpListSharedPtr)
        {
        }

        DisContField2D::DisContField2D(
            const DisContField2D &In, 
            const bool            DeclareCoeffPhysArrays) :
            ExpList2D            (In,DeclareCoeffPhysArrays),
            m_bndCondExpansions  (In.m_bndCondExpansions),
            m_bndConditions      (In.m_bndConditions),
            m_globalBndMat       (In.m_globalBndMat),
            m_trace              (In.m_trace),
            m_traceMap           (In.m_traceMap),
            m_boundaryEdges      (In.m_boundaryEdges),
            m_periodicEdges      (In.m_periodicEdges),
            m_periodicVertices   (In.m_periodicVertices),
            m_perEdgeToExpMap    (In.m_perEdgeToExpMap)
        {
        }

        DisContField2D::DisContField2D(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const SpatialDomains::MeshGraphSharedPtr   &graph2D,
            const std::string                          &variable,
            const bool                                  SetUpJustDG,
            const bool                                  DeclareCoeffPhysArrays)
            : ExpList2D(pSession,graph2D,DeclareCoeffPhysArrays,variable),
              m_bndCondExpansions(),
              m_bndConditions(),
              m_trace(NullExpListSharedPtr)
        {
            SpatialDomains::BoundaryConditions bcs(m_session, graph2D);

            GenerateBoundaryConditionExpansion(graph2D,bcs,variable,
                                               DeclareCoeffPhysArrays);

            if(DeclareCoeffPhysArrays)
            {
                EvaluateBoundaryConditions();
            }

            ApplyGeomInfo();
            
            // Find periodic edges for this variable.
            FindPeriodicEdges(bcs, variable);

            if (SetUpJustDG)
            {
                SetUpDG();
            }
            else
            {
                // Set element edges to point to Robin BC edges if required.
                int i,cnt;
                Array<OneD, int> ElmtID,EdgeID;
                GetBoundaryToElmtMap(ElmtID,EdgeID);

                for(cnt = i = 0; i < m_bndCondExpansions.num_elements(); ++i)
                {
                    MultiRegions::ExpListSharedPtr locExpList;
                    int e;
                    locExpList = m_bndCondExpansions[i];
                    
                    for(e = 0; e < locExpList->GetExpSize(); ++e)
                    {
                        LocalRegions::Expansion2DSharedPtr exp2d
                            = boost::dynamic_pointer_cast<
                                LocalRegions::Expansion2D>((*m_exp)[ElmtID[cnt+e]]);
                        LocalRegions::Expansion1DSharedPtr exp1d
                            = boost::dynamic_pointer_cast<
                                LocalRegions::Expansion1D>(locExpList->GetExp(e));
                        LocalRegions::ExpansionSharedPtr   exp
                            = boost::dynamic_pointer_cast<
                                LocalRegions::Expansion>  (locExpList->GetExp(e));
                        
                        exp2d->SetEdgeExp(EdgeID[cnt+e],exp);
                        exp1d->SetAdjacentElementExp(EdgeID[cnt+e],exp2d);
                    }
                    cnt += m_bndCondExpansions[i]->GetExpSize();
                }
                
                if(m_session->DefinesSolverInfo("PROJECTION"))
                {
                    std::string ProjectStr = m_session->GetSolverInfo("PROJECTION");
                    if((ProjectStr == "MixedCGDG")||(ProjectStr == "Mixed_CG_Discontinuous"))
                    {
                        SetUpDG();
                    }
                    else
                    {
                        SetUpPhysNormals();
                    }
                }
                else
                {
                    SetUpPhysNormals();
                }
            }
        }


        // Copy type constructor which declares new boundary conditions
        // and re-uses mapping info and trace space if possible
        DisContField2D::DisContField2D(
            const DisContField2D                     &In,
            const SpatialDomains::MeshGraphSharedPtr &graph2D,
            const std::string                        &variable,
            const bool                                SetUpJustDG,
            const bool                                DeclareCoeffPhysArrays) :
            ExpList2D(In,DeclareCoeffPhysArrays),
            m_trace(NullExpListSharedPtr)
        {
            // Set up boundary conditions for this variable.
            SpatialDomains::BoundaryConditions bcs(m_session, graph2D);
            GenerateBoundaryConditionExpansion(graph2D,bcs,variable);
            
            if (DeclareCoeffPhysArrays)
            {
                EvaluateBoundaryConditions();
            }
            
            if (!SameTypeOfBoundaryConditions(In))
            {
                // Find periodic edges for this variable.
                FindPeriodicEdges(bcs, variable);
                
                if(SetUpJustDG)
                {
                    SetUpDG();
                }
                else
                {
                    // set elmt edges to point to robin bc edges if required.
                    int i, cnt = 0;
                    Array<OneD, int> ElmtID,EdgeID;
                    GetBoundaryToElmtMap(ElmtID,EdgeID);

                    for(i = 0; i < m_bndCondExpansions.num_elements(); ++i)
                    {
                        MultiRegions::ExpListSharedPtr locExpList;

                        int e;
                        locExpList = m_bndCondExpansions[i];
                        
                        for(e = 0; e < locExpList->GetExpSize(); ++e)
                        {
                            LocalRegions::Expansion2DSharedPtr exp2d
                                = boost::dynamic_pointer_cast<
                                    LocalRegions::Expansion2D>(
                                        (*m_exp)[ElmtID[cnt+e]]);
                            LocalRegions::Expansion1DSharedPtr exp1d
                                = boost::dynamic_pointer_cast<
                                    LocalRegions::Expansion1D>(
                                        locExpList->GetExp(e));
                            LocalRegions::ExpansionSharedPtr   exp
                                = boost::dynamic_pointer_cast<
                                    LocalRegions::Expansion>  (
                                        locExpList->GetExp(e));
                            
                            exp2d->SetEdgeExp(EdgeID[cnt+e],exp);
                            exp1d->SetAdjacentElementExp(EdgeID[cnt+e],exp2d);
                        }
                        cnt += m_bndCondExpansions[i]->GetExpSize();
                    }
                    

                    if(m_session->DefinesSolverInfo("PROJECTION"))
                    {
                        std::string ProjectStr = m_session->GetSolverInfo("PROJECTION");
                        if((ProjectStr == "MixedCGDG")||(ProjectStr == "Mixed_CG_Discontinuous"))
                        {
                            SetUpDG();
                        }
                        else
                        {
                            SetUpPhysNormals();
                        }
                    }
                    else
                    {
                        SetUpPhysNormals();
                    }
                }
            }
            else
            {
                if(SetUpJustDG)
                {
                    m_globalBndMat     = In.m_globalBndMat;
                    m_trace            = In.m_trace;
                    m_traceMap         = In.m_traceMap;
                    m_periodicEdges    = In.m_periodicEdges;
                    m_periodicVertices = In.m_periodicVertices;
                    m_boundaryEdges    = In.m_boundaryEdges;
                    m_perEdgeToExpMap  = In.m_perEdgeToExpMap;
                }
                else 
                {
                    m_globalBndMat     = In.m_globalBndMat;
                    m_trace            = In.m_trace;
                    m_traceMap         = In.m_traceMap;
                    m_periodicEdges    = In.m_periodicEdges;
                    m_periodicVertices = In.m_periodicVertices;
                    m_boundaryEdges    = In.m_boundaryEdges;
                    m_perEdgeToExpMap  = In.m_perEdgeToExpMap;
                    
                    // set elmt edges to point to robin bc edges if required.
                    int i, cnt = 0;
                    Array<OneD, int> ElmtID,EdgeID;
                    GetBoundaryToElmtMap(ElmtID,EdgeID);
					
                    for(i = 0; i < m_bndCondExpansions.num_elements(); ++i)
                    {
                        MultiRegions::ExpListSharedPtr locExpList;

                        int e;
                        locExpList = m_bndCondExpansions[i];
			
                        for(e = 0; e < locExpList->GetExpSize(); ++e)
                        {
                            LocalRegions::Expansion2DSharedPtr exp2d
                                = boost::dynamic_pointer_cast<
                                    LocalRegions::Expansion2D>(
                                        (*m_exp)[ElmtID[cnt+e]]);
                            LocalRegions::Expansion1DSharedPtr exp1d
                                = boost::dynamic_pointer_cast<
                                    LocalRegions::Expansion1D>(
                                        locExpList->GetExp(e));
                            LocalRegions::ExpansionSharedPtr   exp
                                = boost::dynamic_pointer_cast<
                                    LocalRegions::Expansion>  (
                                        locExpList->GetExp(e));
                            
                            exp2d->SetEdgeExp(EdgeID[cnt+e],exp);
                            exp1d->SetAdjacentElementExp(EdgeID[cnt+e],exp2d);
                        }
                        cnt += m_bndCondExpansions[i]->GetExpSize();
                    }
                    
                    SetUpPhysNormals();
                }
            }
        }

        /**
         * @brief Default destructor.
         */
        DisContField2D::~DisContField2D()
        {
            
        }
        
        GlobalLinSysSharedPtr DisContField2D::GetGlobalBndLinSys(
            const GlobalLinSysKey &mkey)
        {
            ASSERTL0(mkey.GetMatrixType() == StdRegions::eHybridDGHelmBndLam,
                     "Routine currently only tested for HybridDGHelmholtz");
            ASSERTL1(mkey.GetGlobalSysSolnType() == 
                         m_traceMap->GetGlobalSysSolnType(),
                     "The local to global map is not set up for the requested "
                     "solution type");

            GlobalLinSysSharedPtr glo_matrix;
            GlobalLinSysMap::iterator matrixIter = m_globalBndMat->find(mkey);

            if(matrixIter == m_globalBndMat->end())
            {
                glo_matrix = GenGlobalBndLinSys(mkey,m_traceMap);
                (*m_globalBndMat)[mkey] = glo_matrix;
            }
            else
            {
                glo_matrix = matrixIter->second;
            }

            return glo_matrix;
        }

        /**
         * @brief Set up all DG member variables and maps.
         */
        void DisContField2D::SetUpDG()
        {
            // Check for multiple calls
            if (m_trace != NullExpListSharedPtr)
            {
                return;
            }
            
            ExpList1DSharedPtr trace;
            SpatialDomains::MeshGraph2DSharedPtr graph2D = 
                boost::dynamic_pointer_cast<SpatialDomains::MeshGraph2D>(
                    m_graph);
                
            // Set up matrix map
            m_globalBndMat = MemoryManager<GlobalLinSysMap>::
                AllocateSharedPtr();
            
            // Set up trace space
            trace = MemoryManager<ExpList1D>::AllocateSharedPtr(
                m_bndCondExpansions, m_bndConditions, *m_exp, 
                graph2D, m_periodicEdges);
            
            m_trace = boost::dynamic_pointer_cast<ExpList>(trace);
            m_traceMap = MemoryManager<AssemblyMapDG>::
                AllocateSharedPtr(m_session, graph2D, trace, *this,
                                  m_bndCondExpansions, m_bndConditions,
                                  m_periodicEdges);
                
            Array<OneD, Array<OneD, StdRegions::StdExpansionSharedPtr> >
                &elmtToTrace = m_traceMap->GetElmtToTrace();
                
            // Scatter trace segments to 2D elements. For each element, we find
            // the trace segment associated to each edge. The element then
            // retains a pointer to the trace space segments, to ensure
            // uniqueness of normals when retrieving from two adjoining elements
            // which do not lie in a plane.
            for (int i = 0; i < m_exp->size(); ++i)
            {
                for (int j = 0; j < (*m_exp)[i]->GetNedges(); ++j)
                {
                    LocalRegions::Expansion2DSharedPtr exp2d =
                        boost::dynamic_pointer_cast<
                            LocalRegions::Expansion2D>((*m_exp)[i]);
                    LocalRegions::Expansion1DSharedPtr exp1d =
                        boost::dynamic_pointer_cast<
                            LocalRegions::Expansion1D>(elmtToTrace[i][j]);
                    LocalRegions::ExpansionSharedPtr exp =
                        boost::dynamic_pointer_cast<
                            LocalRegions::Expansion>  (elmtToTrace[i][j]);
                    exp2d->SetEdgeExp           (j, exp  );
                    exp1d->SetAdjacentElementExp(j, exp2d);
                }
            }
                
            // Set up physical normals
            SetUpPhysNormals();
            
            // Set up information for parallel jobs.
            for (int i = 0; i < m_trace->GetExpSize(); ++i)
            {
                LocalRegions::Expansion1DSharedPtr traceEl = 
                    boost::dynamic_pointer_cast<
                        LocalRegions::Expansion1D>(m_trace->GetExp(i));
                    
                int offset = m_trace->GetPhys_Offset(i);
                    
                if (m_traceMap->GetTraceToUniversalMapUnique(offset) < 0)
                {
                    traceEl->GetLeftAdjacentElementExp()->NegateEdgeNormal(
                        traceEl->GetLeftAdjacentElementEdge());
                }
            }
                
            int cnt, n, e;
                
            // Identify boundary edges
            for(cnt = 0, n = 0; n < m_bndCondExpansions.num_elements(); ++n)
            {
                if (m_bndConditions[n]->GetBoundaryConditionType() != 
                    SpatialDomains::ePeriodic)
                {
                    for(e = 0; e < m_bndCondExpansions[n]->GetExpSize(); ++e)
                    {
                        m_boundaryEdges.insert(m_trace->GetOffset_Elmt_Id(
                            m_traceMap->GetBndCondTraceToGlobalTraceMap(cnt+e)));
                    }
                }
                cnt += m_bndCondExpansions[n]->GetExpSize();
            }
                
            // Set up information for periodic boundary conditions.
            for (n = 0; n < m_exp->size(); ++n)
            {
                for (e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    map<int,int>::iterator it = m_periodicEdges.find(
                        (*m_exp)[n]->GetGeom2D()->GetEid(e));
                        
                    if (it != m_periodicEdges.end())
                    {
                        m_perEdgeToExpMap[it->first] = make_pair(n, e);
                    }
                }
            }
        }
        
        /**
         * For each boundary region, checks that the types and number of
         * boundary expansions in that region match.
         * @param   In          ContField2D to compare with.
         * @return True if boundary conditions match.
         */
        bool DisContField2D::SameTypeOfBoundaryConditions(
            const DisContField2D &In)
        {
            int i;
            bool returnval = true;

            for(i = 0; i < m_bndConditions.num_elements(); ++i)
            {

                // check to see if boundary condition type is the same
                // and there are the same number of boundary
                // conditions in the boundary definition.
                if((m_bndConditions[i]->GetBoundaryConditionType()
                    != In.m_bndConditions[i]->GetBoundaryConditionType())||
                   (m_bndCondExpansions[i]->GetExpSize()
                    != In.m_bndCondExpansions[i]->GetExpSize()))
                {
                    returnval = false;
                    break;
                }
            }

            // Compare with all other processes. Return true only if all
            // processes report having the same boundary conditions.
            int vSame = (returnval?1:0);
            m_comm->GetRowComm()->AllReduce(vSame, LibUtilities::ReduceMin);

            return (vSame == 1);
        }

        /**
         * \brief This function discretises the boundary conditions by setting
         * up a list of one-dimensional boundary expansions.
         *
         * According to their boundary region, the separate segmental boundary
         * expansions are bundled together in an object of the class
         * MultiRegions#ExpList1D.  
         *
         * \param graph2D   A mesh, containing information about the domain and
         *                  the spectral/hp element expansion.
         * \param bcs       An entity containing information about the
         *                  boundaries and boundary conditions.
         * \param variable  An optional parameter to indicate for which variable
         *                  the boundary conditions should be discretised.
         */
        void DisContField2D::GenerateBoundaryConditionExpansion(
            const SpatialDomains::MeshGraphSharedPtr &graph2D,
            const SpatialDomains::BoundaryConditions &bcs,
            const std::string &variable,
            const bool DeclareCoeffPhysArrays)
        {  	
            int i, cnt = 0;
            SpatialDomains::BoundaryConditionShPtr             bc;
            MultiRegions::ExpList1DSharedPtr                   locExpList;
            const SpatialDomains::BoundaryRegionCollection    &bregions = 
                bcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions = 
                bcs.GetBoundaryConditions();

            int nbnd = bregions.size();
            
            // count the number of non-periodic boundary regions
            for(i = 0; i < nbnd; ++i)
            {
                bc = GetBoundaryCondition(bconditions, i, variable);
                
                if (bc->GetBoundaryConditionType() != SpatialDomains::ePeriodic)
                {
                    cnt++;
                }
            }

            m_bndCondExpansions = 
                Array<OneD, MultiRegions::ExpListSharedPtr>(cnt);
            m_bndConditions     = 
                Array<OneD, SpatialDomains::BoundaryConditionShPtr>(cnt);
	    
            cnt = 0;

            // list non-periodic boundaries
            for(i = 0; i < nbnd; ++i)
            {
                bc = GetBoundaryCondition(bconditions, i, variable);
                if(bc->GetBoundaryConditionType() != SpatialDomains::ePeriodic)
                {
                    locExpList = MemoryManager<MultiRegions::ExpList1D>
                        ::AllocateSharedPtr(*(bregions[i]), graph2D, 
                                            DeclareCoeffPhysArrays, variable);
                    

                    // Set up normals on non-Dirichlet boundary conditions
                    if(bc->GetBoundaryConditionType() != 
                           SpatialDomains::eDirichlet)
                    {
                        SetUpPhysNormals();
                    }

                    m_bndCondExpansions[cnt]  = locExpList;
                    m_bndConditions[cnt]      = bc;
                    SpatialDomains::BndUserDefinedType type = 
                        m_bndConditions[cnt++]->GetUserDefined();
                    if (type == SpatialDomains::eI    || 
                        type == SpatialDomains::eCalcBC)
                    {
                        locExpList->SetUpPhysTangents(*m_exp);
                        SetUpPhysNormals();
                    }
                }
            }
        }

        /**
         * @brief Determine the periodic edges and vertices for the given graph.
         * 
         * @param   bcs         Information about the boundary conditions.
         * @param   variable    Specifies the field.
         */
        void DisContField2D::FindPeriodicEdges(
            const SpatialDomains::BoundaryConditions &bcs,
            const std::string                        &variable)
        {
            int i,k;

            const SpatialDomains::BoundaryRegionCollection &bregions
                = bcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions
                = bcs.GetBoundaryConditions();
            
            int region1ID;
            int region2ID;

            SpatialDomains::Composite comp1;
            SpatialDomains::Composite comp2;

            SpatialDomains::SegGeomSharedPtr segmentGeom1;
            SpatialDomains::SegGeomSharedPtr segmentGeom2;

            SpatialDomains::ElementEdgeVectorSharedPtr element1;
            SpatialDomains::ElementEdgeVectorSharedPtr element2;

            StdRegions::Orientation orient1;
            StdRegions::Orientation orient2;

            SpatialDomains::BoundaryConditionShPtr locBCond;

            // This std::map is a check so that the periodic pairs
            // are not treated twice
            map<int, int> doneBndRegions;

            int nbnd = bregions.size();

            for(i = 0; i < nbnd; ++i)
            {
                locBCond = GetBoundaryCondition(bconditions, i, variable);
                if(locBCond->GetBoundaryConditionType()
                   == SpatialDomains::ePeriodic)
                {
                    region1ID = i;
                    region2ID = (boost::static_pointer_cast<
                                 SpatialDomains::PeriodicBoundaryCondition
                                 >(locBCond))->m_connectedBoundaryRegion;

                    if(doneBndRegions.count(region1ID)==0)
                    {
                        ASSERTL0(bregions[region1ID]->size()
                                 == bregions[region2ID]->size(),
                                 "Size of the 2 periodic boundary regions "
                                 "should be equal");


                        map<int,int> periodicVertices;

                        SpatialDomains::BoundaryRegion::iterator bnd1It, bnd2It;
                        for(bnd1It =  bregions[region1ID]->begin(),
                            bnd2It =  bregions[region2ID]->begin();
                            bnd1It != bregions[region1ID]->end();
                            ++bnd1It, ++bnd2It)
                        {
                            comp1 = bnd1It->second;
                            comp2 = bnd2It->second;

                            ASSERTL0(comp1->size() == comp2->size(),
                                     "Size of the 2 periodic composites should "
                                     "be equal");

                            for(k = 0; k < comp1->size(); k++)
                            {
                                if(!(segmentGeom1
                                     = boost::dynamic_pointer_cast<
                                     SpatialDomains::SegGeom>((*comp1)[k]))
                                   || !(segmentGeom2
                                        = boost::dynamic_pointer_cast<
                                        SpatialDomains::SegGeom>((*comp2)[k]))
                                   )
                                {
                                    ASSERTL0(false,"dynamic cast to a SegGeom "
                                             "failed");
                                }

                                element1 = boost::dynamic_pointer_cast<
                                    SpatialDomains::MeshGraph2D>(m_graph)
                                        ->GetElementsFromEdge(segmentGeom1);
                                element2 = boost::dynamic_pointer_cast<
                                    SpatialDomains::MeshGraph2D>(m_graph)
                                        ->GetElementsFromEdge(segmentGeom2);

                                ASSERTL0(element1->size() == 1,
                                         "The periodic boundaries belong to "
                                         "more than one element of the mesh");
                                ASSERTL0(element2->size() == 1,
                                         "The periodic boundaries belong to "
                                         "more than one element of the mesh");

                                orient1 = (boost::dynamic_pointer_cast<
                                           SpatialDomains::Geometry2D>(
                                               (*element1)[0]->m_Element)
                                           )->GetEorient((*element1)[0]
                                                         ->m_EdgeIndx);
                                orient2 = (boost::dynamic_pointer_cast<
                                           SpatialDomains::Geometry2D>(
                                               (*element2)[0]->m_Element)
                                           )->GetEorient((*element2)[0]
                                                         ->m_EdgeIndx);

                                if(orient1!=orient2)
                                {
                                    // Extract the periodic edges.
                                    m_periodicEdges[segmentGeom1->GetEid()]
                                        = segmentGeom2->GetEid();
                                    m_periodicEdges[segmentGeom2->GetEid()]
                                        = segmentGeom1->GetEid();
                                    
                                    // Extract the periodic vertices.
                                    periodicVertices[segmentGeom1->GetVid(0)]
                                        = segmentGeom2->GetVid(0);
                                    periodicVertices[segmentGeom1->GetVid(1)]
                                        = segmentGeom2->GetVid(1);
                                }
                                else
                                {
                                    m_periodicEdges[segmentGeom1->GetEid()]
                                        = -segmentGeom2->GetEid();
                                    m_periodicEdges[segmentGeom2->GetEid()]
                                        = segmentGeom1->GetEid();
                                    
                                    periodicVertices[segmentGeom1->GetVid(0)]
                                        = segmentGeom2->GetVid(1);
                                    periodicVertices[segmentGeom1->GetVid(1)]
                                        = segmentGeom2->GetVid(0);
                                }
                            }
                        }
                        m_periodicVertices.push_back(periodicVertices);
                    }
                    else
                    {
                        ASSERTL0(doneBndRegions[region1ID]==region2ID,
                                 "Boundary regions are not mutually periodic");
                    }
                    doneBndRegions[region2ID] = region1ID;
                }
            }
        }

        bool DisContField2D::IsLeftAdjacentEdge(const int n, const int e)
        {
            set<int>::iterator     it;
            LocalRegions::Expansion1DSharedPtr traceEl = 
                boost::dynamic_pointer_cast<LocalRegions::Expansion1D>(
                    (m_traceMap->GetElmtToTrace())[n][e]);
            
            int offset = m_trace->GetPhys_Offset(traceEl->GetElmtId());
            
            bool fwd = true;
            if (traceEl->GetLeftAdjacentElementEdge () == -1 ||
                traceEl->GetRightAdjacentElementEdge() == -1)
            {
                // Boundary edge (1 connected element). Do nothing in
                // serial.
                //it = m_boundaryEdges.find(elmtToTrace[n][e]->GetElmtId());
                it = m_boundaryEdges.find(traceEl->GetElmtId());
                
                // If the edge does not have a boundary condition set on
                // it, then assume it is a partition edge.
                if (it == m_boundaryEdges.end())
                {
                    fwd = m_traceMap->
                        GetTraceToUniversalMapUnique(offset) > 0;
                }
            }
            else if (traceEl->GetLeftAdjacentElementEdge () != -1 &&
                     traceEl->GetRightAdjacentElementEdge() != -1)
            {
                // Non-boundary edge (2 connected elements).
                fwd = dynamic_cast<Nektar::StdRegions::StdExpansion*>
                    (traceEl->GetLeftAdjacentElementExp().get()) ==
                    (*m_exp)[n].get();
            }
            else
            {
                ASSERTL2(false, "Unconnected trace element!");
            }

            return fwd;
        }
            
        // Construct the two trace vectors of the inner and outer
        // trace solution from the field contained in m_phys, where
        // the Weak dirichlet boundary conditions are listed in the
        // outer part of the vecotr
        void DisContField2D::v_GetFwdBwdTracePhys(
            Array<OneD, NekDouble> &Fwd,
            Array<OneD, NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhys(m_phys,Fwd,Bwd);
        }

        /**
         * @brief This method extracts the "forward" and "backward" trace data
         * from the array @a field and puts the data into output vectors @a Fwd
         * and @a Bwd.
         * 
         * We first define the convention which defines "forwards" and
         * "backwards". First an association is made between the edge of each
         * element and its corresponding edge in the trace space using the
         * mapping #m_traceMap. The element can either be left-adjacent or
         * right-adjacent to this trace edge (see
         * Expansion1D::GetLeftAdjacentElementExp). Boundary edges are always
         * left-adjacent since left-adjacency is populated first.
         * 
         * If the element is left-adjacent we extract the edge trace data from
         * @a field into the forward trace space @a Fwd; otherwise, we place it
         * in the backwards trace space @a Bwd. In this way, we form a unique
         * set of trace normals since these are always extracted from
         * left-adjacent elements.
         *
         * @param field is a NekDouble array which contains the 2D data
         *              from which we wish to extract the backward and forward
         *              orientated trace/edge arrays.
         * @param Fwd   The resulting forwards space.
         * @param Bwd   The resulting backwards space.
         */
        void DisContField2D::v_GetFwdBwdTracePhys(
            const Array<OneD, const NekDouble> &field,
                  Array<OneD,       NekDouble> &Fwd,
                  Array<OneD,       NekDouble> &Bwd)
        {
            // Loop over elements and collect forward expansion
            int nexp = GetExpSize();
            int cnt, n, e, npts, phys_offset;
            Array<OneD,NekDouble> e_tmp;
            map<int,int>::iterator it2;
            boost::unordered_map<int,pair<int,int> >::iterator it3;

            Array<OneD, Array<OneD, StdRegions::StdExpansionSharedPtr> >
                &elmtToTrace = m_traceMap->GetElmtToTrace();
            
            // Zero forward/backward vectors.
            Vmath::Zero(Fwd.num_elements(), Fwd, 1);
            Vmath::Zero(Bwd.num_elements(), Bwd, 1);

            bool fwd = true;
            for(n = 0; n < nexp; ++n)
            {
                phys_offset = GetPhys_Offset(n);

                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    int offset = m_trace->GetPhys_Offset(
                        elmtToTrace[n][e]->GetElmtId());
                    
                    fwd = IsLeftAdjacentEdge(n,e);

                    if (fwd)
                    {
                        (*m_exp)[n]->GetEdgePhysVals(e, elmtToTrace[n][e],
                                                     field + phys_offset,
                                                     e_tmp = Fwd + offset);
                    }
                    else
                    {
                        (*m_exp)[n]->GetEdgePhysVals(e, elmtToTrace[n][e],
                                                     field + phys_offset,
                                                     e_tmp = Bwd + offset);
                    }
                    
                    // Check to see if this edge is periodic.
                    it2 = m_periodicEdges.find(
                        (*m_exp)[n]->GetGeom2D()->GetEid(e));
                    
                    if (it2 != m_periodicEdges.end())
                    {
                        it3 = m_perEdgeToExpMap.find(abs(it2->second));

                        ASSERTL2(fwd, "Periodic edge in non-forward space?");
                        ASSERTL2(it3 != m_perEdgeToExpMap.end(),
                                 "Periodic edge not found!");
                        
                        int offset2 = m_trace->GetPhys_Offset(
                            elmtToTrace[it3->second.first][it3->second.second]->
                                GetElmtId());
                        
                        /*
                         * Copy fwd -> bwd space, reverse if necessary. Note
                         * that for varying polynomial order this condition will
                         * not work (needs some kind of interpolation here).
                         */
                        if (it2->second                       < 0 || 
                            m_periodicEdges[abs(it2->second)] < 0)
                        {
                            Vmath::Reverse(elmtToTrace[n][e]->GetTotPoints(),
                                           &Fwd[offset], 1, &Bwd[offset2], 1);
                        }
                        else
                        {
                            Vmath::Vcopy  (elmtToTrace[n][e]->GetTotPoints(),
                                           &Fwd[offset], 1, &Bwd[offset2], 1);
                        }
                    }
                }
            }
            
            // Fill boundary conditions into missing elements.
            int id1, id2 = 0;
            
            for(cnt = n = 0; n < m_bndCondExpansions.num_elements(); ++n)
            {				
                if (m_bndConditions[n]->GetBoundaryConditionType() == 
                        SpatialDomains::eDirichlet)
                {
                    for(e = 0; e < m_bndCondExpansions[n]->GetExpSize(); ++e)
                    {
                        npts = m_bndCondExpansions[n]->GetExp(e)->GetNumPoints(0);
                        id1  = m_bndCondExpansions[n]->GetPhys_Offset(e);
                        id2  = m_trace->GetPhys_Offset(
                            m_traceMap->GetBndCondTraceToGlobalTraceMap(cnt+e));
                        Vmath::Vcopy(npts,
                            &(m_bndCondExpansions[n]->GetPhys())[id1], 1,
                            &Bwd[id2],                                 1);
                    }
                    
                    cnt += e;
                }
                else if (m_bndConditions[n]->GetBoundaryConditionType() == 
                             SpatialDomains::eNeumann || 
                         m_bndConditions[n]->GetBoundaryConditionType() == 
                             SpatialDomains::eRobin)
                {
                    for(e = 0; e < m_bndCondExpansions[n]->GetExpSize(); ++e)
                    {
                        npts = m_bndCondExpansions[n]->GetExp(e)->GetNumPoints(0);
                        id1  = m_bndCondExpansions[n]->GetPhys_Offset(e);
                        ASSERTL0((m_bndCondExpansions[n]->GetPhys())[id1]==0.0,
                                 "Method not set up for non-zero Neumann "
                                 "boundary condition");
                        id2  = m_trace->GetPhys_Offset(
                            m_traceMap->GetBndCondTraceToGlobalTraceMap(cnt+e));
                        Vmath::Vcopy(npts, &Fwd[id2], 1, &Bwd[id2], 1);
                    }
                    
                    cnt += e;
                }
                else if (m_bndConditions[n]->GetBoundaryConditionType() !=
                             SpatialDomains::ePeriodic)
                {
                    ASSERTL0(false,
                             "Method not set up for this boundary condition.");
                }
            }
            
            // Do parallel exchange for forwards/backwards spaces.
            m_traceMap->UniversalTraceAssemble(Fwd);
            m_traceMap->UniversalTraceAssemble(Bwd);
        }

        void DisContField2D::v_ExtractTracePhys(
            Array<OneD, NekDouble> &outarray)
        {
            ASSERTL1(m_physState == true,
                     "Field must be in physical state to extract trace space.");

            v_ExtractTracePhys(m_phys, outarray);
        }

        /**
         * @brief This method extracts the trace (edges in 2D) from the field @a
         * inarray and puts the values in @a outarray.
         *
         * It assumes the field is C0 continuous so that it can overwrite the
         * edge data when visited by the two adjacent elements.
         *
         * @param inarray   An array containing the 2D data from which we wish
         *                  to extract the edge data.
         * @param outarray  The resulting edge information.
         */
        void DisContField2D::v_ExtractTracePhys(
            const Array<OneD, const NekDouble> &inarray, 
                  Array<OneD,       NekDouble> &outarray)
        {
            // Loop over elemente and collect forward expansion
            int nexp = GetExpSize();
            int n,e,offset,phys_offset;
            Array<OneD,NekDouble> e_tmp;
            Array<OneD, Array<OneD, StdRegions::StdExpansionSharedPtr> >
                &elmtToTrace = m_traceMap->GetElmtToTrace();

            ASSERTL1(outarray.num_elements() >= m_trace->GetNpoints(),
                     "input array is of insufficient length");

            // use m_trace tmp space in element to fill values
            for(n  = 0; n < nexp; ++n)
            {
                phys_offset = GetPhys_Offset(n);
                
                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    offset = m_trace->GetPhys_Offset(
                        elmtToTrace[n][e]->GetElmtId());
                    (*m_exp)[n]->GetEdgePhysVals(e,  elmtToTrace[n][e],
                                                 inarray + phys_offset,
                                                 e_tmp = outarray + offset);
                }
            }
        }

        void DisContField2D::v_AddTraceIntegral(
            const Array<OneD, const NekDouble> &Fx,
            const Array<OneD, const NekDouble> &Fy,
                  Array<OneD,       NekDouble> &outarray)
        {
            int e,n,offset, t_offset;
            Array<OneD, NekDouble> e_outarray;
            Array<OneD, Array<OneD, StdRegions::StdExpansionSharedPtr> >
                &elmtToTrace = m_traceMap->GetElmtToTrace();

            for(n = 0; n < GetExpSize(); ++n)
            {
                offset = GetCoeff_Offset(n);
                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    t_offset = GetTrace()->GetPhys_Offset(
                        elmtToTrace[n][e]->GetElmtId());

                    (*m_exp)[n]->AddEdgeNormBoundaryInt(
                        e,elmtToTrace[n][e],
                        Fx + t_offset,
                        Fy + t_offset,
                        e_outarray = outarray+offset);
                }
            }
        }

        /**
         * @brief Add trace contributions into elemental coefficient spaces.
         * 
         * Given some quantity \f$ \vec{Fn} \f$, which conatins this
         * routine calculates the integral
         * 
         * \f[ 
         * \int_{\Omega^e} \vec{Fn}, \mathrm{d}S
         * \f] 
         * 
         * and adds this to the coefficient space provided by outarray.
         * 
         * @see Expansion2D::AddEdgeNormBoundaryInt
         * 
         * @param Fn        The trace quantities.
         * @param outarray  Resulting 2D coefficient space.
         */
        void DisContField2D::v_AddTraceIntegral(
            const Array<OneD, const NekDouble> &Fn, 
                  Array<OneD,       NekDouble> &outarray)
        {
            int e,n,offset, t_offset;
            Array<OneD, NekDouble> e_outarray;
            Array<OneD, Array<OneD, StdRegions::StdExpansionSharedPtr> >
                &elmtToTrace = m_traceMap->GetElmtToTrace();

            for(n = 0; n < GetExpSize(); ++n)
            {
                offset = GetCoeff_Offset(n);
                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    t_offset = GetTrace()->GetPhys_Offset(
                        elmtToTrace[n][e]->GetElmtId());
                    (*m_exp)[n]->AddEdgeNormBoundaryInt(
                        e, elmtToTrace[n][e], Fn+t_offset,
                        e_outarray = outarray+offset);
                }
            }
        }


        /**
         * @brief Add trace contributions into elemental coefficient spaces.
         * 
         * Given some quantity \f$ \vec{q} \f$, calculate the elemental integral
         * 
         * \f[ 
         * \int_{\Omega^e} \vec{q}, \mathrm{d}S
         * \f] 
         * 
         * and adds this to the coefficient space provided by
         * outarray. The value of q is determined from the routine
         * IsLeftAdjacentEdge() which if true we use Fwd else we use
         * Bwd
         * 
         * @see Expansion2D::AddEdgeNormBoundaryInt
         * 
         * @param Fwd       The trace quantities associated with left (fwd)
         *                  adjancent elmt.
         * @param Bwd       The trace quantities associated with right (bwd)
         *                  adjacent elet.
         * @param outarray  Resulting 2D coefficient space.
         */
        void DisContField2D::v_AddFwdBwdTraceIntegral(
            const Array<OneD, const NekDouble> &Fwd, 
            const Array<OneD, const NekDouble> &Bwd, 
                  Array<OneD,       NekDouble> &outarray)
        {
            int e,n,offset, t_offset;
            Array<OneD, NekDouble> e_outarray;
            Array<OneD, Array<OneD, StdRegions::StdExpansionSharedPtr> >
                &elmtToTrace = m_traceMap->GetElmtToTrace();

            for(n = 0; n < GetExpSize(); ++n)
            {
                offset = GetCoeff_Offset(n);
                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    t_offset = GetTrace()->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());
                    
                    // Evaluate upwind flux less local edge 
                    if(IsLeftAdjacentEdge(n,e))
                    {
                        (*m_exp)[n]->AddEdgeNormBoundaryInt(
                        e, elmtToTrace[n][e], Fwd+t_offset,
                        e_outarray = outarray+offset);
                    }
                    else
                    {
                        (*m_exp)[n]->AddEdgeNormBoundaryInt(
                        e, elmtToTrace[n][e], Bwd+t_offset,
                        e_outarray = outarray+offset);
                    }

                }
            }
        }

        /**
         * @brief Set up a list of element IDs and edge IDs that link to the
         * boundary conditions.
         */
        void DisContField2D::v_GetBoundaryToElmtMap(
            Array<OneD, int> &ElmtID, 
            Array<OneD, int> &EdgeID)
        {
            map<int, int> globalIdMap;
            int i,n;
            int cnt;
            int nbcs = 0;

            SpatialDomains::MeshGraph2DSharedPtr graph2D =
                boost::dynamic_pointer_cast<SpatialDomains::MeshGraph2D>(
                    m_graph);

            // Populate global ID map (takes global geometry ID to local
            // expansion list ID).
            for (i = 0; i < GetExpSize(); ++i)
            {
                globalIdMap[(*m_exp)[i]->GetGeom2D()->GetGlobalID()] = i;
            }
            
            // Determine number of boundary condition expansions.
            for(i = 0; i < m_bndConditions.num_elements(); ++i)
            {
                nbcs += m_bndCondExpansions[i]->GetExpSize();
            }

            // make sure arrays are of sufficient length
            if(ElmtID.num_elements() != nbcs)
            {
                ElmtID = Array<OneD, int>(nbcs);
            }

            if(EdgeID.num_elements() != nbcs)
            {
                EdgeID = Array<OneD, int>(nbcs);
            }

            for(cnt = n = 0; n < m_bndCondExpansions.num_elements(); ++n)
            {
                for(i = 0; i < m_bndCondExpansions[n]->GetExpSize(); ++i, ++cnt)
                {
                    // Use edge to element map from MeshGraph2D.
                    SpatialDomains::ElementEdgeVectorSharedPtr tmp =
                        graph2D->GetElementsFromEdge(
                            m_bndCondExpansions[n]->GetExp(i)->GetGeom1D());

                    ElmtID[cnt] = globalIdMap[(*tmp)[0]->
                                              m_Element->GetGlobalID()];
                    EdgeID[cnt] = (*tmp)[0]->m_EdgeIndx;
                }
            }
        }


        /** 
         * @brief Calculate the \f$ L^2 \f$ error of the \f$ Q_{\rm dir} \f$
         * derivative using the consistent DG evaluation of \f$ Q_{\rm dir} \f$.
         * 
         * The solution provided is of the primative variation at the quadrature
         * points and the derivative is compared to the discrete derivative at
         * these points, which is likely to be undesirable unless using a much
         * higher number of quadrature points than the polynomial order used to
         * evaluate \f$ Q_{\rm dir} \f$.
        */
        NekDouble DisContField2D::L2_DGDeriv(
            const int                           dir,
            const Array<OneD, const NekDouble> &soln)
        {
            int    i,e,ncoeff_edge;
            Array<OneD, const NekDouble> tmp_coeffs;
            Array<OneD, NekDouble> out_d(m_ncoeffs), out_tmp;

            Array<OneD, Array< OneD, StdRegions::StdExpansionSharedPtr> > 
                &elmtToTrace = m_traceMap->GetElmtToTrace();

            StdRegions::Orientation edgedir;

            int     eid,cnt;
            int     LocBndCoeffs = m_traceMap->GetNumLocalBndCoeffs();
            Array<OneD, NekDouble> loc_lambda(LocBndCoeffs), edge_lambda;
            m_traceMap->GlobalToLocalBnd(m_trace->GetCoeffs(),loc_lambda);

            edge_lambda = loc_lambda;
            
            // Calculate Q using standard DG formulation.
            for(i = cnt = 0; i < GetExpSize(); ++i)
            {
                eid = m_offset_elmt_id[i];
                // Probably a better way of setting up lambda than this.
                // Note cannot use PutCoeffsInToElmts since lambda space
                // is mapped during the solve.
                for(e = 0; e < (*m_exp)[eid]->GetNedges(); ++e)
                {
                    edgedir = (*m_exp)[eid]->GetEorient(e);

                    ncoeff_edge = elmtToTrace[eid][e]->GetNcoeffs();
                    elmtToTrace[eid][e]->SetCoeffsToOrientation(
                        edgedir,edge_lambda,edge_lambda);
                    Vmath::Vcopy(ncoeff_edge,edge_lambda,1,
                                 elmtToTrace[eid][e]->UpdateCoeffs(),1);
                    edge_lambda = edge_lambda + ncoeff_edge;
                }

                (*m_exp)[eid]->DGDeriv(dir,
                                       tmp_coeffs=m_coeffs+m_coeff_offset[eid],
                                       elmtToTrace[eid],
                                       out_tmp = out_d+cnt);
                cnt  += (*m_exp)[eid]->GetNcoeffs();
            }
            
            BwdTrans(out_d,m_phys);
            Vmath::Vsub(m_npoints,m_phys,1,soln,1,m_phys,1);

            return L2();
        }

        void DisContField2D::v_HelmSolve(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray,
                const FlagList &flags,
                const StdRegions::ConstFactorMap &factors,
                const StdRegions::VarCoeffMap &varcoeff,
                const Array<OneD, const NekDouble> &dirForcing)
        {
            int i,j,n,cnt,cnt1,nbndry;
            int nexp = GetExpSize();
            StdRegions::StdExpansionSharedPtr BndExp;

            Array<OneD,NekDouble> f(m_ncoeffs);
            DNekVec F(m_ncoeffs,f,eWrapper);
            Array<OneD,NekDouble> e_f, e_l;

            //----------------------------------
            //  Setup RHS Inner product
            //----------------------------------
            IProductWRTBase(inarray,f);
            Vmath::Neg(m_ncoeffs,f,1);

            //----------------------------------
            //  Solve continuous flux System
            //----------------------------------
            int GloBndDofs   = m_traceMap->GetNumGlobalBndCoeffs();
            int NumDirichlet = m_traceMap->GetNumLocalDirBndCoeffs();
            int e_ncoeffs,id;

            // Retrieve block matrix of U^e
            GlobalMatrixKey HDGLamToUKey(StdRegions::eHybridDGLamToU,
                NullAssemblyMapSharedPtr,factors,varcoeff);
            const DNekScalBlkMatSharedPtr &HDGLamToU = GetBlockMatrix(
                HDGLamToUKey);

            // Retrieve global trace space storage, \Lambda, from trace
            // expansion
            Array<OneD,NekDouble> BndSol = m_trace->UpdateCoeffs();

            // Create trace space forcing, F
            Array<OneD,NekDouble> BndRhs(GloBndDofs,0.0);

            // Zero \Lambda
            Vmath::Zero(GloBndDofs,BndSol,1);

            // Retrieve number of local trace space coefficients N_{\lambda},
            // and set up local elemental trace solution \lambda^e.
            int     LocBndCoeffs = m_traceMap->GetNumLocalBndCoeffs();
            Array<OneD, NekDouble> loc_lambda(LocBndCoeffs);
            DNekVec LocLambda(LocBndCoeffs,loc_lambda,eWrapper);

            //----------------------------------
            // Evaluate Trace Forcing vector F
            // Kirby et al, 2010, P23, Step 5.
            //----------------------------------
            // Loop over all expansions in the domain
            for(cnt = cnt1 = n = 0; n < nexp; ++n)
            {
                nbndry = (*m_exp)[m_offset_elmt_id[n]]->NumDGBndryCoeffs();

                e_ncoeffs = (*m_exp)[m_offset_elmt_id[n]]->GetNcoeffs();
                e_f       = f + cnt;
                e_l       = loc_lambda + cnt1;

                // Local trace space \lambda^e
                DNekVec     Floc    (nbndry, e_l, eWrapper);
                // Local forcing f^e
                DNekVec     ElmtFce (e_ncoeffs, e_f, eWrapper);
                // Compute local (U^e)^{\top} f^e
                Floc = Transpose(*(HDGLamToU->GetBlock(n,n)))*ElmtFce;

                cnt   += e_ncoeffs;
                cnt1  += nbndry;
            }

            // Assemble local \lambda_e into global \Lambda
            m_traceMap->AssembleBnd(loc_lambda,BndRhs);

            // Copy Dirichlet boundary conditions and weak forcing into trace
            // space
            cnt = 0;
            for(i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                if(m_bndConditions[i]->GetBoundaryConditionType() == 
                       SpatialDomains::eDirichlet)
                {
                    for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); ++j)
                    {
                        id = m_traceMap->GetBndCondCoeffsToGlobalCoeffsMap(cnt++);
                        BndSol[id] = m_bndCondExpansions[i]->GetCoeffs()[j];
                    }
                }
                else
                {
                    //Add weak boundary condition to trace forcing
                    for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); ++j)
                    {
                        id = m_traceMap->GetBndCondCoeffsToGlobalCoeffsMap(cnt++);
                        BndRhs[id] += m_bndCondExpansions[i]->GetCoeffs()[j];
                    }
                }
            }

            //----------------------------------
            // Solve trace problem: \Lambda = K^{-1} F
            // K is the HybridDGHelmBndLam matrix.
            //----------------------------------
            if(GloBndDofs - NumDirichlet > 0)
            {
                GlobalLinSysKey       key(StdRegions::eHybridDGHelmBndLam,
                                          m_traceMap,factors,varcoeff);
                GlobalLinSysSharedPtr LinSys = GetGlobalBndLinSys(key);
                LinSys->Solve(BndRhs,BndSol,m_traceMap);
            }

            //----------------------------------
            // Internal element solves
            //----------------------------------
            GlobalMatrixKey invHDGhelmkey(StdRegions::eInvHybridDGHelmholtz,
                NullAssemblyMapSharedPtr,factors,varcoeff);
            const DNekScalBlkMatSharedPtr& InvHDGHelm = GetBlockMatrix(
                invHDGhelmkey);
            DNekVec out(m_ncoeffs,outarray,eWrapper);
            Vmath::Zero(m_ncoeffs,outarray,1);

            // get local trace solution from BndSol
            m_traceMap->GlobalToLocalBnd(BndSol,loc_lambda);

            //  out =  u_f + u_lam = (*InvHDGHelm)*f + (LamtoU)*Lam
            out = (*InvHDGHelm)*F + (*HDGLamToU)*LocLambda;
        }


        /**
         * @brief Calculates the result of the multiplication of a global matrix
         * of type specified by @a mkey with a vector given by @a inarray.
         * 
         * @param mkey      Key representing desired matrix multiplication.
         * @param inarray   Input vector.
         * @param outarray  Resulting multiplication.
         */
        void DisContField2D::v_GeneralMatrixOp(
               const GlobalMatrixKey             &gkey,
               const Array<OneD,const NekDouble> &inarray,
               Array<OneD,      NekDouble> &outarray,
               CoeffState coeffstate)
        {
            int     LocBndCoeffs = m_traceMap->GetNumLocalBndCoeffs();
            Array<OneD, NekDouble> loc_lambda(LocBndCoeffs);
            DNekVec LocLambda(LocBndCoeffs,loc_lambda,eWrapper);
            const DNekScalBlkMatSharedPtr& HDGHelm = GetBlockMatrix(gkey);

            m_traceMap->GlobalToLocalBnd(inarray, loc_lambda);
            LocLambda = (*HDGHelm) * LocLambda;
            m_traceMap->AssembleBnd(loc_lambda,outarray);
        }

        /**
         * @brief Search through the edge expansions and identify which ones
         * have Robin/Mixed type boundary conditions.
         * 
         * If a Robin boundary is found then store the edge ID of the boundary
         * condition and the array of points of the physical space boundary
         * condition which are hold the boundary condition primitive variable
         * coefficient at the quatrature points
         *
         * @return A map containing the Robin boundary condition information
         *         using a key of the element ID.
         */
        map<int, RobinBCInfoSharedPtr> DisContField2D::v_GetRobinBCInfo(void)
        {
            int i,cnt;
            map<int, RobinBCInfoSharedPtr> returnval;
            Array<OneD, int> ElmtID,EdgeID;
            GetBoundaryToElmtMap(ElmtID,EdgeID);

            for(cnt = i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                MultiRegions::ExpListSharedPtr locExpList;

                if(m_bndConditions[i]->GetBoundaryConditionType() == 
                       SpatialDomains::eRobin)
                {
                    int e,elmtid;
                    Array<OneD, NekDouble> Array_tmp;

                    locExpList = m_bndCondExpansions[i];

                    for(e = 0; e < locExpList->GetExpSize(); ++e)
                    {
                        RobinBCInfoSharedPtr rInfo = MemoryManager<RobinBCInfo>
                            ::AllocateSharedPtr(
                                EdgeID[cnt+e],
                                Array_tmp = locExpList->GetPhys() + 
                                            locExpList->GetPhys_Offset(e));
                        elmtid = ElmtID[cnt+e];
                        // make link list if necessary
                        if(returnval.count(elmtid) != 0)
                        {
                            rInfo->next = returnval.find(elmtid)->second;
                        }
                        returnval[elmtid] = rInfo;
                    }
                }
                cnt += m_bndCondExpansions[i]->GetExpSize();
            }

            return returnval;
        }

        /**
         * @brief Evaluate HDG post-processing to increase polynomial order of
         * solution.
         * 
         * This function takes the solution (assumed to be one order lower) in
         * physical space, and postprocesses at the current polynomial order by
         * solving the system:
         * 
         * \f[
         * \begin{aligned}
         *   (\nabla w, \nabla u^*) &= (\nabla w, u), \\
         *   \langle \nabla u^*, 1 \rangle &= \langle \nabla u, 1 \rangle
         * \end{aligned}
         * \f]
         * 
         * where \f$ u \f$ corresponds with the current solution as stored
         * inside #m_coeffs.
         * 
         * @param outarray  The resulting field \f$ u^* \f$.
         */
        void  DisContField2D::EvaluateHDGPostProcessing(
            Array<OneD, NekDouble> &outarray)
        {
            int    i,cnt,e,ncoeff_edge;
            Array<OneD, NekDouble> force, out_tmp,qrhs;
            Array<OneD, Array< OneD, StdRegions::StdExpansionSharedPtr> > 
                &elmtToTrace = m_traceMap->GetElmtToTrace();

            StdRegions::Orientation edgedir;

            int     eid,nq_elmt, nm_elmt;
            int     LocBndCoeffs = m_traceMap->GetNumLocalBndCoeffs();
            Array<OneD, NekDouble> loc_lambda(LocBndCoeffs), edge_lambda;
            Array<OneD, NekDouble> tmp_coeffs;
            m_traceMap->GlobalToLocalBnd(m_trace->GetCoeffs(),loc_lambda);

            edge_lambda = loc_lambda;

            // Calculate Q using standard DG formulation.
            for(i = cnt = 0; i < GetExpSize(); ++i)
            {
                eid = m_offset_elmt_id[i];

                nq_elmt = (*m_exp)[eid]->GetTotPoints();
                nm_elmt = (*m_exp)[eid]->GetNcoeffs();
                qrhs  = Array<OneD, NekDouble>(nq_elmt);
                force = Array<OneD, NekDouble>(2*nm_elmt);
                out_tmp = force + nm_elmt;

                // Probably a better way of setting up lambda than this.  Note
                // cannot use PutCoeffsInToElmts since lambda space is mapped
                // during the solve.
                for(e = 0; e < (*m_exp)[eid]->GetNedges(); ++e)
                {
                    edgedir = (*m_exp)[eid]->GetEorient(e);

                    ncoeff_edge = elmtToTrace[eid][e]->GetNcoeffs();
                    elmtToTrace[eid][e]->SetCoeffsToOrientation(
                        edgedir,edge_lambda,edge_lambda);
                    Vmath::Vcopy(ncoeff_edge,edge_lambda,1,
                                 elmtToTrace[eid][e]->UpdateCoeffs(),1);
                    edge_lambda = edge_lambda + ncoeff_edge;
                }

                // (d/dx w, d/dx q_0)
                (*m_exp)[eid]->DGDeriv(
                    0,tmp_coeffs = m_coeffs + m_coeff_offset[eid],
                    elmtToTrace[eid], out_tmp);
                (*m_exp)[eid]->BwdTrans(out_tmp,qrhs);
                (*m_exp)[eid]->IProductWRTDerivBase(0,qrhs,force);

                // + (d/dy w, d/dy q_1)
                (*m_exp)[eid]->DGDeriv(
                    1,tmp_coeffs = m_coeffs + m_coeff_offset[eid],
                    elmtToTrace[eid], out_tmp);
                (*m_exp)[eid]->BwdTrans(out_tmp,qrhs);
                (*m_exp)[eid]->IProductWRTDerivBase(1,qrhs,out_tmp);

                Vmath::Vadd(nm_elmt,force,1,out_tmp,1,force,1);

                // determine force[0] = (1,u)
                (*m_exp)[eid]->BwdTrans(
                    tmp_coeffs = m_coeffs + m_coeff_offset[eid],qrhs);
                force[0] = (*m_exp)[eid]->Integral(qrhs);

                // multiply by inverse Laplacian matrix
                // get matrix inverse
                LocalRegions::MatrixKey  lapkey(
                    StdRegions::eInvLaplacianWithUnityMean,  
                    (*m_exp)[eid]->DetExpansionType(), *(*m_exp)[eid]);
                DNekScalMatSharedPtr lapsys = 
                    boost::dynamic_pointer_cast<LocalRegions::Expansion>(
                        (*m_exp)[eid])->GetLocMatrix(lapkey);
                
                NekVector<NekDouble> in (nm_elmt,force,eWrapper);
                NekVector<NekDouble> out(nm_elmt,
                    tmp_coeffs = outarray + m_coeff_offset[eid],eWrapper);

                out = (*lapsys)*in;
            }
        }

        /**
         * Evaluates the boundary condition expansions, \a bndCondExpansions,
         * given the information provided by \a bndConditions.
         * @param   time        The time at which the boundary conditions
         *                      should be evaluated.
         * @param   bndCondExpansions   List of boundary conditions.
         * @param   bndConditions   Information about the boundary conditions.
         */
        void DisContField2D::v_EvaluateBoundaryConditions(const NekDouble time,
                                                          const NekDouble x2_in,
                                                          const NekDouble x3_in)
        {
            int i;
            int npoints;
            int nbnd = m_bndCondExpansions.num_elements();

            MultiRegions::ExpListSharedPtr locExpList;

            for(i = 0; i < nbnd; ++i)
            {
                if(time == 0.0 || m_bndConditions[i]->GetUserDefined() == 
                   SpatialDomains::eTimeDependent)
                {
                    locExpList = m_bndCondExpansions[i];
                    npoints    = locExpList->GetNpoints();
                    Array<OneD,NekDouble> x0(npoints,0.0);
                    Array<OneD,NekDouble> x1(npoints,0.0);
                    Array<OneD,NekDouble> x2(npoints,0.0);

                    // Homogeneous input case for x2.
                    if(x2_in == NekConstants::kNekUnsetDouble)
                    {
                        locExpList->GetCoords(x0,x1,x2);
                    }
                    else
                    {
                        locExpList->GetCoords(x0,x1,x2);
                        Vmath::Fill(npoints,x2_in,x2,1);
                    }

                    if(m_bndConditions[i]->GetBoundaryConditionType()
                       == SpatialDomains::eDirichlet)
                    {

                        string filebcs = boost::static_pointer_cast<
                            SpatialDomains::DirichletBoundaryCondition>(
                                m_bndConditions[i])->m_filename;
                        
                        if(filebcs != "")
                        {
                             string var = filebcs.substr(
                                 0, filebcs.find_last_of("."));
                             int len = var.length();
                             var = var.substr(len-1,len);

                             cout << "Boundary condition from file:" 
                                  << filebcs << endl;

                             std::vector<SpatialDomains::
                                         FieldDefinitionsSharedPtr> FieldDef;
                             std::vector<std::vector<NekDouble> > FieldData;
                             m_graph->Import(filebcs,FieldDef, FieldData);

                             // copy FieldData into locExpList
                             locExpList->ExtractDataToCoeffs(
                                 FieldDef[0], FieldData[0],
                                 FieldDef[0]->m_fields[0], locExpList->UpdateCoeffs());   
                             locExpList->BwdTrans_IterPerExp(
                                 locExpList->GetCoeffs(), 
                                 locExpList->UpdatePhys());
                             locExpList->FwdTrans_BndConstrained(
                                 locExpList->GetPhys(),
                                 locExpList->UpdateCoeffs());
                        }
                        else
                        {
                            LibUtilities::Equation condition = 
                                boost::static_pointer_cast<
                                    SpatialDomains::DirichletBoundaryCondition
                                >(m_bndConditions[i])->m_dirichletCondition;
                            
                            condition.Evaluate(x0,x1,x2,time, 
                                               locExpList->UpdatePhys());

                            locExpList->FwdTrans_BndConstrained(
                                locExpList->GetPhys(),
                                locExpList->UpdateCoeffs());
                        }
                    }
                    else if(m_bndConditions[i]->GetBoundaryConditionType()
                            == SpatialDomains::eNeumann)
                    {
                        string filebcs = boost::static_pointer_cast<
                            SpatialDomains::NeumannBoundaryCondition>(
                                m_bndConditions[i])->m_filename;
                        
                        if(filebcs != "")
                        {
                             string var = filebcs.substr(
                                 0, filebcs.find_last_of("."));
                             int len=var.length();
                             var = var.substr(len-1,len);

                             cout << "Boundary condition from file: "
                                  << filebcs << endl;

                             std::vector<SpatialDomains::
                                         FieldDefinitionsSharedPtr> FieldDef;
                             std::vector<std::vector<NekDouble> > FieldData;
                             m_graph->Import(filebcs,FieldDef, FieldData);

                             // copy FieldData into locExpList
                             locExpList->ExtractDataToCoeffs(
                                 FieldDef[0], FieldData[0],
                                 FieldDef[0]->m_fields[0], locExpList->UpdateCoeffs());
                             locExpList->BwdTrans_IterPerExp(
                                 locExpList->GetCoeffs(), 
                                 locExpList->UpdatePhys());
                             
                             /*
                             Array<OneD, NekDouble> x(locExpList->GetTotPoints(),0.0);
                             Array<OneD, NekDouble> y(locExpList->GetTotPoints(),0.0);
                             locExpList->GetCoords(x,y);
                             for(int i=0; i< locExpList->GetTotPoints(); i++)
                             {
                                 cout<<i<<"     "<<x[i]<<"    "<<y[i]<<"   "
                                     <<locExpList->GetPhys()[i]<<endl;     
                             } 
                             */
                             
                             locExpList->IProductWRTBase(locExpList->GetPhys(),
                                                    locExpList->UpdateCoeffs());
                        }
                        else
                        {
                            LibUtilities::Equation condition =
                                boost::static_pointer_cast<
                                    SpatialDomains::NeumannBoundaryCondition
                                >(m_bndConditions[i])->m_neumannCondition;
                            condition.Evaluate(x0,x1,x2,time, 
                                               locExpList->UpdatePhys());

                            locExpList->IProductWRTBase(locExpList->GetPhys(),
                                                    locExpList->UpdateCoeffs());
                        }
                    }
                    else if(m_bndConditions[i]->GetBoundaryConditionType()
                            == SpatialDomains::eRobin)
                    {
                        string filebcs  =  boost::static_pointer_cast<
                                                SpatialDomains::RobinBoundaryCondition
                                           >(m_bndConditions[i])->m_filename;
                        if(filebcs != "")
                        {
                            //Never tested!!!
                            string var = filebcs.substr(
                                0, filebcs.find_last_of("."));
                            int len = var.length();
                            var = var.substr(len-1,len);

                            std::vector<SpatialDomains::
                                        FieldDefinitionsSharedPtr> FieldDef;
                            std::vector<std::vector<NekDouble> >   FieldData;

                            m_graph->Import(filebcs,FieldDef, FieldData);

                            // copy FieldData into locExpList
                            locExpList->ExtractDataToCoeffs(
                                FieldDef[0], FieldData[0],
                                FieldDef[0]->m_fields[0],locExpList->UpdateCoeffs());
                            locExpList->BwdTrans_IterPerExp(
                                locExpList->GetCoeffs(), 
                                locExpList->UpdatePhys());
                            locExpList->IProductWRTBase(
                                locExpList->GetPhys(),
                                locExpList->UpdateCoeffs());

                            LibUtilities::Equation coeff = 
                                boost::static_pointer_cast<
                                    SpatialDomains::RobinBoundaryCondition
                                >(m_bndConditions[i])->m_robinPrimitiveCoeff;

                            // Array<OneD,NekDouble> timeArray(npoints, time);
                            // put primitive coefficient into the physical space
                            // storage
                            coeff.Evaluate(x0,x1,x2,time, 
                                           locExpList->UpdatePhys());
                        }
                        else
                        {
                            LibUtilities::Equation condition = 
                                boost::static_pointer_cast<
                                    SpatialDomains::RobinBoundaryCondition
                                >(m_bndConditions[i])->m_robinFunction;
                            LibUtilities::Equation coeff     = 
                                boost::static_pointer_cast<
                                    SpatialDomains::RobinBoundaryCondition
                                >(m_bndConditions[i])->m_robinPrimitiveCoeff;
                            condition.Evaluate(x0,x1,x2,time,
                                               locExpList->UpdatePhys());

                            locExpList->IProductWRTBase(
                                locExpList->GetPhys(),
                                locExpList->UpdateCoeffs());

                            // put primitive coefficient into the physical space
                            // storage
                            coeff.Evaluate(x0,x1,x2,time,
                                           locExpList->UpdatePhys());
                        }
                    }    
                    else
                    {
                        ASSERTL0(false,"This type of BC not implemented yet");
                    }
                }
            }
        }
    } // end of namespace
} //end of namespace
