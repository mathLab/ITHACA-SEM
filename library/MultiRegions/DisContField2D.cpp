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
// Description: Field definition for  2D domain with boundary
// conditions using LDG flux
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/DisContField2D.h>
#include <LocalRegions/SegExp.h>

namespace Nektar
{
    namespace MultiRegions
    {

        DisContField2D::DisContField2D(void):
            ExpList2D(),
            m_bndCondExpansions(),
            m_bndConditions()
        {
        }

        DisContField2D::DisContField2D(const DisContField2D &In, const bool DeclareCoeffPhysArrays):
            ExpList2D(In,DeclareCoeffPhysArrays),
            m_bndCondExpansions   (In.m_bndCondExpansions),
            m_bndConditions       (In.m_bndConditions),
            m_globalBndMat        (In.m_globalBndMat),
            m_trace               (In.m_trace),
            m_traceMap            (In.m_traceMap)
        {
        }


        DisContField2D::DisContField2D(const LibUtilities::SessionReaderSharedPtr &pSession,
                                       const SpatialDomains::MeshGraphSharedPtr &graph2D,
                                       const std::string &variable,
                                       const bool SetUpJustDG,
                                       const bool DeclareCoeffPhysArrays):

            ExpList2D(pSession,graph2D,DeclareCoeffPhysArrays,variable),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            SpatialDomains::BoundaryConditions bcs(m_session, graph2D);

            GenerateBoundaryConditionExpansion(graph2D,bcs,variable,
                                               DeclareCoeffPhysArrays);

            if(DeclareCoeffPhysArrays)
            {
                EvaluateBoundaryConditions();
            }

            ApplyGeomInfo();

            if(SetUpJustDG)
            {
                // Set up matrix map
                m_globalBndMat   = MemoryManager<GlobalLinSysMap>::AllocateSharedPtr();

                map<int,int> periodicEdges;
                vector<map<int,int> > periodicVertices;
                GetPeriodicEdges(graph2D,bcs,variable,periodicVertices,periodicEdges);

                // Set up Trace space
                m_trace = MemoryManager<ExpList1D>::AllocateSharedPtr(m_bndCondExpansions,m_bndConditions,*m_exp,graph2D,periodicEdges);

                // Scatter trace segments to 2D elements. For each element,
                // we find the trace segment associated to each edge. The
                // element then retains a pointer to the trace space segments,
                // to ensure uniqueness of normals when retrieving from two
                // adjoining elements which do not lie in a plane.
                SpatialDomains::Geometry1DSharedPtr ElmtSegGeom;
                SpatialDomains::Geometry1DSharedPtr TraceSegGeom;
                for (int i = 0; i < m_exp->size(); ++i)
                {
                    for (int j = 0; j < (*m_exp)[i]->GetNedges(); ++j)
                    {
                        ElmtSegGeom  = ((*m_exp)[i]->GetGeom2D())->GetEdge(j);
                        for (int k = 0; k < m_trace->GetExpSize(); ++k)
                        {
                            TraceSegGeom = m_trace->GetExp(k)->GetGeom1D();
                            if (TraceSegGeom == ElmtSegGeom)
                            {
                                LocalRegions::Expansion2DSharedPtr exp2d
                                    = boost::dynamic_pointer_cast<LocalRegions::Expansion2D>((*m_exp)[i]);
                                LocalRegions::Expansion1DSharedPtr exp1d
                                    = boost::dynamic_pointer_cast<LocalRegions::Expansion1D>(m_trace->GetExp(k));

                                exp2d->SetEdgeExp(j,exp1d);
                                exp1d->SetAdjacentElementExp(j,exp2d);
                                break;
                            }
                        }
                    }
                }

                SetUpPhysNormals();

                m_traceMap = MemoryManager<LocalToGlobalDGMap>::
                    AllocateSharedPtr(m_session, graph2D,m_trace,*this,
                                      m_bndCondExpansions,m_bndConditions, periodicEdges);
            }
            else
            {
                // set elmt edges to point to robin bc edges if required.
                int i,cnt;
                Array<OneD, int> ElmtID,EdgeID;
                GetBoundaryToElmtMap(ElmtID,EdgeID);

                for(cnt = i = 0; i < m_bndCondExpansions.num_elements(); ++i)
                {
                    MultiRegions::ExpListSharedPtr locExpList;

//                    if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eRobin)
//                    {
                        int e;
                        locExpList = m_bndCondExpansions[i];

                        for(e = 0; e < locExpList->GetExpSize(); ++e)
                        {
                            LocalRegions::Expansion2DSharedPtr exp2d
                                = boost::dynamic_pointer_cast<LocalRegions::Expansion2D>((*m_exp)[ElmtID[cnt+e]]);
                            LocalRegions::Expansion1DSharedPtr exp1d
                                = boost::dynamic_pointer_cast<LocalRegions::Expansion1D>(locExpList->GetExp(e));

                            exp2d->SetEdgeExp(EdgeID[cnt+e],exp1d);
                            exp1d->SetAdjacentElementExp(EdgeID[cnt+e],exp2d);
                        }
//                    }
                    cnt += m_bndCondExpansions[i]->GetExpSize();
                }
				
                SetUpPhysNormals();
            }
        }


        // Copy type constructor which declares new boundary conditions
        // and re-uses mapping info and trace space if possible
        DisContField2D::DisContField2D(const DisContField2D &In,
                                       const SpatialDomains::MeshGraphSharedPtr &graph2D,
                                       const std::string &variable,
                                       const bool SetUpJustDG,
                                       const bool DeclareCoeffPhysArrays):
            ExpList2D(In,DeclareCoeffPhysArrays)
        {
            SpatialDomains::BoundaryConditions bcs(m_session, graph2D);
            // Set up boundary conditions for this variable.
            GenerateBoundaryConditionExpansion(graph2D,bcs,variable);

            if(DeclareCoeffPhysArrays)
            {
                EvaluateBoundaryConditions();
            }

            if(!SameTypeOfBoundaryConditions(In))
            {
                if(SetUpJustDG)
                {
                    // Set up matrix map
                    m_globalBndMat = MemoryManager<GlobalLinSysMap>
                        ::AllocateSharedPtr();
                    map<int,int> periodicEdges;
                    vector<map<int,int> >periodicVertices;
                    GetPeriodicEdges(graph2D,bcs,variable,
                                     periodicVertices,periodicEdges);

                    // Set up Trace space
                    m_trace = MemoryManager<ExpList1D>
                        ::AllocateSharedPtr(m_bndCondExpansions, m_bndConditions,
                                            *m_exp,graph2D, periodicEdges);

                    // Scatter trace segments to 2D elements. For each
                    // element, we find the trace segment associated
                    // to each edge. The element then retains a
                    // pointer to the trace space segments, to ensure
                    // uniqueness of normals when retrieving from two
                    // adjoining elements which do not lie in a plane.
                    SpatialDomains::Geometry1DSharedPtr ElmtSegGeom;
                    SpatialDomains::Geometry1DSharedPtr TraceSegGeom;
                    for (int i = 0; i < m_exp->size(); ++i)
                    {
                        for (int j = 0; j < (*m_exp)[i]->GetNedges(); ++j)
                        {
                            ElmtSegGeom  = ((*m_exp)[i]->GetGeom2D())->GetEdge(j);
                            for (int k = 0; k < m_trace->GetExpSize(); ++k)
                            {
                                TraceSegGeom = m_trace->GetExp(k)->GetGeom1D();
                                if (TraceSegGeom == ElmtSegGeom)
                                {
                                    LocalRegions::Expansion2DSharedPtr exp2d
                                        = boost::dynamic_pointer_cast<LocalRegions::Expansion2D>((*m_exp)[i]);
                                    LocalRegions::Expansion1DSharedPtr exp1d
                                        = boost::dynamic_pointer_cast<LocalRegions::Expansion1D>(m_trace->GetExp(k));

                                    exp2d->SetEdgeExp(j,exp1d);
                                    exp1d->SetAdjacentElementExp(j,exp2d);
                                    break;
                                }
                            }
                        }
                    }

                    SetUpPhysNormals();

                    // Finally set up the trace map between element edges and
                    // trace segment expansions.
                    m_traceMap = MemoryManager<LocalToGlobalDGMap>::
                        AllocateSharedPtr(m_session,graph2D,m_trace,*this,
                                          m_bndCondExpansions,m_bndConditions,
                                          periodicEdges);

                }
                else
                {
                    // set elmt edges to point to robin bc edges if required.
                    int i,cnt;
                    Array<OneD, int> ElmtID,EdgeID;
                    GetBoundaryToElmtMap(ElmtID,EdgeID);

                    for(cnt = i = 0; i < m_bndCondExpansions.num_elements(); ++i)
                    {
                        MultiRegions::ExpListSharedPtr locExpList;

//                        if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eRobin)
//                        {
                            int e;
                            locExpList = m_bndCondExpansions[i];

                            for(e = 0; e < locExpList->GetExpSize(); ++e)
                            {
                                LocalRegions::Expansion2DSharedPtr exp2d
                                    = boost::dynamic_pointer_cast<LocalRegions::Expansion2D>((*m_exp)[ElmtID[cnt+e]]);
                                LocalRegions::Expansion1DSharedPtr exp1d
                                    = boost::dynamic_pointer_cast<LocalRegions::Expansion1D>(locExpList->GetExp(e));

                                exp2d->SetEdgeExp(EdgeID[cnt+e],exp1d);
                                exp1d->SetAdjacentElementExp(EdgeID[cnt+e],exp2d);
                            }
//                        }
                        cnt += m_bndCondExpansions[i]->GetExpSize();
                    }

                    SetUpPhysNormals();
                }
            }
            else
            {
                if(SetUpJustDG)
                {
                    m_globalBndMat = In.m_globalBndMat;
                    m_trace        = In.m_trace;
                    m_traceMap     = In.m_traceMap;
                }
				else 
				{
					m_globalBndMat = In.m_globalBndMat;
                    m_trace        = In.m_trace;
                    m_traceMap     = In.m_traceMap;
					
					// set elmt edges to point to robin bc edges if required.
                    int i,cnt;
                    Array<OneD, int> ElmtID,EdgeID;
                    GetBoundaryToElmtMap(ElmtID,EdgeID);
					
                    for(cnt = i = 0; i < m_bndCondExpansions.num_elements(); ++i)
                    {
                        MultiRegions::ExpListSharedPtr locExpList;
						
						//                        if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eRobin)
						//                        {
						int e;
						locExpList = m_bndCondExpansions[i];
						
						for(e = 0; e < locExpList->GetExpSize(); ++e)
						{
							LocalRegions::Expansion2DSharedPtr exp2d
							= boost::dynamic_pointer_cast<LocalRegions::Expansion2D>((*m_exp)[ElmtID[cnt+e]]);
							LocalRegions::Expansion1DSharedPtr exp1d
							= boost::dynamic_pointer_cast<LocalRegions::Expansion1D>(locExpList->GetExp(e));
							
							exp2d->SetEdgeExp(EdgeID[cnt+e],exp1d);
							exp1d->SetAdjacentElementExp(EdgeID[cnt+e],exp2d);
						}
						//                        }
                        cnt += m_bndCondExpansions[i]->GetExpSize();
                    }

                    SetUpPhysNormals();
				}

            }
        }


        /**
         * For each boundary region, checks that the types and number of
         * boundary expansions in that region match.
         * @param   In          ContField2D to compare with.
         * @returns True if boundary conditions match.
         */
        bool DisContField2D::SameTypeOfBoundaryConditions(const DisContField2D &In)
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
            m_comm->AllReduce(vSame, LibUtilities::ReduceMin);

            return (vSame == 1);
        }


        void DisContField2D::GenerateBoundaryConditionExpansion(const SpatialDomains::MeshGraphSharedPtr &graph2D,
                                                                const SpatialDomains::BoundaryConditions &bcs,
                                                                const std::string &variable,
                                                                const bool DeclareCoeffPhysArrays)
        {  	
            int i,cnt  = 0;
            SpatialDomains::BoundaryConditionShPtr locBCond;
            MultiRegions::ExpList1DSharedPtr       locExpList;
            const SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();

            int nbnd = bregions.size();
            // count the number of non-periodic boundary regions
            for(i = 0; i < nbnd; ++i)
            {
                SpatialDomains::BoundaryConditionShPtr boundaryCondition = GetBoundaryCondition(bconditions, i, variable);
                if( boundaryCondition->GetBoundaryConditionType() != SpatialDomains::ePeriodic )
                {
                    cnt++;
                }              
            }

            m_bndCondExpansions  = Array<OneD,MultiRegions::ExpListSharedPtr>(cnt);
            m_bndConditions      = Array<OneD,SpatialDomains::BoundaryConditionShPtr>(cnt);
	    
            cnt=0;

            // list non-periodic boundaries
            for(i = 0; i < nbnd; ++i)
            {
                locBCond = GetBoundaryCondition(bconditions, i, variable);
                if(locBCond->GetBoundaryConditionType() != SpatialDomains::ePeriodic)
                {
                    locExpList = MemoryManager<MultiRegions::ExpList1D>
                        ::AllocateSharedPtr(*(bregions[i]),
                                            graph2D,
                                            DeclareCoeffPhysArrays);
                    

                    // Set up normals on non-Dirichlet boundary conditions
                    if(locBCond->GetBoundaryConditionType()
                       != SpatialDomains::eDirichlet)
                    {
                        SetUpPhysNormals();
                    }

                    m_bndCondExpansions[cnt]  = locExpList;
                    m_bndConditions[cnt]    = locBCond;
                    SpatialDomains::BndUserDefinedType type = m_bndConditions[cnt++]->GetUserDefined();
                    if((type == SpatialDomains::eI)||(type == SpatialDomains::eCalcBC))
                    {
                        locExpList->SetUpPhysTangents(*m_exp);
                        SetUpPhysNormals();
                    }
                }
            }
        }

        /**
         * @param   graph2D     A mesh containing information about the domain
         *                      and the spectral/hp element expansion.
         * @param   bcs         Information about the boundary conditions.
         * @param   variable    Specifies the field.
         * @param   periodicVerts   Vector of Maps into which the list of
         *                      periodic vertices is placed, one map for each
         *                      boundary region.
         * @param   periodicEdges   Map into which the list of periodic
         *                      edges is placed.
         */
        void DisContField2D::GetPeriodicEdges(
                                              const SpatialDomains::MeshGraphSharedPtr &graph2D,
                                              const SpatialDomains::BoundaryConditions &bcs,
                                              const std::string &variable,
                                              vector<map<int,int> >& periodicVerts,
                                              map<int,int>& periodicEdges)
        {
            ASSERTL0(boost::dynamic_pointer_cast<SpatialDomains::MeshGraph2D>(graph2D),
                     "Expected a MeshGraph2D.");
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

            StdRegions::EdgeOrientation orient1;
            StdRegions::EdgeOrientation orient2;

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

                                // Extract the periodic edges
                                periodicEdges[segmentGeom1->GetEid()]
                                    = segmentGeom2->GetEid();
                                periodicEdges[segmentGeom2->GetEid()]
                                    = segmentGeom1->GetEid();

                                // Extract the periodic vertices
                                element1 = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph2D>(graph2D)
                                    ->GetElementsFromEdge(segmentGeom1);
                                element2 = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph2D>(graph2D)
                                    ->GetElementsFromEdge(segmentGeom2);

                                ASSERTL0(element1->size()==1,
                                         "The periodic boundaries belong to "
                                         "more than one element of the mesh");
                                ASSERTL0(element2->size()==1,
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
                                    periodicVertices[segmentGeom1->GetVid(0)]
                                        = segmentGeom2->GetVid(0);
                                    periodicVertices[segmentGeom1->GetVid(1)]
                                        = segmentGeom2->GetVid(1);
                                }
                                else
                                {
                                    periodicVertices[segmentGeom1->GetVid(0)]
                                        = segmentGeom2->GetVid(1);
                                    periodicVertices[segmentGeom1->GetVid(1)]
                                        = segmentGeom2->GetVid(0);
                                }
                            }
                        }
                        periodicVerts.push_back(periodicVertices);
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


        DisContField2D::~DisContField2D()
        {
        }

        GlobalLinSysSharedPtr DisContField2D::GetGlobalBndLinSys(const GlobalLinSysKey &mkey)
        {
            ASSERTL0(mkey.GetMatrixType() == StdRegions::eHybridDGHelmBndLam,
                     "Routine currently only tested for HybridDGHelmholtz");
            ASSERTL1(mkey.GetGlobalSysSolnType()==m_traceMap->GetGlobalSysSolnType(),
                     "The local to global map is not set up for the requested solution type");

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

        // Construct the two trace vectors of the inner and outer
        // trace solution from the field contained in m_phys, where
        // the Weak dirichlet boundary conditions are listed in the
        // outer part of the vecotr
        void DisContField2D::GetFwdBwdTracePhys(Array<OneD,NekDouble> &Fwd,
                                                Array<OneD,NekDouble> &Bwd)
        {
            GetFwdBwdTracePhys(m_phys,Fwd,Bwd);
        }

        void DisContField2D::GetFwdBwdTracePhys(const Array<OneD,const NekDouble>  &field,
                                                Array<OneD,NekDouble> &Fwd,
                                                Array<OneD,NekDouble> &Bwd)
        {
            // Loop over elements and collect forward expansion
            int nexp = GetExpSize();
            StdRegions::EdgeOrientation edgedir;
            int nquad_e,cnt,n,e,npts,offset, phys_offset;
            Array<OneD,NekDouble> e_tmp;

            Array<OneD, Array<OneD, StdRegions::StdExpansion1DSharedPtr> >
                elmtToTrace = m_traceMap->GetElmtToTrace();

            // zero vectors;
            Vmath::Zero(Fwd.num_elements(),Fwd,1);
            Vmath::Zero(Bwd.num_elements(),Bwd,1);

            for(n  = 0; n < nexp; ++n)
            {
                phys_offset = GetPhys_Offset(n);

                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    nquad_e = (*m_exp)[n]->GetEdgeNumPoints(e);
                    edgedir = (*m_exp)[n]->GetEorient(e);
                    if(edgedir == StdRegions::eForwards)
                    {
                        offset = m_trace->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());
                        (*m_exp)[n]->GetEdgePhysVals(e, elmtToTrace[n][e],
                                                     field + phys_offset,
                                                     e_tmp = Fwd + offset);
                    }
                }
            }

            for(n  = 0; n < nexp; ++n)
            {
                phys_offset = GetPhys_Offset(n);

                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    nquad_e = (*m_exp)[n]->GetEdgeNumPoints(e);
                    edgedir = (*m_exp)[n]->GetEorient(e);
                    if(edgedir == StdRegions::eBackwards)
                    {
                        offset = m_trace->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());
                        (*m_exp)[n]->GetEdgePhysVals(e, elmtToTrace[n][e],
                                                     field + phys_offset,
                                                     e_tmp = Bwd + offset);
                    }
                }
            }
			
            // fill boundary conditions into missing elements
            int id1,id2 = 0;
            cnt = 0;
			
            for(n = 0; n < m_bndCondExpansions.num_elements(); ++n)
            {				
                if(m_bndConditions[n]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    for(e = 0; e < m_bndCondExpansions[n]->GetExpSize(); ++e)
                    {
						
                        npts = m_bndCondExpansions[n]->GetExp(e)->GetNumPoints(0);

                        if(m_traceMap->GetBndExpAdjacentOrient(cnt+e) == eAdjacentEdgeIsForwards)
                        {
                            id1 = m_bndCondExpansions[n]->GetPhys_Offset(e) ;
                            id2 = m_trace->GetPhys_Offset(m_traceMap->GetBndCondTraceToGlobalTraceMap(cnt+e));
                            Vmath::Vcopy(npts,&(m_bndCondExpansions[n]->GetPhys())[id1],1,&Bwd[id2],1);
                        }
                        else
                        {
                            id1 = m_bndCondExpansions[n]->GetPhys_Offset(e) ;
                            id2 = m_trace->GetPhys_Offset(m_traceMap->GetBndCondTraceToGlobalTraceMap(cnt+e));
                            Vmath::Vcopy(npts,&(m_bndCondExpansions[n]->GetPhys())[id1],1,&Fwd[id2],1);
                        }
                    }

                    cnt +=e;
                }
                else if((m_bndConditions[n]->GetBoundaryConditionType() == SpatialDomains::eNeumann)||(m_bndConditions[n]->GetBoundaryConditionType() == SpatialDomains::eRobin))
                {
                    for(e = 0; e < m_bndCondExpansions[n]->GetExpSize(); ++e)
                    {
                        npts = m_bndCondExpansions[n]->GetExp(e)->GetNumPoints(0);

                        if(m_traceMap->GetBndExpAdjacentOrient(cnt+e) == eAdjacentEdgeIsForwards)
                        {
                            id1 = m_bndCondExpansions[n]->GetPhys_Offset(e);
                            ASSERTL0((m_bndCondExpansions[n]->GetPhys())[id1] == 0.0,"method not set up for non-zero Neumann boundary condition");
                            id2 = m_trace->GetPhys_Offset(m_traceMap->GetBndCondTraceToGlobalTraceMap(cnt+e));
                            Vmath::Vcopy(npts,&Fwd[id2],1,&Bwd[id2],1);
                        }
                        else
                        {
                            id1 = m_bndCondExpansions[n]->GetPhys_Offset(e);
                            ASSERTL0((m_bndCondExpansions[n]->GetPhys())[id1] == 0.0,"method not set up for non-zero Neumann boundary condition");
                            id2 = m_trace->GetPhys_Offset(m_traceMap->GetBndCondTraceToGlobalTraceMap(cnt+e));
                            Vmath::Vcopy(npts,&Bwd[id2],1,&Fwd[id2],1);
                        }
                    }

                    cnt +=e;
                }
                else
                {
                    ASSERTL0(false,"method not set up for non-Dirichlet conditions");
                }
            }

        }

        void DisContField2D::ExtractTracePhys(Array<OneD,NekDouble> &outarray)
        {

            ASSERTL1(m_physState == true,
                     "local physical space is not true ");

            ExtractTracePhys(m_phys, outarray);
        }

        void DisContField2D::ExtractTracePhys(const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
        {
            // Loop over elemente and collect forward expansion
            int nexp = GetExpSize();
            int nquad_e,n,e,offset,phys_offset;
            Array<OneD,NekDouble> e_tmp;
            Array<OneD, Array<OneD, StdRegions::StdExpansion1DSharedPtr> >
                elmtToTrace = m_traceMap->GetElmtToTrace();

            ASSERTL1(outarray.num_elements() >= m_trace->GetNpoints(),
                     "input array is of insufficient length");

            // use m_trace tmp space in element to fill values
            for(n  = 0; n < nexp; ++n)
            {
				phys_offset = GetPhys_Offset(n);
				
                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    nquad_e = (*m_exp)[n]->GetEdgeNumPoints(e);
                    offset = m_trace->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());
					(*m_exp)[n]->GetEdgePhysVals(e,  elmtToTrace[n][e],
                                                 inarray + phys_offset,
                                                 e_tmp = outarray + offset);
                }
            }
        }

        /// Note this routine changes m_trace->m_coeffs space;
        void DisContField2D::AddTraceIntegral(const Array<OneD, const NekDouble> &Fx,
                                              const Array<OneD, const NekDouble> &Fy,
                                              Array<OneD, NekDouble> &outarray)
        {
            int e,n,offset, t_offset;
            Array<OneD, NekDouble> e_outarray;
            Array<OneD, Array<OneD, StdRegions::StdExpansion1DSharedPtr> >
                elmtToTrace = m_traceMap->GetElmtToTrace();

            for(n = 0; n < GetExpSize(); ++n)
            {
                offset = GetCoeff_Offset(n);
                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    t_offset = GetTrace()->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());

                    (*m_exp)[n]->AddEdgeNormBoundaryInt(e,elmtToTrace[n][e],
                                                        Fx + t_offset,
                                                        Fy + t_offset,
                                                        e_outarray = outarray+offset);
                }
            }
        }

        /// Note this routine changes m_trace->m_coeffs space;
        void DisContField2D::AddTraceIntegral(const Array<OneD, const NekDouble> &Fn, Array<OneD, NekDouble> &outarray)
        {
			
            int e,n,offset, t_offset;
            Array<OneD, NekDouble> e_outarray;
            Array<OneD, Array<OneD, StdRegions::StdExpansion1DSharedPtr> >
                elmtToTrace = m_traceMap->GetElmtToTrace();

            for(n = 0; n < GetExpSize(); ++n)
            {
                offset = GetCoeff_Offset(n);
                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    t_offset = GetTrace()->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());
                    (*m_exp)[n]->AddEdgeNormBoundaryInt(e,elmtToTrace[n][e],
                                                        Fn + t_offset,
                                                        e_outarray = outarray+offset);
                }
            }
        }

        // Set up a list of element ids and edge ids that link to the
        // boundary conditions
        void DisContField2D::GetBoundaryToElmtMap(Array<OneD, int> &ElmtID, Array<OneD,int> &EdgeID)
        {
            map<int, int> globalIdMap;
            int i,n,id;
            int bid,cnt,Eid;
            int nbcs = 0;

            SpatialDomains::MeshGraph2DSharedPtr graph2D =
                boost::dynamic_pointer_cast<SpatialDomains::MeshGraph2D>(m_graph);

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
                    // Use face to element map from MeshGraph3D.
                    SpatialDomains::ElementEdgeVectorSharedPtr tmp =
                        graph2D->GetElementsFromEdge(
                            m_bndCondExpansions[n]->GetExp(i)->GetGeom1D());

                    ElmtID[cnt] = globalIdMap[(*tmp)[0]->m_Element->GetGlobalID()];
                    EdgeID[cnt] = (*tmp)[0]->m_EdgeIndx;
                }
            }
        }

        /// Note this routine changes m_trace->m_coeffs space;
        void DisContField2D::AddTraceBiIntegral(const Array<OneD, const NekDouble> &Fwd,
                                                const Array<OneD, const NekDouble> &Bwd,
                                                Array<OneD, NekDouble> &outarray)
        {
            int e,n,offset, t_offset;
            Array<OneD, NekDouble> e_outarray;
            Array<OneD, Array<OneD, StdRegions::StdExpansion1DSharedPtr> >
                elmtToTrace = m_traceMap->GetElmtToTrace();

            for(n = 0; n < GetExpSize(); ++n)
            {
                offset = GetCoeff_Offset(n);
                for(e = 0; e < (*m_exp)[n]->GetNedges(); ++e)
                {
                    t_offset = GetTrace()->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());

                    (*m_exp)[n]->AddEdgeNormBoundaryBiInt(e,elmtToTrace[n][e],
                                                          Fwd + t_offset,
                                                          Bwd + t_offset,
                                                          e_outarray = outarray+offset);
                }
            }
        }

        /** Calculate the L2 error of the Q_dir derivative using the
            consistent DG evaluation of Q_dir. The soln provided is of
            the primative variation at the quadrature points and the
            derivative is compared to the discrete derivative at these
            points which is likely to be undesireable unless using a
            much higher number of quadrature points than the
            polynomial order used to evaluate Q_dir
        */

        NekDouble DisContField2D::L2_DGDeriv(const int dir,
                                             const Array<OneD, const NekDouble> &soln)
        {

            int    i,e,ncoeff_edge;
            Array<OneD, const NekDouble> tmp_coeffs;
            Array<OneD, NekDouble> out_d(m_ncoeffs), out_tmp;

            Array<OneD, Array< OneD, StdRegions::StdExpansion1DSharedPtr> > elmtToTrace = m_traceMap->GetElmtToTrace();

            StdRegions::EdgeOrientation edgedir;

            int     eid,cnt;
            int     LocBndCoeffs = m_traceMap->GetNumLocalBndCoeffs();
            Array<OneD, NekDouble> loc_lambda(LocBndCoeffs), edge_lambda;
            m_traceMap->GlobalToLocalBnd(m_trace->GetCoeffs(),loc_lambda);

            edge_lambda = loc_lambda;
            // Calculate Q using standard DG formulation.
            for(i =cnt = 0; i < GetExpSize(); ++i)
            {
                eid = m_offset_elmt_id[i];
                // Probably a better way of setting up lambda than this.
                // Note cannot use PutCoeffsInToElmts since lambda space
                // is mapped during the solve.
                for(e = 0; e < (*m_exp)[eid]->GetNedges(); ++e)
                {
                    edgedir = (*m_exp)[eid]->GetEorient(e);

                    ncoeff_edge = elmtToTrace[eid][e]->GetNcoeffs();
                    elmtToTrace[eid][e]->SetCoeffsToOrientation(edgedir,edge_lambda,edge_lambda);
                    Vmath::Vcopy(ncoeff_edge,edge_lambda,1,elmtToTrace[eid][e]->UpdateCoeffs(),1);
                    edge_lambda = edge_lambda + ncoeff_edge;
                }

                (*m_exp)[eid]->DGDeriv(dir,tmp_coeffs = m_coeffs+m_coeff_offset[eid],
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
            GlobalMatrixKey HDGLamToUKey(StdRegions::eHybridDGLamToU,NullLocalToGlobalBaseMapSharedPtr,factors,varcoeff);
            const DNekScalBlkMatSharedPtr &HDGLamToU = GetBlockMatrix(HDGLamToUKey);

            // Retrieve global trace space storage, \Lambda, from trace expansion
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
                if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
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
            GlobalMatrixKey invHDGhelmkey(StdRegions::eInvHybridDGHelmholtz,NullLocalToGlobalBaseMapSharedPtr,factors,varcoeff);
            const DNekScalBlkMatSharedPtr& InvHDGHelm = GetBlockMatrix(invHDGhelmkey);
            DNekVec out(m_ncoeffs,outarray,eWrapper);
            Vmath::Zero(m_ncoeffs,outarray,1);

            // get local trace solution from BndSol
            m_traceMap->GlobalToLocalBnd(BndSol,loc_lambda);

            //  out =  u_f + u_lam = (*InvHDGHelm)*f + (LamtoU)*Lam
            out = (*InvHDGHelm)*F + (*HDGLamToU)*LocLambda;
        }


        void DisContField2D::v_GeneralMatrixOp(
               const GlobalMatrixKey             &gkey,
               const Array<OneD,const NekDouble> &inarray,
                     Array<OneD,      NekDouble> &outarray,
               bool UseContCoeffs)
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
         * Search through the edge expansions and identify which ones
         * have Robin/Mixed type boundary conditions. If find a Robin
         * boundary then store the edge id of the boundary condition
         * and the array of points of the physical space boundary
         * condition which are hold the boundary condition primitive
         * variable coefficient at the quatrature points
         *
         * \return std map containing the robin boundary condition
         * info using a key of the element id
         *
         * There is a next member to allow for more than one Robin
         * boundary condition per element
         */

        map<int, RobinBCInfoSharedPtr> DisContField2D::GetRobinBCInfo(void)
        {
            int i,cnt;
            map<int, RobinBCInfoSharedPtr> returnval;
            Array<OneD, int> ElmtID,EdgeID;
            GetBoundaryToElmtMap(ElmtID,EdgeID);

            for(cnt = i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                MultiRegions::ExpListSharedPtr locExpList;

                if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eRobin)
                {
                    int e,elmtid;
                    Array<OneD, NekDouble> Array_tmp;

                    locExpList = m_bndCondExpansions[i];

                    for(e = 0; e < locExpList->GetExpSize(); ++e)
                    {
                        RobinBCInfoSharedPtr rInfo = MemoryManager<RobinBCInfo>::AllocateSharedPtr(EdgeID[cnt+e],Array_tmp = locExpList->GetPhys() + locExpList->GetPhys_Offset(e));
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

        //Taking the solution (assumed to be one order lower) in
        //physical space in inarray, postprocess at the current
        //polynomial order by solving:
        //
        // (\Grad w, \Grad u*) = (\Grad w, m_coeffs);
        // <1,\Grad u*> = <1,\Grad m_coeffs>
        //

        void  DisContField2D::EvaluateHDGPostProcessing(Array<OneD, NekDouble> &outarray)
        {

            int    i,cnt,e,ncoeff_edge;
            Array<OneD, NekDouble> force, out_tmp,qrhs;
            Array<OneD, Array< OneD, StdRegions::StdExpansion1DSharedPtr> > elmtToTrace = m_traceMap->GetElmtToTrace();

            StdRegions::EdgeOrientation edgedir;

            int     eid,nq_elmt, nm_elmt;
            int     LocBndCoeffs = m_traceMap->GetNumLocalBndCoeffs();
            Array<OneD, NekDouble> loc_lambda(LocBndCoeffs), edge_lambda, tmp_coeffs;
            m_traceMap->GlobalToLocalBnd(m_trace->GetCoeffs(),loc_lambda);

            edge_lambda = loc_lambda;

            // Calculate Q using standard DG formulation.
            for(i =cnt = 0; i < GetExpSize(); ++i)
            {
                eid = m_offset_elmt_id[i];

                nq_elmt = (*m_exp)[eid]->GetTotPoints();
                nm_elmt = (*m_exp)[eid]->GetNcoeffs();
                qrhs  = Array<OneD, NekDouble>(nq_elmt);
                force = Array<OneD, NekDouble>(2*nm_elmt);
                out_tmp = force + nm_elmt;


                // Probably a better way of setting up lambda than this.
                // Note cannot use PutCoeffsInToElmts since lambda space
                // is mapped during the solve.
                for(e = 0; e < (*m_exp)[eid]->GetNedges(); ++e)
                {
                    edgedir = (*m_exp)[eid]->GetEorient(e);

                    ncoeff_edge = elmtToTrace[eid][e]->GetNcoeffs();
                    elmtToTrace[eid][e]->SetCoeffsToOrientation(edgedir,edge_lambda,edge_lambda);
                    Vmath::Vcopy(ncoeff_edge,edge_lambda,1,
                                 elmtToTrace[eid][e]->UpdateCoeffs(),1);
                    edge_lambda = edge_lambda + ncoeff_edge;
                }

                // (d/dx w, d/dx q_0)
                (*m_exp)[eid]->DGDeriv(0,tmp_coeffs = m_coeffs + m_coeff_offset[eid],
                                       elmtToTrace[eid],
                                       out_tmp);
                (*m_exp)[eid]->BwdTrans(out_tmp,qrhs);
                (*m_exp)[eid]->IProductWRTDerivBase(0,qrhs,force);

                // + (d/dy w, d/dy q_1)
                (*m_exp)[eid]->DGDeriv(1,tmp_coeffs = m_coeffs + m_coeff_offset[eid],
                                       elmtToTrace[eid],
                                       out_tmp);
                (*m_exp)[eid]->BwdTrans(out_tmp,qrhs);
                (*m_exp)[eid]->IProductWRTDerivBase(1,qrhs,out_tmp);

                Vmath::Vadd(nm_elmt,force,1,out_tmp,1,force,1);

                // determine force[0] = (1,u)
                (*m_exp)[eid]->BwdTrans(tmp_coeffs = m_coeffs + m_coeff_offset[eid],qrhs);
                force[0] = (*m_exp)[eid]->Integral(qrhs);

                // multiply by inverse Laplacian matrix
                // get matrix inverse
                LocalRegions::MatrixKey  lapkey(StdRegions::eInvLaplacianWithUnityMean,  (*m_exp)[eid]->DetExpansionType(),  *(*m_exp)[eid]);
                DNekScalMatSharedPtr lapsys = boost::dynamic_pointer_cast<LocalRegions::Expansion>((*m_exp)[eid])->GetLocMatrix(lapkey);

                NekVector<NekDouble> in(nm_elmt,force,eWrapper);
                NekVector<NekDouble> out(nm_elmt,tmp_coeffs = outarray + m_coeff_offset[eid],eWrapper);

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
            int i,j;
            int npoints;
            int nbnd = m_bndCondExpansions.num_elements();

            MultiRegions::ExpListSharedPtr locExpList;

            for(i = 0; i < nbnd; ++i)
            {
                if(time == 0.0 || m_bndConditions[i]->GetUserDefined()==SpatialDomains::eTimeDependent)
                {
                    locExpList = m_bndCondExpansions[i];
                    npoints = locExpList->GetNpoints();
                    Array<OneD,NekDouble> x0(npoints,0.0);
                    Array<OneD,NekDouble> x1(npoints,0.0);
                    Array<OneD,NekDouble> x2(npoints,0.0);

                    if(x2_in == NekConstants::kNekUnsetDouble) //homogeneous input case for x2
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

                        string filebcs  =  boost::static_pointer_cast<
                                                SpatialDomains::DirichletBoundaryCondition
                                            >(m_bndConditions[i])->m_filename;
                        if(filebcs != "")
                        {
                             string var = filebcs.substr(0, filebcs.find_last_of("."));
                             int len=var.length();
                             var = var.substr(len-1,len);

                             cout<<"boundary condition from file:"<<filebcs<<endl;

                             std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef;
                             std::vector<std::vector<NekDouble> > FieldData;
                             m_graph->Import(filebcs,FieldDef, FieldData);

                             // copy FieldData into locExpList
                             locExpList->ExtractDataToCoeffs(FieldDef[0], FieldData[0],
                                                 FieldDef[0]->m_fields[0]);   
                             locExpList->BwdTrans_IterPerExp(locExpList->GetCoeffs(), locExpList->UpdatePhys());
                             locExpList->FwdTrans_BndConstrained(locExpList->GetPhys(),
                                        locExpList->UpdateCoeffs());
                        }
                        else
                        {
                            LibUtilities::Equation  condition = boost::static_pointer_cast<
                                                                   SpatialDomains::DirichletBoundaryCondition
                                                                >(m_bndConditions[i])->m_dirichletCondition;
                            condition.Evaluate(x0,x1,x2,time, locExpList->UpdatePhys());

                            locExpList->FwdTrans_BndConstrained(locExpList->GetPhys(),
                                        locExpList->UpdateCoeffs());
                        }
                    }
                    else if(m_bndConditions[i]->GetBoundaryConditionType()
                            == SpatialDomains::eNeumann)
                    {
                        string filebcs  =  boost::static_pointer_cast<
                                                SpatialDomains::NeumannBoundaryCondition
                                           >(m_bndConditions[i])->m_filename;
                        if(filebcs != "")
                        {
                             string var = filebcs.substr(0, filebcs.find_last_of("."));
                             int len=var.length();
                             var = var.substr(len-1,len);

                             cout<<"boundary condition from file:"<<filebcs<<endl;

                             std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef;
                             std::vector<std::vector<NekDouble> > FieldData;
                             m_graph->Import(filebcs,FieldDef, FieldData);

                             // copy FieldData into locExpList
                             locExpList->ExtractDataToCoeffs(FieldDef[0], FieldData[0],
                                             FieldDef[0]->m_fields[0]);
                             locExpList->BwdTrans_IterPerExp(locExpList->GetCoeffs(), locExpList->UpdatePhys());
                             locExpList->IProductWRTBase(locExpList->GetPhys(),
                                                    locExpList->UpdateCoeffs());
                        }
                        else
                        {
                            LibUtilities::Equation  condition = boost::static_pointer_cast<
                                                                   SpatialDomains::NeumannBoundaryCondition
                                                                >(m_bndConditions[i])->m_neumannCondition;
                            condition.Evaluate(x0,x1,x2,time, locExpList->UpdatePhys());

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
                             string var = filebcs.substr(0, filebcs.find_last_of("."));
                             int len=var.length();
                             var = var.substr(len-1,len);

                             std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef;
                             std::vector<std::vector<NekDouble> > FieldData;

                             m_graph->Import(filebcs,FieldDef, FieldData);

                             // copy FieldData into locExpList
                             locExpList->ExtractDataToCoeffs(FieldDef[0], FieldData[0],
                                                 FieldDef[0]->m_fields[0]);
                             locExpList->BwdTrans_IterPerExp(locExpList->GetCoeffs(), locExpList->UpdatePhys());
                             locExpList->IProductWRTBase(locExpList->GetPhys(),
                                                    locExpList->UpdateCoeffs());


                            LibUtilities::Equation  coeff     = boost::static_pointer_cast<
                                                                   SpatialDomains::RobinBoundaryCondition
                                                                >(m_bndConditions[i])->m_robinPrimitiveCoeff;
                            // Array<OneD,NekDouble> timeArray(npoints, time);
                            // put primitive coefficient into the physical space storage
                            coeff.Evaluate(x0,x1,x2,time, locExpList->UpdatePhys());
                        }
                        else
                        {
                            LibUtilities::Equation  condition = boost::static_pointer_cast<
                                                                   SpatialDomains::RobinBoundaryCondition
                                                                >(m_bndConditions[i])->m_robinFunction;
                            LibUtilities::Equation  coeff     = boost::static_pointer_cast<
                                                                   SpatialDomains::RobinBoundaryCondition
                                                                >(m_bndConditions[i])->m_robinPrimitiveCoeff;
                            condition.Evaluate(x0,x1,x2,time, locExpList->UpdatePhys());

                            locExpList->IProductWRTBase(locExpList->GetPhys(),
                                                    locExpList->UpdateCoeffs());

                            // put primitive coefficient into the physical space storage
                            coeff.Evaluate(x0,x1,x2,time, locExpList->UpdatePhys());

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
