///////////////////////////////////////////////////////////////////////////////
//
// File DisContField1D.cpp
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
// Description: Field definition for 1D domain with boundary
// conditions using LDG-H
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/DisContField1D.h>
#include <StdRegions/StdSegExp.h>


namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class DisContField1D
         * This class augments the list of local expansions inherited from
         * ExpList1D with boundary conditions. Inter-element boundaries are
         * handled using an discontinuous Galerkin scheme.
         */

        /**
         * Constructs an empty expansion list with no boundary conditions.
         */
        DisContField1D::DisContField1D():
            ExpList1D(),
            m_bndCondExpansions(),
            m_bndConditions()
        {
        }


        /**
         * An expansion list for the boundary expansions is generated first for
         * the field. These are subsequently evaluated for time zero. The trace
         * map is then constructed.
         * @param   graph1D     A mesh, containing information about the domain
         *                      and the spectral/hp element expansions.
         * @param   bcs         Information about the enforced boundary
         *                      conditions.
         * @param   variable    The session variable associated with the
         *                      boundary conditions to enforce.
         * @param   solnType    Type of global system to use.
         */
        DisContField1D::DisContField1D(const LibUtilities::SessionReaderSharedPtr &pSession,
                                       const SpatialDomains::MeshGraphSharedPtr &graph1D,
                                       const std::string &variable):
            ExpList1D(pSession,graph1D),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            SpatialDomains::BoundaryConditions bcs(m_session, graph1D);

            GenerateBoundaryConditionExpansion(graph1D,bcs,variable);
            EvaluateBoundaryConditions();
            ApplyGeomInfo();

            map<int,int> periodicVertices;
            GetPeriodicVertices(graph1D,bcs,variable,periodicVertices);

            m_globalBndMat
                        = MemoryManager<GlobalLinSysMap>::AllocateSharedPtr();

            ExpList0DSharedPtr trace = MemoryManager<ExpList0D>::AllocateSharedPtr(
                m_bndCondExpansions,m_bndConditions,*m_exp,graph1D,periodicVertices);
            m_trace = boost::dynamic_pointer_cast<ExpList>(trace);

            m_traceMap = MemoryManager<AssemblyMapDG>::
                AllocateSharedPtr(pSession,graph1D,trace,*this,
                                  m_bndCondExpansions,m_bndConditions,periodicVertices);
            
            tmpBndSol = Array<OneD,NekDouble>(m_traceMap->GetNumLocalBndCoeffs());
            
            // Scatter trace points to 1D elements. For each element, we find
            // the trace point associated to each vertex. The element then
            // retains a pointer to the trace space points, to ensure
            // uniqueness of normals when retrieving from two adjoining
            // elements which do not lie in a plane.
            
            int ElmtPointGeom = 0;
            int TracePointGeom = 0;
            for (int i = 0; i < m_exp->size(); ++i)
            {
                for (int j = 0; j < (*m_exp)[i]->GetNverts(); ++j)
                {
                    ElmtPointGeom  = ((*m_exp)[i]->GetGeom1D())->GetVid(j);
                    
                    for (int k = 0; k < m_trace->GetExpSize(); ++k)
                    {
                        TracePointGeom = m_trace->GetPhys_Offset(k);
                        
                        if (TracePointGeom == ElmtPointGeom)
                        {
                            LocalRegions::Expansion1DSharedPtr exp1d
                                = boost::dynamic_pointer_cast<LocalRegions::Expansion1D>((*m_exp)[i]);
                            LocalRegions::Expansion0DSharedPtr exp0d
                                = boost::dynamic_pointer_cast<LocalRegions::Expansion0D>(m_trace->GetExp(k));
                            
                            exp0d->SetAdjacentElementExp(j,exp1d);
                            break;
                        }
                    }
                }
            }
            
            SetUpPhysNormals();
        }
		
        
        /**
         * New constructor for multidomain computations in arterial network
         * structures.
         * 
         * @param	domain	Subdomain specified in the inputfile from
         *       	      	which the DisContField1D is set up
         */
        DisContField1D::DisContField1D(const LibUtilities::SessionReaderSharedPtr &pSession,
                                       const SpatialDomains::CompositeMap &domain,
                                       const SpatialDomains::MeshGraphSharedPtr &graph1D,
                                       const std::string &variable,
                                       int i):
            ExpList1D(pSession,domain,graph1D,i),
            m_bndCondExpansions(),
            m_bndConditions()
        {			
            //1. a) Read in all boundary conditions specified in inputfile
            SpatialDomains::BoundaryConditions bcs(pSession, graph1D);
            
            //1. b) Generate Boundary Condition Expansions from the read in ones only if the 
            //      BC is in the currently processed subdomain
            GenerateMultiDomainBoundaryConditionExpansion(graph1D,bcs,variable,i);
            
            //1. c) Evaluate the boundary conditions
            EvaluateBoundaryConditions();
            //ApplyGeomInfo();
            
            //2. Set up trace information
            map<int,int> periodicVertices;
            GetPeriodicVertices(graph1D,bcs,variable,periodicVertices);
            
            m_globalBndMat = MemoryManager<GlobalLinSysMap>::AllocateSharedPtr();
            
            ExpList0DSharedPtr traces = MemoryManager<ExpList0D>::AllocateSharedPtr(m_bndCondExpansions,m_bndConditions,*m_exp,graph1D,periodicVertices);
            m_traces = Array<OneD,MultiRegions::ExpListSharedPtr> (domain.size());
            m_traces[i] = boost::dynamic_pointer_cast<ExpList>(traces);

            m_traceMap = MemoryManager<AssemblyMapDG>::
                AllocateSharedPtr(pSession,graph1D,traces,*this,
                                  m_bndCondExpansions,m_bndConditions,periodicVertices);
            
            tmpBndSol = Array<OneD,NekDouble>(m_traceMap->GetNumLocalBndCoeffs());
            
            // Scatter trace points to 1D elements. For each element, we find
            // the trace point associated to each vertex. The element then
            // retains a pointer to the trace space points, to ensure
            // uniqueness of normals when retrieving from two adjoining
            // elements which do not lie in a plane.
            int ElmtPointGeom = 0;
            int TracePointGeom = 0;
            for (int l = 0; l < m_exp->size(); ++l)
            {
                for (int j = 0; j < (*m_exp)[l]->GetNverts(); ++j)
                {
                    ElmtPointGeom  = ((*m_exp)[l]->GetGeom1D())->GetVid(j);
                    
                    for (int k = 0; k < m_traces[i]->GetExpSize(); ++k)
                    {
                        TracePointGeom = m_traces[i]->GetPhys_Offset(k);
			
                        if (TracePointGeom == ElmtPointGeom)
                        {
                            LocalRegions::Expansion1DSharedPtr exp1d
                                = boost::dynamic_pointer_cast<LocalRegions::Expansion1D>((*m_exp)[l]);
                            
                            LocalRegions::Expansion0DSharedPtr exp0d
                                = boost::dynamic_pointer_cast<LocalRegions::Expansion0D>(m_traces[i]->GetExp(k));
                            
                            exp0d->SetAdjacentElementExp(j,exp1d);
                            break;
                        }
                    }
                }
            }
            
            SetUpPhysNormals();
        }
        
        /**
         * Constructs a field as a copy of an existing field.
         * @param   In          Existing DisContField1D object to copy.
         */
        DisContField1D::DisContField1D(const DisContField1D &In):
            ExpList1D(In),
            m_bndCondExpansions(In.m_bndCondExpansions),
            m_bndConditions(In.m_bndConditions),
            m_globalBndMat(In.m_globalBndMat),
            m_trace(In.m_trace),
            m_traceMap(In.m_traceMap)
        {
        }
        
        
        /**
         * Constructs a field as a copy of an existing explist1D field.
         * @param   In          Existing ExpList1D object to copy.
         */
        DisContField1D::DisContField1D(const ExpList1D &In):
            ExpList1D(In)
	{
	}

        /**
         *
         */
        DisContField1D::~DisContField1D()
        {
        }


        /**
         * Generate the boundary condition expansion list
         * @param   graph1D     A mesh, containing information about the domain
         *                      and the spectral/hp element expansions.
         * @param   bcs         Information about the enforced boundary
         *                      conditions.
         * @param   variable    The session variable associated with the
         *                      boundary conditions to enforce.
         */
        void DisContField1D::GenerateBoundaryConditionExpansion(
            const SpatialDomains::MeshGraphSharedPtr &graph1D,
            SpatialDomains::BoundaryConditions &bcs,
            const std::string variable)
        {
            int i;
            int cnt  = 0;

            const SpatialDomains::BoundaryRegionCollection &bregions
                                                = bcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions
                                                = bcs.GetBoundaryConditions();

            int nbnd = bregions.size();
            // count the number of non-periodic boundary points
            for(i = 0; i < nbnd; ++i)
            {
                const SpatialDomains::BoundaryConditionShPtr boundaryCondition =
                        GetBoundaryCondition(bconditions, i, variable);
                if( boundaryCondition->GetBoundaryConditionType() != SpatialDomains::ePeriodic )
                {
                    SpatialDomains::BoundaryRegion::iterator bregionIt;
                    for (bregionIt = bregions[i]->begin(); bregionIt != bregions[i]->end(); bregionIt++)
                    {
                        cnt += bregionIt->second->size();
                    }
                }
            }			

            m_bndCondExpansions
                    = Array<OneD,MultiRegions::ExpListSharedPtr>(cnt);

            m_bndConditions
                    = Array<OneD,SpatialDomains::BoundaryConditionShPtr>(cnt);

            SetBoundaryConditionExpansion(graph1D,bcs,variable,
                                           m_bndCondExpansions,
                                           m_bndConditions);
        }
		
		
        /**
         * Generate the boundary condition expansion list in case of
         * mutidomain solver.
         * @param   graph1D     A mesh, containing information about the domain
         *                      and the spectral/hp element expansions.
         * @param   bcs         Information about the enforced boundary
         *                      conditions.
         * @param   variable    The session variable associated with the
         *                      boundary conditions to enforce
         * @param   subdomain   Currently processed subdomain.
         */
        void DisContField1D::GenerateMultiDomainBoundaryConditionExpansion(
            const SpatialDomains::MeshGraphSharedPtr &graph1D,
            SpatialDomains::BoundaryConditions &bcs,
            const std::string variable,
            int subdomain)
        {
            int i;
            int cnt  = 0;
			
            const SpatialDomains::BoundaryRegionCollection &bregions
			= bcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions
			= bcs.GetBoundaryConditions();
			
            int nbnd = bregions.size();
            // count the number of non-periodic boundary points
            for(i = 0; i < nbnd; ++i)
            {
                const SpatialDomains::BoundaryConditionShPtr boundaryCondition = GetBoundaryCondition(bconditions, i, variable);
                if( boundaryCondition->GetBoundaryConditionType() != SpatialDomains::ePeriodic )
                {
                    SpatialDomains::BoundaryRegion::iterator bregionIt;
                    for (bregionIt = bregions[i]->begin(); bregionIt != bregions[i]->end(); bregionIt++)
                    {
                        cnt += bregionIt->second->size();
                    }
                }
            }			
			
            m_bndCondExpansions = Array<OneD,MultiRegions::ExpListSharedPtr>(2);
			
            m_bndConditions = Array<OneD,SpatialDomains::BoundaryConditionShPtr>(2);
            
            SetMultiDomainBoundaryConditionExpansion(graph1D,bcs,variable,
                                                     m_bndCondExpansions,
                                                     m_bndConditions,
                                                     subdomain);
        }


        /**
         * @param   graph1D     A mesh containing information about the domain
         *                      and the spectral/hp element expansion.
         * @param   bcs         Information about the boundary conditions.
         * @param   variable    Specifies the field.
         * @param   periodicVertices    Map into which the list of periodic
         *                      vertices is placed.
         */
        void DisContField1D::GetPeriodicVertices(
            const SpatialDomains::MeshGraphSharedPtr &graph1D,
            const SpatialDomains::BoundaryConditions &bcs,
            const std::string variable,
            map<int,int>& periodicVertices)
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

            SpatialDomains::VertexComponentSharedPtr vert1;
            SpatialDomains::VertexComponentSharedPtr vert2;

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
                    region2ID = (boost::static_pointer_cast<SpatialDomains::
                                    PeriodicBoundaryCondition>(locBCond))
                                        ->m_connectedBoundaryRegion;

                    if(doneBndRegions.count(region1ID)==0)
                    {
                        ASSERTL0(bregions[region1ID]->size()
                                    == bregions[region2ID]->size(),
                                 "Size of the 2 periodic boundary regions "
                                 "should be equal");

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
                                if(!(vert1 = boost::dynamic_pointer_cast
                                        <SpatialDomains::VertexComponent>(
                                            (*comp1)[k]))||
                                   !(vert2 = boost::dynamic_pointer_cast
                                        <SpatialDomains::VertexComponent>(
                                            (*comp2)[k])))
                                {
                                    ASSERTL0(false,"dynamic cast to a "
                                                   "VertexComponent failed");
                                }

                                // Extract the periodic vertices
                                periodicVertices[vert1->GetVid()]
                                    = vert2->GetVid();
                            }
                        }
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

        /**
         * @param   graph1D     A mesh containing information about the domain
         *                      and the Spectral/hp element expansion.
         * @param   bcs         Information about the boundary conditions.
         * @param   variable    Specifies the field.
         * @param   bndCondExpansions   Array of ExpList1D objects each
         *                      containing a 1D spectral/hp element expansion
         *                      on a single boundary region.
         * @param   bncConditions   Array of BoundaryCondition objects which
         *                      contain information about the boundary
         *                      conditions on the different boundary regions.
         */
        void DisContField1D::SetBoundaryConditionExpansion(
            const SpatialDomains::MeshGraphSharedPtr &graph1D,
            const SpatialDomains::BoundaryConditions &bcs,
            const std::string variable,
            Array<OneD, MultiRegions::ExpListSharedPtr>
                &bndCondExpansions,
            Array<OneD, SpatialDomains
                ::BoundaryConditionShPtr> &bndConditions)
        {
            int i,k;
            int cnt  = 0;

            const SpatialDomains::BoundaryRegionCollection &bregions
                                                = bcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions
                                                = bcs.GetBoundaryConditions();

            MultiRegions::ExpList0DSharedPtr         locPointExp;
            SpatialDomains::BoundaryConditionShPtr   locBCond;
            SpatialDomains::VertexComponentSharedPtr vert;

            int nbnd = bregions.size();

            cnt=0;
            // list Dirichlet boundaries first
            for(i = 0; i < nbnd; ++i)
            {
                locBCond = GetBoundaryCondition(bconditions, i, variable);
                if(locBCond->GetBoundaryConditionType() == SpatialDomains::eDirichlet 
				   || locBCond->GetBoundaryConditionType()== SpatialDomains::eJunction
				   || locBCond->GetBoundaryConditionType()== SpatialDomains::eBifurcation
				   || locBCond->GetBoundaryConditionType()== SpatialDomains::eMerging)
                {
                    SpatialDomains::BoundaryRegion::iterator bregionIt;
                    for (bregionIt = bregions[i]->begin(); bregionIt != bregions[i]->end(); bregionIt++)
                    {
                        for(k = 0; k < bregionIt->second->size(); k++)
                        {
                            if((vert = boost::dynamic_pointer_cast
                                    <SpatialDomains::VertexComponent>(
                                         (*bregionIt->second)[k])))
                            {
                                locPointExp
                                    = MemoryManager<MultiRegions::ExpList0D>
                                                ::AllocateSharedPtr(vert);
                                bndCondExpansions[cnt]  = locPointExp;
                                bndConditions[cnt++]    = locBCond;
                            }
                            else
                            {
                                ASSERTL0(false,
                                         "dynamic cast to a vertex failed");
                            }
                        }
                    }
                }
            } // end if Dirichlet

            // then, list the other (non-periodic) boundaries
            for(i = 0; i < nbnd; ++i)
            {
                locBCond = GetBoundaryCondition(bconditions, i, variable);

                switch(locBCond->GetBoundaryConditionType())
                {
                case SpatialDomains::eNeumann:
                case SpatialDomains::eRobin:
                    {
                        SpatialDomains::BoundaryRegion::iterator bregionIt;
                        for (bregionIt = bregions[i]->begin(); bregionIt != bregions[i]->end(); bregionIt++)
                        {
                            for(k = 0; k < bregionIt->second->size(); k++)
                            {
                                if((vert = boost::dynamic_pointer_cast
                                        <SpatialDomains::VertexComponent>(
                                            (*bregionIt->second)[k])))
                                {
                                    locPointExp
                                        = MemoryManager<MultiRegions::ExpList0D>
                                        ::AllocateSharedPtr(vert);
                                    bndCondExpansions[cnt]  = locPointExp;
                                    bndConditions[cnt++]    = locBCond;
                                }
                                else
                                {
                                    ASSERTL0(false,
                                             "dynamic cast to a vertex failed");
                                }
                            }
                        }
                    }
                case SpatialDomains::eDirichlet: // do nothing for these types
				case SpatialDomains::eJunction: 
				case SpatialDomains::eBifurcation:
				case SpatialDomains::eMerging:
				case SpatialDomains::ePeriodic:
                    break;
                default:
                    ASSERTL0(false,"This type of BC not implemented yet");
                    break;
                }
            }
        }

        /**
         * MultiDomain case
         * @param   graph1D     A mesh containing information about the domain
         *                      and the Spectral/hp element expansion.
         * @param   bcs         Information about the boundary conditions.
         * @param   variable    Specifies the field.
         * @param   bndCondExpansions   Array of ExpList1D objects each
         *                      containing a 1D spectral/hp element expansion
         *                      on a single boundary region.
         * @param   bncConditions   Array of BoundaryCondition objects which
         *                      contain information about the boundary
         *                      conditions on the different boundary regions.
         */
        void DisContField1D::SetMultiDomainBoundaryConditionExpansion(
            const SpatialDomains::MeshGraphSharedPtr &graph1D,
            const SpatialDomains::BoundaryConditions &bcs,
            const std::string variable,Array<OneD, MultiRegions::ExpListSharedPtr>
                &bndCondExpansions,
            Array<OneD, SpatialDomains::BoundaryConditionShPtr> &bndConditions,
            int subdomain)
        {
            int i,k;
            int cnt  = 0;
			
            const SpatialDomains::BoundaryRegionCollection &bregions
				= bcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions
				= bcs.GetBoundaryConditions();
			
            MultiRegions::ExpList0DSharedPtr         locPointExp;
            SpatialDomains::BoundaryConditionShPtr   locBCond;
            SpatialDomains::VertexComponentSharedPtr vert;
            
            // Find the first boundary condition in the current domain
            int firstcondition = subdomain*2;
            int secondcondition = subdomain*2+1;
            cnt = 0;
			
            // list Dirichlet boundaries first
            for(i = firstcondition; i < secondcondition+1; ++i)
            {
                locBCond = GetBoundaryCondition(bconditions, i, variable);
                if(locBCond->GetBoundaryConditionType() == SpatialDomains::eDirichlet 
				   || locBCond->GetBoundaryConditionType()== SpatialDomains::eJunction
				   || locBCond->GetBoundaryConditionType()== SpatialDomains::eBifurcation
				   || locBCond->GetBoundaryConditionType()== SpatialDomains::eMerging)
                {
                    SpatialDomains::BoundaryRegion::iterator bregionIt;
                    for (bregionIt = bregions[i]->begin(); bregionIt != bregions[i]->end(); bregionIt++)
                    {
                        for(k = 0; k < bregionIt->second->size(); k++)
                        {
                            if((vert = boost::dynamic_pointer_cast<SpatialDomains::VertexComponent>((*bregionIt->second)[k])))
                            {
                                locPointExp = MemoryManager<MultiRegions::ExpList0D>::AllocateSharedPtr(vert);
                                bndCondExpansions[cnt]  = locPointExp;
                                bndConditions[cnt++]    = locBCond;
                            }
                            else
                            {
                                ASSERTL0(false,"dynamic cast to a vertex failed");
                            }
                        }
                    }
                }
            } // end if Dirichlet
			
            // then, list the other (non-periodic) boundaries
            for(i = firstcondition; i < secondcondition+1; ++i)
            {
                locBCond = GetBoundaryCondition(bconditions, i, variable);
				
                switch(locBCond->GetBoundaryConditionType())
                {
					case SpatialDomains::eNeumann:
					case SpatialDomains::eRobin:
                    {
                        SpatialDomains::BoundaryRegion::iterator bregionIt;
                        for (bregionIt = bregions[i]->begin(); bregionIt != bregions[i]->end(); bregionIt++)
                        {
                            for(k = 0; k < bregionIt->second->size(); k++)
                            {
                                if((vert = boost::dynamic_pointer_cast<
                                        SpatialDomains::VertexComponent>(
                                            (*bregionIt->second)[k])))
                                {
                                    locPointExp = MemoryManager<MultiRegions::ExpList0D>::AllocateSharedPtr(vert);
                                    bndCondExpansions[cnt]  = locPointExp;
                                    bndConditions[cnt++]    = locBCond;
                                }
                                else
                                {
                                    ASSERTL0(false,"dynamic cast to a vertex failed");
                                }
                            }
                        }
                    }
                    case SpatialDomains::eDirichlet: // do nothing for these types
                    case SpatialDomains::eJunction: 
                    case SpatialDomains::eBifurcation:
                    case SpatialDomains::eMerging:
                    case SpatialDomains::ePeriodic:
                        break;
                    default:
                        ASSERTL0(false,"This type of BC not implemented yet");
                        break;
                }
            }
        }
		

        /**
         *
         */
        GlobalLinSysSharedPtr DisContField1D::GetGlobalBndLinSys(
                    const GlobalLinSysKey &mkey)
        {
            ASSERTL0(mkey.GetMatrixType() == StdRegions::eHybridDGHelmBndLam,
                     "Routine currently only tested for HybridDGHelmholtz");

            ASSERTL1(mkey.GetGlobalSysSolnType() != eDirectFullMatrix,
                     "Full matrix global systems are not supported for HDG "
                     "expansions");

            ASSERTL1(mkey.GetGlobalSysSolnType()
                                        ==m_traceMap->GetGlobalSysSolnType(),
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
         * Generate the forward or backward state for each trace point.
         * @param   Fwd     Forward state.
         * @param   Bwd     Backward state.
         */
        void DisContField1D::v_GetFwdBwdTracePhys(Array<OneD, NekDouble> &Fwd,
                                                  Array<OneD, NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhys(m_phys,Fwd,Bwd);
        }
        
        
        /**
         * Generate the forward or backward state for each trace point.
         * @param   field   field.
         * @param   Fwd     Forward state.
         * @param   Bwd     Backward state.
         */
        void DisContField1D::v_GetFwdBwdTracePhys(const Array<OneD, const NekDouble> &field,
                                                  Array<OneD,       NekDouble> &Fwd,
                                                  Array<OneD,       NekDouble> &Bwd)
        {
            // Counter variables
            int  n, p, i;
            
            // Number of elements
            int nElements = GetExpSize();
            
            // Number of solution points of each element
            int nLocalSolutionPts;
            
            // Initial index of each element
            int phys_offset;
            
            // Index of each interface between adjacent elements
            int interface_offset;
            
            // Index in case of subdomains
            int subdomain_offset;
            
            // Coordinate of each standard element (-1, 1)
            NekDouble vertex_coord;
            
            // Basis shared pointer
            LibUtilities::BasisSharedPtr Basis;

            // Set forward and backard state to zero
            Vmath::Zero(Fwd.num_elements(), Fwd, 1);
            Vmath::Zero(Bwd.num_elements(), Bwd, 1);
			
            // Loop on the elements
            for (n = 0; n < nElements; ++n)
            {
                // Set the offset of each element
                phys_offset = GetPhys_Offset(n);
                
                // Set the number of solution points of each element
                nLocalSolutionPts = (*m_exp)[n]->GetNumPoints(0);
                
                // Temporary vector for interpolation routine
                Array<OneD, NekDouble> tmp(nLocalSolutionPts, 0.0);
                
                // Partition the field vector in local vectors
                Vmath::Vcopy(nLocalSolutionPts, 
                             (&field[phys_offset]), 1,
                             (&tmp[0]), 1);
                
                // Basis definition on each element
                Basis = (*m_exp)[n]->GetBasis(0);

                for (p = 0; p < 2; ++p)
                {
                    // Coordinate vector for interpolation routine
                    Array<OneD, NekDouble> interface_coord(3, 0.0);
                    
                    vertex_coord       = 0.0;
                    interface_offset   = (*m_exp)[n]->GetGeom1D()->GetVid(p);
                    subdomain_offset   = (*m_exp)[0]->GetGeom1D()->GetVid(0);
                    interface_offset  -= subdomain_offset;
                    
                    for (i = 0; i < ((*m_exp)[n]->GetVertexNormal(p)).num_elements(); i++)
                    {
                        vertex_coord += ((*m_exp)[n]->GetVertexNormal(p))[i][0];
                    }
                    
                    // Set the x-coordinate of the standard interface point
                    interface_coord[0] = -1.0;

                    // Implementation for every points except Gauss points
                    if (Basis->GetPointsType() != LibUtilities::eGaussGaussLegendre)
                    {
                        if(vertex_coord >= 0.0)
                        {
                            Fwd[interface_offset] = field[phys_offset+nLocalSolutionPts-1];
                        }			 
                        if(vertex_coord < 0.0) 
                        {
                            Bwd[interface_offset] = field[phys_offset];
                        }
                    }
                    // Implementation for Gauss points 
                    // (it doesn't work for WeakDG)
                    else
                    {
                        StdRegions::StdSegExp StdSeg(Basis->GetBasisKey());
                     
                        if(vertex_coord >= 0.0)
                        {
                            Fwd[interface_offset] = 
                            //(*m_exp)[n]->PhysEvaluate(interface_coord, tmp);
                            StdSeg.PhysEvaluate(interface_coord, tmp);
                        }					 
                        if(vertex_coord < 0.0) 
                        {
                            Bwd[interface_offset] = 
                            //(*m_exp)[n]->PhysEvaluate(interface_coord, tmp);
                            StdSeg.PhysEvaluate(interface_coord, tmp);
                        }
                    }
                }
            }
            
            // Fill boundary conditions into the missing elements
            int id1         = 0;
            int id2         = 0;
            int firstVertex = (*m_exp)[0]->GetGeom1D()->GetVid(0);
            int lastVertex  = (*m_exp)[nElements-1]->GetGeom1D()->GetVid(1);
            Array<OneD, NekDouble>  processed(m_bndCondExpansions.num_elements()+1, -1.0);
            
            // Bug temporary fixed for Periodic boundary conditions
            if(SpatialDomains::ePeriodic)
            {
                int nTracePts      = Fwd.num_elements();
                Fwd[0]             = Fwd[nTracePts-1];
                Bwd[nTracePts-1]   = Bwd[0];
            }
            
            for(n = 0; n < m_bndCondExpansions.num_elements(); ++n)
            {		
                // Check if the current boundary condition belongs to the current subdomain
                if((m_bndCondExpansions[n]->GetVertex()->GetVid() >= firstVertex) 
                   && (m_bndCondExpansions[n]->GetVertex()->GetVid() <= lastVertex) 
                   && (processed[n] != m_bndCondExpansions[n]->GetVertex()->GetVid()))
                {
                    if((m_bndConditions[n]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                       || (m_bndConditions[n]->GetBoundaryConditionType() == SpatialDomains::eJunction)
                       || (m_bndConditions[n]->GetBoundaryConditionType() == SpatialDomains::eBifurcation)
                       || (m_bndConditions[n]->GetBoundaryConditionType() == SpatialDomains::eMerging))
                    {
                        if(m_bndCondExpansions[n]->GetVertex()->GetVid() == lastVertex)
                        {
                            id1            = 0; //GetCoeff_Offset(n)+1;
                            id2            = m_bndCondExpansions[n]->GetVertex()->GetVid() - subdomain_offset;
                            Bwd[id2]       = m_bndCondExpansions[n]->GetCoeff(id1);
                            processed[n+1] = m_bndCondExpansions[n]->GetVertex()->GetVid();
                        }
                        else
                        {
                            id1        = 0; //GetCoeff_Offset(n);
                            id2        = m_bndCondExpansions[n]->GetVertex()->GetVid() - subdomain_offset;
                            Fwd[id2]   = m_bndCondExpansions[n]->GetCoeff(id1);
                        }
                        
                        //Previous working version
                        /*if(m_traceMap->GetBndExpAdjacentOrient(n) == eAdjacentEdgeIsForwards)
						 {
                         id1 = 0; //GetCoeff_Offset(n)+1;
                         id2 = m_bndCondExpansions[n]->GetVertex()->GetVid();
                         Bwd[id2] = m_bndCondExpansions[n]->GetCoeff(id1);
						 }
						 else
						 {
                         id1 = 0; //GetCoeff_Offset(n);
                         id2 = m_bndCondExpansions[n]->GetVertex()->GetVid();
                         Fwd[id2] = m_bndCondExpansions[n]->GetCoeff(id1);
						 }*/
                    }
                    else
                    {
                        ASSERTL0(false, "Method not set up for non-Dirichlet conditions");
                    }
                }	 
            }			 
		}
        
        
	
        void DisContField1D::v_ExtractTracePhys(
            Array<OneD, NekDouble> &outarray)
        {
            ASSERTL1(m_physState == true,"local physical space is not true ");
            v_ExtractTracePhys(m_phys, outarray);
        }
        
        void DisContField1D::v_ExtractTracePhys(
            const Array<OneD, const NekDouble> &inarray, 
                  Array<OneD,       NekDouble> &outarray)
        {
            // Loop over elemente and collect forward expansion
            int nexp = GetExpSize();
            int n,p,offset,phys_offset;
            
            ASSERTL1(outarray.num_elements() >= m_trace->GetExpSize(),
                     "input array is of insufficient length");
            
            for(n  = 0; n < nexp; ++n)
            {
                phys_offset = GetPhys_Offset(n);
		
                for(p = 0; p < (*m_exp)[n]->GetNverts(); ++p)
                {
                    offset = m_trace->GetPhys_Offset(p+n);
                    outarray[offset] = inarray[phys_offset];
                }
            }
        }		 
	
        void DisContField1D::v_AddTraceIntegral(
            const Array<OneD, const NekDouble> &Fn, 
                  Array<OneD,       NekDouble> &outarray)
        {
            int p, n, offset, t_offset;
            NekDouble vertnorm = 0.0;
            
            for (n = 0; n < GetExpSize(); ++n)
            {
                offset = GetCoeff_Offset(n);
		
                for(p = 0; p < 2; ++p)
                {
                    vertnorm = 0.0;
                    for (int i = 0; i<((*m_exp)[n]->GetVertexNormal(p)).num_elements(); i++)
                    {
                        vertnorm += ((*m_exp)[n]->GetVertexNormal(p))[i][0];
                    }
                    
                    t_offset = GetTrace()->GetPhys_Offset(n+p);
                    
                    if (vertnorm >= 0.0) 
                    {
                        outarray[offset+(*m_exp)[n]->GetVertexMap(1)] += Fn[t_offset];
                    }
                    
                    if (vertnorm < 0.0) 
                    {
                        outarray[offset] -= Fn[t_offset];
                    }
                }
            }
        }
	
	
        void DisContField1D::v_HelmSolve(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const FlagList &flags,
            const StdRegions::ConstFactorMap &factors,
            const StdRegions::VarCoeffMap &varcoeff,
            const Array<OneD, const NekDouble> &dirForcing)
        {
            int i,n,cnt,nbndry;
            int nexp = GetExpSize();
            Array<OneD,NekDouble> f(m_ncoeffs);
            DNekVec F(m_ncoeffs,f,eWrapper);
            Array<OneD,NekDouble> e_f, e_l;

            //----------------------------------
            // Setup RHS Inner product
            //----------------------------------
            IProductWRTBase(inarray,f);
            Vmath::Neg(m_ncoeffs,f,1);

            //----------------------------------
            // Solve continuous Boundary System
            //----------------------------------
            int GloBndDofs = m_traceMap->GetNumGlobalBndCoeffs();
            int NumDirBCs  = m_traceMap->GetNumLocalDirBndCoeffs();
            int e_ncoeffs,id;

            GlobalMatrixKey HDGLamToUKey(StdRegions::eHybridDGLamToU,NullAssemblyMapSharedPtr,factors,varcoeff);

            const DNekScalBlkMatSharedPtr &HDGLamToU = GetBlockMatrix(HDGLamToUKey);

			// Retrieve global trace space storage, \Lambda, from trace expansion
            Array<OneD,NekDouble> BndSol = tmpBndSol; //m_trace->UpdateCoeffs();

			
            Array<OneD,NekDouble> BndRhs(GloBndDofs,0.0);
            // Zero trace space
            Vmath::Zero(GloBndDofs,BndSol,1);

            int     LocBndCoeffs = m_traceMap->GetNumLocalBndCoeffs();
            Array<OneD, NekDouble> loc_lambda(LocBndCoeffs);
            DNekVec LocLambda(LocBndCoeffs,loc_lambda,eWrapper);

            //----------------------------------
            // Evaluate Trace Forcing
            //----------------------------------
            // Determing <u_lam,f> terms using HDGLamToU matrix
            for(cnt = n = 0; n < nexp; ++n)
            {
                nbndry = (*m_exp)[n]->NumDGBndryCoeffs();

                e_ncoeffs = (*m_exp)[n]->GetNcoeffs();
                e_f       = f+m_coeff_offset[n];
                e_l       = loc_lambda + cnt;

                // use outarray as tmp space
                DNekVec     Floc    (nbndry, e_l, eWrapper);
                DNekVec     ElmtFce (e_ncoeffs, e_f, eWrapper);
                Floc = Transpose(*(HDGLamToU->GetBlock(n,n)))*ElmtFce;

                cnt += nbndry;
            }

            // Assemble into global operator
            m_traceMap->AssembleBnd(loc_lambda,BndRhs);

            cnt = 0;
            // Copy Dirichlet boundary conditions into trace space
            for(i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    id = m_traceMap->GetBndCondCoeffsToGlobalCoeffsMap(i);
                    BndSol[id] = m_bndCondExpansions[i]->GetCoeff(0);
                }
                else
                {
                    id = m_traceMap->GetBndCondCoeffsToGlobalCoeffsMap(i);
                    BndRhs[id] += m_bndCondExpansions[i]->GetCoeff(0);
                }
            }

            //----------------------------------
            // Solve trace problem
            //----------------------------------
            if(GloBndDofs - NumDirBCs > 0)
            {
                GlobalLinSysKey       key(StdRegions::eHybridDGHelmBndLam,
                                          m_traceMap,factors);
                GlobalLinSysSharedPtr LinSys = GetGlobalBndLinSys(key);
                LinSys->Solve(BndRhs,BndSol,m_traceMap);
            }

            //----------------------------------
            // Internal element solves
            //----------------------------------
            GlobalMatrixKey invHDGhelmkey(StdRegions::eInvHybridDGHelmholtz,
                                            NullAssemblyMapSharedPtr,
                                            factors);

            const DNekScalBlkMatSharedPtr& InvHDGHelm
                                                = GetBlockMatrix(invHDGhelmkey);
            DNekVec out(m_ncoeffs,outarray,eWrapper);
            Vmath::Zero(m_ncoeffs,outarray,1);

            // get local trace solution from BndSol
            m_traceMap->GlobalToLocalBnd(BndSol,loc_lambda);

            //  out =  u_f + u_lam = (*InvHDGHelm)*f + (LamtoU)*Lam
            out = (*InvHDGHelm)*F + (*HDGLamToU)*LocLambda;
        }

        /**
         * Based on the expression \f$g(x,t)\f$ for the boundary conditions,
         * this function evaluates the boundary conditions for all boundaries
         * at time-level \a t.
         * @param   time        The time at which the boundary conditions
         *                      should be evaluated.
         * @param   bndCondExpansions   List of boundary expansions.
         * @param   bndConditions   Information about the boundary conditions.
         */
        void DisContField1D::v_EvaluateBoundaryConditions(const NekDouble time,
                                                          const NekDouble x2_in,
                                                          const NekDouble x3_in)
        {
            int i;

            NekDouble x0;
            NekDouble x1;
            NekDouble x2;
			
            for(i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                m_bndCondExpansions[i]->GetCoords(x0,x1,x2);

                if(x2_in != NekConstants::kNekUnsetDouble && x3_in != NekConstants::kNekUnsetDouble)
                {
                    x1 = x2_in;
                    x2 = x3_in;
                }

                if(m_bndConditions[i]->GetBoundaryConditionType()
                        == SpatialDomains::eDirichlet)
                {
                    m_bndCondExpansions[i]->SetCoeff(
                            (boost::static_pointer_cast<SpatialDomains
                             ::DirichletBoundaryCondition>(m_bndConditions[i])
                             ->m_dirichletCondition).Evaluate(x0,x1,x2,time));
                }
				else if((m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eJunction)||
						(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eBifurcation)||
						(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eMerging))
                {
                   //Do not update this conditions as they will be by the domain-linking Riemann solvers
                }	
				else if(m_bndConditions[i]->GetBoundaryConditionType()
                        == SpatialDomains::eNeumann)
                {
                    m_bndCondExpansions[i]->SetCoeff(
                            (boost::static_pointer_cast<SpatialDomains
                             ::NeumannBoundaryCondition>(m_bndConditions[i])
                             ->m_neumannCondition).Evaluate(x0,x1,x2,time));
                }
                else if(m_bndConditions[i]->GetBoundaryConditionType()
                        == SpatialDomains::eRobin)
                {
                    m_bndCondExpansions[i]->SetCoeff(
                            (boost::static_pointer_cast<SpatialDomains
                             ::RobinBoundaryCondition>(m_bndConditions[i])
                             ->m_robinFunction).Evaluate(x0,x1,x2,time));

                    m_bndCondExpansions[i]->SetPhys(
                            (boost::static_pointer_cast<SpatialDomains
                             ::RobinBoundaryCondition>(m_bndConditions[i])
                             ->m_robinPrimitiveCoeff).Evaluate(x0,x1,x2,time));

                }
                else
                {
                    ASSERTL0(false,"This type of BC not implemented yet");
                }
            }
        }

        // Set up a list of element ids and edge ids that link to the
        // boundary conditions
        void DisContField1D::v_GetBoundaryToElmtMap(
            Array<OneD, int> &ElmtID, Array<OneD,int> &VertID)
        {
            map<int, int> VertGID;
            int i,n,id;
            int bid,cnt,Vid;
            int nbcs = m_bndConditions.num_elements();

            // make sure arrays are of sufficient length
            if(ElmtID.num_elements() != nbcs)
            {
                ElmtID = Array<OneD, int>(nbcs,-1);
            }
            else
            {
                fill(ElmtID.get(), ElmtID.get()+nbcs, -1);
            }

            if(VertID.num_elements() != nbcs)
            {
                VertID = Array<OneD, int>(nbcs);
            }

            // setup map of all global ids along boundary
            for(cnt = n = 0; n < m_bndCondExpansions.num_elements(); ++n)
            {
                Vid =  m_bndCondExpansions[n]->GetGeom()->GetVid();
                VertGID[Vid] = cnt++;
            }

            // Loop over elements and find verts that match;
            for(cnt = n = 0; n < GetExpSize(); ++n)
            {
                for(i = 0; i < (*m_exp)[n]->GetNverts(); ++i)
                {
                    id = (*m_exp)[n]->GetGeom()->GetVid(i);

                    if(VertGID.count(id) > 0)
                    {
                        bid = VertGID.find(id)->second;
                        ASSERTL1(ElmtID[bid] == -1,"Edge already set");
                        ElmtID[bid] = n;
                        VertID[bid] = i;
                        cnt ++;
                    }
                }
            }

            ASSERTL1(cnt == nbcs,"Failed to visit all boundary condtiions");
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
        map<int, RobinBCInfoSharedPtr> DisContField1D::v_GetRobinBCInfo(void)
        {
            int i;
            map<int, RobinBCInfoSharedPtr> returnval;
            Array<OneD, int> ElmtID,VertID;
            GetBoundaryToElmtMap(ElmtID,VertID);

            for(i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                MultiRegions::ExpListSharedPtr locExpList;

                if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eRobin)
                {
                    int elmtid;
                    Array<OneD, NekDouble> Array_tmp;

                    locExpList = m_bndCondExpansions[i];

                    RobinBCInfoSharedPtr rInfo = MemoryManager<RobinBCInfo>::AllocateSharedPtr(VertID[i],Array_tmp = locExpList->GetPhys());

                    elmtid = ElmtID[i];
                    // make link list if necessary (not likely in
                    // 1D but needed in 2D & 3D)
                    if(returnval.count(elmtid) != 0)
                    {
                        rInfo->next = returnval.find(elmtid)->second;
                    }
                    returnval[elmtid] = rInfo;
                }
            }

            return returnval;
        }
    } // end of namespace
} //end of namespace
