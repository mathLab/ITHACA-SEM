//////////////////////////////////////////////////////////////////////////////
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
#include <LibUtilities/Foundations/ManagerAccess.h>

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
        DisContField1D::DisContField1D(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const SpatialDomains::MeshGraphSharedPtr &graph1D,
            const std::string &variable)
            : ExpList1D(pSession,graph1D),
              m_bndCondExpansions(),
              m_bndConditions()
        {
            SpatialDomains::BoundaryConditions bcs(m_session, graph1D);

            SetupBCsTrace(pSession,graph1D,bcs,variable);
        }

        void  DisContField1D::SetupBCsTrace(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const SpatialDomains::MeshGraphSharedPtr   &graph1D,
            const SpatialDomains::BoundaryConditions   &bcs,
            const std::string &variable)
        {
            
            GenerateBoundaryConditionExpansion(graph1D,bcs,variable);
            EvaluateBoundaryConditions();
            ApplyGeomInfo();

            map<int,int> periodicVertices;
            GetPeriodicVertices(graph1D,bcs,variable,periodicVertices);

            m_globalBndMat = MemoryManager<GlobalLinSysMap>::AllocateSharedPtr();

            ExpList0DSharedPtr trace = MemoryManager<ExpList0D>::
                AllocateSharedPtr(
                    m_bndCondExpansions,
                    m_bndConditions,
                    *m_exp,graph1D,
                    periodicVertices);
            
            m_trace = boost::dynamic_pointer_cast<ExpList>(trace);

            m_traceMap = MemoryManager<AssemblyMapDG>::
                AllocateSharedPtr(
                    pSession,graph1D,trace,*this,
                    m_bndCondExpansions,m_bndConditions,periodicVertices);
            
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
                    
                    ElmtPointGeom  = (*m_exp)[i]->GetGeom1D()->GetVid(j);
                    
                    for (int k = 0; k < m_trace->GetExpSize(); ++k)
                    {
                        TracePointGeom = m_trace->GetExp(k)->GetGeom()->GetVid(0);
                        
                        if (TracePointGeom == ElmtPointGeom)
                        {
                            LocalRegions::Expansion1DSharedPtr exp1d
                                = boost::dynamic_pointer_cast
                                    <LocalRegions::Expansion1D>((*m_exp)[i]);
                            LocalRegions::Expansion0DSharedPtr exp0d
                                = boost::dynamic_pointer_cast
                                    <LocalRegions::Expansion0D>
                                        (m_trace->GetExp(k));
                            
                            exp0d->SetAdjacentElementExp(j,exp1d);
                            break;
                        }
                    }
                }
            }
            
            SetUpPhysNormals();

            int cnt, n, e;

            // Identify boundary verts
            for(cnt = 0, n = 0; n < m_bndCondExpansions.num_elements(); ++n)
            {
                if (m_bndConditions[n]->GetBoundaryConditionType() != 
                    SpatialDomains::ePeriodic)
                {
                    for(e = 0; e < m_bndCondExpansions[n]->GetExpSize(); ++e)
                    {
                        m_boundaryVerts.insert(
                            m_traceMap->GetBndCondTraceToGlobalTraceMap(cnt+e));
                    }
                }
                else
                {
                    ASSERTL0(false,"Periodic verts need setting up");
                }
                cnt += m_bndCondExpansions[n]->GetExpSize();
            }

            // Set up left-adjacent edge list.
            m_leftAdjacentVerts.resize(2*((*m_exp).size()));
            
            cnt = 0;
            // count size of trace
            for (cnt = n = 0; n < m_exp->size(); ++n)
            {
                for (e = 0; e < (*m_exp)[n]->GetNverts(); ++e, ++cnt)
                {
                    m_leftAdjacentVerts[cnt] = IsLeftAdjacentVertex(n, e);
                }
            }
        }


        bool DisContField1D::IsLeftAdjacentVertex(const int n, const int e)
        {
            set<int>::iterator it;
            LocalRegions::Expansion0DSharedPtr traceEl = 
                boost::dynamic_pointer_cast<LocalRegions::Expansion0D>(
                    (m_traceMap->GetElmtToTrace())[n][e]);
            
            bool fwd = true;
            if (traceEl->GetLeftAdjacentElementVertex () == -1 ||
                traceEl->GetRightAdjacentElementVertex() == -1)
            {
                // Boundary edge (1 connected element). Do nothing in
                // serial.
                it = m_boundaryVerts.find(traceEl->GetElmtId());
                
                // If the edge does not have a boundary condition set on
                // it, then assume it is a partition edge.
                if (it == m_boundaryVerts.end())
                {
                    ASSERTL0(false,"Parallel method not set up");
#if 0 
                    int traceGeomId = traceEl->GetGeom1D()->GetGlobalID();
                    PeriodicMap::iterator pIt = m_periodicVerts.find(
                        traceGeomId);

                    if (pIt != m_periodicVerts.end() && !pIt->second[0].isLocal)
                    {
                        fwd = traceGeomId == min(traceGeomId,pIt->second[0].id);
                    }
                    else
                    {
                        fwd = m_traceMap->
                            GetTraceToUniversalMapUnique(offset) >= 0;
                    }
#endif
                }
            }
            else if (traceEl->GetLeftAdjacentElementVertex () != -1 &&
                     traceEl->GetRightAdjacentElementVertex() != -1)
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
		
        
        // Given all boundary regions for the whole solution determine
        // which ones (if any) are part of domain and ensure all other
        // conditions are given as UserDefined Dirichlet. 
        SpatialDomains::BoundaryConditionsSharedPtr DisContField1D::GetDomainBCs(
            const SpatialDomains::CompositeMap &domain,
            const SpatialDomains::BoundaryConditions &Allbcs,
            const std::string &variable)
        {            
            SpatialDomains::BoundaryConditionsSharedPtr returnval;
            
            returnval = MemoryManager<SpatialDomains::BoundaryConditions>::AllocateSharedPtr();
            
            SpatialDomains::BoundaryRegion::iterator bregionIt;
            map<int,int> GeometryToRegionsMap;

            SpatialDomains::BoundaryRegionCollection::const_iterator it;

            const SpatialDomains::BoundaryRegionCollection &bregions
                = Allbcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions
                = Allbcs.GetBoundaryConditions();

            // Set up a map of all boundary regions
            for(it = bregions.begin(); it != bregions.end(); ++it)
            {
                SpatialDomains::BoundaryRegion::iterator bregionIt;
                for (bregionIt  = it->second->begin();
                     bregionIt != it->second->end(); bregionIt++)
                {
                    // can assume that all regions only contain one point in 1D
                    // Really do not need lopp above
                    int id = (*(bregionIt->second))[0]->GetGlobalID();
                    GeometryToRegionsMap[id] = it->first;
                }
            }

            SpatialDomains::CompositeMapConstIter domIt; 
            map<int,SpatialDomains::GeometrySharedPtr> EndOfDomain;

            // Now find out which points in domain have only one vertex
            for(domIt = domain.begin(); domIt != domain.end(); ++domIt)
            {
                SpatialDomains::Composite geomvector = domIt->second; 
                for(int i = 0; i < geomvector->size(); ++i)
                {
                    for(int j = 0; j < 2; ++j)
                    {
                        int vid = (*geomvector)[i]->GetVid(j); 
                        if(EndOfDomain.count(vid) == 0)
                        {
                            EndOfDomain[vid] = (*geomvector)[i]->GetVertex(j);
                        }
                        else
                        {
                            EndOfDomain.erase(vid);
                        }
                    }
                }
            }
            ASSERTL1(EndOfDomain.size() == 2,"Did not find two ends of domain");
             
            map<int,SpatialDomains::GeometrySharedPtr>::iterator regIt;
            for(regIt = EndOfDomain.begin(); regIt != EndOfDomain.end(); ++regIt)
            {
                if(GeometryToRegionsMap.count(regIt->first) != 0) // Set up boundary condition up
                {
                    map<int,int>::iterator iter = GeometryToRegionsMap.find(regIt->first);
                    ASSERTL1(iter != GeometryToRegionsMap.end(),"Failied to find GeometryToRegionMap");
                    int regionId = iter->second;
                    SpatialDomains::BoundaryRegionCollection::const_iterator bregionsIter = bregions.find(regionId);
                    ASSERTL1(bregionsIter != bregions.end(),"Failed to find boundary region");
                    SpatialDomains::BoundaryRegionShPtr breg = bregionsIter->second;
                    returnval->AddBoundaryRegions   (regionId,breg);

                    SpatialDomains::BoundaryConditionCollection::const_iterator bconditionsIter = bconditions.find(regionId);
                    ASSERTL1(bconditionsIter != bconditions.end(),"Failed to find boundary collection");
                    SpatialDomains::BoundaryConditionMapShPtr bcond = bconditionsIter->second;
                    returnval->AddBoundaryConditions(regionId,bcond);
                }
                else // Set up an undefined region. 
                {
                    SpatialDomains::BoundaryRegionShPtr breg(MemoryManager<SpatialDomains::BoundaryRegion>::AllocateSharedPtr());
                    
                    // Set up Composite (GemetryVector) to contain vertex and put into bRegion 
                    SpatialDomains::Composite gvec(MemoryManager<SpatialDomains::GeometryVector>::AllocateSharedPtr());
                    gvec->push_back(regIt->second);
                    (*breg)[regIt->first] = gvec;

                    returnval->AddBoundaryRegions(bregions.size()+1,breg);

                    SpatialDomains::BoundaryConditionMapShPtr bCondition = MemoryManager<SpatialDomains::BoundaryConditionMap>::AllocateSharedPtr();

                    // Set up just boundary condition for this variable. 
                    SpatialDomains::BoundaryConditionShPtr notDefinedCondition(MemoryManager<SpatialDomains::NotDefinedBoundaryCondition>::AllocateSharedPtr(m_session, "0"));
                    (*bCondition)[variable] = notDefinedCondition;
                    
                    returnval->AddBoundaryConditions(bregions.size()+1,bCondition);

                }
            }
            
            return returnval; 
        } 
        
        /**
         * Constructor for use in multidomain computations where a
         * domain list can be passed instead of graph1D
         * 
         * @param	domain	Subdomain specified in the inputfile from
         *       	      	which the DisContField1D is set up
         */
        DisContField1D::DisContField1D(
                                       const LibUtilities::SessionReaderSharedPtr &pSession,
                                       const SpatialDomains::MeshGraphSharedPtr &graph1D,
                                       const SpatialDomains::CompositeMap &domain,
                                       const SpatialDomains::BoundaryConditions &Allbcs, 
                                       const std::string &variable,
                                       bool SetToOneSpaceDimension):
            ExpList1D(pSession,graph1D,domain, true,variable,SetToOneSpaceDimension),
            m_bndCondExpansions(),
            m_bndConditions()
        {			
            SpatialDomains::BoundaryConditionsSharedPtr DomBCs = GetDomainBCs(domain,Allbcs,variable); 
            
            SetupBCsTrace(pSession,graph1D,*DomBCs,variable);
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
            m_traceMap(In.m_traceMap),
            m_boundaryVerts(In.m_boundaryVerts),
            m_leftAdjacentVerts(In.m_leftAdjacentVerts)
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
            const SpatialDomains::BoundaryConditions &bcs,
            const std::string variable)
        {
            int cnt  = 0;

            const SpatialDomains::BoundaryRegionCollection &bregions
                                                = bcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions
                                                = bcs.GetBoundaryConditions();
            SpatialDomains::BoundaryRegionCollection::const_iterator it;

            // count the number of non-periodic boundary points
            for (it = bregions.begin(); it != bregions.end(); ++it)
            {
                const SpatialDomains::BoundaryConditionShPtr boundaryCondition =
                    GetBoundaryCondition(bconditions, it->first, variable);
                if (boundaryCondition->GetBoundaryConditionType() !=
                    SpatialDomains::ePeriodic )
                {
                    SpatialDomains::BoundaryRegion::iterator bregionIt;
                    for (bregionIt  = it->second->begin();
                         bregionIt != it->second->end(); bregionIt++)
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
            /*
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
            */
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
            int k;
            int cnt  = 0;

            const SpatialDomains::BoundaryRegionCollection &bregions
                = bcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions
                = bcs.GetBoundaryConditions();
            SpatialDomains::BoundaryRegionCollection::const_iterator it;
            
            MultiRegions::ExpList0DSharedPtr         locPointExp;
            SpatialDomains::BoundaryConditionShPtr   locBCond;
            SpatialDomains::VertexComponentSharedPtr vert;

            cnt = 0;
            // list Dirichlet boundaries first
            for (it = bregions.begin(); it != bregions.end(); ++it)
            {
                locBCond = GetBoundaryCondition(
                    bconditions, it->first, variable);

                if (locBCond->GetBoundaryConditionType() ==
                    SpatialDomains::eDirichlet)

                {
                    SpatialDomains::BoundaryRegion::iterator bregionIt;
                    for (bregionIt  = it->second->begin();
                         bregionIt != it->second->end(); bregionIt++)
                    {
                        for (k = 0; k < bregionIt->second->size(); k++)
                        {
                            if ((vert = boost::dynamic_pointer_cast
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
            for (it = bregions.begin(); it != bregions.end(); ++it)
            {
                locBCond = GetBoundaryCondition(bconditions, it->first, variable);
                
                switch(locBCond->GetBoundaryConditionType())
                {
                case SpatialDomains::eNeumann:
                case SpatialDomains::eRobin:
                case SpatialDomains::eNotDefined: // presume this will be reused as Neuman, Robin or Dirichlet later
                    {
                        SpatialDomains::BoundaryRegion::iterator bregionIt;
                        for (bregionIt  = it->second->begin();
                             bregionIt != it->second->end(); bregionIt++)
                        {
                            for (k = 0; k < bregionIt->second->size(); k++)
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
                    // do nothing for these types
                case SpatialDomains::eDirichlet:
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

            if (matrixIter == m_globalBndMat->end())
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
         * @brief This method extracts the "forward" and "backward" trace data
         * from the array @a field and puts the data into output vectors @a Fwd
         * and @a Bwd.
         * 
         * We first define the convention which defines "forwards" and
         * "backwards". First an association is made between the vertex of each
         * element and its corresponding vertex in the trace space using the
         * mapping #m_traceMap. The element can either be left-adjacent or
         * right-adjacent to this trace edge (see
         * Expansion0D::GetLeftAdjacentElementExp). Boundary edges are never 
         * left-adjacent since elemental left-adjacency is populated first.
         * 
         * If the element is left-adjacent we extract the vertex trace data from
         * @a field into the forward trace space @a Fwd; otherwise, we place it
         * in the backwards trace space @a Bwd. In this way, we form a unique
         * set of trace normals since these are always extracted from
         * left-adjacent elements.
         *
         * @param field is a NekDouble array which contains the 1D data
         *              from which we wish to extract the backward and forward
         *              orientated trace/edge arrays.
         * @param Fwd   The resulting forwards space.
         * @param Bwd   The resulting backwards space.
         */

        void DisContField1D::v_GetFwdBwdTracePhys(
            const Array<OneD, const NekDouble> &field,
            Array<OneD,       NekDouble> &Fwd,
            Array<OneD,       NekDouble> &Bwd)
        {
            // Counter variables
            int  n, v;
            
            // Number of elements
            int nElements = GetExpSize(); 
            
            // Number of solution points of each element
            int nLocalSolutionPts;
            
            // Initial index of each element
            int phys_offset;
            
            Array<OneD, Array<OneD, StdRegions::StdExpansionSharedPtr> >
                &elmtToTrace = m_traceMap->GetElmtToTrace();
            
            // Basis shared pointer
            LibUtilities::BasisSharedPtr Basis;
            
            // Set forward and backard state to zero
            Vmath::Zero(Fwd.num_elements(), Fwd, 1);
            Vmath::Zero(Bwd.num_elements(), Bwd, 1);
			
            int cnt;
            // Loop on the elements
            for (cnt = n = 0; n < nElements; ++n)
            {
                // Set the offset of each element
                phys_offset = GetPhys_Offset(n);
                
                // Set the number of solution points of each element
                nLocalSolutionPts = (*m_exp)[n]->GetNumPoints(0);
                
                Basis = (*m_exp)[n]->GetBasis(0);

                for(v = 0; v < 2; ++v, ++cnt)
                {
                    int offset = m_trace->GetPhys_Offset(elmtToTrace[n][v]->GetElmtId());
                    
                    if (m_leftAdjacentVerts[cnt])
                    {
                        (*m_exp)[n]->GetVertexPhysVals(v, field + phys_offset,
                                                       Fwd[offset]);
                    }
                    else
                    {
                        (*m_exp)[n]->GetVertexPhysVals(v, field + phys_offset,
                                                       Bwd[offset]);
                    }
                }
            }
            
            // Fill boundary conditions into missing elements.
            int id = 0;
            
            for(cnt = n = 0; n < m_bndCondExpansions.num_elements(); ++n)
            {				
                if (m_bndConditions[n]->GetBoundaryConditionType() == 
                        SpatialDomains::eDirichlet)
                {
                    id  = m_trace->GetPhys_Offset(m_traceMap->GetBndCondTraceToGlobalTraceMap(cnt));
                    Bwd[id] = m_bndCondExpansions[n]->GetPhys()[0];
                    cnt++;
                }
                else if (m_bndConditions[n]->GetBoundaryConditionType() == 
                         SpatialDomains::eNeumann || 
                         m_bndConditions[n]->GetBoundaryConditionType() == 
                         SpatialDomains::eRobin)
                {
                    ASSERTL0((m_bndCondExpansions[n]->GetPhys())[0]==0.0,
                             "Method not set up for non-zero Neumann "
                             "boundary condition");
                    id  = m_trace->GetPhys_Offset(m_traceMap->GetBndCondTraceToGlobalTraceMap(cnt));
                    Bwd[id] = Fwd[id];
                    
                    cnt++;
                }
                else if (m_bndConditions[n]->GetBoundaryConditionType() ==
                         SpatialDomains::eNotDefined)
                {
                    // Do nothing
                }
                else if (m_bndConditions[n]->GetBoundaryConditionType() !=
                         SpatialDomains::ePeriodic)
                {
                    ASSERTL0(false,
                             "Method not set up for this boundary condition.");
                }
            }
            
#if 0 
            // Copy any periodic boundary conditions.
            for (n = 0; n < m_periodicFwdCopy.size(); ++n)
            {
                Bwd[m_periodicBwdCopy[n]] = Fwd[m_periodicFwdCopy[n]];
            }
#endif

            // Do parallel exchange for forwards/backwards spaces.
            m_traceMap->UniversalTraceAssemble(Fwd);
            m_traceMap->UniversalTraceAssemble(Bwd);

        }
        
	
        void DisContField1D::v_ExtractTracePhys(
            Array<OneD, NekDouble> &outarray)
        {
            ASSERTL1(m_physState == true,"local physical space is not true ");
            v_ExtractTracePhys(m_phys, outarray);
        }
        
        /**
         * @brief This method extracts the trace (verts in 1D) from the field @a
         * inarray and puts the values in @a outarray.
         *
         * It assumes the field is C0 continuous so that it can overwrite the
         * edge data when visited by the two adjacent elements.
         *
         * @param inarray   An array containing the 1D data from which we wish
         *                  to extract the edge data.
         * @param outarray  The resulting edge information.
         *
         * This will not work for non-boundary expansions
         */
        void DisContField1D::v_ExtractTracePhys(
            const Array<OneD, const NekDouble> &inarray, 
                  Array<OneD,       NekDouble> &outarray)
        {
            // Loop over elemente and collect forward expansion
            int nexp = GetExpSize();
            int n,p,offset,phys_offset;
            
            ASSERTL1(outarray.num_elements() >= m_trace->GetExpSize(),
                "input array is of insufficient length");
            
            for (n  = 0; n < nexp; ++n)
            {
                phys_offset = GetPhys_Offset(n);
		
                for (p = 0; p < (*m_exp)[n]->GetNverts(); ++p)
                {
                    offset = m_trace->GetPhys_Offset(
                                 (m_traceMap->GetElmtToTrace())[n][p]->GetElmtId());
                    (*m_exp)[n]->GetVertexPhysVals(p,  inarray + phys_offset,
                                                   outarray[offset]);
                }
            }
        }		 
	
        void DisContField1D::v_AddTraceIntegral(
            const Array<OneD, const NekDouble> &Fn, 
                  Array<OneD,       NekDouble> &outarray)
        {
            int p,n,offset, t_offset;
            double vertnorm =0.0;


            Array<OneD, Array<OneD, StdRegions::StdExpansionSharedPtr> >
                &elmtToTrace = m_traceMap->GetElmtToTrace();
            
            // Basis shared pointer
            LibUtilities::BasisSharedPtr Basis;
            
            for (n = 0; n < GetExpSize(); ++n)
            {
                // Basis definition on each element
                Basis = (*m_exp)[n]->GetBasis(0);
                
                // Number of coefficients on each element
                int e_ncoeffs = (*m_exp)[n]->GetNcoeffs();
                
                offset = GetCoeff_Offset(n);
                
                // Implementation for every points except Gauss points
                if (Basis->GetBasisType() != LibUtilities::eGauss_Lagrange)
                {
#if 0 // Currently assume p=1 normal is always positive and p=1 is
      // negative since this is the sign-convention of the local
      // segment. Note that we do not put sign in Fn in 1D since there is only one component
                    for(p = 0; p < 2; ++p)
                    {
                       
                        vertnorm = 0.0;
                        for (int i=0; i<((*m_exp)[n]->
                                         GetVertexNormal(p)).num_elements(); i++)
                        {
                            vertnorm += ((*m_exp)[n]->GetVertexNormal(p))[i][0];
                        }
                        
                        t_offset = GetTrace()->GetPhys_Offset(elmtToTrace[n][p]->GetElmtId());
                        
                        if (vertnorm >= 0.0)
                        {
                            outarray[offset+(*m_exp)[n]->GetVertexMap(1)] +=
                            Fn[t_offset];
                        }
                        
                        if (vertnorm < 0.0)
                        {
                            outarray[offset] -= Fn[t_offset];
                        }
                    }
#else
                    t_offset = GetTrace()->GetCoeff_Offset(elmtToTrace[n][0]->GetElmtId());
                    outarray[offset] -= Fn[t_offset];
                    
                    t_offset = GetTrace()->GetCoeff_Offset(elmtToTrace[n][1]->GetElmtId());
                    outarray[offset+(*m_exp)[n]->GetVertexMap(1)] += Fn[t_offset];
#endif

                }
                else
                {
                    DNekMatSharedPtr                     m_Ixm;
                    LibUtilities::BasisSharedPtr BASE;
                    const LibUtilities::PointsKey
                            BS_p(e_ncoeffs,LibUtilities::eGaussGaussLegendre);
                    const LibUtilities::BasisKey
                            BS_k(LibUtilities::eGauss_Lagrange,e_ncoeffs,BS_p);
                    
                    BASE  = LibUtilities::BasisManager()[BS_k];
                    
                    Array<OneD, NekDouble> coords(3, 0.0);
                    
                    int j;
                    
                    for(p = 0; p < 2; ++p)
                    {
                        vertnorm = 0.0;
                        for (int i=0; i<((*m_exp)[n]->
                             GetVertexNormal(p)).num_elements(); i++)
                        {
                            vertnorm += ((*m_exp)[n]->GetVertexNormal(p))[i][0];
                            coords[0] = vertnorm ;
                        }
                        
                        t_offset = GetTrace()->GetPhys_Offset(n+p);
                        
                        if (vertnorm >= 0.0)
                        {
                            m_Ixm = BASE->GetI(coords);
                            
                            
                            for (j = 0; j < e_ncoeffs; j++)
                            {
                                outarray[offset + j]  +=
                                    (m_Ixm->GetPtr())[j] * Fn[t_offset];
                            }
                        }
                        
                        if (vertnorm < 0.0)
                        {
                            m_Ixm = BASE->GetI(coords);
                            
                            for (j = 0; j < e_ncoeffs; j++)
                            {
                                outarray[offset + j] -=
                                    (m_Ixm->GetPtr())[j] * Fn[t_offset];
                            }
                        }
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

            GlobalMatrixKey HDGLamToUKey(
                StdRegions::eHybridDGLamToU,
                NullAssemblyMapSharedPtr,
                factors,
                varcoeff);

            const DNekScalBlkMatSharedPtr &HDGLamToU =
                GetBlockMatrix(HDGLamToUKey);

            // Retrieve global trace space storage, \Lambda, from trace expansion
            Array<OneD,NekDouble> BndSol =  Array<OneD,NekDouble>
                (m_traceMap->GetNumLocalBndCoeffs());
            
			
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
            for (cnt = n = 0; n < nexp; ++n)
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
            for (i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                if (m_bndConditions[i]->GetBoundaryConditionType() ==
                    SpatialDomains::eDirichlet)
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
            if (GloBndDofs - NumDirBCs > 0)
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

            const DNekScalBlkMatSharedPtr& InvHDGHelm =
                GetBlockMatrix(invHDGhelmkey);
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

            Array<OneD, NekDouble> x0(1);
            Array<OneD, NekDouble> x1(1);
            Array<OneD, NekDouble> x2(1);
			
            for (i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                if (time == 0.0 || m_bndConditions[i]->GetUserDefined() ==
                    SpatialDomains::eTimeDependent)
                {
                    m_bndCondExpansions[i]->GetCoords(x0,x1,x2);
                    
                    if (x2_in != NekConstants::kNekUnsetDouble && x3_in !=
                        NekConstants::kNekUnsetDouble)
                    {
                        x1[0] = x2_in;
                        x2[0] = x3_in;
                    }
                    
                    if (m_bndConditions[i]->GetBoundaryConditionType() ==
                        SpatialDomains::eDirichlet)
                    {
                        m_bndCondExpansions[i]->SetCoeff(0,
                            (boost::static_pointer_cast<SpatialDomains
                             ::DirichletBoundaryCondition>(m_bndConditions[i])
                             ->m_dirichletCondition).Evaluate(x0[0],x1[0],x2[0],time));
                        m_bndCondExpansions[i]->SetPhys(0,m_bndCondExpansions[i]->GetCoeff(0));
                    }
                    else if (m_bndConditions[i]->GetBoundaryConditionType()
                            == SpatialDomains::eNeumann)
                    {
                        m_bndCondExpansions[i]->SetCoeff(0,
                            (boost::static_pointer_cast<SpatialDomains
                             ::NeumannBoundaryCondition>(m_bndConditions[i])
                             ->m_neumannCondition).Evaluate(x0[0],x1[0],x2[0],time));
                    }
                    else if (m_bndConditions[i]->GetBoundaryConditionType()
                            == SpatialDomains::eRobin)
                    {
                        m_bndCondExpansions[i]->SetCoeff(0,
                            (boost::static_pointer_cast<SpatialDomains
                             ::RobinBoundaryCondition>(m_bndConditions[i])
                             ->m_robinFunction).Evaluate(x0[0],x1[0],x2[0],time));
                        
                        m_bndCondExpansions[i]->SetPhys(0,
                            (boost::static_pointer_cast<SpatialDomains
                             ::RobinBoundaryCondition>(m_bndConditions[i])
                             ->m_robinPrimitiveCoeff).Evaluate(x0[0],x1[0],x2[0],time));
                        
                    }
                    else if (m_bndConditions[i]->GetBoundaryConditionType()
                             == SpatialDomains::eNotDefined)
                    {
                    }
                    else
                    {
                        ASSERTL0(false,"This type of BC not implemented yet");
                    }
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
            if (ElmtID.num_elements() != nbcs)
            {
                ElmtID = Array<OneD, int>(nbcs,-1);
            }
            else
            {
                fill(ElmtID.get(), ElmtID.get()+nbcs, -1);
            }

            if (VertID.num_elements() != nbcs)
            {
                VertID = Array<OneD, int>(nbcs);
            }

            // setup map of all global ids along boundary
            for (cnt = n = 0; n < m_bndCondExpansions.num_elements(); ++n)
            {
                Vid =  m_bndCondExpansions[n]->GetExp(0)->GetGeom()->GetVid(0);
                VertGID[Vid] = cnt++;
            }

            // Loop over elements and find verts that match;
            for (cnt = n = 0; n < GetExpSize(); ++n)
            {
                for (i = 0; i < (*m_exp)[n]->GetNverts(); ++i)
                {
                    id = (*m_exp)[n]->GetGeom()->GetVid(i);

                    if (VertGID.count(id) > 0)
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

            for (i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                MultiRegions::ExpListSharedPtr locExpList;

                if (m_bndConditions[i]->GetBoundaryConditionType() ==
                   SpatialDomains::eRobin)
                {
                    int elmtid;
                    Array<OneD, NekDouble> Array_tmp;

                    locExpList = m_bndCondExpansions[i];

                    RobinBCInfoSharedPtr rInfo =
                        MemoryManager<RobinBCInfo>::
                            AllocateSharedPtr(
                                VertID[i],Array_tmp = locExpList->GetPhys());

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
