//////////////////////////////////////////////////////////////////////////////
//
// File DisContField.cpp
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

#include <MultiRegions/DisContField.h>
#include <StdRegions/StdSegExp.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <SpatialDomains/MeshGraph.h>
#include <LocalRegions/Expansion0D.h>
#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/HexExp.h>
#include <LocalRegions/TetExp.h>


using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class DisContField
         * This class augments the list of local expansions inherited from
         * ExpList with boundary conditions. Inter-element boundaries are
         * handled using an discontinuous Galerkin scheme.
         */

        /**
         * Constructs an empty expansion list with no boundary conditions.
         */
        DisContField::DisContField():
            ExpList(),
            m_bndCondExpansions(),
            m_bndConditions(),
            m_trace(NullExpListSharedPtr)
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
        DisContField::DisContField(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const SpatialDomains::MeshGraphSharedPtr   &graph,
            const std::string                          &variable,
            const bool                                  SetUpJustDG,
            const bool                                  DeclareCoeffPhysArrays,
            const Collections::ImplementationType       ImpType)
            : ExpList(pSession,graph,DeclareCoeffPhysArrays, variable, ImpType),
              m_trace(NullExpListSharedPtr)
        {
            if (variable.compare("DefaultVar") != 0)
            {
                SpatialDomains::BoundaryConditions bcs(m_session, graph);

                GenerateBoundaryConditionExpansion(graph,bcs,variable,
                                                   DeclareCoeffPhysArrays);
                if(DeclareCoeffPhysArrays)
                {
                    EvaluateBoundaryConditions(0.0, variable);
                }
                ApplyGeomInfo();
                FindPeriodicTraces(bcs,variable);
            }

            if(SetUpJustDG)
            {
                SetUpDG(variable);
            }
            else
            {
                int i,cnt;
                Array<OneD, int> ElmtID, TraceID;
                GetBoundaryToElmtMap(ElmtID, TraceID);

                for(cnt = i = 0; i < m_bndCondExpansions.size(); ++i)
                {
                    MultiRegions::ExpListSharedPtr locExpList;
                    locExpList = m_bndCondExpansions[i];

                    for(int e = 0; e < locExpList->GetExpSize(); ++e)
                    {
                        (*m_exp)[ElmtID[cnt+e]]->SetTraceExp
                            (TraceID[cnt+e], locExpList->GetExp(e));
                        locExpList->GetExp(e)->SetAdjacentElementExp
                            (TraceID[cnt+e], (*m_exp)[ElmtID[cnt+e]]);
                    }
                    cnt += m_bndCondExpansions[i]->GetExpSize();
                }

                if(m_session->DefinesSolverInfo("PROJECTION"))
                {
                    std::string ProjectStr =
                         m_session->GetSolverInfo("PROJECTION");
                    if((ProjectStr == "MixedCGDG") ||
                       (ProjectStr == "Mixed_CG_Discontinuous"))
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

        /**
         * @brief Set up all DG member variables and maps.
         */
        void  DisContField::SetUpDG(const std::string variable)
        {
            // Check for multiple calls
            if (m_trace != NullExpListSharedPtr)
            {
                return;
            }

            // Set up matrix map
            m_globalBndMat = MemoryManager<GlobalLinSysMap>::
                AllocateSharedPtr();

            // Set up trace space
            m_trace = MemoryManager<ExpList>::AllocateSharedPtr
                (m_session, m_bndCondExpansions, m_bndConditions,
                 *m_exp,m_graph);

            PeriodicMap  periodicTraces = (m_expType == e1D)? m_periodicVerts:
                (m_expType == e2D)? m_periodicEdges: m_periodicFaces;

            m_traceMap = MemoryManager<AssemblyMapDG>::
                AllocateSharedPtr(m_session, m_graph, m_trace, *this,
                                  m_bndCondExpansions, m_bndConditions,
                                  periodicTraces, variable);

            if (m_session->DefinesCmdLineArgument("verbose"))
            {
                m_traceMap->PrintStats(std::cout, variable);
            }

            Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >
                &elmtToTrace = m_traceMap->GetElmtToTrace();

            // Scatter trace to elements. For each element, we find
            // the trace  associated to each elemental trace. The element
            // then retains a pointer to the elemental trace, to
            // ensure uniqueness of normals when retrieving from two
            // adjoining elements which do not lie in a plane.

            for (int i = 0; i < m_exp->size(); ++i)
            {
                for (int j = 0; j < (*m_exp)[i]->GetNtraces(); ++j)
                {
                    LocalRegions::ExpansionSharedPtr exp = elmtToTrace[i][j];;

                    exp->SetAdjacentElementExp(j, (*m_exp)[i]);
                    (*m_exp)[i]->SetTraceExp   (j, exp);
                }
            }

            // Set up physical normals
            SetUpPhysNormals();

            int cnt, n;

            // Identify boundary trace
            for(cnt = 0, n = 0; n < m_bndCondExpansions.size(); ++n)
            {
                if (m_bndConditions[n]->GetBoundaryConditionType() !=
                    SpatialDomains::ePeriodic)
                {
                    for(int v = 0; v < m_bndCondExpansions[n]->GetExpSize(); ++v)
                    {
                        m_boundaryTraces.insert
                            (m_traceMap->GetBndCondIDToGlobalTraceID(cnt+v));
                    }
                    cnt += m_bndCondExpansions[n]->GetExpSize();
                }
            }

            // Set up information for periodic boundary conditions.
            std::unordered_map<int,pair<int,int> > perTraceToExpMap;
            for (cnt = n = 0; n < m_exp->size(); ++n)
            {
                for (int v = 0; v < (*m_exp)[n]->GetNtraces(); ++v, ++cnt)
                {
                    auto it = periodicTraces.find
                        ((*m_exp)[n]->GetGeom()->GetTid(v));

                    if (it != periodicTraces.end())
                    {
                        perTraceToExpMap[it->first] = make_pair(n,v);
                    }
                }
            }

            // Set up left-adjacent tracelist.
            m_leftAdjacentTraces.resize(cnt);

            // count size of trace
            for (cnt = n = 0; n < m_exp->size(); ++n)
            {
                for (int v = 0; v < (*m_exp)[n]->GetNtraces(); ++v, ++cnt)
                {
                    m_leftAdjacentTraces[cnt] = IsLeftAdjacentTrace(n, v);
                }
            }


            // Set up mapping to copy Fwd of periodic bcs to Bwd of other edge.
            for (cnt = n = 0; n < m_exp->size(); ++n)
            {
                for (int v = 0; v < (*m_exp)[n]->GetNtraces(); ++v, ++cnt)
                {
                    int GeomId = (*m_exp)[n]->GetGeom()->GetTid(v);

                    // Check to see if this trace is periodic.
                    auto it = periodicTraces.find(GeomId);

                    if (it != periodicTraces.end())
                    {
                        const PeriodicEntity &ent = it->second[0];
                        auto it2 = perTraceToExpMap.find(ent.id);

                        if (it2 == perTraceToExpMap.end())
                        {
                            if (m_session->GetComm()->GetRowComm()->GetSize()
                                > 1 && !ent.isLocal)
                            {
                                continue;
                            }
                            else
                            {
                                ASSERTL1(false, "Periodic trace not found!");
                            }
                        }

                        ASSERTL1(m_leftAdjacentTraces[cnt],
                                 "Periodic trace in non-forward space?");

                        int offset  = m_trace->GetPhys_Offset
                            ((m_traceMap->GetElmtToTrace())
                             [n][v]->GetElmtId());

                        int offset2 = m_trace->GetPhys_Offset
                            ((m_traceMap->GetElmtToTrace())
                             [it2->second.first]
                             [it2->second.second]->GetElmtId());


                        switch(m_expType)
                        {
                        case e1D:
                            {
                                m_periodicFwdCopy.push_back(offset);
                                m_periodicBwdCopy.push_back(offset2);
                            }
                            break;
                        case e2D:
                            {
                                // Calculate relative orientations between trace to
                                // calculate copying map.
                                int nquad = elmtToTrace[n][v]->GetNumPoints(0);

                                vector<int> tmpBwd(nquad);
                                vector<int> tmpFwd(nquad);

                                if (ent.orient == StdRegions::eForwards)
                                {
                                    for (int i = 0; i < nquad; ++i)
                                    {
                                        tmpBwd[i] = offset2 + i;
                                        tmpFwd[i] = offset  + i;
                                    }
                                }
                                else
                                {
                                    for (int i = 0; i < nquad; ++i)
                                    {
                                        tmpBwd[i] = offset2 + i;
                                        tmpFwd[i] = offset  + nquad - i - 1;
                                    }
                                }

                                for (int i = 0; i < nquad; ++i)
                                {
                                    m_periodicFwdCopy.push_back(tmpFwd[i]);
                                    m_periodicBwdCopy.push_back(tmpBwd[i]);
                                }
                            }
                            break;
                        case e3D:
                        {
                            // Calculate relative orientations between faces to
                            // calculate copying map.
                            int nquad1 = elmtToTrace[n][v]->GetNumPoints(0);
                            int nquad2 = elmtToTrace[n][v]->GetNumPoints(1);

                            vector<int> tmpBwd(nquad1*nquad2);
                            vector<int> tmpFwd(nquad1*nquad2);

                            if (ent.orient == StdRegions::eDir1FwdDir2_Dir2FwdDir1 ||
                                ent.orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1 ||
                                ent.orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1 ||
                                ent.orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
                            {
                                for (int i = 0; i < nquad2; ++i)
                                {
                                    for (int j = 0; j < nquad1; ++j)
                                    {
                                        tmpBwd[i*nquad1+j] = offset2 + i*nquad1+j;
                                        tmpFwd[i*nquad1+j] = offset  + j*nquad2+i;
                                    }
                                }
                            }
                            else
                            {
                                for (int i = 0; i < nquad2; ++i)
                                {
                                    for (int j = 0; j < nquad1; ++j)
                                    {
                                        tmpBwd[i*nquad1+j] = offset2 + i*nquad1+j;
                                        tmpFwd[i*nquad1+j] = offset  + i*nquad1+j;
                                    }
                                }
                            }

                            if (ent.orient==StdRegions::eDir1BwdDir1_Dir2FwdDir2 ||
                                ent.orient==StdRegions::eDir1BwdDir1_Dir2BwdDir2 ||
                                ent.orient==StdRegions::eDir1FwdDir2_Dir2BwdDir1 ||
                                ent.orient==StdRegions::eDir1BwdDir2_Dir2BwdDir1)
                            {
                                // Reverse x direction
                                for (int i = 0; i < nquad2; ++i)
                                {
                                    for (int j = 0; j < nquad1/2; ++j)
                                    {
                                        swap(tmpFwd[i*nquad1 + j],
                                             tmpFwd[i*nquad1 + nquad1-j-1]);
                                    }
                                }
                            }

                            if (ent.orient==StdRegions::eDir1FwdDir1_Dir2BwdDir2 ||
                                ent.orient==StdRegions::eDir1BwdDir1_Dir2BwdDir2 ||
                                ent.orient==StdRegions::eDir1BwdDir2_Dir2FwdDir1 ||
                                ent.orient==StdRegions::eDir1BwdDir2_Dir2BwdDir1)
                            {
                                // Reverse y direction
                                for (int j = 0; j < nquad1; ++j)
                                {
                                    for (int i = 0; i < nquad2/2; ++i)
                                    {
                                        swap(tmpFwd[i*nquad1 + j],
                                             tmpFwd[(nquad2-i-1)*nquad1 + j]);
                                    }
                                }
                            }

                            for (int i = 0; i < nquad1*nquad2; ++i)
                            {
                                m_periodicFwdCopy.push_back(tmpFwd[i]);
                                m_periodicBwdCopy.push_back(tmpBwd[i]);
                            }
                        }
                        break;
                        default:
                            ASSERTL1(false,"not set up");
                        }
                    }
                }
            }

            m_locTraceToTraceMap = MemoryManager<LocTraceToTraceMap>::
                AllocateSharedPtr(*this, m_trace, elmtToTrace,
                                  m_leftAdjacentTraces);
        }


        bool DisContField::IsLeftAdjacentTrace(const int n, const int e)
        {
            LocalRegions::ExpansionSharedPtr traceEl =
                m_traceMap->GetElmtToTrace()[n][e];

            PeriodicMap  periodicTraces = (m_expType == e1D)? m_periodicVerts:
                (m_expType == e2D)? m_periodicEdges: m_periodicFaces;

            bool fwd = true;
            if (traceEl->GetLeftAdjacentElementTrace () == -1 ||
                traceEl->GetRightAdjacentElementTrace() == -1)
            {
                // Boundary edge (1 connected element). Do nothing in
                // serial.
                auto it = m_boundaryTraces.find(traceEl->GetElmtId());

                // If the edge does not have a boundary condition set on
                // it, then assume it is a partition edge or periodic.
                if (it == m_boundaryTraces.end())
                {
                    fwd = true;
                }
            }
            else if (traceEl->GetLeftAdjacentElementTrace () != -1 &&
                     traceEl->GetRightAdjacentElementTrace() != -1)
            {
                // Non-boundary edge (2 connected elements).
                fwd = (traceEl->GetLeftAdjacentElementExp().get())
                    == (*m_exp)[n].get();
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
        SpatialDomains::BoundaryConditionsSharedPtr DisContField::GetDomainBCs
             (const SpatialDomains::CompositeMap &domain,
              const SpatialDomains::BoundaryConditions &Allbcs,
              const std::string &variable)
        {
            SpatialDomains::BoundaryConditionsSharedPtr returnval;

            returnval = MemoryManager<SpatialDomains::BoundaryConditions>::AllocateSharedPtr();

            map<int,int> GeometryToRegionsMap;

            const SpatialDomains::BoundaryRegionCollection &bregions
                = Allbcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions
                = Allbcs.GetBoundaryConditions();

            // Set up a map of all boundary regions
            for(auto &it : bregions)
            {
                for (auto &bregionIt : *it.second)
                {
                    // can assume that all regions only contain one point in 1D
                    // Really do not need loop above
                    int id = bregionIt.second->m_geomVec[0]->GetGlobalID();
                    GeometryToRegionsMap[id] = it.first;
                }
            }

            map<int,SpatialDomains::GeometrySharedPtr> EndOfDomain;

            // Now find out which points in domain have only one vertex
            for(auto &domIt : domain)
            {
                SpatialDomains::CompositeSharedPtr geomvector = domIt.second;
                for(int i = 0; i < geomvector->m_geomVec.size(); ++i)
                {
                    for(int j = 0; j < 2; ++j)
                    {
                        int vid = geomvector->m_geomVec[i]->GetVid(j);
                        if(EndOfDomain.count(vid) == 0)
                        {
                            EndOfDomain[vid] = geomvector->m_geomVec[i]->GetVertex(j);
                        }
                        else
                        {
                            EndOfDomain.erase(vid);
                        }
                    }
                }
            }
            ASSERTL1(EndOfDomain.size() == 2,"Did not find two ends of domain");

            int numNewBc = 1;
            for(auto &regIt : EndOfDomain)
            {
                if(GeometryToRegionsMap.count(regIt.first) != 0)
                {
                    // Set up boundary condition up
                    auto iter = GeometryToRegionsMap.find(regIt.first);
                    ASSERTL1(iter != GeometryToRegionsMap.end(),
                             "Failied to find GeometryToRegionMap");

                    int regionId = iter->second;
                    auto bregionsIter = bregions.find(regionId);
                    ASSERTL1(bregionsIter != bregions.end(),
                             "Failed to find boundary region");

                    SpatialDomains::BoundaryRegionShPtr breg =
                        bregionsIter->second;
                    returnval->AddBoundaryRegions(regionId, breg);

                    auto bconditionsIter = bconditions.find(regionId);
                    ASSERTL1(bconditionsIter != bconditions.end(),
                             "Failed to find boundary collection");

                    SpatialDomains::BoundaryConditionMapShPtr bcond =
                        bconditionsIter->second;
                    returnval->AddBoundaryConditions(regionId,bcond);
                }
                else // Set up an undefined region.
                {
                    SpatialDomains::BoundaryRegionShPtr breg(MemoryManager<SpatialDomains::BoundaryRegion>::AllocateSharedPtr());

                    // Set up Composite (GemetryVector) to contain vertex and put into bRegion
                    SpatialDomains::CompositeSharedPtr gvec =
                        MemoryManager<SpatialDomains::Composite>
                        ::AllocateSharedPtr();
                    gvec->m_geomVec.push_back(regIt.second);
                    (*breg)[regIt.first] = gvec;

                    returnval->AddBoundaryRegions(bregions.size()+numNewBc,breg);

                    SpatialDomains::BoundaryConditionMapShPtr bCondition = MemoryManager<SpatialDomains::BoundaryConditionMap>::AllocateSharedPtr();

                    // Set up just boundary condition for this variable.
                    SpatialDomains::BoundaryConditionShPtr notDefinedCondition(MemoryManager<SpatialDomains::NotDefinedBoundaryCondition>::AllocateSharedPtr(m_session, "0"));
                    (*bCondition)[variable] = notDefinedCondition;

                    returnval->AddBoundaryConditions(bregions.size()+numNewBc,bCondition);
		    ++numNewBc;

                }
            }

            return returnval;
        }

        /**
         * Constructor for use in multidomain computations where a
         * domain list can be passed instead of graph1D
         *
         * @param	domain	Subdomain specified in the inputfile from
         *       	      	which the DisContField is set up
         */
        DisContField::DisContField
              (const LibUtilities::SessionReaderSharedPtr &pSession,
               const SpatialDomains::MeshGraphSharedPtr   &graph1D,
               const SpatialDomains::CompositeMap         &domain,
               const SpatialDomains::BoundaryConditions   &Allbcs,
               const std::string                          &variable,
               bool SetToOneSpaceDimension,
               const Collections::ImplementationType ImpType):
                  ExpList(pSession,domain,graph1D,true,variable,
                          SetToOneSpaceDimension, pSession->GetComm(),ImpType)
        {
            if (variable.compare("DefaultVar") != 0)
            {
                SpatialDomains::BoundaryConditionsSharedPtr
                    DomBCs = GetDomainBCs(domain,Allbcs,variable);

                GenerateBoundaryConditionExpansion(m_graph,*DomBCs,variable);
                EvaluateBoundaryConditions(0.0, variable);
                ApplyGeomInfo();
                FindPeriodicTraces(*DomBCs,variable);
            }

            SetUpDG(variable);
        }

        /**
         * Constructs a field as a copy of an existing field.
         * @param   In          Existing DisContField object to copy.
         */
        DisContField::DisContField(const DisContField &In,
                                   const bool DeclareCoeffPhysArrays):
            ExpList(In,DeclareCoeffPhysArrays),
            m_bndCondExpansions(In.m_bndCondExpansions),
            m_bndConditions(In.m_bndConditions),
            m_globalBndMat(In.m_globalBndMat),
            m_traceMap(In.m_traceMap),
            m_boundaryTraces(In.m_boundaryTraces),
            m_periodicVerts(In.m_periodicVerts),
            m_periodicFwdCopy(In.m_periodicFwdCopy),
            m_periodicBwdCopy(In.m_periodicBwdCopy),
            m_leftAdjacentTraces(In.m_leftAdjacentTraces),
            m_locTraceToTraceMap (In.m_locTraceToTraceMap)
        {
            if (In.m_trace)
            {
                m_trace = MemoryManager<ExpList>::AllocateSharedPtr
                    (*In.m_trace,DeclareCoeffPhysArrays);
            }
        }


        /*
         * @brief Copy type constructor which declares new boundary conditions
         * and re-uses mapping info and trace space if possible
         */
        DisContField::DisContField(
            const DisContField                       &In,
            const SpatialDomains::MeshGraphSharedPtr &graph,
            const std::string                        &variable,
            const bool                                SetUpJustDG,
            const bool                                DeclareCoeffPhysArrays):
            ExpList(In,DeclareCoeffPhysArrays)
        {

            m_trace = NullExpListSharedPtr;

            // Set up boundary conditions for this variable.
            // Do not set up BCs if default variable
            if(variable.compare("DefaultVar") != 0)
            {
                SpatialDomains::BoundaryConditions bcs(m_session, graph);
                GenerateBoundaryConditionExpansion(graph, bcs, variable);

                if (DeclareCoeffPhysArrays)
                {
                    EvaluateBoundaryConditions(0.0, variable);
                }

                if (!SameTypeOfBoundaryConditions(In))
                {
                    // Find periodic edges for this variable.
                    FindPeriodicTraces(bcs, variable);

                    if(SetUpJustDG)
                    {
                        SetUpDG(variable);
                    }
                    else
                    {
                        // set elmt edges to point to robin bc edges if required
                        int i, cnt = 0;
                        Array<OneD, int> ElmtID,TraceID;
                        GetBoundaryToElmtMap(ElmtID,TraceID);

                        for (i = 0; i < m_bndCondExpansions.size(); ++i)
                        {
                            MultiRegions::ExpListSharedPtr locExpList;

                            int e;
                            locExpList = m_bndCondExpansions[i];

                            for(e = 0; e < locExpList->GetExpSize(); ++e)
                            {
                                (*m_exp)[ElmtID[cnt+e]]->SetTraceExp
                                    (TraceID[cnt+e], locExpList->GetExp(e));
                                locExpList->GetExp(e)->SetAdjacentElementExp
                                    (TraceID[cnt+e], (*m_exp)[ElmtID[cnt+e]]);
                            }

                            cnt += m_bndCondExpansions[i]->GetExpSize();
                        }


                        if (m_session->DefinesSolverInfo("PROJECTION"))
                        {
                            std::string ProjectStr =
                                m_session->GetSolverInfo("PROJECTION");

                            if ((ProjectStr == "MixedCGDG") ||
                                (ProjectStr == "Mixed_CG_Discontinuous"))
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
                    m_globalBndMat        = In.m_globalBndMat;
                    m_trace               = In.m_trace;
                    m_traceMap            = In.m_traceMap;
                    m_locTraceToTraceMap  = In.m_locTraceToTraceMap;
                    m_periodicVerts       = In.m_periodicVerts;
                    m_periodicEdges       = In.m_periodicEdges;
                    m_periodicFaces       = In.m_periodicFaces;
                    m_periodicFwdCopy     = In.m_periodicFwdCopy;
                    m_periodicBwdCopy     = In.m_periodicBwdCopy;
                    m_boundaryTraces      = In.m_boundaryTraces;
                    m_leftAdjacentTraces  = In.m_leftAdjacentTraces;

                    if (SetUpJustDG == false)
                    {
                        // set elmt edges to point to robin bc edges if required
                        int i, cnt = 0;
                        Array<OneD, int> ElmtID, TraceID;
                        GetBoundaryToElmtMap(ElmtID, TraceID);

                        for (i = 0; i < m_bndCondExpansions.size(); ++i)
                        {
                            MultiRegions::ExpListSharedPtr locExpList;

                            int e;
                            locExpList = m_bndCondExpansions[i];

                            for (e = 0; e < locExpList->GetExpSize(); ++e)
                            {
                                (*m_exp)[ElmtID[cnt+e]]->SetTraceExp
                                    (TraceID[cnt+e], locExpList->GetExp(e));
                                locExpList->GetExp(e)->SetAdjacentElementExp
                                    (TraceID[cnt+e], (*m_exp)[ElmtID[cnt+e]]);
                            }
                            cnt += m_bndCondExpansions[i]->GetExpSize();
                        }

                        SetUpPhysNormals();
                    }
                }
            }
        }

        /**
         * Constructs a field as a copy of an existing explist field.
         * @param   In          Existing ExpList object to copy.
         */
        DisContField::DisContField(const ExpList &In):
            ExpList(In)
	{
	}

        /**
         *
         */
        DisContField::~DisContField()
        {
        }


        /**
         * \brief This function discretises the boundary conditions by setting
         * up a list of one-dimensions lower boundary expansions.
         *
         * According to their boundary region, the separate  boundary
         * expansions are bundled together in an object of the class
         *
         * @param   graph       A mesh, containing information about the domain
         *                      and the spectral/hp element expansions.
         * @param   bcs         Information about the enforced boundary
         *                      conditions.
         * @param   variable    The session variable associated with the
         *                      boundary conditions to enforce.
         * @param DeclareCoeffPhysArrays bool to identify if array
         *                               space should be setup.
         *                               Default is true.
         */
        void DisContField::GenerateBoundaryConditionExpansion(
                                                              const SpatialDomains::MeshGraphSharedPtr &graph,
                                                              const SpatialDomains::BoundaryConditions &bcs,
                                                              const std::string variable,
                                                              const bool DeclareCoeffPhysArrays)
        {
            int cnt  = 0;
            SpatialDomains::BoundaryConditionShPtr            bc;
            MultiRegions::ExpListSharedPtr                    locExpList;
            const SpatialDomains::BoundaryRegionCollection &bregions
                = bcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions
                = bcs.GetBoundaryConditions();

            m_bndCondExpansions =
                Array<OneD, MultiRegions::ExpListSharedPtr>(bregions.size());
            m_bndConditions     =
                Array<OneD, SpatialDomains::BoundaryConditionShPtr>(bregions.size());

            m_bndCondBndWeight = Array<OneD, NekDouble> {bregions.size(), 0.0};

            // count the number of non-periodic boundary points
            for (auto &it : bregions)
            {
                bc = GetBoundaryCondition(bconditions, it.first, variable);

                locExpList = MemoryManager<MultiRegions::ExpList>
                    ::AllocateSharedPtr(m_session, *(it.second), graph,
                                        DeclareCoeffPhysArrays, variable,
                                        false, bc->GetComm());

                m_bndCondExpansions[cnt]  = locExpList;
                m_bndConditions[cnt]      = bc;

                std::string type = m_bndConditions[cnt]->GetUserDefined();

                // Set up normals on non-Dirichlet boundary conditions. Second
                // two conditions ideally should be in local solver setup (when
                // made into factory)
                if(bc->GetBoundaryConditionType() != SpatialDomains::eDirichlet
                   || boost::iequals(type,"I") || boost::iequals(type,"CalcBC"))
                {
                    SetUpPhysNormals();
                }
                cnt ++;
            }
        }


        /**
         * @brief Determine the periodic faces, edges and vertices for the given
         * graph.
         *
         * @param   bcs         Information about the boundary conditions.
         * @param   variable    Specifies the field.
         */
        void DisContField::FindPeriodicTraces(
                                              const SpatialDomains::BoundaryConditions &bcs,
                                              const std::string variable)
        {
            const SpatialDomains::BoundaryRegionCollection &bregions
                = bcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions
                = bcs.GetBoundaryConditions();

            LibUtilities::CommSharedPtr vComm =
                m_session->GetComm()->GetRowComm();


            switch(m_expType)
            {
            case e1D:
                {
                    int i, region1ID, region2ID;

                    SpatialDomains::BoundaryConditionShPtr locBCond;
                    map<int,int> BregionToVertMap;

                    // Construct list of all periodic Region and their
                    // global vertex on this process.
                    for (auto &it : bregions)
                    {
                        locBCond = GetBoundaryCondition(bconditions,
                                                        it.first, variable);

                        if (locBCond->GetBoundaryConditionType()
                            != SpatialDomains::ePeriodic)
                        {
                            continue;
                        }
                        int id = it.second->begin()->second->
                            m_geomVec[0]->GetGlobalID();

                        BregionToVertMap[it.first] = id;
                    }

                    set<int> islocal;

                    int n = vComm->GetSize();
                    int p = vComm->GetRank();

                    Array<OneD, int> nregions(n, 0);
                    nregions[p] = BregionToVertMap.size();
                    vComm->AllReduce(nregions, LibUtilities::ReduceSum);

                    int totRegions = Vmath::Vsum(n, nregions, 1);

                    Array<OneD, int> regOffset(n, 0);

                    for (i = 1; i < n; ++i)
                    {
                        regOffset[i] = regOffset[i-1] + nregions[i-1];
                    }

                    Array<OneD, int> bregmap(totRegions, 0);
                    Array<OneD, int> bregid (totRegions, 0);

                    i = regOffset[p];
                    for (auto &iit : BregionToVertMap)
                    {
                        bregid [i  ] = iit.first;
                        bregmap[i++] = iit.second;
                        islocal.insert(iit.first);
                    }

                    vComm->AllReduce(bregmap, LibUtilities::ReduceSum);
                    vComm->AllReduce(bregid,  LibUtilities::ReduceSum);

                    for (int i = 0; i < totRegions; ++i)
                    {
                        BregionToVertMap[bregid[i]] = bregmap[i];
                    }

                    // Construct list of all periodic pairs local to this process.
                    for (auto &it : bregions)
                    {
                        locBCond = GetBoundaryCondition(bconditions,
                                                        it.first, variable);

                        if (locBCond->GetBoundaryConditionType()
                            != SpatialDomains::ePeriodic)
                        {
                            continue;
                        }

                        // Identify periodic boundary region IDs.
                        region1ID = it.first;
                        region2ID = std::static_pointer_cast<
                            SpatialDomains::PeriodicBoundaryCondition>(
                                                                       locBCond)->m_connectedBoundaryRegion;

                        ASSERTL0(BregionToVertMap.count(region1ID) != 0,
                                 "Cannot determine vertex of region1ID");

                        ASSERTL0(BregionToVertMap.count(region2ID) != 0,
                                 "Cannot determine vertex of region2ID");

                        PeriodicEntity ent(BregionToVertMap[region2ID],
                                           StdRegions::eNoOrientation,
                                           islocal.count(region2ID) != 0);

                        m_periodicVerts[BregionToVertMap[region1ID]].push_back(ent);
                    }
                }
                break;
            case e2D:
                {
                    int region1ID, region2ID, i, j, k, cnt;
                    SpatialDomains::BoundaryConditionShPtr locBCond;

                    SpatialDomains::CompositeOrdering compOrder   =
                        m_graph->GetCompositeOrdering();
                    SpatialDomains::BndRegionOrdering bndRegOrder =
                        m_graph->GetBndRegionOrdering();
                    SpatialDomains::CompositeMap      compMap     =
                        m_graph->GetComposites();

                    // Unique collection of pairs of periodic composites
                    // (i.e. if composites 1 and 2 are periodic then this
                    // map will contain either the pair (1,2) or (2,1) but
                    // not both).
                    map<int,int>                                 perComps;
                    map<int,vector<int>>                         allVerts;
                    set<int>                                     locVerts;
                    map<int, pair<int, StdRegions::Orientation>> allEdges;

                    // Set up a set of all local verts and edges.
                    for(i = 0; i < (*m_exp).size(); ++i)
                    {
                        for(j = 0; j < (*m_exp)[i]->GetNverts(); ++j)
                        {
                            int id = (*m_exp)[i]->GetGeom()->GetVid(j);
                            locVerts.insert(id);
                        }
                    }

                    // Construct list of all periodic pairs local to this process.
                    for (auto &it : bregions)
                    {
                        locBCond = GetBoundaryCondition(
                                                        bconditions, it.first, variable);

                        if (locBCond->GetBoundaryConditionType()
                            != SpatialDomains::ePeriodic)
                        {
                            continue;
                        }

                        // Identify periodic boundary region IDs.
                        region1ID = it.first;
                        region2ID = std::static_pointer_cast<
                            SpatialDomains::PeriodicBoundaryCondition>(
                                                                       locBCond)->m_connectedBoundaryRegion;

                        // From this identify composites. Note that in
                        // serial this will be an empty map.
                        int cId1, cId2;
                        if (vComm->GetSize() == 1)
                        {
                            cId1 = it.second->begin()->first;
                            cId2 = bregions.find(region2ID)->second->begin()->first;
                        }
                        else
                        {
                            cId1 = bndRegOrder.find(region1ID)->second[0];
                            cId2 = bndRegOrder.find(region2ID)->second[0];
                        }

                        ASSERTL0(it.second->size() == 1,
                                 "Boundary region "+boost::lexical_cast<string>(
                                                                                region1ID)+" should only contain 1 composite.");

                        vector<unsigned int> tmpOrder;

                        // Construct set containing all periodic edgesd on
                        // this process
                        SpatialDomains::CompositeSharedPtr
                            c = it.second->begin()->second;

                        for (i = 0; i < c->m_geomVec.size(); ++i)
                        {
                            SpatialDomains::SegGeomSharedPtr segGeom =
                                std::dynamic_pointer_cast<
                                    SpatialDomains::SegGeom>(c->m_geomVec[i]);
                            ASSERTL0(segGeom, "Unable to cast to shared ptr");

                            SpatialDomains::GeometryLinkSharedPtr elmt =
                                m_graph->GetElementsFromEdge(segGeom);
                            ASSERTL0(elmt->size() == 1,
                                     "The periodic boundaries belong to "
                                     "more than one element of the mesh");

                            SpatialDomains::Geometry2DSharedPtr geom =
                                std::dynamic_pointer_cast
                                <SpatialDomains::Geometry2D>(elmt->at(0).first);

                            allEdges[c->m_geomVec[i]->GetGlobalID()] =
                                make_pair(elmt->at(0).second,
                                          geom->GetEorient(elmt->at(0).second));

                            // In serial mesh partitioning will not have occurred so
                            // need to fill composite ordering map manually.
                            if (vComm->GetSize() == 1)
                            {
                                tmpOrder.push_back(c->m_geomVec[i]->GetGlobalID());
                            }

                            vector<int> vertList(2);
                            vertList[0] = segGeom->GetVid(0);
                            vertList[1] = segGeom->GetVid(1);
                            allVerts[c->m_geomVec[i]->GetGlobalID()] = vertList;
                        }

                        if (vComm->GetSize() == 1)
                        {
                            compOrder[it.second->begin()->first] = tmpOrder;
                        }


                        // See if we already have either region1 or
                        // region2 stored in perComps map.
                        if (perComps.count(cId1) == 0)
                        {
                            if (perComps.count(cId2) == 0)
                            {
                                perComps[cId1] = cId2;
                            }
                            else
                            {
                                std::stringstream ss;
                                ss << "Boundary region " << cId2 << " should be "
                                   << "periodic with " << perComps[cId2] << " but "
                                   << "found " << cId1 << " instead!";
                                ASSERTL0(perComps[cId2] == cId1, ss.str());
                            }
                        }
                        else
                        {
                            std::stringstream ss;
                            ss << "Boundary region " << cId1 << " should be "
                               << "periodic with " << perComps[cId1] << " but "
                               << "found " << cId2 << " instead!";
                            ASSERTL0(perComps[cId1] == cId1, ss.str());
                        }
                    }


                    // Process local edge list to obtain relative edge orientations.
                    int              n = vComm->GetSize();
                    int              p = vComm->GetRank();
                    int              totEdges;
                    Array<OneD, int> edgecounts(n,0);
                    Array<OneD, int> edgeoffset(n,0);
                    Array<OneD, int> vertoffset(n,0);

                    edgecounts[p] = allEdges.size();
                    vComm->AllReduce(edgecounts, LibUtilities::ReduceSum);

                    edgeoffset[0] = 0;
                    for (i = 1; i < n; ++i)
                    {
                        edgeoffset[i] = edgeoffset[i-1] + edgecounts[i-1];
                    }

                    totEdges = Vmath::Vsum(n, edgecounts, 1);
                    Array<OneD, int> edgeIds   (totEdges, 0);
                    Array<OneD, int> edgeIdx   (totEdges, 0);
                    Array<OneD, int> edgeOrient(totEdges, 0);
                    Array<OneD, int> edgeVerts (totEdges, 0);

                    auto sIt = allEdges.begin();

                    for (i = 0; sIt != allEdges.end(); ++sIt)
                    {
                        edgeIds   [edgeoffset[p] + i  ] = sIt->first;
                        edgeIdx   [edgeoffset[p] + i  ] = sIt->second.first;
                        edgeOrient[edgeoffset[p] + i  ] = sIt->second.second;
                        edgeVerts [edgeoffset[p] + i++] = allVerts[sIt->first].size();
                    }

                    vComm->AllReduce(edgeIds,    LibUtilities::ReduceSum);
                    vComm->AllReduce(edgeIdx,    LibUtilities::ReduceSum);
                    vComm->AllReduce(edgeOrient, LibUtilities::ReduceSum);
                    vComm->AllReduce(edgeVerts,  LibUtilities::ReduceSum);

                    // Calculate number of vertices on each processor.
                    Array<OneD, int> procVerts(n,0);
                    int nTotVerts;

                    // Note if there are no periodic edges at all calling Vsum will
                    // cause a segfault.
                    if (totEdges > 0)
                    {
                        nTotVerts = Vmath::Vsum(totEdges, edgeVerts, 1);
                    }
                    else
                    {
                        nTotVerts = 0;
                    }

                    for (i = 0; i < n; ++i)
                    {
                        if (edgecounts[i] > 0)
                        {
                            procVerts[i] = Vmath::Vsum(
                                                       edgecounts[i], edgeVerts + edgeoffset[i], 1);
                        }
                        else
                        {
                            procVerts[i] = 0;
                        }
                    }
                    vertoffset[0] = 0;

                    for (i = 1; i < n; ++i)
                    {
                        vertoffset[i] = vertoffset[i-1] + procVerts[i-1];
                    }

                    Array<OneD, int> traceIds(nTotVerts, 0);
                    for (i = 0, sIt = allEdges.begin(); sIt != allEdges.end(); ++sIt)
                    {
                        for (j = 0; j < allVerts[sIt->first].size(); ++j)
                        {
                            traceIds[vertoffset[p] + i++] = allVerts[sIt->first][j];
                        }
                    }

                    vComm->AllReduce(traceIds, LibUtilities::ReduceSum);

                    // For simplicity's sake create a map of edge id -> orientation.
                    map<int, StdRegions::Orientation> orientMap;
                    map<int, vector<int> >            vertMap;

                    for (cnt = i = 0; i < totEdges; ++i)
                    {
                        ASSERTL0(orientMap.count(edgeIds[i]) == 0,
                                 "Already found edge in orientation map!");

                        // Work out relative orientations. To avoid having
                        // to exchange vertex locations, we figure out if
                        // the edges are backwards or forwards orientated
                        // with respect to a forwards orientation that is
                        // CCW. Since our local geometries are
                        // forwards-orientated with respect to the
                        // Cartesian axes, we need to invert the
                        // orientation for the top and left edges of a
                        // quad and the left edge of a triangle.
                        StdRegions::Orientation o =
                            (StdRegions::Orientation)edgeOrient[i];

                        if (edgeIdx[i] > 1)
                        {
                            o = o == StdRegions::eForwards ?
                                StdRegions::eBackwards : StdRegions::eForwards;
                        }

                            orientMap[edgeIds[i]] = o;

                        vector<int> verts(edgeVerts[i]);

                        for (j = 0; j < edgeVerts[i]; ++j)
                        {
                            verts[j] = traceIds[cnt++];
                        }
                        vertMap[edgeIds[i]] = verts;
                    }

                    // Go through list of composites and figure out which
                    // edges are parallel from original ordering in
                    // session file. This includes composites which are
                    // not necessarily on this process.
                    map<int,int> allCompPairs;

                    // Store temporary map of periodic vertices which will hold all
                    // periodic vertices on the entire mesh so that doubly periodic
                    // vertices can be counted properly across partitions. Local
                    // vertices are copied into m_periodicVerts at the end of the
                    // function.
                    PeriodicMap periodicVerts;

                    for (auto &cIt : perComps)
                    {
                        SpatialDomains::CompositeSharedPtr c[2];
                        const int   id1  = cIt.first;
                        const int   id2  = cIt.second;
                        std::string id1s = boost::lexical_cast<string>(id1);
                        std::string id2s = boost::lexical_cast<string>(id2);

                        if (compMap.count(id1) > 0)
                        {
                            c[0] = compMap[id1];
                        }

                        if (compMap.count(id2) > 0)
                        {
                            c[1] = compMap[id2];
                        }

                        ASSERTL0(c[0] || c[1],
                                 "Both composites not found on this process!");

                        // Loop over composite ordering to construct list of all
                        // periodic edges regardless of whether they are on this
                        // process.
                        map<int,int> compPairs;

                        ASSERTL0(compOrder.count(id1) > 0,
                                 "Unable to find composite "+id1s+" in order map.");
                        ASSERTL0(compOrder.count(id2) > 0,
                                 "Unable to find composite "+id2s+" in order map.");
                        ASSERTL0(compOrder[id1].size() == compOrder[id2].size(),
                                 "Periodic composites "+id1s+" and "+id2s+
                                 " should have the same number of elements.");
                        ASSERTL0(compOrder[id1].size() > 0,
                                 "Periodic composites "+id1s+" and "+id2s+
                                 " are empty!");

                        // TODO: Add more checks.
                        for (i = 0; i < compOrder[id1].size(); ++i)
                        {
                            int eId1 = compOrder[id1][i];
                            int eId2 = compOrder[id2][i];

                            ASSERTL0(compPairs.count(eId1) == 0,
                                     "Already paired.");

                            if (compPairs.count(eId2) != 0)
                            {
                                ASSERTL0(compPairs[eId2] == eId1, "Pairing incorrect");
                            }
                            compPairs[eId1] = eId2;
                        }

                        // Construct set of all edges that we have locally on this
                        // processor.
                        set<int> locEdges;
                        for (i = 0; i < 2; ++i)
                        {
                            if (!c[i])
                            {
                                continue;
                            }

                            if (c[i]->m_geomVec.size() > 0)
                            {
                                for (j = 0; j < c[i]->m_geomVec.size(); ++j)
                                {
                                    locEdges.insert(c[i]->m_geomVec[j]->GetGlobalID());
                                }
                            }
                        }

                        // Loop over all edges in the geometry composite.
                        for (auto &pIt : compPairs)
                        {
                            int  ids  [2] = {pIt.first, pIt.second};
                            bool local[2] = {locEdges.count(pIt.first) > 0,
                                locEdges.count(pIt.second) > 0};

                            ASSERTL0(orientMap.count(ids[0]) > 0 &&
                                     orientMap.count(ids[1]) > 0,
                                     "Unable to find edge in orientation map.");

                            allCompPairs[pIt.first ] = pIt.second;
                            allCompPairs[pIt.second] = pIt.first;

                            for (i = 0; i < 2; ++i)
                            {
                                if (!local[i])
                                {
                                    continue;
                                }

                                int other = (i+1) % 2;

                                StdRegions::Orientation o =
                                    orientMap[ids[i]] == orientMap[ids[other]] ?
                                    StdRegions::eBackwards :
                                    StdRegions::eForwards;

                                PeriodicEntity ent(ids  [other], o,
                                                   local[other]);
                                m_periodicEdges[ids[i]].push_back(ent);
                            }

                            for (i = 0; i < 2; ++i)
                            {
                                int other = (i+1) % 2;

                                StdRegions::Orientation o =
                                    orientMap[ids[i]] == orientMap[ids[other]] ?
                                    StdRegions::eBackwards :
                                    StdRegions::eForwards;

                                // Determine periodic vertices.
                                vector<int> perVerts1 = vertMap[ids[i]];
                                vector<int> perVerts2 = vertMap[ids[other]];

                                map<int, pair<int, bool> > tmpMap;
                                if (o == StdRegions::eForwards)
                                {
                                    tmpMap[perVerts1[0]] = make_pair(perVerts2[0],
                                                                     locVerts.count(perVerts2[0]) > 0);
                                    tmpMap[perVerts1[1]] = make_pair(perVerts2[1],
                                                                     locVerts.count(perVerts2[1]) > 0);
                                }
                                else
                                {
                                    tmpMap[perVerts1[0]] = make_pair(perVerts2[1],
                                                                     locVerts.count(perVerts2[1]) > 0);
                                    tmpMap[perVerts1[1]] = make_pair(
                                                                     perVerts2[0], locVerts.count(perVerts2[0]) > 0);
                                }

                                for (auto &mIt : tmpMap)
                                {
                                    // See if this vertex has been recorded already.
                                    PeriodicEntity ent2(mIt.second.first,
                                                        StdRegions::eNoOrientation,
                                                        mIt.second.second);
                                    auto perIt = periodicVerts.find(mIt.first);

                                    if (perIt == periodicVerts.end())
                                    {
                                        periodicVerts[mIt.first].push_back(ent2);
                                        perIt = periodicVerts.find(mIt.first);
                                    }
                                    else
                                    {
                                        bool doAdd = true;
                                        for (j = 0; j < perIt->second.size(); ++j)
                                        {
                                            if (perIt->second[j].id ==
                                                mIt.second.first)
                                            {
                                                doAdd = false;
                                                break;
                                            }
                                        }

                                        if (doAdd)
                                        {
                                            perIt->second.push_back(ent2);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // Distribute the size of the periodic boundary to all 
                    // processes. We assume that all processes who own periodic
                    // faces have the same local copy of allCompPairs
                    int NPairs  = allCompPairs.size();
                    vComm->AllReduce(NPairs, LibUtilities::ReduceMax);

                    // Check that the previous assertion regarding allCompPairs
                    // is correct
                    ASSERTL0(
                        allCompPairs.size() == NPairs || allCompPairs.size() == 0,
                        "Local copy of allCompPairs not the same for all ranks"
                    );

                    // Allocate local vectors that will contain the content of
                    // allCompPairs if process owns faces on periodic boundary
                    Array<OneD, int> first(NPairs, -1);
                    Array<OneD, int> second(NPairs, -1);
                    cnt = 0;
                    for(const auto &it : allCompPairs)
                    {
                        first[cnt]    = it.first;
                        second[cnt++] = it.second;
                    }
                    
                    // Distribute the content in first and second to all processes
                    vComm->AllReduce(first,  LibUtilities::ReduceMax);
                    vComm->AllReduce(second, LibUtilities::ReduceMax);

                    // Check that the MPI Allreduce routine worked
                    ASSERTL0(std::count(first.begin(), first.end(), -1) == 0, 
                        "Distribution of allCompPairs failed");
                    ASSERTL0(std::count(second.begin(), second.end(), -1) == 0, 
                        "Distribution of allCompParis failed")

                    // Put content back in allCompPairs
                    allCompPairs.clear();
                    for(cnt = 0; cnt < NPairs; ++cnt)
                    {
                        allCompPairs[first[cnt]] = second[cnt];
                    }

                    // Search for periodic vertices and edges which are not in
                    // a periodic composite but lie in this process. First, loop
                    // over all information we have from other processors.
                    for (cnt = i = 0; i < totEdges; ++i)
                    {
                        int edgeId    = edgeIds[i];

                        ASSERTL0(allCompPairs.count(edgeId) > 0,
                                 "Unable to find matching periodic edge.");

                        int perEdgeId = allCompPairs[edgeId];

                        for (j = 0; j < edgeVerts[i]; ++j, ++cnt)
                        {
                            int vId = traceIds[cnt];

                            auto perId = periodicVerts.find(vId);

                            if (perId == periodicVerts.end())
                            {
                                // This vertex is not included in the
                                // map. Figure out which vertex it is
                                // supposed to be periodic with. perEdgeId
                                // is the edge ID which is periodic with
                                // edgeId. The logic is much the same as
                                // the loop above.
                                int perVertexId =
                                    orientMap[edgeId] == orientMap[perEdgeId] ?
                                    vertMap[perEdgeId][(j+1)%2]:
                                    vertMap[perEdgeId][j];

                                PeriodicEntity ent(perVertexId,
                                                   StdRegions::eNoOrientation,
                                                   locVerts.count(perVertexId) > 0);

                                periodicVerts[vId].push_back(ent);
                            }
                        }
                    }

                    // Loop over all periodic vertices to complete connectivity
                    // information.
                    for (auto &perIt : periodicVerts)
                    {
                        // Loop over associated vertices.
                        for (i = 0; i < perIt.second.size(); ++i)
                        {
                            auto perIt2 = periodicVerts.find(perIt.second[i].id);
                            ASSERTL0(perIt2 != periodicVerts.end(),
                                     "Couldn't find periodic vertex.");

                            for (j = 0; j < perIt2->second.size(); ++j)
                            {
                                if (perIt2->second[j].id == perIt.first)
                                {
                                    continue;
                                }

                                bool doAdd = true;
                                for (k = 0; k < perIt.second.size(); ++k)
                                {
                                    if (perIt2->second[j].id == perIt.second[k].id)
                                    {
                                        doAdd = false;
                                        break;
                                    }
                                }

                                if (doAdd)
                                {
                                    perIt.second.push_back(perIt2->second[j]);
                                }
                            }
                        }
                    }

                    // Do one final loop over periodic vertices to remove non-local
                    // vertices from map.
                    for (auto &perIt : periodicVerts)
                    {
                        if (locVerts.count(perIt.first) > 0)
                        {
                            m_periodicVerts.insert(perIt);
                        }
                    }
                }
                break;
            case e3D:
                {
                    SpatialDomains::CompositeOrdering compOrder   =
                        m_graph->GetCompositeOrdering();
                    SpatialDomains::BndRegionOrdering bndRegOrder =
                        m_graph->GetBndRegionOrdering();
                    SpatialDomains::CompositeMap      compMap     =
                        m_graph->GetComposites();

                    // perComps: Stores a unique collection of pairs of periodic
                    // composites (i.e. if composites 1 and 2 are periodic then this map
                    // will contain either the pair (1,2) or (2,1) but not both).
                    //
                    // The four maps allVerts, allCoord, allEdges and allOrient map a
                    // periodic face to a vector containing the vertex ids of the face;
                    // their coordinates; the edge ids of the face; and their
                    // orientation within that face respectively.
                    //
                    // Finally the three sets locVerts, locEdges and locFaces store any
                    // vertices, edges and faces that belong to a periodic composite and
                    // lie on this process.
                    map<int,RotPeriodicInfo>                       rotComp;
                    map<int,int>                                   perComps;
                    map<int,vector<int> >                          allVerts;
                    map<int,SpatialDomains::PointGeomVector>       allCoord;
                    map<int,vector<int> >                          allEdges;
                    map<int,vector<StdRegions::Orientation> >      allOrient;
                    set<int>                                       locVerts;
                    set<int>                                       locEdges;
                    set<int>                                       locFaces;

                    int region1ID, region2ID, i, j, k, cnt;
                    SpatialDomains::BoundaryConditionShPtr locBCond;

                    // Set up a set of all local verts and edges.
                    for(i = 0; i < (*m_exp).size(); ++i)
                    {
                        for(j = 0; j < (*m_exp)[i]->GetNverts(); ++j)
                        {
                            int id = (*m_exp)[i]->GetGeom()->GetVid(j);
                            locVerts.insert(id);
                        }

                        for(j = 0; j < (*m_exp)[i]->GetGeom()->GetNumEdges(); ++j)
                        {
                            int id = (*m_exp)[i]->GetGeom()->GetEid(j);
                            locEdges.insert(id);
                        }
                    }

                    // Begin by populating the perComps map. We loop over all periodic
                    // boundary conditions and determine the composite associated with
                    // it, then fill out the all* maps.
                    for (auto &it : bregions)
                    {

                        locBCond = GetBoundaryCondition
                            (bconditions, it.first, variable);

                        if (locBCond->GetBoundaryConditionType()
                            != SpatialDomains::ePeriodic)
                        {
                            continue;
                        }

                        // Identify periodic boundary region IDs.
                        region1ID = it.first;
                        region2ID = std::static_pointer_cast<
                            SpatialDomains::PeriodicBoundaryCondition>
                            (locBCond)->m_connectedBoundaryRegion;

                        // Check the region only contains a single composite.
                        ASSERTL0(it.second->size() == 1,
                                 "Boundary region "+boost::lexical_cast<string>
                                 (region1ID)+" should only contain 1 composite.");

                        // From this identify composites by looking at the original
                        // boundary region ordering. Note that in serial the mesh
                        // partitioner is not run, so this map will be empty and
                        // therefore needs to be populated by using the corresponding
                        // boundary region.
                        int cId1, cId2;
                        if (vComm->GetSize() == 1)
                        {
                            cId1 = it.second->begin()->first;
                            cId2 = bregions.find(region2ID)->second->begin()->first;
                        }
                        else
                        {
                            cId1 = bndRegOrder.find(region1ID)->second[0];
                            cId2 = bndRegOrder.find(region2ID)->second[0];
                        }

                        // check to see if boundary is rotationally aligned
                        if(boost::icontains(locBCond->GetUserDefined(),"Rotated"))
                        {
                            vector<string> tmpstr;

                            boost::split(tmpstr,locBCond->GetUserDefined(),
                                         boost::is_any_of(":"));

                            if(boost::iequals(tmpstr[0],"Rotated"))
                            {
                                ASSERTL1(tmpstr.size() > 2,
                                         "Expected Rotated user defined string to "
                                         "contain direction and rotation angle "
                                         "and optionally a tolerance, "
                                         "i.e. Rotated:dir:PI/2:1e-6");


                                ASSERTL1((tmpstr[1] == "x")||(tmpstr[1] == "y")
                                         ||(tmpstr[1] == "z"), "Rotated Dir is "
                                         "not specified as x,y or z");

                                RotPeriodicInfo RotInfo;
                                RotInfo.m_dir = (tmpstr[1] == "x")? 0:
                                    (tmpstr[1] == "y")? 1:2;

                                LibUtilities::Interpreter strEval;
                                int ExprId = strEval.DefineFunction("", tmpstr[2]);
                                RotInfo.m_angle = strEval.Evaluate(ExprId);

                                if(tmpstr.size() == 4)
                                {
                                    try {
                                        RotInfo.m_tol = boost::lexical_cast
                                            <NekDouble>(tmpstr[3]);
                                    }
                                    catch (...) {
                                        NEKERROR(ErrorUtil::efatal,
                                                 "failed to cast tolerance input "
                                                 "to a double value in Rotated"
                                                 "boundary information");
                                    }
                                }
                                else
                                {
                                    RotInfo.m_tol = 1e-8;
                                }
                                rotComp[cId1] = RotInfo;
                            }
                        }

                        SpatialDomains::CompositeSharedPtr c = it.second->begin()->second;

                        vector<unsigned int> tmpOrder;

                        // store the rotation info of this

                        // From the composite, we now construct the allVerts, allEdges
                        // and allCoord map so that they can be transferred across
                        // processors. We also populate the locFaces set to store a
                        // record of all faces local to this process.
                        for (i = 0; i < c->m_geomVec.size(); ++i)
                        {
                            SpatialDomains::Geometry2DSharedPtr faceGeom =
                                std::dynamic_pointer_cast<
                                    SpatialDomains::Geometry2D>(c->m_geomVec[i]);
                            ASSERTL1(faceGeom, "Unable to cast to shared ptr");

                            // Get geometry ID of this face and store in locFaces.
                            int faceId = c->m_geomVec[i]->GetGlobalID();
                            locFaces.insert(faceId);

                            // In serial, mesh partitioning will not have occurred so
                            // need to fill composite ordering map manually.
                            if (vComm->GetSize() == 1)
                            {
                                tmpOrder.push_back(c->m_geomVec[i]->GetGlobalID());
                            }

                            // Loop over vertices and edges of the face to populate
                            // allVerts, allEdges and allCoord maps.
                            vector<int> vertList, edgeList;
                            SpatialDomains::PointGeomVector coordVec;
                            vector<StdRegions::Orientation> orientVec;
                            for (j = 0; j < faceGeom->GetNumVerts(); ++j)
                            {
                                vertList .push_back(faceGeom->GetVid   (j));
                                edgeList .push_back(faceGeom->GetEid   (j));
                                coordVec .push_back(faceGeom->GetVertex(j));
                                orientVec.push_back(faceGeom->GetEorient(j));
                            }

                            allVerts [faceId] = vertList;
                            allEdges [faceId] = edgeList;
                            allCoord [faceId] = coordVec;
                            allOrient[faceId] = orientVec;
                        }

                        // In serial, record the composite ordering in compOrder for
                        // later in the routine.
                        if (vComm->GetSize() == 1)
                        {
                            compOrder[it.second->begin()->first] = tmpOrder;
                        }

                        // See if we already have either region1 or region2 stored in
                        // perComps map already and do a sanity check to ensure regions
                        // are mutually periodic.
                        if (perComps.count(cId1) == 0)
                        {
                            if (perComps.count(cId2) == 0)
                            {
                                perComps[cId1] = cId2;
                            }
                            else
                            {
                                std::stringstream ss;
                                ss << "Boundary region " << cId2 << " should be "
                                   << "periodic with " << perComps[cId2] << " but "
                                   << "found " << cId1 << " instead!";
                                ASSERTL0(perComps[cId2] == cId1, ss.str());
                            }
                        }
                        else
                        {
                            std::stringstream ss;
                            ss << "Boundary region " << cId1 << " should be "
                               << "periodic with " << perComps[cId1] << " but "
                               << "found " << cId2 << " instead!";
                            ASSERTL0(perComps[cId1] == cId1, ss.str());
                        }
                    }

                    // The next routines process local face lists to exchange vertices,
                    // edges and faces.
                    int              n = vComm->GetSize();
                    int              p = vComm->GetRank();
                    int              totFaces;
                    Array<OneD, int> facecounts(n,0);
                    Array<OneD, int> vertcounts(n,0);
                    Array<OneD, int> faceoffset(n,0);
                    Array<OneD, int> vertoffset(n,0);

                    Array<OneD, int> rotcounts(n,0);
                    Array<OneD, int> rotoffset(n,0);

                    rotcounts[p] = rotComp.size();
                    vComm->AllReduce(rotcounts, LibUtilities::ReduceSum);
                    int totrot  = Vmath::Vsum(n,rotcounts,1);

                    if(totrot)
                    {
                        for (i = 1; i < n ; ++i)
                        {
                            rotoffset[i] = rotoffset[i-1] + rotcounts[i-1];
                        }

                        Array<OneD, int> compid(totrot,0);
                        Array<OneD, int> rotdir(totrot,0);
                        Array<OneD, NekDouble> rotangle(totrot,0.0);
                        Array<OneD, NekDouble> rottol(totrot,0.0);

                        // fill in rotational informaiton
                        auto rIt = rotComp.begin();

                        for(i = 0; rIt != rotComp.end(); ++rIt)
                        {
                            compid  [rotoffset[p] + i  ] = rIt->first;
                            rotdir  [rotoffset[p] + i  ] = rIt->second.m_dir;
                            rotangle[rotoffset[p] + i  ] = rIt->second.m_angle;
                            rottol  [rotoffset[p] + i++] = rIt->second.m_tol;
                        }

                        vComm->AllReduce(compid, LibUtilities::ReduceSum);
                        vComm->AllReduce(rotdir, LibUtilities::ReduceSum);
                        vComm->AllReduce(rotangle, LibUtilities::ReduceSum);
                        vComm->AllReduce(rottol, LibUtilities::ReduceSum);

                        // Fill in full rotational composite list
                        for(i =0; i < totrot; ++i)
                        {
                            RotPeriodicInfo rinfo(rotdir[i],rotangle[i], rottol[i]);

                            rotComp[compid[i]] = rinfo;
                        }
                    }

                    // First exchange the number of faces on each process.
                    facecounts[p] = locFaces.size();
                    vComm->AllReduce(facecounts, LibUtilities::ReduceSum);

                    // Set up an offset map to allow us to distribute face IDs to all
                    // processors.
                    faceoffset[0] = 0;
                    for (i = 1; i < n; ++i)
                    {
                        faceoffset[i] = faceoffset[i-1] + facecounts[i-1];
                    }

                    // Calculate total number of faces.
                    totFaces = Vmath::Vsum(n, facecounts, 1);

                    // faceIds holds face IDs for each periodic face. faceVerts holds
                    // the number of vertices in this face.
                    Array<OneD, int> faceIds  (totFaces, 0);
                    Array<OneD, int> faceVerts(totFaces, 0);

                    // Process p writes IDs of its faces into position faceoffset[p] of
                    // faceIds which allows us to perform an AllReduce to distribute
                    // information amongst processors.
                    auto sIt = locFaces.begin();
                    for (i = 0; sIt != locFaces.end(); ++sIt)
                    {
                        faceIds  [faceoffset[p] + i  ] = *sIt;
                        faceVerts[faceoffset[p] + i++] = allVerts[*sIt].size();
                    }

                    vComm->AllReduce(faceIds,   LibUtilities::ReduceSum);
                    vComm->AllReduce(faceVerts, LibUtilities::ReduceSum);

                    // procVerts holds number of vertices (and also edges since each
                    // face is 2D) on each process.
                    Array<OneD, int> procVerts(n,0);
                    int nTotVerts;

                    // Note if there are no periodic faces at all calling Vsum will
                    // cause a segfault.
                    if (totFaces > 0)
                    {
                        // Calculate number of vertices on each processor.
                        nTotVerts = Vmath::Vsum(totFaces, faceVerts, 1);
                    }
                    else
                    {
                        nTotVerts = 0;
                    }

                    for (i = 0; i < n; ++i)
                    {
                        if (facecounts[i] > 0)
                        {
                            procVerts[i] = Vmath::Vsum
                                (facecounts[i], faceVerts + faceoffset[i], 1);
                        }
                        else
                        {
                            procVerts[i] = 0;
                        }
                    }

                    // vertoffset is defined in the same manner as edgeoffset
                    // beforehand.
                    vertoffset[0] = 0;
                    for (i = 1; i < n; ++i)
                    {
                        vertoffset[i] = vertoffset[i-1] + procVerts[i-1];
                    }

                    // At this point we exchange all vertex IDs, edge IDs and vertex
                    // coordinates for each face. The coordinates are necessary because
                    // we need to calculate relative face orientations between periodic
                    // faces to determined edge and vertex connectivity.
                    Array<OneD, int>       vertIds(nTotVerts,   0);
                    Array<OneD, int>       edgeIds(nTotVerts,   0);
                    Array<OneD, int>       edgeOrt(nTotVerts,   0);
                    Array<OneD, NekDouble> vertX  (nTotVerts, 0.0);
                    Array<OneD, NekDouble> vertY  (nTotVerts, 0.0);
                    Array<OneD, NekDouble> vertZ  (nTotVerts, 0.0);

                    for (cnt = 0, sIt = locFaces.begin();
                         sIt != locFaces.end(); ++sIt)
                    {
                        for (j = 0; j < allVerts[*sIt].size(); ++j)
                        {
                            int vertId = allVerts[*sIt][j];
                            vertIds[vertoffset[p] + cnt  ] = vertId;
                            vertX  [vertoffset[p] + cnt  ] = (*allCoord[*sIt][j])(0);
                            vertY  [vertoffset[p] + cnt  ] = (*allCoord[*sIt][j])(1);
                            vertZ  [vertoffset[p] + cnt  ] = (*allCoord[*sIt][j])(2);
                            edgeIds[vertoffset[p] + cnt  ] = allEdges [*sIt][j];
                            edgeOrt[vertoffset[p] + cnt++] = allOrient[*sIt][j];
                        }
                    }

                    vComm->AllReduce(vertIds, LibUtilities::ReduceSum);
                    vComm->AllReduce(vertX,   LibUtilities::ReduceSum);
                    vComm->AllReduce(vertY,   LibUtilities::ReduceSum);
                    vComm->AllReduce(vertZ,   LibUtilities::ReduceSum);
                    vComm->AllReduce(edgeIds, LibUtilities::ReduceSum);
                    vComm->AllReduce(edgeOrt, LibUtilities::ReduceSum);

                    // Finally now we have all of this information, we construct maps
                    // which make accessing the information easier. These are
                    // conceptually the same as all* maps at the beginning of the
                    // routine, but now hold information for all periodic vertices.
                    map<int, vector<int> >                          vertMap;
                    map<int, vector<int> >                          edgeMap;
                    map<int, SpatialDomains::PointGeomVector>       coordMap;

                    // These final two maps are required for determining the relative
                    // orientation of periodic edges. vCoMap associates vertex IDs with
                    // their coordinates, and eIdMap maps an edge ID to the two vertices
                    // which construct it.
                    map<int, SpatialDomains::PointGeomSharedPtr>    vCoMap;
                    map<int, pair<int, int> >                       eIdMap;

                    for (cnt = i = 0; i < totFaces; ++i)
                    {
                        vector<int> edges(faceVerts[i]);
                        vector<int> verts(faceVerts[i]);
                        SpatialDomains::PointGeomVector coord(faceVerts[i]);

                        // Keep track of cnt to enable correct edge vertices to be
                        // inserted into eIdMap.
                        int tmp = cnt;
                        for (j = 0; j < faceVerts[i]; ++j, ++cnt)
                        {
                            edges[j] = edgeIds[cnt];
                            verts[j] = vertIds[cnt];
                            coord[j]  = MemoryManager<SpatialDomains::PointGeom>
                                ::AllocateSharedPtr(3, verts[j], vertX[cnt],
                                                    vertY[cnt], vertZ[cnt]);
                            vCoMap[vertIds[cnt]] = coord[j];

                            // Try to insert edge into the eIdMap to avoid re-inserting.
                            auto testIns = eIdMap.insert
                                ( make_pair( edgeIds[cnt], make_pair(vertIds[tmp+j],
                                                                     vertIds[tmp+((j+1) % faceVerts[i])])));

                            if (testIns.second == false)
                            {
                                continue;
                            }

                            // If the edge is reversed with respect to the face, then
                            // swap the edges so that we have the original ordering of
                            // the edge in the 3D element. This is necessary to properly
                            // determine edge orientation. Note that the logic relies on
                            // the fact that all edge forward directions are CCW
                            // orientated: we use a tensor product ordering for 2D
                            // elements so need to reverse this for edge IDs 2 and 3.
                            StdRegions::Orientation edgeOrient =
                                static_cast<StdRegions::Orientation>(edgeOrt[cnt]);
                            if (j > 1)
                            {
                                edgeOrient = edgeOrient == StdRegions::eForwards ?
                                    StdRegions::eBackwards : StdRegions::eForwards;
                            }

                            if (edgeOrient == StdRegions::eBackwards)
                            {
                                swap(testIns.first->second.first,
                                     testIns.first->second.second);
                            }
                        }

                        vertMap [faceIds[i]] = verts;
                        edgeMap [faceIds[i]] = edges;
                        coordMap[faceIds[i]] = coord;
                    }

                    // Go through list of composites and figure out which edges are
                    // parallel from original ordering in session file. This includes
                    // composites which are not necessarily on this process.

                    // Store temporary map of periodic vertices which will hold all
                    // periodic vertices on the entire mesh so that doubly periodic
                    // vertices/edges can be counted properly across partitions. Local
                    // vertices/edges are copied into m_periodicVerts and
                    // m_periodicEdges at the end of the function.
                    PeriodicMap periodicVerts, periodicEdges;

                    // Construct two maps which determine how vertices and edges of
                    // faces connect given a specific face orientation. The key of the
                    // map is the number of vertices in the face, used to determine
                    // difference between tris and quads.
                    map<int, map<StdRegions::Orientation, vector<int> > > vmap;
                    map<int, map<StdRegions::Orientation, vector<int> > > emap;

                    map<StdRegions::Orientation, vector<int> > quadVertMap;
                    quadVertMap[StdRegions::eDir1FwdDir1_Dir2FwdDir2] = {0,1,2,3};
                    quadVertMap[StdRegions::eDir1FwdDir1_Dir2BwdDir2] = {3,2,1,0};
                    quadVertMap[StdRegions::eDir1BwdDir1_Dir2FwdDir2] = {1,0,3,2};
                    quadVertMap[StdRegions::eDir1BwdDir1_Dir2BwdDir2] = {2,3,0,1};
                    quadVertMap[StdRegions::eDir1FwdDir2_Dir2FwdDir1] = {0,3,2,1};
                    quadVertMap[StdRegions::eDir1FwdDir2_Dir2BwdDir1] = {1,2,3,0};
                    quadVertMap[StdRegions::eDir1BwdDir2_Dir2FwdDir1] = {3,0,1,2};
                    quadVertMap[StdRegions::eDir1BwdDir2_Dir2BwdDir1] = {2,1,0,3};

                    map<StdRegions::Orientation, vector<int> > quadEdgeMap;
                    quadEdgeMap[StdRegions::eDir1FwdDir1_Dir2FwdDir2] = {0,1,2,3};
                    quadEdgeMap[StdRegions::eDir1FwdDir1_Dir2BwdDir2] = {2,1,0,3};
                    quadEdgeMap[StdRegions::eDir1BwdDir1_Dir2FwdDir2] = {0,3,2,1};
                    quadEdgeMap[StdRegions::eDir1BwdDir1_Dir2BwdDir2] = {2,3,0,1};
                    quadEdgeMap[StdRegions::eDir1FwdDir2_Dir2FwdDir1] = {3,2,1,0};
                    quadEdgeMap[StdRegions::eDir1FwdDir2_Dir2BwdDir1] = {1,2,3,0};
                    quadEdgeMap[StdRegions::eDir1BwdDir2_Dir2FwdDir1] = {3,0,1,2};
                    quadEdgeMap[StdRegions::eDir1BwdDir2_Dir2BwdDir1] = {1,0,3,2};

                    map<StdRegions::Orientation, vector<int> > triVertMap;
                    triVertMap[StdRegions::eDir1FwdDir1_Dir2FwdDir2] = {0,1,2};
                    triVertMap[StdRegions::eDir1BwdDir1_Dir2FwdDir2] = {1,0,2};

                    map<StdRegions::Orientation, vector<int> > triEdgeMap;
                    triEdgeMap[StdRegions::eDir1FwdDir1_Dir2FwdDir2] = {0,1,2};
                    triEdgeMap[StdRegions::eDir1BwdDir1_Dir2FwdDir2] = {0,2,1};

                    vmap[3] = triVertMap;
                    vmap[4] = quadVertMap;
                    emap[3] = triEdgeMap;
                    emap[4] = quadEdgeMap;

                    map<int,int> allCompPairs;

                    // Collect composite id's of each periodic face for use if rotation is required
                    map<int,int> fIdToCompId;

                    // Finally we have enough information to populate the periodic
                    // vertex, edge and face maps. Begin by looping over all pairs of
                    // periodic composites to determine pairs of periodic faces.
                    for (auto &cIt : perComps)
                    {
                        SpatialDomains::CompositeSharedPtr c[2];
                        const int   id1  = cIt.first;
                        const int   id2  = cIt.second;
                        std::string id1s = boost::lexical_cast<string>(id1);
                        std::string id2s = boost::lexical_cast<string>(id2);

                        if (compMap.count(id1) > 0)
                        {
                            c[0] = compMap[id1];
                        }

                        if (compMap.count(id2) > 0)
                        {
                            c[1] = compMap[id2];
                        }

                        ASSERTL0(c[0] || c[1],
                                 "Neither composite not found on this process!");

                        // Loop over composite ordering to construct list of all
                        // periodic faces, regardless of whether they are on this
                        // process.
                        map<int,int> compPairs;


                        ASSERTL0(compOrder.count(id1) > 0,
                                 "Unable to find composite "+id1s+" in order map.");
                        ASSERTL0(compOrder.count(id2) > 0,
                                 "Unable to find composite "+id2s+" in order map.");
                        ASSERTL0(compOrder[id1].size() == compOrder[id2].size(),
                                 "Periodic composites "+id1s+" and "+id2s+
                                 " should have the same number of elements.");
                        ASSERTL0(compOrder[id1].size() > 0,
                                 "Periodic composites "+id1s+" and "+id2s+
                                 " are empty!");

                        // Look up composite ordering to determine pairs.
                        for (i = 0; i < compOrder[id1].size(); ++i)
                        {
                            int eId1 = compOrder[id1][i];
                            int eId2 = compOrder[id2][i];

                            ASSERTL0(compPairs.count(eId1) == 0,
                                     "Already paired.");

                            // Sanity check that the faces are mutually periodic.
                            if (compPairs.count(eId2) != 0)
                            {
                                ASSERTL0(compPairs[eId2] == eId1, "Pairing incorrect");
                            }
                            compPairs[eId1] = eId2;

                            // store  a map of face ids to composite ids
                            fIdToCompId[eId1] = id1;
                            fIdToCompId[eId2] = id2;
                        }

                        // Now that we have all pairs of periodic faces, loop over the
                        // ones local on this process and populate face/edge/vertex
                        // maps.
                        for (auto &pIt : compPairs)
                        {
                            int  ids  [2] = {pIt.first, pIt.second};
                            bool local[2] = {locFaces.count(pIt.first) > 0,
                                locFaces.count(pIt.second) > 0};

                            ASSERTL0(coordMap.count(ids[0]) > 0 &&
                                     coordMap.count(ids[1]) > 0,
                                     "Unable to find face in coordinate map");

                            allCompPairs[pIt.first ] = pIt.second;
                            allCompPairs[pIt.second] = pIt.first;

                            // Loop up coordinates of the faces, check they have the
                            // same number of vertices.
                            SpatialDomains::PointGeomVector tmpVec[2]
                        = { coordMap[ids[0]], coordMap[ids[1]] };

                            ASSERTL0(tmpVec[0].size() == tmpVec[1].size(),
                                     "Two periodic faces have different number "
                                     "of vertices!");

                            // o will store relative orientation of faces. Note that in
                            // some transpose cases (Dir1FwdDir2_Dir2BwdDir1 and
                            // Dir1BwdDir1_Dir2FwdDir1) it seems orientation will be
                            // different going from face1->face2 instead of face2->face1
                            // (check this).
                            StdRegions::Orientation o;
                            bool rotbnd = false;
                            int  dir = 0;
                            NekDouble angle = 0.0;
                            NekDouble sign = 0.0;
                            NekDouble tol = 1e-8;

                            // check to see if perioid boundary is rotated
                            if(rotComp.count(fIdToCompId[pIt.first]))
                            {
                                rotbnd = true;
                                dir   = rotComp[fIdToCompId[pIt.first]].m_dir;
                                angle = rotComp[fIdToCompId[pIt.first]].m_angle;
                                tol   = rotComp[fIdToCompId[pIt.first]].m_tol;
                            }

                            // Record periodic faces.
                            for (i = 0; i < 2; ++i)
                            {
                                if (!local[i])
                                {
                                    continue;
                                }

                                // Reference to the other face.
                                int other = (i+1) % 2;

                                // angle is set up for i = 0 to i = 1
                                sign = (i == 0)? 1.0:-1.0;

                                // Calculate relative face orientation.
                                if (tmpVec[0].size() == 3)
                                {
                                    o = SpatialDomains::TriGeom::GetFaceOrientation
                                        (tmpVec[i], tmpVec[other],
                                         rotbnd, dir, sign*angle, tol);
                                }
                                else
                                {
                                    o = SpatialDomains::QuadGeom::GetFaceOrientation(
                                                                                     tmpVec[i], tmpVec[other],
                                                                                     rotbnd,dir,sign*angle,tol);
                                }

                                // Record face ID, orientation and whether other face is
                                // local.
                                PeriodicEntity ent(ids  [other], o,
                                                   local[other]);
                                m_periodicFaces[ids[i]].push_back(ent);
                            }

                            int nFaceVerts = vertMap[ids[0]].size();

                            // Determine periodic vertices.
                            for (i = 0; i < 2; ++i)
                            {
                                int other = (i+1) % 2;

                                // angle is set up for i = 0 to i = 1
                                sign = (i == 0)? 1.0:-1.0;

                                // Calculate relative face orientation.
                                if (tmpVec[0].size() == 3)
                                {
                                    o = SpatialDomains::TriGeom::GetFaceOrientation(
                                                                                    tmpVec[i], tmpVec[other], rotbnd, dir,
                                                                                    sign*angle, tol);
                                }
                                else
                                {
                                    o = SpatialDomains::QuadGeom::GetFaceOrientation(
                                                                                     tmpVec[i], tmpVec[other], rotbnd, dir,
                                                                                     sign*angle, tol);
                                }

                                if (nFaceVerts == 3)
                                {
                                    ASSERTL0(o == StdRegions::eDir1FwdDir1_Dir2FwdDir2 ||
                                             o == StdRegions::eDir1BwdDir1_Dir2FwdDir2,
                                             "Unsupported face orientation for face "+
                                             boost::lexical_cast<string>(ids[i]));
                                }

                                // Look up vertices for this face.
                                vector<int> per1 = vertMap[ids[i]];
                                vector<int> per2 = vertMap[ids[other]];

                                // tmpMap will hold the pairs of vertices which are
                                // periodic.
                                map<int, pair<int, bool> > tmpMap;

                                // Use vmap to determine which vertices connect given
                                // the orientation o.
                                for (j = 0; j < nFaceVerts; ++j)
                                {
                                    int v = vmap[nFaceVerts][o][j];
                                    tmpMap[per1[j]] = make_pair
                                        (per2[v], locVerts.count(per2[v]) > 0);
                                }

                                // Now loop over tmpMap to associate periodic vertices.
                                for (auto &mIt : tmpMap)
                                {
                                    PeriodicEntity ent2(mIt.second.first,
                                                        StdRegions::eNoOrientation,
                                                        mIt.second.second);

                                    // See if this vertex has been recorded already.
                                    auto perIt = periodicVerts.find(mIt.first);

                                    if (perIt == periodicVerts.end())
                                    {
                                        // Vertex is new - just record this entity as
                                        // usual.
                                        periodicVerts[mIt.first].push_back(ent2);
                                        perIt = periodicVerts.find(mIt.first);
                                    }
                                    else
                                    {
                                        // Vertex is known - loop over the vertices
                                        // inside the record and potentially add vertex
                                        // mIt.second to the list.
                                        for (k = 0; k < perIt->second.size(); ++k)
                                        {
                                            if (perIt->second[k].id == mIt.second.first)
                                            {
                                                break;
                                            }
                                        }

                                        if (k == perIt->second.size())
                                        {
                                            perIt->second.push_back(ent2);
                                        }
                                    }
                                }
                            }

                            // Determine periodic edges. Logic is the same as above,
                            // and perhaps should be condensed to avoid replication.
                            for (i = 0; i < 2; ++i)
                            {
                                int other = (i+1) % 2;

                                // angle is set up for i = 0 to i = 1
                                sign = (i == 0)? 1.0:-1.0;

                                if (tmpVec[0].size() == 3)
                                {
                                    o = SpatialDomains::TriGeom::GetFaceOrientation(
                                                                                    tmpVec[i], tmpVec[other], rotbnd, dir,
                                                                                    sign*angle, tol);
                                }
                                else
                                {
                                    o = SpatialDomains::QuadGeom::GetFaceOrientation(
                                                                                     tmpVec[i], tmpVec[other], rotbnd, dir,
                                                                                     sign*angle, tol);
                                }

                                vector<int> per1 = edgeMap[ids[i]];
                                vector<int> per2 = edgeMap[ids[other]];

                                map<int, pair<int, bool> > tmpMap;

                                for (j = 0; j < nFaceVerts; ++j)
                                {
                                    int e = emap[nFaceVerts][o][j];
                                    tmpMap[per1[j]] = make_pair
                                        (per2[e], locEdges.count(per2[e]) > 0);
                                }

                                for (auto &mIt : tmpMap)
                                {
                                    // Note we assume orientation of edges is forwards -
                                    // this may be reversed later.
                                    PeriodicEntity ent2(mIt.second.first,
                                                        StdRegions::eForwards,
                                                        mIt.second.second);
                                    auto perIt = periodicEdges.find(mIt.first);

                                    if (perIt == periodicEdges.end())
                                    {
                                        periodicEdges[mIt.first].push_back(ent2);
                                        perIt = periodicEdges.find(mIt.first);
                                    }
                                    else
                                    {
                                        for (k = 0; k < perIt->second.size(); ++k)
                                        {
                                            if (perIt->second[k].id == mIt.second.first)
                                            {
                                                break;
                                            }
                                        }

                                        if (k == perIt->second.size())
                                        {
                                            perIt->second.push_back(ent2);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // Make sure that the nubmer of face pairs and the
                    // face Id to composite Id map match in size
                    ASSERTL1(allCompPairs.size() == fIdToCompId.size(),
                            "At this point the size of allCompPairs "
                            "should have been the same as fIdToCompId");

                    // Distribute the size of the periodic boundary to all 
                    // processes. We assume that all processes who own periodic
                    // faces have the same local copy of allCompPairs
                    int NPairs  = allCompPairs.size();
                    vComm->AllReduce(NPairs, LibUtilities::ReduceMax);

                    // Check that the previous assertion regarding allCompPairs
                    // is correct
                    ASSERTL0(
                        allCompPairs.size() == NPairs || allCompPairs.size() == 0,
                        "Local copy of allCompPairs not the same for all ranks"
                    );

                    // Allocate local vectors that will contain the content of
                    // allCompPairs if process owns faces on periodic boundary
                    Array<OneD, int> first(NPairs, -1);
                    Array<OneD, int> second(NPairs, -1);
                    cnt = 0;
                    for(const auto &it : allCompPairs)
                    {
                        first[cnt]    = it.first;
                        second[cnt++] = it.second;
                    }
                    
                    // Distribute the content in first and second to all processes
                    vComm->AllReduce(first,  LibUtilities::ReduceMax);
                    vComm->AllReduce(second, LibUtilities::ReduceMax);

                    // Check that the MPI Allreduce routine worked
                    ASSERTL0(std::count(first.begin(), first.end(), -1) == 0, 
                        "Distribution of allCompPairs failed");
                    ASSERTL0(std::count(second.begin(), second.end(), -1) == 0, 
                        "Distribution of allCompPairs failed")

                    // Put content back in allCompPairs
                    allCompPairs.clear();
                    for(cnt = 0; cnt < NPairs; ++cnt)
                    {
                        allCompPairs[first[cnt]] = second[cnt];
                    }

                    // Store face ID to composite ID map for rotational boundaries
                    if(rotComp.size())
                    {   
                        // Set values to -1
                        std::fill(first.begin(), first.end(), -1);
                        std::fill(second.begin(), second.end(), -1);

                        cnt = 0;
                        for (const auto &pIt : fIdToCompId)
                        {
                            first [cnt  ] = pIt.first;
                            second[cnt++] = pIt.second;
                        }

                        vComm->AllReduce(first,  LibUtilities::ReduceMax);
                        vComm->AllReduce(second, LibUtilities::ReduceMax);

                        // Check that the MPI Allreduce routine worked
                        ASSERTL0(std::count(first.begin(), first.end(), -1) == 0,
                            "Distribution of fIdToCompId failed");
                        ASSERTL0(std::count(second.begin(), second.end(), -1) == 0,
                            "Distribution of fIdToCompId failed")

                        fIdToCompId.clear();
                        for(cnt = 0; cnt < NPairs; ++cnt)
                        {
                            fIdToCompId[first[cnt]] = second[cnt];
                        }
                    }

                    // also will need an edge id to composite id at end of routine
                    map<int,int> eIdToCompId;

                    // Search for periodic vertices and edges which are not
                    // in a periodic composite but lie in this process. First,
                    // loop over all information we have from other
                    // processors.
                    for (cnt = i = 0; i < totFaces; ++i)
                    {
                        bool rotbnd = false;
                        int dir = 0;
                        NekDouble angle = 0.0;
                        NekDouble tol = 1e-8;

                        int faceId    = faceIds[i];

                        ASSERTL0(allCompPairs.count(faceId) > 0,
                                 "Unable to find matching periodic face.");

                        int perFaceId = allCompPairs[faceId];

                        // check to see if periodic boundary is rotated
                        ASSERTL1((rotComp.size() == 0) ||
                                 fIdToCompId.count(faceId) > 0,"Face " +
                                 boost::lexical_cast<string>(faceId) +
                                 " not found in fIdtoCompId map");
                        if(rotComp.count(fIdToCompId[faceId]))
                        {
                            rotbnd = true;
                            dir   = rotComp[fIdToCompId[faceId]].m_dir;
                            angle = rotComp[fIdToCompId[faceId]].m_angle;
                            tol   = rotComp[fIdToCompId[faceId]].m_tol;
                        }

                        for (j = 0; j < faceVerts[i]; ++j, ++cnt)
                        {
                            int vId = vertIds[cnt];

                            auto perId = periodicVerts.find(vId);

                            if (perId == periodicVerts.end())
                            {

                                // This vertex is not included in the
                                // map. Figure out which vertex it is supposed
                                // to be periodic with. perFaceId is the face
                                // ID which is periodic with faceId. The logic
                                // is much the same as the loop above.
                                SpatialDomains::PointGeomVector tmpVec[2]
                            = { coordMap[faceId], coordMap[perFaceId] };

                                int nFaceVerts = tmpVec[0].size();
                                StdRegions::Orientation o = nFaceVerts == 3 ?
                                    SpatialDomains::TriGeom::GetFaceOrientation(
                                                                                tmpVec[0], tmpVec[1], rotbnd, dir, angle, tol):
                                    SpatialDomains::QuadGeom::GetFaceOrientation(
                                                                                 tmpVec[0], tmpVec[1], rotbnd, dir, angle, tol);

                                // Use vmap to determine which vertex of the other face
                                // should be periodic with this one.
                                int perVertexId = vertMap[perFaceId][vmap[nFaceVerts][o][j]];


                                PeriodicEntity ent(perVertexId,
                                                   StdRegions::eNoOrientation,
                                                   locVerts.count(perVertexId) > 0);

                                periodicVerts[vId].push_back(ent);
                            }

                            int eId = edgeIds[cnt];

                            perId = periodicEdges.find(eId);

                            // this map is required at very end to determine rotation of edges.
                            if(rotbnd)
                            {
                                eIdToCompId[eId] = fIdToCompId[faceId];
                            }

                            if (perId == periodicEdges.end())
                            {
                                // This edge is not included in the map. Figure
                                // out which edge it is supposed to be periodic
                                // with. perFaceId is the face ID which is
                                // periodic with faceId. The logic is much the
                                // same as the loop above.
                                SpatialDomains::PointGeomVector tmpVec[2]
                            = { coordMap[faceId], coordMap[perFaceId] };

                                int nFaceEdges = tmpVec[0].size();
                                StdRegions::Orientation o = nFaceEdges == 3 ?
                                    SpatialDomains::TriGeom::GetFaceOrientation(
                                                                                tmpVec[0], tmpVec[1], rotbnd, dir, angle, tol):
                                    SpatialDomains::QuadGeom::GetFaceOrientation(
                                                                                 tmpVec[0], tmpVec[1], rotbnd, dir, angle, tol);

                                // Use emap to determine which edge of the other
                                // face should be periodic with this one.
                                int perEdgeId = edgeMap[perFaceId][emap[nFaceEdges][o][j]];

                                PeriodicEntity ent(perEdgeId,
                                                   StdRegions::eForwards,
                                                   locEdges.count(perEdgeId) > 0);

                                periodicEdges[eId].push_back(ent);


                                // this map is required at very end to
                                // determine rotation of edges.
                                if(rotbnd)
                                {
                                    eIdToCompId[perEdgeId] = fIdToCompId[perFaceId];
                                }
                            }
                        }
                    }

                    // Finally, we must loop over the periodicVerts and periodicEdges
                    // map to complete connectivity information.
                    for (auto &perIt : periodicVerts)
                    {
                        // For each vertex that is periodic with this one...
                        for (i = 0; i < perIt.second.size(); ++i)
                        {
                            // Find the vertex in the periodicVerts map...
                            auto perIt2 = periodicVerts.find(perIt.second[i].id);
                            ASSERTL0(perIt2 != periodicVerts.end(),
                                     "Couldn't find periodic vertex.");

                            // Now search through this vertex's list and make sure that
                            // we have a record of any vertices which aren't in the
                            // original list.
                            for (j = 0; j < perIt2->second.size(); ++j)
                            {
                                if (perIt2->second[j].id == perIt.first)
                                {
                                    continue;
                                }

                                for (k = 0; k < perIt.second.size(); ++k)
                                {
                                    if (perIt2->second[j].id == perIt.second[k].id)
                                    {
                                        break;
                                    }
                                }

                                if (k == perIt.second.size())
                                {
                                    perIt.second.push_back(perIt2->second[j]);
                                }
                            }
                        }
                    }

                    for (auto &perIt : periodicEdges)
                    {
                        for (i = 0; i < perIt.second.size(); ++i)
                        {
                            auto perIt2 = periodicEdges.find(perIt.second[i].id);
                            ASSERTL0(perIt2 != periodicEdges.end(),
                                     "Couldn't find periodic edge.");

                            for (j = 0; j < perIt2->second.size(); ++j)
                            {
                                if (perIt2->second[j].id == perIt.first)
                                {
                                    continue;
                                }

                                for (k = 0; k < perIt.second.size(); ++k)
                                {
                                    if (perIt2->second[j].id == perIt.second[k].id)
                                    {
                                        break;
                                    }
                                }

                                if (k == perIt.second.size())
                                {
                                    perIt.second.push_back(perIt2->second[j]);
                                }
                            }
                        }
                    }

                    // Loop over periodic edges to determine relative edge orientations.
                    for (auto &perIt : periodicEdges)
                    {
                        bool rotbnd = false;
                        int dir = 0;
                        NekDouble angle = 0.0;
                        NekDouble tol = 1e-8;


                        // Find edge coordinates
                        auto eIt = eIdMap.find(perIt.first);
                        SpatialDomains::PointGeom v[2] = {
                            *vCoMap[eIt->second.first],
                            *vCoMap[eIt->second.second]
                        };

                        // check to see if perioid boundary is rotated
                        if(rotComp.count(eIdToCompId[perIt.first]))
                        {
                            rotbnd = true;
                            dir   = rotComp[eIdToCompId[perIt.first]].m_dir;
                            angle = rotComp[eIdToCompId[perIt.first]].m_angle;
                            tol   = rotComp[eIdToCompId[perIt.first]].m_tol;
                        }

                        // Loop over each edge, and construct a vector that takes us
                        // from one vertex to another. Use this to figure out which
                        // vertex maps to which.
                        for (i = 0; i < perIt.second.size(); ++i)
                        {
                            eIt = eIdMap.find(perIt.second[i].id);

                            SpatialDomains::PointGeom w[2] = {
                                *vCoMap[eIt->second.first],
                                *vCoMap[eIt->second.second]
                            };

                            int vMap[2] = {-1,-1};
                            if(rotbnd)
                            {

                                SpatialDomains::PointGeom r;

                                r.Rotate(v[0],dir,angle);

                                if(r.dist(w[0])< tol)
                                {
                                    vMap[0] = 0;
                                }
                                else
                                {
                                    r.Rotate(v[1],dir,angle);
                                    if(r.dist(w[0]) < tol)
                                    {
                                        vMap[0] = 1;
                                    }
                                    else
                                    {
                                        NEKERROR(ErrorUtil::efatal,
                                                 "Unable to align rotationally "
                                                 "periodic edge vertex");
                                    }
                                }
                            }
                            else // translation test
                            {
                                NekDouble cx = 0.5*(w[0](0)-v[0](0)+w[1](0)-v[1](0));
                                NekDouble cy = 0.5*(w[0](1)-v[0](1)+w[1](1)-v[1](1));
                                NekDouble cz = 0.5*(w[0](2)-v[0](2)+w[1](2)-v[1](2));

                                for (j = 0; j < 2; ++j)
                                {
                                    NekDouble x = v[j](0);
                                    NekDouble y = v[j](1);
                                    NekDouble z = v[j](2);
                                    for (k = 0; k < 2; ++k)
                                    {
                                        NekDouble x1 = w[k](0)-cx;
                                        NekDouble y1 = w[k](1)-cy;
                                        NekDouble z1 = w[k](2)-cz;

                                        if (sqrt((x1-x)*(x1-x)+(y1-y)*(y1-y)+(z1-z)*(z1-z))
                                            < 1e-8)
                                        {
                                            vMap[k] = j;
                                            break;
                                        }
                                    }
                                }

                                // Sanity check the map.
                                ASSERTL0(vMap[0] >= 0 && vMap[1] >= 0,
                                         "Unable to align periodic edge vertex.");
                                ASSERTL0((vMap[0] == 0 || vMap[0] == 1) &&
                                         (vMap[1] == 0 || vMap[1] == 1) &&
                                         (vMap[0] != vMap[1]),
                                         "Unable to align periodic edge vertex.");
                            }

                            // If 0 -> 0 then edges are aligned already; otherwise
                            // reverse the orientation.
                            if (vMap[0] != 0)
                            {
                                perIt.second[i].orient = StdRegions::eBackwards;
                            }
                        }
                    }

                    // Do one final loop over periodic vertices/edges to remove
                    // non-local vertices/edges from map.
                    for (auto &perIt : periodicVerts)
                    {
                        if (locVerts.count(perIt.first) > 0)
                        {
                            m_periodicVerts.insert(perIt);
                        }
                    }

                    for (auto &perIt : periodicEdges)
                    {
                        if (locEdges.count(perIt.first) > 0)
                        {
                            m_periodicEdges.insert(perIt);
                        }
                    }
                }
                break;
            default:
                ASSERTL1(false,"Not setup for this expansion");
                break;
            }
        }

        /**
         *
         */
        GlobalLinSysSharedPtr DisContField::GetGlobalBndLinSys(
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
            auto matrixIter = m_globalBndMat->find(mkey);

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
         * For each boundary region, checks that the types and number of
         * boundary expansions in that region match.
         * @param   In          Field to compare with.
         * @return True if boundary conditions match.
         */
        bool DisContField::SameTypeOfBoundaryConditions(
                                                        const DisContField &In)
        {
            int i;
            bool returnval = true;

            for(i = 0; i < m_bndConditions.size(); ++i)
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


        vector<bool> &DisContField::GetNegatedFluxNormal(void)
        {
            if(m_negatedFluxNormal.size() == 0)
            {
                Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >
                    &elmtToTrace = m_traceMap->GetElmtToTrace();

                m_negatedFluxNormal.resize(2*GetExpSize());

                for(int i = 0; i < GetExpSize(); ++i)
                {

                    for(int v = 0; v < 2; ++v)
                    {
                        LocalRegions::Expansion0DSharedPtr vertExp =
                            elmtToTrace[i][v]->as<LocalRegions::Expansion0D>();

                        if(vertExp->GetLeftAdjacentElementExp()->GetGeom()
                           ->GetGlobalID() !=
                           (*m_exp)[i]->GetGeom()->GetGlobalID())
                        {
                            m_negatedFluxNormal[2*i+v] = true;
                        }
                        else
                        {
                            m_negatedFluxNormal[2*i+v] = false;
                        }
                    }
                }

            }

            return m_negatedFluxNormal;
        }



        /**
         * \brief This method extracts the "forward" and "backward" trace
         * data from the array \a field and puts the data into output
         * vectors \a Fwd and \a Bwd.
         *
         * We first define the convention which defines "forwards" and
         * "backwards". First an association is made between the vertex/edge/face of
         * each element and its corresponding vertex/edge/face in the trace space
         * using the mapping #m_traceMap. The element can either be
         * left-adjacent or right-adjacent to this trace face (see
         * Expansion2D::GetLeftAdjacentElementExp). Boundary faces are
         * always left-adjacent since left-adjacency is populated first.
         *
         * If the element is left-adjacent we extract the trace data
         * from \a field into the forward trace space \a Fwd; otherwise,
         * we place it in the backwards trace space \a Bwd. In this way,
         * we form a unique set of trace normals since these are always
         * extracted from left-adjacent elements.
         *
         * \param field is a NekDouble array which contains the fielddata
         * from which we wish to extract the backward and forward
         * orientated trace/face arrays.
         *
         * \return Updates a NekDouble array \a Fwd and \a Bwd
         */
        void DisContField::v_GetFwdBwdTracePhys
           (const Array<OneD, const NekDouble> &field,
            Array<OneD, NekDouble> &Fwd,
            Array<OneD, NekDouble> &Bwd,
            bool FillBnd,           // these should be template params so that compiler can remove them
            bool PutFwdInBwdOnBCs,
            bool DoExchange)
        {
            // Is this zeroing necessary?
            // Zero forward/backward vectors.
            Vmath::Zero(Fwd.size(), Fwd, 1);
            Vmath::Zero(Bwd.size(), Bwd, 1);

            // Basis definition on each element
            LibUtilities::BasisSharedPtr basis = (*m_exp)[0]->GetBasis(0);
            if((basis->GetBasisType() != LibUtilities::eGauss_Lagrange))
            {
                // blocked routine
                Array<OneD, NekDouble> tracevals(m_locTraceToTraceMap->
                                                 GetNLocTracePts());

                m_locTraceToTraceMap->LocTracesFromField(field, tracevals);
                m_locTraceToTraceMap->InterpLocTracesToTrace(0, tracevals, Fwd);

                Array<OneD, NekDouble> invals = tracevals + m_locTraceToTraceMap->
                    GetNFwdLocTracePts();
                m_locTraceToTraceMap->InterpLocTracesToTrace(1, invals, Bwd);
            }
            else
            {
                // Loop over elements and collect forward expansion
                auto nexp = GetExpSize();
                Array<OneD,NekDouble> e_tmp;
                LocalRegions::ExpansionSharedPtr exp;

                Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >
                    &elmtToTrace = m_traceMap->GetElmtToTrace();

                for (int n = 0, cnt = 0; n < nexp; ++n)
                {
                    exp = (*m_exp)[n];
                    auto phys_offset = GetPhys_Offset(n);

                    for(int e = 0; e < exp->GetNtraces(); ++e, ++cnt)
                    {
                        auto offset = m_trace->GetPhys_Offset(
                            elmtToTrace[n][e]->GetElmtId());

                        e_tmp = (m_leftAdjacentTraces[cnt])? Fwd + offset:
                            Bwd + offset;

                        exp->GetTracePhysVals(e, elmtToTrace[n][e],
                                              field + phys_offset, e_tmp);
                    }
                }
            }

            DisContField::v_PeriodicBwdCopy(Fwd,Bwd);

            if (FillBnd)
            {
                    v_FillBwdWithBoundCond(Fwd, Bwd, PutFwdInBwdOnBCs);
            }

            if(DoExchange)
            {
                // Do parallel exchange for forwards/backwards spaces.
                m_traceMap->GetAssemblyCommDG()->PerformExchange(Fwd, Bwd);
            }
        }

        void DisContField::v_FillBwdWithBoundCond(
                                                  const Array<OneD, NekDouble> &Fwd,
                                                  Array<OneD, NekDouble> &Bwd,
                                                  bool PutFwdInBwdOnBCs)
        {
            // Fill boundary conditions into missing elements
            if (PutFwdInBwdOnBCs) // just set Bwd value to be Fwd value on BCs
            {
                // Fill boundary conditions into missing elements
                for (int n = 0, cnt = 0; n < m_bndCondExpansions.size(); ++n)
                {
                    if (m_bndConditions[n]->GetBoundaryConditionType() ==
                        SpatialDomains::eDirichlet)
                    {
                        auto ne = m_bndCondExpansions[n]->GetExpSize();
                        for (int e = 0; e < ne; ++e)
                        {
                            auto npts = m_bndCondExpansions[n]->
                                GetExp(e)->GetTotPoints();
                            auto id2 = m_trace->GetPhys_Offset(m_traceMap->
                                        GetBndCondIDToGlobalTraceID(cnt+e));
                            Vmath::Vcopy(npts, &Fwd[id2], 1, &Bwd[id2], 1);
                        }

                        cnt += ne;
                    }
                    else if (m_bndConditions[n]->GetBoundaryConditionType() ==
                                SpatialDomains::eNeumann ||
                                m_bndConditions[n]->GetBoundaryConditionType() ==
                                SpatialDomains::eRobin)
                    {
                        auto ne = m_bndCondExpansions[n]->GetExpSize();
                        for (int e = 0; e < ne; ++e)
                        {
                            auto npts = m_bndCondExpansions[n]->
                                GetExp(e)->GetTotPoints();
                            auto id2 = m_trace->GetPhys_Offset(m_traceMap->
                                        GetBndCondIDToGlobalTraceID(cnt+e));
                            Vmath::Vcopy(npts, &Fwd[id2], 1, &Bwd[id2], 1);
                        }
                        cnt += ne;
                    }
                    else if (m_bndConditions[n]->GetBoundaryConditionType() !=
                                SpatialDomains::ePeriodic)
                    {
                        ASSERTL0(false,
                                    "Method not set up for this "
                                    "boundary condition.");
                    }
                }
            }
            else
            {
                for (int n = 0, cnt = 0; n < m_bndCondExpansions.size(); ++n)
                {
                    if (m_bndConditions[n]->GetBoundaryConditionType() ==
                        SpatialDomains::eDirichlet)
                    {
                        auto ne = m_bndCondExpansions[n]->GetExpSize();
                        for (int e = 0; e < ne; ++e)
                        {
                            auto npts = m_bndCondExpansions[n]->GetExp(e)
                                ->GetTotPoints();
                            auto id1 = m_bndCondExpansions[n]
                                ->GetPhys_Offset(e);
                            auto id2 = m_trace->GetPhys_Offset(m_traceMap->
                                               GetBndCondIDToGlobalTraceID(cnt+e));
                            Vmath::Vcopy(npts,
                                         &(m_bndCondExpansions[n]->GetPhys())
                                         [id1], 1, &Bwd[id2], 1);
                        }
                        cnt += ne;
                    }
                    else if (m_bndConditions[n]->GetBoundaryConditionType()
                             == SpatialDomains::eNeumann ||
                             m_bndConditions[n]->GetBoundaryConditionType()
                             == SpatialDomains::eRobin)
                    {
                        auto ne = m_bndCondExpansions[n]->GetExpSize();
                        for(int e = 0; e < ne; ++e)
                        {
                            auto npts = m_bndCondExpansions[n]->GetExp(e)
                                ->GetTotPoints();
                            auto id1  = m_bndCondExpansions[n]->
                                GetPhys_Offset(e);
                            ASSERTL0((m_bndCondExpansions[n]->
                                      GetPhys())[id1] == 0.0,
                                     "Method not set up for non-zero "
                                     "Neumann boundary condition");
                            auto id2  = m_trace->GetPhys_Offset(
                                      m_traceMap->GetBndCondIDToGlobalTraceID(cnt+e));
                            Vmath::Vcopy(npts, &Fwd[id2], 1, &Bwd[id2], 1);
                        }
                        
                        cnt += ne;
                    }
                    else if (m_bndConditions[n]->GetBoundaryConditionType()
                             == SpatialDomains::eNotDefined)
                    {
                        // Do nothing
                    }
                    else if (m_bndConditions[n]->GetBoundaryConditionType() !=
                             SpatialDomains::ePeriodic)
                    {
                        NEKERROR(ErrorUtil::efatal,
                                 "Method not set up for this boundary "
                                 "condition.");
                    }
                }
            }
        }

        void DisContField::v_AddTraceQuadPhysToField(
            const Array<OneD, const NekDouble> &Fwd,
            const Array<OneD, const NekDouble> &Bwd,
                  Array<OneD,       NekDouble> &field)
        {
            // Basis definition on each element
            LibUtilities::BasisSharedPtr basis = (*m_exp)[0]->GetBasis(0);
            if (basis->GetBasisType() != LibUtilities::eGauss_Lagrange)
            {
                Array<OneD, NekDouble> edgevals(m_locTraceToTraceMap->
                                               GetNLocTracePts(), 0.0);

                Array<OneD, NekDouble> invals = edgevals +
                    m_locTraceToTraceMap->GetNFwdLocTracePts();
                m_locTraceToTraceMap->RightIPTWLocEdgesToTraceInterpMat(
                                        1, Bwd, invals);

                m_locTraceToTraceMap->RightIPTWLocEdgesToTraceInterpMat(
                                        0, Fwd, edgevals);

                m_locTraceToTraceMap->AddLocTracesToField(edgevals, field);
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                    "v_AddTraceQuadPhysToField not coded for eGauss_Lagrange");
            }
        }

        void DisContField::v_ExtractTracePhys(Array<OneD, NekDouble> &outarray)
        {
            ASSERTL1(m_physState == true,"local physical space is not true ");
            v_ExtractTracePhys(m_phys, outarray);
        }

        /**
         * @brief This method extracts the trace (verts in 1D) from
         * the field @a inarray and puts the values in @a outarray.
         *
         * It assumes the field is C0 continuous so that it can
         * overwrite the edge data when visited by the two adjacent
         * elements.
         *
         * @param inarray   An array containing the 1D data from which we wish
         *                  to extract the edge data.
         * @param outarray  The resulting edge information.
         *
         * This will not work for non-boundary expansions
         */
        void DisContField::v_ExtractTracePhys
           (const Array<OneD, const NekDouble> &inarray,
            Array<OneD,       NekDouble> &outarray)
        {
            LibUtilities::BasisSharedPtr basis = (*m_exp)[0]->GetBasis(0);
            if( (basis->GetBasisType() != LibUtilities::eGauss_Lagrange))
            {
                Vmath::Zero(outarray.size(), outarray, 1);

                Array<OneD, NekDouble> tracevals(
                                    m_locTraceToTraceMap->GetNFwdLocTracePts());
                m_locTraceToTraceMap->FwdLocTracesFromField(inarray,tracevals);
                m_locTraceToTraceMap->InterpLocTracesToTrace(0,tracevals,outarray);
                m_traceMap->GetAssemblyCommDG()->
                    PerformExchange(outarray, outarray);
            }
            else
            {

                // Loop over elemente and collect forward expansion
                int nexp = GetExpSize();
                int n,p,offset,phys_offset;
                Array<OneD,NekDouble> t_tmp;

                Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >
                    &elmtToTrace = m_traceMap->GetElmtToTrace();

                ASSERTL1(outarray.size() >= m_trace->GetNpoints(),
                         "input array is of insufficient length");

                for (n  = 0; n < nexp; ++n)
                {
                    phys_offset = GetPhys_Offset(n);

                    for (p = 0; p < (*m_exp)[n]->GetNtraces(); ++p)
                    {
                        offset = m_trace->GetPhys_Offset
                            (elmtToTrace[n][p]->GetElmtId());
                        (*m_exp)[n]->GetTracePhysVals(p,elmtToTrace[n][p],
                                                    inarray + phys_offset,
                                                    t_tmp = outarray +offset);
                    }
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
         * @see Expansion3D::AddFaceNormBoundaryInt
         *
         * @param Fn        The trace quantities.
         * @param outarray  Resulting 3D coefficient space.
         *
         */
        void DisContField::v_AddTraceIntegral
        (const Array<OneD, const NekDouble> &Fn,
         Array<OneD,       NekDouble> &outarray)
        {
            if(m_expType == e1D)
            {
                int n,offset, t_offset;

                Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >
                    &elmtToTrace = m_traceMap->GetElmtToTrace();

                vector<bool> negatedFluxNormal = GetNegatedFluxNormal();

                for (n = 0; n < GetExpSize(); ++n)
                {
                    // Number of coefficients on each element
                    int e_ncoeffs = (*m_exp)[n]->GetNcoeffs();

                    offset = GetCoeff_Offset(n);

                    // Implementation for every points except Gauss points
                    if ((*m_exp)[n]->GetBasis(0)->GetBasisType() !=
                        LibUtilities::eGauss_Lagrange)
                    {
                        t_offset = GetTrace()->GetCoeff_Offset
                            (elmtToTrace[n][0]->GetElmtId());
                        if(negatedFluxNormal[2*n])
                        {
                            outarray[offset] -= Fn[t_offset];
                        }
                        else
                        {
                            outarray[offset] += Fn[t_offset];
                        }

                        t_offset = GetTrace()->GetCoeff_Offset
                            (elmtToTrace[n][1]->GetElmtId());

                        if(negatedFluxNormal[2*n+1])
                        {
                            outarray[offset+(*m_exp)[n]->GetVertexMap(1)] -=
                                Fn[t_offset];
                        }
                        else
                        {
                            outarray[offset+(*m_exp)[n]->GetVertexMap(1)] +=
                                Fn[t_offset];
                        }

                    }
                    else
                    {
                        int j;
                        static DNekMatSharedPtr   m_Ixm, m_Ixp;
                        static int sav_ncoeffs = 0;
                        if(!m_Ixm.get() || e_ncoeffs != sav_ncoeffs)
                        {
                            LibUtilities::BasisSharedPtr BASE;
                            const LibUtilities::PointsKey
                                BS_p(e_ncoeffs,LibUtilities::eGaussGaussLegendre);
                            const LibUtilities::BasisKey
                                BS_k(LibUtilities::eGauss_Lagrange,e_ncoeffs,BS_p);

                            BASE  = LibUtilities::BasisManager()[BS_k];

                            Array<OneD, NekDouble> coords(1, 0.0);

                            coords[0] = -1.0;
                            m_Ixm = BASE->GetI(coords);

                            coords[0] = 1.0;
                            m_Ixp = BASE->GetI(coords);

                            sav_ncoeffs = e_ncoeffs;
                        }

                        t_offset = GetTrace()->GetCoeff_Offset
                            (elmtToTrace[n][0]->GetElmtId());

                        if(negatedFluxNormal[2*n])
                        {
                            for (j = 0; j < e_ncoeffs; j++)
                            {
                                outarray[offset + j]  -=
                                    (m_Ixm->GetPtr())[j] * Fn[t_offset];
                            }
                        }
                        else
                        {
                            for (j = 0; j < e_ncoeffs; j++)
                            {
                                outarray[offset + j]  +=
                                    (m_Ixm->GetPtr())[j] * Fn[t_offset];
                            }
                        }

                        t_offset = GetTrace()->GetCoeff_Offset
                            (elmtToTrace[n][1]->GetElmtId());

                        if (negatedFluxNormal[2*n+1])
                        {
                            for (j = 0; j < e_ncoeffs; j++)
                            {
                                outarray[offset + j] -=
                                    (m_Ixp->GetPtr())[j] * Fn[t_offset];
                            }
                        }
                        else
                        {
                            for (j = 0; j < e_ncoeffs; j++)
                            {
                                outarray[offset + j] +=
                                    (m_Ixp->GetPtr())[j] * Fn[t_offset];
                            }
                        }
                    }
                }
            }
            else // other dimensions
            {
                // Basis definition on each element
                LibUtilities::BasisSharedPtr basis = (*m_exp)[0]->GetBasis(0);
                if((m_expType != e1D)&&
                   (basis->GetBasisType() != LibUtilities::eGauss_Lagrange))
                {
                    Array<OneD, NekDouble> Fcoeffs(m_trace->GetNcoeffs());
                    m_trace->IProductWRTBase(Fn, Fcoeffs);

                    m_locTraceToTraceMap->AddTraceCoeffsToFieldCoeffs(Fcoeffs,
                                                                      outarray);
                }
                else
                {
                    int e, n, offset, t_offset;
                    Array<OneD, NekDouble> e_outarray;
                    Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >
                        &elmtToTrace = m_traceMap->GetElmtToTrace();

                    for(n = 0; n < GetExpSize(); ++n)
                    {
                        offset = GetCoeff_Offset(n);
                        for(e = 0; e < (*m_exp)[n]->GetNtraces(); ++e)
                        {
                            t_offset = GetTrace()->GetPhys_Offset
                                (elmtToTrace[n][e]->GetElmtId());
                            (*m_exp)[n]->AddEdgeNormBoundaryInt
                                (e, elmtToTrace[n][e],
                                 Fn+t_offset,
                                 e_outarray = outarray+offset);
                        }
                    }
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
         * IsLeftAdjacentTrace() which if true we use Fwd else we use
         * Bwd
         *
         * @param Fwd       The trace quantities associated with left (fwd)
         *                  adjancent elmt.
         * @param Bwd       The trace quantities associated with right (bwd)
         *                  adjacent elet.
         * @param outarray  Resulting  coefficient space.
         */
        void DisContField::v_AddFwdBwdTraceIntegral(
            const Array<OneD, const NekDouble> &Fwd,
            const Array<OneD, const NekDouble> &Bwd,
                  Array<OneD,       NekDouble> &outarray)
        {

            ASSERTL0(m_expType != e1D, "This method is not setup or "
                     "tested for 1D expansion");

            Array<OneD, NekDouble> Coeffs(m_trace->GetNcoeffs());

            m_trace->IProductWRTBase(Fwd,Coeffs);
            m_locTraceToTraceMap->AddTraceCoeffsToFieldCoeffs(0,Coeffs,outarray);
            m_trace->IProductWRTBase(Bwd,Coeffs);
            m_locTraceToTraceMap->AddTraceCoeffsToFieldCoeffs(1,Coeffs,outarray);
        }

        void DisContField::v_HelmSolve
            (const Array<OneD, const NekDouble> &inarray,
             Array<OneD,       NekDouble>       &outarray,
             const StdRegions::ConstFactorMap   &factors,
             const StdRegions::VarCoeffMap      &varcoeff,
             const MultiRegions::VarFactorsMap  &varfactors,
             const Array<OneD, const NekDouble> &dirForcing,
             const bool                          PhysSpaceForcing)
        {
            boost::ignore_unused(varfactors,dirForcing);
            int i,n,cnt,nbndry;
            int nexp = GetExpSize();

            Array<OneD,NekDouble> f(m_ncoeffs);
            DNekVec F(m_ncoeffs,f,eWrapper);
            Array<OneD,NekDouble> e_f, e_l;

            //----------------------------------
            // Setup RHS Inner product if required
            //----------------------------------
            if(PhysSpaceForcing)
            {
                IProductWRTBase(inarray,f);
                Vmath::Neg(m_ncoeffs,f,1);
            }
            else
            {
                Vmath::Smul(m_ncoeffs,-1.0,inarray,1,f,1);
            }

            //----------------------------------
            // Solve continuous Boundary System
            //----------------------------------
            int GloBndDofs = m_traceMap->GetNumGlobalBndCoeffs();
            int NumDirBCs  = m_traceMap->GetNumLocalDirBndCoeffs();
            int e_ncoeffs;

            GlobalMatrixKey HDGLamToUKey(StdRegions::eHybridDGLamToU,
                                         NullAssemblyMapSharedPtr,factors,varcoeff);
            const DNekScalBlkMatSharedPtr &HDGLamToU =
                GetBlockMatrix(HDGLamToUKey);

            // Retrieve number of local trace space coefficients N_{\lambda},
            // and set up local elemental trace solution \lambda^e.
            int     LocBndCoeffs = m_traceMap->GetNumLocalBndCoeffs();
            Array<OneD, NekDouble> bndrhs(LocBndCoeffs,0.0);
            Array<OneD, NekDouble> loclambda(LocBndCoeffs,0.0);
            DNekVec LocLambda(LocBndCoeffs,loclambda,eWrapper);

            //----------------------------------
            // Evaluate Trace Forcing
            // Kirby et al, 2010, P23, Step 5.
            //----------------------------------
            // Determing <u_lam,f> terms using HDGLamToU matrix
            for (cnt = n = 0; n < nexp; ++n)
            {
                nbndry = (*m_exp)[n]->NumDGBndryCoeffs();

                e_ncoeffs = (*m_exp)[n]->GetNcoeffs();
                e_f       = f + m_coeff_offset[n];
                e_l       = bndrhs + cnt;

                // use outarray as tmp space
                DNekVec     Floc    (nbndry, e_l, eWrapper);
                DNekVec     ElmtFce (e_ncoeffs, e_f, eWrapper);
                Floc = Transpose(*(HDGLamToU->GetBlock(n,n)))*ElmtFce;

                cnt += nbndry;
            }

            Array<OneD, const int> bndCondMap =
                m_traceMap->GetBndCondCoeffsToLocalTraceMap();
            Array<OneD, const NekDouble> Sign =
                m_traceMap->GetLocalToGlobalBndSign();

            // Copy Dirichlet boundary conditions and weak forcing
            // into trace space
            int locid;
            cnt = 0;
            for(i = 0; i < m_bndCondExpansions.size(); ++i)
            {
                Array<OneD, const NekDouble> bndcoeffs =
                    m_bndCondExpansions[i]->GetCoeffs();

                if(m_bndConditions[i]->GetBoundaryConditionType() ==
                   SpatialDomains::eDirichlet)
                {
                    for(int j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); ++j)
                    {
                        locid = bndCondMap[cnt + j];
                        loclambda[locid] = Sign[locid]*bndcoeffs[j];
                    }
                }
                else if (m_bndConditions[i]->GetBoundaryConditionType() ==
                         SpatialDomains::eNeumann ||
                         m_bndConditions[i]->GetBoundaryConditionType() ==
                         SpatialDomains::eRobin)
                {
                    //Add weak boundary condition to trace forcing
                    for(int j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); ++j)
                    {
                        locid = bndCondMap[cnt + j];
                        bndrhs[locid] += Sign[locid]*bndcoeffs[j];
                    }
                }

                cnt += (m_bndCondExpansions[i])->GetNcoeffs();
            }

            //----------------------------------
            // Solve trace problem: \Lambda = K^{-1} F
            // K is the HybridDGHelmBndLam matrix.
            //----------------------------------
            if(GloBndDofs - NumDirBCs > 0)
            {
                GlobalLinSysKey       key(StdRegions::eHybridDGHelmBndLam,
                                          m_traceMap,factors,varcoeff);
                GlobalLinSysSharedPtr LinSys = GetGlobalBndLinSys(key);
                LinSys->Solve(bndrhs,loclambda,m_traceMap);

                // For consistency with previous version put global
                // solution into m_trace->m_coeffs
                m_traceMap->LocalToGlobal(loclambda,m_trace->UpdateCoeffs());
            }

            //----------------------------------
            // Internal element solves
            //----------------------------------
            GlobalMatrixKey invHDGhelmkey(StdRegions::eInvHybridDGHelmholtz,
                                          NullAssemblyMapSharedPtr,
                                          factors, varcoeff);

            const DNekScalBlkMatSharedPtr& InvHDGHelm =
                GetBlockMatrix(invHDGhelmkey);
            DNekVec out(m_ncoeffs,outarray,eWrapper);
            Vmath::Zero(m_ncoeffs,outarray,1);

            //  out =  u_f + u_lam = (*InvHDGHelm)*f + (LamtoU)*Lam
            out = (*InvHDGHelm)*F + (*HDGLamToU)*LocLambda;
        }

        /* \brief This function evaluates the boundary conditions at a certain
         * time-level.
         *
         * Based on the boundary condition \f$g(\boldsymbol{x},t)\f$ evaluated
         * at a given time-level \a t, this function transforms the boundary
         * conditions onto the coefficients of the (one-dimensional) boundary
         * expansion. Depending on the type of boundary conditions, these
         * expansion coefficients are calculated in different ways:
         * - <b>Dirichlet boundary conditions</b><BR>
         *   In order to ensure global \f$C^0\f$ continuity of the spectral/hp
         *   approximation, the Dirichlet boundary conditions are projected onto
         *   the boundary expansion by means of a modified \f$C^0\f$ continuous
         *   Galerkin projection. This projection can be viewed as a collocation
         *   projection at the vertices, followed by an \f$L^2\f$ projection on
         *   the interior modes of the edges. The resulting coefficients
         *   \f$\boldsymbol{\hat{u}}^{\mathcal{D}}\f$ will be stored for the
         *   boundary expansion.
         * - <b>Neumann boundary conditions</b>
         *   In the discrete Galerkin formulation of the problem to be solved,
         *   the Neumann boundary conditions appear as the set of surface
         *   integrals: \f[\boldsymbol{\hat{g}}=\int_{\Gamma}
         *   \phi^e_n(\boldsymbol{x})g(\boldsymbol{x})d(\boldsymbol{x})\quad
         *   \forall n \f]
         *   As a result, it are the coefficients \f$\boldsymbol{\hat{g}}\f$
         *   that will be stored in the boundary expansion
         *
         * @param   time        The time at which the boundary conditions
         *                      should be evaluated.
         * @param   bndCondExpansions   List of boundary conditions.
         * @param   bndConditions   Information about the boundary conditions.
         *
         * This will only be undertaken for time dependent
         * boundary conditions unless time == 0.0 which is the
         * case when the method is called from the constructor.
         */
        void DisContField::v_EvaluateBoundaryConditions
                              (const NekDouble   time,
                               const std::string varName,
                               const NekDouble   x2_in,
                               const NekDouble   x3_in)
        {
            int i;
            int npoints;

            MultiRegions::ExpListSharedPtr locExpList;

            for (i = 0; i < m_bndCondExpansions.size(); ++i)
            {
                if (time == 0.0 || m_bndConditions[i]->IsTimeDependent())
                {
                    m_bndCondBndWeight[i] = 1.0;
                    locExpList = m_bndCondExpansions[i];

                    npoints = locExpList->GetNpoints();
                    Array<OneD, NekDouble> x0(npoints, 0.0);
                    Array<OneD, NekDouble> x1(npoints, 0.0);
                    Array<OneD, NekDouble> x2(npoints, 0.0);

                    locExpList->GetCoords(x0, x1, x2);

                    if (x2_in != NekConstants::kNekUnsetDouble && x3_in !=
                        NekConstants::kNekUnsetDouble)
                    {
                        Vmath::Fill(npoints,x2_in,x1,1);
                        Vmath::Fill(npoints,x3_in,x2,1);
                    }
                    else  if(x2_in != NekConstants::kNekUnsetDouble)
                    {
                        Vmath::Fill(npoints,x2_in,x2,1);
                    }

                    // treat 1D expansions separately since we only
                    // require an evaluation at a point rather than
                    // any projections or inner products that are not
                    // available in a PointExp
                    if(m_expType == e1D)
                    {
                        if (m_bndConditions[i]->GetBoundaryConditionType() ==
                            SpatialDomains::eDirichlet)
                        {
                            m_bndCondExpansions[i]->SetCoeff
                                (0,(std::static_pointer_cast<SpatialDomains
                                    ::DirichletBoundaryCondition>
                                    (m_bndConditions[i])
                                    ->m_dirichletCondition).Evaluate
                                 (x0[0],x1[0],x2[0],time));
                            m_bndCondExpansions[i]->SetPhys
                                (0,m_bndCondExpansions[i]->GetCoeff(0));
                        }
                        else if (m_bndConditions[i]->GetBoundaryConditionType()
                                 == SpatialDomains::eNeumann)
                        {
                            m_bndCondExpansions[i]->SetCoeff
                                (0,(std::static_pointer_cast<SpatialDomains
                                    ::NeumannBoundaryCondition>
                                    (m_bndConditions[i])
                                    ->m_neumannCondition).Evaluate
                                 (x0[0],x1[0],x2[0],time));
                        }
                        else if (m_bndConditions[i]->GetBoundaryConditionType()
                                 == SpatialDomains::eRobin)
                        {
                            m_bndCondExpansions[i]->SetCoeff
                                (0,(std::static_pointer_cast<SpatialDomains
                                    ::RobinBoundaryCondition>
                                    (m_bndConditions[i])
                                    ->m_robinFunction).Evaluate
                                 (x0[0],x1[0],x2[0],time));

                        }
                        else if (m_bndConditions[i]->GetBoundaryConditionType()
                                 == SpatialDomains::ePeriodic)
                        {
                            continue;
                        }
                        else if (m_bndConditions[i]->GetBoundaryConditionType()
                                 == SpatialDomains::eNotDefined)
                        {
                        }
                        else
                        {
                            NEKERROR(ErrorUtil::efatal,
                                     "This type of BC not implemented yet");
                        }
                    }
                    else // 2D and 3D versions
                    {
                        if (m_bndConditions[i]->GetBoundaryConditionType()
                            == SpatialDomains::eDirichlet)
                        {
                            SpatialDomains::DirichletBCShPtr bcPtr
                                = std::static_pointer_cast<
                                    SpatialDomains::DirichletBoundaryCondition>
                                (m_bndConditions[i]);

                            Array<OneD, NekDouble> valuesFile(npoints, 1.0),
                                valuesExp(npoints, 1.0);

                            string filebcs = bcPtr->m_filename;
                            string exprbcs = bcPtr->m_expr;

                            if (filebcs != "")
                            {
                                ExtractFileBCs
                                    (filebcs, bcPtr->GetComm(), varName, locExpList);
                                valuesFile = locExpList->GetPhys();
                            }

                            if (exprbcs != "")
                            {
                                LibUtilities::Equation condition =
                                    std::static_pointer_cast<
                                     SpatialDomains::DirichletBoundaryCondition>
                                    (m_bndConditions[i])->m_dirichletCondition;

                                condition.Evaluate(x0, x1, x2, time, valuesExp);
                            }

                            Vmath::Vmul(npoints, valuesExp, 1, valuesFile, 1,
                                        locExpList->UpdatePhys(), 1);

                            locExpList->FwdTrans_BndConstrained
                                (locExpList->GetPhys(),
                                 locExpList->UpdateCoeffs());
                        }
                        else if (m_bndConditions[i]->GetBoundaryConditionType()
                                 == SpatialDomains::eNeumann)
                        {
                            SpatialDomains::NeumannBCShPtr
                                bcPtr = std::static_pointer_cast<
                                    SpatialDomains::NeumannBoundaryCondition>
                                (m_bndConditions[i]);
                            string filebcs  = bcPtr->m_filename;
                            if (filebcs != "")
                            {
                                ExtractFileBCs
                                    (filebcs, bcPtr->GetComm(), varName, locExpList);
                            }
                            else
                            {
                                LibUtilities::Equation condition =
                                    std::static_pointer_cast<
                                     SpatialDomains::NeumannBoundaryCondition>
                                    (m_bndConditions[i])->m_neumannCondition;
                                condition.Evaluate(x0, x1, x2, time,
                                               locExpList->UpdatePhys());
                            }

                            locExpList->IProductWRTBase(locExpList->GetPhys(),
                                                   locExpList->UpdateCoeffs());
                        }
                        else if (m_bndConditions[i]->GetBoundaryConditionType()
                                 == SpatialDomains::eRobin)
                        {
                            SpatialDomains::RobinBCShPtr
                                bcPtr = std::static_pointer_cast<
                                    SpatialDomains::RobinBoundaryCondition>
                                (m_bndConditions[i]);

                            string filebcs = bcPtr->m_filename;

                            if (filebcs != "")
                            {
                                ExtractFileBCs
                                    (filebcs, bcPtr->GetComm(), varName,
                                     locExpList);
                            }
                            else
                            {
                                LibUtilities::Equation condition =
                                    std::static_pointer_cast<
                                        SpatialDomains::RobinBoundaryCondition>
                                    (m_bndConditions[i])->m_robinFunction;
                                condition.Evaluate(x0, x1, x2, time,
                                               locExpList->UpdatePhys());
                            }

                            locExpList->IProductWRTBase
                                (locExpList->GetPhys(),
                                 locExpList->UpdateCoeffs());
                        }
                        else if (m_bndConditions[i]->GetBoundaryConditionType()
                                 == SpatialDomains::ePeriodic)
                        {
                            continue;
                        }
                        else
                        {
                            NEKERROR(ErrorUtil::efatal,
                                     "This type of BC not implemented yet");
                        }
                    }
                }
            }
        }

        /**
         * @brief Fill the weight with m_bndCondBndWeight.
         */
        void DisContField::v_FillBwdWithBwdWeight(
                  Array<OneD,       NekDouble> &weightave,
                  Array<OneD,       NekDouble> &weightjmp)
        {
            int cnt;
            int npts;
            int e = 0;

            // Fill boundary conditions into missing elements
            int id2 = 0;

            for (int n = cnt = 0; n < m_bndCondExpansions.size(); ++n)
            {

                if (m_bndConditions[n]->GetBoundaryConditionType() ==
                        SpatialDomains::eDirichlet)
                {
                    for (e = 0; e < m_bndCondExpansions[n]->GetExpSize(); ++e)
                    {
                        npts = m_bndCondExpansions[n]->
                                GetExp(e)->GetTotPoints();
                        id2 = m_trace->GetPhys_Offset(m_traceMap->
                                        GetBndCondIDToGlobalTraceID(cnt+e));
                        Vmath::Fill(npts,m_bndCondBndWeight[n],
                                    &weightave[id2], 1);
                        Vmath::Fill(npts, 0.0, &weightjmp[id2], 1);

                    }

                    cnt += e;
                }
                else if (m_bndConditions[n]->GetBoundaryConditionType() ==
                             SpatialDomains::eNeumann ||
                         m_bndConditions[n]->GetBoundaryConditionType() ==
                             SpatialDomains::eRobin)
                {
                    for (e = 0; e < m_bndCondExpansions[n]->GetExpSize(); ++e)
                    {
                        npts = m_bndCondExpansions[n]->
                                GetExp(e)->GetTotPoints();
                        id2 = m_trace->GetPhys_Offset(m_traceMap->
                                        GetBndCondIDToGlobalTraceID(cnt+e));
                        Vmath::Fill(npts,
                                    m_bndCondBndWeight[n],
                                    &weightave[id2], 1);
                        Vmath::Fill(npts, 0.0, &weightjmp[id2], 1);
                    }

                    cnt += e;
                }
                else if (m_bndConditions[n]->GetBoundaryConditionType() !=
                             SpatialDomains::ePeriodic)
                {
                    NEKERROR(ErrorUtil::efatal,
                             "Method not set up for this boundary condition.");
                }
            }
        }


        // Set up a list of element ids and trace ids that link to the
        // boundary conditions
        void DisContField::v_GetBoundaryToElmtMap(
                             Array<OneD, int> &ElmtID,
                             Array<OneD,int>  &TraceID)
        {

            if (m_BCtoElmMap.size() == 0)
            {
                switch(m_expType)
                {
                case e1D:
                {
                    map<int, int> VertGID;
                    int i,n,id;
                    int bid,cnt,Vid;
                    int nbcs = m_bndConditions.size();

                    // make sure arrays are of sufficient length
                    m_BCtoElmMap   = Array<OneD, int>(nbcs,-1);
                    m_BCtoTraceMap = Array<OneD, int>(nbcs);

                    // setup map of all global ids along boundary
                    cnt = 0;
                    for (n = 0; n < m_bndCondExpansions.size(); ++n)
                    {
                        Vid =  m_bndCondExpansions[n]->GetExp(0)->
                            GetGeom()->GetVertex(0)->GetVid();
                        VertGID[Vid] = cnt++;
                    }

                    // Loop over elements and find verts that match;
                    LocalRegions::ExpansionSharedPtr exp;
                    for(cnt = n = 0; n < GetExpSize(); ++n)
                    {
                        exp = (*m_exp)[n];
                        for(i = 0; i < exp->GetNverts(); ++i)
                        {
                            id = exp->GetGeom()->GetVid(i);

                            if (VertGID.count(id) > 0)
                            {
                                bid = VertGID.find(id)->second;
                                ASSERTL1(m_BCtoElmMap[bid] == -1,
                                         "Edge already set");
                                m_BCtoElmMap  [bid] = n;
                                m_BCtoTraceMap[bid] = i;
                                cnt ++;
                            }
                        }
                    }
                    ASSERTL1(cnt == nbcs,
                             "Failed to visit all boundary condtiions");
                }
                break;
                case e2D:
                {
                    map<int, int> globalIdMap;
                    int i,n;
                    int cnt;
                    int nbcs = 0;

                    // Populate global ID map (takes global geometry
                    // ID to local expansion list ID).
                    for (i = 0; i < GetExpSize(); ++i)
                    {
                        globalIdMap[(*m_exp)[i]->GetGeom()->GetGlobalID()] = i;
                    }

                    // Determine number of boundary condition expansions.
                    for(i = 0; i < m_bndConditions.size(); ++i)
                    {
                        nbcs += m_bndCondExpansions[i]->GetExpSize();
                    }

                    // Initialize arrays
                    m_BCtoElmMap   = Array<OneD, int>(nbcs);
                    m_BCtoTraceMap = Array<OneD, int>(nbcs);

                    LocalRegions::Expansion1DSharedPtr exp1d;
                    cnt = 0;
                    for (n = 0; n < m_bndCondExpansions.size(); ++n)
                    {
                        for (i = 0; i < m_bndCondExpansions[n]->GetExpSize();
                             ++i, ++cnt)
                        {
                            exp1d = m_bndCondExpansions[n]->GetExp(i)->
                                as<LocalRegions::Expansion1D>();

                            // Use edge to element map from MeshGraph.
                            SpatialDomains::GeometryLinkSharedPtr tmp =
                                m_graph->GetElementsFromEdge(exp1d->GetGeom1D());

                            m_BCtoElmMap[cnt] = globalIdMap[
                                  (*tmp)[0].first->GetGlobalID()];
                            m_BCtoTraceMap[cnt] = (*tmp)[0].second;
                        }
                    }
                }
                break;
                case e3D:
                {
                    map<int,int> globalIdMap;
                    int i, n;
                    int cnt;
                    int nbcs = 0;

                    // Populate global ID map (takes global geometry ID to local
                    // expansion list ID).
                    LocalRegions::Expansion3DSharedPtr exp3d;
                    for (i = 0; i < GetExpSize(); ++i)
                    {
                        globalIdMap[(*m_exp)[i]->GetGeom()->GetGlobalID()] = i;
                    }

                    // Determine number of boundary condition expansions.
                    for(i = 0; i < m_bndConditions.size(); ++i)
                    {
                        nbcs += m_bndCondExpansions[i]->GetExpSize();
                    }

                    // Initialize arrays
                    m_BCtoElmMap = Array<OneD, int>(nbcs);
                    m_BCtoTraceMap = Array<OneD, int>(nbcs);

                    LocalRegions::Expansion2DSharedPtr exp2d;
                    for(cnt = n = 0; n < m_bndCondExpansions.size(); ++n)
                    {
                        for(i = 0; i < m_bndCondExpansions[n]->GetExpSize(); ++i, ++cnt)
                        {
                            exp2d = m_bndCondExpansions[n]->GetExp(i)->
                                as<LocalRegions::Expansion2D>();

                            SpatialDomains::GeometryLinkSharedPtr tmp =
                                m_graph->GetElementsFromFace(exp2d->GetGeom2D());
                            m_BCtoElmMap[cnt] = globalIdMap
                                [tmp->at(0).first->GetGlobalID()];
                            m_BCtoTraceMap[cnt] = tmp->at(0).second;
                        }
                    }
                }
              break;
                default:
                        ASSERTL1(false,"Not setup for this expansion");
                    break;
                }
            }

            ElmtID = m_BCtoElmMap;
            TraceID = m_BCtoTraceMap;
        }

        void DisContField::v_GetBndElmtExpansion(int i,
                            std::shared_ptr<ExpList> &result,
                            const bool DeclareCoeffPhysArrays)
        {
            int n, cnt, nq;
            int offsetOld, offsetNew;
            std::vector<unsigned int> eIDs;

            Array<OneD, int> ElmtID,TraceID;
            GetBoundaryToElmtMap(ElmtID,TraceID);

            // Skip other boundary regions
            for (cnt = n = 0; n < i; ++n)
            {
                cnt += m_bndCondExpansions[n]->GetExpSize();
            }

            // Populate eIDs with information from BoundaryToElmtMap
            for (n = 0; n < m_bndCondExpansions[i]->GetExpSize(); ++n)
            {
                eIDs.push_back(ElmtID[cnt+n]);
            }

            // Create expansion list
            result =
                MemoryManager<ExpList>::AllocateSharedPtr
                    (*this, eIDs, DeclareCoeffPhysArrays);

            // Copy phys and coeffs to new explist
            if( DeclareCoeffPhysArrays)
            {
                Array<OneD, NekDouble> tmp1, tmp2;
                for (n = 0; n < result->GetExpSize(); ++n)
                {
                    nq = GetExp(ElmtID[cnt+n])->GetTotPoints();
                    offsetOld = GetPhys_Offset(ElmtID[cnt+n]);
                    offsetNew = result->GetPhys_Offset(n);
                    Vmath::Vcopy(nq, tmp1 = GetPhys()+ offsetOld, 1,
                                tmp2 = result->UpdatePhys()+ offsetNew, 1);

                    nq = GetExp(ElmtID[cnt+n])->GetNcoeffs();
                    offsetOld = GetCoeff_Offset(ElmtID[cnt+n]);
                    offsetNew = result->GetCoeff_Offset(n);
                    Vmath::Vcopy(nq, tmp1 = GetCoeffs()+ offsetOld, 1,
                                tmp2 = result->UpdateCoeffs()+ offsetNew, 1);
                }
            }
        }

        /**
         * @brief Reset this field, so that geometry information can be updated.
         */
        void DisContField::v_Reset()
        {
            ExpList::v_Reset();

            // Reset boundary condition expansions.
            for (int n = 0; n < m_bndCondExpansions.size(); ++n)
            {
                m_bndCondExpansions[n]->Reset();
            }
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
        map<int, RobinBCInfoSharedPtr> DisContField::v_GetRobinBCInfo(void)
        {
            int i,cnt;
            map<int, RobinBCInfoSharedPtr> returnval;
            Array<OneD, int> ElmtID,TraceID;
            GetBoundaryToElmtMap(ElmtID,TraceID);

            for(cnt = i = 0; i < m_bndCondExpansions.size(); ++i)
            {
                MultiRegions::ExpListSharedPtr locExpList;

                if(m_bndConditions[i]->GetBoundaryConditionType() ==
                       SpatialDomains::eRobin)
                {
                    int e,elmtid;
                    Array<OneD, NekDouble> Array_tmp;

                    locExpList = m_bndCondExpansions[i];

                    int npoints    = locExpList->GetNpoints();
                    Array<OneD, NekDouble> x0(npoints, 0.0);
                    Array<OneD, NekDouble> x1(npoints, 0.0);
                    Array<OneD, NekDouble> x2(npoints, 0.0);
                    Array<OneD, NekDouble> coeffphys(npoints);

                    locExpList->GetCoords(x0, x1, x2);

                    LibUtilities::Equation coeffeqn =
                        std::static_pointer_cast<
                            SpatialDomains::RobinBoundaryCondition>
                        (m_bndConditions[i])->m_robinPrimitiveCoeff;

                    // evalaute coefficient
                    coeffeqn.Evaluate(x0, x1, x2, 0.0, coeffphys);

                    for(e = 0; e < locExpList->GetExpSize(); ++e)
                    {
                        RobinBCInfoSharedPtr rInfo =
                            MemoryManager<RobinBCInfo>
                            ::AllocateSharedPtr(
                                TraceID[cnt+e],
                                Array_tmp = coeffphys +
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
         * @brief Calculate the \f$ L^2 \f$ error of the \f$ Q_{\rm dir} \f$
         * derivative using the consistent DG evaluation of \f$ Q_{\rm dir} \f$.
         *
         * The solution provided is of the primative variation at the quadrature
         * points and the derivative is compared to the discrete derivative at
         * these points, which is likely to be undesirable unless using a much
         * higher number of quadrature points than the polynomial order used to
         * evaluate \f$ Q_{\rm dir} \f$.
        */
        NekDouble DisContField::L2_DGDeriv(
            const int                           dir,
            const Array<OneD, const NekDouble> &soln)
        {
            int    i,e,ncoeff_edge;
            Array<OneD, const NekDouble> tmp_coeffs;
            Array<OneD, NekDouble> out_d(m_ncoeffs), out_tmp;

            Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >
                &elmtToTrace = m_traceMap->GetElmtToTrace();

            StdRegions::Orientation edgedir;

            int     cnt;
            int     LocBndCoeffs = m_traceMap->GetNumLocalBndCoeffs();
            Array<OneD, NekDouble> loc_lambda(LocBndCoeffs), edge_lambda;


            m_traceMap->GlobalToLocalBnd(m_trace->GetCoeffs(),loc_lambda);

            edge_lambda = loc_lambda;

            // Calculate Q using standard DG formulation.
            for(i = cnt = 0; i < GetExpSize(); ++i)
            {
                // Probably a better way of setting up lambda than this.
                // Note cannot use PutCoeffsInToElmts since lambda space
                // is mapped during the solve.
                int nEdges = (*m_exp)[i]->GetGeom()->GetNumEdges();
                Array<OneD, Array<OneD, NekDouble> > edgeCoeffs(nEdges);

                for(e = 0; e < nEdges; ++e)
                {
                    edgedir = (*m_exp)[i]->GetTraceOrient(e);
                    ncoeff_edge = elmtToTrace[i][e]->GetNcoeffs();
                    edgeCoeffs[e] = Array<OneD, NekDouble>(ncoeff_edge);
                    Vmath::Vcopy(ncoeff_edge, edge_lambda, 1, edgeCoeffs[e], 1);
                    elmtToTrace[i][e]->SetCoeffsToOrientation(
                        edgedir, edgeCoeffs[e], edgeCoeffs[e]);
                    edge_lambda = edge_lambda + ncoeff_edge;
                }

                (*m_exp)[i]->DGDeriv(dir,
                                       tmp_coeffs=m_coeffs+m_coeff_offset[i],
                                       elmtToTrace[i],
                                       edgeCoeffs,
                                       out_tmp = out_d+cnt);
                cnt  += (*m_exp)[i]->GetNcoeffs();
            }

            BwdTrans(out_d,m_phys);
            Vmath::Vsub(m_npoints,m_phys,1,soln,1,m_phys,1);
            return L2(m_phys);
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
        void  DisContField::EvaluateHDGPostProcessing(
            Array<OneD, NekDouble> &outarray)
        {
            int    i,cnt,e,ncoeff_trace;
            Array<OneD, NekDouble> force, out_tmp, qrhs, qrhs1;
            Array<OneD, Array< OneD, LocalRegions::ExpansionSharedPtr> >
                &elmtToTrace = m_traceMap->GetElmtToTrace();

            StdRegions::Orientation edgedir;

            int     nq_elmt, nm_elmt;
            int     LocBndCoeffs = m_traceMap->GetNumLocalBndCoeffs();
            Array<OneD, NekDouble> loc_lambda(LocBndCoeffs), trace_lambda;
            Array<OneD, NekDouble> tmp_coeffs;
            m_traceMap->GlobalToLocalBnd(m_trace->GetCoeffs(),loc_lambda);

            trace_lambda = loc_lambda;

            int dim  = (m_expType == e2D)? 2:3;

            int num_points[] = {0,0,0};
            int num_modes [] = {0,0,0};

            // Calculate Q using standard DG formulation.
            for(i = cnt = 0; i < GetExpSize(); ++i)
            {
                nq_elmt = (*m_exp)[i]->GetTotPoints();
                nm_elmt = (*m_exp)[i]->GetNcoeffs();
                qrhs   = Array<OneD, NekDouble>(nq_elmt);
                qrhs1  = Array<OneD, NekDouble>(nq_elmt);
                force  = Array<OneD, NekDouble>(2*nm_elmt);
                out_tmp = force + nm_elmt;
                LocalRegions::ExpansionSharedPtr ppExp;

                for(int j= 0; j < dim; ++j)
                {
                    num_points[j] = (*m_exp)[i]->GetBasis(j)->GetNumPoints();
                    num_modes[j]  = (*m_exp)[i]->GetBasis(j)->GetNumModes();
                }

                // Probably a better way of setting up lambda than this.  Note
                // cannot use PutCoeffsInToElmts since lambda space is mapped
                // during the solve.
                int nTraces = (*m_exp)[i]->GetNtraces();
                Array<OneD, Array<OneD, NekDouble> > traceCoeffs(nTraces);

                for(e = 0; e < (*m_exp)[i]->GetNtraces(); ++e)
                {
                    edgedir = (*m_exp)[i]->GetTraceOrient(e);
                    ncoeff_trace = elmtToTrace[i][e]->GetNcoeffs();
                    traceCoeffs[e] = Array<OneD, NekDouble>(ncoeff_trace);
                    Vmath::Vcopy(ncoeff_trace, trace_lambda, 1, traceCoeffs[e], 1);
                    if(dim == 2)
                    {
                        elmtToTrace[i][e]->SetCoeffsToOrientation(edgedir,
                                                            traceCoeffs[e], traceCoeffs[e]);
                    }
                    else
                    {
                        (*m_exp)[i]->as<LocalRegions::Expansion3D>()->
                            SetFaceToGeomOrientation(e,traceCoeffs[e]);

                    }
                    trace_lambda = trace_lambda + ncoeff_trace;
                }

                //creating orthogonal expansion (checking if we have quads or triangles)
                LibUtilities::ShapeType shape = (*m_exp)[i]->DetShapeType();
                switch(shape)
                {
                    case LibUtilities::eQuadrilateral:
                    {
                        const LibUtilities::PointsKey PkeyQ1(num_points[0],
                                                       LibUtilities::eGaussLobattoLegendre);
                        const LibUtilities::PointsKey PkeyQ2(num_points[1],
                                                      LibUtilities::eGaussLobattoLegendre);
                        LibUtilities::BasisKey  BkeyQ1(LibUtilities::eOrtho_A, num_modes[0],
                                                       PkeyQ1);
                        LibUtilities::BasisKey  BkeyQ2(LibUtilities::eOrtho_A, num_modes[1],
                                                       PkeyQ2);
                        SpatialDomains::QuadGeomSharedPtr qGeom =
                            std::dynamic_pointer_cast<SpatialDomains::QuadGeom>
                            ((*m_exp)[i]->GetGeom());
                        ppExp = MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr
                            (BkeyQ1, BkeyQ2, qGeom);
                    }
                    break;
                    case LibUtilities::eTriangle:
                    {
                        const LibUtilities::PointsKey PkeyT1(num_points[0],
                                                            LibUtilities::eGaussLobattoLegendre);
                        const LibUtilities::PointsKey PkeyT2(num_points[1],
                                                             LibUtilities::eGaussRadauMAlpha1Beta0);
                        LibUtilities::BasisKey  BkeyT1(LibUtilities::eOrtho_A, num_modes[0],
                                                       PkeyT1);
                        LibUtilities::BasisKey  BkeyT2(LibUtilities::eOrtho_B, num_modes[1],
                                                       PkeyT2);
                        SpatialDomains::TriGeomSharedPtr tGeom = std::dynamic_pointer_cast
                            <SpatialDomains::TriGeom>((*m_exp)[i]->GetGeom());
                        ppExp = MemoryManager<LocalRegions::TriExp>::AllocateSharedPtr
                            (BkeyT1, BkeyT2, tGeom);
                    }
                    break;
                    case LibUtilities::eHexahedron:
                    {
                        const LibUtilities::PointsKey PkeyH1(num_points[0],
                                                      LibUtilities::eGaussLobattoLegendre);
                        const LibUtilities::PointsKey PkeyH2(num_points[1],
                                                      LibUtilities::eGaussLobattoLegendre);
                        const LibUtilities::PointsKey PkeyH3(num_points[2],
                                                      LibUtilities::eGaussLobattoLegendre);
                        LibUtilities::BasisKey  BkeyH1(LibUtilities::eOrtho_A,
                                                       num_modes[0], PkeyH1);
                        LibUtilities::BasisKey  BkeyH2(LibUtilities::eOrtho_A,
                                                       num_modes[1], PkeyH2);
                        LibUtilities::BasisKey  BkeyH3(LibUtilities::eOrtho_A,
                                                       num_modes[2], PkeyH3);
                        SpatialDomains::HexGeomSharedPtr hGeom =
                            std::dynamic_pointer_cast<SpatialDomains::HexGeom>
                            ((*m_exp)[i]->GetGeom());
                        ppExp = MemoryManager<LocalRegions::HexExp>::AllocateSharedPtr
                            (BkeyH1, BkeyH2, BkeyH3, hGeom);
                    }
                    break;
                    case LibUtilities::eTetrahedron:
                    {
                        const LibUtilities::PointsKey PkeyT1(num_points[0],
                                                        LibUtilities::eGaussLobattoLegendre);
                        const LibUtilities::PointsKey PkeyT2(num_points[1],
                                                       LibUtilities::eGaussRadauMAlpha1Beta0);
                        const LibUtilities::PointsKey PkeyT3(num_points[2],
                                                       LibUtilities::eGaussRadauMAlpha2Beta0);
                        LibUtilities::BasisKey  BkeyT1(LibUtilities::eOrtho_A,
                                                       num_modes[0], PkeyT1);
                        LibUtilities::BasisKey  BkeyT2(LibUtilities::eOrtho_B,
                                                       num_modes[1], PkeyT2);
                        LibUtilities::BasisKey  BkeyT3(LibUtilities::eOrtho_C,
                                                       num_modes[2], PkeyT3);
                        SpatialDomains::TetGeomSharedPtr tGeom =
                            std::dynamic_pointer_cast<SpatialDomains::TetGeom>
                            ((*m_exp)[i]->GetGeom());
                        ppExp = MemoryManager<LocalRegions::TetExp>::AllocateSharedPtr
                            (BkeyT1, BkeyT2, BkeyT3, tGeom);
                    }
                    break;
                default:
                    NEKERROR(ErrorUtil::efatal, 
                        "Wrong shape type, HDG postprocessing is not "
                        "implemented");
                };


                //DGDeriv
                // (d/dx w, d/dx q_0)
                (*m_exp)[i]->DGDeriv(
                    0,tmp_coeffs = m_coeffs + m_coeff_offset[i],
                    elmtToTrace[i], traceCoeffs, out_tmp);
                (*m_exp)[i]->BwdTrans(out_tmp,qrhs);
                ppExp->IProductWRTDerivBase(0,qrhs,force);

                // + (d/dy w, d/dy q_1)
                (*m_exp)[i]->DGDeriv(
                    1,tmp_coeffs = m_coeffs + m_coeff_offset[i],
                    elmtToTrace[i], traceCoeffs, out_tmp);

                (*m_exp)[i]->BwdTrans(out_tmp,qrhs);
                ppExp->IProductWRTDerivBase(1,qrhs,out_tmp);

                Vmath::Vadd(nm_elmt,force,1,out_tmp,1,force,1);

                // determine force[0] = (1,u)
                (*m_exp)[i]->BwdTrans(
                    tmp_coeffs = m_coeffs + m_coeff_offset[i],qrhs);
                force[0] = (*m_exp)[i]->Integral(qrhs);

                // multiply by inverse Laplacian matrix
                // get matrix inverse
                LocalRegions::MatrixKey  lapkey(StdRegions::eInvLaplacianWithUnityMean,
                                                ppExp->DetShapeType(), *ppExp);
                DNekScalMatSharedPtr lapsys = ppExp->GetLocMatrix(lapkey);

                NekVector<NekDouble> in (nm_elmt,force,eWrapper);
                NekVector<NekDouble> out(nm_elmt);

                out = (*lapsys)*in;

                // Transforming back to modified basis
                Array<OneD, NekDouble> work(nq_elmt);
                ppExp->BwdTrans(out.GetPtr(), work);
                (*m_exp)[i]->FwdTrans(work, tmp_coeffs = outarray + m_coeff_offset[i]);
            }
        }

        void DisContField::v_GetLocTraceFromTracePts(
                const Array<OneD, const NekDouble>  &Fwd,
                const Array<OneD, const NekDouble>  &Bwd,
                Array<OneD,       NekDouble>        &locTraceFwd,
                Array<OneD,       NekDouble>        &locTraceBwd)
        {
            if (NullNekDouble1DArray != locTraceBwd)
            {
                switch(m_expType)
                {
                case e2D:
                    m_locTraceToTraceMap->RightIPTWLocEdgesToTraceInterpMat(
                        1, Bwd, locTraceBwd);
                    break;
                case e3D:
                    m_locTraceToTraceMap->RightIPTWLocFacesToTraceInterpMat(
                        1, Bwd, locTraceBwd);
                    break;
                default:
                    NEKERROR(ErrorUtil::efatal, 
                        "GetLocTraceFromTracePts not defined");
                }
            }

            if (NullNekDouble1DArray != locTraceFwd)
            {
                switch(m_expType)
                {
                case e2D:
                    m_locTraceToTraceMap->RightIPTWLocEdgesToTraceInterpMat(
                        0, Fwd, locTraceFwd);
                    break;
                case e3D:
                    m_locTraceToTraceMap->RightIPTWLocFacesToTraceInterpMat(
                        0, Fwd, locTraceFwd);
                    break;
                default:
                    NEKERROR(ErrorUtil::efatal, 
                        "GetLocTraceFromTracePts not defined");
                }
            }
        }

        void DisContField::v_AddTraceIntegralToOffDiag(
            const Array<OneD, const NekDouble> &FwdFlux, 
            const Array<OneD, const NekDouble> &BwdFlux, 
                  Array<OneD,       NekDouble> &outarray)
        {
            Array<OneD, NekDouble> FCoeffs(m_trace->GetNcoeffs());

            m_trace->IProductWRTBase(FwdFlux,FCoeffs);
            m_locTraceToTraceMap->AddTraceCoeffsToFieldCoeffs(1, FCoeffs, 
                outarray);
            m_trace->IProductWRTBase(BwdFlux,FCoeffs);
            m_locTraceToTraceMap->AddTraceCoeffsToFieldCoeffs(0, FCoeffs, 
                outarray);
        }
    } // end of namespace
} //end of namespace
