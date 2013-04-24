//////////////////////////////////////////////////////////////////////////////
//
// File DisContField3D.cpp
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
// Description: Field definition for 3D domain with boundary
// conditions using LDG flux
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/DisContField3D.h>
#include <LocalRegions/Expansion3D.h>
#include <LocalRegions/Expansion2D.h>
#include <SpatialDomains/MeshGraph3D.h>

#include <boost/assign/std/vector.hpp>
#include <boost/tuple/tuple.hpp>

using namespace boost::assign;

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class DisContField3D
         * Abstraction of a global discontinuous three-dimensional spectral/hp
         * element expansion which approximates the solution of a set of
         * partial differential equations.
         */

        /**
         * @brief Default constructor.
         */
        DisContField3D::DisContField3D() :
            ExpList3D             (),
            m_bndCondExpansions   (),
            m_bndConditions       (),
            m_trace(NullExpListSharedPtr)
        {
        }

        /**
         * @brief Constructs a global discontinuous field based on an input mesh
         * with boundary conditions.
         */
        DisContField3D::DisContField3D(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const SpatialDomains::MeshGraphSharedPtr   &graph3D,
            const std::string                          &variable,
            const bool                                  SetUpJustDG)
            : ExpList3D          (pSession,graph3D,variable),
              m_bndCondExpansions(),
              m_bndConditions    (),
              m_trace(NullExpListSharedPtr)
        {
            SpatialDomains::BoundaryConditions bcs(m_session, graph3D);
            
            GenerateBoundaryConditionExpansion(graph3D,bcs,variable);
            EvaluateBoundaryConditions();
            ApplyGeomInfo();
            
            // Find periodic edges for this variable.
            FindPeriodicFaces(bcs, variable);

            if(SetUpJustDG)
            {
                SetUpDG();
            }
            else
            {
                // Set element edges to point to Robin BC edges if required.
                int i,cnt,f;
                Array<OneD, int> ElmtID, FaceID;
                GetBoundaryToElmtMap(ElmtID, FaceID);

                for(cnt = i = 0; i < m_bndCondExpansions.num_elements(); ++i)
                {
                    MultiRegions::ExpListSharedPtr locExpList;
                    locExpList = m_bndCondExpansions[i];
                    
                    for(f = 0; f < locExpList->GetExpSize(); ++f)
                    {
                        LocalRegions::Expansion3DSharedPtr exp3d
                            = boost::dynamic_pointer_cast<
                                LocalRegions::Expansion3D>((*m_exp)[ElmtID[cnt+f]]);
                        LocalRegions::Expansion2DSharedPtr exp2d
                            = boost::dynamic_pointer_cast<
                                LocalRegions::Expansion2D>(locExpList->GetExp(f));
                        
                        exp3d->SetFaceExp(FaceID[cnt+f],exp2d);
                        exp2d->SetAdjacentElementExp(FaceID[cnt+f],exp3d);
                    }
                    cnt += m_bndCondExpansions[i]->GetExpSize();
                }
            }
        }
        
        /*
         * @brief Copy type constructor which declares new boundary conditions
         * and re-uses mapping info and trace space if possible
         */
        DisContField3D::DisContField3D( 
            const DisContField3D                     &In,
            const SpatialDomains::MeshGraphSharedPtr &graph3D,
            const std::string                        &variable,
            const bool                                SetUpJustDG) 
            : ExpList3D(In), 
              m_trace(NullExpListSharedPtr)
        {
            SpatialDomains::BoundaryConditions bcs(m_session, graph3D);
            
            GenerateBoundaryConditionExpansion(graph3D,bcs,variable);
            EvaluateBoundaryConditions();
            ApplyGeomInfo();
           
            if(!SameTypeOfBoundaryConditions(In))
            {
                // Find periodic edges for this variable.
                FindPeriodicFaces(bcs, variable);
               
                if (SetUpJustDG)
                {
                    SetUpDG(variable);
                }
                else
                {
                    int i,cnt,f;
                    Array<OneD, int> ElmtID,FaceID;
                    GetBoundaryToElmtMap(ElmtID,FaceID);
                   
                    for(cnt = i = 0; i < m_bndCondExpansions.num_elements(); ++i)
                    {
                        MultiRegions::ExpListSharedPtr locExpList;
                        locExpList = m_bndCondExpansions[i];
                       
                        for(f = 0; f < locExpList->GetExpSize(); ++f)
                        {
                            LocalRegions::Expansion3DSharedPtr exp3d
                                = boost::dynamic_pointer_cast<
                                    LocalRegions::Expansion3D>(
                                        (*m_exp)[ElmtID[cnt+f]]);
                            LocalRegions::Expansion2DSharedPtr exp2d
                                = boost::dynamic_pointer_cast<
                                    LocalRegions::Expansion2D>(
                                        locExpList->GetExp(f));
                           
                            exp3d->SetFaceExp(FaceID[cnt+f],exp2d);
                            exp2d->SetAdjacentElementExp(FaceID[cnt+f],exp3d);
                        }
                       
                        cnt += m_bndCondExpansions[i]->GetExpSize();
                    }
                    SetUpPhysNormals();
                }
               
            }
            //else if we have the same boundary condition
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
                   
                    int i,cnt,f;
                    Array<OneD, int> ElmtID,FaceID;
                    GetBoundaryToElmtMap(ElmtID,FaceID);
                   
                    for(cnt = i = 0; i < m_bndCondExpansions.num_elements(); ++i)
                    {
                        MultiRegions::ExpListSharedPtr locExpList;
                        locExpList = m_bndCondExpansions[i];
                       
                        for(f = 0; f < locExpList->GetExpSize(); ++f)
                        {
                            LocalRegions::Expansion3DSharedPtr exp3d
                                = boost::dynamic_pointer_cast<
                                    LocalRegions::Expansion3D>(
                                        (*m_exp)[ElmtID[cnt+f]]);
                            LocalRegions::Expansion2DSharedPtr exp2d
                                = boost::dynamic_pointer_cast<
                                    LocalRegions::Expansion2D>(
                                        locExpList->GetExp(f));
                           
                            exp3d->SetFaceExp(FaceID[cnt+f],exp2d);
                            exp2d->SetAdjacentElementExp(FaceID[cnt+f],exp3d);
                        }
                       
                        cnt += m_bndCondExpansions[i]->GetExpSize();
                    }

                    if(m_session->DefinesSolverInfo("PROJECTION"))
                    {
                        std::string ProjectStr = 
                            m_session->GetSolverInfo("PROJECTION");
                        if (ProjectStr == "MixedCGDG"           ||
                            ProjectStr == "Mixed_CG_Discontinuous")
                        {
                            SetUpDG(variable);
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
        }

        /**
         *
         */
        DisContField3D::DisContField3D(const DisContField3D &In) :
            ExpList3D(In),
            m_bndCondExpansions   (In.m_bndCondExpansions),
            m_bndConditions       (In.m_bndConditions),
            m_globalBndMat        (In.m_globalBndMat),
            m_trace               (In.m_trace),
            m_traceMap            (In.m_traceMap)
        {
        }

        /**
         * @brief Destructor.
         */
        DisContField3D::~DisContField3D()
        {
        }

        GlobalLinSysSharedPtr DisContField3D::GetGlobalBndLinSys(
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
         * @brief Set up all DG member variables and maps.
         */
        void DisContField3D::SetUpDG(const std::string variable)
        {
            if (m_trace != NullExpListSharedPtr)
            {
                return;
            }
            
            ExpList2DSharedPtr trace;
            
            SpatialDomains::MeshGraph3DSharedPtr graph3D = 
                boost::dynamic_pointer_cast<SpatialDomains::MeshGraph3D>(
                    m_graph);
            
            // Set up matrix map
            m_globalBndMat = MemoryManager<GlobalLinSysMap>::
                AllocateSharedPtr();
            
            // Set up Trace space
            bool UseGenSegExp = true;
            trace = MemoryManager<ExpList2D>::AllocateSharedPtr(
                m_bndCondExpansions, m_bndConditions,
                *m_exp,graph3D, m_periodicFaces, UseGenSegExp);

            m_trace    = trace;
            m_traceMap = MemoryManager<AssemblyMapDG>::AllocateSharedPtr(
                m_session,graph3D,trace,*this,m_bndCondExpansions,
                m_bndConditions, m_periodicFaces,variable);

            Array<OneD, Array<OneD, StdRegions::StdExpansionSharedPtr> >
                &elmtToTrace = m_traceMap->GetElmtToTrace();
            
            // Scatter trace segments to 3D elements. For each element, we find
            // the trace segment associated to each edge. The element then
            // retains a pointer to the trace space segments, to ensure
            // uniqueness of normals when retrieving from two adjoining elements
            // which do not lie in a plane.
            for (int i = 0; i < m_exp->size(); ++i)
            {
                for (int j = 0; j < (*m_exp)[i]->GetNfaces(); ++j)
                {
                    LocalRegions::Expansion3DSharedPtr exp3d =
                        boost::dynamic_pointer_cast<
                            LocalRegions::Expansion3D>((*m_exp)[i]);
                    LocalRegions::Expansion2DSharedPtr exp2d =
                        boost::dynamic_pointer_cast<
                            LocalRegions::Expansion2D>(elmtToTrace[i][j]);
                    exp3d->SetFaceExp           (j, exp2d);
                    exp2d->SetAdjacentElementExp(j, exp3d);
                }
            }
            
            // Set up physical normals
            SetUpPhysNormals();
            
            // Set up information for parallel jobs.
            for (int i = 0; i < m_trace->GetExpSize(); ++i)
            {
                LocalRegions::Expansion2DSharedPtr traceEl = 
                    boost::dynamic_pointer_cast<
                        LocalRegions::Expansion2D>(m_trace->GetExp(i));
                    
                int offset = m_trace->GetPhys_Offset(i);
                
                if (m_traceMap->GetTraceToUniversalMapUnique(offset) < 0)
                {
                    traceEl->GetLeftAdjacentElementExp()->NegateFaceNormal(
                        traceEl->GetLeftAdjacentElementFace());
                }
            }
            
            int cnt, n, e;
            
            // Identify boundary faces
            for(cnt = 0, n = 0; n < m_bndCondExpansions.num_elements(); ++n)
            {
                if (m_bndConditions[n]->GetBoundaryConditionType() != 
                    SpatialDomains::ePeriodic)
                {
                    for(e = 0; e < m_bndCondExpansions[n]->GetExpSize(); ++e)
                    {
                        m_boundaryFaces.insert(m_trace->GetOffset_Elmt_Id(
                            m_traceMap->GetBndCondTraceToGlobalTraceMap(cnt+e)));
                    }
                }
                cnt += m_bndCondExpansions[n]->GetExpSize();
            }
                
            // Set up information for periodic boundary conditions.
            boost::unordered_map<int,pair<int,int> > perFaceToExpMap;
            boost::unordered_map<int,pair<int,int> >::iterator it2;
            cnt = 0;
            for (int n = 0; n < m_exp->size(); ++n)
            {
                for (int e = 0; e < (*m_exp)[n]->GetNfaces(); ++e, ++cnt)
                {
                    PeriodicMap::iterator it = m_periodicFaces.find(
                        (*m_exp)[n]->GetGeom3D()->GetFid(e));

                    if (it != m_periodicFaces.end())
                    {
                        perFaceToExpMap[it->first] = make_pair(n, e);
                    }
                }
            }
            
            // Set up left-adjacent face list.
            m_leftAdjacentFaces.resize(cnt);
            cnt = 0;
            for (int i = 0; i < m_exp->size(); ++i)
            {
                for (int j = 0; j < (*m_exp)[i]->GetNfaces(); ++j, ++cnt)
                {
                    m_leftAdjacentFaces[cnt] = IsLeftAdjacentFace(i, j);
                }
            }

            // Set up mapping to copy Fwd of periodic bcs to Bwd of other edge.
            cnt = 0;
            for (int n = 0; n < m_exp->size(); ++n)
            {
                for (int e = 0; e < (*m_exp)[n]->GetNfaces(); ++e, ++cnt)
                {
                    int faceGeomId = (*m_exp)[n]->GetGeom3D()->GetFid(e);
                    int offset = m_trace->GetPhys_Offset(
                        elmtToTrace[n][e]->GetElmtId());

                    // Check to see if this face is periodic.
                    PeriodicMap::iterator it = m_periodicFaces.find(faceGeomId);

                    if (it != m_periodicFaces.end())
                    {
                        const PeriodicEntity &ent = it->second[0];
                        it2 = perFaceToExpMap.find(ent.id);

                        if (it2 == perFaceToExpMap.end())
                        {
                            if (m_session->GetComm()->GetSize() > 1 &&
                                !ent.isLocal)
                            {
                                continue;
                            }
                            else
                            {
                                ASSERTL1(false, "Periodic edge not found!");
                            }
                        }

                        ASSERTL1(m_leftAdjacentFaces[cnt],
                                 "Periodic face in non-forward space?");

                        int offset2 = m_trace->GetPhys_Offset(
                            elmtToTrace[it2->second.first][it2->second.second]->
                                GetElmtId());

                        // Calculate relative orientations between faces to
                        // calculate copying map.
                        int nquad1 = elmtToTrace[n][e]->GetNumPoints(0);
                        int nquad2 = elmtToTrace[n][e]->GetNumPoints(1);

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

                        if (ent.orient == StdRegions::eDir1BwdDir1_Dir2FwdDir2 ||
                            ent.orient == StdRegions::eDir1BwdDir1_Dir2BwdDir2 ||
                            ent.orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1 ||
                            ent.orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
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

                        if (ent.orient == StdRegions::eDir1FwdDir1_Dir2BwdDir2 ||
                            ent.orient == StdRegions::eDir1BwdDir1_Dir2BwdDir2 ||
                            ent.orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1 ||
                            ent.orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
                        {
                            // Reverse y direction
                            for (int j = 0; j < nquad1; ++j)
                            {
                                for (int i = 0; i < nquad2; ++i)
                                {
                                    swap(tmpFwd[offset2 + i*nquad1 + j],
                                         tmpFwd[offset2 + (nquad2-i-1)*nquad1 + j]);
                                }
                            }
                        }

                        for (int i = 0; i < nquad1*nquad2; ++i)
                        {
                            m_periodicFwdCopy.push_back(tmpFwd[i]);
                            m_periodicBwdCopy.push_back(tmpBwd[i]);
                        }
                    }
                }
            }
        }

        /**
         * For each boundary region, checks that the types and number of
         * boundary expansions in that region match.
         * 
         * @param   In          ContField3D to compare with.
         * @return true if boundary conditions match.
         */
        bool DisContField3D::SameTypeOfBoundaryConditions(
            const DisContField3D &In)
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
            int vSame = returnval ? 1 : 0;
            m_comm->AllReduce(vSame, LibUtilities::ReduceMin);

            return (vSame == 1);
        }

        /**
         * According to their boundary region, the separate segmental boundary
         * expansions are bundled together in an object of the class
         * MultiRegions#ExpList2D.
         *
         * \param graph3D A mesh, containing information about the domain and
         * the spectral/hp element expansion.
         * \param bcs An entity containing information about the boundaries and
         * boundary conditions.
         * \param variable An optional parameter to indicate for which variable
         * the boundary conditions should be discretised.
         */
        void DisContField3D::GenerateBoundaryConditionExpansion(
            const SpatialDomains::MeshGraphSharedPtr &graph3D,
            const SpatialDomains::BoundaryConditions &bcs,
            const std::string                        &variable)
        {
            int cnt  = 0;
            MultiRegions::ExpList2DSharedPtr       locExpList;
            SpatialDomains::BoundaryConditionShPtr locBCond;

            const SpatialDomains::BoundaryRegionCollection    &bregions = 
                bcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions = 
                bcs.GetBoundaryConditions();
            SpatialDomains::BoundaryRegionCollection::const_iterator it;

            // count the number of non-periodic boundary regions
            for (it = bregions.begin(); it != bregions.end(); ++it)
            {
                SpatialDomains::BoundaryConditionShPtr boundaryCondition = 
                    GetBoundaryCondition(bconditions, it->first, variable);
                if (boundaryCondition->GetBoundaryConditionType() != 
                        SpatialDomains::ePeriodic)
                {
                    cnt++;
                }
            }

            m_bndCondExpansions = Array<OneD,MultiRegions::ExpListSharedPtr>(cnt);
            m_bndConditions     = Array<OneD,SpatialDomains::BoundaryConditionShPtr>(cnt);

            cnt = 0;

            // list Dirichlet boundaries first
            for (it = bregions.begin(); it != bregions.end(); ++it)
            {
                locBCond = GetBoundaryCondition(
                    bconditions, it->first, variable);

                if(locBCond->GetBoundaryConditionType()
                       != SpatialDomains::ePeriodic)
                {
                    locExpList = MemoryManager<MultiRegions::ExpList2D>
                        ::AllocateSharedPtr(m_session, *(it->second),
                                            graph3D, variable);

                    // Set up normals on non-Dirichlet boundary conditions
                    if(locBCond->GetBoundaryConditionType() != 
                           SpatialDomains::eDirichlet)
                    {
                        SetUpPhysNormals();
                    }

                    m_bndCondExpansions[cnt]  = locExpList;
                    m_bndConditions[cnt++]    = locBCond;
                }
            }
        }

        /**
         * @brief Determine the perioidc edges and vertices for the given graph.
         * 
         * @param   bcs         Information about the boundary conditions.
         * @param   variable    Specifies the field.
         */
        void DisContField3D::FindPeriodicFaces(
            const SpatialDomains::BoundaryConditions &bcs,
            const std::string                        &variable)
        {
            const SpatialDomains::BoundaryRegionCollection &bregions
                = bcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions
                = bcs.GetBoundaryConditions();
            SpatialDomains::MeshGraph3DSharedPtr graph3D
                = boost::dynamic_pointer_cast<
                    SpatialDomains::MeshGraph3D>(m_graph);
            SpatialDomains::BoundaryRegionCollection::const_iterator it;

            LibUtilities::CommSharedPtr     vComm       =
                m_session->GetComm()->GetRowComm();
            LibUtilities::CompositeOrdering compOrder   =
                m_session->GetCompositeOrdering();
            LibUtilities::BndRegionOrdering bndRegOrder =
                m_session->GetBndRegionOrdering();
            SpatialDomains::CompositeMap    compMap     =
                m_graph->GetComposites();

            // Unique collection of pairs of periodic composites (i.e. if
            // composites 1 and 2 are periodic then this map will contain either
            // the pair (1,2) or (2,1) but not both).
            map<int,int>                                   perComps;
            map<int,vector<int> >                          allVerts;
            map<int,SpatialDomains::VertexComponentVector> allCoord;
            map<int,vector<int> >                          allEdges;
            set<int>                                       allFaces;
            set<int>                                       locFaces;
            vector<set<int> >                              locVertEdge(2);

            int region1ID, region2ID, i, j, k, cnt;
            SpatialDomains::BoundaryConditionShPtr locBCond;

            // Construct list of all periodic pairs local to this process.
            for (it = bregions.begin(); it != bregions.end(); ++it)
            {
                locBCond = GetBoundaryCondition(
                    bconditions, it->first, variable);

                if (locBCond->GetBoundaryConditionType()
                        != SpatialDomains::ePeriodic)
                {
                    continue;
                }

                // Identify periodic boundary region IDs.
                region1ID = it->first;
                region2ID = boost::static_pointer_cast<
                    SpatialDomains::PeriodicBoundaryCondition>(
                        locBCond)->m_connectedBoundaryRegion;

                // From this identify composites. Note that in serial this will
                // be an empty map.
                int cId1, cId2;
                if (vComm->GetSize() == 1)
                {
                    cId1 = it->second->begin()->first;
                    cId2 = bregions.find(region2ID)->second->begin()->first;
                }
                else
                {
                    cId1 = bndRegOrder.find(region1ID)->second[0];
                    cId2 = bndRegOrder.find(region2ID)->second[0];
                }

                ASSERTL0(it->second->size() == 1,
                         "Boundary region "+boost::lexical_cast<string>(
                             region1ID)+" should only contain 1 composite.");

                // Construct set containing all periodic edges on this process
                SpatialDomains::Composite c = it->second->begin()->second;

                vector<unsigned int> tmpOrder;

                for (i = 0; i < c->size(); ++i)
                {
                    SpatialDomains::Geometry2DSharedPtr faceGeom =
                        boost::dynamic_pointer_cast<
                            SpatialDomains::Geometry2D>((*c)[i]);
                    ASSERTL0(faceGeom, "Unable to cast to shared ptr");

                    SpatialDomains::ElementFaceVectorSharedPtr elmt =
                        graph3D->GetElementsFromFace(faceGeom);
                    ASSERTL0(elmt->size() == 1,
                             "The periodic boundaries belong to "
                             "more than one element of the mesh");

                    SpatialDomains::Geometry3DSharedPtr geom =
                        boost::dynamic_pointer_cast<SpatialDomains::Geometry3D>(
                            (*elmt)[0]->m_Element);

                    allFaces.insert((*c)[i]->GetGlobalID());

                    // In serial mesh partitioning will not have occurred so
                    // need to fill composite ordering map manually.
                    if (vComm->GetSize() == 1)
                    {
                        tmpOrder.push_back((*c)[i]->GetGlobalID());
                    }

                    vector<int> vertList, edgeList;
                    SpatialDomains::VertexComponentVector tmpVec;
                    for (j = 0; j < faceGeom->GetNumVerts(); ++j)
                    {
                        vertList.push_back(faceGeom->GetVid   (j));
                        edgeList.push_back(faceGeom->GetEid   (j));
                        tmpVec  .push_back(faceGeom->GetVertex(j));
                    }

                    allVerts[(*c)[i]->GetGlobalID()] = vertList;
                    allEdges[(*c)[i]->GetGlobalID()] = edgeList;
                    allCoord[(*c)[i]->GetGlobalID()] = tmpVec;
                }

                if (vComm->GetSize() == 1)
                {
                    compOrder[it->second->begin()->first] = tmpOrder;
                }

                // See if we already have either region1 or region2 stored in
                // perComps map.
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

            // Process local face list to exchange vertices, edges and faces +
            // determine orientations.
            int              n = vComm->GetSize();
            int              p = vComm->GetRank();
            int              totFaces;
            Array<OneD, int> facecounts(n,0);
            Array<OneD, int> vertcounts(n,0);
            Array<OneD, int> faceoffset(n,0);
            Array<OneD, int> vertoffset(n,0);

            facecounts[p] = allFaces.size();
            vComm->AllReduce(facecounts, LibUtilities::ReduceSum);

            faceoffset[0] = 0;
            for (i = 1; i < n; ++i)
            {
                faceoffset[i] = faceoffset[i-1] + facecounts[i-1];
            }

            totFaces = Vmath::Vsum(n, facecounts, 1);
            Array<OneD, int> faceIds   (totFaces, 0);
            Array<OneD, int> faceVerts (totFaces, 0);

            set<int>::iterator sIt;

            for (i = 0, sIt = allFaces.begin(); sIt != allFaces.end(); ++sIt)
            {
                faceIds  [faceoffset[p] + i  ] = *sIt;
                faceVerts[faceoffset[p] + i++] = allVerts[*sIt].size();
            }

            vComm->AllReduce(faceIds,   LibUtilities::ReduceSum);
            vComm->AllReduce(faceVerts, LibUtilities::ReduceSum);

            int nTotVerts = Vmath::Vsum(totFaces, faceVerts, 1);

            // Calculate number of vertices on each processor.
            Array<OneD, int> procVerts(n,0);
            for (i = 0; i < n; ++i)
            {
                procVerts[i] = Vmath::Vsum(
                    facecounts[i], faceVerts + faceoffset[i], 1);
            }

            vertoffset[0] = 0;
            for (i = 1; i < n; ++i)
            {
                vertoffset[i] = vertoffset[i-1] + procVerts[i-1];
            }

            // Exchange x,y,z coordinates of vertices so each process can
            // calculate relative face orientations independently.
            Array<OneD, int>       vertIds(nTotVerts,   0);
            Array<OneD, int>       edgeIds(nTotVerts,   0);
            Array<OneD, NekDouble> vertX  (nTotVerts, 0.0);
            Array<OneD, NekDouble> vertY  (nTotVerts, 0.0);
            Array<OneD, NekDouble> vertZ  (nTotVerts, 0.0);

            for (cnt = 0, sIt = allFaces.begin();
                 sIt != allFaces.end(); ++sIt)
            {
                for (j = 0; j < allVerts[*sIt].size(); ++j)
                {
                    int vertId = allVerts[*sIt][j];
                    vertIds[vertoffset[p] + cnt  ] = vertId;
                    vertX  [vertoffset[p] + cnt  ] = (*allCoord[*sIt][j])(0);
                    vertY  [vertoffset[p] + cnt  ] = (*allCoord[*sIt][j])(1);
                    vertZ  [vertoffset[p] + cnt  ] = (*allCoord[*sIt][j])(2);
                    edgeIds[vertoffset[p] + cnt++] = allEdges[*sIt][j];
                }
            }

            vComm->AllReduce(vertIds, LibUtilities::ReduceSum);
            vComm->AllReduce(edgeIds, LibUtilities::ReduceSum);
            vComm->AllReduce(vertX,   LibUtilities::ReduceSum);
            vComm->AllReduce(vertY,   LibUtilities::ReduceSum);
            vComm->AllReduce(vertZ,   LibUtilities::ReduceSum);

            vector<map<int, vector<int> > >                 vertEdgeMap(2);
            map<int, SpatialDomains::VertexComponentVector> coordMap;

            for (cnt = i = 0; i < totFaces; ++i)
            {
                vector<int> edges(faceVerts[i]);
                vector<int> verts(faceVerts[i]);
                SpatialDomains::VertexComponentVector coord(faceVerts[i]);

                for (j = 0; j < faceVerts[i]; ++j)
                {
                    edges[j] = edgeIds[cnt];
                    verts[j] = vertIds[cnt];
                    coord[j] = MemoryManager<SpatialDomains::VertexComponent>
                        ::AllocateSharedPtr(
                            3, verts[j], vertX[cnt], vertY[cnt], vertZ[cnt++]);
                }

                vertEdgeMap[0][faceIds[i]] = verts;
                vertEdgeMap[1][faceIds[i]] = edges;
                coordMap      [faceIds[i]] = coord;
            }

            // Go through list of composites and figure out which edges are
            // parallel from original ordering in session file. This includes
            // composites which are not necessarily on this process.
            map<int,int>::iterator cIt, pIt;
            map<int,int>::const_iterator oIt;

            // Store temporary map of periodic vertices which will hold all
            // periodic vertices on the entire mesh so that doubly periodic
            // vertices can be counted properly across partitions. Local
            // vertices are copied into m_periodicVerts at the end of the
            // function.
            vector<PeriodicMap> periodicVertEdge(2);

            map<int, map<StdRegions::Orientation, vector<int> > > vmap;
            map<int, map<StdRegions::Orientation, vector<int> > > emap;
            map<int, map<StdRegions::Orientation,
                         vector<StdRegions::Orientation> > > omap;

            map<StdRegions::Orientation, vector<int> > quadVertMap;
            quadVertMap[StdRegions::eDir1FwdDir1_Dir2FwdDir2] += 0,1,2,3;
            quadVertMap[StdRegions::eDir1FwdDir1_Dir2BwdDir2] += 3,2,1,0;
            quadVertMap[StdRegions::eDir1BwdDir1_Dir2FwdDir2] += 1,0,3,2;
            quadVertMap[StdRegions::eDir1BwdDir1_Dir2BwdDir2] += 2,3,0,1;
            quadVertMap[StdRegions::eDir1FwdDir2_Dir2FwdDir1] += 0,3,2,1;
            quadVertMap[StdRegions::eDir1FwdDir2_Dir2BwdDir1] += 3,0,1,2;
            quadVertMap[StdRegions::eDir1BwdDir2_Dir2FwdDir1] += 1,2,3,0;
            quadVertMap[StdRegions::eDir1BwdDir2_Dir2BwdDir1] += 2,1,0,3;

            map<StdRegions::Orientation, vector<int> > quadEdgeMap;
            quadEdgeMap[StdRegions::eDir1FwdDir1_Dir2FwdDir2] += 0,1,2,3;
            quadEdgeMap[StdRegions::eDir1FwdDir1_Dir2BwdDir2] += 2,1,0,3;
            quadEdgeMap[StdRegions::eDir1BwdDir1_Dir2FwdDir2] += 0,3,2,1;
            quadEdgeMap[StdRegions::eDir1BwdDir1_Dir2BwdDir2] += 2,3,0,1;
            quadEdgeMap[StdRegions::eDir1FwdDir2_Dir2FwdDir1] += 3,2,1,0;
            quadEdgeMap[StdRegions::eDir1FwdDir2_Dir2BwdDir1] += 3,0,1,2;
            quadEdgeMap[StdRegions::eDir1BwdDir2_Dir2FwdDir1] += 1,2,3,0;
            quadEdgeMap[StdRegions::eDir1BwdDir2_Dir2BwdDir1] += 1,0,3,2;

            StdRegions::Orientation eF = StdRegions::eForwards;
            StdRegions::Orientation eB = StdRegions::eBackwards;
            
            map<StdRegions::Orientation,
                vector<StdRegions::Orientation> > quadOrientMap;
            quadOrientMap[StdRegions::eDir1FwdDir1_Dir2FwdDir2] += eF,eF,eF,eF;
            quadOrientMap[StdRegions::eDir1FwdDir1_Dir2BwdDir2] += eF,eB,eF,eB;
            quadOrientMap[StdRegions::eDir1BwdDir1_Dir2FwdDir2] += eB,eF,eB,eF;
            quadOrientMap[StdRegions::eDir1BwdDir1_Dir2BwdDir2] += eB,eB,eB,eB;
            quadOrientMap[StdRegions::eDir1FwdDir2_Dir2FwdDir1] += eF,eF,eF,eF;
            quadOrientMap[StdRegions::eDir1FwdDir2_Dir2BwdDir1] += eF,eB,eF,eB;
            quadOrientMap[StdRegions::eDir1BwdDir2_Dir2FwdDir1] += eB,eF,eB,eF;
            quadOrientMap[StdRegions::eDir1BwdDir2_Dir2BwdDir1] += eB,eB,eB,eB;

            map<StdRegions::Orientation, vector<int> > triVertMap;
            triVertMap[StdRegions::eDir1FwdDir1_Dir2FwdDir2] += 0,1,2;
            triVertMap[StdRegions::eDir1BwdDir1_Dir2FwdDir2] += 1,0,2;

            map<StdRegions::Orientation, vector<int> > triEdgeMap;
            triEdgeMap[StdRegions::eDir1FwdDir1_Dir2FwdDir2] += 0,1,2;
            triEdgeMap[StdRegions::eDir1BwdDir1_Dir2FwdDir2] += 0,2,1;

            map<StdRegions::Orientation,
                vector<StdRegions::Orientation> > triOrientMap;
            triOrientMap[StdRegions::eDir1FwdDir1_Dir2FwdDir2] += eF,eF,eF;
            triOrientMap[StdRegions::eDir1BwdDir1_Dir2FwdDir2] += eB,eF,eF;

            vmap[3] = triVertMap;
            vmap[4] = quadVertMap;
            emap[3] = triEdgeMap;
            emap[4] = quadEdgeMap;
            omap[3] = triOrientMap;
            omap[4] = quadOrientMap;

            for (cIt = perComps.begin(); cIt != perComps.end(); ++cIt)
            {
                SpatialDomains::Composite c[2];
                const int   id1  = cIt->first;
                const int   id2  = cIt->second;
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

                ASSERTL0(c[0] || c[1], "Both composites not found on this process!");

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

                // Construct set of all periodic faces that we have locally on
                // this processor.
                for (i = 0; i < 2; ++i)
                {
                    if (!c[i])
                    {
                        continue;
                    }

                    if (c[i]->size() > 0)
                    {
                        for (j = 0; j < c[i]->size(); ++j)
                        {
                            int faceId = c[i]->at(j)->GetGlobalID();
                            locFaces.insert(faceId);
                            for (k = 0; k < vertEdgeMap[0][faceId].size(); ++k)
                            {
                                locVertEdge[0].insert(vertEdgeMap[0][faceId][k]);
                                locVertEdge[1].insert(vertEdgeMap[1][faceId][k]);
                            }
                        }
                    }
                }

                // Loop over all faces in the geometry composite.
                for (pIt = compPairs.begin(); pIt != compPairs.end(); ++pIt)
                {
                    int  ids  [2] = {pIt->first, pIt->second};
                    bool local[2] = {locFaces.count(pIt->first) > 0,
                                     locFaces.count(pIt->second) > 0};

                    ASSERTL0(coordMap.count(ids[0]) > 0 &&
                             coordMap.count(ids[1]) > 0,
                             "Unable to find face in coordinate map");

                    // Calculate relative orientation of faces.
                    SpatialDomains::VertexComponentVector tmpVec1
                        = coordMap[ids[0]];
                    SpatialDomains::VertexComponentVector tmpVec2
                        = coordMap[ids[1]];

                    ASSERTL0(tmpVec1.size() == tmpVec2.size(),
                             "Two periodic faces have different number "
                             "of vertices!");

                    StdRegions::Orientation o;
                    if (tmpVec1.size() == 3)
                    {
                        o = SpatialDomains::TriGeom::GetFaceOrientation(
                            tmpVec1, tmpVec2);
                    }
                    else
                    {
                        o = SpatialDomains::QuadGeom::GetFaceOrientation(
                            tmpVec1, tmpVec2);
                    }

                    for (i = 0; i < 2; ++i)
                    {
                        if (!local[i])
                        {
                            continue;
                        }

                        int other = (i+1) % 2;
                        PeriodicEntity ent(ids  [other], o,
                                           local[other]);
                        m_periodicFaces[ids[i]].push_back(ent);
                    }

                    int nFaceVerts = vertEdgeMap[0][ids[0]].size();

                    for (i = 0; i < 2; ++i)
                    {
                        int other = (i+1) % 2;

                        // Determine periodic vertices.
                        for (j = 0; j < 2; ++j)
                        {
                            vector<int> per1 = vertEdgeMap[j][ids[i]];
                            vector<int> per2 = vertEdgeMap[j][ids[other]];

                            map<int, boost::tuple<
                                int, bool, StdRegions::Orientation> > tmpMap;
                            map<int, boost::tuple<
                                int, bool, StdRegions::Orientation> >::iterator mIt;

                            for (k = 0; k < nFaceVerts; ++k)
                            {
                                int v = j == 0 ? vmap[nFaceVerts][o][k] :
                                    emap[nFaceVerts][o][k];
                                tmpMap[per1[k]] = boost::make_tuple(
                                    per2[v], locVertEdge[j].count(per2[v]) > 0,
                                    j == 0 ? StdRegions::eNoOrientation :
                                       omap[nFaceVerts][o][k]);
                            }

                            for (mIt = tmpMap.begin(); mIt != tmpMap.end(); ++mIt)
                            {
                                // See if this vertex has been recorded already.
                                PeriodicEntity ent2(boost::get<0>(mIt->second),
                                                    boost::get<2>(mIt->second),
                                                    boost::get<1>(mIt->second));
                                PeriodicMap::iterator perIt = periodicVertEdge[j].find(
                                    mIt->first);
                                
                                if (perIt == periodicVertEdge[j].end())
                                {
                                    periodicVertEdge[j][mIt->first].push_back(ent2);
                                    perIt = periodicVertEdge[j].find(mIt->first);
                                }
                                else
                                {
                                    bool doAdd = true;
                                    for (k = 0; k < perIt->second.size(); ++k)
                                    {
                                        if (perIt->second[k].id == boost::get<0>(mIt->second))
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
            }

            // Loop over all periodic vertices to complete connectivity
            // information.
            PeriodicMap::iterator perIt, perIt2;
            for (int l = 0; l < 2; ++l)
            {
                for (perIt  = periodicVertEdge[l].begin();
                     perIt != periodicVertEdge[l].end(); ++perIt)
                {
                    // Loop over associated vertices.
                    for (i = 0; i < perIt->second.size(); ++i)
                    {
                        perIt2 = periodicVertEdge[l].find(perIt->second[i].id);
                        ASSERTL0(perIt2 != periodicVertEdge[l].end(),
                                 "Couldn't find periodic vertex.");
                        
                        for (j = 0; j < perIt2->second.size(); ++j)
                        {
                            if (perIt2->second[j].id == perIt->first)
                            {
                                continue;
                            }
                            
                            bool doAdd = true;
                            for (k = 0; k < perIt->second.size(); ++k)
                            {
                                if (perIt2->second[j].id == perIt->second[k].id)
                                {
                                    doAdd = false;
                                    break;
                                }
                            }
                            
                            if (doAdd)
                            {
                                perIt->second.push_back(perIt2->second[j]);
                            }
                        }
                    }
                }

                // Do one final loop over periodic vertices/edges to remove
                // non-local vertices/edges from map.
                for (perIt  = periodicVertEdge[l].begin();
                     perIt != periodicVertEdge[l].end(); ++perIt)
                {
                    if (locVertEdge[l].count(perIt->first) > 0)
                    {
                        if (l == 0)
                        {
                            m_periodicVerts.insert(*perIt);
                        }
                        else
                        {
                            m_periodicEdges.insert(*perIt);
                        }
                    }
                }
            }

            cout << "Found " << m_periodicFaces.size() << " periodic faces" << endl;
            
            PeriodicMap::iterator asd;
            for (asd = m_periodicFaces.begin(); asd != m_periodicFaces.end(); ++asd)
            {
                cout << asd->first << " <-> " << asd->second[0].id
                     << " " << asd->second[0].orient << " " << asd->second[0].isLocal << endl;
            }

            cout << "RANK: " << vComm->GetRank() << " found " << m_periodicVerts.size() << " periodic verts" << endl;
            for (asd = m_periodicVerts.begin(); asd != m_periodicVerts.end(); ++asd)
            {
                cout << "RANK: " << vComm->GetRank() << " " << asd->first << " <-> ";
                for (i = 0; i < asd->second.size(); ++i)
                {
                    cout << "(" << asd->second[i].id
                         << "," << asd->second[i].isLocal
                         << ") ";
                }
                cout << endl;
            }
            
            cout << "RANK: " << vComm->GetRank() << " found " << m_periodicEdges.size() << " periodic edges" << endl;
            for (asd = m_periodicEdges.begin(); asd != m_periodicEdges.end(); ++asd)
            {
                cout << "RANK: " << vComm->GetRank() << " " << asd->first << " <-> ";
                for (i = 0; i < asd->second.size(); ++i)
                {
                    cout << "(" << asd->second[i].id
                         << "," << asd->second[i].isLocal
                         << ") ";
                }
                cout << endl;
            }
        }

        bool DisContField3D::IsLeftAdjacentFace(const int n, const int e)
        {
            set<int>::iterator it;
            LocalRegions::Expansion2DSharedPtr traceEl = 
                boost::dynamic_pointer_cast<LocalRegions::Expansion2D>(
                    (m_traceMap->GetElmtToTrace())[n][e]);
            
            int offset = m_trace->GetPhys_Offset(traceEl->GetElmtId());
            
            bool fwd = true;
            if (traceEl->GetLeftAdjacentElementFace () == -1 ||
                traceEl->GetRightAdjacentElementFace() == -1)
            {
                // Boundary edge (1 connected element). Do nothing in
                // serial.
                it = m_boundaryFaces.find(traceEl->GetElmtId());
                
                // If the edge does not have a boundary condition set on
                // it, then assume it is a partition edge.
                if (it == m_boundaryFaces.end())
                {
                    fwd = m_traceMap->
                        GetTraceToUniversalMapUnique(offset) > 0;
                }
            }
            else if (traceEl->GetLeftAdjacentElementFace () != -1 &&
                     traceEl->GetRightAdjacentElementFace() != -1)
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

        /**
         * \brief This method extracts the "forward" and "backward" trace
         * data from the array \a field and puts the data into output
         * vectors \a Fwd and \a Bwd.
         * 
         * We first define the convention which defines "forwards" and
         * "backwards". First an association is made between the face of
         * each element and its corresponding face in the trace space
         * using the mapping #m_traceMap. The element can either be
         * left-adjacent or right-adjacent to this trace face (see
         * Expansion2D::GetLeftAdjacentElementExp). Boundary faces are
         * always left-adjacent since left-adjacency is populated first.
         * 
         * If the element is left-adjacent we extract the face trace data
         * from \a field into the forward trace space \a Fwd; otherwise,
         * we place it in the backwards trace space \a Bwd. In this way,
         * we form a unique set of trace normals since these are always
         * extracted from left-adjacent elements.
         *
         * \param field is a NekDouble array which contains the 3D data
         * from which we wish to extract the backward and forward
         * orientated trace/face arrays.
         *
         * \return Updates a NekDouble array \a Fwd and \a Bwd
         */
        void DisContField3D::v_GetFwdBwdTracePhys(Array<OneD, NekDouble> &Fwd,
                                                  Array<OneD, NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhys(m_phys, Fwd, Bwd);
        }
        
        
        void DisContField3D::v_GetFwdBwdTracePhys(
            const Array<OneD, const NekDouble> &field,
                  Array<OneD,       NekDouble> &Fwd,
                  Array<OneD,       NekDouble> &Bwd)
        {
            // Loop over elements and collect forward and backward expansions.
            int nexp = GetExpSize();
            int cnt, n, e, npts, offset, phys_offset;
            Array<OneD,NekDouble> e_tmp;
            
            Array<OneD, Array<OneD, StdRegions::StdExpansionSharedPtr> >
                &elmtToTrace = m_traceMap->GetElmtToTrace();

            set<int>::iterator    it;
            PeriodicMap::iterator it2;
            boost::unordered_map<int,pair<int,int> >::iterator it3;
            
            // Zero vectors.
            Vmath::Zero(Fwd.num_elements(), Fwd, 1);
            Vmath::Zero(Bwd.num_elements(), Bwd, 1);

            bool fwd;

            for(cnt = n = 0; n < nexp; ++n)
            {
                phys_offset = GetPhys_Offset(n);
                for(e = 0; e < (*m_exp)[n]->GetNfaces(); ++e, ++cnt)
                {
                    offset = m_trace->GetPhys_Offset(
                        elmtToTrace[n][e]->GetElmtId());

                    fwd = m_leftAdjacentFaces[cnt];
                    if (fwd)
                    {
                        (*m_exp)[n]->GetFacePhysVals(e, elmtToTrace[n][e],
                                                     field + phys_offset,
                                                     e_tmp = Fwd + offset);
                    }
                    else
                    {
                        (*m_exp)[n]->GetFacePhysVals(e, elmtToTrace[n][e],
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
                if(m_bndConditions[n]->GetBoundaryConditionType() == 
                       SpatialDomains::eDirichlet)
                {
                    for(e = 0; e < m_bndCondExpansions[n]->GetExpSize(); ++e)
                    {
                        npts = m_bndCondExpansions[n]->GetExp(e)->GetTotPoints();
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
                        npts = m_bndCondExpansions[n]->GetExp(e)->GetNumPoints(0)*
                               m_bndCondExpansions[n]->GetExp(e)->GetNumPoints(1);
                        id1  = m_bndCondExpansions[n]->GetPhys_Offset(e);
                        id2  = m_trace->GetPhys_Offset(
                            m_traceMap->GetBndCondTraceToGlobalTraceMap(cnt+e));

                        ASSERTL0((m_bndCondExpansions[n]->GetPhys())[id1] == 0.0,
                                 "method not set up for non-zero Neumann "
                                 "boundary condition");
                        
                        Vmath::Vcopy(npts,&Fwd[id2],1,&Bwd[id2],1);
                    }

                    cnt += e;
                }
                else
                {
                    ASSERTL0(false, "Method only set up for Dirichlet, Neumann "
                             "and Robin conditions.");
                }
            }

            // Copy any periodic boundary conditions.
            for (n = 0; n < m_periodicFwdCopy.size(); ++n)
            {
                Bwd[m_periodicBwdCopy[n]] = Fwd[m_periodicFwdCopy[n]];
            }
            
            // Do parallel exchange for forwards/backwards spaces.
            m_traceMap->UniversalTraceAssemble(Fwd);
            m_traceMap->UniversalTraceAssemble(Bwd);
        }

        void DisContField3D::v_ExtractTracePhys(
            Array<OneD, NekDouble> &outarray)
        {
            ASSERTL1(m_physState == true,
                     "Field is not in physical space.");
            
            v_ExtractTracePhys(m_phys, outarray);
        }

        void DisContField3D::v_ExtractTracePhys(
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
            for(n = 0; n < nexp; ++n)
            {
                phys_offset = GetPhys_Offset(n);
                
                for(e = 0; e < (*m_exp)[n]->GetNfaces(); ++e)
                {
                    offset = m_trace->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());
                    (*m_exp)[n]->GetFacePhysVals(e, elmtToTrace[n][e],
                                                 inarray + phys_offset,
                                                 e_tmp = outarray + offset);
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
        void DisContField3D::v_AddTraceIntegral(
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
                e_outarray = outarray+offset;
                for(e = 0; e < (*m_exp)[n]->GetNfaces(); ++e)
                {
                    t_offset = m_trace->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());
                    (*m_exp)[n]->AddFaceNormBoundaryInt(e,elmtToTrace[n][e],
                                                        Fn + t_offset,
                                                        e_outarray);
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
         * and adds this to the coefficient space provided by
         * outarray. The value of q is determined from the routine
         * IsLeftAdjacentFace() which if true we use Fwd else we use
         * Bwd
         *
         * @see Expansion3D::AddFaceNormBoundaryInt
         * 
         * @param Fwd       The trace quantities associated with left (fwd)
         *                  adjancent elmt.
         * @param Bwd       The trace quantities associated with right (bwd)
         *                  adjacent elet.
         * @param outarray  Resulting 3D coefficient space.
         */
        void DisContField3D::v_AddFwdBwdTraceIntegral(
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
                for(e = 0; e < (*m_exp)[n]->GetNfaces(); ++e)
                {
                    t_offset = m_trace->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());

                    // Evaluate upwind flux less local edge 
                    if(IsLeftAdjacentFace(n,e))
                    {
                        (*m_exp)[n]->AddFaceNormBoundaryInt(
                         e, elmtToTrace[n][e],  Fwd + t_offset,
                         e_outarray = outarray+offset);
                    }
                    else
                    {
                        (*m_exp)[n]->AddFaceNormBoundaryInt(
                         e, elmtToTrace[n][e],  Bwd + t_offset,
                         e_outarray = outarray+offset);
                    }
                }
            }
        }

        /**
         * @brief Set up a list of elemeent IDs and edge IDs that link to the
         * boundary conditions.
         */
        void DisContField3D::v_GetBoundaryToElmtMap(
            Array<OneD, int> &ElmtID,
            Array<OneD, int> &FaceID)
        {
            map<int,int> globalIdMap;
            int i, n;
            int cnt;
            int nbcs = 0;
            
            SpatialDomains::MeshGraph3DSharedPtr graph3D = 
                boost::dynamic_pointer_cast<SpatialDomains::MeshGraph3D>(
                    m_graph);
            
            // Populate global ID map (takes global geometry ID to local
            // expansion list ID).
            for (i = 0; i < GetExpSize(); ++i)
            {
                globalIdMap[(*m_exp)[i]->GetGeom3D()->GetGlobalID()] = i;
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

            if(FaceID.num_elements() != nbcs)
            {
                FaceID = Array<OneD, int>(nbcs);
            }
            
            for(cnt = n = 0; n < m_bndCondExpansions.num_elements(); ++n)
            {
                for(i = 0; i < m_bndCondExpansions[n]->GetExpSize(); ++i, ++cnt)
                {
                    // Use face to element map from MeshGraph3D.
                    SpatialDomains::ElementFaceVectorSharedPtr tmp = 
                        graph3D->GetElementsFromFace(
                            m_bndCondExpansions[n]->GetExp(i)->GetGeom2D());
                    
                    ElmtID[cnt] = globalIdMap[(*tmp)[0]->
                                              m_Element->GetGlobalID()];
                    FaceID[cnt] = (*tmp)[0]->m_FaceIndx;
                }
            }
        }

        /**
         * Solving Helmholtz Equation in 3D
         */
        void DisContField3D::v_HelmSolve(
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
            GlobalMatrixKey HDGLamToUKey(StdRegions::eHybridDGLamToU,NullAssemblyMapSharedPtr,factors,varcoeff);
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
            GlobalMatrixKey invHDGhelmkey(StdRegions::eInvHybridDGHelmholtz,NullAssemblyMapSharedPtr,factors,varcoeff);
            const DNekScalBlkMatSharedPtr& InvHDGHelm = GetBlockMatrix(invHDGhelmkey);
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
        void DisContField3D::v_GeneralMatrixOp(
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
         * Search through the edge expansions and identify which ones have
         * Robin/Mixed type boundary conditions. If find a Robin boundary then
         * store the edge id of the boundary condition and the array of points
         * of the physical space boundary condition which are hold the boundary
         * condition primitive variable coefficient at the quatrature points
         *
         * \return std map containing the robin boundary condition
         * info using a key of the element id
         *
         * There is a next member to allow for more than one Robin
         * boundary condition per element
         */
        map<int, RobinBCInfoSharedPtr> DisContField3D::v_GetRobinBCInfo(void)
        {
            int i,cnt;
            map<int, RobinBCInfoSharedPtr> returnval;
            Array<OneD, int> ElmtID,FaceID;
            GetBoundaryToElmtMap(ElmtID,FaceID);

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
                        RobinBCInfoSharedPtr rInfo = MemoryManager<RobinBCInfo>::AllocateSharedPtr(FaceID[cnt+e],Array_tmp = locExpList->GetPhys() + locExpList->GetPhys_Offset(e));
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
         * \brief This function evaluates the boundary conditions at a certain
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
         */
        void DisContField3D::v_EvaluateBoundaryConditions(const NekDouble time,
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
                    npoints = locExpList->GetNpoints();
                    
                    Array<OneD,NekDouble> x0(npoints,0.0);
                    Array<OneD,NekDouble> x1(npoints,0.0);
                    Array<OneD,NekDouble> x2(npoints,0.0);
                    
                    locExpList->GetCoords(x0,x1,x2);
                    
                    if(m_bndConditions[i]->GetBoundaryConditionType()
                       == SpatialDomains::eDirichlet)
                    {
                    LibUtilities::Equation  condition = boost::static_pointer_cast<
                        SpatialDomains::DirichletBoundaryCondition >(m_bndConditions[i])->m_dirichletCondition;
                    
                    condition.Evaluate(x0,x1,x2,time,locExpList->UpdatePhys());
                    
                    locExpList->FwdTrans_BndConstrained(locExpList->GetPhys(),
                                                        locExpList->UpdateCoeffs());
                    }
                    else if(m_bndConditions[i]->GetBoundaryConditionType()
                            == SpatialDomains::eNeumann)
                    {
                        LibUtilities::Equation  condition = boost::static_pointer_cast<
                        SpatialDomains::NeumannBoundaryCondition
                            >(m_bndConditions[i])->m_neumannCondition;
                        
                        condition.Evaluate(x0,x1,x2,time,locExpList->UpdatePhys());
                        
                        locExpList->IProductWRTBase(locExpList->GetPhys(),
                                                    locExpList->UpdateCoeffs());
                    }
                    else if(m_bndConditions[i]->GetBoundaryConditionType()
                            == SpatialDomains::eRobin)
                    {
                        LibUtilities::Equation  condition = boost::static_pointer_cast<
                        SpatialDomains::RobinBoundaryCondition
                            >(m_bndConditions[i])->m_robinFunction;
                        
                        LibUtilities::Equation coeff     = 
                            boost::static_pointer_cast<
                        SpatialDomains::RobinBoundaryCondition
                            >(m_bndConditions[i])->m_robinPrimitiveCoeff;
                        
                        condition.Evaluate(x0,x1,x2,time,locExpList->UpdatePhys());
                        
                        locExpList->IProductWRTBase(locExpList->GetPhys(),
                                                    locExpList->UpdateCoeffs());
                        
                        // put primitive coefficient into the physical space
                        // storage
                        coeff.Evaluate(x0,x1,x2,time,
                                       locExpList->UpdatePhys());
                        
                    }
                    else
                    {
                        ASSERTL0(false,"This type of BC not implemented yet");
                    }
                }
            }
        }
    } // end of namespace
} // end of namespace
