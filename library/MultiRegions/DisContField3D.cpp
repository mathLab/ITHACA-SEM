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

#include <boost/core/ignore_unused.hpp>

#include <MultiRegions/DisContField3D.h>
#include <LocalRegions/Expansion3D.h>
#include <LocalRegions/Expansion2D.h>
#include <SpatialDomains/MeshGraph.h>
#include <LocalRegions/HexExp.h>
#include <LocalRegions/TetExp.h>
#include <LocalRegions/PrismExp.h>
#include <LibUtilities/Foundations/Interp.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <tuple>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/predicate.hpp>

using namespace std;

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
             m_bndCondBndWeight    (),
             m_bndConditions       (),
             m_trace(NullExpListSharedPtr)
         {
         }

         /**
          * @brief Constructs a global discontinuous field based on an input
          * mesh with boundary conditions.
          */
         DisContField3D::DisContField3D(
             const LibUtilities::SessionReaderSharedPtr &pSession,
             const SpatialDomains::MeshGraphSharedPtr   &graph3D,
             const std::string                          &variable,
             const bool                                  SetUpJustDG,
             const Collections::ImplementationType       ImpType):
             ExpList3D          (pSession, graph3D, variable, ImpType),
               m_bndCondExpansions(),
               m_bndCondBndWeight(),
               m_bndConditions    (),
               m_trace(NullExpListSharedPtr)
         {
             // do not set up BCs if default variable
             if (variable.compare("DefaultVar") != 0)
             {
                 SpatialDomains::BoundaryConditions bcs(m_session, graph3D);

                 GenerateBoundaryConditionExpansion(graph3D,bcs,variable);
                 EvaluateBoundaryConditions(0.0, variable);

                 // Find periodic edges for this variable.
                 FindPeriodicFaces(bcs, variable);
             }

             if(SetUpJustDG)
             {
                 SetUpDG();
                 m_locTraceToTraceMap->
                        TraceLocToElmtLocCoeffMap(*this, m_trace);
             }
             else
             {
                 // Set element edges to point to Robin BC edges if required.
                 int i, cnt, f;
                 Array<OneD, int> ElmtID, FaceID;
                 GetBoundaryToElmtMap(ElmtID, FaceID);

                 for(cnt = i = 0; i < m_bndCondExpansions.size(); ++i)
                 {
                     MultiRegions::ExpListSharedPtr locExpList;
                     locExpList = m_bndCondExpansions[i];

                     for(f = 0; f < locExpList->GetExpSize(); ++f)
                     {
                         LocalRegions::Expansion3DSharedPtr exp3d
                            = (*m_exp)[ElmtID[cnt+f]]->
                                as<LocalRegions::Expansion3D>();
                         LocalRegions::Expansion2DSharedPtr exp2d
                            = locExpList->GetExp(f)->
                                as<LocalRegions::Expansion2D>();

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
             EvaluateBoundaryConditions(0.0, variable);
             ApplyGeomInfo();

             if (!SameTypeOfBoundaryConditions(In))
             {
                 // Find periodic edges for this variable.
                 FindPeriodicFaces(bcs, variable);

                 if (SetUpJustDG)
                 {
                     SetUpDG(variable);
                    m_locTraceToTraceMap->
                            TraceLocToElmtLocCoeffMap(*this, m_trace);
                 }
                 else
                 {
                     int i,cnt,f;
                     Array<OneD, int> ElmtID,FaceID;
                     GetBoundaryToElmtMap(ElmtID,FaceID);

                     for(cnt = i = 0; i < m_bndCondExpansions.size(); ++i)
                     {
                         MultiRegions::ExpListSharedPtr locExpList;
                         locExpList = m_bndCondExpansions[i];

                         for(f = 0; f < locExpList->GetExpSize(); ++f)
                         {
                             LocalRegions::Expansion3DSharedPtr exp3d
                                 = (*m_exp)[ElmtID[cnt+f]]->
                                         as<LocalRegions::Expansion3D>();
                             LocalRegions::Expansion2DSharedPtr exp2d
                                 = locExpList->GetExp(f)->
                                         as<LocalRegions::Expansion2D>();

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
                 m_globalBndMat       = In.m_globalBndMat;
                 m_trace              = In.m_trace;
                 m_traceMap           = In.m_traceMap;
                 m_locTraceToTraceMap = In.m_locTraceToTraceMap;
                 m_periodicVerts      = In.m_periodicVerts;
                 m_periodicEdges      = In.m_periodicEdges;
                 m_periodicFaces      = In.m_periodicFaces;

                 if(SetUpJustDG)
                 {
                 }
                 else
                 {
                     int i,cnt,f;
                     Array<OneD, int> ElmtID,FaceID;
                     GetBoundaryToElmtMap(ElmtID,FaceID);

                     for (cnt = i = 0;
                          i < m_bndCondExpansions.size(); ++i)
                     {
                         MultiRegions::ExpListSharedPtr locExpList;
                         locExpList = m_bndCondExpansions[i];

                         for(f = 0; f < locExpList->GetExpSize(); ++f)
                         {
                             LocalRegions::Expansion3DSharedPtr exp3d
                                 = (*m_exp)[ElmtID[cnt+f]]->
                                         as<LocalRegions::Expansion3D>();
                             LocalRegions::Expansion2DSharedPtr exp2d
                                 = locExpList->GetExp(f)->
                                         as<LocalRegions::Expansion2D>();

                             exp3d->SetFaceExp(FaceID[cnt+f], exp2d);
                             exp2d->SetAdjacentElementExp(FaceID[cnt+f], exp3d);
                         }

                         cnt += m_bndCondExpansions[i]->GetExpSize();
                     }

                     if (m_session->DefinesSolverInfo("PROJECTION"))
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
             m_traceMap            (In.m_traceMap),
             m_locTraceToTraceMap  (In.m_locTraceToTraceMap),
             m_periodicFaces       (In.m_periodicFaces),
             m_periodicEdges       (In.m_periodicEdges),
             m_periodicVerts       (In.m_periodicVerts)
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
             auto matrixIter = m_globalBndMat->find(mkey);

             if (matrixIter == m_globalBndMat->end())
             {
                 glo_matrix = GenGlobalBndLinSys(mkey, m_traceMap);
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

             // Set up matrix map
             m_globalBndMat = MemoryManager<GlobalLinSysMap>::
                 AllocateSharedPtr();

             // Set up Trace space
             bool UseGenSegExp = true;
             trace = MemoryManager<ExpList2D>::AllocateSharedPtr(
                 m_session, m_bndCondExpansions, m_bndConditions,
                 *m_exp, m_graph, m_periodicFaces, UseGenSegExp);

             m_trace    = trace;

             m_traceMap = MemoryManager<AssemblyMapDG>::AllocateSharedPtr(
                 m_session, m_graph, trace, *this, m_bndCondExpansions,
                 m_bndConditions, m_periodicFaces,variable);

             if (m_session->DefinesCmdLineArgument("verbose"))
             {
                 m_traceMap->PrintStats(std::cout, variable);
             }

             Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >
                 &elmtToTrace = m_traceMap->GetElmtToTrace();

             // Scatter trace segments to 3D elements. For each element, we
             // find the trace segment associated to each edge. The element
             // then retains a pointer to the trace space segments, to ensure
             // uniqueness of normals when retrieving from two adjoining
             // elements which do not lie in a plane.
             for (int i = 0; i < m_exp->size(); ++i)
             {
                 for (int j = 0; j < (*m_exp)[i]->GetNfaces(); ++j)
                 {
                     LocalRegions::Expansion3DSharedPtr exp3d =
                             (*m_exp)[i]->as<LocalRegions::Expansion3D>();
                     LocalRegions::Expansion2DSharedPtr exp2d =
                             elmtToTrace[i][j]->as<LocalRegions::Expansion2D>();
                     exp3d->SetFaceExp           (j, exp2d);
                     exp2d->SetAdjacentElementExp(j, exp3d);
                 }
             }

             // Set up physical normals
             SetUpPhysNormals();

             int cnt, n, e;

             // Identify boundary faces
             for(cnt = 0, n = 0; n < m_bndCondExpansions.size(); ++n)
             {
                 if (m_bndConditions[n]->GetBoundaryConditionType() !=
                     SpatialDomains::ePeriodic)
                 {
                     for(e = 0; e < m_bndCondExpansions[n]->GetExpSize(); ++e)
                     {
                         m_boundaryFaces.insert(
                             m_traceMap->GetBndCondIDToGlobalTraceID(cnt+e));
                     }
                     cnt += m_bndCondExpansions[n]->GetExpSize();
                 }
             }

             // Set up information for periodic boundary conditions.
             std::unordered_map<int,pair<int,int> > perFaceToExpMap;
             cnt = 0;
             LocalRegions::Expansion3DSharedPtr exp3d;
             for (int n = 0; n < m_exp->size(); ++n)
             {
                 exp3d = (*m_exp)[n]->as<LocalRegions::Expansion3D>();
                 for (int e = 0; e < exp3d->GetNfaces(); ++e, ++cnt)
                 {
                     auto it = m_periodicFaces.find(
                         exp3d->GetGeom3D()->GetFid(e));

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
                 exp3d = (*m_exp)[n]->as<LocalRegions::Expansion3D>();
                 for (int e = 0; e < exp3d->GetNfaces(); ++e, ++cnt)
                 {
                     int faceGeomId = exp3d->GetGeom3D()->GetFid(e);
                     int offset = m_trace->GetPhys_Offset(
                         elmtToTrace[n][e]->GetElmtId());

                     // Check to see if this face is periodic.
                     auto it = m_periodicFaces.find(faceGeomId);

                     if (it != m_periodicFaces.end())
                     {
                         const PeriodicEntity &ent = it->second[0];
                         auto it2 = perFaceToExpMap.find(ent.id);

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
                            ent.orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1 ||
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
                            ent.orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1 ||
                            ent.orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
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
                }
            }

            m_locTraceToTraceMap = MemoryManager<LocTraceToTraceMap>::
                AllocateSharedPtr(*this, m_trace, elmtToTrace,
                                  m_leftAdjacentFaces);

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
            int cnt = 0;
            MultiRegions::ExpList2DSharedPtr       locExpList;
            SpatialDomains::BoundaryConditionShPtr locBCond;

            const SpatialDomains::BoundaryRegionCollection    &bregions =
                bcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions =
                bcs.GetBoundaryConditions();

            m_bndCondExpansions =
                Array<OneD,MultiRegions::ExpListSharedPtr>(bregions.size());
            m_bndConditions     =
                Array<OneD,SpatialDomains::BoundaryConditionShPtr>(bregions.size());

            m_bndCondBndWeight = Array<OneD, NekDouble> {bregions.size(), 0.0};
            
            // list Dirichlet boundaries first
            for (auto &it : bregions)
            {
                locBCond = GetBoundaryCondition(
                    bconditions, it.first, variable);
                locExpList = MemoryManager<MultiRegions::ExpList2D>
                    ::AllocateSharedPtr(m_session, *(it.second),
                                        graph3D, variable, locBCond->GetComm());

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

        /**
         * @brief Determine the periodic faces, edges and vertices for the given
         * graph.
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

            LibUtilities::CommSharedPtr       vComm       =
                m_session->GetComm()->GetRowComm();
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

                for(j = 0; j < (*m_exp)[i]->GetNedges(); ++j)
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

                // Check the region only contains a single composite.
                ASSERTL0(it.second->size() == 1,
                         "Boundary region "+boost::lexical_cast<string>(
                             region1ID)+" should only contain 1 composite.");

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
                    procVerts[i] = Vmath::Vsum(
                        facecounts[i], faceVerts + faceoffset[i], 1);
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
                        ::AllocateSharedPtr(
                            3, verts[j], vertX[cnt], vertY[cnt], vertZ[cnt]);
                    vCoMap[vertIds[cnt]] = coord[j];

                    // Try to insert edge into the eIdMap to avoid re-inserting.
                    auto testIns = eIdMap.insert(
                        make_pair(
                            edgeIds[cnt],
                            make_pair(vertIds[tmp+j],
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

            // Collect composite ides of each periodic face for use if rotation is required
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
                    bool rotbnd     = false;
                    int  dir        = 0;
                    NekDouble angle = 0.0;
                    NekDouble sign  = 1.0;
                    NekDouble tol   = 1e-8;

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
                            o = SpatialDomains::TriGeom::GetFaceOrientation(
                                   tmpVec[i], tmpVec[other],
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
                            ASSERTL0(
                                o == StdRegions::eDir1FwdDir1_Dir2FwdDir2 ||
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
                            tmpMap[per1[j]] = make_pair(
                                per2[v], locVerts.count(per2[v]) > 0);
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
                            tmpMap[per1[j]] = make_pair(
                                per2[e], locEdges.count(per2[e]) > 0);
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

            Array<OneD, int> pairSizes(n, 0);
            pairSizes[p] = allCompPairs.size();
            vComm->AllReduce(pairSizes, LibUtilities::ReduceSum);

            int totPairSizes = Vmath::Vsum(n, pairSizes, 1);

            Array<OneD, int> pairOffsets(n, 0);
            pairOffsets[0] = 0;

            for (i = 1; i < n; ++i)
            {
                pairOffsets[i] = pairOffsets[i-1] + pairSizes[i-1];
            }


            ASSERTL1(allCompPairs.size() == fIdToCompId.size(),
                     "At this point the size of allCompPairs "
                     "should have been the same as fIdToCompId");

            Array<OneD, int> first (totPairSizes, 0);
            Array<OneD, int> second(totPairSizes, 0);

            cnt = pairOffsets[p];

            for (auto &pIt : allCompPairs)
            {
                first [cnt  ] = pIt.first;
                second[cnt++] = pIt.second;
            }

            vComm->AllReduce(first,  LibUtilities::ReduceSum);
            vComm->AllReduce(second, LibUtilities::ReduceSum);

            allCompPairs.clear();

            for(cnt = 0; cnt < totPairSizes; ++cnt)
            {
                allCompPairs[first[cnt]] = second[cnt];
            }

            // make global list of faces to composite ids if rotComp is non-zero

            if(rotComp.size())
            {
                Vmath::Zero(totPairSizes,first,1);
                Vmath::Zero(totPairSizes,second,1);

                cnt = pairOffsets[p];

                for (auto &pIt : fIdToCompId)
                {
                    first [cnt  ] = pIt.first;
                    second[cnt++] = pIt.second;
                }

                vComm->AllReduce(first,  LibUtilities::ReduceSum);
                vComm->AllReduce(second, LibUtilities::ReduceSum);

                fIdToCompId.clear();

                for(cnt = 0; cnt < totPairSizes; ++cnt)
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
                bool rotbnd     = false;
                int dir         = 0;
                NekDouble angle = 0.0;
                NekDouble tol   = 1e-8;

                int faceId    = faceIds[i];

                ASSERTL0(allCompPairs.count(faceId) > 0,
                         "Unable to find matching periodic face.");

                int perFaceId = allCompPairs[faceId];

                // check to see if periodic boundary is rotated
                ASSERTL1(fIdToCompId.count(faceId) > 0,"Face " +
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
                bool rotbnd     = false;
                int dir         = 0;
                NekDouble angle = 0.0;
                NekDouble tol   = 1e-8;


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
                                ASSERTL0(false,"Unable to align rotationally periodic edge vertex");
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

        bool DisContField3D::IsLeftAdjacentFace(const int n, const int e)
        {
            LocalRegions::Expansion2DSharedPtr traceEl =
                    m_traceMap->GetElmtToTrace()[n][e]->
                         as<LocalRegions::Expansion2D>();

            bool fwd = true;
            if (traceEl->GetLeftAdjacentElementFace () == -1 ||
                traceEl->GetRightAdjacentElementFace() == -1)
            {
                // Boundary edge (1 connected element). Do nothing in
                // serial.
                auto it = m_boundaryFaces.find(traceEl->GetElmtId());

                // If the edge does not have a boundary condition set on
                // it, then assume it is a partition edge.
                if (it == m_boundaryFaces.end())
                {
                    fwd = true; // Partition edge is always fwd
                }
            }
            else if (traceEl->GetLeftAdjacentElementFace () != -1 &&
                     traceEl->GetRightAdjacentElementFace() != -1)
            {
                // Non-boundary edge (2 connected elements).
                fwd = (traceEl->GetLeftAdjacentElementExp().get() == (*m_exp)[n].get() );
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
        void DisContField3D::v_GetFwdBwdTracePhysInterior(
            const Array<OneD, const NekDouble> &field,
                  Array<OneD,       NekDouble> &Fwd,
                  Array<OneD,       NekDouble> &Bwd)
        {
            // Zero vectors.
            Vmath::Zero(Fwd.size(), Fwd, 1);
            Vmath::Zero(Bwd.size(), Bwd, 1);

            Array<OneD, NekDouble> facevals(m_locTraceToTraceMap->
                                            GetNLocTracePts());
            m_locTraceToTraceMap->LocTracesFromField(field,facevals);
            m_locTraceToTraceMap->InterpLocFacesToTrace(0, facevals, Fwd);

            Array<OneD, NekDouble> invals = facevals + m_locTraceToTraceMap->
                                                        GetNFwdLocTracePts();
            m_locTraceToTraceMap->InterpLocFacesToTrace(1, invals, Bwd);

            DisContField3D::v_PeriodicBwdCopy(Fwd, Bwd);
        }

        void DisContField3D::v_AddTraceQuadPhysToField(
            const Array<OneD, const NekDouble> &Fwd,
            const Array<OneD, const NekDouble> &Bwd,
                  Array<OneD,       NekDouble> &field)
        {
            Array<OneD, NekDouble> facevals(m_locTraceToTraceMap->
                                            GetNLocTracePts(), 0.0);

            Array<OneD, NekDouble> invals = facevals + m_locTraceToTraceMap->
                                                    GetNFwdLocTracePts();
            m_locTraceToTraceMap->
                    RightIPTWLocFacesToTraceInterpMat(1, Bwd, invals);
            
            m_locTraceToTraceMap->
                    RightIPTWLocFacesToTraceInterpMat(0, Fwd, facevals);

            m_locTraceToTraceMap->AddLocTracesToField(facevals, field);
        }

        /**
         * @brief Fill the Bwd based on corresponding boundary conditions.
         * Periodic boundary is considered interior traces and is not 
         * treated here.
         */
        void DisContField3D::v_FillBwdWithBound(
            const Array<OneD, const NekDouble> &Fwd,
                  Array<OneD,       NekDouble> &Bwd)
        {
            int cnt = 0;
            int npts = 0;
            int e = 0;
            int id1 = 0;
            int id2 = 0;

            for (int n = 0; n < m_bndCondExpansions.size(); ++n)
            {
                if(m_bndConditions[n]->GetBoundaryConditionType() ==
                       SpatialDomains::eDirichlet)
                {
                    for(e = 0; e < m_bndCondExpansions[n]->GetExpSize(); ++e)
                    {
                        npts = m_bndCondExpansions[n]->GetExp(e)->GetTotPoints();
                        id1  = m_bndCondExpansions[n]->GetPhys_Offset(e);
                        id2  = m_trace->GetPhys_Offset(
                            m_traceMap->GetBndCondIDToGlobalTraceID(cnt+e));
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
                        npts = m_bndCondExpansions[n]->GetExp(e)->GetTotPoints();
                        id1  = m_bndCondExpansions[n]->GetPhys_Offset(e);
                        id2  = m_trace->GetPhys_Offset(
                            m_traceMap->GetBndCondIDToGlobalTraceID(cnt+e));

                        // Turning this off since we can have non-zero
                        //Neumann in mixed CG-DG method
                        //ASSERTL1((m_bndCondExpansions[n]->GetPhys())[id1]
                        //== 0.0, "method not set up for non-zero
                        //Neumann " "boundary condition");

                        Vmath::Vcopy(npts,&Fwd[id2],1,&Bwd[id2],1);
                    }

                    cnt += e;
                }
                else if (m_bndConditions[n]->GetBoundaryConditionType() ==
                             SpatialDomains::ePeriodic)
                {
                    continue;
                }
                else
                {
                    ASSERTL0(false, "Method only set up for Dirichlet, Neumann "
                             "and Robin conditions.");
                }
            }
        }

        /**
         * @brief Fill the Bwd based on corresponding boundary conditions for 
         * derivatives. Periodic boundary is considered interior traces and 
         * is not treated here.
         */
        void DisContField3D::v_FillBwdWithBoundDeriv(
            const int                          Dir,
            const Array<OneD, const NekDouble> &Fwd,
                  Array<OneD,       NekDouble> &Bwd)
        {
            boost::ignore_unused(Dir);
            int cnt = 0;
            int e = 0;
            int npts = 0;
            int id2 = 0;
            
            for(int n = 0; n < m_bndCondExpansions.size(); ++n)
            {
                if(m_bndConditions[n]->GetBoundaryConditionType() == 
                       SpatialDomains::eDirichlet)
                {
                    for(e = 0; e < m_bndCondExpansions[n]->GetExpSize(); ++e)
                    {
                        npts = m_bndCondExpansions[n]->
                                GetExp(e)->GetTotPoints();
                        id2  = m_trace->GetPhys_Offset(
                            m_traceMap->GetBndCondIDToGlobalTraceID(cnt+e));
                        Vmath::Vcopy(npts, &Fwd[id2], 1, &Bwd[id2], 1);
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
                        npts = m_bndCondExpansions[n]->
                                GetExp(e)->GetTotPoints();
                        id2  = m_trace->GetPhys_Offset(
                            m_traceMap->GetBndCondIDToGlobalTraceID(cnt+e));
                        
                        // Turning this off since we can have non-zero
                        //Neumann in mixed CG-DG method
                        
                        Vmath::Vcopy(npts, &Fwd[id2], 1, &Bwd[id2], 1);
                    }

                    cnt += e;
                }
                else if (m_bndConditions[n]->GetBoundaryConditionType() ==
                             SpatialDomains::ePeriodic)
                {
                    continue;
                }
                else
                {
                    ASSERTL0(false, "Method only set up for Dirichlet, Neumann "
                             "and Robin conditions.");
                }
            }
        }

        /**
         * @brief Fill the weight with m_bndCondBndWeight.
         */
        void DisContField3D::v_FillBwdWithBwdWeight(
                  Array<OneD,       NekDouble> &weightave,
                  Array<OneD,       NekDouble> &weightjmp)
        {
            int cnt = 0;
            int e = 0;
            int npts = 0;
            int id2 = 0;
            
            for(int n = 0; n < m_bndCondExpansions.size(); ++n)
            {
                if(m_bndConditions[n]->GetBoundaryConditionType() == 
                       SpatialDomains::eDirichlet)
                {
                    for(e = 0; e < m_bndCondExpansions[n]->GetExpSize(); ++e)
                    {
                        npts = m_bndCondExpansions[n]->
                                GetExp(e)->GetTotPoints();
                        id2  = m_trace->GetPhys_Offset(
                            m_traceMap->GetBndCondIDToGlobalTraceID(cnt+e));
                        Vmath::Fill(npts,
                                    m_bndCondBndWeight[n], 
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
                    for(e = 0; e < m_bndCondExpansions[n]->GetExpSize(); ++e)
                    {
                        npts = m_bndCondExpansions[n]->
                                GetExp(e)->GetTotPoints();
                        id2  = m_trace->GetPhys_Offset(
                            m_traceMap->GetBndCondIDToGlobalTraceID(cnt+e));
                        
                        Vmath::Fill(npts,
                                    m_bndCondBndWeight[n], 
                                    &weightave[id2], 1);
                        Vmath::Fill(npts, 0.0, &weightjmp[id2], 1);
                    }

                    cnt += e;
                }
                else if (m_bndConditions[n]->GetBoundaryConditionType() ==
                             SpatialDomains::ePeriodic)
                {
                    continue;
                }
                else
                {
                    ASSERTL0(false, "Method only set up for Dirichlet, Neumann "
                             "and Robin conditions.");
                }
            }
        }

         const vector<bool> &DisContField3D::v_GetLeftAdjacentFaces(void) const
        {
            return m_leftAdjacentFaces;
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

            Vmath::Zero(outarray.size(), outarray, 1);

            Array<OneD, NekDouble> facevals(m_locTraceToTraceMap->GetNFwdLocTracePts());
            m_locTraceToTraceMap->FwdLocTracesFromField(inarray,facevals);
            m_locTraceToTraceMap->InterpLocFacesToTrace(0,facevals,outarray);

            // gather entries along parallel partitions which have
            // only filled in Fwd part on their own partition
            m_traceMap->GetAssemblyCommDG()->PerformExchange(outarray, outarray);
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

            Array<OneD, NekDouble> Fcoeffs(m_trace->GetNcoeffs());
            m_trace->IProductWRTBase(Fn, Fcoeffs);

            m_locTraceToTraceMap->AddTraceCoeffsToFieldCoeffs(Fcoeffs,
                                                              outarray);
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
            Array<OneD, NekDouble> Coeffs(m_trace->GetNcoeffs());

            m_trace->IProductWRTBase(Fwd,Coeffs);
            m_locTraceToTraceMap->AddTraceCoeffsToFieldCoeffs(0,Coeffs,outarray);
            m_trace->IProductWRTBase(Bwd,Coeffs);
            m_locTraceToTraceMap->AddTraceCoeffsToFieldCoeffs(1,Coeffs,outarray);
        }

        /**
         * @brief Set up a list of elemeent IDs and edge IDs that link to the
         * boundary conditions.
         */
        void DisContField3D::v_GetBoundaryToElmtMap(
            Array<OneD, int> &ElmtID,
            Array<OneD, int> &FaceID)
        {
            if (m_BCtoElmMap.size() == 0)
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
                    exp3d = (*m_exp)[i]->as<LocalRegions::Expansion3D>();
                    globalIdMap[exp3d->GetGeom3D()->GetGlobalID()] = i;
                }

                // Determine number of boundary condition expansions.
                for(i = 0; i < m_bndConditions.size(); ++i)
                {
                    nbcs += m_bndCondExpansions[i]->GetExpSize();
                }

                // Initialize arrays
                m_BCtoElmMap = Array<OneD, int>(nbcs);
                m_BCtoFaceMap = Array<OneD, int>(nbcs);

                LocalRegions::Expansion2DSharedPtr exp2d;
                for(cnt = n = 0; n < m_bndCondExpansions.size(); ++n)
                {
                    for(i = 0; i < m_bndCondExpansions[n]->GetExpSize(); ++i, ++cnt)
                    {
                        exp2d = m_bndCondExpansions[n]->GetExp(i)->
                                            as<LocalRegions::Expansion2D>();

                        SpatialDomains::GeometryLinkSharedPtr tmp =
                            m_graph->GetElementsFromFace(exp2d->GetGeom2D());
                        m_BCtoElmMap[cnt] = globalIdMap[
                            tmp->at(0).first->GetGlobalID()];
                        m_BCtoFaceMap[cnt] = tmp->at(0).second;
                    }
                }
            }
            ElmtID = m_BCtoElmMap;
            FaceID = m_BCtoFaceMap;
        }

        void DisContField3D::v_GetBndElmtExpansion(int i,
                            std::shared_ptr<ExpList> &result,
                            const bool DeclareCoeffPhysArrays)
        {
            int n, cnt, nq;
            int offsetOld, offsetNew;
            std::vector<unsigned int> eIDs;

            Array<OneD, int> ElmtID,EdgeID;
            GetBoundaryToElmtMap(ElmtID,EdgeID);

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
                MemoryManager<ExpList3D>::AllocateSharedPtr
                    (*this, eIDs, DeclareCoeffPhysArrays);

            // Copy phys and coeffs to new explist
            if (DeclareCoeffPhysArrays)
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
        void DisContField3D::v_Reset()
        {
            ExpList::v_Reset();

            // Reset boundary condition expansions.
            for (int n = 0; n < m_bndCondExpansions.size(); ++n)
            {
                m_bndCondBndWeight[n]   =   0.0;
                m_bndCondExpansions[n]->Reset();
            }
        }

        /**
         * Solving Helmholtz Equation in 3D
         */
        void DisContField3D::v_HelmSolve(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray,
                const StdRegions::ConstFactorMap &factors,
                const StdRegions::VarCoeffMap &varcoeff,
                const MultiRegions::VarFactorsMap &varfactors,
                const Array<OneD, const NekDouble> &dirForcing,
                const bool PhysSpaceForcing)
        {
            boost::ignore_unused(varfactors, dirForcing);

            int i,j,n,cnt,cnt1,nbndry;
            int nexp = GetExpSize();
            StdRegions::StdExpansionSharedPtr BndExp;

            Array<OneD,NekDouble> f(m_ncoeffs);
            DNekVec F(m_ncoeffs,f,eWrapper);
            Array<OneD,NekDouble> e_f, e_l;

            //----------------------------------
            //  Setup RHS Inner product
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
            //  Solve continuous flux System
            //----------------------------------
            int GloBndDofs   = m_traceMap->GetNumGlobalBndCoeffs();
            int NumDirichlet = m_traceMap->GetNumLocalDirBndCoeffs();
            int e_ncoeffs;

            // Retrieve block matrix of U^e
            GlobalMatrixKey HDGLamToUKey(StdRegions::eHybridDGLamToU,NullAssemblyMapSharedPtr,factors,varcoeff);
            const DNekScalBlkMatSharedPtr &HDGLamToU = GetBlockMatrix(HDGLamToUKey);

            // Retrieve number of local trace space coefficients N_{\lambda},
            // and set up local elemental trace solution \lambda^e.
            int     LocBndCoeffs = m_traceMap->GetNumLocalBndCoeffs();
            Array<OneD, NekDouble> bndrhs(LocBndCoeffs,0.0);
            Array<OneD, NekDouble> loclambda(LocBndCoeffs,0.0);
            DNekVec LocLambda(LocBndCoeffs,loclambda,eWrapper);

            //----------------------------------
            // Evaluate Trace Forcing vector F
            // Kirby et al, 2010, P23, Step 5.
            //----------------------------------
            // Loop over all expansions in the domain
            for(cnt = cnt1 = n = 0; n < nexp; ++n)
            {
                nbndry = (*m_exp)[n]->NumDGBndryCoeffs();

                e_ncoeffs = (*m_exp)[n]->GetNcoeffs();
                e_f       = f + cnt;
                e_l       = bndrhs + cnt1;

                // Local trace space \lambda^e
                DNekVec     Floc    (nbndry, e_l, eWrapper);
                // Local forcing f^e
                DNekVec     ElmtFce (e_ncoeffs, e_f, eWrapper);
                // Compute local (U^e)^{\top} f^e
                Floc = Transpose(*(HDGLamToU->GetBlock(n,n)))*ElmtFce;

                cnt   += e_ncoeffs;
                cnt1  += nbndry;
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
                    for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); ++j)
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
                    for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); ++j)
                    {
                        locid = bndCondMap[cnt + j];
                        bndrhs[locid] += Sign[locid]*bndcoeffs[j]; 
                    }
                }
                else if (m_bndConditions[i]->GetBoundaryConditionType() ==
                             SpatialDomains::ePeriodic)
                {
                    ASSERTL0(false, "HDG implementation does not support "
                             "periodic boundary conditions at present.");
                }
                cnt += (m_bndCondExpansions[i])->GetNcoeffs();
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
                LinSys->Solve(bndrhs,loclambda,m_traceMap);

                // For consistency with previous version put global
                // solution into m_trace->m_coeffs
                m_traceMap->LocalToGlobal(loclambda,m_trace->UpdateCoeffs());
            }

            //----------------------------------
            // Internal element solves
            //----------------------------------
            GlobalMatrixKey invHDGhelmkey(StdRegions::eInvHybridDGHelmholtz,NullAssemblyMapSharedPtr,factors,varcoeff);
            const DNekScalBlkMatSharedPtr& InvHDGHelm = GetBlockMatrix(invHDGhelmkey);
            DNekVec out(m_ncoeffs,outarray,eWrapper);
            Vmath::Zero(m_ncoeffs,outarray,1);

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
               Array<OneD,      NekDouble> &outarray)
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

            for(cnt = i = 0; i < m_bndCondExpansions.size(); ++i)
            {
                MultiRegions::ExpListSharedPtr locExpList;

                if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eRobin)
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
                            ::AllocateSharedPtr(FaceID[cnt+e],
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
        void  DisContField3D::EvaluateHDGPostProcessing(
            Array<OneD, NekDouble> &outarray)
        {
            int    i,cnt,f,ncoeff_face;
            Array<OneD, NekDouble> force, out_tmp,qrhs,qrhs1;
            Array<OneD, Array< OneD, LocalRegions::ExpansionSharedPtr> >
                &elmtToTrace = m_traceMap->GetElmtToTrace();

            int     nq_elmt, nm_elmt;
            int     LocBndCoeffs = m_traceMap->GetNumLocalBndCoeffs();
            Array<OneD, NekDouble> loc_lambda(LocBndCoeffs), face_lambda;
            Array<OneD, NekDouble> tmp_coeffs;
            m_traceMap->GlobalToLocalBnd(m_trace->GetCoeffs(),loc_lambda);

            face_lambda = loc_lambda;

            // Calculate Q using standard DG formulation.
            for(i = cnt = 0; i < GetExpSize(); ++i)
            {
                LocalRegions::Expansion3DSharedPtr exp =
                        (*m_exp)[i]->as<LocalRegions::Expansion3D>();

                nq_elmt = (*m_exp)[i]->GetTotPoints();
                nm_elmt = (*m_exp)[i]->GetNcoeffs();
                qrhs    = Array<OneD, NekDouble>(nq_elmt);
                qrhs1   = Array<OneD, NekDouble>(nq_elmt);
                force   = Array<OneD, NekDouble>(2*nm_elmt);
                out_tmp = force + nm_elmt;
                LocalRegions::ExpansionSharedPtr ppExp;

                int num_points0 = (*m_exp)[i]->GetBasis(0)->GetNumPoints();
                int num_points1 = (*m_exp)[i]->GetBasis(1)->GetNumPoints();
                int num_points2 = (*m_exp)[i]->GetBasis(2)->GetNumPoints();
                int num_modes0 = (*m_exp)[i]->GetBasis(0)->GetNumModes();
                int num_modes1 = (*m_exp)[i]->GetBasis(1)->GetNumModes();
                int num_modes2 = (*m_exp)[i]->GetBasis(2)->GetNumModes();

                // Probably a better way of setting up lambda than this.  Note
                // cannot use PutCoeffsInToElmts since lambda space is mapped
                // during the solve.
                int nFaces = (*m_exp)[i]->GetNfaces();
                Array<OneD, Array<OneD, NekDouble> > faceCoeffs(nFaces);
                for(f = 0; f < nFaces; ++f)
                {
                    ncoeff_face = elmtToTrace[i][f]->GetNcoeffs();
                    faceCoeffs[f] = Array<OneD, NekDouble>(ncoeff_face);
                    Vmath::Vcopy(ncoeff_face, face_lambda, 1, faceCoeffs[f], 1);
                    exp->SetFaceToGeomOrientation(f, faceCoeffs[f]);
                    face_lambda = face_lambda + ncoeff_face;
                }

                //creating orthogonal expansion (checking if we have quads or triangles)
                LibUtilities::ShapeType shape = (*m_exp)[i]->DetShapeType();
                switch(shape)
                {
                    case LibUtilities::eHexahedron:
                    {
                        const LibUtilities::PointsKey PkeyH1(num_points0,LibUtilities::eGaussLobattoLegendre);
                        const LibUtilities::PointsKey PkeyH2(num_points1,LibUtilities::eGaussLobattoLegendre);
                        const LibUtilities::PointsKey PkeyH3(num_points2,LibUtilities::eGaussLobattoLegendre);
                        LibUtilities::BasisKey  BkeyH1(LibUtilities::eOrtho_A, num_modes0, PkeyH1);
                        LibUtilities::BasisKey  BkeyH2(LibUtilities::eOrtho_A, num_modes1, PkeyH2);
                        LibUtilities::BasisKey  BkeyH3(LibUtilities::eOrtho_A, num_modes2, PkeyH3);
                        SpatialDomains::HexGeomSharedPtr hGeom = std::dynamic_pointer_cast<SpatialDomains::HexGeom>((*m_exp)[i]->GetGeom());
                        ppExp = MemoryManager<LocalRegions::HexExp>::AllocateSharedPtr(BkeyH1, BkeyH2, BkeyH3, hGeom);
                    }
                    break;
                    case LibUtilities::eTetrahedron:
                    {
                        const LibUtilities::PointsKey PkeyT1(num_points0,LibUtilities::eGaussLobattoLegendre);
                        const LibUtilities::PointsKey PkeyT2(num_points1,LibUtilities::eGaussRadauMAlpha1Beta0);
                        const LibUtilities::PointsKey PkeyT3(num_points2,LibUtilities::eGaussRadauMAlpha2Beta0);
                        LibUtilities::BasisKey  BkeyT1(LibUtilities::eOrtho_A, num_modes0, PkeyT1);
                        LibUtilities::BasisKey  BkeyT2(LibUtilities::eOrtho_B, num_modes1, PkeyT2);
                        LibUtilities::BasisKey  BkeyT3(LibUtilities::eOrtho_C, num_modes2, PkeyT3);
                        SpatialDomains::TetGeomSharedPtr tGeom = std::dynamic_pointer_cast<SpatialDomains::TetGeom>((*m_exp)[i]->GetGeom());
                        ppExp = MemoryManager<LocalRegions::TetExp>::AllocateSharedPtr(BkeyT1, BkeyT2, BkeyT3, tGeom);
                    }
                    break;
                    case LibUtilities::ePrism:
                    {
                        const LibUtilities::PointsKey PkeyP1(num_points0,LibUtilities::eGaussLobattoLegendre);
                        const LibUtilities::PointsKey PkeyP2(num_points1,LibUtilities::eGaussLobattoLegendre);
                        const LibUtilities::PointsKey PkeyP3(num_points2,LibUtilities::eGaussRadauMAlpha1Beta0);
                        LibUtilities::BasisKey  BkeyP1(LibUtilities::eOrtho_A, num_modes0, PkeyP1);
                        LibUtilities::BasisKey  BkeyP2(LibUtilities::eOrtho_A, num_modes1, PkeyP2);
                        LibUtilities::BasisKey  BkeyP3(LibUtilities::eOrtho_B, num_modes2, PkeyP3);
                        SpatialDomains::PrismGeomSharedPtr pGeom = std::dynamic_pointer_cast<SpatialDomains::PrismGeom>((*m_exp)[i]->GetGeom());
                        ppExp = MemoryManager<LocalRegions::PrismExp>::AllocateSharedPtr(BkeyP1, BkeyP2, BkeyP3, pGeom);
                    }
                    break;
                    default:
                        ASSERTL0(false, "Wrong shape type, HDG postprocessing is not implemented");
                };


                //DGDeriv
                // (d/dx w, q_0)
                (*m_exp)[i]->DGDeriv(
                    0,tmp_coeffs = m_coeffs + m_coeff_offset[i],
                    elmtToTrace[i], faceCoeffs, out_tmp);
                (*m_exp)[i]->BwdTrans(out_tmp,qrhs);
                ppExp->IProductWRTDerivBase(0,qrhs,force);


                // + (d/dy w, q_1)
                (*m_exp)[i]->DGDeriv(
                    1,tmp_coeffs = m_coeffs + m_coeff_offset[i],
                    elmtToTrace[i], faceCoeffs, out_tmp);
                (*m_exp)[i]->BwdTrans(out_tmp,qrhs);
                ppExp->IProductWRTDerivBase(1,qrhs,out_tmp);

                Vmath::Vadd(nm_elmt,force,1,out_tmp,1,force,1);

                // + (d/dz w, q_2)
                (*m_exp)[i]->DGDeriv(
                    2,tmp_coeffs = m_coeffs + m_coeff_offset[i],
                    elmtToTrace[i], faceCoeffs, out_tmp);
                (*m_exp)[i]->BwdTrans(out_tmp,qrhs);
                ppExp->IProductWRTDerivBase(2,qrhs,out_tmp);

                Vmath::Vadd(nm_elmt,force,1,out_tmp,1,force,1);
                // determine force[0] = (1,u)
                (*m_exp)[i]->BwdTrans(
                    tmp_coeffs = m_coeffs + m_coeff_offset[i],qrhs);
                force[0] = (*m_exp)[i]->Integral(qrhs);

                // multiply by inverse Laplacian matrix
                // get matrix inverse
                LocalRegions::MatrixKey  lapkey(StdRegions::eInvLaplacianWithUnityMean, ppExp->DetShapeType(), *ppExp);
                DNekScalMatSharedPtr lapsys = ppExp->GetLocMatrix(lapkey);

                NekVector<NekDouble> in (nm_elmt, force, eWrapper);
                NekVector<NekDouble> out(nm_elmt);

                out = (*lapsys)*in;

                // Transforming back to modified basis
                Array<OneD, NekDouble> work(nq_elmt);
                ppExp->BwdTrans(out.GetPtr(), work);
                (*m_exp)[i]->FwdTrans(work,
                                tmp_coeffs = outarray + m_coeff_offset[i]);
            }
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
        void DisContField3D::v_EvaluateBoundaryConditions(
            const NekDouble   time,
            const std::string varName,
            const NekDouble   x2_in,
            const NekDouble   x3_in)
        {
            boost::ignore_unused(x2_in, x3_in);

            int i;
            int npoints;
            int nbnd = m_bndCondExpansions.size();
            MultiRegions::ExpListSharedPtr locExpList;

            for (i = 0; i < nbnd; ++i)
            {
                if (time == 0.0 || m_bndConditions[i]->IsTimeDependent())
                {
                    m_bndCondBndWeight[i]   =   1.0;
                    locExpList = m_bndCondExpansions[i];
                    npoints    = locExpList->GetNpoints();

                    Array<OneD, NekDouble> x0(npoints, 0.0);
                    Array<OneD, NekDouble> x1(npoints, 0.0);
                    Array<OneD, NekDouble> x2(npoints, 0.0);
                    Array<OneD, NekDouble> valuesFile(npoints, 1.0), valuesExp(npoints, 1.0);

                    locExpList->GetCoords(x0, x1, x2);

                    if (m_bndConditions[i]->GetBoundaryConditionType()
                        == SpatialDomains::eDirichlet)
                    {
                        SpatialDomains::DirichletBCShPtr bcPtr = std::static_pointer_cast<
                            SpatialDomains::DirichletBoundaryCondition>(
                                m_bndConditions[i]);
                        string filebcs = bcPtr->m_filename;
                        string exprbcs = bcPtr->m_expr;

                        if (filebcs != "")
                        {
                            ExtractFileBCs(filebcs, bcPtr->GetComm(), varName, locExpList);
                            valuesFile = locExpList->GetPhys();
                        }

                        if (exprbcs != "")
                        {
                            LibUtilities::Equation  condition = std::static_pointer_cast<SpatialDomains::
                                    DirichletBoundaryCondition >(
                                    m_bndConditions[i])->m_dirichletCondition;

                            condition.Evaluate(x0, x1, x2, time, valuesExp);
                        }

                        Vmath::Vmul(npoints, valuesExp, 1, valuesFile, 1, locExpList->UpdatePhys(), 1);

                        locExpList->FwdTrans_BndConstrained(
                            locExpList->GetPhys(),
                            locExpList->UpdateCoeffs());
                    }
                    else if (m_bndConditions[i]->GetBoundaryConditionType()
                             == SpatialDomains::eNeumann)
                    {
                        SpatialDomains::NeumannBCShPtr bcPtr = std::static_pointer_cast<
                            SpatialDomains::NeumannBoundaryCondition>(
                                m_bndConditions[i]);
                        string filebcs = bcPtr->m_filename;

                        if (filebcs != "")
                        {
                            ExtractFileBCs(filebcs, bcPtr->GetComm(), varName, locExpList);
                        }
                        else
                        {

                            LibUtilities::Equation condition = std::
                                static_pointer_cast<SpatialDomains::
                                                    NeumannBoundaryCondition>(
                                    m_bndConditions[i])->m_neumannCondition;

                            condition.Evaluate(x0, x1, x2, time,
                                               locExpList->UpdatePhys());

                            locExpList->IProductWRTBase(locExpList->GetPhys(),
                                                        locExpList->UpdateCoeffs());
                        }
                    }
                    else if (m_bndConditions[i]->GetBoundaryConditionType()
                             == SpatialDomains::eRobin)
                    {
                        SpatialDomains::RobinBCShPtr bcPtr = std::static_pointer_cast<
                            SpatialDomains::RobinBoundaryCondition>(
                                m_bndConditions[i]);
                        string filebcs = bcPtr->m_filename;

                        if (filebcs != "")
                        {
                            ExtractFileBCs(filebcs, bcPtr->GetComm(), varName, locExpList);
                        }
                        else
                        {
                            LibUtilities::Equation condition = std::
                                static_pointer_cast<SpatialDomains::
                                                    RobinBoundaryCondition>(
                                    m_bndConditions[i])->m_robinFunction;

                            condition.Evaluate(x0, x1, x2, time,
                                               locExpList->UpdatePhys());

                        }

                        locExpList->IProductWRTBase(locExpList->GetPhys(),
                                                    locExpList->UpdateCoeffs());

                    }
                    else if (m_bndConditions[i]->GetBoundaryConditionType()
                             == SpatialDomains::ePeriodic)
                    {
                        continue;
                    }
                    else
                    {
                        ASSERTL0(false, "This type of BC not implemented yet");
                    }
                }
            }
        }
    } // end of namespace
} // end of namespace
