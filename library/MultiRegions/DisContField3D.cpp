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

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class DisContField3D
         * The class #DisContField3D extends an ExpList3D object through the
         * addition of boundary conditions in a discontinuous galerkin
         * formulation.
         */

        /**
         *
         */
        DisContField3D::DisContField3D() :
            ExpList3D             (),
            m_bndCondExpansions   (),
            m_bndConditions       ()
        {
        }


        /**
         * 
         */
        DisContField3D::DisContField3D(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const SpatialDomains::MeshGraphSharedPtr   &graph3D,
            const std::string                          &variable,
            const bool                                  SetUpJustDG)
        : ExpList3D          (pSession,graph3D),
          m_bndCondExpansions(),
          m_bndConditions    ()
        {
            SpatialDomains::BoundaryConditions bcs(m_session, graph3D);
            
            GenerateBoundaryConditionExpansion(graph3D,bcs,variable);
            EvaluateBoundaryConditions();
            ApplyGeomInfo();
            
            if(SetUpJustDG)
            {
                ExpList2DSharedPtr trace;
                
                // Set up matrix map
                m_globalBndMat = MemoryManager<GlobalLinSysMap>::AllocateSharedPtr();
                map<int,int>          periodicEdges;
                map<int,int>          periodicVertices;
                map<int,PeriodicFace> periodicFaces;
                GetPeriodicFaces(graph3D, bcs, variable,
                                 periodicVertices, periodicEdges, periodicFaces);
                
                // Set up Trace space
                bool UseGenSegExp = true;
                trace = MemoryManager<ExpList2D>
                    ::AllocateSharedPtr(m_bndCondExpansions, m_bndConditions,
                                        *m_exp,graph3D, periodicFaces, UseGenSegExp);
                m_trace         = trace;
                m_periodicFaces = periodicFaces;
                
                // Scatter trace segments to 3D elements. For each element,
                // we find the trace segment associated to each face. The
                // element then retains a pointer to the trace space segments
                SpatialDomains::Geometry2DSharedPtr ElmtFaceGeom;
                SpatialDomains::Geometry2DSharedPtr TraceFaceGeom;
                for (int i = 0; i < m_exp->size(); ++i)
                {
                    for (int j = 0; j < (*m_exp)[i]->GetNfaces(); ++j)
                    {
                        ElmtFaceGeom  = ((*m_exp)[i]->GetGeom3D())->GetFace(j);
                        for (int k = 0; k < m_trace->GetExpSize(); ++k)
                        {
                            TraceFaceGeom = m_trace->GetExp(k)->GetGeom2D();
                            if (TraceFaceGeom == ElmtFaceGeom)
                            {
                                LocalRegions::Expansion3DSharedPtr exp3d
                                    = boost::dynamic_pointer_cast<LocalRegions::Expansion3D>((*m_exp)[i]);
                                LocalRegions::Expansion2DSharedPtr exp2d
                                    = boost::dynamic_pointer_cast<LocalRegions::Expansion2D>(m_trace->GetExp(k));
                                
                                exp3d->SetFaceExp(j, exp2d);
                                exp2d->SetAdjacentElementExp(j, exp3d);
                                break;
                            }
                        }
                    }
                }
                SetUpPhysNormals();
                
                m_traceMap = MemoryManager<AssemblyMapDG>::AllocateSharedPtr(
                    m_session,graph3D,trace,*this,m_bndCondExpansions,
                    m_bndConditions, periodicFaces);

                // Set up information for periodic boundary conditions.
                for (int n = 0; n < m_exp->size(); ++n)
                {
                    for (int e = 0; e < (*m_exp)[n]->GetNfaces(); ++e)
                    {
                        map<int,PeriodicFace>::iterator it = m_periodicFaces.find(
                            (*m_exp)[n]->GetGeom3D()->GetFid(e));
                        
                        if (it != m_periodicFaces.end())
                        {
                            m_perFaceToExpMap[it->first] = make_pair(n, e);
                        }
                    }
                }
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
                            = boost::dynamic_pointer_cast<LocalRegions::Expansion3D>((*m_exp)[ElmtID[cnt+f]]);
                        LocalRegions::Expansion2DSharedPtr exp2d
                            = boost::dynamic_pointer_cast<LocalRegions::Expansion2D>(locExpList->GetExp(f));
                        
                        exp3d->SetFaceExp(FaceID[cnt+f],exp2d);
                        exp2d->SetAdjacentElementExp(FaceID[cnt+f],exp3d);
                    }
                    
                    cnt += m_bndCondExpansions[i]->GetExpSize();
                }
                //normals computation currently not implemented for Prisms, breaks regression tests
                //SetUpPhysNormals();
            }
        }
        
        /*
         * Copy type constructor which declares new boundary conditions
         * and re-uses mapping info and trace space if possible
         */
        DisContField3D::DisContField3D( const DisContField3D &In,
                                        const SpatialDomains::MeshGraphSharedPtr &graph3D,
                                        const std::string &variable,
                                        const bool SetUpJustDG) :
            ExpList3D(In)
       {
           SpatialDomains::BoundaryConditions bcs(m_session, graph3D);
           
           GenerateBoundaryConditionExpansion(graph3D,bcs,variable);
           EvaluateBoundaryConditions();
           ApplyGeomInfo();
           
           if(!SameTypeOfBoundaryConditions(In))
           {
               if(SetUpJustDG)
               {
                   ExpList2DSharedPtr trace;
                   
                   // Set up matrix map
                   m_globalBndMat = MemoryManager<GlobalLinSysMap>::AllocateSharedPtr();
                   map<int,int>          periodicEdges;
                   map<int,int>          periodicVertices;
                   map<int,PeriodicFace> periodicFaces;
                   GetPeriodicFaces(graph3D, bcs, variable,
                                    periodicVertices, periodicEdges, periodicFaces);
                   
                   // Set up Trace space
                   bool UseGenSegExp = true;
                   trace = MemoryManager<ExpList2D>
                       ::AllocateSharedPtr(m_bndCondExpansions, m_bndConditions,
                                           *m_exp,graph3D, periodicFaces, UseGenSegExp);
                   m_trace = boost::dynamic_pointer_cast<ExpList>(trace);
                   
                   // Scatter trace segments to 3D elements. For each element,
                   // we find the trace segment associated to each face. The
                   // element then retains a pointer to the trace space segments
                   SpatialDomains::Geometry2DSharedPtr ElmtFaceGeom;
                   SpatialDomains::Geometry2DSharedPtr TraceFaceGeom;
                   for (int i = 0; i < m_exp->size(); ++i)
                   {
                       for (int j = 0; j < (*m_exp)[i]->GetNfaces(); ++j)
                       {
                           ElmtFaceGeom  = ((*m_exp)[i]->GetGeom3D())->GetFace(j);
                           for (int k = 0; k < m_trace->GetExpSize(); ++k)
                           {
                               TraceFaceGeom = m_trace->GetExp(k)->GetGeom2D();
                               if (TraceFaceGeom == ElmtFaceGeom)
                               {
                                   LocalRegions::Expansion3DSharedPtr exp3d
                                       = boost::dynamic_pointer_cast<LocalRegions::Expansion3D>((*m_exp)[i]);
                                   LocalRegions::Expansion2DSharedPtr exp2d
                                       = boost::dynamic_pointer_cast<LocalRegions::Expansion2D>(m_trace->GetExp(k));
                                   
                                   exp3d->SetFaceExp(j, exp2d);
                                   exp2d->SetAdjacentElementExp(j, exp3d);
                                   break;
                               }
                           }
                       }
                   }
                   SetUpPhysNormals();
                   
                   m_traceMap = MemoryManager<AssemblyMapDG>::AllocateSharedPtr(
                       m_session,graph3D,trace,*this, m_bndCondExpansions,
                       m_bndConditions, periodicFaces);
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
                               = boost::dynamic_pointer_cast<LocalRegions::Expansion3D>((*m_exp)[ElmtID[cnt+f]]);
                           LocalRegions::Expansion2DSharedPtr exp2d
                               = boost::dynamic_pointer_cast<LocalRegions::Expansion2D>(locExpList->GetExp(f));
                           
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
                               = boost::dynamic_pointer_cast<LocalRegions::Expansion3D>((*m_exp)[ElmtID[cnt+f]]);
                           LocalRegions::Expansion2DSharedPtr exp2d
                               = boost::dynamic_pointer_cast<LocalRegions::Expansion2D>(locExpList->GetExp(f));
                           
                           exp3d->SetFaceExp(FaceID[cnt+f],exp2d);
                           exp2d->SetAdjacentElementExp(FaceID[cnt+f],exp3d);
                       }
                       
                       cnt += m_bndCondExpansions[i]->GetExpSize();
                   }
                   SetUpPhysNormals();
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
         *
         */
        DisContField3D::~DisContField3D()
        {
        }

        GlobalLinSysSharedPtr DisContField3D::GetGlobalBndLinSys(const GlobalLinSysKey &mkey)
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
        
        /**
         * According to their boundary region, the separate segmental boundary
         * expansions are bundled together in an object of the class
         * MultiRegions#ExpList1D.
         * The list of expansions of the Dirichlet boundary regions are listed
         * first in the array #m_bndCondExpansions.
         *
         * \param graph2D A mesh, containing information about the domain and
         * the spectral/hp element expansion.
         * \param bcs An entity containing information about the boundaries and
         * boundary conditions.
         * \param variable An optional parameter to indicate for which variable
         * the boundary conditions should be discretised.
         */
        void DisContField3D::GenerateBoundaryConditionExpansion(
            const SpatialDomains::MeshGraphSharedPtr &graph3D,
            const SpatialDomains::BoundaryConditions &bcs,
            const std::string &variable)
        {
            int cnt1  = 0;
            const SpatialDomains::BoundaryRegionCollection    &bregions = 
                bcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions = 
                bcs.GetBoundaryConditions();

            int nbnd = bregions.size();

            // count the number of non-periodic boundary regions
            for(int i = 0; i < nbnd; ++i)
            {
                SpatialDomains::BoundaryConditionShPtr boundaryCondition = 
                    GetBoundaryCondition(bconditions, i, variable);
                if (boundaryCondition->GetBoundaryConditionType() != 
                        SpatialDomains::ePeriodic)
                {
                    cnt1++;
                }
            }

            m_bndCondExpansions  = Array<OneD,MultiRegions::ExpListSharedPtr>(cnt1);
            m_bndConditions      = Array<OneD,SpatialDomains::BoundaryConditionShPtr>(cnt1);

            SetBoundaryConditionExpansion(graph3D,bcs,variable,m_bndCondExpansions,m_bndConditions);
        }


        /**
         * For each boundary region, checks that the types and number of
         * boundary expansions in that region match.
         * @param   In          ContField2D to compare with.
         * @returns True if boundary conditions match.
         */
        bool DisContField3D::SameTypeOfBoundaryConditions(const DisContField3D &In)
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
         * @param   graph3D     A mesh containing information about the domain
         *                      and the spectral/hp element expansions.
         * @param   bcs         Information about the boundary conditions.
         * @param   variable    Specifies the field.
         * @param   bndCondExpansions   Array of ExpList2D objects each
         *                      containing a 2D spectral/hp element expansion
         *                      on a single boundary region.
         * @param   bndConditions   Array of BoundaryCondition objects which
         *                      contain information about the boundary
         *                      conditions on the different boundary regions.
         */
        void DisContField3D::SetBoundaryConditionExpansion(
                        const SpatialDomains::MeshGraphSharedPtr &graph3D,
                        const SpatialDomains::BoundaryConditions &bcs,
                        const std::string variable,
                        Array<OneD, ExpListSharedPtr> &bndCondExpansions,
                        Array<OneD, SpatialDomains::BoundaryConditionShPtr>
                        &bndConditions)
        {
            int i;
            int cnt  = 0;

            const SpatialDomains::BoundaryRegionCollection &bregions
                                        = bcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions
                                        = bcs.GetBoundaryConditions();

            MultiRegions::ExpList2DSharedPtr       locExpList;
            SpatialDomains::BoundaryConditionShPtr locBCond;

            int nbnd = bregions.size();

            cnt=0;
            // list Dirichlet boundaries first
            for(i = 0; i < nbnd; ++i)
            {
                locBCond = GetBoundaryCondition(bconditions, i, variable);

                if(locBCond->GetBoundaryConditionType()
                                        == SpatialDomains::eDirichlet)
                {
                    locExpList = MemoryManager<MultiRegions::ExpList2D>
                                        ::AllocateSharedPtr(m_session,*(bregions[i]),
                                                            graph3D);
                    bndCondExpansions[cnt]  = locExpList;
                    bndConditions[cnt++]    = locBCond;
                } // end if Dirichlet
            }
            // then, list the other (non-periodic) boundaries
            for(i = 0; i < nbnd; ++i)
            {
                locBCond = GetBoundaryCondition(bconditions, i, variable);

                switch(locBCond->GetBoundaryConditionType())
                {
                case SpatialDomains::eNeumann:
                case SpatialDomains::eRobin:
                    {
                        locExpList = MemoryManager<MultiRegions::ExpList2D>
                            ::AllocateSharedPtr(m_session,*(bregions[i]),
                                                graph3D);
                        bndCondExpansions[cnt]  = locExpList;
                        bndConditions[cnt++]    = locBCond;
                        SetUpPhysNormals();
                    }
                    break;
                case SpatialDomains::eDirichlet: // do nothing for these types
                case SpatialDomains::ePeriodic:
                    break;
                default:
                    ASSERTL0(false,"This type of BC not implemented yet");
                    break;
                }
            }
        }
        
        /**
         * @param   graph3D     A mesh containing information about the domain
         *                      and the spectral/hp element expansion.
         * @param   bcs         Information about the boundary conditions.
         * @param   variable    Specifies the field.
         * @param   periodicVertices    Map into which the list of periodic
         *                      vertices is placed.
         * @param   periodicEdges   Map into which the list of periodic edges
         *                      is placed.
         * @param   periodicFaces   Map into which the list of periodic faces
         *                      is placed.
         */
        void DisContField3D::GetPeriodicFaces(
            const SpatialDomains::MeshGraphSharedPtr &graph3D,
            const SpatialDomains::BoundaryConditions &bcs,
            const std::string                        &variable,
            map<int,int>                             &periodicVertices,
            map<int,int>                             &periodicEdges,
            map<int,PeriodicFace>                    &periodicFaces)
        {
            int i,k,l;

            const SpatialDomains::BoundaryRegionCollection &bregions
                                        = bcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions
                                        = bcs.GetBoundaryConditions();

            int                                        region1ID;
            int                                        region2ID;
            StdRegions::Orientation                    orient1;
            StdRegions::Orientation                    orient2;
            SpatialDomains::Composite                  comp1;
            SpatialDomains::Composite                  comp2;
            SpatialDomains::Geometry2DSharedPtr        faceGeom1;
            SpatialDomains::Geometry2DSharedPtr        faceGeom2;
            SpatialDomains::ElementFaceVectorSharedPtr element1;
            SpatialDomains::ElementFaceVectorSharedPtr element2;
            SpatialDomains::BoundaryConditionShPtr     locBCond;

            // This std::map is a check so that the periodic pairs are not
            // treated twice
            map<int, int> doneBndRegions;
            
            int nbnd = bregions.size();

            for(i = 0; i < nbnd; ++i)
            {
                locBCond = GetBoundaryCondition(bconditions, i, variable);
                
                if(locBCond->GetBoundaryConditionType()
                                        != SpatialDomains::ePeriodic)
                {
                    continue;
                }
                
                region1ID = i;
                region2ID = (boost::static_pointer_cast<SpatialDomains
                             ::PeriodicBoundaryCondition>(locBCond))
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
                            if(!(faceGeom1 = boost::dynamic_pointer_cast<
                                 SpatialDomains::Geometry2D>((*comp1)[k]))||
                               !(faceGeom2 = boost::dynamic_pointer_cast<
                                 SpatialDomains::Geometry2D>((*comp2)[k])))
                            {
                                ASSERTL0(false,"dynamic cast to a "
                                               "Geometry2D failed");
                            }
                            
                            element1 = boost::dynamic_pointer_cast<
                                SpatialDomains::MeshGraph3D>(graph3D)
                                ->GetElementsFromFace(faceGeom1);
                            element2 = boost::dynamic_pointer_cast<
                                SpatialDomains::MeshGraph3D>(graph3D)
                                ->GetElementsFromFace(faceGeom2);
                            
                            ASSERTL0(element1->size() == 1,
                                     "The periodic boundaries belong to "
                                     "more than one element of the mesh");
                            ASSERTL0(element2->size() == 1,
                                     "The periodic boundaries belong to "
                                     "more than one element of the mesh");
                            
                            // Obtain face orientation.
                            SpatialDomains::QuadGeomSharedPtr q1, q2;
                            SpatialDomains::TriGeomSharedPtr  t1, t2;
                            StdRegions::Orientation           forient;
                            
                            if ((q1 = boost::dynamic_pointer_cast<
                                 SpatialDomains::QuadGeom>(faceGeom1)) &&
                                (q2 = boost::dynamic_pointer_cast<
                                 SpatialDomains::QuadGeom>(faceGeom2)))
                            {
                                forient = SpatialDomains::QuadGeom::
                                    GetFaceOrientation(*q1, *q2);
                            }
                            else if ((t1 = boost::dynamic_pointer_cast<
                                      SpatialDomains::TriGeom>(faceGeom1)) &&
                                     (t2 = boost::dynamic_pointer_cast<
                                      SpatialDomains::TriGeom>(faceGeom2)))
                            {
                                forient = SpatialDomains::TriGeom::
                                    GetFaceOrientation(*t1, *t2);
                            }
                            else
                            {
                                ASSERTL0(false, "Failed to cast face.");
                            }
                            
                            // Set periodic faces, along with their
                            // orientation to one another.
                            periodicFaces[faceGeom1->GetFid()] = 
                                pair<int, StdRegions::Orientation>(
                                    faceGeom2->GetFid(), forient);
                            periodicFaces[faceGeom2->GetFid()] =
                                pair<int, StdRegions::Orientation>(
                                    faceGeom1->GetFid(), forient);
                            
                            int nVert = faceGeom1->GetNumVerts();

                            // From face orientation, determine periodic
                            // edges and vertices.
                            if (nVert == 3)
                            {
                                ASSERTL0(forient == 
                                         StdRegions::eDir1FwdDir1_Dir2FwdDir2 ||
                                         forient == 
                                         StdRegions::eDir1BwdDir1_Dir2FwdDir2,
                                         "Unrecognised face orientation for "
                                         "triangular face "+
                                         boost::lexical_cast<string>(
                                             faceGeom1->GetGlobalID())+": "+
                                         boost::lexical_cast<string>(
                                             StdRegions::OrientationMap[forient]));
                                
                                int f1 = (*element1)[0]->m_FaceIndx;
                                int f2 = (*element2)[0]->m_FaceIndx;
                                
                                // Vertex/edge maps for fwd/bwd orientation in
                                // a-direction.
                                int vmap[2][3] = {{0,1,2}, {1,0,2}};
                                int emap[2][3] = {{0,1,2}, {0,2,1}};
                                
                                for (l = 0; l < nVert; ++l)
                                {
                                    int fo   = ((int)forient - 5)/2;
                                    int vid1 = faceGeom1->GetVid(l);
                                    int vid2 = faceGeom2->GetVid(vmap[fo][l]);
                                    periodicVertices[vid1] = vid2;
                                    periodicVertices[vid2] = vid1;
                                    
                                    int eid1 = faceGeom1->GetEid(l);
                                    int eid2 = faceGeom2->GetEid(emap[fo][l]);
                                    periodicEdges[eid1] = eid2;
                                    periodicEdges[eid2] = eid1;
                                }
                            }
                            else if (nVert == 4)
                            {
                                // Vertex mapping for all possible face
                                // orientations.
                                int vmap[8][4] = {
                                    {0,1,2,3},{3,2,1,0},{1,0,3,2},{2,3,0,1},
                                    {0,3,2,1},{3,0,1,2},{1,2,3,0},{2,1,0,3}
                                };
                                
                                // Edge mapping for all possible face
                                // orientations.
                                int emap[8][4] = {
                                    {0,1,2,3},{2,1,0,3},{0,3,2,1},{2,3,0,1},
                                    {3,2,1,0},{3,0,1,2},{1,2,3,0},{1,0,3,2}
                                };
                                
                                for (l = 0; l < nVert; ++l)
                                {
                                    int fo   = (int)forient - 5;
                                    int vid1 = faceGeom1->GetVid(l);
                                    int vid2 = faceGeom2->GetVid(vmap[fo][l]);
                                    periodicVertices[vid1] = vid2;
                                    periodicVertices[vid2] = vid1;
                                    
                                    int eid1 = faceGeom1->GetEid(l);
                                    int eid2 = faceGeom2->GetEid(emap[fo][l]);
                                    periodicEdges[eid1] = eid2;
                                    periodicEdges[eid2] = eid1;
                                }
                            }
                            else
                            {
                                ASSERTL0(false, "Unknown number of edges!");
                            }
                        }
                    }
                }
                else
                {
                    ASSERTL0(doneBndRegions[region1ID] == region2ID,
                             "Boundary regions are not mutually periodic");
                }
                doneBndRegions[region2ID] = region1ID;
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

        void DisContField3D::v_GetBoundaryToElmtMap(Array<OneD,int> &ElmtID,
                                                    Array<OneD,int> &FaceID)
        {
            map<int,int> globalIdMap;
            int i,n,id;
            int bid,cnt,Fid;
            int nbcs = 0;
            
            SpatialDomains::MeshGraph3DSharedPtr graph3D = 
                boost::dynamic_pointer_cast<SpatialDomains::MeshGraph3D>(m_graph);
            
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
                    
                    ElmtID[cnt] = globalIdMap[(*tmp)[0]->m_Element->GetGlobalID()];
                    FaceID[cnt] = (*tmp)[0]->m_FaceIndx;
                }
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
        void DisContField3D::v_EvaluateBoundaryConditions(const NekDouble time,
                                                          const NekDouble x2_in,
                                                          const NekDouble x3_in)
        {
            int i,j;
            int npoints;
            int nbnd = m_bndCondExpansions.num_elements();
            MultiRegions::ExpListSharedPtr locExpList;

            for(i = 0; i < nbnd; ++i)
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
                                                           SpatialDomains::DirichletBoundaryCondition
                                                        >(m_bndConditions[i])->m_dirichletCondition;

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

                    condition.Evaluate(x0,x1,x2,time,locExpList->UpdatePhys());

                    locExpList->IProductWRTBase(locExpList->GetPhys(),
                                                locExpList->UpdateCoeffs());

                    /// \todo RobinPrimitiveCoeff forgotten? - PB
                }
                else
                {
                    ASSERTL0(false,"This type of BC not implemented yet");
                }
            }
        }

        const Array<OneD,const MultiRegions::ExpListSharedPtr> & DisContField3D::v_GetBndCondExpansions()
        {
            return m_bndCondExpansions;
        }

        const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>& DisContField3D::v_GetBndConditions()
        {
            return m_bndConditions;
        }

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
            int cnt,n,e,npts,offset, phys_offset;
            Array<OneD,NekDouble> e_tmp;
            
            Array<OneD, Array<OneD, StdRegions::StdExpansionSharedPtr> >
                &elmtToTrace = m_traceMap->GetElmtToTrace();

            set<int>::iterator     it;
            map<int,PeriodicFace>::iterator it2;
            boost::unordered_map<int,pair<int,int> >::iterator it3;
            
            // Zero vectors.
            Vmath::Zero(Fwd.num_elements(), Fwd, 1);
            Vmath::Zero(Bwd.num_elements(), Bwd, 1);
            
            for(n = 0; n < nexp; ++n)
            {
                phys_offset = GetPhys_Offset(n);
                for(e = 0; e < (*m_exp)[n]->GetNfaces(); ++e)
                {
                    LocalRegions::Expansion2DSharedPtr traceEl = 
                        boost::dynamic_pointer_cast<
                            LocalRegions::Expansion2D>(elmtToTrace[n][e]);
                    
                    offset = m_trace->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());
                    bool fwd = true;
                    
                    if (traceEl->GetLeftAdjacentElementFace () == -1 ||
                        traceEl->GetRightAdjacentElementFace() == -1)
                    {
                        // Boundary face (1 connected element) - always put in
                        // forwards space.
                    }
                    else if (traceEl->GetLeftAdjacentElementFace () != -1 &&
                             traceEl->GetRightAdjacentElementFace() != -1)
                    {
                        // Non-boundary face (2 connected elements).
                        fwd = dynamic_cast<Nektar::StdRegions::StdExpansion*>
                                    (traceEl->GetLeftAdjacentElementExp().get()) ==
                              (*m_exp)[n].get();
                    }
                    else
                    {
                        ASSERTL2(false, "Unconnected trace element!");
                    }
                    
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
                    
                    // Check to see if this face is periodic.
                    it2 = m_periodicFaces.find(
                        (*m_exp)[n]->GetGeom3D()->GetFid(e));
                    
                    if (it2 != m_periodicFaces.end())
                    {
                        it3 = m_perFaceToExpMap.find(abs(it2->second.first));

                        ASSERTL2(fwd, "Periodic face in non-forward space?");
                        ASSERTL2(it3 != m_perFaceToExpMap.end(),
                                 "Periodic face not found!");
                        
                        int offset2 = m_trace->GetPhys_Offset(
                            elmtToTrace[it3->second.first][it3->second.second]->
                                GetElmtId());
                        
                        // Extract from 3D element to 2D space. We use the
                        // GetFacePhysVals function since the data will need
                        // reordering depending on relative face orientations.
                        (*m_exp)[n]->GetFacePhysVals(e, elmtToTrace[n][e],
                                                     field + phys_offset,
                                                     e_tmp = Bwd + offset2,
                                                     it2->second.second);
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
                        npts = m_bndCondExpansions[n]->GetExp(e)->GetTotPoints();
                        id1  = m_bndCondExpansions[n]->GetPhys_Offset(e);
                        id2  = m_trace->GetPhys_Offset(m_traceMap->GetBndCondTraceToGlobalTraceMap(cnt+e));
                        Vmath::Vcopy(npts,&(m_bndCondExpansions[n]->GetPhys())[id1],1,&Bwd[id2],1);
                    }

                    cnt += e;
                }
                else if (m_bndConditions[n]->GetBoundaryConditionType() == SpatialDomains::eNeumann || 
                         m_bndConditions[n]->GetBoundaryConditionType() == SpatialDomains::eRobin)
                {
                    for(e = 0; e < m_bndCondExpansions[n]->GetExpSize(); ++e)
                    {
                        npts = m_bndCondExpansions[n]->GetExp(e)->GetNumPoints(0)*
                               m_bndCondExpansions[n]->GetExp(e)->GetNumPoints(1);
                        id1  = m_bndCondExpansions[n]->GetPhys_Offset(e);
                        id2  = m_trace->GetPhys_Offset(m_traceMap->GetBndCondTraceToGlobalTraceMap(cnt+e));

                        ASSERTL0((m_bndCondExpansions[n]->GetPhys())[id1] == 0.0,
                                 "method not set up for non-zero Neumann boundary condition");
                        
                        Vmath::Vcopy(npts,&Fwd[id2],1,&Bwd[id2],1);
                    }

                    cnt += e;
                }
                else
                {
                    ASSERTL0(false, "Method only set up for Dirichlet, Neumann and Robin conditions.");
                }
            }
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
        
        /// Note this routine changes m_trace->m_coeffs space;
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
                for(e = 0; e < (*m_exp)[n]->GetNfaces(); ++e)
                {
                    t_offset = m_trace->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());
                    (*m_exp)[n]->AddFaceNormBoundaryInt(e,elmtToTrace[n][e],
                                                        Fn + t_offset,
                                                        e_outarray = outarray+offset);
                }
            }
        }
    } // end of namespace
} // end of namespace
