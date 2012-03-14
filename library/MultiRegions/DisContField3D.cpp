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
        DisContField3D::DisContField3D( const LibUtilities::SessionReaderSharedPtr &pSession,
                                        const SpatialDomains::MeshGraphSharedPtr &graph3D,
                                        const std::string &variable,
                                        const bool SetUpJustDG) :
            ExpList3D(pSession,graph3D),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            SpatialDomains::BoundaryConditions bcs(m_session, graph3D);
            
            GenerateBoundaryConditionExpansion(graph3D,bcs,variable);
            EvaluateBoundaryConditions();
            ApplyGeomInfo();
            
            if(SetUpJustDG)
            {
                // Set up matrix map
                m_globalBndMat = MemoryManager<GlobalLinSysMap>::AllocateSharedPtr();
                map<int,int> periodicEdges;
                map<int,int> periodicVertices;
                map<int,int> periodicFaces;
                GetPeriodicFaces(graph3D, bcs, variable,
                                 periodicVertices, periodicEdges, periodicFaces);
                
                //ASSERTL0(false, "DisContField3D Constructor needs implementation.");
                
                // Set up Trace space
                bool UseGenSegExp = true;
                m_trace = MemoryManager<ExpList2D>
                    ::AllocateSharedPtr(m_bndCondExpansions, m_bndConditions,
                                        *m_exp,graph3D, periodicFaces, UseGenSegExp);
                
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
                
                m_traceMap = MemoryManager<LocalToGlobalDGMap>::AllocateSharedPtr(
                    m_session,graph3D,m_trace,*this,m_bndCondExpansions,
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
                   // Set up matrix map
                   m_globalBndMat = MemoryManager<GlobalLinSysMap>::AllocateSharedPtr();
                   map<int,int> periodicEdges;
                   map<int,int> periodicVertices;
                   map<int,int> periodicFaces;
                   GetPeriodicFaces(graph3D, bcs, variable,
                                    periodicVertices, periodicEdges, periodicFaces);
                   
                   // Set up Trace space
                   bool UseGenSegExp = true;
                   m_trace = MemoryManager<ExpList2D>
                       ::AllocateSharedPtr(m_bndCondExpansions, m_bndConditions,
                                           *m_exp,graph3D, periodicFaces, UseGenSegExp);
                   
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
                   
                   m_traceMap = MemoryManager<LocalToGlobalDGMap>::AllocateSharedPtr(
                       m_session,graph3D,m_trace,*this, m_bndCondExpansions,
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
            const SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();

            int nbnd = bregions.size();

            // count the number of non-periodic boundary regions
            for(int i = 0; i < nbnd; ++i)
            {
                SpatialDomains::BoundaryConditionShPtr boundaryCondition = GetBoundaryCondition(bconditions, i, variable);
                if( boundaryCondition->GetBoundaryConditionType() != SpatialDomains::ePeriodic )
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
            int vSame = (returnval?1:0);
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
        //TODO: implement
        void DisContField3D::GetPeriodicFaces(
                    const SpatialDomains::MeshGraphSharedPtr &graph3D,
                    const SpatialDomains::BoundaryConditions &bcs,
                    const std::string &variable,
                    map<int,int>& periodicVertices,
                    map<int,int>& periodicEdges,
                    map<int,int>& periodicFaces)
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
                    ASSERTL0(false,"this method needs sorting");
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
                                if(!(segmentGeom1 = boost::dynamic_pointer_cast<
                                        SpatialDomains::SegGeom>((*comp1)[k]))||
                                   !(segmentGeom2 = boost::dynamic_pointer_cast<
                                        SpatialDomains::SegGeom>((*comp2)[k])))
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
//                                 element1 = graph3D.GetElementsFromEdge(segmentGeom1);
//                                 element2 = graph3D.GetElementsFromEdge(segmentGeom2);

                                ASSERTL0(element1->size()==1,
                                         "The periodic boundaries belong to "
                                         "more than one element of the mesh");
                                ASSERTL0(element2->size()==1,
                                         "The periodic boundaries belong to "
                                         "more than one element of the mesh");

                                orient1 = (boost::dynamic_pointer_cast<
                                            SpatialDomains::Geometry2D>(
                                                (*element1)[0]->m_Element))
                                            ->GetEorient(
                                                (*element1)[0]->m_EdgeIndx);
                                orient2 = (boost::dynamic_pointer_cast<
                                            SpatialDomains::Geometry2D>(
                                                (*element2)[0]->m_Element))
                                            ->GetEorient(
                                                (*element2)[0]->m_EdgeIndx);

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
        map<int, RobinBCInfoSharedPtr> DisContField3D::GetRobinBCInfo(void)
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

        map<int, RobinBCInfoSharedPtr> DisContField3D::v_GetRobinBCInfo()
        {
            return GetRobinBCInfo();
        }

        void DisContField3D::GetFwdBwdTracePhys(Array<OneD, NekDouble> &Fwd,
                                                Array<OneD, NekDouble> &Bwd)
        {
            GetFwdBwdTracePhys(m_phys,Fwd,Bwd);
        }
        
        void DisContField3D::GetFwdBwdTracePhys(
            const Array<OneD, const NekDouble> &field,
                  Array<OneD,       NekDouble> &Fwd,
                  Array<OneD,       NekDouble> &Bwd)
        {
            // Loop over elements and collect forward and backward expansions.
            int nexp = GetExpSize();
            int cnt,n,e,npts,offset, phys_offset;
            Array<OneD,NekDouble> e_tmp;
            
            Array<OneD, Array<OneD, StdRegions::StdExpansion2DSharedPtr> >
                elmtToTrace = m_traceMap->GetElmtToFace();
            
            // Zero vectors.
            Vmath::Zero(Fwd.num_elements(),Fwd,1);
            Vmath::Zero(Bwd.num_elements(),Bwd,1);
            
            Array<OneD, const NekDouble> a_tmp;
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
                        fwd = traceEl->GetLeftAdjacentElementExp() == 
                            (*m_exp)[n];
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

        void DisContField3D::ExtractTracePhys()
        {
            ExtractTracePhys(m_trace->UpdatePhys());
        }

        void DisContField3D::ExtractTracePhys(Array<OneD, NekDouble> &outarray)
        {
            ASSERTL1(m_physState == true,
                     "Field is not in physical space.");
            
            ExtractTracePhys(m_phys, outarray);
        }

        void DisContField3D::ExtractTracePhys(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray)
        {
            // Loop over elemente and collect forward expansion
            int nexp = GetExpSize();
            int n,e,offset,phys_offset;
            Array<OneD,NekDouble> e_tmp;
            Array<OneD, Array<OneD, StdRegions::StdExpansion2DSharedPtr> >
                elmtToTrace = m_traceMap->GetElmtToFace();

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
        void DisContField3D::AddTraceIntegral(
            const Array<OneD, const NekDouble> &Fx,
            const Array<OneD, const NekDouble> &Fy,
            const Array<OneD, const NekDouble> &Fz,
                  Array<OneD,       NekDouble> &outarray)
        {
            int e,n,offset, t_offset;
            Array<OneD, NekDouble> e_outarray;
            Array<OneD, Array<OneD, StdRegions::StdExpansion2DSharedPtr> >
                elmtToTrace = m_traceMap->GetElmtToFace();
            
            for(n = 0; n < GetExpSize(); ++n)
            {
                offset = GetCoeff_Offset(n);
                for(e = 0; e < (*m_exp)[n]->GetNfaces(); ++e)
                {
                    t_offset = GetTrace3D()->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());
                    (*m_exp)[n]->AddFaceNormBoundaryInt(e,elmtToTrace[n][e],
                                                        Fx + t_offset,
                                                        Fy + t_offset,
                                                        Fz + t_offset,
                                                        e_outarray = outarray+offset);
                }
            }
        }

        /// Note this routine changes m_trace->m_coeffs space;
        void DisContField3D::AddTraceIntegral(
            const Array<OneD, const NekDouble> &Fn,
                  Array<OneD,       NekDouble> &outarray)
        {
            int e,n,offset, t_offset;
            Array<OneD, NekDouble> e_outarray;
            Array<OneD, Array<OneD, StdRegions::StdExpansion2DSharedPtr> >
                elmtToTrace = m_traceMap->GetElmtToFace();

            for(n = 0; n < GetExpSize(); ++n)
            {
                offset = GetCoeff_Offset(n);
                for(e = 0; e < (*m_exp)[n]->GetNfaces(); ++e)
                {
                    t_offset = GetTrace3D()->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());
                    (*m_exp)[n]->AddFaceNormBoundaryInt(e,elmtToTrace[n][e],
                                                        Fn + t_offset,
                                                        e_outarray = outarray+offset);
                }
            }
        }
    } // end of namespace
} // end of namespace
