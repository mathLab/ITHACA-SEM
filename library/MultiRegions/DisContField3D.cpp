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
         * @todo Implement 3D trace space.
         */
        DisContField3D::DisContField3D( LibUtilities::CommSharedPtr &pComm,
                                        SpatialDomains::MeshGraph3D &graph3D,
                                        const GlobalSysSolnType solnType,
                                        bool SetUpJustDG) :
            ExpList3D(pComm,graph3D),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            ApplyGeomInfo(graph3D);

            if(SetUpJustDG)
            {
                ASSERTL0(false, "DisContField3D Constructor needs trace implementation.");
                // Set up matrix map
                m_globalBndMat = MemoryManager<GlobalLinSysMap>
                                                    ::AllocateSharedPtr();
/*
                map<int,int> periodicEdges;
                map<int,int> periodicVertices;
                map<int,int> periodicFaces;
                GetPeriodicFaces(graph3D,bcs,bcs.GetVariable(bc_loc),
                                 periodicVertices,periodicEdges,periodicFaces);

                // Set up Trace space
                bool UseGenSegExp = true;
                m_trace = MemoryManager<ExpList1D>
                    ::AllocateSharedPtr(m_bndCondExpansions, m_bndConditions,
                                *m_exp,graph2D, periodicEdges, UseGenSegExp);

                m_traceMap = MemoryManager<LocalToGlobalDGMap>::
                    AllocateSharedPtr(graph2D,m_trace,m_exp,solnType,
                                      m_bndCondExpansions,m_bndConditions,
                                      periodicEdges);
*/            }

        }


        /**
         * @todo Implement 3D trace space.
         */
        DisContField3D::DisContField3D( LibUtilities::CommSharedPtr &pComm,
                                        SpatialDomains::MeshGraph3D &graph3D,
                                        SpatialDomains::BoundaryConditions &bcs,
                                        const int bc_loc,
                                        const GlobalSysSolnType solnType,
                                        bool SetUpJustDG) :
            ExpList3D(pComm,graph3D),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph3D,bcs,
                                               bcs.GetVariable(bc_loc));
            EvaluateBoundaryConditions();
            ApplyGeomInfo(graph3D);

            if(SetUpJustDG)
            {
                ASSERTL0(false, "DisContField3D Constructor needs trace implementation.");
                // Set up matrix map
                m_globalBndMat = MemoryManager<GlobalLinSysMap>
                                                    ::AllocateSharedPtr();
                map<int,int> periodicEdges;
                map<int,int> periodicVertices;
                map<int,int> periodicFaces;
                GetPeriodicFaces(graph3D,bcs,bcs.GetVariable(bc_loc),
                                 periodicVertices,periodicEdges,periodicFaces);

/*                // Set up Trace space
                bool UseGenSegExp = true;
                m_trace = MemoryManager<ExpList1D>
                    ::AllocateSharedPtr(m_bndCondExpansions, m_bndConditions,
                                *m_exp,graph2D, periodicEdges, UseGenSegExp);

                m_traceMap = MemoryManager<LocalToGlobalDGMap>::
                    AllocateSharedPtr(graph2D,m_trace,m_exp,solnType,
                                      m_bndCondExpansions,m_bndConditions,
                                      periodicEdges);
*/            }

        }


        /**
         * @todo Implement 3D trace space.
         */
        DisContField3D::DisContField3D( LibUtilities::CommSharedPtr &pComm,
                                        SpatialDomains::MeshGraph3D &graph3D,
                                        SpatialDomains::BoundaryConditions &bcs,
                                        const std::string variable,
                                        const GlobalSysSolnType solnType,
                                        bool SetUpJustDG) :
            ExpList3D(pComm,graph3D),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph3D,bcs,variable);
            EvaluateBoundaryConditions();
            ApplyGeomInfo(graph3D);

            if(SetUpJustDG)
            {
                // Set up matrix map
                m_globalBndMat = MemoryManager<GlobalLinSysMap>
                                                    ::AllocateSharedPtr();
                map<int,int> periodicEdges;
                map<int,int> periodicVertices;
                map<int,int> periodicFaces;
                GetPeriodicFaces(graph3D,bcs,variable,
                                 periodicVertices,periodicEdges,periodicFaces);

                ASSERTL0(false, "DisContField3D Constructor needs implementation.");

                // Set up Trace space
/*                bool UseGenSegExp = true;
                m_trace = MemoryManager<ExpList1D>
                    ::AllocateSharedPtr(m_bndCondExpansions, m_bndConditions,
                                *m_exp,graph2D, periodicEdges, UseGenSegExp);

                m_traceMap = MemoryManager<LocalToGlobalDGMap>::
                    AllocateSharedPtr(graph2D,m_trace,m_exp,solnType,
                                      m_bndCondExpansions,m_bndConditions,
                                      periodicEdges);
*/            }

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
                    SpatialDomains::MeshGraph3D &graph3D,
                    SpatialDomains::BoundaryConditions &bcs,
                    const std::string variable)
        {
            int cnt1  = 0;
            int cnt2  = 0;
            SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();

            int nbnd = bregions.size();

            // count the number of non-periodic boundary regions
            for(int i = 0; i < nbnd; ++i)
            {
                if( ((*(bconditions[i]))[variable])->GetBoundaryConditionType() != SpatialDomains::ePeriodic )
                {
                    cnt1++;
                    if( ((*(bconditions[i]))[variable])->GetBoundaryConditionType() == SpatialDomains::eDirichlet )
                    {
                        cnt2++;
                    }
                }
            }

            m_bndCondExpansions  = Array<OneD,MultiRegions::ExpListSharedPtr>(cnt1);
            m_bndConditions      = Array<OneD,SpatialDomains::BoundaryConditionShPtr>(cnt1);

            SetBoundaryConditionExpansion(graph3D,bcs,variable,m_bndCondExpansions,m_bndConditions);
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
                        SpatialDomains::MeshGraph3D &graph3D,
                        SpatialDomains::BoundaryConditions &bcs,
                        const std::string variable,
                        Array<OneD, ExpListSharedPtr> &bndCondExpansions,
                        Array<OneD, SpatialDomains::BoundaryConditionShPtr>
                        &bndConditions)
        {
            int i;
            int cnt  = 0;

            SpatialDomains::BoundaryRegionCollection &bregions
                                        = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions
                                        = bcs.GetBoundaryConditions();

            MultiRegions::ExpList2DSharedPtr       locExpList;
            SpatialDomains::BoundaryConditionShPtr locBCond;

            int nbnd = bregions.size();

            cnt=0;
            // list Dirichlet boundaries first
            for(i = 0; i < nbnd; ++i)
            {
                locBCond = (*(bconditions[i]))[variable];
                if(locBCond->GetBoundaryConditionType()
                                        == SpatialDomains::eDirichlet)
                {
                    locExpList = MemoryManager<MultiRegions::ExpList2D>
                                        ::AllocateSharedPtr(m_comm,*(bregions[i]),
                                                            graph3D);
                    bndCondExpansions[cnt]  = locExpList;
                    bndConditions[cnt++]    = locBCond;
                } // end if Dirichlet
            }
            // then, list the other (non-periodic) boundaries
            for(i = 0; i < nbnd; ++i)
            {
                locBCond = (*(bconditions[i]))[variable];

                switch(locBCond->GetBoundaryConditionType())
                {
                case SpatialDomains::eNeumann:
                case SpatialDomains::eRobin:
                    {
                        locExpList = MemoryManager<MultiRegions::ExpList2D>
                            ::AllocateSharedPtr(m_comm,*(bregions[i]),
                                                graph3D);
                        bndCondExpansions[cnt]  = locExpList;
                        bndConditions[cnt++]    = locBCond;
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
                    SpatialDomains::MeshGraph3D &graph3D,
                    SpatialDomains::BoundaryConditions &bcs,
                    const std::string variable,
                    map<int,int>& periodicVertices,
                    map<int,int>& periodicEdges,
                    map<int,int>& periodicFaces)
        {
            int i,j,k;

            SpatialDomains::BoundaryRegionCollection &bregions
                                        = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions
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
                locBCond = (*(bconditions[i]))[variable];
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


        // Set up a list of element ids and edge ids that link to the
        // boundary conditions
        void DisContField3D::GetBoundaryToElmtMap(Array<OneD, int> &ElmtID, Array<OneD,int> &FaceID)
        {
            map<int, int> FaceGID;
            int i,n,id;
            int bid,cnt,Fid;
            int nbcs = 0;

            for(i = 0; i < m_bndConditions.num_elements(); ++i)
            {
                nbcs += m_bndCondExpansions[i]->GetExpSize();
            }

            // make sure arrays are of sufficient length
            if(ElmtID.num_elements() != nbcs)
            {
                ElmtID = Array<OneD, int>(nbcs,-1);
            }
            else
            {
                fill(ElmtID.get(), ElmtID.get()+nbcs, -1);
            }

            if(FaceID.num_elements() != nbcs)
            {
                FaceID = Array<OneD, int>(nbcs);
            }

            // setup map of all global ids along boundary
            for(cnt = n = 0; n < m_bndCondExpansions.num_elements(); ++n)
            {
                for(i = 0; i < m_bndCondExpansions[n]->GetExpSize(); ++i)
                {
                    Fid =  m_bndCondExpansions[n]->GetExp(i)->GetGeom2D()->GetFid();
                    FaceGID[Fid] = cnt++;
                }
            }


            // Loop over elements and find edges that match;
            for(cnt = n = 0; n < GetExpSize(); ++n)
            {
                for(i = 0; i < (*m_exp)[n]->GetNfaces(); ++i)
                {
                    id = (*m_exp)[n]->GetGeom3D()->GetFid(i);

                    if(FaceGID.count(id) > 0)
                    {
                        bid = FaceGID.find(id)->second;
                        ASSERTL1(ElmtID[bid] == -1,"Face already set");
                        ElmtID[bid] = n;
                        FaceID[bid] = i;
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
            GetBoundaryToElmtMap(ElmtID,FaceID);
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
                    for(j = 0; j < npoints; j++)
                    {
                        (locExpList->UpdatePhys())[j]
                            = (boost::static_pointer_cast<SpatialDomains
                                ::DirichletBoundaryCondition>(m_bndConditions[i])
                                    ->m_dirichletCondition)
                                        .Evaluate(x0[j],x1[j],x2[j],time);
                    }

                    locExpList->FwdTrans_BndConstrained(locExpList->GetPhys(),
                                                    locExpList->UpdateCoeffs());
                }
                else if(m_bndConditions[i]->GetBoundaryConditionType()
                                        == SpatialDomains::eNeumann)
                {
                    for(j = 0; j < npoints; j++)
                    {
                        (locExpList->UpdatePhys())[j]
                            = (boost::static_pointer_cast<SpatialDomains
                                ::NeumannBoundaryCondition>(m_bndConditions[i])
                                    ->m_neumannCondition)
                                        .Evaluate(x0[j],x1[j],x2[j],time);
                    }

                    locExpList->IProductWRTBase(locExpList->GetPhys(),
                                                locExpList->UpdateCoeffs());
                }
                else if(m_bndConditions[i]->GetBoundaryConditionType()
                                        == SpatialDomains::eRobin)
                {
                    for(j = 0; j < npoints; j++)
                    {
                        (locExpList->UpdatePhys())[j]
                            = (boost::static_pointer_cast<SpatialDomains
                               ::RobinBoundaryCondition>(m_bndConditions[i])
                               ->m_robinFunction).Evaluate(x0[j],x1[j],x2[j],time);
                    }

                    locExpList->IProductWRTBase(locExpList->GetPhys(),
                                                locExpList->UpdateCoeffs());
                }
                else
                {
                    ASSERTL0(false,"This type of BC not implemented yet");
                }
            }
        }
    } // end of namespace
} // end of namespace
