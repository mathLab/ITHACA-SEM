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
        DisContField3D::DisContField3D( SpatialDomains::MeshGraph3D &graph3D,
                                        const GlobalSysSolnType solnType,
                                        bool SetUpJustDG) :
            ExpList3D(graph3D),
            m_numDirBndCondExpansions(0),
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
        DisContField3D::DisContField3D( SpatialDomains::MeshGraph3D &graph3D,
                                        SpatialDomains::BoundaryConditions &bcs,
                                        const int bc_loc,
                                        const GlobalSysSolnType solnType,
                                        bool SetUpJustDG) :
            ExpList3D(graph3D),
            m_numDirBndCondExpansions(0),
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
        DisContField3D::DisContField3D( SpatialDomains::MeshGraph3D &graph3D,
                                        SpatialDomains::BoundaryConditions &bcs,
                                        const std::string variable,
                                        const GlobalSysSolnType solnType,
                                        bool SetUpJustDG) :
            ExpList3D(graph3D),
            m_numDirBndCondExpansions(0),
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

            m_numDirBndCondExpansions = cnt2;
            m_bndCondExpansions  = Array<OneD,MultiRegions::ExpList2DSharedPtr>(cnt1);
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
                        Array<OneD, ExpList2DSharedPtr> &bndCondExpansions,
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
                                        ::AllocateSharedPtr(*(bregions[i]),
                                                            graph3D);
                    bndCondExpansions[cnt]  = locExpList;
                    bndConditions[cnt++]    = locBCond;
                } // end if Dirichlet
            }
            // then, list the other (non-periodic) boundaries
            for(i = 0; i < nbnd; ++i)
            {
                locBCond = (*(bconditions[i]))[variable];
                if(locBCond->GetBoundaryConditionType()
                                        == SpatialDomains::eNeumann)
                {
                    locExpList = MemoryManager<MultiRegions::ExpList2D>
                                        ::AllocateSharedPtr(*(bregions[i]),
                                                            graph3D);
                    bndCondExpansions[cnt]  = locExpList;
                    bndConditions[cnt++]    = locBCond;
                }
                else if((locBCond->GetBoundaryConditionType()
                            != SpatialDomains::eDirichlet)
                        && (locBCond->GetBoundaryConditionType()
                            != SpatialDomains::ePeriodic))
                {
                    ASSERTL0(false,"This type of BC not implemented yet");
                }
            }
        }


        /**
         * @param   time        The time at which the boundary conditions
         *                      should be evaluated.
         * @param   bndCondExpansions   List of boundary conditions.
         * @param   bndConditions   Information about the boundary conditions.
         */
        void DisContField3D::EvaluateBoundaryConditions(const NekDouble time)
        {
            int i,j;
            int npoints;
            int nbnd = m_bndCondExpansions.num_elements();
            MultiRegions::ExpList2DSharedPtr locExpList;

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
                                    ->m_DirichletCondition)
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
                                    ->m_NeumannCondition)
                                        .Evaluate(x0[j],x1[j],x2[j],time);
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


        void DisContField3D::v_EvaluateBoundaryConditions(const NekDouble time)
        {
            EvaluateBoundaryConditions(time);
        }


        const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>
                                        &DisContField3D::v_GetBndConditions()
        {
            return m_bndConditions;
        }

    } // end of namespace
} // end of namespace
