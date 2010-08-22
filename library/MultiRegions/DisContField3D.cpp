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
                MultiRegions::ExpList2DSharedPtr locExpList;

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


        void DisContField3D::v_EvaluateBoundaryConditions(
                const NekDouble time)
        {
            EvaluateBoundaryConditions(time);
        }

    } // end of namespace
} // end of namespace
