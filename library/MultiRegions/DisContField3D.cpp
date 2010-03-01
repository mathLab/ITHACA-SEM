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
//                map<int,int> periodicEdges;
//                map<int,int> periodicVertices;
//                map<int,int> periodicFaces;
//                GetPeriodicFaces(graph3D,bcs,bcs.GetVariable(bc_loc),
//                                 periodicVertices,periodicEdges,periodicFaces);

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
         * 
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
         * 
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
/*        DisContField3D::DisContField3D( SpatialDomains::MeshGraph3D &graph3D,
                                        SpatialDomains::BoundaryConditions &bcs,
                                        const GlobalSysSolnType solnType,
                                        const bool constructMap) :
            ExpList3D(graph3D),
            m_numDirBndCondExpansions(0),
            m_bndCondExpansions(),
            m_bndConditions()
        {
        }
*/

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


    } // end of namespace
} // end of namespace
