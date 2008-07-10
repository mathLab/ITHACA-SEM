///////////////////////////////////////////////////////////////////////////////
//
// File ContField2D.cpp
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
// Description: Field definition for 2D domain with boundary conditions
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ContField2D.h>

namespace Nektar
{
    namespace MultiRegions
    {

        ContField2D::ContField2D(void):
            ContExpList2D(),
            m_bndCondExpansions(),
            m_bndConditions()
        {
        }

        ContField2D::ContField2D(const ContField2D &In):
            ContExpList2D(In),
            m_bndCondExpansions(In.m_bndCondExpansions),
            m_bndConditions(In.m_bndConditions)
        {
        }

        ContField2D::ContField2D(SpatialDomains::MeshGraph2D &graph2D,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const int bc_loc):
            ContExpList2D(graph2D,false)  ,
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph2D,bcs,bcs.GetVariable(bc_loc));
            EvaluateBoundaryConditions();

            map<int,int> periodicEdges;
            map<int,int> periodicVertices;
            GetPeriodicEdges(graph2D,bcs,bcs.GetVariable(bc_loc),periodicVertices,periodicEdges);

            m_locToGloMap = MemoryManager<LocalToGlobalMap2D>::AllocateSharedPtr(m_ncoeffs,*m_exp,
                                                                                 m_bndCondExpansions,
                                                                                 m_bndConditions,
                                                                                 periodicVertices,
                                                                                 periodicEdges);
	    
	    m_contNcoeffs = m_locToGloMap->GetTotGloDofs();
	    m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }

        ContField2D::ContField2D(SpatialDomains::MeshGraph2D &graph2D,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const std::string variable):
            ContExpList2D(graph2D,false)  ,
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph2D,bcs,variable);
            EvaluateBoundaryConditions();

            map<int,int> periodicEdges;
            map<int,int> periodicVertices;
            GetPeriodicEdges(graph2D,bcs,variable,periodicVertices,periodicEdges);

            m_locToGloMap = MemoryManager<LocalToGlobalMap2D>::AllocateSharedPtr(m_ncoeffs,*m_exp,
                                                                                 m_bndCondExpansions,
                                                                                 m_bndConditions,
                                                                                 periodicVertices,
                                                                                 periodicEdges);
	    
	    m_contNcoeffs = m_locToGloMap->GetTotGloDofs();
	    m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }

        ContField2D::ContField2D(const LibUtilities::BasisKey &TriBa, 
                                 const LibUtilities::BasisKey &TriBb, 
                                 const LibUtilities::BasisKey &QuadBa, 
                                 const LibUtilities::BasisKey &QuadBb,
                                 SpatialDomains::MeshGraph2D &graph2D,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const int bc_loc,
                                 const LibUtilities::PointsType TriNb):
            ContExpList2D(TriBa,TriBb,QuadBa,QuadBb,graph2D,TriNb,false)  ,
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph2D,bcs,bcs.GetVariable(bc_loc));
            EvaluateBoundaryConditions();

            map<int,int> periodicEdges;
            map<int,int> periodicVertices;
            GetPeriodicEdges(graph2D,bcs,bcs.GetVariable(bc_loc),periodicVertices,periodicEdges);

            m_locToGloMap = MemoryManager<LocalToGlobalMap2D>::AllocateSharedPtr(m_ncoeffs,*m_exp,
                                                                                 m_bndCondExpansions,
                                                                                 m_bndConditions,
                                                                                 periodicVertices,
                                                                                 periodicEdges);
	    
	    m_contNcoeffs = m_locToGloMap->GetTotGloDofs();
	    m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }

        ContField2D::ContField2D(const LibUtilities::BasisKey &TriBa, 
                                 const LibUtilities::BasisKey &TriBb, 
                                 const LibUtilities::BasisKey &QuadBa, 
                                 const LibUtilities::BasisKey &QuadBb, 
                                 SpatialDomains::MeshGraph2D &graph2D,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const std::string variable,
                                 const LibUtilities::PointsType TriNb):
            ContExpList2D(TriBa,TriBb,QuadBa,QuadBb,graph2D,TriNb,false),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph2D,bcs,variable);
            EvaluateBoundaryConditions();

            map<int,int> periodicEdges;
            map<int,int> periodicVertices;
            GetPeriodicEdges(graph2D,bcs,variable,periodicVertices,periodicEdges);

            m_locToGloMap = MemoryManager<LocalToGlobalMap2D>::AllocateSharedPtr(m_ncoeffs,*m_exp,
                                                                                 m_bndCondExpansions,
                                                                                 m_bndConditions,
                                                                                 periodicVertices,
                                                                                 periodicEdges);
	    
	    m_contNcoeffs = m_locToGloMap->GetTotGloDofs();
	    m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }
        
        ContField2D::~ContField2D()
        {
        }
     
        void ContField2D::GenerateBoundaryConditionExpansion(SpatialDomains::MeshGraph2D &graph2D,
                                                             SpatialDomains::BoundaryConditions &bcs, 
                                                             const std::string variable)
        {
            int i;
            int cnt  = 0;
            
            SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();   
            
            MultiRegions::ExpList1DSharedPtr locExpList;  
            SpatialDomains::BoundaryConditionShPtr locBCond; 

            int nbnd = bregions.size();
            // count the number of non-periodic boundary regions
            for(i = 0; i < nbnd; ++i)
            {
                if( ((*(bconditions[i]))[variable])->GetBoundaryConditionType() != SpatialDomains::ePeriodic )
                {
                    cnt++;
                }
            }
            m_bndCondExpansions  = Array<OneD,MultiRegions::ExpList1DSharedPtr>(cnt);
            m_bndConditions      = Array<OneD,SpatialDomains::BoundaryConditionShPtr>(cnt);
          
            cnt=0;
            // list Dirichlet boundaries first
            for(i = 0; i < nbnd; ++i)
            {  
                locBCond = (*(bconditions[i]))[variable];  
                if(locBCond->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {                   
                    locExpList = MemoryManager<MultiRegions::ExpList1D>::AllocateSharedPtr(*(bregions[i]),graph2D);             
                    m_bndCondExpansions[cnt]  = locExpList;
                    m_bndConditions[cnt++]    = locBCond;
                } // end if Dirichlet
            }
            // then, list the other (non-periodic) boundaries
            for(i = 0; i < nbnd; ++i)
            {        
                locBCond = (*(bconditions[i]))[variable];  
                if(locBCond->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                {                    
                    locExpList = MemoryManager<MultiRegions::ExpList1D>::AllocateSharedPtr(*(bregions[i]),graph2D);
                    m_bndCondExpansions[cnt]  = locExpList;
                    m_bndConditions[cnt++]    = locBCond;
                }  
                else if((locBCond->GetBoundaryConditionType() != SpatialDomains::eDirichlet) && 
                        (locBCond->GetBoundaryConditionType() != SpatialDomains::ePeriodic))
                {
                    ASSERTL0(false,"This type of BC not implemented yet");
                }                  
            }
        }
        
        void ContField2D::GetPeriodicEdges(SpatialDomains::MeshGraph2D &graph2D,
                                           SpatialDomains::BoundaryConditions &bcs, 
                                           const std::string variable,
                                           map<int,int>& periodicVertices,
                                           map<int,int>& periodicEdges)
        {
            int i,j,k;
            
            SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();

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
                if(locBCond->GetBoundaryConditionType() == SpatialDomains::ePeriodic)
                {    
                    region1ID = i;
                    region2ID = (boost::static_pointer_cast<SpatialDomains::PeriodicBoundaryCondition>(locBCond))->m_ConnectedBoundaryRegion;

                    if(doneBndRegions.count(region1ID)==0)
                    {                    
                        ASSERTL0(bregions[region1ID]->size() == bregions[region2ID]->size(),
                                 "Size of the 2 periodic boundary regions should be equal");
                    
                        for(j = 0; j < bregions[region1ID]->size(); j++)
                        {
                            comp1 = (*(bregions[region1ID]))[j];
                            comp2 = (*(bregions[region2ID]))[j];
                            
                            ASSERTL0(comp1->size() == comp2->size(),
                                     "Size of the 2 periodic composites should be equal");
                            
                            for(k = 0; k < comp1->size(); k++)
                            {                                      
                                if(!(segmentGeom1 = boost::dynamic_pointer_cast<SpatialDomains::SegGeom>((*comp1)[k]))||
                                   !(segmentGeom2 = boost::dynamic_pointer_cast<SpatialDomains::SegGeom>((*comp2)[k])))
                                {
                                    ASSERTL0(false,"dynamic cast to a SegGeom failed");
                                } 

                                // Extract the periodic edges
                                periodicEdges[segmentGeom1->GetEid()] = segmentGeom2->GetEid();

                                // Extract the periodic vertices
                                element1 = graph2D.GetElementsFromEdge(segmentGeom1);
                                element2 = graph2D.GetElementsFromEdge(segmentGeom2);

                                ASSERTL0(element1->size()==1,"The periodic boundaries belong to more than one element of the mesh");
                                ASSERTL0(element2->size()==1,"The periodic boundaries belong to more than one element of the mesh");

                                orient1 = (boost::dynamic_pointer_cast<SpatialDomains::Geometry2D>((*element1)[0]->m_Element))->
                                    GetEorient((*element1)[0]->m_EdgeIndx);
                                orient2 = (boost::dynamic_pointer_cast<SpatialDomains::Geometry2D>((*element2)[0]->m_Element))->
                                    GetEorient((*element2)[0]->m_EdgeIndx);

                                if(orient1!=orient2)
                                {
                                    periodicVertices[segmentGeom1->GetVid(0)] = segmentGeom2->GetVid(0);
                                    periodicVertices[segmentGeom1->GetVid(1)] = segmentGeom2->GetVid(1);
                                }
                                else
                                {
                                    periodicVertices[segmentGeom1->GetVid(0)] = segmentGeom2->GetVid(1);
                                    periodicVertices[segmentGeom1->GetVid(1)] = segmentGeom2->GetVid(0);
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

        void ContField2D::EvaluateBoundaryConditions(const NekDouble time)
        {         
            int i,j;
            int npoints;
            int nbnd = m_bndCondExpansions.num_elements();
            MultiRegions::ExpList1DSharedPtr locExpList; 
            
            for(i = 0; i < nbnd; ++i)
            {                 
                locExpList = m_bndCondExpansions[i];                  
                npoints = locExpList->GetPointsTot();
                
                Array<OneD,NekDouble> x0(npoints,0.0);
                Array<OneD,NekDouble> x1(npoints,0.0);
                Array<OneD,NekDouble> x2(npoints,0.0);  
                
                locExpList->GetCoords(x0,x1,x2);

                if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {            
                    for(j = 0; j < npoints; j++)
                    {
                        (locExpList->UpdatePhys())[j] = (boost::static_pointer_cast<SpatialDomains::DirichletBoundaryCondition>(m_bndConditions[i])->
                                                         m_DirichletCondition).Evaluate(x0[j],x1[j],x2[j],time);
                    }
                    
                    locExpList->FwdTrans_BndConstrained(*locExpList);
                }
                else if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                {          
                    for(j = 0; j < npoints; j++)
                    {
                        (locExpList->UpdatePhys())[j] = (boost::static_pointer_cast<SpatialDomains::NeumannBoundaryCondition>(m_bndConditions[i])->
                                                         m_NeumannCondition).Evaluate(x0[j],x1[j],x2[j],time);
                    }

                    locExpList->IProductWRTBase(*locExpList); 
                }
                else
                {
                    ASSERTL0(false,"This type of BC not implemented yet");
                }
            }           
        }
        
        void ContField2D::FwdTrans(const ExpList &In)
        {
            GlobalLinSysKey key(StdRegions::eMass);
            GlobalSolve(key,In);

            m_transState = eContinuous;
            m_physState = false;
        }

        // Solve the helmholtz problem assuming that m_contCoeff vector 
        // contains an intial estimate for solution
        void ContField2D::HelmSolve(const ExpList &In, NekDouble lambda)
        {
            GlobalLinSysKey key(StdRegions::eHelmholtz,lambda);
            // Note -1.0 term necessary to invert forcing function to
            // be consistent with matrix definition
            GlobalSolve(key,In,-1.0);
        }

        void ContField2D::GlobalSolve(const GlobalLinSysKey &key, 
                                      const ExpList &Rhs, NekDouble ScaleForcing)
        {
            int i,j;
            int bndcnt=0;
            int NumDirBcs = m_locToGloMap->GetNumDirichletDofs();
            Array<OneD,NekDouble> sln;
            Array<OneD,NekDouble> init(m_contNcoeffs,0.0);
            Array<OneD,NekDouble> Dir_fce(m_contNcoeffs,0.0);

            //assume m_contCoeffs contains initial estimate
            // Set BCs in m_contCoeffs
            Blas::Dcopy(m_contNcoeffs, m_contCoeffs, 1, init, 1);
            for(i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                if(m_bndConditions[i]->GetBoundaryConditionType() != SpatialDomains::eDirichlet)
                {
                    break;
                }
                for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                {
                    init[m_locToGloMap->GetBndCondMap(bndcnt++)] = (m_bndCondExpansions[i]->GetCoeffs())[j];
                }
            }
            GeneralMatrixOp(key, init, Dir_fce);

            // Set up forcing function
            IProductWRTBase(Rhs);
            
            // apply scaling of forcing term; (typically used to negate helmholtz forcing);
            Vmath::Smul(m_contNcoeffs, ScaleForcing, m_contCoeffs, 1, m_contCoeffs, 1);

            // Forcing function with Dirichlet conditions 
            Vmath::Vsub(m_contNcoeffs, m_contCoeffs, 1,
                        Dir_fce, 1, m_contCoeffs, 1);

            // Forcing function with weak boundary conditions 
            for(i; i < m_bndCondExpansions.num_elements(); ++i)
            {
                for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                {
                    m_contCoeffs[m_locToGloMap->GetBndCondMap(bndcnt++)] +=  
                        (m_bndCondExpansions[i]->GetCoeffs())[j];
                }
            }

            if(m_contNcoeffs - NumDirBcs > 0)
            {

                GlobalLinSysSharedPtr LinSys = GetGlobalLinSys(key);
                       
                sln = m_contCoeffs+NumDirBcs;
                LinSys->Solve(sln,sln,*m_locToGloMap);
            }

            // Recover solution by addinig intial conditons
            Vmath::Zero(NumDirBcs, m_contCoeffs, 1);
            Vmath::Vadd(m_contNcoeffs, init, 1, m_contCoeffs, 1,
                        m_contCoeffs, 1);

            m_transState = eContinuous;
            m_physState = false;
        }

        GlobalLinSysSharedPtr ContField2D::GetGlobalLinSys(const GlobalLinSysKey &mkey) 
        {
            GlobalLinSysSharedPtr glo_matrix;
            GlobalLinSysMap::iterator matrixIter = m_globalMat->find(mkey);

            if(matrixIter == m_globalMat->end())
            {
                glo_matrix = GenGlobalLinSys(mkey,m_locToGloMap);
                (*m_globalMat)[mkey] = glo_matrix;
            }
            else
            {
                glo_matrix = matrixIter->second;
            }

            return glo_matrix;
        }


    } // end of namespace
} //end of namespace
