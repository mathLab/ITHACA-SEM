///////////////////////////////////////////////////////////////////////////////
//
// File ContField1D.cpp
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
// Description: Field definition for 1D domain with boundary conditions
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ContField1D.h>

namespace Nektar
{
    namespace MultiRegions
    {

        ContField1D::ContField1D(void):
            ContExpList1D(),
            m_bndCondExpansions(),
            m_bndConditions()
        {
        }

        ContField1D::ContField1D(const ContField1D &In):
            ContExpList1D(In),
            m_bndCondExpansions(In.m_bndCondExpansions),
            m_bndConditions(In.m_bndConditions)
        {
        }

        ContField1D::ContField1D(SpatialDomains::MeshGraph1D &graph1D,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const int bc_loc):
            ContExpList1D(graph1D,false),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph1D,bcs,bcs.GetVariable(bc_loc));
            EvaluateBoundaryConditions();

            map<int,int> periodicVertices;
            GetPeriodicVertices(graph1D,bcs,bcs.GetVariable(bc_loc),periodicVertices);

            m_locToGloMap = MemoryManager<LocalToGlobalMap1D>::AllocateSharedPtr(m_ncoeffs,*m_exp,
                                                                                 m_bndCondExpansions,
                                                                                 m_bndConditions,
                                                                                 periodicVertices);
	    
	    m_contNcoeffs = m_locToGloMap->GetTotGloDofs();
	    m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }

        ContField1D::ContField1D(SpatialDomains::MeshGraph1D &graph1D,
            SpatialDomains::BoundaryConditions &bcs, 
            const std::string variable):
            ContExpList1D(graph1D,false),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph1D,bcs,variable);
            EvaluateBoundaryConditions();

            map<int,int> periodicVertices;
            GetPeriodicVertices(graph1D,bcs,variable,periodicVertices);

            m_locToGloMap = MemoryManager<LocalToGlobalMap1D>::AllocateSharedPtr(m_ncoeffs,*m_exp,
                                                                                 m_bndCondExpansions,
                                                                                 m_bndConditions,
                                                                                 periodicVertices);
	    
	    m_contNcoeffs = m_locToGloMap->GetTotGloDofs();
	    m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }

        ContField1D::ContField1D(const LibUtilities::BasisKey &Ba, 
                                 const SpatialDomains::MeshGraph1D &graph1D,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const int bc_loc):
            ContExpList1D(Ba,graph1D,false),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph1D,bcs,bcs.GetVariable(bc_loc));
            EvaluateBoundaryConditions();

            map<int,int> periodicVertices;
            GetPeriodicVertices(graph1D,bcs,bcs.GetVariable(bc_loc),periodicVertices);

            m_locToGloMap = MemoryManager<LocalToGlobalMap1D>::AllocateSharedPtr(m_ncoeffs,*m_exp,
                                                                                 m_bndCondExpansions,
                                                                                 m_bndConditions,
                                                                                 periodicVertices);
	    
	    m_contNcoeffs = m_locToGloMap->GetTotGloDofs();
	    m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }

        ContField1D::ContField1D(const LibUtilities::BasisKey &Ba, 
                                 const SpatialDomains::MeshGraph1D &graph1D,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const std::string variable):
            ContExpList1D(Ba,graph1D,false),
            m_bndCondExpansions(), 
            m_bndConditions()         
        {
            GenerateBoundaryConditionExpansion(graph1D,bcs,variable);
            EvaluateBoundaryConditions();

            map<int,int> periodicVertices;
            GetPeriodicVertices(graph1D,bcs,variable,periodicVertices);

            m_locToGloMap = MemoryManager<LocalToGlobalMap1D>::AllocateSharedPtr(m_ncoeffs,*m_exp,
                                                                                 m_bndCondExpansions,
                                                                                 m_bndConditions,
                                                                                 periodicVertices);
	    
	    m_contNcoeffs = m_locToGloMap->GetTotGloDofs();
	    m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }

        ContField1D::~ContField1D()
        {
        }

        void ContField1D::GenerateBoundaryConditionExpansion(const SpatialDomains::MeshGraph1D &graph1D,
                                                             SpatialDomains::BoundaryConditions &bcs, 
                                                             const std::string variable)
        {
            int i,j,k;
            int cnt  = 0;
            
            SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();   
            
            LocalRegions::PointExpSharedPtr          locPointExp;
            SpatialDomains::BoundaryConditionShPtr   locBCond; 
            SpatialDomains::VertexComponentSharedPtr vert;

            int nbnd = bregions.size(); 
            // count the number of non-periodic boundary points
            for(i = 0; i < nbnd; ++i)
            {   
                if( ((*(bconditions[i]))[variable])->GetBoundaryConditionType() != SpatialDomains::ePeriodic )
                {
                    for(j = 0; j < bregions[i]->size(); j++)
                    {
                        cnt += (*bregions[i])[j]->size();
                    } 
                }
            }
                       
            m_bndCondExpansions = Array<OneD,LocalRegions::PointExpSharedPtr>(cnt);
            m_bndConditions     = Array<OneD,SpatialDomains::BoundaryConditionShPtr>(cnt);
            
            cnt=0;
            // list Dirichlet boundaries first
            for(i = 0; i < nbnd; ++i)
            {  
                locBCond = (*(bconditions[i]))[variable];  
                if(locBCond->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {       
                    for(j = 0; j < bregions[i]->size(); j++)
                    {
                        for(k = 0; k < ((*bregions[i])[j])->size(); k++)
                        {
                            if(vert = boost::dynamic_pointer_cast<SpatialDomains::VertexComponent>((*(*bregions[i])[j])[k]))
                            {
                                locPointExp = MemoryManager<LocalRegions::PointExp>::AllocateSharedPtr(vert);
                                m_bndCondExpansions[cnt]  = locPointExp;
                                m_bndConditions[cnt++]    = locBCond;
                            }
                            else
                            {
                                ASSERTL0(false,"dynamic cast to a vertex failed");
                            }
                        }
                    }
                } // end if Dirichlet
            }

            // then, list the other (non-periodic) boundaries
            for(i = 0; i < nbnd; ++i)
            {        
                locBCond = (*(bconditions[i]))[variable];  
                if(locBCond->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                {    
                    for(j = 0; j < bregions[i]->size(); j++)
                    {
                        for(k = 0; k < ((*bregions[i])[j])->size(); k++)
                        {     
                            if(vert = boost::dynamic_pointer_cast<SpatialDomains::VertexComponent>((*(*bregions[i])[j])[k]))
                            {
                                locPointExp = MemoryManager<LocalRegions::PointExp>::AllocateSharedPtr(vert);
                                m_bndCondExpansions[cnt]  = locPointExp;
                                m_bndConditions[cnt++]    = locBCond;
                            }
                            else
                            {
                                ASSERTL0(false,"dynamic cast to a vertex failed");
                            }            
                        }
                    }
                }    
                else if((locBCond->GetBoundaryConditionType() != SpatialDomains::eDirichlet) && 
                        (locBCond->GetBoundaryConditionType() != SpatialDomains::ePeriodic))
                {
                    ASSERTL0(false,"This type of BC not implemented yet");
                }                  
            }
        }

        void ContField1D::GetPeriodicVertices(const SpatialDomains::MeshGraph1D &graph1D,
                                              SpatialDomains::BoundaryConditions &bcs, 
                                              const std::string variable,
                                              map<int,int>& periodicVertices)
        {

            int i,j,k;
            
            SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();

            int region1ID;
            int region2ID;

            SpatialDomains::Composite comp1;
            SpatialDomains::Composite comp2;

            SpatialDomains::VertexComponentSharedPtr vert1;
            SpatialDomains::VertexComponentSharedPtr vert2;
            
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
                                if(!(vert1 = boost::dynamic_pointer_cast<SpatialDomains::VertexComponent>((*comp1)[k]))||
                                   !(vert2 = boost::dynamic_pointer_cast<SpatialDomains::VertexComponent>((*comp2)[k])))
                                {
                                    ASSERTL0(false,"dynamic cast to a VertexComponent failed");
                                } 

                                // Extract the periodic vertices
                                periodicVertices[vert1->GetVid()] = vert2->GetVid();
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

        void ContField1D::EvaluateBoundaryConditions(const NekDouble time)
        {
            int i;

            NekDouble x0;
            NekDouble x1;
            NekDouble x2;
            
            for(i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                m_bndCondExpansions[i]->GetCoords(x0,x1,x2);

                if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                { 
                    m_bndCondExpansions[i]->SetValue((boost::static_pointer_cast<SpatialDomains::DirichletBoundaryCondition>(m_bndConditions[i])->
                                                      m_DirichletCondition).Evaluate(x0,x1,x2,time));
                }
                else if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                { 
                    m_bndCondExpansions[i]->SetValue((boost::static_pointer_cast<SpatialDomains::NeumannBoundaryCondition>(m_bndConditions[i])->
                                                      m_NeumannCondition).Evaluate(x0,x1,x2,time));
                }
                else
                {
                    ASSERTL0(false,"This type of BC not implemented yet");
                }
            }                
        }

        void ContField1D::FwdTrans(const ExpList &In)
        {
            GlobalLinSysKey key(StdRegions::eMass);
            GlobalSolve(key,In);

            m_transState = eContinuous;
            m_physState = false;
        }

        // Solve the helmholtz problem assuming that m_contCoeff vector 
        // contains an intial estimate for solution
        void ContField1D::HelmSolve(const ExpList &In, NekDouble lambda)
        {
            ASSERTL0(lambda != NekUnsetDouble,"Helmholtz constant is not specified");

            GlobalLinSysKey key(StdRegions::eHelmholtz,lambda);
            // Note -1.0 term necessary to invert forcing function to
            // be consistent with matrix definition
            GlobalSolve(key,In,-1.0);
        }

        void ContField1D::GlobalSolve(const GlobalLinSysKey &key, 
                                      const ExpList &Rhs, 
                                      NekDouble ScaleForcing)
        {
            int i;
            int NumDirBcs = m_locToGloMap->GetNumDirichletDofs();
            Array<OneD,NekDouble> sln;
            Array<OneD,NekDouble> init(m_contNcoeffs,0.0);
            Array<OneD,NekDouble> Dir_fce(m_contNcoeffs,0.0);

            //assume m_contCoeffs contains initial estimate
            // Set BCs in m_contCoeffs
            Blas::Dcopy(m_contNcoeffs, m_contCoeffs, 1, init, 1);
            for(i = 0; i < NumDirBcs; ++i)
            {
                init[i] = m_bndCondExpansions[i]->GetValue();
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
            for(i = 0; i < m_bndCondExpansions.num_elements()-NumDirBcs; ++i)
            {
                m_contCoeffs[m_locToGloMap->GetBndCondMap(i+NumDirBcs)] +=  
                    m_bndCondExpansions[i+NumDirBcs]->GetValue();
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

        GlobalLinSysSharedPtr ContField1D::GetGlobalLinSys(const GlobalLinSysKey &mkey) 
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
