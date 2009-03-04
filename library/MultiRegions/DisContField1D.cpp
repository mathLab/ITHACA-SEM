///////////////////////////////////////////////////////////////////////////////
//
// File DisContField1D.cpp
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
// Description: Field definition for 1D domain with boundary
// conditions using LDG-H
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/DisContField1D.h>

namespace Nektar
{
    namespace MultiRegions
    {

        DisContField1D::DisContField1D(void)
        {
        }
        
        // Copy constructor
        DisContField1D::DisContField1D(const DisContField1D &In):
            ExpList1D(In),
            m_bndCondExpansions(In.m_bndCondExpansions),
            m_bndConditions(In.m_bndConditions),
            m_globalBndMat(In.m_globalBndMat),
            m_trace(In.m_trace),
            m_traceMap(In.m_traceMap)
        {
        }

        DisContField1D::DisContField1D(SpatialDomains::MeshGraph1D &graph1D,
            SpatialDomains::BoundaryConditions &bcs, const int bc_loc):
            ExpList1D(graph1D)     
        {
            GenerateBoundaryConditionExpansion(graph1D,bcs,bcs.GetVariable(bc_loc));
            EvaluateBoundaryConditions();

            map<int,int> periodicVertices;
            GetPeriodicVertices(graph1D,bcs,bcs.GetVariable(bc_loc),periodicVertices);

            m_globalBndMat   = MemoryManager<GlobalLinSysMap>::AllocateSharedPtr();

            //GenerateFieldBnd1D(bcs,bcs.GetVariable(bc_loc));

            m_traceMap = MemoryManager<LocalToGlobalDGMap>::AllocateSharedPtr(graph1D,m_exp,m_bndCondExpansions,m_bndConditions);

            m_trace = Array<OneD, NekDouble>(m_traceMap->GetNumLocalBndCoeffs());
        }

        DisContField1D::DisContField1D(SpatialDomains::MeshGraph1D &graph1D,
            SpatialDomains::BoundaryConditions &bcs, 
            const std::string variable):
            ExpList1D(graph1D)  
        {
            GenerateBoundaryConditionExpansion(graph1D,bcs,variable);
            EvaluateBoundaryConditions();
            
            map<int,int> periodicVertices;
            GetPeriodicVertices(graph1D,bcs,variable,periodicVertices);

            m_globalBndMat = MemoryManager<GlobalLinSysMap>::AllocateSharedPtr();

            //GenerateFieldBnd1D(bcs,variable);
            
            m_traceMap = MemoryManager<LocalToGlobalDGMap>::AllocateSharedPtr(graph1D,m_exp,m_bndCondExpansions,m_bndConditions);

            m_trace = Array<OneD, NekDouble>(m_traceMap->GetNumLocalBndCoeffs());
        }

        void DisContField1D::GenerateBoundaryConditionExpansion(const SpatialDomains::MeshGraph1D &graph1D,
                                                             SpatialDomains::BoundaryConditions &bcs, 
                                                             const std::string variable)
        {
            int i,j,k;
            int cnt  = 0;
            
            SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();   
            
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
            
            SetBoundaryConditionExpansion(graph1D,bcs,variable, 
                                           m_bndCondExpansions,
                                           m_bndConditions);
        }
#if 0  
        void DisContField1D::GenerateFieldBnd1D(SpatialDomains::BoundaryConditions &bcs,   const std::string variable)
        {
            int i,nbnd;
            int cnt=0;
            LocalRegions::PointExpSharedPtr  p_exp;
            SpatialDomains::BoundaryConditionShPtr locBCond;
            SpatialDomains::VertexComponentSharedPtr vert;
            SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();

            nbnd = bregions.size();

            // count the number of non-periodic boundary regions
            for(int i = 0; i < nbnd; ++i)
            {
                if( ((*(bconditions[i]))[variable])->GetBoundaryConditionType() != SpatialDomains::ePeriodic )
                {
                    for(j = 0; j < bregions[i]->size(); ++j)
                    {
                        cnt += (*bregions[i])[j]->size();
                    }
                }
            }
            
            m_bndCondExpansions = Array<OneD,LocalRegions::PointExpSharedPtr>(cnt);
            m_bndConditions     = Array<OneD,SpatialDomains::BoundaryConditionShPtr>(cnt);
            
            // Set up matrix map
            m_globalBndMat   = MemoryManager<GlobalLinSysMap>::AllocateSharedPtr();

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
                                Array<OneD,NekDouble> coords(3,0.0);
                                
                                p_exp = MemoryManager<LocalRegions::PointExp>::AllocateSharedPtr(vert);
                                vert->GetCoords(coords);
                                p_exp->SetValue(boost::static_pointer_cast<SpatialDomains::DirichletBoundaryCondition>(locBCond)->m_DirichletCondition.Evaluate(coords[0],coords[1],coords[2]));
                                
                                m_bndCondExpansions[cnt] = p_exp;
                                m_bndConditions[cnt++] = locBCond;
                            }
                            else
                            {
                                ASSERTL0(false,"dynamic cast to a vertex failed");
                            }
                }
            }

            for(i = 0; i < nbnd; ++i)
            {
                locBCond = (*(bconditions[i]))[variable];

                if(locBCond->GetBoundaryConditionType() != SpatialDomains::eDirichlet)
                {
                    if(vert = boost::dynamic_pointer_cast<SpatialDomains::VertexComponent>((*(*bregions[i])[0])[0]))
                    {
                        Array<OneD,NekDouble> coords(3,0.0);

                        p_exp = MemoryManager<LocalRegions::PointExp>::AllocateSharedPtr(vert);
                        vert->GetCoords(coords);

                        if(locBCond->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                        {
                            p_exp->SetValue(boost::static_pointer_cast<SpatialDomains::NeumannBoundaryCondition>(locBCond)
                                            ->m_NeumannCondition.Evaluate(coords[0],coords[1],coords[2]));
                            m_bndCondExpansions[cnt] = p_exp;
                            m_bndConditions[cnt++] = locBCond;
                        }
                        else if(locBCond->GetBoundaryConditionType() == SpatialDomains::eRobin)
                        {
                            boost::shared_ptr<SpatialDomains::RobinBoundaryCondition> robinBC  = 
                                boost::static_pointer_cast<SpatialDomains::RobinBoundaryCondition>(locBCond);
                            p_exp->SetValue(-robinBC->m_a.Evaluate(coords[0],coords[1],coords[2])/
                                            robinBC->m_b.Evaluate(coords[0],coords[1],coords[2]));
                            m_bndCondExpansions[cnt] = p_exp;
                            m_bndConditions[cnt++] = locBCond;
                        }
                        else
                        {
                            ASSERTL0(false,"This type of BC not implemented yet");
                        }
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a vertex failed");
                    }
                }
            }
        }
#endif

        DisContField1D::~DisContField1D()
        {
        }

        
        GlobalLinSysSharedPtr DisContField1D::GetGlobalBndLinSys(const GlobalLinSysKey &mkey) 
        {
            GlobalLinSysSharedPtr glo_matrix;
            GlobalLinSysMap::iterator matrixIter = m_globalBndMat->find(mkey);

            if(matrixIter == m_globalBndMat->end())
            {
                glo_matrix = GenGlobalBndLinSys(mkey,*m_traceMap);
                (*m_globalBndMat)[mkey] = glo_matrix;
            }
            else
            {
                glo_matrix = matrixIter->second;
            }

            return glo_matrix;
        }

        void DisContField1D::HelmSolve(const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,       NekDouble> &outarray,
                                       NekDouble lambda)
        {
            int i,j,n,cnt,cnt1,nbndry;
            int nexp = GetExpSize();
            static DNekScalBlkMatSharedPtr InvHDGHelm;
            Array<OneD, unsigned int> vmap;
            Array<OneD,NekDouble> f(m_ncoeffs);
            Array<OneD,NekDouble> e_f, e_bndsol;
            NekDouble tau = 10;

            //----------------------------------
            // Setup RHS Inner product
            //----------------------------------
            IProductWRTBase(inarray,f);
            Vmath::Neg(m_ncoeffs,f,1);
            
            //----------------------------------
            // Solve continuous Boundary System
            //----------------------------------
            int GloBndDofs = m_traceMap->GetNumGlobalBndCoeffs();
            int NumDirBCs  = m_traceMap->GetNumLocalDirBndCoeffs();
            int e_ncoeffs, loc, id;
            NekDouble bval;

            DNekVec Floc;

            // Zero trace space
            Vmath::Zero(GloBndDofs,m_trace,1);

            // Set up boundary space forcing
            for(cnt = cnt1 = n = 0; n < nexp; ++n)
            {
                nbndry = (*m_exp)[n]->NumDGBndryCoeffs(); 		    

                LocalRegions::MatrixKey Umatkey(StdRegions::eHybridDGLamToU, 
                   (*m_exp)[n]->DetExpansionType(),*((*m_exp)[n]), lambda, tau);

                e_ncoeffs = (*m_exp)[n]->GetNcoeffs();
                e_f       = f+cnt;

                DNekVec      ElmtFce(e_ncoeffs, e_f ,eWrapper);
                DNekScalMat &LamToU = *((*m_exp)[n]->GetLocMatrix(Umatkey)); 

                Floc = Transpose(LamToU)*ElmtFce;

                for(i = 0; i < nbndry; ++i)
                {
                    id  = m_traceMap->GetLocalToGlobalBndMap(i+cnt1);

                    m_trace[id] += Floc[i];
                }

                cnt  += e_ncoeffs;
                cnt1 += nbndry;
            }
            
            
            // Copy Dirichlet boundary conditions into trace space        
            for(i = 0; i < m_bndConditions.num_elements(); ++i)
            {
                if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    id          = m_traceMap->GetBndCondCoeffsToGlobalCoeffsMap(i);
                    m_trace[id] = m_bndCondExpansions[i]->GetValue();
                }
                // Add boundary flux
                else if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                {
                    id           = m_traceMap->GetBndCondCoeffsToGlobalCoeffsMap(i);
                    m_trace[id] += m_bndCondExpansions[i]->GetValue();
                }
            }

            
            // Dirichlet boundary forcing 
            for(cnt = n = 0; n < nexp; ++n)
            {
                nbndry = (*m_exp)[n]->NumDGBndryCoeffs();
                // check to see if element has Dirichlet boundary
                // Probably could use a quicker check here
                if(Vmath::Vmin(nbndry,&(m_traceMap->GetLocalToGlobalBndMap())[cnt],1) < NumDirBCs)
                {
                    // Get BndSys Matrix
                    
                    LocalRegions::MatrixKey Bmatkey(StdRegions::eHybridDGHelmBndLam,  (*m_exp)[n]->DetExpansionType(),*((*m_exp)[n]), lambda, tau);
                    
                    DNekScalMat &BndSys = *((*m_exp)[n]->GetLocMatrix(Bmatkey));               
                    DNekVec vin (nbndry);
                    DNekVec vout(nbndry);

                    // Set up Dirichlet Values

                    for(i = 0 ; i < nbndry; ++i)
                    {
                        id = m_traceMap->GetLocalToGlobalBndMap(cnt+i);
                        if(id < NumDirBCs)
                        {
                            vin[i] = m_trace[id]; 
                        } 
                        else
                        {
                            vin[i] = 0.0;
                        }
                    }

                    vout = BndSys*vin;
                    
                    // Subtract vout from forcing terms 
                    for(i = 0; i < nbndry; ++i)
                    {
                        id = m_traceMap->GetLocalToGlobalBndMap(cnt+i);

                        if(id >= NumDirBCs)
                        {
                            m_trace[id] -= vout[i];
                        }
                    }
                }
                cnt += nbndry;
            }


            if(GloBndDofs - NumDirBCs > 0)
            {
                GlobalLinSysKey key(StdRegions::eHybridDGHelmBndLam,
                                    m_traceMap, lambda,tau,eDirectFullMatrix);
                GlobalLinSysSharedPtr LinSys = GetGlobalBndLinSys(key);
                
                Array<OneD,NekDouble> sln = m_trace+NumDirBCs;
                LinSys->Solve(sln,sln);
            }

            //----------------------------------
            // Setup forcing for local interior solves
            //----------------------------------
            Vmath::Zero(m_ncoeffs,&outarray[0],1);
            for(cnt = cnt1 = i = 0; i < nexp; ++i)
            {
                e_bndsol = outarray + cnt;
                e_f      = f + cnt;

                (*m_exp)[i]->GetBoundaryMap(vmap);
                // put boundary solutions into local space; 
                for(j = 0; j < nbndry; ++j)
                {
                    e_bndsol[vmap[j]] = m_trace[m_traceMap->GetLocalToGlobalBndMap(cnt1+j)];
                }

                (*m_exp)[i]->AddHDGHelmholtzTraceTerms(tau, e_bndsol,e_f);
                cnt  += (*m_exp)[i]->GetNcoeffs();
                cnt1 += nbndry;
            }
            
            // Inverse block diagonal interior solve
            if(!InvHDGHelm.get())
            {
                InvHDGHelm = SetupBlockMatrix(StdRegions::eInvHybridDGHelmholtz, lambda, tau);
            }
            DNekVec in (m_ncoeffs,f,eWrapper);
            DNekVec out(m_ncoeffs,outarray,eWrapper);            
            out = (*InvHDGHelm)*in;            
        }
        
    } // end of namespace
} //end of namespace
