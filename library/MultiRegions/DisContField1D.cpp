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

        DisContField1D::DisContField1D(const DisContField1D &In):
            ExpList1D(In),
            m_bndConstraint(In.m_bndConstraint),
            m_bndTypes(In.m_bndTypes),
            m_globalBndMat(In.m_globalBndMat)
        {
        }

        DisContField1D::DisContField1D(SpatialDomains::MeshGraph1D &graph1D,
            SpatialDomains::BoundaryConditions &bcs, const int bc_loc):
            ExpList1D(graph1D)     
        {
            GenerateField1D(bcs,bcs.GetVariable(bc_loc));
            m_lambdaMap = MemoryManager<LocalToGlobalBndryMap1D>::AllocateSharedPtr(*m_exp,graph1D,m_bndConstraint,m_bndTypes);
        }

        DisContField1D::DisContField1D(SpatialDomains::MeshGraph1D &graph1D,
            SpatialDomains::BoundaryConditions &bcs, 
            const std::string variable):
            ExpList1D(graph1D)  
        {
            GenerateField1D(bcs,variable);
            m_lambdaMap = MemoryManager<LocalToGlobalBndryMap1D>::AllocateSharedPtr(*m_exp,graph1D,m_bndConstraint,m_bndTypes);
        }

        void DisContField1D::GenerateField1D(SpatialDomains::BoundaryConditions &bcs,   const std::string variable)
        {
            int i,nbnd;
            int cnt=0;
            LocalRegions::PointExpSharedPtr  p_exp;
            SpatialDomains::BoundaryConditionShPtr locBCond;
            SpatialDomains::VertexComponentSharedPtr vert;

            SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();

            nbnd = bregions.size();
            m_bndConstraint = Array<OneD,LocalRegions::PointExpSharedPtr>(nbnd);
            m_bndTypes = Array<OneD,SpatialDomains::BoundaryConditionType>(nbnd);

            // Set up matrix map
            m_globalBndMat   = MemoryManager<GlobalLinSysMap>::AllocateSharedPtr();
            // list Dirichlet boundaries first
            for(i = 0; i < nbnd; ++i)
            {
                locBCond = (*(bconditions[i]))[variable];

                if(locBCond->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    if(vert = boost::dynamic_pointer_cast<SpatialDomains::VertexComponent>((*(*bregions[i])[0])[0]))
                    {
                        Array<OneD,NekDouble> coords(3,0.0);

                        p_exp = MemoryManager<LocalRegions::PointExp>::AllocateSharedPtr(vert);
                        vert->GetCoords(coords);
                        p_exp->SetValue(boost::static_pointer_cast<SpatialDomains::DirichletBoundaryCondition>(locBCond)
                                        ->m_DirichletCondition.Evaluate(coords[0],coords[1],coords[2]));

                        m_bndConstraint[cnt] = p_exp;
                        m_bndTypes[cnt++] = SpatialDomains::eDirichlet;
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
                            m_bndConstraint[cnt] = p_exp;
                            m_bndTypes[cnt++] = SpatialDomains::eNeumann;
                        }
                        else if(locBCond->GetBoundaryConditionType() == SpatialDomains::eRobin)
                        {
                            boost::shared_ptr<SpatialDomains::RobinBoundaryCondition> robinBC  = 
                                boost::static_pointer_cast<SpatialDomains::RobinBoundaryCondition>(locBCond);
                            p_exp->SetValue(-robinBC->m_a.Evaluate(coords[0],coords[1],coords[2])/
                                            robinBC->m_b.Evaluate(coords[0],coords[1],coords[2]));
                            m_bndConstraint[cnt] = p_exp;
                            m_bndTypes[cnt++] = SpatialDomains::eRobin;
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

        DisContField1D::~DisContField1D()
        {
        }

        
        GlobalLinSysSharedPtr DisContField1D::GetGlobalBndLinSys(const GlobalLinSysKey &mkey) 
        {
            GlobalLinSysSharedPtr glo_matrix;
            GlobalLinSysMap::iterator matrixIter = m_globalBndMat->find(mkey);

            if(matrixIter == m_globalBndMat->end())
            {
                glo_matrix = GenGlobalBndLinSys(mkey);
                (*m_globalBndMat)[mkey] = glo_matrix;
            }
            else
            {
                glo_matrix = matrixIter->second;
            }

            return glo_matrix;
        }

	GlobalLinSysSharedPtr DisContField1D::GenGlobalBndLinSys(const GlobalLinSysKey &mkey)
	{
            int i,j,n,gid1,gid2,loc_lda,cnt;
            StdRegions::StdExpansionVectorIter def;
            DNekLinSysSharedPtr   linsys;
            GlobalLinSysSharedPtr returnlinsys;

            int totDofs       = (*m_exp).size()+1;
            int NumDirBCs     = m_lambdaMap->GetNumDirichletDofs();
            unsigned int rows = totDofs - NumDirBCs;
            unsigned int cols = totDofs - NumDirBCs;
            NekDouble zero = 0.0; 
            NekDouble helmfactor = mkey.GetFactor1();
            NekDouble tau        = mkey.GetFactor2();

            DNekMat LocBlk;

            DNekMatSharedPtr Gmat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,zero);            
            ASSERTL0(mkey.GetLinSysType() == StdRegions::eUnifiedDGHelmBndSys,
                     "Routine is only set up for UnifiedDGHelmholtz");
            

            // fill global matrix 
            for(n = cnt = 0; n < (*m_exp).size(); ++n)
            {
                // Matrix to Bnd Sys
                LocalRegions::MatrixKey Umatkey(StdRegions::eUnifiedDGHelmBndSys, (*m_exp)[n]->DetShapeType(),*((*m_exp)[n]), helmfactor,tau);
                DNekScalMat &BndSys = *((*m_exp)[n]->GetLocMatrix(Umatkey)); 
                
                loc_lda = BndSys.GetColumns();
                
                for(i = 0; i < loc_lda; ++i)
                {
                    gid1 = m_lambdaMap->GetBndMap(cnt + i);
                    if(gid1 >= NumDirBCs)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2 = m_lambdaMap->GetBndMap(cnt + j);
                            if(gid2 >= NumDirBCs)
                            {
                                (*Gmat)(gid1-NumDirBCs,gid2-NumDirBCs) 
                                    += (BndSys)(i,j);
                            }
                        }		
                    }
                }
                cnt += (*m_exp)[n]->NumBndryCoeffs(); 		    
            }
            
            // Believe that we need a call of the type:
            // linsys=MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmat,eWrapper);
            linsys       = MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmat);
            returnlinsys = MemoryManager<GlobalLinSys>::AllocateSharedPtr(mkey,linsys);
            return returnlinsys;
        }
        
        void DisContField1D::HelmSolve(DisContField1D &Fce, NekDouble lambda)
        {
            int i,j,n,cnt,cnt1,nbndry;
            int nexp = GetExpSize();
            static DNekScalBlkMatSharedPtr InvUDGHelm;
            StdRegions::StdExpMap vmap;
            Array<OneD,NekDouble> f(m_ncoeffs);
            Array<OneD,NekDouble> e_f, e_bndsol;
            NekDouble tau = 10;

            //----------------------------------
            // Setup RHS Inner product
            //----------------------------------
            IProductWRTBase(Fce.GetPhys(),f);
            Vmath::Neg(m_ncoeffs,&f[0],1);
            
            //----------------------------------
            // Solve continuous Boundary System
            //----------------------------------
            int GloBndDofs = m_lambdaMap->GetTotGloBndDofs();
            int NumDirBCs  = m_lambdaMap->GetNumDirichletDofs();
            int e_ncoeffs;
            Array<OneD,NekDouble> BndSol(GloBndDofs,0.0);
            DNekVec Floc;

            // Set up boundary space forcing
            nbndry = (*m_exp)[0]->NumBndryCoeffs(); 		    
            for(cnt = cnt1 = n = 0; n < nexp; ++n)
            {
                LocalRegions::MatrixKey Umatkey(StdRegions::eUnifiedDGLamToU, 
                   (*m_exp)[n]->DetShapeType(),*((*m_exp)[n]), lambda, tau);

                e_ncoeffs = (*m_exp)[n]->GetNcoeffs();
                e_f       = f+cnt;

                DNekVec ElmtFce(e_ncoeffs, e_f ,eWrapper);
                DNekScalMat &LamToU = *((*m_exp)[n]->GetLocMatrix(Umatkey)); 
#if 0
                Floc = LamToU.Transpose()*ElmtFce;
#else
                DNekMat TLamToU = *LamToU.GetOwnedMatrix();
                TLamToU.Transpose();
                Floc = LamToU.Scale()*TLamToU*ElmtFce;
#endif           
                for(i = 0; i < nbndry; ++i)
                {
                    BndSol[m_lambdaMap->GetBndMap(i+cnt1)] += Floc[i];
                }
                cnt  += e_ncoeffs;
                cnt1 += nbndry;
            }
            
            if(GloBndDofs - NumDirBCs > 0)
            {
                GlobalLinSysKey key(StdRegions::eUnifiedDGHelmBndSys,
                                    lambda,tau,eDirectFullMatrix);
                GlobalLinSysSharedPtr LinSys = GetGlobalBndLinSys(key);
                
                Array<OneD,NekDouble> sln = BndSol+NumDirBCs;
                LinSys->Solve(sln,sln,*m_lambdaMap);
            }

            // Set Dirichlet Boundary Conditions into BndSol
            for(i = 0; i < NumDirBCs; ++i)
            {
                if(m_bndTypes[i] == SpatialDomains::eDirichlet)
                {  
                    BndSol[i] = m_bndConstraint[i]->GetValue();
                }
            }

            //----------------------------------
            // Setup forcing for local interior solves
            //----------------------------------
            Vmath::Zero(m_ncoeffs,&m_coeffs[0],1);
            for(cnt = cnt1 = i = 0; i < nexp; ++i)
            {
                e_bndsol = m_coeffs + cnt;
                e_f      = f + cnt;

                (*m_exp)[i]->MapTo(StdRegions::eForwards,vmap);
                // put boundary solutions into local space; 
                for(j = 0; j < nbndry; ++j)
                {
                    e_bndsol[vmap[j]] = BndSol[m_lambdaMap->GetBndMap(cnt1+j)];
                }

                (*m_exp)[i]->UDGHelmholtzBoundaryTerms(tau, e_bndsol,e_f);
                cnt  += (*m_exp)[i]->GetNcoeffs();
                cnt1 += nbndry;
            }
            
            // Inverse block diagonal interior solve
            if(!InvUDGHelm.get())
            {
                InvUDGHelm = SetupBlockMatrix(StdRegions::eInvUnifiedDGHelmholtz, lambda, tau);
            }
            DNekVec in (m_ncoeffs,f,eWrapper);
            DNekVec out(m_ncoeffs,m_coeffs,eWrapper);            
            out = (*InvUDGHelm)*in;            
        }
        
    } // end of namespace
} //end of namespace
