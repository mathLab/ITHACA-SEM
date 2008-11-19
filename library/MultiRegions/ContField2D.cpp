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

            m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_ncoeffs,*m_exp,
                                                                                     m_bndCondExpansions,
                                                                                     m_bndConditions,
                                                                                     periodicVertices,
                                                                                     periodicEdges);
	    
	    m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
	    m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }

        ContField2D::ContField2D(SpatialDomains::MeshGraph2D &graph2D,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const std::string variable):
            ContExpList2D(graph2D,false),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph2D,bcs,variable);
            EvaluateBoundaryConditions();

            map<int,int> periodicEdges;
            map<int,int> periodicVertices;
            GetPeriodicEdges(graph2D,bcs,variable,periodicVertices,periodicEdges);

            m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_ncoeffs,*m_exp,
                                                                                 m_bndCondExpansions,
                                                                                 m_bndConditions,
                                                                                 periodicVertices,
                                                                                 periodicEdges);
	    
	    m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
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

            m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_ncoeffs,*m_exp,
                                                                                 m_bndCondExpansions,
                                                                                 m_bndConditions,
                                                                                 periodicVertices,
                                                                                 periodicEdges);
	    
	    m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
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

            m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_ncoeffs,*m_exp,
                                                                                 m_bndCondExpansions,
                                                                                 m_bndConditions,
                                                                                 periodicVertices,
                                                                                 periodicEdges);
	    
	    m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
	    m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }
        
        ContField2D::~ContField2D()
        {
        }
     
        void ContField2D::GenerateBoundaryConditionExpansion(SpatialDomains::MeshGraph2D &graph2D,
                                                             SpatialDomains::BoundaryConditions &bcs, 
                                                             const std::string variable)
        {
            int cnt  = 0;
            SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();   

            int nbnd = bregions.size();

            // count the number of non-periodic boundary regions
            for(int i = 0; i < nbnd; ++i)
            {
                if( ((*(bconditions[i]))[variable])->GetBoundaryConditionType() != SpatialDomains::ePeriodic )
                {
                    cnt++;
                }
            }
            m_bndCondExpansions  = Array<OneD,MultiRegions::ExpList1DSharedPtr>(cnt);
            m_bndConditions      = Array<OneD,SpatialDomains::BoundaryConditionShPtr>(cnt);
            
            SetBoundaryConditionExpansion(graph2D,bcs,variable,m_bndCondExpansions,m_bndConditions);
        }
        
        void ContField2D::FwdTrans(const ExpList &In)
        {
            Array<OneD,NekDouble> rhs(m_contNcoeffs);
            
            // Inner product of forcing 
            ExpList::IProductWRTBase(In.GetPhys(), rhs);

            // Solve the system
            GlobalLinSysKey key(StdRegions::eMass);
            GlobalSolve(key,rhs,m_contCoeffs);
            
            m_transState = eContinuous;
            m_physState = false;
        }

        void ContField2D::MultiplyByInvMassMatrix(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, bool GlobalArrays, bool ZeroBCs)
        {
            GlobalLinSysKey key(StdRegions::eMass);
            int NumDirBcs = m_locToGloMap->GetNumDirichletBndCoeffs();
            
            if(GlobalArrays)
            {
                if(inarray.data() == outarray.data())
                {
                    Array<OneD,NekDouble> tmp(inarray.num_elements());                    
                    Vmath::Vcopy(inarray.num_elements(),inarray,1,tmp,1);

                    Array<OneD,NekDouble> rhs_H = tmp+NumDirBcs;
                    Array<OneD,NekDouble> sol_H = outarray+NumDirBcs;
                    GlobalLinSysSharedPtr invMass_HH = GetGlobalLinSys(key);
                
                    invMass_HH->Solve(rhs_H,sol_H,*m_locToGloMap);            
                }
                else
                {
                    Array<OneD,NekDouble> rhs_H = inarray+NumDirBcs;
                    Array<OneD,NekDouble> sol_H = outarray+NumDirBcs;
                    GlobalLinSysSharedPtr invMass_HH = GetGlobalLinSys(key);
                
                    invMass_HH->Solve(rhs_H,sol_H,*m_locToGloMap);                    
                }

                if(ZeroBCs)
                {
                    Vmath::Zero(NumDirBcs,outarray,1);
                }

            }
            else
            {
                Array<OneD, NekDouble> globaltmp(m_contNcoeffs);

                if(inarray.data() == outarray.data())
                {
                    Array<OneD,NekDouble> tmp(inarray.num_elements());                    
                    Vmath::Vcopy(inarray.num_elements(),inarray,1,tmp,1);
                    Assemble(tmp,outarray);
                }
                else
                {
                    Assemble(inarray,outarray);
                }
                    
                Array<OneD,NekDouble> rhs_H = outarray+NumDirBcs;
                Array<OneD,NekDouble> sol_H = globaltmp+NumDirBcs;
                GlobalLinSysSharedPtr invMass_HH = GetGlobalLinSys(key);
                
                invMass_HH->Solve(rhs_H,sol_H,*m_locToGloMap);

                if(ZeroBCs)
                {
                    Vmath::Zero(NumDirBcs,globaltmp,1);
                }

                GlobalToLocal(globaltmp,outarray);
            }
        }

        // Solve the helmholtz problem assuming that m_contCoeff vector 
        // contains an intial estimate for solution
        void ContField2D::HelmSolve(const ExpList &In, NekDouble lambda)
        {
            GlobalLinSysKey key(StdRegions::eHelmholtz,lambda);
            Array<OneD,NekDouble> rhs(m_contNcoeffs);
        
            
            // Inner product of forcing 
            ExpList::IProductWRTBase(In.GetPhys(), m_coeffs);
            Assemble(m_coeffs,rhs);

            // Note -1.0 term necessary to invert forcing function to
            // be consistent with matrix definition
            Vmath::Neg(m_contNcoeffs,  rhs, 1);
            
            GlobalSolve(key,rhs,m_contCoeffs);

            m_transState = eContinuous;
            m_physState = false;
        }

        void  ContField2D::LaplaceSolve(const ExpList &In,
                                        const Array<OneD, Array<OneD,NekDouble> >& variablecoeffs,
                                        NekDouble time)
        {
            GlobalLinSysKey key(StdRegions::eLaplacian,time,variablecoeffs);
            Array<OneD,NekDouble> rhs(m_contNcoeffs);        
            
            // Inner product of forcing 
            ExpList::IProductWRTBase(In.GetPhys(), m_coeffs);
            Assemble(m_coeffs,rhs);

            // Note -1.0 term necessary to invert forcing function to
            // be consistent with matrix definition
            Vmath::Neg(m_contNcoeffs,  rhs, 1);
            
            GlobalSolve(key,rhs,m_contCoeffs);

            m_transState = eContinuous;
            m_physState = false;
        }

        // Note inout contains initial guess and final output. 
        void ContField2D::GlobalSolve(const GlobalLinSysKey &key, const Array<OneD, const  NekDouble> &rhs, Array<OneD, NekDouble> &inout)
        {
            int i,j;
            int bndcnt=0;
            int NumDirBcs = m_locToGloMap->GetNumDirichletBndCoeffs();
            Array<OneD,NekDouble> sln;
            Array<OneD,NekDouble> fce(m_contNcoeffs,0.0);
            
            //assume m_contCoeffs contains initial solution
            // Set BCs in inout
            for(i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                if(m_bndConditions[i]->GetBoundaryConditionType() != SpatialDomains::eDirichlet)
                {
                    break;
                }
                for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                {
                    inout[m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsMap(bndcnt++)] = (m_bndCondExpansions[i]->GetCoeffs())[j];
                }
            }

            GeneralMatrixOp(key, inout, fce);            

            // Forcing function with Dirichlet conditions 
            Vmath::Vsub(m_contNcoeffs, rhs, 1, fce, 1, fce, 1);

            // Forcing function with weak boundary conditions 
            for(i; i < m_bndCondExpansions.num_elements(); ++i)
            {
                for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                {
                    fce[m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsMap(bndcnt++)] +=  
                        (m_bndCondExpansions[i]->GetCoeffs())[j];
                }
            }

            if(m_contNcoeffs - NumDirBcs > 0)
            {

                GlobalLinSysSharedPtr LinSys = GetGlobalLinSys(key);
                       
                sln = fce+NumDirBcs;
                LinSys->Solve(sln,sln,*m_locToGloMap);
            }

            // Recover solution by adding
            Vmath::Zero(NumDirBcs, fce, 1);
            Vmath::Vadd(m_contNcoeffs, fce, 1, inout, 1, inout, 1);
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
