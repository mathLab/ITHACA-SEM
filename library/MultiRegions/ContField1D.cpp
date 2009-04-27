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

            m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_ncoeffs,*m_exp,
                                                                                     m_bndCondExpansions,
                                                                                     m_bndConditions,
                                                                                     periodicVertices);
	    
	    m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
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

            m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_ncoeffs,*m_exp,
                                                                                     m_bndCondExpansions,
                                                                                     m_bndConditions,
                                                                                     periodicVertices);
	    
	    m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
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

            m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_ncoeffs,*m_exp,
                                                                                     m_bndCondExpansions,
                                                                                     m_bndConditions,
                                                                                     periodicVertices);
	    
	    m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
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
            
            m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_ncoeffs,*m_exp,
                                                                                     m_bndCondExpansions,
                                                                                     m_bndConditions,
                                                                                     periodicVertices);
	    
	    m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
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
        
        void ContField1D::FwdTrans(const Array<OneD, const NekDouble> &inarray,
                                         Array<OneD,       NekDouble> &outarray,
                                   bool  UseContCoeffs)
        {            
            // Inner product of forcing
            Array<OneD,NekDouble> wsp(m_contNcoeffs);  
            IProductWRTBase(inarray,wsp,true);
            
            // Solve the system
            GlobalLinSysKey key(StdRegions::eMass,m_locToGloMap);

            if(UseContCoeffs)
            {
                GlobalSolve(key,wsp,outarray);
            }
            else
            {
                Array<OneD,NekDouble> tmp(m_contNcoeffs,0.0);
                GlobalSolve(key,wsp,tmp);
                GlobalToLocal(tmp,outarray);
            }
        }

        void ContField1D::MultiplyByInvMassMatrix(const Array<OneD, const NekDouble> &inarray,
                                                        Array<OneD,       NekDouble> &outarray,
                                                  bool  UseContCoeffs)
                                                  
        {
            GlobalLinSysKey key(StdRegions::eMass,m_locToGloMap);
            
            if(UseContCoeffs)
            {
                if(inarray.data() == outarray.data())
                {
                    Array<OneD, NekDouble> tmp(m_contNcoeffs,0.0);   
                    Vmath::Vcopy(m_contNcoeffs,inarray,1,tmp,1);
                    GlobalSolve(key,tmp,outarray);
                }
                else
                {
                    GlobalSolve(key,inarray,outarray);
                }
            }
            else
            {
                Array<OneD, NekDouble> globaltmp(m_contNcoeffs,0.0);

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
                    
                GlobalSolve(key,outarray,globaltmp);
                GlobalToLocal(globaltmp,outarray);
            }
        }

        // Solve the helmholtz problem assuming that m_contCoeff vector 
        // contains an intial estimate for solution
        void ContField1D::HelmSolve(const Array<OneD, const NekDouble> &inarray,
                                          Array<OneD,       NekDouble> &outarray,
                                    NekDouble lambda,
                                    bool UseContCoeffs,
                                    Array<OneD, NekDouble>& dirForcing)
        {
            // Inner product of forcing
            Array<OneD,NekDouble> wsp(m_contNcoeffs);  
            IProductWRTBase(inarray,wsp,true);       
            // Note -1.0 term necessary to invert forcing function to
            // be consistent with matrix definition
            Vmath::Neg(m_contNcoeffs, wsp, 1);

            // Solve the system
            GlobalLinSysKey key(StdRegions::eHelmholtz,m_locToGloMap,lambda);

            if(UseContCoeffs)
            {
                GlobalSolve(key,wsp,outarray,dirForcing);
            }
            else
            {
                Array<OneD,NekDouble> tmp(m_contNcoeffs,0.0);
                GlobalSolve(key,wsp,tmp,dirForcing);
                GlobalToLocal(tmp,outarray);
            }
        }

        void ContField1D::GenerateDirBndCondForcing(const GlobalLinSysKey &key, 
                                                    Array<OneD, NekDouble> &inout, 
                                                    Array<OneD, NekDouble> &outarray)
        {       
            int NumDirBcs = m_locToGloMap->GetNumGlobalDirBndCoeffs();

            for(int i = 0; i < NumDirBcs; ++i)
            {
                inout[i] = m_bndCondExpansions[i]->GetValue();
            }

            GeneralMatrixOp(*(key.GetGlobalMatrixKey()),inout,outarray,true);  
        }

        // Note inout contains initial guess and final output. 
        void ContField1D::GlobalSolve(const GlobalLinSysKey &key, 
                                      const Array<OneD, const  NekDouble> &rhs, 
                                      Array<OneD, NekDouble> &inout,
                                      Array<OneD, NekDouble>& dirForcing)
        {
            int i,j;
            int NumDirBcs = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            Array<OneD,NekDouble> fce(m_contNcoeffs);
                     
            if(dirForcing.num_elements())
            {
                // Use the precomputed dirichlet forcing
                Vmath::Vsub(m_contNcoeffs, rhs, 1, dirForcing, 1, fce, 1);
            }
            else
            {
                // Compute the dirichlet forcing
                GenerateDirBndCondForcing(key,inout,fce);
                Vmath::Vsub(m_contNcoeffs, rhs, 1, fce, 1, fce, 1);
            }

            // Forcing function with weak boundary conditions 
            for(i = 0; i < m_bndCondExpansions.num_elements()-NumDirBcs; ++i)
            {
                fce[m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsMap(i+NumDirBcs)] +=  
                    m_bndCondExpansions[i+NumDirBcs]->GetValue();
            }

            if(m_contNcoeffs - NumDirBcs > 0)
            {
                GlobalLinSysSharedPtr LinSys = GetGlobalLinSys(key);
                                       
                Array<OneD,NekDouble> sln(fce+NumDirBcs);
                LinSys->Solve(sln,sln,*m_locToGloMap);
            }

            // Recover solution by adding
            Vmath::Zero(NumDirBcs, fce, 1);
            Vmath::Vadd(m_contNcoeffs, fce, 1, inout, 1, inout, 1);
        }

        GlobalLinSysSharedPtr ContField1D::GetGlobalLinSys(const GlobalLinSysKey &mkey) 
        {
            GlobalLinSysSharedPtr glo_matrix;
            GlobalLinSysMap::iterator matrixIter = m_globalLinSys->find(mkey);

            if(matrixIter == m_globalLinSys->end())
            {
                glo_matrix = GenGlobalLinSys(mkey,m_locToGloMap);
                (*m_globalLinSys)[mkey] = glo_matrix;
            }
            else
            {
                glo_matrix = matrixIter->second;
            }

            return glo_matrix;
        }

    } // end of namespace
} //end of namespace
