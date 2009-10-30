///////////////////////////////////////////////////////////////////////////////
//
// File ContField3D.cpp
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
// Description: Field definition for 3D domain with boundary conditions
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ContField3D.h>

namespace Nektar
{
  namespace MultiRegions
  {

        ContField3D::ContField3D(void):
            ContExpList3D(),
            m_bndCondExpansions(),
            m_bndConditions()
        {
        }

        ContField3D::ContField3D(const ContField3D &In):
            ContExpList3D(In),
            m_bndCondExpansions(In.m_bndCondExpansions),
            m_bndConditions(In.m_bndConditions)
        {
        }

        ContField3D::ContField3D(SpatialDomains::MeshGraph3D &graph3D,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const int bc_loc,
                                 const GlobalSysSolnType solnType):
            ContExpList3D(graph3D,solnType,false),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph3D,bcs,bcs.GetVariable(bc_loc));
            EvaluateBoundaryConditions();

            map<int,int> periodicFaces;
            map<int,int> periodicEdges;
            map<int,int> periodicVertices;
   
            GetPeriodicFaces(graph3D,bcs,bcs.GetVariable(bc_loc),periodicVertices,periodicEdges,periodicFaces);
            
            m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_ncoeffs,*m_exp,solnType,
                                                                                     m_bndCondExpansions,
                                                                                     m_bndConditions,
                                                                                     periodicVertices,
                                                                                     periodicEdges,
                                                                                     periodicFaces); 
        
            m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
            m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }

        ContField3D::ContField3D(SpatialDomains::MeshGraph3D &graph3D,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const std::string variable,
                                 const GlobalSysSolnType solnType):
            ContExpList3D(graph3D,solnType,false),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph3D, bcs, variable);
            EvaluateBoundaryConditions();

            map<int,int> periodicFaces;
            map<int,int> periodicEdges;
            map<int,int> periodicVertices;
            GetPeriodicFaces(graph3D,bcs,variable,periodicVertices,periodicEdges,periodicFaces);

           m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_ncoeffs,*m_exp,solnType,
                                                                                    m_bndCondExpansions,
                                                                                    m_bndConditions,
                                                                                    periodicVertices,
                                                                                    periodicEdges,
                                                                                    periodicFaces); 
            m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
            m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }

        ContField3D::ContField3D(const LibUtilities::BasisKey &Ba,
                                 const LibUtilities::BasisKey &Bb,
                                 const LibUtilities::BasisKey &Bc,
                                 SpatialDomains::MeshGraph3D &graph3D,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const int bc_loc,
                                 const LibUtilities::PointsType TetNb,
                                 const GlobalSysSolnType solnType):
            ContExpList3D(Ba,Bb,Bc,graph3D,TetNb,solnType,false),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph3D,bcs,bcs.GetVariable(bc_loc));
            EvaluateBoundaryConditions();

            map<int,int> periodicFaces;
            map<int,int> periodicEdges;
            map<int,int> periodicVertices;
            GetPeriodicFaces(graph3D,bcs,bcs.GetVariable(bc_loc),periodicVertices,periodicEdges,periodicFaces);

            m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_ncoeffs,*m_exp,solnType,
                                                                                     m_bndCondExpansions,
                                                                                     m_bndConditions,
                                                                                     periodicVertices,
                                                                                     periodicEdges,
                                                                                     periodicFaces);
        
        m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
        m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }

        ContField3D::ContField3D(const LibUtilities::BasisKey &Ba,
                                 const LibUtilities::BasisKey &Bb,
                                 const LibUtilities::BasisKey &Bc,
                                 SpatialDomains::MeshGraph3D &graph3D,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const std::string variable,
                                 const LibUtilities::PointsType TetNb,
                                 const GlobalSysSolnType solnType):
            ContExpList3D(Ba,Bb,Bc,graph3D,TetNb,solnType,false),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph3D,bcs,variable);
            EvaluateBoundaryConditions();

            map<int,int> periodicFaces;
            map<int,int> periodicEdges;
            map<int,int> periodicVertices;
            GetPeriodicFaces(graph3D,bcs,variable,periodicVertices,periodicEdges,periodicFaces);

            m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_ncoeffs,*m_exp,solnType,
                                                                                     m_bndCondExpansions,
                                                                                     m_bndConditions,
                                                                                     periodicVertices,
                                                                                     periodicEdges,
                                                                                     periodicFaces);       
            m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
            m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }

        ContField3D::~ContField3D()
        {
        }

        void ContField3D::GenerateBoundaryConditionExpansion(SpatialDomains::MeshGraph3D &graph3D,
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

        void ContField3D::FwdTrans(const Array<OneD, const NekDouble> &inarray,
                                         Array<OneD,       NekDouble> &outarray,
                                   bool  UseContCoeffs)
        {            
            // Inner product of forcing
            Array<OneD,NekDouble> wsp(m_contNcoeffs);  
            IProductWRTBase(inarray,wsp,true);
            
            // Solve the system
            GlobalLinSysKey key(StdRegions::eMass,
                                m_locToGloMap,
                                m_locToGloMap->GetGlobalSysSolnType());

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

        void ContField3D::MultiplyByInvMassMatrix(const Array<OneD, const NekDouble> &inarray,
                                                        Array<OneD,       NekDouble> &outarray,
                                                  bool  UseContCoeffs)
                                                  
        {
            GlobalLinSysKey key(StdRegions::eMass,
                                m_locToGloMap,
                                m_locToGloMap->GetGlobalSysSolnType());
            
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
        void ContField3D::HelmSolve(const Array<OneD, const NekDouble> &inarray,
                                          Array<OneD,       NekDouble> &outarray,
                                    NekDouble lambda,
                                    bool      UseContCoeffs,
                                    const Array<OneD, const NekDouble>& dirForcing)
        {
            // Inner product of forcing
            Array<OneD,NekDouble> wsp(m_contNcoeffs);  
            IProductWRTBase(inarray,wsp,true);       
            // Note -1.0 term necessary to invert forcing function to
            // be consistent with matrix definition
            Vmath::Neg(m_contNcoeffs, wsp, 1);

            // Forcing function with weak boundary conditions 
            int i,j;
            int bndcnt=m_locToGloMap->GetNumLocalDirBndCoeffs();
            NekDouble sign;
            for(i = m_numDirBndCondExpansions; i < m_bndCondExpansions.num_elements(); ++i)
            {
                for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                {
                    sign = m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsSign(bndcnt);
                    wsp[m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsMap(bndcnt++)] +=  
                        sign * (m_bndCondExpansions[i]->GetCoeffs())[j];
                }
            }

            // Solve the system
            GlobalLinSysKey key(StdRegions::eHelmholtz,
                                m_locToGloMap,lambda,
                                m_locToGloMap->GetGlobalSysSolnType());

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

        void ContField3D::GenerateDirBndCondForcing(const GlobalLinSysKey &key, 
                                                    Array<OneD, NekDouble> &inout, 
                                                    Array<OneD, NekDouble> &outarray)
        {       
            int bndcnt=0;
            const Array<OneD,const int>& map  = m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsMap();
            NekDouble sign;

            for(int i = 0; i < m_numDirBndCondExpansions; ++i)
            {
                const Array<OneD,const NekDouble>& coeffs = m_bndCondExpansions[i]->GetCoeffs();
                for(int j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); ++j)
                {
                    sign = m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsSign(bndcnt);
                    inout[map[bndcnt++]] = sign * coeffs[j];
                }
            }            
            GeneralMatrixOp(*(key.GetGlobalMatrixKey()),inout,outarray,true);   
        }


        // Note inout contains initial guess and final output. 
        void ContField3D::GlobalSolve(const GlobalLinSysKey &key, 
                                      const Array<OneD, const NekDouble>& rhs, 
                                            Array<OneD,       NekDouble>& inout,
                                      const Array<OneD, const NekDouble>& dirForcing)
        {
            int i,j;
            int bndcnt=0;
            int NumDirBcs = m_locToGloMap->GetNumGlobalDirBndCoeffs();

            // STEP 1: SET THE DIRICHLET DOFS TO THE RIGHT VALUE
            //         IN THE SOLUTION ARRAY       
            const Array<OneD,const int>& map  = m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsMap();
            NekDouble sign;
            for(i = 0; i < m_numDirBndCondExpansions; ++i)
            {
                const Array<OneD,const NekDouble>& coeffs = m_bndCondExpansions[i]->GetCoeffs();
                for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); ++j)
                {
                    sign = m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsSign(bndcnt);
                    inout[map[bndcnt++]] = sign * coeffs[j];
                }
            }    

            // STEP 2: CALCULATE THE HOMOGENEOUS COEFFICIENTS
            if(m_contNcoeffs - NumDirBcs > 0)
            {
                GlobalLinSysSharedPtr LinSys = GetGlobalLinSys(key);
                LinSys->Solve(rhs,inout,m_locToGloMap,dirForcing);
            }
        }

        GlobalLinSysSharedPtr ContField3D::GetGlobalLinSys(const GlobalLinSysKey &mkey)
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
    

  } //end of namespace
} //end of namespace
