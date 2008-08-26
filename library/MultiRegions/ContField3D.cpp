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
                                 const int bc_loc):
            ContExpList3D(graph3D,false),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph3D,bcs,bcs.GetVariable(bc_loc));
            EvaluateBoundaryConditions();

            map<int,int> periodicFaces;
            map<int,int> periodicEdges;
            map<int,int> periodicVertices;
   
            GetPeriodicFaces(graph3D,bcs,bcs.GetVariable(bc_loc),periodicVertices,periodicEdges,periodicFaces);

             //TODO:implement
//             m_locToGloMap = MemoryManager<LocalToGlobalMap3D>::AllocateSharedPtr(m_ncoeffs,*m_exp,                                              m_bndCondExpansions, m_bndConditions, periodicVertices, periodicEdges);
        
            m_contNcoeffs = m_locToGloMap->GetTotGloDofs();
            m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }

        ContField3D::ContField3D(SpatialDomains::MeshGraph3D &graph3D,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const std::string variable):
            ContExpList3D(graph3D,false),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph3D, bcs, variable);
            EvaluateBoundaryConditions();

            map<int,int> periodicFaces;
            map<int,int> periodicEdges;
            map<int,int> periodicVertices;
            GetPeriodicFaces(graph3D,bcs,variable,periodicVertices,periodicEdges,periodicFaces);

            // TODO: implement
/*            m_locToGloMap = MemoryManager<LocalToGlobalMap3D>::AllocateSharedPtr(m_ncoeffs,*m_exp,
                                                                                 m_bndCondExpansions,
                                                                                 m_bndConditions,
                                                                                 periodicVertices,
                                                                                 periodicEdges);   */     
            m_contNcoeffs = m_locToGloMap->GetTotGloDofs();
            m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }

        ContField3D::ContField3D(const LibUtilities::BasisKey &Ba,
                                 const LibUtilities::BasisKey &Bb,
                                 const LibUtilities::BasisKey &Bc,
                                 SpatialDomains::MeshGraph3D &graph3D,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const int bc_loc,
                                 const LibUtilities::PointsType TetNb):
            ContExpList3D(Ba,Bb,Bc,graph3D,TetNb,false),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph3D,bcs,bcs.GetVariable(bc_loc));
            EvaluateBoundaryConditions();

            map<int,int> periodicFaces;
            map<int,int> periodicEdges;
            map<int,int> periodicVertices;
            GetPeriodicFaces(graph3D,bcs,bcs.GetVariable(bc_loc),periodicVertices,periodicEdges,periodicFaces);

//             m_locToGloMap = MemoryManager<LocalToGlobalMap3D>::AllocateSharedPtr(m_ncoeffs,*m_exp,
//                                                                                  m_bndCondExpansions,
//                                                                                  m_bndConditions,
//                                                                                  periodicVertices,
//                                                                                  periodicEdges);
        
        m_contNcoeffs = m_locToGloMap->GetTotGloDofs();
        m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }

        ContField3D::ContField3D(const LibUtilities::BasisKey &Ba,
                                 const LibUtilities::BasisKey &Bb,
                                 const LibUtilities::BasisKey &Bc,
                                 SpatialDomains::MeshGraph3D &graph3D,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const std::string variable,
                                 const LibUtilities::PointsType TetNb):
            ContExpList3D(Ba,Bb,Bc,graph3D,TetNb,false),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph3D,bcs,variable);
            EvaluateBoundaryConditions();

            map<int,int> periodicFaces;
            map<int,int> periodicEdges;
            map<int,int> periodicVertices;
            GetPeriodicFaces(graph3D,bcs,variable,periodicVertices,periodicEdges,periodicFaces);

            // TODO: implement
//             m_locToGloMap = MemoryManager<LocalToGlobalMap3D>::AllocateSharedPtr(m_ncoeffs,*m_exp,
//                                                                                  m_bndCondExpansions,
//                                                                                  m_bndConditions,
//                                                                                  periodicVertices,
//                                                                                  periodicEdges);       
            m_contNcoeffs = m_locToGloMap->GetTotGloDofs();
            m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }

        ContField3D::~ContField3D()
        {
        }

        void ContField3D::GenerateBoundaryConditionExpansion(SpatialDomains::MeshGraph3D &graph3D,
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
            m_bndCondExpansions  = Array<OneD,MultiRegions::ExpList2DSharedPtr>(cnt);
            m_bndConditions      = Array<OneD,SpatialDomains::BoundaryConditionShPtr>(cnt);
            
            SetBoundaryConditionExpansion(graph3D,bcs,variable,m_bndCondExpansions,m_bndConditions);
        }

        void ContField3D::FwdTrans(const ExpList &In)
        {
            GlobalLinSysKey key(StdRegions::eMass);
            GlobalSolve(key,In);

            m_transState = eContinuous;
            m_physState = false;
        }

       // Solve the helmholtz problem assuming that m_contCoeff vector 
        // contains an intial estimate for solution
        void ContField3D::HelmSolve(const ExpList &In, NekDouble lambda)
        {
            GlobalLinSysKey key(StdRegions::eHelmholtz,lambda);
            // Note -1.0 term necessary to invert forcing function to
            // be consistent with matrix definition
            GlobalSolve(key,In,-1.0);
        }

        void ContField3D::GlobalSolve(const GlobalLinSysKey &key,
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

        GlobalLinSysSharedPtr ContField3D::GetGlobalLinSys(const GlobalLinSysKey &mkey)
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
    

  } //end of namespace
} //end of namespace
