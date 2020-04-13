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
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>

using namespace std;

namespace Nektar
{
  namespace MultiRegions
  {

        ContField3D::ContField3D():
            DisContField3D(),
            m_locToGloMap(),
            m_globalMat(),
            m_globalLinSysManager(
                std::bind(
                    &ContField3D::GenGlobalLinSys, this, std::placeholders::_1),
                std::string("GlobalLinSys"))
        {
        }


        /**
         * Given a mesh \a graph2D, containing information about the domain and
         * the spectral/hp element expansion, this constructor fills the list
         * of local expansions #m_exp with the proper expansions, calculates
         * the total number of quadrature points \f$\boldsymbol{x}_i\f$ and
         * local expansion coefficients \f$\hat{u}^e_n\f$ and allocates memory
         * for the arrays #m_coeffs and #m_phys. Furthermore, it constructs the
         * mapping array (contained in #m_locToGloMap) for the transformation
         * between local elemental level and global level, it calculates the
         * total number global expansion coefficients \f$\hat{u}_n\f$ and
         * allocates memory for the array #m_coeffs. The constructor also
         * discretises the boundary conditions, specified by the argument \a
         * bcs, by expressing them in terms of the coefficient of the expansion
         * on the boundary.
         *
         * @param   pSession    Session information.
         * @param   graph3D     A mesh, containing information about the domain
         *                      and the spectral/hp element expansion.
         * @param   variable    The variable associated with this field.
         */
        ContField3D::ContField3D(const LibUtilities::SessionReaderSharedPtr &pSession,
                                 const SpatialDomains::MeshGraphSharedPtr &graph3D,
                                 const std::string &variable,
                                 const bool CheckIfSingularSystem,
                                 const Collections::ImplementationType ImpType):
                DisContField3D(pSession,graph3D,variable,false,ImpType),
                m_globalMat(MemoryManager<GlobalMatrixMap>::AllocateSharedPtr()),
                m_globalLinSysManager(
                    std::bind(&ContField3D::GenGlobalLinSys, this,  std::placeholders::_1),
                    std::string("GlobalLinSys"))
        {
            m_locToGloMap = MemoryManager<AssemblyMapCG>::AllocateSharedPtr(
                m_session,m_ncoeffs,*this,m_bndCondExpansions,m_bndConditions,
                CheckIfSingularSystem, variable,
                m_periodicVerts, m_periodicEdges, m_periodicFaces);

            if (m_session->DefinesCmdLineArgument("verbose"))
            {
                m_locToGloMap->PrintStats(std::cout, variable);
            }
        }


        /**
         * Given a mesh \a graph3D, containing information about the domain and
         * the spectral/hp element expansion, this constructor fills the list
         * of local expansions #m_exp with the proper expansions, calculates
         * the total number of quadrature points \f$\boldsymbol{x}_i\f$ and
         * local expansion coefficients \f$\hat{u}^e_n\f$ and allocates memory
         * for the arrays #m_coeffs and #m_phys. Furthermore, it constructs the
         * mapping array (contained in #m_locToGloMap) for the transformation
         * between local elemental level and global level, it calculates the
         * total number global expansion coefficients \f$\hat{u}_n\f$ and
         * allocates memory for the array #m_coeffs. The constructor also
         * discretises the boundary conditions, specified by the argument \a
         * bcs, by expressing them in terms of the coefficient of the expansion
         * on the boundary.
         *
         * @param   In          Existing ContField2D object used to provide the
         *                      local to global mapping information and
         *                      global solution type.
         * @param   graph3D     A mesh, containing information about the domain
         *                      and the spectral/hp element expansion.
         * @param   bcs         The boundary conditions.
         * @param   bc_loc
         */
        ContField3D::ContField3D(const ContField3D &In,
                                 const SpatialDomains::MeshGraphSharedPtr &graph3D,
                                 const std::string &variable,
                                 const bool CheckIfSingularSystem):
	    DisContField3D(In,graph3D,variable,false),
            m_globalMat   (MemoryManager<GlobalMatrixMap>::AllocateSharedPtr()),
            m_globalLinSysManager(
                std::bind(&ContField3D::GenGlobalLinSys, this,  std::placeholders::_1),
                std::string("GlobalLinSys"))

        {
            if(!SameTypeOfBoundaryConditions(In) || CheckIfSingularSystem)
            {
                SpatialDomains::BoundaryConditions bcs(m_session, graph3D);
                m_locToGloMap = MemoryManager<AssemblyMapCG>::AllocateSharedPtr(
                    m_session,m_ncoeffs,*this,m_bndCondExpansions,m_bndConditions,
                    CheckIfSingularSystem, variable,
                    m_periodicVerts, m_periodicEdges, m_periodicFaces);

                if (m_session->DefinesCmdLineArgument("verbose"))
                {
                    m_locToGloMap->PrintStats(std::cout, variable);
                }
            }
            else
            {
                m_locToGloMap = In.m_locToGloMap;
            }
        }


        ContField3D::ContField3D(const ContField3D &In):
                DisContField3D(In),
                m_locToGloMap(In.m_locToGloMap),
                m_globalMat(In.m_globalMat),
                m_globalLinSysManager(In.m_globalLinSysManager)
        {
        }


        ContField3D::~ContField3D()
        {
        }


        /**
         * Given the coefficients of an expansion, this function evaluates the
         * spectral/hp expansion \f$u^{\delta}(\boldsymbol{x})\f$ at the
         * quadrature points \f$\boldsymbol{x}_i\f$. This operation is
         * evaluated locally by the function ExpList#BwdTrans.
         *
         * The coefficients of the expansion should be contained in the variable
         * #inarray of the ExpList object \a In. The resulting physical values
         * at the quadrature points \f$u^{\delta}(\boldsymbol{x}_i)\f$ are
         * stored in the array #outarray.
         *
         * @param   In          An ExpList, containing the local coefficients
         *                      \f$\hat{u}_n^e\f$ in its array #m_coeffs.
         */
        void ContField3D::v_BwdTrans(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,       NekDouble> &outarray)
        {
            BwdTrans_IterPerExp(inarray,outarray);
        }


        /**
         * The operation is evaluated locally (i.e. with respect to all local
         * expansion modes) by the function ExpList#IProductWRTBase. The inner
         * product with respect to the global expansion modes is than obtained
         * by a global assembly operation.
         *
         * The values of the function \f$f(\boldsymbol{x})\f$ evaluated at the
         * quadrature points \f$\boldsymbol{x}_i\f$ should be contained in the
         * variable #m_phys of the ExpList object \a in. The result is stored
         * in the array #m_coeffs.
         *
         * @param   In          An ExpList, containing the discrete evaluation
         *                      of \f$f(\boldsymbol{x})\f$ at the quadrature
         *                      points in its array #m_phys.
         */
        void ContField3D::v_IProductWRTBase(const Array<OneD, const NekDouble> &inarray,
                                            Array<OneD, NekDouble> &outarray)
        {
            IProductWRTBase_IterPerExp(inarray,outarray);
        }


      void ContField3D::v_FwdTrans(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,       NekDouble> &outarray)
      {
          // Inner product of forcing
          Array<OneD,NekDouble> wsp(m_ncoeffs);
          IProductWRTBase(inarray,wsp);
          
          // Solve the system
          GlobalLinSysKey key(StdRegions::eMass, m_locToGloMap);
          
          GlobalSolve(key,wsp,outarray);
      }


      void ContField3D::v_MultiplyByInvMassMatrix(
                                          const Array<OneD, const NekDouble> &inarray,
                                          Array<OneD,       NekDouble> &outarray)
          
      {
          GlobalLinSysKey key(StdRegions::eMass,m_locToGloMap);
          GlobalSolve(key,inarray,outarray);
      }


      // Note inout contains initial guess and final output.
      void ContField3D::GlobalSolve(
          const GlobalLinSysKey              &key,
          const Array<OneD, const NekDouble> &locrhs,
                Array<OneD,       NekDouble> &inout,
          const Array<OneD, const NekDouble> &dirForcing)
      {
          int NumDirBcs = m_locToGloMap->GetNumGlobalDirBndCoeffs();
          int contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();

          // STEP 1: SET THE DIRICHLET DOFS TO THE RIGHT VALUE
          //         IN THE SOLUTION ARRAY
          v_ImposeDirichletConditions(inout);

          // STEP 2: CALCULATE THE HOMOGENEOUS COEFFICIENTS
          if(contNcoeffs - NumDirBcs > 0)
          {
                GlobalLinSysSharedPtr LinSys = GetGlobalLinSys(key);
                LinSys->Solve(locrhs,inout,m_locToGloMap,dirForcing);
          }
      }

      GlobalLinSysSharedPtr ContField3D::GetGlobalLinSys(const GlobalLinSysKey &mkey)
      {
          return m_globalLinSysManager[mkey];
      }


      GlobalLinSysSharedPtr ContField3D::GenGlobalLinSys(const GlobalLinSysKey &mkey)
      {
          ASSERTL1(mkey.LocToGloMapIsDefined(),
                   "To use method must have a AssemblyMap "
                   "attached to key");
          return ExpList::GenGlobalLinSys(mkey, m_locToGloMap);
      }


      /**
       * Returns the global matrix associated with the given GlobalMatrixKey.
       * If the global matrix has not yet been constructed on this field,
       * it is first constructed using GenGlobalMatrix().
       * @param   mkey        Global matrix key.
       * @returns Assocated global matrix.
       */
      GlobalMatrixSharedPtr ContField3D::GetGlobalMatrix(const GlobalMatrixKey &mkey)
      {
          ASSERTL1(mkey.LocToGloMapIsDefined(),
                   "To use method must have a AssemblyMap "
                   "attached to key");

            GlobalMatrixSharedPtr glo_matrix;
            auto matrixIter = m_globalMat->find(mkey);

            if(matrixIter == m_globalMat->end())
            {
                glo_matrix = GenGlobalMatrix(mkey,m_locToGloMap);
                (*m_globalMat)[mkey] = glo_matrix;
            }
            else
            {
                glo_matrix = matrixIter->second;
            }

            return glo_matrix;
      }


      void ContField3D::v_ImposeDirichletConditions(Array<OneD,NekDouble>& outarray)
      {
            int i,j;
            int bndcnt=0;

            Array<OneD, NekDouble> sign = m_locToGloMap->
                GetBndCondCoeffsToLocalCoeffsSign();
            const Array<OneD, const int> map= m_locToGloMap->
                GetBndCondCoeffsToLocalCoeffsMap();
            
            for(i = 0; i < m_bndCondExpansions.size(); ++i)
            {
                if(m_bndConditions[i]->GetBoundaryConditionType() ==
                   SpatialDomains::eDirichlet)
                {
                    const Array<OneD, NekDouble> bndcoeff =
                        (m_bndCondExpansions[i])->GetCoeffs(); 

                    if(m_locToGloMap->GetSignChange())
                    {
                        for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                        {
                            outarray[map[bndcnt + j]] = sign[bndcnt + j] * bndcoeff[j]; 
                        }
                    }
                    else
                    {
                        for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                        {
                            outarray[map[bndcnt+j]] = bndcoeff[j]; 
                        }
                    }                    
                }
                
                bndcnt += m_bndCondExpansions[i]->GetNcoeffs();
            }

            // communicate local Dirichlet coeffsthat are just
            // touching a dirichlet boundary on another partition
            set<int> &ParallelDirBndSign = m_locToGloMap->GetParallelDirBndSign();

            for (auto &it : ParallelDirBndSign)
            {
                outarray[it] *= -1;
            }
                
            m_locToGloMap->UniversalAbsMaxBnd(outarray);

            for (auto &it : ParallelDirBndSign)
            {
                outarray[it] *= -1;
            }

            // sort local dirichlet coeffs that are just touching a
            // dirichlet boundary
            set<ExtraDirDof> &copyLocalDirDofs = m_locToGloMap->GetCopyLocalDirDofs();

            for (auto &it : copyLocalDirDofs)
            {
                outarray[std::get<0>(it)] =
                    outarray[std::get<1>(it)]*std::get<2>(it);
            }
      }          
      
      void ContField3D::v_FillBndCondFromField(void)
      {
            int bndcnt = 0;

            Array<OneD, NekDouble> sign = m_locToGloMap->
                GetBndCondCoeffsToLocalCoeffsSign();
            const Array<OneD, const int> bndmap= m_locToGloMap->
                GetBndCondCoeffsToLocalCoeffsMap();
            
            for(int i = 0; i < m_bndCondExpansions.size(); ++i)
            {
                Array<OneD, NekDouble>& coeffs = m_bndCondExpansions[i]->UpdateCoeffs();
                
                if(m_locToGloMap->GetSignChange())
                {
                    for(int j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); ++j)
                    {
                        coeffs[j] = sign[bndcnt+j] * m_coeffs[bndmap[bndcnt+j]];
                    }
                }
                else
                {
                    for(int j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); ++j)
                    {
                        coeffs[j] = m_coeffs[bndmap[bndcnt+j]];
                    }
                }

                bndcnt += m_bndCondExpansions[i]->GetNcoeffs();
            }
      }

      void ContField3D::v_FillBndCondFromField(const int nreg)
      {
          int bndcnt = 0;
          
          ASSERTL1(nreg < m_bndCondExpansions.size(),
                   "nreg is out or range since this many boundary "
                   "regions to not exist");
          
          Array<OneD, NekDouble> sign = m_locToGloMap->
              GetBndCondCoeffsToLocalCoeffsSign();
          const Array<OneD, const int> bndmap= m_locToGloMap->
              GetBndCondCoeffsToLocalCoeffsMap();
          
            // Now fill in all other Dirichlet coefficients.
          Array<OneD, NekDouble>& coeffs =
              m_bndCondExpansions[nreg]->UpdateCoeffs();

          for(int j = 0; j < nreg; ++j)
          {
              if(m_bndConditions[j]->GetBoundaryConditionType()
                 == SpatialDomains::ePeriodic)
              {
                  continue;
              }
              bndcnt += m_bndCondExpansions[j]->GetNcoeffs();
          }
          
          if(m_locToGloMap->GetSignChange())
          {
              for(int j = 0; j < (m_bndCondExpansions[nreg])->GetNcoeffs(); ++j)
              {
                  coeffs[j] = sign[bndcnt + j] * m_coeffs[bndmap[bndcnt + j]];
              }
          }
          else
          {
              for(int j = 0; j < (m_bndCondExpansions[nreg])->GetNcoeffs(); ++j)
              {
                  coeffs[j] = m_coeffs[bndmap[bndcnt + j]];
              }
          }
      }

      void ContField3D::v_LocalToGlobal(bool useComm)
      {
          m_locToGloMap->LocalToGlobal(m_coeffs, m_coeffs,useComm);
      }


      void ContField3D::v_LocalToGlobal(
          const Array<OneD, const NekDouble> &inarray,
          Array<OneD,NekDouble> &outarray,
          bool useComm)
      {
          m_locToGloMap->LocalToGlobal(inarray, outarray, useComm);
      }


      void ContField3D::v_GlobalToLocal(void)
      {
          m_locToGloMap->GlobalToLocal(m_coeffs, m_coeffs);
      }


      void ContField3D::v_GlobalToLocal(
          const Array<OneD, const NekDouble> &inarray,
          Array<OneD,NekDouble> &outarray)
      {
          m_locToGloMap->GlobalToLocal(inarray, outarray);
      }


      void ContField3D::v_HelmSolve(
                                    const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD,       NekDouble> &outarray,
                                    const StdRegions::ConstFactorMap &factors,
                                    const StdRegions::VarCoeffMap &varcoeff,
                                    const MultiRegions::VarFactorsMap &varfactors,
                                    const Array<OneD, const NekDouble> &dirForcing,
                                    const bool PhysSpaceForcing)
      {
          int i,j;

          //----------------------------------
          //  Setup RHS Inner product
          //----------------------------------
          // Inner product of forcing
          Array<OneD,NekDouble> wsp(m_ncoeffs);
          if(PhysSpaceForcing)
          {
              IProductWRTBase(inarray,wsp);
              // Note -1.0 term necessary to invert forcing function to
              // be consistent with matrix definition
              Vmath::Neg(m_ncoeffs, wsp, 1);
          }
          else
          {
              Vmath::Smul(m_ncoeffs,-1.0,inarray,1,wsp,1);
          }
          
          Array<OneD, NekDouble> sign = m_locToGloMap->
              GetBndCondCoeffsToLocalCoeffsSign();
          const Array<OneD, const int> map= m_locToGloMap->
              GetBndCondCoeffsToLocalCoeffsMap();
          int bndcnt = 0; 
          // Add weak boundary conditions to forcing
          for(i = 0; i < m_bndCondExpansions.size(); ++i)
          {
              if(m_bndConditions[i]->GetBoundaryConditionType() ==
                 SpatialDomains::eNeumann ||
                 m_bndConditions[i]->GetBoundaryConditionType() ==
                 SpatialDomains::eRobin)
              {
                  const Array<OneD, NekDouble> bndcoeff =
                      (m_bndCondExpansions[i])->GetCoeffs(); 
                  
                  if(m_locToGloMap->GetSignChange())
                  {
                      for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                      {
                          wsp[map[bndcnt + j]] += sign[bndcnt + j] * bndcoeff[j]; 
                      }
                  }
                  else
                  {
                      for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                      {
                          wsp[map[bndcnt + j]] += bndcoeff[j]; 
                      }
                  }
              }
              bndcnt += m_bndCondExpansions[i]->GetNcoeffs();
          }

          // Solve the system
          GlobalLinSysKey key(StdRegions::eHelmholtz, m_locToGloMap,
                              factors,varcoeff,varfactors);
          
          GlobalSolve(key,wsp,outarray,dirForcing);
      }

      void ContField3D::v_GeneralMatrixOp(
          const GlobalMatrixKey             &gkey,
          const Array<OneD,const NekDouble> &inarray,
          Array<OneD,      NekDouble> &outarray)
      {
          GeneralMatrixOp_IterPerExp(gkey,inarray,outarray);
      }

    /**
    * First compute the inner product of forcing function with respect to
    * base, and then solve the system with the linear advection operator.
    * @param   velocity    Array of advection velocities in physical space
    * @param   inarray     Forcing function.
    * @param   outarray    Result.
    * @param   lambda      reaction coefficient
    * @param   dirForcing  Dirichlet Forcing.
    */
    // could combine this with HelmholtzCG.
    void ContField3D::v_LinearAdvectionDiffusionReactionSolve(
        const Array<OneD, Array<OneD, NekDouble> > &velocity,
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &outarray,
        const NekDouble lambda,
        const Array<OneD, const NekDouble> &dirForcing)
    {
        // Inner product of forcing
        Array<OneD,NekDouble> wsp(m_ncoeffs);
        IProductWRTBase(inarray,wsp);
        
        // Note -1.0 term necessary to invert forcing function to
        // be consistent with matrix definition
        Vmath::Neg(m_ncoeffs, wsp, 1);
        
        // Forcing function with weak boundary conditions
        int i,j;
        int bndcnt=0;
        Array<OneD, NekDouble> sign = m_locToGloMap->
            GetBndCondCoeffsToLocalCoeffsSign();
        const Array<OneD, const int> map= m_locToGloMap->
            GetBndCondCoeffsToLocalCoeffsMap();
        // Add weak boundary conditions to forcing
        for (i = 0; i < m_bndCondExpansions.size(); ++i)
        {
            if(m_bndConditions[i]->GetBoundaryConditionType() ==
               SpatialDomains::eNeumann ||
               m_bndConditions[i]->GetBoundaryConditionType() ==
               SpatialDomains::eRobin)
            {
                const Array<OneD, NekDouble> bndcoeff =
                    (m_bndCondExpansions[i])->GetCoeffs(); 
                
                if(m_locToGloMap->GetSignChange())
                {
                    for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                    {
                        wsp[map[bndcnt + j]] += sign[bndcnt + j] * bndcoeff[j]; 
                    }
                }
                else
                {
                    for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                    {
                        wsp[map[bndcnt+j]] += bndcoeff[bndcnt + j]; 
                    }
                }                    
            }
            
            bndcnt += m_bndCondExpansions[i]->GetNcoeffs();
        }
        
        // Solve the system
        StdRegions::ConstFactorMap factors;
        factors[StdRegions::eFactorLambda] = lambda;
        StdRegions::VarCoeffMap varcoeffs;
        varcoeffs[StdRegions::eVarCoeffVelX] = velocity[0];
        varcoeffs[StdRegions::eVarCoeffVelY] = velocity[1];
        varcoeffs[StdRegions::eVarCoeffVelZ] = velocity[2];
        GlobalLinSysKey key(StdRegions::eLinearAdvectionDiffusionReaction,
                            m_locToGloMap,
                            factors,
                            varcoeffs);

        GlobalSolve(key, wsp, outarray, dirForcing);
    }


      int ContField3D::GetGlobalMatrixNnz(const GlobalMatrixKey &gkey)
      {
          ASSERTL1(gkey.LocToGloMapIsDefined(),
                   "To use method must have a AssemblyMap "
                   "attached to key");

          auto matrixIter = m_globalMat->find(gkey);

          if(matrixIter == m_globalMat->end())
          {
              return 0;
          }
          else
          {
              return matrixIter->second->GetNumNonZeroEntries();
          }

          return 0;
      }


      /**
       *
       */
      void ContField3D::v_ClearGlobalLinSysManager(void)
      {
          m_globalLinSysManager.ClearManager("GlobalLinSys");
      }

  } //end of namespace
} //end of namespace
