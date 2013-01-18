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
#include <MultiRegions/AssemblyMap/AssemblyMapCG3D.h>

namespace Nektar
{
  namespace MultiRegions
  {

        ContField3D::ContField3D():
            DisContField3D(),
            m_locToGloMap(),
            m_globalMat(),
            m_globalLinSysManager(
                    boost::bind(&ContField3D::GenGlobalLinSys, this, _1),
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
                                 const std::string &variable):
                DisContField3D(pSession,graph3D,variable,false),
                m_globalMat(MemoryManager<GlobalMatrixMap>::AllocateSharedPtr()),
                m_globalLinSysManager(
                        boost::bind(&ContField3D::GenGlobalLinSys, this, _1),
                        std::string("GlobalLinSys"))
        {
            SpatialDomains::BoundaryConditions bcs(m_session, graph3D);
            
            m_locToGloMap = MemoryManager<AssemblyMapCG3D>::AllocateSharedPtr(
                m_session,m_ncoeffs,*this,m_bndCondExpansions,m_bndConditions,
                m_periodicVertices,m_periodicEdges,m_periodicFaces);
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
         * @param   In          Existing ContField2D object used to provide the
         *                      local to global mapping information and
         *                      global solution type.
         * @param   graph2D     A mesh, containing information about the domain
         *                      and the spectral/hp element expansion.
         * @param   bcs         The boundary conditions.
         * @param   bc_loc
         */
        ContField3D::ContField3D(const ContField3D &In,
                                 const SpatialDomains::MeshGraphSharedPtr &graph3D,
                                 const std::string &variable):
	    DisContField3D(In,graph3D,variable,false),
            m_globalMat   (MemoryManager<GlobalMatrixMap>::AllocateSharedPtr()),
            m_globalLinSysManager(boost::bind(&ContField3D::GenGlobalLinSys, this, _1),
                                  std::string("GlobalLinSys"))

        {
            if(!SameTypeOfBoundaryConditions(In))
            {
                SpatialDomains::BoundaryConditions bcs(m_session, graph3D);

                m_locToGloMap = MemoryManager<AssemblyMapCG3D>::AllocateSharedPtr(
                    m_session,m_ncoeffs,*this,m_bndCondExpansions,m_bndConditions,
                    m_periodicVertices, m_periodicEdges, m_periodicFaces);

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
                                      Array<OneD,       NekDouble> &outarray,
                                     CoeffState coeffstate)
        {
            if(coeffstate == eGlobal)
            {
                bool doGlobalOp = m_globalOptParam->DoGlobalMatOp(
                                                        StdRegions::eBwdTrans);

                if(doGlobalOp)
                {
                    GlobalMatrixKey gkey(StdRegions::eBwdTrans,m_locToGloMap);
                    GlobalMatrixSharedPtr mat = GetGlobalMatrix(gkey);
                    mat->Multiply(inarray,outarray);
                }
                else
                {
                    Array<OneD, NekDouble> wsp(m_ncoeffs);
                    GlobalToLocal(inarray,wsp);
                    BwdTrans_IterPerExp(wsp,outarray);
                }
            }
            else
            {
                BwdTrans_IterPerExp(inarray,outarray);
            }
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
                                            Array<OneD, NekDouble> &outarray,
                                            CoeffState coeffstate)
        {
            if(coeffstate == eGlobal)
            {
                bool doGlobalOp = m_globalOptParam->DoGlobalMatOp(
                                                StdRegions::eIProductWRTBase);

                if(doGlobalOp)
                {
                    GlobalMatrixKey gkey(StdRegions::eIProductWRTBase,
                                         m_locToGloMap);
                    GlobalMatrixSharedPtr mat = GetGlobalMatrix(gkey);
                    mat->Multiply(inarray,outarray);
                    m_locToGloMap->UniversalAssemble(outarray);
                }
                else
                {
                    Array<OneD, NekDouble> wsp(m_ncoeffs);
                    IProductWRTBase_IterPerExp(inarray,wsp);
                    Assemble(wsp,outarray);
                }
            }
            else
            {
                IProductWRTBase_IterPerExp(inarray,outarray);
            }
        }
      
      
      void ContField3D::v_FwdTrans(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,       NekDouble> &outarray,
                                   CoeffState coeffstate)
      {
          // Inner product of forcing
          int contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
          Array<OneD,NekDouble> wsp(contNcoeffs);
          IProductWRTBase(inarray,wsp,eGlobal);
          
          // Solve the system
          GlobalLinSysKey key(StdRegions::eMass, m_locToGloMap);
          
          if(coeffstate == eGlobal)
          {
              GlobalSolve(key,wsp,outarray);
          }
          else
          {
              Array<OneD,NekDouble> tmp(contNcoeffs,0.0);
              GlobalSolve(key,wsp,tmp);
              GlobalToLocal(tmp,outarray);
          }
      }
      
      
      void ContField3D::v_MultiplyByInvMassMatrix(
                                          const Array<OneD, const NekDouble> &inarray,
                                          Array<OneD,       NekDouble> &outarray,
                                          CoeffState coeffstate)
          
      {
          int contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
          GlobalLinSysKey key(StdRegions::eMass, m_locToGloMap);
          
          if(coeffstate == eGlobal)
          {
              if(inarray.data() == outarray.data())
              {
                  Array<OneD, NekDouble> tmp(contNcoeffs,0.0);
                  Vmath::Vcopy(contNcoeffs,inarray,1,tmp,1);
                  GlobalSolve(key,tmp,outarray);
              }
              else
              {
                  GlobalSolve(key,inarray,outarray);
              }
          }
          else
          {
              Array<OneD, NekDouble> globaltmp(contNcoeffs,0.0);
              
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
      
      
      void ContField3D::GenerateDirBndCondForcing(const GlobalLinSysKey &key,
                                                  Array<OneD, NekDouble> &inout,
                                                  Array<OneD, NekDouble> &outarray)
      {
          int bndcnt=0;
          const Array<OneD,const int>& map  = m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsMap();
          NekDouble sign;
          
          for(int i = 0; i < m_bndCondExpansions.num_elements(); ++i)
          {
              if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
              {
                  const Array<OneD,const NekDouble>& coeffs = m_bndCondExpansions[i]->GetCoeffs();
                  for(int j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); ++j)
                  {
                      sign = m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsSign(bndcnt);
                      inout[map[bndcnt++]] = sign * coeffs[j];
                  }
              }
              else
              {
                  bndcnt += m_bndCondExpansions[i]->GetNcoeffs();
              }
            }
          GeneralMatrixOp(key,inout,outarray,eGlobal);
      }
      


      // Note inout contains initial guess and final output.
      void ContField3D::GlobalSolve(
          const GlobalLinSysKey              &key,
          const Array<OneD, const NekDouble> &rhs,
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
              LinSys->Solve(rhs,inout,m_locToGloMap,dirForcing);
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
            GlobalMatrixMap::iterator matrixIter = m_globalMat->find(mkey);
            
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
          int bndcnt      = 0;
          int nDir        = m_locToGloMap->GetNumGlobalDirBndCoeffs();
          int contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
          
          NekDouble sign;
          const Array<OneD,const int> &bndMap = 
              m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsMap();
          
          Array<OneD, NekDouble> tmp(
              m_locToGloMap->GetNumGlobalBndCoeffs(), 0.0);

          // Fill in Dirichlet coefficients that are to be sent to other
          // processors.
          map<int, vector<pair<int, int> > > &extraDirDofs = 
              m_locToGloMap->GetExtraDirDofs();
          map<int, vector<pair<int, int> > >::iterator it;
          for (it = extraDirDofs.begin(); it != extraDirDofs.end(); ++it)
          {
              for (i = 0; i < it->second.size(); ++i)
              {
                  tmp[it->second.at(i).second] = 
                      m_bndCondExpansions[it->first]->GetCoeffs()[
                          it->second.at(i).first];
              }
          }
          m_locToGloMap->UniversalAssembleBnd(tmp);
          
          // Now fill in all other Dirichlet coefficients.
          for(i = 0; i < m_bndCondExpansions.num_elements(); ++i)
          {
              if(m_bndConditions[i]->GetBoundaryConditionType() == 
                     SpatialDomains::eDirichlet)
              {
                  const Array<OneD,const NekDouble>& coeffs = 
                      m_bndCondExpansions[i]->GetCoeffs();
                  for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); ++j)
                  {
                      sign = m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsSign(
                          bndcnt);
                      tmp[bndMap[bndcnt++]] = sign * coeffs[j];
                  }
              }
              else
              {
                  bndcnt += m_bndCondExpansions[i]->GetNcoeffs();
              }
          }
          
          Vmath::Vcopy(nDir, tmp, 1, outarray, 1);
      }          

      void ContField3D::v_LocalToGlobal(void)
      {
          m_locToGloMap->LocalToGlobal(m_coeffs, m_coeffs);
      }

      void ContField3D::v_GlobalToLocal(void)
      {
          m_locToGloMap->GlobalToLocal(m_coeffs, m_coeffs);
      }
      


      void ContField3D::v_HelmSolve(
                                    const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD,       NekDouble> &outarray,
                                    const FlagList &flags,
                                    const StdRegions::ConstFactorMap &factors,
                                    const StdRegions::VarCoeffMap &varcoeff,
                                    const Array<OneD, const NekDouble> &dirForcing)
      {
          // Inner product of forcing
          int contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
          Array<OneD,NekDouble> wsp(contNcoeffs);
          IProductWRTBase(inarray,wsp,eGlobal);
          // Note -1.0 term necessary to invert forcing function to
          // be consistent with matrix definition
          Vmath::Neg(contNcoeffs, wsp, 1);
          
          // Forcing function with weak boundary conditions
          int i,j;
          int bndcnt = 0;
          NekDouble sign;
          Array<OneD, NekDouble> gamma(contNcoeffs, 0.0);
          for(i = 0; i < m_bndCondExpansions.num_elements(); ++i)
          {
              if(m_bndConditions[i]->GetBoundaryConditionType() != SpatialDomains::eDirichlet)
              {
                  for(j = 0; j < (m_bndCondExpansions[i])->GetNcoeffs(); j++)
                  {
                      sign = m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsSign(bndcnt);
                      gamma[m_locToGloMap->GetBndCondCoeffsToGlobalCoeffsMap(bndcnt++)] +=
                          sign * (m_bndCondExpansions[i]->GetCoeffs())[j];
                  }
              }
              else
              {
                  bndcnt += m_bndCondExpansions[i]->GetNcoeffs();
              }
          }
          m_locToGloMap->UniversalAssemble(gamma);
          
          // Add weak boundary conditions to forcing
          Vmath::Vadd(contNcoeffs, wsp, 1, gamma, 1, wsp, 1);
          
          // Solve the system
          GlobalLinSysKey key(StdRegions::eHelmholtz, m_locToGloMap, factors,varcoeff);
          
          if(flags.isSet(eUseGlobal))
          {
              GlobalSolve(key,wsp,outarray,dirForcing);
          }
          else
          {
              Array<OneD,NekDouble> tmp(contNcoeffs, 0.0);
              GlobalSolve(key,wsp,tmp,dirForcing);
              GlobalToLocal(tmp,outarray);
          }
      }
      
      void ContField3D::v_GeneralMatrixOp(
          const GlobalMatrixKey             &gkey,
          const Array<OneD,const NekDouble> &inarray,
                Array<OneD,      NekDouble> &outarray,
          CoeffState                         coeffstate)
      {
          if(coeffstate == eGlobal)
          {
              bool doGlobalOp = m_globalOptParam->DoGlobalMatOp(gkey.GetMatrixType());
              
              if(doGlobalOp)
              {
                  GlobalMatrixSharedPtr mat = GetGlobalMatrix(gkey);
                  mat->Multiply(inarray,outarray);
                  m_locToGloMap->UniversalAssemble(outarray);
              }
              else
              {
                  Array<OneD,NekDouble> tmp1(2*m_ncoeffs);
                  Array<OneD,NekDouble> tmp2(tmp1+m_ncoeffs);
                  GlobalToLocal(inarray,tmp1);
                  GeneralMatrixOp_IterPerExp(gkey,tmp1,tmp2);
                  Assemble(tmp2,outarray);
              }
          }
          else
          {
              GeneralMatrixOp_IterPerExp(gkey,inarray,outarray);
          }
      }
      
      int ContField3D::GetGlobalMatrixNnz(const GlobalMatrixKey &gkey)
      {
          ASSERTL1(gkey.LocToGloMapIsDefined(),
                   "To use method must have a AssemblyMap "
                   "attached to key");
          
          GlobalMatrixMap::iterator matrixIter = m_globalMat->find(gkey);
          
          if(matrixIter == m_globalMat->end())
          {
              return 0;
          }
          else
          {
              return matrixIter->second->GetMatrix()->GetNumNonZeroEntries();
          }
          
          return 0;
      }
      
  } //end of namespace
} //end of namespace
