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
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>

using namespace std;

namespace Nektar
{
  namespace MultiRegions
  {

        ContField3D::ContField3D():
            ContField()
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
        ContField3D::ContField3D
            (const LibUtilities::SessionReaderSharedPtr &pSession,
             const SpatialDomains::MeshGraphSharedPtr &graph3D,
             const std::string &variable,
             const bool DeclareCoeffPhysArrays,
             const bool CheckIfSingularSystem,
             const Collections::ImplementationType ImpType):
                ContField(pSession,graph3D,variable,DeclareCoeffPhysArrays,
                          CheckIfSingularSystem,ImpType)
        {
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
	    ContField(In,graph3D,variable,CheckIfSingularSystem)
        {
        }


        ContField3D::ContField3D(const ContField3D &In,
                                 bool DeclareCoeffPhysArrays):
            ContField(In,DeclareCoeffPhysArrays)
        {
        }
      
      
        ContField3D::~ContField3D()
        {
        }

      
      
      void ContField3D::v_HelmSolve(
                                    const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD,       NekDouble> &outarray,
                                    const FlagList &flags,
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
          for(i = 0; i < m_bndCondExpansions.num_elements(); ++i)
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
      

  } //end of namespace
} //end of namespace
