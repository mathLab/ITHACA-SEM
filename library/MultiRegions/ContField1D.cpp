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
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class ContField1D
         * As opposed to the class #ContExpList1D, the class #ContField1D is
         * able to incorporate the boundary conditions imposed to the problem
         * to be solved. Therefore, the class is equipped with three additional
         * data members:
         * - #m_bndCondExpansions
         * - #m_bndTypes
         * - #m_bndCondEquations
         *
         * The first data structure, #m_bndCondExpansions,
         * contains the point Expansion on the boundary,  #m_bndTypes
         * stores information about the type of boundary condition on the
         * different parts of the boundary while #m_bndCondEquations holds the
         * equation of the imposed boundary conditions.
         *
         * Furthermore, in case of Dirichlet boundary conditions,
         * this class is capable of lifting a known solution satisfying these
         * boundary conditions. If we denote the unknown solution by
         * \f$u^{\mathcal{H}}(\boldsymbol{x})\f$ and the known Dirichlet
         * boundary conditions by \f$u^{\mathcal{D}}(\boldsymbol{x})\f$, the
         * expansion then can be decomposed as
         * \f[ u^{\delta}(\boldsymbol{x}_i)=u^{\mathcal{D}}(\boldsymbol{x}_i)+
         * u^{\mathcal{H}}(\boldsymbol{x}_i)=\sum_{n=0}^{N^{\mathcal{D}}-1}
         * \hat{u}_n^{\mathcal{D}}\Phi_n(\boldsymbol{x}_i)+
         * \sum_{n={N^{\mathcal{D}}}}^{N_{\mathrm{dof}}
         *   -1}\hat{u}_n^{\mathcal{H}}
         * \Phi_n(\boldsymbol{x}_i).\f]
         * This lifting is accomplished by ordering the known global degrees of
         * freedom, prescribed by the Dirichlet boundary conditions, first in
         * the global array \f$\boldsymbol{\hat{u}}\f$, that is,
         * \f[\boldsymbol{\hat{u}}=\left[ \begin{array}{c}
         * \boldsymbol{\hat{u}}^{\mathcal{D}}\\
         * \boldsymbol{\hat{u}}^{\mathcal{H}}
         * \end{array} \right].\f]
         * Such kind of expansions are also referred to as continuoous fields.
         * This class should be used when solving 2D problems using a standard
         * Galerkin approach.
         */

        /**
         * Constructs an empty 1D continuous field.
         */
        ContField1D::ContField1D():
           ContField()
        {
        }


        /**
         * Given a mesh \a graph1D, containing information about the domain and
         * the spectral/hp element expansion, this constructor fills the list
         * of local expansions #m_exp with the proper expansions, calculates
         * the total number of quadrature points \f$\boldsymbol{x}_i\f$ and
         * local expansion coefficients \f$\hat{u}^e_n\f$ and allocates
         * memory for the arrays #m_coeffs and #m_phys. Furthermore, it
         * constructs the mapping array (contained in #m_locToGloMap) for the
         * transformation between local elemental level and global level, it
         * calculates the total number global expansion coefficients
         * \f$\hat{u}_n\f$. 
         * The constructor also discretises the boundary conditions, specified
         * by the argument \a bcs, by expressing them in terms of the
         * coefficient of the expansion on  the boundary.
         *
         * @param   graph1D     A 1D mesh, containing information about the
         *                      domain and the spectral/hp element expansion.
         * @param   bcs         The boundary conditions.
         * @param   variable    An optional parameter to indicate for which
         *                      variable the field should be constructed.
         */
        ContField1D::ContField1D(
                      const LibUtilities::SessionReaderSharedPtr &pSession,
                      const SpatialDomains::MeshGraphSharedPtr &graph1D,
                      const std::string &variable,
                      const bool DeclareCoeffPhysArrays,
                      const bool CheckIfSingularSystem,
                      const Collections::ImplementationType ImpType):
            ContField(pSession,graph1D,variable,DeclareCoeffPhysArrays,
                      CheckIfSingularSystem, ImpType)
        {
        }


        /**
         * Constructs a 1D continuous field as a copy of an existing field
         * including all boundary conditions.
         * @param   In          Existing continuous field to duplicate.
         */
        ContField1D::ContField1D(const ContField1D &In,
                                 bool DeclareCoeffPhysArrays):
            ContField(In,DeclareCoeffPhysArrays)
        {
        }

        /**
         * Constructs a 1D continuous field as a copy of an existing
         * explist 1D field and adding all the boundary conditions.
         *
         * @param In Existing explist1D field .
         */
         ContField1D::ContField1D
         (const LibUtilities::SessionReaderSharedPtr &pSession,
          const ExpList & In):
             ContField(pSession,In)
         {
         }

        /**
         *
         */
        ContField1D::~ContField1D()
        {
        }

#if 0

        void ContField1D::v_ImposeDirichletConditions(Array<OneD,NekDouble>& outarray)
        {
            const Array<OneD, const int> map= m_locToGloMap->
                GetBndCondCoeffsToLocalCoeffsMap();
            
            for(int i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                if(m_bndConditions[i]->GetBoundaryConditionType() ==
                   SpatialDomains::eDirichlet)
                {
                    outarray[map[i]] = m_bndCondExpansions[i]->GetCoeffs(0);
                }
            }
        }

        /**
         * Consider the one dimensional Helmholtz equation,
         * \f[\frac{d^2u}{dx^2}-\lambda u(x) = f(x),\f]
         * supplemented with appropriate boundary conditions (which are
         * contained in the data member #m_bndCondExpansions). Applying a
         * \f$C^0\f$ continuous Galerkin discretisation, this equation leads to
         * the following linear system:
         * \f[\left( \boldsymbol{M}+\lambda\boldsymbol{L}\right)
         * \boldsymbol{\hat{u}}_g=\boldsymbol{\hat{f}}\f]
         * where \f$\boldsymbol{M}\f$ and \f$\boldsymbol{L}\f$ are the mass and
         * Laplacian matrix respectively. This function solves the system above
         * for the global coefficients \f$\boldsymbol{\hat{u}}\f$ by a call to
         * the function #GlobalSolve.
         *
         * The values of the function \f$f(x)\f$ evaluated at the
         * quadrature points \f$\boldsymbol{x}_i\f$ should be contained in the
         * variable #m_phys of the ExpList object \a inarray. The resulting
         * global coefficients \f$\boldsymbol{\hat{u}}_g\f$ are stored in the
         * array #m_coeffs.
         *
         * @param   inarray     Input containing forcing function
         *                      \f$\boldsymbol{f}\f$ at the quadrature points.
         * @param   outarray    Output containing the coefficients
         *                      \f$\boldsymbol{u}_g\f$
         * @param   lambda      Parameter value.
         * @param   Sigma       Coefficients of lambda.
         * @param   varcoeff    Variable diffusivity coefficients.
         * @param   dirForcing  Dirichlet Forcing.
         */
        void ContField1D::v_HelmSolve(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray,
                const FlagList &flags,
                const StdRegions::ConstFactorMap &factors,
                const StdRegions::VarCoeffMap &varcoeff,
                const MultiRegions::VarFactorsMap &varfactors,
                const Array<OneD, const NekDouble> &dirForcing,
                const bool PhysSpaceForcing)
        {
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

            const Array<OneD, const int> map= m_locToGloMap->
                GetBndCondCoeffsToLocalCoeffsMap();
            // Add weak boundary conditions to forcing
            for(int i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                if(m_bndConditions[i]->GetBoundaryConditionType() ==
                   SpatialDomains::eNeumann ||
                   m_bndConditions[i]->GetBoundaryConditionType() ==
                   SpatialDomains::eRobin)
                {
                    wsp[map[i]] += m_bndCondExpansions[i]->GetCoeffs(0); 
                }
            }

            // Solve the system
            GlobalLinSysKey key(StdRegions::eHelmholtz,
                                m_locToGloMap,factors,varcoeff,varfactors);


            GlobalSolve(key,wsp,outarray,dirForcing);
        }
#endif

    } // end of namespace
} //end of namespace
