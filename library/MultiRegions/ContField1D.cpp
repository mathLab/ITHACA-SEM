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
            DisContField1D(),
            m_locToGloMap(),
            m_globalLinSysManager(
                std::bind(&ContField1D::GenGlobalLinSys, this, std::placeholders::_1),
                std::string("GlobalLinSys"))
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
                      const Collections::ImplementationType ImpType):
            DisContField1D(pSession,graph1D,variable,false,ImpType),
            m_locToGloMap(),
            m_globalLinSysManager(
                std::bind(&ContField1D::GenGlobalLinSys, this, std::placeholders::_1),
                std::string("GlobalLinSys"))
        {
            SpatialDomains::BoundaryConditions bcs(pSession, graph1D);

            m_locToGloMap = MemoryManager<AssemblyMapCG>
                ::AllocateSharedPtr(m_session,m_ncoeffs,*this,
                                    m_bndCondExpansions,
                                    m_bndConditions,
                                    false,
                                    variable,
                                    m_periodicVerts);
        }


        /**
         * Constructs a 1D continuous field as a copy of an existing field
         * including all boundary conditions.
         * @param   In          Existing continuous field to duplicate.
         */
        ContField1D::ContField1D(const ContField1D &In):
            DisContField1D(In),
            m_locToGloMap(In.m_locToGloMap),
            m_globalLinSysManager(
                std::bind(&ContField1D::GenGlobalLinSys, this, std::placeholders::_1),
                std::string("GlobalLinSys"))
        {
        }

        /**
         * Constructs a 1D continuous field as a copy of an existing explist 1D field
         * and adding all the boundary conditions.
         * @param   In          Existing explist1D field .
         */
         ContField1D::ContField1D(const LibUtilities::SessionReaderSharedPtr &pSession, const ExpList1D & In):
            DisContField1D(In),
            m_locToGloMap(),
            m_globalLinSysManager(
                std::bind(&ContField1D::GenGlobalLinSys, this, std::placeholders::_1),
                std::string("GlobalLinSys"))
        {
            m_locToGloMap = MemoryManager<AssemblyMapCG>
                ::AllocateSharedPtr(pSession, m_ncoeffs, In);

        }

        /**
         *
         */
        ContField1D::~ContField1D()
        {
        }


        /**
         * Given a function \f$f(x)\f$ defined at the quadrature
         * points, this function determines the unknown global coefficients
         * \f$\boldsymbol{\hat{u}}^{\mathcal{H}}\f$ employing a discrete
         * Galerkin projection from physical space to coefficient
         * space. The operation is evaluated by the function #GlobalSolve using
         * the global mass matrix.
         *
         * The values of the function \f$f(x)\f$ evaluated at the
         * quadrature points \f$x_i\f$ should be contained in the
         * variable #m_phys of the ExpList object \a Sin. The resulting global
         * coefficients \f$\hat{u}_g\f$ are stored in the array #m_coeffs.
         *
         * @param   inarray     Discrete \f$f(x)\f$.
         * @param   outarray    Storage for result.
         */
        void ContField1D::FwdTrans(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,       NekDouble> &outarray)
        {

            // Inner product of forcing
            Array<OneD,NekDouble> wsp(m_ncoeffs);
            IProductWRTBase(inarray,wsp);

            // Solve the system
            GlobalLinSysKey key(StdRegions::eMass, m_locToGloMap);

            GlobalSolve(key,wsp,outarray);
        }


        /**
         * Given the coefficients of an expansion, this function evaluates the
         * spectral/hp expansion \f$u^{\delta}(x)\f$ at the quadrature
         * points \f$x_i\f$. This operation is evaluated locally by the
         * function ExpList#BwdTrans.
         *
         * The coefficients of the expansion should be contained in the
         * variable #m_coeffs of the ExpList object \a In. The resulting
         * physical values at the quadrature points \f$u^{\delta}(x_i)\f$ are
         * stored in the array #m_phys.
         *
         * @param   In          An ExpList, containing the local
         *                      coefficients \f$\hat{u}_n^e\f$ in its array
         *                      #m_coeffs.
         */
        void ContField1D::BwdTrans(
                                const Array<OneD, const NekDouble>  &inarray,
                                Array<OneD,       NekDouble>  &outarray)
        {

            BwdTrans_IterPerExp(inarray,outarray);
        }


        /**
         *
         */
        void ContField1D::MultiplyByInvMassMatrix(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,       NekDouble> &outarray)
        {
            GlobalLinSysKey key(StdRegions::eMass, m_locToGloMap);
            GlobalSolve(key,inarray,outarray);
        }


        /**
         * Given a linear system specified by the key \a key,
         * \f[\boldsymbol{M}\boldsymbol{\hat{u}}_g=\boldsymbol{\hat{f}},\f]
         * this function solves this linear system taking into account the
         * boundary conditions specified in the data member
         * #m_bndCondExpansions.
         * Therefore, it adds an array \f$\boldsymbol{\hat{g}}\f$ which
         * represents the non-zero surface integral resulting from the weak
         * boundary conditions (e.g. Neumann boundary conditions) to the right
         * hand side, that is,
         * \f[\boldsymbol{M}\boldsymbol{\hat{u}}_g=\boldsymbol{\hat{f}}+
         * \boldsymbol{\hat{g}}.\f]
         * Furthermore, it lifts the known degrees of freedom which are
         * prescribed by the Dirichlet boundary conditions. As these known
         * coefficients \f$\boldsymbol{\hat{u}}^{\mathcal{D}}\f$ are numbered
         * first in the global coefficient array \f$\boldsymbol{\hat{u}}_g\f$,
         * the linear system can be decomposed as,
         * \f[\left[\begin{array}{cc}
         * \boldsymbol{M}^{\mathcal{DD}}&\boldsymbol{M}^{\mathcal{DH}}\\
         * \boldsymbol{M}^{\mathcal{HD}}&\boldsymbol{M}^{\mathcal{HH}}
         * \end{array}\right]
         * \left[\begin{array}{c}
         * \boldsymbol{\hat{u}}^{\mathcal{D}}\\
         * \boldsymbol{\hat{u}}^{\mathcal{H}}
         * \end{array}\right]=
         * \left[\begin{array}{c}
         * \boldsymbol{\hat{f}}^{\mathcal{D}}\\
         * \boldsymbol{\hat{f}}^{\mathcal{H}}
         * \end{array}\right]+
         * \left[\begin{array}{c}
         * \boldsymbol{\hat{g}}^{\mathcal{D}}\\
         * \boldsymbol{\hat{g}}^{\mathcal{H}}
         * \end{array}\right]
         * \f]
         * which will then be solved for the unknown coefficients
         * \f$\boldsymbol{\hat{u}}^{\mathcal{H}}\f$ as,
         * \f[
         * \boldsymbol{M}^{\mathcal{HH}}\boldsymbol{\hat{u}}^{\mathcal{H}}
         * = \boldsymbol{\hat{f}}^{\mathcal{H}}
         *   +\boldsymbol{\hat{g}}^{\mathcal{H}}
         *   -\boldsymbol{M}^{\mathcal{HD}}\boldsymbol{\hat{u}}^{\mathcal{D}}\f]
         *
         * @param   key         Specifes the linear system to solve.
         * @param   locrhs         Forcing term \f$\boldsymbol{f}\f$.
         * @param   inout       Solution vector \f$\boldsymbol{\hat{u}}\f$.
         * @param   dirForcing  .
         */
        void ContField1D::GlobalSolve(const GlobalLinSysKey &key,
                                const Array<OneD, const NekDouble>& locrhs,
                                      Array<OneD,       NekDouble>& inout,
                                const Array<OneD, const NekDouble>& dirForcing)
        {
            int NumDirBcs   = m_locToGloMap->GetNumGlobalDirBndCoeffs();
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


        /**
         * The function searches the map #m_globalLinSys to see if the global
         * matrix has been created before. If not, it calls the function
         * #GenglobalLinSys to generate the requested global system.
         *
         * @param   mkey        Key specifying the linear system.
         * @returns             Pointer to the required linear system.
         */
        GlobalLinSysSharedPtr ContField1D::GetGlobalLinSys(
                                const GlobalLinSysKey &mkey)
        {
            return m_globalLinSysManager[mkey];
        }

        GlobalLinSysSharedPtr ContField1D::GenGlobalLinSys(
                                const GlobalLinSysKey &mkey)
        {
            ASSERTL1(mkey.LocToGloMapIsDefined(),
                     "To use method must have a AssemblyMap "
                     "attached to key");
            return ExpList::GenGlobalLinSys(mkey, m_locToGloMap);
        }


        /**
         * The operation is evaluated locally (i.e. with respect to all local
         * expansion modes) by the function ExpList#IProductWRTBase. The inner
         * product with respect to the global expansion modes is than obtained
         * by a global assembly operation.
         *
         * The values of the function \f$f(x)\f$ evaluated at the quadrature
         * points \f$x_i\f$ should be contained in the variable #m_phys of the
         * ExpList object \a in. The result is stored in the array
         * #m_coeffs.
         *
         * @param   In          An ExpList, containing the discrete evaluation
         *                      of \f$f(x)\f$ at the quadrature points in its
         *                      array #m_phys.
         */
        void ContField1D::IProductWRTBase(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,       NekDouble> &outarray)
        {
            IProductWRTBase_IterPerExp(inarray,outarray);
        }


        void ContField1D::v_FwdTrans(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,       NekDouble> &outarray)
        {
            FwdTrans(inarray,outarray);
        }

        void ContField1D::v_MultiplyByInvMassMatrix(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,       NekDouble> &outarray)
        {
            MultiplyByInvMassMatrix(inarray,outarray);
        }

        void ContField1D::v_ImposeDirichletConditions(Array<OneD,NekDouble>& outarray)
        {
            const Array<OneD, const int> map= m_locToGloMap->
                GetBndCondCoeffsToLocalCoeffsMap();
            
            for(int i = 0; i < m_bndCondExpansions.size(); ++i)
            {
                if(m_bndConditions[i]->GetBoundaryConditionType() ==
                   SpatialDomains::eDirichlet)
                {
                    outarray[map[i]] = m_bndCondExpansions[i]->GetCoeffs(0);
                }
            }
        }

        /**
         * This operation is evaluated as:
         * \f{tabbing}
         * \hspace{1cm}  \= Do \= $e=$  $1, N_{\mathrm{el}}$ \\
         * \> \> Do \= $i=$  $0,N_m^e-1$ \\
         * \> \> \> $\boldsymbol{\hat{u}}_g[\mbox{map}[e][i]] =
         * \mbox{sign}[e][i] \cdot \boldsymbol{\hat{u}}^{e}[i]$\\
         * \> \> continue\\
         * \> continue
         * \f}
         * where \a map\f$[e][i]\f$ is the mapping array and \a
         * sign\f$[e][i]\f$ is an array of similar dimensions ensuring the
         * correct modal connectivity between the different elements (both
         * these arrays are contained in the data member #m_locToGloMap). This
         * operation is equivalent to the gather operation
         * \f$\boldsymbol{\hat{u}}_g=\mathcal{A}^{-1}\boldsymbol{\hat{u}}_l\f$,
         * where \f$\mathcal{A}\f$ is the
         * \f$N_{\mathrm{eof}}\times N_{\mathrm{dof}}\f$ permutation matrix.
         *
         */
        void ContField1D::v_LocalToGlobal(
            const Array<OneD, const NekDouble> &inarray,
            Array<OneD,NekDouble> &outarray, bool useComm)
        {
            m_locToGloMap->LocalToGlobal(inarray, outarray, useComm);
        }


        void ContField1D::v_LocalToGlobal(bool useComm)
        {
            m_locToGloMap->LocalToGlobal(m_coeffs,m_coeffs, useComm);
        }

        /**
         * This operation is evaluated as:
         * \f{tabbing}
         * \hspace{1cm}  \= Do \= $e=$  $1, N_{\mathrm{el}}$ \\
         * \> \> Do \= $i=$  $0,N_m^e-1$ \\
         * \> \> \> $\boldsymbol{\hat{u}}^{e}[i] = \mbox{sign}[e][i] \cdot
         * \boldsymbol{\hat{u}}_g[\mbox{map}[e][i]]$ \\
         * \> \> continue \\
         * \> continue
         * \f}
         * where \a map\f$[e][i]\f$ is the mapping array and
         * \a sign\f$[e][i]\f$ is an array of similar dimensions ensuring the
         * correct modal connectivity between the different elements (both
         * these arrays are contained in the data member #m_locToGloMap). This
         * operation is equivalent to the scatter operation
         * \f$\boldsymbol{\hat{u}}_l=\mathcal{A}\boldsymbol{\hat{u}}_g\f$, where
         * \f$\mathcal{A}\f$ is the
         * \f$N_{\mathrm{eof}}\times N_{\mathrm{dof}}\f$ permutation matrix.
         *
         */
        void ContField1D::v_GlobalToLocal(
            const Array<OneD, const NekDouble> &inarray,
            Array<OneD,NekDouble> &outarray)
        {
            m_locToGloMap->GlobalToLocal(inarray, outarray);
        }

        void ContField1D::v_GlobalToLocal(void)
        {
            m_locToGloMap->GlobalToLocal(m_coeffs,m_coeffs);
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
            int i;
            for(i = 0; i < m_bndCondExpansions.size(); ++i)
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

        const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>&
                                ContField1D::v_GetBndConditions()
        {
            return GetBndConditions();
        }

        void ContField1D::v_BwdTrans(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,       NekDouble> &outarray)
        {
            BwdTrans(inarray,outarray);
        }

        void ContField1D::v_IProductWRTBase(
                                const Array<OneD, const NekDouble>   &inarray,
                                Array<OneD,       NekDouble>   &outarray)
                                        {
            IProductWRTBase(inarray,outarray);
        }

        /**
         * This is equivalent to the operation:
         * \f[\boldsymbol{M\hat{u}}_g\f]
         * where \f$\boldsymbol{M}\f$ is the global matrix of type specified by
         * \a mkey. After scattering the global array \a inarray to local
         * level, this operation is evaluated locally by the function
         * ExpList#GeneralMatrixOp. The global result is then obtained by a
         * global assembly procedure.
         *
         * @param   mkey        This key uniquely defines the type matrix
         *                      required for the operation.
         * @param   inarray     The vector \f$\boldsymbol{\hat{u}}_g\f$ of size
         *                      \f$N_{\mathrm{dof}}\f$.
         * @param   outarray    The resulting vector of size
         *                      \f$N_{\mathrm{dof}}\f$.
         */
        void ContField1D::v_GeneralMatrixOp(
                                const GlobalMatrixKey                &gkey,
                                const Array<OneD,const NekDouble>    &inarray,
                                Array<OneD,      NekDouble>    &outarray)
        {
            GeneralMatrixOp_IterPerExp(gkey,inarray,outarray);
        }



        /**
         * Reset the GlobalLinSys Manager
         */
        void ContField1D::v_ClearGlobalLinSysManager(void)
        {
            m_globalLinSysManager.ClearManager("GlobalLinSys");
        }

    } // end of namespace
} //end of namespace
