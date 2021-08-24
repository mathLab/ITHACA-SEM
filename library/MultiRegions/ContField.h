///////////////////////////////////////////////////////////////////////////////
//
// File ContField.h
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
// Description: Field definition in tow-dimensions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_CONTFIELD2D_H
#define NEKTAR_LIBS_MULTIREGIONS_CONTFIELD2D_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <SpatialDomains/Conditions.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/DisContField.h>
#include <MultiRegions/GlobalMatrix.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>


namespace Nektar
{
    namespace MultiRegions
    {
        /// This class is the abstraction of a global continuous two-
        /// dimensional spectral/hp element expansion which approximates the
        /// solution of a set of partial differential equations.
        class ContField: public DisContField
        {
        public:
            /// The default constructor.
            MULTI_REGIONS_EXPORT ContField();

            /// This constructor sets up global continuous field based on an
            /// input mesh and boundary conditions.
            MULTI_REGIONS_EXPORT ContField(
                        const LibUtilities::SessionReaderSharedPtr &pSession,
                        const SpatialDomains::MeshGraphSharedPtr &graph2D,
                        const std::string &variable  = "DefaultVar",
                        const bool DeclareCoeffPhysArrays = true,
                        const bool CheckIfSingularSystem = false,
                        const Collections::ImplementationType ImpType
                                             = Collections::eNoImpType);

            /// Construct a global continuous field with solution type based on
            /// another field but using a separate input mesh and boundary
            /// conditions.
            MULTI_REGIONS_EXPORT ContField(const ContField &In,
                        const SpatialDomains::MeshGraphSharedPtr &graph2D,
                        const std::string &variable,
                        const bool DeclareCoeffPhysArrays = true,
                        const bool CheckIfSingularSystem = false);

            /// The copy constructor.
            MULTI_REGIONS_EXPORT ContField
                (const ContField &In, bool DeclareCoeffPhysArrays = true);

            /// Copy constructor.
            MULTI_REGIONS_EXPORT ContField
                (const LibUtilities::SessionReaderSharedPtr &pSession,
                 const ExpList & In);

            /// The default destructor.
            MULTI_REGIONS_EXPORT virtual ~ContField();

            /// Assembles the global coefficients \f$\boldsymbol{\hat{u}}_g\f$
            /// from the local coefficients \f$\boldsymbol{\hat{u}}_l\f$.
            inline void Assemble();

            /// Assembles the global coefficients \f$\boldsymbol{\hat{u}}_g\f$
            /// from the local coefficients \f$\boldsymbol{\hat{u}}_l\f$.
            inline void Assemble(
                            const Array<OneD, const NekDouble> &inarray,
                                  Array<OneD,NekDouble> &outarray) const;

            /// Returns the map from local to global level.
            inline const AssemblyMapCGSharedPtr& GetLocalToGlobalMap()
                const;

            /// Calculates the inner product of a function
            /// \f$f(\boldsymbol{x})\f$ with respect to all <em>global</em>
            /// expansion modes \f$\phi_n^e(\boldsymbol{x})\f$.
            inline void IProductWRTBase(
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD, NekDouble> &outarray);

            /// Performs the global forward transformation of a function
            /// \f$f(\boldsymbol{x})\f$, subject to the boundary conditions
            /// specified.
            MULTI_REGIONS_EXPORT void FwdTrans(const Array<OneD,
                                               const NekDouble> &inarray,
                                               Array<OneD,NekDouble> &outarray);

            /// Performs the backward transformation of the spectral/hp
            /// element expansion.
            inline void BwdTrans(
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,       NekDouble> &outarray);

            /// Multiply a solution by the inverse mass matrix.
            MULTI_REGIONS_EXPORT void MultiplyByInvMassMatrix(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,  NekDouble> &outarray);

            /// Solves the two-dimensional Laplace equation, subject to the
            /// boundary conditions specified.
            MULTI_REGIONS_EXPORT void LaplaceSolve(
                               const Array<OneD, const NekDouble> &inarray,
                               Array<OneD,       NekDouble> &outarray,
                               const Array<OneD, const NekDouble> &dirForcing
                               = NullNekDouble1DArray,
                               const Array<OneD,       Array<OneD,NekDouble> >&
                               variablecoeffs = NullNekDoubleArrayOfArray,
                               NekDouble time = 0.0);

            /// Compute the eigenvalues of the linear advection operator.
            MULTI_REGIONS_EXPORT void LinearAdvectionEigs(const NekDouble ax,
                                     const NekDouble ay,
                                     Array<OneD, NekDouble> &Real,
                                     Array<OneD, NekDouble> &Imag,
                                     Array<OneD, NekDouble> &Evecs
                                                          = NullNekDouble1DArray);

            /// Returns the boundary conditions expansion.
            inline const Array<OneD,const MultiRegions::ExpListSharedPtr>&
                                                        GetBndCondExpansions();

            /// Returns the boundary conditions.
            inline const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>&
                GetBndConditions();

            inline int GetGlobalMatrixNnz(const GlobalMatrixKey &gkey);

        protected:

            //private:
            /// (A shared pointer to) the object which contains all the
            /// required information for the transformation from local to
            /// global degrees of freedom.
            AssemblyMapCGSharedPtr m_locToGloMap;

            /// (A shared pointer to) a list which collects all the global
            /// matrices being assembled, such that they should be constructed
            /// only once.
            GlobalMatrixMapShPtr   m_globalMat;

            /// A manager which collects all the global
            /// linear systems being assembled, such that they should be
            /// constructed only once.
            LibUtilities::NekManager<GlobalLinSysKey, GlobalLinSys> m_globalLinSysManager;

            /// Solves the linear system specified by the key \a key.
            MULTI_REGIONS_EXPORT void GlobalSolve(const GlobalLinSysKey &key,
                             const Array<OneD, const  NekDouble> &rhs,
                                   Array<OneD,        NekDouble> &inout,
                             const Array<OneD, const NekDouble> &dirForcing
                                                        = NullNekDouble1DArray);

            /// Returns the global matrix specified by \a mkey.
            MULTI_REGIONS_EXPORT GlobalMatrixSharedPtr
                GetGlobalMatrix(const GlobalMatrixKey &mkey);

            /// Returns the linear system specified by the key \a mkey.
            MULTI_REGIONS_EXPORT GlobalLinSysSharedPtr
                GetGlobalLinSys(const GlobalLinSysKey &mkey);

            MULTI_REGIONS_EXPORT GlobalLinSysSharedPtr
                GenGlobalLinSys(const GlobalLinSysKey &mkey);

            /// Impose the Dirichlet Boundary Conditions on outarray
            MULTI_REGIONS_EXPORT virtual void
                v_ImposeDirichletConditions(Array<OneD,NekDouble>& outarray);

            MULTI_REGIONS_EXPORT virtual void
                v_FillBndCondFromField();

            MULTI_REGIONS_EXPORT virtual void
                v_FillBndCondFromField(const int nreg);

            /// Gathers the global coefficients \f$\boldsymbol{\hat{u}}_g\f$
            /// from the local coefficients \f$\boldsymbol{\hat{u}}_l\f$.
            MULTI_REGIONS_EXPORT virtual void v_LocalToGlobal(
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray,
                bool useComm);

            MULTI_REGIONS_EXPORT virtual void v_LocalToGlobal(bool useComm);

            /// Scatters from the global coefficients
            /// \f$\boldsymbol{\hat{u}}_g\f$ to the local coefficients
            /// \f$\boldsymbol{\hat{u}}_l\f$.
            MULTI_REGIONS_EXPORT virtual void v_GlobalToLocal(
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray);

            MULTI_REGIONS_EXPORT virtual void v_GlobalToLocal(void);

            // /// Template method virtual forwarder for FwdTrans().
            // MULTI_REGIONS_EXPORT virtual void v_BwdTrans(
            //                     const Array<OneD, const NekDouble> &inarray,
            //                     Array<OneD,       NekDouble> &outarray);


            /// Template method virtual forwarder for FwdTrans().
            MULTI_REGIONS_EXPORT virtual void v_FwdTrans(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,       NekDouble> &outarray);

            /// Template method virtual forwarded for SmoothField().
            MULTI_REGIONS_EXPORT virtual void v_SmoothField(
                                      Array<OneD,NekDouble> &field);

            /// Template method virtual forwarder for MultiplyByInvMassMatrix().
            MULTI_REGIONS_EXPORT virtual void v_MultiplyByInvMassMatrix(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,       NekDouble> &outarray);

            /// Solves the two-dimensional Helmholtz equation, subject to the
            /// boundary conditions specified.
            MULTI_REGIONS_EXPORT virtual void v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                    const StdRegions::ConstFactorMap &factors,
                    const StdRegions::VarCoeffMap &varcoeff,
                    const MultiRegions::VarFactorsMap &varfactors,
                    const Array<OneD, const NekDouble> &dirForcing,
                    const bool PhysSpaceForcing);

            /// Calculates the result of the multiplication of a global
            /// matrix of type specified by \a mkey with a vector given by \a
            /// inarray.
            virtual void v_GeneralMatrixOp(
                   const GlobalMatrixKey             &gkey,
                   const Array<OneD,const NekDouble> &inarray,
                   Array<OneD,      NekDouble> &outarray);

            // Solve the linear advection problem assuming that m_coeffs
            // vector contains an intial estimate for solution
            MULTI_REGIONS_EXPORT virtual void
                v_LinearAdvectionDiffusionReactionSolve
                     (const Array<OneD, Array<OneD, NekDouble> > &velocity,
                      const Array<OneD, const NekDouble> &inarray,
                      Array<OneD, NekDouble> &outarray,
                      const NekDouble lambda,
                      const Array<OneD, const NekDouble>&
                      dirForcing = NullNekDouble1DArray);

            // Solve the linear advection problem assuming that m_coeff
            // vector contains an intial estimate for solution
            MULTI_REGIONS_EXPORT void
                v_LinearAdvectionReactionSolve
                        (const Array<OneD, Array<OneD, NekDouble> > &velocity,
                         const Array<OneD, const NekDouble> &inarray,
                         Array<OneD, NekDouble> &outarray,
                         const NekDouble lambda,
                         const Array<OneD, const NekDouble>&
                         dirForcing = NullNekDouble1DArray);

            /// Template method virtual forwarder for GetBndConditions().
            MULTI_REGIONS_EXPORT virtual const Array<OneD,const SpatialDomains
                                ::BoundaryConditionShPtr>& v_GetBndConditions();
            MULTI_REGIONS_EXPORT virtual void v_ClearGlobalLinSysManager(void);
        };

        typedef std::shared_ptr<ContField>      ContFieldSharedPtr;

        /**
         * This operation is evaluated as:
         * \f{tabbing}
         * \hspace{1cm}  \= Do \= $e=$  $1, N_{\mathrm{el}}$ \\
         * \> \> Do \= $i=$  $0,N_m^e-1$ \\
         * \> \> \> $\boldsymbol{\hat{u}}_g[\mbox{map}[e][i]] =
         * \boldsymbol{\hat{u}}_g[\mbox{map}[e][i]]+\mbox{sign}[e][i] \cdot
         * \boldsymbol{\hat{u}}^{e}[i]$\\
         * \> \> continue\\
         * \> continue
         * \f}
         * where \a map\f$[e][i]\f$ is the mapping array and \a
         * sign\f$[e][i]\f$ is an array of similar dimensions ensuring the
         * correct modal connectivity between the different elements (both
         * these arrays are contained in the data member #m_locToGloMap). This
         * operation is equivalent to the gather operation
         * \f$\boldsymbol{\hat{u}}_g=\mathcal{A}^{T}\boldsymbol{\hat{u}}_l\f$,
         * where \f$\mathcal{A}\f$ is the
         * \f$N_{\mathrm{eof}}\times N_{\mathrm{dof}}\f$ permutation matrix.
         *
         * @note    The array #m_coeffs should be filled with the local
         *          coefficients \f$\boldsymbol{\hat{u}}_l\f$ and that the
         *          resulting global coefficients \f$\boldsymbol{\hat{u}}_g\f$
         *          will be stored in #m_coeffs.
         */
        inline void ContField::Assemble()
        {
            m_locToGloMap->Assemble(m_coeffs,m_coeffs);
        }

        /**
         * This operation is evaluated as:
         * \f{tabbing}
         * \hspace{1cm}  \= Do \= $e=$  $1, N_{\mathrm{el}}$ \\
         * \> \> Do \= $i=$  $0,N_m^e-1$ \\
         * \> \> \> $\boldsymbol{\hat{u}}_g[\mbox{map}[e][i]] =
         * \boldsymbol{\hat{u}}_g[\mbox{map}[e][i]]+\mbox{sign}[e][i] \cdot
         * \boldsymbol{\hat{u}}^{e}[i]$\\
        *  \> \> continue\\
         * \> continue
         * \f}
         * where \a map\f$[e][i]\f$ is the mapping array and \a
         * sign\f$[e][i]\f$ is an array of similar dimensions ensuring the
         * correct modal connectivity between the different elements (both
         * these arrays are contained in the data member #m_locToGloMap). This
         * operation is equivalent to the gather operation
         * \f$\boldsymbol{\hat{u}}_g=\mathcal{A}^{T}\boldsymbol{\hat{u}}_l\f$,
         * where \f$\mathcal{A}\f$ is the
         * \f$N_{\mathrm{eof}}\times N_{\mathrm{dof}}\f$ permutation matrix.
         *
         * @param   inarray     An array of size \f$N_\mathrm{eof}\f$
         *                      containing the local degrees of freedom
         *                      \f$\boldsymbol{x}_l\f$.
         * @param   outarray    The resulting global degrees of freedom
         *                      \f$\boldsymbol{x}_g\f$ will be stored in this
         *                      array of size \f$N_\mathrm{dof}\f$.
         */
        inline void ContField::Assemble(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,NekDouble> &outarray) const
        {
            m_locToGloMap->Assemble(inarray,outarray);
        }


        inline const AssemblyMapCGSharedPtr&
            ContField::GetLocalToGlobalMap() const
        {
            return  m_locToGloMap;
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
        inline void ContField::IProductWRTBase(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD, NekDouble> &outarray)

        {
            IProductWRTBase_IterPerExp(inarray,outarray);
        }

        /**
         * Given the coefficients of an expansion, this function evaluates the
         * spectral/hp expansion \f$u^{\delta}(\boldsymbol{x})\f$ at the
         * quadrature points \f$\boldsymbol{x}_i\f$. This operation is
         * evaluated locally by the function ExpList#BwdTrans.
         *
         * The coefficients of the expansion should be contained in the variable
         * #m_coeffs of the ExpList object \a In. The resulting physical values
         * at the quadrature points \f$u^{\delta}(\boldsymbol{x}_i)\f$ are
         * stored in the array #m_phys.
         *
         * @param   In          An ExpList, containing the local coefficients
         *                      \f$\hat{u}_n^e\f$ in its array #m_coeffs.
         */
        inline void ContField::BwdTrans(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,       NekDouble> &outarray)
        {
            BwdTrans_IterPerExp(inarray,outarray);
        }

        inline const Array<OneD,const MultiRegions::ExpListSharedPtr>&
                                            ContField::GetBndCondExpansions()
        {
            return m_bndCondExpansions;
        }

        inline const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>&
                                                ContField::GetBndConditions()
        {
            return m_bndConditions;
        }

        inline int ContField::GetGlobalMatrixNnz(const GlobalMatrixKey &gkey)
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

#endif // MULTIERGIONS_CONTFIELD2D_H
