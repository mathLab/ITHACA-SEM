///////////////////////////////////////////////////////////////////////////////
//
// File ContField1D.h
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
// Description: Field definition in one-dimension
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_CONTFIELD1D_H
#define NEKTAR_LIBS_MULTIREGIONS_CONTFIELD1D_H


#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/DisContField1D.h>
#include <MultiRegions/GlobalMatrix.h>
#include <MultiRegions/GlobalLinSys.h>
#include <SpatialDomains/Conditions.h>
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /// Abstraction of a global continuous one-dimensional spectral/hp
        /// element expansion which approximates the solution of a set of
        /// partial differential equations.
        class ContField1D: public DisContField1D
        {
        public:
            /// Default constructor.
            MULTI_REGIONS_EXPORT ContField1D();

            /// Set up global continuous field based on an input mesh and
            /// boundary conditions.
            MULTI_REGIONS_EXPORT ContField1D(
                        const LibUtilities::SessionReaderSharedPtr &pSession,
                        const SpatialDomains::MeshGraphSharedPtr &graph1D,
                        const std::string &variable  = "DefaultVar",
                        const Collections::ImplementationType ImpType
                                  = Collections::eNoImpType);

            /// Copy constructor.
            MULTI_REGIONS_EXPORT ContField1D(const ContField1D &In);

            /// Copy constructor.
            MULTI_REGIONS_EXPORT ContField1D( const LibUtilities::SessionReaderSharedPtr &pSession,
                                              const ExpList1D & In);
            /// Destructor
            MULTI_REGIONS_EXPORT virtual ~ContField1D();

            /// Perform global forward transformation of a function \f$f(x)\f$,
            //  subject to the boundary conditions specified.
            MULTI_REGIONS_EXPORT void FwdTrans(const Array<OneD,
                                               const NekDouble> &inarray,
                                               Array<OneD,NekDouble> &outarray);

            /// This function performs the backward transformation of the
            /// spectral/hp element expansion.
            MULTI_REGIONS_EXPORT void BwdTrans(const Array<OneD,
                                               const NekDouble> &inarray,
                                               Array<OneD,NekDouble> &outarray);

            ///
            MULTI_REGIONS_EXPORT void MultiplyByInvMassMatrix(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,       NekDouble> &outarray);

            /// Return the boundary conditions expansion.
            // inline
            MULTI_REGIONS_EXPORT const Array<OneD,const MultiRegions::ExpListSharedPtr>&
                                                     GetBndCondExpansions();

            // inline
            MULTI_REGIONS_EXPORT const Array<OneD,const SpatialDomains
                                ::BoundaryConditionShPtr>& GetBndConditions();


            /// Assembles the global coefficients \f$\boldsymbol{\hat{u}}_g\f$
            /// from the local coefficients \f$\boldsymbol{\hat{u}}_l\f$.
            // inline
            MULTI_REGIONS_EXPORT void Assemble();

            /// Assembles the global coefficients \f$\boldsymbol{\hat{u}}_g\f$
            /// from the local coefficients \f$\boldsymbol{\hat{u}}_l\f$.
            // inline
            MULTI_REGIONS_EXPORT void Assemble(const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,NekDouble> &outarray);

            /// Returns the map from local to global level.
            // inline
            MULTI_REGIONS_EXPORT const AssemblyMapCGSharedPtr& GetLocalToGlobalMap() const;

            /// Calculates the inner product of a function \f$f(x)\f$ with
            /// respect to all <em>global</em> expansion modes
            /// \f$\phi_n^e(x)\f$.
            MULTI_REGIONS_EXPORT void IProductWRTBase(const Array<OneD, const NekDouble> &inarray,
                                                      Array<OneD, NekDouble> &outarray);

            /// Calculates the result of the multiplication of a global matrix
            /// of type specified by \a mkey with a vector given by \a inarray.
            MULTI_REGIONS_EXPORT void GeneralMatrixOp(const GlobalMatrixKey             &gkey,
                                const Array<OneD,const NekDouble> &inarray,
                                Array<OneD,      NekDouble> &outarray);   

        protected:
            /// (A shared pointer to) the object which contains all the required
            /// information for the transformation from local to global degrees
            /// of freedom.
            AssemblyMapCGSharedPtr m_locToGloMap;

            /// matrices being assembled, such that they should be constructed
            /// only once.
            GlobalMatrixMapShPtr            m_globalMat;

            /// A manager which collects all the global
            /// linear systems being assembled, such that they should be
            /// constructed only once.
            LibUtilities::NekManager<GlobalLinSysKey, GlobalLinSys> m_globalLinSysManager;

        private:
            /// Returns the linear system specified by \a mkey.
            GlobalLinSysSharedPtr GetGlobalLinSys(const GlobalLinSysKey &mkey);

            GlobalLinSysSharedPtr GenGlobalLinSys(const GlobalLinSysKey &mkey);

            /// Solve the linear system specified by the key \a key.
            void GlobalSolve(   const GlobalLinSysKey &key,
                                const Array<OneD, const NekDouble> &rhs,
                                      Array<OneD,       NekDouble> &inout,
                                const Array<OneD, const NekDouble> &dirForcing
                                                        = NullNekDouble1DArray);

            /// Perform a forward transform
            virtual void v_FwdTrans(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,       NekDouble> &outarray);

            virtual void v_MultiplyByInvMassMatrix(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,       NekDouble> &outarray);

            /// Impose the Dirichlet Boundary Conditions on outarray 
            MULTI_REGIONS_EXPORT virtual void v_ImposeDirichletConditions(Array<OneD,NekDouble>& outarray);

            /// Scatters from the global coefficients
            /// \f$\boldsymbol{\hat{u}}_g\f$ to the local coefficients
            /// \f$\boldsymbol{\hat{u}}_l\f$.
            // inline
            MULTI_REGIONS_EXPORT virtual void v_GlobalToLocal(
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray);

            MULTI_REGIONS_EXPORT virtual void v_GlobalToLocal(void);

            /// Gathers the global coefficients \f$\boldsymbol{\hat{u}}_g\f$
            /// from the local coefficients \f$\boldsymbol{\hat{u}}_l\f$.
            // inline
            MULTI_REGIONS_EXPORT virtual void v_LocalToGlobal(
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray,
                bool useComm);

            MULTI_REGIONS_EXPORT virtual void v_LocalToGlobal(bool useComm);

            virtual void v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                    const StdRegions::ConstFactorMap &factors,
                    const StdRegions::VarCoeffMap &varcoeff,
                    const MultiRegions::VarFactorsMap &varfactors,
                    const Array<OneD, const NekDouble> &dirForcing,
                    const bool PhysSpaceForcing);

            virtual const Array<OneD,const SpatialDomains
                                ::BoundaryConditionShPtr>& v_GetBndConditions();

            virtual void v_BwdTrans(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,       NekDouble> &outarray);

            virtual void v_IProductWRTBase(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,       NekDouble> &outarray);

            /// Calculates the result of the multiplication of a global matrix
            /// of type specified by \a mkey with a vector given by \a inarray.
            virtual void v_GeneralMatrixOp(
                                const GlobalMatrixKey             &gkey,
                                const Array<OneD,const NekDouble> &inarray,
                                Array<OneD,      NekDouble> &outarray);

            virtual void v_ClearGlobalLinSysManager(void);

        };
        typedef std::shared_ptr<ContField1D>      ContField1DSharedPtr;

        // Inline implementations follow

        inline const Array<OneD,const MultiRegions::ExpListSharedPtr>&
                                ContField1D::GetBndCondExpansions()
        {
            return m_bndCondExpansions;
        }

        inline const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>&
                                ContField1D::GetBndConditions()
        {
            return m_bndConditions;
        }


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
         */
        inline void ContField1D::Assemble()
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
         * @param   inarray     An array of size \f$N_\mathrm{eof}\f$
         *                      containing the local degrees of freedom
         *                      \f$\boldsymbol{x}_l\f$.
         * @param   outarray    The resulting global degrees of freedom
         *                      \f$\boldsymbol{x}_g\f$ will be stored in this
         *                      array of size \f$N_\mathrm{dof}\f$.
         */
        inline void ContField1D::Assemble(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,NekDouble> &outarray)
        {
            m_locToGloMap->Assemble(inarray,outarray);
        }

        inline const AssemblyMapCGSharedPtr&
                                    ContField1D::GetLocalToGlobalMap() const
        {
            return  m_locToGloMap;
        }

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_CONTSOLNFIELD1D_H
