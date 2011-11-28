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
// Description: Field definition in one-dimension
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_CONTFIELD1D_H
#define NEKTAR_LIBS_MULTIREGIONS_CONTFIELD1D_H

#include <LibUtilities/Communication/Comm.h>

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/DisContField1D.h>
#include <MultiRegions/LocalToGlobalC0ContMap.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList0D.h>
#include <LocalRegions/PointExp.h>
#include <SpatialDomains/MeshGraph1D.h>
#include <SpatialDomains/Conditions.h>
#include <MultiRegions/MultiRegionsDeclspec.h>

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
                        const std::string &variable);

            /// Copy constructor.
            MULTI_REGIONS_EXPORT ContField1D(const ContField1D &In);

            /// Copy constructor.
            MULTI_REGIONS_EXPORT ContField1D( const LibUtilities::SessionReaderSharedPtr &pSession,
                                              const ExpList1D & In);
            /// Destructor
            MULTI_REGIONS_EXPORT ~ContField1D();

            /// Perform global forward transformation of a function \f$f(x)\f$,
            //  subject to the boundary conditions specified.
            MULTI_REGIONS_EXPORT void FwdTrans(      const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,      NekDouble> &outarray,
                                bool  UseContCoeffs = false);

            /// This function performs the backward transformation of the
            /// spectral/hp element expansion.
            MULTI_REGIONS_EXPORT void BwdTrans(      const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                bool  UseContCoeffs = false);

            ///
            MULTI_REGIONS_EXPORT void MultiplyByInvMassMatrix(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                bool  UseContCoeffs = false);

            /// Return the boundary conditions expansion.
            // inline
            MULTI_REGIONS_EXPORT const Array<OneD,const MultiRegions::ExpListSharedPtr>&
                                                     GetBndCondExpansions();

            // inline
            MULTI_REGIONS_EXPORT const Array<OneD,const SpatialDomains
                                ::BoundaryConditionShPtr>& GetBndConditions();

            /// Returns the total number of global degrees of freedom
            /// \f$N_{\mathrm{dof}}\f$.
            // inline
            inline int GetContNcoeffs();

            /// Returns (a reference to) the array \f$\boldsymbol{\hat{u}}_g\f$
            /// (implemented as #m_contCoeffs) containing all global expansion
            /// coefficients.
            // inline
            MULTI_REGIONS_EXPORT Array<OneD, NekDouble> &UpdateContCoeffs();

            /// Returns (a reference to) the array \f$\boldsymbol{\hat{u}}_g\f$
            /// (implemented as #m_contCoeffs) containing all global expansion
            /// coefficients.
            // inline
            MULTI_REGIONS_EXPORT const Array<OneD, const NekDouble> &GetContCoeffs() const;

            /// Scatters from the global coefficients
            /// \f$\boldsymbol{\hat{u}}_g\f$ to the local coefficients
            /// \f$\boldsymbol{\hat{u}}_l\f$.
            // inline
            MULTI_REGIONS_EXPORT void GlobalToLocal();

            /// Scatters from the global coefficients
            /// \f$\boldsymbol{\hat{u}}_g\f$ to the local coefficients
            /// \f$\boldsymbol{\hat{u}}_l\f$.
            // inline
            MULTI_REGIONS_EXPORT void GlobalToLocal( const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,NekDouble> &outarray);

            /// Gathers the global coefficients \f$\boldsymbol{\hat{u}}_g\f$
            /// from the local coefficients \f$\boldsymbol{\hat{u}}_l\f$.
            // inline
            MULTI_REGIONS_EXPORT void LocalToGlobal();

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
            MULTI_REGIONS_EXPORT const LocalToGlobalC0ContMapSharedPtr& GetLocalToGlobalMap() const;

            /// Calculates the inner product of a function \f$f(x)\f$ with
            /// respect to all <em>global</em> expansion modes
            /// \f$\phi_n^e(x)\f$.
            MULTI_REGIONS_EXPORT void IProductWRTBase(const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD, NekDouble> &outarray,
                                bool  UseContCoeffs = false);

            /// Calculates the result of the multiplication of a global matrix
            /// of type specified by \a mkey with a vector given by \a inarray.
            MULTI_REGIONS_EXPORT void GeneralMatrixOp(const GlobalMatrixKey             &gkey,
                                const Array<OneD,const NekDouble> &inarray,
                                      Array<OneD,      NekDouble> &outarray,
                                bool  UseContCoeffs = false);

        protected:
            /// (A shared pointer to) the object which contains all the required
            /// information for the transformation from local to global degrees
            /// of freedom.
            LocalToGlobalC0ContMapSharedPtr m_locToGloMap;

            /// The total number of global degrees of freedom.
            /// #m_contNcoeffs\f$=N_{\mathrm{dof}}\f$
            int                             m_contNcoeffs;

            /// The array of length #m_ncoeffs\f$=N_{\mathrm{dof}}\f$ containing
            /// the global expansion coefficients.
            Array<OneD, NekDouble>          m_contCoeffs;

            /// (A shared pointer to) a list which collects all the global
            /// matrices being assembled, such that they should be constructed
            /// only once.
            GlobalMatrixMapShPtr            m_globalMat;

            /// (A shared pointer to) a list which collects all the global
            /// matrices being assembled, such that they should be constructed
            /// only once.
            GlobalLinSysMapShPtr            m_globalLinSys;

        private:
			
			MULTI_REGIONS_EXPORT virtual int v_GetContNcoeffs() const;
			
			MULTI_REGIONS_EXPORT virtual void v_SetContCoeffsArray(Array<OneD, NekDouble> &inarray);
			
            /// Returns the linear system specified by \a mkey.
            GlobalLinSysSharedPtr GetGlobalLinSys(const GlobalLinSysKey &mkey);

            /// Solve the linear system specified by the key \a key.
            void GlobalSolve(   const GlobalLinSysKey &key,
                                const Array<OneD, const NekDouble> &rhs,
                                      Array<OneD,       NekDouble> &inout,
                                const Array<OneD, const NekDouble> &dirForcing
                                                        = NullNekDouble1DArray);

            /// Perform a forward transform
            virtual void v_FwdTrans(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                bool  UseContCoeffs);

            virtual void v_MultiplyByInvMassMatrix(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                bool  UseContCoeffs);

            virtual void v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff);

            virtual void v_HelmSolveCG(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff,
                          bool UseContCoeffs,
                    const Array<OneD, const NekDouble> &dirForcing);

            virtual const Array<OneD,const SpatialDomains
                                ::BoundaryConditionShPtr>& v_GetBndConditions();

            virtual const Array<OneD, const NekDouble> &v_GetContCoeffs() const;

            virtual void v_BwdTrans(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                bool  UseContCoeffs);

            virtual void v_IProductWRTBase(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                bool  UseContCoeffs);

            /// Calculates the result of the multiplication of a global matrix
            /// of type specified by \a mkey with a vector given by \a inarray.
            virtual void v_GeneralMatrixOp(
                                const GlobalMatrixKey             &gkey,
                                const Array<OneD,const NekDouble> &inarray,
                                      Array<OneD,      NekDouble> &outarray,
                                bool  UseContCoeffs);

        };
        typedef boost::shared_ptr<ContField1D>      ContField1DSharedPtr;

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
         * @return  #m_contNcoeffs, the total number of global degrees of
         * freedom.
         */
        inline int ContField1D::GetContNcoeffs()
        {
            return m_contNcoeffs;
        }

        /**
         * If one wants to get hold of the underlying data without modifying
         * them, rather use the function #GetContCoeffs instead.
         *
         * @return (A reference to) the array #m_contCoeffs.
         */
        inline Array<OneD, NekDouble> &ContField1D::UpdateContCoeffs()
        {
            m_transState = eContinuous;
            return m_contCoeffs;
        }

        /**
         * As the function returns a constant reference to a
         * <em>const Array</em>, it is not possible to modify the underlying
         * data of the array #m_contCoeffs. In order to do so, use the function
         * #UpdateContCoeffs instead.
         *
         * \return (A reference to) the array #m_contCoeffs.
         */
        inline const Array<OneD, const NekDouble>&
                                ContField1D::GetContCoeffs() const
        {
            return m_contCoeffs;
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
         * @note The array #m_contCoeffs should be filled with the global
         * coefficients \f$\boldsymbol{\hat{u}}_g\f$ and that the resulting
         * local coefficients \f$\boldsymbol{\hat{u}}_l\f$ will be stored in
         * #m_coeffs.
         */
        inline void ContField1D::GlobalToLocal()
        {
            m_locToGloMap->GlobalToLocal(m_contCoeffs,m_coeffs);
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
         * where \a map\f$[e][i]\f$ is the mapping array and \a
         * sign\f$[e][i]\f$ is an array of similar dimensions ensuring the
         * correct modal connectivity between the different elements (both
         * these arrays are contained in the data member #m_locToGloMap). This
         * operation is equivalent to the scatter operation
         * \f$\boldsymbol{\hat{u}}_l=\mathcal{A}\boldsymbol{\hat{u}}_g\f$, where
         * \f$\mathcal{A}\f$ is the
         * \f$N_{\mathrm{eof}}\times N_{\mathrm{dof}}\f$ permutation matrix.
         *
         * @param   inarray     An array of size \f$N_\mathrm{dof}\f$
         *                      containing the global degrees of freedom
         *                      \f$\boldsymbol{x}_g\f$.
         * @param   outarray    The resulting local degrees of freedom
         *                      \f$\boldsymbol{x}_l\f$ will be stored in this
         *                      array of size \f$N_\mathrm{eof}\f$.
         */
        inline void ContField1D::GlobalToLocal(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,NekDouble> &outarray)
        {
            m_locToGloMap->GlobalToLocal(inarray,outarray);
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
         * @note    The array #m_coeffs should be filled with the local
         *          coefficients \f$\boldsymbol{\hat{u}}_l\f$ and that the
         *          resulting global coefficients \f$\boldsymbol{\hat{u}}_g\f$
         *          will be stored in #m_contCoeffs.
         */
        inline void ContField1D::LocalToGlobal()
        {
            m_locToGloMap->LocalToGlobal(m_coeffs,m_contCoeffs);
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
         * @note    The array #m_coeffs should be filled with the local
         * coefficients \f$\boldsymbol{\hat{u}}_l\f$ and that the resulting
         * global coefficients \f$\boldsymbol{\hat{u}}_g\f$ will be stored in
         * #m_contCoeffs.
         */
        inline void ContField1D::Assemble()
        {
            m_locToGloMap->Assemble(m_coeffs,m_contCoeffs);
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

        inline const LocalToGlobalC0ContMapSharedPtr&
                                    ContField1D::GetLocalToGlobalMap() const
        {
            return  m_locToGloMap;
        }

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_CONTSOLNFIELD1D_H
