///////////////////////////////////////////////////////////////////////////////
//
// File ContField3D.h
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
// Description: Field definition in three-dimensions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_CONTFIELD3D_H
#define NEKTAR_LIBS_MULTIREGIONS_CONTFIELD3D_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <LibUtilities/Communication/Comm.h>

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/DisContField3D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/LocalToGlobalC0ContMap.h>
#include <MultiRegions/GlobalLinSys.h>

#include <SpatialDomains/MeshGraph3D.h>
#include <SpatialDomains/Conditions.h>

namespace Nektar
{
    namespace MultiRegions
    {
        class ContField3D: public DisContField3D
        {
        public:
            MULTI_REGIONS_EXPORT ContField3D();

            /// Construct a global continuous field.
            MULTI_REGIONS_EXPORT ContField3D(
                        const LibUtilities::SessionReaderSharedPtr &pSession,
                        const SpatialDomains::MeshGraphSharedPtr &graph3D,
                        const std::string &variable);

            /// Construct a global continuous field with solution type based on
            /// another field but using a separate input mesh and boundary
            /// conditions.
            MULTI_REGIONS_EXPORT ContField3D(const ContField3D &In,
                        const SpatialDomains::MeshGraphSharedPtr &graph3D,
                        const std::string &variable);

            MULTI_REGIONS_EXPORT ContField3D(const ContField3D &In);

            MULTI_REGIONS_EXPORT ~ContField3D();

            inline const Array<OneD,const MultiRegions::ExpListSharedPtr>& GetBndCondExpansions()
            {
                return m_bndCondExpansions;
            }

            /// This function return the boundary conditions expansion.
            inline const Array<OneD,const MultiRegions::ExpListSharedPtr>
                    &GetBndCondExp();

            MULTI_REGIONS_EXPORT void GenerateDirBndCondForcing(
                    const GlobalLinSysKey &key,
                    Array<OneD, NekDouble> &inout,
                    Array<OneD, NekDouble> &outarray);

            /// Returns (a reference to) the array \f$\boldsymbol{\hat{u}}_g\f$
            /// (implemented as #m_contCoeffs) containing all global expansion
            /// coefficients.
            inline Array<OneD, NekDouble> &UpdateContCoeffs();

            /// Returns (a reference to) the array \f$\boldsymbol{\hat{u}}_g\f$
            /// (implemented as #m_contCoeffs) containing all global expansion
            /// coefficients.
            inline const Array<OneD, const NekDouble> &GetContCoeffs() const;

            /// Returns the total number of global degrees of freedom
            /// \f$N_{\mathrm{dof}}\f$.
            inline int GetContNcoeffs();

            inline void GlobalToLocal();

            inline void GlobalToLocal(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,NekDouble> &outarray);

            inline void LocalToGlobal();

            inline void Assemble();

            inline void Assemble(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,NekDouble> &outarray);

            inline const LocalToGlobalC0ContMapSharedPtr& GetLocalToGlobalMap()
                                                                        const;

            MULTI_REGIONS_EXPORT int GetGlobalMatrixNnz(const GlobalMatrixKey &gkey);


        protected:
            LocalToGlobalC0ContMapSharedPtr m_locToGloMap;
            int                             m_contNcoeffs;
            Array<OneD, NekDouble>          m_contCoeffs;

            /// (A shared pointer to) a list which collects all the global
            /// matrices being assembled, such that they should be constructed
            /// only once.
            GlobalMatrixMapShPtr            m_globalMat;

            /// (A shared pointer to) a list which collects all the global
            /// linear system being assembled, such that they should be
            /// constructed only once.
            GlobalLinSysMapShPtr            m_globalLinSys;

            /// Performs the backward transformation of the spectral/hp
            /// element expansion.
            virtual void v_BwdTrans(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          bool UseContCoeffs = false);

            /// Calculates the inner product of a function
            /// \f$f(\boldsymbol{x})\f$ with respect to all <em>global</em>
            /// expansion modes \f$\phi_n^e(\boldsymbol{x})\f$.
            virtual void v_IProductWRTBase(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          bool  UseContCoeffs = false);

            virtual void v_FwdTrans(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          bool  UseContCoeffs);

        private:
            GlobalLinSysSharedPtr GetGlobalLinSys(const GlobalLinSysKey &mkey);

            /// Returns the global matrix specified by \a mkey.
            GlobalMatrixSharedPtr GetGlobalMatrix(const GlobalMatrixKey &mkey);


            void GlobalSolve(const GlobalLinSysKey &key,
                    const Array<OneD, const NekDouble> &rhs,
                          Array<OneD,       NekDouble> &inout,
                    const Array<OneD, const NekDouble> &dirForcing
                                                     = NullNekDouble1DArray);

            virtual void v_MultiplyByInvMassMatrix(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          bool  UseContCoeffs);

            virtual void v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                    const FlagList &flags,
                    const StdRegions::ConstFactorMap &factors,
                    const StdRegions::VarCoeffMap &varcoeff,
                    const Array<OneD, const NekDouble> &dirForcing);

            virtual void v_GeneralMatrixOp(
                    const GlobalMatrixKey             &gkey,
                    const Array<OneD,const NekDouble> &inarray,
                          Array<OneD,      NekDouble> &outarray,
                          bool  UseContCoeffs);


        };
        typedef boost::shared_ptr<ContField3D>      ContField3DSharedPtr;

        inline const Array<OneD,const MultiRegions::ExpListSharedPtr>
                &ContField3D::GetBndCondExp()
        {
            return m_bndCondExpansions;
        }

        /**
         * If one wants to get hold of the underlying data without modifying
         * them, rather use the function #GetContCoeffs instead.
         *
         * @return  (A reference to) the array #m_contCoeffs.
         */
        inline Array<OneD, NekDouble> &ContField3D::UpdateContCoeffs()
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
        inline const Array<OneD, const NekDouble>
                                            &ContField3D::GetContCoeffs() const
        {
            return m_contCoeffs;
        }

        /**
         * @return  #m_contNcoeffs, the total number of global degrees of
         * freedom.
         */
        inline int ContField3D::GetContNcoeffs()
        {
            return m_contNcoeffs;
        }

        inline void ContField3D::GlobalToLocal()
        {
            m_locToGloMap->GlobalToLocal(m_contCoeffs, m_coeffs);
        }

        inline void ContField3D::GlobalToLocal(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray)
        {
            m_locToGloMap->GlobalToLocal(inarray, outarray);
        }

        inline void ContField3D::LocalToGlobal()
        {
            m_locToGloMap->LocalToGlobal(m_coeffs, m_contCoeffs);
        }

        inline void ContField3D::Assemble()
        {
            m_locToGloMap->Assemble(m_coeffs, m_contCoeffs);
        }

        inline void ContField3D::Assemble(
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray)
        {
            m_locToGloMap->Assemble(inarray, outarray);
        }

        inline const LocalToGlobalC0ContMapSharedPtr&
                ContField3D::GetLocalToGlobalMap() const
        {
            return  m_locToGloMap;
        }

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_CONTFIELD3D_H
