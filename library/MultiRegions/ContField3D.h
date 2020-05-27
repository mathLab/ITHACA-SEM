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
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/DisContField3D.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/GlobalMatrix.h>
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

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
                        const std::string &variable  = "DefaultVar",
                        const bool CheckIfSingularSystem = false,
                        const Collections::ImplementationType ImpType
                        = Collections::eNoImpType);

            /// Construct a global continuous field with solution type based on
            /// another field but using a separate input mesh and boundary
            /// conditions.
            MULTI_REGIONS_EXPORT ContField3D(const ContField3D &In,
                        const SpatialDomains::MeshGraphSharedPtr &graph3D,
                        const std::string &variable,
                        const bool CheckIfSingularSystem = false);

            MULTI_REGIONS_EXPORT ContField3D(const ContField3D &In);

            MULTI_REGIONS_EXPORT virtual ~ContField3D();

            inline const Array<OneD,const MultiRegions::ExpListSharedPtr>& GetBndCondExpansions()
            {
                return m_bndCondExpansions;
            }

            /// This function return the boundary conditions expansion.
            inline const Array<OneD,const MultiRegions::ExpListSharedPtr>
                    &GetBndCondExp();

            inline void Assemble();

            inline void Assemble(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,NekDouble> &outarray);

            inline const AssemblyMapCGSharedPtr& GetLocalToGlobalMap()
                                                                        const;

            MULTI_REGIONS_EXPORT int GetGlobalMatrixNnz(const GlobalMatrixKey &gkey);


        protected:
            AssemblyMapCGSharedPtr m_locToGloMap;

            /// (A shared pointer to) a list which collects all the global
            /// matrices being assembled, such that they should be constructed
            /// only once.
            GlobalMatrixMapShPtr            m_globalMat;

            /// (A shared pointer to) a list which collects all the global
            /// linear system being assembled, such that they should be
            /// constructed only once.
            LibUtilities::NekManager<GlobalLinSysKey, GlobalLinSys> m_globalLinSysManager;

            /// Performs the backward transformation of the spectral/hp
            /// element expansion.
            virtual void v_BwdTrans(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD,       NekDouble> &outarray);

            /// Calculates the inner product of a function
            /// \f$f(\boldsymbol{x})\f$ with respect to all <em>global</em>
            /// expansion modes \f$\phi_n^e(\boldsymbol{x})\f$.
            virtual void v_IProductWRTBase(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD,       NekDouble> &outarray);

            virtual void v_FwdTrans(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD,       NekDouble> &outarray);

            
        private:
            GlobalLinSysSharedPtr GetGlobalLinSys(const GlobalLinSysKey &mkey);

            GlobalLinSysSharedPtr GenGlobalLinSys(const GlobalLinSysKey &mkey);

            /// Returns the global matrix specified by \a mkey.
            GlobalMatrixSharedPtr GetGlobalMatrix(const GlobalMatrixKey &mkey);


            void GlobalSolve(const GlobalLinSysKey &key,
                    const Array<OneD, const NekDouble> &rhs,
                          Array<OneD,       NekDouble> &inout,
                    const Array<OneD, const NekDouble> &dirForcing
                                                     = NullNekDouble1DArray);

            /// Impose the Dirichlet Boundary Conditions on outarray 
            virtual void v_ImposeDirichletConditions(Array<OneD,NekDouble>& outarray);

            virtual void v_FillBndCondFromField();

            virtual void v_FillBndCondFromField(const int nreg);

            virtual void v_LocalToGlobal(bool useComm);


            virtual void v_LocalToGlobal(
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray,
                bool useComm);


            virtual void v_GlobalToLocal(void);


            virtual void v_GlobalToLocal(
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray);


            virtual void v_MultiplyByInvMassMatrix(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD,       NekDouble> &outarray);

            virtual void v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                    const StdRegions::ConstFactorMap &factors,
                    const StdRegions::VarCoeffMap &varcoeff,
                    const MultiRegions::VarFactorsMap &varfactors,
                    const Array<OneD, const NekDouble> &dirForcing,
                    const bool PhysSpaceForcing);
            virtual void v_GeneralMatrixOp(
                    const GlobalMatrixKey             &gkey,
                    const Array<OneD,const NekDouble> &inarray,
                    Array<OneD,      NekDouble> &outarray);

            // Solve the linear advection problem assuming that m_coeffs
            // vector contains an intial estimate for solution
            MULTI_REGIONS_EXPORT virtual void v_LinearAdvectionDiffusionReactionSolve(
                    const Array<OneD, Array<OneD, NekDouble> > &velocity,
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray,
                    const NekDouble lambda,
                    const Array<OneD, const NekDouble> &dirForcing = NullNekDouble1DArray);
            
            virtual void v_ClearGlobalLinSysManager(void);

        };
        typedef std::shared_ptr<ContField3D>      ContField3DSharedPtr;

        inline const Array<OneD,const MultiRegions::ExpListSharedPtr>
                &ContField3D::GetBndCondExp()
        {
            return m_bndCondExpansions;
        }


        inline void ContField3D::Assemble()
        {
            m_locToGloMap->Assemble(m_coeffs, m_coeffs);
        }

        inline void ContField3D::Assemble(
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray)
        {
            m_locToGloMap->Assemble(inarray, outarray);
        }

        inline const AssemblyMapCGSharedPtr&
                ContField3D::GetLocalToGlobalMap() const
        {
            return  m_locToGloMap;
        }


    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_CONTFIELD3D_H
