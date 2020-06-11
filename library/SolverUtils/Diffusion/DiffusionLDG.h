///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLDG.h
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
// Description: LDG diffusion class.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DIFFUSIONWEAKDG
#define NEKTAR_SOLVERUTILS_DIFFUSIONWEAKDG

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Diffusion/Diffusion.h>

namespace Nektar
{
    namespace SolverUtils
    {
        class DiffusionLDG : public Diffusion
        {
        public:
            static DiffusionSharedPtr create(std::string diffType)
            {
                boost::ignore_unused(diffType);
                return DiffusionSharedPtr(new DiffusionLDG());
            }

            static std::string type;

        protected:
            DiffusionLDG();

            std::string                            m_shockCaptureType;

            /// Coefficient of penalty term
            NekDouble                                         m_C11;

            Array<OneD, Array<OneD, NekDouble> >              m_traceNormals;
            LibUtilities::SessionReaderSharedPtr              m_session;

            virtual void v_InitObject(
                LibUtilities::SessionReaderSharedPtr               pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>        pFields);

            virtual void v_Diffuse(
                const std::size_t                                 nConvective,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray,
                const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >        &pBwd);

            virtual void v_DiffuseCoeffs(
                const std::size_t                                 nConvective,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray,
                const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >        &pBwd);

            virtual void v_DiffuseCalculateDerivative(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                TensorOfArray3D<NekDouble>                        &qfields,
                const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >        &pBwd);

            virtual void v_DiffuseVolumeFlux(
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble>>           &inarray,
                TensorOfArray3D<NekDouble>                          &qfields,
                TensorOfArray3D<NekDouble>                          &VolumeFlux,
                Array< OneD, int >                               &nonZeroIndex);

            virtual void v_DiffuseTraceFlux(
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble>>           &inarray,
                TensorOfArray3D<NekDouble>                          &qfields,
                TensorOfArray3D<NekDouble>                          &VolumeFlux,
                Array<OneD, Array<OneD, NekDouble> >                &TraceFlux,
                const Array<OneD, Array<OneD, NekDouble>>           &pFwd,
                const Array<OneD, Array<OneD, NekDouble>>           &pBwd,
                Array< OneD, int >                               &nonZeroIndex);

            void NumFluxforScalar(
                const Array<OneD, MultiRegions::ExpListSharedPtr>       &fields,
                const Array<OneD, Array<OneD, NekDouble> >              &ufield,
                TensorOfArray3D<NekDouble>                              &uflux,
                const Array<OneD, Array<OneD, NekDouble> >              &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >              &pBwd);

            void ApplyScalarBCs(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const std::size_t                                  var,
                const Array<OneD, const NekDouble>                &ufield,
                const Array<OneD, const NekDouble>                &Fwd,
                const Array<OneD, const NekDouble>                &Bwd,
                      Array<OneD,       NekDouble>                &penaltyflux);

            void NumFluxforVector(
                const Array<OneD, MultiRegions::ExpListSharedPtr>       &fields,
                const Array<OneD, Array<OneD, NekDouble> >              &ufield,
                TensorOfArray3D<NekDouble>                              &qfield,
                Array<OneD, Array<OneD, NekDouble> >                    &qflux);

            void ApplyVectorBCs(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const std::size_t                                  var,
                const std::size_t                                  dir,
                const Array<OneD, const NekDouble>                &qfield,
                const Array<OneD, const NekDouble>                &qFwd,
                const Array<OneD, const NekDouble>                &qBwd,
                      Array<OneD,       NekDouble>                &penaltyflux);
        };
    }
}

#endif
