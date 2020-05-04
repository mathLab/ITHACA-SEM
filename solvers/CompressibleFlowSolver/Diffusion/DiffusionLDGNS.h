///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLDGNS.h
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

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_DIFFUSIONLDGNS
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_DIFFUSIONLDGNS

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Diffusion/Diffusion.h>
#include <CompressibleFlowSolver/Misc/EquationOfState.h>
#include <LocalRegions/Expansion3D.h>
#include <LocalRegions/Expansion2D.h>


using namespace Nektar::SolverUtils;

namespace Nektar
{
    class DiffusionLDGNS : public Diffusion
    {
    public:
        static DiffusionSharedPtr create(std::string diffType)
        {
            boost::ignore_unused(diffType);
            return DiffusionSharedPtr(new DiffusionLDGNS());
        }

        static std::string type;

    protected:
        DiffusionLDGNS();

        /// Penalty coefficient for LDGNS
        NekDouble                            m_C11;

        /// h scaling for penalty term
        Array<OneD, NekDouble>               m_traceOneOverH;

        Array<OneD, Array<OneD, NekDouble> > m_traceVel;
        Array<OneD, Array<OneD, NekDouble> > m_traceNormals;
        LibUtilities::SessionReaderSharedPtr m_session;
        NekDouble                            m_Twall;
        /// Equation of system for computing temperature
        EquationOfStateSharedPtr             m_eos;

        TensorOfArray3D<NekDouble> m_viscTensor;

        Array<OneD, Array<OneD, NekDouble> > m_homoDerivs;

        std::size_t                                  m_spaceDim;
        std::size_t                                  m_diffDim;

        virtual void v_InitObject(
            LibUtilities::SessionReaderSharedPtr               pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr>        pFields);

        virtual void v_Diffuse(
            const std::size_t                                  nConvective,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> > &pFwd,
            const Array<OneD, Array<OneD, NekDouble> > &pBwd);
        virtual void v_DiffuseCoeffs(
            const std::size_t                                 nConvective,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd);

        virtual void v_DiffuseCalculateDerivative(
            const Array<OneD, MultiRegions::ExpListSharedPtr>&fields,
            const Array<OneD, Array<OneD, NekDouble> >       &inarray,
            TensorOfArray3D<NekDouble>                       &qfields,
            const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >       &pBwd);

        virtual void v_DiffuseVolumeFlux(
            const Array<OneD, MultiRegions::ExpListSharedPtr>&fields,
            const Array<OneD, Array<OneD, NekDouble>>        &inarray,
            TensorOfArray3D<NekDouble>                       &qfields,
            TensorOfArray3D<NekDouble>                       &VolumeFlux,
            Array< OneD, int >                               &nonZeroIndex);
        virtual void v_DiffuseTraceFlux(
            const Array<OneD, MultiRegions::ExpListSharedPtr>&fields,
            const Array<OneD, Array<OneD, NekDouble>>        &inarray,
            TensorOfArray3D<NekDouble>                       &qfields,
            TensorOfArray3D<NekDouble>                       &VolumeFlux,
            Array<OneD, Array<OneD, NekDouble> >             &TraceFlux,
            const Array<OneD, Array<OneD, NekDouble>>        &pFwd,
            const Array<OneD, Array<OneD, NekDouble>>        &pBwd,
            Array< OneD, int >                               &nonZeroIndex);

        void NumericalFluxO1(
            const Array<OneD, MultiRegions::ExpListSharedPtr>  &fields,
            const Array<OneD, Array<OneD, NekDouble> >         &inarray,
            TensorOfArray3D<NekDouble>                         &numericalFluxO1,
            const Array<OneD, Array<OneD, NekDouble> >         &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >         &pBwd);

        void ApplyBCsO1(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd,
                  Array<OneD, Array<OneD, NekDouble> >        &flux01);

        void NumericalFluxO2(
            const Array<OneD, MultiRegions::ExpListSharedPtr>       &fields,
            TensorOfArray3D<NekDouble>                              &qfield,
            Array<OneD, Array<OneD, NekDouble> >                    &qflux,
            const Array<OneD, Array<OneD, NekDouble> >              &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >              &pBwd);

        void ApplyBCsO2(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const std::size_t                                  var,
            const std::size_t                                  dir,
            const Array<OneD, const NekDouble>                &qfield,
            const Array<OneD, const NekDouble>                &qFwd,
            const Array<OneD, const NekDouble>                &qBwd,
                  Array<OneD,       NekDouble>                &penaltyflux);

        virtual void v_SetHomoDerivs(
            Array<OneD, Array<OneD, NekDouble> > &deriv)
        {
            m_homoDerivs = deriv;
        }

        virtual TensorOfArray3D<NekDouble> &v_GetFluxTensor()
        {
            return m_viscTensor;
        }
    };

    typedef std::shared_ptr<DiffusionLDGNS> DiffusionLDGNSSharedPtr;
}

#endif

