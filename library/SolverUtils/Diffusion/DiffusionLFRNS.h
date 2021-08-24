///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLFRNS.h
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
// Description: LFRNS diffusion class.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DIFFUSIONLFRNS
#define NEKTAR_SOLVERUTILS_DIFFUSIONLFRNS

#include <SolverUtils/Diffusion/Diffusion.h>

namespace Nektar
{
    namespace SolverUtils
    {
        class DiffusionLFRNS : public Diffusion
        {
        public:
            static DiffusionSharedPtr create(std::string diffType)
            {
                return DiffusionSharedPtr(new DiffusionLFRNS(diffType));
            }

            static std::string                   type[];

            Array<OneD, NekDouble>               m_jac;
            Array<OneD, Array<OneD, NekDouble> > m_gmat;

            Array<OneD, Array<OneD, NekDouble> > m_Q2D_e0;
            Array<OneD, Array<OneD, NekDouble> > m_Q2D_e1;
            Array<OneD, Array<OneD, NekDouble> > m_Q2D_e2;
            Array<OneD, Array<OneD, NekDouble> > m_Q2D_e3;

            Array<OneD, Array<OneD, NekDouble> > m_dGL_xi1;
            Array<OneD, Array<OneD, NekDouble> > m_dGR_xi1;
            Array<OneD, Array<OneD, NekDouble> > m_dGL_xi2;
            Array<OneD, Array<OneD, NekDouble> > m_dGR_xi2;
            Array<OneD, Array<OneD, NekDouble> > m_dGL_xi3;
            Array<OneD, Array<OneD, NekDouble> > m_dGR_xi3;
            DNekMatSharedPtr                     m_Ixm;
            DNekMatSharedPtr                     m_Ixp;

        protected:
            DiffusionLFRNS(std::string diffType);

            Array<OneD, Array<OneD, NekDouble> > m_traceVel;
            Array<OneD, Array<OneD, NekDouble> > m_traceNormals;
            LibUtilities::SessionReaderSharedPtr m_session;
            NekDouble                            m_gamma;
            NekDouble                            m_gasConstant;
            NekDouble                            m_Twall;
            std::string                          m_ViscosityType;
            NekDouble                            m_mu;
            NekDouble                            m_thermalConductivity;
            NekDouble                            m_rhoInf;
            NekDouble                            m_pInf;


            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_IF1;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_DU1;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_DFC1;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_BD1;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_D1;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_DD1;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_viscTensor;
            Array<OneD, Array<OneD, NekDouble> >               m_viscFlux;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_DFC2;
            Array<OneD, Array<OneD, NekDouble> >               m_divFD;
            Array<OneD, Array<OneD, NekDouble> >               m_divFC;

            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_tmp1;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_tmp2;

            Array<OneD, Array<OneD, NekDouble> > m_homoDerivs;

            int                                  m_spaceDim;
            int                                  m_diffDim;

            std::string m_diffType;

            virtual void v_InitObject(
                LibUtilities::SessionReaderSharedPtr               pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>        pFields);

            virtual void v_SetupMetrics(
                LibUtilities::SessionReaderSharedPtr               pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>        pFields);

            virtual void v_SetupCFunctions(
                LibUtilities::SessionReaderSharedPtr               pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>        pFields);

            virtual void v_Diffuse(
                const std::size_t                                  nConvective,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >        &pBwd);

            virtual void v_NumericalFluxO1(
                const Array<OneD, MultiRegions::ExpListSharedPtr>      &fields,
                const Array<OneD, Array<OneD, NekDouble> >             &inarray,
                      Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                                            &numericalFluxO1);

            virtual void v_WeakPenaltyO1(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >      &penaltyfluxO1);

            virtual void v_NumericalFluxO2(
                const Array<OneD, MultiRegions::ExpListSharedPtr>       &fields,
                const Array<OneD, Array<OneD, NekDouble> >              &ufield,
                      Array<OneD, Array<OneD, Array<OneD, NekDouble> > >&qfield,
                      Array<OneD, Array<OneD, NekDouble> >              &qflux);

            virtual void v_WeakPenaltyO2(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const int                                          var,
                const int                                          dir,
                const Array<OneD, const NekDouble>                &qfield,
                      Array<OneD,       NekDouble>                &penaltyflux);

            virtual void v_DerCFlux_1D(
                const int                                     nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>&fields,
                const Array<OneD, const NekDouble>               &flux,
                const Array<OneD, const NekDouble>               &iFlux,
                      Array<OneD,       NekDouble>               &derCFlux);

            virtual void v_DerCFlux_2D(
                const int                                     nConvectiveFields,
                const int                                        direction,
                const Array<OneD, MultiRegions::ExpListSharedPtr>&fields,
                const Array<OneD, const NekDouble>               &flux,
                const Array<OneD,       NekDouble>               &iFlux,
                      Array<OneD,       NekDouble>               &derCFlux);

            virtual void v_DivCFlux_2D(
                const int                                      nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>&fields,
                const Array<OneD, const NekDouble>               &fluxX1,
                const Array<OneD, const NekDouble>               &fluxX2,
                const Array<OneD, const NekDouble>               &numericalFlux,
                      Array<OneD,       NekDouble>               &divCFlux);

            virtual void v_DivCFlux_2D_Gauss(
                const int                                     nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>&fields,
                const Array<OneD, const NekDouble>               &fluxX1,
                const Array<OneD, const NekDouble>               &fluxX2,
                const Array<OneD, const NekDouble>               &numericalFlux,
                      Array<OneD,       NekDouble>               &divCFlux);

            virtual void v_FluxVec(
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &fluxvector)
            {
                fluxvector = m_viscTensor;
            };

            virtual void v_SetHomoDerivs(
                                         Array<OneD, Array<OneD, NekDouble> > &deriv)
            {
                m_homoDerivs = deriv;
            }

            virtual Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &v_GetFluxTensor()
            {
                return m_viscTensor;
            }
        };

        typedef std::shared_ptr<DiffusionLFRNS> DiffusionLFRNSSharedPtr;
    }
}

#endif

