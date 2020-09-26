///////////////////////////////////////////////////////////////////////////////
//
// File: AdvectionFR.h
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
// Description: FR advection class.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_ADVECTIONFR
#define NEKTAR_SOLVERUTILS_ADVECTIONFR

#include <SolverUtils/Advection/Advection.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
    namespace SolverUtils
    {
        class AdvectionFR : public Advection
        {
        public:
            static AdvectionSharedPtr create(std::string advType)
            {
                return AdvectionSharedPtr(new AdvectionFR(advType));
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
            AdvectionFR(std::string advType);

            Array<OneD, Array<OneD, NekDouble> >               m_traceNormals;

            std::string m_advType;

            virtual void v_InitObject(
                LibUtilities::SessionReaderSharedPtr              pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>       pFields);

            virtual void v_Advect(
                const int nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &advVel,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const NekDouble                                   &time,
                const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayofArray);

            virtual void v_SetupMetrics(
                LibUtilities::SessionReaderSharedPtr              pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>       pFields);

            virtual void v_SetupCFunctions(
                LibUtilities::SessionReaderSharedPtr              pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>       pFields);

            virtual void v_DivCFlux_1D(
                const int nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, const NekDouble> &fluxX1,
                const Array<OneD, const NekDouble> &numericalFlux,
                      Array<OneD,       NekDouble> &divCFlux);

            virtual void v_DivCFlux_2D(
                const int nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, const NekDouble> &fluxX1,
                const Array<OneD, const NekDouble> &fluxX2,
                const Array<OneD, const NekDouble> &numericalFlux,
                      Array<OneD,       NekDouble> &divCFlux);

            virtual void v_DivCFlux_2D_Gauss(
                const int nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, const NekDouble> &fluxX1,
                const Array<OneD, const NekDouble> &fluxX2,
                const Array<OneD, const NekDouble> &numericalFlux,
                Array<OneD,       NekDouble> &divCFlux);

            virtual void v_DivCFlux_3D(
                const int nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, const NekDouble> &fluxX1,
                const Array<OneD, const NekDouble> &fluxX2,
                const Array<OneD, const NekDouble> &fluxX3,
                const Array<OneD, const NekDouble> &numericalFlux,
                      Array<OneD,       NekDouble> &divCFlux);
        };
    }
}

#endif
