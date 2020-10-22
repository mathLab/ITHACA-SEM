///////////////////////////////////////////////////////////////////////////////
//
// File PulseWavePressureArea.h
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
// Description: PulseWavePressureArea header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_PULSEWAVEPRESSUREAREA_H
#define NEKTAR_PULSEWAVEPRESSUREAREA_H

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>

namespace Nektar
{

class PulseWavePressureArea;

typedef std::shared_ptr<PulseWavePressureArea> PulseWavePressureAreaSharedPtr;

static PulseWavePressureAreaSharedPtr NullPulseWavePressureAreaSharedPtr;

typedef LibUtilities::NekFactory<std::string, PulseWavePressureArea,
                                  Array<OneD, MultiRegions::ExpListSharedPtr> &,
                                  const LibUtilities::SessionReaderSharedPtr &>
    PressureAreaFactory;
PressureAreaFactory &GetPressureAreaFactory();

class PulseWavePressureArea
{
    public:
        PulseWavePressureArea(Array<OneD, MultiRegions::ExpListSharedPtr>
                &pVessel, const LibUtilities::SessionReaderSharedPtr &pSession);

        virtual ~PulseWavePressureArea();

        inline void GetPressure(NekDouble &P, const NekDouble &beta,
                const NekDouble &A, const NekDouble &A0, const NekDouble &dAUdx,
                      const NekDouble &gamma = 0, const NekDouble &alpha = 0.5);

        inline void GetC(NekDouble &c, const NekDouble &beta, const NekDouble &A,
                             const NekDouble &A0, const NekDouble &alpha = 0.5);

        inline void GetW1(NekDouble &W1, const NekDouble &u,
                 const NekDouble &beta, const NekDouble &A, const NekDouble &A0,
                                                  const NekDouble &alpha = 0.5);

        inline void GetW2(NekDouble &W2, const NekDouble &u,
                 const NekDouble &beta, const NekDouble &A, const NekDouble &A0,
                                                  const NekDouble &alpha = 0.5);

        inline void GetAFromChars(NekDouble &A, const NekDouble &W1,
                const NekDouble &W2, const NekDouble &beta, const NekDouble &A0,
                                                  const NekDouble &alpha = 0.5);

        inline void GetUFromChars(NekDouble &u, const NekDouble &W1,
                                                           const NekDouble &W2);

        inline void GetCharIntegral(NekDouble &I, const NekDouble &beta,
         const NekDouble &A, const NekDouble &A0, const NekDouble &alpha = 0.5);

        inline void GetJacobianInverse(NekMatrix<NekDouble> &invJ,
             const Array<OneD, NekDouble> &Au, const Array<OneD, NekDouble> &uu,
           const Array<OneD, NekDouble> &beta, const Array<OneD, NekDouble> &A0,
                  const Array<OneD, NekDouble> &alpha, const std::string &type);

    protected:
        virtual void v_GetPressure(NekDouble &P, const NekDouble &beta,
                const NekDouble &A, const NekDouble &A0, const NekDouble &dAUdx,
                  const NekDouble &gamma = 0, const NekDouble &alpha = 0.5) = 0;

        virtual void v_GetC(NekDouble &c, const NekDouble &beta,
                                        const NekDouble &A, const NekDouble &A0,
                                              const NekDouble &alpha = 0.5) = 0;

        virtual void v_GetW1(NekDouble &W1, const NekDouble &u,
                 const NekDouble &beta, const NekDouble &A, const NekDouble &A0,
                                              const NekDouble &alpha = 0.5) = 0;

        virtual void v_GetW2(NekDouble &W2, const NekDouble &u,
                 const NekDouble &beta, const NekDouble &A, const NekDouble &A0,
                                              const NekDouble &alpha = 0.5) = 0;

        virtual void v_GetAFromChars(NekDouble &A, const NekDouble &W1,
                const NekDouble &W2, const NekDouble &beta, const NekDouble &A0,
                                              const NekDouble &alpha = 0.5) = 0;

        virtual void v_GetUFromChars(NekDouble &u, const NekDouble &W1,
                                                       const NekDouble &W2) = 0;

        virtual void v_GetCharIntegral(NekDouble &I, const NekDouble &beta,
                                        const NekDouble &A, const NekDouble &A0,
                                              const NekDouble &alpha = 0.5) = 0;

        virtual void v_GetJacobianInverse(NekMatrix<NekDouble> &invJ,
             const Array<OneD, NekDouble> &Au, const Array<OneD, NekDouble> &uu,
           const Array<OneD, NekDouble> &beta, const Array<OneD, NekDouble> &A0,
              const Array<OneD, NekDouble> &alpha, const std::string &type) = 0;

        Array<OneD, MultiRegions::ExpListSharedPtr> m_vessels;
        LibUtilities::SessionReaderSharedPtr m_session;

        NekDouble m_PExt;
        NekDouble m_rho;

    private:
};

inline void PulseWavePressureArea::GetPressure(NekDouble &P,
                 const NekDouble &beta, const NekDouble &A, const NekDouble &A0,
         const NekDouble &dAUdx, const NekDouble &gamma, const NekDouble &alpha)
{
    v_GetPressure(P, beta, A, A0, dAUdx, gamma, alpha);
}

inline void PulseWavePressureArea::GetC(NekDouble &c, const NekDouble &beta,
                const NekDouble &A, const NekDouble &A0, const NekDouble &alpha)
{
    v_GetC(c, beta, A, A0, alpha);
}

inline void PulseWavePressureArea::GetW1(NekDouble &W1, const NekDouble &u,
                 const NekDouble &beta, const NekDouble &A, const NekDouble &A0,
                                                         const NekDouble &alpha)
{
    v_GetW1(W1, u, beta, A, A0, alpha);
}

inline void PulseWavePressureArea::GetW2(NekDouble &W2, const NekDouble &u,
                 const NekDouble &beta, const NekDouble &A, const NekDouble &A0,
                                                         const NekDouble &alpha)
{
    v_GetW2(W2, u, beta, A, A0, alpha);
}

inline void PulseWavePressureArea::GetAFromChars(NekDouble &A,
                                       const NekDouble &W1, const NekDouble &W2,
             const NekDouble &beta, const NekDouble &A0, const NekDouble &alpha)
{
    v_GetAFromChars(A, W1, W2, beta, A0, alpha);
}

inline void PulseWavePressureArea::GetUFromChars(NekDouble &u,
                                       const NekDouble &W1, const NekDouble &W2)
{
    v_GetUFromChars(u, W1, W2);
}

inline void PulseWavePressureArea::GetCharIntegral(NekDouble &I,
                 const NekDouble &beta, const NekDouble &A, const NekDouble &A0,
                                                         const NekDouble &alpha)
{
    v_GetCharIntegral(I, beta, A, A0, alpha);
}

inline void PulseWavePressureArea::GetJacobianInverse(
                                                     NekMatrix<NekDouble> &invJ,
             const Array<OneD, NekDouble> &Au, const Array<OneD, NekDouble> &uu,
           const Array<OneD, NekDouble> &beta, const Array<OneD, NekDouble> &A0,
                   const Array<OneD, NekDouble> &alpha, const std::string &type)
{
    v_GetJacobianInverse(invJ, Au, uu, beta, A0, alpha, type);
}

} // namespace Nektar
#endif
