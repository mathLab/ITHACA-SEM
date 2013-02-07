///////////////////////////////////////////////////////////////////////////////
//
// File NekFFTW.cpp
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
// Description: Wrapper around FFTW library
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/FFT/NekFFTW.h>
#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
    namespace LibUtilities
    {
        string NekFFTW::className
            = GetNektarFFTFactory().RegisterCreatorFunction("NekFFTW",
                                                            NekFFTW::create);

        NekFFTW::NekFFTW(int N)
                : NektarFFT(N)
        {
            m_wsp = Array<OneD, NekDouble>(m_N);
            phys = Array<OneD,NekDouble>(m_N);
            coef = Array<OneD,NekDouble>(m_N);

            plan_forward  = fftw_plan_r2r_1d(m_N, &phys[0], &coef[0],
                                             FFTW_R2HC, FFTW_ESTIMATE);
            plan_backward = fftw_plan_r2r_1d(m_N, &coef[0], &phys[0],
                                             FFTW_HC2R, FFTW_ESTIMATE);

            m_FFTW_w = Array<OneD,NekDouble>(m_N);
            m_FFTW_w_inv = Array<OneD,NekDouble>(m_N);

            m_FFTW_w[0] = 1.0/(NekDouble)m_N;
            m_FFTW_w[1] = 0.0;

            m_FFTW_w_inv[0] = m_N;
            m_FFTW_w_inv[1] = 0.0;

            for(int i=2;i<m_N;i++)
            {
                m_FFTW_w[i]     = m_FFTW_w[0];
                m_FFTW_w_inv[i] = m_FFTW_w_inv[0];
            }
        }

        // Distructor
        NekFFTW::~NekFFTW()
        {

        }

        // Forward transformation
        void NekFFTW::v_FFTFwdTrans(
                Array<OneD,NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray)
        {
            Vmath::Vcopy(m_N, inarray, 1, phys, 1);

            fftw_execute(plan_forward);

            Reshuffle_FFTW2Nek(coef);

            Vmath::Vcopy(m_N, coef, 1, outarray, 1);
        }

        // Backward transformation
        void NekFFTW::v_FFTBwdTrans(
                Array<OneD,NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray)
        {
            Vmath::Vcopy(m_N, inarray, 1, coef, 1);

            Reshuffle_Nek2FFTW(coef);

            fftw_execute(plan_backward);

            Vmath::Vcopy(m_N, phys, 1, outarray, 1);
        }

        // Reshuffle FFTW2Nek
        void NekFFTW::Reshuffle_FFTW2Nek(Array<OneD,NekDouble> &coef)
        {
            int halfN = m_N/2;

            m_wsp[1] = coef[halfN];

            Vmath::Vcopy(halfN, coef, 1, m_wsp, 2);

            for(int i = 0; i < (halfN - 1); i++)
            {
                m_wsp[(m_N-1)-2*i] = coef[halfN+1+i];
            }

            Vmath::Vmul(m_N, m_wsp, 1, m_FFTW_w, 1, coef, 1);

            return;
        }

        // Reshuffle Nek2FFTW
        void NekFFTW::Reshuffle_Nek2FFTW(Array<OneD,NekDouble> &coef)
        {
            int halfN = m_N/2;

            Vmath::Vmul(m_N, coef, 1, m_FFTW_w_inv, 1, coef, 1);

            m_wsp[halfN] = coef[1];

            Vmath::Vcopy(halfN, coef, 2, m_wsp, 1);

            for(int i = 0; i < (halfN-1); i++)
            {
                m_wsp[halfN+1+i] = coef[(m_N-1)-2*i];
            }

            Vmath::Smul(m_N, 1.0/m_N, m_wsp, 1, coef, 1);

            return;
        }
    }
}
