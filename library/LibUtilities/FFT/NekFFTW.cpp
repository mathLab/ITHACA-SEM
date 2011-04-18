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

namespace Nektar
{
    namespace LibUtilities
    {
		string NekFFTW::className = NektarFFTFactory::RegisterCreatorFunction("NekFFTW", NekFFTW::create);
		
		NekFFTW::NekFFTW(int N):
		NektarFFT(N)
		//NekFFTW::NekFFTW(int N)
		{
			//m_N = N;
			
			phys = Array<OneD,NekDouble>(m_N);
			coef = Array<OneD,NekDouble>(m_N);
			
			plan_forward  = fftw_plan_r2r_1d(m_N,&phys[0],&coef[0],FFTW_R2HC,FFTW_ESTIMATE);
			plan_backward = fftw_plan_r2r_1d(m_N,&coef[0],&phys[0],FFTW_HC2R,FFTW_ESTIMATE);
			
			m_FFTW_w = Array<OneD,NekDouble>(m_N);
			
			m_FFTW_w[0] = 1.0/(double)m_N;
			m_FFTW_w[1] = m_FFTW_w[0];
			
			for(int i=2;i<m_N;i++)
			{
				m_FFTW_w[i]=2.0*m_FFTW_w[0];
			}
		}
		
		// Distructor
		NekFFTW::~NekFFTW()
		{
			
		}
		//=================================================================================
		// Forward transformation
		void NekFFTW::v_FFTFwdTrans(Array<OneD,NekDouble> &phys, Array<OneD,NekDouble> &coef)
		//void NekFFTW::FFTFwdTrans(Array<OneD,NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
		{
			for(int i=0;i< m_N; i++)
			{
				phys[i] = inarray[i];
			}
			
			fftw_execute(plan_forward);
			
			Reshuffle_FFTW2Nek(coef);
			
			for(int i=0;i< m_N; i++)
			{
				outarray[i] = coef[i];
			}
		}
		
		//=================================================================================
		// Backward transformation
		void NekFFTW::v_FFTBwdTrans(Array<OneD,NekDouble> &coef, Array<OneD,NekDouble> &phys)
		//void NekFFTW::FFTBwdTrans(Array<OneD,NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
		{
			for(int i=0;i< m_N; i++)
			{
				coef[i] = inarray[i];
			}
			
			Reshuffle_Nek2FFTW(coef);
			
			fftw_execute(plan_backward);
			
			for(int i=0;i< m_N; i++)
			{
				outarray[i] = phys[i];
			}
		}
		
		//==================================================================================
		// Reshuffle FFTW2Nek
		void NekFFTW::Reshuffle_FFTW2Nek(Array<OneD,NekDouble> &coef)
		{
			Array<OneD,NekDouble> tmp = Array<OneD,NekDouble>(m_N);
			
			tmp[1]=coef[m_N/2];
			
			for(int i=0;i<m_N/2;i++)
			{
				tmp[2*i]=coef[i];
			}
			
			for(int i=0;i<(m_N/2-1);i++)
			{
				tmp[(m_N-1)-2*i]=coef[m_N/2+1+i];
			}
			
			for(int i=0;i<m_N;i++)
			{
				coef[i]=tmp[i]*m_FFTW_w[i];
			}
			
			if(m_N>2)
			{
				int sign = -1;
				for(int i=2;i<(m_N-1);i=i+2)
				{
					
					coef[i]=coef[i]*sign;
					sign=sign*(-1);
					coef[i+1]=coef[i+1]*sign;
				}
			}
			return;
		}
		
		//===================================================================================
		// Reshuffle Nek2FFTW
		void NekFFTW::Reshuffle_Nek2FFTW(Array<OneD,NekDouble> &coef)
		{
			Array<OneD,NekDouble> tmp = Array<OneD,NekDouble>(m_N);
			
			for(int i=0;i<m_N;i++)
			{
				coef[i]=coef[i]/m_FFTW_w[i];
			}
			
			if(m_N>2)
			{
				int sign = -1;
				for(int i=2;i<(m_N-1);i=i+2)
				{
					
					coef[i]=coef[i]*sign;
					sign=sign*(-1);
					coef[i+1]=coef[i+1]*sign;
				}
			}
			
			tmp[m_N/2]=coef[1];
			
			for(int i=0;i<m_N/2;i++)
			{
				tmp[i]=coef[2*i];
			}
			
			for(int i=0;i<(m_N/2-1);i++)
			{
				tmp[m_N/2+1+i]=coef[(m_N-1)-2*i];
			}
			
			for(int i=0;i<m_N;i++)
			{
				coef[i]=tmp[i]/m_N;
			}
			return;			
		}
		//===================================================================================
	}//end namespace LibUtilities
}//end of namespace Nektar
