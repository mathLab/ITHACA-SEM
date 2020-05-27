///////////////////////////////////////////////////////////////////////////////
//
// File NekFFTW.h
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
// Description: Header file for the wrapper around FFTW library
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILIITIES_FFT_NEKFFTW_H
#define NEKTAR_LIB_UTILIITIES_FFT_NEKFFTW_H

#include <LibUtilities/FFT/NektarFFT.h>

#include <LibUtilities/BasicUtils/NekManager.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicConst/NektarUnivConsts.hpp>

#include <fftw3.h>

namespace Nektar
{
    template <typename Dim, typename DataType>
    class Array;

	namespace LibUtilities
	{
		class NekFFTW;
		
		// A shared pointer to the NekFFTW object
		typedef std::shared_ptr<NekFFTW>  NekFFTWSharedPtr;
		
		class NekFFTW: public NektarFFT
		{
		public:
			
			/// Creates an instance of this class
			static NektarFFTSharedPtr create(int N) 
			{
				return MemoryManager<NekFFTW>::AllocateSharedPtr(N);
			}
			
			/// Name of class
			static std::string className;
			            
			// constructor (initialisation of the FFTW planes and fill up the m_FFTW_w vector)
			NekFFTW(int N);
			
			// Distructor
			virtual ~NekFFTW();
			
			
			
			virtual void v_FFTFwdTrans(Array<OneD,NekDouble> &inarray, Array<OneD,NekDouble> &outarray);
			
			virtual void v_FFTBwdTrans(Array<OneD,NekDouble> &inarray, Array<OneD,NekDouble> &outarray);
			
			
			
		protected:
			
			Array<OneD,NekDouble> m_FFTW_w;  // weights to convert arrays form Nektar++ to FFTW format
			Array<OneD,NekDouble> m_FFTW_w_inv; // weights to convert arrays from FFTW to Nektar++ format

			Array<OneD,NekDouble> phys;
			Array<OneD,NekDouble> coef;

			Array<OneD,NekDouble> m_wsp;     // Workspace area for transforms

			fftw_plan plan_backward;         // plan to execute a backward FFT in FFTW
			fftw_plan plan_forward;          // plan to execute a forward FFT in FFTW
			/**
			 * Reshuffling routines to put the coefficients in Nektar++/FFTW format.
			 * The routines take as an input the number of points N, the vector of coeffcients
			 * and the vector containing the weights of the numerical integration.
			 * The routines modify directly the coef vector.
			 */
			void Reshuffle_FFTW2Nek(Array<OneD,NekDouble> &coef);
			
			void Reshuffle_Nek2FFTW(Array<OneD,NekDouble> &coef);
			
			

		private:
		};

	}//end namespace LibUtilities
}//end of namespace Nektar
#endif //NEKTAR_LIB_UTILIITIES_FFT_NEKFFTW_H
