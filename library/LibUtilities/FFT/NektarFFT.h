///////////////////////////////////////////////////////////////////////////////
//
// File NektarFFT.h
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
// Description: Header file for the Fast Fourier Transform class in Nektar++
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILIITIES_FFT_NEKTARFFT_H
#define NEKTAR_LIB_UTILIITIES_FFT_NEKTARFFT_H

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

namespace Nektar
{
    template <typename Dim, typename DataType>
    class Array;
    
	namespace LibUtilities
	{
		/**
		 * The NektarFFT class is a virtual class to manage the use of the FFT to do Fwd/Bwd transformations
		 * and convolutions. The function here defined will link to a proper implementation of the FFT algorithm.
		 * Depending on the user definition the functions can link to a class which is a wrapper around the FFTW 
		 * library or to a specific FFT implementation.
		 */
		class NektarFFT;
		
		// A shared pointer to the NektarFFT object
		typedef std::shared_ptr<NektarFFT>  NektarFFTSharedPtr;
		
		/// Datatype of the NekFactory used to instantiate classes derived from
		/// the NektarFFT class.
		typedef LibUtilities::NekFactory< std::string, NektarFFT, int> NektarFFTFactory;
		
		LIB_UTILITIES_EXPORT NektarFFTFactory& GetNektarFFTFactory();

		class NektarFFT
		{
		public:
			
			/// Initialises NektarFFT class members.
			LIB_UTILITIES_EXPORT NektarFFT(int N);
			
			// Distructor
			LIB_UTILITIES_EXPORT  ~NektarFFT();
			
			/**
			 * m_N is the dimension of the Fourier transform.
			 * It means that the coefficient vector and the vector of the variable in physical
			 * space have size m_N. It is becasue everything is managed just with real data.
			 */
			int m_N;
		
			/**
			 * Forward transformation to pass from physical to coefficient space using the FFT.
			 * This method will take the place of the Matrix-Vector multiplication
			 * input:
			 * N         = number of Fourier points
			 * inarray   = vector in physical space (length N)
			 * output:
			 * outarray  = vector in coefficient space (length N)
			 */
			LIB_UTILITIES_EXPORT void FFTFwdTrans(Array<OneD,NekDouble> &phy, Array<OneD,NekDouble> &coef);
			
			/**
			 * Backward transformation to pass from coefficient to physical space using the FFT.
			 * This method will take the place of the Matrix-Vector multiplication
			 * input:
			 * N          = number of Fourier points
			 * inarrray   = vector in coefficient space (length N)
			 * output:
			 * outarray   = vector in physical space (length N)
			 */
			LIB_UTILITIES_EXPORT void FFTBwdTrans(Array<OneD,NekDouble> &coef, Array<OneD,NekDouble> &phys);
			
		protected:
			
			
			virtual void v_FFTFwdTrans(Array<OneD,NekDouble> &phys, Array<OneD,NekDouble> &coef);
						
			virtual void v_FFTBwdTrans(Array<OneD,NekDouble> &coef, Array<OneD,NekDouble> &phys);
			
		private:
			
		};
	}//end namespace LibUtilities
}//end of namespace Nektar
#endif //NEKTAR_LIB_UTILIITIES_FFT_NEKTARFFT_H
