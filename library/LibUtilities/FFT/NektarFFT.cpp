///////////////////////////////////////////////////////////////////////////////
//
// File NektarFFT.cpp
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
// Description: Fast Fourier Transform base class in Nektar++
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/FFT/NektarFFT.h>

namespace Nektar
{
	namespace LibUtilities
	{
		/**
		 * @class NektarFFT
		 *
		 * This class is a base class for all FFT implementation. It provides
		 * the underlying generic functionality and interface for perform a Fourier transform.
		 *
		 * To perform the FFT with a differebt algorithm, create a derived class from this class
		 * and reimplement the virtual functions to provide custom implementation
		 * of the algorithm.
		 *
		 */
		
		/**
		 * This constructor is protected as the objects of this class are never
		 * instantiated directly.
		 */
		NektarFFT::NektarFFT(int N)
		{
			m_N = N;
		}
		
		NektarFFT::~NektarFFT()
		{
			
		}
		
		NektarFFTFactory& GetNektarFFTFactory()
		{
                    static NektarFFTFactory instance;
                    return instance;
		}

		/**
		 * This allows initialisation of the class which cannot be completed
		 * during object construction (such as setting of initial conditions).
		 *
		 * Public interface routine to virtual function implementation.
		 */
		
		void NektarFFT::FFTFwdTrans(Array<OneD,NekDouble> &phys, Array<OneD,NekDouble> &coef)
		{
			v_FFTFwdTrans(phys,coef);
		}
		
		void NektarFFT::FFTBwdTrans(Array<OneD,NekDouble> &coef, Array<OneD,NekDouble> &phys)
		{
			v_FFTBwdTrans(coef,phys);
		}
		
		void NektarFFT::v_FFTFwdTrans(Array<OneD,NekDouble> &phys, Array<OneD,NekDouble> &coef)
		{
            boost::ignore_unused(phys, coef);
		}
		
		void NektarFFT::v_FFTBwdTrans(Array<OneD,NekDouble> &coef, Array<OneD,NekDouble> &phys)
		{
            boost::ignore_unused(coef, phys);
		}
		
	}//end namespace LibUtilities
}//end of namespace Nektar
