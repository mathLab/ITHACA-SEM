///////////////////////////////////////////////////////////////////////////////
//
// File kernel.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Scientific Computing and Imaging Institute,
// University of Utah (USA) and Department of Aeronautics, Imperial
// College London (UK).
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILIITIES_KERNEL_KERNEL_H
#define NEKTAR_LIB_UTILIITIES_KERNEL_KERNEL_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicConst/NektarUnivConsts.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

namespace Nektar
{
	namespace LibUtilities
	{
		class Kernel
		{
		public:
			/**
			 * \brief The default constructor
			 */
			LIB_UTILITIES_EXPORT Kernel();

			LIB_UTILITIES_EXPORT Kernel(int order);

			/**
			 * \brief The default destructor
			 */
			//~Kernel();

			/**
			 * \brief This funciton updates the bspline to the appropriate order.
			 */
			LIB_UTILITIES_EXPORT void UpdateKernelBspline();

			/**
			 * \brief This funciton updates the kernel coefficients.
			 */
			LIB_UTILITIES_EXPORT void UpdateKernelCoeffs();

			/**
			 * \brief This funciton updates the kernel breaks.
			 * \param h represents the mesh spacing
			 */
			LIB_UTILITIES_EXPORT void UpdateKernelBreaks(NekDouble h);

			/**
			 * \brief This funciton returns a 2D array representing the bspline
			 *  of the appropriate order.
			 */
			Array<TwoD,NekDouble> GetKernelBspline()
			{
				return b_spline;
			}

			/**
			 * \brief This funciton returns a 1D array representing the kernel coefficients
			 */
			Array<OneD,NekDouble> GetKernelCoeffs()
			{
				return k_coeffs;
			}

			/**
			 * \brief This funciton returns a 1D array representing the kernel breaks
			 */
			Array<OneD,NekDouble> GetKernelBreaks()
			{
				return k_breaks;
			}

			/**
			 * \brief This funciton sets the k_order variable
			 */
			void UpdateKernelOrder(int order)
			{
				k_order = order;
			}

			/**
			 * \brief This funciton sets the k_ncoeffs variable
			 */
			void UpdateKernelNumOfCoeffs()
			{
				k_ncoeffs = 2*(k_order-1)+1;
			}

			/**
			 * \brief This funciton sets the kernel width size
			 */
			void UpdateKernelWidth()
			{
				k_width = 3*(k_order-1)+1;
			}

			/**
			 * \brief This funciton returns the order of the kernel
			 */
			int GetKernelOrder()
			{
				return k_order;
			}

			/**
			 * \brief This funciton returns the number of kernel coefficients
			 */
			int GetKernelNumeOfCoeffs()
			{
				return k_ncoeffs;
			}

			/**
			 * \brief This funciton returns the size of the kernel width
			 */
			int GetKernelWidth()
			{
				return k_width;
			}

			/**
			 * \brief This funciton moves the center of the kernel to the
			   \param x_value.
			   \param outarray is used to store the result
			 */
			LIB_UTILITIES_EXPORT void MoveKernelCenter(NekDouble x_value, Array<OneD,NekDouble> &outarray);

			/**
			 * \brief This funciton calculates the mesh breaks under the kernel support
			   \param inarray contains the local kernel breaks
			   \param h is the mesh spacing
			   \param outarray contains the coordinate of the mesh breaks under the kernel support
			 */
			LIB_UTILITIES_EXPORT void FindMeshUnderKernel(Array<OneD,NekDouble> &inarray, NekDouble h,
									 Array<OneD,NekDouble> &outarray);

			/**
			 * \brief This funciton evaluates the kernel at input values
			   \param inarray.
			   \param h is the mesh spacing.
			   \param outarray contains the kernel values
			 */
			LIB_UTILITIES_EXPORT void EvaluateKernel(Array<OneD,NekDouble> inarray,NekDouble h,
								Array<OneD,NekDouble> &outarray);

			/**
			 * \brief This function evaluates the bspline at input values
			   \param inarray               input values.
			   \param h                     the mesh spacing.
			   \param offset
			   \param outarray              contains the bspline values.
			 */
			LIB_UTILITIES_EXPORT void EvaluateBspline(Array<OneD,NekDouble> inarray, NekDouble h,
								 NekDouble offset, Array<OneD,NekDouble> &outarray);


			/**
			 * \brief This funciton performs the ordered merge of \param inarray1
			 * and \param inarray2 and puts the result in \param outarray
			 */
			LIB_UTILITIES_EXPORT void Sort(Array<OneD,NekDouble> &inarray1, Array<OneD,NekDouble> &inarray2,
					  Array<OneD,NekDouble> &outarray);

		protected:

			int k_order;						/**< bsplines are obtained by convolving the characteristic
													 fucntion with itself (k_order - 1) times */
			int k_ncoeffs;						/**< Represents the number of kernel coefficients */
			int k_width;						/**< Represents the width of the kernel */
			NekDouble k_center;					/**< holds the center of the kernel */
			Array<TwoD,NekDouble> b_spline;		/**< 2D array representing the bspline */
			Array<OneD,NekDouble> k_coeffs;     /**< 1D array representing the kernel coefficients */
			Array<OneD,NekDouble> k_breaks;     /**< 1D array representing the kernel breaks */

			/**
			 * \brief This funciton evaluates the piecewise bspline polynomial
			   \param interval \param interval at point \param x_value.
			 */
			NekDouble EvaluateBsplinePoly(NekDouble x_value,int interval);

		private:

		};

		typedef std::shared_ptr<Kernel>      KernelSharedPtr;

	}//end of namespace
}// end of namespace

#endif //NEKTAR_LIB_UTILIITIES_KERNEL_KERNEL_H
