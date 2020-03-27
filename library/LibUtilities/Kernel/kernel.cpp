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

#include <LibUtilities/Kernel/kernel.h>

namespace Nektar
{
	namespace LibUtilities
	{
		Kernel::Kernel(int order)
		{
			UpdateKernelOrder(order);
			UpdateKernelNumOfCoeffs();
			UpdateKernelWidth();
			UpdateKernelBspline();
			UpdateKernelCoeffs();
			k_center = 0.0;
		}

		void Kernel::UpdateKernelBspline()
		{
			Array<TwoD,NekDouble> out_final(k_order,k_order);

			int i,j;

			if(k_order == 1.0)
			{
				b_spline[0][0] = 1.0;

			}else if(k_order == 2)
			{
				NekDouble out[2][2] = {{1.0,0.0},
									{-1.0,1.0}};
				for(i = 0; i < k_order; i++)
				{
					for(j = 0; j < k_order; j++)
					{
						out_final[i][j] = out[i][j];
					}
				}

			}else if(k_order == 3)
			{
				NekDouble out[3][3] = {{0.5, 0.0, 0.0},
									{-1.0, 1.0, 0.5},
									{0.5, -1.0, 0.5}};
				for(i = 0; i < k_order; i++)
				{
					for(j = 0; j < k_order; j++)
					{
						out_final[i][j] = out[i][j];
					}
				}

			}else if(k_order == 4)
			{
				NekDouble out[4][4] = {{1.0/6.0, 0, 0, 0},
									{-0.5, 0.5, 0.5, 1.0/6.0},
									{0.5, -1.0, 0, 2.0/3},
									{-1.0/6.0, 0.5, -0.5, 1.0/6.0}};
				for(i = 0; i < k_order; i++)
				{
					for(j = 0; j < k_order; j++)
					{
						out_final[i][j] = out[i][j];
					}
				}

			}else if(k_order == 5)
			{
				NekDouble out[5][5] = {{1.0 / 24.0, 0, 0, 0, 0},
									{-1.0/6, 1.0/6, 0.25, 1.0/6, 1.0/24},
									{0.25, -0.5, -0.25, 0.5, 11.0/24},
									{-1.0/6, 0.5, -0.25, -0.5, 11.0/24},
									{1.0/24, -1.0/6.0, 0.25, -1.0/6, 1.0/24}};
				for(i = 0; i < k_order; i++)
				{
					for(j = 0; j < k_order; j++)
					{
						out_final[i][j] = out[i][j];
					}
				}

			}else if(k_order == 6)
			{
				NekDouble out[6][6] = {{1.0/1.020, 0, 0, 0, 0, 0},
									{-1.0/24, 1.0/24, 1.0/12, 1.0/12, 1.0/24, 1.0/120},
									{1.0/12, -1.0/6, -1.0/6, 1.0/6, 5.0/12, 13.0/60},
									{-1.0/12, 0.25, 0, -1.0/2, 0, 11.0/20},
									{1.0/24, -1.0/6, 1.0/6, 1.0/6, -5.0/12, 13.0/60},
									{-1.0/120, 1.0/24, -1.0/12, 1.0/12, -1.0/24, 1.0/120}};
				for(i = 0; i < k_order; i++)
				{
					for(j = 0; j < k_order; j++)
					{
						out_final[i][j] = out[i][j];
					}
				}

			}else if(k_order == 7)
			{
				NekDouble out[7][7] = {{1.0/720, 0, 0, 0, 0, 0, 0},
									{-1.0/120, 1.0/120, 1.0/48, 1.0/36, 1.0/48, 1.0/120, 1.0/720},
									{1.0/48, -1.0/24, -1.0/16, 1.0/36, 3.0/16, 5.0/24, 19.0/240},
									{-1.0/36, 1.0/1.02, 1.0/24, -2.0/9, -5.0/24, 1.0/3, 151.0/360},
									{1.0/48, -1.0/12, 1.0/24, 2.0/9, -5.0/24, -1.0/3, 151.0/360},
									{-1.0/120, 1.0/24, -1.0/16, -1.0/36, 3.0/16, -5.0/24, 19.0/240},
									{1.0/720, -1.0/120, 1.0/48, -1.0/36, 1.0/48, -1.0/120, 1.0/720}};
				for(i = 0; i < k_order; i++)
				{
					for(j = 0; j < k_order; j++)
					{
						out_final[i][j] = out[i][j];
					}
				}

			}else if(k_order == 8)
			{
				NekDouble out[8][8] = {{1.0/5040, 0, 0, 0, 0, 0, 0, 0},
									{-1.0/720, 1.0/720, 1.0/240, 1.0/144, 1.0/144, 1.0/240, 1.0/720, 1.0/5040 },
									{1.0/240, -1.0/120, -1.0/60, 0, 1.0/18, 1.0/10, 7.0/90, 1.0/42},
									{-1.0/144, 1.0/48, 1.0/48, -1.0/16, -19.0/144, 1.0/16, 49.0/144, 397.0/1680},
									{1.0/144, -1.0/36, 0, 1.0/9, 0, -1.0/3, 0, 151.0/315},
									{-1.0/240, 1.0/48, -1.0/48, -1.0/16, 19.0/144, 1.0/16, -49.0/144, 397.0/1680},
									{ 1.0/720, -1.0/120, 1.0/60, 0, -1.0/18, 1.0/10, -7.0/90, 1.0/42},
									{ -1.0/5040, 1.0/720, -1.0/240, 1.0/144, -1.0/144, 1.0/240, -1.0/720, 1.0/5040}};
				for(i = 0; i < k_order; i++)
				{
					for(j = 0; j < k_order; j++)
					{
						out_final[i][j] = out[i][j];
					}
				}

			}

			b_spline = out_final;

			ASSERTL0(k_order <= 8, "Order is not supported");


		}

		void Kernel::UpdateKernelCoeffs()
		{
			int i;
			Array<OneD,NekDouble> out_final(k_ncoeffs);

			if (k_order == 2)
			{
				NekDouble out[3] = {-1.0/12, 7.0/6, -1.0/12};
				for(i = 0; i < k_ncoeffs; i++)
				{
					out_final[i] = out[i];
				}
			}else if (k_order == 3)
			{
				NekDouble out[5] = {37.0/1920, -97.0/480, 437.0/320, -97.0/480, 37.0/1920};
				for(i = 0; i < k_ncoeffs; i++)
				{
					out_final[i] = out[i];
				}
			}else if (k_order == 4)
			{
				NekDouble out[7] = {-41.0/7560, 311.0/5040, -919.0/2520, 12223.0/7560, -919.0/2520, 311.0/5040, -41.0/7560};
				for(i = 0; i < k_ncoeffs; i++)
				{
					out_final[i] = out[i];
				}
			}else if (k_order == 5)
			{
				NekDouble out[9] = {153617.0/92897280, -35411.0/1658880, 3153959.0/23224320, -6803459.0/11612160, 18017975.0/9289728, -6803459.0/11612160, 3153959.0/23224320, -35411.0/1658880, 153617.0/92897280};
				for(i = 0; i < k_ncoeffs; i++)
				{
					out_final[i] = out[i];
				}
			}else if (k_order == 6)
			{
				NekDouble out[11] = {-4201.0/7983360, 30773.0/3991680, -20813.0/380160, 2825.0/11088, -1179649.0/1330560, 1569217.0/665280, -1179649.0/1330560, 2825.0/11088, -20813.0/380160, 30773.0/3991680, -4201.0/7983360};
				for(i = 0; i < k_ncoeffs; i++)
				{
					out_final[i] = out[i];
				}
			}else if (k_order == 7)
			{
				NekDouble out[13] = {13154671847.0/76517631590400.0, -18073154507.0/6376469299200.0, 287360344573.0/12752938598400.0, -2217732343517.0/19129407897600.0, 1240941746699.0/2833986355200.0, -275386671493.0/212548976640.0, 2648644782397.0/910924185600.0, -275386671493.0/212548976640.0, 1240941746699.0/2833986355200.0, -2217732343517.0/19129407897600.0, 287360344573.0/12752938598400.0, -18073154507.0/6376469299200.0, 13154671847.0/76517631590400.0};
				for(i = 0; i < k_ncoeffs; i++)
				{
					out_final[i] = out[i];
				}
			}else if (k_order == 8)
			{
				NekDouble out[15] = {-800993.0/14010796800.0, 73587167.0/70053984000.0, -651305719.0/70053984000.0, 3714581677.0/70053984000.0, -3085236289.0/14010796800.0, 1426328231.0/2001542400.0, -43268401973.0/23351328000.0, 42401344373.0/11675664000.0, -43268401973.0/23351328000.0, 1426328231.0/2001542400.0, -3085236289.0/14010796800.0, 3714581677.0/70053984000.0, -651305719.0/70053984000.0, 73587167.0/70053984000.0, -800993.0/14010796800.0};
				for(i = 0; i < k_ncoeffs; i++)
				{
					out_final[i] = out[i];
				}
			}

			k_coeffs = out_final;
			ASSERTL0(k_order <= 8, "Order is not supported");

		}

		void Kernel::UpdateKernelBreaks(NekDouble h)
		{
			int i;
			Array<OneD,NekDouble> temp(k_width+1);
			temp[0] = -(k_width/2.0)*h; // it needs to get scaled by h
			for(i = 1; i < k_width+1; i++)
			{
				temp[i] = temp[i-1]+h;
			}
			k_breaks = temp;
		}

		void Kernel::MoveKernelCenter(NekDouble x_value, Array<OneD,NekDouble> &outarray)
		{
			int i;
			for(i = 0; i < k_width+1; i++)
			{
				outarray[i] = k_breaks[i] + x_value;
			}

			// Update the center of the kernel
			k_center = x_value;
		}

		void Kernel::FindMeshUnderKernel(Array<OneD,NekDouble> &inarray, NekDouble h,
										 Array<OneD,NekDouble> &outarray)
		{
			int j;
			NekDouble first = ceil(inarray[0]/h)*h;
			int index = k_width;
			NekDouble last = floor(inarray[index]/h)*h;
			int count = (int)((last-first)/h)+1; // number of mesh breaks under the kernel support
			Array<OneD,NekDouble> mesh_breaks(count);
			mesh_breaks[0] = first;
			for(j = 1; j < count; j++)
			{
				mesh_breaks[j] = mesh_breaks[j-1]+h;
			}
			outarray = mesh_breaks;

		}

		void Kernel::EvaluateKernel(Array<OneD,NekDouble> inarray,NekDouble h,
									Array<OneD,NekDouble> &outarray)
		{
			int gamma,i;
			int degree = k_order - 1;
			int nvalues = inarray.size();
			Array<OneD,NekDouble> bs_values(nvalues);

			for(i = 0; i < nvalues; i++ )
			{
				outarray[i] = 0.0;
			}

			for(gamma = -degree; gamma <= degree; gamma++)
			{
				int cIndex = gamma+degree;

				// Evaluate the bSpline values
				EvaluateBspline(inarray,h,k_center+(gamma*h),bs_values);

				for(i = 0; i < nvalues; i++ )
				{
					outarray[i] += k_coeffs[cIndex]*bs_values[i];
				}
			}

		}

		void Kernel::EvaluateBspline(Array<OneD,NekDouble> inarray, NekDouble h,
									 NekDouble offset, Array<OneD,NekDouble> &outarray)
		{
			int i;
			NekDouble min_value = -k_order/2.0;
			NekDouble max_value = k_order/2.0;

			int nvalues = inarray.size();

			// Make a copy for further modifications
			Array<OneD,NekDouble> inarray_cp(nvalues);

			for(i = 0; i < nvalues; i++)
			{
				inarray_cp[i] = inarray[i] - offset;
				inarray_cp[i] = inarray_cp[i]/h;
				int interval = (int)floor(inarray_cp[i] - min_value); // determines to which interval of the bspline the value belongs

				if(inarray_cp[i] >= min_value && inarray_cp[i] <= max_value)
				{
					if(interval >= k_order)
					{
						interval -= 1;
					}
					NekDouble shift = min_value + interval;
					inarray_cp[i] -= shift;
					outarray[i] = EvaluateBsplinePoly(inarray_cp[i],interval);
				}else
				{
					outarray[i] = 0.0;
				}

			}
		}

		NekDouble Kernel::EvaluateBsplinePoly(NekDouble x_value,int interval)
		{
			int i;
			int deg = k_order - 1;
			NekDouble poly_value = b_spline[interval][0];

			for(i = 0; i < deg; i++)
			{
				poly_value = poly_value*x_value + b_spline[interval][i+1];

			}

			return poly_value;
		}

		void Kernel::Sort(Array<OneD,NekDouble> &inarray1, Array<OneD,NekDouble> &inarray2,
							   Array<OneD,NekDouble> &outarray)
		{
			int j;
			int kIndex = 0; // Keeps track of the kernel breaks
			int mIndex = 0; // Keeps track of the mesh breaks
			for(j = 0; j < outarray.size(); j++)
			{
				if(mIndex >= inarray2.size())
				{
					outarray[j] = inarray1[kIndex];
					kIndex++;
				}else if(kIndex >= inarray1.size())
				{
					outarray[j] = inarray2[mIndex];
					mIndex++;

				}else if(inarray1[kIndex] < inarray2[mIndex])
				{
					outarray[j] = inarray1[kIndex];
					kIndex++;
				}else
				{
					outarray[j] = inarray2[mIndex];
					mIndex++;
				}

			}
		}

	}// end of namespace
}// end of namespace
