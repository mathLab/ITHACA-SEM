///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansion2D.cpp
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
// Description: Daughter of StdExpansion. This class contains routine
// which are common to 2D expansion. Typically this inolves physiocal
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdExpansion2D.h>

#ifdef max
#undef max
#endif

namespace Nektar
{
    namespace StdRegions
    {

        StdExpansion2D::StdExpansion2D()
        {
        }

        StdExpansion2D::StdExpansion2D(int numcoeffs, 
				       const LibUtilities::BasisKey &Ba,
                                       const LibUtilities::BasisKey &Bb):
	    StdExpansion(numcoeffs,2, Ba, Bb)
        {
        }

        StdExpansion2D::StdExpansion2D(const StdExpansion2D &T):
                StdExpansion(T)
        {
        }

        StdExpansion2D::~StdExpansion2D()
        {
        }

        //----------------------------
        // Differentiation Methods
        //----------------------------

        void StdExpansion2D::PhysTensorDeriv(ConstNekDoubleSharedArray inarray,
					     NekDoubleSharedArray &outarray_d0, 
					     NekDoubleSharedArray &outarray_d1)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();
            DNekMatSharedPtr D0, D1;
            NekDoubleSharedArray wsp = GetDoubleTmpSpace(nquad0 * nquad1);

            // copy inarray to wsp in case inarray is used as outarray 
	    Vmath::Vcopy(nquad0*nquad1, &inarray[0], 1, &wsp[0], 1);

            D0 = ExpPointsProperties(0)->GetD();
            D1 = ExpPointsProperties(1)->GetD();

            if (outarray_d0) // calculate du/dx_0
            {
                Blas::Dgemm('T', 'N', nquad0, nquad1, nquad0, 1.0,
                            &(D0->GetPtr())[0], nquad0, &wsp[0], nquad0, 0.0,
                            &outarray_d0[0], nquad0);
            }

            // calculate du/dx_1
            if (outarray_d1)
            {
                Blas:: Dgemm('N', 'N', nquad0, nquad1, nquad1, 1.0, &wsp[0], nquad0,
                         &(D1->GetPtr())[0], nquad1, 0.0, &outarray_d1[0], nquad0);
            }

        }

        NekDouble StdExpansion2D::PhysEvaluate2D(ConstNekDoubleSharedArray coords)
        {
            NekDouble val;
            int i;
            int nq0 = m_base[0]->GetNumPoints();
            int nq1 = m_base[1]->GetNumPoints();
            NekDoubleSharedArray wsp1 = GetDoubleTmpSpace(nq1);

            DNekMatSharedPtr I;

            ASSERTL2(coord[0] < -1, "coord[0] < -1");
            ASSERTL2(coord[0] > 1, "coord[0] >  1");
            ASSERTL2(coord[1] < -1, "coord[1] < -1");
            ASSERTL2(coord[1] > 1, "coord[1] >  1");

            // interpolate first coordinate direction
            I = ExpPointsProperties(0)->GetI(coords);
            for (i = 0; i < nq1;++i)
	    {
                wsp1[i] = Blas::Ddot(nq0, &(I->GetPtr())[0], 1, 
				     &m_phys[i * nq0], 1);
	    }

            // interpolate in second coordinate direction
            I = ExpPointsProperties(1)->GetI(coords+1);

            val = Blas::Ddot(nq1, &(I->GetPtr())[0], 1, &wsp1[0], 1);

            return val;
        }

        //////////////////////////////
        /// Integration Methods
        //////////////////////////////

        NekDouble StdExpansion2D::Integral(ConstNekDoubleSharedArray inarray, 
					   ConstNekDoubleSharedArray w0,
					   ConstNekDoubleSharedArray w1)
        {
            int i;
            NekDouble Int = 0.0;
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();
            NekDoubleSharedArray tmp = GetDoubleTmpSpace(nquad0 * nquad1);

            // multiply by integration constants
            for (i = 0; i < nquad1; ++i)
            {
                Vmath::Vmul(nquad0, &inarray[0] + i*nquad0, 1, w0.get(),
                            1, &tmp[0] + i*nquad0, 1);
            }

            for (i = 0; i < nquad0; ++i)
            {
                Vmath::Vmul(nquad1, &tmp[0]+ i, nquad0, w1.get(), 1,
                            &tmp[0] + i, nquad0);
            }

            Int = Vmath::Vsum(nquad0 * nquad1, &tmp[0], 1);

            return Int;
        }

    } //end namespace
} //end namespace


/**
* $Log: StdExpansion2D.cpp,v $
* Revision 1.10  2007/04/05 15:20:11  sherwin
* Updated 2D stuff to comply with SharedArray philosophy
*
* Revision 1.9  2007/03/29 19:35:09  bnelson
* Replaced boost::shared_array with SharedArray
*
* Revision 1.8  2007/03/20 16:58:42  sherwin
* Update to use NekDoubleSharedArray storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
*
* Revision 1.7  2007/03/14 21:24:09  sherwin
* Update for working version of MultiRegions up to ExpList1D
*
* Revision 1.6  2007/03/05 19:28:50  bcarmo
* StdExpansion2D.cpp modified according to StdExpansion1D. Compiles.
*
* Revision 1.6  2007/03/05 22:35:21  bcarmo
* Version with StdExpansion2D compiling
*
* Revision 1.5  2007/01/20 22:35:21  sherwin
* Version with StdExpansion compiling
*
* Revision 1.4  2007/01/18 18:44:45  bnelson
* Updates to compile on Visual Studio 2005.
*
* Revision 1.3  2007/01/17 16:05:35  pvos
* updated doxygen documentation
*
* Revision 1.2  2006/06/01 14:46:16  kirby
* *** empty log message ***
*
* Revision 1.1  2006/05/04 18:58:31  kirby
* *** empty log message ***
*
* Revision 1.16  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.15  2006/04/01 21:59:27  sherwin
* Sorted new definition of ASSERT
*
* Revision 1.14  2006/03/21 09:21:32  sherwin
* Introduced NekMemoryManager
*
* Revision 1.13  2006/02/27 23:47:23  sherwin
*
* Standard coding update upto compilation of StdHexExp.cpp
*
* Revision 1.12  2006/02/26 23:37:29  sherwin
*
* Updates and compiling checks upto StdExpansions1D
*
**/


