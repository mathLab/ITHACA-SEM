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

        StdExpansion2D::StdExpansion2D(const BasisKey &Ba, const BasisKey &Bb, 
            int numcoeffs, double *coeffs, double *phys, 
            bool spaceowner):
        StdExpansion(2,Ba,Bb,BasisKey(),numcoeffs,coeffs,phys,spaceowner)
        {
        }

        StdExpansion2D::StdExpansion2D(const StdExpansion2D &T):
        StdExpansion(T)
        {
        }

        StdExpansion2D::~StdExpansion2D()
        { 
        }

        void StdExpansion2D::TensorDeriv(const double *inarray, 
            double *outarray_d0,double *outarray_d1)
        {
            int    nquad0 = m_base[0]->GetPointsOrder();
            int    nquad1 = m_base[1]->GetPointsOrder();
            const double *D0,*D1;
            BstShrDArray wsp = GetDoubleTmpSpace(nquad0*nquad1);
            double *tmp = wsp.get();

            // check to see if either calling array is inarray
            if((outarray_d0 == inarray)||(outarray_d1 == inarray))
            { 
                Vmath::Vcopy(nquad0*nquad1,inarray,1,tmp,1);
            }
            else
            {
                tmp =  (double *)inarray;
            }

            BasisManagerSingleton::Instance().GetD(m_base[0],D0);
            BasisManagerSingleton::Instance().GetD(m_base[1],D1);

            if(outarray_d0) // calculate du/dx_0
            {
                Blas::Dgemm('T','N',nquad0,nquad1,nquad0,1.0,D0,nquad0,tmp,
                    nquad0,0.0,outarray_d0,nquad0);
            }

            // calculate du/dx_1
            if(outarray_d1)
            {
                Blas:: Dgemm('N','N',nquad0,nquad1,nquad1,1.0,tmp,nquad0,D1,
                    nquad1,0.0,outarray_d1,nquad0);
            }

        }

        void StdExpansion2D::TensorDeriv(double *outarray_d0, 
            double *outarray_d1)
        {
            TensorDeriv(this->m_phys,outarray_d0,outarray_d1);
        }

        double StdExpansion2D::PhysEvaluate(const double *coords)
        {
            double  val;
            int i;
            int nq0 = m_base[0]->GetPointsOrder();
            int nq1 = m_base[1]->GetPointsOrder();
            BstShrDArray tmp  = GetDoubleTmpSpace(std::max(nq0,nq1));
            BstShrDArray wsp1 = GetDoubleTmpSpace(nq1);
            double  *tmp1 = wsp1.get();

            ASSERTL2(coord[0] < -1,"coord[0] < -1");
            ASSERTL2(coord[0] >  1,"coord[0] >  1");
            ASSERTL2(coord[1] < -1,"coord[1] < -1");
            ASSERTL2(coord[1] >  1,"coord[1] >  1");

            m_base[0]->GetInterpVec(coords[0],tmp.get());

            // interpolate first coordinate direction
            for(i = 0; i < nq1;++i)
            {
                tmp1[i] =  Blas::Ddot(nq0,tmp.get(),1,m_phys+i*nq0,1);
            }

            // interpolate in second coordinate direction 
            m_base[1]->GetInterpVec(coords[1],tmp.get());

            val = Blas::Ddot(nq1,tmp.get(),1,tmp1,1);

            return val;    
        }


        //////////////////////////////
        /// Integration Methods
        //////////////////////////////

        double StdExpansion2D::Integral(const double *inarray, const double *w0, 
            const double* w1)
        {
            int    i;
            double Int = 0.0;
            int    nquad0 = m_base[0]->GetPointsOrder();
            int    nquad1 = m_base[1]->GetPointsOrder();
            BstShrDArray tmp  = GetDoubleTmpSpace(nquad0*nquad1);

            // multiply by integration constants 
            for(i = 0; i < nquad1; ++i)
            {
                Vmath::Vmul(nquad0,(double*)inarray+i*nquad0,1,(double*)w0,1,
                    tmp.get()+i*nquad0,1);
            }

            for(i = 0; i < nquad0; ++i)
            {
                Vmath::Vmul(nquad1,tmp.get()+i,nquad0,(double*)w1,1,tmp.get()+i, nquad0);
            }

            Int = Vmath::Vsum(nquad0*nquad1,tmp.get(),1);

            return Int;
        }

    }//end namespace
}//end namespace


/** 
* $Log: StdExpansion2D.cpp,v $
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


