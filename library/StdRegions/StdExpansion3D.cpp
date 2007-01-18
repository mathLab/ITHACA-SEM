///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansioneD.cpp
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
// which are common to 3D expansion. Typically this inolves physiocal
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdExpansion3D.h>

#ifdef max
#undef max
#endif

namespace Nektar
{
  namespace StdRegions
  {


    StdExpansion3D::StdExpansion3D()
    {
    }
  
    StdExpansion3D::StdExpansion3D(const BasisKey &Ba, const BasisKey &Bb, 
				   const BasisKey &Bc, int numcoeffs, 
				   double *coeffs, double *phys, 
				   bool spaceowner):
      StdExpansion(3,Ba,Bb,Bc,numcoeffs,coeffs,phys,spaceowner)
    {
    }
  
    StdExpansion3D::StdExpansion3D(const StdExpansion3D &T):
      StdExpansion(T)
    {
    }
  
    StdExpansion3D::~StdExpansion3D()
    { 
    }
  
    void StdExpansion3D::TensorDeriv(const double *inarray, 
				     double *outarray_d0, 
				     double *outarray_d1, 
				     double *outarray_d2)
    {
      int    i;
      int    nquad0 = m_base[0]->GetPointsOrder();
      int    nquad1 = m_base[1]->GetPointsOrder();
      int    nquad2 = m_base[2]->GetPointsOrder();
      const double *D0,*D1,*D2;
      BstShrDArray wsp;
      double * tmp;
     
      // check to see if either calling array is inarray
      if((outarray_d0 == inarray)||(outarray_d1 == inarray)||
	 (outarray_d2 == inarray))
      { 
	wsp = GetDoubleTmpSpace(nquad0*nquad1*nquad2);
	tmp = wsp.get();
	Vmath::Vcopy(nquad0*nquad1*nquad2,inarray,1,tmp,1);
      }
      else
      {
	tmp = (double *) inarray;
      }
      
      BasisManagerSingleton::Instance().GetD(m_base[0],D0);
      BasisManagerSingleton::Instance().GetD(m_base[1],D1);
      BasisManagerSingleton::Instance().GetD(m_base[2],D2);

      // calculate du/dx_0
      if(outarray_d0)
      {
	for(i=0; i < nquad2; ++i)
	{
	  Blas::Dgemm('T','N',nquad0,nquad1,nquad0,1.0,D0,nquad0,
		      tmp+i*nquad0*nquad1, nquad0,0.0,
		      outarray_d0+i*nquad0*nquad1,nquad0);
	}
      }

      // calculate du/dx_1
      if(outarray_d1)
      {
	for(i=0; i < nquad2; ++i)
	{
	  Blas:: Dgemm('N','N',nquad0,nquad1,nquad1,1.0,
		       tmp+i*nquad0*nquad1,nquad0,D1,nquad1,0.0,
		       outarray_d1+i*nquad0*nquad1,nquad0);
	}
      }
    
      // calculate du/dx_2
      if(outarray_d2)
      {
	for(i=0; i < nquad0*nquad1; ++i)
	{
	  Blas:: Dgemv('T',nquad2,nquad2,1.0,D2,nquad2,
		       tmp+i,nquad0*nquad1,0.0,
		       outarray_d1+i,nquad0*nquad1);
	}
      }    
    }
  
    void StdExpansion3D::TensorDeriv(double *outarray_d0, 
				    double *outarray_d1,
				    double *outarray_d2)
    {
      TensorDeriv(this->m_phys,outarray_d0,outarray_d1,outarray_d2);
    }

  
    double StdExpansion3D::PhysEvaluate(const double *coords)
    {
      double  val;
      int i;
      int nq1 = m_base[0]->GetPointsOrder();
      int nq2 = m_base[1]->GetPointsOrder();
      int nq3 = m_base[2]->GetPointsOrder();


      BstShrDArray Ivec = GetDoubleTmpSpace(std::max(std::max(nq1,nq2),nq3));
      BstShrDArray wsp  = GetDoubleTmpSpace(nq2*nq3);
      double *tmp = wsp.get();
      BstShrDArray wsp1 = GetDoubleTmpSpace(nq3);
      double *tmp1 = wsp1.get();
    
      ASSERTL2(coord[0] < -1,"coord[0] < -1");
      ASSERTL2(coord[0] >  1,"coord[0] >  1");
      ASSERTL2(coord[1] < -1,"coord[1] < -1");
      ASSERTL2(coord[1] >  1,"coord[1] >  1");
      ASSERTL2(coord[2] < -1,"coord[2] < -1");
      ASSERTL2(coord[2] >  1,"coord[2] >  1");
      
      // interpolate first coordinate direction
      m_base[0]->GetInterpVec(coords[0],Ivec.get());
      for(i = 0; i < nq2*nq3;++i)
      {
	tmp[i] =  Blas::Ddot(nq1,Ivec.get(),1,m_phys+i*nq1,1);
      }
    
      // interpolate in second coordinate direction 
      m_base[1]->GetInterpVec(coords[1],Ivec.get());
      for(i =0; i < nq3; ++i)
      {
	tmp1[i] = Blas::Ddot(nq2,tmp+i*nq2,1,Ivec.get(),1);
      }
    
      // interpolate in third coordinate direction 
      m_base[2]->GetInterpVec(coords[2],Ivec.get());
      val = Blas::Ddot(nq3,tmp1,1,Ivec.get(),1);
    
      return val;    
    }

  }//end namespace
}//end namespace

/** 
 * $Log: StdExpansion3D.cpp,v $
 * Revision 1.3  2007/01/17 16:05:40  pvos
 * updated doxygen documentation
 *
 * Revision 1.2  2006/06/01 14:46:16  kirby
 * *** empty log message ***
 *
 * Revision 1.1  2006/05/04 18:58:31  kirby
 * *** empty log message ***
 *
 * Revision 1.9  2006/04/01 21:59:27  sherwin
 * Sorted new definition of ASSERT
 *
 * Revision 1.8  2006/03/21 09:21:32  sherwin
 * Introduced NekMemoryManager
 *
 * Revision 1.7  2006/02/27 23:47:23  sherwin
 *
 * Standard coding update upto compilation of StdHexExp.cpp
 *
 **/ 



