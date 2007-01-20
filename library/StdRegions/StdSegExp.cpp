///////////////////////////////////////////////////////////////////////////////
//
// File StdSegExp.cpp
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
// Description: Routines within Standard Segment Expansions
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdSegExp.h>


namespace Nektar
{
  namespace StdRegions
  {

    StdMatrix StdSegExp::s_elmtmats;
    
    StdSegExp::StdSegExp(const BasisKey &Ba):
      StdExpansion1D(Ba,Ba.GetNumBasis(),NULL,NULL,true)
    {    
    }

    StdSegExp::StdSegExp(const BasisKey &Ba,double *coeffs, double *phys):
      StdExpansion1D(Ba,Ba.GetNumBasis(),coeffs,phys,false)
    {
    }

    StdSegExp::StdSegExp(const StdSegExp &T):StdExpansion1D(T)
    {
    }
  

    StdSegExp::~StdSegExp()
    {    
    }
    
  
    //----------------------------
    // Integration Methods
    //----------------------------

    double StdSegExp::Integral(const double *inarray)
    {
      double Int = 0.0;
      const double *z,*w0;
      int    nquad0 = m_base[0]->GetPointsOrder();
      BstShrDArray tmp  = GetDoubleTmpSpace(nquad0);
    
      PointsManager()[m_base[0].m_basiskey.m_pointskey]->GetZW(z,w0);

      // multiply by integration constants 
      Vmath::Vmul(nquad0,(double*)inarray,1,(double*)w0,1,tmp.get(),1);
    
      Int = Vmath::Vsum(nquad0,tmp.get(),1);
      
      return Int;
    }

  void StdSegExp::IProductWRTBase(const double *base, const double * inarray, 
				 double * outarray, int coll_check)
  {
    int    nquad = m_base[0]->GetPointsOrder();
    const double *z,*w;
    BstShrDArray tmp  = GetDoubleTmpSpace(nquad);

    
    ASSERTL2(m_base[0]->GetAlpha() == 0.0,"Basis[0] has non-zero alpha weight");
    ASSERTL2(m_base[0]->GetBeta() == 0.0,"Basis[1] has non-zero beta weight");

    PointsManager()[m_base[0].m_basiskey..m_pointskey]->GetZW(z,w0);

    Vmath::Vmul(nquad,(double*)inarray,1,(double*)w,1,tmp.get(),1);

    if(coll_check&&m_base[0]->Collocation())
    {
      Vmath::Vcopy(nquad,tmp.get(),1,outarray,1);
    }
    else
    {
      Blas::Dgemv('T',nquad,m_base[0]->GetNumBasis(),1.0,base,nquad,
		  tmp.get(),1,0.0,outarray,1);
    }
    
  }
  
  void StdSegExp::IProductWRTBase(const double * inarray, double * outarray)
  {
    IProductWRTBase(m_base[0]->GetBdata(),inarray,outarray,1);
  }

  void StdSegExp::FillMode(const int mode, double *outarray)
  {
    int   nquad = m_base[0]->GetPointsOrder();
    const double * base  = m_base[0]->GetBdata();

    ASSERTL2(modes <= m_ncoeffs , 
	     "calling argument mode is larger than total expansion order");

    Vmath::Vcopy(nquad,(double *)base+mode*nquad,1, outarray,1);
  }
    
  void StdSegExp::GenMassMatrix(double * outarray)
  {
    StdExpansion::GenerateMassMatrix(outarray);

    // For Fourier basis set the imaginary component of mean mode
    // to have a unit diagonal component in mass matrix 
    if(m_base[0]->GetBasisType() == eFourier)
    {
      outarray[m_base[0]->GetNumBasis()+1] = 1.0;
    }
  }
 
  void StdSegExp::GenLapMatrix(double * outarray)
  {
    int    i;
    int   nquad = m_base[0]->GetPointsOrder();
    const double * dbase  = m_base[0]->GetDbdata();
    const double *z,*w;
    BstShrDArray tmp  = GetDoubleTmpSpace(nquad);
    
    BasisManagerSingleton::Instance().GetZW(m_base[0],z,w);
    
    for(i = 0; i < m_base[0]->GetNumBasis(); ++i)
    {
      Vmath::Vcopy(nquad,(double *)dbase+i*nquad,1, tmp.get(),1);
      IProductWRTBase(m_base[0]->GetDbdata(), tmp.get(),
		      outarray+i*m_base[0]->GetNumBasis(),0);
    }

  }

  StdMatContainer * StdSegExp::GetMassMatrix() 
  {
    StdMatContainer * mat;
    mat = s_elmtmats.GetLocalMass(this);
    return mat;
  }

  StdMatContainer * StdSegExp::GetLapMatrix() 
  {
    StdMatContainer * mat;
    mat = s_elmtmats.GetLocalLap(this);
    return mat;
  }

  //----------------------------
  // Differentiation Methods
  //-----------------------------

  inline void StdSegExp::Deriv(double * outarray)
  {
    TensorDeriv(outarray);
  }  

  void StdSegExp::Deriv(const double *inarray, double * outarray)
  {
    TensorDeriv(inarray,outarray);
  }

  //----------------------------
  // Evaluation Methods
  //----------------------------
   
    void StdSegExp::BwdTrans(double * outarray)
    {
    int           nquad = m_base[0]->GetPointsOrder();
    const double *base  = m_base[0]->GetBdata();
    
    if(m_base[0]->Collocation())
    {
      Vmath::Vcopy(nquad,m_coeffs,1,outarray,1);
    }
    else
    {
      Blas::Dgemv('N',nquad,m_base[0]->GetNumBasis(),1.0,base,nquad,
		  m_coeffs,1,0.0,outarray,1);
    }
  }

  void StdSegExp::FwdTrans(const double *inarray)
  {
    StdMatContainer *M;

    if(m_base[0]->Collocation())
    {
      Vmath::Vcopy(GetNcoeffs(),inarray,1,m_coeffs,1);
    }
    else{
      IProductWRTBase(inarray,m_coeffs);
      M = GetMassMatrix();
      M->Solve(m_coeffs,1);
    }
  }
 
  double StdSegExp::Evaluate(const double *Lcoord)
  {
    return PhysEvaluate(Lcoord);
  }
    
    
  void StdSegExp::MapTo(EdgeOrientation dir, StdExpMap& Map)
  {

    if(Map.GetLen() < 2)
    {
      Map.SetMapMemory(2);
    }
    
    switch(m_base[0]->GetBasisType())
    {
    case eGLL_Lagrange:
    {
      int order = m_base[0]->GetNumBasis();
      if(dir == eForwards)
      {
	Map[0] = 0;
	Map[1] = order-1;
      }
      else
      {
	Map[0] = order-1;
	Map[1] = 0;
      }
    }
    break;
    case eModified_A:
      
      if(dir == eForwards)
      {
	Map[0] = 0;
	Map[1] = 1;
      }
      else
      {
	Map[0] = 1;
	Map[1] = 0;
      }
      break;
    default:
      ASSERTL0(0,"Mapping not defined for this expansion");
    }
  }    
  
  void StdSegExp::SetInvInfo(StdMatContainer *mat, MatrixType Mform)
  {
      mat->SetLda(m_ncoeffs);
      mat->SetMatForm(eSymmetric_Positive);
      
      if(GeoFacType() == eRegular)
      {
	  switch(Mform)
	  {
	  case eMassMatrix:
	    {
		// default setting  for this matrix 
		mat->SetMatForm(eSymmetric_Positive);
		
		switch(m_base[0]->GetBasisType())
		{
		case eOrtho_A: case eLegendre: 
		if(m_base[0]->ExactIprodInt()) // diagonal matrix 
		{
		    mat->SetMatForm(eSymmetric_Positive_Banded);
		    mat->SetBwidth(1);
		}	
		break;
		case eModified_A:
		    // Banded matrix. Note only makes sense to use banded storage
		    // when rank > 2*bandwidth-1
		    if((m_base[0]->ExactIprodInt())&&(m_base[0]->GetNumBasis()>7))
		    { 
			mat->SetMatForm(eSymmetric_Positive_Banded);
			mat->SetBwidth(4);
		    }  
		    break;
		case eFourier:
		    mat->SetMatForm(eSymmetric_Positive_Banded);
		    mat->SetBwidth(1);
		    break;	
		case eGLL_Lagrange:
		    // diagonal matrix 
		    if(m_base[0]->GetPointsOrder() == m_base[0]->GetNumBasis()) 
		    {
			mat->SetMatForm(eSymmetric_Positive_Banded);
			mat->SetBwidth(1);
		    }
		    break;
		default:
		    break;
		}
		break;
	    case eLapMatrix:
		mat->SetMatForm(eSymmetric);	
		
		break;
	    default:
		ASSERTL0(false, "MatrixType not known");
                   break;
	    
	    }
	    
	    
	    }
	}
  }
  
  /** \brief Define an integer mapping vector to map to different types
   */
  
  //  void StdSegExp::Mapto(int *map,  StdExpansion2D  &2DExp, int edge, enum Dir){
  //  
  }//end namespace
}//end namespace

/** 
 * $Log: StdSegExp.cpp,v $
 * Revision 1.7  2007/01/15 15:07:20  pvos
 * updating doxygen documentation
 *
 * Revision 1.6  2006/12/10 19:00:54  sherwin
 * Modifications to handle nodal expansions
 *
 * Revision 1.5  2006/08/16 23:34:42  jfrazier
 * *** empty log message ***
 *
 * Revision 1.4  2006/08/05 19:03:48  sherwin
 * Update to make the multiregions 2D expansion in connected regions work
 *
 * Revision 1.3  2006/06/06 15:25:21  jfrazier
 * Removed unreferenced variables and replaced ASSERTL0(false, ....) with
 * NEKERROR.
 *
 * Revision 1.2  2006/06/01 13:43:20  kirby
 * *** empty log message ***
 *
 * Revision 1.1  2006/05/04 18:58:33  kirby
 * *** empty log message ***
 *
 * Revision 1.52  2006/04/01 21:59:27  sherwin
 * Sorted new definition of ASSERT
 *
 * Revision 1.51  2006/03/21 09:21:32  sherwin
 * Introduced NekMemoryManager
 *
 * Revision 1.50  2006/03/05 22:11:03  sherwin
 *
 * Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
 *
 * Revision 1.49  2006/03/01 08:25:04  sherwin
 *
 * First compiling version of StdRegions
 *
 * Revision 1.48  2006/02/26 23:37:30  sherwin
 *
 * Updates and compiling checks upto StdExpansions1D
 *
 **/ 




