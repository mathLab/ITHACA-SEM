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
    
    /// \brief Constructor using BasisKey class for quadrature points and
    /// order definition.
    ///
    /// Inputs:\n
    ///	
    ///	- \a Ba: BasisKey class definition containsing order and
    /// quadrature points.
    
    StdSegExp::StdSegExp(const BasisKey &Ba):
      StdExpansion1D(Ba,Ba.GetBasisOrder(),NULL,NULL,true)
    {    
    }


    /// \brief Constructor using BasisKey class for quadrature points and
    /// order definition where _coeffs and _phys are all set.
    ///      
    /// Inputs:\n
    ///
    ///  - \a Ba: BasisKey definition containsing order and
    /// quadrature points.
    /// - \a coeffs list of expansions coefficient to be set in m_coeffs
    /// - \a phys: list of physical values to be set in m_phys
    StdSegExp::StdSegExp(const BasisKey &Ba,double *coeffs, double *phys):
      StdExpansion1D(Ba,Ba.GetBasisOrder(),coeffs,phys,false)
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
    
    /// \brief Integrate the physical point list \a inarray over region
    /// and return the value
    ///
    /// Inputs:\n
    ///
    /// - \a inarray: definition of function to be returned at quadrature point 
    /// of expansion. 
    ///
    /// Outputs:\n
    ///
    /// - returns \f$\int^1_{-1} u(\xi_1)d \xi_1 \f$ where \f$inarray[i]
    /// = u(\xi_{1i}) \f$
    double StdSegExp::Integral(const double *inarray)
    {
      double Int = 0.0;
      const double *z,*w0;
      int    nquad0 = m_base[0]->GetPointsOrder();
      BstShrDArray tmp  = GetDoubleTmpSpace(nquad0);
    
      BasisManagerSingleton::Instance().GetZW(m_base[0],z,w0);

    // multiply by integration constants 
      Vmath::Vmul(nquad0,(double*)inarray,1,(double*)w0,1,tmp.get(),1);
    
      Int = Vmath::Vsum(nquad0,tmp.get(),1);
      
      return Int;
    }

  /**
     \brief  Inner product of \a inarray over region with respect to
     expansion basis \a base and return in \a outarray 

     Calculate \f$ I[p] = \int^{1}_{-1} \phi_p(\xi_1) u(\xi_1) d\xi_1
     = \sum_{i=0}^{nq-1} \phi_p(\xi_{1i}) u(\xi_{1i}) w_i \f$ where
     \f$ outarray[p] = I[p], inarray[i] = u(\xi_{1i}), base[p*nq+i] =
     \phi_p(\xi_{1i}) \f$.

     Inputs: \n 

     - \a base: an array definiing the local basis for the inner
     product usually passed from Basis->get_bdata() or
     Basis->get_Dbdata()
     - \a inarray: physical point array of function to be integrated
     \f$ u(\xi_1) \f$
     - \a coll_check: Flag to identify when a Basis->collocation()
       call should be performed to see if this is a GLL_Lagrange basis
       with a collocation property. (should be set to 0 if taking the
       inner product with respect to the derivative of basis)
       
     Output: \n

     - \a outarray: array of coefficients representing the inner
       product of function with ever  mode in the exapnsion

  **/
    
  void StdSegExp::IProductWRTBase(const double *base, const double * inarray, 
				 double * outarray, int coll_check)
  {
    int    nquad = m_base[0]->GetPointsOrder();
    const double *z,*w;
    BstShrDArray tmp  = GetDoubleTmpSpace(nquad);

    
    ASSERTL2(m_base[0]->GetAlpha() == 0.0,"Basis[0] has non-zero alpha weight");
    ASSERTL2(m_base[0]->GetBeta() == 0.0,"Basis[1] has non-zero beta weight");

    BasisManagerSingleton::Instance().GetZW(m_base[0],z,w);
    
    Vmath::Vmul(nquad,(double*)inarray,1,(double*)w,1,tmp.get(),1);

    if(coll_check&&m_base[0]->Collocation())
    {
      Vmath::Vcopy(nquad,tmp.get(),1,outarray,1);
    }
    else
    {
      Blas::Dgemv('T',nquad,m_base[0]->GetBasisOrder(),1.0,base,nquad,
		  tmp.get(),1,0.0,outarray,1);
    }
    
  }
  
  /** \brief  Inner product of \a inarray over region with respect to the 
      expansion basis (this)->_Base[0] and return in \a outarray 
      
      Wrapper call to StdSegExp::IProduct_WRT_B
      
      Input:\n
      
      - \a inarray: array of function evaluated at the physical
      collocation points
      
      Output:\n

      - \a outarray: array of inner product with respect to each
      basis over region
  */
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
    
  /** \brief Generate local mass matrix \f$ {\bf M}[i][j] =
      \int^1_{-1} \phi_i(\xi_1) \phi_j(\xi_1) d\xi_1 \f$ in standard
      region and store in \a outarray  

      Input:\n

      None - all values set up in expansion class

      Output:\n
      
      - \a outarray: Local mass matrix in standard region. Matrix is
        in row major format and is of stored as \a
        outarray[j*(this)->_ncoeffs + i]

  */
  void StdSegExp::GenMassMatrix(double * outarray)
  {
    StdExpansion::GenerateMassMatrix(outarray);

    // For Fourier basis set the imaginary component of mean mode
    // to have a unit diagonal component in mass matrix 
    if(m_base[0]->GetBasisType() == eFourier)
    {
      outarray[m_base[0]->GetBasisOrder()+1] = 1.0;
    }
  }

  /** \brief Generate local weak Laplacian matrix \f$ {\bf L}[i][j] =
      \int^1_{-1} \frac{d \phi_i(\xi_1)}{d \xi_1}\frac{d
      \phi_j(\xi_1)}{d \xi_1} d\xi_1 \f$ in standard region and store in
      \a outarray  

      Input:\n

      None - all values set up in expansion class

      Output:\n
      
      - \a outarray: Local mass matrix in standard region. Matrix is
        in row major format and is of stored as \a
        outarray[j*(this)->_ncoeffs + i]
  */

  void StdSegExp::GenLapMatrix(double * outarray)
  {
    int    i;
    int   nquad = m_base[0]->GetPointsOrder();
    const double * dbase  = m_base[0]->GetDbdata();
    const double *z,*w;
    BstShrDArray tmp  = GetDoubleTmpSpace(nquad);
    
    BasisManagerSingleton::Instance().GetZW(m_base[0],z,w);
    
    for(i = 0; i < m_base[0]->GetBasisOrder(); ++i)
    {
      Vmath::Vcopy(nquad,(double *)dbase+i*nquad,1, tmp.get(),1);
      IProductWRTBase(m_base[0]->GetDbdata(), tmp.get(),
		      outarray+i*m_base[0]->GetBasisOrder(),0);
    }

  }

  /** \brief Get the mass matrix attached to this expansion by using
      the StdMatrix manager _ElmtMats and return the standard Matrix
      container */
  StdMatContainer * StdSegExp::GetMassMatrix() 
  {
    StdMatContainer * mat;
    mat = s_elmtmats.GetLocalMass(this);
    return mat;
  }

  /** \brief Get the weak Laplacian matrix attached to this
      expansion by using the StdMatrix manager _ElmtMats and return
      the standard Matrix container */
  StdMatContainer * StdSegExp::GetLapMatrix() 
  {
    StdMatContainer * mat;
    mat = s_elmtmats.GetLocalLap(this);
    return mat;
  }

  //----------------------------
  // Differentiation Methods
  //-----------------------------
  
  /** \brief Evaluate the derivative \f$ d/d{\xi_1} \f$ at the
      physical quadrature opoints in the expansion (i.e. (this)->_phys)
      and return in \a outarray. 
      
      This is a wrapper function around StdExpansion1D::TensorDeriv
      
      Input:\n

      - \a (this)->_phys: array of function evaluated at the
        quadrature points

      Output: \n

      - \a outarray: array of the derivative \f$
        du/d_{\xi_1}|_{\xi_{1i}} \f$
  */
  inline void StdSegExp::Deriv(double * outarray)
  {
    TensorDeriv(outarray);
  }  

  /** \brief Evaluate the derivative \f$ d/d{\xi_1} \f$ at the
      physical quadrature points given by \a inarray and return in \a
      outarray.

      This is a wrapper around StdExpansion1D::Tensor_Deriv

      Input:\n

      - \a inarray: array of function evaluated at the quadrature
        points

      Output: \n

      - \a outarray: array of the derivative \f$
        du/d_{\xi_1}|_{\xi_{1i}} \f$

*/
  void StdSegExp::Deriv(const double *inarray, double * outarray)
  {
    TensorDeriv(inarray,outarray);
  }

  //----------------------------
  // Evaluation Methods
  //----------------------------
  
  
  /** \brief Backward transform from coefficient space (stored in
      (this)->_coeffs) and evaluate at the physical quadrature points \a
      outarray 

      Operation can be evaluated as \f$ u(\xi_{1i}) =
      \sum_{p=0}^{order-1} \hat{u}_p \phi_p(\xi_{1i}) \f$ or
      equivalently \f$ {\bf u} = {\bf B}^T {\bf \hat{u}} \f$ where
      \f${\bf B}[i][j] = \phi_i(\xi_{1j}), \mbox{\_coeffs}[p] = {\bf
      \hat{u}}[p] \f$
      

      Inputs:\n

      - (this)->_coeffs: Expansion coefficients

      Output: \n

      - \a outarray: Array of physical quadrature points of function
  */
  
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
      Blas::Dgemv('N',nquad,m_base[0]->GetBasisOrder(),1.0,base,nquad,
		  m_coeffs,1,0.0,outarray,1);
    }
  }

  
  /** \brief Forward transform from physical quadrature space
      stored in \a inarray and evaluate the expansion coefficients and
      store in \a (this)->_coeffs  

      Perform a forward transform using a Galerkin projection by
      taking the inner product of the physical points and multiplying
      by the inverse of the mass matrix using the Solve method of the
      standard matrix container holding the local mass matrix, i.e.
      \f$ {\bf \hat{u}} = {\bf M}^{-1} {\bf I} \f$ where \f$ {\bf I}[p] =
      \int^1_{-1} \phi_p(\xi_1) u(\xi_1) d\xi_1 \f$

      Inputs:\n
      
      - \a inarray: array of physical quadrature points to be transformed

      Outputs:\n

      - (this)->_coeffs: updated array of expansion coefficients. 

  */ 
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
  
  /** \brief Single Point Evaluation: \f$ u(x) = \sum_p \phi_p(x)
      \hat{u}_p = \sum_p h_p(x) u(x_p)\f$  */
  double StdSegExp::Evaluate(const double *Lcoord)
  {
    return PhysEvaluate(Lcoord);
  }
    
    
  void StdSegExp::MapTo(EdgeOrientation dir, StdExpMap& Map)
  {

    if(Map.GetLen() < 2)
    {
      Map.SetMap(2);
    }
    
    switch(m_base[0]->GetBasisType())
    {
    case eGLL_Lagrange:
    {
      int order = m_base[0]->GetBasisOrder();
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
  
  /** \brief Define an integer mapping vector to map to different types
   */
  
  //  void StdSegExp::Mapto(int *map,  StdExpansion2D  &2DExp, int edge, enum Dir){
  //  
  }//end namespace
}//end namespace

/** 
 * $Log: StdSegExp.cpp,v $
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




