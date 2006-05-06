///////////////////////////////////////////////////////////////////////////////
//
// File $Source: /usr/sci/projects/Nektar/cvs/Nektar++/libs/LocalRegions/TriExp.cpp,v $ 
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
// Description: 
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/TriExp.h>
#include <stdio.h>

namespace Nektar
{
  namespace LocalRegions 
  {
  
    TriExp::TriExp(const StdRegions::BasisKey &Ba,  
		   const StdRegions::BasisKey &Bb, 
		   SpatialDomains::TriGeom *geom):
      StdRegions::StdTriExp(Ba,Bb)
    {
      m_geom = geom;
    }
  
    TriExp::TriExp(const StdRegions::BasisKey &Ba, 
		   const StdRegions::BasisKey &Bb, 
		   double *coeffs, double *phys, SpatialDomains::TriGeom *geom)
      :StdRegions::StdTriExp(Ba,Bb,coeffs,phys)
    {
      m_geom = geom;
    }
  
 
    TriExp::TriExp(const TriExp &T):StdRegions::StdTriExp(T)
    {
      m_geom = T.m_geom;
    }
  
 
    TriExp::~TriExp()
    {
    }

  
    MetricRelatedInfo *TriExp::GenGeoFac()
    {
      MetricRelatedInfo *minfo;
      SpatialDomains::GeoFac *Xgfac;
      const double **gmat, *odata;
      double *ndata;
      
    
      if((Xgfac = m_geom->GetXGeoFac()) == NULL)  // define geometric version 
      {
	m_geom->SetXGeoFac(Xgfac = m_geom->GenXGeoFac());
      }

      int coordim = m_geom->GetCoordDim();
      StdRegions::GeomType gtype = Xgfac->GetGtype();
    
      minfo = new MetricRelatedInfo (gtype,2,coordim);
    
      // interp gfac to 
      if(gtype == StdRegions::eDeformed)
      {
	const StdRegions::BasisKey *CBasis0;
	const StdRegions::BasisKey *CBasis1;
	int coordim = m_geom->GetCoordDim();

	CBasis0 = m_geom->GetBasis(0,0); // this assumes all goembasis are same
	CBasis1 = m_geom->GetBasis(0,1);
      
	// basis are different distributions
	if(!(m_base[0]->SamePoints(*CBasis0))||
	   !(m_base[1]->SamePoints(*CBasis1)))
	{
	  int i;
	  int nq = m_base[0]->GetPointsOrder()*m_base[1]->GetPointsOrder();
	  double *ndata;
	  const double **gmat, *odata;
	  
	  // interpolate Geometric data
	  ndata = new double [2*coordim*nq];	
	  gmat  = Xgfac->GetGmat();
	
	  for(i = 0; i < 2*coordim; ++i)
	  {
	    Interp2D(CBasis0,CBasis1, gmat[i], m_base[0],
		     m_base[1], ndata + i*nq);
	  }
	
	  minfo->ResetGmat(ndata,nq,2,coordim);
	  
	  // interpolate Jacobian
	  ndata = new double [nq];	
	  odata = Xgfac->GetJac();
	  
	  Interp2D(CBasis0,CBasis1,odata,m_base[0],m_base[1],ndata);
	  
	  minfo->ResetJac(ndata);
	}
      }
      else
      {
	// interpolate Geometric data
	ndata = new double [2*coordim];	
	gmat  = Xgfac->GetGmat();

	Blas::Dcopy(2*coordim,gmat[0],1,ndata,1);
	minfo->ResetGmat(ndata,2,2,coordim);
      
	// interpolate Jacobian
	ndata = new double [1];	
	odata = Xgfac->GetJac();
	
	ndata[0] = odata[0];
	minfo->ResetJac(ndata);
      }
      
      return minfo;
    }


    //----------------------------
    // Integration Methods
    //----------------------------
    
    /** \brief Integrate the physical point list \a inarray over region
	and return the value
	
	Inputs:\n
	
	- \a inarray: definition of function to be returned at quadrature point 
	of expansion. 
	
	Outputs:\n
	
	- returns \f$\int^1_{-1}\int^1_{-1} u(\xi_1, \xi_2) J[i,j] d
        \xi_1 d \xi_2 \f$ where \f$inarray[i,j] = u(\xi_{1i},\xi_{2j})
        \f$ and \f$ J[i,j] \f$ is the Jacobian evaluated at the
        quadrature point.
    */
    
    
    double TriExp::Integral(const double *inarray)
    {
      
      int    nquad0 = m_base[0]->GetPointsOrder();
      int    nquad1 = m_base[1]->GetPointsOrder();
      const double *jac = m_minfo->GetJac();
      double ival, *tmp = new double [nquad0*nquad1];

      // multiply inarray with Jacobian
      if(m_minfo->GetGtype() == StdRegions::eDeformed)
      {
	Vmath::Vmul(nquad0*nquad1,jac,1,(double*)inarray,1,tmp,1);
      }
      else
      {
	Vmath::Smul(nquad0*nquad1,(double) jac[0],(double*)inarray,1,tmp,1);
      }
    
      // call StdQuadExp version;
      ival = StdTriExp::Integral(tmp);
      delete[] tmp;
      return ival; 
    }
    
    
    /** 
	\brief Calculate the inner product of inarray with respect to
	the basis B=base0*base1 and put into outarray:
	
	\f$ 
	\begin{array}{rcl}
	I_{pq} = (\phi_q \phi_q, u) & = & \sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1}
	\phi_p(\xi_{0,i}) \phi_q(\xi_{1,j}) w^0_i w^1_j u(\xi_{0,i} \xi_{1,j}) 
	J_{i,j}\\
	& = & \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i})
	\sum_{j=0}^{nq_1} \phi_q(\xi_{1,j}) \tilde{u}_{i,j}  J_{i,j}
	\end{array}
	\f$ 
	
	where
	
	\f$  \tilde{u}_{i,j} = w^0_i w^1_j u(\xi_{0,i},\xi_{1,j}) \f$
	
	which can be implemented as
	
	\f$  f_{qi} = \sum_{j=0}^{nq_1} \phi_q(\xi_{1,j}) \tilde{u}_{i,j} = 
	{\bf B_1 U}  \f$
	\f$  I_{pq} = \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i}) f_{qi} = 
	{\bf B_0 F}  \f$
    **/
    
    void TriExp::IProductWRTBase(const double *base0, const double *base1,
				 const double * inarray, double * outarray)
    {

      int    nquad0 = m_base[0]->GetPointsOrder();
      int    nquad1 = m_base[1]->GetPointsOrder();
      const double *jac = m_minfo->GetJac();
      double *tmp = new double [nquad0*nquad1];
      
      // multiply inarray with Jacobian
      if(m_minfo->GetGtype() == StdRegions::eDeformed)
      {
	Vmath::Vmul(nquad0*nquad1,jac,1,(double*)inarray,1,tmp,1);
      }
      else
      {
	Vmath::Smul(nquad0*nquad1,jac[0],(double*)inarray,1,tmp,1);
      }
    
      StdTriExp::IProductWRTBase(base0,base1,tmp,outarray);
    }

    /** \brief  Inner product of \a inarray over region with respect to the 
	expansion basis (this)->_Base[0] and return in \a outarray 
	
	Wrapper call to TriExp::IProductWRTBase
	
	Input:\n
	
	- \a inarray: array of function evaluated at the physical
	collocation points
	
	Output:\n
	
	- \a outarray: array of inner product with respect to each
	basis over region
	
    */
    
    
    void TriExp::IProductWRTBase(const double * inarray, double * outarray)
    {
      IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(),
		      inarray,outarray);
    }

    /** \brief Get the mass matrix attached to this expansion by using
	the StdMatrix manager _ElmtMats and return the standard Matrix
	container */
    
    
    StdRegions::StdMatContainer * TriExp::GetMassMatrix() 
    {
      StdRegions::StdMatContainer * tmp;
      tmp = s_elmtmats.GetLocalMass(this);
      return tmp;
    }
  
    /** \brief Get the weak Laplacian matrix attached to this
	expansion by using the StdMatrix manager _ElmtMats and return
	the standard Matrix container */
    
    StdRegions::StdMatContainer * TriExp::GetLapMatrix() 
    {
      StdRegions::StdMatContainer * tmp;
      tmp = s_elmtmats.GetLocalLap(this);
      return tmp;
    }


    ///////////////////////////////
    /// Differentiation Methods
    ///////////////////////////////
    
    /** 
	\brief Calculate the deritive of the physical points 
    **/
    void TriExp::Deriv(const int n, double **outarray)
    {
      Deriv(n, this->m_phys, outarray);
    }
    

    /** 
	\brief Calculate the derivative of the physical points 
	
	Basic call uses Tensor_Deriv funcion defined under StdExpansion2D.
    **/
    void TriExp::Deriv(const int n, const double *inarray, double **outarray)
    {
      int    i;
      int    nquad0 = m_base[0]->GetPointsOrder();
      int    nquad1 = m_base[0]->GetPointsOrder();
      const double **gmat = m_minfo->GetGmat();
      double *tmp = new double [nquad0*nquad1];
    
      if(m_geom)
      {
	ASSERTL2(n <= m_geom->GeCoordDim(),
		 "value of n is larger than the number of coordinates");
      }
      
      StdTriExp::Deriv(inarray, tmp, outarray[n-1]);
     
      if(m_minfo->GetGtype() == StdRegions::eDeformed)
      {
	for(i = 0; i < n; ++i)
	{
	  Vmath::Vmul  (nquad0*nquad1,gmat[2*i+1],1,outarray[n-1],1,
			outarray[i],1);
	  Vmath::Vvtvp (nquad0*nquad1,gmat[2*i]  ,1,tmp          ,1, 
			outarray[i],1,outarray[i],1);
	}
      }
      else
      {
	for(i = 0; i < n; ++i)
	{
	  Vmath::Vmul  (nquad0*nquad1,gmat[2*i+1],1,outarray[n-1],1,
			outarray[i],1);
	  Vmath::Vvtvp (nquad0*nquad1,gmat[2*i]  ,1,tmp          ,1,
			outarray[i],1,outarray[i],1);
	}
      }
    }
    

    /** \brief Forward transform from physical quadrature space
	stored in \a inarray and evaluate the expansion coefficients and
	store in \a (this)->m_coeffs  
	
	Inputs:\n
	
	- \a inarray: array of physical quadrature points to be transformed
	
	Outputs:\n
	
      - (this)->_coeffs: updated array of expansion coefficients. 
      
    */ 

  
    void TriExp::FwdTrans(const double *inarray)
    {
      StdRegions::StdMatContainer *M;

      IProductWRTBase(inarray,m_coeffs);
      M = GetMassMatrix();
      M->Solve(m_coeffs,1);
    }
   
    
    void TriExp::GetCoords(double **coords)
    {
      int  i;
      const double *x;
      const StdRegions::BasisKey *CBasis0;
      const StdRegions::BasisKey *CBasis1;

      // get physical points defined in Geom
      m_geom->FillGeom();

      for(i = 0; i < m_geom->GetCoordDim(); ++i)
      {
	CBasis0 = m_geom->GetBasis(i,0); 
	CBasis1 = m_geom->GetBasis(i,1);
      
	if((m_base[0]->SamePoints(*CBasis0))&&
	   (m_base[1]->SamePoints(*CBasis1)))
	{
	  x = m_geom->GetPhys(i);
	  Blas::Dcopy(m_base[0]->GetPointsOrder()*m_base[1]->GetPointsOrder(),
		      x,1,coords[i],1);
	}
	else // Interpolate to Expansion point distribution
	{
	  Interp2D(CBasis0, CBasis1,m_geom->GetPhys(i),
		   m_base[0],m_base[1],coords[i]);
	}
      }

    }
  

    // get the coordinates "coords" at the local coordinates "Lcoords"
    void TriExp::GetCoord(const double *Lcoords, double *coords)
    {
      int  i;
      
      ASSERTL1(Lcoords[0] >= -1.0 && Lcoords[1] <= 1.0 && 
	       Lcoords[1] >= -1.0 && Lcoords[1]  <=1.0,
	       "Local coordinates are not in region [-1,1]");
      
      m_geom->FillGeom();

      for(i = 0; i < m_geom->GetCoordDim(); ++i)
      {
	coords[i] = m_geom->GetCoord(i,Lcoords);
      }
    }

  
    void TriExp::WriteToFile(FILE *outfile)
    {
      int i,j;
      double *coords[3];
      int  nquad0 = m_base[0]->GetPointsOrder();
      int  nquad1 = m_base[1]->GetPointsOrder();
      int  gdim   = m_geom->GetCoordDim();
      

      coords[0] = new double [3*nquad0*nquad1];
      coords[1] = coords[0] + nquad0*nquad1;
      coords[2] = coords[1] + nquad0*nquad1;

      GetCoords(coords);
    
      fprintf(outfile,"Variables = x");

      if(gdim == 2)
      {
	fprintf(outfile,", y");
      }
      else if (gdim == 3)
      {
	fprintf(outfile,", y, z");
      }
      fprintf(outfile,", v\n");
    
      fprintf(outfile,"Zone, I=%d, J=%d, F=Point\n",nquad0,nquad1);
      for(i = 0; i < nquad0*nquad1; ++i)
      {
	for(j = 0; j < gdim; ++j)
	{
	  fprintf(outfile,"%lf ",coords[j][i]);
	}
	fprintf(outfile,"%lf \n",m_phys[i]);
      }
    }
    
    double TriExp::Evaluate(const double *coord)
    {
      double Lcoord[2];
    
      ASSERTL0(m_geom,"_geom not defined");
      m_geom->GetLocCoords(Lcoord,coord);
      
      return StdTriExp::Evaluate(Lcoord);
    }
    
  }//end of namespace
}//end of namespace

/** 
 *    $Log: TriExp.cpp,v $
 *    Revision 1.1  2006/05/04 18:58:46  kirby
 *    *** empty log message ***
 *
 *    Revision 1.17  2006/03/13 19:47:54  sherwin
 *
 *    Fixed bug related to constructor of GeoFac and also makde arguments to GeoFac all consts
 *
 *    Revision 1.16  2006/03/13 18:20:33  sherwin
 *
 *    Fixed error in ResetGmat
 *
 *    Revision 1.15  2006/03/12 21:59:48  sherwin
 *
 *    compiling version of LocalRegions
 *
 *    Revision 1.14  2006/03/12 07:43:32  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/

