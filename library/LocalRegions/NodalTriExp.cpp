///////////////////////////////////////////////////////////////////////////////
//
// File NodalTriExp.cpp
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
// Description: NodalTriExp routines
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/NodalTriExp.h>

namespace Nektar
{
  namespace LocalRegions 
  {
      
      NodalTriExp::NodalTriExp(const StdRegions::BasisKey &Ba,  
			       const StdRegions::BasisKey &Bb, 
			       StdRegions::NodalBasisType Ntype, 
			       SpatialDomains::TriGeomSharedPtr geom):
	  StdRegions::StdNodalTriExp(Ba,Bb,Ntype)
      {
	  m_geom = geom;
      }
      
      NodalTriExp::NodalTriExp(const StdRegions::BasisKey &Ba, 
			       const StdRegions::BasisKey &Bb, 
			       StdRegions::NodalBasisType Ntype,
			       double *coeffs, double *phys, 
			       SpatialDomains::TriGeomSharedPtr geom):
	  StdRegions::StdNodalTriExp(Ba,Bb,Ntype,coeffs,phys)
      {
	  m_geom = geom;
      }
      
	
      NodalTriExp::NodalTriExp(const NodalTriExp &T):StdRegions::StdNodalTriExp(T)
      {
	  m_geom = T.m_geom;
      }
      
	
      NodalTriExp::~NodalTriExp()
      {
      }
      
      MetricRelatedInfoSharedPtr NodalTriExp::GenGeoFac()
      {
	  MetricRelatedInfoSharedPtr minfo;
	  SpatialDomains::GeoFacSharedPtr Xgfac;
	  const double **gmat, *odata;
	  double *ndata;
	  
	  
	  if((Xgfac = m_geom->GetXGeoFac()).get() == NULL)  // define geometric version 
	  {
	      m_geom->SetXGeoFac(Xgfac = m_geom->GenXGeoFac());
	  }
	  
	  int coordim = m_geom->GetCoordDim();
	  StdRegions::GeomType gtype = Xgfac->GetGtype();
	  
	  minfo.reset(new MetricRelatedInfo (gtype,2,coordim));
	  
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
	
	
      double NodalTriExp::Integral(const double *inarray)
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
	  ival = StdNodalTriExp::Integral(tmp);
	  delete[] tmp;
	  return ival; 
      }
      

      void NodalTriExp::IProductWRTBase(const double *base0, 
					const double *base1,
					const double *inarray, 
					double *outarray)
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
	  
	  StdNodalTriExp::IProductWRTBase(base0,base1,tmp,outarray);
      }

      void NodalTriExp::IProductWRTBase(const double * inarray, double * outarray)
      {
	  IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(),
			  inarray,outarray);
      }
      

      /** \brief Get the mass matrix attached to this expansion by using
	  the StdMatrix manager _ElmtMats and return the standard Matrix
	  container */

      StdRegions::StdMatContainer * NodalTriExp::GetMassMatrix() 
      {
	  StdRegions::StdMatContainer * tmp;
	  tmp = s_elmtmats.GetLocalMass(this);
	  return tmp;
      }
	

	
      /** \brief Get the weak Laplacian matrix attached to this
	  expansion by using the StdMatrix manager _ElmtMats and return
	  the standard Matrix container */
      
      StdRegions::StdMatContainer * NodalTriExp::GetLapMatrix() 
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
      void NodalTriExp::Deriv(const int n, double **outarray)
      {
	  Deriv(n, this->m_phys, outarray);
      }
      

      /** 
	  \brief Calculate the derivative of the physical points 
	  
	  Basic call uses Tensor_Deriv funcion defined under StdExpansion2D.
      **/
      void NodalTriExp::Deriv(const int n, const double *inarray, 
			      double **outarray)
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
	  
	  StdNodalTriExp::Deriv(inarray, tmp, outarray[n-1]);
	  
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
      
      
      void NodalTriExp::FwdTrans(const double *inarray)
      {
	  StdRegions::StdMatContainer *M;
	  IProductWRTBase(inarray,m_coeffs);
	  M = GetMassMatrix();
	  M->Solve(m_coeffs,1);
      }
	
	
      void NodalTriExp::GetCoords(double **coords)
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
		  Blas::Dcopy(m_base[0]->GetPointsOrder()*
			      m_base[1]->GetPointsOrder(),
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
      void NodalTriExp::GetCoord(const double *Lcoords, double *coords)
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
      
      
      void NodalTriExp::WriteToFile(FILE *outfile)
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
      
      void NodalTriExp::WriteToFile(std::ofstream &outfile, const int dumpVar)
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
	  
	  if(dumpVar)
	  {
	      outfile << "Variables = x";
	      
	      if(gdim == 2)
	      {
		  outfile << ", y";
	      }
	      else if (gdim == 3)
	      {
		  outfile << ", y, z";
	      }
	      outfile << ", v\n" << std::endl;
	  }
	  
	  outfile << "Zone, I=" << nquad0 << ", J=" << 
	      nquad1 <<", F=Point" << std::endl;
	  
	  for(i = 0; i < nquad0*nquad1; ++i)
	  {
	      for(j = 0; j < gdim; ++j)
	      {
		  outfile << coords[j][i] << " ";
	      }
	      outfile << m_phys[i] << std::endl;
	  }
	  delete[] coords[0];
      }
      
      
      double NodalTriExp::Evaluate(const double *coord)
      {
	  double Lcoord[2];
	  
	  ASSERTL0(m_geom,"_geom not defined");
	  m_geom->GetLocCoords(Lcoord,coord);
	  
	  return StdNodalTriExp::Evaluate(Lcoord);
      }

  }//end of namespace
}//end of namespace

/** 
 *    $Log: NodalTriExp.cpp,v $
 *    Revision 1.1  2006/05/04 18:58:45  kirby
 *    *** empty log message ***
 *
 *    Revision 1.3  2006/03/12 07:43:32  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/
