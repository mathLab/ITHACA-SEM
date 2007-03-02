///////////////////////////////////////////////////////////////////////////////
//
// File SegExp.cpp
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
// Description: SegExp routines
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/LocalRegions.hpp>
#include <LocalRegions/SegExp.h>

namespace Nektar
{
    namespace LocalRegions 
    {
	
	SegExp::SegExp(const LibUtilities::BasisKey &Ba, 
		       SpatialDomains::SegGeomSharedPtr geom): 
	    StdRegions::StdSegExp(Ba)
	{
	    m_geom = geom;    


            m_matrixManager.RegisterCreator(MatrixKey(StdRegions::eMassMatrix,
		StdRegions::eNoShapeType,*this),
	        boost::bind(&SegExp::CreateMatrix, this, _1));

            m_linSysManager.RegisterCreator(LinSysKey(StdRegions::eMassMatrix,
		StdRegions::eNoShapeType,*this),
		boost::bind(&SegExp::CreateLinSys, this, _1));
	}
	
  
	SegExp::SegExp(const SegExp &S):StdRegions::StdSegExp(S)
	{
	    m_geom   = S.m_geom;
	    m_minfo  = S.m_minfo;
	}
	
	// by default the StdExpansion1D destructor will be called
	
	SegExp::~SegExp()
	{
	}

	// interpolate and possibly generate geometric factors. 
	MetricRelatedInfoSharedPtr SegExp::GenMinfo()
	{
	    MetricRelatedInfoSharedPtr minfo;
	    SpatialDomains::GeoFacSharedPtr Xgfac;
	    double *ndata;
	    const double **gmat, *odata;
	    
	    if((Xgfac = m_geom->GetXGeoFac()).get() == NULL) //define geometric version 
	    {
		Xgfac = m_geom->GenXGeoFac();
		m_geom->SetXGeoFac(Xgfac);
	    }
	    
	    int coordim = m_geom->GetCoordim();
	    StdRegions::GeomType gtype = Xgfac->GetGtype();
	    
	    minfo.reset(new MetricRelatedInfo (gtype,1,coordim));
	    
	    // interp gfac to 
	    if(gtype == StdRegions::eDeformed)
	    {
		LibUtilities::BasisSharedPtr CBasis0;
		
		CBasis0 = m_geom->GetBasis(0,0); // this assumes all geombasis are same
		
		// basis are different distributions
		if(!(m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey())))
		{
		    int i, nq = m_base[0]->GetNumPoints();
		    
		    // interpolate Geometric data
		    ndata = MemoryManager::AllocateArray<double>(coordim*nq);	
		    gmat  = Xgfac->GetGmat();
		    
		    for(i = 0; i < 2*coordim; ++i)
		    {
			Interp1D(CBasis0->GetBasisKey(), gmat[i], 
				 m_base[0]->GetBasisKey(), 
				 ndata + i*nq);
		    }
		    
		    minfo->ResetGmat(ndata,nq,1,coordim);
		    
		    // interpolate Jacobian
		    ndata = MemoryManager::AllocateArray<double>(nq);	
		    odata = Xgfac->GetJac();
		    
		    Interp1D(CBasis0->GetBasisKey(),odata,
			     m_base[0]->GetBasisKey(), ndata);
		    
		    minfo->ResetJac(ndata);
		    
		    ErrorUtil::Error(ErrorUtil::ewarning,__FILE__,__LINE__,
			 "Need to check/debug routine for deformed elements");
		}
	    }
	    else    // regular geometry
	    {
		// interpolate Geometric data
		ndata = new double [coordim];	
		gmat  = Xgfac->GetGmat();
		
		Blas::Dcopy(coordim,gmat[0],1,ndata,1);
		minfo->ResetGmat(ndata,1,1,coordim);
		
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
	
	- \a inarray: definition of function to be returned at
	quadrature point of expansion.
	
	Outputs:\n
	
	- returns \f$\int^1_{-1} u(\xi_1)d \xi_1 \f$ where \f$inarray[i]
        = u(\xi_{1i}) \f$
	*/
	
	double SegExp::Integral(const double *inarray)
	{
	    
	    int    nquad0 = m_base[0]->GetNumPoints();
	    const double *jac = m_minfo->GetJac();
	    double  ival;
	    BstShrDArray tmp   = GetDoubleTmpSpace(nquad0);
	    
	    // multiply inarray with Jacobian
	    
	    if(m_minfo->GetGtype() == StdRegions::eDeformed)
	    {
		Vmath::Vmul(nquad0,jac,1,(double*)inarray,1,&tmp[0],1);
	    }
	    else
	    {
		Vmath::Smul(nquad0,(double)jac[0],(double*)inarray,1,&tmp[0],1);
	    }
	    
	    // call StdSegExp version;
	    ival = StdSegExp::Integral(&tmp[0]);
	    
	    return ival;
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
	
	
	void SegExp::IProductWRTBase(const double *base, const double *inarray,
				     double * outarray, int coll_check)
	{
	    int    nquad0 = m_base[0]->GetNumPoints();
	    const double *jac = m_minfo->GetJac();
	    BstShrDArray tmp   = GetDoubleTmpSpace(nquad0);
	    
	    
	    // multiply inarray with Jacobian
	    
	    if(m_minfo->GetGtype() == StdRegions::eDeformed)
	    {
		Vmath::Vmul(nquad0,jac,1,(double*)inarray,1,&tmp[0],1);
	    }
	    else
	    {
		Vmath::Smul(nquad0,jac[0],(double*)inarray,1,&tmp[0],1);
	    }
	    
	    StdSegExp::IProductWRTBase(base,&tmp[0],outarray,coll_check);
	}
	
	
	/** \brief  Inner product of \a inarray over region with respect to 
	    the expansion basis (this)->_Base[0] and return in \a outarray 
	    
	    Wrapper call to SegExp::IProduct_WRT_B
	    
	    Input:\n
	    
	    - \a inarray: array of function evaluated at the physical
	    collocation points
	    
	    Output:\n
	    
	    - \a outarray: array of inner product with respect to each
	    basis over region
	*/
	
	
	void SegExp::IProductWRTBase(const double * inarray, double * outarray)
	{
	    IProductWRTBase(m_base[0]->GetBdata(),inarray,outarray,1);
	}
	
	
	//----------------------------
	// Differentiation Methods
	//-----------------------------
	
	/** \brief Evaluate the derivative \f$ d/d{\xi_1} \f$ at the
	    physical quadrature opoints in the expansion (i.e. (this)->_phys)
	    and return in \a outarray. 
	    
	    This is a wrapper function around SegExp::Deriv(dim,inarray,outarray)
	    
	    Input:\n
	    
	    - \a n: number of derivatives to be evaluated where \f$ n \leq  dim\f$
	
	    - \a (this)->_phys: array of function evaluated at the
	    quadrature points
	    
	    Output: \n
	    
	    - \a outarray: array of the derivatives \f$
	    du/d_{\xi_1}|_{\xi_{1i}} d\xi_1/dx, 
	    du/d_{\xi_1}|_{\xi_{1i}} d\xi_1/dy, 
	    du/d_{\xi_1}|_{\xi_{1i}} d\xi_1/dz, 
	    \f$ depending on value of \a dim
	*/
	
	
	void SegExp::PhysDeriv(const int n, double **outarray)
	{
	    PhysDeriv(n,&m_phys[0],outarray);
	}
	
	/** \brief Evaluate the derivative \f$ d/d{\xi_1} \f$ at the
	    physical quadrature points given by \a inarray and return in \a
	    outarray.
	    
	    This is a wrapper around StdExpansion1D::Tensor_Deriv
	    
	    Input:\n
	    
	    - \a n: number of derivatives to be evaluated where \f$ n \leq  dim\f$
	    
	    - \a inarray: array of function evaluated at the quadrature points
	    
	    Output: \n
	    
	    - \a outarray: array of the derivatives \f$
	    du/d_{\xi_1}|_{\xi_{1i}} d\xi_1/dx, 
	    du/d_{\xi_1}|_{\xi_{1i}} d\xi_1/dy, 
	    du/d_{\xi_1}|_{\xi_{1i}} d\xi_1/dz, 
	    \f$ depending on value of \a dim
	*/
	
	
	void SegExp::PhysDeriv(const int n, const double *inarray, double ** outarray)
	{
	    int    i;
	    int    nquad0 = m_base[0]->GetNumPoints();
	    const double **gmat = m_minfo->GetGmat();
	    
	    if(m_geom)
	    {
		ASSERTL2(n <= m_geom->GetCoordim(),
			 "value of n is larger than the number of coordinates");
	    }
	    
	    StdExpansion1D::PhysTensorDeriv(inarray,outarray[n-1]);
	    
	    if(m_minfo->GetGtype() == StdRegions::eDeformed)
	    {
		for(i = 0; i < n; ++i)
		{
		    Vmath::Vmul(nquad0,gmat[i],1,outarray[n-1],1,outarray[i],1);
		}
	    }
	    else
	    {
		for(i = 0; i < n; ++i)
		{
		    Vmath::Smul(nquad0,gmat[i][0],outarray[n-1],1,outarray[i],1);
		}
	    }
	}
	
	//----------------------------
	// Evaluation Methods
	//----------------------------
	
	
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
	    
	    - (this)->m_coeffs: updated array of expansion coefficients. 
	    
	*/ 
	
	
	// need to sort out family of matrices 
	void SegExp::FwdTrans(const double *inarray)
	{
	    if(m_base[0]->Collocation())
	    {
		Vmath::Vcopy(GetNcoeffs(),inarray,1,&m_coeffs[0],1);
	    }
	    else
	    {
		IProductWRTBase(inarray,&m_coeffs[0]);

		LinSysKey   masskey(StdRegions::eMassMatrix,DetShapeType(),
				    *this);
		DNekLinSysSharedPtr matsys = m_linSysManager[masskey];

		DNekVec   v(m_ncoeffs,m_coeffs,eWrapper);
		matsys->Solve(v,v);
	    }
	}
	
	void SegExp::GetCoords(double **coords)
	{
	    int  i;
	    BstShrDArray x;
	    DNekMatSharedPtr   I;
	    LibUtilities::BasisSharedPtr CBasis;
	    
	    ASSERTL0(m_geom, "m_geom not define");
	    
	    // get physical points defined in Geom
	    m_geom->FillGeom();
	    
	    for(i = 0; i < m_geom->GetCoordim(); ++i)
	    {
		CBasis = m_geom->GetBasis(i,0);
		
		if(m_base[0]->GetBasisKey().SamePoints(CBasis->GetBasisKey()))
		{
		    x = m_geom->GetPhys(i);
		    Blas::Dcopy(m_base[0]->GetNumPoints(),&x[0],1,coords[i],1);
		}
		else // Interpolate to Expansion point distribution
		{
		    I = LibUtilities::PointsManager()[CBasis->GetPointsKey()]
			->GetI(m_base[0]->GetPointsKey());
		    
		    for(i = 0; i < m_geom->GetCoordim(); ++i)
		    {
			x = m_geom->GetPhys(i);
			Blas::Dgemv('T',CBasis->GetNumPoints(),
				    m_base[0]->GetNumPoints(),1.0,
				    &((*I).GetPtr())[0],
				    CBasis->GetNumPoints(),&x[0],1,
				    0.0,coords[i],1);
		    }
		}
	    }
	}
	
	void SegExp::GetCoord(const double *Lcoords, double *coords)
	{
	    int  i;
	    
	    ASSERTL1(Lcoords[0] >= -1.0&& Lcoords[1] <= 1.0,
		     "Local coordinates are not in region [-1,1]");
	    
	    m_geom->FillGeom();
	    
	    for(i = 0; i < m_geom->GetCoordim(); ++i)
	    {
		coords[i] = m_geom->GetCoord(i,Lcoords);
	    }
	}
    
	
	void SegExp::WriteToFile(FILE *outfile)
	{
	    int i,j;
	    double *coords[3];
	    int  nquad = m_base[0]->GetNumPoints();
	    
	    ASSERTL0(m_geom,"_geom not defined");
	    
	    int  coordim  = m_geom->GetCoordim();
	    
	    
	    coords[0] = new double [3*nquad];
	    coords[1] = coords[0] + nquad;
	    coords[2] = coords[1] + nquad;
	    
	    GetCoords(coords);
	    
	    fprintf(outfile,"Variables = x");
	    if(coordim == 2)
	    {
		fprintf(outfile,", y");
	    }
	    else if (coordim == 3)
	    {
		fprintf(outfile,", y, z");
	    }
	    fprintf(outfile,", v\n");
	    
	    fprintf(outfile,"Zone, I=%d, F=Point\n",nquad);
	    
	    for(i = 0; i < nquad; ++i)
	    {
		for(j = 0; j < coordim; ++j)
		{
		    fprintf(outfile,"%lf ",coords[j][i]);
		}
		fprintf(outfile,"%lf \n",m_phys[i]);
	    }
	    
	    delete[] coords[0];
	}
	
	void SegExp::WriteToFile(std::ofstream &outfile, const int dumpVar)
	{
	    int i,j;
	    double *coords[3];
	    int     nquad = m_base[0]->GetNumPoints();
	    
	    ASSERTL0(m_geom,"m_geom not defined");
	    
	    int     coordim  = m_geom->GetCoordim();
	    
	    
	    BstShrDArray tmp   = GetDoubleTmpSpace(3*nquad);
	    
	    
	    coords[0] = &tmp[0];
	    coords[1] = coords[0] + nquad;
	    coords[2] = coords[1] + nquad;
	    
	    GetCoords(coords);
	    
	    if(dumpVar)
	    {
		outfile << "Variables = x";
		
		if(coordim == 2)
		{
		    outfile << ", y";
		}
		else if (coordim == 3)
		{
		    outfile << ", y, z";
		}
		outfile << ", v\n" << std::endl;
	    }
	    
	    outfile << "Zone, I=" << nquad <<", F=Point" << std::endl;
	    
	    for(i = 0; i < nquad; ++i)
	    {
		for(j = 0; j < coordim; ++j)
		{
		    outfile << coords[j][i] << " ";
		}
		outfile << m_phys[i] << std::endl;
	    }
	}
	
	double SegExp::Evaluate(const double *coord)
	{
	    double val;
	    double Lcoord;
	    
	    ASSERTL0(m_geom,"_geom not defined");
	    m_geom->GetLocCoords(&Lcoord,coord);
	    
	    return val = StdSegExp::Evaluate(&Lcoord);
	}
	
        DNekMatSharedPtr SegExp::CreateMatrix(const MatrixKey &mkey)
        {
            DNekMatSharedPtr returnval;

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eMassMatrix:
                returnval = GenMassMatrix();
                break;
            default:
                NEKERROR(ErrorUtil::efatal, "Matrix creation not defined");
                break;
            }

            return returnval;
        }


        DNekLinSysSharedPtr SegExp::CreateLinSys(const LinSysKey &mkey)
        {
            DNekLinSysSharedPtr returnval;

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eMassMatrix:
                returnval = MemoryManager::AllocateSharedPtr<DNekLinSys> (m_matrixManager[mkey]);
                break;
            default:
                NEKERROR(ErrorUtil::efatal, "Linear System creation not defined");
                break;
            }

            return returnval;
        }


    } // end of namespace    
}//end of namespace

//
// $Log: SegExp.cpp,v $
// Revision 1.6  2006/06/01 15:27:27  sherwin
// Modifications to account for LibUtilities reshuffle
//
// Revision 1.5  2006/06/01 14:15:58  sherwin
// Added typdef of boost wrappers and made GeoFac a boost shared pointer.
//
// Revision 1.4  2006/05/30 14:00:03  sherwin
// Updates to make MultiRegions and its Demos work
//
// Revision 1.3  2006/05/29 17:05:49  sherwin
// Modified to put shared_ptr around geom definitions
//
// Revision 1.2  2006/05/06 20:36:16  sherwin
// Modifications to get LocalRegions/Project1D working
//
// Revision 1.1  2006/05/04 18:58:46  kirby
// *** empty log message ***
//
// Revision 1.38  2006/03/13 19:47:54  sherwin
//
// Fixed bug related to constructor of GeoFac and also makde arguments to GeoFac all consts
//
// Revision 1.37  2006/03/13 18:20:33  sherwin
//
// Fixed error in ResetGmat
//
// Revision 1.36  2006/03/12 21:59:48  sherwin
//
// compiling version of LocalRegions
//
//
