///////////////////////////////////////////////////////////////////////////////
//
// File QuadExp.cpp 
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

#include <LocalRegions/LocalRegions.hpp>
#include <stdio.h>

#include <LocalRegions/QuadExp.h>


namespace Nektar
{
    namespace LocalRegions 
    {
	
	
	QuadExp::QuadExp(const LibUtilities::BasisKey &Ba, 
			 const LibUtilities::BasisKey &Bb, 
			 const SpatialDomains::QuadGeomSharedPtr &geom):
	    StdRegions::StdQuadExp(Ba,Bb)
	{
            m_geom = geom;

            m_matrixManager.RegisterCreator(MatrixKey(StdRegions::eMassMatrix,
                StdRegions::eNoShapeType,*this),
                boost::bind(&QuadExp::CreateMatrix, this, _1));

            m_linSysManager.RegisterCreator(LinSysKey(StdRegions::eMassMatrix,
                StdRegions::eNoShapeType,*this),
                boost::bind(&QuadExp::CreateLinSys, this, _1));

            GenMetricInfo();
	}
	
	
	QuadExp::QuadExp(const QuadExp &T):StdRegions::StdQuadExp(T)
	{
	    m_geom       = T.m_geom;
            m_metricinfo = T.m_metricinfo;
	}
	
	// by default the StdQuadExp destructor will be called
	
	QuadExp::~QuadExp()
	{
	}
	
	void QuadExp::GenMetricInfo()
	{
	    SpatialDomains::GeomFactorsSharedPtr Xgfac;
            
            Xgfac = m_geom->GetGeomFactors();

            if(Xgfac->GetGtype() != SpatialDomains::eDeformed)
            {
                m_metricinfo = Xgfac;
            }
            else
            {
		LibUtilities::BasisSharedPtr CBasis0;
		LibUtilities::BasisSharedPtr CBasis1;

                Array<OneD,NekDouble> ndata;
                ConstArray<TwoD,NekDouble> gmat;
                ConstArray<OneD,NekDouble> odata;
                int coordim = m_geom->GetCoordim();

		CBasis0 = m_geom->GetBasis(0,0); // this assumes all goembasis are same
		CBasis1 = m_geom->GetBasis(0,1);
                
                m_metricinfo = MemoryManager<SpatialDomains::GeomFactors>::AllocateSharedPtr(); 
                
		// basis are different distributions
		if(!(m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))||
		   !(m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey())))
		{
		    int i;
		    int nq = m_base[0]->GetNumPoints()*
                        m_base[1]->GetNumPoints();
                    
                    ndata = Array<OneD,NekDouble>(coordim*nq);
                    gmat  = Xgfac->GetGmat();
                    
                    // interpolate Geometric data
				    
		    for(i = 0; i < 2*coordim; ++i)
		    {
			Interp2D(CBasis0->GetBasisKey(),
                                 CBasis1->GetBasisKey(), 
                                 &gmat[i][0], 
                                 m_base[0]->GetBasisKey(),
				 m_base[1]->GetBasisKey(),
                                 &ndata[0] + i*nq);
		    }

                    m_metricinfo->ResetGmat(ndata,nq,2,coordim);

		    // interpolate Jacobian
                    ndata = Array<OneD,NekDouble>(nq);	
                    odata = Xgfac->GetJac();
		    
		    Interp2D(CBasis0->GetBasisKey(),
                             CBasis1->GetBasisKey(),
                             &odata[0],
                             m_base[0]->GetBasisKey(),
			     m_base[1]->GetBasisKey(), 
                             &ndata[0]);
		    
                    m_metricinfo->ResetJac(nq,ndata);

                    NEKERROR(ErrorUtil::ewarning,
                        "Need to check/debug routine for deformed elements");

		}
                else // Same data can be used 
                {
		    int nq = m_base[0]->GetNumPoints()*
                        m_base[1]->GetNumPoints();
                    
                    // interpolate Geometric data
                    ndata = Array<OneD,NekDouble>(2*nq*coordim); 
                    gmat  = Xgfac->GetGmat();
                    Blas::Dcopy(2*coordim*nq,&gmat[0][0],1,&ndata[0],1);
                    
                    m_metricinfo->ResetGmat(ndata,nq,2,coordim);
		
                    // Copy Jacobian
                    ndata = Array<OneD,NekDouble>(nq);	
                    
                    odata = Xgfac->GetJac();
                    Blas::Dcopy(nq,&odata[0],1,&ndata[0],1);
                    
                    m_metricinfo->ResetJac(nq,ndata);
 
                }
            }
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
            \xi_1 d \xi_2 \f$ where \f$inarray[i,j] =
            u(\xi_{1i},\xi_{2j}) \f$ and \f$ J[i,j] \f$ is the
            Jacobian evaluated at the quadrature point.
	*/
        NekDouble QuadExp::Integral(const ConstArray<OneD,NekDouble> &inarray)
	{
	    int    nquad0 = m_base[0]->GetNumPoints();
	    int    nquad1 = m_base[1]->GetNumPoints();
	    ConstArray<OneD,NekDouble> jac = m_metricinfo->GetJac();
	    NekDouble ival;
            Array<OneD,NekDouble> tmp   = Array<OneD,NekDouble>(nquad0*nquad1);
	    
	    // multiply inarray with Jacobian
	    
	    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
	    {
		Vmath::Vmul(nquad0*nquad1,&jac[0],1,(NekDouble*)&inarray[0],1,
                            &tmp[0],1);
	    }
	    else
	    {
		Vmath::Smul(nquad0*nquad1,(NekDouble) jac[0],
                            (NekDouble*)&inarray[0],1,&tmp[0],1);
	    }
	    
	    // call StdQuadExp version;
	    ival = StdQuadExp::Integral(tmp);
	    
	    return  ival;
	}
	
	/** 
	    \brief Calculate the inner product of inarray with respect to
	    the basis B=base0*base1 and put into outarray:
	    
	    \f$ \begin{array}{rcl} I_{pq} = (\phi_q \phi_q, u) & = &
	    \sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1} \phi_p(\xi_{0,i})
	    \phi_q(\xi_{1,j}) w^0_i w^1_j u(\xi_{0,i} \xi_{1,j})
	    J_{i,j}\\ & = & \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i})
	    \sum_{j=0}^{nq_1} \phi_q(\xi_{1,j}) \tilde{u}_{i,j}
	    J_{i,j} \end{array} \f$
	    
	    where
	    
	    \f$  \tilde{u}_{i,j} = w^0_i w^1_j u(\xi_{0,i},\xi_{1,j}) \f$
	    
	    which can be implemented as
	    
	    \f$  f_{qi} = \sum_{j=0}^{nq_1} \phi_q(\xi_{1,j}) \tilde{u}_{i,j} = 
	    {\bf B_1 U}  \f$
	  \f$  I_{pq} = \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i}) f_{qi} = 
	  {\bf B_0 F}  \f$
	**/
        void QuadExp::IProductWRTBase(const ConstArray<OneD,NekDouble> &base0, 
                                     const ConstArray<OneD,NekDouble> &base1, 
                                     const ConstArray<OneD,NekDouble> &inarray,
                                     Array<OneD,NekDouble> &outarray, 
                                     const int coll_check)
	{
	    int    nquad0 = m_base[0]->GetNumPoints();
	    int    nquad1 = m_base[1]->GetNumPoints();
	    ConstArray<OneD,NekDouble> jac = m_metricinfo->GetJac();
            Array<OneD,NekDouble> tmp = Array<OneD,NekDouble>(nquad0*nquad1);

	    // multiply inarray with Jacobian
	    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
	    {
		Vmath::Vmul(nquad0*nquad1,&jac[0],1,(NekDouble*)&inarray[0],1,&tmp[0],1);
	    }
	    else
	    {
		Vmath::Smul(nquad0*nquad1,jac[0],(NekDouble*)&inarray[0],1,&tmp[0],1);
	    }
	    
	    StdQuadExp::IProductWRTBase(base0,base1,tmp,outarray,coll_check);
	    
        }

	
	/** \brief  Inner product of \a inarray over region with respect to the 
	    expansion basis (this)->_Base[0] and return in \a outarray 
	    
	    Wrapper call to QuadExp::IProduct_WRT_B
	    
	    Input:\n
	    
	    - \a inarray: array of function evaluated at the physical
	    collocation points
	    
	    Output:\n
	    
	    - \a outarray: array of inner product with respect to each
	    basis over region
	    
	*/	
	void QuadExp::IProductWRTBase(const ConstArray<OneD,NekDouble> &inarray, 
                                      Array<OneD,NekDouble> &outarray)
	{
	    IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(),
			    inarray,outarray,1);
	}
	
	
	
	///////////////////////////////
	/// Differentiation Methods
	///////////////////////////////
	
	/** 
	    \brief Calculate the deritive of the physical points 
	    
	    For quadrilateral region can use the Tensor_Deriv function
	    defined under StdExpansion.
	    
	**/
        void QuadExp::PhysDeriv(const ConstArray<OneD,NekDouble> & inarray,
                                Array<OneD,NekDouble> &out_d0,
                                Array<OneD,NekDouble> &out_d1,
                                Array<OneD,NekDouble> &out_d2)
        {
	    int    nquad0 = m_base[0]->GetNumPoints();
	    int    nquad1 = m_base[0]->GetNumPoints();
	    ConstArray<TwoD,NekDouble> gmat = m_metricinfo->GetGmat();
            Array<OneD,NekDouble> Diff0 = Array<OneD,NekDouble>(nquad0*nquad1);
            Array<OneD,NekDouble> Diff1 = Array<OneD,NekDouble>(nquad0*nquad1);
	    
// 	    if(m_geom)
// 	    {
// 		ASSERTL2(n <= m_geom->GetCoordDim(),
// 			 "value of n is larger than the number of coordinates");
// 	    }
	    
	    StdExpansion2D::PhysTensorDeriv(inarray, Diff0, Diff1);
	    
	    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
	    {
                if(out_d0.num_elements())
		{
		    Vmath::Vmul  (nquad0*nquad1,&gmat[0][0],1,&Diff0[0],1, &out_d0[0], 1);
		    Vmath::Vvtvp (nquad0*nquad1,&gmat[1][0],1,&Diff1[0],1, &out_d0[0], 1,
                                  &out_d0[0],1);
		}

                if(out_d1.num_elements())
		{
		    Vmath::Vmul  (nquad0*nquad1,&gmat[2][0],1,&Diff0[0],1, &out_d1[0], 1);
		    Vmath::Vvtvp (nquad0*nquad1,&gmat[3][0],1,&Diff1[0],1, &out_d1[0], 1,
                                  &out_d1[0],1);
		}

                if(out_d2.num_elements())
		{
		    Vmath::Vmul  (nquad0*nquad1,&gmat[4][0],1,&Diff0[0],1, &out_d2[0], 1);
		    Vmath::Vvtvp (nquad0*nquad1,&gmat[5][0],1,&Diff1[0],1, &out_d2[0], 1,
                                  &out_d2[0],1);
		}
	    }
	    else // regular geometry
	    {
                if(out_d0.num_elements())
		{
		    Vmath::Smul  (nquad0*nquad1,gmat[0][0],&Diff0[0],1, &out_d0[0], 1);
                    Blas::Daxpy (nquad0*nquad1,gmat[1][0],&Diff1[0],1, &out_d0[0], 1);
// 		    Vmath::Svtvp (nquad0*nquad1,gmat[1][0],&Diff1[0],1, &out_d0[0], 1,
//                                   &out_d0[0],1);
		}

                if(out_d1.num_elements())
		{
		    Vmath::Smul  (nquad0*nquad1,gmat[2][0],&Diff0[0],1, &out_d1[0], 1);
                    Blas::Daxpy (nquad0*nquad1,gmat[3][0],&Diff1[0],1, &out_d1[0], 1);
		}

                if(out_d2.num_elements())
		{
		    Vmath::Smul  (nquad0*nquad1,gmat[4][0],&Diff0[0],1, &out_d2[0], 1);
                    Blas::Daxpy (nquad0*nquad1,gmat[5][0],&Diff1[0],1, &out_d2[0], 1);
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
        void QuadExp::FwdTrans(const ConstArray<OneD,NekDouble> & inarray, 
                               Array<OneD,NekDouble> &outarray)
	{
	    if((m_base[0]->Collocation())&&(m_base[1]->Collocation()))
	    {
		Vmath::Vcopy(GetNcoeffs(),&inarray[0],1,&m_coeffs[0],1);
	    }
	    else
	    {
		IProductWRTBase(inarray,outarray);

                LinSysKey  masskey(StdRegions::eMassMatrix,DetShapeType(),*this);
                DNekLinSysSharedPtr matsys = m_linSysManager[masskey];

                DNekVec   v(m_ncoeffs,outarray,eWrapper);
                matsys->Solve(v,v);
	    }
	}
	

        void QuadExp::GetCoords(Array<OneD,NekDouble> &coords_0,
                               Array<OneD,NekDouble> &coords_1,
                               Array<OneD,NekDouble> &coords_2)
	{
	    LibUtilities::BasisSharedPtr CBasis0;
	    LibUtilities::BasisSharedPtr CBasis1;
            Array<OneD,NekDouble>  x;

            ASSERTL0(m_geom, "m_geom not define");
	    
	    // get physical points defined in Geom
	    m_geom->FillGeom();

            switch(m_geom->GetCoordim())
            {
            case 3:
                ASSERTL0(coords_2.num_elements() != 0, 
                         "output coords_2 is not defined");

		CBasis0 = m_geom->GetBasis(2,0); 
		CBasis1 = m_geom->GetBasis(2,1);
		
		if((m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))&&
		   (m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey())))
		{
		    x = m_geom->UpdatePhys(2);
		    Blas::Dcopy(m_base[0]->GetNumPoints()*
				m_base[1]->GetNumPoints(),
				&x[0],1,&coords_2[0],1);
		}
		else // Interpolate to Expansion point distribution
		{
		    Interp2D(CBasis0->GetBasisKey(), CBasis1->GetBasisKey(),&(m_geom->UpdatePhys(2))[0],
			     m_base[0]->GetBasisKey(),m_base[1]->GetBasisKey(),&coords_2[0]);
		}
            case 2:
                ASSERTL0(coords_1.num_elements(), 
                         "output coords_1 is not defined");

		CBasis0 = m_geom->GetBasis(1,0); 
		CBasis1 = m_geom->GetBasis(1,1);
		
		if((m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))&&
		   (m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey())))
		{
		    x = m_geom->UpdatePhys(1);
		    Blas::Dcopy(m_base[0]->GetNumPoints()*
				m_base[1]->GetNumPoints(),
				&x[0],1,&coords_1[0],1);
		}
		else // Interpolate to Expansion point distribution
		{
		    Interp2D(CBasis0->GetBasisKey(), CBasis1->GetBasisKey(), &(m_geom->UpdatePhys(1))[0],
			     m_base[0]->GetBasisKey(),m_base[1]->GetBasisKey(),&coords_1[0]);
		}
            case 1:
                ASSERTL0(coords_0.num_elements(), 
                         "output coords_0 is not defined");

		CBasis0 = m_geom->GetBasis(0,0); 
		CBasis1 = m_geom->GetBasis(0,1);
		
		if((m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))&&
		   (m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey())))
		{
		    x = m_geom->UpdatePhys(0);
		    Blas::Dcopy(m_base[0]->GetNumPoints()*
				m_base[1]->GetNumPoints(),
				&x[0],1,&coords_0[0],1);
		}
		else // Interpolate to Expansion point distribution
		{
		    Interp2D(CBasis0->GetBasisKey(), CBasis1->GetBasisKey(), &(m_geom->UpdatePhys(0))[0],
			     m_base[0]->GetBasisKey(),m_base[1]->GetBasisKey(),&coords_0[0]);
		}
                break;
            default:
                ASSERTL0(false,"Number of dimensions are greater than 2");
                break;
	    }
	}
	
	// get the coordinates "coords" at the local coordinates "Lcoords"
	
        void QuadExp::GetCoord(const ConstArray<OneD,NekDouble> &Lcoords, 
                               Array<OneD,NekDouble> &coords)
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
	
	void QuadExp::WriteToFile(FILE *outfile)
	{
	    int i,j;
            Array<OneD,NekDouble> coords[3];
            int  nquad0 = m_base[0]->GetNumPoints();
            int  nquad1 = m_base[1]->GetNumPoints();

            ASSERTL0(m_geom,"_geom not defined");
            
	    int  coordim   = m_geom->GetCoordDim();
	    
            coords[0] = Array<OneD,NekDouble>(nquad0*nquad1);
            coords[1] = Array<OneD,NekDouble>(nquad0*nquad1);
            coords[2] = Array<OneD,NekDouble>(nquad0*nquad1);
	    
	    std::fprintf(outfile,"Variables = x");
	    if(coordim == 2)
	    {
                GetCoords(coords[0],coords[1]);
		fprintf(outfile,", y");
	    }
	    else if (coordim == 3)
	    {
                GetCoords(coords[0],coords[1],coords[2]);
		fprintf(outfile,", y, z");
	    }

	    fprintf(outfile,", v\n");
	    
	    fprintf(outfile,"Zone, I=%d, J=%d, F=Point\n",nquad0,nquad1);
	    for(i = 0; i < nquad0*nquad1; ++i)
	    {
		for(j = 0; j < coordim; ++j)
		{
		    fprintf(outfile,"%lf ",coords[j][i]);
		}
		fprintf(outfile,"%lf \n",m_phys[i]);
	    }
	}
	
	
	void QuadExp::WriteToFile(std::ofstream &outfile, const int dumpVar)
	{
	    int i,j;
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();
            Array<OneD,NekDouble> coords[3];

            ASSERTL0(m_geom,"m_geom not defined");

            int     coordim  = m_geom->GetCoordim();

            coords[0] = Array<OneD,NekDouble>(nquad0*nquad1);
            coords[1] = Array<OneD,NekDouble>(nquad0*nquad1);
            coords[2] = Array<OneD,NekDouble>(nquad0*nquad1);
	    
	    GetCoords(coords[0],coords[1],coords[2]);
      
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
    
	    outfile << "Zone, I=" << nquad0 << ", J=" << 
		nquad1 <<", F=Point" << std::endl;
	    
	    for(i = 0; i < nquad0*nquad1; ++i)
	    {
		for(j = 0; j < coordim; ++j)
		{
		    outfile << coords[j][i] << " ";
		}
		outfile << m_phys[i] << std::endl;
	    }
	}
	

	NekDouble QuadExp::PhysEvaluate(const ConstArray<OneD,NekDouble> &coord)
	{
            Array<OneD,NekDouble> Lcoord = Array<OneD,NekDouble>(2);
	    
	    ASSERTL0(m_geom,"m_geom not defined");
	    m_geom->GetLocCoords(coord,Lcoord);
	    
	    return StdQuadExp::PhysEvaluate(Lcoord);
	}
	
        DNekMatSharedPtr QuadExp::CreateMatrix(const MatrixKey &mkey)
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


        DNekLinSysSharedPtr QuadExp::CreateLinSys(const LinSysKey &mkey)
        {
            DNekLinSysSharedPtr returnval;

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eMassMatrix:
                returnval = MemoryManager<DNekLinSys>::AllocateSharedPtr (m_matrixManager[mkey]);
                break;
            default:
                NEKERROR(ErrorUtil::efatal, "Linear System creation not defined");
                break;
            }

            return returnval;
        }


    }//end of namespace
}//end of namespace

/** 
 *    $Log: QuadExp.cpp,v $
 *    Revision 1.14  2007/06/07 15:54:19  pvos
 *    Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
 *    Also made corrections to various ASSERTL2 calls
 *
 *    Revision 1.13  2007/06/06 11:29:31  pvos
 *    Changed ErrorUtil::Error into NEKERROR (modifications in ErrorUtil.hpp caused compiler errors)
 *
 *    Revision 1.12  2007/06/01 17:08:07  pvos
 *    Modification to make LocalRegions/Project2D run correctly (PART1)
 *
 *    Revision 1.11  2007/05/31 19:13:12  pvos
 *    Updated NodalTriExp + LocalRegions/Project2D + some other modifications
 *
 *    Revision 1.10  2007/05/31 11:38:16  pvos
 *    Updated QuadExp and TriExp
 *
 *    Revision 1.9  2007/04/26 15:00:16  sherwin
 *    SJS compiling working version using SHaredArrays
 *
 *    Revision 1.8  2006/08/05 19:03:47  sherwin
 *    Update to make the multiregions 2D expansion in connected regions work
 *
 *    Revision 1.7  2006/06/02 18:48:39  sherwin
 *    Modifications to make ProjectLoc2D run bit there are bus errors for order > 3
 *
 *    Revision 1.6  2006/06/01 15:27:27  sherwin
 *    Modifications to account for LibUtilities reshuffle
 *
 *    Revision 1.5  2006/06/01 14:15:57  sherwin
 *    Added typdef of boost wrappers and made GeoFac a boost shared pointer.
 *
 *    Revision 1.4  2006/05/30 14:00:03  sherwin
 *    Updates to make MultiRegions and its Demos work
 *
 *    Revision 1.3  2006/05/29 17:05:49  sherwin
 *    Modified to put shared_ptr around geom definitions
 *
 *    Revision 1.2  2006/05/06 20:36:16  sherwin
 *    Modifications to get LocalRegions/Project1D working
 *
 *    Revision 1.1  2006/05/04 18:58:46  kirby
 *    *** empty log message ***
 *
 *    Revision 1.23  2006/03/13 19:47:54  sherwin
 *
 *    Fixed bug related to constructor of GeoFac and also makde arguments to GeoFac all consts
 *
 *    Revision 1.22  2006/03/13 18:20:33  sherwin
 *
 *    Fixed error in ResetGmat
 *
 *    Revision 1.21  2006/03/12 21:59:48  sherwin
 *
 *    compiling version of LocalRegions
 *
 *    Revision 1.20  2006/03/12 07:43:32  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/
