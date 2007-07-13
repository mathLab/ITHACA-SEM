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

        // constructor
        SegExp::SegExp(const LibUtilities::BasisKey &Ba, 
                       const SpatialDomains::SegGeomSharedPtr &geom):
        StdRegions::StdSegExp(Ba)
        {
            m_geom = geom;    

            for(int i = 0; i < StdRegions::SIZE_MatrixType; ++i)
            {
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                StdRegions::eNoShapeType,*this),
                boost::bind(&SegExp::CreateMatrix, this, _1));
            }

            GenMetricInfo();
        }

        SegExp::SegExp(const LibUtilities::BasisKey &Ba):
            StdRegions::StdSegExp(Ba)
        {

            for(int i = 0; i < StdRegions::SIZE_MatrixType; ++i)
            {
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                StdRegions::eNoShapeType,*this),
                boost::bind(&SegExp::CreateMatrix, this, _1));
            }

            // Set up unit geometric factors. 
            m_metricinfo = MemoryManager<SpatialDomains::GeomFactors>::AllocateSharedPtr(); 
            int coordim = 1;
            Array<OneD,NekDouble> ndata = Array<OneD,NekDouble>(coordim,0.0); 
            ndata[0] = 1.0;
            m_metricinfo->ResetGmat(ndata,1,1,coordim);
            m_metricinfo->ResetJac(1,ndata);
        }

        // copy constructor
        SegExp::SegExp(const SegExp &S):
        StdRegions::StdSegExp(S)
        {
            m_geom        = S.m_geom;
            m_metricinfo  = S.m_metricinfo;
            
        }

        // by default the StdExpansion1D destructor will be called

        SegExp::~SegExp()
        {
        }

        // interpolate and possibly generate geometric factors. 
        void SegExp::GenMetricInfo()
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
                Array<OneD,NekDouble> ndata;
                ConstArray<OneD,NekDouble> odata;
                ConstArray<TwoD,NekDouble> gmat;
                int coordim = m_geom->GetCoordim();
                int  nq = m_base[0]->GetNumPoints();

                // assume all directiosn of geombasis are same
                CBasis0 = m_geom->GetBasis(0,0); 

                m_metricinfo = MemoryManager<SpatialDomains::GeomFactors>::AllocateSharedPtr(); 

                // check to see if basis are different distributions
                if(!(m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey())))
                {
                    int i;

                    // interpolate Geometric data
                    ndata = Array<OneD,NekDouble>(coordim*nq);	
                    gmat  = Xgfac->GetGmat();

                    for(i = 0; i < 2*coordim; ++i)
                    {
                        Interp1D(CBasis0->GetBasisKey(), &gmat[i][0], 
                                 m_base[0]->GetBasisKey(), 
                                 &ndata[0] + i*nq);
                    }

                    m_metricinfo->ResetGmat(ndata,nq,1,coordim);

                    // interpolate Jacobian
                    ndata = Array<OneD,NekDouble>(nq);	
                    odata = Xgfac->GetJac();

                    Interp1D(CBasis0->GetBasisKey(),odata,
                        m_base[0]->GetBasisKey(), ndata);

                    m_metricinfo->ResetJac(nq,ndata);

                    NEKERROR(ErrorUtil::ewarning,
                        "Need to check/debug routine for deformed elements");
                }
                else  // Same data can be used 
                {
                    // Copy Geometric data
                    ndata = Array<OneD,NekDouble>(coordim*nq); 
                    gmat  = Xgfac->GetGmat();
                    Blas::Dcopy(nq*coordim,&gmat[0][0],1,&ndata[0],1);

                    m_metricinfo->ResetGmat(ndata,nq,1,coordim);

                    // Copy Jacobian
                    ndata = Array<OneD,NekDouble>(nq);	

                    odata = Xgfac->GetJac();
                    Blas::Dcopy(nq,&odata[0],1,&ndata[0],1);

                    m_metricinfo->ResetJac(1,ndata);
                }   

            }
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

        NekDouble SegExp::Integral(const ConstArray<OneD,NekDouble>&  inarray)
        {

            int    nquad0 = m_base[0]->GetNumPoints();
            ConstArray<OneD,NekDouble> jac = m_metricinfo->GetJac();
            NekDouble  ival;
            Array<OneD,NekDouble> tmp   = Array<OneD,NekDouble>(nquad0);

            // multiply inarray with Jacobian

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0,&jac[0],1,(NekDouble *)&inarray[0],1,&tmp[0],1);
            }
            else
            {
                Vmath::Smul(nquad0,(NekDouble)jac[0],(NekDouble *)&inarray[0],1,&tmp[0],1);
            }

            // call StdSegExp version;
            ival = StdSegExp::Integral(tmp);

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


        void SegExp::IProductWRTBase(const ConstArray<OneD,NekDouble>& base, 
                                     const ConstArray<OneD,NekDouble>& inarray,
                                     Array<OneD,NekDouble> &outarray, 
                                     const int coll_check)
        {
            int        nquad0 = m_base[0]->GetNumPoints();
            ConstArray<OneD,NekDouble> jac = m_metricinfo->GetJac();
            Array<OneD,NekDouble> tmp = Array<OneD,NekDouble>(nquad0);


            // multiply inarray with Jacobian

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0,&jac[0],1,(NekDouble *)&inarray[0],1,&tmp[0],1);
            }
            else
            {
                Vmath::Smul(nquad0,jac[0],(NekDouble *)&inarray[0],1,&tmp[0],1);
            }

            StdSegExp::IProductWRTBase(base,tmp,outarray,coll_check);
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


        void SegExp::IProductWRTBase(const ConstArray<OneD,NekDouble>&  inarray, 
            Array<OneD,NekDouble> &outarray)
        {
            IProductWRTBase(m_base[0]->GetBdata(),inarray,outarray,1);
        }


        //----------------------------
        // Differentiation Methods
        //-----------------------------

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


        void SegExp::PhysDeriv(const ConstArray<OneD,NekDouble>& inarray,
                               Array<OneD,NekDouble> &out_d0,
                               Array<OneD,NekDouble> &out_d1,
                               Array<OneD,NekDouble> &out_d2)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            ConstArray<TwoD,NekDouble>  gmat = m_metricinfo->GetGmat();
            Array<OneD,NekDouble> Diff = Array<OneD,NekDouble>(nquad0);

//             if(m_geom)
//             {
//                 ASSERTL2(n <= m_geom->GetCoordim(),
//                     "value of n is larger than the number of coordinates");
//             }

            StdExpansion1D::PhysTensorDeriv(inarray,Diff);

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                if(out_d0.num_elements())
                {
                    Vmath::Vmul(nquad0,&gmat[0][0],1,&Diff[0],1,
                                &out_d0[0],1);
                }

                if(out_d1.num_elements())
                {
                    Vmath::Vmul(nquad0,&gmat[1][0],1,&Diff[0],1,
                        &out_d1[0],1);
                }

                if(out_d2.num_elements())
                {
                    Vmath::Vmul(nquad0,&gmat[2][0],1,&Diff[0],1,
                        &out_d2[0],1);
                }
            }
            else 
            {
                if(out_d0.num_elements())
                {
                    Vmath::Smul(nquad0,gmat[0][0],&Diff[0],1,
                        &out_d0[0],1);
                }

                if(out_d1.num_elements())
                {
                    Vmath::Smul(nquad0,gmat[1][0],&Diff[0],1,
                        &out_d1[0],1);
                }

                if(out_d2.num_elements())
                {
                    Vmath::Smul(nquad0,gmat[2][0],&Diff[0],1,
                        &out_d2[0],1);
                }
            } 
        } 

        //----------------------------
        // Evaluation Methods
        //----------------------------


        /** \brief Forward transform from physical quadrature space
        stored in \a inarray and evaluate the expansion coefficients and
        store in \a outarray

        Perform a forward transform using a Galerkin projection by
        taking the inner product of the physical points and multiplying
        by the inverse of the mass matrix using the Solve method of the
        standard matrix container holding the local mass matrix, i.e.
        \f$ {\bf \hat{u}} = {\bf M}^{-1} {\bf I} \f$ where \f$ {\bf I}[p] =
        \int^1_{-1} \phi_p(\xi_1) u(\xi_1) d\xi_1 \f$

        Inputs:\n

        - \a inarray: array of physical quadrature points to be transformed

        Outputs:\n

        - \a outarray: updated array of expansion coefficients. 

        */ 


        // need to sort out family of matrices 
        void SegExp::FwdTrans(const ConstArray<OneD,NekDouble>& inarray, 
                              Array<OneD,NekDouble> &outarray)
        {
            if(m_base[0]->Collocation())
            {
                Vmath::Vcopy(GetNcoeffs(),&inarray[0],1,&outarray[0],1);
            }
            else 
            {
                IProductWRTBase(inarray,outarray);

		// get Mass matrix inverse
		MatrixKey             masskey(StdRegions::eInvMass, DetShapeType(),*this);
		DNekScalMatSharedPtr  matsys = m_matrixManager[masskey];

                // copy inarray in case inarray == outarray
                DNekVec in (m_ncoeffs,outarray);
                DNekVec out(m_ncoeffs,outarray,eWrapper);
                
                out = (*matsys)*in;
            }
        }

        void SegExp::GetCoords(Array<OneD, NekDouble> &coords_0,
            Array<OneD, NekDouble> &coords_1,
            Array<OneD, NekDouble> &coords_2)
        { 
            Array<OneD,NekDouble>  x;

            LibUtilities::BasisSharedPtr CBasis; 

            ASSERTL0(m_geom, "m_geom not define");

            // get physical points defined in Geom
            m_geom->FillGeom();

            switch(m_geom->GetCoordim())
            {
            case 3:
                ASSERTL0(coords_2.num_elements() != 0, "output coords_2 is not defined");
                CBasis = m_geom->GetBasis(2,0);

                if(m_base[0]->GetBasisKey().SamePoints(CBasis->GetBasisKey()))
                {
                    x = m_geom->UpdatePhys(2);
                    Blas::Dcopy(m_base[0]->GetNumPoints(),&x[0],1,&coords_2[0],1);
                }
                else // Interpolate to Expansion point distribution
                {
                    Interp1D(CBasis->GetBasisKey(),&(m_geom->UpdatePhys(2))[0],
                             m_base[0]->GetBasisKey(),&coords_2[0]);
                }
            case 2:
                ASSERTL0(coords_1.num_elements() != 0, 
                    "output coords_1 is not defined");
                CBasis = m_geom->GetBasis(1,0);

                if(m_base[0]->GetBasisKey().SamePoints(CBasis->GetBasisKey()))
                {
                    x = m_geom->UpdatePhys(1);
                    Blas::Dcopy(m_base[0]->GetNumPoints(),&x[0],1,&coords_1[0],1);
                }
                else // Interpolate to Expansion point distribution
                {
                    Interp1D(CBasis->GetBasisKey(),&(m_geom->UpdatePhys(1))[0],
                             m_base[0]->GetBasisKey(),&coords_1[0]);
                }
            case 1:
                ASSERTL0(coords_0.num_elements() != 0, 
                    "output coords_2 is not defined");
                CBasis = m_geom->GetBasis(0,0);

                if(m_base[0]->GetBasisKey().SamePoints(CBasis->GetBasisKey()))
                {
                    x = m_geom->UpdatePhys(0);
                    Blas::Dcopy(m_base[0]->GetNumPoints(),&x[0],1,&coords_0[0],1);
                }
                else // Interpolate to Expansion point distribution
                {
                    Interp1D(CBasis->GetBasisKey(),&(m_geom->UpdatePhys(0))[0],
                             m_base[0]->GetBasisKey(),&coords_0[0]);
                }
                break;
            default:
                ASSERTL0(false,"Number of dimensions are greater than 2");
                break;
            }
        }


        void SegExp::GetCoord(const ConstArray<OneD,NekDouble>& Lcoords, 
                              Array<OneD,NekDouble> &coords)
        {
            int  i;

            ASSERTL1(Lcoords[0] >= -1.0&& Lcoords[0] <= 1.0,
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
            Array<OneD,NekDouble> coords[3];
            int  nquad = m_base[0]->GetNumPoints();

            ASSERTL0(m_geom,"_geom not defined");

            int  coordim  = m_geom->GetCoordim();

            coords[0] = Array<OneD,NekDouble>(nquad);
            coords[1] = Array<OneD,NekDouble>(nquad);
            coords[2] = Array<OneD,NekDouble>(nquad);

            std::fprintf(outfile,"Variables = x");
            if(coordim == 2)
            {
                GetCoords(coords[0],coords[1]);
                fprintf(outfile,", y");
            }
            else if (coordim == 3)
            {
                GetCoords(coords[0],coords[1], coords[2]);
                fprintf(outfile,", y, z");
            }
            else
            {
                GetCoords(coords[0]);
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
        }

        void SegExp::WriteToFile(std::ofstream &outfile, const int dumpVar)
        {
            int i,j;
            int     nquad = m_base[0]->GetNumPoints();

            Array<OneD,NekDouble> coords[3];

            ASSERTL0(m_geom,"m_geom not defined");

            int     coordim  = m_geom->GetCoordim();

            coords[0] = Array<OneD,NekDouble>(nquad);
            coords[1] = Array<OneD,NekDouble>(nquad);
            coords[2] = Array<OneD,NekDouble>(nquad);

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

        DNekMatSharedPtr SegExp::GetStdMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            // Need to check if matrix exists in stdMatrixManager.
            // If not then make a local expansion with standard metric info
            // and generate matrix. Otherwise direct call is OK. 
            if(!StdMatManagerAlreadyCreated(mkey))
            {
                SegExpSharedPtr tmp = MemoryManager<SegExp>::AllocateSharedPtr(m_base[0]->GetBasisKey());

                return tmp->StdSegExp::GetStdMatrix(mkey);                
            }
            else
            {
                return StdSegExp::GetStdMatrix(mkey);
            }

        }



        NekDouble SegExp::PhysEvaluate(const ConstArray<OneD,NekDouble>& coord)
        {
            Array<OneD,NekDouble> Lcoord = Array<OneD,NekDouble>(1);

            ASSERTL0(m_geom,"_geom not defined");
            m_geom->GetLocCoords(coord,Lcoord);

            return StdSegExp::PhysEvaluate(Lcoord);
        }

        DNekScalMatSharedPtr SegExp::CreateMatrix(const MatrixKey &mkey)
        {
            DNekScalMatSharedPtr returnval;

            ASSERTL2(m_metricinfo->GetGtype == SpatialDomains::eNoGeomType,"Geometric information is not set up");

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eMass:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        DNekMatSharedPtr mat = GenMatrix(StdRegions::eMass);
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(1.0,mat);
                    }
                    else
                    {
                        DNekMatSharedPtr mat = GetStdMatrix(*mkey.GetStdMatKey());

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(
                                    (m_metricinfo->GetJac())[0],mat);
                    }
                }
                break;
            case StdRegions::eInvMass:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        DNekMatSharedPtr mat = GenMatrix(StdRegions::eMass);
                        ASSERTL0(false,"Need a matrix inverse routine");

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(1.0,mat);
                    }
                    else
                    {

                        DNekMatSharedPtr mat = GetStdMatrix(*mkey.GetStdMatKey());

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(
                                      1.0/(m_metricinfo->GetJac())[0],mat );
                    }
                }
                break;
            case StdRegions::eLaplacian:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        DNekMatSharedPtr mat = GenMatrix(StdRegions::eLaplacian);

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(1.0,mat);
                    }
                    else
                    {
                        DNekMatSharedPtr mat = GetStdMatrix(*mkey.GetStdMatKey());
                        NekDouble  gfac = m_metricinfo->GetGmat()[0][0];
                        NekDouble  jac  = m_metricinfo->GetJac()[0];
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(
                                                                      gfac*gfac*jac,mat);
                    }
                }
                break;
            case StdRegions::eHelmholtz:
                {

                    NekDouble factor = mkey.GetScaleFactor();
                    MatrixKey masskey(StdRegions::eMass,
                                      mkey.GetShapeType(), *this);    
                    DNekScalMatSharedPtr mass = m_matrixManager[masskey];
                    MatrixKey lapkey(StdRegions::eLaplacian,
                                     mkey.GetShapeType(), *this);
                    DNekScalMatSharedPtr lap = m_matrixManager[lapkey];

                    DNekMatSharedPtr helm = MemoryManager<DNekMat>::AllocateSharedPtr(lap->GetRows(),lap->GetColumns());

                    (*helm) = (*lap) + 1.0/factor*(*mass);
                    // Even better:          helm  = lap + 1.0/factor*mass;
                    
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,helm);            
                }
                break;
            default:
                NEKERROR(ErrorUtil::efatal, "Matrix creation not defined");
                break;
            }

            return returnval;
        }

    } // end of namespace    
 }//end of namespace

//
// $Log: SegExp.cpp,v $
// Revision 1.23  2007/07/13 09:02:22  sherwin
// Mods for Helmholtz solver
//
// Revision 1.22  2007/07/12 12:53:01  sherwin
// Updated to have a helmholtz matrix
//
// Revision 1.21  2007/07/11 19:26:04  sherwin
// update for new Manager structure
//
// Revision 1.20  2007/07/11 06:36:23  sherwin
// Updates with MatrixManager update
//
// Revision 1.19  2007/07/10 17:17:26  sherwin
// Introduced Scaled Matrices into the MatrixManager
//
// Revision 1.18  2007/06/07 15:54:19  pvos
// Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
// Also made corrections to various ASSERTL2 calls
//
// Revision 1.17  2007/06/06 11:29:31  pvos
// Changed ErrorUtil::Error into NEKERROR (modifications in ErrorUtil.hpp caused compiler errors)
//
// Revision 1.16  2007/05/30 15:59:05  sherwin
// Fixed bug with new definition of GetLocCoords
//
// Revision 1.15  2007/05/28 08:35:25  sherwin
// Updated for localregions up to Project1D
//
// Revision 1.14  2007/05/27 16:10:28  bnelson
// Update to new Array type.
//
// Revision 1.13  2007/04/26 15:00:16  sherwin
// SJS compiling working version using SHaredArrays
//
// Revision 1.12  2007/04/08 03:33:30  jfrazier
// Minor reformatting and fixing SharedArray usage.
//
// Revision 1.11  2007/04/04 21:49:24  sherwin
// Update for SharedArray
//
// Revision 1.10  2007/03/25 15:48:21  sherwin
// UPdate LocalRegions to take new NekDouble and shared_array formats. Added new Demos
//
// Revision 1.9  2007/03/20 09:13:37  kirby
// new geomfactor routines; update for metricinfo; update style
//
// Revision 1.8  2007/03/14 21:24:07  sherwin
// Update for working version of MultiRegions up to ExpList1D
//
// Revision 1.7  2007/03/02 12:01:54  sherwin
// Update for working version of LocalRegions/Project1D
//
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
