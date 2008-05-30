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
#include <LocalRegions/LocalRegions.h>
#include <LocalRegions/LocalRegions.hpp>
#include <LocalRegions/SegExp.h>

namespace Nektar
{
    namespace LocalRegions 
    {

        // constructor
        SegExp::SegExp(const LibUtilities::BasisKey &Ba, 
                       const SpatialDomains::Geometry1DSharedPtr &geom):
            StdRegions::StdSegExp(Ba),
            m_geom(geom),
            m_metricinfo(MemoryManager<SpatialDomains::GeomFactors>::AllocateSharedPtr()),
            m_matrixManager(std::string("SegExpMatrix")),
            m_staticCondMatrixManager(std::string("SegExpStaticCondMatrix"))
        {         
            for(int i = 0; i < StdRegions::SIZE_MatrixType; ++i)
            {
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,  StdRegions::eNoExpansionType,*this),   
                                                boost::bind(&SegExp::CreateMatrix, this, _1));
                m_staticCondMatrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i, StdRegions::eNoExpansionType,*this),   
                                                          boost::bind(&SegExp::CreateStaticCondMatrix, this, _1));
            }
            GenMetricInfo();
        }

        SegExp::SegExp(const LibUtilities::BasisKey &Ba):
            StdRegions::StdSegExp(Ba),
            m_geom(),
            m_metricinfo(MemoryManager<SpatialDomains::GeomFactors>::AllocateSharedPtr()),
            m_matrixManager(std::string("SegExpMatrix")),
            m_staticCondMatrixManager(std::string("SegExpStaticCondMatrix"))
        {
            for(int i = 0; i < StdRegions::SIZE_MatrixType; ++i)
            {
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,StdRegions::eNoExpansionType,*this),  
                                                boost::bind(&SegExp::CreateMatrix, this, _1));
                m_staticCondMatrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,StdRegions::eNoExpansionType,*this),  
                                                          boost::bind(&SegExp::CreateStaticCondMatrix, this, _1));
            }

            // Set up unit geometric factors. ; 
            Array<OneD,NekDouble> ndata(1,1.0); 
            m_metricinfo->ResetGmat(ndata,1,1,1);
            m_metricinfo->ResetJac(1,ndata);
        }

        // copy constructor
        SegExp::SegExp(const SegExp &S):
            StdRegions::StdSegExp(S),
            m_geom(S.m_geom),
            m_metricinfo(S.m_metricinfo),
            m_matrixManager(std::string("SegExpMatrix")),
            m_staticCondMatrixManager(std::string("SegExpStaticCondMatrix"))
        {      
            for(int i = 0; i < StdRegions::SIZE_MatrixType; ++i)
            {
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,StdRegions::eNoExpansionType,*this),  
                                                boost::bind(&SegExp::CreateMatrix, this, _1));
                m_staticCondMatrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,StdRegions::eNoExpansionType,*this),  
                                                          boost::bind(&SegExp::CreateStaticCondMatrix, this, _1));
            }       
        }

        // by default the StdExpansion1D destructor will be called

        SegExp::~SegExp()
        {
        }

        // interpolate and possibly generate geometric factors. 
        void SegExp::GenMetricInfo()
        {
            SpatialDomains::GeomFactorsSharedPtr Xgfac = m_geom->GetGeomFactors();

            if( Xgfac->GetGtype() != SpatialDomains::eDeformed)
            {
                m_metricinfo = Xgfac;
            }
            else
            {                
                int coordim = m_geom->GetCoordim();
                int  nq = m_base[0]->GetNumPoints();
                Array<OneD,NekDouble> ndata(coordim*nq); 
                Array<OneD, const NekDouble> odata = Xgfac->GetJac();
                Array<TwoD, const NekDouble> gmat  = Xgfac->GetGmat();;

                // assume all directiosn of geombasis are same
                LibUtilities::BasisSharedPtr CBasis0 = m_geom->GetBasis(0,0); 

                // check to see if basis are different distributions
                if(!(m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey())))
                {
                    // interpolate Geometric data   
                    for(int i = 0; i < 2*coordim; ++i)
                    {
                        Interp1D(CBasis0->GetBasisKey(), &gmat[i][0], 
                                 m_base[0]->GetBasisKey(), 
                                 &ndata[0] + i*nq);
                    }

                    m_metricinfo->ResetGmat(ndata,nq,1,coordim);

                    // interpolate Jacobian
                    ndata = Array<OneD,NekDouble>(nq);   
                     
                    Interp1D(CBasis0->GetBasisKey(),odata,
                        m_base[0]->GetBasisKey(), ndata);

                    m_metricinfo->ResetJac(nq,ndata);                    
                    m_metricinfo->ResetNormals(Xgfac->GetNormals());

                    NEKERROR(ErrorUtil::ewarning,
                        "Need to check/debug routine for deformed elements");
                }
                else  // Same data can be used 
                {
                    // Copy Geometric data
                    Blas::Dcopy(nq*coordim,&gmat[0][0],1,&ndata[0],1);

                    m_metricinfo->ResetGmat(ndata,nq,1,coordim);

                    // Copy Jacobian
                    ndata = Array<OneD,NekDouble>(nq);    
                    Blas::Dcopy(nq,&odata[0],1,&ndata[0],1);

                    m_metricinfo->ResetJac(1,ndata);                    
                    m_metricinfo->ResetNormals(Xgfac->GetNormals());
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

        NekDouble SegExp::Integral(const Array<OneD, const NekDouble>&  inarray)
        {

            int    nquad0 = m_base[0]->GetNumPoints();
            Array<OneD, const NekDouble> jac = m_metricinfo->GetJac();
            NekDouble  ival;
            Array<OneD,NekDouble> tmp(nquad0);

            // multiply inarray with Jacobian
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0, jac, 1, inarray, 1, tmp,1);
            }
            else
            {
                Vmath::Smul(nquad0, jac[0], inarray, 1, tmp, 1);
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


        void SegExp::IProductWRTBase(const Array<OneD, const NekDouble>& base, 
                                     const Array<OneD, const NekDouble>& inarray,
                                     Array<OneD,NekDouble> &outarray, 
                                     const int coll_check)
        {
            int        nquad0 = m_base[0]->GetNumPoints();
            Array<OneD, const NekDouble> jac = m_metricinfo->GetJac();
            Array<OneD,NekDouble> tmp(nquad0);


            // multiply inarray with Jacobian
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0, jac, 1, inarray, 1, tmp, 1);
            }
            else
            {
                Vmath::Smul(nquad0, jac[0], inarray, 1, tmp, 1);
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


        void SegExp::IProductWRTBase(const Array<OneD, const NekDouble>&  inarray, 
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


        void SegExp::PhysDeriv(const Array<OneD, const NekDouble>& inarray,
                               Array<OneD,NekDouble> &out_d0,
                               Array<OneD,NekDouble> &out_d1,
                               Array<OneD,NekDouble> &out_d2)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            Array<TwoD, const NekDouble>  gmat = m_metricinfo->GetGmat();
            Array<OneD,NekDouble> diff(nquad0);

            StdExpansion1D::PhysTensorDeriv(inarray,diff);

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                if(out_d0.num_elements())
                {
                    Vmath::Vmul(nquad0,&gmat[0][0],1,&diff[0],1,
                                &out_d0[0],1);
                }

                if(out_d1.num_elements())
                {
                    Vmath::Vmul(nquad0,&gmat[1][0],1,&diff[0],1,
                                &out_d1[0],1);
                }

                if(out_d2.num_elements())
                {
                    Vmath::Vmul(nquad0,&gmat[2][0],1,&diff[0],1,
                        &out_d2[0],1);
                }
            }
            else 
            {
                if(out_d0.num_elements())
                {
                    Vmath::Smul(nquad0, gmat[0][0], diff, 1,
                        out_d0, 1);
                }

                if(out_d1.num_elements())
                {
                    Vmath::Smul(nquad0, gmat[1][0], diff, 1,
                        out_d1, 1);
                }

                if(out_d2.num_elements())
                {
                    Vmath::Smul(nquad0, gmat[2][0], diff, 1,
                        out_d2, 1);
                }
            } 
        } 

        //----------------------------
        // Evaluation Methods
        //----------------------------


        void SegExp::SetCoeffsToOrientation(StdRegions::EdgeOrientation dir)
        {

            if(dir == StdRegions::eBackwards)
            {
                Array<OneD,NekDouble> outarray(m_ncoeffs);
                ReverseCoeffsAndSign(m_coeffs,outarray);
                Vmath::Vcopy(m_ncoeffs,outarray,1,m_coeffs,1);
            }
        }

        // Reverse the coefficients in a boundary interior expansion
        // this routine is of use when we need the segment
        // coefficients corresponding to a expansion in the reverse
        // coordinate direction
        void SegExp::ReverseCoeffsAndSign(const Array<OneD,NekDouble> &inarray,
                                          Array<OneD,NekDouble> &outarray)
        {

            int m;
            NekDouble sgn = 1;

            ASSERTL1(&inarray[0] != &outarray[0],"inarray and outarray can not be the same");
            switch(GetBasisType(0))
            {
            case LibUtilities::eModified_A:
                
                //Swap vertices
                outarray[0] = inarray[1]; 
                outarray[1] = inarray[0];
                
                // negate odd modes
                for(m = 2; m < m_ncoeffs; ++m)
                {
                    outarray[m] = sgn*inarray[m];
                    sgn = -sgn;
                }
                break;
            case LibUtilities::eGLL_Lagrange:

                for(m = 0; m < m_ncoeffs; ++m)
                {
                    outarray[m_ncoeffs-1-m] = inarray[m];
                }
                
            default:
                ASSERTL0(false,"This basis is not allowed in this method");
                break;
            }
        }

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
        void SegExp::FwdTrans(const Array<OneD, const NekDouble>& inarray, 
                              Array<OneD,NekDouble> &outarray)
        {
            if(m_base[0]->Collocation())
            {
                Vmath::Vcopy(m_ncoeffs, inarray, 1, outarray, 1);
            }
            else 
            {
                IProductWRTBase(inarray,outarray);

                // get Mass matrix inverse
                MatrixKey             masskey(StdRegions::eInvMass, DetExpansionType(),*this);
                DNekScalMatSharedPtr& matsys = m_matrixManager[masskey];
                
                // copy inarray in case inarray == outarray
                NekVector<const NekDouble> in(m_ncoeffs,outarray,eWrapper);
                NekVector<NekDouble> out(m_ncoeffs,outarray,eWrapper);
                
                out = (*matsys)*in;
            }
        }

        void SegExp::FwdTrans_BndConstrained(const Array<OneD, const NekDouble>& inarray, 
                                             Array<OneD, NekDouble> &outarray)
        {
            if(m_base[0]->Collocation())
            {
                Vmath::Vcopy(m_ncoeffs, inarray, 1, outarray, 1);
            }
            else
            {
                int nInteriorDofs = m_ncoeffs-2;
                int offset;

                switch(m_base[0]->GetBasisType())
                {
                case LibUtilities::eGLL_Lagrange:
                    {
                        offset = 1;
                    }
                    break;
                case LibUtilities::eModified_A:
                    {
                        offset = 2;
                    }
                    break;
                default:
                    ASSERTL0(false,"This type of FwdTrans is not defined for this expansion type");
                }    

                fill(outarray.get(), outarray.get()+m_ncoeffs, 0.0 );

                outarray[GetVertexMap(0)] = inarray[0];
                outarray[GetVertexMap(1)] = inarray[m_base[0]->GetNumPoints()-1];

                Array<OneD, NekDouble> tmp0(m_ncoeffs); //ideally, we would like to have tmp0 to be replaced by outarray (currently MassMatrixOp does not allow aliasing)
                Array<OneD, NekDouble> tmp1(m_ncoeffs);

                MassMatrixOp(outarray,tmp0);
                IProductWRTBase(inarray,tmp1);

                Vmath::Vsub(m_ncoeffs, tmp1, 1, tmp0, 1, tmp1, 1);
                
                // get Mass matrix inverse (only of interior DOF)
                MatrixKey             masskey(StdRegions::eMass, DetExpansionType(),*this);
                DNekScalMatSharedPtr  matsys = (m_staticCondMatrixManager[masskey])->GetBlock(1,1);
                
                Blas::Dgemv('N',nInteriorDofs,nInteriorDofs, matsys->Scale(), &((matsys->GetOwnedMatrix())->GetPtr())[0],
                            nInteriorDofs,tmp1.get()+offset,1,0.0,outarray.get()+offset,1);       
            }
        }

        void SegExp::GetCoords(Array<OneD, NekDouble> &coords_0,
            Array<OneD, NekDouble> &coords_1,
            Array<OneD, NekDouble> &coords_2)
        { 
            Array<OneD,NekDouble>  x;

            LibUtilities::BasisSharedPtr CBasis; 
            ASSERTL0(m_geom, "m_geom not defined");

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
                    Blas::Dcopy(m_base[0]->GetNumPoints(), x, 1, coords_2, 1);
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
                    Blas::Dcopy(m_base[0]->GetNumPoints(), x, 1, coords_1, 1);
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
                    Blas::Dcopy(m_base[0]->GetNumPoints(), x, 1, coords_0, 1);
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


        void SegExp::GetCoord(const Array<OneD, const NekDouble>& Lcoords, 
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

        void SegExp::WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar)
        {
            if(format==eTecplot)
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
            else if(format==eGmsh)
            {  
                int i,j;
                int     nquad = m_base[0]->GetNumPoints();
                
                Array<OneD,NekDouble> coords[3];
                
                ASSERTL0(m_geom,"m_geom not defined");             
                coords[0] = Array<OneD,NekDouble>(nquad,0.0);
                coords[1] = Array<OneD,NekDouble>(nquad,0.0);
                coords[2] = Array<OneD,NekDouble>(nquad,0.0);
                
                GetCoords(coords[0],coords[1],coords[2]);

                if(dumpVar)
                {
                    outfile<<"View.Type = 2;"<<endl;
                    outfile<<"View \" \" {"<<endl;
                }
 
                for(i = 0; i < nquad; ++i)
                {
                    outfile << "SP(" << coords[0][i] << ", " << coords[1][i] << ", " << coords[2][i] << ")";
                    outfile << "{" << m_phys[i] << "};" << endl;
                }

                if(dumpVar)
                { 
                    outfile << "};" << endl;
                }
            }
            else
            {
                ASSERTL0(false, "Output routine not implemented for requested type of output");
            }
        }

        DNekMatSharedPtr SegExp::CreateStdMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            LibUtilities::BasisKey bkey = m_base[0]->GetBasisKey();
            StdRegions::StdSegExpSharedPtr tmp = MemoryManager<StdSegExp>::AllocateSharedPtr(bkey);
            
            return tmp->GetStdMatrix(mkey); 
        }

        NekDouble SegExp::PhysEvaluate(const Array<OneD, const NekDouble>& coord)
        {
            Array<OneD,NekDouble> Lcoord = Array<OneD,NekDouble>(1);

            ASSERTL0(m_geom,"_geom not defined");
            m_geom->GetLocCoords(coord,Lcoord);

            return StdSegExp::PhysEvaluate(Lcoord);
        }

        DNekScalMatSharedPtr SegExp::CreateMatrix(const MatrixKey &mkey)
        {
            DNekScalMatSharedPtr returnval;
            NekDouble fac;

            ASSERTL2(m_metricinfo->GetGtype == SpatialDomains::eNoGeomType,"Geometric information is not set up");

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eMass:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        fac = 1.0;
                        goto UseLocRegionsMatrix;
                    }
                    else
                    {
                        fac = (m_metricinfo->GetJac())[0];
                        goto UseStdRegionsMatrix;
                    }
                }
                break;
            case StdRegions::eInvMass:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        StdRegions::StdMatrixKey masskey(StdRegions::eMass,DetExpansionType(),
                                                         *this);
                        DNekMatSharedPtr mat = GenMatrix(masskey);
                        mat->Invert();

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        fac = 1.0/(m_metricinfo->GetJac())[0];
                        goto UseStdRegionsMatrix;
                    }
                }
                break;
            case StdRegions::eWeakDeriv0:
            case StdRegions::eWeakDeriv1:
            case StdRegions::eWeakDeriv2:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        fac = 1.0; 
                        goto UseLocRegionsMatrix;
                    }
                    else
                    {
                        int dir;
                        switch(mkey.GetMatrixType())
                        {
                        case StdRegions::eWeakDeriv0:
                            dir = 0;
                            break;
                        case StdRegions::eWeakDeriv1:
                            ASSERTL1(m_geom->GetCoordim() >= 2,"Cannot call eWeakDeriv2 in a coordinate system which is not at least two-dimensional");
                            dir = 1;
                            break;
                        case StdRegions::eWeakDeriv2:
                            ASSERTL1(m_geom->GetCoordim() == 3,"Cannot call eWeakDeriv2 in a coordinate system which is not three-dimensional");
                            dir = 2;
                            break;
                        }
                        fac = m_metricinfo->GetGmat()[dir][0]*
                            m_metricinfo->GetJac()[0];
                        goto UseStdRegionsMatrix;
                    }
                }
                break;
            case StdRegions::eLaplacian:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        fac = 1.0;
                        goto UseLocRegionsMatrix;
                    }
                    else
                    {
                        int coordim = m_geom->GetCoordim();
                        fac = 0.0;
                        for(int i = 0; i < coordim; ++i)
                        {
                            fac += m_metricinfo->GetGmat()[i][0]*
                                m_metricinfo->GetGmat()[i][0];
                        }
                        fac *= m_metricinfo->GetJac()[0];
                        goto UseStdRegionsMatrix;
                    }
                }
                break;
            case StdRegions::eHelmholtz:
                {
                    NekDouble factor = mkey.GetScaleFactor();
                    MatrixKey masskey(StdRegions::eMass,
                                      mkey.GetExpansionType(), *this);    
                    DNekScalMat &MassMat = *(this->m_matrixManager[masskey]);
                    MatrixKey lapkey(StdRegions::eLaplacian,
                                     mkey.GetExpansionType(), *this);
                    DNekScalMat &LapMat = *(this->m_matrixManager[lapkey]);

                    int rows = LapMat.GetRows();
                    int cols = LapMat.GetColumns();

                    DNekMatSharedPtr helm = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);
                    
                    NekDouble one = 1.0;
                    (*helm) = LapMat + factor*MassMat;
                    
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,helm);            
                }
                break;
            case StdRegions::eUnifiedDGHelmholtz:
            case StdRegions::eUnifiedDGLamToU:
            case StdRegions::eUnifiedDGLamToQ0:
            case StdRegions::eUnifiedDGHelmBndSys:
                {
                    NekDouble one    = 1.0;
                    
                    DNekMatSharedPtr mat = GenMatrix(*mkey.GetStdMatKey());
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                }
            break;
            case StdRegions::eInvUnifiedDGHelmholtz:
                {
                    NekDouble one = 1.0;

                    StdRegions::StdMatrixKey hkey(StdRegions::eUnifiedDGHelmholtz,
                                                  DetExpansionType(),*this,
                                                  mkey.GetConstant(0),
                                                  mkey.GetConstant(1));
                    DNekMatSharedPtr mat = GenMatrix(hkey);

                    mat->Invert();
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                }
                break;
            UseLocRegionsMatrix:
                {
                    DNekMatSharedPtr mat = GenMatrix(*mkey.GetStdMatKey());
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(fac,mat);
                }
                break;
            UseStdRegionsMatrix:
                {
                    DNekMatSharedPtr mat = GetStdMatrix(*mkey.GetStdMatKey());
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(fac,mat);
                }
                break;
            default:
                NEKERROR(ErrorUtil::efatal, "Matrix creation not defined");
                break;
            }

            return returnval;
        }
        
        DNekMatSharedPtr SegExp::GenMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            
            DNekMatSharedPtr returnval;

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eUnifiedDGHelmBndSys:
                {
                    int nbndry = NumBndryCoeffs();
                    int nquad  = GetNumPoints(0);
                    int i,j,k;
                    Array<OneD,NekDouble> b(nbndry);
                    const Array<OneD, const NekDouble> &Basis = m_base[0]->GetBdata();
                    NekDouble lambdaval = mkey.GetConstant(0);
                    NekDouble tau       = mkey.GetConstant(1);

                    Array<OneD,unsigned int> bmap(NumBndryCoeffs());
                    GetBoundaryMap(bmap);
                    
                    // declare matrix space
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry,
                                                                          nbndry);
                    DNekMat &BndMat = *returnval;
                    
                    // Matrix to map Lambda to U
                    LocalRegions::MatrixKey Umatkey(StdRegions::eUnifiedDGLamToU, DetExpansionType(),*this, lambdaval,tau);
                    DNekScalMat &LamToU = *GetLocMatrix(Umatkey); 
                
                    // Matrix to map Lambda to Q
                    LocalRegions::MatrixKey Qmatkey(StdRegions::eUnifiedDGLamToQ0,
                                            DetExpansionType(),*this, lambdaval,tau);
                    DNekScalMat &LamToQ = *GetLocMatrix(Qmatkey); 
                    

#if 0  // Equivalent symmetric form over all elements
                    // Mass matrix 
                    LocalRegions::MatrixKey masskey(StdRegions::eMass,
                                                    DetShapeType(),*this);
                    DNekScalMat &Mass = *GetLocMatrix(masskey); 

#if 0 // What I would like to call
                    BndMat = Transpose(LamToQ)*Mass*LamToQ + lambdaval*Transpose(LamToU)*Mass*LamToU;
#else
                    DNekMat QTransMat = *LamToQ.GetOwnedMatrix();
                    QTransMat.Transpose();

                    DNekMat UTransMat = *LamToU.GetOwnedMatrix();
                    UTransMat.Transpose();

                    BndMat = LamToQ.Scale()*QTransMat*Mass*LamToQ + 
                        lambdaval*LamToU.Scale()*UTransMat*Mass*LamToU; 
#endif
                    
                    // Add boundary terms
                    // upper terms
                    
                    for(i = 0; i < nbndry; ++i)
                    {
                        b[i] = 0.0;
                        for(k = 0; k < nbndry; ++k)
                        {
                            b[i] += LamToU(bmap[k],i)*Basis[(bmap[k]+1)*nquad-1];
                        }
                    }
                    
                    for(i = 0; i < nbndry; ++i)
                    {
                        for(j = 0; j < nbndry; ++j)
                        {
                            BndMat(i,j) += tau*(b[i]*b[j] - 
                                                b[i]*Basis[(j+1)*nquad-1] - 
                                                Basis[(i+1)*nquad-1]*b[j] +
                                                Basis[(i+1)*nquad-1]*
                                                Basis[(j+1)*nquad-1]);
                        }
                    }
                    
                    // lower terms
                    for(i = 0; i < nbndry; ++i)
                    {
                        b[i] = 0.0;
                        for(k = 0; k < nbndry; ++k)
                        {
                            b[i] += LamToU(bmap[k],i)*Basis[bmap[k]*nquad];
                        }
                    }
                    
                    for(i = 0; i < nbndry; ++i)
                    {
                        for(j = 0; j < nbndry; ++j)
                        {
                            BndMat(i,j) -= tau*(b[i]*b[j] - 
                                                b[i]*Basis[j*nquad] - 
                                                Basis[i*nquad]*b[j] +
                                                Basis[i*nquad]*Basis[j*nquad]);
                        }
                    }
#else // use flux form directly 
                    for(j = 0; j < nbndry; ++j)
                    {
                        BndMat(0,j) = -LamToQ(bmap[0],j) + tau*LamToU(bmap[0],j);
                    }

                    for(j = 0; j < nbndry; ++j)
                    {
                        BndMat(1,j) =  LamToQ(bmap[1],j) - tau*LamToU(bmap[1],j);
                    }

#endif
                    //  cout << BndMat << endl;
                }
                break;
            default:
                returnval = StdSegExp::GenMatrix(mkey);
                break;
            }

            return returnval;
        }

        void SegExp::AddNormBoundaryInt(const int dir, 
                                        Array<OneD, const NekDouble> &inarray,
                                        Array<OneD,NekDouble> &outarray) 
       {

            int k;
            int nbndry = NumBndryCoeffs();
            int nquad  = GetNumPoints(0);
            const Array<OneD, const NekDouble> &Basis  = m_base[0]->GetBdata();
            StdRegions::StdExpMap vmap;

            MapTo(StdRegions::eForwards,vmap);

            // add G (\lambda - ulam) = G x F term (can
            // assume G is diagonal since one of the basis
            // is zero at boundary otherwise)
            for(k = 0; k < nbndry; ++k)
            {
                outarray[vmap[k]] += (Basis[(vmap[k]+1)*nquad-1]*Basis[(vmap[k]+1)*nquad-1] - Basis[vmap[k]*nquad]*Basis[vmap[k]*nquad])*inarray[vmap[k]];
            }
        }

        void SegExp::AddUDGHelmholtzBoundaryTerms(const NekDouble tau, 
                                    const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD,NekDouble> &outarray,
                                    bool MatrixTerms)
        {
            int i,j,k,n;
            int nbndry = NumBndryCoeffs();
            int nquad  = GetNumPoints(0);
            NekDouble  val, val1;
            StdRegions::StdExpMap vmap;
            
            ASSERTL0(&inarray[0] != &outarray[0],"Input and output arrays use the same memory");

            const Array<OneD, const NekDouble> &Dbasis = m_base[0]->GetDbdata();
            const Array<OneD, const NekDouble> &Basis  = m_base[0]->GetBdata();
            
            MatrixKey    masskey(StdRegions::eInvMass, DetExpansionType(),*this);
            DNekScalMat  &invMass = *m_matrixManager[masskey];
            Array<TwoD, const NekDouble>  gmat = m_metricinfo->GetGmat();
            NekDouble rx0,rx1;
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                rx0 = gmat[0][0];
                rx1 = gmat[0][nquad-1];
            }
            else
            {
                rx0 = rx1 = gmat[0][0];
            }

            MapTo(StdRegions::eForwards,vmap);


            if(MatrixTerms == true) //term which arise in matrix formulations but not rhs boundary terms. 
            {
                // Add -D^T M^{-1}G operation =-<n phi_i, d\phi_j/dx>
                for(i = 0; i < nbndry; ++i)
                {
                    for(j = 0; j < m_ncoeffs; ++j)
                    {
                        outarray[vmap[i]] -= Basis[(vmap[i]+1)*nquad-1]*Dbasis[(j+1)*nquad-1]*rx1*inarray[j];
                        outarray[vmap[i]] += Basis[vmap[i]*nquad]*Dbasis[j*nquad]*rx0*inarray[j];
                    }
                }
            }
            
            //Add -E^T M^{-1}D_i^e = -< d\phi_i/dx, n  phi_j>
            for(i = 0; i < m_ncoeffs; ++i)
            {
                for(j = 0; j < nbndry; ++j)
                {
                    outarray[i] -= Dbasis[(i+1)*nquad-1]*rx1*Basis[(vmap[j]+1)*nquad-1]*inarray[vmap[j]];
                    outarray[i] += Dbasis[i*nquad]*rx0*Basis[vmap[j]*nquad]*inarray[vmap[j]];
                }
            }
            
            // Add F = \tau <phi_i,phi_j> (note phi_i is zero if phi_j is non-zero)
            for(i = 0; i < nbndry; ++i)
            {
                outarray[vmap[i]] += tau*Basis[(vmap[i]+1)*nquad-1]*Basis[(vmap[i]+1)*nquad-1]*inarray[vmap[i]];
                outarray[vmap[i]] += tau*Basis[vmap[i]*nquad]*Basis[vmap[i]*nquad]*inarray[vmap[i]];
            }

            // Add E M^{-1} G term 
            for(i = 0; i < nbndry; ++i)
            {
                for(n = 0; n < nbndry; ++n)
                {
                    // evaluate M^{-1} G
                    val1 = 0.0;
                    for(k = 0; k < nbndry; ++k)
                    {
                        val = 0.0;
                        for(j = 0; j < nbndry; ++j)
                        {
                            val += (Basis[(vmap[k]+1)*nquad-1]*Basis[(vmap[j]+1)*nquad-1] - Basis[vmap[k]*nquad]*Basis[vmap[j]*nquad])*inarray[vmap[j]];
                        }
                        
                        val1 += invMass(vmap[n],vmap[k])*val; 
                    }

                    outarray[vmap[i]] += (Basis[(vmap[i]+1)*nquad-1]*Basis[(vmap[n]+1)*nquad-1] - Basis[vmap[i]*nquad]*Basis[vmap[n]*nquad])*val1; 
                }
            }
        }
        
        DNekScalBlkMatSharedPtr SegExp::CreateStaticCondMatrix(const MatrixKey &mkey)
        {
            DNekScalBlkMatSharedPtr returnval;
            
            ASSERTL2(m_metricinfo->GetGtype == SpatialDomains::eNoGeomType,"Geometric information is not set up");

            // set up block matrix system
            int nbdry = NumBndryCoeffs();
            int nint = m_ncoeffs - nbdry;
            Array<OneD, unsigned int> exp_size(2);
            exp_size[0] = nbdry;
            exp_size[1] = nint;
            returnval = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(exp_size,exp_size); //Really need a constructor which takes Arrays
            NekDouble factor = 1.0;

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eHelmholtz: // special case since Helmholtz not defined in StdRegions

                // use Deformed case for both regular and deformed geometries 
                factor = 1.0;
                goto UseLocRegionsMatrix;
                break;
            default:
                if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                {
                    factor = 1.0;
                    goto UseLocRegionsMatrix;
                }
                else
                {
                    DNekScalMatSharedPtr mat = GetLocMatrix(mkey);
                    factor = mat->Scale();
                    goto UseStdRegionsMatrix;
                }
                break;
            UseStdRegionsMatrix:
                {
                    NekDouble            invfactor = 1.0/factor;
                    NekDouble            one = 1.0;
                    DNekBlkMatSharedPtr& mat = GetStdStaticCondMatrix(*(mkey.GetStdMatKey()));                    
                    DNekScalMatSharedPtr Atmp;
                    DNekMatSharedPtr     Asubmat;

                    returnval->SetBlock(0,0,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,Asubmat = mat->GetBlock(0,0)));
                    returnval->SetBlock(0,1,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Asubmat = mat->GetBlock(0,1)));
                    returnval->SetBlock(1,0,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,Asubmat = mat->GetBlock(1,0)));
                    returnval->SetBlock(1,1,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(invfactor,Asubmat = mat->GetBlock(1,1)));
                }
                break;
            UseLocRegionsMatrix:
                {
                    int i,j;
                    int cnt = 0;
                    NekDouble            invfactor = 1.0/factor;
                    NekDouble            one = 1.0;
                    DNekScalMat &mat = *GetLocMatrix(mkey);
                    DNekMatSharedPtr A = MemoryManager<DNekMat>::AllocateSharedPtr(nbdry,nbdry);
                    DNekMatSharedPtr B = MemoryManager<DNekMat>::AllocateSharedPtr(nbdry,nint);
                    DNekMatSharedPtr C = MemoryManager<DNekMat>::AllocateSharedPtr(nint,nbdry);
                    DNekMatSharedPtr D = MemoryManager<DNekMat>::AllocateSharedPtr(nint,nint);

                    Array<OneD,unsigned int> bmap(nbdry);
                    Array<OneD,unsigned int> imap(nint);
                    GetBoundaryMap(bmap);
                    GetInteriorMap(imap);

                    for(i = 0; i < nbdry; ++i)
                    {
                        for(j = 0; j < nbdry; ++j)
                        {
                            (*A)(i,j) = mat(bmap[i],bmap[j]);
                        }
                        
                        for(j = 0; j < nint; ++j)
                        {
                            (*B)(i,j) = mat(bmap[i],imap[j]);
                        }
                    }
                    
                    for(i = 0; i < nint; ++i)
                    {
                        for(j = 0; j < nbdry; ++j)
                        {
                            (*C)(i,j) = mat(imap[i],bmap[j]);
                        }
                        
                        for(j = 0; j < nint; ++j)
                        {
                            (*D)(i,j) = mat(imap[i],imap[j]);
                        }
                    }
                    
                    // Calculate static condensed system 
                    if(nint)
                    {
                        D->Invert();
                        (*B) = (*B)*(*D);
                        (*A) = (*A) - (*B)*(*C);
                    }
                    
                    DNekScalMatSharedPtr     Atmp;

                    returnval->SetBlock(0,0,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,A));
                    returnval->SetBlock(0,1,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,B));
                    returnval->SetBlock(1,0,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,C));
                    returnval->SetBlock(1,1,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(invfactor,D));
                }
            }


            
            return returnval;
        }

    } // end of namespace    
 }//end of namespace

//
// $Log: SegExp.cpp,v $
// Revision 1.43  2008/05/29 21:33:37  pvos
// Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
//
// Revision 1.42  2008/05/29 01:02:13  bnelson
// Added precompiled header support.
//
// Revision 1.41  2008/05/14 18:06:50  sherwin
// mods to fix Seggeom to Geometry1D casting
//
// Revision 1.40  2008/05/10 18:27:33  sherwin
// Modifications necessary for QuadExp Unified DG Solver
//
// Revision 1.39  2008/04/06 05:59:05  bnelson
// Changed ConstArray to Array<const>
//
// Revision 1.38  2008/04/02 22:19:26  pvos
// Update for 2D local to global mapping
//
// Revision 1.37  2008/03/18 14:12:53  pvos
// Update for nodal triangular helmholtz solver
//
// Revision 1.36  2008/03/12 15:24:29  pvos
// Clean up of the code
//
// Revision 1.35  2008/02/28 10:04:11  sherwin
// Modes for UDG codes
//
// Revision 1.34  2008/01/21 19:59:32  sherwin
// Updated to take SegGeoms instead of EdgeComponents
//
// Revision 1.33  2007/12/17 13:04:30  sherwin
// Modified GenMatrix to take a StdMatrixKey and removed m_constant from MatrixKey
//
// Revision 1.32  2007/12/06 22:49:09  pvos
// 2D Helmholtz solver updates
//
// Revision 1.31  2007/11/20 16:28:46  sherwin
// Added terms for UDG Helmholtz solver
//
// Revision 1.30  2007/10/04 12:10:03  sherwin
// Update for working version of static condensation in Helmholtz1D and put lambda coefficient on the mass matrix rather than the Laplacian operator.
//
// Revision 1.29  2007/10/03 11:37:50  sherwin
// Updates relating to static condensation implementation
//
// Revision 1.28  2007/08/10 03:38:08  jfrazier
// Updated with new rev of NekManager.
//
// Revision 1.27  2007/07/30 21:00:06  sherwin
// Temporary fix in CreateMatrix
//
// Revision 1.26  2007/07/28 05:09:33  sherwin
// Fixed version with updated MemoryManager
//
// Revision 1.25  2007/07/20 00:45:51  bnelson
// Replaced boost::shared_ptr with Nektar::ptr
//
// Revision 1.24  2007/07/13 15:22:11  sherwin
// Update for Helmholtz (working without bcs )
//
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
