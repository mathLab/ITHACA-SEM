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

#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/SegExp.h>
#include <LibUtilities/Foundations/Interp.h>


namespace Nektar
{
    namespace LocalRegions 
    {

        /**
         * @class SegExp
         * Defines a Segment local expansion.
         */

        /// Constructor using BasisKey class for quadrature points and
        /// order definition.
        /**
         * @param   Ba          Basis key of segment expansion.
         * @param   geom        Description of geometry.
         */
        SegExp::SegExp( const LibUtilities::BasisKey &Ba,
                        const SpatialDomains::Geometry1DSharedPtr &geom):
            StdExpansion(Ba.GetNumModes(), 1, Ba),
            StdExpansion1D(Ba.GetNumModes(), Ba),
            StdRegions::StdSegExp(Ba),
            m_geom(geom),
            m_metricinfo(m_geom->GetGeomFactors(m_base)),
            m_matrixManager(
                    boost::bind(&SegExp::CreateMatrix, this, _1),
                    std::string("SegExpMatrix")),
            m_staticCondMatrixManager(
                    boost::bind(&SegExp::CreateStaticCondMatrix, this, _1),
                    std::string("SegExpStaticCondMatrix"))
        {
        }


        /// Copy Constructor
        /**
         * @param   S           Existing segment to duplicate.
         */
        SegExp::SegExp(const SegExp &S):
            StdRegions::StdSegExp(S),
            m_geom(S.m_geom),
            m_metricinfo(S.m_metricinfo),
            m_matrixManager(S.m_matrixManager),
            m_staticCondMatrixManager(S.m_staticCondMatrixManager)
        {
        }


        /**
         *
         */
        SegExp::~SegExp()
        {
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

        NekDouble SegExp::v_Integral(
                const Array<OneD, const NekDouble>&  inarray)
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
            ival = StdSegExp::v_Integral(tmp);
            //ival = StdSegExp::Integral(tmp);
            return ival;
        }


        //-----------------------------
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
        void SegExp::v_PhysDeriv(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,NekDouble> &out_d0,
                      Array<OneD,NekDouble> &out_d1,
                      Array<OneD,NekDouble> &out_d2)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            Array<TwoD, const NekDouble>  gmat = m_metricinfo->GetGmat();
            Array<OneD,NekDouble> diff(nquad0);

            //StdExpansion1D::PhysTensorDeriv(inarray,diff);
            PhysTensorDeriv(inarray,diff);
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

        /**
        *\brief Evaluate the derivative along a line:
        * \f$ d/ds=\frac{spacedim}{||tangent||}d/d{\xi}  \f$.
        * The derivative is calculated performing
        *the product \f$ du/d{s}=\nabla u \cdot tangent \f$.
        *\param inarray function to derive
        *\param out_ds result of the derivative operation 
        **/
        void SegExp::v_PhysDeriv_s(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,NekDouble> &out_ds)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int     coordim  = m_geom->GetCoordim();
            Array<OneD, NekDouble> diff (nquad0);
            //this operation is needed if you put out_ds==inarray
            Vmath::Zero(nquad0,out_ds,1);
            switch(coordim)
            {
            case 2:

                    Array<OneD, Array<OneD, NekDouble> > tangents;
                    tangents = Array<OneD, Array<OneD, NekDouble> >(coordim);
                    for(int k=0; k<coordim; ++k)
                    {
                        tangents[k]= Array<OneD, NekDouble>(nquad0); 
                    }
                    tangents = GetMetricInfo()->GetEdgeTangent();
                    ASSERTL0(tangents!=NullNekDoubleArrayofArray, 
                            "tangent vectors do not exist: check if a boundary region is defined as I ");
                    //diff= dU/de
                    Array<OneD,NekDouble> diff(nquad0);

                    PhysTensorDeriv(inarray,diff);

                    //get dS/de= (Jac)^-1
                    Array<OneD, NekDouble> Jac = m_metricinfo->GetJac();
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                         //calculate the derivative as (dU/de)*(Jac)^-1
                         Vmath::Vdiv(nquad0,diff,1,Jac ,1,out_ds,1);
                    }
                    else
                    {
                          NekDouble invJac = 1/Jac[0];
                          Vmath::Smul(nquad0, invJac,diff,1,out_ds,1);
                    }

            }
        }

        /**
           *\brief Evaluate the derivative normal to a line: 
           * \f$ d/dn=\frac{spacedim}{||normal||}d/d{\xi}  \f$.
           * The derivative is calculated performing
           *the product \f$ du/d{s}=\nabla u \cdot normal \f$.
           *\param inarray function to derive
           *\param out_dn result of the derivative operation 
        **/
        void SegExp::v_PhysDeriv_n(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD, NekDouble>& out_dn)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            Array<TwoD, const NekDouble>  gmat = m_metricinfo->GetGmat();
            int     coordim  = m_geom->GetCoordim();
            Array<OneD, NekDouble> out_dn_tmp(nquad0,0.0);
            switch(coordim)
            {
            case 2:

                    Array<OneD, NekDouble> inarray_d0(nquad0);
                    Array<OneD, NekDouble> inarray_d1(nquad0);

                    v_PhysDeriv(inarray,inarray_d0,inarray_d1);
                    Array<OneD, Array<OneD, NekDouble> > normals;
                    normals = Array<OneD, Array<OneD, NekDouble> >(coordim);
cout<<"der_n"<<endl;
                    for(int k=0; k<coordim; ++k)
                    {
                        normals[k]= Array<OneD, NekDouble>(nquad0);
                    }
// @TODO: this routine no longer makes sense, since normals are not unique on an edge
//                    normals = GetMetricInfo()->GetNormal();
                    for(int i=0; i<nquad0; i++)
                    {
cout<<"nx= "<<normals[0][i]<<"  ny="<<normals[1][i]<<endl;
                    }
                    ASSERTL0(normals!=NullNekDoubleArrayofArray, 
                            "normal vectors do not exist: check if a boundary region is defined as I ");
                    // \nabla u \cdot normal
                    Vmath::Vmul(nquad0,normals[0],1,inarray_d0,1,out_dn_tmp,1);
                    Vmath::Vadd(nquad0,out_dn_tmp,1,out_dn,1,out_dn,1);
                    Vmath::Zero(nquad0,out_dn_tmp,1);
                    Vmath::Vmul(nquad0,normals[1],1,inarray_d1,1,out_dn_tmp,1);
                    Vmath::Vadd(nquad0,out_dn_tmp,1,out_dn,1,out_dn,1);

                    for(int i=0; i<nquad0; i++)
                    {
cout<<"deps/dx ="<<inarray_d0[i]<<"  deps/dy="<<inarray_d1[i]<<endl;
                    }


            }

        }
        void SegExp::v_PhysDeriv(const int dir, 
                               const Array<OneD, const NekDouble>& inarray,
                               Array<OneD, NekDouble> &outarray)
        {
            switch(dir)
            {
            case 0:
                {
                    PhysDeriv(inarray, outarray, NullNekDouble1DArray, NullNekDouble1DArray);
                }
                break;
            case 1:
                {
                    PhysDeriv(inarray, NullNekDouble1DArray, outarray, NullNekDouble1DArray);   
                }
                break;
            case 2:
                {
                    PhysDeriv(inarray, NullNekDouble1DArray, NullNekDouble1DArray, outarray);   
                }
                break;
            default:
                {
                    ASSERTL1(false,"input dir is out of range");
                }
                break;
            }
        }


        //-----------------------------
        // Transforms
        //-----------------------------

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
        void SegExp::v_FwdTrans(
                const Array<OneD, const NekDouble>& inarray, 
                      Array<OneD,NekDouble> &outarray)
        {
            if(m_base[0]->Collocation())
            {
                Vmath::Vcopy(m_ncoeffs, inarray, 1, outarray, 1);
            }
            else 
            {
                v_IProductWRTBase(inarray,outarray);

                // get Mass matrix inverse
                MatrixKey             masskey(StdRegions::eInvMass, DetExpansionType(),*this);
                DNekScalMatSharedPtr  matsys = m_matrixManager[masskey];

                // copy inarray in case inarray == outarray
                NekVector<NekDouble> in(m_ncoeffs,outarray,eCopy);
                NekVector<NekDouble> out(m_ncoeffs,outarray,eWrapper);

                out = (*matsys)*in;
            }
        }

        void SegExp::v_FwdTrans_BndConstrained(
                const Array<OneD, const NekDouble>& inarray,
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
                case LibUtilities::eGauss_Lagrange:
                case LibUtilities::eGLL_Lagrange:
                    {
                        offset = 1;
                    }
                    break;
                case LibUtilities::eModified_A:
                case LibUtilities::eModified_B:
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

                if(m_ncoeffs>2)
                {
                    //  ideally, we would like to have tmp0 to be replaced by
                    //  outarray (currently MassMatrixOp does not allow aliasing)
                    Array<OneD, NekDouble> tmp0(m_ncoeffs); 
                    Array<OneD, NekDouble> tmp1(m_ncoeffs);

                    StdRegions::StdMatrixKey  stdmasskey(StdRegions::eMass,DetExpansionType(),*this);
                    MassMatrixOp(outarray,tmp0,stdmasskey);
                    v_IProductWRTBase(inarray,tmp1);

                    Vmath::Vsub(m_ncoeffs, tmp1, 1, tmp0, 1, tmp1, 1);

                    // get Mass matrix inverse (only of interior DOF)
                    MatrixKey             masskey(StdRegions::eMass, DetExpansionType(),*this);
                    DNekScalMatSharedPtr  matsys = (m_staticCondMatrixManager[masskey])->GetBlock(1,1);

                    Blas::Dgemv('N',nInteriorDofs,nInteriorDofs, matsys->Scale(), 
                                &((matsys->GetOwnedMatrix())->GetPtr())[0],
                                nInteriorDofs,tmp1.get()+offset,1,0.0,outarray.get()+offset,1);
                }
            }
        }


        //-----------------------------
        // Inner product functions
        //-----------------------------

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
        void SegExp::v_IProductWRTBase(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD, NekDouble> &outarray)
        {
            v_IProductWRTBase(m_base[0]->GetBdata(),inarray,outarray,1);
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
        void SegExp::v_IProductWRTBase(
                const Array<OneD, const NekDouble>& base,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD, NekDouble> &outarray,
                      int coll_check)
        {
            int   nquad0 = m_base[0]->GetNumPoints();
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
            StdSegExp::v_IProductWRTBase(base,tmp,outarray,coll_check);
        }


        void SegExp::v_IProductWRTDerivBase(
                const int dir, 
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD, NekDouble> & outarray)
        {
            int    nquad = m_base[0]->GetNumPoints();
            const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();

            Array<OneD, NekDouble> tmp1(nquad);

            switch(dir)
            {
            case 0:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul(nquad,gmat[0],1,inarray,1,tmp1,1);
                    }
                    else
                    {
                        Vmath::Smul(nquad, gmat[0][0], inarray, 1, tmp1, 1);
                    }
                }
                break;
            case 1:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul(nquad,gmat[1],1,inarray,1,tmp1,1);
                    }
                    else
                    {
                        Vmath::Smul(nquad, gmat[1][0], inarray, 1, tmp1, 1);
                    }
                }
                break;
            case 2:
                {
                    ASSERTL1(m_geom->GetCoordim() == 3,"input dir is out of range");
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul(nquad,gmat[2],1,inarray,1,tmp1,1);
                    }
                    else
                    {
                        Vmath::Smul(nquad, gmat[2][0], inarray, 1, tmp1, 1);
                    }
                }
                break;
            default:
                {
                    ASSERTL1(false,"input dir is out of range");
                }
                break;
            }
            v_IProductWRTBase(m_base[0]->GetDbdata(),tmp1,outarray,1);
        }

        void SegExp::v_NormVectorIProductWRTBase(
            const Array<OneD, const NekDouble> &Fx,
            const Array<OneD, const NekDouble> &Fy,
                  Array<OneD,       NekDouble> &outarray)
        {
            int nq = m_base[0]->GetNumPoints();
            Array<OneD, NekDouble > Fn(nq);
//            cout << "I am segment " << GetGeom()->GetGlobalID() << endl;
//            cout << "I want edge " << GetLeftAdjacentElementEdge() << endl;
// @TODO: This routine no longer makes sense as a normal is not unique to an edge
            const Array<OneD, const Array<OneD, NekDouble> >
                 &normals = GetLeftAdjacentElementExp()->GetEdgeNormal(GetLeftAdjacentElementEdge());
            Vmath::Vmul (nq,&Fx[0],1,&normals[0][0], 1,&Fn[0],1);
            Vmath::Vvtvp(nq,&Fy[0],1,&normals[1][0],1,&Fn[0],1,&Fn[0],1);

            v_IProductWRTBase(Fn,outarray);
        }

        //-----------------------------
        // Evaluation functions
        //-----------------------------

        NekDouble SegExp::v_PhysEvaluate(const Array<OneD, const NekDouble>& coord)
        {
            return PhysEvaluate(coord,m_phys);
        }

        NekDouble SegExp::v_PhysEvaluate(
                const Array<OneD, const NekDouble>& coord,
                const Array<OneD, const NekDouble> &physvals)
        {
            Array<OneD,NekDouble> Lcoord = Array<OneD,NekDouble>(1);

            ASSERTL0(m_geom,"m_geom not defined");
            m_geom->GetLocCoords(coord,Lcoord);

            return StdSegExp::v_PhysEvaluate(Lcoord, physvals);
        }


        /// Returns the locations of the quadrature points in up to
        /// three-dimensions.
        void SegExp::v_GetCoords(
                    Array<OneD, NekDouble> &coords_0,
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
                    LibUtilities::Interp1D(CBasis->GetBasisKey(),&(m_geom->UpdatePhys(2))[0],
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
                else // LibUtilities::Interpolate to Expansion point distribution
                {
                    LibUtilities::Interp1D(CBasis->GetBasisKey(),&(m_geom->UpdatePhys(1))[0],
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
                    LibUtilities::Interp1D(CBasis->GetBasisKey(),&(m_geom->UpdatePhys(0))[0],
                             m_base[0]->GetBasisKey(),&coords_0[0]);
                }
                break;
            default:
                ASSERTL0(false,"Number of dimensions are greater than 2");
                break;
            }
        }


        void SegExp::v_GetCoord(
                const Array<OneD, const NekDouble>& Lcoords, 
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

        //-----------------------------
        // Helper functions
        //-----------------------------

        /// Writes out the physical space data to file.
        void SegExp::v_WriteToFile(
                    std::ofstream &outfile,
                    OutputFormat format,
                    const bool dumpVar,
                    std::string var)
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
                    outfile << ", "<< var << std::endl << std::endl;
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
                int i;
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
                    outfile << "SP(" << coords[0][i] << ", ";
                    outfile << coords[1][i] << ", " << coords[2][i] << ")";
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

        void SegExp::v_SetCoeffsToOrientation(StdRegions::Orientation dir)
        {
            v_SetCoeffsToOrientation(dir,m_coeffs,m_coeffs);
        }

        void SegExp::v_SetCoeffsToOrientation(
                StdRegions::Orientation dir,
                Array<OneD, const NekDouble> &inarray,
                Array<OneD, NekDouble> &outarray)
        {

            if(dir == StdRegions::eBackwards)
            {
                if(&inarray[0] != &outarray[0])
                {
                    Array<OneD,NekDouble> intmp (inarray);
                    ReverseCoeffsAndSign(intmp,outarray);
                }
                else
                {
                    ReverseCoeffsAndSign(inarray,outarray);
                }
            }
        }

        StdRegions::Orientation SegExp::v_GetPorient(int point)
        {
            return m_geom->GetPorient(point);
        }


        int SegExp::v_GetCoordim()
        {
            return m_geom->GetCoordim();
        }

        /// Returns a pointer to the GeomFactors object describing the
        /// metric information for the segment.
        const SpatialDomains::GeomFactorsSharedPtr& SegExp::v_GetMetricInfo() const
        {
            return m_metricinfo;
        }


        const Array<OneD, const NekDouble>& SegExp::v_GetPhysNormals(void)
        {
            NEKERROR(ErrorUtil::efatal, "Got to SegExp");
            return NullNekDouble1DArray;
        }

        SpatialDomains::GeomType SegExp::v_MetricInfoType()
        {
            return m_metricinfo->GetGtype();
        }



        /// Returns a pointer to a Geometry object describing the
        /// geometry of the segment.
        const SpatialDomains::GeometrySharedPtr SegExp::v_GetGeom() const
        {
            return m_geom;
        }

        /// Returns a pointer to a Geometry1D object describing the
        /// geometry of the segment.
        const SpatialDomains::Geometry1DSharedPtr& SegExp::v_GetGeom1D() const
        {
            return m_geom;
        }

        /** \brief Virtual function to evaluate the discrete \f$ L_\infty\f$
            error \f$ |\epsilon|_\infty = \max |u - u_{exact}|\f$ where \f$
            u_{exact}\f$ is given by the array \a sol.

            The full function is defined in StdExpansion::Linf

            Input:

            - \a _phys: takes the physical value space array as
            approximate solution

            - \a sol: array of solution function  at physical quadrature points

            output:

            - returns the \f$ L_\infty \f$ error as a NekDouble.
        */
        NekDouble SegExp::v_Linf(const Array<OneD, const NekDouble>& sol)
        {
            return Linf(sol);
        }

        /** \brief Virtual function to evaluate the \f$ L_\infty \f$ norm of
            the function defined at the physical points \a (this)->_phys.

            The full function is defined in StdExpansion::Linf

            Input:

            - \a _phys: uses the physical value space array as discrete
            function to be evaulated.

            output:

            - returns the \f$ L_\infty \f$  as a NekDouble.
        */
        NekDouble SegExp::v_Linf()
        {
            return Linf();
        }

        /** \brief Virtual function to evaluate the \f$ L_2\f$, \f$ |
            \epsilon |_{2} = \left [ \int^1_{-1} [u - u_{exact}]^2 dx
            \right]^{1/2} d\xi_1 \f$ where \f$ u_{exact}\f$ is given by the
            array sol.

            The full function is defined in StdExpansion::L2

            Input:

            - \a _phys: takes the physical value space array as
            approximate solution
            - \a sol: array of solution function  at physical quadrature points

            output:

            - returns the \f$ L_2 \f$ error as a NekDouble.
        */
        NekDouble SegExp::v_L2(const Array<OneD, const NekDouble>& sol)
        {
            return StdExpansion::L2(sol);
        }

        /** \brief Virtual function to evaluate the \f$ L_2\f$ norm of the
            function defined at the physical points \a (this)->_phys.

            The full function is defined in StdExpansion::L2

            Input:

            - \a _phys: uses the physical value space array as discrete
            function to be evaulated.

            output:

            - returns the \f$ L_2 \f$  as a NekDouble.
        */
        NekDouble SegExp::v_L2()
        {
            return StdExpansion::L2();
        }



        int SegExp::v_GetNumPoints(const int dir) const
        {
            return GetNumPoints(dir);
        }


        int SegExp::v_GetNcoeffs(void) const
        {
            return m_ncoeffs;
        }

        const LibUtilities::BasisSharedPtr& SegExp::v_GetBasis(int dir) const
        {
            return GetBasis(dir);
        }


        int SegExp::v_NumBndryCoeffs() const
        {
            return 2;
        }


        int SegExp::v_NumDGBndryCoeffs() const
        {
            return 2;
        }


        void SegExp::v_SetUpPhysTangents(
                const StdRegions::StdExpansionSharedPtr &exp2D,
                const int edge)
        {
             GetMetricInfo()->ComputeEdgeTangents(exp2D->GetGeom(),
                                edge, GetBasis(0)->GetPointsKey());
        }


        /// Unpack data from input file assuming it comes from the same expansion type
        void SegExp::v_ExtractDataToCoeffs(
                const NekDouble *data, 
                const std::vector<unsigned int > &nummodes, 
                const int mode_offset,
                NekDouble *coeffs)
        {
            switch(m_base[0]->GetBasisType())
            { 
            case LibUtilities::eModified_A:
                {
                    int fillorder = min((int) nummodes[mode_offset],m_ncoeffs);

                    Vmath::Zero(m_ncoeffs,coeffs,1);
                    Vmath::Vcopy(fillorder,&data[0],1,&coeffs[0],1);
                }
                break;
            case LibUtilities::eGLL_Lagrange:
                {
                    // Assume that input is also Gll_Lagrange but no way to check;
                    LibUtilities::PointsKey p0(nummodes[mode_offset],LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::Interp1D(p0,data, m_base[0]->GetPointsKey(), coeffs);
                }
                break;
            case LibUtilities::eGauss_Lagrange:
                {
                    // Assume that input is also Gauss_Lagrange but no way to check;
                    LibUtilities::PointsKey p0(nummodes[mode_offset],LibUtilities::eGaussGaussLegendre);
                    LibUtilities::Interp1D(p0,data, m_base[0]->GetPointsKey(), coeffs);
                }
                break;
            default:
                ASSERTL0(false,"basis is either not set up or not hierarchicial");
            }
        }

        void SegExp::v_ComputeVertexNormal(const int vertex)
        {
            int i;
            const SpatialDomains::GeomFactorsSharedPtr &geomFactors = GetGeom()->GetMetricInfo();
            SpatialDomains::GeomType type = geomFactors->GetGtype();
            const Array<TwoD, const NekDouble> &gmat = geomFactors->GetGmat();
            int nqe = m_base[0]->GetNumPoints();
            int vCoordDim = GetCoordim();

            m_vertexNormals[vertex] = Array<OneD, Array<OneD, NekDouble> >(vCoordDim);
            Array<OneD, Array<OneD, NekDouble> > &normal = m_vertexNormals[vertex];
            for (i = 0; i < vCoordDim; ++i)
            {
                normal[i] = Array<OneD, NekDouble>(nqe);
            }

            // Regular geometry case
            if((type == SpatialDomains::eRegular)||(type == SpatialDomains::eMovingRegular))
            {
                NekDouble vert;
                // Set up normals
                switch(vertex)
                {
                    case 0:
                        for(i = 0; i < vCoordDim; ++i)
                        {
                            Vmath::Fill(nqe,-gmat[i][0],normal[i],1);
                        }
                        break;
                    case 1:
                        for(i = 0; i < vCoordDim; ++i)
                        {
                            Vmath::Fill(nqe,gmat[i][0],normal[i],1);
                        }
                        break;
                    default:
                        ASSERTL0(false,"point is out of range (point < 2)");
                }

                // normalise
                vert = 0.0;
                for(i =0 ; i < vCoordDim; ++i)
                {
                    vert += normal[i][0]*normal[i][0];
                }
                vert = 1.0/sqrt(vert);
                for (i = 0; i < vCoordDim; ++i)
                {
                    Vmath::Smul(nqe,vert,normal[i],1,normal[i],1);
                }
            }
        }

        //-----------------------------
        // Operator creation functions
        //-----------------------------


        void SegExp::v_LaplacianMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdRegions::StdMatrixKey     &mkey)
        {
            int    nquad = m_base[0]->GetNumPoints();
            const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();

            Array<OneD,NekDouble> physValues(nquad);
            Array<OneD,NekDouble> dPhysValuesdx(nquad);

            BwdTrans(inarray,physValues);

            // Laplacian matrix operation
            switch(m_geom->GetCoordim())
            {
            case 1:
                {
                    PhysDeriv(physValues,dPhysValuesdx);

                    // multiply with the proper geometric factors
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul(nquad,&gmat[0][0],1,dPhysValuesdx.get(),1,dPhysValuesdx.get(),1);
                    }
                    else
                    {
                        Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.get(), 1);
                    }
                }
                break;
            case 2:
                {
                    Array<OneD,NekDouble> dPhysValuesdy(nquad);

                    PhysDeriv(physValues,dPhysValuesdx,dPhysValuesdy);

                    // multiply with the proper geometric factors
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul (nquad,&gmat[0][0],1,dPhysValuesdx.get(),1,dPhysValuesdx.get(),1);
                        Vmath::Vvtvp(nquad,&gmat[1][0],1,dPhysValuesdy.get(),1,dPhysValuesdx.get(),1,
                                     dPhysValuesdx.get(),1);
                    }
                    else
                    {
                        Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.get(), 1);
                        Blas::Daxpy(nquad, gmat[1][0], dPhysValuesdy.get(), 1, dPhysValuesdx.get(), 1);
                    }
                }
                break;
            case 3:
                {
                    Array<OneD,NekDouble> dPhysValuesdy(nquad);
                    Array<OneD,NekDouble> dPhysValuesdz(nquad);

                    PhysDeriv(physValues,dPhysValuesdx,dPhysValuesdy,dPhysValuesdz); 

                    // multiply with the proper geometric factors
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul (nquad,&gmat[0][0],1,dPhysValuesdx.get(),1,dPhysValuesdx.get(),1);
                        Vmath::Vvtvp(nquad,&gmat[1][0],1,dPhysValuesdy.get(),1,dPhysValuesdx.get(),1,
                                     dPhysValuesdx.get(),1);
                        Vmath::Vvtvp(nquad,&gmat[2][0],1,dPhysValuesdz.get(),1,dPhysValuesdx.get(),1,
                                     dPhysValuesdx.get(),1);
                    }
                    else
                    {
                        Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.get(), 1);
                        Blas::Daxpy(nquad, gmat[1][0], dPhysValuesdy.get(), 1, dPhysValuesdx.get(), 1);
                        Blas::Daxpy(nquad, gmat[2][0], dPhysValuesdz.get(), 1, dPhysValuesdx.get(), 1);
                    }  
                }
                break;
            default:
                ASSERTL0(false,"Wrong number of dimensions");
                break;
            }

            v_IProductWRTBase(m_base[0]->GetDbdata(),dPhysValuesdx,outarray,1);
        }


        void SegExp::v_HelmholtzMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdRegions::StdMatrixKey     &mkey)
        {
            int    nquad = m_base[0]->GetNumPoints();
            const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();
            const NekDouble lambda = mkey.GetConstFactor(StdRegions::eFactorLambda);

            Array<OneD,NekDouble> physValues(nquad);
            Array<OneD,NekDouble> dPhysValuesdx(nquad);
            Array<OneD,NekDouble> wsp(m_ncoeffs);

            BwdTrans(inarray,physValues);

            // mass matrix operation
            v_IProductWRTBase((m_base[0]->GetBdata()),physValues,wsp,1);

            // Laplacian matrix operation
            switch(m_geom->GetCoordim())
            {
            case 1:
                {
                    PhysDeriv(physValues,dPhysValuesdx);

                    // multiply with the proper geometric factors
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul(nquad,&gmat[0][0],1,dPhysValuesdx.get(),1,dPhysValuesdx.get(),1);
                    }
                    else
                    {
                        Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.get(), 1);
                    }
                }
                break;
            case 2:
                {
                    Array<OneD,NekDouble> dPhysValuesdy(nquad);

                    PhysDeriv(physValues,dPhysValuesdx,dPhysValuesdy);

                    // multiply with the proper geometric factors
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul (nquad,&gmat[0][0],1,dPhysValuesdx.get(),1,dPhysValuesdx.get(),1);
                        Vmath::Vvtvp(nquad,&gmat[1][0],1,dPhysValuesdy.get(),1,dPhysValuesdx.get(),1,
                                     dPhysValuesdx.get(),1);
                    }
                    else
                    {
                        Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.get(), 1);
                        Blas::Daxpy(nquad, gmat[1][0], dPhysValuesdy.get(), 1, dPhysValuesdx.get(), 1);
                    }
                }
                break;
            case 3:
                {
                    Array<OneD,NekDouble> dPhysValuesdy(nquad);
                    Array<OneD,NekDouble> dPhysValuesdz(nquad);

                    PhysDeriv(physValues,dPhysValuesdx,dPhysValuesdy,dPhysValuesdz); 

                    // multiply with the proper geometric factors
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul (nquad,&gmat[0][0],1,dPhysValuesdx.get(),1,dPhysValuesdx.get(),1);
                        Vmath::Vvtvp(nquad,&gmat[1][0],1,dPhysValuesdy.get(),1,dPhysValuesdx.get(),1,
                                     dPhysValuesdx.get(),1);
                        Vmath::Vvtvp(nquad,&gmat[2][0],1,dPhysValuesdz.get(),1,dPhysValuesdx.get(),1,
                                     dPhysValuesdx.get(),1);
                    }
                    else
                    {
                        Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.get(), 1);
                        Blas::Daxpy(nquad, gmat[1][0], dPhysValuesdy.get(), 1, dPhysValuesdx.get(), 1);
                        Blas::Daxpy(nquad, gmat[2][0], dPhysValuesdz.get(), 1, dPhysValuesdx.get(), 1);
                    }  
                }
                break;
            default:
                ASSERTL0(false,"Wrong number of dimensions");
                break;
            }

            v_IProductWRTBase(m_base[0]->GetDbdata(),dPhysValuesdx,outarray,1);
            Blas::Daxpy(m_ncoeffs, lambda, wsp.get(), 1, outarray.get(), 1);
        }

        //-----------------------------
        // Matrix creation functions
        //-----------------------------


        DNekScalBlkMatSharedPtr SegExp::v_GetLocStaticCondMatrix(
                const MatrixKey &mkey)
        {
            return m_staticCondMatrixManager[mkey];
        }

        DNekScalMatSharedPtr SegExp::v_GetLocMatrix(const MatrixKey &mkey)
        {
            return m_matrixManager[mkey];
        }


        DNekMatSharedPtr SegExp::v_CreateStdMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            LibUtilities::BasisKey bkey = m_base[0]->GetBasisKey();
            StdRegions::StdSegExpSharedPtr tmp = 
                        MemoryManager<StdSegExp>::AllocateSharedPtr(bkey);

            return tmp->GetStdMatrix(mkey); 
        }


        DNekScalMatSharedPtr SegExp::CreateMatrix(const MatrixKey &mkey)
        {
            DNekScalMatSharedPtr returnval;
            NekDouble fac;

            ASSERTL2(m_metricinfo->GetGtype() != 
                     SpatialDomains::eNoGeomType,"Geometric information is not set up");

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eMass:
                {
                    if((m_metricinfo->GetGtype() == SpatialDomains::eDeformed)||
                       (mkey.GetNVarCoeff()))
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
                    if((m_metricinfo->GetGtype() == SpatialDomains::eDeformed)||
                       (mkey.GetNVarCoeff()))
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
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed || 
                       mkey.GetNVarCoeff())
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
                                ASSERTL1(m_geom->GetCoordim() >= 2,
                                         "Cannot call eWeakDeriv2 in a "
                                         "coordinate system which is not at "
                                         "least two-dimensional");
                                dir = 1;
                                break;
                            case StdRegions::eWeakDeriv2:
                                ASSERTL1(m_geom->GetCoordim() == 3,
                                         "Cannot call eWeakDeriv2 in a "
                                         "coordinate system which is not "
                                         "three-dimensional");
                                dir = 2;
                                break;
                            default:
                                break;
                        }

                        MatrixKey deriv0key(StdRegions::eWeakDeriv0,
                                            mkey.GetExpansionType(), *this);  

                        DNekMatSharedPtr WeakDerivStd = GetStdMatrix(deriv0key);
                        fac = m_metricinfo->GetGmat()[dir][0]*
                            m_metricinfo->GetJac()[0];

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(fac,WeakDerivStd);
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
                    NekDouble factor = mkey.GetConstFactor(StdRegions::eFactorLambda);
                    MatrixKey masskey(StdRegions::eMass,
                                      mkey.GetExpansionType(), *this);    
                    DNekScalMat &MassMat = *(this->m_matrixManager[masskey]);
                    MatrixKey lapkey(StdRegions::eLaplacian,
                                     mkey.GetExpansionType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekScalMat &LapMat = *(this->m_matrixManager[lapkey]);

                    int rows = LapMat.GetRows();
                    int cols = LapMat.GetColumns();

                    DNekMatSharedPtr helm = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);

                    NekDouble one = 1.0;
                    (*helm) = LapMat + factor*MassMat;

                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,helm);
                }
                break;
            case StdRegions::eHybridDGHelmholtz:
            case StdRegions::eHybridDGLamToU:
            case StdRegions::eHybridDGLamToQ0:
            case StdRegions::eHybridDGHelmBndLam:
                {
                    NekDouble one    = 1.0;

                    DNekMatSharedPtr mat = GenMatrix(mkey);
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                }
                break;
            case StdRegions::eInvHybridDGHelmholtz:
                {
                    NekDouble one = 1.0;

//                    StdRegions::StdMatrixKey hkey(StdRegions::eHybridDGHelmholtz,
//                                                  DetExpansionType(),*this,
//                                                  mkey.GetConstant(0),
//                                                  mkey.GetConstant(1));
                    MatrixKey hkey(StdRegions::eHybridDGHelmholtz, DetExpansionType(),
                                   *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekMatSharedPtr mat = GenMatrix(hkey);

                    mat->Invert();
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                }
                break;
            UseLocRegionsMatrix:
                {
                    DNekMatSharedPtr mat = GenMatrix(mkey);
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(fac,mat);
                }
                break;
            UseStdRegionsMatrix:
                {
                    DNekMatSharedPtr mat = GetStdMatrix(mkey);
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(fac,mat);
                }
                break;
            default:
                NEKERROR(ErrorUtil::efatal, "Matrix creation not defined");
                break;
            }

            return returnval;
        }

        DNekMatSharedPtr SegExp::v_GenMatrix(const StdRegions::StdMatrixKey &mkey)
        {

            DNekMatSharedPtr returnval;

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eHybridDGHelmholtz:
            case StdRegions::eHybridDGLamToU:
            case StdRegions::eHybridDGLamToQ0:
            case StdRegions::eHybridDGLamToQ1:
            case StdRegions::eHybridDGLamToQ2:
            case StdRegions::eHybridDGHelmBndLam:
                returnval = Expansion1D::v_GenMatrix(mkey);
                break;
            default:
                returnval = StdSegExp::v_GenMatrix(mkey);
                break;
            }

            return returnval;
        }


        DNekScalBlkMatSharedPtr SegExp::CreateStaticCondMatrix(const MatrixKey &mkey)
        {
            DNekScalBlkMatSharedPtr returnval;

            ASSERTL2(m_metricinfo->GetGtype() != SpatialDomains::eNoGeomType,"Geometric information is not set up");

            // set up block matrix system
            int nbdry = NumBndryCoeffs();
            int nint = m_ncoeffs - nbdry;
            Array<OneD, unsigned int> exp_size(2);
            exp_size[0] = nbdry;
            exp_size[1] = nint;

            /// \todo Really need a constructor which takes Arrays
            returnval = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(exp_size,exp_size);
            NekDouble factor = 1.0;

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eLaplacian: 
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
                    DNekBlkMatSharedPtr  mat = GetStdStaticCondMatrix(mkey);
                    DNekScalMatSharedPtr Atmp;
                    DNekMatSharedPtr     Asubmat;

                    returnval->SetBlock(0,0,Atmp = 
                        MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,Asubmat = mat->GetBlock(0,0)));
                    returnval->SetBlock(0,1,Atmp =
                        MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Asubmat = mat->GetBlock(0,1)));
                    returnval->SetBlock(1,0,Atmp =
                        MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,Asubmat = mat->GetBlock(1,0)));
                    returnval->SetBlock(1,1,Atmp =
                        MemoryManager<DNekScalMat>::AllocateSharedPtr(invfactor,Asubmat = mat->GetBlock(1,1)));
                }
                break;
            UseLocRegionsMatrix:
                {
                    int i,j;
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


        //-----------------------------
        // Private methods
        //-----------------------------


        ///  Reverse the coefficients in a boundary interior expansion
        ///  this routine is of use when we need the segment
        ///  coefficients corresponding to a expansion in the reverse
        ///  coordinate direction
        void SegExp::ReverseCoeffsAndSign(
                const Array<OneD,NekDouble> &inarray,
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
                break;
            default:
                ASSERTL0(false,"This basis is not allowed in this method");
                break;
            }
        }

        /* \brief Mass inversion product from \a inarray to \a outarray
         *
         *     Multiply by the inverse of the mass matrix
         *     \f$ {\bf \hat{u}} = {\bf M}^{-1} {\bf I} \f$
         *
        **/
        void SegExp::MultiplyByElmtInvMass(
                const Array<OneD, const NekDouble>& inarray, 
                      Array<OneD,NekDouble> &outarray)
        {
            // get Mass matrix inverse
            MatrixKey             masskey(StdRegions::eInvMass, DetExpansionType(),*this);
            DNekScalMatSharedPtr  matsys = m_matrixManager[masskey];

            NekVector<NekDouble> in(m_ncoeffs,inarray,eCopy);
            NekVector<NekDouble> out(m_ncoeffs,outarray,eWrapper);

            out = (*matsys)*in;
        }


    } // end of namespace
}//end of namespace

