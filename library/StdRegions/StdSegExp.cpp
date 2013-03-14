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

        /** \brief Default constructor */

        StdSegExp::StdSegExp()
        {
        }


        /** \brief Constructor using BasisKey class for quadrature points and 
         *  order definition
         *
         *  \param Ba BasisKey class definition containing order and quadrature 
         *  points.
         */

        StdSegExp::StdSegExp(const LibUtilities::BasisKey &Ba):
        StdExpansion(Ba.GetNumModes(), 1, Ba),
        StdExpansion1D(Ba.GetNumModes(),Ba)
        {
        }


        /** \brief Copy Constructor */

        StdSegExp::StdSegExp(const StdSegExp &T):
                StdExpansion(T),
                StdExpansion1D(T)
        {
        }


        StdSegExp::~StdSegExp()
        {
        }

        /** \brief Return Shape of region, using  ShapeType enum list.
         *  i.e. Segment
         */
        LibUtilities::ShapeType StdSegExp::v_DetShapeType() const
        {
            return LibUtilities::eSegment;
        }

        bool StdSegExp::v_IsBoundaryInteriorExpansion()
        {

            bool returnval = false;

            if(m_base[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                returnval = true;
            }

            if(m_base[0]->GetBasisType() == LibUtilities::eGLL_Lagrange)
            {
                returnval = true;
            }

            return returnval;
        }




        //---------------------------------------------------------------------
        // Integration Methods
        //---------------------------------------------------------------------

        /** \brief Integrate the physical point list \a inarray over region 
         *  and return the value
         *
         *  \param inarray definition of function to be integrated evauluated at 
         *  quadrature point of expansion. 
         *  \return returns \f$\int^1_{-1} u(\xi_1)d \xi_1 \f$ where \f$inarray[i]
         *  = u(\xi_{1i}) \f$
         */
        NekDouble StdSegExp::v_Integral(const Array<OneD, const NekDouble>& inarray )
        {
            NekDouble Int = 0.0;
            int    nquad0 = m_base[0]->GetNumPoints();
            Array<OneD, NekDouble> tmp(nquad0);
            Array<OneD, const NekDouble> z  = m_base[0]->GetZ();
            Array<OneD, const NekDouble> w0 = m_base[0]->GetW();

            // multiply by integration constants 
            Vmath::Vmul(nquad0, inarray, 1, w0, 1, tmp, 1);

            Int = Vmath::Vsum(nquad0, tmp, 1);

            return Int;
        }




        //---------------------------------------------------------------------
        // Differentiation Methods
        //---------------------------------------------------------------------


        /** \brief Evaluate the derivative \f$ d/d{\xi_1} \f$ at the physical 
         *  quadrature points given by \a inarray and return in \a outarray.
         *
         *  This is a wrapper around StdExpansion1D::Tensor_Deriv
         *  \param inarray array of a function evaluated at the quadrature points
         *  \param  outarray the resulting array of the derivative \f$
         *  du/d_{\xi_1}|_{\xi_{1i}} \f$ will be stored in the array \a outarra 
         */

        void StdSegExp::v_PhysDeriv(const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> &out_d0,
                Array<OneD, NekDouble> &out_d1,
                Array<OneD, NekDouble> &out_d2)
        {
            PhysTensorDeriv(inarray,out_d0);
        }


        void StdSegExp::v_PhysDeriv(const int dir,
                const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray)
        {
            ASSERTL1(dir==0,"input dir is out of range");
            PhysTensorDeriv(inarray,outarray);
            // PhysDeriv(inarray, outarray);
        }

        void StdSegExp::v_StdPhysDeriv(
                const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> &out_d0,
                Array<OneD, NekDouble> &out_d1,
                Array<OneD, NekDouble> &out_d2)
        {
            PhysTensorDeriv(inarray,out_d0);
            // PhysDeriv(inarray, out_d0);
        }

        void StdSegExp::v_StdPhysDeriv(
                const int dir,
                const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray)
        {
            ASSERTL1(dir==0,"input dir is out of range");
            PhysTensorDeriv(inarray,outarray);
            // PhysDeriv(inarray, outarray);
        }



        //---------------------------------------------------------------------
        // Transforms
        //---------------------------------------------------------------------

        /** \brief Backward transform from coefficient space given
         *  in \a inarray and evaluate at the physical quadrature
         *  points \a outarray
         *
         *  Operation can be evaluated as \f$ u(\xi_{1i}) =
         *  \sum_{p=0}^{order-1} \hat{u}_p \phi_p(\xi_{1i}) \f$ or equivalently 
         *  \f$ {\bf u} = {\bf B}^T {\bf \hat{u}} \f$ where
         *  \f${\bf B}[i][j] = \phi_i(\xi_{1j}), \mbox{\_coeffs}[p] = {\bf
         *  \hat{u}}[p] \f$
         *
         *  The function takes the coefficient array \a inarray as
         *  input for the transformation
         *
         *  \param inarray: the coeffficients of the expansion 
         *
         *  \param outarray: the resulting array of the values of the function at 
         *  the physical quadrature points will be stored in the array \a outarray
         */

        void StdSegExp::v_BwdTrans(
                const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray)
        {
            int  nquad = m_base[0]->GetNumPoints();

            if(m_base[0]->Collocation())
            {
                Vmath::Vcopy(nquad, inarray, 1, outarray, 1);
            }
            else
            {

#ifdef NEKTAR_USING_DIRECT_BLAS_CALLS

                Blas::Dgemv('N',nquad,m_base[0]->GetNumModes(),1.0, (m_base[0]->GetBdata()).get(),
                            nquad,&inarray[0],1,0.0,&outarray[0],1);

#else //NEKTAR_USING_DIRECT_BLAS_CALLS

                NekVector<NekDouble> in(m_ncoeffs,inarray,eWrapper);
                NekVector<NekDouble> out(nquad,outarray,eWrapper);
                NekMatrix<NekDouble> B(nquad,m_ncoeffs,m_base[0]->GetBdata(),eWrapper);
                out = B * in;

#endif //NEKTAR_USING_DIRECT_BLAS_CALLS 

            }
        }

        /** \brief Forward transform from physical quadrature space stored in 
         *  \a inarray and evaluate the expansion coefficients and store in 
         *  \a outarray
         *
         *  Perform a forward transform using a Galerkin projection by taking the 
         *  inner product of the physical points and multiplying by the inverse of
         *  the mass matrix using the Solve method of the standard matrix 
         *  container holding the local mass matrix, i.e. \f$ {\bf \hat{u}} = 
         *  {\bf M}^{-1} {\bf I} \f$ where \f$ {\bf I}[p] =  \int^1_{-1} 
         *  \phi_p(\xi_1) u(\xi_1) d\xi_1 \f$
         *
         *  This function stores the expansion coefficients calculated by the 
         *  transformation in the coefficient space array \a outarray
         *
         *  \param inarray: array of physical quadrature points to be transformed
         *
         *  \param outarray: the coeffficients of the expansion 
         */ 

        void StdSegExp::v_FwdTrans(const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray)
        {
            if(m_base[0]->Collocation())
            {
                Vmath::Vcopy(m_ncoeffs, inarray, 1, outarray, 1);
            }
            else
            {
                v_IProductWRTBase(inarray,outarray);

                // get Mass matrix inverse
                StdMatrixKey      masskey(eInvMass,v_DetShapeType(),*this);
                DNekMatSharedPtr  matsys = GetStdMatrix(masskey);

                NekVector<NekDouble> in(m_ncoeffs,outarray,eCopy);
                NekVector<NekDouble> out(m_ncoeffs,outarray,eWrapper);

                out = (*matsys)*in;
            }
        }

        void StdSegExp::v_FwdTrans_BndConstrained(
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
                    // ideally, we would like to have tmp0 to be replaced by
                    // outarray (currently MassMatrixOp does not allow aliasing)
                    Array<OneD, NekDouble> tmp0(m_ncoeffs);
                    Array<OneD, NekDouble> tmp1(m_ncoeffs);

                    StdMatrixKey      masskey(eMass,v_DetShapeType(),*this);
                    MassMatrixOp(outarray,tmp0,masskey);
                    v_IProductWRTBase(inarray,tmp1);

                    Vmath::Vsub(m_ncoeffs, tmp1, 1, tmp0, 1, tmp1, 1);

                    // get Mass matrix inverse (only of interior DOF)
                    DNekMatSharedPtr  matsys = (m_stdStaticCondMatrixManager[masskey])->GetBlock(1,1);

                    Blas::Dgemv('N',nInteriorDofs,nInteriorDofs,1.0, &(matsys->GetPtr())[0],
                                nInteriorDofs,tmp1.get()+offset,1,0.0,outarray.get()+offset,1);
                }
            }

        }


        void StdSegExp::v_BwdTrans_SumFac(const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray)
        {
            v_BwdTrans(inarray, outarray);
        }



        //---------------------------------------------------------------------
        // Inner product functions
        //---------------------------------------------------------------------



        /** \brief  Inner product of \a inarray over region with respect to 
         *  expansion basis \a base and return in \a outarray 
         *
         *  Calculate \f$ I[p] = \int^{1}_{-1} \phi_p(\xi_1) u(\xi_1) d\xi_1
         *  = \sum_{i=0}^{nq-1} \phi_p(\xi_{1i}) u(\xi_{1i}) w_i \f$ where
         *  \f$ outarray[p] = I[p], inarray[i] = u(\xi_{1i}), base[p*nq+i] =
         *  \phi_p(\xi_{1i}) \f$.
         *
         *  \param  base an array defining the local basis for the inner product 
         *  usually passed from Basis->GetBdata() or Basis->GetDbdata()
         *  \param inarray: physical point array of function to be integrated
         *  \f$ u(\xi_1) \f$
         *  \param coll_check flag to identify when a Basis->Collocation() call 
         *  should be performed to see if this is a GLL_Lagrange basis with a 
         *  collocation property. (should be set to 0 if taking the inner  
         *  product with respect to the derivative of basis)
         *  \param outarray  the values of the inner product with respect to 
         *  each basis over region will be stored in the array \a outarray as
         *  output of the function
         */

        void StdSegExp::v_IProductWRTBase(
            const Array<OneD, const NekDouble> &base,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            int coll_check)
        {
            int    nquad = m_base[0]->GetNumPoints();
            Array<OneD, NekDouble> tmp(nquad);
            Array<OneD, const NekDouble> z =  m_base[0]->GetZ();
            Array<OneD, const NekDouble> w =  m_base[0]->GetW();

            Vmath::Vmul(nquad, inarray, 1, w, 1, tmp, 1);
            
            /* Comment below was a bug for collocated basis
            if(coll_check&&m_base[0]->Collocation())
            {
                Vmath::Vcopy(nquad, tmp, 1, outarray, 1);
            }
            else
            {
                Blas::Dgemv('T',nquad,m_ncoeffs,1.0,base.get(),nquad,
                            &tmp[0],1,0.0,outarray.get(),1);
            }*/
            
            // Correct implementation
            Blas::Dgemv('T',nquad,m_ncoeffs,1.0,base.get(),nquad,
                        &tmp[0],1,0.0,outarray.get(),1);
        }

        /** \brief Inner product of \a inarray over region with respect to the
         *  expansion basis (this)->m_base[0] and return in \a outarray
         *
         *  Wrapper call to StdSegExp::IProductWRTBase
         *  \param inarray array of function values evaluated at the physical
         *  collocation points
         *  \param outarray  the values of the inner product with respect to 
         *  each basis over region will be stored in the array \a outarray as
         *  output of the function
         */
        void StdSegExp::v_IProductWRTBase(const Array<OneD, const NekDouble>& inarray, 
                Array<OneD, NekDouble> &outarray)
        {
            v_IProductWRTBase(m_base[0]->GetBdata(),inarray,outarray,1);
        }

        void StdSegExp::v_IProductWRTDerivBase(
                const int dir,
                const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> & outarray)
      {
            ASSERTL1(dir >= 0 && dir < 1,"input dir is out of range");
            v_IProductWRTBase(m_base[0]->GetDbdata(),inarray,outarray,1);
        }

        void StdSegExp::v_IProductWRTBase_SumFac(
                const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray)
        {
            v_IProductWRTBase(m_base[0]->GetBdata(),inarray,outarray,1);
        }





        //----------------------------
        // Evaluation
        //----------------------------

        void StdSegExp::v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
        {
            int   nquad = m_base[0]->GetNumPoints();
            const NekDouble * base  = m_base[0]->GetBdata().get();

            ASSERTL2(mode <= m_ncoeffs,
             "calling argument mode is larger than total expansion order");

            Vmath::Vcopy(nquad,(NekDouble *)base+mode*nquad,1, &outarray[0],1);
        }

        NekDouble StdSegExp::v_PhysEvaluate(
                const Array<OneD, const NekDouble>& coords)
        {
            return  StdExpansion1D::v_PhysEvaluate(coords, m_phys);
        }

        NekDouble StdSegExp::v_PhysEvaluate(
                const Array<OneD, const NekDouble>& coords,
                const Array<OneD, const NekDouble>& physvals)
        {
            return  StdExpansion1D::v_PhysEvaluate(coords, physvals);
        }

        void StdSegExp::v_LaplacianMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
        {
            int    nquad = m_base[0]->GetNumPoints();

            Array<OneD,NekDouble> physValues(nquad);
            Array<OneD,NekDouble> dPhysValuesdx(nquad);

            v_BwdTrans(inarray,physValues);

            // Laplacian matrix operation
            v_PhysDeriv(physValues,dPhysValuesdx);
            v_IProductWRTBase(m_base[0]->GetDbdata(),dPhysValuesdx,outarray,1);
        }


        void StdSegExp::v_HelmholtzMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
        {
            int    nquad = m_base[0]->GetNumPoints();

            Array<OneD,NekDouble> physValues(nquad);
            Array<OneD,NekDouble> dPhysValuesdx(nquad);
            Array<OneD,NekDouble> wsp(m_ncoeffs);

            v_BwdTrans(inarray,physValues);

            // mass matrix operation
            v_IProductWRTBase((m_base[0]->GetBdata()),physValues,wsp,1);

            // Laplacian matrix operation
            v_PhysDeriv(physValues,dPhysValuesdx);
            v_IProductWRTBase(m_base[0]->GetDbdata(),dPhysValuesdx,outarray,1);
            Blas::Daxpy(m_ncoeffs, mkey.GetConstFactor(eFactorLambda), wsp.get(), 1, outarray.get(), 1);
        }


        void StdSegExp::v_GetCoords(
                Array<OneD, NekDouble> &coords_0,
                Array<OneD, NekDouble> &coords_1,
                Array<OneD, NekDouble> &coords_2)
        {
            Blas::Dcopy(GetNumPoints(0),(m_base[0]->GetZ()).get(),
                        1,&coords_0[0],1);
        }





        //---------------------------------------------------------------------
        // Helper functions
        //---------------------------------------------------------------------


        void StdSegExp::v_WriteToFile(
                std::ofstream &outfile,
                OutputFormat format,
                const bool dumpVar,
                std::string var)
        {
            if(format==eTecplot)
            {
                int i;
                int     nquad = m_base[0]->GetNumPoints();
                Array<OneD, const NekDouble> z = m_base[0]->GetZ();

                if(dumpVar)
                {
                    outfile << "Variables = z";   
                    outfile << ", "<< var << std::endl << std::endl;
                }

                outfile << "Zone, I=" << nquad <<", F=Point" << std::endl;

                for(i = 0; i < nquad; ++i)
                {
                    outfile << z[i] << " " << m_phys[i] << std::endl;
                }
            }
            else if(format==eGmsh)
            {
                int i;
                int     nquad = m_base[0]->GetNumPoints();
                Array<OneD, const NekDouble> z = m_base[0]->GetZ();

                if(dumpVar)
                {
                    outfile<<"View.Type = 2;"<<endl;
                    outfile<<"View \" \" {"<<endl;
                }
 
                for(i = 0; i < nquad; ++i)
                {
                    outfile << "SP(" << z[i] <<",  0.0, 0.0){" << m_phys[i] << "};" << endl;
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

        int StdSegExp::v_GetNverts() const
        {
            return 2;
        } 

        int StdSegExp::v_NumBndryCoeffs() const
        {
            return 2;
        } 

        int StdSegExp::v_NumDGBndryCoeffs() const
        {
            return 2;
        } 

        int StdSegExp::v_CalcNumberOfCoefficients(
                const std::vector<unsigned int> &nummodes,
                int &modes_offset)
        {
            int nmodes = nummodes[modes_offset];
            modes_offset += 1;

            return nmodes;
        }

        //---------------------------------------------------------------------
        // Wrapper functions
        //---------------------------------------------------------------------

        DNekMatSharedPtr StdSegExp::v_GenMatrix(const StdMatrixKey &mkey)
        {
            DNekMatSharedPtr Mat;
            MatrixType mattype;

            switch(mattype = mkey.GetMatrixType())
            {
            case eFwdTrans:
                {
                    Mat = MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs,m_ncoeffs);
                    StdMatrixKey iprodkey(eIProductWRTBase,v_DetShapeType(),*this);
                    DNekMat &Iprod = *GetStdMatrix(iprodkey);
                    StdMatrixKey imasskey(eInvMass,v_DetShapeType(),*this);
                    DNekMat &Imass = *GetStdMatrix(imasskey);

                    (*Mat) = Imass*Iprod;
					
					
					
                }
                break;
            default:
                {
                    Mat = StdExpansion::CreateGeneralMatrix(mkey);

                    if(mattype ==  eMass)
                    {
                        // For Fourier basis set the imaginary component
                        // of mean mode to have a unit diagonal component
                        // in mass matrix
                        if(m_base[0]->GetBasisType() == LibUtilities::eFourier)
                        {
                            (*Mat)(1,1) = 1.0;
                        }
                    }
                }
                break;
            }

            return Mat;
        }


        DNekMatSharedPtr StdSegExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
        {
            return v_GenMatrix(mkey);
        }

        //---------------------------------------------------------------------
        // Mappings
        //---------------------------------------------------------------------


        void StdSegExp::v_GetBoundaryMap(Array<OneD, unsigned int>& outarray)
        {
            if(outarray.num_elements() != NumBndryCoeffs())
            {
                outarray = Array<OneD, unsigned int>(NumBndryCoeffs());
            }
            const LibUtilities::BasisType Btype = GetBasisType(0);
            int nummodes = m_base[0]->GetNumModes();

            outarray[0] = 0;

            switch(Btype)
            {
            case LibUtilities::eGLL_Lagrange:
            case LibUtilities::eGauss_Lagrange:
            case LibUtilities::eChebyshev:
            case LibUtilities::eFourier:
                outarray[1]= nummodes-1;
                break;
            case LibUtilities::eModified_A:
            case LibUtilities::eModified_B:
                outarray[1] = 1;
                break;
            default:
                ASSERTL0(0,"Mapping array is not defined for this expansion");
                break;
            }
        }

        void StdSegExp::v_GetInteriorMap(Array<OneD, unsigned int>& outarray)
        {
            int i;
            if(outarray.num_elements()!=GetNcoeffs()-NumBndryCoeffs())
            {
                outarray = Array<OneD, unsigned int>(GetNcoeffs()-NumBndryCoeffs());
            }
            const LibUtilities::BasisType Btype = GetBasisType(0);

            switch(Btype)
            {
            case LibUtilities::eGLL_Lagrange:
            case LibUtilities::eGauss_Lagrange:
            case LibUtilities::eChebyshev:
            case LibUtilities::eFourier:
                for(i = 0 ; i < GetNcoeffs()-2;i++)
                {
                    outarray[i] = i+1;
                }
                break;
            case LibUtilities::eModified_A:
            case LibUtilities::eModified_B:
                for(i = 0 ; i < GetNcoeffs()-2;i++)
                {
                    outarray[i] = i+2;
                }
                break;
            default:
                ASSERTL0(0,"Mapping array is not defined for this expansion");
                break;
            }
        }

        int StdSegExp::v_GetVertexMap(int localVertexId)
        {
            ASSERTL0((localVertexId==0)||(localVertexId==1),"local vertex id"
                     "must be between 0 or 1");

            int localDOF = localVertexId;

            if( (m_base[0]->GetBasisType()==LibUtilities::eGLL_Lagrange) &&
                (localVertexId==1) )
            {
                localDOF = m_base[0]->GetNumModes()-1;
            }
            return localDOF;
        }


    }//end namespace
}//end namespace

