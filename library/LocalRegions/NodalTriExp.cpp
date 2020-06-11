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
#include <LibUtilities/Foundations/Interp.h>

using namespace std;

namespace Nektar
{
    namespace LocalRegions
    {
        NodalTriExp::NodalTriExp(const LibUtilities::BasisKey &Ba,
                                 const LibUtilities::BasisKey &Bb,
                                 const LibUtilities::PointsType Ntype,
                                 const SpatialDomains::TriGeomSharedPtr &geom):
            StdExpansion  (LibUtilities::StdTriData::getNumberOfCoefficients(Ba.GetNumModes(),(Bb.GetNumModes())),2,Ba,Bb),
            StdExpansion2D(LibUtilities::StdTriData::getNumberOfCoefficients(Ba.GetNumModes(),(Bb.GetNumModes())),Ba,Bb),
            StdNodalTriExp(Ba,Bb,Ntype),
            Expansion     (geom),
            Expansion2D   (geom),
            m_matrixManager(
                std::bind(&NodalTriExp::CreateMatrix, this, std::placeholders::_1),
                std::string("NodalTriExpMatrix")),
            m_staticCondMatrixManager(
                std::bind(&NodalTriExp::CreateStaticCondMatrix, this, std::placeholders::_1),
                std::string("NodalTriExpStaticCondMatrix"))
        {
        }

        NodalTriExp::NodalTriExp(const NodalTriExp &T):
            StdExpansion(T),
            StdExpansion2D(T),
            StdRegions::StdNodalTriExp(T),
            Expansion   (T),
            Expansion2D (T),
            m_matrixManager(T.m_matrixManager),
            m_staticCondMatrixManager(T.m_staticCondMatrixManager)
        {
        }

        NodalTriExp::~NodalTriExp()
        {
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


        NekDouble NodalTriExp::Integral(const Array<OneD, const NekDouble> &inarray)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            Array<OneD, const NekDouble> jac = m_metricinfo->GetJac(GetPointsKeys());
            NekDouble ival;
            Array<OneD,NekDouble> tmp(nquad0*nquad1);

            // multiply inarray with Jacobian
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0*nquad1, jac, 1, inarray, 1,tmp, 1);
            }
            else
            {
                Vmath::Smul(nquad0*nquad1, jac[0], inarray, 1, tmp, 1);
            }

            // call StdQuadExp version;
            ival = StdNodalTriExp::v_Integral(tmp);
            return ival;
        }


        void NodalTriExp::IProductWRTBase_SumFac(const Array<OneD, const NekDouble>& inarray,
                                                 Array<OneD, NekDouble> &outarray,
                                                 bool multiplybyweights)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    order1 = m_base[1]->GetNumModes();

            if(multiplybyweights)
            {
                Array<OneD,NekDouble> tmp(nquad0*nquad1+nquad0*order1);
                Array<OneD,NekDouble> wsp(tmp+nquad0*nquad1);

                MultiplyByQuadratureMetric(inarray,tmp);
                StdTriExp::IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),m_base[1]->GetBdata(),tmp,outarray,wsp);
                NodalToModalTranspose(outarray,outarray);
            }
            else
            {
                Array<OneD,NekDouble> wsp(nquad0*order1);

                StdTriExp::IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),m_base[1]->GetBdata(),inarray,outarray,wsp);
                NodalToModalTranspose(outarray,outarray);
            }
        }

        void NodalTriExp::IProductWRTBase_MatOp(const Array<OneD, const NekDouble>& inarray,
                                   Array<OneD, NekDouble> &outarray)
        {
            int nq = GetTotPoints();
            MatrixKey      iprodmatkey(StdRegions::eIProductWRTBase,DetShapeType(),*this);
            DNekScalMatSharedPtr iprodmat = m_matrixManager[iprodmatkey];

            Blas::Dgemv('N',m_ncoeffs,nq,iprodmat->Scale(),(iprodmat->GetOwnedMatrix())->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
        }

        void NodalTriExp::IProductWRTDerivBase_SumFac(const int dir,
                                                      const Array<OneD, const NekDouble>& inarray,
                                                      Array<OneD, NekDouble> & outarray)
        {
            ASSERTL1((dir==0)||(dir==1)||(dir==2),"Invalid direction.");
            ASSERTL1((dir==2)?(m_geom->GetCoordim()==3):true,"Invalid direction.");

            int    i;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nqtot  = nquad0*nquad1;
            int    wspsize = max(nqtot,m_ncoeffs);

            const Array<TwoD, const NekDouble>& df =
                                m_metricinfo->GetDerivFactors(GetPointsKeys());

            Array<OneD, NekDouble> tmp0 (6*wspsize);
            Array<OneD, NekDouble> tmp1 (tmp0 +   wspsize);
            Array<OneD, NekDouble> tmp2 (tmp0 + 2*wspsize);
            Array<OneD, NekDouble> tmp3 (tmp0 + 3*wspsize);
            Array<OneD, NekDouble> gfac0(tmp0 + 4*wspsize);
            Array<OneD, NekDouble> gfac1(tmp0 + 5*wspsize);

            const Array<OneD, const NekDouble>& z0 = m_base[0]->GetZ();
            const Array<OneD, const NekDouble>& z1 = m_base[1]->GetZ();

            // set up geometric factor: 2/(1-z1)
            for(i = 0; i < nquad1; ++i)
            {
                gfac0[i] = 2.0/(1-z1[i]);
            }
            for(i = 0; i < nquad0; ++i)
            {
                gfac1[i] = 0.5*(1+z0[i]);
            }

            for(i = 0; i < nquad1; ++i)
            {
                Vmath::Smul(nquad0,gfac0[i],&inarray[0]+i*nquad0,1,&tmp0[0]+i*nquad0,1);
            }

            for(i = 0; i < nquad1; ++i)
            {
                Vmath::Vmul(nquad0,&gfac1[0],1,&tmp0[0]+i*nquad0,1,&tmp1[0]+i*nquad0,1);
            }

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nqtot,&df[2*dir][0],  1,&tmp0[0],   1,&tmp0[0],1);
                Vmath::Vmul(nqtot,&df[2*dir+1][0],1,&tmp1[0],   1,&tmp1[0],1);
                Vmath::Vmul(nqtot,&df[2*dir+1][0],1,&inarray[0],1,&tmp2[0],1);
            }
            else
            {
                Vmath::Smul(nqtot, df[2*dir][0],   tmp0,    1, tmp0, 1);
                Vmath::Smul(nqtot, df[2*dir+1][0], tmp1,    1, tmp1, 1);
                Vmath::Smul(nqtot, df[2*dir+1][0], inarray, 1, tmp2, 1);
            }
            Vmath::Vadd(nqtot, tmp0, 1, tmp1, 1, tmp1, 1);

            MultiplyByQuadratureMetric(tmp1,tmp1);
            MultiplyByQuadratureMetric(tmp2,tmp2);

            IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),m_base[1]->GetBdata() ,tmp1,tmp3    ,tmp0);
            IProductWRTBase_SumFacKernel(m_base[0]->GetBdata() ,m_base[1]->GetDbdata(),tmp2,outarray,tmp0);
            Vmath::Vadd(m_ncoeffs, tmp3, 1, outarray, 1, outarray, 1);

            NodalToModalTranspose(outarray,outarray);
        }

        void NodalTriExp::IProductWRTDerivBase_MatOp(const int dir,
                                                     const Array<OneD, const NekDouble>& inarray,
                                                     Array<OneD, NekDouble> &outarray)
        {
            int nq = GetTotPoints();
            StdRegions::MatrixType mtype = StdRegions::eIProductWRTDerivBase0;

            switch(dir)
            {
            case 0:
                {
                    mtype = StdRegions::eIProductWRTDerivBase0;
                }
                break;
            case 1:
                {
                    mtype = StdRegions::eIProductWRTDerivBase1;
                }
                break;
            case 2:
                {
                    mtype = StdRegions::eIProductWRTDerivBase2;
                }
                break;
            default:
                {
                    ASSERTL1(false,"input dir is out of range");
                }
                break;
            }

            MatrixKey      iprodmatkey(mtype,DetShapeType(),*this);
            DNekScalMatSharedPtr iprodmat = m_matrixManager[iprodmatkey];

            Blas::Dgemv('N',m_ncoeffs,nq,iprodmat->Scale(),(iprodmat->GetOwnedMatrix())->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);

        }

        ///////////////////////////////
        /// Differentiation Methods
        ///////////////////////////////

        /**
            \brief Calculate the deritive of the physical points
        **/
        void NodalTriExp::PhysDeriv(const Array<OneD, const NekDouble> & inarray,
                                    Array<OneD,NekDouble> &out_d0,
                                    Array<OneD,NekDouble> &out_d1,
                                    Array<OneD,NekDouble> &out_d2)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int     nqtot = nquad0*nquad1;
            const Array<TwoD, const NekDouble>& df
                            = m_metricinfo->GetDerivFactors(GetPointsKeys());

            Array<OneD,NekDouble> diff0(2*nqtot);
            Array<OneD,NekDouble> diff1(diff0+nqtot);

            StdNodalTriExp::v_PhysDeriv(inarray, diff0, diff1);

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                if(out_d0.size())
                {
                    Vmath::Vmul  (nqtot,df[0],1,diff0,1, out_d0, 1);
                    Vmath::Vvtvp (nqtot,df[1],1,diff1,1, out_d0, 1, out_d0,1);
                }

                if(out_d1.size())
                {
                    Vmath::Vmul  (nqtot,df[2],1,diff0,1, out_d1, 1);
                    Vmath::Vvtvp (nqtot,df[3],1,diff1,1, out_d1, 1, out_d1,1);
                }

                if(out_d2.size())
                {
                    Vmath::Vmul  (nqtot,df[4],1,diff0,1, out_d2, 1);
                    Vmath::Vvtvp (nqtot,df[5],1,diff1,1, out_d2, 1, out_d2,1);
                }
            }
            else // regular geometry
            {
                if(out_d0.size())
                {
                    Vmath::Smul (nqtot, df[0][0], diff0, 1, out_d0, 1);
                    Blas::Daxpy (nqtot, df[1][0], diff1, 1, out_d0, 1);
                }

                if(out_d1.size())
                {
                    Vmath::Smul (nqtot, df[2][0], diff0, 1, out_d1, 1);
                    Blas::Daxpy (nqtot, df[3][0], diff1, 1, out_d1, 1);
                }

                if(out_d2.size())
                {
                    Vmath::Smul (nqtot, df[4][0], diff0, 1, out_d2, 1);
                    Blas::Daxpy (nqtot, df[5][0], diff1, 1, out_d2, 1);
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
        void NodalTriExp::FwdTrans(const Array<OneD, const NekDouble> & inarray,
                                   Array<OneD,NekDouble> &outarray)
        {
            IProductWRTBase(inarray,outarray);

            // get Mass matrix inverse
            MatrixKey  masskey(StdRegions::eInvMass, DetShapeType(),*this,StdRegions::NullConstFactorMap,StdRegions::NullVarCoeffMap,m_nodalPointsKey.GetPointsType());
            DNekScalMatSharedPtr  matsys = m_matrixManager[masskey];

            // copy inarray in case inarray == outarray
            NekVector<NekDouble> in(m_ncoeffs,outarray,eCopy);
            NekVector<NekDouble> out(m_ncoeffs,outarray,eWrapper);

            out = (*matsys)*in;
        }

        void NodalTriExp::GeneralMatrixOp_MatOp(const Array<OneD, const NekDouble> &inarray,
                                                Array<OneD,NekDouble> &outarray,
                                                const StdRegions::StdMatrixKey &mkey)
        {
            DNekScalMatSharedPtr   mat = GetLocMatrix(mkey);

            if(inarray.get() == outarray.get())
            {
                Array<OneD,NekDouble> tmp(m_ncoeffs);
                Vmath::Vcopy(m_ncoeffs,inarray.get(),1,tmp.get(),1);

                Blas::Dgemv('N',m_ncoeffs,m_ncoeffs,mat->Scale(),(mat->GetOwnedMatrix())->GetPtr().get(),
                            m_ncoeffs, tmp.get(), 1, 0.0, outarray.get(), 1);
            }
            else
            {
                Blas::Dgemv('N',m_ncoeffs,m_ncoeffs,mat->Scale(),(mat->GetOwnedMatrix())->GetPtr().get(),
                            m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
            }
        }

        void NodalTriExp::GetCoords(Array<OneD,NekDouble> &coords_0,
                                    Array<OneD,NekDouble> &coords_1,
                                    Array<OneD,NekDouble> &coords_2)
        {
            Expansion::v_GetCoords(coords_0, coords_1, coords_2);
        }

        // get the coordinates "coords" at the local coordinates "Lcoords"
        void NodalTriExp::GetCoord(const Array<OneD, const NekDouble> &Lcoords,
                                   Array<OneD,NekDouble> &coords)
        {
            int  i;

            ASSERTL1(Lcoords[0] >= -1.0 && Lcoords[1] <= 1.0 &&
                     Lcoords[1] >= -1.0 && Lcoords[1]  <=1.0,
                     "Local coordinates are not in region [-1,1]");

            m_geom->FillGeom();

            for(i = 0; i < m_geom->GetCoordim(); ++i)
            {
                coords[i] = m_geom->GetCoord(i,Lcoords);
            }
        }

        DNekMatSharedPtr NodalTriExp::CreateStdMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
            LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
            LibUtilities::PointsType ntype = m_nodalPointsKey.GetPointsType();
            StdRegions::StdNodalTriExpSharedPtr tmp = MemoryManager<StdNodalTriExp>::
                AllocateSharedPtr(bkey0,bkey1,ntype);

            return tmp->GetStdMatrix(mkey);
        }

        NekDouble NodalTriExp::PhysEvaluate(
            const Array<OneD, const NekDouble> &coord,
            const Array<OneD, const NekDouble> &physvals)

        {
            Array<OneD,NekDouble> Lcoord = Array<OneD,NekDouble>(2);

            ASSERTL0(m_geom,"m_geom not defined");
            m_geom->GetLocCoords(coord,Lcoord);

            return StdNodalTriExp::v_PhysEvaluate(Lcoord, physvals);
        }

        DNekScalMatSharedPtr NodalTriExp::CreateMatrix(const MatrixKey &mkey)
        {
            DNekScalMatSharedPtr returnval;
            LibUtilities::PointsKeyVector ptsKeys = GetPointsKeys();

            ASSERTL2(m_metricinfo->GetGtype() != SpatialDomains::eNoGeomType,"Geometric information is not set up");

            StdRegions::MatrixType mtype = mkey.GetMatrixType();

            switch(mtype)
            {
            case StdRegions::eMass:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble jac = (m_metricinfo->GetJac(ptsKeys))[0];
                        DNekMatSharedPtr mat = GetStdMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(jac,mat);
                    }
                }
                break;
            case StdRegions::eInvMass:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        StdRegions::StdMatrixKey masskey(StdRegions::eMass,DetShapeType(),
                                                         *this);
                        DNekMatSharedPtr mat = GenMatrix(masskey);
                        mat->Invert();

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble fac = 1.0/(m_metricinfo->GetJac(ptsKeys))[0];
                        DNekMatSharedPtr mat = GetStdMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(fac,mat);
                    }
                }
                break;
            case StdRegions::eLaplacian:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(mkey);

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        ASSERTL1(m_geom->GetCoordim() == 2,"Standard Region Laplacian is only set up for Quads in two-dimensional");
                        MatrixKey lap00key(StdRegions::eLaplacian00,
                                           mkey.GetShapeType(), *this);
                        MatrixKey lap01key(StdRegions::eLaplacian01,
                                           mkey.GetShapeType(), *this);
                        MatrixKey lap11key(StdRegions::eLaplacian11,
                                           mkey.GetShapeType(), *this);

                        DNekMat &lap00 = *GetStdMatrix(lap00key);
                        DNekMat &lap01 = *GetStdMatrix(lap01key);
                        DNekMat &lap11 = *GetStdMatrix(lap11key);

                        NekDouble jac = (m_metricinfo->GetJac(ptsKeys))[0];
                        Array<TwoD, const NekDouble> gmat =
                                                m_metricinfo->GetGmat(ptsKeys);

                        int rows = lap00.GetRows();
                        int cols = lap00.GetColumns();

                        DNekMatSharedPtr lap = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);

                        (*lap) = gmat[0][0] * lap00 +
                                 gmat[1][0] * (lap01 + Transpose(lap01)) +
                                 gmat[3][0] * lap11;

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(jac,lap);
                    }
                }
                break;
            case StdRegions::eHelmholtz:
                {
                    NekDouble factor = mkey.GetConstFactor(StdRegions::eFactorLambda);
                    MatrixKey masskey(StdRegions::eMass,
                                      mkey.GetShapeType(), *this);
                    DNekScalMat &MassMat = *(this->m_matrixManager[masskey]);
                    MatrixKey lapkey(StdRegions::eLaplacian,
                                     mkey.GetShapeType(), *this);
                    DNekScalMat &LapMat = *(this->m_matrixManager[lapkey]);

                    int rows = LapMat.GetRows();
                    int cols = LapMat.GetColumns();

                    DNekMatSharedPtr helm = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);

                    NekDouble one = 1.0;
                    (*helm) = LapMat + factor*MassMat;

                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,helm);
                }
                break;
            default:
                NEKERROR(ErrorUtil::efatal, "Matrix creation not defined");
                break;
            }

            return returnval;
        }

        DNekScalBlkMatSharedPtr NodalTriExp::CreateStaticCondMatrix(const MatrixKey &mkey)
        {
            DNekScalBlkMatSharedPtr returnval;

            ASSERTL2(m_metricinfo->GetGtype() != SpatialDomains::eNoGeomType,"Geometric information is not set up");

            // set up block matrix system
            unsigned int nbdry = NumBndryCoeffs();
            unsigned int nint = (unsigned int)(m_ncoeffs - nbdry);
            unsigned int exp_size[] = {nbdry,nint};
            unsigned int nblks = 2;
            returnval = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nblks,nblks,exp_size,exp_size); //Really need a constructor which takes Arrays
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

                    returnval->SetBlock(0,0,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,Asubmat = mat->GetBlock(0,0)));
                    returnval->SetBlock(0,1,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Asubmat = mat->GetBlock(0,1)));
                    returnval->SetBlock(1,0,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,Asubmat = mat->GetBlock(1,0)));
                    returnval->SetBlock(1,1,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(invfactor,Asubmat = mat->GetBlock(1,1)));
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


        StdRegions::StdExpansionSharedPtr NodalTriExp::v_GetStdExp(void) const
        {

            return MemoryManager<StdRegions::StdNodalTriExp>
                    ::AllocateSharedPtr(m_base[0]->GetBasisKey(),
                                        m_base[1]->GetBasisKey(),
                                        m_nodalPointsKey.GetPointsType());
        }

        StdRegions::StdExpansionSharedPtr NodalTriExp::v_GetLinStdExp(void) const
        {
            LibUtilities::BasisKey bkey0(m_base[0]->GetBasisType(),
                           2, m_base[0]->GetPointsKey());
            LibUtilities::BasisKey bkey1(m_base[1]->GetBasisType(),
                           2, m_base[1]->GetPointsKey());

            return MemoryManager<StdRegions::StdNodalTriExp>
                ::AllocateSharedPtr( bkey0, bkey1, m_nodalPointsKey.GetPointsType());
        }

        StdRegions::Orientation NodalTriExp::v_GetEorient(int edge)
        {
            return GetGeom2D()->GetEorient(edge);
        }

        DNekMatSharedPtr NodalTriExp::v_GenMatrix(const StdRegions::StdMatrixKey &mkey)
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
                returnval = Expansion2D::v_GenMatrix(mkey);
                break;
            default:
                returnval = StdNodalTriExp::v_GenMatrix(mkey);
                break;
            }
            return returnval;
        }

        void NodalTriExp::v_ComputeEdgeNormal(const int edge)
        {
            int i;
            const SpatialDomains::GeomFactorsSharedPtr & geomFactors = GetGeom()->GetMetricInfo();
            const SpatialDomains::GeomType type = geomFactors->GetGtype();

            LibUtilities::PointsKeyVector ptsKeys = GetPointsKeys();
            const Array<TwoD, const NekDouble> & df = geomFactors->GetDerivFactors(ptsKeys);
            const Array<OneD, const NekDouble> & jac  = geomFactors->GetJac(ptsKeys);
            int nqe = m_base[0]->GetNumPoints();
            int dim = GetCoordim();

            m_edgeNormals[edge] = Array<OneD, Array<OneD, NekDouble> >(dim);
            Array<OneD, Array<OneD, NekDouble> > &normal = m_edgeNormals[edge];
            for (i = 0; i < dim; ++i)
            {
                normal[i] = Array<OneD, NekDouble>(nqe);
            }

            size_t nqb = nqe;
            size_t nbnd= edge;
            m_elmtBndNormDirElmtLen[nbnd] = Array<OneD, NekDouble> {nqb, 0.0};
            Array<OneD, NekDouble> &length = m_elmtBndNormDirElmtLen[nbnd];

            // Regular geometry case
            if((type == SpatialDomains::eRegular)||(type == SpatialDomains::eMovingRegular))
            {
                NekDouble fac;
                // Set up normals
                switch(edge)
                {
                case 0:
                    for(i = 0; i < GetCoordim(); ++i)
                    {
                        Vmath::Fill(nqe,-df[2*i+1][0],normal[i],1);
                    }
                    break;
                case 1:
                    for(i = 0; i < GetCoordim(); ++i)
                    {
                        Vmath::Fill(nqe,df[2*i+1][0] + df[2*i][0],normal[i],1);
                    }
                        break;
                case 2:
                    for(i = 0; i < GetCoordim(); ++i)
                    {
                        Vmath::Fill(nqe,-df[2*i][0],normal[i],1);
                    }
                    break;
                default:
                    ASSERTL0(false,"Edge is out of range (edge < 3)");
                }

                // normalise
                fac = 0.0;
                for(i =0 ; i < GetCoordim(); ++i)
                {
                    fac += normal[i][0]*normal[i][0];
                }
                fac = 1.0/sqrt(fac);

                Vmath::Fill(nqb, fac, length, 1);

                for (i = 0; i < GetCoordim(); ++i)
                {
                    Vmath::Smul(nqe,fac,normal[i],1,normal[i],1);
                }
            }
            else   // Set up deformed normals
            {
                int j;

                int nquad0 = ptsKeys[0].GetNumPoints();
                int nquad1 = ptsKeys[1].GetNumPoints();

                LibUtilities::PointsKey from_key;

                Array<OneD,NekDouble> normals(GetCoordim()*max(nquad0,nquad1),0.0);
                Array<OneD,NekDouble> edgejac(GetCoordim()*max(nquad0,nquad1),0.0);

                // Extract Jacobian along edges and recover local
                // derivates (dx/dr) for polynomial interpolation by
                // multiplying m_gmat by jacobian
                switch(edge)
                {
                case 0:
                    for(j = 0; j < nquad0; ++j)
                    {
                        edgejac[j] = jac[j];
                        for(i = 0; i < GetCoordim(); ++i)
                        {
                            normals[i*nquad0+j] = -df[2*i+1][j]*edgejac[j];
                        }
                    }
                    from_key = ptsKeys[0];
                    break;
                case 1:
                    for(j = 0; j < nquad1; ++j)
                    {
                        edgejac[j] = jac[nquad0*j+nquad0-1];
                        for(i = 0; i < GetCoordim(); ++i)
                        {
                            normals[i*nquad1+j] = (df[2*i][nquad0*j + nquad0-1] +  df[2*i+1][nquad0*j + nquad0-1])*edgejac[j];
                        }
                    }
                    from_key = ptsKeys[1];
                    break;
                case 2:
                    for(j = 0; j < nquad1; ++j)
                    {
                        edgejac[j] = jac[nquad0*j];
                        for(i = 0; i < GetCoordim(); ++i)
                        {
                            normals[i*nquad1+j] = -df[2*i][nquad0*j]*edgejac[j];
                        }
                    }
                    from_key = ptsKeys[1];
                    break;
                default:
                    ASSERTL0(false,"edge is out of range (edge < 3)");

                }

                int nq  = from_key.GetNumPoints();
                Array<OneD,NekDouble> work(nqe,0.0);

                // interpolate Jacobian and invert
                LibUtilities::Interp1D(from_key,jac,m_base[0]->GetPointsKey(),work);
                Vmath::Sdiv(nq,1.0,&work[0],1,&work[0],1);

                // interpolate
                for(i = 0; i < GetCoordim(); ++i)
                {
                    LibUtilities::Interp1D(from_key,&normals[i*nq],m_base[0]->GetPointsKey(),&normal[i][0]);
                    Vmath::Vmul(nqe,work,1,normal[i],1,normal[i],1);
                }

                //normalise normal vectors
                Vmath::Zero(nqe,work,1);
                for(i = 0; i < GetCoordim(); ++i)
                {
                    Vmath::Vvtvp(nqe,normal[i],1, normal[i],1,work,1,work,1);
                }

                Vmath::Vsqrt(nqe,work,1,work,1);
                Vmath::Sdiv(nqe,1.0,work,1,work,1);

                Vmath::Vcopy(nqb, work, 1, length, 1);

                for(i = 0; i < GetCoordim(); ++i)
                {
                    Vmath::Vmul(nqe,normal[i],1,work,1,normal[i],1);
                }

                // Reverse direction so that points are in
                // anticlockwise direction if edge >=2
                if(edge >= 2)
                {
                    for(i = 0; i < GetCoordim(); ++i)
                    {
                        Vmath::Reverse(nqe,normal[i],1, normal[i],1);
                    }
                }
            }
        }

    }//end of namespace
}//end of namespace
