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
// Description: Expansion for quadrilateral elements.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/Foundations/InterpCoeff.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/Expansion3D.h>
#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/Foundations/Interp.h>
#include <LocalRegions/SegExp.h>

using namespace std;

namespace Nektar
{
    namespace LocalRegions
    {
        QuadExp::QuadExp(const LibUtilities::BasisKey &Ba,
                         const LibUtilities::BasisKey &Bb,
                         const SpatialDomains::QuadGeomSharedPtr &geom):
             StdExpansion  (Ba.GetNumModes()*Bb.GetNumModes(),2,Ba,Bb),
             StdExpansion2D(Ba.GetNumModes()*Bb.GetNumModes(),Ba,Bb),
             StdQuadExp    (Ba,Bb),
             Expansion     (geom),
             Expansion2D   (geom),
             m_matrixManager(
                 std::bind(&QuadExp::CreateMatrix, this, std::placeholders::_1),
                 std::string("QuadExpMatrix")),
             m_staticCondMatrixManager(
                 std::bind(&QuadExp::CreateStaticCondMatrix, this, std::placeholders::_1),
                 std::string("QuadExpStaticCondMatrix"))
        {
        }


        QuadExp::QuadExp(const QuadExp &T):
            StdExpansion(T),
            StdExpansion2D(T),
            StdQuadExp(T),
            Expansion   (T),
            Expansion2D (T),
            m_matrixManager(T.m_matrixManager),
            m_staticCondMatrixManager(T.m_staticCondMatrixManager)
        {
        }


        QuadExp::~QuadExp()
        {
        }


        NekDouble QuadExp::v_Integral(
            const Array<OneD,
            const NekDouble> &inarray)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            Array<OneD, const NekDouble> jac = m_metricinfo->GetJac(GetPointsKeys());
            NekDouble ival;
            Array<OneD,NekDouble> tmp(nquad0*nquad1);

            // multiply inarray with Jacobian
            if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0*nquad1, jac, 1, inarray, 1, tmp, 1);
            }
            else
            {
                Vmath::Smul(nquad0*nquad1, jac[0], inarray, 1, tmp, 1);
            }

            // call StdQuadExp version;
            ival = StdQuadExp::v_Integral(tmp);
            return  ival;
        }

        void QuadExp::v_PhysDeriv(
            const Array<OneD, const NekDouble> & inarray,
                  Array<OneD,NekDouble> &out_d0,
                  Array<OneD,NekDouble> &out_d1,
                  Array<OneD,NekDouble> &out_d2)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int     nqtot = nquad0*nquad1;
            const Array<TwoD, const NekDouble>& df = m_metricinfo->GetDerivFactors(GetPointsKeys());
            Array<OneD,NekDouble> diff0(2*nqtot);
            Array<OneD,NekDouble> diff1(diff0+nqtot);

            StdQuadExp::v_PhysDeriv(inarray, diff0, diff1);

            if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                if (out_d0.size())
                {
                    Vmath::Vmul  (nqtot, df[0], 1, diff0, 1, out_d0, 1);
                    Vmath::Vvtvp (nqtot, df[1], 1, diff1, 1, out_d0, 1,
                    					 out_d0,1);
                }

                if(out_d1.size())
                {
                    Vmath::Vmul  (nqtot,df[2],1,diff0,1, out_d1, 1);
                    Vmath::Vvtvp (nqtot,df[3],1,diff1,1, out_d1, 1, out_d1,1);
                }

                if (out_d2.size())
                {
                    Vmath::Vmul  (nqtot,df[4],1,diff0,1, out_d2, 1);
                    Vmath::Vvtvp (nqtot,df[5],1,diff1,1, out_d2, 1, out_d2,1);
                }
            }
            else // regular geometry
            {
                if (out_d0.size())
                {
                    Vmath::Smul (nqtot, df[0][0], diff0, 1, out_d0, 1);
                    Blas::Daxpy (nqtot, df[1][0], diff1, 1, out_d0, 1);
                }

                if (out_d1.size())
                {
                    Vmath::Smul (nqtot, df[2][0], diff0, 1, out_d1, 1);
                    Blas::Daxpy (nqtot, df[3][0], diff1, 1, out_d1, 1);
                }

                if (out_d2.size())
                {
                    Vmath::Smul (nqtot, df[4][0], diff0, 1, out_d2, 1);
                    Blas::Daxpy (nqtot, df[5][0], diff1, 1, out_d2, 1);
                }
            }
        }


        void QuadExp::v_PhysDeriv(
            const int dir,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD, NekDouble> &outarray)
        {
            switch (dir)
            {
            case 0:
                {
                    v_PhysDeriv(inarray, outarray, NullNekDouble1DArray,
                                NullNekDouble1DArray);
                }
                break;
            case 1:
                {
                    v_PhysDeriv(inarray, NullNekDouble1DArray, outarray,
                                NullNekDouble1DArray);
                }
                break;
            case 2:
                {
                    v_PhysDeriv(inarray, NullNekDouble1DArray,
                                NullNekDouble1DArray, outarray);
                }
                break;
            default:
                {
                    ASSERTL1(false,"input dir is out of range");
                }
                break;
            }
        }


        void QuadExp::v_PhysDirectionalDeriv(
            const Array<OneD, const NekDouble> & inarray,
            const Array<OneD, const NekDouble>& direction,
                  Array<OneD,NekDouble> &out)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nqtot = nquad0*nquad1;

            const Array<TwoD, const NekDouble>& df = m_metricinfo->GetDerivFactors(GetPointsKeys());

            Array<OneD,NekDouble> diff0(2*nqtot);
            Array<OneD,NekDouble> diff1(diff0+nqtot);

            StdQuadExp::v_PhysDeriv(inarray, diff0, diff1);

            if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Array<OneD, Array<OneD, NekDouble> > tangmat(2);

                // d/dx_v^s = v_x*ds/dx + v_y*ds/dy + v_z*dx/dz
                for (int i=0; i< 2; ++i)
                {
                    tangmat[i] = Array<OneD, NekDouble>(nqtot,0.0);
                    for (int k=0; k<(m_geom->GetCoordim()); ++k)
                    {
                        Vmath::Vvtvp(nqtot,
                                     &df[2*k+i][0], 1,
                                     &direction[k*nqtot], 1,
                                     &tangmat[i][0], 1,
                                     &tangmat[i][0], 1);
                    }
                }

                /// D_v = d/dx_v^s + d/dx_v^r
                if (out.size())
                {
                    Vmath::Vmul  (nqtot,
                                  &tangmat[0][0], 1,
                                  &diff0[0], 1,
                                  &out[0], 1);
                    Vmath::Vvtvp (nqtot,
                                  &tangmat[1][0], 1,
                                  &diff1[0], 1,
                                  &out[0], 1,
                                  &out[0], 1);
                }

            }
            else
            {
                ASSERTL1(m_metricinfo->GetGtype() ==
                         SpatialDomains::eDeformed,"Wrong route");
            }
        }


        void QuadExp::v_FwdTrans(
            const Array<OneD, const NekDouble> & inarray,
            Array<OneD,NekDouble> &outarray)
        {
            if ((m_base[0]->Collocation())&&(m_base[1]->Collocation()))
            {
                Vmath::Vcopy(m_ncoeffs, inarray, 1, outarray, 1);
            }
            else
            {
                IProductWRTBase(inarray,outarray);

                // get Mass matrix inverse
                MatrixKey             masskey(StdRegions::eInvMass,
                                              DetShapeType(),*this);
                DNekScalMatSharedPtr  matsys = m_matrixManager[masskey];

                // copy inarray in case inarray == outarray
                NekVector<NekDouble> in(m_ncoeffs,outarray,eCopy);
                NekVector<NekDouble> out(m_ncoeffs,outarray,eWrapper);

                out = (*matsys)*in;
            }
        }


        void QuadExp::v_FwdTrans_BndConstrained(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD, NekDouble> &outarray)
        {
            if ((m_base[0]->Collocation())&&(m_base[1]->Collocation()))
            {
                Vmath::Vcopy(m_ncoeffs, inarray, 1, outarray, 1);
            }
            else
            {
                int i,j;
                int npoints[2] = {m_base[0]->GetNumPoints(),
                                  m_base[1]->GetNumPoints()};
                int nmodes[2]  = {m_base[0]->GetNumModes(),
                                  m_base[1]->GetNumModes()};

                fill(outarray.get(), outarray.get()+m_ncoeffs, 0.0 );

                if(nmodes[0] == 1 && nmodes[1] == 1)
                {
                    outarray[0] = inarray[0];
                    return;
                }

                Array<OneD, NekDouble> physEdge[4];
                Array<OneD, NekDouble> coeffEdge[4];
                StdRegions::Orientation orient[4];
                for (i = 0; i < 4; i++)
                {
                    physEdge[i]  = Array<OneD, NekDouble>(npoints[i%2]);
                    coeffEdge[i] = Array<OneD, NekDouble>(nmodes[i%2]);
                    orient[i]    = GetEorient(i);
                }

                for (i = 0; i < npoints[0]; i++)
                {
                    physEdge[0][i] = inarray[i];
                    physEdge[2][i] = inarray[npoints[0]*(npoints[1]-1)+i];
                }

                for (i = 0; i < npoints[1]; i++)
                {
                    physEdge[1][i] =
                        inarray[npoints[0]-1+i*npoints[0]];
                    physEdge[3][i] =
                        inarray[i*npoints[0]];
                }

                for (i = 0; i < 4; i++)
                {
                    if ( orient[i] == StdRegions::eBackwards )
                    {
                        reverse((physEdge[i]).get(),
                                (physEdge[i]).get() + npoints[i%2] );
                    }
                }

                SegExpSharedPtr segexp[4];
                for (i = 0; i < 4; i++)
                {
                    segexp[i] = MemoryManager<LocalRegions::SegExp>::
                        AllocateSharedPtr(
                            m_base[i%2]->GetBasisKey(),GetGeom2D()->GetEdge(i));
                }

                Array<OneD, unsigned int> mapArray;
                Array<OneD, int>          signArray;
                NekDouble sign;

                for (i = 0; i < 4; i++)
                {
                    segexp[i%2]->FwdTrans_BndConstrained(
                        physEdge[i],coeffEdge[i]);

                    GetEdgeToElementMap(i,orient[i],mapArray,signArray);
                    for (j=0; j < nmodes[i%2]; j++)
                    {
                        sign = (NekDouble) signArray[j];
                        outarray[ mapArray[j] ] = sign * coeffEdge[i][j];
                    }
                }

                int nBoundaryDofs = NumBndryCoeffs();
                int nInteriorDofs = m_ncoeffs - nBoundaryDofs;

                if (nInteriorDofs > 0) {
                    Array<OneD, NekDouble> tmp0(m_ncoeffs);
                    Array<OneD, NekDouble> tmp1(m_ncoeffs);

                    StdRegions::StdMatrixKey
                        stdmasskey(StdRegions::eMass,DetShapeType(),*this);
                    MassMatrixOp(outarray,tmp0,stdmasskey);
                    IProductWRTBase(inarray,tmp1);

                    Vmath::Vsub(m_ncoeffs, tmp1, 1, tmp0, 1, tmp1, 1);

                    // get Mass matrix inverse (only of interior DOF)
                    // use block (1,1) of the static condensed system
                    // note: this block alreay contains the inverse matrix
                    MatrixKey
                        masskey(StdRegions::eMass,DetShapeType(),*this);
                    DNekScalMatSharedPtr matsys =
                        (m_staticCondMatrixManager[masskey])->GetBlock(1,1);

                    Array<OneD, NekDouble> rhs(nInteriorDofs);
                    Array<OneD, NekDouble> result(nInteriorDofs);

                    GetInteriorMap(mapArray);

                    for (i = 0; i < nInteriorDofs; i++)
                    {
                        rhs[i] = tmp1[ mapArray[i] ];
                    }

                    Blas::Dgemv('N', nInteriorDofs, nInteriorDofs,
                                matsys->Scale(),
                                &((matsys->GetOwnedMatrix())->GetPtr())[0],
                                nInteriorDofs,rhs.get(),1,0.0,result.get(),1);

                    for (i = 0; i < nInteriorDofs; i++)
                    {
                        outarray[ mapArray[i] ] = result[i];
                    }
                }
            }

        }


        void QuadExp::v_IProductWRTBase(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD, NekDouble> &outarray)
        {
            if (m_base[0]->Collocation() && m_base[1]->Collocation())
            {
                MultiplyByQuadratureMetric(inarray,outarray);
            }
            else
            {
                IProductWRTBase_SumFac(inarray,outarray);
            }
        }


        void QuadExp::v_IProductWRTDerivBase(
            const int dir,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD, NekDouble> & outarray)
        {
            IProductWRTDerivBase_SumFac(dir,inarray,outarray);
        }

        
        void QuadExp::v_IProductWRTBase_SumFac(
            const Array<OneD, const NekDouble>& inarray,
            Array<OneD, NekDouble> &outarray,
                                               bool multiplybyweights)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    order0 = m_base[0]->GetNumModes();

            if(multiplybyweights)
            {
                Array<OneD,NekDouble> tmp(nquad0*nquad1+nquad1*order0);
                Array<OneD,NekDouble> wsp(tmp+nquad0*nquad1);

                MultiplyByQuadratureMetric(inarray,tmp);
                StdQuadExp::IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                                         m_base[1]->GetBdata(),
                                                         tmp,outarray,wsp,true,true);
            }
            else
            {
                Array<OneD,NekDouble> wsp(nquad1*order0);

                StdQuadExp::IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                                         m_base[1]->GetBdata(),
                                                         inarray,outarray,wsp,true,true);
            }
        }


        void QuadExp::v_IProductWRTBase_MatOp(
             const Array<OneD, const NekDouble>& inarray,
                   Array<OneD, NekDouble> &outarray)
        {
            int nq = GetTotPoints();
            MatrixKey
                iprodmatkey(StdRegions::eIProductWRTBase,DetShapeType(),*this);
            DNekScalMatSharedPtr iprodmat = m_matrixManager[iprodmatkey];

            Blas::Dgemv('N',m_ncoeffs,nq,iprodmat->Scale(),
                        (iprodmat->GetOwnedMatrix())->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);

        }


        void QuadExp::v_IProductWRTDerivBase_SumFac(
            const int dir,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD, NekDouble> & outarray)
        {
            ASSERTL1((dir==0) || (dir==1) || (dir==2),
                     "Invalid direction.");
            ASSERTL1((dir==2) ? (m_geom->GetCoordim() ==3):true,
                     "Invalid direction.");

            int    nquad0  = m_base[0]->GetNumPoints();
            int    nquad1  = m_base[1]->GetNumPoints();
            int    nqtot   = nquad0*nquad1;
            int    nmodes0 = m_base[0]->GetNumModes();

            const Array<TwoD, const NekDouble>& df = m_metricinfo->GetDerivFactors(GetPointsKeys());

            Array<OneD, NekDouble> tmp1(2*nqtot+m_ncoeffs+nmodes0*nquad1);
            Array<OneD, NekDouble> tmp2(tmp1 +   nqtot);
            Array<OneD, NekDouble> tmp3(tmp1 + 2*nqtot);
            Array<OneD, NekDouble> tmp4(tmp1 + 2*nqtot+m_ncoeffs);

            if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nqtot,
                            &df[2*dir][0], 1,
                            inarray.get(), 1,
                            tmp1.get(), 1);
                Vmath::Vmul(nqtot,
                            &df[2*dir+1][0], 1,
                            inarray.get(), 1,
                            tmp2.get(),1);
            }
            else
            {
                Vmath::Smul(nqtot,
                            df[2*dir][0], inarray.get(), 1,
                            tmp1.get(), 1);
                Vmath::Smul(nqtot,
                            df[2*dir+1][0], inarray.get(), 1,
                            tmp2.get(), 1);
            }

            MultiplyByQuadratureMetric(tmp1,tmp1);
            MultiplyByQuadratureMetric(tmp2,tmp2);

            IProductWRTBase_SumFacKernel(
                m_base[0]->GetDbdata(), m_base[1]->GetBdata(),
                tmp1, tmp3, tmp4, false, true);
            IProductWRTBase_SumFacKernel(
                m_base[0]->GetBdata() , m_base[1]->GetDbdata(),
                tmp2, outarray, tmp4, true, false);
            Vmath::Vadd(m_ncoeffs, tmp3, 1, outarray, 1, outarray, 1);
        }


        void QuadExp::v_IProductWRTDerivBase_MatOp(
            const int dir,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD, NekDouble> &outarray)
        {
            int nq = GetTotPoints();
            StdRegions::MatrixType mtype = StdRegions::eIProductWRTDerivBase0;

            switch (dir)
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

            Blas::Dgemv('N', m_ncoeffs, nq, iprodmat->Scale(),
                        (iprodmat->GetOwnedMatrix())->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
        }


        void QuadExp::v_NormVectorIProductWRTBase(
            const Array<OneD, const NekDouble> &Fx,
            const Array<OneD, const NekDouble> &Fy,
            const Array<OneD, const NekDouble> &Fz,
                  Array<OneD,       NekDouble> &outarray)
        {
            int nq = m_base[0]->GetNumPoints()*m_base[1]->GetNumPoints();
            Array<OneD, NekDouble> Fn(nq);

            const Array<OneD, const Array<OneD, NekDouble> > &normals =
                GetLeftAdjacentElementExp()->GetFaceNormal(
                    GetLeftAdjacentElementFace());

            if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vvtvvtp(nq,&normals[0][0],1,&Fx[0],1,
                                  &normals[1][0],1,&Fy[0],1,&Fn[0],1);
                Vmath::Vvtvp  (nq,&normals[2][0],1,&Fz[0],1,&Fn[0],1,&Fn[0],1);
            }
            else
            {
                Vmath::Svtsvtp(nq,normals[0][0],&Fx[0],1,
                                  normals[1][0],&Fy[0],1,&Fn[0],1);
                Vmath::Svtvp  (nq,normals[2][0],&Fz[0],1,&Fn[0],1,&Fn[0],1);
            }

            IProductWRTBase(Fn,outarray);
        }

        void QuadExp::v_NormVectorIProductWRTBase(
            const Array<OneD, const Array<OneD, NekDouble> > &Fvec,
                  Array<OneD,       NekDouble>               &outarray)
        {
            NormVectorIProductWRTBase(Fvec[0], Fvec[1], Fvec[2], outarray);
        }

        StdRegions::StdExpansionSharedPtr QuadExp::v_GetStdExp(void) const
        {
            return MemoryManager<StdRegions::StdQuadExp>
                ::AllocateSharedPtr(m_base[0]->GetBasisKey(),
                                    m_base[1]->GetBasisKey());
        }


        StdRegions::StdExpansionSharedPtr QuadExp::v_GetLinStdExp(void) const
        {
            LibUtilities::BasisKey bkey0(m_base[0]->GetBasisType(),
                           2, m_base[0]->GetPointsKey());
            LibUtilities::BasisKey bkey1(m_base[1]->GetBasisType(),
                           2, m_base[1]->GetPointsKey());

            return MemoryManager<StdRegions::StdQuadExp>
                ::AllocateSharedPtr( bkey0, bkey1);
        }

        void QuadExp::v_GetCoords(
            Array<OneD, NekDouble> &coords_0,
            Array<OneD, NekDouble> &coords_1,
            Array<OneD, NekDouble> &coords_2)
        {
            Expansion::v_GetCoords(coords_0, coords_1, coords_2);
        }


        void QuadExp::v_GetCoord(const Array<OneD, const NekDouble> &Lcoords,
                               Array<OneD,NekDouble> &coords)
        {
            int  i;

            ASSERTL1(Lcoords[0] >= -1.0 && Lcoords[1] <= 1.0 &&
                     Lcoords[1] >= -1.0 && Lcoords[1]  <=1.0,
                     "Local coordinates are not in region [-1,1]");

            m_geom->FillGeom();
            for (i = 0; i < m_geom->GetCoordim(); ++i)
            {
                coords[i] = m_geom->GetCoord(i,Lcoords);
            }
        }



        /**
         * Given the local cartesian coordinate \a Lcoord evaluate the
         * value of physvals at this point by calling through to the
         * StdExpansion method
         */
        NekDouble QuadExp::v_StdPhysEvaluate(
            const Array<OneD, const NekDouble> &Lcoord,
            const Array<OneD, const NekDouble> &physvals)
        {
            // Evaluate point in local (eta) coordinates.
            return StdQuadExp::v_PhysEvaluate(Lcoord,physvals);
        }

        NekDouble QuadExp::v_PhysEvaluate(
            const Array<OneD, const NekDouble> &coord,
            const Array<OneD, const NekDouble> &physvals)
        {
            Array<OneD,NekDouble> Lcoord = Array<OneD, NekDouble>(2);

            ASSERTL0(m_geom,"m_geom not defined");
            m_geom->GetLocCoords(coord,Lcoord);

            return StdQuadExp::v_PhysEvaluate(Lcoord, physvals);
        }


        // Get edge values from the 2D Phys space along an edge
        // following a counter clockwise edge convention for definition
        // of edgedir, Note that point distribution is given by QuadExp.
        void QuadExp::v_GetEdgePhysVals(
            const int edge,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,NekDouble> &outarray)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();

            StdRegions::Orientation edgedir = GetEorient(edge);
            switch(edge)
            {
            case 0:
                if (edgedir == StdRegions::eForwards)
                {
                    Vmath::Vcopy(nquad0,&(inarray[0]),1,&(outarray[0]),1);
                }
                else
                {
                    Vmath::Vcopy(nquad0,&(inarray[0])+(nquad0-1),-1,
                                 &(outarray[0]),1);
                }
                break;
            case 1:
                if (edgedir == StdRegions::eForwards)
                {
                    Vmath::Vcopy(nquad1,&(inarray[0])+(nquad0-1),nquad0,
                                 &(outarray[0]),1);
                }
                else
                {
                    Vmath::Vcopy(nquad1,&(inarray[0])+(nquad0*nquad1-1),
                                 -nquad0, &(outarray[0]),1);
                }
                break;
            case 2:
                if (edgedir == StdRegions::eForwards)
                {
                    Vmath::Vcopy(nquad0,&(inarray[0])+(nquad0*nquad1-1),-1,
                                 &(outarray[0]),1);
                }
                else
                {
                    Vmath::Vcopy(nquad0,&(inarray[0])+nquad0*(nquad1-1),1,
                                 &(outarray[0]),1);
                }
                break;
            case 3:
                if (edgedir == StdRegions::eForwards)
                {
                    Vmath::Vcopy(nquad1,&(inarray[0]) + nquad0*(nquad1-1),
                                 -nquad0,&(outarray[0]),1);
                }
                else
                {
                    Vmath::Vcopy(nquad1,&(inarray[0]),nquad0,
                                 &(outarray[0]),1);
                }
                break;
            default:
                ASSERTL0(false,"edge value (< 3) is out of range");
                break;
            }
        }


        void QuadExp::v_GetTracePhysVals(
             const int edge,
             const StdRegions::StdExpansionSharedPtr &EdgeExp,
             const Array<OneD, const NekDouble> &inarray,
             Array<OneD,NekDouble> &outarray,
             StdRegions::Orientation  orient)
        {
            boost::ignore_unused(orient);
            v_GetEdgePhysVals(edge,EdgeExp,inarray,outarray);
        }


        void QuadExp::v_GetEdgePhysVals(
            const int edge,
            const StdRegions::StdExpansionSharedPtr &EdgeExp,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,NekDouble> &outarray)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();

            // Implementation for all the basis except Gauss points
            if (m_base[0]->GetPointsType() !=
                LibUtilities::eGaussGaussLegendre &&
                m_base[1]->GetPointsType() !=
                LibUtilities::eGaussGaussLegendre)
            {
                switch (edge)
                {
                    case 0:
                        Vmath::Vcopy(nquad0,&(inarray[0]),1,&(outarray[0]),1);
                        break;
                    case 1:
                        Vmath::Vcopy(nquad1,&(inarray[0])+(nquad0-1),
                                     nquad0,&(outarray[0]),1);
                        break;
                    case 2:
                        Vmath::Vcopy(nquad0,&(inarray[0])+nquad0*(nquad1-1),1,
                                     &(outarray[0]),1);
                        break;
                    case 3:
                        Vmath::Vcopy(nquad1,&(inarray[0]),nquad0,
                                     &(outarray[0]),1);
                        break;
                    default:
                        ASSERTL0(false,"edge value (< 3) is out of range");
                        break;
                }
            }
            else
            {
                QuadExp::v_GetEdgeInterpVals(edge, inarray, outarray);
            }

            // Interpolate if required
            if (m_base[edge%2]->GetPointsKey() !=
                EdgeExp->GetBasis(0)->GetPointsKey())
            {
                Array<OneD,NekDouble> outtmp(max(nquad0,nquad1));

                outtmp = outarray;

                LibUtilities::Interp1D(
                    m_base[edge%2]->GetPointsKey(), outtmp,
                    EdgeExp->GetBasis(0)->GetPointsKey(), outarray);
            }

            //Reverse data if necessary
            if(GetEorient(edge) == StdRegions::eBackwards)
            {
                Vmath::Reverse(EdgeExp->GetNumPoints(0),&outarray[0], 1,
                               &outarray[0], 1);
            }
        }

        void QuadExp::v_GetEdgeInterpVals(const int edge,
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD,NekDouble> &outarray)
        {
             int i;
             int nq0 = m_base[0]->GetNumPoints();
             int nq1 = m_base[1]->GetNumPoints();

             StdRegions::ConstFactorMap factors;
             factors[StdRegions::eFactorGaussEdge] = edge;

             StdRegions::StdMatrixKey key(
                 StdRegions::eInterpGauss,
                 DetShapeType(),*this,factors);

             DNekScalMatSharedPtr mat_gauss = m_matrixManager[key];

             switch (edge)
             {
                 case 0:
                 {
                     for (i = 0; i < nq0; i++)
                     {
                         outarray[i] = Blas::Ddot(
                             nq1, mat_gauss->GetOwnedMatrix()->GetPtr().get(),
                             1, &inarray[i], nq0);
                     }
                     break;
                 }
                 case 1:
                 {
                     for (i = 0; i < nq1; i++)
                     {
                         outarray[i] =  Blas::Ddot(
                             nq0, mat_gauss->GetOwnedMatrix()->GetPtr().get(),
                             1, &inarray[i * nq0], 1);
                     }
                     break;
                 }
                 case 2:
                 {
                     for (i = 0; i < nq0; i++)
                     {
                         outarray[i] = Blas::Ddot(
                             nq1, mat_gauss->GetOwnedMatrix()->GetPtr().get(),
                             1, &inarray[i], nq0);
                     }
                     break;
                 }
                 case 3:
                 {
                     for (i = 0; i < nq1; i++)
                     {
                         outarray[i] =  Blas::Ddot(
                             nq0, mat_gauss->GetOwnedMatrix()->GetPtr().get(),
                             1, &inarray[i * nq0], 1);
                     }
                     break;
                 }
                 default:
                     ASSERTL0(false, "edge value (< 3) is out of range");
                     break;
             }
        }


        void QuadExp::v_GetEdgePhysMap(
            const int                edge,
            Array<OneD, int>        &outarray)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();

            // Get points in Cartesian orientation
            switch (edge)
            {
                case 0:
                    outarray = Array<OneD, int>(nquad0);
                    for (int i = 0; i < nquad0; ++i)
                    {
                        outarray[i] = i;
                    }
                    break;
                case 1:
                    outarray = Array<OneD, int>(nquad1);
                    for (int i = 0; i < nquad1; ++i)
                    {
                        outarray[i] = (nquad0-1) + i*nquad0;
                    }
                    break;
                case 2:
                    outarray = Array<OneD, int>(nquad0);
                    for (int i = 0; i < nquad0; ++i)
                    {
                        outarray[i] = i + nquad0*(nquad1-1);
                    }
                    break;
                case 3:
                    outarray = Array<OneD, int>(nquad1);
                    for (int i = 0; i < nquad1; ++i)
                    {
                        outarray[i] = i*nquad0;
                    }
                    break;
                default:
                    ASSERTL0(false, "edge value (< 3) is out of range");
                    break;
            }

        }




        void QuadExp::v_GetEdgeQFactors(
                const int edge,
                Array<OneD, NekDouble> &outarray)
        {
            int i;
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();

            LibUtilities::PointsKeyVector ptsKeys = GetPointsKeys();
            const Array<OneD, const NekDouble>& jac = m_metricinfo->GetJac(ptsKeys);
            const Array<TwoD, const NekDouble>& df  = m_metricinfo->GetDerivFactors(ptsKeys);

            Array<OneD, NekDouble> j (max(nquad0, nquad1), 0.0);
            Array<OneD, NekDouble> g0(max(nquad0, nquad1), 0.0);
            Array<OneD, NekDouble> g1(max(nquad0, nquad1), 0.0);
            Array<OneD, NekDouble> g2(max(nquad0, nquad1), 0.0);
            Array<OneD, NekDouble> g3(max(nquad0, nquad1), 0.0);

            if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                // Implementation for all the basis except Gauss points
                if (m_base[0]->GetPointsType()
                    != LibUtilities::eGaussGaussLegendre
                    && m_base[1]->GetPointsType() !=
                    LibUtilities::eGaussGaussLegendre)
                {
                    switch (edge)
                    {
                        case 0:
                            Vmath::Vcopy(nquad0, &(df[1][0]),
                                         1, &(g1[0]), 1);
                            Vmath::Vcopy(nquad0, &(df[3][0]),
                                         1, &(g3[0]), 1);
                            Vmath::Vcopy(nquad0, &(jac[0]),1, &(j[0]),  1);

                            for (i = 0; i < nquad0; ++i)
                            {
                                outarray[i] = j[i]*sqrt(g1[i]*g1[i]
                                                                 + g3[i]*g3[i]);
                            }
                            break;
                        case 1:
                            Vmath::Vcopy(nquad1,
                                         &(df[0][0])+(nquad0-1), nquad0,
                                         &(g0[0]), 1);

                            Vmath::Vcopy(nquad1,
                                         &(df[2][0])+(nquad0-1), nquad0,
                                         &(g2[0]), 1);

                            Vmath::Vcopy(nquad1,
                                         &(jac[0])+(nquad0-1), nquad0,
                                         &(j[0]), 1);

                            for (i = 0; i < nquad1; ++i)
                            {
                                outarray[i] = j[i]*sqrt(g0[i]*g0[i] +
                                                                   g2[i]*g2[i]);
                            }
                            break;
                        case 2:

                            Vmath::Vcopy(nquad0,
                                         &(df[1][0])+(nquad0*(nquad1-1)), 1,
                                         &(g1[0]), 1);

                            Vmath::Vcopy(nquad0,
                                         &(df[3][0])+(nquad0*(nquad1-1)), 1,
                                         &(g3[0]), 1);

                            Vmath::Vcopy(nquad0,
                                         &(jac[0])+(nquad0*(nquad1-1)), 1,
                                         &(j[0]), 1);

                            for (i = 0; i < nquad0; ++i)
                            {
                                outarray[i] =
                                    j[i]*sqrt(g1[i]*g1[i]+ g3[i]*g3[i]);
                            }
                            break;
                        case 3:

                            Vmath::Vcopy(nquad1, &(df[0][0]), nquad0,&(g0[0]), 1);
                            Vmath::Vcopy(nquad1, &(df[2][0]), nquad0,&(g2[0]), 1);
                            Vmath::Vcopy(nquad1, &(jac[0]), nquad0, &(j[0]), 1);

                            for (i = 0; i < nquad1; ++i)
                            {
                                outarray[i] = j[i]*sqrt(g0[i]*g0[i] +
                                                        g2[i]*g2[i]);
                            }
                            break;
                        default:
                            ASSERTL0(false,"edge value (< 3) is out of range");
                            break;
                    }
                }
                else
                {
                    int nqtot =  nquad0 * nquad1;
                    Array<OneD, NekDouble> tmp_gmat0(nqtot, 0.0);
                    Array<OneD, NekDouble> tmp_gmat1(nqtot, 0.0);
                    Array<OneD, NekDouble> tmp_gmat2(nqtot, 0.0);
                    Array<OneD, NekDouble> tmp_gmat3(nqtot, 0.0);
                    Array<OneD, NekDouble> g0_edge(max(nquad0, nquad1), 0.0);
                    Array<OneD, NekDouble> g1_edge(max(nquad0, nquad1), 0.0);
                    Array<OneD, NekDouble> g2_edge(max(nquad0, nquad1), 0.0);
                    Array<OneD, NekDouble> g3_edge(max(nquad0, nquad1), 0.0);
                    Array<OneD, NekDouble> jac_edge(max(nquad0, nquad1), 0.0);

                    switch (edge)
                    {
                        case 0:
                            Vmath::Vmul(nqtot,&(df[1][0]),1,&jac[0],1,
                                        &(tmp_gmat1[0]),1);
                            Vmath::Vmul(nqtot,&(df[3][0]),1,&jac[0],1,
                                        &(tmp_gmat3[0]),1);
                            QuadExp::v_GetEdgeInterpVals(
                                edge, tmp_gmat1, g1_edge);
                            QuadExp::v_GetEdgeInterpVals(
                                edge, tmp_gmat3, g3_edge);

                            for (i = 0; i < nquad0; ++i)
                            {
                                outarray[i] = sqrt(g1_edge[i]*g1_edge[i] +
                                                   g3_edge[i]*g3_edge[i]);
                            }
                            break;

                        case 1:
                            Vmath::Vmul(nqtot,
                                        &(df[0][0]), 1,
                                        &jac[0], 1,
                                        &(tmp_gmat0[0]), 1);
                            Vmath::Vmul(nqtot,
                                        &(df[2][0]), 1,
                                        &jac[0], 1,
                                        &(tmp_gmat2[0]),
                                        1);
                            QuadExp::v_GetEdgeInterpVals(
                                edge, tmp_gmat0, g0_edge);
                            QuadExp::v_GetEdgeInterpVals(
                                edge, tmp_gmat2, g2_edge);

                            for (i = 0; i < nquad1; ++i)
                            {
                                outarray[i] = sqrt(g0_edge[i]*g0_edge[i]
                                                   + g2_edge[i]*g2_edge[i]);
                            }

                            break;
                        case 2:

                            Vmath::Vmul(nqtot,
                                        &(df[1][0]), 1,
                                        &jac[0], 1,
                                        &(tmp_gmat1[0]), 1);
                            Vmath::Vmul(nqtot,
                                        &(df[3][0]), 1,
                                        &jac[0], 1,
                                        &(tmp_gmat3[0]),1);
                            QuadExp::v_GetEdgeInterpVals(
                                edge, tmp_gmat1, g1_edge);
                            QuadExp::v_GetEdgeInterpVals(
                                edge, tmp_gmat3, g3_edge);


                            for (i = 0; i < nquad0; ++i)
                            {
                                outarray[i] = sqrt(g1_edge[i]*g1_edge[i]
                                                   + g3_edge[i]*g3_edge[i]);
                            }

                            Vmath::Reverse(nquad0,&outarray[0],1,&outarray[0],1);

                            break;
                        case 3:
                            Vmath::Vmul(nqtot,
                                        &(df[0][0]), 1,
                                        &jac[0], 1,
                                        &(tmp_gmat0[0]), 1);
                            Vmath::Vmul(nqtot,
                                        &(df[2][0]),1,
                                        &jac[0], 1,
                                        &(tmp_gmat2[0]),1);
                            QuadExp::v_GetEdgeInterpVals(
                                edge, tmp_gmat0, g0_edge);
                            QuadExp::v_GetEdgeInterpVals(
                                edge, tmp_gmat2, g2_edge);


                            for (i = 0; i < nquad1; ++i)
                            {
                                outarray[i] = sqrt(g0_edge[i]*g0_edge[i] +
                                                   g2_edge[i]*g2_edge[i]);
                            }

                            Vmath::Reverse(nquad1,
                                           &outarray[0], 1,
                                           &outarray[0], 1);

                            break;
                        default:
                            ASSERTL0(false,"edge value (< 3) is out of range");
                            break;
                    }
                }
            }
            else
            {

                switch (edge)
                {
                    case 0:



                        for (i = 0; i < nquad0; ++i)
                        {
                            outarray[i] = jac[0]*sqrt(df[1][0]*df[1][0] +
                                                      df[3][0]*df[3][0]);
                        }
                        break;
                    case 1:
                        for (i = 0; i < nquad1; ++i)
                        {
                            outarray[i] = jac[0]*sqrt(df[0][0]*df[0][0] +
                                                      df[2][0]*df[2][0]);
                        }
                        break;
                    case 2:
                        for (i = 0; i < nquad0; ++i)
                        {
                            outarray[i] = jac[0]*sqrt(df[1][0]*df[1][0] +
                                                      df[3][0]*df[3][0]);
                        }
                        break;
                    case 3:
                        for (i = 0; i < nquad1; ++i)
                        {
                            outarray[i] = jac[0]*sqrt(df[0][0]*df[0][0] +
                                                      df[2][0]*df[2][0]);
                        }
                        break;
                    default:
                        ASSERTL0(false,"edge value (< 3) is out of range");
                        break;
                }
            }
        }


        void QuadExp::v_ComputeEdgeNormal(const int edge)
        {
            int i;
            const SpatialDomains::GeomFactorsSharedPtr & geomFactors =
            GetGeom()->GetMetricInfo();
            SpatialDomains::GeomType type = geomFactors->GetGtype();

            LibUtilities::PointsKeyVector ptsKeys = GetPointsKeys();
            for(i = 0; i < ptsKeys.size(); ++i)
            {
                // Need at least 2 points for computing normals
                if (ptsKeys[i].GetNumPoints() == 1)
                {
                    LibUtilities::PointsKey pKey(2, ptsKeys[i].GetPointsType());
                    ptsKeys[i] = pKey;
                }
            }

            const Array<TwoD, const NekDouble> & df = geomFactors->GetDerivFactors(ptsKeys);
            const Array<OneD, const NekDouble> & jac  = geomFactors->GetJac(ptsKeys);
            int nqe;
            if (edge == 0 || edge == 2)
            {
                nqe = m_base[0]->GetNumPoints();
            }
            else
            {
                nqe = m_base[1]->GetNumPoints();
            }
            int vCoordDim = GetCoordim();

            m_edgeNormals[edge] = Array<OneD, Array<OneD, NekDouble> >
                                                                    (vCoordDim);
            Array<OneD, Array<OneD, NekDouble> > &normal = m_edgeNormals[edge];
            for (i = 0; i < vCoordDim; ++i)
            {
                normal[i] = Array<OneD, NekDouble>(nqe);
            }

            size_t nqb = nqe;
            size_t nbnd= edge;
            m_elmtBndNormDirElmtLen[nbnd] = Array<OneD, NekDouble> {nqb, 0.0};
            Array<OneD, NekDouble> &length = m_elmtBndNormDirElmtLen[nbnd];

            // Regular geometry case
            if ((type == SpatialDomains::eRegular)||
               (type == SpatialDomains::eMovingRegular))
            {
                NekDouble fac;
                // Set up normals
                switch (edge)
                {
                    case 0:
                        for (i = 0; i < vCoordDim; ++i)
                        {
                            Vmath::Fill(nqe, -df[2*i+1][0], normal[i], 1);
                        }
                        break;
                    case 1:
                        for (i = 0; i < vCoordDim; ++i)
                        {
                            Vmath::Fill(nqe, df[2*i][0], normal[i], 1);
                        }
                        break;
                    case 2:
                        for (i = 0; i < vCoordDim; ++i)
                        {
                            Vmath::Fill(nqe, df[2*i+1][0], normal[i], 1);
                        }
                        break;
                    case 3:
                        for (i = 0; i < vCoordDim; ++i)
                        {
                            Vmath::Fill(nqe, -df[2*i][0], normal[i], 1);
                        }
                        break;
                    default:
                        ASSERTL0(false, "edge is out of range (edge < 4)");
                }

                // normalise
                fac = 0.0;
                for (i =0 ; i < vCoordDim; ++i)
                {
                    fac += normal[i][0]*normal[i][0];
                }
                fac = 1.0/sqrt(fac);

                Vmath::Fill(nqb, fac, length, 1);

                for (i = 0; i < vCoordDim; ++i)
                {
                    Vmath::Smul(nqe, fac, normal[i], 1,normal[i], 1);
                }
            }
            else   // Set up deformed normals
            {
                int j;

                int nquad0 = ptsKeys[0].GetNumPoints();
                int nquad1 = ptsKeys[1].GetNumPoints();

                LibUtilities::PointsKey from_key;

                Array<OneD,NekDouble> normals(vCoordDim*max(nquad0,nquad1),0.0);
                Array<OneD,NekDouble> edgejac(vCoordDim*max(nquad0,nquad1),0.0);

                // Extract Jacobian along edges and recover local
                // derivates (dx/dr) for polynomial interpolation by
                // multiplying m_gmat by jacobian

                // Implementation for all the basis except Gauss points
                if (m_base[0]->GetPointsType() !=
                   LibUtilities::eGaussGaussLegendre
                   && m_base[1]->GetPointsType() !=
                   LibUtilities::eGaussGaussLegendre)
                {
                    switch (edge)
                    {
                        case 0:
                            for (j = 0; j < nquad0; ++j)
                            {
                                edgejac[j] = jac[j];
                                for (i = 0; i < vCoordDim; ++i)
                                {
                                    normals[i*nquad0+j] =
                                        -df[2*i+1][j]*edgejac[j];
                                }
                            }
                            from_key = ptsKeys[0];
                            break;
                        case 1:
                            for (j = 0; j < nquad1; ++j)
                            {
                                edgejac[j] = jac[nquad0*j+nquad0-1];
                                for (i = 0; i < vCoordDim; ++i)
                                {
                                    normals[i*nquad1+j]  =
                                        df[2*i][nquad0*j + nquad0-1]
                                        *edgejac[j];
                                }
                            }
                            from_key = ptsKeys[1];
                            break;
                        case 2:
                            for (j = 0; j < nquad0; ++j)
                            {
                                edgejac[j] = jac[nquad0*(nquad1-1)+j];
                                for (i = 0; i < vCoordDim; ++i)
                                {
                                    normals[i*nquad0+j] =
                                        (df[2*i+1][nquad0*(nquad1-1)+j])
                                        *edgejac[j];
                                }
                            }
                            from_key = ptsKeys[0];
                            break;
                        case 3:
                            for (j = 0; j < nquad1; ++j)
                            {
                                edgejac[j] = jac[nquad0*j];
                                for (i = 0; i < vCoordDim; ++i)
                                {
                                    normals[i*nquad1+j] =
                                        -df[2*i][nquad0*j]*edgejac[j];
                                }
                            }
                            from_key = ptsKeys[1];
                            break;
                        default:
                            ASSERTL0(false,"edge is out of range (edge < 3)");
                    }
                }
                else
                {
                    int nqtot =  nquad0 * nquad1;
                    Array<OneD,  NekDouble> tmp_gmat(nqtot, 0.0);
                    Array<OneD,  NekDouble> tmp_gmat_edge(nqe, 0.0);

                    switch (edge)
                    {
                        case 0:
                            for (j = 0; j < nquad0; ++j)
                            {
                                for (i = 0; i < vCoordDim; ++i)
                                {
                                    Vmath::Vmul(nqtot,
                                                &(df[2*i+1][0]), 1,
                                                &jac[0], 1,
                                                &(tmp_gmat[0]), 1);
                                    QuadExp::v_GetEdgeInterpVals(
                                        edge, tmp_gmat, tmp_gmat_edge);
                                    normals[i*nquad0+j] = -tmp_gmat_edge[j];
                                }
                            }
                            from_key = ptsKeys[0];
                            break;
                        case 1:
                            for (j = 0; j < nquad1; ++j)
                            {
                                for (i = 0; i < vCoordDim; ++i)
                                {
                                    Vmath::Vmul(nqtot,
                                                &(df[2*i][0]), 1,
                                                &jac[0], 1,
                                                &(tmp_gmat[0]), 1);
                                    QuadExp::v_GetEdgeInterpVals(
                                        edge, tmp_gmat, tmp_gmat_edge);
                                    normals[i*nquad1+j]  = tmp_gmat_edge[j];
                                }
                            }
                            from_key = ptsKeys[1];
                            break;
                        case 2:
                            for (j = 0; j < nquad0; ++j)
                            {
                                for (i = 0; i < vCoordDim; ++i)
                                {
                                    Vmath::Vmul(nqtot,
                                                &(df[2*i+1][0]), 1,
                                                &jac[0], 1,
                                                &(tmp_gmat[0]), 1);
                                    QuadExp::v_GetEdgeInterpVals(
                                        edge, tmp_gmat, tmp_gmat_edge);
                                    normals[i*nquad0+j] = tmp_gmat_edge[j];
                                }
                            }
                            from_key = ptsKeys[0];
                            break;
                        case 3:
                            for (j = 0; j < nquad1; ++j)
                            {
                                for (i = 0; i < vCoordDim; ++i)
                                {
                                    Vmath::Vmul(nqtot,
                                                &(df[2*i][0]), 1,
                                                &jac[0], 1,
                                                &(tmp_gmat[0]) ,1);
                                    QuadExp::v_GetEdgeInterpVals(
                                        edge, tmp_gmat, tmp_gmat_edge);
                                    normals[i*nquad1+j] = -tmp_gmat_edge[j];
                                }
                            }
                            from_key = ptsKeys[1];
                            break;
                        default:
                            ASSERTL0(false,"edge is out of range (edge < 3)");
                    }
                }

                int nq  = from_key.GetNumPoints();
                Array<OneD,NekDouble> work(nqe,0.0);

                // interpolate Jacobian and invert
                LibUtilities::Interp1D(
                    from_key,jac, m_base[0]->GetPointsKey(), work);
                Vmath::Sdiv(nqe,1.0,&work[0],1,&work[0],1);

                // interpolate
                for (i = 0; i < GetCoordim(); ++i)
                {
                    LibUtilities::Interp1D(
                        from_key,&normals[i*nq],
                        m_base[0]->GetPointsKey(),
                        &normal[i][0]);
                    Vmath::Vmul(nqe, work, 1, normal[i], 1, normal[i], 1);
                }

                //normalise normal vectors
                Vmath::Zero(nqe,work,1);
                for (i = 0; i < GetCoordim(); ++i)
                {
                    Vmath::Vvtvp(nqe,
                                 normal[i], 1,
                                 normal[i],1 ,
                                 work, 1,
                                 work, 1);
                }

                Vmath::Vsqrt(nqe,work,1,work,1);
                Vmath::Sdiv(nqe,1.0,work,1,work,1);

                Vmath::Vcopy(nqb, work, 1, length, 1);

                for (i = 0; i < GetCoordim(); ++i)
                {
                    Vmath::Vmul(nqe, normal[i], 1, work, 1, normal[i], 1);
                }
            }
            if (GetGeom()->GetEorient(edge) == StdRegions::eBackwards)
            {
                for (i = 0; i < vCoordDim; ++i)
                {
                    if (geomFactors->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Reverse(nqe, normal[i], 1, normal[i],1);
                    }
                }
            }
        }

        const SpatialDomains::GeomFactorsSharedPtr& QuadExp::v_GetMetricInfo() const
        {
            return m_metricinfo;
        }


        int QuadExp::v_GetCoordim()
        {
            return m_geom->GetCoordim();
        }


        void QuadExp::v_ExtractDataToCoeffs(
            const NekDouble *data,
            const std::vector<unsigned int > &nummodes,
            int mode_offset,
            NekDouble *coeffs,
            std::vector<LibUtilities::BasisType> &fromType)
        {
            int data_order0 = nummodes[mode_offset];
            int fillorder0  = std::min(m_base[0]->GetNumModes(),data_order0);

            int data_order1 = nummodes[mode_offset + 1];
            int order1      = m_base[1]->GetNumModes();
            int fillorder1  = min(order1,data_order1);

            // Check if same basis
            if (fromType[0] != m_base[0]->GetBasisType() ||
                fromType[1] != m_base[1]->GetBasisType())
            {
                // Construct a quad with the appropriate basis type at our
                // quadrature points, and one more to do a forwards
                // transform. We can then copy the output to coeffs.
                StdRegions::StdQuadExp tmpQuad(
                    LibUtilities::BasisKey(
                        fromType[0], data_order0, m_base[0]->GetPointsKey()),
                    LibUtilities::BasisKey(
                        fromType[1], data_order1, m_base[1]->GetPointsKey()));
                StdRegions::StdQuadExp tmpQuad2(m_base[0]->GetBasisKey(),
                                                m_base[1]->GetBasisKey());

                Array<OneD, const NekDouble> tmpData(tmpQuad.GetNcoeffs(), data);
                Array<OneD, NekDouble> tmpBwd(tmpQuad2.GetTotPoints());
                Array<OneD, NekDouble> tmpOut(tmpQuad2.GetNcoeffs());

                tmpQuad.BwdTrans(tmpData, tmpBwd);
                tmpQuad2.FwdTrans(tmpBwd, tmpOut);
                Vmath::Vcopy(tmpOut.size(), &tmpOut[0], 1, coeffs, 1);

                return;
            }

            switch (m_base[0]->GetBasisType())
            {
                case LibUtilities::eModified_A:
                {
                    int i;
                    int cnt = 0;
                    int cnt1 = 0;

                    ASSERTL1(m_base[1]->GetBasisType() ==
                            LibUtilities::eModified_A,
                            "Extraction routine not set up for this basis");

                    Vmath::Zero(m_ncoeffs,coeffs,1);
                    for (i = 0; i < fillorder0; ++i)
                    {
                        Vmath::Vcopy(fillorder1, data + cnt, 1, coeffs +cnt1, 1);
                        cnt  += data_order1;
                        cnt1 += order1;
                    }
                }
                    break;
                case LibUtilities::eGLL_Lagrange:
                {
                    LibUtilities::PointsKey
                        p0(nummodes[0], LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::PointsKey
                        p1(nummodes[1], LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::PointsKey t0(
                        m_base[0]->GetNumModes(),
                        LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::PointsKey t1(
                        m_base[1]->GetNumModes(),
                        LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::Interp2D(p0, p1, data, t0, t1, coeffs);
                }
                    break;
                case LibUtilities::eGauss_Lagrange:
                {
                    // Assume that input is also Gll_Lagrange but no way to check;
                    LibUtilities::PointsKey
                        p0(nummodes[0],LibUtilities::eGaussGaussLegendre);
                    LibUtilities::PointsKey
                        p1(nummodes[1],LibUtilities::eGaussGaussLegendre);
                    LibUtilities::PointsKey t0(
                        m_base[0]->GetNumModes(),
                        LibUtilities::eGaussGaussLegendre);
                    LibUtilities::PointsKey t1(
                        m_base[1]->GetNumModes(),
                        LibUtilities::eGaussGaussLegendre);
                    LibUtilities::Interp2D(p0, p1, data, t0, t1, coeffs);
                }
                    break;
                default:
                    ASSERTL0(false,
                    "basis is either not set up or not hierarchicial");
            }
        }


        StdRegions::Orientation QuadExp::v_GetEorient(int edge)
        {
            return m_geom->GetEorient(edge);
        }


        const LibUtilities::BasisSharedPtr& QuadExp::v_GetBasis(int dir) const
        {
            ASSERTL1(dir >= 0 &&dir <= 1, "input dir is out of range");
            return m_base[dir];
        }


        int QuadExp::v_GetNumPoints(const int dir) const
        {
            return GetNumPoints(dir);
        }

        DNekMatSharedPtr QuadExp::v_GenMatrix(
            const StdRegions::StdMatrixKey &mkey)
        {
            DNekMatSharedPtr returnval;
            switch (mkey.GetMatrixType())
            {
                case StdRegions::eHybridDGHelmholtz:
                case StdRegions::eHybridDGLamToU:
                case StdRegions::eHybridDGLamToQ0:
                case StdRegions::eHybridDGLamToQ1:
                case StdRegions::eHybridDGLamToQ2:
                case StdRegions::eHybridDGHelmBndLam:
				case StdRegions::eInvLaplacianWithUnityMean:
                    returnval = Expansion2D::v_GenMatrix(mkey);
                    break;
                default:
                    returnval = StdQuadExp::v_GenMatrix(mkey);
            }
            return returnval;
        }

        DNekMatSharedPtr QuadExp::v_CreateStdMatrix(
            const StdRegions::StdMatrixKey &mkey)
        {
            LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
            LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
            StdRegions::StdQuadExpSharedPtr tmp =
                MemoryManager<StdQuadExp>::AllocateSharedPtr(bkey0,bkey1);
            return tmp->GetStdMatrix(mkey);
        }


        DNekScalMatSharedPtr QuadExp::CreateMatrix(const MatrixKey &mkey)
        {
            DNekScalMatSharedPtr returnval;
            LibUtilities::PointsKeyVector ptsKeys = GetPointsKeys();

            ASSERTL2(m_metricinfo->GetGtype() != SpatialDomains::eNoGeomType,
                     "Geometric information is not set up");

            switch (mkey.GetMatrixType())
            {
                case StdRegions::eMass:
                {
                    if ((m_metricinfo->GetGtype() ==
                         SpatialDomains::eDeformed) || (mkey.GetNVarCoeff()))
                    {
                        NekDouble        one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>::
                            AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble        jac = (m_metricinfo->GetJac(ptsKeys))[0];
                        DNekMatSharedPtr mat = GetStdMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>::
                            AllocateSharedPtr(jac,mat);
                    }
                }
                    break;
                case StdRegions::eInvMass:
                {
                    if ((m_metricinfo->GetGtype() ==
                        SpatialDomains::eDeformed) || (mkey.GetNVarCoeff()))
                    {
                        NekDouble one = 1.0;
                        StdRegions::StdMatrixKey masskey(
                            StdRegions::eMass, DetShapeType(), *this);
                        DNekMatSharedPtr mat = GenMatrix(masskey);
                        mat->Invert();

                        returnval = MemoryManager<DNekScalMat>::
                            AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble fac = 1.0/(m_metricinfo->GetJac(ptsKeys))[0];
                        DNekMatSharedPtr mat = GetStdMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>::
                            AllocateSharedPtr(fac,mat);
                    }
                }
                    break;
                case StdRegions::eWeakDeriv0:
                case StdRegions::eWeakDeriv1:
                case StdRegions::eWeakDeriv2:
                {
                    if((m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                       || (mkey.GetNVarCoeff()))
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(mkey);

                        returnval = MemoryManager<DNekScalMat>::
                            AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble jac = (m_metricinfo->GetJac(ptsKeys))[0];
                        Array<TwoD, const NekDouble> df =
                            m_metricinfo->GetDerivFactors(ptsKeys);
                        int dir = 0;

                        switch(mkey.GetMatrixType())
                        {
                            case StdRegions::eWeakDeriv0:
                                dir = 0;
                                break;
                            case StdRegions::eWeakDeriv1:
                                dir = 1;
                                break;
                            case StdRegions::eWeakDeriv2:
                                dir = 2;
                                break;
                            default:
                                break;
                        }

                        MatrixKey deriv0key(StdRegions::eWeakDeriv0,
                                            mkey.GetShapeType(), *this);
                        MatrixKey deriv1key(StdRegions::eWeakDeriv1,
                                            mkey.GetShapeType(), *this);

                        DNekMat &deriv0 = *GetStdMatrix(deriv0key);
                        DNekMat &deriv1 = *GetStdMatrix(deriv1key);

                        int rows = deriv0.GetRows();
                        int cols = deriv1.GetColumns();

                        DNekMatSharedPtr WeakDeriv = MemoryManager<DNekMat>::
                            AllocateSharedPtr(rows,cols);
                        (*WeakDeriv) = df[2*dir][0]*deriv0 +
                                       df[2*dir+1][0]*deriv1;
                        returnval = MemoryManager<DNekScalMat>::
                            AllocateSharedPtr(jac,WeakDeriv);
                    }
                }
                    break;
                case StdRegions::eLaplacian:
                {
                    if( (m_metricinfo->GetGtype() ==
                         SpatialDomains::eDeformed) || (mkey.GetNVarCoeff() > 0)
                       || (mkey.ConstFactorExists
                           (StdRegions::eFactorSVVCutoffRatio)))
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(mkey);

                        returnval = MemoryManager<DNekScalMat>::
                            AllocateSharedPtr(one,mat);
                    }
                    else
                    {
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
                        Array<TwoD, const NekDouble>
                            gmat = m_metricinfo->GetGmat(ptsKeys);

                        int rows = lap00.GetRows();
                        int cols = lap00.GetColumns();

                        DNekMatSharedPtr lap =
                            MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);

                        (*lap) = gmat[0][0] * lap00 +
                                 gmat[1][0] * (lap01 + Transpose(lap01)) +
                                 gmat[3][0] * lap11;

                        returnval = MemoryManager<DNekScalMat>::
                            AllocateSharedPtr(jac,lap);
                    }
                }
                    break;
                case StdRegions::eInvLaplacianWithUnityMean:
                {
                    DNekMatSharedPtr mat = GenMatrix(mkey);
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(1.0,mat);
                }
                    break;
                case StdRegions::eHelmholtz:
                {
                    NekDouble lambda =
                        mkey.GetConstFactor(StdRegions::eFactorLambda);

                    MatrixKey masskey(mkey, StdRegions::eMass);
                    DNekScalMat &MassMat = *(this->m_matrixManager[masskey]);

                    MatrixKey lapkey(mkey, StdRegions::eLaplacian);
                    DNekScalMat &LapMat = *(this->m_matrixManager[lapkey]);

                    int rows = LapMat.GetRows();
                    int cols = LapMat.GetColumns();

                    DNekMatSharedPtr helm = MemoryManager<DNekMat>::
                        AllocateSharedPtr(rows,cols);

                    NekDouble one = 1.0;
                    (*helm) = LapMat + lambda*MassMat;

                    returnval =
                        MemoryManager<DNekScalMat>::AllocateSharedPtr(one,helm);
                }
                    break;
                case StdRegions::eIProductWRTBase:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>::
                            AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble jac = (m_metricinfo->GetJac(ptsKeys))[0];
                        DNekMatSharedPtr mat = GetStdMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>::
                            AllocateSharedPtr(jac,mat);
                    }
                }
                    break;
                case StdRegions::eIProductWRTDerivBase0:
                case StdRegions::eIProductWRTDerivBase1:
                case StdRegions::eIProductWRTDerivBase2:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>::
                            AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble jac = (m_metricinfo->GetJac(ptsKeys))[0];
                        const Array<TwoD, const NekDouble>& df =
                                                        m_metricinfo->GetDerivFactors(ptsKeys);
                        int dir = 0;

                        switch(mkey.GetMatrixType())
                        {
                            case StdRegions::eIProductWRTDerivBase0:
                                dir = 0;
                                break;
                            case StdRegions::eIProductWRTDerivBase1:
                                dir = 1;
                                break;
                            case StdRegions::eIProductWRTDerivBase2:
                                dir = 2;
                                break;
                            default:
                                break;
                        }

                        MatrixKey iProdDeriv0Key(
                            StdRegions::eIProductWRTDerivBase0,
                            mkey.GetShapeType(), *this);
                        MatrixKey iProdDeriv1Key(
                            StdRegions::eIProductWRTDerivBase1,
                            mkey.GetShapeType(), *this);

                        DNekMat &stdiprod0 = *GetStdMatrix(iProdDeriv0Key);
                        DNekMat &stdiprod1 = *GetStdMatrix(iProdDeriv0Key);

                        int rows = stdiprod0.GetRows();
                        int cols = stdiprod1.GetColumns();

                        DNekMatSharedPtr mat = MemoryManager<DNekMat>::
                            AllocateSharedPtr(rows,cols);
                        (*mat) = df[2*dir][0]*stdiprod0 +
                                 df[2*dir+1][0]*stdiprod1;

                        returnval = MemoryManager<DNekScalMat>::
                            AllocateSharedPtr(jac,mat);
                    }
                }
                    break;
                case StdRegions::eInvHybridDGHelmholtz:
                {
                    NekDouble one = 1.0;

                    MatrixKey hkey(StdRegions::eHybridDGHelmholtz,
                                   DetShapeType(), *this,
                                   mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekMatSharedPtr mat = GenMatrix(hkey);

                    mat->Invert();
                    returnval =
                         MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                }
                    break;
                case StdRegions::eInterpGauss:
                {
                    DNekMatSharedPtr m_Ix;
                    Array<OneD, NekDouble> coords(1, 0.0);
                    StdRegions::ConstFactorMap factors = mkey.GetConstFactors();
                    int edge = (int)factors[StdRegions::eFactorGaussEdge];

                    coords[0] = (edge == 0 || edge == 3) ? -1.0 : 1.0;

                    m_Ix = m_base[(edge + 1) % 2]->GetI(coords);
                    returnval =
                        MemoryManager<DNekScalMat>::AllocateSharedPtr(1.0,m_Ix);
                }
                    break;
                case StdRegions::ePreconLinearSpace:
                {
                    NekDouble one = 1.0;
                    MatrixKey helmkey(
                        StdRegions::eHelmholtz, mkey.GetShapeType(), *this,
                        mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekScalBlkMatSharedPtr helmStatCond =
                        GetLocStaticCondMatrix(helmkey);
                    DNekScalMatSharedPtr A =helmStatCond->GetBlock(0,0);
                    DNekMatSharedPtr R=BuildVertexMatrix(A);

                    returnval =
                        MemoryManager<DNekScalMat>::AllocateSharedPtr(one, R);
                }
                break;
            default:
                {
                    NekDouble        one = 1.0;
                    DNekMatSharedPtr mat = GenMatrix(mkey);

                    returnval =
                        MemoryManager<DNekScalMat>::AllocateSharedPtr(one, mat);
                }
                break;
            }

            return returnval;
        }


        DNekScalBlkMatSharedPtr QuadExp::CreateStaticCondMatrix(
            const MatrixKey &mkey)
        {
            DNekScalBlkMatSharedPtr returnval;

            ASSERTL2(m_metricinfo->GetGtype()
                     != SpatialDomains::eNoGeomType,
                     "Geometric information is not set up");

            // set up block matrix system
            unsigned int nbdry = NumBndryCoeffs();
            unsigned int nint = (unsigned int)(m_ncoeffs - nbdry);
            unsigned int exp_size[] = {nbdry,nint};
            unsigned int nblks = 2;
            returnval = MemoryManager<DNekScalBlkMat>::
                AllocateSharedPtr(nblks,nblks,exp_size,exp_size);
                //Really need a constructor which takes Arrays
            NekDouble factor = 1.0;

            switch (mkey.GetMatrixType())
            {
                // this can only use stdregions statically condensed system
                // for mass matrix
                case StdRegions::eMass:
                    if ((m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                        ||(mkey.GetNVarCoeff()))
                    {
                        factor = 1.0;
                        goto UseLocRegionsMatrix;
                    }
                    else
                    {
                        factor = (m_metricinfo->GetJac(GetPointsKeys()))[0];
                        goto UseStdRegionsMatrix;
                    }
                    break;
                default: // use Deformed case for both
                        // regular and deformed geometries
                    factor = 1.0;
                    goto UseLocRegionsMatrix;
                    break;
                UseStdRegionsMatrix:
                {
                    NekDouble            invfactor = 1.0/factor;
                    NekDouble            one = 1.0;
                    DNekBlkMatSharedPtr  mat = GetStdStaticCondMatrix(mkey);
                    DNekScalMatSharedPtr Atmp;
                    DNekMatSharedPtr     Asubmat;

                    returnval->SetBlock(0,0,Atmp = MemoryManager<DNekScalMat>::
                        AllocateSharedPtr(factor,Asubmat = mat->GetBlock(0,0)));
                    returnval->SetBlock(0,1,Atmp = MemoryManager<DNekScalMat>::
                        AllocateSharedPtr(one,Asubmat = mat->GetBlock(0,1)));
                    returnval->SetBlock(1,0,Atmp = MemoryManager<DNekScalMat>::
                        AllocateSharedPtr(factor,Asubmat = mat->GetBlock(1,0)));
                    returnval->SetBlock(1,1,Atmp = MemoryManager<DNekScalMat>::
                        AllocateSharedPtr(invfactor,Asubmat = mat->GetBlock(1,1)));
                }
                    break;
                UseLocRegionsMatrix:
                {
                    int i,j;
                    NekDouble            invfactor = 1.0/factor;
                    NekDouble            one = 1.0;
                    DNekScalMat &mat = *GetLocMatrix(mkey);
                    DNekMatSharedPtr A = MemoryManager<DNekMat>::
                        AllocateSharedPtr(nbdry,nbdry);
                    DNekMatSharedPtr B = MemoryManager<DNekMat>::
                        AllocateSharedPtr(nbdry,nint);
                    DNekMatSharedPtr C = MemoryManager<DNekMat>::
                        AllocateSharedPtr(nint,nbdry);
                    DNekMatSharedPtr D = MemoryManager<DNekMat>::
                        AllocateSharedPtr(nint,nint);

                    Array<OneD,unsigned int> bmap(nbdry);
                    Array<OneD,unsigned int> imap(nint);
                    GetBoundaryMap(bmap);
                    GetInteriorMap(imap);

                    for (i = 0; i < nbdry; ++i)
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

                    for (i = 0; i < nint; ++i)
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

                    returnval->SetBlock(0,0,Atmp = MemoryManager<DNekScalMat>::
                                                AllocateSharedPtr(factor, A));
                    returnval->SetBlock(0,1,Atmp = MemoryManager<DNekScalMat>::
                                                AllocateSharedPtr(one, B));
                    returnval->SetBlock(1,0,Atmp = MemoryManager<DNekScalMat>::
                                                AllocateSharedPtr(factor, C));
                    returnval->SetBlock(1,1,Atmp = MemoryManager<DNekScalMat>::
                                                AllocateSharedPtr(invfactor, D));

                }
            }
            return returnval;
        }


        DNekScalMatSharedPtr QuadExp::v_GetLocMatrix(const MatrixKey &mkey)
        {
            return m_matrixManager[mkey];
        }


        DNekScalBlkMatSharedPtr QuadExp::v_GetLocStaticCondMatrix(
                                                          const MatrixKey &mkey)
        {
            return m_staticCondMatrixManager[mkey];
        }

        void QuadExp::v_DropLocStaticCondMatrix(const MatrixKey &mkey)
        {
            m_staticCondMatrixManager.DeleteObject(mkey);
        }


        void QuadExp::v_MassMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,NekDouble> &outarray,
            const StdRegions::StdMatrixKey &mkey)
        {
            StdExpansion::MassMatrixOp_MatFree(inarray, outarray, mkey);
        }


        void QuadExp::v_LaplacianMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,NekDouble> &outarray,
            const StdRegions::StdMatrixKey &mkey)
        {
            QuadExp::LaplacianMatrixOp_MatFree(inarray, outarray, mkey);
        }


        void QuadExp::v_LaplacianMatrixOp(
            const int k1,
            const int k2,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,NekDouble> &outarray,
            const StdRegions::StdMatrixKey &mkey)
        {
            StdExpansion::LaplacianMatrixOp_MatFree(
                k1, k2, inarray, outarray, mkey);
        }


        void QuadExp::v_WeakDerivMatrixOp(
            const int i,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,NekDouble> &outarray,
            const StdRegions::StdMatrixKey &mkey)
        {
            StdExpansion::WeakDerivMatrixOp_MatFree(i, inarray, outarray, mkey);
        }


        void QuadExp::v_WeakDirectionalDerivMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,NekDouble> &outarray,
            const StdRegions::StdMatrixKey &mkey)
        {
            StdExpansion::WeakDirectionalDerivMatrixOp_MatFree(
                inarray, outarray, mkey);
        }


        void QuadExp::v_MassLevelCurvatureMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,NekDouble> &outarray,
            const StdRegions::StdMatrixKey &mkey)
        {
            StdExpansion::MassLevelCurvatureMatrixOp_MatFree(
                inarray, outarray, mkey);
        }


        void QuadExp::v_HelmholtzMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,NekDouble> &outarray,
            const StdRegions::StdMatrixKey &mkey)
        {
            QuadExp::HelmholtzMatrixOp_MatFree(inarray, outarray, mkey);
        }


        void QuadExp::v_GeneralMatrixOp_MatOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,NekDouble> &outarray,
            const StdRegions::StdMatrixKey &mkey)
        {
            MatrixKey newkey(mkey);
            DNekScalMatSharedPtr   mat = GetLocMatrix(newkey);

            if (inarray.get() == outarray.get())
            {
                Array<OneD,NekDouble> tmp(m_ncoeffs);
                Vmath::Vcopy(m_ncoeffs,inarray.get(),1,tmp.get(),1);

                Blas::Dgemv('N',m_ncoeffs,m_ncoeffs, mat->Scale(),
                            (mat->GetOwnedMatrix())->GetPtr().get(), m_ncoeffs,
                            tmp.get(), 1, 0.0, outarray.get(), 1);
            }
            else
            {
                Blas::Dgemv('N',m_ncoeffs,m_ncoeffs,mat->Scale(),
                            (mat->GetOwnedMatrix())->GetPtr().get(), m_ncoeffs,
                            inarray.get(), 1, 0.0, outarray.get(), 1);
            }
        }

        void QuadExp::v_ReduceOrderCoeffs(
            int                                 numMin,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray)
        {
            int n_coeffs = inarray.size();

            Array<OneD, NekDouble> coeff    (n_coeffs);
            Array<OneD, NekDouble> coeff_tmp(n_coeffs, 0.0);
            Array<OneD, NekDouble> tmp, tmp2;

            int nmodes0 = m_base[0]->GetNumModes();
            int nmodes1 = m_base[1]->GetNumModes();
            int numMax  = nmodes0;

            Vmath::Vcopy(n_coeffs,inarray,1,coeff_tmp,1);

            const LibUtilities::PointsKey Pkey0(
                nmodes0, LibUtilities::eGaussLobattoLegendre);
            const LibUtilities::PointsKey Pkey1(
                nmodes1, LibUtilities::eGaussLobattoLegendre);
            LibUtilities::BasisKey b0(
                m_base[0]->GetBasisType(), nmodes0, Pkey0);
            LibUtilities::BasisKey b1(
                m_base[1]->GetBasisType(), nmodes1, Pkey1);
            LibUtilities::BasisKey bortho0(
                LibUtilities::eOrtho_A,    nmodes0, Pkey0);
            LibUtilities::BasisKey bortho1(
                LibUtilities::eOrtho_A,    nmodes1, Pkey1);

            LibUtilities::InterpCoeff2D(
                b0, b1, coeff_tmp, bortho0, bortho1, coeff);

            Vmath::Zero(n_coeffs, coeff_tmp, 1);

            int cnt = 0;
            for (int i = 0; i < numMin+1; ++i)
            {
                Vmath::Vcopy(numMin,
                             tmp  = coeff+cnt,1,
                             tmp2 = coeff_tmp+cnt,1);

                cnt = i*numMax;
            }

            LibUtilities::InterpCoeff2D(
                bortho0, bortho1, coeff_tmp,
                b0,      b1,      outarray);
        }

        void QuadExp::v_LaplacianMatrixOp_MatFree_Kernel(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          Array<OneD,       NekDouble> &wsp)
        {
            if (m_metrics.count(eMetricLaplacian00) == 0)
            {
                ComputeLaplacianMetric();
            }

            int       nquad0  = m_base[0]->GetNumPoints();
            int       nquad1  = m_base[1]->GetNumPoints();
            int       nqtot   = nquad0*nquad1;
            int       nmodes0 = m_base[0]->GetNumModes();
            int       nmodes1 = m_base[1]->GetNumModes();
            int       wspsize = max(max(max(nqtot,m_ncoeffs),nquad1*nmodes0),nquad0*nmodes1);

            ASSERTL1(wsp.size() >= 3*wspsize,
                     "Workspace is of insufficient size.");

            const Array<OneD, const NekDouble>& base0  = m_base[0]->GetBdata();
            const Array<OneD, const NekDouble>& base1  = m_base[1]->GetBdata();
            const Array<OneD, const NekDouble>& dbase0 = m_base[0]->GetDbdata();
            const Array<OneD, const NekDouble>& dbase1 = m_base[1]->GetDbdata();
            const Array<OneD, const NekDouble>& metric00 = m_metrics[eMetricLaplacian00];
            const Array<OneD, const NekDouble>& metric01 = m_metrics[eMetricLaplacian01];
            const Array<OneD, const NekDouble>& metric11 = m_metrics[eMetricLaplacian11];

            // Allocate temporary storage
            Array<OneD,NekDouble> wsp0(wsp);
            Array<OneD,NekDouble> wsp1(wsp+wspsize);
            Array<OneD,NekDouble> wsp2(wsp+2*wspsize);

            StdExpansion2D::PhysTensorDeriv(inarray,wsp1,wsp2);

            // wsp0 = k = g0 * wsp1 + g1 * wsp2 = g0 * du_dxi1 + g1 * du_dxi2
            // wsp2 = l = g1 * wsp1 + g2 * wsp2 = g0 * du_dxi1 + g1 * du_dxi2
            // where g0, g1 and g2 are the metric terms set up in the GeomFactors class
            // especially for this purpose
            Vmath::Vvtvvtp(nqtot,&metric00[0],1,&wsp1[0],1,&metric01[0],1,&wsp2[0],1,&wsp0[0],1);
            Vmath::Vvtvvtp(nqtot,&metric01[0],1,&wsp1[0],1,&metric11[0],1,&wsp2[0],1,&wsp2[0],1);

            // outarray = m = (D_xi1 * B)^T * k
            // wsp1     = n = (D_xi2 * B)^T * l
            IProductWRTBase_SumFacKernel(dbase0,base1,wsp0,outarray,wsp1,false,true);
            IProductWRTBase_SumFacKernel(base0,dbase1,wsp2,wsp1,    wsp0,true,false);

            // outarray = outarray + wsp1
            //          = L * u_hat
            Vmath::Vadd(m_ncoeffs,wsp1.get(),1,outarray.get(),1,outarray.get(),1);
        }

        void QuadExp::v_ComputeLaplacianMetric()
        {
            if (m_metrics.count(eMetricQuadrature) == 0)
            {
                ComputeQuadratureMetric();
            }

            const SpatialDomains::GeomType type = m_metricinfo->GetGtype();
            const unsigned int nqtot = GetTotPoints();
            const unsigned int dim = 2;
            const MetricType m[3][3] = { {eMetricLaplacian00, eMetricLaplacian01, eMetricLaplacian02},
                                       {eMetricLaplacian01, eMetricLaplacian11, eMetricLaplacian12},
                                       {eMetricLaplacian02, eMetricLaplacian12, eMetricLaplacian22}
            };

            const Array<TwoD, const NekDouble> gmat =
                                    m_metricinfo->GetGmat(GetPointsKeys());
            for (unsigned int i = 0; i < dim; ++i)
            {
                for (unsigned int j = i; j < dim; ++j)
                {
                    m_metrics[m[i][j]] = Array<OneD, NekDouble>(nqtot);
                    if (type == SpatialDomains::eDeformed)
                    {
                        Vmath::Vcopy(nqtot, &gmat[i*dim+j][0], 1,
                                     &m_metrics[m[i][j]][0], 1);
                    }
                    else
                    {
                        Vmath::Fill(nqtot, gmat[i*dim+j][0],
                                    &m_metrics[m[i][j]][0], 1);
                    }
                    MultiplyByQuadratureMetric(m_metrics[m[i][j]],
                                               m_metrics[m[i][j]]);

                }
            }
        }

        void QuadExp::v_SVVLaplacianFilter(
                    Array<OneD, NekDouble> &array,
                    const StdRegions::StdMatrixKey &mkey)
        {
            int nq = GetTotPoints();

            // Calculate sqrt of the Jacobian
            Array<OneD, const NekDouble> jac =
                                    m_metricinfo->GetJac(GetPointsKeys());
            Array<OneD, NekDouble> sqrt_jac(nq);
            if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vsqrt(nq,jac,1,sqrt_jac,1);
            }
            else
            {
                Vmath::Fill(nq,sqrt(jac[0]),sqrt_jac,1);
            }

            // Multiply array by sqrt(Jac)
            Vmath::Vmul(nq,sqrt_jac,1,array,1,array,1);

            // Apply std region filter
            StdQuadExp::v_SVVLaplacianFilter( array, mkey);

            // Divide by sqrt(Jac)
            Vmath::Vdiv(nq,array,1,sqrt_jac,1,array,1);
        }

    }//end of namespace
}//end of namespace
