///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansion.cpp
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
// Description: Definition of methods in class StdExpansion which is
// the base class to all expansion shapes
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <StdRegions/StdExpansion.h>
#include <LibUtilities/Foundations/ManagerAccess.h>  // for BasisManager, etc

namespace Nektar
{
    namespace StdRegions
    {
        StdExpansion::StdExpansion(void):
            m_elmt_id(0),
            m_ncoeffs(0)
        {
        }

        StdExpansion::StdExpansion(const int numcoeffs, const int numbases,
            const LibUtilities::BasisKey &Ba,
            const LibUtilities::BasisKey &Bb,
            const LibUtilities::BasisKey &Bc):
            m_base(numbases),
            m_elmt_id(0),
            m_ncoeffs(numcoeffs),
            m_stdMatrixManager(
                std::bind(&StdExpansion::CreateStdMatrix, this, std::placeholders::_1),
                std::string("StdExpansionStdMatrix")),
            m_stdStaticCondMatrixManager(
                std::bind(&StdExpansion::CreateStdStaticCondMatrix, this, std::placeholders::_1),
                std::string("StdExpansionStdStaticCondMatrix")),
            m_IndexMapManager(
                std::bind(&StdExpansion::CreateIndexMap,this, std::placeholders::_1),
                std::string("StdExpansionIndexMap"))
        {
            switch(m_base.size())
            {
            case 3:
                ASSERTL2(Bc!=LibUtilities::NullBasisKey,
                    "NULL Basis attempting to be used.");
                m_base[2] = LibUtilities::BasisManager()[Bc];
                /* Falls through. */
            case 2:
                ASSERTL2(Bb!=LibUtilities::NullBasisKey,
                    "NULL Basis attempting to be used.");
                m_base[1] = LibUtilities::BasisManager()[Bb];
                /* Falls through. */
            case 1:
                ASSERTL2(Ba!=LibUtilities::NullBasisKey,
                    "NULL Basis attempting to be used.");
                m_base[0] = LibUtilities::BasisManager()[Ba];
                break;
            default:
                break;
//                ASSERTL0(false, "numbases incorrectly specified");
            };

        } //end constructor


        StdExpansion::StdExpansion(const StdExpansion &T):
            std::enable_shared_from_this<StdExpansion>(T),
            m_base(T.m_base),
            m_elmt_id(T.m_elmt_id),
            m_ncoeffs(T.m_ncoeffs),
            m_stdMatrixManager(T.m_stdMatrixManager),
            m_stdStaticCondMatrixManager(T.m_stdStaticCondMatrixManager),
            m_IndexMapManager(T.m_IndexMapManager)
        {
        }

        StdExpansion::~StdExpansion()
        {
        }

        NekDouble StdExpansion::Linf(const Array<OneD, const NekDouble>& phys,
                                     const Array<OneD, const NekDouble>& sol)
        {
            NekDouble  val;
            int     ntot = GetTotPoints();
            Array<OneD, NekDouble>  wsp(ntot);

            if(sol ==  NullNekDouble1DArray)
            {
                Vmath::Vabs(ntot, phys, 1, wsp, 1);
            }
            else
            {
                Vmath::Vsub(ntot, sol, 1, phys, 1, wsp, 1);
                Vmath::Vabs(ntot, wsp, 1, wsp, 1);
            }

            val = Vmath::Vamax(ntot, wsp, 1);

            return  val;
        }

        NekDouble StdExpansion::L2(const Array<OneD, const NekDouble>& phys,
                                   const Array<OneD, const NekDouble>& sol)
        {
            NekDouble  val;
            int     ntot = GetTotPoints();
            Array<OneD, NekDouble> wsp(ntot);

            if (sol.size() == 0)
            {
                Vmath::Vmul(ntot, phys, 1, phys, 1, wsp, 1);
            }
            else
            {
                Vmath::Vsub(ntot, sol, 1, phys, 1, wsp, 1);
                Vmath::Vmul(ntot, wsp, 1, wsp, 1, wsp, 1);
            }

            val = v_Integral(wsp);

            // if val too small, sqrt returns nan.
            if (fabs(val) < NekConstants::kNekSqrtTol*NekConstants::kNekSqrtTol)
            {
                return 0.0;
            }
            else
            {
                return sqrt(val);
            }
        }

        NekDouble StdExpansion::H1(const Array<OneD, const NekDouble>& phys,
                                   const Array<OneD, const NekDouble>& sol)
        {
            int         i;
            NekDouble  val;
            int     ntot = GetTotPoints();
            int     coordim = v_GetCoordim();
            Array<OneD, NekDouble> wsp(3*ntot);
            Array<OneD, NekDouble> wsp_deriv = wsp + ntot;
            Array<OneD, NekDouble> sum = wsp_deriv + ntot;

            if(sol ==  NullNekDouble1DArray)
            {
                Vmath::Vcopy(ntot,phys, 1, wsp, 1);
                Vmath::Vmul(ntot, phys, 1, phys, 1, sum, 1);
            }
            else
            {
                Vmath::Vsub(ntot, sol, 1, phys, 1, wsp, 1);
                Vmath::Vmul(ntot, wsp, 1, wsp, 1, sum, 1);
            }


            for(i = 0; i < coordim; ++i)
            {
                v_PhysDeriv(i,wsp,wsp_deriv);
                Vmath::Vvtvp(ntot,wsp_deriv,1,wsp_deriv,1,sum,1,sum,1);
            }

            val = sqrt(v_Integral(sum));

            return val;
        }


        DNekBlkMatSharedPtr StdExpansion::CreateStdStaticCondMatrix(const StdMatrixKey &mkey)
        {
            DNekBlkMatSharedPtr returnval;

            DNekMatSharedPtr  mat = GetStdMatrix(mkey);
            int nbdry = NumBndryCoeffs(); // also checks to see if this is a boundary interior decomposed expansion
            int nint = m_ncoeffs - nbdry;
            DNekMatSharedPtr A = MemoryManager<DNekMat>::AllocateSharedPtr(nbdry,nbdry);
            DNekMatSharedPtr B = MemoryManager<DNekMat>::AllocateSharedPtr(nbdry,nint);
            DNekMatSharedPtr C = MemoryManager<DNekMat>::AllocateSharedPtr(nint,nbdry);
            DNekMatSharedPtr D = MemoryManager<DNekMat>::AllocateSharedPtr(nint,nint);

            int i,j;

            Array<OneD,unsigned int> bmap(nbdry);
            Array<OneD,unsigned int> imap(nint);
            GetBoundaryMap(bmap);
            GetInteriorMap(imap);

            for(i = 0; i < nbdry; ++i)
            {
                for(j = 0; j < nbdry; ++j)
                {
                    (*A)(i,j) = (*mat)(bmap[i],bmap[j]);
                }

                for(j = 0; j < nint; ++j)
                {
                    (*B)(i,j) = (*mat)(bmap[i],imap[j]);
                }
            }

            for(i = 0; i < nint; ++i)
            {
                for(j = 0; j < nbdry; ++j)
                {
                    (*C)(i,j) = (*mat)(imap[i],bmap[j]);
                }

                for(j = 0; j < nint; ++j)
                {
                    (*D)(i,j) = (*mat)(imap[i],imap[j]);
                }
            }

            // Calculate static condensed system
            if(nint)
            {
                D->Invert();
                (*B) = (*B)*(*D);
                (*A) = (*A) - (*B)*(*C);
            }

            // set up block matrix system
            Array<OneD, unsigned int> exp_size(2);
            exp_size[0] = nbdry;
            exp_size[1] = nint;
            returnval = MemoryManager<DNekBlkMat>::AllocateSharedPtr(exp_size,exp_size);

            returnval->SetBlock(0,0,A);
            returnval->SetBlock(0,1,B);
            returnval->SetBlock(1,0,C);
            returnval->SetBlock(1,1,D);

            return returnval;
        }

        IndexMapValuesSharedPtr StdExpansion::CreateIndexMap(const IndexMapKey &ikey)
        {
            IndexMapValuesSharedPtr returnval;

            IndexMapType itype = ikey.GetIndexMapType();

            int entity = ikey.GetIndexEntity();

            Orientation orient = ikey.GetIndexOrientation();

            Array<OneD,unsigned int>     map;
            Array<OneD,int>             sign;

            switch(itype)
            {
                case eEdgeToElement:
                    {
                    v_GetEdgeToElementMap(entity,orient,map,sign);
                }
                break;
                case eFaceToElement:
                {
                    v_GetFaceToElementMap(entity,orient,map,sign);
                }
                break;
                case eEdgeInterior:
                {
                    v_GetEdgeInteriorMap(entity,orient,map,sign);
                }
                break;
                case eFaceInterior:
                {
                    v_GetFaceInteriorMap(entity,orient,map,sign);
                }
                break;
                case eBoundary:
                {
                    ASSERTL0(false,"Boundary Index Map not implemented yet.");
                }
                break;
                case eVertex:
                {
                    ASSERTL0(false,"Vertex Index Map not implemented yet.");
                }
                break;
                default:
                {
                    ASSERTL0(false,"The Index Map you are requiring is not between the possible options.");
                }
            }

            returnval = MemoryManager<IndexMapValues>::AllocateSharedPtr(map.size());

            for(int i = 0; i < map.size(); i++)
            {
                (*returnval)[i].index =  map[i];
                (*returnval)[i].sign  =  sign[i];
            }

            return returnval;
        }

        DNekMatSharedPtr StdExpansion::CreateGeneralMatrix(const StdMatrixKey &mkey)
        {
            int     i;
            DNekMatSharedPtr  returnval;

            switch(mkey.GetMatrixType())
            {
            case eInvMass:
                {
                    StdMatrixKey masskey(eMass,mkey.GetShapeType(),*this,NullConstFactorMap,NullVarCoeffMap,mkey.GetNodalPointsType());
                    DNekMatSharedPtr mmat = GetStdMatrix(masskey);

                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(*mmat); //Populate standard mass matrix.
                    returnval->Invert();
                }
                break;
            case eInvNBasisTrans:
                {
                    StdMatrixKey tmpkey(eNBasisTrans,mkey.GetShapeType(),*this,NullConstFactorMap,NullVarCoeffMap,mkey.GetNodalPointsType());
                    DNekMatSharedPtr tmpmat = GetStdMatrix(tmpkey);
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(*tmpmat); //Populate  matrix.
                    returnval->Invert();
                }
                break;
            case eBwdTrans:
                {
                    int nq = GetTotPoints();
                    Array<OneD, NekDouble> tmpin(m_ncoeffs);
                    Array<OneD, NekDouble> tmpout(nq);

                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(nq,m_ncoeffs);

                    for(int i=0; i<m_ncoeffs; ++i)
                    {
                        Vmath::Zero(m_ncoeffs, tmpin, 1);
                        tmpin[i] = 1.0;

                        BwdTrans_SumFac(tmpin,tmpout);

                        Vmath::Vcopy(nq,tmpout.get(),1,
                                     returnval->GetRawPtr()+i*nq,1);
                    }
                }
                break;
            case eIProductWRTBase:
                {
                    int nq = GetTotPoints();
                    Array<OneD, NekDouble> tmpin(nq);
                    Array<OneD, NekDouble> tmpout(m_ncoeffs);

                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs,nq);

                    for(i=0; i < nq; ++i)
                    {
                        Vmath::Zero(nq, tmpin, 1);
                        tmpin[i] = 1.0;

                        IProductWRTBase_SumFac(tmpin,tmpout);

                        Vmath::Vcopy(m_ncoeffs,tmpout.get(),1,
                                     returnval->GetRawPtr()+i*m_ncoeffs,1);
                    }
                }
                break;
            case eIProductWRTDerivBase0:
                {
                    int nq = GetTotPoints();
                    Array<OneD, NekDouble> tmpin(nq);
                    Array<OneD, NekDouble> tmpout(m_ncoeffs);

                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(nq,m_ncoeffs);

                    for(i=0; i < nq; ++i)
                    {
                        Vmath::Zero(nq, tmpin, 1);
                        tmpin[i] = 1.0;

                        IProductWRTDerivBase_SumFac(0,tmpin,tmpout);

                        Vmath::Vcopy(m_ncoeffs,tmpout.get(),1,
                                     returnval->GetRawPtr()+i*m_ncoeffs,1);
                    }
                }
                break;
            case eIProductWRTDerivBase1:
                {
                    int nq = GetTotPoints();
                    Array<OneD, NekDouble> tmpin(nq);
                    Array<OneD, NekDouble> tmpout(m_ncoeffs);

                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(nq,m_ncoeffs);

                    for(i=0; i < nq; ++i)
                    {
                        Vmath::Zero(nq, tmpin, 1);
                        tmpin[i] = 1.0;

                        IProductWRTDerivBase_SumFac(1,tmpin,tmpout);

                        Vmath::Vcopy(m_ncoeffs,tmpout.get(),1,
                                     returnval->GetRawPtr()+i*m_ncoeffs,1);
                    }
                }
                break;
            case eIProductWRTDerivBase2:
                {
                    int nq = GetTotPoints();
                    Array<OneD, NekDouble> tmpin(nq);
                    Array<OneD, NekDouble> tmpout(m_ncoeffs);

                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(nq,m_ncoeffs);

                    for(i=0; i < nq; ++i)
                    {
                        Vmath::Zero(nq, tmpin, 1);
                        tmpin[i] = 1.0;

                        IProductWRTDerivBase_SumFac(2,tmpin,tmpout);

                        Vmath::Vcopy(m_ncoeffs,tmpout.get(),1,
                                     returnval->GetRawPtr()+i*m_ncoeffs,1);
                    }
                }
                break;
            case eEquiSpacedToCoeffs:
                {
                    // check to see if equispaced basis
                    int nummodes = m_base[0]->GetNumModes();
                    bool equispaced = true;
                    for(int i = 1; i < m_base.size(); ++i)
                    {
                        if(m_base[i]->GetNumModes() != nummodes)
                        {
                            equispaced = false;
                        }
                    }

                    ASSERTL0(equispaced,
                             "Currently need to have same num modes in all "
                             "directionmodes to use EquiSpacedToCoeff method");

                    int ntot = GetTotPoints();
                    Array<OneD, NekDouble>               qmode(ntot);
                    Array<OneD, NekDouble>               emode(m_ncoeffs);

                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(
                                                        m_ncoeffs,m_ncoeffs);
                    for(int i = 0; i < m_ncoeffs; ++i)
                    {
                        // Get mode at quadrature points
                        FillMode(i,qmode);

                        // interpolate to equi spaced
                        PhysInterpToSimplexEquiSpaced(qmode,emode,nummodes);

                        // fill matrix
                        Vmath::Vcopy(m_ncoeffs, &emode[0], 1,
                                     returnval->GetRawPtr() + i*m_ncoeffs, 1);
                    }
                    // invert matrix
                    returnval->Invert();

                }
                break;
            case eMass:
            case eHelmholtz:
            case eLaplacian:
            case eLaplacian00:
            case eLaplacian01:
            case eLaplacian02:
            case eLaplacian11:
            case eLaplacian12:
            case eLaplacian22:
            case eWeakDeriv0:
            case eWeakDeriv1:
            case eWeakDeriv2:
            case eWeakDirectionalDeriv:
            case eMassLevelCurvature:
            case eLinearAdvectionReaction:
            case eLinearAdvectionDiffusionReaction:
                {
                    Array<OneD, NekDouble> tmp(m_ncoeffs);
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs,m_ncoeffs);
                    DNekMat &Mat = *returnval;

                    for(i=0; i < m_ncoeffs; ++i)
                    {
                        Vmath::Zero(m_ncoeffs, tmp, 1);
                        tmp[i] = 1.0;

                        GeneralMatrixOp_MatFree(tmp,tmp,mkey);

                        Vmath::Vcopy(m_ncoeffs,&tmp[0],1,
                                     &(Mat.GetPtr())[0]+i*m_ncoeffs,1);
                    }
                }
                break;
            default:
                {
                    NEKERROR(ErrorUtil::efatal, "This type of matrix can not be created using a general approach");
                }
                break;
            }

            return returnval;
        }

        void StdExpansion::GeneralMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                           Array<OneD,NekDouble> &outarray,
                                           const StdMatrixKey &mkey)
        {
            switch(mkey.GetMatrixType())
            {
            case eMass:
                MassMatrixOp(inarray,outarray,mkey);
                break;
            case eWeakDeriv0:
                WeakDerivMatrixOp(0,inarray,outarray,mkey);
                break;
            case eWeakDeriv1:
                WeakDerivMatrixOp(1,inarray,outarray,mkey);
                break;
            case eWeakDeriv2:
                WeakDerivMatrixOp(2,inarray,outarray,mkey);
                break;
            case eWeakDirectionalDeriv:
                WeakDirectionalDerivMatrixOp(inarray,outarray,mkey);
                break;
            case eMassLevelCurvature:
                MassLevelCurvatureMatrixOp(inarray,outarray,mkey);
                break;
            case eLinearAdvectionReaction:
                LinearAdvectionDiffusionReactionMatrixOp(inarray,outarray,mkey,false);
                break;
            case eLinearAdvectionDiffusionReaction:
                LinearAdvectionDiffusionReactionMatrixOp(inarray,outarray,mkey);
                break;
            case eLaplacian:
                LaplacianMatrixOp(inarray,outarray,mkey);
                break;
            case eLaplacian00:
                LaplacianMatrixOp(0,0,inarray,outarray,mkey);
                break;
            case eLaplacian01:
                LaplacianMatrixOp(0,1,inarray,outarray,mkey);
                break;
            case eLaplacian02:
                LaplacianMatrixOp(0,2,inarray,outarray,mkey);
                break;
            case eLaplacian10:
                LaplacianMatrixOp(1,0,inarray,outarray,mkey);
                break;
            case eLaplacian11:
                LaplacianMatrixOp(1,1,inarray,outarray,mkey);
                break;
            case eLaplacian12:
                LaplacianMatrixOp(1,2,inarray,outarray,mkey);
                break;
            case eLaplacian20:
                LaplacianMatrixOp(2,0,inarray,outarray,mkey);
                break;
            case eLaplacian21:
                LaplacianMatrixOp(2,1,inarray,outarray,mkey);
                break;
            case eLaplacian22:
                LaplacianMatrixOp(2,2,inarray,outarray,mkey);
                break;
            case eHelmholtz:
                HelmholtzMatrixOp(inarray,outarray,mkey);
                break;
            default:
                NEKERROR(ErrorUtil::efatal, "This matrix does not have an operator");
                break;
            }
        }

        void StdExpansion::GeneralMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                   Array<OneD,NekDouble> &outarray,
                                                   const StdMatrixKey &mkey)
        {
            switch(mkey.GetMatrixType())
            {
            case eMass:
                MassMatrixOp_MatFree(inarray,outarray,mkey);
                break;
            case eWeakDeriv0:
                WeakDerivMatrixOp_MatFree(0,inarray,outarray,mkey);
                break;
            case eWeakDeriv1:
                WeakDerivMatrixOp_MatFree(1,inarray,outarray,mkey);
                break;
            case eWeakDeriv2:
                WeakDerivMatrixOp_MatFree(2,inarray,outarray,mkey);
                break;
            case eWeakDirectionalDeriv:
                WeakDirectionalDerivMatrixOp_MatFree(inarray,outarray,mkey);
                break;
            case eMassLevelCurvature:
                MassLevelCurvatureMatrixOp_MatFree(inarray,outarray,mkey);
                break;
            case eLinearAdvectionReaction:
                LinearAdvectionDiffusionReactionMatrixOp_MatFree(inarray,outarray,mkey,false);
                break;
            case eLinearAdvectionDiffusionReaction:
                LinearAdvectionDiffusionReactionMatrixOp_MatFree(inarray,outarray,mkey);
                break;
            case eLaplacian:
                LaplacianMatrixOp_MatFree(inarray,outarray,mkey);
                break;
            case eLaplacian00:
                LaplacianMatrixOp_MatFree(0,0,inarray,outarray,mkey);
                break;
            case eLaplacian01:
                LaplacianMatrixOp_MatFree(0,1,inarray,outarray,mkey);
                break;
            case eLaplacian02:
                LaplacianMatrixOp_MatFree(0,2,inarray,outarray,mkey);
                break;
            case eLaplacian10:
                LaplacianMatrixOp_MatFree(1,0,inarray,outarray,mkey);
                break;
            case eLaplacian11:
                LaplacianMatrixOp_MatFree(1,1,inarray,outarray,mkey);
                break;
            case eLaplacian12:
                LaplacianMatrixOp_MatFree(1,2,inarray,outarray,mkey);
                break;
            case eLaplacian20:
                LaplacianMatrixOp_MatFree(2,0,inarray,outarray,mkey);
                break;
            case eLaplacian21:
                LaplacianMatrixOp_MatFree(2,1,inarray,outarray,mkey);
                break;
            case eLaplacian22:
                LaplacianMatrixOp_MatFree(2,2,inarray,outarray,mkey);
                break;
            case eHelmholtz:
                HelmholtzMatrixOp_MatFree(inarray,outarray,mkey);
                break;
            default:
                NEKERROR(ErrorUtil::efatal, "This matrix does not have an operator");
                break;
            }
        }

        void StdExpansion::MassMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                Array<OneD,NekDouble> &outarray,
                                                const StdMatrixKey &mkey)
        {
            int nq = GetTotPoints();
            Array<OneD, NekDouble> tmp(nq);

            v_BwdTrans(inarray,tmp);

            if(mkey.HasVarCoeff(eVarCoeffMass))
            {
                Vmath::Vmul(nq, mkey.GetVarCoeff(eVarCoeffMass), 1, tmp, 1, tmp, 1);
            }

            v_IProductWRTBase(tmp, outarray);
        }

        void StdExpansion::LaplacianMatrixOp_MatFree(const int k1, const int k2,
                                                     const Array<OneD, const NekDouble> &inarray,
                                                     Array<OneD,NekDouble> &outarray,
                                                     const StdMatrixKey &mkey)
        {
            ASSERTL1(k1 >= 0 && k1 < GetCoordim(),"invalid first  argument");
            ASSERTL1(k2 >= 0 && k2 < GetCoordim(),"invalid second argument");

            int nq = GetTotPoints();
            Array<OneD, NekDouble> tmp(nq);
            Array<OneD, NekDouble> dtmp(nq);
            VarCoeffType varcoefftypes[3][3]
                = { {eVarCoeffD00, eVarCoeffD01, eVarCoeffD02},
                    {eVarCoeffD01, eVarCoeffD11, eVarCoeffD12},
                    {eVarCoeffD02, eVarCoeffD12, eVarCoeffD22}
            };

            v_BwdTrans(inarray,tmp);
            v_PhysDeriv(k2,tmp,dtmp);
            if (mkey.GetNVarCoeff()&&
                (!mkey.ConstFactorExists(eFactorSVVDiffCoeff)))
            {
                if (k1 == k2)
                {
                    // By default, k1 == k2 has \sigma = 1 (diagonal entries)
                    if(mkey.HasVarCoeff(varcoefftypes[k1][k1]))
                    {
                        Vmath::Vmul(nq, mkey.GetVarCoeff(varcoefftypes[k1][k1]), 1, dtmp, 1, dtmp, 1);
                    }
                    v_IProductWRTDerivBase(k1, dtmp, outarray);
                }
                else
                {
                    // By default, k1 != k2 has \sigma = 0 (off-diagonal entries)
                    if(mkey.HasVarCoeff(varcoefftypes[k1][k2]))
                    {
                        Vmath::Vmul(nq, mkey.GetVarCoeff(varcoefftypes[k1][k2]), 1, dtmp, 1, dtmp, 1);
                        v_IProductWRTDerivBase(k1, dtmp, outarray);
                    }
                    else
                    {
                        Vmath::Zero(GetNcoeffs(), outarray, 1);
                    }
                }
            }
            else
            {
                // Multiply by svv tensor
                if(mkey.ConstFactorExists(eFactorSVVDiffCoeff))
                {
                    Vmath::Vcopy(nq, dtmp, 1, tmp, 1);
                    SVVLaplacianFilter(dtmp,mkey);
                    Vmath::Vadd(nq, tmp, 1, dtmp, 1, dtmp, 1);
                }
                v_IProductWRTDerivBase(k1, dtmp, outarray);
            }
        }

        void StdExpansion::LaplacianMatrixOp_MatFree_GenericImpl(const Array<OneD, const NekDouble> &inarray,
                                                                 Array<OneD,NekDouble> &outarray,
                                                                 const StdMatrixKey &mkey)
        {
            const int dim = GetCoordim();

            int i,j;

            Array<OneD,NekDouble> store(m_ncoeffs);
            Array<OneD,NekDouble> store2(m_ncoeffs,0.0);

            if(mkey.GetNVarCoeff() == 0||mkey.ConstFactorExists(eFactorSVVDiffCoeff))
            {
                // just call diagonal matrix form of laplcian operator
                for(i = 0; i < dim; ++i)
                {
                    LaplacianMatrixOp(i,i,inarray,store,mkey);
                    Vmath::Vadd(m_ncoeffs, store, 1, store2, 1, store2, 1);
                }
            }
            else
            {
                const MatrixType mtype[3][3]
                    = {{eLaplacian00,eLaplacian01,eLaplacian02},
                       {eLaplacian01,eLaplacian11,eLaplacian12},
                       {eLaplacian02,eLaplacian12,eLaplacian22}};
                StdMatrixKeySharedPtr mkeyij;

                for(i = 0; i < dim; i++)
                {
                    for(j = 0; j < dim; j++)
                    {
                        mkeyij = MemoryManager<StdMatrixKey>::AllocateSharedPtr(mkey,mtype[i][j]);
                        LaplacianMatrixOp(i,j,inarray,store,*mkeyij);
                        Vmath::Vadd(m_ncoeffs, store, 1, store2, 1, store2, 1);
                    }
                }
            }

            Vmath::Vcopy(m_ncoeffs,store2.get(),1,outarray.get(),1);
        }

        void StdExpansion::WeakDerivMatrixOp_MatFree(const int k1,
                                                     const Array<OneD, const NekDouble> &inarray,
                                                     Array<OneD,NekDouble> &outarray,
                                                     const StdMatrixKey &mkey)
        {
            Array<OneD, NekDouble> tmp(GetTotPoints());
            int nq = GetTotPoints();

            v_BwdTrans(inarray,tmp);
            v_PhysDeriv(k1,tmp,tmp);

            VarCoeffType keys[] = {eVarCoeffD00, eVarCoeffD11, eVarCoeffD22};
            if(mkey.HasVarCoeff(keys[k1]))
            {
                Vmath::Vmul(nq, &(mkey.GetVarCoeff(keys[k1]))[0], 1, &tmp[0], 1, &tmp[0], 1);
            }

            v_IProductWRTBase(tmp, outarray);
        }

        void StdExpansion::WeakDirectionalDerivMatrixOp_MatFree(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,NekDouble> &outarray,
            const StdMatrixKey &mkey)
        {
            int nq = GetTotPoints();

            Array<OneD, NekDouble> tmp(nq), Dtmp(nq);
            Array<OneD, NekDouble> Mtmp(nq), Mout(m_ncoeffs);

            v_BwdTrans(inarray,tmp);
            v_PhysDirectionalDeriv(tmp, mkey.GetVarCoeff(eVarCoeffMF), Dtmp);

            v_IProductWRTBase(Dtmp, outarray);

            // Compte M_{div tv}
            Vmath::Vmul(nq, &(mkey.GetVarCoeff(eVarCoeffMFDiv))[0], 1,
                            &tmp[0],                                1,
                            &Mtmp[0],                               1);

            v_IProductWRTBase(Mtmp, Mout);

            // Add D_tv + M_{div tv}
            Vmath::Vadd(m_ncoeffs, &Mout[0],     1,
                                   &outarray[0], 1,
                                   &outarray[0], 1);
        }

        void StdExpansion::MassLevelCurvatureMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                              Array<OneD,NekDouble> &outarray,
                                                              const StdMatrixKey &mkey)
        {
            boost::ignore_unused(inarray, outarray, mkey);
            ///@todo fix this
            //          int nqtot = GetTotPoints();
            //          int matrixid = mkey.GetMatrixID();
            //
            //          NekDouble checkweight=0.0;
            //          Array<OneD, NekDouble> tmp(nqtot), tan(nqtot), dtan0(nqtot), dtan1(nqtot), weight(nqtot,0.0);
            //
            //          int gmatnumber = (mkey.GetVariableCoefficient(1)).size();
            //
            //          v_BwdTrans(inarray,tmp);
            //
            //          // weight = \grad \cdot tanvec
            //          for(int k = 0; k < GetCoordim(); ++k)
            //          {
            //              Vmath::Vcopy(nqtot, &(mkey.GetVariableCoefficient(0))[k*nqtot],
            //                              1, &tan[0], 1);
            //
            //              // For Regular mesh ...
            //              if(gmatnumber==1)
            //              {
            //                  // D_{/xi} and D_{/eta}
            //                  v_PhysDeriv(0,tan,dtan0);
            //                  v_PhysDeriv(1,tan,dtan1);
            //
            //                  // d v / d x_i = (d \xi / d x_i)*( d v / d \xi ) + (d \eta / d x_i)*( d v / d \eta )
            //                  Vmath::Svtvp(nqtot,(mkey.GetVariableCoefficient(2*k+1))[0],&dtan0[0],1,&weight[0],1,&weight[0],1);
            //                  Vmath::Svtvp(nqtot,(mkey.GetVariableCoefficient(2*k+2))[0],&dtan1[0],1,&weight[0],1,&weight[0],1);
            //              }
            //
            //              // For Curved mesh ...
            //              else if(gmatnumber==nqtot)
            //              {
            //                  // D_{x} and D_{y}
            //                  v_PhysDeriv(k,tan,dtan0);
            //                  Vmath::Vadd(nqtot,&dtan0[0],1,&weight[0],1,&weight[0],1);
            //              }
            //
            //              else
            //              {
            //                  ASSERTL1( ((gmatnumber=1) || (gmatnumber==nqtot) ), "Gmat is not in a right size");
            //              }
            //          }
            //
            //          Vmath::Vmul(nqtot, &weight[0], 1, &tmp[0], 1, &tmp[0], 1);
            //          v_IProductWRTBase(tmp, outarray);
        }

        void StdExpansion::LinearAdvectionDiffusionReactionMatrixOp_MatFree( const Array<OneD, const NekDouble> &inarray,
                                                                             Array<OneD,NekDouble> &outarray,
                                                                             const StdMatrixKey &mkey,
                                                                             bool addDiffusionTerm)
        {

            int i;
            int ndir = mkey.GetNVarCoeff(); // assume num.r consts corresponds to directions
            ASSERTL0(ndir,"Must define at least one advection velocity");

            NekDouble   lambda = mkey.GetConstFactor(eFactorLambda);
            int         totpts = GetTotPoints();
            Array<OneD, NekDouble> tmp(3*totpts);
            Array<OneD, NekDouble> tmp_deriv = tmp + totpts;
            Array<OneD, NekDouble> tmp_adv   = tmp_deriv + totpts;


            ASSERTL1(ndir <= GetCoordim(),"Number of constants is larger than coordinate dimensions");

            v_BwdTrans(inarray,tmp);

            VarCoeffType varcoefftypes[] = {eVarCoeffVelX, eVarCoeffVelY, eVarCoeffVelZ};

            //calculate u dx + v dy + ..
            Vmath::Zero(totpts,tmp_adv,1);
            for(i = 0; i < ndir; ++i)
            {
                v_PhysDeriv(i,tmp,tmp_deriv);
                Vmath::Vvtvp(totpts,mkey.GetVarCoeff(varcoefftypes[i]),1,tmp_deriv,1,tmp_adv,1,tmp_adv,1);
            }

            if(lambda) // add -lambda*u
            {
                Vmath::Svtvp(totpts,-lambda,tmp,1,tmp_adv,1,tmp_adv,1);
            }


            if(addDiffusionTerm)
            {
                Array<OneD, NekDouble> lap(m_ncoeffs);
                StdMatrixKey mkeylap(eLaplacian,DetShapeType(),*this,
                                     mkey.GetConstFactors(),
                                     mkey.GetVarCoeffs(),
                                     mkey.GetNodalPointsType());
                LaplacianMatrixOp(inarray,lap,mkeylap);

                v_IProductWRTBase(tmp_adv, outarray);
                // Lap v - u.grad v + lambda*u
                // => (grad u, grad v) + u.grad v - lambda*u
                Vmath::Vadd(m_ncoeffs,lap,1,outarray,1,outarray,1);
            }
            else
            {
                v_IProductWRTBase(tmp_adv, outarray);
            }

        }


        void StdExpansion::HelmholtzMatrixOp_MatFree_GenericImpl(const Array<OneD, const NekDouble> &inarray,
                                                                 Array<OneD,NekDouble> &outarray,
                                                                 const StdMatrixKey &mkey)
        {
            NekDouble lambda = mkey.GetConstFactor(eFactorLambda);
            Array<OneD,NekDouble> tmp(m_ncoeffs);
            StdMatrixKey mkeymass(eMass,DetShapeType(),*this);
            StdMatrixKey mkeylap(eLaplacian,DetShapeType(),*this,
                                 mkey.GetConstFactors(),
                                 mkey.GetVarCoeffs(),
                                 mkey.GetNodalPointsType());

            MassMatrixOp(inarray,tmp,mkeymass);
            LaplacianMatrixOp(inarray,outarray,mkeylap);

            Blas::Daxpy(m_ncoeffs, lambda, tmp, 1, outarray, 1);
        }

        void StdExpansion::BwdTrans_MatOp(const Array<OneD, const NekDouble>& inarray,
                                          Array<OneD, NekDouble> &outarray)
        {
            int nq = GetTotPoints();
            StdMatrixKey      bwdtransmatkey(eBwdTrans,DetShapeType(),*this);
            DNekMatSharedPtr  bwdtransmat = GetStdMatrix(bwdtransmatkey);

            Blas::Dgemv('N',nq,m_ncoeffs,1.0,bwdtransmat->GetPtr().get(),
                        nq, inarray.get(), 1, 0.0, outarray.get(), 1);
        }

        // VIRTUAL INLINE FUNCTIONS FROM HEADER FILE
        void StdExpansion::SetUpPhysNormals(const int edge)
        {
            v_SetUpPhysNormals(edge);
        }

        NekDouble StdExpansion::StdPhysEvaluate(const Array<OneD, const NekDouble> &Lcoord,
                                                const Array<OneD, const NekDouble> &physvals)
        {
            return v_StdPhysEvaluate(Lcoord,physvals);
        }

        int StdExpansion::v_GetElmtId(void)
        {
            return m_elmt_id;
        }

        const Array<OneD, const NekDouble>& StdExpansion::v_GetPhysNormals(void)
        {
            NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
            return NullNekDouble1DArray;
        }


        void StdExpansion::v_SetPhysNormals(Array<OneD, const NekDouble> &normal)
        {
            boost::ignore_unused(normal);
            NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
        }

        void StdExpansion::v_SetUpPhysNormals(const int edge)
        {
            boost::ignore_unused(edge);
            NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
        }

        int StdExpansion::v_CalcNumberOfCoefficients(const std::vector<unsigned int>  &nummodes, int &modes_offset)
        {
            boost::ignore_unused(nummodes, modes_offset);
            NEKERROR(ErrorUtil::efatal, "This function is not defined for this class");
            return 0;
        }

        void StdExpansion::v_NormVectorIProductWRTBase(const Array<OneD, const NekDouble> &Fx, Array< OneD, NekDouble> &outarray)
        {
            boost::ignore_unused(Fx, outarray);
            NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
        }

        void StdExpansion::v_NormVectorIProductWRTBase(const Array<OneD, const NekDouble> &Fx, const Array<OneD, const NekDouble> &Fy, Array< OneD, NekDouble> &outarray)
        {
            boost::ignore_unused(Fx, Fy, outarray);
            NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
        }

        void StdExpansion::v_NormVectorIProductWRTBase(const Array<OneD, const NekDouble> &Fx,
                                                       const Array<OneD, const NekDouble> &Fy,
                                                       const Array<OneD, const NekDouble> &Fz,
                                                       Array< OneD, NekDouble> &outarray)
        {
            boost::ignore_unused(Fx, Fy, Fz, outarray);
            NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
        }

        void StdExpansion::v_NormVectorIProductWRTBase(const Array<OneD, const Array<OneD, NekDouble> > &Fvec, Array< OneD, NekDouble> &outarray)
        {
            boost::ignore_unused(Fvec, outarray);
            NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
        }


        DNekScalBlkMatSharedPtr StdExpansion::v_GetLocStaticCondMatrix(const LocalRegions::MatrixKey &mkey)
        {
            boost::ignore_unused(mkey);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
            return NullDNekScalBlkMatSharedPtr;
        }

        void StdExpansion::v_DropLocStaticCondMatrix(const LocalRegions::MatrixKey &mkey)
        {
            boost::ignore_unused(mkey);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
        }

        StdRegions::Orientation StdExpansion::v_GetForient(int face)

        {
            boost::ignore_unused(face);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for three-dimensional  LocalRegions");
            return eDir1FwdDir1_Dir2FwdDir2;
        }

        StdRegions::Orientation StdExpansion::v_GetEorient(int edge)
        {
            boost::ignore_unused(edge);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for two-dimensional  LocalRegions");
            return eForwards;
        }

        void StdExpansion::v_SetCoeffsToOrientation(StdRegions::Orientation dir,
                                                    Array<OneD, const NekDouble> &inarray,
                                                    Array<OneD, NekDouble> &outarray)
        {
            boost::ignore_unused(dir, inarray, outarray);
            NEKERROR(ErrorUtil::efatal, "This function is not defined for this shape");
        }

        void StdExpansion::v_SetCoeffsToOrientation(
            Array<OneD, NekDouble> &coeffs,
            StdRegions::Orientation dir)
        {
            boost::ignore_unused(coeffs, dir);
            NEKERROR(ErrorUtil::efatal, "This function is not defined for this shape");
        }


        NekDouble StdExpansion::v_StdPhysEvaluate(const Array<OneD, const NekDouble> &Lcoord,
                                                  const Array<OneD, const NekDouble> &physvals)

        {
            boost::ignore_unused(Lcoord, physvals);
            NEKERROR(ErrorUtil::efatal, "This function is not defined for this shape");
            return 0;
        }

        void StdExpansion::v_LocCoordToLocCollapsed(const Array<OneD, const NekDouble>& xi,Array<OneD, NekDouble>& eta)
        {
            boost::ignore_unused(xi, eta);
            NEKERROR(ErrorUtil::efatal, "This function is not defined for this shape");
        }

        int StdExpansion::v_GetNedges() const
        {
            ASSERTL0(false, "This function is needs defining for this shape");
            return 0;
        }

        int StdExpansion::v_GetNfaces() const
        {
            ASSERTL0(false, "This function is needs defining for this shape");
            return 0;
        }

        int StdExpansion::v_NumBndryCoeffs() const
        {
            ASSERTL0(false, "This function is needs defining for this shape");
            return 0;
        }

        int StdExpansion::v_NumDGBndryCoeffs() const
        {
            ASSERTL0(false, "This function is needs defining for this shape");
            return 0;
        }

        int StdExpansion::v_GetEdgeNcoeffs(const int i) const
        {
            boost::ignore_unused(i);
            ASSERTL0(false, "This function is not valid or not defined");
            return 0;
        }

        int StdExpansion::v_GetTotalEdgeIntNcoeffs() const
        {
            ASSERTL0(false, "This function is not valid or not defined");
            return 0;
        }

        int StdExpansion::v_GetEdgeNumPoints(const int i) const
        {
            boost::ignore_unused(i);
            ASSERTL0(false, "This function is not valid or not defined");
            return 0;
        }

        int StdExpansion::v_DetCartesianDirOfEdge(const int edge)
        {
            boost::ignore_unused(edge);
            ASSERTL0(false, "This function is not valid or not defined");
            return 0;
        }

        const LibUtilities::BasisKey StdExpansion::v_DetEdgeBasisKey(const int i) const
        {
            boost::ignore_unused(i);
            ASSERTL0(false, "This function is not valid or not defined");
            return LibUtilities::NullBasisKey;
        }

        const LibUtilities::BasisKey StdExpansion::v_DetFaceBasisKey(const int i, const int k) const
        {
            boost::ignore_unused(i, k);
            ASSERTL0(false, "This function is not valid or not defined");
            return LibUtilities::NullBasisKey;
        }

        int StdExpansion::v_GetFaceNumPoints(const int i) const
        {
            boost::ignore_unused(i);
            ASSERTL0(false, "This function is not valid or not defined");
            return 0;
        }

        int StdExpansion::v_GetFaceNcoeffs(const int i) const
        {
            boost::ignore_unused(i);
            ASSERTL0(false, "This function is not valid or not defined");
            return 0;
        }

        int StdExpansion::v_GetFaceIntNcoeffs(const int i) const
        {
            boost::ignore_unused(i);
            ASSERTL0(false, "This function is not valid or not defined");
            return 0;
        }

        int StdExpansion::v_GetTotalFaceIntNcoeffs() const
        {
            ASSERTL0(false, "This function is not valid or not defined");
            return 0;
        }

        int StdExpansion::v_GetTraceNcoeffs(const int i) const
        {
            boost::ignore_unused(i);
            ASSERTL0(false, "This function is not valid or not defined");
            return 0;
        }


        LibUtilities::PointsKey StdExpansion::v_GetFacePointsKey(const int i, const int j) const
        {
            boost::ignore_unused(i, j);
            ASSERTL0(false, "This function is not valid or not defined");
            return LibUtilities::NullPointsKey;
        }

        LibUtilities::BasisType StdExpansion::v_GetEdgeBasisType(const int i) const
        {
            boost::ignore_unused(i);
            ASSERTL0(false, "This function is not valid or not defined");

            return LibUtilities::eNoBasisType;
        }

        const LibUtilities::PointsKey StdExpansion::v_GetNodalPointsKey() const
        {
            ASSERTL0(false, "This function is not valid or not defined");

            return LibUtilities::NullPointsKey;
        }

        LibUtilities::ShapeType StdExpansion::v_DetShapeType() const
        {
            ASSERTL0(false, "This expansion does not have a shape type defined");
            return LibUtilities::eNoShapeType;
        }

        std::shared_ptr<StdExpansion>
        StdExpansion::v_GetStdExp(void) const
        {
            ASSERTL0(false,"This method is not defined for this expansion");
            StdExpansionSharedPtr returnval;
            return returnval;
        }

        std::shared_ptr<StdExpansion>
        StdExpansion::v_GetLinStdExp(void) const
        {
            ASSERTL0(false,"This method is not defined for this expansion");
            StdExpansionSharedPtr returnval;
            return returnval;
        }

        int StdExpansion::v_GetShapeDimension() const
        {
            ASSERTL0(false, "This function is not valid or not defined");
            return 0;
        }

        bool StdExpansion::v_IsBoundaryInteriorExpansion()
        {
            ASSERTL0(false,"This function has not been defined for this expansion");
            return false;
        }


        bool StdExpansion::v_IsNodalNonTensorialExp()
        {
            return false;
        }

        void  StdExpansion::v_IProductWRTDerivBase (const int dir,
                                                    const Array<OneD, const NekDouble>& inarray,
                                                    Array<OneD, NekDouble> &outarray)
        {
            boost::ignore_unused(dir, inarray, outarray);
            NEKERROR(ErrorUtil::efatal, "This method has not been defined");
        }


        /**
         *
         */
        void  StdExpansion::v_IProductWRTDirectionalDerivBase (
            const Array<OneD, const NekDouble>& direction,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD, NekDouble> &outarray)
        {
            boost::ignore_unused(direction, inarray, outarray);
            NEKERROR(ErrorUtil::efatal, "This method has not been defined");
        }


        /**
         *
         */
        void StdExpansion::v_FwdTrans_BndConstrained(const Array<OneD, const NekDouble>& inarray,
                                                     Array<OneD, NekDouble> &outarray)
        {
            boost::ignore_unused(inarray, outarray);
            NEKERROR(ErrorUtil::efatal, "This method has not been defined");
        }


        /**
         * @brief Integrates the specified function over the domain.
         * @see StdRegions#StdExpansion#Integral.
         */
        NekDouble StdExpansion::v_Integral(const Array<OneD, const NekDouble>& inarray )
        {
            boost::ignore_unused(inarray);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                     "local expansions");
            return 0;
        }


        void StdExpansion::v_AddRobinMassMatrix(const int edgeid, const Array<OneD, const NekDouble > &primCoeffs, DNekMatSharedPtr &inoutmat)
        {
            boost::ignore_unused(edgeid, primCoeffs, inoutmat);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                     "specific element types");
        }

        void StdExpansion::v_AddRobinEdgeContribution(const int edgeid,
                                        const Array<OneD, const NekDouble > &primCoeffs,
                                        const Array<OneD, NekDouble> &incoeffs,
                                        Array<OneD, NekDouble> &coeffs)
        {
            boost::ignore_unused(edgeid, primCoeffs, incoeffs, coeffs);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                     "specific element types");
        }

        /**
         * @brief Calculate the derivative of the physical points
         * @see StdRegions#StdExpansion#PhysDeriv
         */
        void StdExpansion::v_PhysDeriv (const Array<OneD, const NekDouble>& inarray,
                                        Array<OneD, NekDouble> &out_d1,
                                        Array<OneD, NekDouble> &out_d2,
                                        Array<OneD, NekDouble> &out_d3)
        {
            boost::ignore_unused(inarray, out_d1, out_d2, out_d3);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                     "local expansions");
        }

        void StdExpansion::v_PhysDeriv_s(const Array<OneD, const NekDouble>& inarray,
                                         Array<OneD, NekDouble> &out_ds)
        {
            boost::ignore_unused(inarray, out_ds);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                     "local expansions");
        }
        void StdExpansion::v_PhysDeriv_n(const Array<OneD, const NekDouble>& inarray,
                                         Array<OneD, NekDouble>& out_dn)
        {
            boost::ignore_unused(inarray, out_dn);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                     "local expansions");
        }

        /**
         * @brief Calculate the derivative of the physical points in a
         * given direction
         * @see StdRegions#StdExpansion#PhysDeriv
         */
        void StdExpansion::v_PhysDeriv(const int dir,
                                       const Array<OneD, const NekDouble>& inarray,
                                       Array<OneD, NekDouble> &out_d0)

        {
            boost::ignore_unused(dir, inarray, out_d0);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                     "specific element types");
        }

        /**
         * @brief Physical derivative along a direction vector.
         * @see StdRegions#StdExpansion#PhysDirectionalDeriv
         */
        void StdExpansion::v_PhysDirectionalDeriv(const Array<OneD, const NekDouble>& inarray,
                                                  const Array<OneD, const NekDouble>& direction,
                                                  Array<OneD, NekDouble> &outarray)
        {
            boost::ignore_unused(inarray, direction, outarray);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                     "specific element types");
        }

        void StdExpansion::v_StdPhysDeriv (const Array<OneD, const NekDouble>& inarray,
                                           Array<OneD, NekDouble> &out_d1,
                                           Array<OneD, NekDouble> &out_d2,
                                           Array<OneD, NekDouble> &out_d3)
        {
            boost::ignore_unused(inarray, out_d1, out_d2, out_d3);
            NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
        }

        void   StdExpansion::v_StdPhysDeriv (const int dir,
                                             const Array<OneD, const NekDouble>& inarray,
                                             Array<OneD, NekDouble> &outarray)
        {
            boost::ignore_unused(dir, inarray, outarray);
            NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
        }


        NekDouble StdExpansion::v_PhysEvaluate(const Array<OneD, const NekDouble>& coords, const Array<OneD, const NekDouble>& physvals)
        {
            boost::ignore_unused(coords, physvals);
            NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
            return 0;
        }


        NekDouble StdExpansion::v_PhysEvaluate(const Array<OneD, DNekMatSharedPtr > & I, const Array<OneD, const NekDouble>& physvals)
        {
            boost::ignore_unused(I, physvals);
            NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
            return 0;
        }


        void StdExpansion::v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
        {
            boost::ignore_unused(mode, outarray);
            NEKERROR(ErrorUtil::efatal, "This function has not "
                     "been defined for this shape");
        }

        DNekMatSharedPtr StdExpansion::v_GenMatrix(const StdMatrixKey &mkey)
        {
            boost::ignore_unused(mkey);
            NEKERROR(ErrorUtil::efatal, "This function has not "
                     "been defined for this element");
            DNekMatSharedPtr returnval;
            return returnval;
        }

        DNekMatSharedPtr StdExpansion::v_CreateStdMatrix(const StdMatrixKey &mkey)
        {
            boost::ignore_unused(mkey);
            NEKERROR(ErrorUtil::efatal, "This function has not "
                     "been defined for this element");
            DNekMatSharedPtr returnval;
            return returnval;
        }

        void StdExpansion::v_GetCoords(Array<OneD, NekDouble> &coords_0,
                                       Array<OneD, NekDouble> &coords_1,
                                       Array<OneD, NekDouble> &coords_2)
        {
            boost::ignore_unused(coords_0, coords_1, coords_2);
            NEKERROR(ErrorUtil::efatal, "Write coordinate definition method");
        }

        void StdExpansion::v_GetCoord(const Array<OneD, const NekDouble>& Lcoord,
                                      Array<OneD, NekDouble> &coord)
        {
            boost::ignore_unused(Lcoord, coord);
            NEKERROR(ErrorUtil::efatal, "Write coordinate definition method");
        }

            int StdExpansion::v_GetCoordim(void)
            {
                NEKERROR(ErrorUtil::efatal, "Write method");
                return -1;
            }

            void StdExpansion::v_GetBoundaryMap(Array<OneD, unsigned int>& outarray)
            {
                boost::ignore_unused(outarray);
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
            }

            void StdExpansion::v_GetInteriorMap(Array<OneD, unsigned int>& outarray)
            {
                boost::ignore_unused(outarray);
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
            }

            int StdExpansion::v_GetVertexMap(const int localVertexId,
                                         bool useCoeffPacking)
            {
                boost::ignore_unused(localVertexId, useCoeffPacking);
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
                return 0;
            }

            void StdExpansion::v_GetEdgeInteriorMap(const int eid, const Orientation edgeOrient,
                                              Array<OneD, unsigned int> &maparray,
                                              Array<OneD, int> &signarray)
            {
                boost::ignore_unused(eid, edgeOrient, maparray, signarray);
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
            }

            void StdExpansion::v_GetFaceNumModes(
                                              const int fid,
                                              const Orientation faceOrient,
                                              int &numModes0,
                                              int &numModes1)
            {
                boost::ignore_unused(fid, faceOrient, numModes0, numModes1);
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
            }

            void StdExpansion::v_GetFaceInteriorMap(const int fid, const Orientation faceOrient,
                                              Array<OneD, unsigned int> &maparray,
                                              Array<OneD, int> &signarray)
            {
                boost::ignore_unused(fid, faceOrient, maparray, signarray);
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
            }

            void StdExpansion::v_GetEdgeToElementMap(
                const int                  eid,
                const Orientation          edgeOrient,
                Array<OneD, unsigned int>& maparray,
                Array<OneD, int>&          signarray,
                int                        P)
            {
                boost::ignore_unused(eid, edgeOrient, maparray, signarray, P);
                NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
            }

            void StdExpansion::v_GetFaceToElementMap(const int fid, const Orientation faceOrient,
                                                     Array<OneD, unsigned int> &maparray,
                                                     Array<OneD, int> &signarray,
                                                     int nummodesA, int nummodesB)
            {
                boost::ignore_unused(fid, faceOrient, maparray, signarray,
                                     nummodesA, nummodesB);
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
            }

            void StdExpansion::v_GetEdgePhysVals(const int edge, const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
            {
                boost::ignore_unused(edge, inarray, outarray);
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape or library" );
            }

            void StdExpansion::v_GetEdgePhysVals(const int edge,  const std::shared_ptr<StdExpansion>  &EdgeExp, const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
            {
                boost::ignore_unused(edge, EdgeExp, inarray, outarray);
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape or library" );
            }

        void StdExpansion::v_GetTracePhysVals(const int edge,  const std::shared_ptr<StdExpansion>  &EdgeExp, const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray, StdRegions::Orientation  orient)
            {
                boost::ignore_unused(edge, EdgeExp, inarray, outarray, orient);
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape or library" );
            }

            void StdExpansion::v_GetVertexPhysVals(const int vertex, const Array<OneD, const NekDouble> &inarray, NekDouble &outarray)
            {
                boost::ignore_unused(vertex, inarray, outarray);
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape or library" );
            }

            void StdExpansion::v_AddVertexPhysVals(
                const int                 vertex,
                const NekDouble           &inarray,
                Array<OneD, NekDouble>   &outarray)
            {
                boost::ignore_unused(vertex, inarray, outarray);
                NEKERROR(ErrorUtil::efatal, 
                    "Method does not exist for this shape or library" );
            }

            void StdExpansion::v_GetEdgeInterpVals(const int edge,const Array<OneD, const NekDouble> &inarray,Array<OneD,NekDouble> &outarray)
            {
                boost::ignore_unused(edge, inarray, outarray);
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape or library" );
            }

            void StdExpansion::v_GetEdgeQFactors(
                    const int edge,
                    Array<OneD, NekDouble> &outarray)
            {
                boost::ignore_unused(edge, outarray);
                NEKERROR(ErrorUtil::efatal,
                     "Method does not exist for this shape or library");
            }

            void StdExpansion::v_GetFacePhysVals( const int                                face,
                                             const std::shared_ptr<StdExpansion>   &FaceExp,
                const Array<OneD, const NekDouble>      &inarray,
                      Array<OneD,       NekDouble>      &outarray,
                StdRegions::Orientation                  orient)
            {
                boost::ignore_unused(face, FaceExp, inarray, outarray, orient);
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape or library" );
            }

            void StdExpansion::v_GetEdgePhysMap(
                const int  edge,
                Array<OneD, int>   &outarray)
            {
                boost::ignore_unused(edge, outarray);
                NEKERROR(ErrorUtil::efatal,
                     "Method does not exist for this shape or library" );
            }

            void StdExpansion::v_GetFacePhysMap(const int  face,
                                                Array<OneD, int>   &outarray)
            {
                boost::ignore_unused(face, outarray);
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape or library" );
            }

            void StdExpansion::v_MultiplyByQuadratureMetric(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray)
            {
                boost::ignore_unused(inarray, outarray);
                v_MultiplyByStdQuadratureMetric(inarray,outarray);
            }

            void StdExpansion::v_DivideByQuadratureMetric(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray)
            {
                v_DivideByStdQuadratureMetric(inarray, outarray);
            }

            void StdExpansion::v_MultiplyByStdQuadratureMetric(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray)
            {
                boost::ignore_unused(inarray, outarray);
                NEKERROR(ErrorUtil::efatal,
                    "Method does not exist for this shape or library");
            }

            void StdExpansion::v_DivideByStdQuadratureMetric(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray)
            {
                boost::ignore_unused(inarray,outarray);
                NEKERROR(ErrorUtil::efatal,
                         "v_DivideByStdQuadratureMetric does not exist for this"
                         " shape or library");
            }

            void StdExpansion::v_BwdTrans_SumFac(const Array<OneD, const NekDouble>& inarray,
                                           Array<OneD, NekDouble> &outarray)
            {
                boost::ignore_unused(inarray, outarray);
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
            }

            void StdExpansion::v_IProductWRTBase_SumFac(const Array<OneD, const NekDouble>& inarray,
                                                        Array<OneD, NekDouble> &outarray,
                                                        bool multiplybyweights)
            {
                boost::ignore_unused(inarray, outarray, multiplybyweights);
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
            }

            /**
             *
             */
            void StdExpansion::v_IProductWRTDirectionalDerivBase_SumFac(
                const Array<OneD, const NekDouble>& direction,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD, NekDouble> &outarray)
            {
                boost::ignore_unused(direction, inarray, outarray);
                NEKERROR(ErrorUtil::efatal,
                         "Method does not exist for this shape" );
            }

            void StdExpansion::v_IProductWRTDerivBase_SumFac(const int dir,
                                                       const Array<OneD, const NekDouble>& inarray,
                                                       Array<OneD, NekDouble> &outarray)
            {
                boost::ignore_unused(dir, inarray, outarray);
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
            }

            void StdExpansion::v_MassMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                        Array<OneD,NekDouble> &outarray,
                                        const StdMatrixKey &mkey)
            {
                // If this function is not reimplemented on shape level, the function
                // below will be called
                MassMatrixOp_MatFree(inarray,outarray,mkey);
            }

            void StdExpansion::v_LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             const StdMatrixKey &mkey)
            {
                // If this function is not reimplemented on shape level, the function
                // below will be called
                LaplacianMatrixOp_MatFree(inarray,outarray,mkey);
            }


             void StdExpansion::v_SVVLaplacianFilter(Array<OneD,NekDouble> &array,
                                             const StdMatrixKey &mkey)
             {
                 boost::ignore_unused(array, mkey);
                 ASSERTL0(false, "This function is not defined in StdExpansion.");
             }

            void StdExpansion::v_ExponentialFilter(
                                          Array<OneD, NekDouble> &array,
                                    const NekDouble        alpha,
                                    const NekDouble        exponent,
                                    const NekDouble        cutoff)
             {
                 boost::ignore_unused(array, alpha, exponent, cutoff);
                 ASSERTL0(false, "This function is not defined in StdExpansion.");
             }

            void StdExpansion::v_ReduceOrderCoeffs(int numMin,
                                                   const Array<OneD, const NekDouble> &inarray,
                                                   Array<OneD, NekDouble> &outarray)
            {
                boost::ignore_unused(numMin, inarray, outarray);
                ASSERTL0(false, "This function is not defined in StdExpansion.");
            }

            void StdExpansion::v_LaplacianMatrixOp(const int k1, const int k2,
                                             const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             const StdMatrixKey &mkey)
            {
                // If this function is not reimplemented on shape level, the function
                // below will be called
                LaplacianMatrixOp_MatFree(k1,k2,inarray,outarray,mkey);
            }

            void StdExpansion::v_WeakDerivMatrixOp(const int i,
                                             const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             const StdMatrixKey &mkey)
            {
                // If this function is not reimplemented on shape level, the function
                // below will be called
                WeakDerivMatrixOp_MatFree(i,inarray,outarray,mkey);

            }

        void StdExpansion::v_WeakDirectionalDerivMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                                          Array<OneD,NekDouble> &outarray,
                                                          const StdMatrixKey &mkey)
        {
            // If this function is not reimplemented on shape level, the function
            // below will be called
            WeakDirectionalDerivMatrixOp_MatFree(inarray,outarray,mkey);

        }

        void StdExpansion::v_MassLevelCurvatureMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                                        Array<OneD,NekDouble> &outarray,
                                                        const StdMatrixKey &mkey)
        {
            // If this function is not reimplemented on shape level, the function
            // below will be called
            MassLevelCurvatureMatrixOp_MatFree(inarray,outarray,mkey);
        }

        void StdExpansion::v_LinearAdvectionDiffusionReactionMatrixOp(const Array<OneD,
                                                                      const NekDouble> &inarray,
                                                                      Array<OneD,NekDouble> &outarray,
                                                                      const StdMatrixKey &mkey, bool addDiffusionTerm)
        {
            // If this function is not reimplemented on shape level, the function
            // below will be called
            LinearAdvectionDiffusionReactionMatrixOp_MatFree(inarray,outarray,mkey,addDiffusionTerm);

        }

        void StdExpansion::v_HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                               Array<OneD,NekDouble> &outarray,
                                               const StdMatrixKey &mkey)
        {
            // If this function is not reimplemented on shape level, the function
            // below will be called
            HelmholtzMatrixOp_MatFree(inarray,outarray,mkey);
        }

        void StdExpansion::v_LaplacianMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                       Array<OneD,NekDouble> &outarray,
                                                       const StdMatrixKey &mkey)
        {
            // If this function is not reimplemented on shape level, the function
            // below will be called
            LaplacianMatrixOp_MatFree_GenericImpl(inarray,outarray,mkey);
        }

        void StdExpansion::v_LaplacianMatrixOp_MatFree_Kernel(
                            const Array<OneD, const NekDouble> &inarray,
                                  Array<OneD,       NekDouble> &outarray,
                                  Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(inarray, outarray, wsp);
            ASSERTL0(false, "Not implemented.");
        }

        void StdExpansion::v_HelmholtzMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                       Array<OneD,NekDouble> &outarray,
                                                       const StdMatrixKey &mkey)
        {
            // If this function is not reimplemented on shape level, the function
            // below will be called
            HelmholtzMatrixOp_MatFree_GenericImpl(inarray,outarray,mkey);
        }

        const NormalVector & StdExpansion::v_GetEdgeNormal(const int edge) const
        {
            boost::ignore_unused(edge);
            ASSERTL0(false, "Cannot get edge normals for this expansion.");
            static NormalVector result;
            return result;
        }

        void StdExpansion::v_ComputeEdgeNormal(const int edge)
        {
            boost::ignore_unused(edge);
            ASSERTL0(false, "Cannot compute edge normal for this expansion.");
        }

        void StdExpansion::v_ComputeFaceNormal(const int face)
        {
            boost::ignore_unused(face);
            ASSERTL0(false, "Cannot compute face normal for this expansion.");
        }

        void StdExpansion::v_ComputeVertexNormal(const int vertex)
        {
            boost::ignore_unused(vertex);
            ASSERTL0(false, "Cannot compute vertex normal for this expansion.");
        }

        const NormalVector & StdExpansion::v_GetFaceNormal(const int face) const
        {
            boost::ignore_unused(face);
            ASSERTL0(false, "Cannot get face normals for this expansion.");
            static NormalVector result;
            return result;
        }

        const NormalVector & StdExpansion::v_GetVertexNormal(const int vertex) const
        {
            boost::ignore_unused(vertex);
            ASSERTL0(false, "Cannot get vertex normals for this expansion.");
            static NormalVector result;
            return result;
        }

        const NormalVector & StdExpansion::v_GetSurfaceNormal(const int id) const
        {
            boost::ignore_unused(id);
            ASSERTL0(false, "Cannot get face normals for this expansion.");
            static NormalVector result;
            return result;
        }

        Array<OneD, unsigned int>
        StdExpansion::v_GetEdgeInverseBoundaryMap(int eid)
        {
            boost::ignore_unused(eid);
            ASSERTL0(false, "Not implemented.");
            Array<OneD, unsigned int> noinversemap(1);
            return noinversemap;
        }

        Array<OneD, unsigned int>
        StdExpansion::v_GetFaceInverseBoundaryMap(int fid,
                                                  StdRegions::Orientation faceOrient,
                                                  int P1,
                                                  int P2)
        {
            boost::ignore_unused(fid, faceOrient, P1, P2);
            ASSERTL0(false, "Not implemented.");
            Array<OneD, unsigned int> noinversemap(1);
            return noinversemap;
        }

        void StdExpansion::v_GetInverseBoundaryMaps(
                    Array<OneD, unsigned int> &vmap,
                    Array<OneD, Array<OneD, unsigned int> > &emap,
                    Array<OneD, Array<OneD, unsigned int> > &fmap)
        {
            boost::ignore_unused(vmap, emap, fmap);
            ASSERTL0(false, "Not implemented.");
        }

        DNekMatSharedPtr
        StdExpansion::v_BuildInverseTransformationMatrix(
            const DNekScalMatSharedPtr & m_transformationmatrix)
        {
            boost::ignore_unused(m_transformationmatrix);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
            return NullDNekMatSharedPtr;
        }

        void StdExpansion::PhysInterpToSimplexEquiSpaced(
            const Array<OneD, const NekDouble> &inarray,
            Array<OneD, NekDouble>       &outarray,
            int npset)
        {
            LibUtilities::ShapeType shape = DetShapeType();
            DNekMatSharedPtr  intmat;

            int nqtot = GetTotPoints();
            int np = 0;
            if(npset == -1) // use values from basis num points()
            {
                int nqbase;
                for(int i = 0; i < m_base.size(); ++i)
                {
                    nqbase = m_base[i]->GetNumPoints();
                    np     = std::max(np,nqbase);
                }

                StdMatrixKey Ikey(ePhysInterpToEquiSpaced, shape, *this);
                intmat = GetStdMatrix(Ikey);
            }
            else
            {
                np = npset;

                ConstFactorMap cmap;
                cmap[eFactorConst] = np;
                StdMatrixKey Ikey(ePhysInterpToEquiSpaced, shape, *this, cmap);
                intmat = GetStdMatrix(Ikey);

            }

            NekVector<NekDouble> in (nqtot,inarray,eWrapper);
            NekVector<NekDouble> out(LibUtilities::GetNumberOfCoefficients(shape,np,np,np),outarray,eWrapper);
            out = (*intmat) * in;
        }

        void StdExpansion::v_GetSimplexEquiSpacedConnectivity(
            Array<OneD, int> &conn,
            bool              standard)
        {
            boost::ignore_unused(conn, standard);
            ASSERTL0(false, "Not implemented.");
        }

        void StdExpansion::EquiSpacedToCoeffs(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD, NekDouble>       &outarray)
        {
            LibUtilities::ShapeType shape = DetShapeType();

            // inarray has to be consistent with NumModes definition
            // There is also a check in GetStdMatrix to see if all
            // modes are of the same size
            ConstFactorMap cmap;

            cmap[eFactorConst] = m_base[0]->GetNumModes();
            StdMatrixKey      Ikey(eEquiSpacedToCoeffs, shape, *this,cmap);
            DNekMatSharedPtr  intmat = GetStdMatrix(Ikey);

            NekVector<NekDouble> in (m_ncoeffs, inarray, eWrapper);
            NekVector<NekDouble> out(m_ncoeffs, outarray,eWrapper);
            out = (*intmat) * in;
        }
        
    }//end namespace
}//end namespace
