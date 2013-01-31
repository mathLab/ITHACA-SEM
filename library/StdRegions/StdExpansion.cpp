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
// Description: Definition of methods in class StdExpansion which is
// the base class to all expansion shapes
//
///////////////////////////////////////////////////////////////////////////////


#include <StdRegions/StdExpansion.h>
#include <LibUtilities/Foundations/ManagerAccess.h>  // for BasisManager, etc

namespace Nektar
{
    namespace StdRegions
    {

        /** define list of number of vertices corresponding to each ShapeType */
        const int g_shapenverts[SIZE_ExpansionType] = {0,2,3,4,4,5,6,8};

        /** define list of number of edges corresponding to each ShapeType */
        const int g_shapenedges[SIZE_ExpansionType] = {0,1,3,4,6,8,9,12};

        /** define list of number of faces corresponding to each ShapeType */
        const int g_shapenfaces[SIZE_ExpansionType] = {0,0,0,0,4,5,5,6};

        StdExpansion::StdExpansion(void):
            m_elmt_id(0),
            m_numbases(0),
            m_ncoeffs(0)
        {
        }

        StdExpansion::StdExpansion(const int numcoeffs, const int numbases,
            const LibUtilities::BasisKey &Ba,
            const LibUtilities::BasisKey &Bb,
            const LibUtilities::BasisKey &Bc):
            m_elmt_id(0),
            m_numbases(numbases),
            m_base(m_numbases),
            m_ncoeffs(numcoeffs),
            m_coeffs(m_ncoeffs,0.0),
            m_stdMatrixManager(
                    boost::bind(&StdExpansion::CreateStdMatrix, this, _1),
                    std::string("StdExpansionStdMatrix")),
            m_stdStaticCondMatrixManager(
                    boost::bind(&StdExpansion::CreateStdStaticCondMatrix, this, _1),
                    std::string("StdExpansionStdStaticCondMatrix")),
		    m_IndexMapManager(
					boost::bind(&StdExpansion::CreateIndexMap,this, _1),
					std::string("StdExpansionIndexMap"))
        {
            switch(m_numbases)
            {
            case 3:
                ASSERTL2(Bc!=LibUtilities::NullBasisKey,
                    "NULL Basis attempting to be used.");
                m_base[2] = LibUtilities::BasisManager()[Bc];

            case 2:
                ASSERTL2(Bb!=LibUtilities::NullBasisKey,
                    "NULL Basis attempting to be used.");

                m_base[1] = LibUtilities::BasisManager()[Bb];
            case 1:
                ASSERTL2(Ba!=LibUtilities::NullBasisKey,
                    "NULL Basis attempting to be used.");
                m_base[0] = LibUtilities::BasisManager()[Ba];
                break;
            default:
                ASSERTL0(false, "numbases incorrectly specified");
            };

            //allocate memory for phys
            m_phys = Array<OneD, NekDouble>(GetTotPoints());

        } //end constructor


        StdExpansion::StdExpansion(const StdExpansion &T):
            m_elmt_id(T.m_elmt_id),
            m_numbases(T.m_numbases),
            m_base(T.m_base),
            m_ncoeffs(T.m_ncoeffs),
            m_coeffs(m_ncoeffs),
            m_phys((T.m_phys).num_elements()),
            m_stdMatrixManager(T.m_stdMatrixManager),
            m_stdStaticCondMatrixManager(T.m_stdStaticCondMatrixManager),
		    m_IndexMapManager(T.m_IndexMapManager)
        {
            //CopyArray(T.m_base, m_base);
            CopyArray(T.m_coeffs, m_coeffs);
            CopyArray(T.m_phys, m_phys);
        }

        StdExpansion::~StdExpansion()
        {
        }

        NekDouble StdExpansion::Linf(const Array<OneD, const NekDouble>& sol)
        {
            NekDouble  val;
            int     ntot = GetTotPoints();
            Array<OneD, NekDouble>  wsp(ntot);

            Vmath::Vsub(ntot, sol, 1, m_phys, 1, wsp, 1);
            Vmath::Vabs(ntot, wsp, 1, wsp, 1);
            val = Vmath::Vamax(ntot, wsp, 1);

            return  val;
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
			
            returnval = MemoryManager<IndexMapValues>::AllocateSharedPtr(map.num_elements());
			
            for(int i = 0; i < map.num_elements(); i++)
            {
                (*returnval)[i].index =  map[i];
                (*returnval)[i].sign  =  sign[i];
            }
			
            return returnval;
        }

        NekDouble StdExpansion::Linf()
        {
            return Vmath::Vamax(GetTotPoints(), m_phys, 1);
        }

        NekDouble StdExpansion::L2(const Array<OneD, const NekDouble>& sol)
        {
            NekDouble  val;
            int     ntot = GetTotPoints();
            Array<OneD, NekDouble> wsp(ntot);

            Vmath::Vsub(ntot, sol, 1, m_phys, 1, wsp, 1);
            Vmath::Vmul(ntot, wsp, 1, wsp, 1, wsp, 1);

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

        NekDouble StdExpansion::L2()
        {       	
            NekDouble  val;
            int     ntot = GetTotPoints();
            Array<OneD, NekDouble> wsp(ntot);

            Vmath::Vmul(ntot, m_phys, 1, m_phys, 1, wsp, 1);
            val   = sqrt(v_Integral(wsp));

            return val;
        }

        NekDouble StdExpansion::H1(const Array<OneD, const NekDouble>& sol)
        {
            int         i;
            NekDouble  val;
            int     ntot = GetTotPoints();
            int     coordim = v_GetCoordim();
            Array<OneD, NekDouble> wsp(3*ntot);
            Array<OneD, NekDouble> wsp_deriv = wsp + ntot;
            Array<OneD, NekDouble> sum = wsp_deriv + ntot;

            Vmath::Vsub(ntot, sol, 1, m_phys, 1, wsp, 1);
            Vmath::Vmul(ntot, wsp, 1, wsp, 1, sum, 1);
            for(i = 0; i < coordim; ++i)
            {
                v_PhysDeriv(i,wsp,wsp_deriv);
                Vmath::Vvtvp(ntot,wsp_deriv,1,wsp_deriv,1,sum,1,sum,1);
            }

            val = sqrt(v_Integral(sum));

            return val;
        }

        NekDouble StdExpansion::H1()
        {
            int i;
            NekDouble  val;
            int     ntot = GetTotPoints();
            int     coordim = v_GetCoordim();
            Array<OneD, NekDouble> wsp_deriv(2*ntot);
            Array<OneD, NekDouble> sum = wsp_deriv + ntot;

            Vmath::Vmul(ntot, m_phys, 1, m_phys, 1, sum, 1);

            for(i = 0; i < coordim; ++i)
            {
                v_PhysDeriv(i,m_phys,wsp_deriv);
                Vmath::Vvtvp(ntot,wsp_deriv,1,wsp_deriv,1,sum,1,sum,1);
            }

            val = sqrt(v_Integral(sum));

            return val;
        }


        DNekMatSharedPtr StdExpansion::CreateGeneralMatrix(const StdMatrixKey &mkey)
        {
            int     i;
            DNekMatSharedPtr  returnval;

            switch(mkey.GetMatrixType())
            {
            case eInvMass:
                {
                    StdMatrixKey masskey(eMass,mkey.GetExpansionType(),*this,NullConstFactorMap,NullVarCoeffMap,mkey.GetNodalPointsType());
                    DNekMatSharedPtr mmat = GetStdMatrix(masskey);
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(*mmat); //Populate standard mass matrix.
                    returnval->Invert();
                }
                break;
            case eInvNBasisTrans:
                {
                    StdMatrixKey tmpkey(eNBasisTrans,mkey.GetExpansionType(),*this,NullConstFactorMap,NullVarCoeffMap,mkey.GetNodalPointsType());
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
            if (mkey.GetNVarCoeff())
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
            Vmath::Vcopy(m_ncoeffs,store2.get(),1,outarray.get(),1);
        }

        void StdExpansion::WeakDerivMatrixOp_MatFree(const int k1,
                                                           const Array<OneD, const NekDouble> &inarray,
                                                           Array<OneD,NekDouble> &outarray,
                                                           const StdMatrixKey &mkey)
        {
            // ASSERTL1(k1 >= 0 && k1 < ExpansionTypeDimMap[v_DetExpansionType()],"invalid first  argument");
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

        void StdExpansion::WeakDirectionalDerivMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                                Array<OneD,NekDouble> &outarray,
                                                                const StdMatrixKey &mkey)
        {
            int nq = GetTotPoints();
//            int varsize = ((mkey.GetVariableCoefficient(0)).num_elements())/dim;
            Array<OneD, NekDouble> tmp(nq);

             v_BwdTrans(inarray,tmp);
            // For Deformed mesh ==============
//            if (varsize==nq)
//            {
//                v_PhysDirectionalDeriv(tmp,mkey.GetVariableCoefficient(0),tmp);
//            }
//
//            // For Regular mesh ==========
//            else
//            {
//                ASSERTL0(false, "Wrong route");
//            }

            v_IProductWRTBase(tmp, outarray);
        }

        void StdExpansion::MassLevelCurvatureMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                               Array<OneD,NekDouble> &outarray,
                                                               const StdMatrixKey &mkey)
      {
///@todo fix this
//          int nqtot = GetTotPoints();
//          int matrixid = mkey.GetMatrixID();
//
//          NekDouble checkweight=0.0;
//          Array<OneD, NekDouble> tmp(nqtot), tan(nqtot), dtan0(nqtot), dtan1(nqtot), weight(nqtot,0.0);
//
//          int gmatnumber = (mkey.GetVariableCoefficient(1)).num_elements();
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

            VarCoeffType varcoefftypes[] = {eVarCoeffVelX, eVarCoeffVelY};

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
                StdMatrixKey mkeylap(eLaplacian,DetExpansionType(),*this);
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
            StdMatrixKey mkeymass(eMass,DetExpansionType(),*this);
            StdMatrixKey mkeylap(eLaplacian,DetExpansionType(),*this);

            MassMatrixOp(inarray,tmp,mkeymass);
            LaplacianMatrixOp(inarray,outarray,mkeylap);

            Blas::Daxpy(m_ncoeffs, lambda, tmp, 1, outarray, 1);
        }

        void StdExpansion::BwdTrans_MatOp(const Array<OneD, const NekDouble>& inarray,
                                          Array<OneD, NekDouble> &outarray)
        {
            int nq = GetTotPoints();
            StdMatrixKey      bwdtransmatkey(eBwdTrans,DetExpansionType(),*this);
            DNekMatSharedPtr  bwdtransmat = GetStdMatrix(bwdtransmatkey);

            Blas::Dgemv('N',nq,m_ncoeffs,1.0,bwdtransmat->GetPtr().get(),
                        nq, inarray.get(), 1, 0.0, outarray.get(), 1);
        }

        //   I/O routine
        void StdExpansion::WriteCoeffsToFile(std::ofstream &outfile)
        {
            int i;
            for(i=0; i<m_ncoeffs; ++i)
            {
                outfile << m_coeffs[i] << std::endl;
            }
        }

        void StdExpansion::WriteTecplotZone(std::ofstream &outfile)
        {
            int i,j;

            int coordim   = GetCoordim();
            int totpoints = GetTotPoints();

            Array<OneD,NekDouble> coords[3];

            coords[0] = Array<OneD,NekDouble>(totpoints);
            coords[1] = Array<OneD,NekDouble>(totpoints);
            coords[2] = Array<OneD,NekDouble>(totpoints);

            GetCoords(coords[0],coords[1],coords[2]);

            switch(DetExpansionType())
            {
                case eSegment:
                    outfile << "Zone, I=" << GetNumPoints(0) << ", F=Block" << std::endl;
                    break;
                case eTriangle: 
                case eQuadrilateral:
                    outfile << "Zone, I=" << GetNumPoints(0) << ", J=" << GetNumPoints(1) <<", F=Block" << std::endl;
                    break;
                case eTetrahedron: 
                case ePrism: 
                case ePyramid: 
                case eHexahedron:
                    outfile << "Zone, I=" << GetNumPoints(0) << ", J=" << GetNumPoints(1) << ", K="<< GetNumPoints(2) << ", F=Block" << std::endl;
                    break;
                default:
                    ASSERTL0(false, "Unsupported expansion type.");
                    break;
            }

            for(j = 0; j < coordim; ++j)
            {
                for(i = 0; i < totpoints; ++i)
                {
                    outfile << coords[j][i] << " ";
                }
                outfile << std::endl;
            }

        }

        void StdExpansion::WriteTecplotField(std::ofstream &outfile)
        {
            int i;

            int totpoints = GetTotPoints();

            // printing the fields of that zone
            for(i = 0; i < totpoints; ++i)
            {
                outfile << m_phys[i] << " ";
            }
            outfile << std::endl;
        }


        // VIRTUAL INLINE FUNCTIONS FROM HEADER FILE
            void StdExpansion::SetUpPhysNormals(const int edge)
            {
                v_SetUpPhysNormals(edge);
            }

	    void StdExpansion::SetUpPhysTangents(const boost::shared_ptr<StdExpansion> &exp2d, const int edge)
	    {
	    	v_SetUpPhysTangents(exp2d, edge);
	    }
            void StdExpansion::AddEdgeNormBoundaryInt(const int edge,
                                                boost::shared_ptr<StdExpansion>    &EdgeExp,
                                                const Array<OneD, const NekDouble> &Fx,
                                                const Array<OneD, const NekDouble> &Fy,
                                                Array<OneD, NekDouble> &outarray)
            {
                v_AddEdgeNormBoundaryInt(edge,EdgeExp,Fx,Fy,outarray);
            }

            void StdExpansion::AddEdgeNormBoundaryInt(const int edge,
                                                boost::shared_ptr<StdExpansion>    &EdgeExp,
                                                const Array<OneD, const NekDouble> &Fn,
                                                Array<OneD, NekDouble> &outarray)
            {
                v_AddEdgeNormBoundaryInt(edge,EdgeExp,Fn,outarray);
            }

            void StdExpansion::AddFaceNormBoundaryInt(const int face,
                                                boost::shared_ptr<StdExpansion>    &FaceExp,
                                                const Array<OneD, const NekDouble> &Fn,
                                                Array<OneD, NekDouble> &outarray)
            {
                v_AddFaceNormBoundaryInt(face,FaceExp,Fn,outarray);
            }

            const Array<OneD, const NekDouble>& StdExpansion::v_GetPhysNormals(void)
            {
                NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
                return NullNekDouble1DArray;
            }


            void StdExpansion::v_SetPhysNormals(Array<OneD, const NekDouble> &normal)
            {
                NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
            }

            void StdExpansion::v_SetUpPhysNormals(const int edge)
            {
                NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
            }

            void StdExpansion::v_SetUpPhysTangents(const boost::shared_ptr<StdExpansion> &exp2d, const int edge)
	    {
                NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
	    }

        int StdExpansion::v_CalcNumberOfCoefficients(const std::vector<unsigned int>  &nummodes, int &modes_offset)
            {
                NEKERROR(ErrorUtil::efatal, "This function is not defined for this class");
                return 0;
            }
            
        void StdExpansion::v_ExtractDataToCoeffs(const NekDouble *data, 
                                   const std::vector<unsigned int > &nummodes, 
                                   const int nmode_offset,
                                   NekDouble *coeffs)
        {
                NEKERROR(ErrorUtil::efatal, "This function is not defined for this class");
        }



            void StdExpansion::v_NormVectorIProductWRTBase(const Array<OneD, const NekDouble> &Fx, const Array<OneD, const NekDouble> &Fy, Array< OneD, NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
            }

        void StdExpansion::v_NormVectorIProductWRTBase(
                              const Array<OneD, const NekDouble> &Fx, 
                              const Array<OneD, const NekDouble> &Fy, 
                              const Array<OneD, const NekDouble> &Fz, 
                              Array< OneD, NekDouble> &outarray)
        {
            NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
        }

            DNekScalBlkMatSharedPtr StdExpansion::v_GetLocStaticCondMatrix(const LocalRegions::MatrixKey &mkey)
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
                return NullDNekScalBlkMatSharedPtr;
            }


            StdRegions::Orientation StdExpansion::v_GetFaceOrient(int face)

            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for three-dimensional  LocalRegions");
                return eDir1FwdDir1_Dir2FwdDir2;
            }

            StdRegions::Orientation StdExpansion::v_GetEorient(int edge)
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for two-dimensional  LocalRegions");
                return eForwards;
            }
		
			StdRegions::Orientation StdExpansion::v_GetPorient(int point)
			{
				NEKERROR(ErrorUtil::efatal, "This function is only valid for one-dimensional  LocalRegions");
				return eFwd;
			}


            StdRegions::Orientation StdExpansion::v_GetCartesianEorient(int edge)
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for two-dimensional  LocalRegions");
                return eForwards;
            }


        void StdExpansion::v_SetCoeffsToOrientation(StdRegions::Orientation dir,
                                                    Array<OneD, const NekDouble> &inarray,
                                                    Array<OneD, NekDouble> &outarray)
        {
            NEKERROR(ErrorUtil::efatal, "This function is not defined for this shape");
        }

        void StdExpansion::v_SetCoeffsToOrientation(StdRegions::Orientation dir)
        {
            NEKERROR(ErrorUtil::efatal, "This function is not defined for this shape");
        }
        

            void StdExpansion::v_AddEdgeNormBoundaryInt(const int edge,
                                                        boost::shared_ptr<StdExpansion>    &EdgeExp,
                                                        const Array<OneD, const NekDouble> &Fx,
                                                        const Array<OneD, const NekDouble> &Fy,
                                                        Array<OneD, NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal, "This function is not defined for this shape");
            }

             void StdExpansion::v_AddEdgeNormBoundaryInt(const int edge,
                                                  boost::shared_ptr<StdExpansion>    &EdgeExp,
                                                  const Array<OneD, const NekDouble> &Fn,
                                                  Array<OneD, NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal, "This function is not defined for this shape");
            }

             void StdExpansion::v_AddFaceNormBoundaryInt(const int face,
                                                  boost::shared_ptr<StdExpansion>    &FaceExp,
                                                  const Array<OneD, const NekDouble> &Fn,
                                                  Array<OneD, NekDouble> &outarray)
            {
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
                ASSERTL0(false, "This function is not valid or not defined");
                return 0;
            }

            int StdExpansion::v_DetCartesianDirOfEdge(const int edge)
            {
                ASSERTL0(false, "This function is not valid or not defined");
                return 0;
            }

            const LibUtilities::BasisKey StdExpansion::v_DetEdgeBasisKey(const int i) const
            {
                ASSERTL0(false, "This function is not valid or not defined");
                return LibUtilities::NullBasisKey;
            }

            const LibUtilities::BasisKey StdExpansion::v_DetFaceBasisKey(const int i, const int k) const
            {
                ASSERTL0(false, "This function is not valid or not defined");
                return LibUtilities::NullBasisKey;
            }

            int StdExpansion::v_GetFaceNumPoints(const int i) const
            {
                ASSERTL0(false, "This function is not valid or not defined");
                return 0;
            }

            int StdExpansion::v_GetFaceNcoeffs(const int i) const
            {
                ASSERTL0(false, "This function is not valid or not defined");
                return 0;
            }

            int StdExpansion::v_GetFaceIntNcoeffs(const int i) const
            {
                ASSERTL0(false, "This function is not valid or not defined");
                return 0;
            }

            int StdExpansion::v_GetTotalFaceIntNcoeffs() const
            {
                ASSERTL0(false, "This function is not valid or not defined");
                return 0;
            }        

            LibUtilities::PointsKey StdExpansion::v_GetFacePointsKey(const int i, const int j) const
            {
                ASSERTL0(false, "This function is not valid or not defined");
                return LibUtilities::NullPointsKey;
            }

            LibUtilities::BasisType StdExpansion::v_GetEdgeBasisType(const int i) const
            {
                ASSERTL0(false, "This function is not valid or not defined");

                return LibUtilities::eNoBasisType;
            }

            ExpansionType StdExpansion::v_DetExpansionType() const
            {
                ASSERTL0(false, "This expansion does not have a shape type defined");
                return eNoExpansionType;
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

            void  StdExpansion::v_IProductWRTDerivBase (const int dir,
                                                   const Array<OneD, const NekDouble>& inarray,
                                                   Array<OneD, NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal, "This method has not been defined");
            }

            /**
             *
             */
            void StdExpansion::v_FwdTrans_BndConstrained(const Array<OneD, const NekDouble>& inarray,
                                                   Array<OneD, NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal, "This method has not been defined");
            }


            /**
             * @brief Integrates the specified function over the domain.
             * @see StdRegions#StdExpansion#Integral.
             */
            NekDouble StdExpansion::v_Integral(const Array<OneD, const NekDouble>& inarray )
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                         "local expansions");
                return 0;
            }


            void StdExpansion::v_AddRobinMassMatrix(const int edgeid, const Array<OneD, const NekDouble > &primCoeffs, DNekMatSharedPtr &inoutmat)
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                         "specific element types");
            }

            void StdExpansion::v_AddRobinEdgeContribution(const int edgeid, const Array<OneD, const NekDouble > &primCoeffs, Array<OneD, NekDouble> &coeffs)
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                         "specific element types");
            }

            void StdExpansion::v_DGDeriv(const int dir,
                                         const Array<OneD, const NekDouble>& inarray,
                                         Array<OneD, boost::shared_ptr< StdExpansion > > &EdgeExp,
                                         Array<OneD, NekDouble> &outarray)
            {
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
                NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                         "local expansions");
            }

            void StdExpansion::v_PhysDeriv_s(const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> &out_ds)
            {
                    NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                            "local expansions");
            }
            void StdExpansion::v_PhysDeriv_n(const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble>& out_dn)
            {
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
                NEKERROR(ErrorUtil::efatal, "This function is only valid for "
                         "specific element types");
            }

            void StdExpansion::v_StdPhysDeriv (const Array<OneD, const NekDouble>& inarray,
                                         Array<OneD, NekDouble> &out_d1,
                                         Array<OneD, NekDouble> &out_d2,
                                         Array<OneD, NekDouble> &out_d3)
            {
                NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
            }

            void   StdExpansion::v_StdPhysDeriv (const int dir,
                                           const Array<OneD, const NekDouble>& inarray,
                                           Array<OneD, NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
            }

            NekDouble StdExpansion::v_PhysEvaluate(const Array<OneD, const NekDouble>& coords)
            {
                NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
                return 0;
            }


        NekDouble StdExpansion::v_PhysEvaluate(const Array<OneD, const NekDouble>& coords, const Array<OneD, const NekDouble>& physvals)
            {
                NEKERROR(ErrorUtil::efatal, "Method does not exist for this shape");
                return 0;
            }


            void StdExpansion::v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal, "This function has not "
                         "been defined for this shape");
            }

            DNekMatSharedPtr StdExpansion::v_GenMatrix(const StdMatrixKey &mkey)
            {
                NEKERROR(ErrorUtil::efatal, "This function has not "
                         "been defined for this element");
                DNekMatSharedPtr returnval;
                return returnval;
            }

            DNekMatSharedPtr StdExpansion::v_CreateStdMatrix(const StdMatrixKey &mkey)
            {
                NEKERROR(ErrorUtil::efatal, "This function has not "
                         "been defined for this element");
                DNekMatSharedPtr returnval;
                return returnval;
            }

            void StdExpansion::v_GetCoords(Array<OneD, NekDouble> &coords_0,
                                     Array<OneD, NekDouble> &coords_1,
                                     Array<OneD, NekDouble> &coords_2)
            {
                NEKERROR(ErrorUtil::efatal, "Write coordinate definition method");
            }

            void StdExpansion::v_GetCoord(const Array<OneD, const NekDouble>& Lcoord,
                                    Array<OneD, NekDouble> &coord)
            {
                NEKERROR(ErrorUtil::efatal, "Write coordinate definition method");
            }

            int StdExpansion::v_GetCoordim(void)
            {
                NEKERROR(ErrorUtil::efatal, "Write method");
                return -1;
            }

            void StdExpansion::v_GetBoundaryMap(Array<OneD, unsigned int>& outarray)
            {
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
            }

            void StdExpansion::v_GetInteriorMap(Array<OneD, unsigned int>& outarray)
            {
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
            }

            int StdExpansion::v_GetVertexMap(const int localVertexId)
            {
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
                return 0;
            }

            void StdExpansion::v_GetEdgeInteriorMap(const int eid, const Orientation edgeOrient,
                                              Array<OneD, unsigned int> &maparray,
                                              Array<OneD, int> &signarray)
            {
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
            }

            void StdExpansion::v_GetFaceInteriorMap(const int fid, const Orientation faceOrient,
                                              Array<OneD, unsigned int> &maparray,
                                              Array<OneD, int> &signarray)
            {
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
            }

            void StdExpansion::v_GetEdgeToElementMap(const int eid, const Orientation edgeOrient,
                                               Array<OneD, unsigned int> &maparray,
                                               Array<OneD, int> &signarray)
            {
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
            }

            void StdExpansion::v_GetFaceToElementMap(const int fid, const Orientation faceOrient,
                                                     Array<OneD, unsigned int> &maparray,
                                                     Array<OneD, int> &signarray,
                                                     int nummodesA, int nummodesB)
            {
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
            }


            void StdExpansion::v_GetEdgePhysVals(const int edge, const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape or library" );
            }

            void StdExpansion::v_GetEdgePhysVals(const int edge,  const boost::shared_ptr<StdExpansion>  &EdgeExp, const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape or library" );
            }
        
            void StdExpansion::v_GetEdgeQFactors(
                    const int edge,  
                    Array<OneD, NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal,
                     "Method does not exist for this shape or library");
            }

            void StdExpansion::v_GetFacePhysVals(
                const int                                face,
                const boost::shared_ptr<StdExpansion>   &FaceExp,
                const Array<OneD, const NekDouble>      &inarray,
                      Array<OneD,       NekDouble>      &outarray,
                StdRegions::Orientation                  orient)
            {
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape or library" );
            }

            void StdExpansion::v_WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar, std::string var)
            {
                NEKERROR(ErrorUtil::efatal, "WriteToFile: Write method");
            }

            void StdExpansion::v_ReadFromFile(std::ifstream &infile, OutputFormat format, const bool dumpVar)
            {
                NEKERROR(ErrorUtil::efatal, "ReadFromFile: Write method");
            }

            const  boost::shared_ptr<SpatialDomains::GeomFactors>& StdExpansion::v_GetMetricInfo() const
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
                return SpatialDomains::NullGeomFactorsSharedPtr;

            }

            const boost::shared_ptr<SpatialDomains::Geometry> StdExpansion::v_GetGeom() const
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");

                return SpatialDomains::NullGeometrySharedPtr;
            }

            const boost::shared_ptr<SpatialDomains::Geometry1D>& StdExpansion::v_GetGeom1D() const
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");

                return SpatialDomains::NullGeometry1DSharedPtr;
            }

            const boost::shared_ptr<SpatialDomains::Geometry2D>& StdExpansion::v_GetGeom2D() const
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");

                return SpatialDomains::NullGeometry2DSharedPtr;
            }

            const boost::shared_ptr<SpatialDomains::Geometry3D>& StdExpansion::v_GetGeom3D() const
            {
                NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");

                return SpatialDomains::NullGeometry3DSharedPtr;
            }

            void StdExpansion::v_BwdTrans_SumFac(const Array<OneD, const NekDouble>& inarray,
                                           Array<OneD, NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
            }

            void StdExpansion::v_IProductWRTBase_SumFac(const Array<OneD, const NekDouble>& inarray,
                                                  Array<OneD, NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
            }

            void StdExpansion::v_IProductWRTDerivBase_SumFac(const int dir,
                                                       const Array<OneD, const NekDouble>& inarray,
                                                       Array<OneD, NekDouble> &outarray)
            {
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
                ASSERTL0(false, "Cannot get edge normals for this expansion.");
                static NormalVector result;
                return result;
            }

            void StdExpansion::v_ComputeEdgeNormal(const int edge)
            {
                ASSERTL0(false, "Cannot compute edge normal for this expansion.");
            }

            void StdExpansion::v_ComputeFaceNormal(const int face)
            {
                ASSERTL0(false, "Cannot compute face normal for this expansion.");
            }

            void StdExpansion::v_NegateEdgeNormal(const int edge)
            {
                ASSERTL0(false, "Not implemented.");
            }
        
            void StdExpansion::v_NegateFaceNormal(const int face)
            {
                ASSERTL0(false, "Not implemented.");
            }

		
			void StdExpansion::v_ComputeVertexNormal(const int vertex)
			{
				ASSERTL0(false, "Cannot compute vertex normal for this expansion.");
			}

            const NormalVector & StdExpansion::v_GetFaceNormal(const int face) const
            {
                ASSERTL0(false, "Cannot get face normals for this expansion.");
                static NormalVector result;
                return result;
            }
		
			const NormalVector & StdExpansion::v_GetVertexNormal(const int vertex) const
			{
				ASSERTL0(false, "Cannot get vertex normals for this expansion.");
                static NormalVector result;
                return result;
			}	
		
            const NormalVector & StdExpansion::v_GetSurfaceNormal() const
            {
                ASSERTL0(false, "Cannot get face normals for this expansion.");
                static NormalVector result;
                return result;
            }
    }//end namespace
}//end namespace
