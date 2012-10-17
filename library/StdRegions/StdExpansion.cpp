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
            MatrixType mtype = mkey.GetMatrixType();

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
					ASSERTL0(false,"Face to Element Index Map not implemented yet.");
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
            int dim = 3;
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
            case eTriangle: case eQuadrilateral:

                outfile << "Zone, I=" << GetNumPoints(0) << ", J=" << GetNumPoints(1) <<", F=Block" << std::endl;
                break;
            case eTetrahedron: case ePrism: case ePyramid: case eHexahedron:
                outfile << "Zone, I=" << GetNumPoints(0) << ", J=" << GetNumPoints(1) << ", K="<< GetNumPoints(2) << ", F=Block" << std::endl;
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

            void StdExpansion::AddEdgeNormBoundaryBiInt(const int edge,
                                                boost::shared_ptr<StdExpansion>    &EdgeExp,
                                                const Array<OneD, const NekDouble> &Fwd,
                                                const Array<OneD, const NekDouble> &Bwd,
                                                Array<OneD, NekDouble> &outarray)
            {
                v_AddEdgeNormBoundaryBiInt(edge,EdgeExp,Fwd,Bwd,outarray);
            }

            void StdExpansion::AddNormTraceInt(const int dir,
                                         Array<OneD, const NekDouble> &inarray,
                                         Array<OneD,NekDouble> &outarray)
            {
                v_AddNormTraceInt(dir,inarray,outarray);
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
            
        void StdExpansion::v_ExtractDataToCoeffs(const std::vector<NekDouble> &data, 
                                   const int offset, 
                                   const std::vector<unsigned int > &nummodes, 
                                   const int nmode_offset,
                                   Array<OneD, NekDouble> &coeffs)
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


            void StdExpansion::v_AddHDGHelmholtzTraceTerms(const NekDouble tau,
                                                     const Array<OneD, const NekDouble> &inarray,
                                                     Array<OneD,NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal, "This function is not defined for this shape");
            }


            void StdExpansion::v_AddHDGHelmholtzTraceTerms(const NekDouble tau,
                                                     const Array<OneD, const NekDouble> &inarray,
                                                     Array<OneD, boost::shared_ptr< StdExpansion1D > > &edgeExp,
                                                     Array<OneD,NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal, "This function is not defined for this shape");
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

            void StdExpansion::v_AddEdgeNormBoundaryBiInt(const int edge,
                                                    boost::shared_ptr<StdExpansion>    &EdgeExp,
                                                    const Array<OneD, const NekDouble> &Fwd,
                                                    const Array<OneD, const NekDouble> &Bwd,
                                                    Array<OneD, NekDouble> &outarray)
            {
                NEKERROR(ErrorUtil::efatal, "v_AddEdgeNormBoundaryBiInt is not defined for this shape");
            }

           void StdExpansion::v_AddNormTraceInt(const int dir,
                                           Array<OneD, const NekDouble> &inarray,
                                           Array<OneD,NekDouble> &outarray)
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

/**
* $Log: StdExpansion.cpp,v $
* Revision 1.94  2010/01/03 19:39:09  cantwell
* Added FldToVtk converter.
* Added XmlToVtk converter.
*
* Revision 1.93  2009/12/18 00:11:03  bnelson
* Fixed windows compiler warnings.
*
* Revision 1.92  2009/12/14 18:03:18  cbiotto
* Adding functions for tecplot file
*
* Revision 1.91  2009/11/13 16:17:46  sehunchun
* *** empty log message ***
*
* Revision 1.90  2009/11/11 18:43:58  sehunchun
* *** empty log message ***
*
* Revision 1.89  2009/11/10 19:01:37  sehunchun
* Update related to Variable coefficients of HDG2D Solver
*
* Revision 1.88  2009/11/06 21:42:16  sherwin
* Added call to DGDeriv function
*
* Revision 1.87  2009/11/02 19:15:43  cantwell
* Moved ContField1D to inherit from DisContField1D.
* Moved ContField3D to inherit from DisContField3D.
* Incorporated GenExpList1D functionality into ExpList1D.
* Tidied up and added documentation to various classes.
* Moved Namespace documentation and introductions to separate files along with
* doxygen configuration.
* Added option to use system ZLIB library instead of libboost_zlib on UNIX.
* Added extra search paths to FindMetis.cmake and FindNektar++.cmake.
* Updated Linux compiling instructions.
* Updated regDemo to use Helmholtz2D-g when built as debug.
*
* Revision 1.86  2009/10/30 14:01:19  pvos
* Multi-level static condensation updates
*
* Revision 1.85  2009/10/25 18:53:39  sherwin
* Added H1 norm definition
*
* Revision 1.84  2009/07/09 21:43:29  sehunchun
* Mass Matrix multiplicatin with variablecoefficient
*
* Revision 1.83  2009/07/02 13:27:51  sehunchun
* Unnecessary error message for 2D restriction is deleted
*
* Revision 1.82  2009/04/27 09:20:21  pvos
* Fixed small bug
*
* Revision 1.81  2009/04/03 14:57:34  sherwin
* Linear Advection matrices added, corrected unsigned int intialisation
*
* Revision 1.80  2008/11/24 10:31:14  pvos
* Changed name from _PartitionedOp to _MatFree
*
* Revision 1.79  2008/11/19 16:02:47  pvos
* Added functionality for variable Laplacian coeffcients
*
* Revision 1.78  2008/11/07 11:37:31  pvos
* Minor correction
*
* Revision 1.77  2008/11/05 16:08:15  pvos
* Added elemental optimisation functionality
*
* Revision 1.76  2008/09/09 14:12:51  sherwin
* Removed Interp1D/2D/3D and put the into LibUtilities
*
* Revision 1.75  2008/08/14 22:09:50  sherwin
* Modifications to remove HDG routines from StdRegions and removed StdExpMap
*
* Revision 1.74  2008/07/29 22:21:15  sherwin
* A bunch of mods for DG advection and separaring the GetGeom calls into GetGeom1D ...
*
* Revision 1.73  2008/07/19 21:12:54  sherwin
* Removed MapTo function and made orientation convention anticlockwise in UDG routines
*
* Revision 1.72  2008/07/12 16:30:07  sherwin
* Added an new member m_elmt_id so that there is an element number for use later in lists
*
* Revision 1.71  2008/07/02 14:08:56  pvos
* Implementation of HelmholtzMatOp and LapMatOp on shape level
*
* Revision 1.70  2008/05/30 00:33:49  delisi
* Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
*
* Revision 1.69  2008/05/10 18:27:33  sherwin
* Modifications necessary for QuadExp Unified DG Solver
*
* Revision 1.68  2008/04/06 06:04:14  bnelson
* Changed ConstArray to Array<const>
*
* Revision 1.67  2008/04/02 22:18:10  pvos
* Update for 2D local to global mapping
*
* Revision 1.66  2008/03/18 14:15:45  pvos
* Update for nodal triangular helmholtz solver
*
* Revision 1.65  2008/03/12 15:25:09  pvos
* Clean up of the code
*
* Revision 1.63  2008/02/29 19:15:19  sherwin
* Update for UDG stuff
*
* Revision 1.62  2008/02/16 05:59:14  ehan
* Added interpolation 3D.
*
* Revision 1.61  2008/01/23 09:09:46  sherwin
* Updates for Hybrized DG
*
* Revision 1.60  2007/12/17 13:03:45  sherwin
* Modified StdMatrixKey to contain a list of constants and GenMatrix to take a StdMatrixKey
*
* Revision 1.59  2007/12/06 22:44:46  pvos
* 2D Helmholtz solver updates
*
* Revision 1.58  2007/11/29 21:40:20  sherwin
* updates for MultiRegions and DG solver
*
* Revision 1.57  2007/11/08 16:55:12  pvos
* Updates towards 2D helmholtz solver
*
* Revision 1.56  2007/10/15 20:38:32  ehan
* Tested standard mass matrix
*
* Revision 1.55  2007/10/04 12:10:04  sherwin
* Update for working version of static condensation in Helmholtz1D and put lambda coefficient on the mass matrix rather than the Laplacian operator.
*
* Revision 1.54  2007/10/03 11:37:51  sherwin
* Updates relating to static condensation implementation
*
* Revision 1.53  2007/09/27 12:55:57  pvos
* Column major Blas calls corrections
*
* Revision 1.52  2007/09/25 14:25:56  pvos
* Update for helmholtz1D with different expansion orders
*
* Revision 1.51  2007/08/29 23:26:48  jfrazier
* Created non-static manager that shares data across instances.
*
* Revision 1.50  2007/07/27 16:56:50  jfrazier
* Changed manager to static.
*
* Revision 1.49  2007/07/27 00:22:53  bnelson
* Memory manager now accepts non-const parameters to the allocate methods.
*
* Revision 1.48  2007/07/22 23:04:25  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.47  2007/07/16 18:28:43  sherwin
* Modification to introduce non-zero Dirichlet boundary conditions into the Helmholtz1D Demo
*
* Revision 1.46  2007/07/15 19:28:28  bnelson
* *** empty log message ***
*
* Revision 1.45  2007/07/13 15:20:19  kirby
* *** empty log message ***
*
* Revision 1.43  2007/07/13 09:02:25  sherwin
* Mods for Helmholtz solver
*
* Revision 1.42  2007/07/12 12:55:14  sherwin
* Simplified Matrix Generation
*
* Revision 1.41  2007/07/10 20:41:52  kirby
* more fixes
*
* Revision 1.40  2007/07/10 19:27:58  kirby
* Update for new matrix structures
*
* Revision 1.39  2007/07/09 15:19:14  sherwin
* Introduced an InvMassMatrix and replaced the StdLinSysManager call with a StdMatrixManager call to the inverse matrix
*
* Revision 1.38  2007/06/07 15:54:19  pvos
* Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
* Also made corrections to various ASSERTL2 calls
*
* Revision 1.37  2007/05/30 23:56:54  sherwin
* Silly errors
*
* Revision 1.36  2007/05/30 20:49:12  sherwin
* Updates to do with LocalRegions and SpatialDomains
*
* Revision 1.35  2007/05/23 15:12:45  pvos
* removed some obsolete lines
*
* Revision 1.34  2007/05/15 05:18:23  bnelson
* Updated to use the new Array object.
*
* Revision 1.33  2007/04/26 15:00:17  sherwin
* SJS compiling working version using SHaredArrays
*
* Revision 1.32  2007/04/18 16:09:12  pvos
* Added some new Tensor Operations routines
*
* Revision 1.31  2007/04/10 14:00:45  sherwin
* Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
*
* Revision 1.30  2007/04/08 03:36:57  jfrazier
* Updated to use SharedArray consistently and minor reformatting.
*
* Revision 1.29  2007/03/31 00:40:02  bnelson
* *** empty log message ***
*
* Revision 1.28  2007/03/29 19:35:08  bnelson
* Replaced boost::shared_array with SharedArray
*
* Revision 1.27  2007/03/25 15:48:22  sherwin
* UPdate LocalRegions to take new NekDouble and shared_array formats. Added new Demos
*
* Revision 1.26  2007/03/21 20:56:42  sherwin
* Update to change BasisSharedVector to boost::shared_array<BasisSharedPtr> and removed tthe Vector definitions in GetCoords and PhysDeriv
*
* Revision 1.25  2007/03/20 16:58:42  sherwin
* Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
*
* Revision 1.24  2007/03/14 21:24:09  sherwin
* Update for working version of MultiRegions up to ExpList1D
*
* Revision 1.23  2007/03/02 12:01:51  sherwin
* Update for working version of LocalRegions/Project1D
*
* Revision 1.22  2007/02/28 19:05:11  sherwin
* Moved key definitions to their own files to make things more transparent
*
* Revision 1.21  2007/02/24 09:07:25  sherwin
* Working version of stdMatrixManager and stdLinSysMatrix
*
* Revision 1.20  2007/02/23 19:26:07  jfrazier
* General bug fix and formatting.
*
* Revision 1.19  2007/02/22 22:02:27  sherwin
* Update with executing StdMatManager
*
* Revision 1.18  2007/02/22 18:11:31  sherwin
* Version with some create functions introduced for StdMatManagers
*
* Revision 1.17  2007/02/21 22:55:16  sherwin
* First integration of StdMatrixManagers
*
* Revision 1.16  2007/02/17 04:03:22  jfrazier
* Added NekManager for holding matrices.  Need to finish the create function.
*
* Revision 1.15  2007/02/14 16:35:50  pvos
* Corrected an error in the code
*
* Revision 1.14  2007/02/13 09:52:27  sherwin
* Updates to fix mass matrix inverse issues
*
* Revision 1.13  2007/02/07 12:51:52  sherwin
* Compiling version of Project1D
*
* Revision 1.12  2007/02/06 02:23:28  jfrazier
* Minor cleanup.
*
* Revision 1.11  2007/01/30 20:01:35  sherwin
* Update for first compiling Project1D routine
*
* Revision 1.10  2007/01/29 15:04:53  sherwin
* StdBasis.h moved to LibUtilities. Other minor mods
*
* Revision 1.9  2007/01/28 18:34:18  sherwin
* More modifications to make Demo Project1D compile
*
* Revision 1.8  2007/01/23 23:20:20  sherwin
* New version after Jan 07 update
*
* Revision 1.7  2007/01/20 22:35:20  sherwin
* Version with StdExpansion compiling
*
* Revision 1.6  2007/01/15 11:08:37  pvos
* Updating doxygen documentation
*
* Revision 1.5  2006/12/10 19:00:54  sherwin
* Modifications to handle nodal expansions
*
* Revision 1.4  2006/08/05 19:03:48  sherwin
* Update to make the multiregions 2D expansion in connected regions work
*
* Revision 1.3  2006/06/01 14:46:16  kirby
* *** empty log message ***
*
* Revision 1.2  2006/05/29 19:03:08  sherwin
* Modifications to wrap geometric information in shared_ptr
*
* Revision 1.1  2006/05/04 18:58:31  kirby
* *** empty log message ***
*
* Revision 1.54  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.53  2006/04/01 21:59:26  sherwin
* Sorted new definition of ASSERT
*
* Revision 1.52  2006/03/21 09:21:31  sherwin
* Introduced NekMemoryManager
*
* Revision 1.51  2006/03/12 14:20:44  sherwin
*
* First compiling version of SpatialDomains and associated modifications
*
* Revision 1.50  2006/03/06 12:39:59  sherwin
*
* Added NekConstants class for all constants in this library
*
* Revision 1.49  2006/03/05 22:11:02  sherwin
*
* Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
*
* Revision 1.48  2006/03/03 23:04:54  sherwin
*
* Corrected Mistake in StdBasis.cpp to do with eModified_B
*
* Revision 1.47  2006/03/02 16:20:20  sherwin
*
* Introduced method GetPointsTot
*
* Revision 1.46  2006/03/01 08:25:03  sherwin
*
* First compiling version of StdRegions
*
* Revision 1.45  2006/02/26 23:37:29  sherwin
*
* Updates and compiling checks upto StdExpansions1D
*
* Revision 1.44  2006/02/26 21:23:20  bnelson
* Fixed a variety of compiler errors caused by updates to the coding standard.
*
* Revision 1.43  2006/02/15 08:06:36  sherwin
*
* Put files into coding standard (although they do not compile)
*
* Revision 1.42  2006/02/12 21:51:42  sherwin
*
* Added licence
*
* Revision 1.41  2006/02/10 16:44:10  sherwin
*
* Updated to comply with coding standard
*
**/
