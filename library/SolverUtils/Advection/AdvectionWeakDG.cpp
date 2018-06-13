///////////////////////////////////////////////////////////////////////////////
//
// File: AdvectionWeakDG.cpp
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
// Description: Weak DG advection class.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Advection/AdvectionWeakDG.h>
#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/Expansion3D.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/Expansion0D.h>
#include <iostream>
#include <iomanip>
#include <math.h>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string AdvectionWeakDG::type = GetAdvectionFactory().
            RegisterCreatorFunction("WeakDG", AdvectionWeakDG::create);

        AdvectionWeakDG::AdvectionWeakDG()
        {
        }

        /**
         * @brief Initialise AdvectionWeakDG objects and store them before
         * starting the time-stepping.
         *
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         */
        void AdvectionWeakDG::v_InitObject(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            Advection::v_InitObject(pSession, pFields);
        }

        /**
         * @brief Compute the advection term at each time-step using the
         * Discontinuous Galerkin approach (DG).
         *
         * @param nConvectiveFields   Number of fields.
         * @param fields              Pointer to fields.
         * @param advVel              Advection velocities.
         * @param inarray             Solution at the previous time-step.
         * @param outarray            Advection term to be passed at the
         *                            time integration class.
         */
        void AdvectionWeakDG::v_Advect(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &advVel,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const NekDouble                                   &time,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            int nPointsTot      = fields[0]->GetTotPoints();
            int nCoeffs         = fields[0]->GetNcoeffs();
            int nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();
            int i, j;
            
            Array<OneD, Array<OneD, NekDouble> > tmp(nConvectiveFields);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > fluxvector(
                nConvectiveFields);

            // Allocate storage for flux vector F(u).
            for (i = 0; i < nConvectiveFields; ++i)
            {
                fluxvector[i] =
                    Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
                for (j = 0; j < m_spaceDim; ++j)
                {
                    fluxvector[i][j] = Array<OneD, NekDouble>(nPointsTot);
                }
            }

            ASSERTL1(m_riemann,
                     "Riemann solver must be provided for AdvectionWeakDG.");

            m_fluxVector(inarray, fluxvector);

            // Get the advection part (without numerical flux)
            for(i = 0; i < nConvectiveFields; ++i)
            {
                tmp[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);

                fields[i]->IProductWRTDerivBase(fluxvector[i],tmp[i]);
            }

            // Store forwards/backwards space along trace space
            Array<OneD, Array<OneD, NekDouble> > Fwd    (nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> > Bwd    (nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> > numflux(nConvectiveFields);

            if (pFwd == NullNekDoubleArrayofArray ||
                pBwd == NullNekDoubleArrayofArray)
            {
                for(i = 0; i < nConvectiveFields; ++i)
                {
                    Fwd[i]     = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                    Bwd[i]     = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                    numflux[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                    fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
                }
            }
            else
            {
                for(i = 0; i < nConvectiveFields; ++i)
                {
                    Fwd[i]     = pFwd[i];
                    Bwd[i]     = pBwd[i];
                    numflux[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                }
            }

            m_riemann->Solve(m_spaceDim, Fwd, Bwd, numflux);

  /*           

            if(FALSE)
            {

                Array<OneD, unsigned int> n_blks(nTracePointsTot);
                for(int i=0;i<nTracePointsTot;i++)
                {
                    n_blks[i]    = nConvectiveFields;
                }
                DNekBlkMatSharedPtr FJac = MemoryManager<DNekBlkMat>
                    ::AllocateSharedPtr(n_blks, n_blks, eDIAGONAL);
                DNekBlkMatSharedPtr BJac = MemoryManager<DNekBlkMat>
                    ::AllocateSharedPtr(n_blks, n_blks, eDIAGONAL);
                m_riemann->CalcFluxJacobian(m_spaceDim, Fwd, Bwd, FJac,BJac);




                Array<OneD, NekDouble> PointFwd(nConvectiveFields,0.0),PointBwd(nConvectiveFields,0.0);
                Array<OneD, NekDouble> PntFluxFwd(nConvectiveFields,0.0);
                Array<OneD, NekDouble> PntFluxBwd(nConvectiveFields,0.0);
                Array<OneD, NekDouble> PointFlux(nConvectiveFields,0.0);
                int cnt=0;
                for(int i=0; i<nTracePointsTot;  i++)
                {
                    for(int j=0; j<nConvectiveFields;j++)
                    {
                        PointFwd[j] = Fwd[j][i];
                        PointBwd[j] = Bwd[j][i];
                    }
                    NekVector<NekDouble> VectFluxFwd(nConvectiveFields,PntFluxFwd,eWrapper);
                    NekVector<NekDouble> VectFluxBwd(nConvectiveFields,PntFluxBwd,eWrapper);

                    DNekMat &MF = (*FJac->GetBlock(i,i));
                    NekVector<NekDouble> VectFwd(nConvectiveFields,PointFwd,eWrapper);
                    VectFluxFwd = MF * VectFwd;


                    DNekMat &MB = (*BJac->GetBlock(i,i));
                    NekVector<NekDouble> VectBwd(nConvectiveFields,PointBwd,eWrapper);
                    VectFluxBwd = MB * VectBwd;

                    NekDouble error=0.0;
                    for(int j=0;j<nConvectiveFields;j++)
                    {
                        PointFlux[j] = PntFluxFwd[j]+PntFluxBwd[j];
                        error += abs(PointFlux[j]-numflux[j][i]);
                    }

                    if(error>1.0E-7)
                    {
                        cnt++;
                        std::cout   <<std::scientific<<std::setw(12)<<std::setprecision(5)
                                <<"abs(PointFlux[0]-numflux[0][i])   =   "<<abs(PointFlux[0]-numflux[0][i])<<"    "<<PointFlux[0]<<"    "<<numflux[0][i]<<std::endl
                                <<"abs(PointFlux[1]-numflux[1][i])   =   "<<abs(PointFlux[1]-numflux[1][i])<<"    "<<PointFlux[1]<<"    "<<numflux[1][i]<<std::endl
                                <<"abs(PointFlux[2]-numflux[2][i])   =   "<<abs(PointFlux[2]-numflux[2][i])<<"    "<<PointFlux[2]<<"    "<<numflux[2][i]<<std::endl
                                <<"abs(PointFlux[3]-numflux[3][i])   =   "<<abs(PointFlux[3]-numflux[3][i])<<"    "<<PointFlux[3]<<"    "<<numflux[3][i]<<std::endl;
                    }
                }
                
                std::cout   <<"cnt= "<<cnt<<std::endl;
                
                int j = 0;
            }

 */

            // Evaulate <\phi, \hat{F}\cdot n> - OutField[i]
            for(i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Neg                      (nCoeffs, tmp[i], 1);
                fields[i]->AddTraceIntegral     (numflux[i], tmp[i]);
                fields[i]->MultiplyByElmtInvMass(tmp[i], tmp[i]);
                fields[i]->BwdTrans             (tmp[i], outarray[i]);
            }
        }

        
        void AdvectionWeakDG::v_AddVolumJac2Mat( const int nConvectiveFields,
                                        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                                        const   Array<OneD, const  Array<OneD, DNekMatSharedPtr> >&ElmtJac,
                                        const int nDirctn, 
                                        Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray)
        {
            MultiRegions::ExpListSharedPtr explist = pFields[0];
                std::shared_ptr<LocalRegions::ExpansionVector> pexp = explist->GetExp();
            int ntotElmt            = (*pexp).size();
            int nElmtPnt,nElmtCoef;

            NekDouble tmp;
            DNekMatSharedPtr        tmpGmtx,ElmtMat;

            Array<OneD, DNekMatSharedPtr>  mtxPerVar(ntotElmt);
            Array<OneD, Array<OneD, NekDouble> > JacArray(ntotElmt);
            Array<OneD, int > elmtpnts(ntotElmt);
            Array<OneD, int > elmtcoef(ntotElmt);
            for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
            {
                nElmtCoef           = (*pexp)[nelmt]->GetNcoeffs();
                nElmtPnt            = (*pexp)[nelmt]->GetTotPoints();
                elmtpnts[nelmt]     =   nElmtPnt;
                elmtcoef[nelmt]     =   nElmtCoef;
                mtxPerVar[nelmt]    =MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nElmtCoef, nElmtPnt);
                JacArray[nelmt]    =Array<OneD, NekDouble>(nElmtPnt,0.0);
            }

            for(int m = 0; m < nConvectiveFields; m++)
            {
                for(int n = 0; n < nConvectiveFields; n++)
                {
                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        nElmtPnt            = elmtpnts[nelmt];
                        for(int npnt = 0; npnt < nElmtPnt; npnt++)
                        {
                            JacArray[nelmt][npnt]   =   (*(ElmtJac[nelmt][npnt]))(m,n);
                        }
                    }

                    // explist->GetMatIpwrtdbWeightBwd(JacArray,nDirctn,mtxPerVar);
                    explist->GetMatIpwrtDeriveBase(JacArray,nDirctn,mtxPerVar);

                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        nElmtCoef       = elmtcoef[nelmt];
                        nElmtPnt        = elmtpnts[nelmt];

                        tmpGmtx         = gmtxarray[m][n]->GetBlock(nelmt,nelmt);
                        ElmtMat         = mtxPerVar[nelmt];

                        for(int ncl = 0; ncl < nElmtPnt; ncl++)
                        {
                            for(int nrw = 0; nrw < nElmtCoef; nrw++)
                            {
                                tmp   =   (*tmpGmtx)(nrw,ncl) + (*ElmtMat)(nrw,ncl);
                                tmpGmtx->SetValue(nrw,ncl,tmp);
                            }
                        }
                    }
                }
            }
        }


        

        void AdvectionWeakDG::v_AddTraceJac2Mat(
            const int                                          nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const Array<OneD, DNekBlkMatSharedPtr>            &TraceJac,
            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray)
        {
            MultiRegions::ExpListSharedPtr explist = pFields[0]->GetTrace();
            std::shared_ptr<LocalRegions::ExpansionVector> pexp= explist->GetExp();
            int ntotElmt            = (*pexp).size();
            int nElmtPnt,nElmtCoef,noffset,pntoffset;

            NekDouble tmp;
            DNekMatSharedPtr        tmpGmtx,ElmtMat ;
            
            // Array<OneD,unsigned int> map;
            // Array<OneD,int> sign;

            Array<OneD, DNekMatSharedPtr>  mtxPerVar(ntotElmt);
            Array<OneD, Array<OneD, NekDouble> > JacFwd(ntotElmt);
            Array<OneD, Array<OneD, NekDouble> > JacBwd(ntotElmt);
            Array<OneD, int > elmtpnts(ntotElmt);
            Array<OneD, int > elmtcoef(ntotElmt);

            for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
            {
                nElmtCoef           = (*pexp)[nelmt]->GetNcoeffs();
                nElmtPnt            = (*pexp)[nelmt]->GetTotPoints();
                elmtpnts[nelmt]     =   nElmtPnt;
                elmtcoef[nelmt]     =   nElmtCoef;
                mtxPerVar[nelmt]    =MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nElmtCoef, nElmtPnt,0.0);
                JacFwd[nelmt]     =Array<OneD, NekDouble>(nElmtPnt,0.0);
                JacBwd[nelmt]     =Array<OneD, NekDouble>(nElmtPnt,0.0);
            }

            StdRegions::Orientation orient;
            Array<OneD, int > LAdjExpid(ntotElmt);
            Array<OneD, int > RAdjExpid(ntotElmt);
            Array<OneD, bool> RAdjflag(ntotElmt,false);

            Array<OneD, Array<OneD,unsigned int> > elmtLeftMap(ntotElmt);
            Array<OneD, Array<OneD,int> > elmtLeftSign(ntotElmt);
            Array<OneD, Array<OneD,unsigned int> > elmtRightMap(ntotElmt);
            Array<OneD, Array<OneD,int> > elmtRightSign(ntotElmt);
            // std::shared_ptr<LocalRegions::Expansion2DSharedPtr>    LAdjExp,RAdjExp;
            // int LAdjExpid, RAdjExpid;

            int ntmp=0;

            switch(explist->GetGraph()->GetSpaceDimension())
            {
                case 1:
                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        LocalRegions::Expansion0DSharedPtr traceEl =
                            explist->GetExp(nelmt)->as<LocalRegions::Expansion0D>();
                        LocalRegions::Expansion1DSharedPtr  LAdjExp      =   traceEl->GetLeftAdjacentElementExp();
                        LocalRegions::Expansion1DSharedPtr  RAdjExp      =   traceEl->GetRightAdjacentElementExp();

                        int LAdjBndid    =   traceEl->GetLeftAdjacentElementVertex();
                        int RAdjBndid    =   traceEl->GetRightAdjacentElementVertex();

                        orient  =   LAdjExp->v_GetEorient(LAdjBndid);
                        // LAdjExp->GetFaceToElementMap(LAdjBndid,orient,elmtLeftMap[nelmt],elmtLeftSign[nelmt]);
                        LAdjExpid[nelmt]    =   LAdjExp->GetElmtId();

                        if (RAdjBndid>-1) 
                        {
                            RAdjflag[nelmt] = true;
                            orient  =   RAdjExp->v_GetEorient(RAdjBndid);
                            // RAdjExp->GetFaceToElementMap(RAdjBndid,orient,elmtRightMap[nelmt],elmtRightSign[nelmt]);
                            RAdjExpid[nelmt]    =   RAdjExp->GetElmtId();
                        }
                    }
                    break;
                case 2:
                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        LocalRegions::Expansion1DSharedPtr traceEl =
                            explist->GetExp(nelmt)->as<LocalRegions::Expansion1D>();
                        LocalRegions::Expansion2DSharedPtr  LAdjExp      =   traceEl->GetLeftAdjacentElementExp();
                        LocalRegions::Expansion2DSharedPtr  RAdjExp      =   traceEl->GetRightAdjacentElementExp();

                        int LAdjBndid    =   traceEl->GetLeftAdjacentElementEdge();
                        int RAdjBndid    =   traceEl->GetRightAdjacentElementEdge();

                        ntmp    =   LAdjExp->GetElmtId();
                        LAdjExpid[nelmt]    =   ntmp;
                        orient  =   LAdjExp->v_GetEorient(LAdjBndid);
                        LAdjExp->GetEdgeToElementMap(LAdjBndid,orient,elmtLeftMap[nelmt],elmtLeftSign[nelmt]);
                        // std::cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<std::endl<<"Ntrace=    "<<nelmt<<std::endl<<"LAdjExp=    "<<ntmp<<endl;
                        // for(int i = 0; i < elmtLeftMap[nelmt].num_elements(); i++)
                        // {
                        //     std::cout<<"elmtLeftMap=    "<<elmtLeftMap[nelmt][i]<<" elmtLeftSign=    "<<elmtLeftSign[nelmt][i]<<endl;
                        // }
                        if (RAdjBndid>-1) 
                        {
                            RAdjflag[nelmt] = true;
                            ntmp    =   RAdjExp->GetElmtId();
                            RAdjExpid[nelmt]    =   ntmp;
                            orient  =   RAdjExp->v_GetEorient(RAdjBndid);
                            RAdjExp->GetEdgeToElementMap(RAdjBndid,orient,elmtRightMap[nelmt],elmtRightSign[nelmt]);
                            // std::cout<<"RAdjExp=    "<<ntmp<<endl;
                            // for(int i = 0; i < elmtRightMap[nelmt].num_elements(); i++)
                            // {
                            //     std::cout<<"elmtRightMap=    "<<elmtRightMap[nelmt][i]<<" elmtRightSign=    "<<elmtRightSign[nelmt][i]<<endl;
                            // }
                        }
                    }
                    break;
                case 3:
                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        LocalRegions::Expansion2DSharedPtr traceEl =
                            explist->GetExp(nelmt)->as<LocalRegions::Expansion2D>();
                        LocalRegions::Expansion3DSharedPtr  LAdjExp      =   traceEl->GetLeftAdjacentElementExp();
                        LocalRegions::Expansion3DSharedPtr  RAdjExp      =   traceEl->GetRightAdjacentElementExp();

                        int LAdjBndid    =   traceEl->GetLeftAdjacentElementFace();
                        int RAdjBndid    =   traceEl->GetRightAdjacentElementFace();

                        orient  =   LAdjExp->v_GetEorient(LAdjBndid);
                        LAdjExp->GetFaceToElementMap(LAdjBndid,orient,elmtLeftMap[nelmt],elmtLeftSign[nelmt]);
                        LAdjExpid[nelmt]    =   LAdjExp->GetElmtId();

                        if (RAdjBndid>-1) 
                        {
                            RAdjflag[nelmt] = true;
                            orient  =   RAdjExp->v_GetEorient(RAdjBndid);
                            RAdjExp->GetFaceToElementMap(RAdjBndid,orient,elmtRightMap[nelmt],elmtRightSign[nelmt]);
                            RAdjExpid[nelmt]    =   RAdjExp->GetElmtId();
                        }
                    }
                    break;
                default:
                    ASSERTL0(false,"Trace SpaceDimension>2")
                    break;
            }

            for(int m = 0; m < nConvectiveFields; m++)
            {
                for(int n = 0; n < nConvectiveFields; n++)
                {
                    // std::cout<<std::endl<<"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<std::endl<<"m=    "<<m<<"  n=    "<<n<<endl;
                    
                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        nElmtPnt            = elmtpnts[nelmt];
                        noffset             =   explist->GetPhys_Offset(nelmt);
                        for(int npnt = 0; npnt < nElmtPnt; npnt++)
                        {
                            pntoffset = noffset+npnt;
                            // cout<<"pntoffset="<<pntoffset<<endl;
                            JacFwd[nelmt][npnt]   =   (*(TraceJac[0]->GetBlock(pntoffset,pntoffset)))(m,n);
                            JacBwd[nelmt][npnt]   =   (*(TraceJac[1]->GetBlock(pntoffset,pntoffset)))(m,n);
                        }
                    }

                    // explist->GetMatIpwrtbWeightBwd(JacFwd,mtxPerVar);
                    explist->GetMatIpwrtBase(JacFwd,mtxPerVar);

                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        nElmtCoef       = elmtcoef[nelmt];
                        nElmtPnt        = elmtpnts[nelmt];

                        ElmtMat         = mtxPerVar[nelmt];
                        // std::cout   <<std::endl<<"*********************************"<<std::endl<<"element :   "<<nelmt<<std::endl;
                        // std::cout   <<(*ElmtMat)<<endl;
                       
                        tmpGmtx         = gmtxarray[m][n]->GetBlock(LAdjExpid[nelmt],LAdjExpid[nelmt]);

                        for(int ncl = 0; ncl < nElmtPnt; ncl++)
                        {
                            int nclAdjExp = elmtLeftMap[nelmt][ncl];

                            for(int nrw = 0; nrw < nElmtCoef; nrw++)
                            {
                                int nrwAdjExp = elmtLeftMap[nelmt][nrw];
                                tmp   =   (*tmpGmtx)(nclAdjExp,nrwAdjExp);
                                tmp   +=  elmtLeftSign[nelmt][ncl]*elmtLeftSign[nelmt][nrw]*(*ElmtMat)(nrw,ncl);
                                tmpGmtx->SetValue(nclAdjExp,nrwAdjExp,tmp);
                            }
                        }
                    }
                    
                    // explist->GetMatIpwrtbWeightBwd(JacBwd,mtxPerVar);
                    explist->GetMatIpwrtBase(JacBwd,mtxPerVar);

                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        if(RAdjflag[nelmt])
                        {
                            nElmtCoef       = elmtcoef[nelmt];
                            nElmtPnt        = elmtpnts[nelmt];

                            ElmtMat         = mtxPerVar[nelmt];

                            // std::cout   <<std::endl<<"*********************************"<<std::endl<<"element :   "<<nelmt<<std::endl;
                            // std::cout   <<(*ElmtMat)<<endl;

                            tmpGmtx         = gmtxarray[m][n]->GetBlock(RAdjExpid[nelmt],RAdjExpid[nelmt]);
                            for(int ncl = 0; ncl < nElmtPnt; ncl++)
                            {
                                int nclAdjExp = elmtRightMap[nelmt][ncl];

                                for(int nrw = 0; nrw < nElmtCoef; nrw++)
                                {
                                    int nrwAdjExp = elmtRightMap[nelmt][nrw];
                                    tmp   =   (*tmpGmtx)(nclAdjExp,nrwAdjExp);
                                    tmp   -=  elmtRightSign[nelmt][ncl]*elmtRightSign[nelmt][nrw]*(*ElmtMat)(nrw,ncl);
                                    tmpGmtx->SetValue(nclAdjExp,nrwAdjExp,tmp);
                                }
                            }
                        }
                    }

                    
                }
            }
        }
    }//end of namespace SolverUtils
}//end of namespace Nektar
