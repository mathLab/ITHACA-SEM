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
#include <MultiRegions/AssemblyMap/LocTraceToTraceMap.h>
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
            
            // NumericalFlux(nConvectiveFields, fields, advVel, inarray, numflux,time,pFwd,pBwd);

            // Evaulate <\phi, \hat{F}\cdot n> - OutField[i]
            for(i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Neg                      (nCoeffs, tmp[i], 1);
                fields[i]->AddTraceIntegral     (numflux[i], tmp[i]);
                fields[i]->MultiplyByElmtInvMass(tmp[i], tmp[i]);
                fields[i]->BwdTrans             (tmp[i], outarray[i]);
            }
        }


        void AdvectionWeakDG::v_Advect_coeff(
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
            for(int i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Fill(outarray[i].num_elements(),0.0,outarray[i],1);
                fields[i]->IProductWRTDerivBase(fluxvector[i],outarray[i]);
            }

            // Store forwards/backwards space along trace space
            Array<OneD, Array<OneD, NekDouble> > numflux(nConvectiveFields);
            for(i = 0; i < nConvectiveFields; ++i)
            {
                numflux[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
            }

            NumericalFlux(nConvectiveFields, fields, advVel, inarray, numflux,time,pFwd,pBwd);

            // Evaulate <\phi, \hat{F}\cdot n> - OutField[i]
            for(i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Neg                      (nCoeffs, outarray[i], 1);
                fields[i]->AddTraceIntegral     (numflux[i], outarray[i]);
                fields[i]->MultiplyByElmtInvMass(outarray[i], outarray[i]);
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

                    explist->GetMatIpwrtdbWeightBwd(JacArray,nDirctn,mtxPerVar);
                    // explist->GetMatIpwrtDeriveBase(JacArray,nDirctn,mtxPerVar);

                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        nElmtCoef       = elmtcoef[nelmt];
                        nElmtPnt        = elmtpnts[nelmt];

                        tmpGmtx         = gmtxarray[m][n]->GetBlock(nelmt,nelmt);
                        ElmtMat         = mtxPerVar[nelmt];

                        for(int ncl = 0; ncl < nElmtCoef; ncl++)
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
            MultiRegions::ExpListSharedPtr tracelist = pFields[0]->GetTrace();
            std::shared_ptr<LocalRegions::ExpansionVector> traceExp= tracelist->GetExp();
            int ntotTrac            = (*traceExp).size();
            int nTracPnt,nTracCoef,noffset,pntoffset;


            MultiRegions::ExpListSharedPtr fieldlist = pFields[0]->GetTrace();
            std::shared_ptr<LocalRegions::ExpansionVector> fieldExp= tracelist->GetExp();
            int ntotElmt            = (*traceExp).size();
            int nElmtPnt,nElmtCoef,nElmtoffset,Elmtpntoffset;

            NekDouble tmp;
            DNekMatSharedPtr        tmpGmtx,ElmtMat ;
            
            Array<OneD, DNekMatSharedPtr>  mtxPerVar(ntotTrac);
            Array<OneD, Array<OneD, NekDouble> > JacFwd(ntotTrac);
            Array<OneD, Array<OneD, NekDouble> > JacBwd(ntotTrac);
            Array<OneD, int > elmtpnts(ntotTrac);
            Array<OneD, int > elmtcoef(ntotTrac);

            for(int  nelmt = 0; nelmt < ntotTrac; nelmt++)
            {
                nTracCoef           = (*traceExp)[nelmt]->GetNcoeffs();
                nTracPnt            = (*traceExp)[nelmt]->GetTotPoints();
                elmtpnts[nelmt]     =   nTracPnt;
                elmtcoef[nelmt]     =   nTracCoef;
                // mtxPerVar[nelmt]    =MemoryManager<DNekMat>
                //                     ::AllocateSharedPtr(nTracCoef, nTracPnt,0.0);
                mtxPerVar[nelmt]    =MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nTracCoef, nTracCoef,0.0);
                JacFwd[nelmt]     =Array<OneD, NekDouble>(nTracPnt,0.0);
                JacBwd[nelmt]     =Array<OneD, NekDouble>(nTracPnt,0.0);
            }

            StdRegions::Orientation orient;
            Array<OneD, int > LAdjExpid(ntotTrac);
            Array<OneD, int > RAdjExpid(ntotTrac);
            Array<OneD, bool> LAdjflag(ntotTrac,false);
            Array<OneD, bool> RAdjflag(ntotTrac,false);

            Array<OneD, Array<OneD,unsigned int > > elmtLeftMap(ntotTrac);
            Array<OneD, Array<OneD,int          > > elmtLeftSign(ntotTrac);
            Array<OneD, Array<OneD,unsigned int > > elmtRightMap(ntotTrac);
            Array<OneD, Array<OneD,int          > > elmtRightSign(ntotTrac);

            // const Array<OneD, const Array<OneD, int >> LRAdjExpid;
            // const Array<OneD, const Array<OneD, bool>> LRAdjflag;

            // const Array<OneD, const Array<OneD, Array<OneD, int > > > elmtLRMap;
            // const Array<OneD, const Array<OneD, Array<OneD, int > > > elmtLRSign;

            int ntmp=0;

            switch(tracelist->GetGraph()->GetSpaceDimension())
            {
                case 1:
                {
                    ASSERTL0(false,"GetSpaceDimension==1 not coded yet");
                    break;
                }
                case 2:
                case 3:
                {
                    const MultiRegions::LocTraceToTraceMapSharedPtr locTraceToTraceMap = pFields[0]->GetlocTraceToTraceMap();
                    
                    const Array<OneD, const Array<OneD, int >> LRAdjExpid  =   locTraceToTraceMap->GetLeftRightAdjacentExpId();
                    // for(int lr=0;lr<2;lr++)
                    // {

                    // }
                    const Array<OneD, const Array<OneD, bool>> LRAdjflag   =   locTraceToTraceMap->GetLeftRightAdjacentExpFlag();

                    const Array<OneD, const Array<OneD, Array<OneD, int > > > elmtLRMap   =   locTraceToTraceMap->GetTraceceffToLeftRightExpcoeffMap();
                    const Array<OneD, const Array<OneD, Array<OneD, int > > > elmtLRSign  =   locTraceToTraceMap->GetTraceceffToLeftRightExpcoeffSign();

                    for(int  nelmt = 0; nelmt < ntotTrac; nelmt++)
                    {
                        LocalRegions::Expansion1DSharedPtr traceEl =
                            tracelist->GetExp(nelmt)->as<LocalRegions::Expansion1D>();
                        LocalRegions::Expansion2DSharedPtr  LAdjExp      =   traceEl->GetLeftAdjacentElementExp();
                        LocalRegions::Expansion2DSharedPtr  RAdjExp      =   traceEl->GetRightAdjacentElementExp();

                        int ntracecoeff     =   traceEl->GetNcoeffs();
                        elmtLeftMap     [nelmt] =   Array<OneD,unsigned int > (ntracecoeff);
                        elmtLeftSign    [nelmt] =   Array<OneD,int          > (ntracecoeff);
                        elmtRightMap    [nelmt] =   Array<OneD,unsigned int > (ntracecoeff);
                        elmtRightSign   [nelmt] =   Array<OneD,int          > (ntracecoeff);
                    
                        LAdjExpid[nelmt]    =   LRAdjExpid[0][nelmt];
                        RAdjExpid[nelmt]    =   LRAdjExpid[1][nelmt];

                        LAdjflag[nelmt]    =   LRAdjflag[0][nelmt];
                        RAdjflag[nelmt]    =   LRAdjflag[1][nelmt];

                        for(int i = 0; i < elmtLRMap[0][nelmt].num_elements(); i++)
                        {
                            elmtLeftMap[nelmt][i]   =   elmtLRMap[0][nelmt][i];
                            elmtRightMap[nelmt][i]  =   elmtLRMap[1][nelmt][i];
                            elmtLeftSign[nelmt][i]   =   elmtLRSign[0][nelmt][i];
                            elmtRightSign[nelmt][i]  =   elmtLRSign[1][nelmt][i];
                        }
                    }
                    break;
                }
                default:
                {
                    ASSERTL0(false,"Trace SpaceDimension>2");
                    break;
                }
            }

            DNekMatSharedPtr FtmpMat,BtmpMat;
            for(int m = 0; m < nConvectiveFields; m++)
            {
                for(int n = 0; n < nConvectiveFields; n++)
                {
                    for(int  nelmt = 0; nelmt < ntotTrac; nelmt++)
                    {
                        nTracPnt            = elmtpnts[nelmt];
                        noffset             =   tracelist->GetPhys_Offset(nelmt);
                        for(int npnt = 0; npnt < nTracPnt; npnt++)
                        {
                            NekDouble ftmp,btmp;
                            pntoffset = noffset+npnt;
                            // cout<<"pntoffset="<<pntoffset<<endl;
                            FtmpMat                 = TraceJac[0]->GetBlock(pntoffset,pntoffset);
                            ftmp                    =   (*FtmpMat)(m,n);
                            JacFwd[nelmt][npnt]     =   ftmp;

                            BtmpMat                 = TraceJac[1]->GetBlock(pntoffset,pntoffset);
                            btmp                    =   (*BtmpMat)(m,n);
                            JacBwd[nelmt][npnt]     =   btmp;
                        }
                    }

                    tracelist->GetMatIpwrtbWeightBwd(JacFwd,mtxPerVar);
                    // tracelist->GetMatIpwrtBase(JacFwd,mtxPerVar);

                    for(int  nelmt = 0; nelmt < ntotTrac; nelmt++)
                    {
                        if(LAdjflag[nelmt])
                        {
                            nTracCoef       = elmtcoef[nelmt];
                            nTracPnt        = elmtpnts[nelmt];

                            ElmtMat         = mtxPerVar[nelmt];
                        
                            tmpGmtx         = gmtxarray[m][n]->GetBlock(LAdjExpid[nelmt],LAdjExpid[nelmt]);

                            for(int ncl = 0; ncl < nTracCoef; ncl++)
                            {
                                int nclAdjExp = elmtLeftMap[nelmt][ncl];

                                for(int nrw = 0; nrw < nTracCoef; nrw++)
                                {
                                    int nrwAdjExp = elmtLeftMap[nelmt][nrw];
                                    tmp   =   (*tmpGmtx)(nrwAdjExp,nclAdjExp);
                                    // tmp   -=  elmtLeftSign[nelmt][ncl]*elmtLeftSign[nelmt][nrw]*(*ElmtMat)(nrw,ncl);
                                    tmp   -=  (*ElmtMat)(nrw,ncl);
                                    tmpGmtx->SetValue(nrwAdjExp,nclAdjExp,tmp);
                                }
                            }
                        }
                    }

                                        
                    
                    tracelist->GetMatIpwrtbWeightBwd(JacBwd,mtxPerVar);
                    // tracelist->GetMatIpwrtBase(JacBwd,mtxPerVar);

                    for(int  nelmt = 0; nelmt < ntotTrac; nelmt++)
                    {
                        if(RAdjflag[nelmt])
                        {
                            nTracCoef       = elmtcoef[nelmt];
                            nTracPnt        = elmtpnts[nelmt];

                            ElmtMat         = mtxPerVar[nelmt];

                            tmpGmtx         = gmtxarray[m][n]->GetBlock(RAdjExpid[nelmt],RAdjExpid[nelmt]);

                            for(int ncl = 0; ncl < nTracCoef; ncl++)
                            {
                                int nclAdjExp = elmtRightMap[nelmt][ncl];

                                for(int nrw = 0; nrw < nTracCoef; nrw++)
                                {
                                    int nrwAdjExp = elmtRightMap[nelmt][nrw];
                                    tmp   =   (*tmpGmtx)(nrwAdjExp,nclAdjExp);
                                    // tmp   +=  elmtRightSign[nelmt][ncl]*elmtRightSign[nelmt][nrw]*(*ElmtMat)(nrw,ncl);
                                    tmp   +=  (*ElmtMat)(nrw,ncl);
                                    tmpGmtx->SetValue(nrwAdjExp,nclAdjExp,tmp);
                                }
                            }
                        }
                    }
                }
            }
        }


        void AdvectionWeakDG::v_AddTraceJac2Mat_new(
            const int                                           nConvectiveFields,
            const int                                           nSpaceDim,
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &pFields,
            const Array<OneD, DNekBlkMatSharedPtr>              &TraceJacCons,
            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >      &gmtxarray,
            const Array<OneD, DNekBlkMatSharedPtr>              &TraceJacGrad)
        {
            MultiRegions::ExpListSharedPtr tracelist = pFields[0]->GetTrace();
            std::shared_ptr<LocalRegions::ExpansionVector> traceExp= tracelist->GetExp();
            int ntotTrac            = (*traceExp).size();
            int nTracPnt,nTracCoef,noffset,pntoffset;

            NekDouble tmp;
            DNekMatSharedPtr                    tmpGmtx,ElmtMat ;
            Array<OneD, DNekBlkMatSharedPtr>    TraceJac;
            Array<OneD, DNekMatSharedPtr>       mtxPerVar;


            Array<OneD, DNekMatSharedPtr>  TraceJacFwd(ntotTrac);
            Array<OneD, DNekMatSharedPtr>  TraceJacBwd(ntotTrac);
            Array<OneD, Array<OneD, NekDouble> > JacFwd(ntotTrac);
            Array<OneD, Array<OneD, NekDouble> > JacBwd(ntotTrac);
            Array<OneD, int > elmtpnts(ntotTrac);
            Array<OneD, int > elmtcoef(ntotTrac);

            for(int  nelmt = 0; nelmt < ntotTrac; nelmt++)
            {
                nTracCoef           = (*traceExp)[nelmt]->GetNcoeffs();
                nTracPnt            = (*traceExp)[nelmt]->GetTotPoints();
                elmtpnts[nelmt]     =   nTracPnt;
                elmtcoef[nelmt]     =   nTracCoef;
                TraceJacFwd[nelmt]    =MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nTracCoef, nTracPnt,0.0);

                TraceJacBwd[nelmt]    =MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nTracCoef, nTracPnt,0.0);

                JacFwd[nelmt]     =Array<OneD, NekDouble>(nTracPnt,0.0);
                JacBwd[nelmt]     =Array<OneD, NekDouble>(nTracPnt,0.0);
            }

            StdRegions::Orientation orient;
            Array<OneD, int > LAdjExpid(ntotTrac);
            Array<OneD, int > RAdjExpid(ntotTrac);
            Array<OneD, bool> LAdjflag(ntotTrac,false);
            Array<OneD, bool> RAdjflag(ntotTrac,false);

            Array<OneD, Array<OneD,unsigned int > > elmtLeftMap(ntotTrac);
            Array<OneD, Array<OneD,int          > > elmtLeftSign(ntotTrac);
            Array<OneD, Array<OneD,unsigned int > > elmtRightMap(ntotTrac);
            Array<OneD, Array<OneD,int          > > elmtRightSign(ntotTrac);

            // const Array<OneD, const Array<OneD, int >> LRAdjExpid;
            // const Array<OneD, const Array<OneD, bool>> LRAdjflag;

            // const Array<OneD, const Array<OneD, Array<OneD, int > > > elmtLRMap;
            // const Array<OneD, const Array<OneD, Array<OneD, int > > > elmtLRSign;

            int ntmp=0;

            switch(tracelist->GetGraph()->GetSpaceDimension())
            {
                case 1:
                {
                    ASSERTL0(false,"GetSpaceDimension==1 not coded yet");
                    break;
                }
                case 2:
                case 3:
                {
                    const MultiRegions::LocTraceToTraceMapSharedPtr locTraceToTraceMap = pFields[0]->GetlocTraceToTraceMap();
                    
                    const Array<OneD, const Array<OneD, int >> LRAdjExpid  =   locTraceToTraceMap->GetLeftRightAdjacentExpId();
                    // for(int lr=0;lr<2;lr++)
                    // {

                    // }
                    const Array<OneD, const Array<OneD, bool>> LRAdjflag   =   locTraceToTraceMap->GetLeftRightAdjacentExpFlag();

                    const Array<OneD, const Array<OneD, Array<OneD, int > > > elmtLRMap   =   locTraceToTraceMap->GetTraceceffToLeftRightExpcoeffMap();
                    const Array<OneD, const Array<OneD, Array<OneD, int > > > elmtLRSign  =   locTraceToTraceMap->GetTraceceffToLeftRightExpcoeffSign();

                    for(int  nelmt = 0; nelmt < ntotTrac; nelmt++)
                    {
                        LocalRegions::Expansion1DSharedPtr traceEl =
                            tracelist->GetExp(nelmt)->as<LocalRegions::Expansion1D>();
                        LocalRegions::Expansion2DSharedPtr  LAdjExp      =   traceEl->GetLeftAdjacentElementExp();
                        LocalRegions::Expansion2DSharedPtr  RAdjExp      =   traceEl->GetRightAdjacentElementExp();

                        int ntracecoeff     =   traceEl->GetNcoeffs();
                        elmtLeftMap     [nelmt] =   Array<OneD,unsigned int > (ntracecoeff);
                        elmtLeftSign    [nelmt] =   Array<OneD,int          > (ntracecoeff);
                        elmtRightMap    [nelmt] =   Array<OneD,unsigned int > (ntracecoeff);
                        elmtRightSign   [nelmt] =   Array<OneD,int          > (ntracecoeff);
                    
                        LAdjExpid[nelmt]    =   LRAdjExpid[0][nelmt];
                        RAdjExpid[nelmt]    =   LRAdjExpid[1][nelmt];

                        LAdjflag[nelmt]    =   LRAdjflag[0][nelmt];
                        RAdjflag[nelmt]    =   LRAdjflag[1][nelmt];

                        for(int i = 0; i < elmtLRMap[0][nelmt].num_elements(); i++)
                        {
                            elmtLeftMap[nelmt][i]   =   elmtLRMap[0][nelmt][i];
                            elmtRightMap[nelmt][i]  =   elmtLRMap[1][nelmt][i];
                            elmtLeftSign[nelmt][i]   =   elmtLRSign[0][nelmt][i];
                            elmtRightSign[nelmt][i]  =   elmtLRSign[1][nelmt][i];
                        }
                    }
                    break;
                }
                default:
                {
                    ASSERTL0(false,"Trace SpaceDimension>2");
                    break;
                }
            }


            MultiRegions::ExpListSharedPtr fieldlist = pFields[0]->GetTrace();
            std::shared_ptr<LocalRegions::ExpansionVector> fieldExp= tracelist->GetExp();
            int ntotElmt            = (*traceExp).size();
            int nElmtPnt,nElmtCoef,nElmtoffset,Elmtpntoffset;

            Array<OneD, DNekMatSharedPtr>  ElmtJacGrad(ntotElmt);
            Array<OneD, DNekMatSharedPtr>  ElmtJacCons(ntotElmt);

            for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
            {
                nElmtPnt            = (*fieldExp)[nelmt]->GetNcoeffs();
                nElmtCoef           = (*fieldExp)[nelmt]->GetTotPoints();
                
                ElmtJacGrad[nelmt]    =MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nElmtCoef, nElmtPnt,0.0);
                ElmtJacCons[nelmt]    =MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nElmtCoef, nElmtPnt,0.0);
            }


            for(int m = 0; m < nConvectiveFields; m++)
            {
                for(int n = 0; n < nConvectiveFields; n++)
                {
                    TraceJac = TraceJacGrad;
                    // ElmtJacCons to set 0
                    // for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    // {
                    //     (*ElmtJacCons[nelmt]) =  0.0;
                    // }
                    for(int ndir = 0; ndir < nSpaceDim; ndir++)
                    {
                        // ElmtJacGrad to set 0
                        // for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                        // {
                        //     (*ElmtJacGrad[nelmt]) =  0.0;
                        // }
                        int ngrad = n*nSpaceDim+ndir;
                        
                        CalcJacobTraceInteg(pFields,m,ngrad,TraceJac,TraceJacFwd,TraceJacBwd);
                        
                        
                        pFields[0]->AddTraceJacToElmtJac(TraceJacFwd,TraceJacBwd,ElmtJacGrad);


                        // pFields[0]->RightIPTPhysDeriv(ndir,ElmtJacGrad,ElmtJacGrad);

                        for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                        {
                            (*ElmtJacCons[nelmt]) = (*ElmtJacCons[nelmt])+ (*ElmtJacGrad[nelmt]); 
                        }

                    }


                    TraceJac = TraceJacCons;

                    CalcJacobTraceInteg(pFields,m,n,TraceJac,TraceJacFwd,TraceJacBwd);
                        
                        
                    pFields[0]->AddTraceJacToElmtJac(TraceJacFwd, TraceJacBwd,ElmtJacCons);


                    // pFields[0]->RightIPBwdMatrix(ndir,ElmtJacCons,ElmtJacCons);


                    // add ElmtJacCons to gmtxarray[m][n]


                }
            }
        }


        void AdvectionWeakDG::CalcJacobTraceInteg(
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &pFields,
            const int                                         m,
            const int                                         n,
            const Array<OneD, const DNekBlkMatSharedPtr>    & TraceJac,
            Array<OneD, DNekMatSharedPtr>                   & TraceJacFwd,
            Array<OneD, DNekMatSharedPtr>                   & TraceJacBwd)
        {
            MultiRegions::ExpListSharedPtr tracelist = pFields[0]->GetTrace();
            std::shared_ptr<LocalRegions::ExpansionVector> traceExp= tracelist->GetExp();
            int ntotTrac            = (*traceExp).size();
            int nTracPnt,nTracCoef,noffset,pntoffset;

            Array<OneD, int > elmtpnts(ntotTrac);

            Array<OneD, Array<OneD, NekDouble> > JacFwd(ntotTrac);
            Array<OneD, Array<OneD, NekDouble> > JacBwd(ntotTrac);


            for(int  nelmt = 0; nelmt < ntotTrac; nelmt++)
            {
                nTracCoef           = (*traceExp)[nelmt]->GetNcoeffs();
                nTracPnt            = 
                elmtpnts[nelmt]     =   nTracPnt;

                JacFwd[nelmt]     =Array<OneD, NekDouble>(nTracPnt,0.0);
                JacBwd[nelmt]     =Array<OneD, NekDouble>(nTracPnt,0.0);
            }



            DNekMatSharedPtr FtmpMat,BtmpMat;
            NekDouble ftmp,btmp;
            for(int  nelmt = 0; nelmt < ntotTrac; nelmt++)
            {
                nTracPnt            =   elmtpnts[nelmt];
                noffset             =   tracelist->GetPhys_Offset(nelmt);
                for(int npnt = 0; npnt < nTracPnt; npnt++)
                {
                    pntoffset = noffset+npnt;
                    FtmpMat                 = TraceJac[0]->GetBlock(pntoffset,pntoffset);
                    ftmp                    =   (*FtmpMat)(m,n);
                    JacFwd[nelmt][npnt]     =   ftmp;

                    BtmpMat                 = TraceJac[1]->GetBlock(pntoffset,pntoffset);
                    btmp                    =   (*BtmpMat)(m,n);
                    JacBwd[nelmt][npnt]     =   btmp;
                }
            }

            // tracelist->GetMatIpwrtbWeightBwd(JacFwd,mtxPerVar);
            tracelist->GetMatIpwrtBase(JacFwd,TraceJacFwd);
            
            // tracelist->GetMatIpwrtbWeightBwd(JacBwd,mtxPerVar);
            tracelist->GetMatIpwrtBase(JacBwd,TraceJacBwd);

        }

        void AdvectionWeakDG::v_NumCalRiemFluxJac(
            const int                                          nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &AdvVel,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd,
            DNekBlkMatSharedPtr &FJac,
            DNekBlkMatSharedPtr &BJac)
        {
            int nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();

            ASSERTL1(m_riemann,
                     "Riemann solver must be provided for AdvectionWeakDG.");

            // Store forwards/backwards space along trace space
            Array<OneD, Array<OneD, NekDouble> > Fwd    (nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> > Bwd    (nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> > flux   (nConvectiveFields);

            if (pFwd == NullNekDoubleArrayofArray ||
                pBwd == NullNekDoubleArrayofArray)
            {
                for(int i = 0; i < nConvectiveFields; ++i)
                {
                    Fwd[i]      = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                    Bwd[i]      = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                    flux[i]     = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                    fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
                }
            }
            else
            {
                for(int i = 0; i < nConvectiveFields; ++i)
                {
                    Fwd[i]      = pFwd[i];
                    Bwd[i]      = pBwd[i];
                    flux[i]     = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                }
            }

            m_riemann->Solve(m_spaceDim, Fwd, Bwd, flux);

            NekDouble eps   =   1.0E-5;

            int nFields = Fwd   .num_elements();
            int nPts    = Fwd[0].num_elements();


            Array<OneD,       Array<OneD, NekDouble> >  plusflux(nFields),plusfield(nFields);
            Array<OneD,       Array<OneD, NekDouble> >  Jacvect(nFields);
            
            // estimate the magnitude of each flow variables 
            Array<OneD, NekDouble>  magnitdEstimat(nFields,0.0);
            for(int i = 0; i < nFields; i++)
            {
                for(int j=0;j<nPts;j++)
                {
                    magnitdEstimat[i]   +=   Fwd[i][j]*Fwd[i][j] ;
                }
            }
            for(int i = 2; i < nFields-1; i++)
            {
                magnitdEstimat[1]   +=   magnitdEstimat[i] ;
            }
            for(int i = 2; i < nFields-1; i++)
            {
                magnitdEstimat[i]   =   magnitdEstimat[1] ;
            }
            NekDouble ototpnts = 1.0/(nPts);
            for(int i = 0; i < nFields; i++)
            {
                magnitdEstimat[i] = sqrt(magnitdEstimat[i]*ototpnts);
            }

            // Allocate the Jacobian matrix
            for(int i=0;i<nPts;i++)
            {
                DNekMatSharedPtr    tmpMat;
                tmpMat = MemoryManager<DNekMat>
                        ::AllocateSharedPtr(nFields, nFields);
                FJac->SetBlock(i,i,tmpMat);

                tmpMat = MemoryManager<DNekMat>
                        ::AllocateSharedPtr(nFields, nFields);
                BJac->SetBlock(i,i,tmpMat);
            }

            // Allocate temporary variables
            for(int i = 0; i < nFields; i++)
            {
                Jacvect[i]      =    Array<OneD, NekDouble>(nPts,0.0);
                plusflux[i]     =    Array<OneD, NekDouble>(nPts,0.0);
                plusfield[i]    =    Array<OneD, NekDouble>(nPts,0.0);
            }

            DNekMatSharedPtr tmpMat;


            for(int i = 0; i < nFields; i++)
            {
                Vmath::Vcopy(nPts, Fwd[i],1,plusfield[i],1);
            }

            // Fwd Jacobian
            for(int i = 0; i < nFields; i++)
            {
                NekDouble epsvar = eps*magnitdEstimat[i];
                NekDouble oepsvar   =   1.0/epsvar;
                for(int j = 0; j < nPts; j++)
                {
                    plusfield[i][j] =    Fwd[i][j]  +   epsvar;
                }


                

                m_riemann->Solve(m_spaceDim, plusfield, Bwd, plusflux);

                for (int n = 0; n < nFields; n++)
                {
                    Vmath::Vsub(nPts,&plusflux[n][0],1,&flux[n][0],1,&Jacvect[n][0],1);
                    Vmath::Smul(nPts, oepsvar ,&Jacvect[n][0],1,&Jacvect[n][0],1);
                }
                for(int j = 0; j < nPts; j++)
                {
                    tmpMat  =   FJac->GetBlock(j,j);
                    for (int n = 0; n < nFields; n++)
                    {
                        (*tmpMat)(n,i) = Jacvect[n][j];
                    }
                }
                
                for(int j = 0; j < nPts; j++)
                {
                    plusfield[i][j] =    Fwd[i][j];
                }
            }

            for(int i = 0; i < nFields; i++)
            {
                Vmath::Vcopy(nPts, Bwd[i],1,plusfield[i],1);
            }

            for(int i = 0; i < nFields; i++)
            {
                NekDouble epsvar = eps*magnitdEstimat[i];
                NekDouble oepsvar   =   1.0/epsvar;

                for(int j = 0; j < nPts; j++)
                {
                    plusfield[i][j] =    Bwd[i][j]  +   epsvar;
                }

                m_riemann->Solve(m_spaceDim, Fwd, plusfield, plusflux);

                for (int n = 0; n < nFields; n++)
                {
                    Vmath::Vsub(nPts,&plusflux[n][0],1,&flux[n][0],1,&Jacvect[n][0],1);
                    Vmath::Smul(nPts, oepsvar ,&Jacvect[n][0],1,&Jacvect[n][0],1);
                }
                for(int j = 0; j < nPts; j++)
                {
                    tmpMat  =   BJac->GetBlock(j,j);
                    for (int n = 0; n < nFields; n++)
                    {
                        (*tmpMat)(n,i) = Jacvect[n][j];
                    }
                }
                
                for(int j = 0; j < nPts; j++)
                {
                    plusfield[i][j] =    Bwd[i][j];
                }
            }
        }
    }//end of namespace SolverUtils
}//end of namespace Nektar
