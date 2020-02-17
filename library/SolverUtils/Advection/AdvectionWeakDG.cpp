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
#include <iostream>
#include <iomanip>

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

#ifdef CFS_DEBUGMODE
        pSession->LoadParameter("DebugVolTraceSwitch",                 m_DebugVolTraceSwitch      ,    0);
#endif
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
            int nCoeffs         = fields[0]->GetNcoeffs();
            
            Array<OneD, Array<OneD, NekDouble> > tmp(nConvectiveFields);
            for(int i = 0; i < nConvectiveFields; ++i)
            {
                tmp[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
            }

            AdvectionWeakDG::v_Advect_coeff(nConvectiveFields,fields,advVel,inarray,tmp,time,pFwd,pBwd);
            
            for(int i = 0; i < nConvectiveFields; ++i)
            {
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
            
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > fluxvector(nConvectiveFields);
            // Allocate storage for flux vector F(u).
            for (i = 0; i < nConvectiveFields; ++i)
            {
                fluxvector[i] =
                    Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
                for (j = 0; j < m_spaceDim; ++j)
                {
                    fluxvector[i][j] = Array<OneD, NekDouble>(nPointsTot,0.0);
                }
            }
#ifdef CFS_DEBUGMODE
        if(2!=m_DebugVolTraceSwitch)
        {
#endif
            v_AdvectVolumeFlux(nConvectiveFields,fields,advVel,inarray,fluxvector,time);
#ifdef CFS_DEBUGMODE
        }
#endif
            
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
#ifdef CFS_DEBUGMODE
        if(1!=m_DebugVolTraceSwitch)
        {
#endif
            v_AdvectTraceFlux(nConvectiveFields, fields, advVel, inarray, numflux,time,pFwd,pBwd);
#ifdef CFS_DEBUGMODE
        }
#endif
            // Evaulate <\phi, \hat{F}\cdot n> - OutField[i]
            for(i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Neg                      (nCoeffs, outarray[i], 1);
                fields[i]->AddTraceIntegral     (numflux[i], outarray[i]);
                fields[i]->MultiplyByElmtInvMass(outarray[i], outarray[i]);
            }
        }

           /**
         * @brief Compute the advection term at each time-step using the
         * Discontinuous Galerkin approach (DG).
         *
         * @param nConvectiveFields   Number of fields.
         * @param fields              Pointer to fields.
         * @param advVel              Advection velocities.
         * @param inarray             Solution at the previous time-step.
         * @param TraceFlux         Advection Trace flux
         *                            time integration class.
         */
        void AdvectionWeakDG::v_AdvectTraceFlux(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble>>         &advVel,
            const Array<OneD, Array<OneD, NekDouble>>         &inarray,
                  Array<OneD, Array<OneD, NekDouble>>         &TraceFlux, 
            const NekDouble                                   &time,
            const Array<OneD, Array<OneD, NekDouble>>         &pFwd,
            const Array<OneD, Array<OneD, NekDouble>>         &pBwd)
        {
            int nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();

            ASSERTL1(m_riemann, "Riemann solver must be provided for AdvectionWeakDG.");

            // Store forwards/backwards space along trace space
            Array<OneD, Array<OneD, NekDouble>> Fwd(nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble>> Bwd(nConvectiveFields);
            // Array<OneD, Array<OneD, NekDouble> > numflux(nConvectiveFields);

            if (pFwd == NullNekDoubleArrayofArray || pBwd == NullNekDoubleArrayofArray)
            {
                for (int i = 0; i < nConvectiveFields; ++i)
                {
                    Fwd[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                    Bwd[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                    // numflux[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                    fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
                }
            }
            else
            {
                for (int i = 0; i < nConvectiveFields; ++i)
                {
                    Fwd[i] = pFwd[i];
                    Bwd[i] = pBwd[i];
                    // numflux[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                }
            }

            m_riemann->Solve(m_spaceDim, Fwd, Bwd, TraceFlux);
        }

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF

        void AdvectionWeakDG::v_AddVolumJacToMat( 
            const Array<OneD, MultiRegions::ExpListSharedPtr>                       &pFields,
            const int                                                               &nConvectiveFields,
            const Array<OneD, const Array<OneD,  Array<OneD, 
                Array<OneD,  Array<OneD,  NekDouble> > > > >                        &ElmtJacArray,
            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >                          &gmtxarray)
        {
            MultiRegions::ExpListSharedPtr explist = pFields[0];
                std::shared_ptr<LocalRegions::ExpansionVector> pexp = explist->GetExp();
            int ntotElmt            = (*pexp).size();
            int nElmtPnt,nElmtCoef;

            NekDouble tmp;
            DNekMatSharedPtr        tmpGmtx,ElmtMat;

            Array<OneD,NekDouble> GMat_data;
            Array<OneD,NekDouble> Elmt_data;

            Array<OneD, DNekMatSharedPtr>  mtxPerVar(ntotElmt);
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
            }

            Array<OneD, DNekMatSharedPtr>  mtxPerVarCoeff(ntotElmt);
            for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
            {
                nElmtCoef   =   elmtcoef[nelmt];
                mtxPerVarCoeff[nelmt]    =MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nElmtCoef, nElmtCoef);
            }

            for(int m = 0; m < nConvectiveFields; m++)
            {
                for(int n = 0; n < nConvectiveFields; n++)
                {
                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        (*mtxPerVarCoeff[nelmt])   =    0.0;
                        (*mtxPerVar[nelmt])        =    0.0;
                    }

                    explist->GetMatIpwrtDeriveBase(ElmtJacArray[m][n],mtxPerVar);
                    //TODO:: To reuse mtxPerVar
                    explist->AddRightIPTBaseMatrix(mtxPerVar,mtxPerVarCoeff);

                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        nElmtCoef       = elmtcoef[nelmt];
                        nElmtPnt        = elmtpnts[nelmt];

                        tmpGmtx         = gmtxarray[m][n]->GetBlock(nelmt,nelmt);
                        ElmtMat         = mtxPerVarCoeff[nelmt];

                        GMat_data       = tmpGmtx->GetPtr();
                        Elmt_data       = ElmtMat->GetPtr();

                        Vmath::Vadd(nElmtCoef*nElmtCoef,GMat_data,1,Elmt_data,1,GMat_data,1);
                    }
                }
            }
        }

        void AdvectionWeakDG::v_AddVolumJacToMat( 
            const Array<OneD, MultiRegions::ExpListSharedPtr>                       &pFields,
            const int                                                               &nConvectiveFields,
            const Array<OneD, const Array<OneD,  Array<OneD, 
                Array<OneD,  Array<OneD,  NekDouble> > > > >                        &ElmtJacArray,
            Array<OneD, Array<OneD, SNekBlkMatSharedPtr> >                          &gmtxarray)
        {
            MultiRegions::ExpListSharedPtr explist = pFields[0];
                std::shared_ptr<LocalRegions::ExpansionVector> pexp = explist->GetExp();
            int ntotElmt            = (*pexp).size();
            int nElmtPnt,nElmtCoef;

            NekDouble tmp;
            SNekMatSharedPtr        tmpGmtx;
            DNekMatSharedPtr        ElmtMat;

            Array<OneD,NekSingle> GMat_data;
            Array<OneD,NekDouble> Elmt_data;
            Array<OneD,NekSingle> Elmt_dataSingle;

            Array<OneD, DNekMatSharedPtr>  mtxPerVar(ntotElmt);
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
            }

            Array<OneD, DNekMatSharedPtr>  mtxPerVarCoeff(ntotElmt);
            for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
            {
                nElmtCoef   =   elmtcoef[nelmt];
                mtxPerVarCoeff[nelmt]    =MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nElmtCoef, nElmtCoef);
            }

            for(int m = 0; m < nConvectiveFields; m++)
            {
                for(int n = 0; n < nConvectiveFields; n++)
                {
                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        (*mtxPerVarCoeff[nelmt])   =    0.0;
                        (*mtxPerVar[nelmt])        =    0.0;
                    }

                    explist->GetMatIpwrtDeriveBase(ElmtJacArray[m][n],mtxPerVar);
                    //TODO:: To reuse mtxPerVar
                    explist->AddRightIPTBaseMatrix(mtxPerVar,mtxPerVarCoeff);

                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        nElmtCoef       = elmtcoef[nelmt];
                        nElmtPnt        = elmtpnts[nelmt];
                        int ntotDofs    = nElmtCoef*nElmtCoef;

                        if(Elmt_dataSingle.num_elements()<ntotDofs)
                        {
                            Elmt_dataSingle = Array<OneD, NekSingle> (ntotDofs);
                        }

                        tmpGmtx         = gmtxarray[m][n]->GetBlock(nelmt,nelmt);
                        ElmtMat         = mtxPerVarCoeff[nelmt];

                        GMat_data       = tmpGmtx->GetPtr();
                        Elmt_data       = ElmtMat->GetPtr();

                        for(int i=0;i<ntotDofs;i++)
                        {
                            Elmt_dataSingle[i]  =   NekSingle( Elmt_data[i] );
                        }

                        Vmath::Vadd(ntotDofs,GMat_data,1,Elmt_dataSingle,1,GMat_data,1);
                    }
                }
            }
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
#endif
    }//end of namespace SolverUtils
}//end of namespace Nektar
