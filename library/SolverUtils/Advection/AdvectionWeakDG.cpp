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

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Advection/AdvectionWeakDG.h>
#include <iostream>
#include <iomanip>

#include <LibUtilities/BasicUtils/Timer.h>

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
            LibUtilities::Timer timer1;
            timer1.Start();
            
            boost::ignore_unused(advVel, time);
            size_t nCoeffs         = fields[0]->GetNcoeffs();

            Array<OneD, Array<OneD, NekDouble> > tmp{size_t(nConvectiveFields)};
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                tmp[i] = Array<OneD, NekDouble> {nCoeffs, 0.0};
            }

            AdvectionWeakDG::v_AdvectCoeffs(
                nConvectiveFields, fields, advVel, inarray, tmp, time,
                pFwd, pBwd);

            // why was this broken in many loops over convective fields?
            // this is terrible for locality
            LibUtilities::Timer timer;
            timer.Start();
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                fields[i]->BwdTrans(tmp[i], outarray[i]);
            }
            timer.Stop();
            timer.AccumulateRegion("AdvWeakDG:_BwdTrans",1);
            timer1.Stop();
            timer1.AccumulateRegion("AdvWeakDG:All");
        }

        void AdvectionWeakDG::v_AdvectCoeffs(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &advVel,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const NekDouble                                   &time,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            LibUtilities::Timer timer1;
            timer1.Start();
            size_t nPointsTot      = fields[0]->GetTotPoints();
            size_t nCoeffs         = fields[0]->GetNcoeffs();
            size_t nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();

            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                fluxvector(nConvectiveFields);
            // Allocate storage for flux vector F(u).
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                fluxvector[i] =
                    Array<OneD, Array<OneD, NekDouble> > {size_t(m_spaceDim)};
                for (int j = 0; j < m_spaceDim; ++j)
                {
                    fluxvector[i][j] = Array<OneD, NekDouble>(nPointsTot,0.0);
                }
            }

            LibUtilities::Timer timer;
            timer.Start();
            v_AdvectVolumeFlux(nConvectiveFields, fields, advVel, inarray,
                                fluxvector, time);
            timer.Stop();
            timer.AccumulateRegion("AdvWeakDG:_fluxVector",1);


            timer.Start();
            // Get the advection part (without numerical flux)
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Fill(outarray[i].size(), 0.0, outarray[i], 1);
                fields[i]->IProductWRTDerivBase(fluxvector[i], outarray[i]);
            }
            timer.Stop();
            timer.AccumulateRegion("AdvWeakDG:_IProductWRTDerivBase",1);

            Array<OneD, Array<OneD, NekDouble> >
                numflux{size_t(nConvectiveFields)};

            for (int i = 0; i < nConvectiveFields; ++i)
            {
                numflux[i] = Array<OneD, NekDouble> {nTracePointsTot, 0.0};
            }

            v_AdvectTraceFlux(nConvectiveFields, fields, advVel, inarray,
                numflux, time, pFwd, pBwd);

            // Evaulate <\phi, \hat{F}\cdot n> - OutField[i]
            for(int i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Neg                      (nCoeffs, outarray[i], 1);

                timer.Start();
                fields[i]->AddTraceIntegral     (numflux[i], outarray[i]);
                timer.Stop();
                timer.AccumulateRegion("AdvWeakDG:_AddTraceIntegral",1);

                timer.Start();
                fields[i]->MultiplyByElmtInvMass(outarray[i], outarray[i]);
                timer.Stop();
                timer.AccumulateRegion("AdvWeakDG:_MultiplyByElmtInvMass",1);
            }
            timer1.Stop();
            timer1.AccumulateRegion("AdvWeakDG: Coeff All");
        }
        
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
            boost::ignore_unused(advVel, time);
            int nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();

            ASSERTL1(m_riemann,
                "Riemann solver must be provided for AdvectionWeakDG.");

            // Store forwards/backwards space along trace space
            Array<OneD, Array<OneD, NekDouble>> Fwd(nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble>> Bwd(nConvectiveFields);

            if (pFwd == NullNekDoubleArrayofArray ||
                pBwd == NullNekDoubleArrayofArray)
            {
                for (int i = 0; i < nConvectiveFields; ++i)
                {
                    Fwd[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                    Bwd[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                    fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
                }
            }
            else
            {
                for (int i = 0; i < nConvectiveFields; ++i)
                {
                    Fwd[i] = pFwd[i];
                    Bwd[i] = pBwd[i];
                }
            }

            LibUtilities::Timer timer;
            timer.Start();
            m_riemann->Solve(m_spaceDim, Fwd, Bwd, TraceFlux);
            timer.Stop();
            timer.AccumulateRegion("AdvWeakDG:_Riemann",1);

        }

        void AdvectionWeakDG::v_AddVolumJacToMat( 
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const int &nConvectiveFields,
            const TensorOfArray5D<NekDouble> &ElmtJacArray,
            Array<OneD, Array<OneD, SNekBlkMatSharedPtr>> &gmtxarray)
        {
            MultiRegions::ExpListSharedPtr explist = pFields[0];
                std::shared_ptr<LocalRegions::ExpansionVector> pexp = explist->GetExp();
            int nTotElmt            = (*pexp).size();
            int nElmtPnt,nElmtCoef;

            SNekMatSharedPtr        tmpGmtx;
            DNekMatSharedPtr        ElmtMat;

            Array<OneD,NekSingle> GMat_data;
            Array<OneD,NekDouble> Elmt_data;
            Array<OneD,NekSingle> Elmt_dataSingle;

            Array<OneD, DNekMatSharedPtr>  mtxPerVar(nTotElmt);
            Array<OneD, int > elmtpnts(nTotElmt);
            Array<OneD, int > elmtcoef(nTotElmt);
            for(int nelmt = 0; nelmt < nTotElmt; nelmt++)
            {
                nElmtCoef           = (*pexp)[nelmt]->GetNcoeffs();
                nElmtPnt            = (*pexp)[nelmt]->GetTotPoints();
                elmtpnts[nelmt]     =   nElmtPnt;
                elmtcoef[nelmt]     =   nElmtCoef;
                mtxPerVar[nelmt]    =MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nElmtCoef, nElmtPnt);
            }

            Array<OneD, DNekMatSharedPtr>  mtxPerVarCoeff(nTotElmt);
            for(int nelmt = 0; nelmt < nTotElmt; nelmt++)
            {
                nElmtCoef   =   elmtcoef[nelmt];
                mtxPerVarCoeff[nelmt]    =MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nElmtCoef, nElmtCoef);
            }

            for(int m = 0; m < nConvectiveFields; m++)
            {
                for(int n = 0; n < nConvectiveFields; n++)
                {
                    for(int nelmt = 0; nelmt < nTotElmt; nelmt++)
                    {
                        (*mtxPerVarCoeff[nelmt])   =    0.0;
                        (*mtxPerVar[nelmt])        =    0.0;
                    }

                    explist->GetMatIpwrtDeriveBase(ElmtJacArray[m][n], 
                        mtxPerVar);
                    // Memory can be saved by reusing mtxPerVar
                    explist->AddRightIPTBaseMatrix(mtxPerVar,mtxPerVarCoeff);

                    for(int nelmt = 0; nelmt < nTotElmt; nelmt++)
                    {
                        nElmtCoef = elmtcoef[nelmt];
                        nElmtPnt = elmtpnts[nelmt];
                        int ntotDofs = nElmtCoef*nElmtCoef;

                        if(Elmt_dataSingle.size()<ntotDofs)
                        {
                            Elmt_dataSingle = Array<OneD, NekSingle> (ntotDofs);
                        }

                        tmpGmtx = gmtxarray[m][n]->GetBlock(nelmt,nelmt);
                        ElmtMat = mtxPerVarCoeff[nelmt];

                        GMat_data = tmpGmtx->GetPtr();
                        Elmt_data = ElmtMat->GetPtr();

                        for(int i=0;i<ntotDofs;i++)
                        {
                            Elmt_dataSingle[i]  =   NekSingle( Elmt_data[i] );
                        }

                        Vmath::Vadd(ntotDofs,GMat_data,1,Elmt_dataSingle,1,GMat_data,1);
                    }
                }
            }
        }

    }//end of namespace SolverUtils
}//end of namespace Nektar
