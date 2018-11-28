///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionIP.cpp
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
// Description: IP diffusion class.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Diffusion/DiffusionIP.h>
#include <iostream>
#include <iomanip>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string DiffusionIP::type = GetDiffusionFactory().
            RegisterCreatorFunction("InteriorPenalty", DiffusionIP::create);

        DiffusionIP::DiffusionIP()
        {
        }

        void DiffusionIP::v_InitObject(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            m_session = pSession;

            m_session->LoadSolverInfo("ShockCaptureType",
                                  m_shockCaptureType,    "Off");	

            m_session->LoadParameter("IPDebugParameter",
                                  m_IPDebugParameter,   0.5);	

            m_session->LoadParameter("IPPenaltyFactor2",
                                  m_IPPenaltyFactor2,   1.0/12.0);	

            // Setting up the normals
            int i;
            int nDim = pFields[0]->GetCoordim(0);
            int nVariable = pFields.num_elements();
            int nTracePts = pFields[0]->GetTrace()->GetTotPoints();
            
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for(i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts,0.0);
            }
            pFields[0]->GetTrace()->GetNormals(m_traceNormals);
            m_traceNormDirctnElmtLength =   Array<OneD, NekDouble> (nTracePts,0.0);
            pFields[0]->GetTrace()->GetElmtNormalLength(m_traceNormDirctnElmtLength);

            // TODO:: to check parallel case
            Array<OneD, NekDouble> lengthstmp(nTracePts,0.0);
            pFields[0]->PeriodicBwdCopy(m_traceNormDirctnElmtLength,lengthstmp);
            Vmath::Vadd(nTracePts,lengthstmp,1,m_traceNormDirctnElmtLength,1,m_traceNormDirctnElmtLength,1);

            m_tracBwdWeight  =   Array<OneD, NekDouble> (nTracePts,0.0);
            pFields[0]->GetBwdWeight(m_tracBwdWeight);
            Array<OneD, NekDouble> tmpBwdWeight(nTracePts,0.0);
            for(int i =1; i<nVariable;i++)
            {
                pFields[i]->GetBwdWeight(tmpBwdWeight);
                Vmath::Vsub(nTracePts,tmpBwdWeight,1,m_tracBwdWeight,1,tmpBwdWeight,1);
                Vmath::Vabs(nTracePts,tmpBwdWeight,1,tmpBwdWeight,1);
                NekDouble norm = 0.0;
                for(int j = 0; j<nTracePts; j++)
                {
                    norm += tmpBwdWeight[j];
                }
                ASSERTL0(norm<1.0E-11,"different BWD for different variable not coded yet");
            }
        }
        
        void DiffusionIP::v_Diffuse(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {

            int nCoeffs   = fields[0]->GetNcoeffs();
            Array<OneD, Array<OneD, NekDouble> > tmp(nConvectiveFields);
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                tmp[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
            }
            DiffusionIP::v_Diffuse_coeff(nConvectiveFields,fields,inarray,tmp,pFwd,pBwd);
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                fields[i]->BwdTrans             (tmp[i], outarray[i]);
            }
        }

        void DiffusionIP::v_Diffuse_coeff(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            const NekDouble PenaltyFactor2  =  m_IPPenaltyFactor2;
            int i, j;
            int nDim      = fields[0]->GetCoordim(0);
            int nPts      = fields[0]->GetTotPoints();
            int nCoeffs   = fields[0]->GetNcoeffs();
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();

            Array<OneD, Array<OneD, NekDouble> > tmparray2D;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > tmparray3D;
            tmparray2D = NullNekDoubleArrayofArray;
            tmparray3D = NullNekDoubleArrayofArrayofArray;
            Array<OneD, NekDouble>    Fwd   =   NullNekDouble1DArray;
            Array<OneD, NekDouble>    Bwd   =   NullNekDouble1DArray;
            Fwd = Array<OneD, NekDouble>(nTracePts,0.0);
            Bwd = Array<OneD, NekDouble>(nTracePts,0.0);

            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > elmtFlux(nDim);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > qfield(nDim);
            for (j = 0; j < nDim; ++j)
            {
                qfield[j]       = Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                elmtFlux[j]     = Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    qfield[j][i] = Array<OneD, NekDouble>(nPts, 0.0);
                }
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    elmtFlux[j][i]   = Array<OneD, NekDouble>(nPts, 0.0);
                }
            }

            Array<OneD, Array<OneD, NekDouble> >    vFwd(nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> >    vBwd(nConvectiveFields);
            

            if (pFwd == NullNekDoubleArrayofArray ||
                pBwd == NullNekDoubleArrayofArray)
            {
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    vFwd[i]    =   Array<OneD, NekDouble>(nTracePts,0.0);
                    vBwd[i]    =   Array<OneD, NekDouble>(nTracePts,0.0);
                }
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    fields[i]->GetFwdBwdTracePhys(inarray[i], vFwd[i], vBwd[i]);
                }
            }
            else
            {
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    vFwd[i]    =   pFwd[i];
                    vBwd[i]    =   pBwd[i];
                }
            }

            physFieldDeriv(nConvectiveFields,fields,inarray,vFwd,vBwd,qfield);

            // m_FunctorDerivBndCond(inarray,qfield,m_time,vFwd,tmparray3D);

            Array<OneD, NekDouble> muvar        =   NullNekDouble1DArray;
            Array<OneD, NekDouble> MuVarTrace   =   NullNekDouble1DArray;
            if (m_ArtificialDiffusionVector)
            {

                MuVarTrace  =   Array<OneD, NekDouble>(nTracePts, 0.0);
                muvar       =   Array<OneD, NekDouble>(nPts, 0.0);
                GetAVmu(fields,inarray,muvar,MuVarTrace);
            }

            Array<OneD, int > nonZeroIndex;
            // TODO: qfield AND elmtFlux share storage????
            m_FunctorDiffusionfluxCons(nConvectiveFields,nDim,inarray,qfield, elmtFlux,nonZeroIndex,tmparray2D,muvar);
            
            Array<OneD, Array<OneD, NekDouble> > tmpFluxIprdct(nDim);
            // volume intergration: the nonZeroIndex indicates which flux is nonzero
            for(i = 0; i < nonZeroIndex.num_elements(); ++i)
            {
                int j = nonZeroIndex[i];
                for (int k = 0; k < nDim; ++k)
                {
                    tmpFluxIprdct[k] = elmtFlux[k][j];
                }
                fields[j]->IProductWRTDerivBase(tmpFluxIprdct,outarray[j]);
                Vmath::Neg                      (nCoeffs, outarray[j], 1);
            }
            // release qfield, elmtFlux and muvar;
            muvar   =   NullNekDouble1DArray;
            for (j = 0; j < nDim; ++j)
            {
                elmtFlux[j]     = NullNekDoubleArrayofArray;
            }

            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > traceflux(1);
            traceflux[0]    =   Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
            for (int j = 0; j < nConvectiveFields; ++j)
            {
                traceflux[0][j]   = Array<OneD, NekDouble>(nTracePts, 0.0);
            }

            Array<OneD, Array<OneD, NekDouble> >    solution_Aver(nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> >    solution_jump(nConvectiveFields);
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                solution_Aver[i]    =   Array<OneD, NekDouble>(nTracePts,0.0);
                solution_jump[i]    =   Array<OneD, NekDouble>(nTracePts,0.0);
            }

            CalTraceNumFlux(nConvectiveFields,nDim,nPts,nTracePts,PenaltyFactor2,
                            fields,inarray,qfield,vFwd,vBwd,MuVarTrace,nonZeroIndex,traceflux,solution_Aver,solution_jump);

            if(abs(m_IPDebugParameter)>1.0E-12)
            {
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > traceSymflux(nDim);
                for (int nd = 0; nd < nDim; ++nd)
                {
                    traceSymflux[nd]    = Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
                    for (int j = 0; j < nConvectiveFields; ++j)
                    {
                        traceSymflux[nd][j]   = Array<OneD, NekDouble>(nTracePts, 0.0);
                    }
                }

                Array<OneD, int > nonZeroIndexSymm;
                CalTraceSymFlux(nConvectiveFields,nDim,fields,solution_Aver,solution_jump,
                            nonZeroIndexSymm,traceSymflux);

                for (int i = 0; i < nConvectiveFields; ++i)
                {
                    solution_Aver[i]    =   NullNekDouble1DArray;
                    solution_jump[i]    =   NullNekDouble1DArray;
                }


                AddSymmFluxIntegral(nConvectiveFields,nDim,nPts,nTracePts,fields,nonZeroIndexSymm,traceSymflux,outarray);
            }

            for (int i = 0; i < nConvectiveFields; ++i)
            {
                solution_Aver[i]    =   NullNekDouble1DArray;
                solution_jump[i]    =   NullNekDouble1DArray;
            }

            for(i = 0; i < nonZeroIndex.num_elements(); ++i)
            {
                int j = nonZeroIndex[i];

                fields[j]->AddTraceIntegral     (traceflux[0][j], outarray[j]);
                fields[j]->SetPhysState         (false);
                fields[j]->MultiplyByElmtInvMass(outarray[j], outarray[j]);
            }
        }

        void DiffusionIP::v_CalTraceNumFlux(
            const int                                                           nConvectiveFields,
            const int                                                           nDim,
            const int                                                           nPts,
            const int                                                           nTracePts,
            const NekDouble                                                     PenaltyFactor2,
            const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
            const Array<OneD, Array<OneD, NekDouble> >                          &inarray,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &qfield,
            const Array<OneD, Array<OneD, NekDouble> >                          &vFwd,
            const Array<OneD, Array<OneD, NekDouble> >                          &vBwd,
            const Array<OneD, NekDouble >                                       &MuVarTrace,
                  Array<OneD, int >                                             &nonZeroIndexflux,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &traceflux,
                  Array<OneD, Array<OneD, NekDouble> >                          &solution_Aver,
                  Array<OneD, Array<OneD, NekDouble> >                          &solution_jump)
        {
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    numDeriv(nDim);
            for (int nd = 0; nd < nDim; ++nd)
            {
                numDeriv[nd]     =   Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
                for (int i = 0; i < nConvectiveFields; ++i)
                {
                    numDeriv[nd][i]    = Array<OneD, NekDouble>(nTracePts,0.0);
                }
            }

            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  tmparray3D = NullNekDoubleArrayofArrayofArray;
            m_FunctorDerivBndCond(inarray,qfield,m_time,vFwd,tmparray3D);

            Array<OneD, NekDouble> Fwd(nTracePts,0.0);
            Array<OneD, NekDouble> Bwd(nTracePts,0.0);

            if(abs(PenaltyFactor2)>1.0E-12)
            {
                Add2ndDeriv2Trace(nConvectiveFields,nDim,nPts,nTracePts,PenaltyFactor2,fields,qfield,numDeriv);
            }

            for (int nd = 0; nd < nDim; ++nd)
            {
                for (int i = 0; i < nConvectiveFields; ++i)
                {
                    Vmath::Fill(nTracePts, 0.0, Bwd,1);
                    fields[i]->GetFwdBwdTracePhysDeriv(nd,qfield[nd][i], Fwd, Bwd);
                    for (int nt = 0; nt < nTracePts; ++nt)
                    {
                        numDeriv[nd][i][nt] +=   0.5*( Fwd[nt] + Bwd[nt] );  
                    }
                }
            }
        
            ConsVarAveJump(nConvectiveFields,nTracePts,vFwd,vBwd,solution_Aver,solution_jump);

            Array<OneD, NekDouble>  jumpTmp         =   Fwd;
            Array<OneD, NekDouble>  PenaltyFactor   =   Bwd;
            Fwd     =   NullNekDouble1DArray;
            Bwd     =   NullNekDouble1DArray;
            GetPenaltyFactor(fields,PenaltyFactor);
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Vmul(nTracePts,solution_jump[i],1, PenaltyFactor,1,jumpTmp,1);
                Vmath::Vdiv(nTracePts,jumpTmp,1, m_traceNormDirctnElmtLength,1,jumpTmp,1);
                for (int nd = 0; nd < nDim; ++nd)
                {
                    Vmath::Vvtvp(nTracePts, m_traceNormals[nd],1,jumpTmp,1, numDeriv[nd][i],1, numDeriv[nd][i],1);
                }
            }
            jumpTmp         =   NullNekDouble1DArray;
            PenaltyFactor   =   NullNekDouble1DArray;

            // Calculate normal viscous flux
            m_FunctorDiffusionfluxCons(nConvectiveFields,nDim,solution_Aver,numDeriv,traceflux,nonZeroIndexflux,m_traceNormals,MuVarTrace);
        }


        void DiffusionIP::CalTraceSymFlux(
            const int                                                           nConvectiveFields,
            const int                                                           nDim,
            const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
            const Array<OneD, Array<OneD, NekDouble> >                          &solution_Aver,
                  Array<OneD, Array<OneD, NekDouble> >                          &solution_jump,
                  Array<OneD, int >                                             &nonZeroIndexsymm,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &traceSymflux)
        {
            int nTracePts = solution_jump[nConvectiveFields-1].num_elements();
// //Debug
// for (int j = 0; j < nConvectiveFields; ++j)
// {
//     for (int k = 0; k < nTracePts; ++k)
//     {
//         cout<< "solution_jump["<<j<<"]["<<k<<"]= "<<solution_jump[j][k]<<endl;
//     }
// }
            
// //Debug
// for (int j = 0; j < nConvectiveFields; ++j)
// {
//     for (int k = 0; k < nTracePts; ++k)
//     {
//         cout<< "solution_jump["<<j<<"]["<<k<<"]= "<<solution_jump[j][k]<<endl;
//     }
// }
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Smul(nTracePts,m_IPDebugParameter,solution_jump[i],1,solution_jump[i],1);
            }
            
            m_FunctorSymmetricfluxCons(nConvectiveFields,nDim,solution_Aver,solution_jump,traceSymflux,nonZeroIndexsymm,m_traceNormals);

            
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                MultiRegions::ExpListSharedPtr tracelist = fields[i]->GetTrace();
                for(int nd=0;nd<nDim;nd++)
                {
                    tracelist->MultiplyByQuadratureMetric(traceSymflux[nd][i],traceSymflux[nd][i]);
                }
            }
        }

        void DiffusionIP::Add2ndDeriv2Trace(
            const int                                                           nConvectiveFields,
            const int                                                           nDim,
            const int                                                           nPts,
            const int                                                           nTracePts,
            const NekDouble                                                     PenaltyFactor2,
            const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &qfield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &numDeriv)
        {
            if(PenaltyFactor2>0.0)
            {
                const NekDouble PenaltyFactor2 = 1.0/12.0;
                Array<OneD, NekDouble> Fwd(nTracePts,0.0);
                Array<OneD, NekDouble> Bwd(nTracePts,0.0);

                Array<OneD, Array<OneD, NekDouble> > elmt2ndDerv(nDim);
                for(int nd1=0; nd1<nDim; nd1++)
                {
                    elmt2ndDerv[nd1]    =   Array<OneD, NekDouble>(nPts,0.0);
                }

                Array<OneD, Array<OneD, NekDouble> > qtmp(3);
                for(int nd=0; nd<3; nd++)
                {
                    qtmp[nd]    =   NullNekDouble1DArray;
                }
                for(int nd2=0; nd2<nDim; nd2++)
                {
                    qtmp[nd2]    =   elmt2ndDerv[nd2];
                }

                // the derivatives are assumed to be exchangable  
                for(int nd1=0; nd1<nDim; nd1++)
                {
                    for(int i = 0; i < nConvectiveFields; ++i)
                    {
                        fields[i]->PhysDeriv(qfield[nd1][i], qtmp[0], qtmp[1], qtmp[2]);

                        // for(int nd2=0; nd2<nDim; nd2++)
                        // {
                        //     Vmath::Fill(nTracePts, 0.0, Bwd,1);
                        //     fields[i]->GetFwdBwdTracePhysNoBndFill(elmt2ndDerv[nd2], Fwd, Bwd);
                        //     Vmath::Vsub(nTracePts,Bwd,1,Fwd,1,Fwd,1);
                        //     Vmath::Smul(nTracePts,PenaltyFactor2,Fwd,1,Fwd,1);
                        //     Vmath::Vmul(nTracePts,m_traceNormDirctnElmtLength,1,Fwd,1,Fwd,1);
                        //     // Vmath::Vvtvp(nTracePts,m_traceNormals[nd2],1,Fwd,1,numDeriv[nd1][i],1,numDeriv[nd1][i],1);
                        //     Vmath::Vvtvp(nTracePts,m_traceNormals[nd1],1,Fwd,1,numDeriv[nd2][i],1,numDeriv[nd2][i],1);
                        // }

                        for(int nd2=nd1; nd2<nDim; nd2++)
                        {
                            Vmath::Fill(nTracePts, 0.0, Bwd,1);
                            fields[i]->GetFwdBwdTracePhysDeriv(nd2,elmt2ndDerv[nd2], Fwd, Bwd);
                            Vmath::Vsub(nTracePts,Bwd,1,Fwd,1,Fwd,1);
                            Vmath::Smul(nTracePts,PenaltyFactor2,Fwd,1,Fwd,1);
                            Vmath::Vmul(nTracePts,m_traceNormDirctnElmtLength,1,Fwd,1,Fwd,1);
                            Vmath::Vvtvp(nTracePts,m_traceNormals[nd2],1,Fwd,1,numDeriv[nd1][i],1,numDeriv[nd1][i],1);
                            if(nd2!=nd1)
                            {
                                Vmath::Vvtvp(nTracePts,m_traceNormals[nd1],1,Fwd,1,numDeriv[nd2][i],1,numDeriv[nd2][i],1);
                            }
                        }
                    }
                }
            }
        }

        void DiffusionIP::AddSymmFluxIntegral(
            const int                                                           nConvectiveFields,
            const int                                                           nDim,
            const int                                                           nPts,
            const int                                                           nTracePts,
            const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
            const Array<OneD, const int >                                       &nonZeroIndex,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &tracflux,
                  Array<OneD, Array<OneD, NekDouble> >                          &outarray)
        {
            int nCoeffs =   outarray[nConvectiveFields-1].num_elements();
            Array<OneD, NekDouble > tmpCoeff(nCoeffs,0.0);
            Array<OneD, Array<OneD, NekDouble> > tmpfield(nDim);
            for(int i = 0;i<nDim;i++)
            {
                tmpfield[i]    =   Array<OneD, NekDouble>(nPts,0.0);
            }
            int nv = 0;
            for(int j=0;j<nonZeroIndex.num_elements();j++)
            {
                nv  =   nonZeroIndex[j];
                for(int nd=0;nd<nDim;nd++)
                {
                    Vmath::Zero(nPts,tmpfield[nd],1);

                    fields[nv]->AddTraceQuadPhysToField(tracflux[nd][nv],tracflux[nd][nv],tmpfield[nd]);
                    fields[nv]->DividByQuadratureMetric(tmpfield[nd],tmpfield[nd]);
                }
                fields[nv]->IProductWRTDerivBase(tmpfield,tmpCoeff);
                Vmath::Vadd(nCoeffs,tmpCoeff,1,outarray[nv],1,outarray[nv],1);
            }
        }
        
        void DiffusionIP::v_physFieldDeriv(
            const int                                                   nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>           &fields,
            const Array<OneD, Array<OneD, NekDouble> >                  &inarray,
            const Array<OneD, Array<OneD, NekDouble> >                  &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >                  &pBwd,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    &qfield)
        {
            int nDim      = fields[0]->GetCoordim(0);

            Array<OneD, Array<OneD, NekDouble> > qtmp(3);
            for(int nd=0; nd<3; nd++)
            {
                qtmp[nd]    =   NullNekDouble1DArray;
            }
            for(int i = 0; i < nConvectiveFields; ++i)
            {
                for(int nd=0; nd<nDim; nd++)
                {
                    qtmp[nd]    =   qfield[nd][i];
                }
                fields[i]->PhysDeriv(inarray[i], qtmp[0], qtmp[1], qtmp[2]);
            }
        }
        

        void DiffusionIP::GetPenaltyFactor(
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                  Array<OneD, NekDouble >                       factor)
        {
            MultiRegions::ExpListSharedPtr tracelist = fields[0]->GetTrace();
            std::shared_ptr<LocalRegions::ExpansionVector> traceExp= tracelist->GetExp();
            int ntotTrac            = (*traceExp).size();
            int nTracPnt,noffset;
            
            const MultiRegions::LocTraceToTraceMapSharedPtr locTraceToTraceMap = fields[0]->GetlocTraceToTraceMap();
            
            const Array<OneD, const Array<OneD, int >> LRAdjExpid  =   locTraceToTraceMap->GetLeftRightAdjacentExpId();
            const Array<OneD, const Array<OneD, bool>> LRAdjflag   =   locTraceToTraceMap->GetLeftRightAdjacentExpFlag();



            std::shared_ptr<LocalRegions::ExpansionVector> fieldExp= fields[0]->GetExp();

            Array<OneD, NekDouble > factorFwdBwd(2,0.0);

            NekDouble spaceDim    =   NekDouble( fields[0]->GetCoordim(0) );

            for(int ntrace = 0; ntrace < ntotTrac; ++ntrace)
            {
                noffset     = tracelist->GetPhys_Offset(ntrace);
                nTracPnt    = tracelist->GetTotPoints(ntrace);

                factorFwdBwd[0] =   0.0;
                factorFwdBwd[1] =   0.0;
                
                for(int  nlr = 0; nlr < 2; nlr++)
                {
                    if(LRAdjflag[nlr][ntrace])
                    {
                        int numModes        = fields[0]->GetNcoeffs(LRAdjExpid[nlr][ntrace]);  
                        NekDouble numModesdir     = pow(NekDouble(numModes),(1.0/spaceDim));
                        factorFwdBwd[nlr]   =   1.0 * numModesdir * (numModesdir + 1.0);
                    }
                }

                for(int np = 0; np < nTracPnt; ++np)
                {
                    factor[noffset+np]    =   max(factorFwdBwd[0],factorFwdBwd[1]);
//debug
factor[noffset+np]    =   4.0;
                }
            }
        }

        void DiffusionIP::v_ConsVarAveJump(
            const int                                           nConvectiveFields,
            const int                                           npnts,
            const Array<OneD, const Array<OneD, NekDouble> >    &vFwd,
            const Array<OneD, const Array<OneD, NekDouble> >    &vBwd,
                  Array<OneD,       Array<OneD, NekDouble> >    &aver,
                  Array<OneD,       Array<OneD, NekDouble> >    &jump)
        {
            // TODO:    Direchlet boundary are not that accurate using this solution_Aver
            ConsVarAve(nConvectiveFields,npnts,vFwd,vBwd,aver);

            m_SpecialBndTreat(nConvectiveFields,aver);

            // note: here the jump is 2.0*(aver-vFwd) 
            //       because Viscous wall use a symmetry boundary not the target one   
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Vsub(npnts,aver[i],1,vFwd[i],1,jump[i],1);
                Vmath::Smul(npnts,2.0,jump[i],1,jump[i],1);
            }
        }
        
        void DiffusionIP::ConsVarAve(
            const int                                           nConvectiveFields,
            const int                                           npnts,
            const Array<OneD, const Array<OneD, NekDouble> >    &vFwd,
            const Array<OneD, const Array<OneD, NekDouble> >    &vBwd,
                  Array<OneD,       Array<OneD, NekDouble> >    &aver)
        {
            NekDouble LinternalEngy =0.0;
            NekDouble RinternalEngy =0.0;
            NekDouble AinternalEngy =0.0;

            Array<OneD, NekDouble> Fweight (npnts,1.0);
            Array<OneD, NekDouble> Bweight;

            Bweight = m_tracBwdWeight;

            Vmath::Vsub(npnts,Fweight,1,Bweight,1,Fweight,1);

            for (int i = 0; i < nConvectiveFields-1; ++i)
            {
                Vmath::Vmul (npnts,Fweight,1,vFwd[i],1,aver[i],1);
                Vmath::Vvtvp(npnts,Bweight,1,vBwd[i],1,aver[i],1,aver[i],1);
            }
            
            int nengy = nConvectiveFields-1;
            int nvelst    = 1;
            int nveled    = nengy;
            for (int nt = 0; nt < npnts; ++nt)
            {
                LinternalEngy =0.0;
                for(int j=nvelst;j<nveled;j++)
                {
                    LinternalEngy += vFwd[j][nt]*vFwd[j][nt];
                }
                LinternalEngy *= -0.5/vFwd[0][nt];
                LinternalEngy += vFwd[nengy][nt];

                RinternalEngy =0.0;
                for(int j=nvelst;j<nveled;j++)
                {
                    RinternalEngy += vBwd[j][nt]*vBwd[j][nt];
                }
                RinternalEngy *= -0.5/vBwd[0][nt];
                RinternalEngy += vBwd[nengy][nt];

                AinternalEngy =0.0;
                aver[nengy][nt] = Fweight[nt]*LinternalEngy + Bweight[nt]*RinternalEngy;
                for(int j=nvelst;j<nveled;j++)
                {
                    AinternalEngy += aver[j][nt]*aver[j][nt];
                }
                aver[nengy][nt] += AinternalEngy*(0.5/aver[0][nt]);
            }
        }

        void DiffusionIP::v_AddVolumDerivJac2Mat( 
            const int                                               nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>       &pFields,
            const Array<OneD, const Array<OneD, DNekMatSharedPtr> > &ElmtJac,
            const int                                               nfluxDir, 
            const int                                               nDervDir, 
                  Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >    &gmtxarray)
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
                    explist->GetMatIpwrtDeriveBase(JacArray,nfluxDir,mtxPerVar);
                    explist->AddRightIPTPhysDerivBase(nDervDir,mtxPerVar,mtxPerVar);

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

        void DiffusionIP::v_Diffuse_coeffOld(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            const NekDouble PenaltyFactor2 = 1.0/12.0;

            int i, j;
            int nDim      = fields[0]->GetCoordim(0);
            int nPts      = fields[0]->GetTotPoints();
            int nCoeffs   = fields[0]->GetNcoeffs();
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();

            Array<OneD, Array<OneD, NekDouble> > tmparray2D;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > tmparray3D;
            tmparray2D = NullNekDoubleArrayofArray;
            tmparray3D = NullNekDoubleArrayofArrayofArray;
            Array<OneD, NekDouble>    Fwd   =   NullNekDouble1DArray;
            Array<OneD, NekDouble>    Bwd   =   NullNekDouble1DArray;
            Fwd = Array<OneD, NekDouble>(nTracePts,0.0);
            Bwd = Array<OneD, NekDouble>(nTracePts,0.0);

            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > elmtFlux(nDim);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > qfield(nDim);
            for (j = 0; j < nDim; ++j)
            {
                qfield[j]       = Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                elmtFlux[j]     = Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    qfield[j][i] = Array<OneD, NekDouble>(nPts, 0.0);
                }
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    elmtFlux[j][i]   = Array<OneD, NekDouble>(nPts, 0.0);
                }
            }

            Array<OneD, Array<OneD, NekDouble> > qtmp(3);
            for(int nd=0; nd<3; nd++)
            {
                qtmp[nd]    =   NullNekDouble1DArray;
            }
            for(int i = 0; i < nConvectiveFields; ++i)
            {
                for(int nd=0; nd<nDim; nd++)
                {
                    qtmp[nd]    =   qfield[nd][i];
                }
                fields[i]->PhysDeriv(inarray[i], qtmp[0], qtmp[1], qtmp[2]);
            }

            m_FunctorDerivBndCond(inarray,qfield,m_time,pFwd,tmparray3D);

            Array<OneD, NekDouble> muvar        =   NullNekDouble1DArray;
            Array<OneD, NekDouble> MuVarTrace   =   NullNekDouble1DArray;
            if (m_ArtificialDiffusionVector)
            {
                MuVarTrace  =   Array<OneD, NekDouble>(nTracePts, 0.0);
                muvar       =   Array<OneD, NekDouble>(nPts, 0.0);
                m_ArtificialDiffusionVector(inarray, muvar);

                // for(int i=0;i<nPts;i++)
                // {
                //     if(muvar[i]>0.0)
                //     {
                //         cout << "   muvar["<<i<<"]= "<<muvar[i]<<endl;
                //     }
                // }

                // BwdMuvar is left to be 0.0 according to DiffusionLDG.cpp
                fields[0]->GetFwdBwdTracePhysNoBndFill(muvar,Fwd,Bwd);

                for(int k = 0; k < nTracePts; ++k)
                {
                    MuVarTrace[k] = 0.5 * (Fwd[k] + Bwd[k]) ;
                    // if(MuVarTrace[k]>0.0)
                    // {
                    //     cout << "   MuVarTrace["<<k<<"]= "<<MuVarTrace[k]<<endl;
                    // }
                }
            }
            
            // int ndebug=0;
            // cin >>ndebug;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    numDeriv(nDim);
            for (int nd = 0; nd < nDim; ++nd)
            {
                numDeriv[nd]     =   Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
                for (int i = 0; i < nConvectiveFields; ++i)
                {
                    numDeriv[nd][i]    = Array<OneD, NekDouble>(nTracePts,0.0);
                }
            }
            
            Array<OneD, Array<OneD, NekDouble> > elmt2ndDerv(nDim);
            for(int nd1=0; nd1<nDim; nd1++)
            {
                elmt2ndDerv[nd1]    =   Array<OneD, NekDouble>(nPts,0.0);
            }

            for(int nd2=0; nd2<nDim; nd2++)
            {
                qtmp[nd2]    =   elmt2ndDerv[nd2];
            }

            // the derivatives are assumed to be exchangable  
            for(int nd1=0; nd1<nDim; nd1++)
            {
                for(int i = 0; i < nConvectiveFields; ++i)
                {
                    fields[i]->PhysDeriv(qfield[nd1][i], qtmp[0], qtmp[1], qtmp[2]);

                    // for(int nd2=0; nd2<nDim; nd2++)
                    // {
                    //     Vmath::Fill(nTracePts, 0.0, Bwd,1);
                    //     fields[i]->GetFwdBwdTracePhysNoBndFill(elmt2ndDerv[nd2], Fwd, Bwd);
                    //     Vmath::Vsub(nTracePts,Bwd,1,Fwd,1,Fwd,1);
                    //     Vmath::Smul(nTracePts,PenaltyFactor2,Fwd,1,Fwd,1);
                    //     Vmath::Vmul(nTracePts,m_traceNormDirctnElmtLength,1,Fwd,1,Fwd,1);
                    //     // Vmath::Vvtvp(nTracePts,m_traceNormals[nd2],1,Fwd,1,numDeriv[nd1][i],1,numDeriv[nd1][i],1);
                    //     Vmath::Vvtvp(nTracePts,m_traceNormals[nd1],1,Fwd,1,numDeriv[nd2][i],1,numDeriv[nd2][i],1);
                    // }

                    for(int nd2=nd1; nd2<nDim; nd2++)
                    {
                        Vmath::Fill(nTracePts, 0.0, Bwd,1);
                        fields[i]->GetFwdBwdTracePhysNoBndFill(elmt2ndDerv[nd2], Fwd, Bwd);
                        Vmath::Vsub(nTracePts,Bwd,1,Fwd,1,Fwd,1);
                        Vmath::Smul(nTracePts,PenaltyFactor2,Fwd,1,Fwd,1);
                        Vmath::Vmul(nTracePts,m_traceNormDirctnElmtLength,1,Fwd,1,Fwd,1);
                        Vmath::Vvtvp(nTracePts,m_traceNormals[nd2],1,Fwd,1,numDeriv[nd1][i],1,numDeriv[nd1][i],1);
                        if(nd2!=nd1)
                        {
                            Vmath::Vvtvp(nTracePts,m_traceNormals[nd1],1,Fwd,1,numDeriv[nd2][i],1,numDeriv[nd2][i],1);
                        }
                    }
                }
            }

            for (int nd = 0; nd < nDim; ++nd)
            {
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    Vmath::Fill(nTracePts, 0.0, Bwd,1);
                    fields[i]->GetFwdBwdTracePhysDeriv(nd,qfield[nd][i], Fwd, Bwd);
                    for (int nt = 0; nt < nTracePts; ++nt)
                    {
                        numDeriv[nd][i][nt] +=   0.5*( Fwd[nt] + Bwd[nt] );  
                    }
                }
            }

            Array<OneD, int > nonZeroIndex;
            // TODO: qfield AND elmtFlux share storage????
            m_FunctorDiffusionfluxCons(nConvectiveFields,nDim,inarray,qfield, elmtFlux,nonZeroIndex,tmparray2D,muvar);
            
            Array<OneD, Array<OneD, NekDouble> > tmpFluxIprdct(nDim);
            // volume intergration: the nonZeroIndex indicates which flux is nonzero
            for(i = 0; i < nonZeroIndex.num_elements(); ++i)
            {
                int j = nonZeroIndex[i];
                for (int k = 0; k < nDim; ++k)
                {
                    tmpFluxIprdct[k] = elmtFlux[k][j];
                }
                fields[j]->IProductWRTDerivBase(tmpFluxIprdct,outarray[j]);
                Vmath::Neg                      (nCoeffs, outarray[j], 1);
            }
            // release qfield, elmtFlux and muvar;
            muvar   =   NullNekDouble1DArray;
            for (j = 0; j < nDim; ++j)
            {
                qfield[j]       = NullNekDoubleArrayofArray;
                elmtFlux[j]     = NullNekDoubleArrayofArray;
            }
            
            Array<OneD, Array<OneD, NekDouble> >    solution_jump(nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> >    solution_Aver(nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> >    vFwd(nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> >    vBwd(nConvectiveFields);
            for (i = 0; i < nConvectiveFields; ++i)
            {
                solution_jump[i]    =   Array<OneD, NekDouble>(nTracePts,0.0);
                solution_Aver[i]    =   Array<OneD, NekDouble>(nTracePts,0.0);
            }

            if (pFwd == NullNekDoubleArrayofArray ||
                pBwd == NullNekDoubleArrayofArray)
            {
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    vFwd[i]    =   Array<OneD, NekDouble>(nTracePts,0.0);
                    vBwd[i]    =   Array<OneD, NekDouble>(nTracePts,0.0);
                }
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    fields[i]->GetFwdBwdTracePhys(inarray[i], vFwd[i], vBwd[i]);
                }
            }
            else
            {
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    vFwd[i]    =   pFwd[i];
                    vBwd[i]    =   pBwd[i];
                }
            }

            // TODO:    Direchlet boundary are not that accurate using this solution_Aver
            ConsVarAve(nConvectiveFields,nTracePts,vFwd,vBwd,solution_Aver);

            m_SpecialBndTreat(nConvectiveFields,solution_Aver);

            // note: here the jump is 2.0*(aver-vFwd) 
            //       because Viscous wall use a symmetry boundary not the target one   
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Vsub(nTracePts,solution_Aver[i],1,vFwd[i],1,solution_jump[i],1);
                Vmath::Smul(nTracePts,2.0,solution_jump[i],1,solution_jump[i],1);
            }
            
            Array<OneD, NekDouble>  jumpTmp         =   Fwd;
            Array<OneD, NekDouble>  PenaltyFactor   =   Bwd;
            Fwd     =   NullNekDouble1DArray;
            Bwd     =   NullNekDouble1DArray;
            GetPenaltyFactor(fields,PenaltyFactor);
            for (i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Vmul(nTracePts,solution_jump[i],1, PenaltyFactor,1,jumpTmp,1);
                Vmath::Vdiv(nTracePts,jumpTmp,1, m_traceNormDirctnElmtLength,1,jumpTmp,1);
                for (int nd = 0; nd < nDim; ++nd)
                {
                    Vmath::Vvtvp(nTracePts, m_traceNormals[nd],1,jumpTmp,1, numDeriv[nd][i],1, numDeriv[nd][i],1);
                }
            }
            jumpTmp         =   NullNekDouble1DArray;
            PenaltyFactor   =   NullNekDouble1DArray;

            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > traceflux(1);
            traceflux[0]    =   Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
            for (int j = 0; j < nConvectiveFields; ++j)
            {
                traceflux[0][j]   = Array<OneD, NekDouble>(nTracePts, 0.0);
            }
            // Calculate normal viscous flux
            m_FunctorDiffusionfluxCons(nConvectiveFields,nDim,solution_Aver,numDeriv,traceflux,nonZeroIndex,m_traceNormals,MuVarTrace);

            for(i = 0; i < nonZeroIndex.num_elements(); ++i)
            {
                int j = nonZeroIndex[i];
                fields[j]->AddTraceIntegral     (traceflux[0][j], outarray[j]);
                fields[j]->SetPhysState         (false);
                fields[j]->MultiplyByElmtInvMass(outarray[j], outarray[j]);
            }
        }
    }
}
