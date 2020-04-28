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

            m_session->LoadParameter("IPSymmFtluxCoeff",
                                  m_IPSymmFtluxCoeff,   1.0);	//1.0：SIPG; -1.0: NIPG; 0.0:IIPG

            m_session->LoadParameter("IP2ndDervCoeff",
                                  m_IP2ndDervCoeff,   0.0); // 1.0/12.0	

            m_session->LoadParameter("IPPenaltyCoeff",
                                  m_IPPenaltyCoeff,   4.0); 

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
            m_traceAver = Array<OneD, Array<OneD, NekDouble> >(nVariable);
            m_traceJump = Array<OneD, Array<OneD, NekDouble> >(nVariable);
            for(i = 0; i < nVariable; ++i)
            {
                m_traceAver[i] = Array<OneD, NekDouble> (nTracePts,0.0);
                m_traceJump[i] = Array<OneD, NekDouble> (nTracePts,0.0);
            }

            pFields[0]->GetTrace()->GetNormals(m_traceNormals);
            Array<OneD, NekDouble>  lengthFwd(nTracePts,0.0);
            Array<OneD, NekDouble>  lengthBwd(nTracePts,0.0);
            pFields[0]->GetTrace()->GetElmtNormalLength(lengthFwd,lengthBwd);

            const MultiRegions::AssemblyMapDGSharedPtr TraceMap=pFields[0]->GetTraceMap();
            pFields[0]->PeriodicBwdCopy(lengthFwd,lengthBwd);
            TraceMap->UniversalTraceAssemble(lengthFwd);
            TraceMap->UniversalTraceAssemble(lengthBwd);

            // for(int i=0;i<nTracePts;i++)
            // {
            //     if(lengthBwd[i]<0.0)
            //     {
            //         lengthBwd[i] = lengthFwd[i];
            //     }
            // }
            Array<OneD, NekDouble>  lengthsum(nTracePts,0.0);
            Array<OneD, NekDouble>  lengthmul(nTracePts,0.0);
            Vmath::Vadd(nTracePts,lengthBwd,1,lengthFwd,1,lengthsum,1);
            Vmath::Vmul(nTracePts,lengthBwd,1,lengthFwd,1,lengthmul,1);
            Vmath::Vdiv(nTracePts,lengthsum,1,lengthmul,1,lengthFwd,1);
            m_traceNormDirctnElmtLength = lengthsum;
            m_oIPPenaltyLength =   lengthFwd;
            Vmath::Smul(nTracePts,0.25,m_oIPPenaltyLength,1,m_oIPPenaltyLength,1);

            m_tracBwdWeightAver  =   Array<OneD, NekDouble> (nTracePts,0.0);
            m_tracBwdWeightJump  =   Array<OneD, NekDouble> (nTracePts,0.0);
            pFields[0]->GetBwdWeight(m_tracBwdWeightAver,m_tracBwdWeightJump);
            Array<OneD, NekDouble> tmpBwdWeight(nTracePts,0.0);
            Array<OneD, NekDouble> tmpBwdWeightJump(nTracePts,0.0);
            for(int i =1; i<nVariable;i++)
            {
                pFields[i]->GetBwdWeight(tmpBwdWeight,tmpBwdWeightJump);
                Vmath::Vsub(nTracePts,tmpBwdWeight,1,m_tracBwdWeightAver,1,tmpBwdWeight,1);
                Vmath::Vabs(nTracePts,tmpBwdWeight,1,tmpBwdWeight,1);
                NekDouble norm = 0.0;
                for(int j = 0; j<nTracePts; j++)
                {
                    norm += tmpBwdWeight[j];
                }
                ASSERTL0(norm<1.0E-11,"different BWD for different variable not coded yet");
            }

            m_MuVarTrace   =   NullNekDouble1DArray;
            if (m_ArtificialDiffusionVector)
            {
                m_MuVarTrace  =   Array<OneD, NekDouble>(nTracePts, 0.0);
            }

#ifdef CFS_DEBUGMODE
            m_session->LoadParameter("DebugVolTraceSwitch",                 m_DebugVolTraceSwitch      ,    0);
            m_session->LoadParameter("DebugIP_DDGSwitch",                   m_DebugIP_DDGSwitch      ,    0);
#endif
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
            int i, j;
            int nDim      = fields[0]->GetCoordim(0);
            int nPts      = fields[0]->GetTotPoints();
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();

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

            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > qfield(nDim);
            for (j = 0; j < nDim; ++j)
            {
                qfield[j]       = Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    qfield[j][i] = Array<OneD, NekDouble>(nPts, 0.0);
                }
            }
            DiffuseCalculateDerivative(nConvectiveFields,fields,inarray,qfield,vFwd,vBwd);

            Array<OneD, int > nonZeroIndex;
            Diffuse_coeff(nConvectiveFields, fields, inarray, outarray, vFwd, vBwd,qfield,nonZeroIndex);

            for(i = 0; i < nonZeroIndex.num_elements(); ++i)
            {
                int j = nonZeroIndex[i];

                fields[j]->MultiplyByElmtInvMass(outarray[j], outarray[j]);
            }
        }

        void DiffusionIP::v_Diffuse_coeff(
            const int                                                   nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>           &fields,
            const Array<OneD, Array<OneD, NekDouble> >                  &inarray,
            Array<OneD, Array<OneD, NekDouble> >                        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >                  &vFwd,
            const Array<OneD, Array<OneD, NekDouble> >                  &vBwd,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >          &qfield,
            Array< OneD, int >                                          &nonZeroIndex)
        {
            int i, j;
            int nDim      = fields[0]->GetCoordim(0);
            int nPts      = fields[0]->GetTotPoints();
            int nCoeffs   = fields[0]->GetNcoeffs();
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();

            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > elmtFlux(nDim);
            for (j = 0; j < nDim; ++j)
            {
                elmtFlux[j]     = Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    elmtFlux[j][i]   = Array<OneD, NekDouble>(nPts, 0.0);
                }
            }
#ifdef CFS_DEBUGMODE
            if(2!=m_DebugVolTraceSwitch)
            {
#endif
            DiffuseVolumeFlux(nConvectiveFields,fields,inarray,qfield,elmtFlux,nonZeroIndex);
#ifdef CFS_DEBUGMODE
            }
#endif

            //TODO: TO GET TRACE QFIELD FIRST AND RELEASE qfield. AddDiffusionSymmFluxToCoeff DON'T NEED qfield
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

            Array<OneD, Array<OneD, NekDouble > > Traceflux(nConvectiveFields);
            for (int j = 0; j < nConvectiveFields; ++j)
            {
                Traceflux[j]   = Array<OneD, NekDouble>(nTracePts, 0.0);
            }
#ifdef CFS_DEBUGMODE
            if(1!=m_DebugVolTraceSwitch)
            {
#endif
            DiffuseTraceFlux(nConvectiveFields,fields,inarray,qfield,elmtFlux,Traceflux,vFwd,vBwd,nonZeroIndex);
#ifdef CFS_DEBUGMODE
            }
#endif
            // release qfield, elmtFlux and muvar;
            for (j = 0; j < nDim; ++j)
            {
                elmtFlux[j]     = NullNekDoubleArrayofArray;
            }

            for(i = 0; i < nonZeroIndex.num_elements(); ++i)
            {
                int j = nonZeroIndex[i];

                fields[j]->AddTraceIntegral     (Traceflux[j], outarray[j]);
                fields[j]->SetPhysState         (false);
            }
#ifdef CFS_DEBUGMODE
            if(1!=m_DebugVolTraceSwitch)
            {
#endif
            AddDiffusionSymmFluxToCoeff(nConvectiveFields, fields, inarray,qfield,elmtFlux, outarray, vFwd, vBwd);
#ifdef CFS_DEBUGMODE
            }
#endif
        }

        void DiffusionIP::v_DiffuseCalculateDerivative(
            const int                                                   nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>           &fields,
            const Array<OneD, Array<OneD, NekDouble> >                  &inarray,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    &qfield,
            const Array<OneD, Array<OneD, NekDouble> >                  &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >                  &pBwd)
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

        void DiffusionIP::v_DiffuseVolumeFlux(
            const int                                           nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
            const Array<OneD, Array<OneD, NekDouble>>           &inarray,
            Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &qfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
            Array< OneD, int >                                  &nonZeroIndex) 
        {
            int nDim      = fields[0]->GetCoordim(0);
            int nPts      = fields[0]->GetTotPoints();

            Array<OneD, NekDouble> muvar        =   NullNekDouble1DArray;
            if (m_ArtificialDiffusionVector)
            {
                muvar       =   Array<OneD, NekDouble>(nPts, 0.0);
                GetAVmu(fields,inarray,muvar,m_MuVarTrace);
            }

            Array<OneD, Array<OneD, NekDouble> > tmparray2D = NullNekDoubleArrayofArray;

            // TODO: qfield AND elmtFlux share storage????
            m_FunctorDiffusionfluxCons(nConvectiveFields,nDim,inarray,qfield, VolumeFlux,nonZeroIndex,tmparray2D,muvar);
        }
            
        void DiffusionIP::v_DiffuseTraceFlux(
            const int                                           nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
            const Array<OneD, Array<OneD, NekDouble>>           &inarray,
            Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &qfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
            Array<OneD, Array<OneD, NekDouble> >                &TraceFlux,
            const Array<OneD, Array<OneD, NekDouble>>           &pFwd,
            const Array<OneD, Array<OneD, NekDouble>>           &pBwd,
            Array< OneD, int >                                  &nonZeroIndex)
        {
            int nDim      = fields[0]->GetCoordim(0);
            int nPts      = fields[0]->GetTotPoints();
            // int nCoeffs   = fields[0]->GetNcoeffs();
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();

            Array<OneD, Array<OneD, Array<OneD, NekDouble > > > traceflux3D(1);
            traceflux3D[0]  =   TraceFlux;

            const MultiRegions::AssemblyMapDGSharedPtr  TraceMap=fields[0]->GetTraceMap();
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    qBwd(nDim);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    qFwd(nDim);
            for (int nd = 0; nd < nDim; ++nd)
            {
                qBwd[nd]     =   Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
                qFwd[nd]     =   Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
                for (int i = 0; i < nConvectiveFields; ++i)
                {
                    qBwd[nd][i]    = Array<OneD, NekDouble>(nTracePts,0.0);
                    qFwd[nd][i]    = Array<OneD, NekDouble>(nTracePts,0.0);

                    fields[i]->GetFwdBwdTracePhysDeriv_serial(nd,qfield[nd][i], qFwd[nd][i], qBwd[nd][i]);
                    TraceMap->UniversalTraceAssemble(qBwd[nd][i]);
                    TraceMap->UniversalTraceAssemble(qFwd[nd][i]);
                }
            }

            CalTraceNumFlux_ReduceComm(
                nConvectiveFields, nDim, nPts, nTracePts, m_IP2ndDervCoeff,
                fields, inarray, qfield, pFwd, pBwd, qFwd, qBwd, m_MuVarTrace,
                nonZeroIndex, traceflux3D, m_traceAver, m_traceJump);
            
            ApplyFluxBndConds(nConvectiveFields,fields,TraceFlux);

        }

        void DiffusionIP::v_DiffuseTraceFlux(
            const int                                                       nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>               &fields,
            const Array<OneD, Array<OneD, NekDouble>>                       &inarray,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >  &qfield,
            Array<OneD, Array<OneD, NekDouble> >                            &TraceFlux,
            const Array<OneD, Array<OneD, NekDouble>>                       &pFwd,
            const Array<OneD, Array<OneD, NekDouble>>                       &pBwd,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >  &qFwd,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >  &qBwd,
            const Array<OneD, NekDouble>                                    &MuAVTrace,
            Array< OneD, int >                                              &nonZeroIndex  ,
            const Array<OneD, Array<OneD, NekDouble>>                       &Aver          ,
            const Array<OneD, Array<OneD, NekDouble>>                       &Jump          )
        {
            int nDim      = fields[0]->GetCoordim(0);
            int nPts      = fields[0]->GetTotPoints();
            // int nCoeffs   = fields[0]->GetNcoeffs();
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();

            Array<OneD, Array<OneD, Array<OneD, NekDouble > > > traceflux3D(1);
            traceflux3D[0]  =   TraceFlux;

            Array<OneD, Array<OneD, NekDouble> >           pAver;
            Array<OneD, Array<OneD, NekDouble> >           pJump;
            if((Aver.num_elements()&&Jump.num_elements()))
            {
                pAver = Aver;
                pJump = Jump;
            }
            else
            {
                pAver   =   Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
                pJump   =   Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
                for(int i = 0;i<nConvectiveFields;i++)
                {
                    pAver[i]    =   Array<OneD, NekDouble> (nTracePts,0.0);
                    pJump[i]    =   Array<OneD, NekDouble> (nTracePts,0.0);
                }
                // ConsVarAveJump(nConvectiveFields,nTracePts,pFwd,pBwd,pAver,pJump);
            }

            CalTraceNumFlux_ReduceComm(
                nConvectiveFields, nDim, nPts, nTracePts, m_IP2ndDervCoeff,
                fields, inarray, qfield, pFwd, pBwd, qFwd, qBwd, m_MuVarTrace,
                nonZeroIndex, traceflux3D, pAver, pJump);
                
            ApplyFluxBndConds(nConvectiveFields,fields,TraceFlux);
        }

        void DiffusionIP::v_AddDiffusionSymmFluxToCoeff(
            const int                                           nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
            const Array<OneD, Array<OneD, NekDouble> >          &inarray,
            Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &qfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
            Array<OneD, Array<OneD, NekDouble> >                &outarray,
            const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >          &pBwd)
        {
            if(abs(m_IPSymmFtluxCoeff)>1.0E-12)
            {
                int nDim      = fields[0]->GetCoordim(0);
                int nPts      = fields[0]->GetTotPoints();
                int nTracePts = fields[0]->GetTrace()->GetTotPoints();
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > traceSymflux(nDim);
                for (int nd = 0; nd < nDim; ++nd)
                {
                    traceSymflux[nd]    = Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
                    for (int j = 0; j < nConvectiveFields; ++j)
                    {
                        traceSymflux[nd][j]   = Array<OneD, NekDouble>(nTracePts, 0.0);
                    }
                }
                Array< OneD, int >  nonZeroIndex;
                DiffuseTraceSymmFlux(nConvectiveFields,fields,inarray,qfield,VolumeFlux,
                                        traceSymflux,pFwd,pBwd,nonZeroIndex,m_traceAver,m_traceJump);
              
                AddSymmFluxIntegralToCoeff(nConvectiveFields,nDim,nPts,nTracePts,fields,nonZeroIndex,traceSymflux,outarray);
            }
        }

        // void DiffusionIP::v_AddDiffusionSymmFluxToPhys(
        //     const int                                           nConvectiveFields,
        //     const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
        //     const Array<OneD, Array<OneD, NekDouble> >          &inarray,
        //     Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &qfield,
        //     Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
        //     Array<OneD, Array<OneD, NekDouble> >                &outarray,
        //     const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
        //     const Array<OneD, Array<OneD, NekDouble> >          &pBwd)
        // {
            
        //     if(abs(m_IPSymmFtluxCoeff)>1.0E-12)
        //     {
        //         int nDim      = fields[0]->GetCoordim(0);
        //         int nPts      = fields[0]->GetTotPoints();
        //         int nTracePts = fields[0]->GetTrace()->GetTotPoints();
        //         Array<OneD, Array<OneD, Array<OneD, NekDouble> > > traceSymflux(nDim);
        //         for (int nd = 0; nd < nDim; ++nd)
        //         {
        //             traceSymflux[nd]    = Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
        //             for (int j = 0; j < nConvectiveFields; ++j)
        //             {
        //                 traceSymflux[nd][j]   = Array<OneD, NekDouble>(nTracePts, 0.0);
        //             }
        //         }
        //         Array< OneD, int >  nonZeroIndex;
        //         DiffuseTraceSymmFlux(nConvectiveFields,fields,inarray,qfield,VolumeFlux,traceSymflux,pFwd,pBwd,nonZeroIndex);

        //         AddSymmFluxIntegralToPhys(nConvectiveFields,nDim,nPts,nTracePts,fields,nonZeroIndex,traceSymflux,outarray);
        //     }
        // }
        
        void DiffusionIP::v_DiffuseTraceSymmFlux(
            const int                                                   nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>           &fields,
            const Array<OneD, Array<OneD, NekDouble>>                   &inarray,
            const Array<OneD,Array<OneD, Array<OneD, NekDouble> > >     &qfield,
            const Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    &VolumeFlux,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >          &SymmFlux,
            const Array<OneD, Array<OneD, NekDouble>>                   &pFwd,
            const Array<OneD, Array<OneD, NekDouble>>                   &pBwd,
            Array< OneD, int >                                          &nonZeroIndex,
            Array<OneD, Array<OneD, NekDouble> >                        &solution_Aver,
            Array<OneD, Array<OneD, NekDouble> >                        &solution_jump)
        {
            int nDim      = fields[0]->GetCoordim(0);
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();

            Array<OneD, Array<OneD, NekDouble> >                pAver;
            Array<OneD, Array<OneD, NekDouble> >                pjump;
            if(0==solution_jump.num_elements())
            {
                pAver = Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                pjump = Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                for(int m=0;m<nConvectiveFields;m++)
                {
                    pAver[m] = Array<OneD, NekDouble>(nTracePts);
                    pjump[m] = Array<OneD, NekDouble>(nTracePts);
                }
            }
            else
            {
                pAver   =   solution_Aver;
                pjump   =   solution_jump;
            }

            CalTraceSymFlux(nConvectiveFields,nDim,fields,pAver,pjump,
                        nonZeroIndex,SymmFlux);
            for (int nd = 0; nd < nDim; ++nd)
            {
                for (int j = 0; j < nonZeroIndex.num_elements(); ++j)
                {
                    int i = nonZeroIndex[j];
                    Vmath::Smul(nTracePts,-0.5*m_IPSymmFtluxCoeff,SymmFlux[nd][i],1,SymmFlux[nd][i],1);
                }
            }
        }

        void DiffusionIP::CalTraceSymFlux(
            const int                                                           nConvectiveFields,
            const int                                                           nDim,
            const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
            const Array<OneD, Array<OneD, NekDouble> >                          &solution_Aver,
            const Array<OneD, Array<OneD, NekDouble> >                          &solution_jump,
                  Array<OneD, int >                                             &nonZeroIndexsymm,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &traceSymflux)
        {
            int nTracePts = solution_jump[nConvectiveFields-1].num_elements();

            m_FunctorSymmetricfluxCons(nConvectiveFields,nDim,solution_Aver,solution_jump,traceSymflux,nonZeroIndexsymm,m_traceNormals);
        }

        void DiffusionIP::AddSymmFluxIntegralToCoeff(
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
                MultiRegions::ExpListSharedPtr tracelist = fields[nv]->GetTrace();
                for(int nd=0;nd<nDim;nd++)
                {
                    Vmath::Zero(nPts,tmpfield[nd],1);

                    tracelist->MultiplyByQuadratureMetric(tracflux[nd][nv],tracflux[nd][nv]);

                    fields[nv]->AddTraceQuadPhysToField(tracflux[nd][nv],tracflux[nd][nv],tmpfield[nd]);
                    fields[nv]->DividByQuadratureMetric(tmpfield[nd],tmpfield[nd]);
                }
                fields[nv]->IProductWRTDerivBase(tmpfield,tmpCoeff);
                Vmath::Vadd(nCoeffs,tmpCoeff,1,outarray[nv],1,outarray[nv],1);
            }
        }

        void DiffusionIP::AddSymmFluxIntegralToPhys(
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
            Array<OneD, NekDouble > tmpPhysi(nPts,0.0);
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
                fields[nv]->BwdTrans            (tmpCoeff,tmpPhysi);
                Vmath::Vadd(nPts,tmpPhysi,1,outarray[nv],1,outarray[nv],1);
            }
        }
        
        void DiffusionIP::GetPenaltyFactor(
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                  Array<OneD, NekDouble >                       &factor)
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

            int ntmp,numModes;

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
                        numModes    =   0;  
                        for(int nd=0;nd<spaceDim;nd++)
                        {
                            ntmp        = fields[0]->GetExp(LRAdjExpid[nlr][ntrace])->GetBasisNumModes(nd);  
                            numModes    = max(ntmp,numModes);
                        }
                        factorFwdBwd[nlr]   = (numModes)*(numModes);
                    }
                }

                for(int np = 0; np < nTracPnt; ++np)
                {
                    factor[noffset+np]    =   max(factorFwdBwd[0],factorFwdBwd[1]);
                }
            }
        }

        void DiffusionIP::GetPenaltyFactor_const(
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                  Array<OneD, NekDouble >                       &factor)
        {
            Vmath::Fill(factor.num_elements(),m_IPPenaltyCoeff,factor,1);
        }

        void DiffusionIP::v_ConsVarAveJump(
            const int                                           nConvectiveFields,
            const int                                           npnts,
            const Array<OneD, const Array<OneD, NekDouble> >    &vFwd,
            const Array<OneD, const Array<OneD, NekDouble> >    &vBwd,
                  Array<OneD,       Array<OneD, NekDouble> >    &aver,
                  Array<OneD,       Array<OneD, NekDouble> >    &jump)
        {
            ConsVarAve(nConvectiveFields,npnts,vFwd,vBwd,aver);

            m_SpecialBndTreat(nConvectiveFields,aver);

            // note: here the jump is 2.0*(aver-vFwd) 
            //       because Viscous wall use a symmetry value as the Bwd, not the target one   
            Array<OneD, NekDouble> tmpF (npnts,0.0);
            Array<OneD, NekDouble> tmpB (npnts,0.0);

            Array<OneD, NekDouble> Fweight (npnts,2.0);
            Array<OneD, NekDouble> Bweight;
            Bweight = m_tracBwdWeightJump;
            Vmath::Vsub(npnts,Fweight,1,Bweight,1,Fweight,1);

            for (int i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Vsub(npnts,aver[i],1,vFwd[i],1,tmpF,1);
                Vmath::Vsub(npnts,vBwd[i],1,aver[i],1,tmpB,1);

                Vmath::Vmul(npnts,tmpF,1,Fweight,1,tmpF,1);
                Vmath::Vmul(npnts,tmpB,1,Bweight,1,tmpB,1);
                Vmath::Vadd(npnts,tmpF,1,tmpB,1,jump[i],1);
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

            Bweight = m_tracBwdWeightAver;

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
       

        /*1*
     * @brief aplly Neuman boundary conditions on flux 
     *        Currently only consider WallAdiabatic
     *
     */
        void DiffusionIP::ApplyFluxBndConds(
            const int                                               nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>       &fields,
            Array<OneD,       Array<OneD, NekDouble> >              &flux)
        {            
            int ndens       = 0;
            int nengy       = nConvectiveFields-1;
            int nvelst      = ndens + 1;
            int nveled      = nengy;
            
            int cnt;
            int j, e;
            int id2;

            int nBndEdgePts, nBndEdges, nBndRegions;

            int nLengthArray    =0;

            // Compute boundary conditions  for Energy
            cnt = 0;
            nBndRegions = fields[nengy]->
            GetBndCondExpansions().num_elements();
            for (j = 0; j < nBndRegions; ++j)
            {
                if (fields[nengy]->GetBndConditions()[j]->
                    GetBoundaryConditionType() ==
                    SpatialDomains::ePeriodic)
                {
                    continue;
                }

                nBndEdges = fields[nengy]->
                GetBndCondExpansions()[j]->GetExpSize();
                for (e = 0; e < nBndEdges; ++e)
                {
                    nBndEdgePts = fields[nengy]->
                    GetBndCondExpansions()[j]->GetExp(e)->GetTotPoints();

                    id2 = fields[0]->GetTrace()->
                    GetPhys_Offset(fields[0]->GetTraceMap()->
                                GetBndCondTraceToGlobalTraceMap(cnt++));

                    // Imposing Temperature Twall at the wall 
                    if (boost::iequals(fields[nengy]->GetBndConditions()[j]->
                        GetUserDefined(),"WallAdiabatic"))
                    {
                        Vmath::Zero(nBndEdgePts, &flux[nengy][id2], 1);
                    }                    
                }
            }
        }
        
        void DiffusionIP::CalTraceNumFlux_ReduceComm(
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
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &qFwd,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &qBwd,
            const Array<OneD, NekDouble >                                       &MuVarTrace,
                  Array<OneD, int >                                             &nonZeroIndexflux,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &traceflux,
                  Array<OneD, Array<OneD, NekDouble> >                          &solution_Aver,
                  Array<OneD, Array<OneD, NekDouble> >                          &solution_jump)
        {
            const MultiRegions::AssemblyMapDGSharedPtr                      TraceMap=fields[0]->GetTraceMap();

            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    numDerivBwd(nDim);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    numDerivFwd(nDim);
            for (int nd = 0; nd < nDim; ++nd)
            {
                numDerivBwd[nd]     =   Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
                numDerivFwd[nd]     =   Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
                for (int i = 0; i < nConvectiveFields; ++i)
                {
                    numDerivBwd[nd][i]    = Array<OneD, NekDouble>(nTracePts,0.0);
                    numDerivFwd[nd][i]    = Array<OneD, NekDouble>(nTracePts,0.0);
                }
            }

            if(abs(PenaltyFactor2)>1.0E-12)
            {
                AddSecondDerivTOTrace_ReduceComm(nConvectiveFields,nDim,nPts,nTracePts,PenaltyFactor2,fields,qfield,numDerivFwd,numDerivBwd);
                for (int nd = 0; nd < nDim; ++nd)
                {
                    for (int i = 0; i < nConvectiveFields; ++i)
                    {
                        Vmath::Svtvp(nTracePts,0.5,qBwd[nd][i],1,numDerivBwd[nd][i],1,numDerivBwd[nd][i],1);
                        Vmath::Svtvp(nTracePts,0.5,qFwd[nd][i],1,numDerivFwd[nd][i],1,numDerivFwd[nd][i],1);
                        TraceMap->UniversalTraceAssemble(numDerivBwd[nd][i]);
                        TraceMap->UniversalTraceAssemble(numDerivFwd[nd][i]);
                        Vmath::Vadd(nTracePts,numDerivFwd[nd][i],1,numDerivBwd[nd][i],1,numDerivFwd[nd][i],1);
                        numDerivBwd[nd][i]    = NullNekDouble1DArray;
                    }
                }
            }
            else
            {
                for (int nd = 0; nd < nDim; ++nd)
                {
                    for (int i = 0; i < nConvectiveFields; ++i)
                    {
                        Vmath::Svtvp(nTracePts,0.5,qBwd[nd][i],1,numDerivBwd[nd][i],1,numDerivBwd[nd][i],1);
                        Vmath::Svtvp(nTracePts,0.5,qFwd[nd][i],1,numDerivFwd[nd][i],1,numDerivFwd[nd][i],1);
                        Vmath::Vadd(nTracePts,numDerivFwd[nd][i],1,numDerivBwd[nd][i],1,numDerivFwd[nd][i],1);
                    }
                }
            }

            for (int nd = 0; nd < nDim; ++nd)
            {
                for (int i = 0; i < nConvectiveFields; ++i)
                {
                    numDerivBwd[nd][i]    = NullNekDouble1DArray;
                }
            }

            ConsVarAveJump(nConvectiveFields,nTracePts,vFwd,vBwd,solution_Aver,solution_jump);

            Array<OneD, NekDouble> jumpTmp(nTracePts,0.0);
            Array<OneD, NekDouble> PenaltyFactor(nTracePts,0.0);

            // GetPenaltyFactor_const(fields,PenaltyFactor);
            GetPenaltyFactor(fields,PenaltyFactor);

            Vmath::Vmul(nTracePts,PenaltyFactor,1, m_oIPPenaltyLength,1,PenaltyFactor,1);
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Vmul(nTracePts,solution_jump[i],1, PenaltyFactor,1,jumpTmp,1);
                for (int nd = 0; nd < nDim; ++nd)
                {
                    Vmath::Vvtvp(nTracePts, m_traceNormals[nd],1,jumpTmp,1, numDerivFwd[nd][i],1, numDerivFwd[nd][i],1);
                }
            }
            jumpTmp         =   NullNekDouble1DArray;
            PenaltyFactor   =   NullNekDouble1DArray;
            // Calculate normal viscous flux
            m_FunctorDiffusionfluxCons(nConvectiveFields,nDim,solution_Aver,numDerivFwd,traceflux,nonZeroIndexflux,m_traceNormals,MuVarTrace);
        }
        
        void DiffusionIP::AddSecondDerivTOTrace_ReduceComm(
            const int                                                           nConvectiveFields,
            const int                                                           nDim,
            const int                                                           nPts,
            const int                                                           nTracePts,
            const NekDouble                                                     PenaltyFactor2,
            const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &qfield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &numDerivFwd,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &numDerivBwd)
        {
            Array<OneD, NekDouble> Fwd(nTracePts,0.0);
            Array<OneD, NekDouble> Bwd(nTracePts,0.0);
            Array<OneD,NekDouble>  tmp(nTracePts,0.0);

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

            Vmath::Smul(nTracePts,PenaltyFactor2,m_traceNormDirctnElmtLength,1,tmp,1);
            // the derivatives are assumed to be exchangable  
            for(int nd1=0; nd1<nDim; nd1++)
            {
                for(int i = 0; i < nConvectiveFields; ++i)
                {
                    fields[i]->PhysDeriv(qfield[nd1][i], qtmp[0], qtmp[1], qtmp[2]);

                    for(int nd2=nd1; nd2<nDim; nd2++)
                    {
                        Vmath::Zero(nTracePts,Bwd,1);
                        fields[i]->GetFwdBwdTracePhysDeriv_serial(nd2,elmt2ndDerv[nd2], Fwd, Bwd);
                        Vmath::Vmul(nTracePts,tmp,1,Bwd,1,Bwd,1);
                        Vmath::Vvtvp(nTracePts,m_traceNormals[nd2],1,Bwd,1,numDerivBwd[nd1][i],1,numDerivBwd[nd1][i],1);
                        Vmath::Vmul(nTracePts,tmp,1,Fwd,1,Fwd,1);
                        Vmath::Vvtvm(nTracePts,m_traceNormals[nd2],1,Fwd,1,numDerivFwd[nd1][i],1,numDerivFwd[nd1][i],1);
                        Vmath::Neg(nTracePts,numDerivFwd[nd1][i],1);

                        if(nd2!=nd1)
                        {
                            Vmath::Vvtvp(nTracePts,m_traceNormals[nd1],1,Bwd,1,numDerivBwd[nd2][i],1,numDerivBwd[nd2][i],1);
                            Vmath::Vvtvm(nTracePts,m_traceNormals[nd1],1,Fwd,1,numDerivFwd[nd2][i],1,numDerivFwd[nd2][i],1);
                            Vmath::Neg(nTracePts,numDerivFwd[nd2][i],1);
                        }
                    }
                }
            }
        }

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
        void DiffusionIP::v_MinusVolumDerivJacToMat( 
            const int                                                   nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>           &pFields,
            const Array<OneD, const Array<OneD,  Array<OneD, 
                Array<OneD,  Array<OneD,  NekDouble> > > > >            &ElmtJacArray,
            const int                                                   nDervDir, 
            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >              &gmtxarray)
        {
            MultiRegions::ExpListSharedPtr explist = pFields[0];
                std::shared_ptr<LocalRegions::ExpansionVector> pexp = explist->GetExp();
            int ntotElmt            = (*pexp).size();
            int nElmtPnt,nElmtCoef;

            NekDouble tmp;
            DNekMatSharedPtr        tmpGmtx,ElmtMat;

            Array<OneD, NekDouble>  GlobMat_data;
            Array<OneD, NekDouble>  ElmtMat_data;

            Array<OneD, DNekMatSharedPtr>  mtxPerVar(ntotElmt);
            Array<OneD, DNekMatSharedPtr>  mtxPerVarCoeff(ntotElmt);
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
                (*mtxPerVar[nelmt])    = 0.0;       
                mtxPerVarCoeff[nelmt]    =MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nElmtCoef, nElmtCoef);
                (*mtxPerVarCoeff[nelmt])   =   0.0;
            }

            for(int m = 0; m < nConvectiveFields; m++)
            {
                for(int n = 0; n < nConvectiveFields; n++)
                {

                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        (*mtxPerVarCoeff[nelmt])   =   0.0;
                        (*mtxPerVar[nelmt])   =   0.0;
                    }
                    explist->GetMatIpwrtDeriveBase(ElmtJacArray[m][n],mtxPerVar);
                    //TODO: To check whether it is ok to reuse ElmtJacQuad as output
                    explist->AddRightIPTPhysDerivBase(nDervDir,mtxPerVar,mtxPerVarCoeff);

                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        nElmtCoef       = elmtcoef[nelmt];
                        nElmtPnt        = elmtpnts[nelmt];

                        tmpGmtx         = gmtxarray[m][n]->GetBlock(nelmt,nelmt);
                        ElmtMat         = mtxPerVarCoeff[nelmt];

                        GlobMat_data    = tmpGmtx->GetPtr();
                        ElmtMat_data    = ElmtMat->GetPtr();

                        Vmath::Vsub(nElmtCoef*nElmtCoef,GlobMat_data,1,ElmtMat_data,1,GlobMat_data,1);
                    }
                }
            }
        }

        void DiffusionIP::v_MinusVolumDerivJacToMat( 
            const int                                                   nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>           &pFields,
            const Array<OneD, const Array<OneD,  Array<OneD, 
                Array<OneD,  Array<OneD,  NekDouble> > > > >            &ElmtJacArray,
            const int                                                   nDervDir, 
            Array<OneD, Array<OneD, SNekBlkMatSharedPtr> >              &gmtxarray)
        {
            MultiRegions::ExpListSharedPtr explist = pFields[0];
                std::shared_ptr<LocalRegions::ExpansionVector> pexp = explist->GetExp();
            int ntotElmt            = (*pexp).size();
            int nElmtPnt,nElmtCoef;

            NekDouble tmp;
            SNekMatSharedPtr        tmpGmtx;
            DNekMatSharedPtr        ElmtMat;

            Array<OneD, NekSingle>  GlobMat_data;
            Array<OneD, NekDouble>  ElmtMat_data;
            Array<OneD,NekSingle> Elmt_dataSingle;

            Array<OneD, DNekMatSharedPtr>  mtxPerVar(ntotElmt);
            Array<OneD, DNekMatSharedPtr>  mtxPerVarCoeff(ntotElmt);
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
                (*mtxPerVar[nelmt])    = 0.0;       
                mtxPerVarCoeff[nelmt]    =MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nElmtCoef, nElmtCoef);
                (*mtxPerVarCoeff[nelmt])   =   0.0;
            }

            for(int m = 0; m < nConvectiveFields; m++)
            {
                for(int n = 0; n < nConvectiveFields; n++)
                {

                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        (*mtxPerVarCoeff[nelmt])   =   0.0;
                        (*mtxPerVar[nelmt])   =   0.0;
                    }
                    explist->GetMatIpwrtDeriveBase(ElmtJacArray[m][n],mtxPerVar);
                    //TODO: To check whether it is ok to reuse ElmtJacQuad as output
                    explist->AddRightIPTPhysDerivBase(nDervDir,mtxPerVar,mtxPerVarCoeff);

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

                        GlobMat_data    = tmpGmtx->GetPtr();
                        ElmtMat_data    = ElmtMat->GetPtr();

                        for(int i=0;i<ntotDofs;i++)
                        {
                            Elmt_dataSingle[i]  =   NekSingle( ElmtMat_data[i] );
                        }

                        Vmath::Vsub(ntotDofs,GlobMat_data,1,Elmt_dataSingle,1,GlobMat_data,1);
                    }
                }
            }
        }

        // void DiffusionIP::v_DiffuseTraceFlux(
        //     const int                                           nConvectiveFields,
        //     const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
        //     const Array<OneD, Array<OneD, NekDouble>>           &inarray,
        //     Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &qfield,
        //     Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
        //     Array<OneD, Array<OneD, NekDouble> >                &TraceFlux,
        //     const Array<OneD, Array<OneD, NekDouble>>           &pFwd,
        //     const Array<OneD, Array<OneD, NekDouble>>           &pBwd,
        //     Array< OneD, int >                                  &nonZeroIndex)
        // {
        //     int nDim      = fields[0]->GetCoordim(0);
        //     int nPts      = fields[0]->GetTotPoints();
        //     // int nCoeffs   = fields[0]->GetNcoeffs();
        //     int nTracePts = fields[0]->GetTrace()->GetTotPoints();

        //     Array<OneD, Array<OneD, Array<OneD, NekDouble > > > traceflux3D(1);
        //     traceflux3D[0]  =   TraceFlux;

        //     const MultiRegions::AssemblyMapDGSharedPtr  TraceMap=fields[0]->GetTraceMap();
        //     Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    qBwd(nDim);
        //     Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    qFwd(nDim);
        //     for (int nd = 0; nd < nDim; ++nd)
        //     {
        //         qBwd[nd]     =   Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
        //         qFwd[nd]     =   Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
        //         for (int i = 0; i < nConvectiveFields; ++i)
        //         {
        //             qBwd[nd][i]    = Array<OneD, NekDouble>(nTracePts,0.0);
        //             qFwd[nd][i]    = Array<OneD, NekDouble>(nTracePts,0.0);

        //             fields[i]->GetFwdBwdTracePhysDeriv_serial(nd,VolumeFlux[nd][i], qFwd[nd][i], qBwd[nd][i]);
        //             TraceMap->UniversalTraceAssemble(qBwd[nd][i]);
        //             TraceMap->UniversalTraceAssemble(qFwd[nd][i]);
        //         }
        //     }

        //     CalTraceNumFlux_ReduceComm_Flux(
        //         nConvectiveFields, nDim, nPts, nTracePts, m_IP2ndDervCoeff,
        //         fields, inarray, VolumeFlux, pFwd, pBwd, qFwd, qBwd, m_MuVarTrace,
        //         nonZeroIndex, traceflux3D, m_traceAver, m_traceJump);
        // }

        // void DiffusionIP::v_CalDiffusionSymmFlux(
        //     const int                                           nConvectiveFields,
        //     const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
        //     const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
        //     const Array<OneD, Array<OneD, NekDouble> >          &pBwd,
        //     const Array<OneD, Array<OneD, NekDouble> >          &jump,
        //     Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &flux)
        // {
        //     int nDim      = fields[0]->GetCoordim(0);
        //     int nTracePts = fields[0]->GetTrace()->GetTotPoints();
        //     Array< OneD, int >  nonZeroIndex;
        //     // DiffuseTraceSymmFlux_IP(nConvectiveFields,fields,flux,inarray,traceJumpTmp,nonZeroIndex);
        //     CalTraceSymFlux(nConvectiveFields,nDim,fields,pFwd,flux,nonZeroIndex,flux);
        //     Array<OneD, Array<OneD, Array<OneD, NekDouble> > > tmpFlux(nDim);
        //     for (int nd = 0; nd < nDim; ++nd)
        //     {
        //         tmpFlux[nd] = Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
        //         for (int i = 0; i < nConvectiveFields; ++i)
        //         {
        //             tmpFlux[nd][i]   = Array<OneD, NekDouble>(nTracePts,0.0);
        //         }
        //     }
        //     CalTraceSymFlux(nConvectiveFields,nDim,fields,pBwd,tmpFlux,nonZeroIndex,flux);

        //     for (int nd = 0; nd < nDim; ++nd)
        //     {
        //         for (int j = 0; j < nonZeroIndex.num_elements(); ++j)
        //         {
        //             int i = nonZeroIndex[j];
        //             Vmath::Vadd(nTracePts,tmpFlux[nd][i],1,flux[nd][i],1,flux[nd][i],1);
        //             Vmath::Smul(nTracePts,0.5,flux[nd][i],1,flux[nd][i],1);
        //         }
        //     }
        // }

        // void DiffusionIP::AddPenaltyFluxToTraceFlux(
        //     const int                                           nConvectiveFields,
        //     const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
        //     const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
        //     const Array<OneD, Array<OneD, NekDouble> >          &pBwd,
        //     const Array<OneD, Array<OneD, NekDouble> >          &jump,
        //     Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &flux)
        // {
        //     int nDim      = fields[0]->GetCoordim(0);
        //     int nTracePts = fields[0]->GetTrace()->GetTotPoints();
        //     Array< OneD, int >  nonZeroIndex;
        //     // DiffuseTraceSymmFlux_IP(nConvectiveFields,fields,flux,inarray,traceJumpTmp,nonZeroIndex);
        //     CalTraceSymFlux(nConvectiveFields,nDim,fields,pFwd,flux,nonZeroIndex,flux);
        //     Array<OneD, Array<OneD, Array<OneD, NekDouble> > > tmpFlux(nDim);
        //     for (int nd = 0; nd < nDim; ++nd)
        //     {
        //         tmpFlux[nd] = Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
        //         for (int i = 0; i < nConvectiveFields; ++i)
        //         {
        //             tmpFlux[nd][i]   = Array<OneD, NekDouble>(nTracePts,0.0);
        //         }
        //     }
        //     CalTraceSymFlux(nConvectiveFields,nDim,fields,pBwd,tmpFlux,nonZeroIndex,flux);

        //     for (int nd = 0; nd < nDim; ++nd)
        //     {
        //         for (int j = 0; j < nonZeroIndex.num_elements(); ++j)
        //         {
        //             int i = nonZeroIndex[j];
        //             Vmath::Vadd(nTracePts,tmpFlux[nd][i],1,flux[nd][i],1,flux[nd][i],1);
        //             Vmath::Smul(nTracePts,0.5,flux[nd][i],1,flux[nd][i],1);
        //         }
        //     }
        // }

        // void DiffusionIP::v_AddDiffusionSymmFluxToCoeff(
        //     const int                                           nConvectiveFields,
        //     const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
        //     const Array<OneD, Array<OneD, NekDouble> >          &symflux,
        //     Array<OneD, Array<OneD, NekDouble> >                &outarray,
        //     const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
        //     const Array<OneD, Array<OneD, NekDouble> >          &pBwd)
        // {
        //     if(abs(m_IPSymmFtluxCoeff)>1.0E-12)
        //     {
        //         AddSymmFluxIntegralToCoeff(nConvectiveFields,nDim,nPts,nTracePts,fields,nonZeroIndex,traceSymfluxFwd,outarray);
        //     }
        // }

        // //TODO::WIP
        // void DiffusionIP::CalTraceNumFlux_ReduceComm_Flux(
        //     const int                                                           nConvectiveFields,
        //     const int                                                           nDim,
        //     const int                                                           nPts,
        //     const int                                                           nTracePts,
        //     const NekDouble                                                     PenaltyFactor2,
        //     const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
        //     const Array<OneD, Array<OneD, NekDouble> >                          &inarray,
        //     const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &VolFlux,
        //     const Array<OneD, Array<OneD, NekDouble> >                          &vFwd,
        //     const Array<OneD, Array<OneD, NekDouble> >                          &vBwd,
        //     const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &FluxFwd,
        //     const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &FluxBwd,
        //     const Array<OneD, NekDouble >                                       &MuVarTrace,
        //           Array<OneD, int >                                             &nonZeroIndexflux,
        //           Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &traceflux,
        //           Array<OneD, Array<OneD, NekDouble> >                          &solution_Aver,
        //           Array<OneD, Array<OneD, NekDouble> >                          &solution_jump)
        // {
        //     const MultiRegions::AssemblyMapDGSharedPtr                      TraceMap=fields[0]->GetTraceMap();

        //     Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    numDerivFwd(nDim);
        //     for (int nd = 0; nd < nDim; ++nd)
        //     {
        //         numDerivFwd[nd]     =   Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
        //         for (int i = 0; i < nConvectiveFields; ++i)
        //         {
        //             numDerivFwd[nd][i]    = Array<OneD, NekDouble>(nTracePts,0.0);
        //         }
        //     }

        //     for (int nd = 0; nd < nDim; ++nd)
        //     {
        //         for (int i = 0; i < nConvectiveFields; ++i)
        //         {
        //             Vmath::Vadd(nTracePts,FluxFwd[nd][i],1,FluxBwd[nd][i],1,numDerivFwd[nd][i],1);
        //             Vmath::Vvtvp(nTracePts,m_traceNormals[nd],1,numDerivFwd[nd][i],1,traceflux[0][i],1,traceflux[0][i],1);
        //         }
        //     }

        //     for (int i = 0; i < nConvectiveFields; ++i)
        //     {
        //         Vmath::Smul(nTracePts,0.5,traceflux[0][i],1,traceflux[0][i],1);
        //     }

        //     ConsVarAveJump(nConvectiveFields,nTracePts,vFwd,vBwd,solution_Aver,solution_jump);

        //     Array<OneD, NekDouble> mu(nTracePts,0.0);

        //     m_CalcViscosity(solution_Aver,mu);

        //     Array<OneD, NekDouble> PenaltyFactor(nTracePts,0.0);
        //     // GetPenaltyFactor(fields,PenaltyFactor);
        //     GetPenaltyFactor_const(fields,PenaltyFactor);

        //     Vmath::Vmul(nTracePts,PenaltyFactor,1, m_oIPPenaltyLength,1,PenaltyFactor,1);
        //     Vmath::Vmul(nTracePts,PenaltyFactor,1, mu,1,PenaltyFactor,1);
        //     mu         =   NullNekDouble1DArray;
        //     Array<OneD, Array<OneD, NekDouble> > PenaltyFlux (nConvectiveFields);
        //     for (int i = 0; i < nConvectiveFields; ++i)
        //     {
        //         PenaltyFlux[i]  =   Array<OneD, NekDouble> (nTracePts,0.0);
        //         Vmath::Vmul(nTracePts,solution_jump[i],1, PenaltyFactor,1,PenaltyFlux[i],1);
        //         // Vmath::Vdiv(nTracePts,PenaltyFlux[i],1, solution_Aver[0],1,PenaltyFlux[i],1);
        //     }
        //     PenaltyFactor   =   NullNekDouble1DArray;

        //     for(int j=0;j<nConvectiveFields;j++)
        //     {
        //         Vmath::Vadd(nTracePts,&PenaltyFlux[j][0],1,&traceflux[0][j][0],1,&traceflux[0][j][0],1);
        //     }

        //     ApplyFluxBndConds(nConvectiveFields,fields,traceflux[0]);

        //     int n_nonZero = nConvectiveFields;

        //     nonZeroIndexflux = Array< OneD, int > (n_nonZero,0);
        //     for(int i=1;i<n_nonZero+1; i++)
        //     {
        //         nonZeroIndexflux[n_nonZero-i] =   nConvectiveFields-i;
        //     }
        // }
#endif

    }
}
