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

            // Setting up the normals
            int i;
            int nDim = pFields[0]->GetCoordim(0);
            int nTracePts = pFields[0]->GetTrace()->GetTotPoints();
            
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for(i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts);
            }
            pFields[0]->GetTrace()->GetNormals(m_traceNormals);
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
            int nBndEdgePts, i, j, k, e;
            int nDim      = fields[0]->GetCoordim(0);
            int nPts      = fields[0]->GetTotPoints();
            int nCoeffs   = fields[0]->GetNcoeffs();
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();

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

            Array< Array<OneD, NekDouble> > qtmp(3);
            for(int nd=0; nd<3, nd++)
            {
                qtmp[nd]    =   NullNekDouble1DArray;
            }
            for(int i = 0; i < nConvectiveFields; ++i)
            {
                for(int nd=0; nd<nDim, nd++)
                {
                    qtmp[nd]    =   qfield[nd][i];
                }
                m_fields[i]->PhysDeriv(inarray[i], qtmp[0], qtmp[1], qtmp[2]);
            }

            Array<OneD, NekDouble> muvar        =   NullNekDouble1DArray;
            Array<OneD, NekDouble> FwdMuVar     =   NullNekDouble1DArray;
            Array<OneD, NekDouble> BwdMuVar     =   NullNekDouble1DArray;
            Array<OneD, NekDouble> MuVarTrace   =   NullNekDouble1DArray;
            if (m_ArtificialDiffusionVector)
            {
                muvar   =   Array<OneD, NekDouble>(nPts, 0.0);
                m_ArtificialDiffusionVector(inarray, muvar);

                FwdMuVar =  Array<OneD, NekDouble>(nTracePts, 0.0);
                BwdMuVar =  Array<OneD, NekDouble>(nTracePts, 0.0);

                // TODO: CHECK !!
                fields[0]->GetFwdBwdTracePhysInterior(muvar,FwdMuVar,BwdMuVar);
                fields[0]->FillBwdWITHBoundZero(FwdMuVar,BwdMuVar);

                for(k = 0; k < nTracePts; ++k)
                {
                    FwdMuVar[k] = 0.5 * (FwdMuVar[k] + BwdMuVar[k]) ;
                }
                MuVarTrace  =   FwdMuVar;
                FwdMuVar    =   NullNekDouble1DArray;
                BwdMuVar    =   NullNekDouble1DArray;
            }
            
            Array<OneD, NekDouble>    Fwd   =   NullNekDouble1DArray;
            Array<OneD, NekDouble>    Bwd   =   NullNekDouble1DArray;
            Bwd = Array<OneD, NekDouble>(nTracePts,0.0);
            
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    numDeriv(nDim);
            for (int nd = 0; nd < nDim; ++nd)
            {
                numDeriv[nd]     =   Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    numDeriv[nd][i]    = Array<OneD, NekDouble>(nTracePts,0.0);
                    Fwd =  numDeriv[nd][i]; 
                    Vmath::Fill(nTracePts, 1.0, Bwd,1);
                    fields[i]->GetFwdBwdTracePhysDeriv(qfield[nd][i], Fwd, Bwd);
                    for (int nt = 0; nt < nTracePts; ++nt)
                    {
                        numDeriv[nd][i][nt] =   0.5*( Fwd[nt] + Bwd[nt] );  
                    }
                }
            }
            // release storage
            Bwd   =   NullNekDouble1DArray;
            Array< OneD, int > nonZeroIndex;
            DiffusionFlux(nConvectiveFields,nDim,inarray,qfield, elmtFlux,nonZeroIndex,muvar);
            // volume intergration 
            for(i = 0; i < nonZeroIndex.num_elements(); ++i)
            {
                int j = nonZeroIndex[i];
                fields[j]->IProductWRTDerivBase(elmtFlux[j],outarray[j]);
                Vmath::Neg                      (nCoeffs, outarray[j], 1);
            }
            // release qfield and muvar;
            muvar   =   NullNekDouble1DArray;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > qfield(nDim);
            for (j = 0; j < nDim; ++j)
            {
                qfield[j]       = NullNekDoubleArrayofArray;
                elmtFlux[j]     = NullNekDoubleArrayofArray;
            }
            
            Array<OneD, Array<OneD, NekDouble> >    solution_jump(nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> >    solution_Aver(nConvectiveFields);
            for (i = 0; i < nConvectiveFields; ++i)
            {
                solution_jump[i]    =   Array<OneD, NekDouble>(nTracePts,0.0);
                solution_Aver[i]    =   Array<OneD, NekDouble>(nTracePts,0.0);
            }

            // Careful Fwd and solution_jump, solution_Aver and Bwd share same storage
            if (pFwd == NullNekDoubleArrayofArray ||
                pBwd == NullNekDoubleArrayofArray)
            {
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    Fwd       =    solution_jump[i];  
                    Bwd       =    solution_Aver[i]
                    fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd, Bwd);
                    for (int nt = 0; nt < nTracePts; ++nt)
                    {
                        solution_jump[i][nt]   =   Fwd[nt] - Bwd[nt];  
                        solution_Aver[i][nt]   =   0.5*solution_jump[nt] + Bwd[nt];  
                    }
                }
            }
            else
            {
                Fwd   =   NullNekDouble1DArray;
                Bwd   =   NullNekDouble1DArray;
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    Fwd  = pFwd[i];
                    Bwd  = pBwd[i];
                    for (int nt = 0; nt < nTracePts; ++nt)
                    {
                        solution_jump[i][nt]   =   Fwd[nt] - Bwd[nt];  
                        solution_Aver[i][nt]   =   0.5*solution_jump[nt] + Bwd[nt];  
                    }
                }
            }

            Array<OneD, NekDouble>  PenaltyFactor(nTracePts,0.0);
            GetPenaltyFactor(PenaltyFactor);
            for (i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Vmul(nTracePts,solution_jump[i],1, PenaltyFactor,1,solution_jump[i],1);
            }

            Array<OneD, NekDouble> ElmtLength;
            ElmtLength  =   PenaltyFactor; 
            PenaltyFactor = NullNekDouble1DArray;
            Vmath::Fill(nTracePts, 1.0, ElmtLength,1);
            // getElmtLength(ElmtLength);
            for (i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Vdiv(nTracePts,solution_jump[i],1, ElmtLength,1,solution_jump[i],1);
            }

            for (int nd = 0; nd < nDim; ++nd)
            {
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    Vmath::Vvtvp(nTracePts, m_traceNormals[nd],1,solution_jump[i],1, numDeriv[nd][i],1, numDeriv[nd][i],1);
                }
            }

            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > traceflux(1);
            traceflux[0]    =   Array<OneD, Array<OneD, NekDouble> > traceflux(nConvectiveFields);
            for (int j = 0; j < nConvectiveFields; ++j)
            {
                traceflux[0][j]   = Array<OneD, NekDouble>(nTracePts, 0.0);
            }
            // Calculate normal viscous flux
            DiffusionFlux(nConvectiveFields,nDim,solution_Aver,numDeriv,traceflux,nonZeroIndex,MuVarTrace,m_traceNormals);
            for(i = 0; i < nonZeroIndex.num_elements(); ++i)
            {
                int j = nonZeroIndex[i];
                fields[j]->AddTraceIntegral     (traceflux[j], outarray[j]);
                fields[j]->MultiplyByElmtInvMass(outarray[j], outarray[j]);
            }
        }
    }
}
