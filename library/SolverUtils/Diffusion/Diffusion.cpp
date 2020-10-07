///////////////////////////////////////////////////////////////////////////////
//
// File: Diffusion.cpp
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
// Description: Abstract base class for diffusion.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Diffusion/Diffusion.h>

namespace Nektar
{
    namespace SolverUtils
    {
        DiffusionFactory& GetDiffusionFactory()
        {
            static DiffusionFactory instance;
            return instance;
        }

        void Diffusion::InitObject(
            const LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr>       pFields)
        {
            v_InitObject(pSession, pFields);

            // Div curl storage
            int nPts = pFields[0]->GetTotPoints();
            m_divVel        = Array<OneD, NekDouble>(nPts, 0.0);
            m_divVelSquare  = Array<OneD, NekDouble>(nPts, 0.0);
            m_curlVelSquare = Array<OneD, NekDouble>(nPts, 0.0);
        }

        void Diffusion::Diffuse(
            const std::size_t                                 nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            v_Diffuse(nConvectiveFields, fields, inarray, outarray, pFwd, pBwd);
        }

        /**
         * @brief Similar with Diffusion::Diffuse(): calculate diffusion flux
         * The difference is in the outarray:
         *  it is the coefficients of basis for DiffuseCoeffs()
         *  it is the physics on quadrature points for Diffuse()
         */
        void Diffusion::DiffuseCoeffs(
            const std::size_t                                 nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            v_DiffuseCoeffs(nConvectiveFields, fields, inarray,
                            outarray, pFwd, pBwd);
        }

        void Diffusion::Diffuse(
            const std::size_t                                 nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, Array<OneD, NekDouble> >              &outarray,
            NekDouble                                         time,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            m_time  =    time;
            v_Diffuse(nConvectiveFields, fields, inarray, outarray, pFwd, pBwd);
        }

        void Diffusion::DiffuseCoeffs(
            const std::size_t                             nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, Array<OneD, NekDouble> >              &outarray,
            NekDouble                                         time,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            m_time  =    time;
            v_DiffuseCoeffs(nConvectiveFields, fields, inarray, outarray,
                            pFwd, pBwd);
        }
        void Diffusion::GetAVmu(
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
            const Array<OneD, Array<OneD, NekDouble> >          &inarray,
                  Array<OneD, NekDouble >                       &muvar,
                  Array<OneD, NekDouble >                       &MuVarTrace)
        {
            size_t nTracePts = fields[0]->GetTrace()->GetTotPoints();

            Array<OneD, NekDouble> Fwd{nTracePts, 0.0};
            Array<OneD, NekDouble> Bwd{nTracePts, 0.0};

            m_ArtificialDiffusionVector(inarray, muvar);

            // BwdMuvar is left to be 0.0 according to DiffusionLDG.cpp
            fields[0]->GetFwdBwdTracePhys(muvar, Fwd, Bwd, false);

            for (int k = 0; k < nTracePts; ++k)
            {
                MuVarTrace[k] = 0.5 * (Fwd[k] + Bwd[k]) ;
            }
        }

        void Diffusion::AddSymmFluxIntegralToOffDiag(
            const int                                          nConvectiveFields,
            const int                                             nDim,
            const int                                             nPts,
            const int                                             nTracePts,
            const Array<OneD, MultiRegions::ExpListSharedPtr>     &fields,
            const Array<OneD, const int >                         &nonZeroIndex,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    &Fwdflux,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    &Bwdflux,
            Array<OneD, Array<OneD, NekDouble> >                  &outarray)
        {
            boost::ignore_unused(nTracePts);
            
            int nCoeffs =   outarray[nConvectiveFields-1].size();
            Array<OneD, NekDouble > tmpCoeff(nCoeffs,0.0);
            Array<OneD, Array<OneD, NekDouble> > tmpfield(nDim);

            for(int i = 0;i<nDim;i++)
            {
                tmpfield[i]    =   Array<OneD, NekDouble>(nPts,0.0);
            }
            int nv = 0;

            for(int j=0;j<nonZeroIndex.size();j++)
            {
                nv  =   nonZeroIndex[j];
                MultiRegions::ExpListSharedPtr tracelist = fields[nv]->GetTrace();
                for(int nd=0;nd<nDim;nd++)
                {
                    Vmath::Zero(nPts,tmpfield[nd],1);

                    tracelist->MultiplyByQuadratureMetric(Fwdflux[nd][nv],
                                                          Fwdflux[nd][nv]);
                    tracelist->MultiplyByQuadratureMetric(Bwdflux[nd][nv],
                                                          Bwdflux[nd][nv]);

                    fields[nv]->AddTraceQuadPhysToOffDiag(Fwdflux[nd][nv],
                                                Bwdflux[nd][nv],tmpfield[nd]);
                    fields[nv]->DivideByQuadratureMetric(tmpfield[nd],tmpfield[nd]);
                }
                fields[nv]->IProductWRTDerivBase(tmpfield,tmpCoeff);
                Vmath::Vadd(nCoeffs,tmpCoeff,1,outarray[nv],1,outarray[nv],1);
            }
        }
        
        void Diffusion::v_Diffuse(
            const std::size_t                             nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                    Array<OneD, Array<OneD, NekDouble> >      &outarray,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            boost::ignore_unused(nConvectiveFields, fields, inarray, outarray,
                                    pFwd, pBwd);
            ASSERTL0(false, "v_Diffuse not defined");
        }

        void Diffusion::v_DiffuseCoeffs(
            const std::size_t                                 nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            boost::ignore_unused(nConvectiveFields, fields, inarray,
                                    outarray, pFwd, pBwd);
            ASSERTL0(false, "v_DiffuseCoeffs not defined");
        }

        void Diffusion::v_DiffuseCoeffs(
                const std::size_t                            nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray,
                const Array<OneD, Array<OneD, NekDouble> >        &vFwd,
                const Array<OneD, Array<OneD, NekDouble> >        &vBwd,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &qfield,
                Array< OneD, int >                                 &nonZeroIndex)
        {
            boost::ignore_unused(nConvectiveFields, fields, inarray,
                                 outarray, vFwd, vBwd,qfield,nonZeroIndex);
            ASSERTL0(false,"v_DiffuseCoeffs not defined");
        }

        // Diffusion Calculate the physical derivatives
        void Diffusion::DiffuseCalculateDerivative(
              const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
              const Array<OneD, Array<OneD, NekDouble>>         &inarray,
              TensorOfArray3D<NekDouble>                        &qfields,
              const Array<OneD, Array<OneD, NekDouble>>         &pFwd,
              const Array<OneD, Array<OneD, NekDouble>>         &pBwd)
        {
            v_DiffuseCalculateDerivative(fields, inarray, qfields, pFwd, pBwd);
        }


        const Array<OneD, const Array<OneD, NekDouble> >
        &Diffusion::v_GetTraceNormal()
        {
            ASSERTL0(false,"v_GetTraceNormal not defined");
            return NullNekDoubleArrayofArray;
        }

        void Diffusion::v_ConsVarAveJump(
                                         const std::size_t                             nConvectiveFields,
                                         const size_t                                        npnts,
                                         const Array<OneD, const Array<OneD, NekDouble> >    &vFwd,
                                         const Array<OneD, const Array<OneD, NekDouble> >    &vBwd,
                                         Array<OneD,       Array<OneD, NekDouble> >    &aver,
                                         Array<OneD,       Array<OneD, NekDouble> >    &jump)
        {
            boost::ignore_unused(nConvectiveFields, npnts, vFwd, vBwd,
                                 aver, jump);
            ASSERTL0(false, "v_ConsVarAveJump not defined");
        }

        void Diffusion::v_DiffuseCalculateDerivative(
                                                     const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                                                     const Array<OneD, Array<OneD, NekDouble>>         &inarray,
                                                     TensorOfArray3D<NekDouble>                       &qfields,
                                                     const Array<OneD, Array<OneD, NekDouble>>         &pFwd,
                                                     const Array<OneD, Array<OneD, NekDouble>>         &pBwd)
        {
            boost::ignore_unused(fields, inarray, qfields,
                                 pFwd, pBwd);
            ASSERTL0(false, "Not defined for function DiffuseVolumeFLux.");
        }

        void Diffusion::v_DiffuseVolumeFlux(
                                            const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                                            const Array<OneD, Array<OneD, NekDouble>>           &inarray,
                                            TensorOfArray3D<NekDouble>                          &qfields,
                                            TensorOfArray3D<NekDouble>                          &VolumeFlux,
                                            Array< OneD, int >                                  &nonZeroIndex)
        {
            boost::ignore_unused(fields, inarray, qfields, VolumeFlux,
                                 nonZeroIndex);
            ASSERTL0(false, "Not defined for function DiffuseVolumeFLux.");
        }

        void Diffusion::v_DiffuseTraceFlux(
              const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
              const Array<OneD, Array<OneD, NekDouble>>           &inarray,
              TensorOfArray3D<NekDouble>                          &qfields,
              TensorOfArray3D<NekDouble>                          &VolumeFlux,
              Array<OneD, Array<OneD, NekDouble> >                &TraceFlux,
              const Array<OneD, Array<OneD, NekDouble>>           &pFwd,
              const Array<OneD, Array<OneD, NekDouble>>           &pBwd,
              Array< OneD, int >                                  &nonZeroIndex)
        {
            boost::ignore_unused(fields, inarray, qfields, VolumeFlux,
                                 TraceFlux, pFwd, pBwd, nonZeroIndex);
            ASSERTL0(false, "Not defined function DiffuseTraceFLux.");
        }

        /**
         * @brief Compute primary variables
         *
         */
        void Diffusion::v_GetPrimVar(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble>>         &inarray,
                  Array<OneD, Array<OneD, NekDouble>>         &primVar)
        {

            int nDim = fields[0]->GetCoordim(0);
            int nPts = fields[0]->GetTotPoints();
            for(int i = 0; i < nDim; ++i)
            {
                primVar[i] = Array<OneD, NekDouble>(nPts, 0.0);
                Vmath::Vdiv(nPts, inarray[i+1], 1, inarray[0], 1, primVar[i], 1);
            }
        }

        /**
         * @brief Compute and store scalars needed for shock sensor, i.e.
         * divergence and curl squared
         * @param dimensions, points, primary variables derivative
         *
         */
        void Diffusion::GetDivCurl(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &pVarDer)
        {
            int nDim = fields[0]->GetCoordim(0);
            int nPts = fields[0]->GetTotPoints();
            // Compute and store scalars needed for shock sensor
            Array<OneD, NekDouble>              tmp3(nPts, 0.0);
            Array<OneD, NekDouble>              tmp4(nPts, 0.0);

            // div vel
            for (int j = 0; j < nDim; ++j)
            {
                Vmath::Vadd(nPts, m_divVel, 1, pVarDer[j][j], 1, m_divVel, 1);
            }
            // (div vel)**2
            Vmath::Vmul(nPts, m_divVel, 1, m_divVel, 1, m_divVelSquare, 1);

            // curl vel
            if ( nDim > 2 )
            {
                // curl[0] 3/2
                Vmath::Vadd(nPts, tmp4, 1, pVarDer[2][1], 1, tmp4, 1);
                // curl[0]-2/3
                Vmath::Vcopy(nPts, pVarDer[1][2], 1, tmp3, 1);
                Vmath::Neg(nPts, tmp3, 1);
                Vmath::Vadd(nPts, tmp4, 1, tmp3, 1, tmp4, 1);
                // square curl[0]
                Vmath::Vmul(nPts, tmp4, 1, tmp4, 1, tmp3, 1);
                Vmath::Vadd(nPts, m_curlVelSquare, 1, tmp3, 1, m_curlVelSquare, 1);

                tmp4 = Array<OneD, NekDouble>(nPts, 0.0);
                // curl[1] 3/1
                Vmath::Vadd(nPts, tmp4, 1, pVarDer[2][0], 1, tmp4, 1);
                // curl[1]-1/3
                Vmath::Vcopy(nPts, pVarDer[0][2], 1, tmp3, 1);
                Vmath::Neg(nPts, tmp3, 1);
                Vmath::Vadd(nPts, tmp4, 1, tmp3, 1, tmp4, 1);
                // square curl[1]
                Vmath::Vmul(nPts, tmp4, 1, tmp4, 1, tmp3, 1);
                Vmath::Vadd(nPts, m_curlVelSquare, 1, tmp3, 1, m_curlVelSquare, 1);

                tmp4 = Array<OneD, NekDouble>(nPts, 0.0);
                // curl[2] 1/2
                Vmath::Vadd(nPts, tmp4, 1, pVarDer[0][1], 1, tmp4, 1);
                // curl[2]-2/1
                Vmath::Vcopy(nPts, pVarDer[1][0], 1, tmp3, 1);
                Vmath::Neg(nPts, tmp3, 1);
                Vmath::Vadd(nPts, tmp4, 1, tmp3, 1, tmp4, 1);
                // square curl[2]
                Vmath::Vmul(nPts, tmp4, 1, tmp4, 1, tmp3, 1);
                Vmath::Vadd(nPts, m_curlVelSquare, 1, tmp3, 1, m_curlVelSquare, 1);

            }
            else if ( nDim > 1 )
            {
                // curl[2] 1/2
                Vmath::Vadd(nPts, tmp4, 1, pVarDer[0][1], 1, tmp4, 1);
                // curl[2]-2/1
                Vmath::Vcopy(nPts, pVarDer[1][0], 1, tmp3, 1);
                Vmath::Neg(nPts, tmp3, 1);
                Vmath::Vadd(nPts, tmp4, 1, tmp3, 1, tmp4, 1);
                // square curl[2]
                Vmath::Vmul(nPts, tmp4, 1, tmp4, 1, tmp3, 1);
                Vmath::Vadd(nPts, m_curlVelSquare, 1, tmp3, 1, m_curlVelSquare, 1);
            }
            else
            {
                Vmath::Fill(nPts, 0.0, m_curlVelSquare, 1);
            }
        }
    }
}
