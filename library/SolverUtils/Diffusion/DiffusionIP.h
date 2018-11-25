///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionIP.h
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

#ifndef NEKTAR_SOLVERUTILS_DIFFUSIONWEAKDG
#define NEKTAR_SOLVERUTILS_DIFFUSIONWEAKDG

#include <SolverUtils/Diffusion/Diffusion.h>

namespace Nektar
{
    namespace SolverUtils
    {
        class DiffusionIP : public Diffusion
        {
        public:
            static DiffusionSharedPtr create(std::string diffType)
            {
                return DiffusionSharedPtr(new DiffusionIP());
            }
            
            static std::string type;

            void ConsVarAve(
                const int                                           nConvectiveFields,
                const int                                           npnts,
                const Array<OneD, const Array<OneD, NekDouble> >    &vFwd,
                const Array<OneD, const Array<OneD, NekDouble> >    &vBwd,
                      Array<OneD,       Array<OneD, NekDouble> >    &aver); 
            void CalTraceSymFlux(
                const int                                                           nConvectiveFields,
                const int                                                           nDim,
                const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
                const Array<OneD, Array<OneD, NekDouble> >                          &solution_Aver,
                      Array<OneD, Array<OneD, NekDouble> >                          &solution_jump,
                      Array<OneD, int >                                             &nonZeroIndexsymm,
                      Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &traceSymflux);

        protected:
            DiffusionIP();
   		
            std::string                             m_shockCaptureType;
            NekDouble                               m_IPDebugParameter;
	         
            Array<OneD, Array<OneD, NekDouble> >              m_traceNormals;
            Array<OneD, NekDouble>                            m_tracBwdWeight;
            Array<OneD, NekDouble>                            m_traceNormDirctnElmtLength;
            LibUtilities::SessionReaderSharedPtr              m_session;

            void GetPenaltyFactor(
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                Array<OneD, NekDouble >                             factor); 

            void AddSymmFluxIntegral(
                const int                                                           nConvectiveFields,
                const int                                                           nDim,
                const int                                                           nPts,
                const int                                                           nTracePts,
                const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
                const Array<OneD, const int >                                       &nonZeroIndex,
                      Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &tracflux,
                      Array<OneD, Array<OneD, NekDouble> >                          &outarray);

            void Add2ndDeriv2Trace(
                const int                                                           nConvectiveFields,
                const int                                                           nDim,
                const int                                                           nPts,
                const int                                                           nTracePts,
                const NekDouble                                                     PenaltyFactor2,
                const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
                const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &qfield,
                      Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &numDeriv);

            virtual void v_InitObject(
                LibUtilities::SessionReaderSharedPtr               pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>        pFields);
            
            virtual void v_Diffuse(
                const int                                          nConvective,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayofArray);

            virtual void v_Diffuse_coeff(
                const int                                          nConvective,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayofArray);

            virtual void v_Diffuse_coeffOld(
                const int                                          nConvective,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayofArray);
            virtual void v_CalTraceNumFlux(
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
                      Array<OneD, Array<OneD, NekDouble> >                          &solution_jump);
            
            virtual void v_physFieldDeriv(
                const int                                                   nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>           &fields,
                const Array<OneD, Array<OneD, NekDouble> >                  &inarray,
                const Array<OneD, Array<OneD, NekDouble> >                  &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >                  &pBwd,
                      Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    &qfield);
            virtual void v_AddVolumDerivJac2Mat( 
                const int                                               nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>       &pFields,
                const Array<OneD, const Array<OneD, DNekMatSharedPtr> > &ElmtJac,
                const int                                               nfluxDir, 
                const int                                               nDervDir, 
                      Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >    &gmtxarray);

            virtual void v_ConsVarAveJump(
                const int                                           nConvectiveFields,
                const int                                           npnts,
                const Array<OneD, const Array<OneD, NekDouble> >    &vFwd,
                const Array<OneD, const Array<OneD, NekDouble> >    &vBwd,
                      Array<OneD,       Array<OneD, NekDouble> >    &aver,
                      Array<OneD,       Array<OneD, NekDouble> >    &jump);

            virtual const Array<OneD, const Array<OneD, NekDouble> > &v_GetTraceNormal()
            {
                return m_traceNormals;
            }
            
        }; 
    }
}
    
#endif
