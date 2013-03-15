///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLFR.h
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
// Description: LFR diffusion class.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DIFFUSIONLFR
#define NEKTAR_SOLVERUTILS_DIFFUSIONLFR

#include <SolverUtils/Diffusion/Diffusion.h>

namespace Nektar
{
    namespace SolverUtils
    {
        class DiffusionLFR : public Diffusion
        {
        public:
            static DiffusionSharedPtr create(std::string diffType)
            {
                return DiffusionSharedPtr(new DiffusionLFR(diffType));
            }
            
            static std::string                   type[];
            
            Array<OneD, Array<OneD, NekDouble> > m_Q2D_e0; 
            Array<OneD, Array<OneD, NekDouble> > m_Q2D_e1; 
            Array<OneD, Array<OneD, NekDouble> > m_Q2D_e2; 
            Array<OneD, Array<OneD, NekDouble> > m_Q2D_e3; 
            
            Array<OneD, Array<OneD, NekDouble> > m_dGL_xi1;                  
            Array<OneD, Array<OneD, NekDouble> > m_dGR_xi1;
            Array<OneD, Array<OneD, NekDouble> > m_dGL_xi2;                  
            Array<OneD, Array<OneD, NekDouble> > m_dGR_xi2;
            Array<OneD, Array<OneD, NekDouble> > m_dGL_xi3;                  
            Array<OneD, Array<OneD, NekDouble> > m_dGR_xi3;
            DNekMatSharedPtr                     m_Ixm;
            DNekMatSharedPtr                     m_Ixp;
            
        protected:
            DiffusionLFR(std::string diffType);
            
            Array<OneD, Array<OneD, NekDouble> >              m_traceNormals;
            Array<OneD, Array<OneD, Array<OneD,NekDouble> > > m_tanbasis;
            LibUtilities::SessionReaderSharedPtr              m_session;
            
            std::string m_diffType;
            
            virtual void v_InitObject(
                LibUtilities::SessionReaderSharedPtr               pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>        pFields);
            
            virtual void v_SetupMetrics(
                LibUtilities::SessionReaderSharedPtr               pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>        pFields);
            
            virtual void v_SetupCFunctions(
                LibUtilities::SessionReaderSharedPtr               pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>        pFields);
            
            virtual void v_SetupInterpolationMatrices(
                LibUtilities::SessionReaderSharedPtr               pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>        pFields);
            
            virtual void v_Diffuse(
                const int                                          nConvective,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray);
                        
            virtual void v_NumFluxforScalar(
                const Array<OneD, MultiRegions::ExpListSharedPtr>       &fields,
                const Array<OneD, Array<OneD, NekDouble> >              &ufield,
                      Array<OneD, Array<OneD, Array<OneD, NekDouble> > >&uflux);
            
            virtual void v_WeakPenaltyforScalar(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const int                                          var,
                const Array<OneD, const NekDouble>                &ufield,
                      Array<OneD,       NekDouble>                &penaltyflux);
            
            virtual void v_NumFluxforVector(
                const Array<OneD, MultiRegions::ExpListSharedPtr>       &fields,
                const Array<OneD, Array<OneD, NekDouble> >              &ufield,
                      Array<OneD, Array<OneD, Array<OneD, NekDouble> > >&qfield,
                      Array<OneD, Array<OneD, NekDouble> >              &qflux);
            
            virtual void v_WeakPenaltyforVector(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const int                                          var,
                const int                                          dir,
                const Array<OneD, const NekDouble>                &qfield,
                      Array<OneD,       NekDouble>                &penaltyflux,
                NekDouble                                          C11);
            
            virtual void v_DerCFlux_1D(
                const int                                     nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>&fields,
                const Array<OneD, const NekDouble>               &flux, 
                const Array<OneD, const NekDouble>               &iFlux,
                      Array<OneD,       NekDouble>               &derCFlux);
            
            virtual void v_DerCFlux_2D(
                const int                                     nConvectiveFields,
                const int                                        direction,
                const Array<OneD, MultiRegions::ExpListSharedPtr>&fields,
                const Array<OneD, const NekDouble>               &flux, 
                const Array<OneD,       NekDouble>               &iFlux,
                      Array<OneD,       NekDouble>               &derCFlux);
            
            virtual void v_DivCFlux_2D(
                const int                                      nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>&fields,
                const Array<OneD, const NekDouble>               &fluxX1, 
                const Array<OneD, const NekDouble>               &fluxX2, 
                const Array<OneD, const NekDouble>               &numericalFlux,
                      Array<OneD,       NekDouble>               &divCFlux);
        }; 
    }
}
    
#endif
