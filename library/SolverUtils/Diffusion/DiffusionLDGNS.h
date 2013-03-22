///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLDGNS.h
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
// Description: LDG diffusion class.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DIFFUSIONLDGNS
#define NEKTAR_SOLVERUTILS_DIFFUSIONLDGNS

#include <SolverUtils/Diffusion/Diffusion.h>

namespace Nektar
{
    namespace SolverUtils
    {
        class DiffusionLDGNS : public Diffusion
        {
        public:
            static DiffusionSharedPtr create(std::string diffType)
            {
                return DiffusionSharedPtr(new DiffusionLDGNS());
            }
            
            static std::string type;
            
        protected:
            DiffusionLDGNS();
            
            Array<OneD, Array<OneD, NekDouble> > m_traceNormals;
            LibUtilities::SessionReaderSharedPtr m_session;
            NekDouble                            m_gamma;
            NekDouble                            m_gasConstant;
            NekDouble                            m_Twall;
            std::string                          m_ViscosityType;
            NekDouble                            m_mu;
            NekDouble                            m_thermalConductivity;
            NekDouble                            m_rhoInf;
            NekDouble                            m_pInf;

            
            virtual void v_InitObject(
                LibUtilities::SessionReaderSharedPtr               pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>        pFields);
            
            virtual void v_Diffuse(
                const int                                          nConvective,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray);
            
            virtual void v_NumericalFluxO1(
                const Array<OneD, MultiRegions::ExpListSharedPtr>      &fields,
                const Array<OneD, Array<OneD, NekDouble> >             &inarray,
                      Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                                            &numericalFluxO1);
            
            virtual void v_WeakPenaltyO1(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >      &penaltyfluxO1);
            
            virtual void v_NumericalFluxO2(
                const Array<OneD, MultiRegions::ExpListSharedPtr>       &fields,
                const Array<OneD, Array<OneD, NekDouble> >              &ufield,
                      Array<OneD, Array<OneD, Array<OneD, NekDouble> > >&qfield,
                      Array<OneD, Array<OneD, NekDouble> >              &qflux);
            
            virtual void v_WeakPenaltyO2(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const int                                          var,
                const int                                          dir,
                const Array<OneD, const NekDouble>                &qfield,
                      Array<OneD,       NekDouble>                &penaltyflux);
        }; 
    }
}

#endif

