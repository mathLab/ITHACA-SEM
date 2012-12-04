///////////////////////////////////////////////////////////////////////////////
//
// File: Advection.cpp
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
// Description: Abstract base class for advection.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Advection/Advection.h>

namespace Nektar
{
    namespace SolverUtils
    {
        AdvectionFactory& GetAdvectionFactory()
        {
            typedef Loki::SingletonHolder<AdvectionFactory,
            Loki::CreateUsingNew,
            Loki::NoDestroy > Type;
            return Type::Instance();
        }
        
        void Advection::InitObject(
            const LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr>       pFields)
        {
            v_InitObject(pSession, pFields);
        }
        
        void Advection::SetupMetrics(
            const LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr>       pFields)
        {
            v_SetupMetrics(pSession, pFields);
        }
        
        void Advection::SetupCFunctions(
            const LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr>       pFields)
        {
            v_SetupCFunctions(pSession, pFields);
        }
        
        void Advection::SetupInterpolationMatrices(
            const LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr>       pFields)
        {
            v_SetupInterpolationMatrices(pSession, pFields);
        }
        
        void Advection::Advect(
            const int nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &advVel,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {
            v_Advect(nConvectiveFields, fields, advVel, inarray, outarray);
        }
        
        void Advection::DivCFlux_1D(
                const int nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, const NekDouble> &fluxX1, 
                const Array<OneD, const NekDouble> &numericalFlux,
                      Array<OneD,       NekDouble> &divCFlux)
        {
            v_DivCFlux_1D(nConvectiveFields, fields, fluxX1, numericalFlux, divCFlux);
        }
        
        void Advection::DivCFlux_2D(
                const int nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, const NekDouble> &fluxX1, 
                const Array<OneD, const NekDouble> &fluxX2, 
                const Array<OneD, const NekDouble> &numericalFlux,
                      Array<OneD,       NekDouble> &divCFlux)
        {
            v_DivCFlux_2D(nConvectiveFields, fields, fluxX1, fluxX2, numericalFlux, divCFlux);
        }
        
        void Advection::DivCFlux_3D(
                const int nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, const NekDouble> &fluxX1, 
                const Array<OneD, const NekDouble> &fluxX2,
                const Array<OneD, const NekDouble> &fluxX3, 
                const Array<OneD, const NekDouble> &numericalFlux,
                      Array<OneD,       NekDouble> &divCFlux)
        {
            v_DivCFlux_3D(nConvectiveFields, fields, fluxX1, fluxX2, fluxX3, numericalFlux, divCFlux);
        }
    }
}
