///////////////////////////////////////////////////////////////////////////////
//
// File: AdvectionFR.h
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
// Description: FR advection class.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_ADVECTIONFR
#define NEKTAR_SOLVERUTILS_ADVECTIONFR

#include <SolverUtils/Advection.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Communication/Comm.h>
#include <SpatialDomains/SpatialData.h>


namespace Nektar
{
    namespace SolverUtils
    {
        class AdvectionFR : public Advection
        {
        public:
            static AdvectionSharedPtr create()
            {
                return AdvectionSharedPtr(new AdvectionFR());
            }
            
            static std::string                              type;
            
            std::string                                     m_FRSchemeRecovered;    ///< FR scheme recovered (DG, SD, HU)
            Array<OneD, NekDouble>                          m_dGL;                  ///< Left derivatives of the correction functions
            Array<OneD, NekDouble>                          m_dGR;                  ///< Right derivatives of the correction functions
            
        protected:
            AdvectionFR();
            
            Array<OneD, Array<OneD, NekDouble> > m_traceNormals;
            
            virtual void v_InitObject(
                LibUtilities::SessionReaderSharedPtr              pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>       pFields);
            
            virtual void v_Advect(
                const int                                          nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &advVel,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray);
        }; 
    }
}
    
#endif
