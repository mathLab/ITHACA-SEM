#pragma once

///////////////////////////////////////////////////////////////////////////////
//
// File: AdamsBashforthOrder3TimeIntegrationScheme.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2018 Division of Applied Mathematics, Brown University (USA),
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
// Description: Header file of time integration scheme wrappers
//
///////////////////////////////////////////////////////////////////////////////

#define LUE LIB_UTILITIES_EXPORT

///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

///////////////////////////////////////////////////////////////////////////////

namespace Nektar {
  namespace LibUtilities {

    class AdamsBashforthOrder3TimeIntegrationScheme : public TimeIntegrationScheme
    {
    public:
  
      AdamsBashforthOrder3TimeIntegrationScheme() : TimeIntegrationScheme() 
      {
          m_integration_phases = TimeIntegrationSchemeDataVector( 2 );
          m_integration_phases[ 0 ] = TimeIntegrationSchemeDataSharedPtr( new TimeIntegrationSchemeData( this ) );
          m_integration_phases[ 1 ] = TimeIntegrationSchemeDataSharedPtr( new TimeIntegrationSchemeData( this ) );

          ForwardEulerTimeIntegrationScheme::SetupSchemeData(         m_integration_phases[0] );
          AdamsBashforthOrder3TimeIntegrationScheme::SetupSchemeData( m_integration_phases[1] );
      }

      virtual ~AdamsBashforthOrder3TimeIntegrationScheme()
      {
      }

      /////////////

      static TimeIntegrationSchemeSharedPtr create()
      {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<AdamsBashforthOrder3TimeIntegrationScheme>::AllocateSharedPtr();
        return p;
      }

      static std::string className;

      //////////////

      LUE virtual
          TimeIntegrationMethod GetIntegrationMethod() const { return TimeIntegrationMethod::eAdamsBashforthOrder3; }

      //////////////

      LUE
      static
      void SetupSchemeData( TimeIntegrationSchemeDataSharedPtr & phase )
      {
        phase->m_method = TimeIntegrationMethod::eAdamsBashforthOrder3;
        phase->m_schemeType = eExplicit;

        phase->m_numsteps  = 4;
        phase->m_numstages = 1;

        phase->m_A = Array<OneD, Array<TwoD,NekDouble> >(1);
        phase->m_B = Array<OneD, Array<TwoD,NekDouble> >(1);

        phase->m_A[0] = Array<TwoD,NekDouble>( phase->m_numstages, phase->m_numstages, 0.0 );
        phase->m_B[0] = Array<TwoD,NekDouble>( phase->m_numsteps,  phase->m_numstages, 0.0 );
        phase->m_U    = Array<TwoD,NekDouble>( phase->m_numstages, phase->m_numsteps,  0.0 );
        phase->m_V    = Array<TwoD,NekDouble>( phase->m_numsteps,  phase->m_numsteps,  0.0 );
                    
        phase->m_B[0][1][0] = 1.0;

        phase->m_U[0][0] = 1.0;
        phase->m_U[0][1] = 23.0/12.0;
        phase->m_U[0][2] = -4.0/3.0;
        phase->m_U[0][3] = 5.0/12.0;

        phase->m_V[0][0] = 1.0;
        phase->m_V[0][1] = 23.0/12.0;
        phase->m_V[0][2] = -4.0/3.0;
        phase->m_V[0][3] = 5.0/12.0;
        phase->m_V[2][1] = 1.0;
        phase->m_V[3][2] = 1.0;

        phase->m_numMultiStepValues = 1;
        phase->m_numMultiStepDerivs = 3;
        phase->m_timeLevelOffset = Array<OneD,unsigned int>( phase->m_numsteps );
        phase->m_timeLevelOffset[0] = 0;
        phase->m_timeLevelOffset[1] = 1;
        phase->m_timeLevelOffset[2] = 2;
        phase->m_timeLevelOffset[3] = 3;

        phase->m_firstStageEqualsOldSolution = phase->CheckIfFirstStageEqualsOldSolution( phase->m_A, phase->m_B, phase->m_U, phase->m_V );
        phase->m_lastStageEqualsNewSolution  = phase->CheckIfLastStageEqualsNewSolution(  phase->m_A, phase->m_B, phase->m_U, phase->m_V );

        ASSERTL1( phase->VerifyIntegrationSchemeType( phase->m_schemeType, phase->m_A, phase->m_B, phase->m_U, phase->m_V ),
                  "Time integration scheme coefficients do not match its type" );
      }

    }; // end class AdamsBashforthOrder3TimeIntegrationScheme

  } // end namespace LibUtilities
} // end namespace Nektar
