#pragma once

///////////////////////////////////////////////////////////////////////////////
//
// File: DIRKOrder2TimeIntegrator.h
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

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

///////////////////////////////////////////////////////////////////////////////

#define LUE LIB_UTILITIES_EXPORT

///////////////////////////////////////////////////////////////////////////////

namespace Nektar {
  namespace LibUtilities {

    class DIRKOrder2TimeIntegrationScheme : public TimeIntegrationScheme
    {
    public:
  
      DIRKOrder2TimeIntegrationScheme() : TimeIntegrationScheme() 
      {
          m_integration_phases = TimeIntegrationSchemeDataVector( 1 );
          m_integration_phases[ 0 ] = TimeIntegrationSchemeDataSharedPtr( new TimeIntegrationSchemeData( this ) );

          DIRKOrder2TimeIntegrationScheme::SetupSchemeData( m_integration_phases[0] );
      }

      virtual ~DIRKOrder2TimeIntegrationScheme()
      {
      }

      /////////////

      static TimeIntegrationSchemeSharedPtr create()
      {
        std::cout << "DIRKOrder2TimeIntegrationScheme::create()\n";

        TimeIntegrationSchemeSharedPtr p = MemoryManager<DIRKOrder2TimeIntegrationScheme>::AllocateSharedPtr();
        return p;
      }

      static std::string className;

      //////////////

      LUE virtual
          TimeIntegrationMethod GetIntegrationMethod() const { return TimeIntegrationMethod::eIMEXdirk_3_4_3; }

      //////////////

      LUE
      static
      void SetupSchemeData( TimeIntegrationSchemeDataSharedPtr & phase )
      {
        phase->m_method = TimeIntegrationMethod::eIMEXdirk_3_4_3;
        phase->m_schemeType = eDiagonallyImplicit;

        phase->m_numsteps  = 1;
        phase->m_numstages = 1;

        phase->m_A = Array<OneD, Array<TwoD,NekDouble> >(1);
        phase->m_B = Array<OneD, Array<TwoD,NekDouble> >(1);

        phase->m_A[0] = Array<TwoD,NekDouble>( phase->m_numstages, phase->m_numstages, 0.0 );
        phase->m_B[0] = Array<TwoD,NekDouble>( phase->m_numsteps,  phase->m_numstages, 0.0 );
        phase->m_U    = Array<TwoD,NekDouble>( phase->m_numstages, phase->m_numsteps,  1.0 );
        phase->m_V    = Array<TwoD,NekDouble>( phase->m_numsteps,  phase->m_numsteps,  1.0 );
                    
        NekDouble lambda = (2.0-sqrt(2.0))/2.0;

        phase->m_A[0][0][0] = lambda;
        phase->m_A[0][1][0] = 1.0 - lambda;
        phase->m_A[0][1][1] = lambda;

        phase->m_B[0][0][0] = 1.0 - lambda;
        phase->m_B[0][0][1] = lambda;

        phase->m_numMultiStepValues = 1;
        phase->m_numMultiStepDerivs = 0;
        phase->m_timeLevelOffset = Array<OneD,unsigned int>( phase->m_numsteps );
        phase->m_timeLevelOffset[0] = 0;

        phase->m_firstStageEqualsOldSolution = phase->CheckIfFirstStageEqualsOldSolution( phase->m_A, phase->m_B, phase->m_U, phase->m_V );
        phase->m_lastStageEqualsNewSolution  = phase->CheckIfLastStageEqualsNewSolution(  phase->m_A, phase->m_B, phase->m_U, phase->m_V );

        ASSERTL1( phase->VerifyIntegrationSchemeType( phase->m_schemeType, phase->m_A, phase->m_B, phase->m_U, phase->m_V ),
                  "Time integration scheme coefficients do not match its type" );
      }
      
    }; // end class DIRKOrder2TimeIntegrator

  } // end namespace LibUtilities
} // end namespace Nektar
