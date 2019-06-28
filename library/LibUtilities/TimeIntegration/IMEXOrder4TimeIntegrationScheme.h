#pragma once

///////////////////////////////////////////////////////////////////////////////
//
// File: IMEXOrder4TimeIntegrationScheme.h
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

    class IMEXOrder4TimeIntegrationScheme : public TimeIntegrationScheme
    {
    public:
  
      IMEXOrder4TimeIntegrationScheme() : TimeIntegrationScheme() 
      {
          std::cout << "IMEXOrder4TimeIntegrationScheme Construtor: this is " << this << "\n";

          m_integration_phases = TimeIntegrationSchemeDataVector( 4 );
          m_integration_phases[ 0 ] = TimeIntegrationSchemeDataSharedPtr( new TimeIntegrationSchemeData( this ) );
          m_integration_phases[ 1 ] = TimeIntegrationSchemeDataSharedPtr( new TimeIntegrationSchemeData( this ) );
          m_integration_phases[ 2 ] = TimeIntegrationSchemeDataSharedPtr( new TimeIntegrationSchemeData( this ) );
          m_integration_phases[ 3 ] = TimeIntegrationSchemeDataSharedPtr( new TimeIntegrationSchemeData( this ) );

          IMEXdirk_2_3_3TimeIntegrationScheme::SetupSchemeData( m_integration_phases[0] );
          IMEXdirk_2_3_3TimeIntegrationScheme::SetupSchemeData( m_integration_phases[1] );
          IMEXOrder3TimeIntegrationScheme::SetupSchemeData(     m_integration_phases[2] );
          IMEXOrder4TimeIntegrationScheme::SetupSchemeData(     m_integration_phases[3] );

          std::cout << "done with IMEXOrder4TimeIntegrationScheme constructor\n";
      }

      virtual ~IMEXOrder4TimeIntegrationScheme()
      {
      }

      /////////////

      static TimeIntegrationSchemeSharedPtr create()
      {
        std::cout << "IMEXOrder4TimeIntegrationScheme::create()\n";
        TimeIntegrationSchemeSharedPtr p = MemoryManager<IMEXOrder4TimeIntegrationScheme>::AllocateSharedPtr();
        return p;
      }

      static std::string className;

      //////////////

      LUE virtual
          TimeIntegrationMethod GetIntegrationMethod() const { return TimeIntegrationMethod::eIMEXOrder4; }

      //////////////

      LUE
      static
      void SetupSchemeData( TimeIntegrationSchemeDataSharedPtr & phase )
      {
          
        phase->m_method = TimeIntegrationMethod::eIMEXOrder4;
        phase->m_schemeType = eIMEX;

        phase->m_numsteps  = 8;
        phase->m_numstages = 1;

        phase->m_numMultiStepValues = 4;
        phase->m_numMultiStepDerivs = 4;

        phase->m_timeLevelOffset = Array<OneD,unsigned int>( phase->m_numsteps );

        phase->m_timeLevelOffset[0] = 0;
        phase->m_timeLevelOffset[1] = 1;
        phase->m_timeLevelOffset[2] = 2;
        phase->m_timeLevelOffset[3] = 3;
        phase->m_timeLevelOffset[4] = 0;
        phase->m_timeLevelOffset[5] = 1;
        phase->m_timeLevelOffset[6] = 2;
        phase->m_timeLevelOffset[7] = 3;

        phase->m_A = Array<OneD, Array<TwoD,NekDouble> >(2);
        phase->m_B = Array<OneD, Array<TwoD,NekDouble> >(2);

        NekDouble twentyfifth = 1.0/25.0;

        phase->m_A[0] = Array<TwoD,NekDouble>( phase->m_numstages, phase->m_numstages, 12*twentyfifth );
        phase->m_B[0] = Array<TwoD,NekDouble>( phase->m_numsteps,  phase->m_numstages, 0.0 );
        phase->m_A[1] = Array<TwoD,NekDouble>( phase->m_numstages, phase->m_numstages, 0.0 );
        phase->m_B[1] = Array<TwoD,NekDouble>( phase->m_numsteps,  phase->m_numstages, 0.0 );
        phase->m_U    = Array<TwoD,NekDouble>( phase->m_numstages, phase->m_numsteps,  48*twentyfifth );
        phase->m_V    = Array<TwoD,NekDouble>( phase->m_numsteps,  phase->m_numsteps,  0.0 );
                    
        phase->m_B[0][0][0] = 12*twentyfifth;
        phase->m_B[1][4][0] = 1.0;
        phase->m_U[0][1] = -36*twentyfifth;
        phase->m_U[0][2] =  16*twentyfifth;
        phase->m_U[0][3] = -3*twentyfifth;
        phase->m_U[0][5] = -72*twentyfifth;
        phase->m_U[0][7] = -12*twentyfifth;

        phase->m_V[0][0] =  48*twentyfifth;
        phase->m_V[0][1] = -36*twentyfifth;
        phase->m_V[0][2] =  16*twentyfifth;
        phase->m_V[0][3] = -3*twentyfifth;
        phase->m_V[0][4] =  48*twentyfifth;
        phase->m_V[0][5] = -72*twentyfifth;
        phase->m_V[0][6] =  48*twentyfifth;
        phase->m_V[0][7] = -12*twentyfifth;
        phase->m_V[1][0] =  1.0;
        phase->m_V[2][1] =  1.0;
        phase->m_V[4][3] =  1.0;
        phase->m_V[5][4] =  1.0;
        phase->m_V[6][5] =  1.0;
        phase->m_V[7][6] =  1.0;

        phase->m_firstStageEqualsOldSolution = phase->CheckIfFirstStageEqualsOldSolution( phase->m_A, phase->m_B, phase->m_U, phase->m_V );
        phase->m_lastStageEqualsNewSolution  = phase->CheckIfLastStageEqualsNewSolution(  phase->m_A, phase->m_B, phase->m_U, phase->m_V );

        ASSERTL1( phase->VerifyIntegrationSchemeType( phase->GetIntegrationSchemeType(), phase->m_A, phase->m_B, phase->m_U, phase->m_V ),
                  "Time integration scheme coefficients do not match its type" );
      }

    }; // end class IMEXOrder4TimeIntegrationScheme

  } // end namespace LibUtilities
} // end namespace Nektar
