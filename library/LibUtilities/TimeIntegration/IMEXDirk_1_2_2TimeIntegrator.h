#pragma once

///////////////////////////////////////////////////////////////////////////////
//
// File: IMEXDirk_1_2_2TimeIntegrator.h
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

#include <LibUtilities/TimeIntegration/TimeIntegratorBase.h>

///////////////////////////////////////////////////////////////////////////////

#define LUE LIB_UTILITIES_EXPORT

///////////////////////////////////////////////////////////////////////////////

namespace Nektar {
  namespace LibUtilities {

    class IMEXDirk_1_2_2TimeIntegrator : public TimeIntegratorBase
    {
    public:
  
      virtual ~IMEXDirk_1_2_2TimeIntegrator()
      {
      }

      /////////////

      static TimeIntegratorSharedPtr create()
      {
        TimeIntegratorSharedPtr p = MemoryManager<IMEXDirk_1_2_2TimeIntegrator>::AllocateSharedPtr();
        p->InitObject();
        return p;
      }

      static std::string className; // Will be set in TimeIntegratorBase.cpp during program start up.

      virtual void v_InitObject()
      {
        int steps   = 1;
        m_intScheme = vector<TimeIntegrationSchemeSharedPtr>( steps );

        TimeIntegrationSchemeKey IntKey0( eIMEXdirk_1_2_2 );
        m_intScheme[0] = GetTimeIntegrationSchemeManager()[ IntKey0 ];
      }

      //////////////

      // Replaces (from TimeIntegrationScheme.h): return m_schemeKey.GetIntegrationMethod();
      LUE TimeIntegrationMethod GetIntegrationMethod() const { return TimeIntegrationMethod::eIMEXdirk_1_2_2; }
      // LUE unsigned int          GetIntegrationSteps() const  { return m_intSteps; }

      //////////////

      LUE
      static
      void SetupScheme( TimeIntegrationScheme * scheme )
      {
        scheme->m_numsteps  = 1;
        scheme->m_numstages = 2;

        scheme->m_A = Array<OneD, Array<TwoD,NekDouble> >(2);
        scheme->m_B = Array<OneD, Array<TwoD,NekDouble> >(2);

        scheme->m_A[0] = Array<TwoD,NekDouble>( scheme->m_numstages, scheme->m_numstages, 0.0 );
        scheme->m_B[0] = Array<TwoD,NekDouble>( scheme->m_numsteps,  scheme->m_numstages, 0.0 );
        scheme->m_A[1] = Array<TwoD,NekDouble>( scheme->m_numstages, scheme->m_numstages, 0.0 );
        scheme->m_B[1] = Array<TwoD,NekDouble>( scheme->m_numsteps,  scheme->m_numstages, 0.0 );
        scheme->m_U    = Array<TwoD,NekDouble>( scheme->m_numstages, scheme->m_numsteps,  1.0 );
        scheme->m_V    = Array<TwoD,NekDouble>( scheme->m_numsteps,  scheme->m_numsteps,  1.0 );

        scheme->m_A[0][1][1] = 1.0/2.0;

        scheme->m_B[0][0][1] = 1.0;

        scheme->m_A[1][1][0] = 1.0/2.0;

        scheme->m_B[1][0][1] = 1.0;

        scheme->m_schemeType = eIMEX;

        scheme->m_numMultiStepValues = 1;
        scheme->m_numMultiStepDerivs = 0;
        scheme->m_timeLevelOffset = Array<OneD,unsigned int>( scheme->m_numsteps );
        scheme->m_timeLevelOffset[0] = 0;

        scheme->m_firstStageEqualsOldSolution = scheme->CheckIfFirstStageEqualsOldSolution( scheme->m_A, scheme->m_B, scheme->m_U, scheme->m_V );
        scheme->m_lastStageEqualsNewSolution  = scheme->CheckIfLastStageEqualsNewSolution(  scheme->m_A, scheme->m_B, scheme->m_U, scheme->m_V );

        ASSERTL1( scheme->VerifyIntegrationSchemeType( scheme->m_schemeType, scheme->m_A, scheme->m_B, scheme->m_U, scheme->m_V ),
                  "Time integration scheme coefficients do not match its type" );
      }
      
    }; // end class IMEXDirk_1_2_2TimeIntegrator

  } // end namespace LibUtilities
} // end namespace Nektar
