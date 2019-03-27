#pragma once

///////////////////////////////////////////////////////////////////////////////
//
// File: IMEXOrder3TimeIntegrator.h
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

#include <LibUtilities/TimeIntegration/TimeIntegratorBase.h>

///////////////////////////////////////////////////////////////////////////////

namespace Nektar {
  namespace LibUtilities {

    class IMEXOrder3TimeIntegrator : public TimeIntegratorBase
    {
    public:
  
      virtual ~IMEXOrder3TimeIntegrator()
      {
      }

      /////////////

      static TimeIntegratorSharedPtr create()
      {
        TimeIntegratorSharedPtr p = MemoryManager<IMEXOrder3TimeIntegrator>::AllocateSharedPtr();
        p->InitObject();
        return p;
      }

      static std::string className; // Will be set to "IMEXOrder3" in TimeIntegratorBase.cpp during program start up.

      virtual void v_InitObject()
      {
        int steps   = 3;
        m_intScheme = vector<TimeIntegrationSchemeSharedPtr>( steps );

        TimeIntegrationSchemeKey IntKey0( eIMEXdirk_3_4_3 );
        TimeIntegrationSchemeKey IntKey1( eIMEXdirk_3_4_3 );
        TimeIntegrationSchemeKey IntKey2( eIMEXOrder3 );
        m_intScheme[0] = GetTimeIntegrationSchemeManager()[ IntKey0 ];
        m_intScheme[1] = GetTimeIntegrationSchemeManager()[ IntKey1 ];
        m_intScheme[2] = GetTimeIntegrationSchemeManager()[ IntKey2 ];
      }

      //////////////


      // Replaces (from TimeIntegrationScheme.h): return m_schemeKey.GetIntegrationMethod();
      LUE virtual
          TimeIntegrationMethod GetIntegrationMethod() const { return TimeIntegrationMethod::eIMEXOrder3; }
      // LUE unsigned int          GetIntegrationSteps() const  { return m_intSteps; }

      //////////////

      LUE
      static
      void SetupScheme( TimeIntegrationScheme * scheme )
      {
        cout << "SetupScheme for eIMEXOrder3\n";

        NekDouble eleventh = 1.0/11.0;

        scheme->m_numsteps  = 6;
        scheme->m_numstages = 1;

        scheme->m_A = Array<OneD, Array<TwoD,NekDouble> >(2);
        scheme->m_B = Array<OneD, Array<TwoD,NekDouble> >(2);

        scheme->m_A[0] = Array<TwoD,NekDouble>( scheme->m_numstages, scheme->m_numstages, 6*eleventh );
        scheme->m_B[0] = Array<TwoD,NekDouble>( scheme->m_numsteps,  scheme->m_numstages, 0.0 );
        scheme->m_A[1] = Array<TwoD,NekDouble>( scheme->m_numstages, scheme->m_numstages, 0.0 );
        scheme->m_B[1] = Array<TwoD,NekDouble>( scheme->m_numsteps,  scheme->m_numstages, 0.0 );
        scheme->m_U    = Array<TwoD,NekDouble>( scheme->m_numstages, scheme->m_numsteps,  18*eleventh );
        scheme->m_V    = Array<TwoD,NekDouble>( scheme->m_numsteps,  scheme->m_numsteps,  0.0 );
                    
        scheme->m_B[0][0][0] = 6*eleventh;
        scheme->m_B[1][3][0] = 1.0;
        scheme->m_U[0][1] = -9*eleventh;
        scheme->m_U[0][2] =  2*eleventh;
        scheme->m_U[0][4] = -18*eleventh;
        scheme->m_U[0][5] =  6*eleventh;

        scheme->m_V[0][0] =  18*eleventh;
        scheme->m_V[0][1] = -9*eleventh;
        scheme->m_V[0][2] =  2*eleventh;
        scheme->m_V[0][3] =  18*eleventh;
        scheme->m_V[0][4] = -18*eleventh;
        scheme->m_V[0][5] =  6*eleventh;
        scheme->m_V[1][0] =  1.0;
        scheme->m_V[2][1] =  1.0;
        scheme->m_V[4][3] =  1.0;
        scheme->m_V[5][4] =  1.0;

        scheme->m_schemeType = eIMEX;
        scheme->m_numMultiStepValues = 3;
        scheme->m_numMultiStepDerivs = 3;
        scheme->m_timeLevelOffset = Array<OneD,unsigned int>( scheme->m_numsteps );
        scheme->m_timeLevelOffset[0] = 0;
        scheme->m_timeLevelOffset[1] = 1;
        scheme->m_timeLevelOffset[2] = 2;
        scheme->m_timeLevelOffset[3] = 0;
        scheme->m_timeLevelOffset[4] = 1;
        scheme->m_timeLevelOffset[5] = 2;

        scheme->m_firstStageEqualsOldSolution = scheme->CheckIfFirstStageEqualsOldSolution( scheme->m_A, scheme->m_B, scheme->m_U, scheme->m_V );
        scheme->m_lastStageEqualsNewSolution  = scheme->CheckIfLastStageEqualsNewSolution(  scheme->m_A, scheme->m_B, scheme->m_U, scheme->m_V );

        ASSERTL1( scheme->VerifyIntegrationSchemeType( scheme->m_schemeType, scheme->m_A, scheme->m_B, scheme->m_U, scheme->m_V ),
                  "Time integration scheme coefficients do not match its type" );
      }

    }; // end class IMEXOrder3TimeIntegrator

  } // end namespace LibUtilities
} // end namespace Nektar
