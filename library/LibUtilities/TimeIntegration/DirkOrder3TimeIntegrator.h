#pragma once

///////////////////////////////////////////////////////////////////////////////
//
// File: DirkOrder3TimeIntegrator.h
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

    class DirkOrder3TimeIntegrator : public TimeIntegratorBase
    {
    public:
  
      virtual ~DirkOrder3TimeIntegrator()
      {
      }

      /////////////

      static TimeIntegratorSharedPtr create()
      {
        TimeIntegratorSharedPtr p = MemoryManager<DirkOrder3TimeIntegrator>::AllocateSharedPtr();
        p->InitObject();
        return p;
      }

      static std::string className; // Will be set in TimeIntegratorBase.cpp during program start up.

      virtual void v_InitObject()
      {
        TimeIntegrationSchemeKey IntKey0( eDIRKOrder3 );
        int steps   = 1;
        m_intScheme = vector<TimeIntegrationSchemeSharedPtr>( steps );

        m_intScheme[0] = GetTimeIntegrationSchemeManager()[ IntKey0 ]; // Dd: not sure why this is used or if it is needed in this new version...
      }

      //////////////

      // Replaces (from TimeIntegrationScheme.h): return m_schemeKey.GetIntegrationMethod();
      LUE virtual
          TimeIntegrationMethod GetIntegrationMethod() const     { return TimeIntegrationMethod::eIMEXdirk_3_4_3; }
      // LUE unsigned int          GetIntegrationSteps() const      { return m_intSteps; }

      //////////////

      LUE
      static
      void SetupScheme( TimeIntegrationScheme * scheme )
      {
        scheme->m_numsteps  = 1;
        scheme->m_numstages = 3;

        scheme->m_A = Array<OneD, Array<TwoD,NekDouble> >(1);
        scheme->m_B = Array<OneD, Array<TwoD,NekDouble> >(1);

        scheme->m_A[0] = Array<TwoD,NekDouble>( scheme->m_numstages, scheme->m_numstages, 0.0 );
        scheme->m_B[0] = Array<TwoD,NekDouble>( scheme->m_numsteps,  scheme->m_numstages, 0.0 );
        scheme->m_U    = Array<TwoD,NekDouble>( scheme->m_numstages, scheme->m_numsteps,  1.0 );
        scheme->m_V    = Array<TwoD,NekDouble>( scheme->m_numsteps,  scheme->m_numsteps,  1.0 );
                    
        NekDouble lambda = 0.4358665215;

        scheme->m_A[0][0][0] = lambda;
        scheme->m_A[0][1][0] = 0.5  * (1.0 - lambda);
        scheme->m_A[0][2][0] = 0.25 * (-6.0*lambda*lambda + 16.0*lambda - 1.0);
        scheme->m_A[0][1][1] = lambda;
        scheme->m_A[0][2][1] = 0.25 * ( 6.0*lambda*lambda - 20.0*lambda + 5.0);
        scheme->m_A[0][2][2] = lambda;

        scheme->m_B[0][0][0] = 0.25 * (-6.0*lambda*lambda + 16.0*lambda - 1.0);
        scheme->m_B[0][0][1] = 0.25 * ( 6.0*lambda*lambda - 20.0*lambda + 5.0);
        scheme->m_B[0][0][2] = lambda;

        scheme->m_schemeType = eDiagonallyImplicit;

        scheme->m_numMultiStepValues = 1;
        scheme->m_numMultiStepDerivs = 0;
        scheme->m_timeLevelOffset = Array<OneD,unsigned int>( scheme->m_numsteps );
        scheme->m_timeLevelOffset[0] = 0;

        scheme->m_firstStageEqualsOldSolution = scheme->CheckIfFirstStageEqualsOldSolution( scheme->m_A, scheme->m_B, scheme->m_U, scheme->m_V );
        scheme->m_lastStageEqualsNewSolution  = scheme->CheckIfLastStageEqualsNewSolution(  scheme->m_A, scheme->m_B, scheme->m_U, scheme->m_V );

        ASSERTL1( scheme->VerifyIntegrationSchemeType( scheme->m_schemeType, scheme->m_A, scheme->m_B, scheme->m_U, scheme->m_V ),
                  "Time integration scheme coefficients do not match its type" );
      }
      
    }; // end class DirkOrder3TimeIntegrator

  } // end namespace LibUtilities
} // end namespace Nektar
