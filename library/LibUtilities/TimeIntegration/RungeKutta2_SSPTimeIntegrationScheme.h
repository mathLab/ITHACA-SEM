#pragma once

///////////////////////////////////////////////////////////////////////////////
//
// File: RungeKutta2_SSPTimeIntegrationScheme.h
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

#include <LibUtilities/TimeIntegration/RungeKutta2_ImprovedEulerTimeIntegrationScheme.h>

///////////////////////////////////////////////////////////////////////////////

#define LUE LIB_UTILITIES_EXPORT

///////////////////////////////////////////////////////////////////////////////

namespace Nektar {
  namespace LibUtilities {

    class RungeKutta2_SSPTimeIntegrationScheme : public TimeIntegrationScheme
    {
    public:
  
      RungeKutta2_SSPTimeIntegrationScheme() : TimeIntegrationScheme() 
      {
          m_integration_phases = TimeIntegrationSchemeDataVector( 1 );
          m_integration_phases[ 0 ] = TimeIntegrationSchemeDataSharedPtr( new TimeIntegrationSchemeData( this ) );

          RungeKutta2_SSPTimeIntegrationScheme::SetupSchemeData( m_integration_phases[0] );
      }

      virtual ~RungeKutta2_SSPTimeIntegrationScheme()
      {
      }

      /////////////

      static TimeIntegrationSchemeSharedPtr create()
      {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<RungeKutta2_SSPTimeIntegrationScheme>::AllocateSharedPtr();
        return p;
      }

      static std::string className;

      //////////////

      LUE virtual
          TimeIntegrationMethod GetIntegrationMethod() const { return TimeIntegrationMethod::eRungeKutta2_SSP; }

      //////////////

      LUE
      static
      void SetupSchemeData( TimeIntegrationSchemeDataSharedPtr & phase )
      {
        // FIXME: Dd: Do it like this?  Note, this is never (currently) called anyway because of the switch statement in TimeIntegrationScheme::TimeIntegrationScheme()...
        // FIXME: are these the same things?
        RungeKutta2_ImprovedEulerTimeIntegrationScheme::SetupSchemeData( phase );
      }

    }; // end class RungeKutta2_SSPTimeIntegrator

  } // end namespace LibUtilities
} // end namespace Nektar
