#pragma once

///////////////////////////////////////////////////////////////////////////////
//
// File: BDFImplicitOrder1TimeIntegrationScheme.h
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

#include <LibUtilities/TimeIntegration/BackwardEulerTimeIntegrationScheme.h>

///////////////////////////////////////////////////////////////////////////////

namespace Nektar {
  namespace LibUtilities {

    class BDFImplicitOrder1TimeIntegrationScheme : public TimeIntegrationScheme
    {
    public:
  
      BDFImplicitOrder1TimeIntegrationScheme() : TimeIntegrationScheme() 
      {
          m_integration_phases = TimeIntegrationSchemeDataVector( 1 );
          m_integration_phases[ 0 ] = TimeIntegrationSchemeDataSharedPtr( new TimeIntegrationSchemeData( this ) );

          BDFImplicitOrder1TimeIntegrationScheme::SetupSchemeData( m_integration_phases[0] );
      }

      virtual ~BDFImplicitOrder1TimeIntegrationScheme()
      {
      }

      /////////////

      static TimeIntegrationSchemeSharedPtr create()
      {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<BDFImplicitOrder1TimeIntegrationScheme>::AllocateSharedPtr();
        return p;
      }

      static std::string className;

      //////////////

      LUE virtual
          TimeIntegrationMethod GetIntegrationMethod() const { return TimeIntegrationMethod::eBDFImplicitOrder1; }

      //////////////

      LUE
      static
      void SetupSchemeData( TimeIntegrationSchemeDataSharedPtr & phase )
      {
        // FIXME: Dd? Correct way to do this?
        BackwardEulerTimeIntegrationScheme::SetupSchemeData( phase );
      }

    }; // end class BDFImplicitOrder1TimeIntegrationScheme

  } // end namespace LibUtilities
} // end namespace Nektar
