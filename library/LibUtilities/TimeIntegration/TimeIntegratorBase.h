#pragma once

///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationWrapper.h
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

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

#include <string>

///////////////////////////////////////////////////////////////////////////////

#define LUE LIB_UTILITIES_EXPORT

///////////////////////////////////////////////////////////////////////////////

namespace Nektar {
  namespace LibUtilities {

    ///////////////////////////////////////////////////////////////////////////////////
    // Provide using code the ability to get the TimeIntegrator Factory so that it can
    // create TimeIntegrators.

    class TimeIntegratorBase; // Forward declaration

    /// Datatype of the NekFactory used to instantiate classes derived from the EquationSystem class.
    typedef NekFactory < std::string, TimeIntegratorBase > TimeIntegratorFactory;

    // Allows a code to create a TimeIntegrator. Usually used like this:
    //
    //    LibUtilities::TimeIntegratorSharedPtr timeIntegrator = LibUtilities::GetTimeIntegratorFactory().CreateInstance( "IMEXOrder1" );
    //
    LUE TimeIntegratorFactory & GetTimeIntegratorFactory();

    typedef std::shared_ptr<TimeIntegratorBase> TimeIntegratorSharedPtr;
  
    ///////////////////////////////////////////////////////////////////////////////////

    class TimeIntegratorBase
    {
    public:

      LUE
      TimeIntegratorBase() {}

      LUE
      virtual ~TimeIntegratorBase() {}

      LUE
      inline void InitObject() { v_InitObject(); }

      // Dd: This was originally named InitializeScheme()... however, I believe that is a misnomer, so have renamed it...  Need to verify.
      LUE
      TimeIntegrationSolutionSharedPtr InitializeIntegrator( const NekDouble                                 timestep,
                                                                   TimeIntegrationScheme::ConstDoubleArray & y_0,
                                                             const NekDouble                                 time,
                                                             const TimeIntegrationSchemeOperators          & op )
      {
        cout << "here: " << m_intScheme.size() << "\n";
        int steps = GetIntegrationSteps();
        return m_intScheme[ steps - 1 ]->InitializeScheme( timestep, y_0, time, op );
      }

      LUE
      TimeIntegrationScheme::ConstDoubleArray &
      TimeIntegrate( const int                                timestep,
                     const NekDouble                          delta_t,
                           TimeIntegrationSolutionSharedPtr & solvector,
                     const TimeIntegrationSchemeOperators   & op)
      {
        int steps = GetIntegrationSteps();
        return m_intScheme[ min( timestep, steps - 1 ) ]->TimeIntegrate( delta_t, solvector, op );
      }

      virtual TimeIntegrationMethod GetIntegrationMethod() const = 0;
              unsigned int          GetIntegrationSteps() const { return m_intScheme.size(); } //{ return m_intSteps; }

    protected:
      // TimeIntegrationMethod                       m_method;    // Dd: not sure if the integrator or the scheme should hold the method... or both?  Thought just the scheme but it is needed (for now)
      // int                                         m_intSteps;  // Dd: use size of m_intScheme instead of explicit?
      std::vector<TimeIntegrationSchemeSharedPtr> m_intScheme;

      virtual void v_InitObject() = 0;

    };

  } // end namespace LibUtilities
} // end namespace Nektar
