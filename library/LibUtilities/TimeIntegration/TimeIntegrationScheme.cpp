///////////////////////////////////////////////////////////////////////////////
//
// File TimeIntegrationScheme.cpp
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
// Description: implementation of time integration key class
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/TimeIntegration/TimeIntegrationSolution.h>

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeData.h>

#include <LibUtilities/TimeIntegration/AdamsBashforthOrder2TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/AdamsBashforthOrder3TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/AdamsMoultonOrder1TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/AdamsMoultonOrder2TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/BDFImplicitOrder1TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/BDFImplicitOrder2TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/BackwardEulerTimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/ClassicalRungeKutta4TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/CNABTimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/DIRKOrder2TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/DIRKOrder3TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/ForwardEulerTimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/IMEXdirk_1_2_1TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/IMEXdirk_1_2_2TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/IMEXdirk_2_2_2TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/IMEXdirk_2_3_2TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/IMEXdirk_2_3_3TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/IMEXdirk_3_4_3TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/IMEXdirk_4_4_3TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/IMEXGearTimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/IMEXOrder1TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/IMEXOrder2TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/IMEXOrder3TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/MCNABTimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/RungeKutta2TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/RungeKutta2_ImprovedEulerTimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/RungeKutta2_SSPTimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/RungeKutta3_SSPTimeIntegrationScheme.h>

#include <LibUtilities/BasicConst/NektarUnivConsts.hpp>

#include <algorithm>
#include <iostream>

#include <math.h>

using namespace std;

namespace Nektar
{
    namespace LibUtilities 
    {  
        // !!! Always make sure that this matches TimeIntegrationMethod enum... !!!
        //
        const char* const TimeIntegrationScheme::TimeIntegrationMethodMap[ 33 ] =
          {
            "NoTimeIntegrationMethod",
            "AdamsBashforthOrder1",
            "AdamsBashforthOrder2",
            "AdamsBashforthOrder3",
            "AdamsMoultonOrder1",
            "AdamsMoultonOrder2",
            "BDFImplicitOrder1",
            "BDFImplicitOrder2",
            "ClassicalRungeKutta4",
            "RungeKutta4",
            "RungeKutta3_SSP",
            "RungeKutta2_ImprovedEuler",
            "RungeKutta2_SSP",
            "ForwardEuler",
            "BackwardEuler",
            "IMEXOrder1",
            "IMEXOrder2",
            "IMEXOrder3",
            "Midpoint",
            "RungeKutta2",
            "DIRKOrder2",
            "DIRKOrder3",
            "CNAB",
            "IMEXGear",
            "MCNAB",
            "IMEXdirk_1_1_1",
            "IMEXdirk_1_2_1",
            "IMEXdirk_1_2_2",
            "IMEXdirk_2_2_2",
            "IMEXdirk_2_3_2",
            "IMEXdirk_2_3_3",
            "IMEXdirk_3_4_3",
            "IMEXdirk_4_4_3",
          };

      TimeIntegrationMethod TimeIntegrationScheme::methodFromName( const string & name )
      {
        if(      name == "AdamsBashforthOrder1" ) {      return eAdamsBashforthOrder1; }
        else if( name == "AdamsBashforthOrder2" ) {      return eAdamsBashforthOrder2; }
        else if( name == "AdamsBashforthOrder3" ) {      return eAdamsBashforthOrder3; }
        else if( name == "AdamsMoultonOrder1" ) {        return eAdamsMoultonOrder1; }
        else if( name == "AdamsMoultonOrder2" ) {        return eAdamsMoultonOrder2; }
        else if( name == "BDFImplicitOrder1" ) {         return eBDFImplicitOrder1; }
        else if( name == "BDFImplicitOrder2" ) {         return eBDFImplicitOrder2; }
        else if( name == "ClassicalRungeKutta4" ) {      return eClassicalRungeKutta4; }
        else if( name == "RungeKutta4" ) {               return eRungeKutta4; }
        else if( name == "RungeKutta3_SSP" ) {           return eRungeKutta3_SSP; }
        else if( name == "RungeKutta2_ImprovedEuler" ) { return eRungeKutta2_ImprovedEuler; }
        else if( name == "RungeKutta2_SSP" ) {           return eRungeKutta2_SSP; }
        else if( name == "ForwardEuler" ) {              return eForwardEuler; }
        else if( name == "BackwardEuler" ) {             return eBackwardEuler; }
        else if( name == "IMEXOrder1" ) {                return eIMEXOrder1; }
        else if( name == "IMEXOrder2" ) {                return eIMEXOrder2; }
        else if( name == "IMEXOrder3" ) {                return eIMEXOrder3; }
        else if( name == "Midpoint" ) {                  return eMidpoint; }
        else if( name == "RungeKutta2" ) {               return eRungeKutta2; }
        else if( name == "DIRKOrder2" ) {                return eDIRKOrder2; }
        else if( name == "DIRKOrder3" ) {                return eDIRKOrder3; }
        else if( name == "CNAB" ) {                      return eCNAB; }
        else if( name == "IMEXGear" ) {                  return eIMEXGear; }
        else if( name == "MCNAB" ) {                     return eMCNAB; }
        else if( name == "IMEXdirk_1_1_1" ) {            return eIMEXdirk_1_1_1; }
        else if( name == "IMEXdirk_1_2_1" ) {            return eIMEXdirk_1_2_1; }
        else if( name == "IMEXdirk_1_2_2" ) {            return eIMEXdirk_1_2_2; }
        else if( name == "IMEXdirk_2_2_2" ) {            return eIMEXdirk_2_2_2; }
        else if( name == "IMEXdirk_2_3_2" ) {            return eIMEXdirk_2_3_2; }
        else if( name == "IMEXdirk_2_3_3" ) {            return eIMEXdirk_2_3_3; }
        else if( name == "IMEXdirk_3_4_3" ) {            return eIMEXdirk_3_4_3; }
        else if( name == "IMEXdirk_4_4_3" ) {            return eIMEXdirk_4_4_3; }
        else {
          string msg = "'" + name + "' is not a known TimeIntegrationMethod. (Check spelling?)";
          NEKERROR( ErrorUtil::efatal, msg );
          // return eNoTimeIntegrationMethod;
        }
        return eNoTimeIntegrationMethod; // Note: This line should never be reached, but it does remove a compiler warning.
      }

      string TimeIntegrationScheme::nameFromMethod( const TimeIntegrationMethod method )
      {
        return TimeIntegrationMethodMap[ method ];
      }

      TimeIntegrationSchemeFactory & GetTimeIntegrationSchemeFactory()
      {
        static TimeIntegrationSchemeFactory instance;
        return instance;
      }

      // FIME: Dd: A "Manager" gives you back a "shared" instance... In theory we want separate versions of each of the data items...
      //           So the manager should be removed (perhaps replaced with a Factory, but not sure...)
        // TimeIntegrationSchemeManagerT & GetTimeIntegrationSchemeManager()
        // {
        //     static TimeIntegrationSchemeManagerT instance;
        //     instance.RegisterGlobalCreator( TimeIntegrationScheme::Create );
        //     return instance;
        // }
        

        // bool operator==(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs)
        // {
        //     return (lhs.m_method == rhs.m_method);
        // }
        
        // bool operator<(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs)
        // {
        //     return (lhs.m_method < rhs.m_method);
        // }
        
        // bool TimeIntegrationSchemeKey::opLess::operator()(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs) const
        // {
        //     return (lhs.m_method < rhs.m_method);
        // }

        std::ostream& operator<<( std::ostream& os, const TimeIntegrationSchemeSharedPtr& rhs )
        {
          os << *rhs.get();
          return os;
        }

        std::ostream& operator<<( std::ostream& os, const TimeIntegrationScheme& rhs )
        {
          os << "Time Integration Scheme: " << TimeIntegrationScheme::nameFromMethod( rhs.GetIntegrationMethod() ) << ".\n";
          os << "        Has " << rhs.m_integration_phases.size() << " phases.\n";
          for( int i = 0; i < rhs.m_integration_phases.size(); i++ )
          {
            os << "            - " << TimeIntegrationScheme::nameFromMethod( rhs.m_integration_phases[i]->GetIntegrationMethod() ) << "\n";
          }
          return os;
        }

      TimeIntegrationScheme::ConstDoubleArray &
      TimeIntegrationScheme::TimeIntegrate( const int                                timestep,
                                            const NekDouble                          delta_t,
                                                  TimeIntegrationSolutionSharedPtr & solvector,
                                            const TimeIntegrationSchemeOperators   & op)
      {
        int phases = GetNumIntegrationPhases();
        TimeIntegrationSchemeDataSharedPtr & data = m_integration_phases[ std::min( timestep, phases - 1 ) ];
        return data->TimeIntegrate( delta_t, solvector, op );
      }

      TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr
      TimeIntegrationScheme::
      InitializeScheme( const NekDouble                                 deltaT,
                              TimeIntegrationScheme::ConstDoubleArray & y_0,
                        const NekDouble                                 time,
                        const TimeIntegrationSchemeOperators          & op )
      {
        std::cout << "here: " << m_integration_phases.size() << "\n";
        return m_integration_phases.back()->InitializeData( deltaT, y_0, time, op );
      }

      // Register all the schemes with the Time Integration Scheme Facatory...
      //
#     define REGISTER(x) \
         TimeIntegrationMethod junk##x = GetTimeIntegrationSchemeFactory().RegisterCreatorFunction( e##x, x##TimeIntegrationScheme::create )

//       string x##TimeIntegrationScheme::className = GetTimeIntegrationSchemeFactory().RegisterCreatorFunction( TimeIntegrationMethodMap[ e##x ], x##TimeIntegrationScheme::create )

      REGISTER( AdamsBashforthOrder2 );
      REGISTER( AdamsBashforthOrder3 );
      REGISTER( AdamsMoultonOrder1 );
      REGISTER( AdamsMoultonOrder2 );
      REGISTER( BDFImplicitOrder1 );
      REGISTER( BDFImplicitOrder2 );
      REGISTER( BackwardEuler );
      REGISTER( ClassicalRungeKutta4 );
      REGISTER( CNAB );
      REGISTER( DIRKOrder2 );
      REGISTER( DIRKOrder3 );
      REGISTER( ForwardEuler );
      REGISTER( IMEXdirk_1_2_1 );
      REGISTER( IMEXdirk_1_2_2 );
      REGISTER( IMEXdirk_2_2_2 );
      REGISTER( IMEXdirk_2_3_2 );
      REGISTER( IMEXdirk_2_3_3 );
      REGISTER( IMEXdirk_3_4_3 );
      REGISTER( IMEXdirk_4_4_3 );
      REGISTER( IMEXGear );
      REGISTER( IMEXOrder1 );
      REGISTER( IMEXOrder2 );
      REGISTER( IMEXOrder3 );
      REGISTER( MCNAB );
      REGISTER( RungeKutta2 );
      REGISTER( RungeKutta2_ImprovedEuler );
      REGISTER( RungeKutta2_SSP );
      REGISTER( RungeKutta3_SSP );

    } // end namespace LibUtilities
} // end namespace NekTar
