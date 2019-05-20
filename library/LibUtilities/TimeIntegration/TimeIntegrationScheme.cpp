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
#include <LibUtilities/TimeIntegration/AdamsMoultonOrder2TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/BDFImplicitOrder1TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/BDFImplicitOrder2TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/BackwardEulerTimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/ClassicalRungeKutta4TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/CNABTimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/DIRKOrder2TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/DIRKOrder3TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/ForwardEulerTimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/IMEXdirk_1_2_2TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/IMEXdirk_2_2_2TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/IMEXdirk_2_3_2TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/IMEXdirk_3_4_3TimeIntegrationScheme.h>
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
          os << "Time Integration Scheme: " << TimeIntegrationMethodMap[ rhs.GetIntegrationMethod() ] << ".\n";
          os << "        Has " << rhs.m_integration_phases.size() << " phases.\n";
          for( int i = 0; i < rhs.m_integration_phases.size(); i++ )
          {
            os << "            - " << TimeIntegrationMethodMap[ rhs.m_integration_phases[i]->GetIntegrationSchemeType() ] << "\n";
          }
              
                

          return os;
        }

// FIXME: removing manager for now...
        // // This is a class static function.
        // TimeIntegrationSchemeSharedPtr TimeIntegrationScheme::Create( const TimeIntegrationSchemeKey & key )
        // {
        //     TimeIntegrationSchemeSharedPtr returnval( MemoryManager<TimeIntegrationScheme>::AllocateSharedPtr( key ) );
        //     return returnval;
        // }

//             switch( key.GetIntegrationMethod() )
//             {
//             case eForwardEuler:
//               //            case eAdamsBashforthOrder1:  ForwardEulerTimeIntegrator::SetupScheme( this );         break;
// #if 0
//             case eAdamsBashforthOrder2:  AdamsBashforthOrder2TimeIntegrationScheme::SetupScheme( this ); break;
//             case eAdamsBashforthOrder3:  AdamsBashforthOrder3TimeIntegrator::SetupScheme( this ); break;
//             case eBackwardEuler:
//             case eBDFImplicitOrder1:
//             case eAdamsMoultonOrder1:    BackwardEulerTimeIntegrator::SetupScheme( this );        break;
//             case eIMEXOrder1:            IMEXOrder1TimeIntegrator::SetupScheme( this );           break;
//             case eIMEXOrder2:            IMEXOrder2TimeIntegrator::SetupScheme( this );           break;
//             case eIMEXOrder3:            IMEXOrder3TimeIntegrator::SetupScheme( this );           break;
//             case eAdamsMoultonOrder2:    AdamsMoultonOrder2TimeIntegrator::SetupScheme( this );   break;
//             case eBDFImplicitOrder2:     BDFImplicitOrder2TimeIntegrator::SetupScheme( this );    break;
//             case eMidpoint:
//             case eRungeKutta2:           RungeKutta2TimeIntegrator::SetupScheme( this );          break;
//             case eRungeKutta2_ImprovedEuler:
//             case eRungeKutta2_SSP:       RungeKutta2_ImprovedEulerTimeIntegrator::SetupScheme( this ); break;
//             case eRungeKutta3_SSP:       RungeKutta3_SSPTimeIntegrator::SetupScheme( this );      break;
//             case eClassicalRungeKutta4:
//             case eRungeKutta4:           ClassicalRungeKutta4TimeIntegrator::SetupScheme( this ); break;
//             case eDIRKOrder2:            DirkOrder2TimeIntegrator::SetupScheme( this );           break;
//             case eDIRKOrder3:            DirkOrder3TimeIntegrator::SetupScheme( this );           break;
//             case eIMEXdirk_2_3_2:        IMEXDirk_2_3_2TimeIntegrator::SetupScheme( this );       break;
//             case eIMEXdirk_3_4_3:        IMEXDirk_3_4_3TimeIntegrator::SetupScheme( this );       break;
//             case eCNAB:                  CNABTimeIntegrator::SetupScheme( this );                 break;
//             case eIMEXGear:              IMEXGearTimeIntegrator::SetupScheme( this );             break;
//             case eMCNAB:                 MCNABTimeIntegrator::SetupScheme( this );                break;
//             case eIMEXdirk_2_2_2:        IMEXDirk_2_2_2TimeIntegrator::SetupScheme( this );       break;
//             case eIMEXdirk_2_3_3:
//                 {
//                   cout << "error with eIMEXdirk_2_3_3\n";
//                   throw "not done yet";
//                 }
//                 break;
//             case eIMEXdirk_1_1_1:
//                 {
//                   cout << "error with eIMEXdirk_1_1_1\n";
//                   throw "not done yet";
//                 }
//                 break;
//             case eIMEXdirk_1_2_1:
//                 {
//                   cout << "error with eIMEXdirk_1_2_1\n";
//                   throw "not done yet";
//                 }
//                 break;
//             case eIMEXdirk_1_2_2:        IMEXDirk_1_2_2TimeIntegrator::SetupScheme( this );       break;
//             case eIMEXdirk_4_4_3:
//                 {
//                   cout << "error with eIMEXdirk_4_4_3\n";
//                   throw "not done yet";
//                 }
//                 break;
// #endif
//             default:
//                 {
//                     NEKERROR( ErrorUtil::efatal, "Invalid Time Integration Scheme" );
//                 }
//             }
            
//             // m_firstStageEqualsOldSolution = CheckIfFirstStageEqualsOldSolution(m_A,m_B,m_U,m_V);
//             // m_lastStageEqualsNewSolution  = CheckIfLastStageEqualsNewSolution(m_A,m_B,m_U,m_V);

//             // ASSERTL1(VerifyIntegrationSchemeType(m_schemeType,m_A,m_B,m_U,m_V),
//             //          "Time integration scheme coefficients do not match its type");
//         }


        // bool TimeIntegrationScheme::
        // VerifyIntegrationSchemeType( TimeIntegrationSchemeType type,
        //                              const Array<OneD, const Array<TwoD, NekDouble> >& A,
        //                              const Array<OneD, const Array<TwoD, NekDouble> >& B,
        //                              const Array<TwoD, const NekDouble>& U,
        //                              const Array<TwoD, const NekDouble>& V ) const
        // {
        //     int i;
        //     int j;
        //     int m;
        //     int  IMEXdim = A.num_elements();
        //     int  dim     = A[0].GetRows();

        //     Array<OneD, TimeIntegrationSchemeType> vertype( IMEXdim, eExplicit );

        //     for(m = 0; m < IMEXdim; m++)
        //     {
        //         for(i = 0; i < dim; i++)
        //         {
        //             if( fabs(A[m][i][i]) > NekConstants::kNekZeroTol )
        //             {
        //                 vertype[m] = eDiagonallyImplicit;
        //             }
        //         }
                
        //         for(i = 0; i < dim; i++)
        //         {
        //             for(j = i+1; j < dim; j++)
        //             {
        //                 if( fabs(A[m][i][j]) > NekConstants::kNekZeroTol )
        //                 {
        //                     vertype[m] = eImplicit;
        //                     ASSERTL1(false,"Fully Implicit schemes cannnot be handled by the TimeIntegrationScheme class");
        //                 }
        //             }
        //         }
        //     }

        //     if(IMEXdim == 2)
        //     {
        //         ASSERTL1(B.num_elements()==2,"Coefficient Matrix B should have an implicit and explicit part for IMEX schemes");
        //         if((vertype[0] == eDiagonallyImplicit) &&
        //            (vertype[1] == eExplicit))
        //         {
        //             vertype[0] = eIMEX;
        //         }
        //         else
        //         {
        //             ASSERTL1(false,"This is not a proper IMEX scheme");
        //         }
        //     }

        //     return (vertype[0] == type);
        // }

        // TimeIntegrationSolutionSharedPtr 
        // TimeIntegrationScheme::InitializeScheme( const NekDouble                        timestep,
        //                                                ConstDoubleArray               & y_0,
        //                                          const NekDouble                        time,
        //                                          const TimeIntegrationSchemeOperators & op)
        // {
        //   cout << "THIS IS " << this << "\n";
        //   cout << "aaaaaaaaa: " << m_schemeKey << "\n";

        //     // create a TimeIntegrationSolution object based upon the
        //     // initial value. Initialise all other multi-step values
        //     // and derivatives to zero
        //     TimeIntegrationSolutionSharedPtr y_out = 
        //         MemoryManager<TimeIntegrationSolution>::AllocateSharedPtr( m_schemeKey, y_0, time, timestep );

        //   cout << "AAAAAAAAA\n";
        //     if( GetIntegrationSchemeType() == eExplicit)
        //     {
        //         // ensure initial solution is in correct space
        //         op.DoProjection(y_0,y_out->UpdateSolution(),time);
        //     }
        //   cout << "aaaaaaaaa\n";

        //     // calculate the initial derivative, if is part of the
        //     // solution vector of the current scheme
        //     if( m_numMultiStepDerivs > 0 )
        //     {
        //         if(m_timeLevelOffset[m_numMultiStepValues] == 0)
        //         {
        //             int i;
        //             int nvar    = y_0.num_elements();
        //             int npoints = y_0[0].num_elements();
        //             DoubleArray f_y_0(nvar);
        //             for(i = 0; i < nvar; i++)
        //             {
        //                 f_y_0[i] = Array<OneD,NekDouble>(npoints);
        //             }
        //             // calculate the derivative of the initial value
        //             op.DoOdeRhs(y_0,f_y_0,time);
                    
        //             // multiply by the step size
        //             for(i = 0; i < nvar; i++)
        //             {
        //                 Blas::Dscal(npoints,timestep,f_y_0[i].get(),1);
        //             }
        //             y_out->SetDerivative(0,f_y_0,timestep);
        //         }
        //     }
           
        //     cout << "bbbbbbbbb\n";
        //     return y_out;
        // }
        
      TimeIntegrationSchemeType
      TimeIntegrationScheme::GetIntegrationSchemeType() const { return m_integration_phases.back()->GetIntegrationSchemeType(); }

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
         string x##TimeIntegrationScheme::className = GetTimeIntegrationSchemeFactory().RegisterCreatorFunction( TimeIntegrationMethodMap[ e##x ], x##TimeIntegrationScheme::create )

      REGISTER( AdamsBashforthOrder2 );
      REGISTER( AdamsBashforthOrder3 );
      REGISTER( AdamsMoultonOrder2 );
      REGISTER( BDFImplicitOrder1 );
      REGISTER( BDFImplicitOrder2 );
      REGISTER( BackwardEuler );
      REGISTER( ClassicalRungeKutta4 );
      REGISTER( CNAB );
      REGISTER( DIRKOrder2 );
      REGISTER( DIRKOrder3 );
      REGISTER( ForwardEuler );
      REGISTER( IMEXdirk_1_2_2 );
      REGISTER( IMEXdirk_2_2_2 );
      REGISTER( IMEXdirk_2_3_2 );
      REGISTER( IMEXdirk_3_4_3 );
      REGISTER( IMEXGear );
      REGISTER( IMEXOrder1 );
      REGISTER( IMEXOrder2 );
      REGISTER( IMEXOrder3 );
      REGISTER( MCNAB );
      REGISTER( RungeKutta2 );
      REGISTER( RungeKutta2_ImprovedEuler );
      REGISTER( RungeKutta2_SSP );
      REGISTER( RungeKutta3_SSP );

      //    string AdamsBashforthOrder2TimeIntegrationScheme::className = GTIF().RegisterCreatorFunction( TimeIntegrationMethodMap[ eAdamsBashforthOrder2 ], AdamsBashforthOrder2TimeIntegrationScheme::create );
      //    string IMEXOrder1TimeIntegrationScheme::className           = GTIF().RegisterCreatorFunction( TimeIntegrationMethodMap[ eIMEXOrder1 ],           IMEXOrder1TimeIntegrationScheme::create );
      //    string IMEXOrder2TimeIntegrationScheme::className           = GTIF().RegisterCreatorFunction( TimeIntegrationMethodMap[ eIMEXOrder2 ],           IMEXOrder2TimeIntegrationScheme::create );


      //    string AdamsBashforthOrder3TimeIntegrator::className = GTIF().RegisterCreatorFunction( "AdamsBashforthOrder3", AdamsBashforthOrder3TimeIntegrator::create );
      //    string AdamsMoultonOrder2TimeIntegrator::className   = GTIF().RegisterCreatorFunction( "AdamsMoultonOrder2",   AdamsMoultonOrder2TimeIntegrator::create );
      //    string BackwardEulerTimeIntegrator::className        = GTIF().RegisterCreatorFunction( "BackwardEuler",        BackwardEulerTimeIntegrator::create );
      //    string BDFImplicitOrder1TimeIntegrator::className    = GTIF().RegisterCreatorFunction( "BDFImplicitOrder1",    BDFImplicitOrder1TimeIntegrator::create );
      //    string BDFImplicitOrder2TimeIntegrator::className    = GTIF().RegisterCreatorFunction( "BDFImplicitOrder2",    BDFImplicitOrder2TimeIntegrator::create );
      //    string ClassicalRungeKutta4TimeIntegrator::className = GTIF().RegisterCreatorFunction( "ClassicalRungeKutta4", ClassicalRungeKutta4TimeIntegrator::create );
      //    string CNABTimeIntegrator::className                 = GTIF().RegisterCreatorFunction( "CNAB",                 CNABTimeIntegrator::create );
      //    string DirkOrder2TimeIntegrator::className           = GTIF().RegisterCreatorFunction( "DIRKOrder2",           DirkOrder2TimeIntegrator::create );
      //    string DirkOrder3TimeIntegrator::className           = GTIF().RegisterCreatorFunction( "DIRKOrder3",           DirkOrder3TimeIntegrator::create );
      //    string ForwardEulerTimeIntegrator::className         = GTIF().RegisterCreatorFunction( "ForwardEuler",         ForwardEulerTimeIntegrator::create );
      //    string IMEXDirk_1_2_2TimeIntegrator::className       = GTIF().RegisterCreatorFunction( "IMEXdirk_1_2_2",       IMEXDirk_1_2_2TimeIntegrator::create );
      //    string IMEXDirk_2_2_2TimeIntegrator::className       = GTIF().RegisterCreatorFunction( "IMEXdirk_2_2_2",       IMEXDirk_2_2_2TimeIntegrator::create );
      //    string IMEXDirk_2_3_2TimeIntegrator::className       = GTIF().RegisterCreatorFunction( "IMEXdirk_2_3_2",       IMEXDirk_2_3_2TimeIntegrator::create );
      //    string IMEXDirk_3_4_3TimeIntegrator::className       = GTIF().RegisterCreatorFunction( "IMEXdirk_3_4_3",       IMEXDirk_3_4_3TimeIntegrator::create );
      //    string IMEXGearTimeIntegrator::className             = GTIF().RegisterCreatorFunction( "IMEXGear",             IMEXGearTimeIntegrator::create );
      //    string IMEXOrder1TimeIntegrator::className           = GTIF().RegisterCreatorFunction( "IMEXOrder1",           IMEXOrder1TimeIntegrator::create );
      //    string IMEXOrder2TimeIntegrator::className           = GTIF().RegisterCreatorFunction( "IMEXOrder2",           IMEXOrder2TimeIntegrator::create );
      //    string IMEXOrder3TimeIntegrator::className           = GTIF().RegisterCreatorFunction( "IMEXOrder3",           IMEXOrder3TimeIntegrator::create );
      //    string MCNABTimeIntegrator::className                = GTIF().RegisterCreatorFunction( "MCNAB",                MCNABTimeIntegrator::create );
      //    string RungeKutta2TimeIntegrator::className          = GTIF().RegisterCreatorFunction( "RungeKutta2",          RungeKutta2TimeIntegrator::create );
      //    string RungeKutta2_ImprovedEulerTimeIntegrator::className = GTIF().RegisterCreatorFunction( "RungeKutta2_ImprovedEuler", RungeKutta2_ImprovedEulerTimeIntegrator::create );
      //    string RungeKutta2_SSPTimeIntegrator::className      = GTIF().RegisterCreatorFunction( "RungeKutta2_SSP", RungeKutta2_SSPTimeIntegrator::create );
      //    string RungeKutta3_SSPTimeIntegrator::className      = GTIF().RegisterCreatorFunction( "RungeKutta3_SSP", RungeKutta3_SSPTimeIntegrator::create );


    } // end namespace LibUtilities
} // end namespace NekTar
