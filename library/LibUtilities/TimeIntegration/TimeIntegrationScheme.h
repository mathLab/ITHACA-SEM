#pragma once

///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationScheme.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Header file of time integration scheme base class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_FOUNDATIONS_TIMEINTEGRATIONSCHEME_H
#define NEKTAR_LIB_UTILITIES_FOUNDATIONS_TIMEINTEGRATIONSCHEME_H

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
//#include <LibUtilities/BasicUtils/NekManager.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

#include <string>

///////////////////////////////////////////////////////////////////////////////

#define LUE LIB_UTILITIES_EXPORT

///////////////////////////////////////////////////////////////////////////////

namespace Nektar
{
    namespace LibUtilities
    {
        // Forward declaration of some of the classes in this file...
        class TimeIntegrationScheme;
        class TimeIntegrationSchemeData;
        class TimeIntegrationSolution;

        // typedefs


      // FIXME: Why a shared pointer?  when does the scheme go out of scope... who should delete it?


        typedef std::shared_ptr<TimeIntegrationScheme> TimeIntegrationSchemeSharedPtr;


        // =========================================================================
        // ==== ENUM LIST OF ALL SUPPORTED INTEGRATION SCHEMES
        // =========================================================================
        enum TimeIntegrationMethod
        {
            eNoTimeIntegrationMethod,
            eAdamsBashforthOrder1,            //!< Adams-Bashforth Forward multi-step scheme of order 1
            eAdamsBashforthOrder2,            //!< Adams-Bashforth Forward multi-step scheme of order 2
            eAdamsBashforthOrder3,            //!< Adams-Bashforth Forward multi-step scheme of order 3
            eAdamsBashforthOrder4,            //!< Adams-Bashforth Forward multi-step scheme of order 4
            eAdamsMoultonOrder1,              //!< Adams-Moulton Forward multi-step scheme of order 1
            eAdamsMoultonOrder2,              //!< Adams-Moulton Forward multi-step scheme of order 2
            eBDFImplicitOrder1,               //!< BDF multi-step scheme of order 1 (implicit)
            eBDFImplicitOrder2,               //!< BDF multi-step scheme of order 2 (implicit)
            eClassicalRungeKutta4,            //!< Runge-Kutta multi-stage scheme 4th order explicit (old name)
            eRungeKutta4,                     //!< Classical RungeKutta4 method (new name for eClassicalRungeKutta4)
            eRungeKutta5,                     //!< RungeKutta5 method
            eRungeKutta3_SSP,                 //!< Nonlinear SSP RungeKutta3 explicit
            eRungeKutta2_ImprovedEuler,       //!< Improved RungeKutta2 explicit (old name meaning Heun's method)
            eRungeKutta2_SSP,                 //!< Nonlinear SSP RungeKutta2 explicit (surrogate for eRungeKutta2_ImprovedEuler)
            eForwardEuler,                    //!< Forward Euler scheme
            eBackwardEuler,                   //!< Backward Euler scheme
            eIMEXOrder1,                      //!< IMEX 1st order scheme using Euler Backwards/Euler Forwards
            eIMEXOrder2,                      //!< IMEX 2nd order scheme using Backward Different Formula & Extrapolation
            eIMEXOrder3,                      //!< IMEX 3rd order scheme using Backward Different Formula & Extrapolation
            eIMEXOrder4,                      //!< IMEX 4th order scheme using Backward Different Formula & Extrapolation
            eMidpoint,                        //!< midpoint method (old name)
            eRungeKutta2,                     //!< Classical RungeKutta2 method (new name for eMidpoint)
            eDIRKOrder2,                      //!< Diagonally Implicit Runge Kutta scheme of order 2
            eDIRKOrder3,                      //!< Diagonally Implicit Runge Kutta scheme of order 3
            eCNAB,                            //!< Crank-Nicolson/Adams-Bashforth Order 2 (CNAB)
            eIMEXGear,                        //!< IMEX Gear Order 2
            eMCNAB,                           //!< Modified Crank-Nicolson/Adams-Bashforth Order 2 (MCNAB)
            eIMEXdirk_1_1_1,                  //!< Forward-Backward Euler IMEX DIRK(1,1,1)
            eIMEXdirk_1_2_1,                  //!< Forward-Backward Euler IMEX DIRK(1,2,1)
            eIMEXdirk_1_2_2,                  //!< Implicit-Explicit Midpoint IMEX DIRK(1,2,2)
            eIMEXdirk_2_2_2,                  //!< L-stable, two stage, second order IMEX DIRK(2,2,2)
            eIMEXdirk_2_3_2,                  //!< L-stable, three stage, third order IMEX DIRK(3,4,3)
            eIMEXdirk_2_3_3,                  //!< L-stable, two stage, third order IMEX DIRK(2,3,3)
            eIMEXdirk_3_4_3,                  //!< L-stable, three stage, third order IMEX DIRK(3,4,3)
            eIMEXdirk_4_4_3,                  //!< L-stable, four stage, third order IMEX DIRK(4,4,3)
            SIZE_TimeIntegrationMethod        //!< Length of enum list
        };

      ///////////////////////////////////////////////////////////////////////////////////
      // Provide using code the ability to get the TimeIntegrator Factory so that it can
      // create TimeIntegrators.

      /// Datatype of the NekFactory used to instantiate classes derived from the EquationSystem class.
      typedef NekFactory < std::string, TimeIntegrationScheme > TimeIntegrationSchemeFactory;


      // Allows a code to create a TimeIntegrator. Usually used like this:
      //
      //    LibUtilities::TimeIntegrationSchemeSharedPtr timeIntegrationScheme = LibUtilities::GetTimeIntegrationSchemeFactory().CreateInstance( "IMEXOrder1" );
      //
      LUE TimeIntegrationSchemeFactory & GetTimeIntegrationSchemeFactory();

      //////////////////////////////////////////////////////////////////////////////////////

        enum TimeIntegrationSchemeType
        {
            eNoTimeIntegrationSchemeType,
            eExplicit,              //!< Formally explicit scheme
            eDiagonallyImplicit,    //!< Diagonally implicit scheme (e.g. the DIRK schemes)
            eIMEX,                  //!< Implicit Explicit General Linear Method
            eImplicit,              //!< Fully implicit scheme
        };

        const char* const TimeIntegrationSchemeTypeMap[] = 
        {
            "NoTimeIntegrationSchemeType",
            "Explicit",
            "DiagonallyImplicit",
            "IMEX",
            "Implicit"
        };

        // =====================================================================

        // =====================================================================
        // ==== DEFINITION OF THE CLASS TimeIntegrationSchemeOperators
        // =====================================================================
        class TimeIntegrationSchemeOperators
        {
        public:
            typedef const Array<OneD, const Array<OneD, NekDouble> > InArrayType;
            typedef       Array<OneD,       Array<OneD, NekDouble> > OutArrayType;
            
            typedef std::function< void (InArrayType&, OutArrayType&, const NekDouble) >                  FunctorType1;
            typedef std::function< void (InArrayType&, OutArrayType&, const NekDouble, const NekDouble) > FunctorType2;

            typedef const FunctorType1& ConstFunctorType1Ref;
            typedef const FunctorType2& ConstFunctorType2Ref;

            typedef Array<OneD, FunctorType1> FunctorType1Array;
            typedef Array<OneD, FunctorType2> FunctorType2Array;
            
            TimeIntegrationSchemeOperators(void):
            m_functors1(4),
            m_functors2(1)
            {
            }

            template<typename FuncPointerT, typename ObjectPointerT> 
            void DefineOdeRhs( FuncPointerT func, ObjectPointerT obj )
            {
                m_functors1[0] = std::bind( func, obj, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 );
            }

            template<typename FuncPointerT, typename ObjectPointerT> 
            void DefineOdeExplicitRhs( FuncPointerT func, ObjectPointerT obj )
            {
                m_functors1[1] = std::bind( func, obj, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 );
            }

            template<typename FuncPointerT, typename ObjectPointerT> 
            void DefineOdeImplicitRhs( FuncPointerT func, ObjectPointerT obj )
            {
                m_functors1[2] = std::bind( func, obj, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 );
            }

            template<typename FuncPointerT, typename ObjectPointerT> 
            void DefineProjection( FuncPointerT func, ObjectPointerT obj )
            {
                m_functors1[3] = std::bind( func, obj, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 );
            }

            template<typename FuncPointerT, typename ObjectPointerT> 
            void DefineImplicitSolve(FuncPointerT func, ObjectPointerT obj)
            {
                m_functors2[0] = std::bind( func, obj, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4 );
            }

            inline void DoOdeRhs(       InArrayType  & inarray, 
                                        OutArrayType & outarray, 
                                  const NekDouble      time ) const
            {
                ASSERTL1(m_functors1[0],"OdeRhs should be defined for this time integration scheme");
                m_functors1[0]( inarray, outarray, time );
            }
            
            inline void DoOdeExplicitRhs(InArrayType     &inarray, 
                                         OutArrayType    &outarray, 
                                         const NekDouble time) const
            {
                ASSERTL1(m_functors1[1],"OdeExplicitRhs should be defined for this time integration scheme");
                m_functors1[1](inarray,outarray,time);
            }
            
            inline void DoOdeImplicitRhs(InArrayType     &inarray, 
                                        OutArrayType    &outarray, 
                                        const NekDouble time) const
            {
                ASSERTL1(m_functors1[2],"OdeImplictRhs should be defined for this time integration scheme");
                m_functors1[2](inarray,outarray,time);
            }

            inline void DoProjection(InArrayType     &inarray, 
                                     OutArrayType    &outarray, 
                                     const NekDouble time) const
            {
                ASSERTL1(m_functors1[3],"Projection operation should be defined for this time integration scheme");
                m_functors1[3](inarray,outarray,time);
            }
            
            inline void DoImplicitSolve(       InArrayType  & inarray, 
                                               OutArrayType & outarray, 
                                         const NekDouble      time, 
                                         const NekDouble      lambda ) const
            {
                ASSERTL1(m_functors2[0],"ImplicitSolve should be defined for this time integration scheme");
                m_functors2[0]( inarray, outarray, time, lambda );
            }


        protected:
            FunctorType1Array m_functors1;
            FunctorType2Array m_functors2;

        private:

        }; // end class TimeIntegrationSchemeOperators()

        // =====================================================================
        // ==== DEFINITION OF THE CLASS TimeIntegrationScheme
        // =====================================================================

        // FIXME: Dd: rename to TimeIntegrationSchemeBase (<- Add BASE)????

        class TimeIntegrationScheme
        {
        public:

          // typedefs
          typedef std::shared_ptr<TimeIntegrationSchemeData>                  TimeIntegrationSchemeDataSharedPtr;
          typedef std::vector<TimeIntegrationSchemeDataSharedPtr>             TimeIntegrationSchemeDataVector;
          typedef std::vector<TimeIntegrationSchemeDataSharedPtr>::iterator   TimeIntegrationSchemeDataVectorIter; 

          typedef std::shared_ptr<TimeIntegrationSolution>                TimeIntegrationSolutionSharedPtr;
          typedef std::vector<TimeIntegrationSolutionSharedPtr>           TimeIntegrationSolutionVector; 
          typedef std::vector<TimeIntegrationSolutionSharedPtr>::iterator TimeIntegrationSolutionVectorIter; 


                typedef const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > > ConstTripleArray;
                typedef       Array<OneD,       Array<OneD, Array<OneD, NekDouble> > > TripleArray;
                typedef const Array<OneD, const Array<OneD, NekDouble> >               ConstDoubleArray;
                typedef       Array<OneD,       Array<OneD, NekDouble> >               DoubleArray;
                typedef const Array<OneD, const NekDouble >                            ConstSingleArray;
                typedef       Array<OneD,       NekDouble >                            SingleArray;
                typedef std::function< void (ConstDoubleArray&, DoubleArray&, const NekDouble) >                  FunctorType1;
                typedef std::function< void (ConstDoubleArray&, DoubleArray&, const NekDouble, const NekDouble) > FunctorType2;


          static TimeIntegrationMethod methodFromName( const std::string & name );
          static std::string           nameFromMethod( const TimeIntegrationMethod method );

          unsigned int GetNumIntegrationPhases() const { return m_integration_phases.size(); }

          LUE
          TimeIntegrationSolutionSharedPtr InitializeScheme( const NekDouble                                 deltaT,
                                                                   TimeIntegrationScheme::ConstDoubleArray & y_0,
                                                             const NekDouble                                 time,
                                                             const TimeIntegrationSchemeOperators          & op );
          LUE
          TimeIntegrationScheme::ConstDoubleArray &
          TimeIntegrate( const int                                timestep,
                         const NekDouble                          delta_t,
                               TimeIntegrationSolutionSharedPtr & solvector,
                         const TimeIntegrationSchemeOperators   & op);

          // The "Method" of the "master scheme" is actually the method of the last phase (sub-scheme)...
          // See child classes for implementation.
          virtual TimeIntegrationMethod GetIntegrationMethod() const = 0;

          // FIXME: Dd: this should be protected...
          // NOTE: FIXME: It seems to me that it doesn't make sense to ask a (multi-phase) Scheme what its type is... as
          //              each phase can have a different type... Ask Chris/Mike about this.
          //              For now, returning type of last phase... 
          TimeIntegrationSchemeType GetIntegrationSchemeType() const;

        protected:// <- Dd: Now going to use as a base class, so igore -> Dd: don't think anything is inheriting from this, so don't need protected...
          // private:

          friend class TimeIntegrationSolution;
          friend class TimeIntegrationSchemeData;

            virtual ~TimeIntegrationScheme() {}
 
          friend std::ostream& operator<<( std::ostream& os, const TimeIntegrationScheme& rhs );
          friend std::ostream& operator<<( std::ostream& os, const TimeIntegrationSchemeSharedPtr& rhs );

          friend std::ostream& operator<<( std::ostream& os, const TimeIntegrationSchemeData & rhs );

          //friend std::ostream& operator<<( std::ostream& os, const TimeIntegrationScheme& rhs );
          //          friend std::ostream& operator<<( std::ostream& os, const TimeIntegrationSchemeSharedPtr& rhs );
          
            /**
             * \brief Explicit integration of an ODE.
             *
             * This function explicitely perfroms a signle integration step of the ODE system:
             * \f[
             * \frac{d\boldsymbol{y}}{dt}=\boldsymbol{f}(t,\boldsymbol{y})
             * \f]
             *
             * \param timestep The size of the timestep, i.e. \f$\Delta t\f$.
             * \param f an object of the class FuncType, where FuncType should have a method FuncType::ODEforcing
             *       to evaluate the right hand side \f$f(t,\boldsymbol{y})\f$ of the ODE.
             * \param y on input: the vectors \f$\boldsymbol{y}^{[n-1]}\f$ and \f$t^{[n-1]}\f$ (which corresponds to the 
             *    solution at the old time level)
             * \param y on output: the vectors \f$\boldsymbol{y}^{[n]}\f$ and \f$t^{[n]}\f$ (which corresponds to the 
             *    solution at the old new level)
             * \return The actual solution \f$\boldsymbol{y}^{n}\f$ at the new time level 
             *    (which in fact is also embedded in the argument y).
             */
            LUE ConstDoubleArray& TimeIntegrate( const NekDouble                          timestep,          
                                                       TimeIntegrationSolutionSharedPtr & y      ,
                                                 const TimeIntegrationSchemeOperators   & op );

          // TimeIntegrationSchemeKey  m_schemeKey; 
          //TimeIntegrationSchemeType m_schemeType;
          //unsigned int              m_numsteps;   //< Number of steps in multi-step component. 
          // unsigned int              m_numstages;  //< Number of stages in multi-stage component.

        public: // FIXME m_integration_phases should not be public...
            TimeIntegrationSchemeDataVector m_integration_phases; // Was called m_intSchemes
        protected: // FIXME

            bool m_firstStageEqualsOldSolution; //< Optimisation-flag 
            bool m_lastStageEqualsNewSolution;  //< Optimisation-flag

            // Array<OneD, Array<TwoD,NekDouble> > m_A;
            // Array<OneD, Array<TwoD,NekDouble> > m_B;
            // Array<TwoD,NekDouble>               m_U;
            // Array<TwoD,NekDouble>               m_V;

          // bool m_initialised;   /// bool to identify if array has been initialised 
            // int  m_nvar;          /// The number of variables in integration scheme. 
            // int  m_npoints;       /// The size of inner data which is stored for reuse. 
            // DoubleArray m_Y;      /// Array containing the stage values 
            // DoubleArray m_tmp;    /// explicit right hand side of each stage equation

            // TripleArray m_F;      /// Array corresponding to the stage Derivatives 
            // TripleArray m_F_IMEX; /// Used to store the Explicit stage derivative of IMEX schemes

            // NekDouble   m_T;     ///  Time at the different stages

          //            template <typename> friend class Nektar::MemoryManager;

          //            friend TimeIntegrationSchemeManagerT & GetTimeIntegrationSchemeManager();

          // LUE static TimeIntegrationSchemeSharedPtr Create( const TimeIntegrationSchemeKey & key );

            // This should never be used directly... only used by child classes...
          TimeIntegrationScheme() {}
            
            // These should never be called
          // TimeIntegrationScheme() : m_method( eNoTimeIntegrationMethod ) { NEKERROR(ErrorUtil::efatal,"Default Constructor for the TimeIntegrationScheme class should not be called"); }
          TimeIntegrationScheme( const TimeIntegrationScheme & in ) {
              boost::ignore_unused(in);
              NEKERROR(ErrorUtil::efatal,"Copy Constructor for the TimeIntegrationScheme class should not be called");
          }

            inline int GetFirstDim(ConstTripleArray &y) const
            {
                return y[0].num_elements();
            }

            inline int GetSecondDim(ConstTripleArray &y) const
            {
                return y[0][0].num_elements();
            }

            LUE bool CheckTimeIntegrateArguments(const NekDouble timestep,
                                                   ConstTripleArray               &y_old  ,
                                                   ConstSingleArray               &t_old  ,
                                                   TripleArray                    &y_new  ,
                                                   SingleArray                    &t_new  ,
                                             const TimeIntegrationSchemeOperators &op) const;

          // !!! Always make sure that this matches TimeIntegrationMethod enum... !!!
          //
          static const char* const TimeIntegrationMethodMap[ 36 ];


        }; // end class TimeIntegrationScheme


      LUE std::ostream& operator<<( std::ostream& os, const TimeIntegrationScheme& rhs );
      LUE std::ostream& operator<<( std::ostream& os, const TimeIntegrationSchemeSharedPtr& rhs );
        
        // =========================================================================

    } // end of namespace LibUtilities
} // end of namespace Nektar

#endif
