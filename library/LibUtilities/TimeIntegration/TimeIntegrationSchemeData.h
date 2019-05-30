#pragma once

///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationSchemeData.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2019 Division of Applied Mathematics, Brown University (USA),
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
// Description: Header file of time integration scheme data class
//
// The TimeIntegrationSchemeData class should only be used by the 
// TimeIntegrationScheme class.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

///////////////////////////////////////////////////////////////////////////////

#define LUE LIB_UTILITIES_EXPORT

///////////////////////////////////////////////////////////////////////////////

namespace Nektar
{
    namespace LibUtilities
    {
        // =====================================================================
        // ==== DEFINITION OF THE CLASS TimeIntegrationScheme
        // =====================================================================
        class TimeIntegrationSchemeData
        {
        public:
          
            TimeIntegrationSchemeData( const TimeIntegrationScheme * parent ) : m_parent( parent ), m_initialised( false ) {}
            ~TimeIntegrationSchemeData() {}

        private: // Note: All functions are private as they should only be used by their friend
                 // classes TimeIntegrationScheme and TimeIntegrationSolution

            friend std::ostream& operator<<( std::ostream& os, const TimeIntegrationScheme& rhs );
            friend std::ostream& operator<<( std::ostream& os, const TimeIntegrationSchemeSharedPtr& rhs );

            friend std::ostream& operator<<( std::ostream& os, const TimeIntegrationSchemeData & rhs );
            friend std::ostream& operator<<( std::ostream& os, const TimeIntegrationScheme::TimeIntegrationSchemeDataSharedPtr& rhs );

            friend class TimeIntegrationScheme;
            friend class TimeIntegrationSolution;

            friend class AdamsBashforthOrder2TimeIntegrationScheme;
            friend class AdamsBashforthOrder3TimeIntegrationScheme;
            friend class AdamsMoultonOrder1TimeIntegrationScheme;
            friend class AdamsMoultonOrder2TimeIntegrationScheme;
            friend class BackwardEulerTimeIntegrationScheme;
            friend class BDFImplicitOrder2TimeIntegrationScheme;
            friend class ClassicalRungeKutta4TimeIntegrationScheme;
            friend class CNABTimeIntegrationScheme;
            friend class DIRKOrder2TimeIntegrationScheme;
            friend class DIRKOrder3TimeIntegrationScheme;
            friend class ForwardEulerTimeIntegrationScheme;
            friend class IMEXdirk_1_2_1TimeIntegrationScheme;
            friend class IMEXdirk_1_2_2TimeIntegrationScheme;
            friend class IMEXdirk_2_2_2TimeIntegrationScheme;
            friend class IMEXdirk_2_3_2TimeIntegrationScheme;
            friend class IMEXdirk_2_3_3TimeIntegrationScheme;
            friend class IMEXdirk_3_4_3TimeIntegrationScheme;
            friend class IMEXdirk_4_4_3TimeIntegrationScheme;
            friend class IMEXGearTimeIntegrationScheme;
            friend class IMEXOrder1TimeIntegrationScheme;
            friend class IMEXOrder2TimeIntegrationScheme;
            friend class IMEXOrder3TimeIntegrationScheme;
            friend class MCNABTimeIntegrationScheme;
            friend class RungeKutta2TimeIntegrationScheme;
            friend class RungeKutta2_ImprovedEulerTimeIntegrationScheme;
            friend class RungeKutta2_SSPTimeIntegrationScheme;
            friend class RungeKutta3_SSPTimeIntegrationScheme;

            //////////////////////////////////////////////////////////////


            // FIXME: Dd: move to the right spot:
            const TimeIntegrationScheme * m_parent;

            inline TimeIntegrationMethod     GetIntegrationMethod() const     { return m_method; }
            inline TimeIntegrationSchemeType GetIntegrationSchemeType() const { return m_schemeType; }

            inline NekDouble A(const unsigned int i, const unsigned int j) const { return m_A[0][i][j]; }
            inline NekDouble B(const unsigned int i, const unsigned int j) const { return m_B[0][i][j]; }
            inline NekDouble U(const unsigned int i, const unsigned int j) const { return m_U[i][j]; }
            inline NekDouble V(const unsigned int i, const unsigned int j) const { return m_V[i][j]; }

            inline NekDouble A_IMEX(const unsigned int i, const unsigned int j) const { return m_A[1][i][j]; }
            inline NekDouble B_IMEX(const unsigned int i, const unsigned int j) const { return m_B[1][i][j]; }

            inline unsigned int GetNstages(void) const { return m_numstages; }

            inline unsigned int GetNmultiStepValues() const { return m_numMultiStepValues; }
            inline unsigned int GetNmultiStepDerivs() const { return m_numMultiStepDerivs; }
            inline const Array<OneD, const unsigned int>& GetTimeLevelOffset() const { return m_timeLevelOffset; }

            typedef TimeIntegrationScheme::SingleArray                        SingleArray;
            typedef TimeIntegrationScheme::ConstSingleArray                   ConstSingleArray;
            typedef TimeIntegrationScheme::DoubleArray                        DoubleArray;
            typedef TimeIntegrationScheme::ConstDoubleArray                   ConstDoubleArray;
            typedef TimeIntegrationScheme::TripleArray                        TripleArray;
            typedef TimeIntegrationScheme::ConstTripleArray                   ConstTripleArray;
            typedef TimeIntegrationScheme::TimeIntegrationSchemeDataSharedPtr TimeIntegrationSchemeDataSharedPtr;
            typedef TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr   TimeIntegrationSolutionSharedPtr;

            /**
             * \brief This function initialises the time integration
             * scheme
             *
             * Given the solution at the initial time level
             * \f$\boldsymbol{y}(t_0)\f$, this function generates the
             * vectors \f$\boldsymbol{y}^{[n]}\f$ and \f$t^{[n]}\f$
             * needed to evaluate the time integration scheme
             * formulated as a General Linear Method. These vectors
             * are embedded in an object of the class
             * TimeIntegrationSolution. This class is the abstraction
             * of the input and output vectors of the General Linear
             * Method.
             *
             * For single-step methods, this function is trivial as it
             * just wraps a TimeIntegrationSolution object around the
             * given initial values and initial time.  However, for
             * multistep methods, actual time stepping is being done
             * to evaluate the necessary parameters at multiple time
             * levels needed to start the actual integration.
             *
             * \param timestep The size of the timestep, i.e. \f$\Delta t\f$.
             * \param time on input: the initial time, i.e. \f$t_0\f$.
             * \param time on output: the new time-level after initialisation, in general this yields
             *       \f$t = t_0 + (r-1)\Delta t\f$ where \f$r\f$ is the number of steps of the multi-step method.
             * \param nsteps on output: he number of initialisation steps required. In general this corresponds to \f$r-1\f$
             *       where \f$r\f$ is the number of steps of the multi-step method.
             * \param f an object of the class FuncType, where FuncType should have a method FuncType::ODEforcing
             *       to evaluate the right hand side \f$f(t,\boldsymbol{y})\f$ of the ODE.
             * \param y_0 the initial value \f$\boldsymbol{y}(t_0)\f$
             * \return An object of the class TimeIntegrationSolution which represents the vectors \f$\boldsymbol{y}^{[n]}\f$ and \f$t^{[n]}\f$
             *  that can be used to start the actual integration.
             */

            LUE TimeIntegrationSolutionSharedPtr InitializeData( const NekDouble                        deltaT,
                                                                       ConstDoubleArray               & y_0,
                                                                 const NekDouble                        time,
                                                                 const TimeIntegrationSchemeOperators & op );


            /**
             * \brief Explicit integration of an ODE.
             *
             * This function explicitely perfroms a single integration step of the ODE system:
             * \f[
             * \frac{d\boldsymbol{y}}{dt}=\boldsymbol{f}(t,\boldsymbol{y})
             * \f]
             *
             * \param deltaT The size of the timestep, i.e. \f$\Delta t\f$.
             * \param f an object of the class FuncType, where FuncType should have a method FuncType::ODEforcing
             *       to evaluate the right hand side \f$f(t,\boldsymbol{y})\f$ of the ODE.
             * \param y on input: the vectors \f$\boldsymbol{y}^{[n-1]}\f$ and \f$t^{[n-1]}\f$ (which corresponds to the 
             *    solution at the old time level)
             * \param y on output: the vectors \f$\boldsymbol{y}^{[n]}\f$ and \f$t^{[n]}\f$ (which corresponds to the 
             *    solution at the old new level)
             * \return The actual solution \f$\boldsymbol{y}^{n}\f$ at the new time level 
             *    (which in fact is also embedded in the argument y).
             */
            LUE ConstDoubleArray& TimeIntegrate( const NekDouble                          deltaT,
                                                       TimeIntegrationSolutionSharedPtr & y,
                                                 const TimeIntegrationSchemeOperators   & op );

            TimeIntegrationMethod      m_method;
            TimeIntegrationSchemeType  m_schemeType;

            unsigned int m_numMultiStepValues; // number of entries in input and output vector that correspond
                                               // to VALUES at previous time levels
            unsigned int m_numMultiStepDerivs; // number of entries in input and output vector that correspond
                                               // to DERIVATIVES at previous time levels
            Array<OneD,unsigned int> m_timeLevelOffset; // denotes to which time-level the entries in both 
                                                        // input and output vector correspond, e.g.
                                                        //     INPUT VECTOR --------> m_inputTimeLevelOffset
                                                        //    _            _               _ _
                                                        //   | u^n          |             | 0 | 
                                                        //   | u^{n-1}      |             | 1 | 
                                                        //   | u^{n-2}      |  ----->     | 2 | 
                                                        //   | dt f(u^{n-1})|             | 1 | 
                                                        //   | dt f(u^{n-2})|             | 2 | 
                                                        //    -            -               - -

            unsigned int m_numsteps;  //< Number of steps in this integration phase
            unsigned int m_numstages; //< Number of stages in multi-stage component. 

            bool m_firstStageEqualsOldSolution; //< Optimisation-flag 
            bool m_lastStageEqualsNewSolution;  //< Optimisation-flag

            Array<OneD, Array<TwoD,NekDouble> > m_A;
            Array<OneD, Array<TwoD,NekDouble> > m_B;
            Array<TwoD,NekDouble>               m_U;
            Array<TwoD,NekDouble>               m_V;

            bool m_initialised;   /// bool to identify if array has been initialised 
            int  m_nvar;          /// The number of variables in integration scheme. 
            int  m_npoints;       /// The size of inner data which is stored for reuse. 
            DoubleArray m_Y;      /// Array containing the stage values 
            DoubleArray m_tmp;    /// explicit right hand side of each stage equation

            TripleArray m_F;      /// Array corresponding to the stage Derivatives 
            TripleArray m_F_IMEX; /// Used to store the Explicit stage derivative of IMEX schemes

            NekDouble   m_T;     ///  Time at the different stages

          // template <typename> friend class Nektar::MemoryManager;

          // LUE friend TimeIntegrationSchemeManagerT & GetTimeIntegrationSchemeManager();

          // LUE static TimeIntegrationSchemeDataSharedPtr Create( const TimeIntegrationSchemeKey & key );

          // TimeIntegrationSchemeData( const TimeIntegrationSchemeKey & key );
            
            // These should never be called
            // TimeIntegrationSchemeData() : m_schemeKey( NullTimeIntegrationSchemeKey )
            // {
            //     NEKERROR(ErrorUtil::efatal,"Default Constructor for the TimeIntegrationScheme class should not be called");
            // }
            
            // TimeIntegrationSchemeData( const TimeIntegrationScheme & in ) : m_schemeKey( NullTimeIntegrationSchemeKey )
            // {
            //     NEKERROR(ErrorUtil::efatal,"Copy Constructor for the TimeIntegrationScheme class should not be called");
            // }

          static
            LUE bool VerifyIntegrationSchemeType(       TimeIntegrationSchemeType                    type,
                                                  const Array<OneD, const Array<TwoD, NekDouble> > & A,
                                                  const Array<OneD, const Array<TwoD, NekDouble> > & B,
                                                  const Array<TwoD, const NekDouble>               & U,
                                                  const Array<TwoD, const NekDouble>               & V );

            LUE void TimeIntegrate( const NekDouble                        deltaT,
                                          ConstTripleArray               & y_old,
                                          ConstSingleArray               & t_old,
                                          TripleArray                    & y_new,
                                          SingleArray                    & t_new,
                                    const TimeIntegrationSchemeOperators & op );


            inline int GetFirstDim(  ConstTripleArray & y ) const { return y[0].num_elements(); }
            inline int GetSecondDim( ConstTripleArray & y ) const { return y[0][0].num_elements(); }

            LUE bool CheckTimeIntegrateArguments( // const NekDouble timestep,
                                                         ConstTripleArray               & y_old,
                                                         ConstSingleArray               & t_old,
                                                         TripleArray                    & y_new,
                                                         SingleArray                    & t_new, 
                                                   const TimeIntegrationSchemeOperators & op ) const;

            LUE bool CheckIfFirstStageEqualsOldSolution( const Array<OneD, const Array<TwoD, NekDouble> >& A,
                                                         const Array<OneD, const Array<TwoD, NekDouble> >& B,
                                                         const Array<TwoD, const NekDouble>& U,
                                                         const Array<TwoD, const NekDouble>& V) const;

            LUE bool CheckIfLastStageEqualsNewSolution( const Array<OneD, const Array<TwoD, NekDouble> >& A,
                                                        const Array<OneD, const Array<TwoD, NekDouble> >& B,
                                                        const Array<TwoD, const NekDouble>& U,
                                                        const Array<TwoD, const NekDouble>& V) const;

        }; // end class TimeIntegrationSchemeData

        LUE std::ostream& operator<<(std::ostream& os, const TimeIntegrationSchemeData& rhs);
        LUE std::ostream& operator<<(std::ostream& os, const TimeIntegrationSchemeData::TimeIntegrationSchemeDataSharedPtr& rhs);
        
        // =========================================================================

    } // end of namespace LibUtilities
} // end of namespace Nektar

