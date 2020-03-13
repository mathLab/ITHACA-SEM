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
// Description: Header file of time integration scheme base class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_FOUNDATIONS_TIMEINTEGRATIONSCHEME_H
#define NEKTAR_LIB_UTILITIES_FOUNDATIONS_TIMEINTEGRATIONSCHEME_H

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/NekManager.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

namespace Nektar
{
    namespace LibUtilities
    {
        // Forward declaration of all classes in this file
        class TimeIntegrationScheme;
        class TimeIntegrationSolution;

        // typedefs
        typedef boost::shared_ptr<TimeIntegrationScheme>                TimeIntegrationSchemeSharedPtr;
        typedef std::vector<TimeIntegrationSchemeSharedPtr>             TimeIntegrationSchemeVector; 
        typedef std::vector<TimeIntegrationSchemeSharedPtr>::iterator   TimeIntegrationSchemeVectorIter; 
        typedef boost::shared_ptr<TimeIntegrationSolution>              TimeIntegrationSolutionSharedPtr;
        typedef std::vector<TimeIntegrationSolutionSharedPtr>           TimeIntegrationSolutionVector; 
        typedef std::vector<TimeIntegrationSolutionSharedPtr>::iterator TimeIntegrationSolutionVectorIter; 

        // =========================================================================
        // ==== ENUM LIST OF ALL SUPPORTED INTEGRATION SCHEMES
        // =========================================================================
        enum TimeIntegrationMethod
        {
            eNoTimeIntegrationMethod,
            eAdamsBashforthOrder1,            //!< Adams-Bashforth Forward multi-step scheme of order 1
            eAdamsBashforthOrder2,            //!< Adams-Bashforth Forward multi-step scheme of order 2
            eAdamsBashforthOrder3,            //!< Adams-Bashforth Forward multi-step scheme of order 3
            eAdamsMoultonOrder1,              //!< Adams-Moulton Forward multi-step scheme of order 1
            eAdamsMoultonOrder2,              //!< Adams-Moulton Forward multi-step scheme of order 2
            eBDFImplicitOrder1,               //!< BDF multi-step scheme of order 1 (implicit)
            eBDFImplicitOrder2,               //!< BDF multi-step scheme of order 2 (implicit)
            eClassicalRungeKutta4,            //!< Runge-Kutta multi-stage scheme 4th order explicit (old name)
            eRungeKutta4,                     //!< Classical RungeKutta4 method (new name for eClassicalRungeKutta4)
            eRungeKutta3_SSP,                 //!< Nonlinear SSP RungeKutta3 explicit
            eRungeKutta2_ImprovedEuler,       //!< Improved RungeKutta2 explicit (old name meaning Heun's method)
            eRungeKutta2_SSP,                 //!< Nonlinear SSP RungeKutta2 explicit (surrogate for eRungeKutta2_ImprovedEuler)
            eForwardEuler,                    //!< Forward Euler scheme
            eBackwardEuler,                   //!< Backward Euler scheme
            eIMEXOrder1,                      //!< IMEX 1st order scheme using Euler Backwards/Euler Forwards
            eIMEXOrder2,                      //!< IMEX 2nd order scheme using Backward Different Formula & Extrapolation
            eIMEXOrder3,                      //!< IMEX 3rd order scheme using Backward Different Formula & Extrapolation
            eMidpoint,                        //!< midpoint method (old name)
            eRungeKutta2,                     //!< Classical RungeKutta2 method (new name for eMidpoint)
            eDIRKOrder2,                      //!< Diagonally Implicit Runge Kutta scheme of order 3
            eDIRKOrder3,                      //!< Diagonally Implicit Runge Kutta scheme of order 3
            eCNAB,		                      //!< Crank-Nicolson/Adams-Bashforth Order 2 (CNAB)
            eIMEXGear,		                  //!< IMEX Gear Order 2
            eMCNAB,		                      //!< Modified Crank-Nicolson/Adams-Bashforth Order 2 (MCNAB)
            eIMEXdirk_1_1_1,		      	  //!< Forward-Backward Euler IMEX DIRK(1,1,1)
            eIMEXdirk_1_2_1,		      	  //!< Forward-Backward Euler IMEX DIRK(1,2,1)
            eIMEXdirk_1_2_2,		      	  //!< Implicit-Explicit Midpoint IMEX DIRK(1,2,2)
            eIMEXdirk_2_2_2,		          //!< L-stable, two stage, second order IMEX DIRK(2,2,2)
            eIMEXdirk_2_3_2,                  //!< L-stable, three stage, third order IMEX DIRK(3,4,3)
            eIMEXdirk_2_3_3,		      	  //!< L-stable, two stage, third order IMEX DIRK(2,3,3)
            eIMEXdirk_3_4_3,                  //!< L-stable, three stage, third order IMEX DIRK(3,4,3)
            eIMEXdirk_4_4_3,		      	  //!< L-stable, four stage, third order IMEX DIRK(4,4,3)
            SIZE_TimeIntegrationMethod        //!< Length of enum list
        };

        const char* const TimeIntegrationMethodMap[] = 
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
            
            typedef boost::function< void (InArrayType&, OutArrayType&, const NekDouble) >                  FunctorType1;
            typedef boost::function< void (InArrayType&, OutArrayType&, const NekDouble, const NekDouble) > FunctorType2;

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
                void DefineOdeRhs(FuncPointerT func, ObjectPointerT obj)
            {
                m_functors1[0] =  boost::bind(func, obj, _1, _2, _3);
            }

            template<typename FuncPointerT, typename ObjectPointerT> 
                void DefineOdeExplicitRhs(FuncPointerT func, ObjectPointerT obj)
            {
                m_functors1[1] =  boost::bind(func, obj, _1, _2, _3);
            }

            template<typename FuncPointerT, typename ObjectPointerT> 
                void DefineOdeImplicitRhs(FuncPointerT func, ObjectPointerT obj)
            {
                m_functors1[2] =  boost::bind(func, obj, _1, _2, _3);
            }

            template<typename FuncPointerT, typename ObjectPointerT> 
                void DefineProjection(FuncPointerT func, ObjectPointerT obj)
            {
                m_functors1[3] =  boost::bind(func, obj, _1, _2, _3);
            }

            template<typename FuncPointerT, typename ObjectPointerT> 
                void DefineImplicitSolve(FuncPointerT func, ObjectPointerT obj)
            {
                m_functors2[0] =  boost::bind(func, obj, _1, _2, _3, _4);
            }

            
            inline void DoOdeRhs(InArrayType     &inarray, 
                                 OutArrayType    &outarray, 
                                 const NekDouble time) const
            {
                ASSERTL1(!(m_functors1[0].empty()),"OdeRhs should be defined for this time integration scheme");
                m_functors1[0](inarray,outarray,time);
            }
            
            inline void DoOdeExplicitRhs(InArrayType     &inarray, 
                                         OutArrayType    &outarray, 
                                         const NekDouble time) const
            {
                ASSERTL1(!(m_functors1[1].empty()),"OdeExplicitRhs should be defined for this time integration scheme");
                m_functors1[1](inarray,outarray,time);
            }
            
            inline void DoOdeImplicitRhs(InArrayType     &inarray, 
                                        OutArrayType    &outarray, 
                                        const NekDouble time) const
            {
                ASSERTL1(!(m_functors1[2].empty()),"OdeImplictRhs should be defined for this time integration scheme");
                m_functors1[2](inarray,outarray,time);
            }

            inline void DoProjection(InArrayType     &inarray, 
                                     OutArrayType    &outarray, 
                                     const NekDouble time) const
            {
                ASSERTL1(!(m_functors1[3].empty()),"Projection operation should be defined for this time integration scheme");
                m_functors1[3](inarray,outarray,time);
            }
            
            inline void DoImplicitSolve(InArrayType     &inarray, 
                                        OutArrayType    &outarray, 
                                        const NekDouble time, 
                                        const NekDouble lambda) const
            {
                ASSERTL1(!(m_functors2[0].empty()),"ImplicitSolve should be defined for this time integration scheme");
                m_functors2[0](inarray,outarray,time,lambda);
            }


        protected:
            FunctorType1Array m_functors1;
            FunctorType2Array m_functors2;

        private:

        };

        // =========================================================================
        // ==== DEFINITION OF THE CLASS TimeIntegrationSchemeKey
        // =========================================================================
        class TimeIntegrationSchemeKey
        {
        public:
            
            struct opLess
            {
                LIB_UTILITIES_EXPORT bool operator()(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs) const;
            };
            
            TimeIntegrationSchemeKey(const TimeIntegrationMethod &method): 
            m_method(method)
            {
            }
            
            virtual ~TimeIntegrationSchemeKey()
            {
            }

            TimeIntegrationSchemeKey(const TimeIntegrationSchemeKey &key)
            {
                *this = key; // defer to assignment operator
            }

            TimeIntegrationSchemeKey& operator=(const TimeIntegrationSchemeKey &key)
            {
                m_method  = key.m_method;
                
                return *this;
            }

            inline TimeIntegrationMethod GetIntegrationMethod() const
            {
                return m_method;
            }

            inline bool operator==(const TimeIntegrationSchemeKey &key)
            {
                return (m_method == key.m_method);
            }

            inline bool operator== (const TimeIntegrationSchemeKey *y)
            {
                return (*this == *y);
            }

            inline bool operator != (const TimeIntegrationSchemeKey& y)
            {
                return (!(*this == y));
            }

            inline bool operator != (const TimeIntegrationSchemeKey *y)
            {
                return (!(*this == *y));
            }

            LIB_UTILITIES_EXPORT friend bool operator==(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs);
            LIB_UTILITIES_EXPORT friend bool operator<(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs);
            LIB_UTILITIES_EXPORT friend bool opLess::operator()(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs) const;

        protected:
            TimeIntegrationMethod m_method;  //!< integration method

        private:
            // This should never be called
            TimeIntegrationSchemeKey()
            {
                NEKERROR(ErrorUtil::efatal,"Default Constructor for the TimeIntegrationSchemeKey class should not be called");
            }
        };

        static const TimeIntegrationSchemeKey NullTimeIntegrationSchemeKey(eNoTimeIntegrationMethod);
        LIB_UTILITIES_EXPORT bool operator==(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs);
        LIB_UTILITIES_EXPORT bool operator<(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs);
        LIB_UTILITIES_EXPORT std::ostream& operator<<(std::ostream& os, const TimeIntegrationSchemeKey& rhs);

        // =====================================================================


        // =====================================================================
        // ==== DEFINITION OF THE FUNCTION TimeIntegrationSolution
        // =====================================================================
        // This is the interface the user can use to get hold of the different
        // time integration schemes. It returns you the NekManager 
        // singleton which manages all the time integration schemes...
        //
        typedef NekManager<TimeIntegrationSchemeKey, 
                           TimeIntegrationScheme, 
                           TimeIntegrationSchemeKey::opLess> TimeIntegrationSchemeManagerT;
        LIB_UTILITIES_EXPORT TimeIntegrationSchemeManagerT &TimeIntegrationSchemeManager(void);
        // =====================================================================

        // =====================================================================
        // ==== DEFINITION OF THE CLASS TimeIntegrationScheme
        // =====================================================================
        class TimeIntegrationScheme
        {
        public:

                typedef const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > > ConstTripleArray;
                typedef       Array<OneD,       Array<OneD, Array<OneD, NekDouble> > > TripleArray;
                typedef const Array<OneD, const Array<OneD, NekDouble> >               ConstDoubleArray;
                typedef       Array<OneD,       Array<OneD, NekDouble> >               DoubleArray;
                typedef const Array<OneD, const NekDouble >                            ConstSingleArray;
                typedef       Array<OneD,       NekDouble >                            SingleArray;
                typedef boost::function< void (ConstDoubleArray&, DoubleArray&, const NekDouble) >                  FunctorType1;
                typedef boost::function< void (ConstDoubleArray&, DoubleArray&, const NekDouble, const NekDouble) > FunctorType2;

        public:

            virtual ~TimeIntegrationScheme()
            {
            }

            inline const TimeIntegrationSchemeKey& GetIntegrationSchemeKey() const
            {
                return m_schemeKey;
            }

            inline TimeIntegrationMethod GetIntegrationMethod() const
            {
                return m_schemeKey.GetIntegrationMethod();
            }

            inline TimeIntegrationSchemeType GetIntegrationSchemeType() const
            {
                return m_schemeType;
            }

            inline NekDouble A(const unsigned int i, const unsigned int j) const
            {
                return m_A[0][i][j];
            }

            inline NekDouble B(const unsigned int i, const unsigned int j) const
            {
                return m_B[0][i][j];
            }

            inline NekDouble U(const unsigned int i, const unsigned int j) const
            {
                return m_U[i][j];
            }

            inline NekDouble V(const unsigned int i, const unsigned int j) const
            {
                return m_V[i][j];
            }

            inline NekDouble A_IMEX(const unsigned int i, const unsigned int j) const
            {
                return m_A[1][i][j];
            }

            inline NekDouble B_IMEX(const unsigned int i, const unsigned int j) const
            {
                return m_B[1][i][j];
            }

            inline unsigned int GetNsteps(void) const
            {
                return m_numsteps;
            }

            inline unsigned int GetNstages(void) const
            {
                return m_numstages;
            }

            inline unsigned int GetNmultiStepValues(void) const
            {
                return m_numMultiStepValues;
            }

            inline unsigned int GetNmultiStepDerivs(void) const
            {
                return m_numMultiStepDerivs;
            }

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
            LIB_UTILITIES_EXPORT TimeIntegrationSolutionSharedPtr InitializeScheme(const NekDouble timestep,
                                                                    ConstDoubleArray               &y_0    ,
                                                              const NekDouble                      time    ,
                                                              const TimeIntegrationSchemeOperators &op     );

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
            LIB_UTILITIES_EXPORT ConstDoubleArray& TimeIntegrate(const NekDouble   timestep,          
                                                  TimeIntegrationSolutionSharedPtr &y      ,
                                            const TimeIntegrationSchemeOperators   &op     );

            inline const Array<OneD, const unsigned int>& GetTimeLevelOffset()
            {
                return m_timeLevelOffset;
            }


        protected:
            TimeIntegrationSchemeKey  m_schemeKey; 
            TimeIntegrationSchemeType m_schemeType;
            unsigned int              m_numsteps;   //< Number of steps in multi-step component. 
            unsigned int              m_numstages;  //< Number of stages in multi-stage component. 

            bool m_firstStageEqualsOldSolution; //< Optimisation-flag 
            bool m_lastStageEqualsNewSolution;  //< Optimisation-flag

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

            Array<OneD, Array<TwoD,NekDouble> > m_A;
            Array<OneD, Array<TwoD,NekDouble> > m_B;
            Array<TwoD,NekDouble>               m_U;
            Array<TwoD,NekDouble>               m_V;

        private: 
            bool m_initialised;   /// bool to identify if array has been initialised 
            int  m_nvar;          /// The number of variables in integration scheme. 
            int  m_npoints;       /// The size of inner data which is stored for reuse. 
            DoubleArray m_Y;      /// Array containing the stage values 
            DoubleArray m_tmp;    /// explicit right hand side of each stage equation

            TripleArray m_F;      /// Array corresponding to the stage Derivatives 
            TripleArray m_F_IMEX; /// Used to store the Explicit stage derivative of IMEX schemes

            NekDouble   m_T;     ///  Time at the different stages

            template <typename> friend class Nektar::MemoryManager;
            LIB_UTILITIES_EXPORT friend TimeIntegrationSchemeManagerT &TimeIntegrationSchemeManager(void);

            LIB_UTILITIES_EXPORT static boost::shared_ptr<TimeIntegrationScheme> Create(const TimeIntegrationSchemeKey &key);

            TimeIntegrationScheme(const TimeIntegrationSchemeKey &key);
            
            // These should never be called
            TimeIntegrationScheme():m_schemeKey(NullTimeIntegrationSchemeKey)
            {
                NEKERROR(ErrorUtil::efatal,"Default Constructor for the TimeIntegrationScheme class should not be called");
            }
            
            TimeIntegrationScheme(const TimeIntegrationScheme &in):m_schemeKey(NullTimeIntegrationSchemeKey)
            {
                NEKERROR(ErrorUtil::efatal,"Copy Constructor for the TimeIntegrationScheme class should not be called");
            }

            LIB_UTILITIES_EXPORT bool VerifyIntegrationSchemeType(TimeIntegrationSchemeType type,
                                             const Array<OneD, const Array<TwoD, NekDouble> >& A,
                                             const Array<OneD, const Array<TwoD, NekDouble> >& B,
                                             const Array<TwoD, const NekDouble>& U,
                                             const Array<TwoD, const NekDouble>& V) const;

            LIB_UTILITIES_EXPORT void TimeIntegrate(const NekDouble timestep,
                                     ConstTripleArray               &y_old  ,
                                     ConstSingleArray               &t_old  ,
                                     TripleArray                    &y_new  ,
                                     SingleArray                    &t_new  ,
                               const TimeIntegrationSchemeOperators &op     );


            inline int GetFirstDim(ConstTripleArray &y) const
            {
                return y[0].num_elements();
            }

            inline int GetSecondDim(ConstTripleArray &y) const
            {
                return y[0][0].num_elements();
            }

            LIB_UTILITIES_EXPORT bool CheckTimeIntegrateArguments(const NekDouble timestep,
                                                   ConstTripleArray               &y_old  ,
                                                   ConstSingleArray               &t_old  ,
                                                   TripleArray                    &y_new  ,
                                                   SingleArray                    &t_new  ,
                                             const TimeIntegrationSchemeOperators &op) const;

            LIB_UTILITIES_EXPORT bool CheckIfFirstStageEqualsOldSolution(const Array<OneD, const Array<TwoD, NekDouble> >& A,
                                                    const Array<OneD, const Array<TwoD, NekDouble> >& B,
                                                    const Array<TwoD, const NekDouble>& U,
                                                    const Array<TwoD, const NekDouble>& V) const;
            
            LIB_UTILITIES_EXPORT bool CheckIfLastStageEqualsNewSolution(const Array<OneD, const Array<TwoD, NekDouble> >& A,
                                                   const Array<OneD, const Array<TwoD, NekDouble> >& B,
                                                   const Array<TwoD, const NekDouble>& U,
                                                   const Array<TwoD, const NekDouble>& V) const;



        };


        // =========================================================================
        // ==== DEFINITION OF THE CLASS TimeIntegrationSolution
        // =========================================================================
        class TimeIntegrationSolution
        {
        public:
            typedef Array<OneD, Array<OneD, Array<OneD, NekDouble> > > TripleArray;
            typedef Array<OneD, Array<OneD, NekDouble> >               DoubleArray;
            
            // Constructor for single step methods
            LIB_UTILITIES_EXPORT TimeIntegrationSolution(const TimeIntegrationSchemeKey &key, 
                                    const DoubleArray& y, 
                                    const NekDouble time, 
                                    const NekDouble timestep);
            // Constructor for multi-step methods
            LIB_UTILITIES_EXPORT TimeIntegrationSolution(const TimeIntegrationSchemeKey &key, 
                                    const TripleArray& y, 
                                    const Array<OneD, NekDouble>& t);
            LIB_UTILITIES_EXPORT TimeIntegrationSolution(const TimeIntegrationSchemeKey &key,
                                    unsigned int nvar,
                                    unsigned int npoints);
            LIB_UTILITIES_EXPORT TimeIntegrationSolution(const TimeIntegrationSchemeKey &key);

            inline const TimeIntegrationSchemeSharedPtr& GetIntegrationScheme() const
            {
                return m_scheme;
            }

            inline const TimeIntegrationSchemeKey& GetIntegrationSchemeKey() const
            {
                return m_scheme->GetIntegrationSchemeKey();
            }
            
            inline TimeIntegrationMethod GetIntegrationMethod() const
            {
                return m_scheme->GetIntegrationMethod();
            }

            inline const TripleArray& GetSolutionVector() const
            {
                return m_solVector;
            }

            inline TripleArray& UpdateSolutionVector()
            {
                return m_solVector;
            }

            inline const DoubleArray& GetSolution() const
            {
                return m_solVector[0];
            }

            inline DoubleArray& UpdateSolution()
            {
                return m_solVector[0];
            }

            inline const Array<OneD, const NekDouble>& GetTimeVector() const
            {
                return m_t;
            }

            inline Array<OneD, NekDouble>& UpdateTimeVector()
            {
                return m_t;
            }

            inline NekDouble GetTime() const
            {
                return m_t[0];
            }

            inline int GetNsteps()
            {
                return m_scheme->GetNsteps();
            }

            inline int GetFirstDim() const
            {
                return m_solVector[0].num_elements();
            }

            inline int GetSecondDim() const
            {
                return m_solVector[0][0].num_elements();
            }

            // Return the number of entries in the solution vector that correspond
            // to (multi-step) values
            inline unsigned int GetNvalues(void) const
            {
                return m_scheme->GetNmultiStepValues();
            }

            // Return the number of entries in the solution vector that correspond
            // to (multi-step) derivatives
            inline unsigned int GetNderivs(void) const
            {
                return m_scheme->GetNmultiStepDerivs();
            }

            // Returns an array which indicates to which time-level the entries in the
            // solution vector correspond
            inline const Array<OneD, const unsigned int>& GetTimeLevelOffset()
            {
                return m_scheme->GetTimeLevelOffset();
            }

            // returns the entry in the solution vector which corresponds to 
            // the (multi-step) value at the time-level with specified offset
            inline DoubleArray& GetValue(const unsigned int timeLevelOffset)
            {
                int nMultiStepVals = m_scheme->GetNmultiStepValues();
                const Array<OneD, const unsigned int>& offsetvec = GetTimeLevelOffset();

                for(int i = 0; i < nMultiStepVals; i++)
                {
                    if( timeLevelOffset == offsetvec[i] )
                    {
                        return m_solVector[i];
                    }
                }
                ASSERTL1(false,"The solution vector of this scheme does not contain a value at the requested time-level");
                return m_solVector[0];
            }

            // returns the entry in the solution vector which corresponds to 
            // the (multi-step) derivative at the time-level with specified offset
            inline DoubleArray& GetDerivative(const unsigned int timeLevelOffset)
            {
                int nMultiStepVals = m_scheme->GetNmultiStepValues();
                int size           = m_scheme->GetNsteps();
                const Array<OneD, const unsigned int>& offsetvec = GetTimeLevelOffset();

                for(int i = nMultiStepVals; i < size; i++)
                {
                    if( timeLevelOffset == offsetvec[i] )
                    {
                        return m_solVector[i];
                    }
                }
                ASSERTL1(false,"The solution vector of this scheme does not contain a derivative at the requested time-level");
                return m_solVector[0];
            }

            // returns the time associated with the (multi-step) value at the time-level with the 
            // given offset
            inline NekDouble GetValueTime(const unsigned int timeLevelOffset)
            {
                int nMultiStepVals = m_scheme->GetNmultiStepValues();
                const Array<OneD, const unsigned int>& offsetvec = GetTimeLevelOffset();

                for(int i = 0; i < nMultiStepVals; i++)
                {
                    if( timeLevelOffset == offsetvec[i] )
                    {
                        return m_t[i];
                    }
                }
                ASSERTL1(false,"The solution vector of this scheme does not contain a value at the requested time-level");
                return m_t[0];
            }

            // sets the (multi-step) value and time in the solution
            // vector which corresponds to
            // the value at the time-level with specified offset
            inline void SetValue(const unsigned int timeLevelOffset, const DoubleArray& y, const NekDouble t)
            {
                int nMultiStepVals = m_scheme->GetNmultiStepValues();
                const Array<OneD, const unsigned int>& offsetvec = GetTimeLevelOffset();

                for(int i = 0; i < nMultiStepVals; i++)
                {
                    if( timeLevelOffset == offsetvec[i] )
                    {
                        m_solVector[i] = y;
                        m_t[i] = t;
                        return;
                    }
                }
            }

            // sets the (multi-step) derivative and time in the
            // solution vector which corresponds to
            // the derivative at the time-level with specified offset
            inline void SetDerivative(const unsigned int timeLevelOffset, const DoubleArray& y, const NekDouble timestep)
            {
                int nMultiStepVals = m_scheme->GetNmultiStepValues();
                int size           = m_scheme->GetNsteps();
                const Array<OneD, const unsigned int>& offsetvec = GetTimeLevelOffset();

                for(int i = nMultiStepVals; i < size; i++)
                {
                    if( timeLevelOffset == offsetvec[i] )
                    {
                        m_solVector[i] = y;
                        m_t[i] = timestep;
                        return;
                    }
                }
            }

            // sets the soln Vector 
            inline void SetSolVector(const  int Offset, const DoubleArray& y)
            {
                m_solVector[Offset] = y;
            }

            // Rotate the solution vector 
            // (i.e. updating without calculating/inserting new values)
            inline void RotateSolutionVector()
            {
                int i;
                int nMultiStepVals = m_scheme->GetNmultiStepValues();
                int size           = m_scheme->GetNsteps();
                for(i = (nMultiStepVals-1); i > 0; i--)
                {
                    m_solVector[i] = m_solVector[i-1];
                }

                for(i = (size-1); i > nMultiStepVals; i--)
                {
                    m_solVector[i] = m_solVector[i-1];
                }
            }


        private:
            TimeIntegrationSchemeSharedPtr m_scheme;
            TripleArray m_solVector;
            Array<OneD,NekDouble> m_t;
        };
        // =========================================================================


        LIB_UTILITIES_EXPORT std::ostream& operator<<(std::ostream& os, const TimeIntegrationScheme& rhs);
        LIB_UTILITIES_EXPORT std::ostream& operator<<(std::ostream& os, const TimeIntegrationSchemeSharedPtr& rhs);
        
        // =========================================================================

    }; // end of namespace
} // end of namespace 

#endif //NEKTAR_LIB_UTILITIES_FOUNDATIONS_TIMEINTEGRATIONSCHEME_H
