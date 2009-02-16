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

#include <math.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/NekManager.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicConst/NektarUnivConsts.hpp>


namespace Nektar
{
    namespace LibUtilities
    {
        // Forward declaration of all classes in this file
        class TimeIntegrationSchemeKey;
        class TimeIntegrationScheme;
        class TimeIntegrationSolution;
        class TimeIntegrationSchemeOperators;

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
            eAdamsMoultonOrder1,              //!< Adams-Moulton Forward multi-step scheme of order 1
            eAdamsMoultonOrder2,              //!< Adams-Moulton Forward multi-step scheme of order 2
            eClassicalRungeKutta4,            //!< Runge-Kutta multi-stage scheme
            eForwardEuler,                    //!< Forward euler scheme
            eBackwardEuler,                   //!< Backward euler scheme
            eMidpoint,                        //!< midpoint method
            eDIRKOrder2,                      //!< Diagonally Implicit Runge Kutta scheme of order 3
            eDIRKOrder3,                      //!< Diagonally Implicit Runge Kutta scheme of order 3
            eIMEXdirk_3_4_3,                  //!< L-stable, three stage, third order IMEX DIRK(3,4,3)
            SIZE_TimeIntegrationMethod        //!< Length of enum list
        };

        const char* const TimeIntegrationMethodMap[] = 
        {
            "NoTimeIntegrationMethod",
            "AdamsBashforthOrder1",            
            "AdamsBashforthOrder2",          
            "AdamsMoultonOrder1",           
            "AdamsMoultonOrder2",     
            "ClassicalRungeKutta4",           
            "ForwardEuler",                    
            "BackwardEuler",                  
            "Midpoint",                       
            "DIRKOrder2",                     
            "DIRKOrder3",                     
            "IMEXdirk_3_4_3"                
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

        // =========================================================================



        // =========================================================================
        // ==== DEFINITION OF THE CLASS TimeIntegrationSchemeOperators
        // =========================================================================
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
            m_functors1(5),
            m_functors2(1)
            {
            }

            template<typename FuncPointerT, typename ObjectPointerT> 
                void DefineOdeLhs(FuncPointerT func, ObjectPointerT obj)
            {
                m_functors1[0] =  boost::bind(func, obj, _1, _2, _3);
            }

            template<typename FuncPointerT, typename ObjectPointerT> 
                void DefineOdeLhsSolve(FuncPointerT func, ObjectPointerT obj)
            {
                m_functors1[1] =  boost::bind(func, obj, _1, _2, _3);
            }

            template<typename FuncPointerT, typename ObjectPointerT> 
                void DefineOdeRhs(FuncPointerT func, ObjectPointerT obj)
            {
                m_functors1[2] =  boost::bind(func, obj, _1, _2, _3);
            }

            template<typename FuncPointerT, typename ObjectPointerT> 
                void DefineOdeExplicitRhs(FuncPointerT func, ObjectPointerT obj)
            {
                m_functors1[3] =  boost::bind(func, obj, _1, _2, _3);
            }

            template<typename FuncPointerT, typename ObjectPointerT> 
                void DefineOdeImplicitRhs(FuncPointerT func, ObjectPointerT obj)
            {
                m_functors1[4] =  boost::bind(func, obj, _1, _2, _3);
            }

            template<typename FuncPointerT, typename ObjectPointerT> 
                void DefineImplicitSolve(FuncPointerT func, ObjectPointerT obj)
            {
                m_functors2[0] =  boost::bind(func, obj, _1, _2, _3, _4);
            }

            inline void DoOdeLhs(InArrayType     &inarray, 
                                 OutArrayType    &outarray, 
                                 const NekDouble time) const
            {
                ASSERTL1(!(m_functors1[0].empty()),"OdeLhs should be defined for this time integration scheme");
                m_functors1[0](inarray,outarray,time);
            }

            inline void DoOdeLhsSolve(InArrayType     &inarray, 
                                      OutArrayType    &outarray, 
                                      const NekDouble time) const
            {
                ASSERTL1(!(m_functors1[1].empty()),"OdeLhsSolve should be defined for this time integration scheme");
                m_functors1[1](inarray,outarray,time);
            }
            
            inline void DoOdeRhs(InArrayType     &inarray, 
                                 OutArrayType    &outarray, 
                                 const NekDouble time) const
            {
                ASSERTL1(!(m_functors1[2].empty()),"OdeRhs should be defined for this time integration scheme");
                m_functors1[2](inarray,outarray,time);
            }
            
            inline void DoOdeExplicitRhs(InArrayType     &inarray, 
                                         OutArrayType    &outarray, 
                                         const NekDouble time) const
            {
                ASSERTL1(!(m_functors1[3].empty()),"OdeExplicitRhs should be defined for this time integration scheme");
                m_functors1[3](inarray,outarray,time);
            }
            
            inline void DoOdeImplicitRhs(InArrayType     &inarray, 
                                        OutArrayType    &outarray, 
                                        const NekDouble time) const
            {
                ASSERTL1(!(m_functors1[4].empty()),"OdeImplictRhs should be defined for this time integration scheme");
                m_functors1[4](inarray,outarray,time);
            }
            
            inline void DoImplicitSolve(InArrayType     &inarray, 
                                        OutArrayType    &outarray, 
                                        const NekDouble lambda,
                                        const NekDouble time) const
            {
                ASSERTL1(!(m_functors2[0].empty()),"ImplicitSolve should be defined for this time integration scheme");
                m_functors2[0](inarray,outarray,lambda,time);
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
                bool operator()(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs) const;
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

            friend bool operator==(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs);
            friend bool operator<(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs);
            friend bool opLess::operator()(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs) const;

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
        bool operator==(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs);
        bool operator<(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs);
        std::ostream& operator<<(std::ostream& os, const TimeIntegrationSchemeKey& rhs);

        // =========================================================================


        // =========================================================================
        // ==== DEFINITION OF THE FUNCTION TimeIntegrationSolution
        // =========================================================================
        // This is the interface the user can use to get hold of the different
        // time integration schemes. It returns you the NekManager singleton which
        // manages all the time integration schemes...
        typedef NekManager<TimeIntegrationSchemeKey, 
                           TimeIntegrationScheme, 
                           TimeIntegrationSchemeKey::opLess> TimeIntegrationSchemeManagerT;
        TimeIntegrationSchemeManagerT &TimeIntegrationSchemeManager(void);
        // =========================================================================

        // =========================================================================
        // ==== DEFINITION OF THE CLASS TimeIntegrationSolution
        // =========================================================================
        class TimeIntegrationSolution
        {
        public:
            typedef Array<OneD, Array<OneD, Array<OneD, NekDouble> > > TripleArray;
            typedef Array<OneD, Array<OneD, NekDouble> >               DoubleArray;
            
            // Constructor for single step methods
            TimeIntegrationSolution(TimeIntegrationMethod method, 
                                    const DoubleArray& y, 
                                    const DoubleArray& My, 
                                    NekDouble t);
            // Constructor for multi-step methods
            TimeIntegrationSolution(TimeIntegrationMethod method, 
                                    const DoubleArray& y, 
                                    const TripleArray& My, 
                                    const Array<OneD, NekDouble>& t);
            TimeIntegrationSolution(TimeIntegrationMethod method,
                                    unsigned int nsteps,
                                    unsigned int nvar,
                                    unsigned int npoints);
            
            inline TimeIntegrationMethod GetIntegrationMethod() const
            {
                return m_method;
            }

            inline const TripleArray& GetSolutionVector()
            {
                return m_solVector;
            }

            inline TripleArray& UpdateSolutionVector()
            {
                return m_solVector;
            }

            inline const DoubleArray& GetSolution()
            {
                return m_sol;
            }

            inline DoubleArray& UpdateSolution()
            {
                return m_sol;
            }

            inline const Array<OneD, const NekDouble>& GetTimeVector()
            {
                return m_t;
            }

            inline Array<OneD, NekDouble>& UpdateTimeVector()
            {
                return m_t;
            }

            inline NekDouble GetTime()
            {
                return m_t[0];
            }

            inline int GetFirstDim() const
            {
                return m_solVector[0].num_elements();
            }

            inline int GetSecondDim() const
            {
                return m_solVector[0][0].num_elements();
            }

        private:
            TimeIntegrationMethod m_method;
            DoubleArray m_sol;
            TripleArray m_solVector;
            Array<OneD,NekDouble> m_t;
        };
        // =========================================================================


        // =========================================================================
        // ==== DEFINITION OF THE CLASS TimeIntegrationScheme
        // =========================================================================
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

            /**
             * \brief This function initialises the time integration scheme
             *
             * Given the solution at the initial time level \f$\boldsymbol{y}(t_0)\f$, this function
             * generates the vectors \f$\boldsymbol{y}^{[n]}\f$ and \f$t^{[n]}\f$
             * needed to evaluate the time integration scheme formulated 
             * as a General Linear Method. These vectors are embedded in an object of the
             * class TimeIntegrationSolution. This class is the abstraction of the
             * input and output vectors of the General Linear Method.
             *
             * For single-step methods, this function is trivial as it just wraps a
             * TimeIntegrationSolution object around the given initial values and initial time.
             * However, for multistep methods, actual time stepping is being done
             * to evaluate the necessary parameters at multiple time levels needed to start
             * the actual integration.
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
            TimeIntegrationSolutionSharedPtr InitializeScheme(const NekDouble                      timestep,
                                                                    NekDouble                      &time   ,
                                                                    int                            &nsteps ,
                                                                    ConstDoubleArray               &y_0    ,
                                                              const TimeIntegrationSchemeOperators &op     ) const;

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
            ConstDoubleArray& TimeIntegrate(const NekDouble                        timestep,          
                                                  TimeIntegrationSolutionSharedPtr &y      ,
                                            const TimeIntegrationSchemeOperators   &op     ) const;


        protected:
            TimeIntegrationSchemeKey  m_schemeKey; 
            TimeIntegrationSchemeType m_schemeType;
            unsigned int              m_numsteps;   //< Number of steps in multi-step component. 
            unsigned int              m_numstages;  //< Number of stages in multi-stage component. 

            bool m_firstStageEqualsOldSolution;  //< Optimisation-flag 
            bool m_lastStageEqualsNewSolution;   //< Optimisation-flag

            Array<OneD, Array<TwoD,NekDouble> >   m_A;
            Array<OneD, Array<TwoD,NekDouble> >   m_B;
            Array<TwoD,NekDouble>    m_U;
            Array<TwoD,NekDouble>    m_V;

        private: 
            
            template <typename> friend class Nektar::MemoryManager;
            friend TimeIntegrationSchemeManagerT &TimeIntegrationSchemeManager(void);

            static boost::shared_ptr<TimeIntegrationScheme> Create(const TimeIntegrationSchemeKey &key);

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

            bool VerifyIntegrationSchemeType(TimeIntegrationSchemeType type,
                                             const Array<OneD, const Array<TwoD, NekDouble> >& A,
                                             const Array<OneD, const Array<TwoD, NekDouble> >& B,
                                             const Array<TwoD, const NekDouble>& U,
                                             const Array<TwoD, const NekDouble>& V) const;

            void TimeIntegrate(const NekDouble                      timestep,      
                                     ConstTripleArray               &y_old  ,
                                     ConstDoubleArray               &sol_old,
                                     ConstSingleArray               &t_old  ,
                                     TripleArray                    &y_new  ,
                                     DoubleArray                    &sol_new,
                                     SingleArray                    &t_new  ,
                               const TimeIntegrationSchemeOperators &op     ) const;


            inline int GetFirstDim(ConstTripleArray &y) const
            {
                return y[0].num_elements();
            }

            inline int GetSecondDim(ConstTripleArray &y) const
            {
                return y[0][0].num_elements();
            }

            bool CheckTimeIntegrateArguments(const NekDouble                      timestep,      
                                                   ConstTripleArray               &y_old  ,
                                                   ConstDoubleArray               &sol_old,
                                                   ConstSingleArray               &t_old  ,
                                                   TripleArray                    &y_new  ,
                                                   DoubleArray                    &sol_new,
                                                   SingleArray                    &t_new  ,
                                             const TimeIntegrationSchemeOperators &op) const;

            bool CheckIfFirstStageEqualsOldSolution(const Array<OneD, const Array<TwoD, NekDouble> >& A,
                                                    const Array<OneD, const Array<TwoD, NekDouble> >& B,
                                                    const Array<TwoD, const NekDouble>& U,
                                                    const Array<TwoD, const NekDouble>& V) const;
            
            bool CheckIfLastStageEqualsNewSolution(const Array<OneD, const Array<TwoD, NekDouble> >& A,
                                                   const Array<OneD, const Array<TwoD, NekDouble> >& B,
                                                   const Array<TwoD, const NekDouble>& U,
                                                   const Array<TwoD, const NekDouble>& V) const;



        };

        std::ostream& operator<<(std::ostream& os, const TimeIntegrationScheme& rhs);
        std::ostream& operator<<(std::ostream& os, const TimeIntegrationSchemeSharedPtr& rhs);
        
        // =========================================================================

    }; // end of namespace
} // end of namespace 

#endif //NEKTAR_LIB_UTILITIES_FOUNDATIONS_TIMEINTEGRATIONSCHEME_H
