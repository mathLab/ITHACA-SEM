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


namespace Nektar
{
    namespace LibUtilities
    {
        // Forward declaration of all classes in this file
        class TimeIntegrationSchemeKey;
        class TimeIntegrationScheme;
        class TimeIntegrationSolution;

        // typedefs
        typedef boost::shared_ptr<TimeIntegrationScheme> TimeIntegrationSchemeSharedPtr;
        typedef std::vector< TimeIntegrationSchemeSharedPtr > TimeIntegrationSchemeVector; 
        typedef std::vector< TimeIntegrationSchemeSharedPtr >::iterator TimeIntegrationSchemeVectorIter; 
        typedef boost::shared_ptr<TimeIntegrationSolution> TimeIntegrationSolutionSharedPtr;
        typedef std::vector< TimeIntegrationSolutionSharedPtr > TimeIntegrationSolutionVector; 
        typedef std::vector< TimeIntegrationSolutionSharedPtr >::iterator TimeIntegrationSolutionVectorIter; 

        // =========================================================================
        // ==== ENUM LIST OF ALL SUPPORTED INTEGRATION SCHEMES
        // =========================================================================
        enum TimeIntegrationType
        {
            eNoTimeIntegrationType,
            eAdamsBashforthOrder1,            //!< Adams-Bashforth Forward multi-step scheme of order 1
            eAdamsBashforthOrder2,            //!< Adams-Bashforth Forward multi-step scheme of order 2
            eAdamsMoultonOrder1,              //!< Adams-Moulton Forward multi-step scheme of order 1
            eAdamsMoultonOrder2,              //!< Adams-Moulton Forward multi-step scheme of order 2
            eClassicalRungeKutta4,            //!< Runge-Kutta multi-stage scheme
            eForwardEuler,                    //!< Forward euler scheme
            eBackwardEuler,                   //!< Backward euler scheme
            eMidpoint,                        //!< midpoint method
            SIZE_TimeIntegrationType          //!< Length of enum list
        };

        const char* const TimeIntegrationTypeMap[] = 
        {
            "NoTimeIntegrationType",
            "1st order Adams-Bashforth",
            "2nd order Adams-Bashforth",
            "1st order Adams-Moulton",
            "2nd order Adams-Moulton",
            "Classical Runge-Kutta 4",
            "Forward Euler",
            "Backward Euler"
            "Midpoint method"
        };
        // =========================================================================


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
            
            TimeIntegrationSchemeKey(const TimeIntegrationType &integrationtype): 
            m_integrationtype(integrationtype)
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
                m_integrationtype  = key.m_integrationtype;
                
                return *this;
            }

            inline TimeIntegrationType GetIntegrationSchemeType() const
            {
                return m_integrationtype;
            }

            inline bool operator==(const TimeIntegrationSchemeKey &key)
            {
                return (m_integrationtype == key.m_integrationtype);
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
            TimeIntegrationType m_integrationtype;  //!< Type of the integration scheme

        private:
            // This should never be called
            TimeIntegrationSchemeKey()
            {
                NEKERROR(ErrorUtil::efatal,"Default Constructor for the TimeIntegrationSchemeKey class should not be called");
            }
        };

        static const TimeIntegrationSchemeKey NullTimeIntegrationSchemeKey(eNoTimeIntegrationType);
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
        typedef NekManager<TimeIntegrationSchemeKey, TimeIntegrationScheme, TimeIntegrationSchemeKey::opLess> TimeIntegrationSchemeManagerT;
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
            TimeIntegrationSolution(TimeIntegrationType schemeType, const DoubleArray& y, NekDouble t);
            TimeIntegrationSolution(TimeIntegrationType schemeType, const TripleArray& y, const Array<OneD, NekDouble>& t);
            TimeIntegrationSolution(TimeIntegrationType schemeType, 
                                    unsigned int nsteps,
                                    unsigned int nvar,
                                    unsigned int npoints);
            
            inline TimeIntegrationType GetIntegrationSchemeType() const
            {
                return m_schemeType;
            }

            inline const TripleArray& GetSolutionVector()
            {
                return m_y;
            }

            inline TripleArray& UpdateSolutionVector()
            {
                return m_y;
            }

            inline const DoubleArray& GetSolution()
            {
                return m_y[0];
            }

            inline DoubleArray& UpdateSolution()
            {
                return m_y[0];
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
    
            inline void SetSolutionVector(const TripleArray& y_new)
            {
                ASSERTL1(m_y.num_elements()==y_new.num_elements(),"Dimensions do not match");
                ASSERTL1(m_y[0].num_elements()==y_new[0].num_elements(),"Dimensions do not match");
                ASSERTL1(m_y[0][0].num_elements()==y_new[0][0].num_elements(),"Dimensions do not match");

                m_y = y_new;
            }

        private:
            TimeIntegrationType m_schemeType;
            TripleArray m_y;
            Array<OneD,NekDouble> m_t;
        };
        // =========================================================================


        // =========================================================================
        // ==== DEFINITION OF THE CLASS TimeIntegrationScheme
        // =========================================================================
        class TimeIntegrationScheme
        {
        public:

                typedef Array<OneD, Array<OneD, Array<OneD, NekDouble> > > TripleArray;
                typedef Array<OneD, Array<OneD, NekDouble> >               DoubleArray;

        public:

            static boost::shared_ptr<TimeIntegrationScheme> Create(const TimeIntegrationSchemeKey &key);

            virtual ~TimeIntegrationScheme()
            {
            }

            inline TimeIntegrationType GetIntegrationSchemeType() const
            {
                return m_schemeKey.GetIntegrationSchemeType();
            }

            inline const Array<TwoD, const NekDouble>& GetA(void) const
            {
                return m_A;
            }  

            inline const Array<TwoD, const NekDouble>& GetB(void) const
            {
                return m_B;
            }       

            inline const Array<TwoD, const NekDouble>& GetU(void) const
            {
                return m_U;
            }       

            inline const Array<TwoD, const NekDouble>& GetV(void) const
            {
                return m_V;
            }

            inline unsigned int GetNsteps(void) const
            {
                return m_numsteps;
            }

            inline unsigned int GetNstages(void) const
            {
                return m_numstages;
            }

            bool IsExplicitMethod() const;

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
            template<typename FuncType>
                TimeIntegrationSolutionSharedPtr InitializeScheme(NekDouble timestep,NekDouble& time,int& nsteps,FuncType& f,const DoubleArray& y_0)
            {
                TimeIntegrationSolutionSharedPtr y_out;
                TimeIntegrationType type = m_schemeKey.GetIntegrationSchemeType();
                
                switch(type)
                {
                case eForwardEuler:
                case eBackwardEuler:
                case eAdamsBashforthOrder1:
                case eAdamsMoultonOrder1:
                case eMidpoint:
                case eClassicalRungeKutta4:
                    {
                        nsteps = 0;
                        y_out = MemoryManager<TimeIntegrationSolution>::AllocateSharedPtr(type,y_0,time);
                    }
                break;
                case eAdamsBashforthOrder2:
                    {
                        // To initialise, we do a 1st order Adams Bashforth step (=Forward euler) to calculate the
                        // solution at the first time level
                        TimeIntegrationSchemeKey       IntKey(eAdamsBashforthOrder1);
                        TimeIntegrationSchemeSharedPtr IntScheme = TimeIntegrationSchemeManager()[IntKey];
                        TimeIntegrationSolutionSharedPtr u = IntScheme->InitializeScheme(timestep,time,nsteps,f,y_0);
                        IntScheme->ExplicitIntegration(timestep,f,u);
                        
                        int i;
                        int nvar = y_0.num_elements();
                        int npoints = y_0[0].num_elements();
                        DoubleArray f_y_0(nvar);
                        for(i = 0; i < nvar; i++)
                        {
                            f_y_0[i] = Array<OneD,NekDouble>(npoints);
                        }
                        f.ODEforcing(y_0,f_y_0,time);

                        for(i = 0; i < nvar; i++)
                        {
                            Blas::Dscal(npoints,timestep,f_y_0[i].get(),1);
                        }

                        TripleArray y(2);
                        y[0] = u->GetSolution();
                        y[1] = f_y_0;

                        Array<OneD,NekDouble> t(2);
                        t[0] = u->GetTime();
                        t[1] = timestep;

                        time = u->GetTime();
                        nsteps = 1;
                        y_out = MemoryManager<TimeIntegrationSolution>::AllocateSharedPtr(type,y,t);
                    }
                    break;
                default:
                    {
                        NEKERROR(ErrorUtil::efatal,"Methods need implementation for specified integration scheme type.");
                    }
                }
                
            return y_out;
            }


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
            // Explicitly integrate for one time step the system 
            //              du/dt = f(u)
            // where we pass a class "InClass" which must have the methods
            //    InClass.Forcing(const Array<OneD, Array<OneD, NekDouble> > u) to evaluate f(0)    
            template<typename FuncType>
                const DoubleArray& ExplicitIntegration(NekDouble timestep,                                     
                                                      FuncType& f,
                                                      TimeIntegrationSolutionSharedPtr& y)
            {
                ASSERTL1(IsExplicitMethod()==true,"Implicit integration scheme cannot be used for explit integration");

                TimeIntegrationType type = y->GetIntegrationSchemeType();
                int nvar    = (y->GetSolutionVector())[0].num_elements();
                int npoints = (y->GetSolutionVector())[0][0].num_elements();
                
                ASSERTL1(type == m_schemeKey.GetIntegrationSchemeType(),
                         "Input and output argument are constructed for a different type of time integration scheme.");
                
                TimeIntegrationSolutionSharedPtr y_new = MemoryManager<TimeIntegrationSolution>::AllocateSharedPtr(type,m_numsteps,nvar,npoints);  
                
                ExplicitIntegration(timestep,f,
                                    y->GetSolutionVector(),
                                    y->GetTimeVector(),
                                    y_new->UpdateSolutionVector(),
                                    y_new->UpdateTimeVector()); 

                y = y_new;

                return y->GetSolution();
            }

        protected:
            TimeIntegrationSchemeKey m_schemeKey; 
            unsigned int             m_numsteps;  //< Number of steps in multi-step component. 
            unsigned int             m_numstages; //< Number of stages in multi-stage component. 

            Array<TwoD,NekDouble>    m_A;
            Array<TwoD,NekDouble>    m_B;
            Array<TwoD,NekDouble>    m_U;
            Array<TwoD,NekDouble>    m_V;

        private: 
            
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

            template<typename FuncType>
            void ExplicitIntegration(NekDouble timestep,                       
                                     FuncType& f,
                                     const TripleArray& y_old,
                                     const Array<OneD,const NekDouble>& t_old,
                                     TripleArray& y_new,
                                     Array<OneD, NekDouble>& t_new)
            {
                ASSERTL1(y_old.num_elements()==m_numsteps,"Arguments not appropriate for this method.");                
                ASSERTL1(IsExplicitMethod()==true,"Implicit integration scheme cannot be used for explit integration");

                unsigned int i,j,k;    
                unsigned int nvar      = y_old[0].num_elements();
                unsigned int npoints   = y_old[0][0].num_elements();
                
                // First, we are going to calculate the various stage values and stage derivatives
                // (this is the multi-stage part of the method)
                // - Y corresponds to the stage values
                // - F corresponds to the stage derivatives
                // - T corresponds to the time at the different stages
                TripleArray Y(m_numstages);
                TripleArray F(m_numstages);
                Array<OneD,NekDouble> T(m_numstages,0.0);

                // Allocate memory for the arrays Y and F   
                for(i = 0; i < m_numstages; ++i)
                {    
                    Y[i] = DoubleArray(nvar);
                    F[i] = DoubleArray(nvar);
                    for(j = 0; j < nvar; j++)
                    {
                        Y[i][j] = Array<OneD, NekDouble>(npoints,0.0);
                        F[i][j] = Array<OneD, NekDouble>(npoints);
                    }
                }
    
                // The loop below calculates the stage values and derivatives
                for( i = 0; i < m_numstages; i++)
                {
                    // The stage values Y are a linear combination of:
                    // 1: the stage derivatives
                    for( j = 0; j < i; j++ )
                    {
                        for(k = 0; k < nvar; k++)
                        {
                            Vmath::Svtvp(npoints,timestep*m_A[i][j],F[j][k],1,
                                         Y[i][k],1,Y[i][k],1);
                        }

                        T[i] += m_A[i][j]*timestep;
                    }

                    // 2: the imported multi-step solution of the previous time level
                    for( j = 0; j < m_numsteps; j++)
                    {
                        for(k = 0; k < nvar; k++)
                        {
                            Vmath::Svtvp(npoints,m_U[i][j],y_old[j][k],1,
                                         Y[i][k],1,Y[i][k],1);
                        }
                        T[i] += m_U[i][j]*t_old[j];
                    }

                    // The stage derivative can now be calculated based
                    // on the stage value using the forcing function
                    f.ODEforcing(Y[i],F[i],T[i]);
                }

                // Next, the solution at the new time level will be calculated.
                // For multi-step methods, this includes updating the values 
                // of the auxiliary paremeters

                // Make sure that the y_new is of same dimenension as y_old
                // and that its values are set to zero  
                if(y_new.num_elements() != m_numsteps)
                {
                    y_new = TripleArray(m_numsteps);
                    for(i = 0; i < m_numsteps; ++i)
                    { 
                        y_new[i] = DoubleArray(nvar);
                        for(j = 0; j < nvar; j++)
                        {
                            y_new[i][j] = Array<OneD,NekDouble>(npoints,0.0);
                        } 
                    }
                }
                else
                {
                    for(i = 0; i < m_numsteps; ++i)
                    { 
                        if(y_new[i].num_elements() != nvar)
                        {
                            y_new[i] = DoubleArray(nvar);
                            for(j = 0; j < nvar; j++)
                            {
                                y_new[i][j] = Array<OneD,NekDouble>(npoints,0.0);
                            } 
                        }
                        else
                        {
                            for(j = 0; j < nvar; j++)
                            {
                                if(y_new[i][j].num_elements() != npoints)
                                {
                                    y_new[i][j] = Array<OneD,NekDouble>(npoints,0.0);
                                }
                                else
                                {
                                    Vmath::Zero(npoints,y_new[i][j],1);
                                }
                            } 
                        }
                    }
                }


                if(t_new.num_elements() != m_numsteps)
                {
                    t_new = Array<OneD,NekDouble>(m_numsteps,0.0);
                }
                else
                {
                    Vmath::Zero(m_numsteps,t_new,1);
                }
                
                // The loop below calculates the solution at the new time level
                for( i = 0; i < m_numsteps; i++)
                {
                    // The solution at the new time level is a linear combination of:
                    // 1: the stage derivatives
                    for( j = 0; j < m_numstages; j++ )
                    {
                        for(k = 0; k < nvar; k++)
                        {
                            Vmath::Svtvp(npoints,timestep*m_B[i][j],F[j][k],1,
                                         y_new[i][k],1,y_new[i][k],1);
                        }
                        t_new[i] += m_B[i][j]*timestep;

                    }
        
                    // 2: the imported multi-step solution of the previous time level
                    for( j = 0; j < m_numsteps; j++ )
                    {
                        for(k = 0; k < nvar; k++)
                        {
                            Vmath::Svtvp(npoints,m_V[i][j],y_old[j][k],1,
                                         y_new[i][k],1,y_new[i][k],1);
                        }
                        t_new[i] += m_V[i][j]*t_old[j];
                    }
                }
            }

        };

        std::ostream& operator<<(std::ostream& os, const TimeIntegrationScheme& rhs);
        
        // =========================================================================

    }; // end of namespace
} // end of namespace 

#endif //NEKTAR_LIB_UTILITIES_FOUNDATIONS_TIMEINTEGRATIONSCHEME_H
