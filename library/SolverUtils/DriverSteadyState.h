///////////////////////////////////////////////////////////////////////////////
//
// File AdvectionTerm.h
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
// Description: Driver class for the stability solver
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DRIVERSTEADYSTATE_H
#define NEKTAR_SOLVERUTILS_DRIVERSTEADYSTATE_H

#include <SolverUtils/Driver.h>

////////////////////
#include <SolverUtils/DriverArnoldi.h>
#include <SolverUtils/DriverModifiedArnoldi.h>
///////////////////

namespace Nektar
{
    namespace SolverUtils
    {
        /// Base class for the development of solvers.
        //         class DriverSteadyState: public Driver//, public DriverArnoldi
        class DriverSteadyState: public DriverModifiedArnoldi
        {
        public:
            friend class MemoryManager<DriverSteadyState>;
            
            /// Creates an instance of this class
            static DriverSharedPtr create(const LibUtilities::SessionReaderSharedPtr& pSession) {
                DriverSharedPtr p = MemoryManager<DriverSteadyState>::AllocateSharedPtr(pSession);
                p->InitObject();
                return p;
            }
            
            ///Name of the class
            static std::string className;
            
            void ConvergenceHistory(const Array<OneD, const Array<OneD, NekDouble> > &qBar1,
                                    const Array<OneD, const Array<OneD, NekDouble> > &q0,
                                    NekDouble &MaxNormDiff_q_qBar,
                                    NekDouble &MaxNormDiff_q1_q0);
            
            void ComputeSFD(const int i,
                              const Array<OneD, const Array<OneD, NekDouble> > &q0,
                              const Array<OneD, const Array<OneD, NekDouble> > &qBar0,
                              Array<OneD, Array<OneD, NekDouble> > &q1,
                              Array<OneD, Array<OneD, NekDouble> > &qBar1);
            
            void EvalEV_ScalarSFD(const NekDouble &X_input,
                                  const NekDouble &Delta_input,
                                  const complex<NekDouble> &alpha,
                                  NekDouble &MaxEV);
            
            void GradientDescentMethod(const complex<NekDouble> &alpha,
                                       NekDouble &X_output,
                                       NekDouble &Delta_output);
            
            void ReadEVfile(const int &KrylovSubspaceDim, NekDouble &growthEV, NekDouble &frequencyEV);
            
            void ComputeOptimization();
            
            void SetSFDOperator(const NekDouble X_input, const NekDouble Delta_input);
            
            
        protected:
            /// Constructor
            SOLVER_UTILS_EXPORT DriverSteadyState(const LibUtilities::SessionReaderSharedPtr pSession);
            
            /// Destructor
            SOLVER_UTILS_EXPORT virtual ~DriverSteadyState();
            
            /// Second-stage initialisation
            SOLVER_UTILS_EXPORT virtual void v_InitObject(ostream &out = cout);
            
            /// Virtual function for solve implementation.
            SOLVER_UTILS_EXPORT virtual void v_Execute(ostream &out = cout);
            
            static std::string driverLookupId;
            
        private:
            int NumVar_SFD;  
            
            //Definition of the SFD parameters:
            NekDouble m_Delta;
            NekDouble m_X;
            NekDouble m_dt; 
            
            //For implementation of the exact solution of the filters equation
            NekDouble M11;
            NekDouble M12;
            NekDouble M21;
            NekDouble M22;
            
            int m_stepCounter;
            int m_Check;
            NekDouble TOL;
            
            
            int m_infosteps;
            int m_checksteps;
            
            int MPIrank;
            
            NekDouble Diff_q_qBar;
            NekDouble Diff_q1_q0;
            
            NekDouble Min_MaxNormDiff_q_qBar;
            
            //For coupling SFD and Arnoldi
            bool OptumiumParametersFound;
            NekDouble PartialTOL;
            NekDouble PartialTOL_init;
            NekDouble TimeToRestart;
            int m_NonConvergingStepsCounter;
            NekDouble ParametersTOL;
            NekDouble UpdateCoefficient;
            
            std::ofstream m_file;
        };
    }   
} //end of namespace


#endif //NEKTAR_SOLVERUTILS_DRIVERSTEADYSTATE_H

