///////////////////////////////////////////////////////////////////////////////
//
// File DriverSteadyState.cpp
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
// Description: Incompressible Navier Stokes solver
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/DriverSteadyState.h>


namespace Nektar
{
    namespace SolverUtils
    {
        string DriverSteadyState::className = GetDriverFactory().RegisterCreatorFunction("SteadyState", DriverSteadyState::create);
        string DriverSteadyState::driverLookupId = LibUtilities::SessionReader::RegisterEnumValue("Driver","SteadyState",0);
        
        /**
         * 
         */
        DriverSteadyState::DriverSteadyState(const LibUtilities::SessionReaderSharedPtr pSession)
        : Driver(pSession)
        {
        }
        
        
        /**
         * 
         */
        DriverSteadyState:: ~DriverSteadyState()
        {
        }
        
        
        /**
         * 
         */
        void DriverSteadyState::v_InitObject(ostream &out)
        {
            Driver::v_InitObject(out);
        }
        
        
        void DriverSteadyState::v_Execute(ostream &out)
        
        {
            //With a loop over "DoSolve", this Driver implements the "encaplulated" Selective Frequency Damping method
            //to find the steady state of a flow above the critical Reynolds number.
            
            m_equ[0]->PrintSummary(out);
            m_equ[0]->DoInitialise();
            
            // - SFD Routine -
            NumElmVelocity = m_equ[0]->GetNumElmVelocity();
            
            Array<OneD, Array<OneD, NekDouble> > q0(NumElmVelocity);
            Array<OneD, Array<OneD, NekDouble> > q1(NumElmVelocity);
            Array<OneD, Array<OneD, NekDouble> > qBar0(NumElmVelocity);
            Array<OneD, Array<OneD, NekDouble> > qBar1(NumElmVelocity);
            
            NekDouble TOL(0);
            m_n=0;
            m_Check=0;			
            
            m_session->LoadParameter("FilterWidth", m_Delta0, 1);
            m_session->LoadParameter("ControlCoeff",m_X0, 1);
            m_session->LoadParameter("TOL", TOL, 1.0e-08);
            m_session->LoadParameter("IO_InfoSteps", m_infosteps, 1000);
            m_session->LoadParameter("IO_CheckSteps", m_checksteps, 100000);
            
            m_dt = m_equ[0]->GetTimeStep();
            
            //The time-step is applied here to clarify the matrix F expression
            m_X = m_X0*m_dt;
            m_Delta = m_Delta0/m_dt;
            
            //Exact solution of the Filters equation (ici on a X_tilde et Delta_tile!!!)
            c1 = 1.0/(1.0 + m_X*m_Delta);
            F11 = c1*(1.0 + m_X*m_Delta*exp(-(m_X + 1.0/m_Delta)));
            F12 = c1*(m_X*m_Delta*(1.0 - exp(-(m_X + 1.0/m_Delta))));
            F21 = c1*(1.0 - exp(-(m_X + 1.0/m_Delta)));
            F22 = c1*(m_X*m_Delta + exp(-(m_X + 1.0/m_Delta)));
            
            cout << "------------------ SFD Parameters ------------------" << endl;
            cout << "\tX = " << m_X0 << endl;
            cout << "\tDelta = " << m_Delta0 << endl;
            cout << "----------------------------------------------------" << endl;
            
            m_equ[0]->SetStepsToOne(); //m_steps is set to 1. Then "m_equ[0]->DoSolve()" will run for only one time step			
            ofstream m_file("ConvergenceHistory.txt", ios::out | ios::trunc);
            
            for(int i = 0; i < NumElmVelocity; ++i)
            {
                q0[i] = Array<OneD, NekDouble> (m_equ[0]->GetTotPoints(), 0.0); //q0 is initialised
                qBar0[i] = Array<OneD, NekDouble> (m_equ[0]->GetTotPoints(), 0.0);
                m_equ[0]->CopyFromPhysField(i, qBar0[i]); //qBar0 is initially set at beiing equal to the initial conditions provided in the input file
            }
            
            MaxNormDiff_q_qBar = 1.0;
            MaxNormDiff_q1_q0 = 1.0;
            Min_MaxNormDiff_q_qBar = MaxNormDiff_q_qBar;
            
            
            while (max(MaxNormDiff_q_qBar, MaxNormDiff_q1_q0) > TOL)
            {
                //First order Splitting with exact resolution of the filters equation
                m_equ[0]->DoSolve();
                
                for(int i = 0; i < NumElmVelocity; ++i)
                {
                    m_equ[0]->CopyFromPhysField(i, q0[i]);
                    
                    //Apply the linear operator F to the outcome of the solver
                    ExactFilters(i, q0, qBar0, q1, qBar1);
                    
                    qBar0[i] = qBar1[i];
                    m_equ[0]->CopyToPhysField(i, q1[i]);     
                }
                
                
                if(m_infosteps && !((m_n+1)%m_infosteps))
                {
                    ConvergenceHistory(qBar1, q0, MaxNormDiff_q_qBar, MaxNormDiff_q1_q0);
                }
                
                if(m_checksteps && m_n&&(!((m_n+1)%m_checksteps)))
                {                      
                    m_Check++;                    
                    m_equ[0]->Checkpoint_Output(m_Check);                       
                }
                
                m_n++;
            }
            
            cout << "\nFINAL Filter Width: Delta = " << m_Delta << "; FINAL Control Coeff: X = " << m_X << "\n" << endl;
            m_Check++;
            m_equ[0]->Checkpoint_Output(m_Check);
            
            m_file.close();
            // - End SFD Routine -
            
            m_equ[0]->Output();
            
            // Evaluate and output computation time and solution accuracy.
            // The specific format of the error output is essential for the
            // regression tests to work.
            // Evaluate L2 Error
            for(int i = 0; i < m_equ[0]->GetNvariables(); ++i)
            {
                NekDouble vL2Error = m_equ[0]->L2Error(i,false);
                NekDouble vLinfError = m_equ[0]->LinfError(i);
                if (m_comm->GetRank() == 0)
                {
                    out << "L 2 error (variable " << m_equ[0]->GetVariable(i) << ") : " << vL2Error << endl;
                    out << "L inf error (variable " << m_equ[0]->GetVariable(i) << ") : " << vLinfError << endl;
                }
            }
        }
        
        
        void DriverSteadyState::ExactFilters(const int i,
                                             const Array<OneD, const Array<OneD, NekDouble> > &q0,
                                             const Array<OneD, const Array<OneD, NekDouble> > &qBar0,
                                             Array<OneD, Array<OneD, NekDouble> > &q1,
                                             Array<OneD, Array<OneD, NekDouble> > &qBar1)
        {
            q1[i] = Array<OneD, NekDouble> (m_equ[0]->GetTotPoints(),0.0);
            qBar1[i] = Array<OneD, NekDouble> (m_equ[0]->GetTotPoints(),0.0);
            
            //Exact solution of the Filters equation            
            Vmath::Svtvp(q1[i].num_elements(), F11, q0[i], 1, q1[i], 1, q1[i], 1 );    
            Vmath::Svtvp(q1[i].num_elements(), F12, qBar0[i], 1, q1[i], 1, q1[i], 1 ); 
            
            Vmath::Svtvp(qBar1[i].num_elements(), F21, q0[i], 1, qBar1[i], 1, qBar1[i], 1 );    
            Vmath::Svtvp(qBar1[i].num_elements(), F22, qBar0[i], 1, qBar1[i], 1, qBar1[i], 1 );             
        }
        
        
        void DriverSteadyState::ConvergenceHistory(const Array<OneD, const Array<OneD, NekDouble> > &qBar1,
                                                   const Array<OneD, const Array<OneD, NekDouble> > &q0,
                                                   NekDouble &MaxNormDiff_q_qBar,
                                                   NekDouble &MaxNormDiff_q1_q0)
        {
            //This routine evaluates |q-qBar|_L2 and save the value in "ConvergenceHistory.txt"
            
            Array<OneD, NekDouble > NormDiff_q_qBar(NumElmVelocity, 1.0);
            Array<OneD, NekDouble > NormDiff_q1_q0(NumElmVelocity, 1.0);
            
            MaxNormDiff_q_qBar=0.0;
            MaxNormDiff_q1_q0=0.0;
            
            //Norm Calculation
            for(int i = 0; i < NumElmVelocity; ++i)
            {
                //To check convergence of SFD
                //NormDiff_q_qBar[i] = m_equ[0]->L2Error(i, qBar1[i], false);
                NormDiff_q_qBar[i] = m_equ[0]->LinfError(i, qBar1[i]);
                
                //To check convergence of Navier-Stokes
                NormDiff_q1_q0[i] = m_equ[0]->LinfError(i, q0[i]);
                
                if (MaxNormDiff_q_qBar < NormDiff_q_qBar[i])
                {
                    MaxNormDiff_q_qBar = NormDiff_q_qBar[i];
                }
                
                if (MaxNormDiff_q1_q0 < NormDiff_q1_q0[i])
                {
                    MaxNormDiff_q1_q0 = NormDiff_q1_q0[i];
                }
            }
            
            #if NEKTAR_USE_MPI   
            MPI_Comm_rank(MPI_COMM_WORLD,&MPIrank);
            if (MPIrank==0)
            {
                cout << "SFD (MPI) - Step: " << m_n+1 << "; Time: " << m_equ[0]->GetFinalTime() <<  "; |q-qBar|inf = " << MaxNormDiff_q_qBar << "; |q1-q0|inf = " << MaxNormDiff_q1_q0 << ";\t for X = " << m_X0 <<" and Delta = " << m_Delta0 <<endl;
                std::ofstream m_file( "ConvergenceHistory.txt", std::ios_base::app); 
                m_file << m_n+1 << "\t" << m_equ[0]->GetFinalTime() << "\t" << MaxNormDiff_q_qBar << "\t" << MaxNormDiff_q1_q0 << endl;
                m_file.close();
            }
            #else
            cout << "SFD - Step: " << m_n+1 << "; Time: " << m_equ[0]->GetFinalTime() <<  "; |q-qBar|inf = " << MaxNormDiff_q_qBar << "; |q1-q0|inf = " << MaxNormDiff_q1_q0 << ";\t for X = " << m_X0 <<" and Delta = " << m_Delta0 <<endl;
            std::ofstream m_file( "ConvergenceHistory.txt", std::ios_base::app); 
            m_file << m_n+1 << "\t" << m_equ[0]->GetFinalTime() << "\t" << MaxNormDiff_q_qBar << "\t" << MaxNormDiff_q1_q0 << endl;
            m_file.close();
            #endif              
            
        }
    }
}

