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
#include <SolverUtils/AdvectionSystem.h>

////////////////////////
#include <SolverUtils/DriverModifiedArnoldi.h>
#include <SolverUtils/DriverArnoldi.h>
///////////////////////


namespace Nektar
{
    namespace SolverUtils
    {
        string DriverSteadyState::className = GetDriverFactory().RegisterCreatorFunction("SteadyState", DriverSteadyState::create);
        string DriverSteadyState::driverLookupId = LibUtilities::SessionReader::RegisterEnumValue("Driver","SteadyState",0);
        
        /**
         * 
         */
        //         DriverSteadyState::DriverSteadyState(const LibUtilities::SessionReaderSharedPtr pSession)
        //         : Driver(pSession)
        //         {
            //         }
            
            DriverSteadyState::DriverSteadyState(const LibUtilities::SessionReaderSharedPtr pSession)
            : DriverModifiedArnoldi(pSession)
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
                DriverModifiedArnoldi::v_InitObject(out);
            }
            
            
            void DriverSteadyState::v_Execute(ostream &out)
            
            {
                //With a loop over "DoSolve", this Driver implements the "encaplulated" Selective Frequency Damping method
                //to find the steady state of a flow above the critical Reynolds number.
                
                m_equ[m_nequ - 1]->PrintSummary(out);
                m_equ[m_nequ - 1]->DoInitialise();
                
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////
                AdvectionSystemSharedPtr A = boost::dynamic_pointer_cast<AdvectionSystem>(m_equ[0]);
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////
                
                // - SFD Routine -
                // Compressible case
                NumVar_SFD = m_equ[m_nequ - 1]->UpdateFields()[0]->GetCoordim(0);            
                if (m_session->GetSolverInfo("EqType") == "EulerCFE" || 
                    m_session->GetSolverInfo("EqType") == "NavierStokesCFE")
                {
                    NumVar_SFD += 2; //Number of variables for the compressible equations
                }
                
                Array<OneD, Array<OneD, NekDouble> > q0(NumVar_SFD);
                Array<OneD, Array<OneD, NekDouble> > q1(NumVar_SFD);
                Array<OneD, Array<OneD, NekDouble> > qBar0(NumVar_SFD);
                Array<OneD, Array<OneD, NekDouble> > qBar1(NumVar_SFD);
                Array<OneD, Array<OneD, NekDouble> > partialSteadyFlow(NumVar_SFD); //Only for coupling with SFD
                
                NekDouble TOL(0);
                m_n=0;
                m_Check=0;          
                
                m_session->LoadParameter("FilterWidth", m_Delta0, 1);
                m_session->LoadParameter("ControlCoeff",m_X0, 1);
                m_session->LoadParameter("TOL", TOL, 1.0e-08);
                //             m_session->LoadParameter("TOL", TOL, 1.0e-04);
                m_session->LoadParameter("IO_InfoSteps", m_infosteps, 1000);
                m_session->LoadParameter("IO_CheckSteps", m_checksteps, 100000);
                
                m_dt = m_equ[m_nequ - 1]->GetTimeStep();
                
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
                
                m_equ[m_nequ - 1]->SetStepsToOne(); //m_steps is set to 1. Then "m_equ[m_nequ - 1]->DoSolve()" will run for only one time step      
                ofstream m_file("ConvergenceHistory.txt", ios::out | ios::trunc);
                
                for(int i = 0; i < NumVar_SFD; ++i)
                {
                    q0[i] = Array<OneD, NekDouble> (m_equ[m_nequ - 1]->GetTotPoints(), 0.0); //q0 is initialised
                    qBar0[i] = Array<OneD, NekDouble> (m_equ[m_nequ - 1]->GetTotPoints(), 0.0);
                    m_equ[m_nequ - 1]->CopyFromPhysField(i, qBar0[i]); //qBar0 is initially set at beiing equal to the initial conditions provided in the input file
                    
                    partialSteadyFlow[i] = Array<OneD, NekDouble> (m_equ[m_nequ - 1]->GetTotPoints(), 0.0); //Only for coupling with SFD
                }
                
                Diff_q_qBar = 1.0;
                Diff_q1_q0 = 1.0;
                MinDiff = 1.0;
                int OscillationCounter(0);
                
                //while (max(Diff_q_qBar, Diff_q1_q0) > TOL)
                while (min(Diff_q_qBar, Diff_q1_q0) > TOL)
                {
                    //First order Splitting with exact resolution of the filters equation
                    m_equ[m_nequ - 1]->DoSolve();
                    
                    for(int i = 0; i < NumVar_SFD; ++i)
                    {
                        m_equ[m_nequ - 1]->CopyFromPhysField(i, q0[i]);
                        
                        //Apply the linear operator F to the outcome of the solver
                        ExactFilters(i, q0, qBar0, q1, qBar1);
                        
                        qBar0[i] = qBar1[i];
                        m_equ[m_nequ - 1]->CopyToPhysField(i, q1[i]);     
                    }
                    
                    if(m_infosteps && !((m_n+1)%m_infosteps))
                    {
                        ConvergenceHistory(qBar1, q1, Diff_q_qBar, Diff_q1_q0);
                        
                        //////////////////////////////////////////////////////////////////////////////////////////
                        if (m_EvolutionOperator == eOptimizedSteadyState)
                        {
                            ///cout << "We UPDATE the base floooooow !!!!" << endl;
                            
                            
                            ///We add this block for defining when we run stability-analysis and optimisation of the parameters
                            if(MinDiff > Diff_q1_q0) //The norm is decreasing (the flow is getting closer to its steady-state)
                            {
                                MinDiff = Diff_q1_q0;
                                
                                for(int i = 0; i < NumVar_SFD; ++i)
                                {
                                    Vmath::Vcopy(q1[i].num_elements(), q1[i], 1, partialSteadyFlow[i], 1);
                                }
                                
                                if (MinDiff < 10.0*TOL)
                                {
                                    cout << "\n\t We compute optimization because the flow is partially converged ! \n" << endl;
                                    // -----> UpdateBaseFlow(partialSteadyFlow)
                                    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                    A->GetAdvObject()->SetBaseFlow(partialSteadyFlow);
                                    m_equ[m_nequ - 1]->Checkpoint_Output(010101010); //We save the flow field into a .chk file 
                                    ////////////////////////////////////////////////////////////////////////////////////////////////////////////                                    
                                    DriverModifiedArnoldi::v_Execute(out);
                                    ComputeOptimization();
                                }
                            }
                            else  //The norm is increasing
                            {
                                OscillationCounter += 1;
                                
                                if (Diff_q1_q0 > 100.0*MinDiff || OscillationCounter > 25)
                                {
                                    cout << "\n\t We compute optimization because the flow NOT converging towards SS ! \n" << endl;
                                    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                    A->GetAdvObject()->SetBaseFlow(partialSteadyFlow);
                                    m_equ[m_nequ - 1]->Checkpoint_Output(010101010); //We save the flow field into a .chk file 
                                    ////////////////////////////////////////////////////////////////////////////////////////////////////////////                                       
                                    DriverModifiedArnoldi::v_Execute(out);
                                    ComputeOptimization();
                                    OscillationCounter = 0;
                                }
                            }
                        }
                        //////////////////////////////////////////////////////////////////////////////////////////
                        
                    }
                    
                    if(m_checksteps && m_n&&(!((m_n+1)%m_checksteps)))
                    {                      
                        m_Check++;                    
                        m_equ[m_nequ - 1]->Checkpoint_Output(m_Check);                       
                    }
                    
                    m_n++;
                }
                
                
                m_Check++;
                m_equ[m_nequ - 1]->Checkpoint_Output(m_Check);
                
                m_file.close();
                // - End SFD Routine -
                
                m_equ[m_nequ - 1]->Output();
                
                // Evaluate and output computation time and solution accuracy.
                // The specific format of the error output is essential for the
                // regression tests to work.
                // Evaluate L2 Error
                for(int i = 0; i < NumVar_SFD; ++i)
                {
                    NekDouble vL2Error = m_equ[m_nequ - 1]->L2Error(i,false);
                    NekDouble vLinfError = m_equ[m_nequ - 1]->LinfError(i);
                    if (m_comm->GetRank() == 0)
                    {
                        out << "L 2 error (variable " << m_equ[m_nequ - 1]->GetVariable(i) << ") : " << vL2Error << endl;
                        out << "L inf error (variable " << m_equ[m_nequ - 1]->GetVariable(i) << ") : " << vLinfError << endl;
                    }
                }
            }
            
            
            
            void DriverSteadyState::ComputeOptimization()
            {            
                NekDouble growthEV(0.0);
                NekDouble frequencyEV(0.0);
                
                ///m_kdim is the dimension of Krylov subspace
                ReadEVfile(m_kdim, growthEV, frequencyEV);
                
                cout << "growthEV = " << growthEV << endl;
                cout << "frequencyEV = " << frequencyEV << endl;
                cout << "Norm EV = " << exp(growthEV) << endl;
                
                complex<NekDouble> ApproxEV = polar(exp(growthEV), frequencyEV);   
                
                NekDouble X_opt;
                NekDouble Delta_opt;
                X_opt = m_X/m_dt;
                Delta_opt = m_Delta*m_dt;
                GradientDescentMethod(ApproxEV, X_opt, Delta_opt);
            }
            
            
            void DriverSteadyState::ExactFilters(const int i,
                                                 const Array<OneD, const Array<OneD, NekDouble> > &q0,
                                                 const Array<OneD, const Array<OneD, NekDouble> > &qBar0,
                                                 Array<OneD, Array<OneD, NekDouble> > &q1,
                                                 Array<OneD, Array<OneD, NekDouble> > &qBar1)
            {
                q1[i] = Array<OneD, NekDouble> (m_equ[m_nequ - 1]->GetTotPoints(),0.0);
                qBar1[i] = Array<OneD, NekDouble> (m_equ[m_nequ - 1]->GetTotPoints(),0.0);
                
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
                Array<OneD, NekDouble > NormDiff_q_qBar(NumVar_SFD, 1.0);
                Array<OneD, NekDouble > NormDiff_q1_q0(NumVar_SFD, 1.0);
                
                MaxNormDiff_q_qBar=0.0;
                MaxNormDiff_q1_q0=0.0;
                
                //Norm Calculation
                for(int i = 0; i < NumVar_SFD; ++i)
                {
                    //To check convergence of SFD
                    //NormDiff_q_qBar[i] = m_equ[m_nequ - 1]->L2Error(i, qBar1[i], false);
                    NormDiff_q_qBar[i] = m_equ[m_nequ - 1]->LinfError(i, qBar1[i]);
                    
                    //To check convergence of Navier-Stokes
                    NormDiff_q1_q0[i] = m_equ[m_nequ - 1]->LinfError(i, q0[i]);
                    
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
                    cout << "SFD (MPI) - Step: " << m_n+1 << "; Time: " << m_equ[m_nequ - 1]->GetFinalTime() <<  "; |q-qBar|inf = " << MaxNormDiff_q_qBar << "; |q1-q0|inf = " << MaxNormDiff_q1_q0 << ";\t for X = " << m_X0 <<" and Delta = " << m_Delta0 <<endl;
                    std::ofstream m_file( "ConvergenceHistory.txt", std::ios_base::app); 
                    m_file << m_n+1 << "\t" << m_equ[m_nequ - 1]->GetFinalTime() << "\t" << MaxNormDiff_q_qBar << "\t" << MaxNormDiff_q1_q0 << endl;
                    m_file.close();
                }
                #else
                cout << "SFD - Step: " << m_n+1 << "; Time: " << m_equ[m_nequ - 1]->GetFinalTime() <<  "; |q-qBar|inf = " << MaxNormDiff_q_qBar << "; |q1-q0|inf = " << MaxNormDiff_q1_q0 << ";\t for X = " << m_X0 <<" and Delta = " << m_Delta0 <<endl;
                std::ofstream m_file( "ConvergenceHistory.txt", std::ios_base::app); 
                m_file << m_n+1 << "\t" << m_equ[m_nequ - 1]->GetFinalTime() << "\t" << MaxNormDiff_q_qBar << "\t" << MaxNormDiff_q1_q0 << endl;
                m_file.close();
                #endif              
            }
            
            
            
            void DriverSteadyState::GradientDescentMethod(const complex<NekDouble> &alpha,
                                                          NekDouble &X_output,
                                                          NekDouble &Delta_output)
            {
                //This routine implements a gradient descent method to find the parameters X end Delta which give the minimum 
                //eigenlavue of the SFD problem applied to the scalar case.
                
                bool OptParmFound = false;
                bool Descending = true;
                
                NekDouble X_input = 1.0;
                NekDouble Delta_input = 1.0;
                
                NekDouble F0(0.0);
                NekDouble F0x(0.0);
                NekDouble F0y(0.0);
                NekDouble F1(0.0);
                NekDouble dx = 0.00000001;
                NekDouble dirx(0.0);
                NekDouble diry(0.0);
                NekDouble s(0.0);
                NekDouble CurrentMin = 1.0;
                
                while (OptParmFound == false)
                {                
                    Descending = true;
                    
                    EvalEV_ScalarSFD(X_input, Delta_input, alpha, F0);
                    EvalEV_ScalarSFD(X_input + dx, Delta_input, alpha, F0x);
                    EvalEV_ScalarSFD(X_input, Delta_input + dx, alpha, F0y);
                    
                    dirx =  (F0 - F0x);
                    diry =  (F0 - F0y);
                    
                    s = abs(0.000001/dirx);
                    X_output = X_input + s*dirx;
                    Delta_output = Delta_input + s*diry;
                    
                    F1 = F0;
                    
                    while (Descending == true)
                    {
                        CurrentMin = F1;
                        
                        X_input = X_output;
                        Delta_input = Delta_output;
                        
                        EvalEV_ScalarSFD(X_output, Delta_output, alpha, F1);
                        
                        if (F1 > CurrentMin)
                        {
                            Descending = false;
                        }
                        else
                        {
                            s = s + s*0.01;  
                            X_output = X_input + s*dirx;
                            Delta_output = Delta_input + s*diry;
                        }
                    }
                    
                    if (abs(F0-F1) < dx)
                    {
                        EvalEV_ScalarSFD(X_output, Delta_output, alpha, F1);
                        cout << "\n \t The optimum paramters are: X_opt = " << X_output << " and Delta_opt = " << Delta_output << "\n" << endl;
                        cout << "\t The minimum EV is: " << F1 << "\n" << endl;
                        OptParmFound = true; 
                    }
                }            
            }
            
            
            
            void DriverSteadyState::EvalEV_ScalarSFD(const NekDouble &X_input,
                                                     const NekDouble &Delta_input,
                                                     const complex<NekDouble> &alpha,
                                                     NekDouble &MaxEV)
            {
                //The routine evaluates the maximum eigenvalue of the SFD system when applied to the 1D problem u(n+1) = alpha*u(n)
                //This will help to triger the optimum parameters for the SFD applied to the flow equations
                
                NekDouble A11 = ( 1.0 + X_input*Delta_input*exp(-(X_input + 1.0/Delta_input)) )/(1.0 + X_input*Delta_input);
                NekDouble A12 = ( X_input*Delta_input - X_input*Delta_input*exp(-(X_input + 1.0/Delta_input)) )/(1.0 + X_input*Delta_input);
                NekDouble A21 = ( 1.0 - 1.0*exp(-(X_input + 1.0/Delta_input)) )/(1 + X_input*Delta_input);
                NekDouble A22 = ( X_input*Delta_input + 1.0*exp(-(X_input + 1.0/Delta_input)) )/(1.0 + X_input*Delta_input);
                
                complex<NekDouble> B11 = alpha;
                NekDouble B12 = 0.0;
                NekDouble B21 = 0.0;
                NekDouble B22 = 1.0;
                
                complex<NekDouble> a = A11*B11 + A12*B21;
                NekDouble b = A11*B12 + A12*B22;
                complex<NekDouble> c = A21*B11 + A22*B21;
                NekDouble d = A21*B12 + A22*B22;
                
                complex<NekDouble> delt = sqrt((a-d)*(a-d) + 4.0*b*c);
                complex<NekDouble> lambda_1 = 0.5*(a+d + delt);
                complex<NekDouble> lambda_2 = 0.5*(a+d - delt);
                
                NekDouble NORM_1 = abs(lambda_1);
                NekDouble NORM_2 = abs(lambda_2);
                
                MaxEV = max(NORM_1, NORM_2);                
            }
            
            
            void DriverSteadyState::ReadEVfile(const int &KrylovSubspaceDim, NekDouble &growthEV, NekDouble &frequencyEV)
            {           
                std::string EVfileName = m_session->GetSessionName() +  ".evl";
                std::ifstream EVfile(EVfileName.c_str());
                
                int NumLinesInFile(0);
                NekDouble NonReleventNumber(0.0);
                
                if(EVfile)
                {                
                    std::string line; 
                    
                    //This block count the total number of lines of the .evl file
                    while(getline(EVfile, line)) //Tant qu'on n'est pas à la fin, on lit
                {
                    NumLinesInFile += 1;
                }                
                EVfile.close();
                
                std::ifstream EVfile(EVfileName.c_str()); //go back to the beginning of the file
                
                //We now want to go to the the line where the most unstable eigenlavue was written
                for(int i = 0; i < (NumLinesInFile - KrylovSubspaceDim); ++i) 
                {
                    std::getline(EVfile, line);
                }
                
                //Then we read this line by skipping the first three values written
                EVfile >> NonReleventNumber;
                EVfile >> NonReleventNumber;
                EVfile >> NonReleventNumber;
                
                //The growth rate and the frequency of the EV are at the 4th and 5th colums of the .evl file 
                EVfile >> growthEV;
                EVfile >> frequencyEV;                
                }
                else
                {
                    cout << "An error occured when openning the .evl file" << endl;
                }
                EVfile.close();
            }
            
    }
}

