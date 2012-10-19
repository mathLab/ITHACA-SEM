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
			m_equ[0]->PrintSummary(out);
			m_equ[0]->DoInitialise();
			
			//--------------------------------------------------------- SFD Routine ---------------------------------------------------------
			int NumElmVelocity(0);
			NumElmVelocity = m_equ[0]->GetNumElmVelocity();
			cout << "NumElmVelocity = " << NumElmVelocity << endl;
			
			Array<OneD, Array<OneD, NekDouble> > q0(NumElmVelocity);
			Array<OneD, Array<OneD, NekDouble> > q1(NumElmVelocity);
			Array<OneD, Array<OneD, NekDouble> > qBar0(NumElmVelocity);
			Array<OneD, Array<OneD, NekDouble> > qBar1(NumElmVelocity);
			
			Array<OneD, NekDouble > NormDiff_q_qBar(NumElmVelocity, 1.0);
			Array<OneD, Array<OneD, NekDouble> > Diff_q_qBar(NumElmVelocity);
			
			//Definition of the SFD parameters:
			NekDouble Delta;
			NekDouble X;
			NekDouble cst1;
			NekDouble cst2;
			NekDouble cst3;
			NekDouble cst4;
			NekDouble cst5;
			NekDouble TOL(0);
			int Check(0);
			int n(0);
			NekDouble MinNormDiff_q_qBar(1000);
			NekDouble MaxNormDiff_q_qBar(0);
			NekDouble First_MinNormDiff_q_qBar(0);
			bool Growing=false;
			bool Shrinking=false;
			int Oscillation(0);
			
			m_session->LoadParameter("FilterWidth", Delta);
			m_session->LoadParameter("ControlCoeff", X);
			m_session->LoadParameter("TOL", TOL);
			
			int infosteps(0);
			int checksteps(0);
			m_session->LoadParameter("IO_InfoSteps", infosteps);
			m_session->LoadParameter("IO_CheckSteps", checksteps);
			
			NekDouble dt;
			dt = m_equ[0]->GetTimeStep();
			
			cst1=X*dt;
			cst2=1.0/(1.0 + cst1);
			cst3=dt/Delta;
			cst4=cst2*cst3;
			cst5=1.0/(1.0 + cst3*(1.0-cst1*cst2));
			
			cout << "------------------ SFD Parameters ------------------" << endl;
			cout << "\tDelta = " << Delta << endl;
			cout << "\tX = " << X << endl;
			cout << "----------------------------------------------------" << endl;
			
			m_equ[0]->SetStepsToOne();			
			
			for(int i = 0; i < NumElmVelocity; ++i)
			{
				q0[i] = Array<OneD, NekDouble> (m_equ[0]->GetTotPoints(), 0.0);
				qBar0[i] = Array<OneD, NekDouble> (m_equ[0]->GetTotPoints(), 0.0);
			}
			
		while (NormDiff_q_qBar[0] > TOL)
			{
				m_equ[0]->DoSolve();
				
				for(int i = 0; i < NumElmVelocity; ++i)
				{						
					m_equ[0]->CopyFromPhysField(i, q0[i]);
					
					q1[i] = Array<OneD, NekDouble> (m_equ[0]->GetTotPoints(),0.0);
					qBar1[i] = Array<OneD, NekDouble> (m_equ[0]->GetTotPoints(),0.0);
					
					Diff_q_qBar[i] = Array<OneD, NekDouble> (m_equ[0]->GetTotPoints(),0.0);
					
					Vmath::Smul(qBar1[i].num_elements(), cst4, q0[i], 1, qBar1[i], 1);
					Vmath::Vadd(qBar1[i].num_elements(), qBar0[i], 1, qBar1[i], 1, qBar1[i], 1);
					Vmath::Smul(qBar1[i].num_elements(), cst5, qBar1[i], 1, qBar1[i], 1);
					
					Vmath::Smul(q1[i].num_elements(), cst1, qBar1[i], 1, q1[i], 1);
					Vmath::Vadd(q1[i].num_elements(), q0[i], 1, q1[i], 1, q1[i], 1);
					Vmath::Smul(q1[i].num_elements(), cst2, q1[i], 1, q1[i], 1);
					
					qBar0[i] = qBar1[i];
					
					m_equ[0]->CopyToPhysField(i, q1[i]);
				}
				
				
				if(infosteps && !((n+1)%infosteps))
				{
					//Norm Calculation
					Vmath::Vsub(Diff_q_qBar[0].num_elements(), q1[0], 1, qBar1[0], 1, Diff_q_qBar[0], 1);
					Vmath::Smul(Diff_q_qBar[0].num_elements(), cst1, Diff_q_qBar[0], 1, Diff_q_qBar[0], 1);
					
					//for(int i = 0; i < m_velocity.num_elements(); ++i)
					for(int i = 0; i < 1; ++i)
					{
						NormDiff_q_qBar[i] = 0.0;
						for(int j = 0; j < Diff_q_qBar[i].num_elements(); ++j)
						{
							NormDiff_q_qBar[i] += Diff_q_qBar[i][j]*Diff_q_qBar[i][j];
						}
						NormDiff_q_qBar[i]=sqrt(NormDiff_q_qBar[i]);
					}					
					cout << "SFD - Step: " << n+1 << "; Time: " << m_equ[0]->GetFinalTime() <<  "; |q-qBar|L2 = " << NormDiff_q_qBar[0] <<endl;
					
					std::ofstream file( "Zfichier.txt", std::ios_base::app ); 
					file << n << "\t" << NormDiff_q_qBar[0] << endl;
					file.close();
					
					
					if (Shrinking==false && NormDiff_q_qBar[0] > MaxNormDiff_q_qBar)
					{
						MaxNormDiff_q_qBar=NormDiff_q_qBar[0];
						Growing=true;
						if (MaxNormDiff_q_qBar < First_MinNormDiff_q_qBar)
						{
							Oscillation = 0;
							//cout << "\nOscillation REINITIALIZED !!!\n" << endl; 
							First_MinNormDiff_q_qBar=0;
						}
					}
					
					if (NormDiff_q_qBar[0] < MaxNormDiff_q_qBar && Growing==true)
					{
						Growing=false;
						MaxNormDiff_q_qBar=0;
					}
					
					if (Growing==false && NormDiff_q_qBar[0] < MinNormDiff_q_qBar)
					{
						MinNormDiff_q_qBar=NormDiff_q_qBar[0];
						Shrinking=true;
						if (Oscillation==0)
						{
							First_MinNormDiff_q_qBar=MinNormDiff_q_qBar;
							//cout << "\nNew Oscillation MIN = " << First_MinNormDiff_q_qBar << "\n" << endl; 
						}
					}
					
					if (NormDiff_q_qBar[0] > MinNormDiff_q_qBar && Shrinking==true)
					{
						Shrinking=false;
						MinNormDiff_q_qBar=1000;
						Oscillation=Oscillation+1;
						//cout << "\nOscillation = " << Oscillation << "\n" << endl; 
					}
					
					if (Oscillation==25)
					{					
						Delta = Delta + 0.25;
						
						X = 0.99*(1.0/Delta);
						
						cst1=X*dt;
						cst2=1.0/(1.0 + cst1);
						cst3=dt/Delta;
						cst4=cst2*cst3;
						cst5=1.0/(1.0 + cst3*(1.0-cst1*cst2));
						
						cout << "\nNew Filter Width: Delta = " << Delta << "; New Control Coeff: X = " << X << "\n" << endl;
						
						Oscillation=0;
						
						//We reinitialize the problem:
						m_equ[0]->DoInitialise();
						/*m_time=0.0;
						//DoInitialise();
						SetBoundaryConditions(0.0);
						SetInitialConditions(0.0);*/
					}
				}
				
				if(checksteps && n&&(!((n+1)%checksteps)))
				{
					Check++;
					m_equ[0]->Checkpoint_Output(Check);
				}
				n++;
			}
			//--------------------------------------------------------- End SFD Routine ---------------------------------------------------------

			
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
    }
}

