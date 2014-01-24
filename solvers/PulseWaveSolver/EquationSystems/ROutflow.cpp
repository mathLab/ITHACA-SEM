///////////////////////////////////////////////////////////////////////////////
//
// File CommMpi.cpp
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
// Description: MPI communication implementation
//
///////////////////////////////////////////////////////////////////////////////

#include <PulseWaveSolver/EquationSystems/ROutflow.h>

namespace Nektar
{

    std::string Routflow::className
    = GetFlowFactory().RegisterCreatorFunction(
        "Resistance",
        Routflow::create,
        "Resistive outflow boundary condition");

    /**
     *
     */
    Routflow::Routflow(Array<OneD, MultiRegions::ExpListSharedPtr> pVessel, const LibUtilities::SessionReaderSharedPtr pSession)
      : PulseWaveFlow(pVessel,pSession)
    {
    }

    /**
     *
     */
    Routflow::~Routflow()
    {

    }

    void Routflow::v_DoBoundary(int omega, int n)
    { 
	/*NekDouble Q, A_r, u_r;
	NekDouble A_u, u_u;
        NekDouble R_t, A_l, u_l, u_0, c_0, c_l;

        Array<OneD, MultiRegions::ExpListSharedPtr> vessel(2);

        vessel[0] = m_vessels[2*omega];
        vessel[1] = m_vessels[2*omega+1];

        NekDouble RT=((vessel[0]->GetBndCondExpansions())[n])->GetCoeffs()[0];
        NekDouble pout = m_pout;
        int nq = vessel[0]->GetTotPoints(); 
                
        // Get the values of all variables needed for the Riemann problem
        A_l = m_fields[0]->GetCoeffs()[1];
        u_l = m_fields[1]->GetCoeffs()[1];

        // Call the R RiemannSolver
        R_RiemannSolver(RT,A_l,u_l,m_A_0[omega][nq-1],
                        m_beta[omega][nq-1],pout,A_u,u_u);
			
        // Calculates the new boundary conditions
        A_r=A_l;
        u_r=2*u_u-u_l;
                        
        // Store the updated values in the boundary condition
                        
        (vessel[0]->UpdateBndCondExpansion(n))->UpdatePhys()[0] = A_r;
        (vessel[1]->UpdateBndCondExpansion(n))->UpdatePhys()[0] = u_r;*/

    }

    void Routflow::R_RiemannSolver(NekDouble R,NekDouble A_l,NekDouble u_l,NekDouble A_0, 
                                               NekDouble beta, NekDouble pout,
                                               NekDouble &A_u,NekDouble &u_u)
    {		
        /*NekDouble W1 = 0.0;
        NekDouble c_l = 0.0;
        NekDouble pext = m_pext;
        NekDouble A_calc = 0.0;
        NekDouble fa = 0.0;
        NekDouble dfa = 0.0;
        NekDouble delta_A_calc = 0.0;
        NekDouble rho = m_rho;
        
        int proceed = 1;
        int iter = 0;
        int MAX_ITER = 200;
        
        // Tolerances for the algorithm
        NekDouble Tol = 1.0e-10;
        
        // Calculate the wave speed
        c_l = sqrt(beta/(2*rho))*sqrt(sqrt(A_l));
	
        // Riemann invariant \f$W_1(Al,ul)\f$
        W1 = u_l + 4*c_l;
        
        // Newton Iteration (Area only)
        A_calc = A_l;
        while ((proceed) && (iter < MAX_ITER))
        {	
            iter =iter+1;
	    
            fa = R*W1*A_calc-4*R*sqrt(beta/(2*rho))*A_calc*sqrt(sqrt(A_calc))-pext-beta*(sqrt(A_calc)-sqrt(A_0))+pout;
            dfa = R*W1-5*R*sqrt(beta/(2*rho))*sqrt(sqrt(A_calc))-beta/(2*sqrt(A_calc));
            delta_A_calc = fa/dfa;
            A_calc = A_calc - delta_A_calc;
            
            if (sqrt(delta_A_calc*delta_A_calc) < Tol)
                proceed = 0;
        }
        
	// Obtain u_u and A_u
	//u_u = W1 - 4*sqrt(beta/(2*rho))*(sqrt(sqrt(A_calc))); 
        u_u=(pext+beta*(sqrt(A_calc)-sqrt(A_0))-pout)/(R*A_calc);
        A_u = A_calc;
        */
    }


}
