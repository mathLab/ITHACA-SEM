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
// Description: ROuflow class
//
///////////////////////////////////////////////////////////////////////////////

#include <PulseWaveSolver/EquationSystems/TerminalOutflow.h>

using namespace std;

namespace Nektar
{

    std::string TerminalOutflow::className
    = GetBoundaryFactory().RegisterCreatorFunction(
        "Terminal",
        TerminalOutflow::create,
        "Terminal outflow boundary condition");

    /**
     *
     */
    TerminalOutflow::TerminalOutflow(Array<OneD, MultiRegions::ExpListSharedPtr> pVessel, 
                                     const LibUtilities::SessionReaderSharedPtr pSession,
                                     PulseWavePressureAreaSharedPtr pressureArea)
        : PulseWaveBoundary(pVessel,pSession,pressureArea)
    {
    }

    /**
     *
     */
    TerminalOutflow::~TerminalOutflow()
    {

    }

    void TerminalOutflow::v_DoBoundary(
        const Array<OneD,const Array<OneD, NekDouble> > &inarray,
        Array<OneD, Array<OneD, NekDouble> > &A_0,
        Array<OneD, Array<OneD, NekDouble> > &beta,
        const NekDouble time,
        int omega,int offset,int n)
    { 
	NekDouble A_r, u_r;
        NekDouble RT, A_l, u_l, u_0, c_0, c_l;

        Array<OneD, MultiRegions::ExpListSharedPtr> vessel(2);

        vessel[0] = m_vessels[2*omega];
        vessel[1] = m_vessels[2*omega+1];

        /* Find the terminal resistance boundary condition and
         * calculate the reflection. We assume A_r = A_l and
         * apply the reflection in u_r after paper
         * "Computational Modelling of 1D blood flow"*/
                        
        // Note: The R_t value is contained in A in the inputfile
        RT = (vessel[0]->UpdateBndCondExpansion(n))->GetCoeffs()[0];

        ASSERTL0((-1<=RT && RT<=1),
                 "RT must be comprised between -1 and 1");
        int nq = vessel[0]->GetTotPoints(); 
                        
        // Get the left values A_l and u_l needed for Eq. 37
        A_l = inarray[0][offset+nq-1];
        u_l = inarray[1][offset+nq-1];
                        
        // Get the values at initial state u_0, c_0
        u_0 = 0.0; //for all vessels start from initial condition 0
        c_0 = sqrt(beta[omega][nq-1]/(2*m_rho))*sqrt(sqrt(A_0[omega][nq-1])); 	
                        
        // Calculate the boundary values
        A_r = A_l;
        c_l = sqrt(beta[omega][nq-1]/(2*m_rho))*sqrt(sqrt(A_l));
        u_r = (1-RT)*((u_l-u_0) + 4*(c_l-c_0)) - u_l;
                        
        // Store the new values in the boundary condition
        (vessel[0]->UpdateBndCondExpansion(n))->UpdatePhys()[0] = A_r;
        (vessel[1]->UpdateBndCondExpansion(n))->UpdatePhys()[0] = u_r;
    }

}
