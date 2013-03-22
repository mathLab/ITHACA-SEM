///////////////////////////////////////////////////////////////////////////////
//
// File: AverageSolver.cpp
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
// Description: Average Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/AverageSolver.h>

using namespace std;

namespace Nektar
{
    std::string AverageSolver::solverName =
        SolverUtils::GetRiemannSolverFactory().RegisterCreatorFunction(
            "Average",
            AverageSolver::create,
            "Average Riemann solver");

    AverageSolver::AverageSolver() : CompressibleSolver()
    {
        m_pointSolve = false;
    }

    /**
     * @brief Average Riemann solver.
     *
     * @param rhoL      Density left state.
     * @param rhoR      Density right state.  
     * @param rhouL     x-momentum component left state.  
     * @param rhouR     x-momentum component right state.  
     * @param rhovL     y-momentum component left state.  
     * @param rhovR     y-momentum component right state.  
     * @param rhowL     z-momentum component left state.  
     * @param rhowR     z-momentum component right state.
     * @param EL        Energy left state.  
     * @param ER        Energy right state. 
     * @param rhof      Computed Riemann flux for density.
     * @param rhouf     Computed Riemann flux for x-momentum component 
     * @param rhovf     Computed Riemann flux for y-momentum component 
     * @param rhowf     Computed Riemann flux for z-momentum component 
     * @param Ef        Computed Riemann flux for energy.
     */
    void AverageSolver::v_ArraySolve(
        const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
        const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
              Array<OneD,       Array<OneD, NekDouble> > &flux)
    {
        static NekDouble gamma = m_params["gamma"]();
        
        int expDim = Fwd.num_elements()-2;
        int i, j;
        
        for (j = 0; j < Fwd[0].num_elements(); ++j)
        {
            NekDouble tmp1 = 0.0, tmp2 = 0.0;
            Array<OneD, NekDouble> Ufwd(expDim);
            Array<OneD, NekDouble> Ubwd(expDim);
            
            for (i = 0; i < expDim; ++i)
            {
                Ufwd[i] = Fwd[i+1][j]/Fwd[0][j];
                Ubwd[i] = Bwd[i+1][j]/Bwd[0][j];
                tmp1   += Ufwd[i]*Fwd[i+1][j];
                tmp2   += Ubwd[i]*Bwd[i+1][j];
            }
            
            NekDouble Pfwd = (gamma - 1.0) * (Fwd[expDim+1][j] - 0.5 * tmp1);
            NekDouble Pbwd = (gamma - 1.0) * (Bwd[expDim+1][j] - 0.5 * tmp2);
            
            // Compute the average flux
            flux[0][j] = 0.5 * (Fwd[1][j] + Bwd[1][j]);
            flux[expDim+1][j] = 0.5 * (Ufwd[0] * (Fwd[expDim+1][j] + Pfwd) + 
                                       Ubwd[0] * (Bwd[expDim+1][j] + Pbwd));
            
            for (i = 0; i < expDim; ++i)
            {
                flux[i+1][j] = 0.5 * (Fwd[0][j] * Ufwd[0] * Ufwd[i] + 
                                      Bwd[0][j] * Ubwd[0] * Ubwd[i]);
            }

            // Add in pressure contribution to u field
            flux[1][j] += 0.5 * (Pfwd + Pbwd);
        }
    }
}
