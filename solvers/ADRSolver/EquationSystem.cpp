///////////////////////////////////////////////////////////////////////////////
//
// File EquationSystem.cpp
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
// Description: Main wrapper class for Advection Diffusion Reaction Solver
//
///////////////////////////////////////////////////////////////////////////////

#include <string>
using std::string;

#include <ADRSolver/EquationSystem.h>

namespace Nektar
{
    /**
     * @class EquationSystem
     *
     * This class is a base class for all solver implementations. It provides
     * the underlying generic functionality and interface for solving equations.
     *
     * To solve a steady-state equation, create a derived class from this class
     * and reimplement the virtual functions to provide custom implementation
     * for the problem.
     *
     * To solve unsteady problems, derive from the UnsteadySystem class instead
     * which provides general time integration.
     */

    /**
     * This constructor is protected as the objects of this class are never
     * instantiated directly.
     * @param   pSession        The session reader holding problem parameters.
     */
    EquationSystem::EquationSystem(SessionReaderSharedPtr& pSession)
        : ADRBase(pSession->GetFilename(), true),
          m_session(pSession)
    {
        ZeroPhysFields();
    }


    /**
     *
     */
    EquationSystem::~EquationSystem()
    {

    }


    /**
     * This allows initialisation of the solver which cannot be completed
     * during object construction (such as setting of initial conditions).
     *
     * Public interface routine to virtual function implementation.
     */
    void EquationSystem::DoInitialise()
    {
        v_DoInitialise();
    }


    /**
     * Performs the actual solve.
     *
     * Public interface routine to virtual function implementation.
     */
    void EquationSystem::DoSolve()
    {
        v_DoSolve();
    }


    /**
     * Prints a summary of variables and problem parameters.
     *
     * Public interface routine to virtual function implementation.
     *
     * @param   out             The ostream object to write to.
     */
    void EquationSystem::PrintSummary(std::ostream &out)
    {
        out << "=======================================================================" << endl;
        out << "\tEquation Type   : " << m_session->GetSolverInfo("EQTYPE") << endl;
        ADRBase::SessionSummary(out);

        v_PrintSummary(out);

        out << "=======================================================================" << endl;
    }


    /**
     * Evaluates a physical function at each quadrature point in the domain.
     *
     * @param   pArray          The array into which to write the values.
     * @param   pEqn            The equation to evaluate.
     */
    void EquationSystem::EvaluateFunction(Array<OneD, NekDouble>& pArray,
                              SpatialDomains::ConstUserDefinedEqnShPtr pEqn)
    {
        int nq = m_fields[0]->GetNpoints();

        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);

        // get the coordinates (assuming all fields have the same
        // discretisation)
        m_fields[0]->GetCoords(x0,x1,x2);

        if (pArray.num_elements() != nq)
        {
            pArray = Array<OneD, NekDouble>(nq);
        }
        for(int i = 0; i < nq; i++)
        {
            pArray[i] = pEqn->Evaluate(x0[i],x1[i],x2[i]);
        }

    }


    /**
     * If boundary conditions are time-dependent, they will be evaluated at the
     * time specified.
     *
     * @param   time            The time at which to evaluate the BCs
     */
    void EquationSystem::SetBoundaryConditions(NekDouble time)
    {
      int nvariables = m_fields.num_elements();
      for (int i = 0; i < nvariables; ++i)
      {
          m_fields[i]->EvaluateBoundaryConditions(time);
      }
    }


    /**
     * By default, nothing needs initialising at the EquationSystem level.
     */
    void EquationSystem::v_DoInitialise()
    {

    }


    /**
     *
     */
    void EquationSystem::v_DoSolve()
    {
        ASSERTL0(false, "Solve not defined for this equation.");
    }


    /**
     * By default, there are no further parameters to display.
     */
    void EquationSystem::v_PrintSummary(std::ostream &out)
    {

    }

}
