///////////////////////////////////////////////////////////////////////////////
//
// File Bidomain.cpp
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
// Description: Bidomain cardiac electrophysiology homogenised model.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <CardiacEPSolver/EquationSystems/Bidomain.h>

namespace Nektar
{
    /**
     * @class Bidomain
     *
     * Base model of cardiac electrophysiology of the form
     * \f{align*}{
     *     \frac{\partial u}{\partial t} = \nabla^2 u + J_{ion},
     * \f}
     * where the reaction term, \f$J_{ion}\f$ is defined by a specific cell
     * model.
     *
     * This implementation, at present, treats the reaction terms explicitly
     * and the diffusive element implicitly.
     */

    /**
     * Registers the class with the Factory.
     */
    string Bidomain::className
            = GetEquationSystemFactory().RegisterCreatorFunction(
                "Bidomain",
                Bidomain::create,
                "Bidomain model of cardiac electrophysiology.");


    /**
     *
     */
    Bidomain::Bidomain(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession)
    {
    }

    void Bidomain::v_InitObject()
    {
        UnsteadySystem::v_InitObject();
        m_session->LoadParameter("chi",       m_chi, 1400);
        m_session->LoadParameter("cm",        m_cm, 1.0);
        m_session->LoadParameter("sigmai",    m_sigmai, 1.75);
        m_session->LoadParameter("sigmae",    m_sigmae, 7.0);

        std::string vCellModel;
        m_session->LoadSolverInfo("CELLMODEL", vCellModel, "");
	
        ASSERTL0(vCellModel != "", "Cell Model not specified.");

        m_cell = GetCellModelFactory().CreateInstance(vCellModel, m_session, m_fields[0]);

        m_intVariables.push_back(0);
        m_intVariables.push_back(1);

        if (!m_explicitDiffusion)
        {
            m_ode.DefineImplicitSolve (&Bidomain::DoImplicitSolve, this);
        }
        m_ode.DefineOdeRhs(&Bidomain::DoOdeRhs, this);
    }


    /**
     *
     */
    Bidomain::~Bidomain()
    {

    }


    /**
     * @param   inarray         Input array.
     * @param   outarray        Output array.
     * @param   time            Current simulation time.
     * @param   lambda          Timestep.
     */
    void Bidomain::DoImplicitSolve(
            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                  Array<OneD, Array<OneD, NekDouble> >&outarray,
            const NekDouble time,
            const NekDouble lambda)
    {
        int nvariables  = inarray.num_elements();
        int ncoeffs     = inarray[0].num_elements();
        int nq          = m_fields[0]->GetNpoints();

        Array<OneD, NekDouble> grad0(nq), grad1(nq), grad2(nq), grad3(nq);
        Array<OneD, NekDouble> ggrad0(nq), ggrad1(nq), ggrad2(nq), ggrad3(nq), ttmp(nq);
        NekDouble ratio1=(-1.0*m_sigmai)/(m_sigmai+m_sigmae);
        NekDouble ratio2=(m_chi*m_cm)/(m_sigmai);
	
        // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
        // inarray = input: \hat{rhs} -> output: \hat{Y}
        // outarray = output: nabla^2 \hat{Y}
        // where \hat = modal coeffs
        for (int i = 0; i < nvariables; ++i)
        {
            if (i == 0) {
                StdRegions::ConstFactorMap factors;
                factors[StdRegions::eFactorLambda] = (1.0/lambda)*ratio2;
                m_fields[i]->PhysDeriv(inarray[1],ggrad0,ggrad1,ggrad2);
                // Take second derivative
                m_fields[i]->PhysDeriv(0,ggrad0,ggrad0);
                m_fields[i]->PhysDeriv(1,ggrad1,ggrad1);
                m_fields[i]->PhysDeriv(2,ggrad2,ggrad2);
                // and sum terms
                Vmath::Vadd(nq, &ggrad0[0], 1, &ggrad1[0], 1, &ggrad3[0], 1);
                Vmath::Vadd(nq, &ggrad2[0], 1, &ggrad3[0], 1, &ggrad3[0], 1);
                Vmath::Smul(nq, -1.0, &ggrad3[0], 1, &ggrad3[0], 1);

                // Multiply 1.0/timestep/lambda
                Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], &inarray[i][0], 1, &ttmp[0], 1);
                Vmath::Vadd(nq, ggrad3, 1, ttmp, 1, m_fields[i]->UpdatePhys(), 1);
                // Solve a system of equations with Helmholtz solver and transform
                // back into physical space.
                m_fields[i]->HelmSolve(m_fields[i]->GetPhys(), m_fields[i]->UpdateCoeffs(),NullFlagList,factors);
                m_fields[i]->BwdTrans( m_fields[i]->GetCoeffs(), m_fields[i]->UpdatePhys());
                m_fields[i]->SetPhysState(true);
                // Copy the solution vector (required as m_fields must be set).
                outarray[i] = m_fields[i]->UpdatePhys();
            }
            else if (i == 1) {
                StdRegions::ConstFactorMap factors;
                factors[StdRegions::eFactorLambda] = 0.0;
                m_fields[i]->PhysDeriv(m_fields[0]->UpdatePhys(),grad0,grad1,grad2);
                // Take second derivative
                m_fields[i]->PhysDeriv(0,grad0,grad0);
                m_fields[i]->PhysDeriv(1,grad1,grad1);
                m_fields[i]->PhysDeriv(2,grad2,grad2);
                // and sum terms
                Vmath::Vadd(nq, &grad0[0], 1, &grad1[0], 1, &grad3[0], 1);
                Vmath::Vadd(nq, &grad2[0], 1, &grad3[0], 1, &grad3[0], 1);
                Vmath::Smul(nq, ratio1, grad3, 1, grad3, 1);
                // Now solve Poisson problem for \phi_e
                m_fields[i]->SetPhys(grad3);
                m_fields[i]->HelmSolve(m_fields[i]->GetPhys(), m_fields[i]->UpdateCoeffs(),NullFlagList,factors);
                m_fields[i]->BwdTrans( m_fields[i]->GetCoeffs(), m_fields[i]->UpdatePhys());
                m_fields[i]->SetPhysState(true);
                // Copy the solution vector (required as m_fields must be set).
                outarray[i] = m_fields[i]->UpdatePhys();
            }
            else {
                Vmath::Vcopy(nq, &inarray[i][0], 1, &outarray[i][0], 1);
            }
        }
    }


    void Bidomain::DoOdeRhs(
            const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                  Array<OneD,        Array<OneD, NekDouble> >&outarray,
            const NekDouble time)
    {
        m_cell->TimeIntegrate(inarray, outarray, time);
    }


    void Bidomain::v_SetInitialConditions(NekDouble initialtime,
                        bool dumpInitialConditions)
    {
        int nq = m_fields[0]->GetNpoints();
        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);

        m_fields[0]->GetCoords(x0,x1,x2);
        //NekDouble xmax = Vmath::Vmax(nq,x0,1);
        NekDouble unew=100.0, m_beta = 1.0, uinit=0.0, f, fd,Tol=0.00000001;
        // get the coordinates (assuming all fields have the same discretisation)vinit=0.0, 
        m_fields[0]->GetCoords(x0,x1,x2);
        // Get the initial constant u and v
        for(int j = 0; j<1000; j++)
        {
            f = uinit*uinit*uinit + 3.0*uinit + 6.0*m_beta;
            fd = 3.0*uinit*uinit + 3.0;
            unew = uinit - f/fd;
            if(abs(unew-uinit) < Tol)
            {
                break;
            }
            uinit = unew;
        }
        //vinit = uinit - (1.0/3.0)*uinit*uinit*uinit;
        m_uinit = uinit; // m_vinit = vinit;
        for(int i = 0 ; i < m_fields.num_elements(); i++)
        {
            LibUtilities::EquationSharedPtr ifunc
            = m_session->GetFunction("InitialConditions", i);
            for(int j = 0; j < nq; j++)
            {
                if( x0[j] > 3.5 )
                {
                    (m_fields[0]->UpdatePhys())[j] = m_uinit;
                }
                (m_fields[i]->UpdatePhys())[j]
                          = ifunc->Evaluate(x0[j],x1[j],x2[j],initialtime);
            }
            m_fields[i]->SetPhysState(true);
            m_fields[i]->FwdTrans_IterPerExp(m_fields[i]->GetPhys(),
                    m_fields[i]->UpdateCoeffs());
            cout << "\tField "<< m_session->GetVariable(i)
                                 <<": " << ifunc->GetExpression() << endl;
        }
        if(dumpInitialConditions)
        {
            std::string outname = m_sessionName +"_0.chk";

            // dump initial conditions to file
            WriteFld(outname);
        }
        m_cell->Initialise();
    }


    /**
     *
     */
    void Bidomain::v_PrintSummary(std::ostream &out)
    {
        UnsteadySystem::v_PrintSummary(out);
        out << "\tChi       	  : " << m_chi << endl;
        out << "\tCm       	      : " << m_cm << endl;
        out << "\tI-Conductivity  : " << m_sigmai << endl;
        out << "\tE-Conductivity  : " << m_sigmae << endl;
        m_cell->PrintSummary(out);
    }

}
