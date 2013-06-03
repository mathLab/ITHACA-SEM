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
#include <CardiacEPSolver/Filters/FilterCheckpointCellModel.h>

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
            "Bidomain model of cardiac electrophysiology with 3D diffusion.");


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
        m_session->LoadParameter("Chi",       m_chi);
        m_session->LoadParameter("Cm",        m_capMembrane);

        std::string vCellModel;
        m_session->LoadSolverInfo("CELLMODEL", vCellModel, "");

        ASSERTL0(vCellModel != "", "Cell Model not specified.");

        m_cell = GetCellModelFactory().CreateInstance(vCellModel, m_session, m_fields[0]);
        m_intVariables.push_back(0);
        m_intVariables.push_back(1);       

        // Load variable coefficients
        StdRegions::VarCoeffType varCoeffEnum[3] = {
                StdRegions::eVarCoeffD00,
                StdRegions::eVarCoeffD11,
                StdRegions::eVarCoeffD22
        };
        std::string varName[3] = {
                "AnisotropicConductivityX", 
                "AnisotropicConductivityY", 
                "AnisotropicConductivityZ"
        };


	if (m_session->DefinesFunction("IntracellularConductivity") && m_session->DefinesFunction("ExtracellularConductivity"))
	{
            for (int i = 0; i < m_spacedim; ++i)
            	{
                    int nq = m_fields[0]->GetNpoints();
                    Array<OneD,NekDouble> x0(nq);
                    Array<OneD,NekDouble> x1(nq);
                    Array<OneD,NekDouble> x2(nq);

                    // get the coordinates
                    m_fields[0]->GetCoords(x0,x1,x2);
		    tmp1 = Array<OneD, const Array<OneD, NekDouble> >(nq);
		    tmp2 = Array<OneD, const Array<OneD, NekDouble> >(nq);
		    tmp3 = Array<OneD, const Array<OneD, NekDouble> >(nq);
  		    tmp1[i] = Array<OneD, NekDouble>(nq);
  		    tmp2[i] = Array<OneD, NekDouble>(nq);
  		    tmp3[i] = Array<OneD, NekDouble>(nq);

                    LibUtilities::EquationSharedPtr ifunc1
                            = m_session->GetFunction("IntracellularConductivity", varName[i]);
                    LibUtilities::EquationSharedPtr ifunc2
                            = m_session->GetFunction("ExtracellularConductivity", varName[i]);
                    for(int j = 0; j < nq; j++)
                    {
                        tmp1[i][j] = ifunc1->Evaluate(x0[j],x1[j],x2[j],0.0);
                        tmp2[i][j] = ifunc2->Evaluate(x0[j],x1[j],x2[j],0.0);
                    }
		    Vmath::Vadd(nq, tmp1[i], 1, tmp2[i], 1, tmp3[i], 1);
                    m_vardiffi[varCoeffEnum[i]] = tmp1[i];
                    m_vardiffie[varCoeffEnum[i]] = tmp3[i];
            }
        } 

        
        if (m_session->DefinesParameter("StimulusDuration"))
        {
            ASSERTL0(m_session->DefinesFunction("Stimulus", "u"),
                    "Stimulus function not defined.");
            m_session->LoadParameter("StimulusDuration", m_stimDuration);
        }
        else
        {
            m_stimDuration = 0;
        }


        // Search through the loaded filters and pass the cell model to any
        // CheckpointCellModel filters loaded.
        int k = 0;
        const LibUtilities::FilterMap& f = m_session->GetFilters();
        LibUtilities::FilterMap::const_iterator x;
        for (x = f.begin(); x != f.end(); ++x, ++k)
        {
            if (x->first == "CheckpointCellModel")
            {
                boost::shared_ptr<FilterCheckpointCellModel> c
                    = boost::dynamic_pointer_cast<FilterCheckpointCellModel>(
                                                                m_filters[k]);
                c->SetCellModel(m_cell);
            }
        }

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
        int nq          = m_fields[0]->GetNpoints();

        Array<OneD, NekDouble> grad0(nq), grad1(nq), grad2(nq), grad(nq);
        Array<OneD, NekDouble> ggrad0(nq), ggrad1(nq), ggrad2(nq), ggrad(nq), temp(nq);

        // We solve ( \sigma\nabla^2 - HHlambda ) Y[i] = rhs [i]
        // inarray = input: \hat{rhs} -> output: \hat{Y}
        // outarray = output: nabla^2 \hat{Y}
        // where \hat = modal coeffs
        for (int i = 0; i < nvariables; ++i)
        {
            // Only apply diffusion to first variable.
            if (i > 1) {
                Vmath::Vcopy(nq, &inarray[i][0], 1, &outarray[i][0], 1);
                continue;
            }
            if (i == 0) {
                StdRegions::ConstFactorMap factors;
                factors[StdRegions::eFactorLambda] = (1.0/lambda)*(m_capMembrane*m_chi);
                if (m_spacedim==1) {
                // Take first partial derivative
                m_fields[i]->PhysDeriv(inarray[1],ggrad0);
                // Take second partial derivative
                m_fields[i]->PhysDeriv(0,ggrad0,ggrad0);
                // Multiply by Intracellular-Conductivity
		if (m_session->DefinesFunction("IntracellularConductivity") && m_session->DefinesFunction("ExtracellularConductivity"))
		{
                Vmath::Smul(nq, m_session->GetParameter("sigmaix"), ggrad0, 1, ggrad0, 1);
                }
                // Add partial derivatives together
                Vmath::Vcopy(nq, ggrad0, 1, ggrad, 1);
                Vmath::Smul(nq, -1.0, ggrad, 1, ggrad, 1);
                // Multiply 1.0/timestep/lambda
                Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], inarray[i], 1, temp, 1);
                Vmath::Vadd(nq, ggrad, 1, temp, 1, m_fields[i]->UpdatePhys(), 1);
                // Solve a system of equations with Helmholtz solver and transform
                // back into physical space.
                m_fields[i]->HelmSolve(m_fields[i]->GetPhys(), m_fields[i]->UpdateCoeffs(),NullFlagList,factors);
                m_fields[i]->BwdTrans( m_fields[i]->GetCoeffs(), m_fields[i]->UpdatePhys());
                m_fields[i]->SetPhysState(true);
                // Copy the solution vector (required as m_fields must be set).
                outarray[i] = m_fields[i]->GetPhys();
                }
                
                if (m_spacedim==2) {
                // Take first partial derivative
                m_fields[i]->PhysDeriv(inarray[1],ggrad0,ggrad1);
                // Take second partial derivative
                m_fields[i]->PhysDeriv(0,ggrad0,ggrad0);
                m_fields[i]->PhysDeriv(1,ggrad1,ggrad1);
                // Multiply by Intracellular-Conductivity 
		if (m_session->DefinesFunction("IntracellularConductivity") && m_session->DefinesFunction("ExtracellularConductivity"))
		{                 
                Vmath::Smul(nq, m_session->GetParameter("sigmaix"), ggrad0, 1, ggrad0, 1);   
                Vmath::Smul(nq, m_session->GetParameter("sigmaiy"), ggrad1, 1, ggrad1, 1); 
                }
                // Add partial derivatives together
                Vmath::Vadd(nq, ggrad0, 1, ggrad1, 1, ggrad, 1);
                Vmath::Smul(nq, -1.0, ggrad, 1, ggrad, 1); 
                // Multiply 1.0/timestep/lambda
                Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], inarray[i], 1, temp, 1);
                Vmath::Vadd(nq, ggrad, 1, temp, 1, m_fields[i]->UpdatePhys(), 1);
                // Solve a system of equations with Helmholtz solver and transform
                // back into physical space.
                m_fields[i]->HelmSolve(m_fields[i]->GetPhys(), m_fields[i]->UpdateCoeffs(),NullFlagList,factors,m_vardiffi);
                m_fields[i]->BwdTrans( m_fields[i]->GetCoeffs(), m_fields[i]->UpdatePhys());
                m_fields[i]->SetPhysState(true);
                // Copy the solution vector (required as m_fields must be set).
                outarray[i] = m_fields[i]->GetPhys();                                                                          
                }
                
                if (m_spacedim==3) {
                // Take first partial derivative
                m_fields[i]->PhysDeriv(inarray[1],ggrad0,ggrad1,ggrad2);
                // Take second partial derivative
                m_fields[i]->PhysDeriv(0,ggrad0,ggrad0);
                m_fields[i]->PhysDeriv(1,ggrad1,ggrad1);                
                m_fields[i]->PhysDeriv(2,ggrad2,ggrad2); 
                // Multiply by Intracellular-Conductivity 
		if (m_session->DefinesFunction("IntracellularConductivity") && m_session->DefinesFunction("ExtracellularConductivity"))
		{                 
                Vmath::Smul(nq, m_session->GetParameter("sigmaix"), ggrad0, 1, ggrad0, 1);  
                Vmath::Smul(nq, m_session->GetParameter("sigmaiy"), ggrad1, 1, ggrad1, 1);
                Vmath::Smul(nq, m_session->GetParameter("sigmaiz"), ggrad2, 1, ggrad2, 1);
                }
                // Add partial derivatives together
                Vmath::Vadd(nq, ggrad0, 1, ggrad1, 1, ggrad, 1);
                Vmath::Vadd(nq, ggrad2, 1, ggrad, 1, ggrad, 1);
                Vmath::Smul(nq, -1.0, ggrad, 1, ggrad, 1);
                // Multiply 1.0/timestep/lambda
                Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], inarray[i], 1, temp, 1);
                Vmath::Vadd(nq, ggrad, 1, temp, 1, m_fields[i]->UpdatePhys(), 1);
                // Solve a system of equations with Helmholtz solver and transform
                // back into physical space.
                m_fields[i]->HelmSolve(m_fields[i]->GetPhys(), m_fields[i]->UpdateCoeffs(),NullFlagList,factors,m_vardiffi);
                m_fields[i]->BwdTrans( m_fields[i]->GetCoeffs(), m_fields[i]->UpdatePhys());
                m_fields[i]->SetPhysState(true);
                // Copy the solution vector (required as m_fields must be set).
                outarray[i] = m_fields[i]->GetPhys();                                         
                }                                

            }
            if (i == 1) {
                StdRegions::ConstFactorMap factors;
                factors[StdRegions::eFactorLambda] = 0.0;
                if (m_spacedim==1) {
                // Take first partial derivative
                m_fields[i]->PhysDeriv(m_fields[0]->UpdatePhys(),grad0);
                // Take second derivative
                m_fields[i]->PhysDeriv(0,grad0,grad0);  
                // Multiply by Intracellular-Conductivity
		if (m_session->DefinesFunction("IntracellularConductivity") && m_session->DefinesFunction("ExtracellularConductivity"))
		{                
                Vmath::Smul(nq, m_session->GetParameter("sigmaix"), grad0, 1, grad0, 1);   
                }
                // and sum terms
                Vmath::Vcopy(nq, grad0, 1, grad, 1);
                Vmath::Smul(nq, (-1.0*m_session->GetParameter("sigmaix"))/(m_session->GetParameter("sigmaix")+m_session->GetParameter("sigmaix")), grad, 1, grad, 1);  
                // Now solve Poisson problem for \phi_e
                m_fields[i]->SetPhys(grad);
                m_fields[i]->HelmSolve(m_fields[i]->GetPhys(), m_fields[i]->UpdateCoeffs(), NullFlagList, factors);
                m_fields[i]->BwdTrans( m_fields[i]->GetCoeffs(), m_fields[i]->UpdatePhys());
                m_fields[i]->SetPhysState(true);
                // Copy the solution vector (required as m_fields must be set).
                outarray[i] = m_fields[i]->GetPhys();                                                         
                }
                
                if (m_spacedim==2) {
                // Take first partial derivative
                m_fields[i]->PhysDeriv(m_fields[0]->UpdatePhys(),grad0,grad1);
                // Take second derivative
                m_fields[i]->PhysDeriv(0,grad0,grad0);
                m_fields[i]->PhysDeriv(1,grad1,grad1);  
                // Multiply by Intracellular-Conductivity
		if (m_session->DefinesFunction("IntracellularConductivity") && m_session->DefinesFunction("ExtracellularConductivity"))
		{                
                Vmath::Smul(nq, m_session->GetParameter("sigmaix"), grad0, 1, grad0, 1);
                Vmath::Smul(nq, m_session->GetParameter("sigmaiy"), grad1, 1, grad1, 1); 
                }
                // and sum terms
                Vmath::Vadd(nq, grad0, 1, grad1, 1, grad, 1);
                Vmath::Smul(nq, -1.0, grad, 1, grad, 1);   
                // Now solve Poisson problem for \phi_e
                m_fields[i]->SetPhys(grad);
                m_fields[i]->HelmSolve(m_fields[i]->GetPhys(), m_fields[i]->UpdateCoeffs(), NullFlagList, factors, m_vardiffie);
                m_fields[i]->BwdTrans( m_fields[i]->GetCoeffs(), m_fields[i]->UpdatePhys());
                m_fields[i]->SetPhysState(true);
                // Copy the solution vector (required as m_fields must be set).
                outarray[i] = m_fields[i]->GetPhys();                                                        
                }
                
                if (m_spacedim==3) {
                // Take first partial derivative
                m_fields[i]->PhysDeriv(m_fields[0]->UpdatePhys(),grad0,grad1,grad2);
                // Take second derivative
                m_fields[i]->PhysDeriv(0,grad0,grad0);
                m_fields[i]->PhysDeriv(1,grad1,grad1);
                m_fields[i]->PhysDeriv(2,grad2,grad2); 
                // Multiply by Intracellular-Conductivity
		if (m_session->DefinesFunction("IntracellularConductivity") && m_session->DefinesFunction("ExtracellularConductivity"))
		{                
                Vmath::Smul(nq, m_session->GetParameter("sigmaix"), grad0, 1, grad0, 1);
                Vmath::Smul(nq, m_session->GetParameter("sigmaiy"), grad1, 1, grad1, 1);
                Vmath::Smul(nq, m_session->GetParameter("sigmaiz"), grad2, 1, grad2, 1);       
                }
                // and sum terms
                Vmath::Vadd(nq, grad0, 1, grad1, 1, grad, 1);
                Vmath::Vadd(nq, grad2, 1, grad, 1, grad, 1);
                Vmath::Smul(nq, -1.0, grad, 1, grad, 1);  
                // Now solve Poisson problem for \phi_e
                m_fields[i]->SetPhys(grad);
                m_fields[i]->HelmSolve(m_fields[i]->GetPhys(), m_fields[i]->UpdateCoeffs(), NullFlagList, factors, m_vardiffie);
                m_fields[i]->BwdTrans( m_fields[i]->GetCoeffs(), m_fields[i]->UpdatePhys());
                m_fields[i]->SetPhysState(true);
                // Copy the solution vector (required as m_fields must be set).
                outarray[i] = m_fields[i]->GetPhys();                                                      
                }
                                                
            }
        }
    }


    void Bidomain::DoOdeRhs(
            const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
            Array<OneD,        Array<OneD, NekDouble> >&outarray,
            const NekDouble time)
    {
        int nq = m_fields[0]->GetNpoints();
        m_cell->TimeIntegrate(inarray, outarray, time);
        if (m_stimDuration > 0 && time < m_stimDuration)
        {
            Array<OneD,NekDouble> x0(nq);
            Array<OneD,NekDouble> x1(nq);
            Array<OneD,NekDouble> x2(nq);
            Array<OneD,NekDouble> result(nq);

            // get the coordinates
            m_fields[0]->GetCoords(x0,x1,x2);

            LibUtilities::EquationSharedPtr ifunc
                    = m_session->GetFunction("Stimulus", "u");
            ifunc->Evaluate(x0,x1,x2,time, result);

            Vmath::Vadd(nq, outarray[0], 1, result, 1, outarray[0], 1);
        }        
        Vmath::Smul(nq, 1.0/m_capMembrane, outarray[0], 1, outarray[0], 1);
    }


    void Bidomain::v_SetInitialConditions(NekDouble initialtime,
            bool dumpInitialConditions)
    {
        EquationSystem::v_SetInitialConditions(initialtime, dumpInitialConditions);
        m_cell->Initialise();
    }

    /**
     *
     */
    void Bidomain::v_PrintSummary(std::ostream &out)
    {
        UnsteadySystem::v_PrintSummary(out);
        out << "\tChi       	: " << m_chi << endl;
        out << "\tCm       	: " << m_capMembrane << endl;
        if (m_session->DefinesFunction("IntracellularConductivity", "AnisotropicConductivityX") &&
            m_session->GetFunctionType("IntracellularConductivity", "AnisotropicConductivityX") == LibUtilities::eFunctionTypeExpression)
        {
            out << "\tIntra-Diffusivity-x   : "
                << m_session->GetFunction("IntracellularConductivity", "AnisotropicConductivityX")->GetExpression()
                << endl;
        }
        if (m_session->DefinesFunction("IntracellularConductivity", "AnisotropicConductivityY") &&
            m_session->GetFunctionType("IntracellularConductivity", "AnisotropicConductivityY") == LibUtilities::eFunctionTypeExpression)
        {
            out << "\tIntra-Diffusivity-y   : "
                << m_session->GetFunction("IntracellularConductivity", "AnisotropicConductivityY")->GetExpression()
                << endl;
        }
        if (m_session->DefinesFunction("IntracellularConductivity", "AnisotropicConductivityZ") &&
            m_session->GetFunctionType("IntracellularConductivity", "AnisotropicConductivityZ") == LibUtilities::eFunctionTypeExpression)
        {
            out << "\tIntra-Diffusivity-z   : "
                << m_session->GetFunction("IntracellularConductivity", "AnisotropicConductivityZ")->GetExpression()
                << endl;
        }
        if (m_session->DefinesFunction("ExtracellularConductivity", "AnisotropicConductivityX") &&
            m_session->GetFunctionType("ExtracellularConductivity", "AnisotropicConductivityX") == LibUtilities::eFunctionTypeExpression)
        {
            out << "\tExtra-Diffusivity-x   : "
                << m_session->GetFunction("ExtracellularConductivity", "AnisotropicConductivityX")->GetExpression()
                << endl;
        }
        if (m_session->DefinesFunction("ExtracellularConductivity", "AnisotropicConductivityY") &&
            m_session->GetFunctionType("ExtracellularConductivity", "AnisotropicConductivityY") == LibUtilities::eFunctionTypeExpression)
        {
            out << "\tExtra-Diffusivity-y   : "
                << m_session->GetFunction("ExtracellularConductivity", "AnisotropicConductivityY")->GetExpression()
                << endl;
        }
        if (m_session->DefinesFunction("ExtracellularConductivity", "AnisotropicConductivityZ") &&
            m_session->GetFunctionType("ExtracellularConductivity", "AnisotropicConductivityZ") == LibUtilities::eFunctionTypeExpression)
        {
            out << "\tExtra-Diffusivity-z   : "
                << m_session->GetFunction("ExtracellularConductivity", "AnisotropicConductivityZ")->GetExpression()
                << endl;
        }        
        m_cell->PrintSummary(out);
    }

}
