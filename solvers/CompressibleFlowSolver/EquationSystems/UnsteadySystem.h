///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadySystem.h
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
// Description: Generic timestepping for unsteady compressible flow solvers
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_UNSTEADYSYSTEM_H
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_UNSTEADYSYSTEM_H

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <SolverUtils/EquationSystem.h>

namespace Nektar
{     

  enum UpwindType
  {           
    eNotSet,             ///< flux not defined
    eAverage,            ///< averaged (or centred) flux
    eLLF,                ///< local Lax-Friedrich flux
    eHLL,                ///< Harten-Lax-Leer flux
    eHLLC,               ///< Harten-Lax-Leer Contact corrected flux
    eAUSM,               ///< Advection-Upstream-Splitting-Method flux
    eAUSMPlus,           ///< Advection-Upstream-Splitting-Method flux with reduced pressure oscillations
    eAUSMPlusUp,         ///< Advection-Upstream-Splitting-Method flux for low Mach flows   
    eAUSMPlusUpAllSpeed, ///< Advection-Upstream-Splitting-Method flux for all Mach flows
    eExact,              ///< Exact flux
    SIZE_UpwindType      ///< Length of enum list
  };
  
  const char* const UpwindTypeMap[] =
    {
      "NoSet",
      "Average",
      "LF",
      "HLL",
      "HLLC",
      "AUSM",
      "AUSMPlus",
      "AUSMPlusUp",
      "AUSMPlusUpAllSpeed",
      "Exact"
    };
  
  
  class UnsteadySystem: public SolverUtils::EquationSystem
  {
  public:           
    /// Destructor
    virtual ~UnsteadySystem();

    /// GetTimeStep.
    NekDouble GetTimeStep(const Array<OneD, Array<OneD,NekDouble> > physarray, const Array<OneD,int> ExpOrder, const Array<OneD,NekDouble> CFLDG, NekDouble timeCFL)
    {
      NekDouble TimeStep = v_GetTimeStep(physarray, ExpOrder, CFLDG, timeCFL);
      
      return TimeStep;
    }
    
  protected:
    ///< numerical upwind flux selector
    UpwindType                                      m_upwindType;       
    /// Number of time steps between outputting status information.
    int                                             m_infosteps;
    /// The time integration method to use.
    LibUtilities::TimeIntegrationMethod             m_timeIntMethod;
    /// The time integration scheme operators to use.
    LibUtilities::TimeIntegrationSchemeOperators    m_ode;
    /// Indicates if explicit or implicit treatment of diffusion is used.
    bool                                            m_explicitDiffusion;
    /// Indicates if explicit or implicit treatment of advection is used.
    bool                                            m_explicitAdvection;
    /// Indicates if variables are primitive or conservative
    bool                                            m_primitive;
    /// Cp/Cv ratio gamma
    NekDouble m_gamma;
    /// Gas Constant
    NekDouble m_GasConstant;
    
    /// Initialises UnsteadySystem class members.
    UnsteadySystem(const LibUtilities::SessionReaderSharedPtr& pSession);
    
    virtual void v_InitObject();

    /// Solves an unsteady problem.
    virtual void v_DoSolve();
    
    /// Sets up initial conditions.
    virtual void v_DoInitialise();

    /// Print a summary of time stepping parameters.
    virtual void v_PrintSummary(std::ostream &out);
    
    ///
    virtual void v_NumericalFlux(
				 Array<OneD, Array<OneD, NekDouble> > &physfield,
				 Array<OneD, Array<OneD, NekDouble> > &numflux);
    
    ///
    virtual void v_NumericalFlux(
				 Array<OneD, Array<OneD, NekDouble> > &physfield,
				 Array<OneD, Array<OneD, NekDouble> > &numfluxX,
				 Array<OneD, Array<OneD, NekDouble> > &numfluxY );
    
    virtual void v_GetFluxVector(const int i, const int j,
				 Array<OneD, Array<OneD, NekDouble> > &physfield,
				 Array<OneD, Array<OneD, NekDouble> > &flux);

    virtual void v_SetInitialConditions(NekDouble initialtime = 0.0,
					bool dumpInitialConditions = true) = 0;

    virtual NekDouble v_GetTimeStep(const Array<OneD, Array<OneD,NekDouble> > physarray, 
			       const Array<OneD,int> ExpOrder, 
			       const Array<OneD,NekDouble> CFLDG, NekDouble timeCFL) = 0;

    virtual void v_EvaluateExactSolution(unsigned int field,
    					 Array<OneD, NekDouble> &outfield,
					 const NekDouble time)
    {

    }
    
  private:

    ///< CFL parameter selection
    NekDouble m_cfl;

    inline NekDouble GetStabilityLimit(int n)
    {
        NekDouble CFLDG[21] = {2,6,11.8424,19.1569,27.8419,37.8247,49.0518,61.4815,75.0797,89.8181,105.67,122.62,140.64,159.73,179.85,201.01,223.18,246.36,270.53,295.69,321.83}; //CFLDG 1D [0-20]

        if(n<=20)
        {
            return CFLDG[n];
        }
        else
        {
            ASSERTL0(false,"illegal modes dimension (CFL DG)");
            return 0;
        }
    }

    inline Array<OneD,NekDouble> GetStabilityLimitVector(const Array<OneD,int> &ExpOrder)
    {
        int i;
        Array<OneD,NekDouble> returnval(m_fields[0]->GetExpSize(),0.0);
        for(i =0; i<m_fields[0]->GetExpSize(); i++)
        {
            returnval[i] = GetStabilityLimit(ExpOrder[i]);
        }
        return returnval;
    }

  };

}

#endif
