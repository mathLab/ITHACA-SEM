///////////////////////////////////////////////////////////////////////////////
//
// File ShallowWaterSystem.h
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
// Description: Generic timestepping for shallow water solvers
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_SHALLOWWATERSOLVER_EQUATIONSYSTEMS_SHALLOWWATERSYSTEM_H
#define NEKTAR_SOLVERS_SHALLOWWATERSOLVER_EQUATIONSYSTEMS_SHALLOWWATERSYSTEM_H

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <SolverUtils/EquationSystem.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{
  
    enum UpwindType
  {           
    eNotSet,             ///< flux not defined
    eAverage,            ///< averaged (or centred) flux
    eHLL,                ///< Harten-Lax-Leer flux
    eHLLC,               ///< Harten-Lax-Leer Contact wave flux
    SIZE_UpwindType      ///< Length of enum list
  };
  
  const char* const UpwindTypeMap[] =
    {
      "NoSet",
      "Average",
      "HLL",
      "HLLC",
    };

    /// Base class for unsteady solvers.
    class ShallowWaterSystem : public EquationSystem
    {
    public:
        /// Destructor
        virtual ~ShallowWaterSystem();

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
	/// Acceleration of gravity 
	NekDouble                                       m_g;
	/// Still water depth
	Array<OneD, NekDouble>                          m_depth;
	/// Coriolis force     
	Array<OneD, NekDouble>                          m_coriolis;
	
        /// Initialises UnsteadySystem class members.
        ShallowWaterSystem(const LibUtilities::SessionReaderSharedPtr& pSession);

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

     
        /// Evaulate flux = m_fields*ivel for i th component of Vu for
        /// direction j
        virtual void v_GetFluxVector(const int i, const int j,
                    Array<OneD, Array<OneD, NekDouble> > &physfield,
                    Array<OneD, Array<OneD, NekDouble> > &flux);
	
	void PrimitiveToConservative()
	{
	  v_PrimitiveToConservative();
	}
	virtual void v_PrimitiveToConservative();

	void ConservativeToPrimitive()
	{
	  v_ConservativeToPrimitive();
	}
	virtual void v_ConservativeToPrimitive();
	
	
    private:
	/// 
	void EvaluateWaterDepth(void);
	
	///
	void EvaluateCoriolis(void);

    };
}

#endif
