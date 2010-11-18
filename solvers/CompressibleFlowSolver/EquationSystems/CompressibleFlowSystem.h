///////////////////////////////////////////////////////////////////////////////
//
// File CompressibleFlowSystem.h
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
// Description: Generic timestepping for compressible flow solvers
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_COMPRESSIBLEFLOWSYSTEM_H
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_COMPRESSIBLEFLOWSYSTEM_H

#include <ADRSolver/EquationSystem.h>

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
  
  
  class CompressibleFlowSystem: public EquationSystem
  {
  public:           
    /// Destructor
    virtual ~CompressibleFlowSystem();
    
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
    
    /// Initialises UnsteadySystem class members.
    CompressibleFlowSystem(SessionReaderSharedPtr& pSession);
    
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

    virtual void v_SetInitialConditions(NekDouble initialtime = 0.0,
					bool dumpInitialConditions = true)
    {
      ASSERTL0(false, "This function is not implemented for this equation.");
    }
    
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

  };

}

#endif


    /* /\** */
/*      * Constructor. */
/*      * /param  */
/*      *  */
/*      * */
/*      **\/ */
/*     EulerEquations(string &fileStringName); */
    
/*     void GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> >&physfield,  */
/* 		       Array<OneD, Array<OneD, NekDouble> >&flux); */
    
/*     void NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, */
/* 		       Array<OneD, Array<OneD, NekDouble> > &numfluxX, */
/* 		       Array<OneD, Array<OneD, NekDouble> > &numfluxY); */
    
/*     void ODElhs(const Array<OneD, const  Array<OneD, NekDouble> >&inarray,  */
/* 		Array<OneD,        Array<OneD, NekDouble> >&outarray,  */
/* 		const NekDouble time); */
    
/*     void ODElhsSolve(const Array<OneD, const  Array<OneD, NekDouble> >&inarray,  */
/* 		     Array<OneD,        Array<OneD, NekDouble> >&outarray,  */
/* 		     const NekDouble time); */
    
/*     void ODErhs(const Array<OneD, const  Array<OneD, NekDouble> >&inarray,  */
/* 		          Array<OneD,        Array<OneD, NekDouble> >&outarray,  */
/* 		const NekDouble time); */
    
/*     void ODEdirkSolve(const Array<OneD, const  Array<OneD, NekDouble> >&inarray,  */
/* 		      Array<OneD,        Array<OneD, NekDouble> >&outarray,  */
/* 		      const NekDouble lambda, */
/* 		      const NekDouble time); */
    
/*     void ExplicitlyIntegrateAdvection(int nsteps); */
    
    
/*     void GeneralTimeIntegration(int nsteps,  */
/* 				LibUtilities::TimeIntegrationMethod IntMethod, */
/* 				LibUtilities::TimeIntegrationSchemeOperators ode); */
    
/*     void Summary(std::ostream &out); */
    
/*     void SetIsenTropicVortex(void); */
    
/*     void GetExactIsenTropicVortex(Array<OneD, NekDouble> &outarray, int field); */
    
/*     void SetInitialRinglebFlow(void); */
/*     void SetBoundaryRinglebFlow(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray); */
/*     void GetExactRinglebFlow(Array<OneD, NekDouble> &outarray, int field); */
    
/*     Array<OneD, MultiRegions::ExpListSharedPtr> Derived_field; */

/*     UpwindType GetUpwindType(void) */
/*     { */
/*       return m_upwindType; */
/*     }; */
    
/*     ProblemType m_problemType;   ///< euler problem selector */
    
/*     ProblemType GetProblemType(void) */
/*     { */
/*       return m_problemType; */
/*     }; */

/*   protected: */

/*   private:  */

/*     UpwindType m_upwindType;     ///< numerical upwind flux selector */

/*     int m_infosteps;  ///< dump info to stdout at steps time */

/*     int m_timemax;    ///< period of the simulation */
    
/*     int m_checktime;  ///< check time for the dump file */

/*     NekDouble m_gamma;  */
    
/*     void GetPressure(Array<OneD, Array<OneD, NekDouble> > &physfield, */
/* 		     Array<OneD, NekDouble> &pressure); */
    
/*     void GetSoundSpeed(Array<OneD, Array<OneD, NekDouble> > &physfield, */
/* 		       Array<OneD, NekDouble> &pressure, */
/* 		       Array<OneD, NekDouble> &soundspeed); */
    
/*     void GetMach(Array<OneD, Array<OneD, NekDouble> > &physfield, */
/* 		 Array<OneD, NekDouble> &soundspeed, */
/* 		 Array<OneD, NekDouble> &mach); */
    
/*     void SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble> > &physarray, NekDouble time); */
    
/*     void WallBoundary(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray); */

/*     void SymmetryBoundary(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray); */
    
/*     void RiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL, */
/* 			  NekDouble rhoR, NekDouble rhouR, NekDouble rhovR, NekDouble ER, */
/* 			  NekDouble &rhoflux, NekDouble &rhouflux, NekDouble &rhovflux, NekDouble &Eflux); */
    
/*     void HLLRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL, */
/* 			  NekDouble rhoR, NekDouble rhouR, NekDouble rhovR, NekDouble ER, */
/* 			  NekDouble &rhoflux, NekDouble &rhouflux, NekDouble &rhovflux, NekDouble &Eflux); */
    
/*     void HLLCRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL, */
/* 			   NekDouble rhoR, NekDouble rhouR, NekDouble rhovR, NekDouble ER, */
/* 			   NekDouble &rhoflux, NekDouble &rhouflux, NekDouble &rhovflux, NekDouble &Eflux); */
    
/*     void LFRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL, */
/* 			 NekDouble rhoR, NekDouble rhouR, NekDouble rhovR, NekDouble ER, */
/* 			 NekDouble &rhoflux, NekDouble &rhouflux, NekDouble &rhovflux, NekDouble &Eflux); */
	    
/*     void AverageRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL, */
/* 			      NekDouble rhoR, NekDouble rhouR, NekDouble rhovR, NekDouble ER, */
/* 			      NekDouble &rhoflux, NekDouble &rhouflux, NekDouble &rhovflux, NekDouble &Eflux); */
    
/*     void ExactRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL, */
/* 			    NekDouble rhoR, NekDouble rhouR, NekDouble rhovR, NekDouble ER, */
/* 			    NekDouble &rhoflux, NekDouble &rhouflux, NekDouble &rhovflux, NekDouble &Eflux); */

/*     void AUSMRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL, */
/* 			   NekDouble rhoR, NekDouble rhouR, NekDouble rhovR, NekDouble ER, */
/* 			   NekDouble &rhoflux, NekDouble &rhouflux, NekDouble &rhovflux, NekDouble &Eflux); */
    
/*     NekDouble M1Function(int A, NekDouble M); */
/*     NekDouble M2Function(int A, NekDouble M); */
/*     NekDouble M4Function(int A, NekDouble beta, NekDouble M); */
/*     NekDouble P5Function(int A, NekDouble alpha, NekDouble M); */
    
/*     void ExtractWall(int nchk, Array<OneD, Array<OneD, NekDouble> > &physarray); */

/*     // CFL checking functions */
/*     void GetMinLength(Array<OneD, NekDouble> &MinLength); */
/*     void GetTimeStep(const NekDouble CFL,Array<OneD, NekDouble> &MinLengt,Array<OneD, Array<OneD, NekDouble> > inarray,NekDouble &TimeStep); */
    
/*     inline NekDouble GetCFLNumber(int n) */
/*     { */
/*       NekDouble CFLDG[21] = {2,6,11.8424,19.1569,27.8419,37.8247,49.0518,61.4815,75.0797,89.8181,105.67,122.62,140.64,159.73,179.85,201.01,223.18,246.36,270.53,295.69,321.83}; //CFLDG 1D [0-20] */
      
/*       if(n<=20) */
/* 	return CFLDG[n]; */
/*       else */
/* 	ASSERTL0(false,"illegal modes dimension (CFL DG)"); */
/*     } */
    
/*     virtual void v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield,  */
/* 				 Array<OneD, Array<OneD, NekDouble> > &flux)  */
/*     {  */
/*       GetFluxVector(i,physfield,flux); */
/*     } */
    
/*     virtual void v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,  */
/* 				 Array<OneD, Array<OneD, NekDouble> > &numfluxX,  */
/* 				 Array<OneD, Array<OneD, NekDouble> > &numfluxY ) */
/*     { */
/*       NumericalFlux(physfield, numfluxX, numfluxY);  */
/*     } */
    
/*   }; */
  
/*   typedef boost::shared_ptr<EulerEquations> EulerEquationsSharedPtr; */
  
/* } //end of namespace */

/* #endif //NEKTAR_SOLVERS_EULEREQUATIONS_EULEREQUATIONS_H */

/* /\** */
/* * $Log: EulerEquations.h,v $ */
/* * Revision 1.7  2009/09/14 16:09:28  cbiotto */
/* * *** empty log message *** */
/* * */
/* * Revision 1.6  2009/08/20 10:27:52  cbiotto */
/* * Subsonic and smooth supersonic Euler. Adding numerical fluxes, boundary conditions, */
/* * initial conditions, CFL check. */
/* * */
/* * Revision 1.5  2009/07/17 09:39:51  cbiotto */
/* * Ringleb flow and CFL check */
/* * */
/* * Revision 1.4  2009/06/29 07:47:33  claes */
/* * bug fix in WallBoundary */
/* * */
/* * Revision 1.3  2009/03/17 12:32:01  claes */
/* * Updates to get the EulerSolver to work with the latest version of the TimeIntegrationScheme */
/* * */
/* * Revision 1.2  2009/02/02 16:10:16  claes */
/* * Update to make SWE, Euler and Boussinesq solvers up to date with the time integrator scheme. Linear and classical Boussinsq solver working */
/* * */
/* * Revision 1.1  2009/01/13 10:59:32  pvos */
/* * added new solvers file */
/* * */
/* * Revision 1.1  2008/11/17 08:42:06  claes */
/* * Initial commit of restructured Euler Solver */
/* * */
/* * Revision 1.1  2008/08/22 09:48:23  pvos */
/* * Added Claes' AdvectionDiffusionReaction, ShallowWater and Euler solver */
/* * */
/* **\/ */
