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
// Description: Auxiliary functions for the compressible flow system
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_COMPRESSIBLEFLOWSYSTEM_H
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_COMPRESSIBLEFLOWSYSTEM_H

#include <CompressibleFlowSolver/EquationSystems/UnsteadySystem.h>

namespace Nektar
{  
  /**
   * 
   * 
   **/
  class CompressibleFlowSystem: public UnsteadySystem
  {
  public:

      friend class MemoryManager<CompressibleFlowSystem>;

    /// Creates an instance of this class
    static EquationSystemSharedPtr create(
            LibUtilities::CommSharedPtr& pComm,
            LibUtilities::SessionReaderSharedPtr& pSession)
    {
      return MemoryManager<CompressibleFlowSystem>::AllocateSharedPtr(pComm, pSession);
    }
    /// Name of class
    static std::string className;
    
    virtual ~CompressibleFlowSystem();

  protected:

    CompressibleFlowSystem(
            LibUtilities::CommSharedPtr& pComm,
            LibUtilities::SessionReaderSharedPtr& pSession);

    virtual void v_InitObject();
    
    /// Print a summary of time stepping parameters.
    virtual void v_PrintSummary(std::ostream &out);

    virtual void v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, 
				 Array<OneD, Array<OneD, NekDouble> > &flux) 
    {
      switch(m_expdim)
	{
	case 1:
	  ASSERTL0(false,"1D not implemented for Compressible Flow Equations");
	  break;
	case 2:
	  GetFluxVector2D(i,physfield,flux);
	  break;
	case 3:
	  ASSERTL0(false,"3D not implemented for Compressible Flow Equations");
	  break;
	default:
	  ASSERTL0(false,"Illegal dimension");
	}
    }
    
    
    virtual void v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
				 Array<OneD, Array<OneD, NekDouble> > &numfluxX,
				 Array<OneD, Array<OneD, NekDouble> > &numfluxY)
    {
      switch(m_expdim)
	{
	case 1:
	  ASSERTL0(false,"1D not implemented for Compressible Flow Equations");
	  break;
	case 2:
	  NumericalFlux2D(physfield,numfluxX,numfluxY);
	  break;
	case 3:
	  ASSERTL0(false,"3D not implemented for Compressible Flow Equations");
	  break;
	default:
	  ASSERTL0(false,"Illegal dimension");
	}      
    }

    void WallBoundary2D(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray);

    void SymmetryBoundary(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray);

    virtual NekDouble v_GetTimeStep(const Array<OneD, Array<OneD,NekDouble> > physarray, 
			       const Array<OneD,int> ExpOrder, 
			       const Array<OneD,NekDouble> CFLDG, NekDouble timeCFL);

    virtual void v_SetInitialConditions(NekDouble initialtime = 0.0,
					bool dumpInitialConditions = true)
    {
    }

    void GetVelocityVector(const Array<OneD, Array<OneD, NekDouble> > &physfield,
			   Array<OneD, Array<OneD, NekDouble> > &velocity);
    
    void GetSoundSpeed(const Array<OneD, Array<OneD, NekDouble> > &physfield,
		       Array<OneD, NekDouble> &pressure,
		       Array<OneD, NekDouble> &soundspeed);
    
    void GetMach(Array<OneD, Array<OneD, NekDouble> > &physfield,
		 Array<OneD, NekDouble> &soundspeed,
		 Array<OneD, NekDouble> &mach);

    void GetTemperature(Array<OneD, Array<OneD, NekDouble> > &physfield,
			Array<OneD, NekDouble> &pressure,
			Array<OneD, NekDouble> &temperature);
    
    private:
    
    void GetFluxVector2D(const int i, const Array<OneD, const Array<OneD, NekDouble> > &physfield, 
			 Array<OneD, Array<OneD, NekDouble> > &flux);

    void NumericalFlux2D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
			 Array<OneD, Array<OneD, NekDouble> > &numfluxX, 
			 Array<OneD, Array<OneD, NekDouble> > &numfluxY);
    
    void GetPressure(const Array<OneD, const Array<OneD, NekDouble> > &physfield,
 		     Array<OneD, NekDouble> &pressure); 
    
    void RiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
		       NekDouble rhoR, NekDouble rhouR, NekDouble rhovR, NekDouble ER,
		       NekDouble &rhoflux, NekDouble &rhouflux, NekDouble &rhovflux, NekDouble &Eflux);
    
    void HLLRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
			  NekDouble rhoR, NekDouble rhouR, NekDouble rhovR, NekDouble ER,
			  NekDouble &rhoflux, NekDouble &rhouflux, NekDouble &rhovflux, NekDouble &Eflux);
    
    void HLLCRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
			   NekDouble rhoR, NekDouble rhouR, NekDouble rhovR, NekDouble ER,
			   NekDouble &rhoflux, NekDouble &rhouflux, NekDouble &rhovflux, NekDouble &Eflux);
    
    void LFRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
			 NekDouble rhoR, NekDouble rhouR, NekDouble rhovR, NekDouble ER,
			 NekDouble &rhoflux, NekDouble &rhouflux, NekDouble &rhovflux, NekDouble &Eflux);
	    
    void AverageRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
			      NekDouble rhoR, NekDouble rhouR, NekDouble rhovR, NekDouble ER,
			      NekDouble &rhoflux, NekDouble &rhouflux, NekDouble &rhovflux, NekDouble &Eflux);
    
    void ExactRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
			    NekDouble rhoR, NekDouble rhouR, NekDouble rhovR, NekDouble ER,
			    NekDouble &rhoflux, NekDouble &rhouflux, NekDouble &rhovflux, NekDouble &Eflux);

    void AUSMRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
			   NekDouble rhoR, NekDouble rhouR, NekDouble rhovR, NekDouble ER,
			   NekDouble &rhoflux, NekDouble &rhouflux, NekDouble &rhovflux, NekDouble &Eflux);
    
    NekDouble M1Function(int A, NekDouble M);
    NekDouble M2Function(int A, NekDouble M);
    NekDouble M4Function(int A, NekDouble beta, NekDouble M);
    NekDouble P5Function(int A, NekDouble alpha, NekDouble M);

    Array<OneD,NekDouble> GetStdVelocity(const Array<OneD, Array<OneD,NekDouble> > inarray);
    
  };
}
#endif
