///////////////////////////////////////////////////////////////////////////////
//
// File EulerCFE.h
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
// Description: Euler equations in conservative variables without artificial diffusion
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_EULERCFE_H
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_EULERCFE_H

#include <CompressibleFlowSolver/EquationSystems/CompressibleFlowSystem.h>

namespace Nektar
{  

  enum ProblemType
  {           
    eGeneral,          ///< No problem defined - Default Inital data
    eIsentropicVortex, ///< Isentropic Vortex
    eSubsonicCylinder, ///< Subsonic Cylinder
    eRinglebFlow,      ///< Ringleb Flow
    SIZE_ProblemType   ///< Length of enum list
  };
  
  const char* const ProblemTypeMap[] =
    {
      "General",
      "IsentropicVortex",
      "SubsonicCylinder",
      "RinglebFlow"
    };
  
  /**
   * 
   * 
   **/
  class EulerCFE : public CompressibleFlowSystem
  {
  public:
    /// Creates an instance of this class
    static EquationSystemSharedPtr create(SessionReaderSharedPtr& pSession) 
    {
      return MemoryManager<EulerCFE>::AllocateSharedPtr(pSession);
    }
    /// Name of class
    static std::string className;
    
    EulerCFE(SessionReaderSharedPtr& pSession);
    
    virtual ~EulerCFE();

    ///< problem type selector
    ProblemType                                     m_problemType;   
    
  protected:

    /// Print a summary of time stepping parameters.
    virtual void v_PrintSummary(std::ostream &out);

    void DoOdeRhs(const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
		  Array<OneD,  Array<OneD, NekDouble> > &outarray,
		  const NekDouble time);
    
    void DoOdeProjection(const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
			 Array<OneD,  Array<OneD, NekDouble> > &outarray,
			 const NekDouble time);

    virtual void v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, 
				 Array<OneD, Array<OneD, NekDouble> > &flux) 
    {
      switch(m_expdim)
	{
	case 1:
	  ASSERTL0(false,"1D not implemented for Shallow Water Equations");
	  break;
	case 2:
	  GetFluxVector2D(i,physfield,flux);
	  break;
	case 3:
	  ASSERTL0(false,"3D not implemented for Shallow Water Equations");
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
	  ASSERTL0(false,"1D not implemented for Shallow Water Equations");
	  break;
	case 2:
	  NumericalFlux2D(physfield,numfluxX,numfluxY);
	  break;
	case 3:
	  ASSERTL0(false,"3D not implemented for Shallow Water Equations");
	  break;
	default:
	  ASSERTL0(false,"Illegal dimension");
	}      
    }

    virtual void v_SetInitialConditions(NekDouble initialtime = 0.0,
					bool dumpInitialConditions = true);
    
    private:
    
    void GetFluxVector2D(const int i, const Array<OneD, const Array<OneD, NekDouble> > &physfield, 
			 Array<OneD, Array<OneD, NekDouble> > &flux);

    void NumericalFlux2D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
			 Array<OneD, Array<OneD, NekDouble> > &numfluxX, 
			 Array<OneD, Array<OneD, NekDouble> > &numfluxY);

    void SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble> > &physarray, NekDouble time);

    void WallBoundary2D(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray);
    
    void GetPressure(const Array<OneD, const Array<OneD, NekDouble> > &physfield,
 		     Array<OneD, NekDouble> &pressure); 
    
    void SymmetryBoundary(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray);
    
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

    // Isentropic Vortex Test Case
    void SetIsenTropicVortex(NekDouble initialtime);
    void SetBoundaryIsentropicVortex(int bcRegion, NekDouble time, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray);
    void GetExactIsenTropicVortex(Array<OneD, NekDouble> &outarray, int field);
    
  };
}
#endif
