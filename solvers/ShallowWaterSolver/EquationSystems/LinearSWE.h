///////////////////////////////////////////////////////////////////////////////
//
// File LinearSWE.h
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
// Description: Linearized Shallow water equations in primitive variables
//              valid for a constant water depth
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_SHALLOWWATERSOLVER_EQUATIONSYSTEMS_LINEARSWE_H
#define NEKTAR_SOLVERS_SHALLOWWATERSOLVER_EQUATIONSYSTEMS_LINEARSWE_H

#include <ShallowWaterSolver/EquationSystems/ShallowWaterSystem.h>

namespace Nektar
{     
  

  /**
   * 
   * 
   **/
  class LinearSWE : public ShallowWaterSystem
  {
  public:
      friend class MemoryManager<LinearSWE>;

    /// Creates an instance of this class
    static EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession)
    {
        EquationSystemSharedPtr p = MemoryManager<LinearSWE>::AllocateSharedPtr(pSession);
        p->InitObject();
        return p;
    }
    /// Name of class
    static std::string className;
    
    virtual ~LinearSWE();

  protected:
    LinearSWE(const LibUtilities::SessionReaderSharedPtr& pSession);
    
    virtual void v_InitObject();
    
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

    void v_PrimitiveToConservative()
    {
      m_primitive = true;
    }
      
    void v_ConservativeToPrimitive()
    {
      m_primitive = true;
    }
  
    private:

	void GetFluxVector1D(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, 
			     Array<OneD, Array<OneD, NekDouble> > &flux);
	
	void GetFluxVector2D(const int i, const Array<OneD, const Array<OneD, NekDouble> > &physfield, 
			     Array<OneD, Array<OneD, NekDouble> > &flux);

	void NumericalFlux1D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
			     Array<OneD, Array<OneD, NekDouble> > &numfluxX);
	
	void NumericalFlux2D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
			     Array<OneD, Array<OneD, NekDouble> > &numfluxX, 
			     Array<OneD, Array<OneD, NekDouble> > &numfluxY);

	void AverageLinearNumericalFlux2D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
					  Array<OneD, Array<OneD, NekDouble> > &numfluxX, 
					  Array<OneD, Array<OneD, NekDouble> > &numfluxY);

	void RiemannLinearNumericalFlux2D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
					  Array<OneD, Array<OneD, NekDouble> > &numfluxX, 
					  Array<OneD, Array<OneD, NekDouble> > &numfluxY);

	void RiemannSolver(NekDouble hL,NekDouble uL,NekDouble vL,NekDouble hR,NekDouble uR, 
			   NekDouble vR, NekDouble &hflux, NekDouble &huflux,NekDouble &hvflux );

	void SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble> > &physarray, NekDouble time);

	void WallBoundary1D(int bcRegion, Array<OneD, Array<OneD, NekDouble> > &physarray);

	void WallBoundary2D(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray);

	void AddCoriolis(const Array<OneD, const Array<OneD, NekDouble> > &physarray,
			       Array<OneD,       Array<OneD, NekDouble> > &outarray);
    };
}

#endif

