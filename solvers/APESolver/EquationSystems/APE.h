///////////////////////////////////////////////////////////////////////////////
//
// File APE.h
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
// Description: Acoustic perturbation equations in conservative variables
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_APESOLVER_EQUATIONSYSTEMS_APE_H
#define NEKTAR_SOLVERS_APESOLVER_EQUATIONSYSTEMS_APE_H

#include <SolverUtils/EquationSystem.h>
#include <APESolver/EquationSystems/APESystem.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{     
  

  /**
   * 
   * 
   **/
 
  class APE : public APESystem
  {
  public:
	  friend class MemoryManager<APE>;
    
	/// Creates an instance of this class
    static EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession)
    {
		EquationSystemSharedPtr p = MemoryManager<APE>::AllocateSharedPtr(pSession);
		p->InitObject();
		return p;
    }
    /// Name of class
    static std::string className;
    
    virtual ~APE();
    
  protected:
	  
	APE(const LibUtilities::SessionReaderSharedPtr& pSession);
	  
	virtual void v_InitObject();
   
	Array<OneD, Array<OneD, NekDouble> > basefield;

	  
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
	  ASSERTL0(false,"1D not implemented for Acoustic perturbation equations");
	  break;
	case 2:
	  GetFluxVector2D(i,physfield,flux);
	  break;
	case 3:
	  ASSERTL0(false,"3D not implemented for Acoustic perturbation equations");
	  break;
	default:
	  ASSERTL0(false,"Illegal dimension");
	}
    }
    
	  
    virtual void v_PrintSummary(std::ostream &out); 

    virtual void v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
				 Array<OneD, Array<OneD, NekDouble> > &numfluxX,
				 Array<OneD, Array<OneD, NekDouble> > &numfluxY)
    {
      switch(m_expdim)
	{
	case 1:
	  ASSERTL0(false,"1D not implemented for Acoustic perturbation equations");
	  break;
	case 2:
	  NumericalFlux2D(physfield,numfluxX,numfluxY);
	  break;
	case 3:
	  ASSERTL0(false,"3D not implemented for Acoustic perturbation equations");
	  break;
	default:
	  ASSERTL0(false,"Illegal dimension");
	}
      
    }
      
    virtual void v_PrimitiveToConservative( );
    
    virtual void v_ConservativeToPrimitive( );

    void InitialiseBaseFlowAnalytical(Array<OneD, Array<OneD, NekDouble> > &base, const NekDouble time);
	  
    void GetSource(Array<OneD, NekDouble> &source, const NekDouble time);


    void AddSource(const Array< OneD, Array< OneD, NekDouble > > &inarray, Array< OneD, Array< OneD, NekDouble > > &outarray);


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


	void RiemannSolverUpwind(NekDouble hL,NekDouble uL,NekDouble vL,NekDouble hR,NekDouble uR, NekDouble vR, NekDouble P0, NekDouble U0, NekDouble V0,
			    NekDouble &hflux, NekDouble &huflux,NekDouble &hvflux );

	void SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble> > &physarray, NekDouble time);

	void WallBoundary1D(int bcRegion, Array<OneD, Array<OneD, NekDouble> > &physarray);

	void WallBoundary2D(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray);

	void ConservativeToPrimitive(const Array<OneD, const Array<OneD, NekDouble> >&physin,
				     Array<OneD,       Array<OneD, NekDouble> >&physout);
	void PrimitiveToConservative(const Array<OneD, const Array<OneD, NekDouble> >&physin,
				     Array<OneD,       Array<OneD, NekDouble> >&physout);
	

  };
}

#endif 

