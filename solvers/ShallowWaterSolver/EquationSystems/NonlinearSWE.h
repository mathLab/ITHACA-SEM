///////////////////////////////////////////////////////////////////////////////
//
// File NonlinearSWE.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Nonlinear Shallow water equations in conservative variables
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_SHALLOWWATERSOLVER_EQUATIONSYSTEMS_NONLINEARSWE_H
#define NEKTAR_SOLVERS_SHALLOWWATERSOLVER_EQUATIONSYSTEMS_NONLINEARSWE_H

#include <ShallowWaterSolver/EquationSystems/ShallowWaterSystem.h>

namespace Nektar
{     
  

  /**
   * 
   * 
   **/
 
  class NonlinearSWE : public ShallowWaterSystem
  {
  public:
      friend class MemoryManager<NonlinearSWE>;

    /// Creates an instance of this class
      static SolverUtils::EquationSystemSharedPtr create(
          const LibUtilities::SessionReaderSharedPtr& pSession,
          const SpatialDomains::MeshGraphSharedPtr& pGraph)
    {
      SolverUtils::EquationSystemSharedPtr p = MemoryManager<NonlinearSWE>
          ::AllocateSharedPtr(pSession, pGraph);
      p->InitObject();
      return p;
    }
    /// Name of class
    static std::string className;
    
    virtual ~NonlinearSWE();

  protected:

    NonlinearSWE(const LibUtilities::SessionReaderSharedPtr& pSession,
                 const SpatialDomains::MeshGraphSharedPtr& pGraph);

    virtual void v_InitObject();
    
    void DoOdeRhs(const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
		  Array<OneD,  Array<OneD, NekDouble> > &outarray,
		  const NekDouble time);
    
    void DoOdeProjection(const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
			 Array<OneD,  Array<OneD, NekDouble> > &outarray,
			 const NekDouble time);

    void GetFluxVector(
     const Array<OneD, const Array<OneD, NekDouble> > &physfield, 
     Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux);

    virtual void v_GenerateSummary(SolverUtils::SummaryList& s);
    
    virtual void v_PrimitiveToConservative( );
    
    virtual void v_ConservativeToPrimitive( );


    private:

	void NumericalFlux1D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
			     Array<OneD, Array<OneD, NekDouble> > &numfluxX);
	
	void NumericalFlux2D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
			     Array<OneD, Array<OneD, NekDouble> > &numfluxX, 
			     Array<OneD, Array<OneD, NekDouble> > &numfluxY);


	void SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble> > &physarray, NekDouble time);

        void WallBoundary2D(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &Fwd, Array<OneD, Array<OneD, NekDouble> > &physarray);
	void WallBoundary(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &Fwd, Array<OneD, Array<OneD, NekDouble> > &physarray);

	void AddCoriolis( const Array<OneD,  const Array<OneD, NekDouble> > &physarray,
			       Array<OneD,       Array<OneD, NekDouble> > &outarray);

	void AddVariableDepth(const Array<OneD, const Array<OneD, NekDouble> > &physarray,
			      Array<OneD, Array<OneD, NekDouble> > &outarray);
	
	void ConservativeToPrimitive(const Array<OneD, const Array<OneD, NekDouble> >&physin,
				     Array<OneD,       Array<OneD, NekDouble> >&physout);
	void PrimitiveToConservative(const Array<OneD, const Array<OneD, NekDouble> >&physin,
				     Array<OneD,       Array<OneD, NekDouble> >&physout);

	void GetVelocityVector(
			       const Array<OneD, Array<OneD, NekDouble> > &physfield,
			       Array<OneD, Array<OneD, NekDouble> > &velocity);
  };
 
}

#endif 

