///////////////////////////////////////////////////////////////////////////////
//
// File ShallowWaterEquations.h
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
// Description: Shallow Water Equations definition in two-dimensions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_SHALLOWWATEREQUATIONS_SHALLOWWATEREQUATIONS_H
#define NEKTAR_SOLVERS_SHALLOWWATEREQUATIONS_SHALLOWWATEREQUATIONS_H

#include <MultiRegions/DisContField2D.h>
#include <Auxiliary/ADRBase.h>

namespace Nektar
{     

  /**
   * 
   * 
   **/
  
  class ShallowWaterEquations: public ADRBase
  {
  public:           
    
    /**
     * Default constructor. 
     * 
     **/ 
    ShallowWaterEquations();
    
    
    /**
     * Constructor.
     * /param 
     * 
     *
     **/
    ShallowWaterEquations(string &fileStringName);
    
    
    void PrimitiveToConservative(const Array<OneD, const Array<OneD, NekDouble> >&physin,
				       Array<OneD,       Array<OneD, NekDouble> >&physout);

    void ConservativeToPrimitive(const Array<OneD, const Array<OneD, NekDouble> >&physin,
				       Array<OneD,       Array<OneD, NekDouble> >&physout);
    
    void SWEAdvectionNonConservativeForm(const Array<OneD, const Array<OneD, NekDouble> >&physarray,
					       Array<OneD,       Array<OneD, NekDouble> >&outarray);
    
   void GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> >&physfield, 
		       Array<OneD, Array<OneD, NekDouble> >&flux)
    {
      switch(m_expdim)
	{
	case 1:
	  GetFluxVector1D(i,physfield,flux);
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
    
    void NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
		       Array<OneD, Array<OneD, NekDouble> > &numfluxX,
		       Array<OneD, Array<OneD, NekDouble> > &numfluxY)
    {
      switch(m_expdim)
	{
	case 1:
	  NumericalFlux1D(physfield,numfluxX);
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
    
    //void ODEforcing(const Array<OneD, const  Array<OneD, NekDouble> >&inarray, 
    //Array<OneD, Array<OneD, NekDouble> >&outarray, NekDouble time);

    void ODElhs(const Array<OneD, const  Array<OneD, NekDouble> >&inarray, 
	              Array<OneD,        Array<OneD, NekDouble> >&outarray, 
		const NekDouble time);
    
    void ODElhsSolve(const Array<OneD, const  Array<OneD, NekDouble> >&inarray, 
		           Array<OneD,        Array<OneD, NekDouble> >&outarray, 
		     const NekDouble time);
    
    void ODErhs(const Array<OneD, const  Array<OneD, NekDouble> >&inarray, 
		      Array<OneD,        Array<OneD, NekDouble> >&outarray, 
		const NekDouble time);
    
    void ODEdirkSolve(const Array<OneD, const  Array<OneD, NekDouble> >&inarray, 
		            Array<OneD,        Array<OneD, NekDouble> >&outarray, 
		      const NekDouble lambda,
		      const NekDouble time);
    
    void ExplicitlyIntegrateAdvection(int nsteps);
    
    void Summary(std::ostream &out);
    
    enum UpwindType
    {           ///< flux not defined
      eNotSet,  ///< averaged (or centred) flux
      eAverage, ///< simple upwind flux
      eUpwind,  ///< local Lax-Friedrich flux
      eLLF,     ///< Harten-Lax-Leer Contact corrected flux
      eHLLC,    ///< Roe flux
      eRoe,    
    };
    
    enum LinearType
    {
      eLinear, 
      eNonLinear
    };
    
    enum VariableType
    {
      ePrimitive, 
      eConservative
    };


    void SetGradientBoundary(Array<OneD, NekDouble> &inarray, NekDouble time, int field_0);
    void GradientWallBoundary2D(int bcRegion, int cnt, Array<OneD, NekDouble> &physarray, int field_0);
    void SetDivergenceBoundary(Array<OneD, Array<OneD, NekDouble> > &inarray, NekDouble time, int field_0, int field_1);
    void DivergenceWallBoundary2D(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray, int field_0, int field_1);
    
    void GradientFluxTerms(const Array<OneD, NekDouble> &in, 
			   Array<OneD, Array<OneD, NekDouble> > &out,
			   int field_0);
    
    void NumericalFluxGradient(const Array <OneD, NekDouble> &physarray,
			       Array<OneD, NekDouble> &outX, 
			       Array<OneD, NekDouble> &outY,
			       int field_0);
    
    void DivergenceFluxTerms(const Array<OneD, const Array<OneD, NekDouble> >&in, 
			     Array<OneD, NekDouble> &out,
			     int field_0, int field_1);
    
    void NumericalFluxDivergence(const Array <OneD, const Array<OneD, NekDouble> > &physarray,
				 Array <OneD, NekDouble> &outX, 
				 Array<OneD, NekDouble> &outY,
				 int field_0, int field_1);
    
    
  protected:
    int m_infosteps;  ///< dump info to stdout at steps time
    NekDouble m_g;    ///< Acceleration of gravity
    Array<OneD, NekDouble>  m_coriolis;
    //	Array<OneD, NekDouble>  m_depth;
    Array<OneD, NekDouble>  m_friction;
    NekDouble m_depth; // constant water depth
  
    enum LinearType m_linearType; ///< linear or nonlinear
    enum VariableType m_variableType; ///< primitive or conservative
    
    void AddCoriolis(Array<OneD, Array<OneD, NekDouble> > &physarray, 
		     Array<OneD, Array<OneD, NekDouble> > &outarray);
    
    void SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble> > &physarray, NekDouble time);
    
  private: 
    
    void EvaluateCoriolis(void);
    
    
    void WallBoundary1D(int bcRegion, Array<OneD, Array<OneD, NekDouble> > &physarray);
    void WallBoundary2D(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray);
    
    
    void RiemannSolver(NekDouble hL,NekDouble huL,NekDouble hvL,NekDouble hR,NekDouble huR, 
		       NekDouble hvR, NekDouble &hflux, NekDouble &huflux,NekDouble &hvflux );
    
    void GetFluxVector1D(const int i, Array<OneD, Array<OneD, NekDouble> >&physfield, 
			 Array<OneD, Array<OneD, NekDouble> >&flux);
    
    void GetFluxVector2D(const int i, Array<OneD, Array<OneD, NekDouble> >&physfield, 
			 Array<OneD, Array<OneD, NekDouble> >&flux);
    
    void NumericalFlux1D(Array<OneD, Array<OneD, NekDouble> > &physfield,
			 Array<OneD, Array<OneD, NekDouble> > &numfluxX);
    
    
    void NumericalFlux2D(Array<OneD, Array<OneD, NekDouble> > &physfield,
			 Array<OneD, Array<OneD, NekDouble> > &numfluxX,
			 Array<OneD, Array<OneD, NekDouble> > &numfluxY);		    
    void LinearNumericalFlux2D(Array<OneD, Array<OneD, NekDouble> > &physfield,
			       Array<OneD, Array<OneD, NekDouble> > &numfluxX,
			       Array<OneD, Array<OneD, NekDouble> > &numfluxY);
    void NonlinearNumericalFlux2D(Array<OneD, Array<OneD, NekDouble> > &physfield,
				  Array<OneD, Array<OneD, NekDouble> > &numfluxX,
				  Array<OneD, Array<OneD, NekDouble> > &numfluxY);
    
    virtual void v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, 
				 Array<OneD, Array<OneD, NekDouble> > &flux) 
    { 
      GetFluxVector(i,physfield,flux);
    }
    
    virtual void v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, 
				 Array<OneD, Array<OneD, NekDouble> > &numfluxX, 
				 Array<OneD, Array<OneD, NekDouble> > &numfluxY )
    {
      NumericalFlux(physfield, numfluxX, numfluxY); 
    }
    
  };
  
  typedef boost::shared_ptr<ShallowWaterEquations> ShallowWaterEquationsSharedPtr;
  
} //end of namespace

#endif //NEKTAR_SOLVERS_SHALLOWWATEREQUATIONS_SHALLOWWATEREQUATIONS_H

/**
* $Log: ShallowWaterEquations.h,v $
* Revision 1.2  2009/02/02 16:10:16  claes
* Update to make SWE, Euler and Boussinesq solvers up to date with the time integrator scheme. Linear and classical Boussinsq solver working
*
* Revision 1.1  2009/01/13 10:59:32  pvos
* added new solvers file
*
* Revision 1.1  2008/11/17 08:37:00  claes
* Restructured Shallow Water Solver. Only 2D DG is properly working
*
* Revision 1.2  2008/09/15 14:54:15  claes
* Small changes associated with the BoussinesqSolver
*
* Revision 1.1  2008/08/22 09:48:23  pvos
* Added Claes' AdvectionDiffusionReaction, ShallowWater and Euler solver
*
**/
