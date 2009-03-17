///////////////////////////////////////////////////////////////////////////////
//
// File EulerEquations.h
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
// Description: Euler Equations definition in two-dimensions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_EULEREQUATIONS_EULEREQUATIONS_H
#define NEKTAR_SOLVERS_EULEREQUATIONS_EULEREQUATIONS_H

#include <MultiRegions/DisContField2D.h>
#include <Auxiliary/ADRBase.h>

namespace Nektar
{     
    /**
     * 
     * 
     **/
    
    class EulerEquations: public ADRBase
    {
    public:           

        /**
         * Default constructor. 
         * 
         **/ 
        EulerEquations();


        /**
         * Constructor.
         * /param 
         * 
         *
         **/
        EulerEquations(string &fileStringName);
	
	void GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> >&physfield, 
			   Array<OneD, Array<OneD, NekDouble> >&flux);
	
	void NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &numfluxX,
			   Array<OneD, Array<OneD, NekDouble> > &numfluxY);
        
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


	void GeneralTimeIntegration(int nsteps, 
				    LibUtilities::TimeIntegrationMethod IntMethod,
				    LibUtilities::TimeIntegrationSchemeOperators ode);
	
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
	
	void SetIsenTropicVortex(void);
	
	void GetExactIsenTropicVortex(Array<OneD, NekDouble> &outarray, int field);

    protected:

    private: 
        int m_infosteps;  ///< dump info to stdout at steps time

	NekDouble m_gamma; 

	void GetPressure(Array<OneD, Array<OneD, NekDouble> > &physfield,
			 Array<OneD, NekDouble> &pressure);


	void SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble> > &physarray, NekDouble time);
	
	void WallBoundary(int bcRegion, Array<OneD, Array<OneD, NekDouble> > &physarray);


	void RiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
			   NekDouble rhoR, NekDouble rhouR, NekDouble rhovR, NekDouble ER,
			   NekDouble &rhoflux, NekDouble &rhouflux, NekDouble &rhovflux, NekDouble &Eflux);
	
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
    
    typedef boost::shared_ptr<EulerEquations> EulerEquationsSharedPtr;
    
} //end of namespace

#endif //NEKTAR_SOLVERS_EULEREQUATIONS_EULEREQUATIONS_H

/**
* $Log: EulerEquations.h,v $
* Revision 1.2  2009/02/02 16:10:16  claes
* Update to make SWE, Euler and Boussinesq solvers up to date with the time integrator scheme. Linear and classical Boussinsq solver working
*
* Revision 1.1  2009/01/13 10:59:32  pvos
* added new solvers file
*
* Revision 1.1  2008/11/17 08:42:06  claes
* Initial commit of restructured Euler Solver
*
* Revision 1.1  2008/08/22 09:48:23  pvos
* Added Claes' AdvectionDiffusionReaction, ShallowWater and Euler solver
*
**/
