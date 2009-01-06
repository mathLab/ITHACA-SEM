///////////////////////////////////////////////////////////////////////////////
//
// File ADRBase.h
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
// Description: Basic class for AdvectionDiffusionReaction class,
// Euler Class, ShallowWaterEquations and BoussinesqEquations
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_AUXILIARY_ADRBASE_H
#define NEKTAR_SOLVERS_AUXILIARY_ADRBASE_H

#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>

#include <MultiRegions/DisContField1D.h>
#include <MultiRegions/DisContField2D.h>

namespace Nektar
{     
    /**
     * \brief This class is the base class for the development of solvers.
     *
     * It is basically a class handling vector valued fields
     **/
    
    class ADRBase
    {
    public:           
      
        /**
         * Default constructor. 
	 **/ 
        ADRBase();


        /**
         * Constructor.
         **/
        ADRBase(string &fileStringName,
                bool UseInputFileForProjectionType = false, 
                bool UseContinuousField = false);

        void SetADRBase(SpatialDomains::MeshGraphSharedPtr &graph,
                        int nvariables);
                
        void SetInitialConditions(NekDouble initialtime = 0.0);
        void SetPhysForcingFunctions(Array<OneD, MultiRegions::ExpListSharedPtr> &force);
            
	void EvaluateExactSolution(int field, Array<OneD, NekDouble > &outfield);

        void EvaluateUserDefinedEqn(Array<OneD, Array<OneD, NekDouble> > &outfield);
      
     	void WeakAdvectionGreensDivergenceForm(const Array<OneD, Array<OneD, NekDouble> > &F, 
					       Array<OneD, NekDouble> &outarray);

        void WeakAdvectionDivergenceForm(const Array<OneD, Array<OneD, NekDouble> > &F, 
					 Array<OneD, NekDouble> &outarray);
        
        void WeakAdvectionNonConservativeForm(const Array<OneD, Array<OneD, NekDouble> > &V, 
					      const Array<OneD, const NekDouble> &u, 
					      Array<OneD, NekDouble> &outarray);
        
        void WeakDGAdvection(const Array<OneD, Array<OneD, NekDouble> >& InField, 
			     Array<OneD, Array<OneD, NekDouble> >& OutField, 
			     bool NumericalFluxIncludesNormal = true, 
			     bool InFieldIsInPhysSpace = false); 

	NekDouble L2Error(int field);
	
        void Output     (void);
	
        void Checkpoint_Output(const int n);
	
	void Summary          (std::ostream &out);
	void SessionSummary   (std::ostream &out);
	void TimeParamSummary (std::ostream &out);
        
        inline Array<OneD, MultiRegions::ExpListSharedPtr> &UpdateFields(void)
        {
            return m_fields; 
        }

	inline int GetNcoeffs(void)
        {
            return m_fields[0]->GetNcoeffs();
        }
        inline int GetNvariables(void)
        {
            return m_fields.num_elements();
        }
      
        inline const std::string &GetVariable(unsigned int i)
        {
            return m_boundaryConditions->GetVariable(i);
        }


        inline int GetTraceTotPoints(void)
        {
            return GetTraceNpoints();
        }
            
        inline int GetTraceNpoints(void)
        {
	  switch(m_expdim)
	    {
	    case 1:
	      // can't have two &GetTrace in ExpList.h hmm...
	      //return m_fields[0]->GetTrace().num_elements();
	      break;
	    case 2:
	    case 3:
	      return m_fields[0]->GetTrace()->GetNpoints();
	      break;
	    default:
	      ASSERTL0(false,"illegal expansion dimension");
	    }
        }

        inline int GetTotPoints(void)
        {
            return m_fields[0]->GetNpoints();
        }

        inline int GetNpoints(void)
        {
            return m_fields[0]->GetNpoints();
        }
      
	inline int GetSteps(void)
        {
            return m_steps;
        }

        inline NekDouble GetParameter(const std::string &parmName)
        {
            return m_boundaryConditions->GetParameter(parmName);
        }

        void ZeroPhysFields(void);

	enum ProjectionType
        {
            eGalerkin, 
            eDiscontinuousGalerkin
        };
	
	//--------------------------------------------------------------------------
        // virtual functions wrappers
	
        void GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> >&physfield, 
			   Array<OneD, Array<OneD, NekDouble> >&flux)
        {
	  v_GetFluxVector(i,physfield, flux);
        }
	
	void GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> >&physfield, 
			   Array<OneD, Array<OneD, NekDouble> >&fluxX, 
			   Array<OneD, Array<OneD, NekDouble> > &fluxY)
        {
	  v_GetFluxVector(i,physfield, fluxX, fluxY);
        }

	void NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, 
                           Array<OneD, Array<OneD, NekDouble> > &numflux)
        {
            v_NumericalFlux(physfield, numflux);
        }

	void NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, 
                           Array<OneD, Array<OneD, NekDouble> > &numfluxX,
			   Array<OneD, Array<OneD, NekDouble> > &numfluxY)
        {
	  v_NumericalFlux(physfield, numfluxX, numfluxY);
        }
	//--------------------------------------------------------------------------

    protected:
        Array<OneD, MultiRegions::ExpListSharedPtr> m_fields; ///< Array holding all dependent variables
        SpatialDomains::BoundaryConditionsSharedPtr m_boundaryConditions;
        std::string m_sessionName;   ///< Name of the sessions
        NekDouble m_time;            ///< continous time
        NekDouble m_timestep;        ///< time step size
        int m_steps;                 ///< number of steps to be taken during the simulation
        int m_checksteps;            ///< number of steps between dumping check files 
        int m_spacedim;              ///< Dimension of the space (> expansion dim)
	int m_expdim;                ///< Dimension of the expansion

        enum ProjectionType m_projectionType; ///< Type of projection, i.e. Galerkin or DG 

        Array<OneD, Array<OneD, NekDouble> > m_traceNormals; ///< Array holding the forward normals 
        
    private: 
        
        virtual void v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> >&physfield, 
				     Array<OneD, Array<OneD, NekDouble> >&flux)
        {
	  ASSERTL0(false,"v_GetFluxVector: This function is not valid for the Base class");
        }
        
	virtual void v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> >&physfield, 
				     Array<OneD, Array<OneD, NekDouble> >&fluxX, 
				     Array<OneD, Array<OneD, NekDouble> > &fluxY)
        {
	  ASSERTL0(false,"v_GetFluxVector: This function is not valid for the Base class");
	  
        }

        virtual void v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, 
				     Array<OneD, Array<OneD, NekDouble> > &numflux)
        {
	  ASSERTL0(false,"v_NumericalFlux: This function is not valid for the Base class");
        }

	virtual void v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, 
				     Array<OneD, Array<OneD, NekDouble> > &numfluxX,
				     Array<OneD, Array<OneD, NekDouble> > &numfluxY )
        {
	  ASSERTL0(false,"v_NumericalFlux: This function is not valid for the Base class");
        }

    };
    
    typedef boost::shared_ptr<ADRBase> ADRBaseSharedPtr;
    
} //end of namespace

#endif //NEKTAR_SOLVERS_AUXILIARY_ADRBASE_H

/**
* $Log: ADRBase.h,v $
* Revision 1.5  2008/11/17 08:10:07  claes
* Removed functions that were no longer used after the solver library was restructured
*
* Revision 1.4  2008/10/31 10:50:10  pvos
* Restructured directory and CMakeFiles
*
* Revision 1.3  2008/10/29 22:51:07  sherwin
* Updates for const correctness and ODEforcing
*
* Revision 1.2  2008/10/19 15:59:20  sherwin
* Added Summary method
*
* Revision 1.1  2008/10/16 15:25:45  sherwin
* Working verion of restructured AdvectionDiffusionReactionSolver
*

**/
