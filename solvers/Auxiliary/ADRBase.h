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
                
        void SetInitialConditions(const int eqntype = 1, NekDouble initialtime = 0.0);

        void SetPhysForcingFunctions(Array<OneD, MultiRegions::ExpListSharedPtr> &force);
            
        void EvaluateExactSolution(int field, Array<OneD, NekDouble > &outfield, 
                                   const NekDouble time, const int eqntype = 0);
        
        void EvaluateUserDefinedEqn(Array<OneD, Array<OneD, NekDouble> > &outfield);
        
        void AdvectionNonConservativeForm(const Array<OneD, Array<OneD, NekDouble> > &V, 
                                          const Array<OneD, const NekDouble> &u, Array<OneD, NekDouble> &outarray, 
                                          Array<OneD, NekDouble> &wk = NullNekDouble1DArray);
            
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
			     bool InFieldIsInPhysSpace = false, int nvariables = 0); 
				 
        void WeakDGDiffusion3D(const Array<OneD, Array<OneD, NekDouble> >& InField, 
                             Array<OneD, Array<OneD, NekDouble> >& OutField,
                             bool NumericalFluxIncludesNormal = true, 
                             bool InFieldIsInPhysSpace = false);

        void WeakDGDiffusion(const Array<OneD, Array<OneD, NekDouble> >& InField, 
			     Array<OneD, Array<OneD, NekDouble> >& OutField,
			     bool NumericalFluxIncludesNormal = true, 
			     bool InFieldIsInPhysSpace = false);

        NekDouble L2Error(int field, const int eqntype = 0, 
			  const Array<OneD,NekDouble> &exactsoln = NullNekDouble1DArray);
	
        void Output     (void);
	
        void Checkpoint_Output(const int n);

        void WriteFld(std::ofstream &outfile);
        
	void Array_Output(const int n, std::string name, 
			  const Array<OneD, const NekDouble>&inarray, bool IsInPhysicalSpace);
	
	void Summary          (std::ostream &out);
	void SessionSummary   (std::ostream &out);
	void TimeParamSummary (std::ostream &out);
	void Unitlength(Array<OneD, Array<OneD, NekDouble> > &array);
        
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
		
	virtual void GetFluxVector(const int i, const int j, 
                                   Array<OneD, Array<OneD, NekDouble> > &physfield, 
                                   Array<OneD, Array<OneD, NekDouble> > &flux)
        {
            v_GetFluxVector(i,j,physfield,flux);
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

	void NumFluxforScalar(Array<OneD, Array<OneD, NekDouble> > &ufield, 
                              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
        {
            v_NumFluxforScalar(ufield, uflux);
        }
	
	void NumFluxforVector(Array<OneD, Array<OneD, NekDouble> > &ufield,
			      Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
			      Array<OneD, Array<OneD, NekDouble> >  &qflux)
        {
	  v_NumFluxforVector(ufield,qfield, qflux);
        }

	void NumFluxforDiff(Array<OneD, Array<OneD, NekDouble> > &ufield, 
			    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
	{
	  v_NumFluxforDiff(ufield, uflux);
	}
	
	void NumFluxforDiff(Array<OneD, Array<OneD, NekDouble> > &ufield,
			    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
			    Array<OneD, Array<OneD, NekDouble> >  &qflux)
	{
	  v_NumFluxforDiff(ufield,qfield, qflux);
	}

	NekDouble AdvectionSphere(const NekDouble x0j, const NekDouble x1j,
				  const NekDouble x2j, const NekDouble time)
	{
	  return v_AdvectionSphere(x0j, x1j, x2j, time);
	}

	NekDouble Morphogenesis(const int field, const NekDouble x0j, const NekDouble x1j, 
				const NekDouble x2j, const NekDouble time)
	{
	  return v_Morphogenesis(field, x0j, x1j, x2j, time);
	}
        
    protected:
        Array<OneD, MultiRegions::ExpListSharedPtr> m_fields; ///< Array holding all dependent variables
        SpatialDomains::BoundaryConditionsSharedPtr m_boundaryConditions;
        SpatialDomains::MeshGraphSharedPtr          m_graph;
        std::string m_sessionName;   ///< Name of the sessions
        NekDouble m_time;            ///< continous time
        NekDouble m_timestep;        ///< time step size
        int m_steps;                 ///< number of steps to be taken during the simulation
        int m_checksteps;            ///< number of steps between dumping check files 
        int m_spacedim;              ///< Dimension of the space (> expansion dim)
	int m_expdim;                ///< Dimension of the expansion

        enum ProjectionType m_projectionType; ///< Type of projection, i.e. Galerkin or DG 

        Array<OneD, Array<OneD, NekDouble> > m_traceNormals; ///< Array holding the forward normals

        Array<OneD, Array<OneD, Array<OneD,NekDouble> > > m_curvature; // 1 by nvariable by nq
        Array<OneD, Array<OneD, Array<OneD,NekDouble> > > m_tanbasis; // 2 by m_spacedim by nq 
	
        int NoCaseStringCompare(const string & s1, const string& s2) ;

        // Here for consistency purposes with old version
	int nocase_cmp(const string & s1, const string& s2)
        {
            return NoCaseStringCompare(s1,s2);
        }
	    
    private: 
        
        virtual void v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> >&physfield, 
				     Array<OneD, Array<OneD, NekDouble> >&flux)
        {
            ASSERTL0(false,"v_GetFluxVector: This function is not valid for the Base class");
        }
        
        virtual void v_GetFluxVector(const int i, const int j, Array<OneD, Array<OneD, NekDouble> >&physfield, 
				     Array<OneD, Array<OneD, NekDouble> >&flux)
        {
            ASSERTL0(false,"v_GetqFluxVector: This function is not valid for the Base class");
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
	
	virtual void v_NumFluxforScalar(Array<OneD, Array<OneD, NekDouble> > &ufield, 
                                        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
        {
	  ASSERTL0(false,"v_NumFluxforScalar: This function is not valid for the Base class");
        }
       
        virtual void v_NumFluxforVector(Array<OneD, Array<OneD, NekDouble> > &ufield,
					Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
					Array<OneD, Array<OneD, NekDouble > >  &qflux)
        {
	  ASSERTL0(false,"v_NumFluxforVector: This function is not valid for the Base class");
        }    

	virtual void v_NumFluxforDiff(Array<OneD, Array<OneD, NekDouble> > &ufield, 
				      Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
        {
	  ASSERTL0(false,"v_NumFluxforDiffu: This function is not valid for the Base class");
        }

	 virtual void v_NumFluxforDiff(Array<OneD, Array<OneD, NekDouble> > &ufield,
				       Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
				       Array<OneD, Array<OneD, NekDouble > >  &qflux)
        {
	  ASSERTL0(false,"v_NumFluxforDiffq: This function is not valid for the Base class");
        }

	 virtual NekDouble v_AdvectionSphere(const NekDouble x0j, const NekDouble x1j,
					     const NekDouble x2j, const NekDouble time)
	{
	  ASSERTL0(false,"v_AdvectionSphere: This function is not valid for the Base class");
	}

	virtual NekDouble v_Morphogenesis(const int field, const NekDouble x0j, const NekDouble x1j, 
					  const NekDouble x2j, const NekDouble time)
	{
	  ASSERTL0(false,"v_Morphogenesis: This function is not valid for the Base class");
	}

    };
    
    typedef boost::shared_ptr<ADRBase> ADRBaseSharedPtr;
    
} //end of namespace

#endif //NEKTAR_SOLVERS_AUXILIARY_ADRBASE_H

/**
* $Log: ADRBase.h,v $
* Revision 1.9  2009/07/09 21:29:13  sehunchun
* VS: ----------------------------------------------------------------------
* Add SetUpSurfaceNormal function..
*
* Revision 1.8  2009/07/02 15:57:36  sehunchun
* "ReadBoundaryCondition" options with extenstion to 2D geometry imbedded in 3D
*
* Revision 1.7  2009/04/29 20:45:55  sherwin
* Update for new definition of enum
*
* Revision 1.6  2009/04/27 21:37:14  sherwin
* Updated to dump .fld and .chk file in compressed coefficient format
*
* Revision 1.5  2009/03/10 23:37:14  claes
* Updated the ShallowWaterSolver to work with the general timestepping scheme
*
* Revision 1.4  2009/02/28 22:00:38  sehunchun
*  Explicit Diffusion solver is added
*
* Revision 1.3  2009/02/03 14:33:44  pvos
* Modifications for solvers with time-dependent dirichlet BC's
*
* Revision 1.2  2009/02/02 16:10:16  claes
* Update to make SWE, Euler and Boussinesq solvers up to date with the time integrator scheme. Linear and classical Boussinsq solver working
*
* Revision 1.1  2009/01/13 10:59:32  pvos
* added new solvers file
*
* Revision 1.6  2009/01/06 21:11:03  sherwin
* Updates for Virtual ExpList calls
*
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
