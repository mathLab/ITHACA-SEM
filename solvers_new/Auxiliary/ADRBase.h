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
// Euler Class and ShallowWaterEquations
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
     * It is basically a class handling vector valued fields where every field is
     * a DisContField2D class 
     */
    
    class ADRBase
    {
    public:           
      
        /**
         * Default constructor. 
         * 
         */ 
        ADRBase();


        /**
         * Constructor.
         * /param 
         * 
         *
         */
        ADRBase(string &fileStringName,
                bool UseInputFileForProjectionType = false, 
                bool UseContinuousField = false);

        void SetADRBase(SpatialDomains::MeshGraphSharedPtr &graph,
                        int nvariables);
                

        void SetInitialConditions(NekDouble initialtime = 0.0);
        
        NekDouble L2Error(int field);

        void EvaluateExactSolution(int field, Array<OneD, NekDouble > &outfield);

        void EvaluateUserDefinedEqn(Array<OneD, Array<OneD, NekDouble> > &outfield);
      
        void FwdTrans(const ADRBase &In);
      
        void BwdTrans(void);
     
        void BwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, int field_no);

        void BwdTrans(Array<OneD, Array<OneD, NekDouble> > &outarray);
    
        /**
         * 
         * \return ...
         */
        int GetNcoeffs(void)
        {
            return m_fields[0]->GetNcoeffs();
        }
      
        int GetNpoints(void)
        {
            return m_fields[0]->GetTrace()->GetNpoints();
        }

        int GetTracePointsTot(void)
        {
            return m_fields[0]->GetTrace()->GetNpoints();
        }

        int GetPointsTot(void)
        {
            return m_fields[0]->GetPointsTot();
        }
      
        inline int GetNvariables(void) const
        {
            return m_fields.num_elements();
        }
      
        inline NekDouble GetTime(void)
        {
            return m_time;
        }
      
        inline void SetTime(NekDouble time)
        {
            m_time = time;
        }

        inline NekDouble GetTimeStep(void)
        {
            return m_timestep;
        }

        inline void SetTimeStep(NekDouble dt)
        {
            m_timestep = dt;
        }
      
        inline int GetSteps(void)
        {
            return m_steps;
        }

        inline void SetSteps(int steps)
        {
            m_steps = steps;
        }

        inline int GetChecksteps(void)
        {
            return m_checksteps;
        }

        inline void SetChecksteps(int check)
        {
            m_checksteps = check;
        }

        void GetCoords(Array<OneD, NekDouble>& x0, Array<OneD, NekDouble>& x1, Array<OneD, NekDouble>& x2=NullNekDouble1DArray);
      
        const Array<OneD, const NekDouble> &GetPhys(int field_no = -1);

        void GetPhys(Array<OneD, Array<OneD, NekDouble> >&outarray);


        //-----------------------------
        /**
         *
         */
        void SetPhys(Array<OneD, Array<OneD, NekDouble> >&inarray, int field_no = -1);
        //-----------------------------
      
        void SetPhys(Array<OneD, NekDouble> &inarray, int field_no = -1);
      
        void WriteToFile(std::ofstream &out, OutputFormat format, int field_no);
      
        NekDouble L2(const MultiRegions::ExpList2D &In, int field_no);
      
        void ExtractTracePhys(Array<OneD, NekDouble> &out, int field_no);

        void GetFwdBwdTracePhys(Array<OneD, NekDouble> &Fwd, Array<OneD,NekDouble> &Bwd, int field_no);
      
        void GetTraceNormals(Array<OneD, Array<OneD, NekDouble> > &Normals);

        void UpwindTrace(const Array<OneD, Array<OneD, NekDouble> > &Vel, 
                         const Array<OneD, const NekDouble> &Fwd, 
                         const Array<OneD, const NekDouble> &Bwd, 
                         Array<OneD, NekDouble> &Upwind); 
      
        /*  inline const boost::shared_ptr<Array<OneD, Array< OneD, NekDouble> > >GetTraceInnerNormals(void) const */
        /* 	{ */
        /* 	  return m_fields[0]->GetTraceInnerNormals();  */
        /* 	} */

     
        void IProductWRTDerivBase(const int dir, const Array< OneD, const NekDouble > &in,
                                  Array< OneD, NekDouble > &out, int field_no);
      
        void IProductWRTBase(const Array< OneD, const NekDouble > &in,
                             Array< OneD, NekDouble > &out, int field_no);

        void MultiplyByElmtInvMass(const Array< OneD, const NekDouble > &in,
                                   Array< OneD, NekDouble > &out, int field_no = 0);
        void AddTraceIntegral(Array<OneD, const NekDouble> &Fx, 
                              Array<OneD, const NekDouble> &Fy, 
                              Array<OneD, NekDouble> &outarray,
                              int field_no);
      
        Array<OneD, NekDouble> &UpdateCoeffs(int field_no);
	   
        const Array<OneD, const NekDouble> &GetCoeffs(int field_no);


        void WeakAdvectionGreensDivergenceForm(const Array<OneD, Array<OneD, NekDouble> > &F, Array<OneD, NekDouble> &outarray);

        void WeakAdvectionDivergenceForm(const Array<OneD, Array<OneD, NekDouble> > &F, Array<OneD, NekDouble> &outarray);
        
        void WeakAdvectionNonConservativeForm(const Array<OneD, Array<OneD, NekDouble> > &V, const Array<OneD, const NekDouble> &u, Array<OneD, NekDouble> &outarray);
        
        void WeakDGAdvection(const Array<OneD, Array<OneD, NekDouble> >& InField, Array<OneD, Array<OneD, NekDouble> >& OutField );

        void Output     (void);
        void Checkpoint_Output(const int n);

        // virtual functions wrappers
        void GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> >&physfield, Array<OneD, Array<OneD, NekDouble> >&flux)
        {
            v_GetFluxVector(i,physfield, flux);
        }
        
        void NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, 
                           Array<OneD, Array<OneD, NekDouble> > &numflux)
        {
            v_NumericalFlux(physfield, numflux);
        }

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

        enum ProjectionType
        {
            eGalerkin, 
            eDiscontinuousGalerkin
        };

    protected:
        Array<OneD, MultiRegions::ExpListSharedPtr> m_fields; ///< Array holding all dependent variables

        SpatialDomains::BoundaryConditionsSharedPtr m_boundaryConditions;
        std::string m_sessionName; ///< Name of the sessions
        NekDouble m_time;         ///< continous time
        NekDouble m_timestep;     ///< time step size
        int m_steps;          ///< number of steps to be taken during the simulation
        int m_checksteps;     ///< number of steps between dumping check files 
        int m_spacedim;        ///< Dimension of the space (> expansion dim). 

        enum ProjectionType m_projectionType; ///< Type of projection, i.e. Galerkin or DG 

        Array<OneD, Array<OneD, NekDouble> > m_traceNormals;
        
    private: 
        
        virtual void v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> >&physfield, Array<OneD, Array<OneD, NekDouble> >&flux)
        {
            ASSERTL0(false,"This function is not valid for the Base class");
        }
        
        virtual void v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &numflux)
        {
            ASSERTL0(false,"This function is not valid for the Base class");
        }
        
    };
    
    typedef boost::shared_ptr<ADRBase> ADRBaseSharedPtr;
    
} //end of namespace

#endif //NEKTAR_SOLVERS_AUXILIARY_ADRBASE_H

/**
* $Log: ADRBase.h,v $
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
