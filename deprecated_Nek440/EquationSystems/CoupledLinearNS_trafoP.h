///////////////////////////////////////////////////////////////////////////////
//
// File CoupledLinearNS.h
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
// Description: Coupled Stokes solver scheme header
//
///////////////////////////////////////////////////////////////////////////////


#include "./CoupledLocalToGlobalC0ContMap.h"
#include "./IncNavierStokes.h"
#include "./CoupledLinearNS.h"
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList2D.h>
#include <boost/shared_ptr.hpp>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include "../Eigen/Dense"
//#include <MultiRegions/GlobalLinSysDirectStaticCond.h>

namespace Nektar
{     

  
    class CoupledLinearNS_trafoP: public CoupledLinearNS
    {
    public:
        friend class MemoryManager<CoupledLinearNS_trafoP>;
        
        /// Creates an instance of this class
        static SolverUtils::EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        {
            SolverUtils::EquationSystemSharedPtr p = MemoryManager<CoupledLinearNS_trafoP>::AllocateSharedPtr(pSession);
            p->InitObject();
            return p;
        }
        /// Name of class
        static std::string className;        
        
        /**
         *  Generate the linearised Navier Stokes solver based on the
         *  static condensation of the interior velocity space and
         *  pressure modes.
         */
        void SetUpCoupledMatrix(const NekDouble lambda = 0.0, 
                                const Array< OneD, Array<OneD, NekDouble>  > &Advfield = NullNekDoubleArrayofArray, 
                                bool IsLinearNSEquation = true);
        
        
        const SpatialDomains::ExpansionMap &GenPressureExp(const SpatialDomains::ExpansionMap &VelExp);
        
        void Solve(void);
        
        /**
         *   Solve the coupled linear Navier-Stokes solve using matrix
         *   systems set up at construction.  The solution is stored
         *   in #m_velocity and #m_pressure.
         */
        void SolveLinearNS(const Array<OneD, Array<OneD, NekDouble> > &forcing);
        
        void SolveLinearNS(const Array<OneD, Array<OneD, NekDouble> > &forcing,
                           Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                           MultiRegions::ExpListSharedPtr &pressure,
                           const int HomogeneousMode = 0);
        
        void SolveUnsteadyStokesSystem(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                                       Array<OneD, Array<OneD, NekDouble> > &outarray, 
                                       const NekDouble time,
                                       const NekDouble a_iixDt);
        
        void EvaluateAdvection(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                               Array<OneD, Array<OneD, NekDouble> > &outarray, 
                               const NekDouble time);
        
        void SolveSteadyNavierStokes(void);
        
        void Continuation(void);
        
        /*void EvaluateNewtonRHS(Array<OneD, Array<OneD, NekDouble> > &Velocity,
         *								 Array<OneD, Array<OneD, NekDouble> > &PreviousForcing,
         *								 Array<OneD, Array<OneD, NekDouble> > &outarray);*/
        
        void EvaluateNewtonRHS(Array<OneD, Array<OneD, NekDouble> > &Velocity,
                               Array<OneD, Array<OneD, NekDouble> > &outarray);
                
        void InfNorm(Array<OneD, Array<OneD, NekDouble> > &inarray,
                     Array<OneD, NekDouble> &outarray);
        
        void L2Norm(Array<OneD, Array<OneD, NekDouble> > &inarray,
                    Array<OneD, NekDouble> &outarray);
        
        void DefineForcingTerm(void);
        Array<OneD, Array<OneD, NekDouble> > m_ForcingTerm;
        Array<OneD, Array<OneD, NekDouble> > m_ForcingTerm_Coeffs;
        
        Array<OneD, CoupledLocalToGlobalC0ContMapSharedPtr> m_locToGloMap;

        CoupledLinearNS_trafoP(const LibUtilities::SessionReaderSharedPtr &pSesssion);
        
        void DoInitialiseAdv(Array<OneD, NekDouble> myAdvField_x, Array<OneD, NekDouble> myAdvField_y);
        Eigen::MatrixXd DoTrafo(Array<OneD, Array<OneD, NekDouble> > snapshot_x_collection, Array<OneD, Array<OneD, NekDouble> > snapshot_y_collection, Array<OneD, NekDouble> param_vector);
	Array<OneD, Array<OneD, NekDouble> > DoSolve_at_param(Array<OneD, NekDouble> init_snapshot_x, Array<OneD, NekDouble> init_snapshot_y, NekDouble parameter);

	Eigen::VectorXd curr_f_bnd;
	Eigen::VectorXd curr_f_p;
	Eigen::VectorXd curr_f_int;

	NekDouble Get_m_kinvis(void);
	void Set_m_kinvis(NekDouble);
	int use_Newton;
	int snapshot_computation_plot_rel_errors;
	int debug_mode;
	Array<OneD, Array<OneD, NekDouble> > myAdvField;

    protected:

        
        virtual void v_InitObject();
        
    private:
        /// Id to identify when single mode is mean mode (i.e. beta=0);
        bool m_zeroMode;
        
        int m_counter;
        bool m_initialStep;
        NekDouble   m_tol;        // Tolerence
        int m_maxIt;           // Max number of iteration
        int m_Restart;    // 0=Stokes solution as init guess; 1=Restart.cont as init guess
        int m_MatrixSetUpStep; 
        NekDouble m_kinvisMin;
        NekDouble m_kinvisStep;
        NekDouble m_KinvisPercentage;
        
        
        
        
        Array<OneD, CoupledSolverMatrices> m_mat;
        
        
        /**
         *  Generate the linearised Navier Stokes solver based on the
         *  static condensation of the interior velocity space and
         *  pressure modes. This call also allows for a Fourier mode
         *  to be specified, however if HomogeneneousMode= 0 then can
         *  be used for a standared 2D and hopefully 3D (in the
         *  future).
         */
        void SetUpCoupledMatrix(const NekDouble lambda, 
                                const Array< OneD, Array<OneD, NekDouble> > &Advfield, 
                                bool       IsLinearNSEquation,
                                const int  HomogeneousMode,
                                CoupledSolverMatrices &mat,
                                CoupledLocalToGlobalC0ContMapSharedPtr &locToGloMap,
                                const NekDouble lambda_imag = NekConstants::kNekUnsetDouble);
        
        virtual void v_GenerateSummary(SolverUtils::SummaryList& s);
        
        virtual void v_DoInitialise(void);
        
        virtual void v_DoSolve(void);
        
        virtual bool v_NegatedOp(void);
        
        virtual void v_TransCoeffToPhys(void);
        
        virtual void v_TransPhysToCoeff(void);
        
        virtual void v_Output(void);
        
        virtual int v_GetForceDimension(void);
    };
    
    
    
} //end of namespace


/**
 * $Log:$
 *
 **/

