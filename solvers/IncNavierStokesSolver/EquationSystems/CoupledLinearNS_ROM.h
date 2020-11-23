///////////////////////////////////////////////////////////////////////////////
//
// File CoupledLinearNS_ROM.h
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
// Description: Coupled Stokes solver scheme header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COUPLEDSTOKESSCHEME_H
#define NEKTAR_SOLVERS_COUPLEDSTOKESSCHEME_H

#include <IncNavierStokesSolver/EquationSystems/CoupledLocalToGlobalC0ContMap.h>
#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include "ThirdParty/Eigen/Dense"
//#include <MultiRegions/GlobalLinSysDirectStaticCond.h>

namespace Nektar
{     
    
    typedef struct coupledSolverMatrices
    {
        /**  \brief Boundary-interior Laplacian plus linearised convective
         *             terms pre-multiplying Cinv: 
         *             \f$ m\_BCinv[n,m] = (\nabla \phi^b_n, \nu \nabla
         *             \phi^i_j) + (\phi^b_n,{\bf U \cdot \nabla} \phi^i_j) +
         *             (\phi^b_n \nabla^T {\bf U} \phi^i_j)  m_Cinv[j,m]\f$  
         */
        DNekScalBlkMatSharedPtr m_BCinv;
        
        /** \brief Interior-boundary Laplacian plus linearised convective terms
         * \f$ m\_Btilde^T[n,m] = (\nabla \phi^i_n, \nu \nabla \phi^b_m) +
         * (\phi^i_n,{\bf U \cdot \nabla} \phi^b_m) + (\phi^i_n \nabla^T
         * {\bf U} \phi^b_m) \f$ */
        DNekScalBlkMatSharedPtr m_Btilde;
        
        /** \brief Interior-Interior Laplaican plus linearised convective
         *   terms inverted, i.e. the inverse of 
         *   \f$ m\_C[n,m] = (\nabla \phi^i_n, \nu \nabla
         *   \phi^i_m) + (\phi^i_n,{\bf U \cdot \nabla} \phi^i_m) +
         *   (\phi^i_n \nabla^T {\bf U} \phi^i_m),\f$ */
        DNekScalBlkMatSharedPtr  m_Cinv; 
        
        /** \brief Inner product of pressure system with divergence of the
         *   boundary velocity space  
         *   \f$ m\_D\_{bnd}[n,m] = (\psi_n,\nabla \phi^b_m),\f$ 
         */
        DNekScalBlkMatSharedPtr  m_D_bnd; 
        
        /** \brief Inner product of pressure system with divergence of the
         *  interior velocity space  
         *  \f$ m\_D\_{int}[n,m] = (\psi_j,\nabla \phi^i_m) \f$
         */
        DNekScalBlkMatSharedPtr  m_D_int; 
        
        MultiRegions::GlobalLinSysSharedPtr m_CoupledBndSys;
    } CoupledSolverMatrices;
    
    class CoupledLinearNS_ROM: public IncNavierStokes
    {
    public:
        friend class MemoryManager<CoupledLinearNS_ROM>;
        
        /// Creates an instance of this class
        static SolverUtils::EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr &pGraph)
        {
            SolverUtils::EquationSystemSharedPtr p = MemoryManager<CoupledLinearNS_ROM>::AllocateSharedPtr(pSession, pGraph);
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
        
        
        const SpatialDomains::ExpansionInfoMap
            &GenPressureExp(const SpatialDomains::ExpansionInfoMap &VelExp);
        
        void Solve(void);
        
        /**
         *   Solve the coupled linear Navier-Stokes solve using matrix
         *   systems set up at construction.  The solution is stored
         *   in #m_velocity and #m_pressure.
         */
        void SolveLinearNS(const Array<OneD, Array<OneD, NekDouble> > &forcing);
        
        void SolveLinearNS(Array<OneD, Array<OneD, NekDouble> > &forcing,
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
        
        
	// ROM functions     
        void load_snapshots(void);
        void load_session_parameters(void);
        void ROM_offline_phase(void);
        void ROM_online_phase(void);
        Eigen::MatrixXd DoTrafo(void);
        Eigen::MatrixXd DoTrafo_1p_kinvis(void);
        Eigen::MatrixXd DoTrafo_1p_geo(void);
        void DoInitialiseAdv(Array<OneD, NekDouble> myAdvField_x, Array<OneD, NekDouble> myAdvField_y);
        void Set_m_kinvis(NekDouble);
        void setDBC(Eigen::MatrixXd);
        void set_f_bnd_dbc();
        void gen_phys_base_vecs();
        void gen_proj_adv_terms();
        void gen_reference_matrices();
	Eigen::MatrixXd gen_adv_mats_proj_x(Array<OneD, double> curr_PhysBaseVec_x, int use_Newton);
        Eigen::MatrixXd gen_adv_mats_proj_y(Array<OneD, double> curr_PhysBaseVec_y, int use_Newton);
        Eigen::VectorXd remove_rows(Eigen::VectorXd the_vector, std::set<int> elements_to_be_removed);
        Eigen::MatrixXd remove_cols_and_rows(Eigen::MatrixXd the_matrix, std::set<int> elements_to_be_removed);
        Eigen::MatrixXd gen_no_advection_matrix_pressure();
        Eigen::MatrixXd gen_no_advection_matrix_ABCD();
        Eigen::MatrixXd project_onto_basis( Array<OneD, NekDouble> snapshot_x, Array<OneD, NekDouble> snapshot_y);
        Eigen::MatrixXd reproject_from_basis( Eigen::MatrixXd curr_xy_proj );
        Eigen::MatrixXd gen_affine_mat_proj(double current_nu);
        Eigen::VectorXd gen_affine_vec_proj(double current_nu, int current_index);
        Eigen::VectorXd reconstruct_solution_w_dbc(Eigen::VectorXd reprojected_solve);
        void compute_snapshots_kinvis(void);
        Array<OneD, Array<OneD, NekDouble> > DoSolve_at_param(Array<OneD, NekDouble> init_snapshot_x, Array<OneD, NekDouble> init_snapshot_y, NekDouble parameter);
        Array<OneD, Array<OneD, NekDouble> > trafo_current_para(Array<OneD, NekDouble> snapshot_x, Array<OneD, NekDouble> snapshot_y, Array<OneD, NekDouble> parameter_of_interest, Eigen::VectorXd & ref_f_bnd, Eigen::VectorXd & ref_f_p, Eigen::VectorXd & ref_f_int);
	int get_curr_elem_pos(int);   
	double Geo_T(double w, int elemT, int index);
	void do_geo_trafo();
	void write_curr_field(std::string filename);
	
    
    
            
//	int get_curr_elem_pos(int curr_elem);




        
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
        
    protected:
        CoupledLinearNS_ROM(const LibUtilities::SessionReaderSharedPtr &pSesssion,
                        const SpatialDomains::MeshGraphSharedPtr &pGraph);
        
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
        
        // ROM variables
        bool load_snapshot_data_from_files;
        int number_of_snapshots;
        int Nmax;
        int parameter_space_dimension;
        Array<OneD, int> parameter_types;
        double POD_tolerance;
	Array<OneD, NekDouble> param_vector;   
	Array<OneD, Array<OneD, NekDouble> > snapshot_x_collection;
	Array<OneD, Array<OneD, NekDouble> > snapshot_y_collection;
	Array<OneD, Array<OneD, NekDouble> > snapshot_x_collection_VV;
	Array<OneD, Array<OneD, NekDouble> > snapshot_y_collection_VV;
	int f_bnd_size;
	int f_p_size;
	int f_int_size;
	Eigen::VectorXd curr_f_bnd;
	Eigen::VectorXd curr_f_p;
	Eigen::VectorXd curr_f_int;	
	Eigen::MatrixXd collect_f_all;
	Array<OneD, Array<OneD, NekDouble> > myAdvField;
	bool ROM_started;
	bool ongoing_snapshot_computation;	
	bool debug_mode;
	int RBsize;
	int no_dbc_in_loc;
	int no_not_dbc_in_loc;
	std::set<int> elem_loc_dbc;   // works for all globally connected scenarios
	std::set<int> elem_not_loc_dbc;
        Eigen::MatrixXd PODmodes;
	int M_truth_size;               // works for all globally connected scenarios
	int M_truth_size_without_DBC;   // works for all globally connected scenarios
	int nBndDofs;        
	Eigen::VectorXd f_bnd_dbc;
	Eigen::VectorXd f_bnd_dbc_full_size;
        Eigen::MatrixXd RB;        
       	Eigen::MatrixXd eigen_phys_basis_x;
	Eigen::MatrixXd eigen_phys_basis_y;
	Array<OneD, Array<OneD, double> > orth_PhysBaseVec_x;
	Array<OneD, Array<OneD, double> > orth_PhysBaseVec_y;
        Array<OneD, Array<OneD, double> > PhysBaseVec_x;
	Array<OneD, Array<OneD, double> > PhysBaseVec_y;
	Array<OneD, Eigen::MatrixXd> adv_mats_proj_x;
	Array<OneD, Eigen::MatrixXd> adv_mats_proj_y;
	Array<OneD, Eigen::VectorXd> adv_vec_proj_x;
	Array<OneD, Eigen::VectorXd> adv_vec_proj_y;
       	Array<OneD, Eigen::VectorXd> adv_vec_proj_x_newton;
	Array<OneD, Eigen::VectorXd> adv_vec_proj_y_newton;
	Array<OneD, Eigen::MatrixXd> adv_vec_proj_x_newton_RB;
	Array<OneD, Eigen::MatrixXd> adv_vec_proj_y_newton_RB;
	bool use_Newton;        
	Eigen::MatrixXd the_const_one;
	Eigen::MatrixXd the_ABCD_one;
	Eigen::MatrixXd the_const_one_simplified;
	Eigen::MatrixXd the_ABCD_one_simplified;
	Eigen::MatrixXd the_const_one_proj;
	Eigen::MatrixXd the_ABCD_one_proj;
	Eigen::VectorXd the_ABCD_one_rhs;
	Eigen::VectorXd the_const_one_rhs;
	Eigen::VectorXd the_ABCD_one_rhs_simplified;
	Eigen::VectorXd the_const_one_rhs_simplified;
	Eigen::VectorXd the_ABCD_one_rhs_proj;
	Eigen::VectorXd the_const_one_rhs_proj;        
       	Eigen::MatrixXd curr_xy_projected;
	Array<OneD, Array<OneD, NekDouble> > myAdvField_Newton;
	Array<OneD, Array<OneD, NekDouble> > general_param_vector;
	Array<OneD, Array<OneD, NekDouble> > fine_general_param_vector;
	bool do_trafo_check;
	bool load_cO_snapshot_data_from_files;
	int qoi_dof;
	int fine_grid_dir0;
	int fine_grid_dir1;
	Array<OneD, std::set<int> > elements_trafo;
	int number_elem_trafo;




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

#endif //COUPLEDSTOKESSCHEME_H
