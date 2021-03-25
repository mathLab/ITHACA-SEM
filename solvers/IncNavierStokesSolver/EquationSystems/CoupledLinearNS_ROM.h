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
	
	
	inline double Geo_T(double w, int elemT, int index)
	{
		if (Geo_trafo_load == 1)
			return Geo_T_no1(w,  elemT,  index);
		else if (Geo_trafo_load == 2)
			return Geo_T_no2(w,  elemT,  index); 
		else if (Geo_trafo_load == 3)
			return Geo_T_no3(w,  elemT,  index); 
		else if (Geo_trafo_load == 4)
			return Geo_T_no4(w,  elemT,  index); 
		else if (Geo_trafo_load == 5)
			return Geo_T_no5(w,  elemT,  index);
		return -1; 
	}
	
	inline double Geo_T_no1(double w, int elemT, int index)
	{
		Eigen::Matrix2d T;
	if (elemT == 0) 
	{
		T << 1, 0, 0, 1/w;
	}
	else if (elemT == 1) 
	{
		T << 1, 0, 0, 2/(3-w);
	}
	else if (elemT == 2) 
	{
		T << 1, 0, -(1-w)/2, 1;
	}
	else if (elemT == 3) 
	{
		T << 1, 0, -(w-1)/2, 1;
	}
	else if (elemT == 4) 
	{
		T << 1, 0, 0, 1;
	}

//	cout << T << endl;

	if (index == 0) 
	{
		return 1/(T(0,0)*T(1,1) - T(0,1)*T(1,0)); // 1/det
	}
	else if (index == 1) 
	{
		return T(0,0);
	}
	else if (index == 2) 
	{
		return T(0,1);
	}
	else if (index == 3) 
	{
		return T(1,0);
	}
	else if (index == 4) 
	{
		return T(1,1);
	}

	return 0;


    }
    
    inline Array<OneD, std::set<int> > set_elem_trafo(int number_elem_trafo)
    {
    		Array<OneD, std::set<int> > elements_trafo = Array<OneD, std::set<int> > (number_elem_trafo);
		elements_trafo[0].insert(32); 
		elements_trafo[0].insert(30);
		elements_trafo[0].insert(11);
		elements_trafo[0].insert(10);
		elements_trafo[0].insert(9);
		elements_trafo[0].insert(8);
		elements_trafo[0].insert(17);
		elements_trafo[0].insert(16);
		elements_trafo[0].insert(15);
		elements_trafo[0].insert(14);
		elements_trafo[0].insert(25);
		elements_trafo[0].insert(24);
		elements_trafo[0].insert(23);
		elements_trafo[0].insert(22);

		elements_trafo[1].insert(12);
		elements_trafo[1].insert(28);
		elements_trafo[1].insert(34);
		elements_trafo[1].insert(13);
		elements_trafo[1].insert(21);
		elements_trafo[1].insert(20);
		elements_trafo[1].insert(18);
		elements_trafo[1].insert(19);
		elements_trafo[1].insert(26);
		elements_trafo[1].insert(27);

		elements_trafo[2].insert(29);
		elements_trafo[2].insert(31);

		elements_trafo[3].insert(33);
		elements_trafo[3].insert(35);

		elements_trafo[4].insert(0);
		elements_trafo[4].insert(1);
		elements_trafo[4].insert(2);
		elements_trafo[4].insert(3);
		elements_trafo[4].insert(4);
		elements_trafo[4].insert(5);
		elements_trafo[4].insert(6);
		elements_trafo[4].insert(7);
		
		return elements_trafo;
    }

	inline double Geo_T_no2(double w, int elemT, int index)
	{
		Eigen::Matrix2d T;
		T << 1, 0, 0, 1;
		if (index == 0) 
		{
			return 1/(T(0,0)*T(1,1) - T(0,1)*T(1,0)); // 1/det
		}
		else if (index == 1) 
		{
			return T(0,0);
		}
		else if (index == 2) 
		{
			return T(0,1);
		}
		else if (index == 3) 
		{
			return T(1,0);
		}
		else if (index == 4) 
		{
			return T(1,1);
		}
		return 0;
	    }

	inline double Geo_T_no3(double w, int elemT, int index)
	{
		// placeholder
		Eigen::Matrix2d T;
		T << 1, 0, 0, 1;
		if (index == 0) 
		{
			return 1/(T(0,0)*T(1,1) - T(0,1)*T(1,0)); // 1/det
		}
		else if (index == 1) 
		{
			return T(0,0);
		}
		else if (index == 2) 
		{
			return T(0,1);
		}
		else if (index == 3) 
		{
			return T(1,0);
		}
		else if (index == 4) 
		{
			return T(1,1);
		}
		return 0;
	    }
	    
	inline double Geo_T_no4(double w, int elemT, int index)
	{
		// placeholder
		Eigen::Matrix2d T;
		T << 1, 0, 0, 1;
		if (index == 0) 
		{
			return 1/(T(0,0)*T(1,1) - T(0,1)*T(1,0)); // 1/det
		}
		else if (index == 1) 
		{
			return T(0,0);
		}
		else if (index == 2) 
		{
			return T(0,1);
		}
		else if (index == 3) 
		{
			return T(1,0);
		}
		else if (index == 4) 
		{
			return T(1,1);
		}
		return 0;
	    }
	    
    	inline double Geo_T_no5(double w, int elemT, int index)
	{
		// placeholder
		Eigen::Matrix2d T;
		T << 1, 0, 0, 1;
		if (index == 0) 
		{
			return 1/(T(0,0)*T(1,1) - T(0,1)*T(1,0)); // 1/det
		}
		else if (index == 1) 
		{
			return T(0,0);
		}
		else if (index == 2) 
		{
			return T(0,1);
		}
		else if (index == 3) 
		{
			return T(1,0);
		}
		else if (index == 4) 
		{
			return T(1,1);
		}
		return 0;
	    }
    
    inline Array<OneD, std::set<int> > set_elem_trafo_no2(int number_elem_trafo)
    {
    		Array<OneD, std::set<int> > elements_trafo = Array<OneD, std::set<int> > (number_elem_trafo);
    		for (int i=0; i<36; ++i)
    		{
    			elements_trafo[0].insert(i); 
		}
		return elements_trafo;
    }


    inline Array<OneD, std::set<int> > set_elem_trafo_no3(int number_elem_trafo)
    {
    		Array<OneD, std::set<int> > elements_trafo = Array<OneD, std::set<int> > (number_elem_trafo);

		// placeholder

		return elements_trafo;
    }

    inline Array<OneD, std::set<int> > set_elem_trafo_no4(int number_elem_trafo)
    {
    		Array<OneD, std::set<int> > elements_trafo = Array<OneD, std::set<int> > (number_elem_trafo);

		// placeholder

		return elements_trafo;
    }


    inline Array<OneD, std::set<int> > set_elem_trafo_no5(int number_elem_trafo)
    {
    		Array<OneD, std::set<int> > elements_trafo = Array<OneD, std::set<int> > (number_elem_trafo);

		// placeholder

		return elements_trafo;
    }


    
	void write_curr_field(std::string filename);
	Array<OneD, Array<OneD, Eigen::MatrixXd > > gen_adv_mats_proj_x_geo(Array<OneD, double> curr_PhysBaseVec_x, Array<OneD, Array<OneD, Eigen::VectorXd > > &adv_vec_proj_x_2d);
	Array<OneD, Array<OneD, Eigen::MatrixXd > > gen_adv_mats_proj_y_geo(Array<OneD, double> curr_PhysBaseVec_y, Array<OneD, Array<OneD, Eigen::VectorXd > > &adv_vec_proj_y_2d);
        Eigen::MatrixXd adv_geo_mat_projector(Array<OneD, Array<OneD, Eigen::MatrixXd > > Ah_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > B_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > C_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > D_elem, int curr_elem_trafo, int deriv_index, Eigen::VectorXd &adv_vec_proj);	
        void gen_proj_adv_terms_geo();
	void gen_reference_matrices_geo();
     Eigen::MatrixXd press_geo_mat_projector(Array<OneD, Array<OneD, Eigen::MatrixXd > > Dbnd_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > Dint_elem, int curr_elem_trafo, int deriv_index, Eigen::VectorXd &press_vec_proj);
         Eigen::MatrixXd ABCD_geo_mat_projector(Array<OneD, Array<OneD, Eigen::MatrixXd > > A_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > B_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > C_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > D_elem, int curr_elem_trafo, int deriv_index, Eigen::VectorXd &ABCD_vec_proj);
             Eigen::MatrixXd gen_affine_mat_proj_geo(double current_nu, double w);
    Eigen::VectorXd gen_affine_vec_proj_geo(double current_nu, double w, int current_index);
    void recover_snapshot_loop(Eigen::VectorXd reconstruct_solution, Array<OneD, double> & field_x, Array<OneD, double> & field_y);
    double Linfnorm_ITHACA( Array< OneD, NekDouble > component_x, Array< OneD, NekDouble > component_y );
    double L2norm_ITHACA( Array< OneD, NekDouble > component_x, Array< OneD, NekDouble > component_y );
    double Linfnorm_abs_error_ITHACA( Array< OneD, NekDouble > component1_x, Array< OneD, NekDouble > component1_y, Array< OneD, NekDouble > component2_x, Array< OneD, NekDouble > component2_y );
    double L2norm_abs_error_ITHACA( Array< OneD, NekDouble > component1_x, Array< OneD, NekDouble > component1_y, Array< OneD, NekDouble > component2_x, Array< OneD, NekDouble > component2_y );
    void online_snapshot_check_with_smaller_basis(int);
    double lagrange_interp(double curr_param, int curr_index, int sparse_poly_approx_dimension);
    void compute_sparse_poly_approx();
            
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
        int Nmax_VV;
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
       	Array<OneD, Array<OneD, Array<OneD, Eigen::VectorXd > > > adv_vec_proj_x_2d;
	Array<OneD, Array<OneD, Array<OneD, Eigen::VectorXd > > > adv_vec_proj_y_2d;
	Array<OneD, Array<OneD, Array<OneD, Eigen::MatrixXd > > > adv_mats_proj_x_2d;
	Array<OneD, Array<OneD, Array<OneD, Eigen::MatrixXd > > > adv_mats_proj_y_2d;	
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
	double do_trafo_check_relative_error;
	bool load_cO_snapshot_data_from_files;
	bool compute_smaller_model_errs;
	int qoi_dof;
	int fine_grid_dir0;
	int fine_grid_dir1;
	Array<OneD, std::set<int> > elements_trafo;
	int number_elem_trafo;
	bool use_fine_grid_VV;
	Array<OneD, Array<OneD, Eigen::MatrixXd > > the_const_one_proj_2d;
	Array<OneD, Array<OneD, Eigen::MatrixXd > > the_ABCD_one_proj_2d;
	Array<OneD, Array<OneD, Eigen::VectorXd > > the_ABCD_one_rhs_proj_2d;
	Array<OneD, Array<OneD, Eigen::VectorXd > > the_const_one_rhs_proj_2d;
	bool use_fine_grid_VV_random;
	bool use_fine_grid_VV_and_load_ref;
	bool use_sparse_poly;
	int max_sparse_poly_approx_dimension;
	int Geo_trafo_load;
	bool replace_snapshot_with_transformed;

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
