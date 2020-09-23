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
#include "../Eigen/Dense"
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
//#include <MultiRegions/GlobalLinSysDirectStaticCond.h>

namespace Nektar
{     
    
      
    class CoupledLinearNS_TT: public CoupledLinearNS
    {
    public:
        friend class MemoryManager<CoupledLinearNS_TT>;
        
        /// Creates an instance of this class
        static SolverUtils::EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        {
            SolverUtils::EquationSystemSharedPtr p = MemoryManager<CoupledLinearNS_TT>::AllocateSharedPtr(pSession);
            p->InitObject();
            return p;
        }
        /// Name of class
        static std::string className;      

        void Solve(void);

        void SolveLinearNS(const Array<OneD, Array<OneD, NekDouble> > &forcing);
        
        void SolveLinearNS(const Array<OneD, Array<OneD, NekDouble> > &forcing,
                           Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                           MultiRegions::ExpListSharedPtr &pressure,
                           const int HomogeneousMode = 0);

        void SetUpCoupledMatrix(const NekDouble lambda = 0.0, 
                                const Array< OneD, Array<OneD, NekDouble>  > &Advfield = NullNekDoubleArrayofArray, 
                                bool IsLinearNSEquation = true);

        Array<OneD, CoupledLocalToGlobalC0ContMapSharedPtr> m_locToGloMap;

	void DefineRBspace(Eigen::MatrixXd RB_via_POD);

	Eigen::VectorXd trafoSnapshot(Eigen::VectorXd RB_via_POD, double kInvis);

	double Geo_T(double w, int elemT, int index); // array of matrices not suitable since there is no way to be symbolic

	void trafoSnapshot_simple(Eigen::MatrixXd RB_via_POD);

        Eigen::MatrixXd DoTrafo(Array<OneD, Array<OneD, NekDouble> > snapshot_x_collection, Array<OneD, Array<OneD, NekDouble> > snapshot_y_collection, Array<OneD, NekDouble> param_vector);

	void load_snapshots(int number_of_snapshots);
	void load_snapshots_geometry_params(int );
	void load_snapshots_geometry_params_conv_Oseen(int );
	void compute_snapshots(int number_of_snapshots);
	void compute_snapshots_geometry_params();
	void do_geo_trafo();
	void write_curr_field(std::string filename);
        
	int parameter_space_dimension;
	bool load_cO_snapshot_data_from_files;
	bool do_trafo_check;
	double POD_tolerance;
	double start_param_dir0;
	double end_param_dir0;
	double start_param_dir1;
	double end_param_dir1;
	int use_fine_grid_VV;
	int use_fine_grid_VV_and_load_ref;
	int use_fine_grid_VV_random;
	int use_sparse_poly;

        Array<OneD, Array<OneD, NekDouble> > m_ForcingTerm;
        Array<OneD, Array<OneD, NekDouble> > m_ForcingTerm_Coeffs;

        Eigen::MatrixXd PODmodes;
	Eigen::MatrixXd collect_f_all;
        Eigen::MatrixXd RB;
        Eigen::MatrixXd M_RB;
        Eigen::MatrixXd M_collect_f_all;
	Eigen::MatrixXd Get_no_advection_matrix(void);
	Eigen::MatrixXd Get_no_advection_matrix_ABCD(void);
	Eigen::MatrixXd Get_no_advection_matrix_pressure(void);
	Eigen::MatrixXd Get_advection_matrix(void);
	Eigen::MatrixXd Get_complete_matrix(void);
	Eigen::MatrixXd curr_xy_projected;

	Array<OneD, Array<OneD, double> > PhysBaseVec_x;
	Array<OneD, Array<OneD, double> > PhysBaseVec_y;
	Eigen::MatrixXd eigen_phys_basis_x;
	Eigen::MatrixXd eigen_phys_basis_y;
	Array<OneD, Array<OneD, double> > orth_PhysBaseVec_x;
	Array<OneD, Array<OneD, double> > orth_PhysBaseVec_y;
	Array<OneD, std::set<int> > elements_trafo;
	int number_elem_trafo;
	int qoi_dof;
	int use_LocROM;
	int only_single_cluster;
	int which_single_cluster;
	int load_predef_cluster;
	int fine_grid_dir0;
	int fine_grid_dir1;
	double L2norm_ITHACA( Array< OneD, NekDouble > component_x, Array< OneD, NekDouble > component_y );
	double Linfnorm_ITHACA( Array< OneD, NekDouble > component_x, Array< OneD, NekDouble > component_y );
	double L2norm_abs_error_ITHACA( Array< OneD, NekDouble > component1_x, Array< OneD, NekDouble > component1_y, Array< OneD, NekDouble > component2_x, Array< OneD, NekDouble > component2_y );
	double Linfnorm_abs_error_ITHACA( Array< OneD, NekDouble > component1_x, Array< OneD, NekDouble > component1_y, Array< OneD, NekDouble > component2_x, Array< OneD, NekDouble > component2_y );
	void k_means_ITHACA(int no_clusters, Array<OneD, std::set<int> > &clusters, double &CVT_energy);
	void evaluate_local_clusters(Array<OneD, std::set<int> > optimal_clusters);
        void run_local_ROM_offline(Eigen::MatrixXd collect_f_all);
        void run_local_ROM_offline_add_transition(Eigen::MatrixXd , Eigen::MatrixXd, int);
	void run_local_ROM_online(std::set<int>, int );
	void associate_VV_to_clusters(Array<OneD, std::set<int> >);

//	Array<OneD, Eigen::Matrix2d > elements_trafo_matrix; // put this as a function or find a way with symbolic computation
        void gen_phys_base_vecs();

	Array<OneD, Eigen::MatrixXd> adv_mats_proj_x;
	Array<OneD, Eigen::MatrixXd> adv_mats_proj_y;
	Array<OneD, Eigen::VectorXd> adv_vec_proj_x;
	Array<OneD, Eigen::VectorXd> adv_vec_proj_y;
	Array<OneD, Eigen::VectorXd> adv_vec_proj_x_newton;
	Array<OneD, Eigen::VectorXd> adv_vec_proj_y_newton;
	Array<OneD, Eigen::MatrixXd> adv_vec_proj_x_newton_RB;
	Array<OneD, Eigen::MatrixXd> adv_vec_proj_y_newton_RB;
	Array<OneD, Array<OneD, Array<OneD, Eigen::MatrixXd > > > adv_mats_proj_x_2d;
	Array<OneD, Array<OneD, Array<OneD, Eigen::MatrixXd > > > adv_mats_proj_y_2d;
	Array<OneD, Array<OneD, Array<OneD, Eigen::VectorXd > > > adv_vec_proj_x_2d;
	Array<OneD, Array<OneD, Array<OneD, Eigen::VectorXd > > > adv_vec_proj_y_2d;


	void gen_proj_adv_terms();
	void gen_proj_adv_terms_2d();
	void online_snapshot_check_with_smaller_basis(int);
	void online_snapshot_check_with_smaller_basis_VV(int);
	Eigen::MatrixXd gen_no_advection_matrix_pressure();
	Eigen::MatrixXd gen_no_advection_matrix_ABCD();

	Eigen::MatrixXd gen_adv_mats_proj_x(Array<OneD, double>, int);
	Eigen::MatrixXd gen_adv_mats_proj_y(Array<OneD, double>, int);
	Array<OneD, Array<OneD, Eigen::MatrixXd > > gen_adv_mats_proj_x_2d(Array<OneD, double>, Array<OneD, Array<OneD, Eigen::VectorXd > > &adv_vec_proj_x_2d);
	Array<OneD, Array<OneD, Eigen::MatrixXd > > gen_adv_mats_proj_y_2d(Array<OneD, double>, Array<OneD, Array<OneD, Eigen::VectorXd > > &adv_vec_proj_y_2d);

	double recover_snapshot_data(Eigen::VectorXd, int);
	void recover_snapshot_loop(Eigen::VectorXd, Array<OneD, double> &, Array<OneD, double> &);

	Eigen::MatrixXd ANN_POD_coeffs;
	void offline_phase();
	void online_phase();
	void compute_sparse_poly_approx();
	void compute_ANN_approx();
	Array<OneD, NekDouble> param_point;
	Array<OneD, Array<OneD, NekDouble> > general_param_vector;
	Array<OneD, Array<OneD, NekDouble> > fine_general_param_vector;
	int find_closest_snapshot_location(Array<OneD, NekDouble>, Array<OneD, Array<OneD, NekDouble> >);
	int find_closest_snapshot_location_linf(Array<OneD, NekDouble>, Array<OneD, Array<OneD, NekDouble> >);
	int find_closest_snapshot_location_l1(Array<OneD, NekDouble>, Array<OneD, Array<OneD, NekDouble> >);
	Array<OneD, NekDouble> param_vector;
	int Nmax;
	int RBsize;
	int number_of_snapshots_dir0;
	int number_of_snapshots_dir1;
	int globally_connected;
	int use_Newton;
	int use_non_unique_up_to_two;
	bool debug_mode;
	int use_overlap_p_space;
	int write_ROM_field;
	int snapshot_computation_plot_rel_errors;
	int compute_smaller_model_errs;

        Eigen::MatrixXd MtM;
        Eigen::MatrixXd Mtrafo;
        Eigen::MatrixXd RB_A;
        Eigen::MatrixXd RB_A_adv;
        Eigen::MatrixXd RB_A_no_adv;
        Eigen::MatrixXd RB_B;
        Eigen::MatrixXd RB_B_adv;
        Eigen::MatrixXd RB_B_no_adv;
        Eigen::MatrixXd RB_C;
        Eigen::MatrixXd RB_C_adv;
        Eigen::MatrixXd RB_C_no_adv;
        Eigen::MatrixXd RB_D;
        Eigen::MatrixXd RB_D_adv;
        Eigen::MatrixXd RB_D_no_adv;
        Eigen::MatrixXd RB_Dbnd;
        Eigen::MatrixXd RB_Dint;

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
	Eigen::MatrixXd gen_affine_mat_proj(double);
	Eigen::VectorXd gen_affine_vec_proj(double, int);
	Eigen::MatrixXd gen_affine_mat_proj_2d(double, double);
	Eigen::VectorXd gen_affine_vec_proj_2d(double, double, int);

	Eigen::MatrixXd reproject_from_basis( Eigen::MatrixXd curr_xy_proj );

	Array<OneD, Array<OneD, Eigen::MatrixXd > > the_const_one_proj_2d;
	Array<OneD, Array<OneD, Eigen::MatrixXd > > the_ABCD_one_proj_2d;
	Array<OneD, Array<OneD, Eigen::VectorXd > > the_ABCD_one_rhs_proj_2d;
	Array<OneD, Array<OneD, Eigen::VectorXd > > the_const_one_rhs_proj_2d;

	int no_dbc_in_loc;
	int no_not_dbc_in_loc;
	std::set<int> elem_loc_dbc;   // works for all globally connected scenarios
	std::set<int> elem_not_loc_dbc;
	int M_no_dbc_in_loc;
	int M_no_not_dbc_in_loc;
	std::set<int> M_elem_loc_dbc;
	std::set<int> M_elem_not_loc_dbc;
	Eigen::VectorXd f_bnd_dbc;
	Eigen::VectorXd f_bnd_dbc_full_size;
	Eigen::VectorXd M_f_bnd_dbc;
	Eigen::VectorXd M_f_bnd_dbc_full_size;
	int M_truth_size;               // works for all globally connected scenarios
	int M_truth_size_without_DBC;   // works for all globally connected scenarios
	int nBndDofs;
	int f_bnd_size;
	int f_p_size;
	int f_int_size;

	Eigen::VectorXd curr_f_bnd;
	Eigen::VectorXd curr_f_p;
	Eigen::VectorXd curr_f_int;

	Array<OneD, Array<OneD, NekDouble> > snapshot_x_collection;
	Array<OneD, Array<OneD, NekDouble> > snapshot_y_collection;
	Array<OneD, Array<OneD, NekDouble> > snapshot_x_collection_VV;
	Array<OneD, Array<OneD, NekDouble> > snapshot_y_collection_VV;

        
	CoupledLinearNS_TT(const LibUtilities::SessionReaderSharedPtr &pSesssion);

	double ref_param_nu;
	int ref_param_index;
	Eigen::MatrixXd adv_geo_mat_projector(Array<OneD, Array<OneD, Eigen::MatrixXd > > Ah_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > B_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > C_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > D_elem, int, int, Eigen::VectorXd &adv_vec_proj );
	Eigen::MatrixXd press_geo_mat_projector(Array<OneD, Array<OneD, Eigen::MatrixXd > >, Array<OneD, Array<OneD, Eigen::MatrixXd > >, int curr_elem_trafo, int deriv_index, Eigen::VectorXd &press_vec_proj);
	Eigen::MatrixXd ABCD_geo_mat_projector(Array<OneD, Array<OneD, Eigen::MatrixXd > > A_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > B_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > C_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > D_elem, int curr_elem_trafo, int deriv_index, Eigen::VectorXd &ABCD_vec_proj);

	void gen_reference_matrices();
	void gen_reference_matrices_2d();
	Eigen::VectorXd reconstruct_solution_w_dbc(Eigen::VectorXd);
        void setDBC(Eigen::MatrixXd collect_f_all);
	void setDBC_M(Eigen::MatrixXd collect_f_all);
	Eigen::MatrixXd project_onto_basis(Array<OneD, NekDouble> snapshot_x, Array<OneD, NekDouble> snapshot_y);
	Array<OneD, Array<OneD, NekDouble> > trafo_current_para(Array<OneD, NekDouble>, Array<OneD, NekDouble>, Array<OneD, NekDouble>, Eigen::VectorXd &, Eigen::VectorXd &, Eigen::VectorXd &);
	int get_curr_elem_pos(int);

        void set_MtM();
        void DoInitialiseAdv(Array<OneD, NekDouble> myAdvField_x, Array<OneD, NekDouble> myAdvField_y);
        Eigen::MatrixXd remove_cols_and_rows(Eigen::MatrixXd the_matrix, std::set<int> elements_to_be_removed);
        Eigen::VectorXd remove_rows(Eigen::VectorXd the_vector, std::set<int> elements_to_be_removed);

	NekDouble Get_m_kinvis(void);
	double lagrange_interp(double curr_param, int curr_index, int sparse_poly_approx_dimension);
	void sparse_approx_VV(int sparse_poly_approx_dimension, double& max, double& mean);
	void Set_m_kinvis(NekDouble);
	Array<OneD, Array<OneD, NekDouble> > myAdvField_Newton;

	int max_sparse_poly_approx_dimension;
	bool use_ANN;
	bool use_ANN_local;
	Eigen::MatrixXd train_data_x;
	Eigen::MatrixXd train_data_y;

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

        void SetUpCoupledMatrix(const NekDouble lambda, 
                                const Array< OneD, Array<OneD, NekDouble> > &Advfield, 
                                bool       IsLinearNSEquation,
                                const int  HomogeneousMode,
                                CoupledSolverMatrices &mat,
                                CoupledLocalToGlobalC0ContMapSharedPtr &locToGloMap,
                                const NekDouble lambda_imag = NekConstants::kNekUnsetDouble);

        virtual void v_DoSolve(void);        

        virtual void v_DoInitialise(void);

        virtual void v_Output(void);

    };
    
    
    
} //end of namespace


/**
 * $Log:$
 *
 **/

