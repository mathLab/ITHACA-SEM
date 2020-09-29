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
        friend class Nektar::SolverUtils::EquationSystem;
        //friend class CoupledLinearNS_trafoP;
        
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

	void trafoSnapshot_simple(Eigen::MatrixXd RB_via_POD);

        Eigen::MatrixXd DoTrafo(Array<OneD, Array<OneD, NekDouble> > snapshot_x_collection, Array<OneD, Array<OneD, NekDouble> > snapshot_y_collection, Array<OneD, NekDouble> param_vector);

	void load_snapshots(int number_of_snapshots);
	void compute_snapshots(int number_of_snapshots);
        
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
        void gen_phys_base_vecs();

	Array<OneD, Eigen::MatrixXd> adv_mats_proj_x;
	Array<OneD, Eigen::MatrixXd> adv_mats_proj_y;
	Array<OneD, Eigen::VectorXd> adv_vec_proj_x;
	Array<OneD, Eigen::VectorXd> adv_vec_proj_y;
	Array<OneD, Eigen::VectorXd> adv_vec_proj_x_newton;
	Array<OneD, Eigen::VectorXd> adv_vec_proj_y_newton;
	Array<OneD, Eigen::MatrixXd> adv_vec_proj_x_newton_RB;
	Array<OneD, Eigen::MatrixXd> adv_vec_proj_y_newton_RB;

	void gen_proj_adv_terms();

	void recover_snapshot_data(Eigen::VectorXd, int);

	void offline_phase();
	void online_phase();
	Array<OneD, NekDouble> param_vector;
	int Nmax;
	int RBsize;
	int globally_connected;
	int use_Newton;
	int debug_mode;
	int write_ROM_field;
	int write_SEM_field;
	int snapshot_computation_plot_rel_errors;

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

	Eigen::VectorXd curr_f_bnd;
	Eigen::VectorXd curr_f_p;
	Eigen::VectorXd curr_f_int;

	Array<OneD, Array<OneD, NekDouble> > snapshot_x_collection;
	Array<OneD, Array<OneD, NekDouble> > snapshot_y_collection;
        
	CoupledLinearNS_TT(const LibUtilities::SessionReaderSharedPtr &pSesssion);

	double ref_param_nu;
	int ref_param_index;


	void gen_reference_matrices();
	Eigen::VectorXd reconstruct_solution_w_dbc(Eigen::VectorXd);
        void setDBC(Eigen::MatrixXd collect_f_all);
	void setDBC_M(Eigen::MatrixXd collect_f_all);
	Eigen::MatrixXd project_onto_basis(Array<OneD, NekDouble> snapshot_x, Array<OneD, NekDouble> snapshot_y);

        void set_MtM();
        void DoInitialiseAdv(Array<OneD, NekDouble> myAdvField_x, Array<OneD, NekDouble> myAdvField_y);
        Eigen::MatrixXd remove_cols_and_rows(Eigen::MatrixXd the_matrix, std::set<int> elements_to_be_removed);
        Eigen::VectorXd remove_rows(Eigen::VectorXd the_vector, std::set<int> elements_to_be_removed);

	NekDouble Get_m_kinvis(void);
	void Set_m_kinvis(NekDouble);
	
	
	Eigen::VectorXd gen_affine_vec(double, double, Eigen::VectorXd);
	
	
	std::vector< Array<OneD, double> > reproject_back(Eigen::VectorXd);
	double FarrelOutput(Eigen::VectorXd);
	Eigen::VectorXd reconstruct_solution_w_different_dbc(Eigen::VectorXd, double);
	void error_analysis(int, double, double, std::ofstream &);
	Eigen::MatrixXd load_collect_f_all();
	
	std::vector<int> flipperMap;
	int max_dimension;
	int continuation_from_files;
	int compare_accuracy_mode;
	int no_offline_files;
	Array<OneD, NekDouble> param_vector2;
	NekDouble offline_average_time, online_average_time;
	unsigned int online_no_solves;
	std::vector<Eigen::VectorXd> solve_affine;
	bool create_error_file;
	
	

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

