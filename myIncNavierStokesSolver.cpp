///////////////////////////////////////////////////////////////////////////////
//
// File IncNavierStokesSolver.cpp
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
// Description: Incompressible Navier Stokes solver
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Driver.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include "EquationSystems/VelocityCorrectionScheme.h"
#include "Eigen/Dense"
#include "Eigen/Core"
#include "EquationSystems/CoupledLinearNS_TT.h"
#include "EquationSystems/CoupledLinearNS_trafoP.h"

//#include <MultiRegions/ExpList2D.h>

using namespace std;
using namespace Nektar;
using namespace Nektar::SolverUtils;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr session;
    string vDriverModule;
    DriverSharedPtr drv;
  
    try
    {
        // Create session reader.
        session = LibUtilities::SessionReader::CreateInstance(argc, argv);
        session->LoadSolverInfo("Driver", vDriverModule, "Standard");
        // Create driver
        drv = GetDriverFactory().CreateInstance(vDriverModule, session); 

	// offline data acquisition
	int load_snapshot_data_from_files = session->GetParameter("load_snapshot_data_from_files");
	int number_of_snapshots = session->GetParameter("number_of_snapshots");
	int snapshots_to_be_collected_aka_Nmax = number_of_snapshots;  
	int Nmax = snapshots_to_be_collected_aka_Nmax;
	Array<OneD, NekDouble> param_vector(snapshots_to_be_collected_aka_Nmax);
	// = [0.1 0.5 1 10];
	param_vector[0] = 0.1; // should also wander to the xml file
	param_vector[1] = 0.5;
	param_vector[2] = 1;
	param_vector[3] = 10;
	
//	for (int i = 0; i<Nmax; i++)
//	{
//		cout << "param_vector[i] " << param_vector[i] << endl;   
//	}

	CoupledLinearNS_TT CLNS(session);
	CLNS.InitObject();


	CLNS.load_snapshots(number_of_snapshots);
	Array<OneD, Array<OneD, NekDouble> > snapshot_x_collection = CLNS.snapshot_x_collection;
	Array<OneD, Array<OneD, NekDouble> > snapshot_y_collection = CLNS.snapshot_y_collection;



	// trafo the snapshots and see if they have really converged
	// so I need to capture the bnd / p / int as well as the local / global / phys
	// Domanda: how to change the k_invis for an already declared class -- Risposta: use Getter and Setter functions (to be implemented - done)

	// Q: how to know the bnd / p / int sizes beforehand
	CLNS.DoInitialise(); 
//	cout << "after CLNS.DoInitialise(); " << endl;
	CLNS.DoSolve();
//	cout << "after CLNS.DoSolve(); " << endl;
//	Eigen::MatrixXd collect_f_all = CLNS.DoTrafo(CLNS.snapshot_x_collection, CLNS.snapshot_y_collection, param_vector);
	session->SetSolverInfo("SolverType", "CoupledLinearisedNS_trafoP");
	CoupledLinearNS_trafoP babyCLNS_trafo(session);
	babyCLNS_trafo.InitObject();
	Eigen::MatrixXd collect_f_all = babyCLNS_trafo.DoTrafo(snapshot_x_collection, snapshot_y_collection, param_vector);

	Eigen::BDCSVD<Eigen::MatrixXd> svd_collect_f_all(collect_f_all, Eigen::ComputeThinU);
	cout << "svd_collect_f_all.singularValues() " << svd_collect_f_all.singularValues() << endl << endl;
	Eigen::MatrixXd collect_f_all_PODmodes = svd_collect_f_all.matrixU();
	// here probably limit to something like 99.99 percent of PODenergy, this will set RBsize
	int RBsize = Nmax; // the case of using all

	CLNS.setDBC(collect_f_all);
	// identifiy the D-BC dofs --- can do this by checking the c_f_bnd dofs -> needs a func as well
	int no_dbc_in_loc = CLNS.no_dbc_in_loc;
	int no_not_dbc_in_loc = CLNS.no_not_dbc_in_loc;
	std::set<int> elem_loc_dbc = CLNS.elem_loc_dbc;
	set<int> elem_not_loc_dbc = CLNS.elem_not_loc_dbc;

	Array<OneD, MultiRegions::ExpListSharedPtr> m_fields = CLNS.UpdateFields();
        int  nel  = m_fields[0]->GetNumElmts(); // number of spectral elements

	CLNS.set_MtM();
	Eigen::MatrixXd MtM = CLNS.MtM;

	Eigen::MatrixXd c_f_all_PODmodes_wo_dbc(collect_f_all_PODmodes.rows() - no_dbc_in_loc, collect_f_all_PODmodes.cols());
	Eigen::VectorXd f_bnd_dbc(no_dbc_in_loc);
	Eigen::VectorXd f_bnd_dbc_full_size(collect_f_all_PODmodes.rows());
	int counter_all = 0;
	int counter_dbc = 0;
	for (int index=0; index < collect_f_all_PODmodes.rows(); ++index)
	{
		if (!elem_loc_dbc.count(index))
		{
			f_bnd_dbc_full_size(index) = 0;
			c_f_all_PODmodes_wo_dbc.row(counter_all) = collect_f_all_PODmodes.row(index);
			counter_all++;
		}
		else
		{
			f_bnd_dbc_full_size(index) = collect_f_all(index,0);
			f_bnd_dbc(counter_dbc) = collect_f_all(index,0);
			counter_dbc++;
		}
	}

	Array<OneD, double> PhysBase_zero(CLNS.GetNpoints(), 0.0);

	Array<OneD, Array<OneD, double> > PhysBaseVec_x(RBsize); // or do I want it as an Eigen-matrix
	Array<OneD, Array<OneD, double> > PhysBaseVec_y(RBsize);

	
	for (int curr_trafo_iter=0; curr_trafo_iter < RBsize; curr_trafo_iter++)
	{
		Eigen::VectorXd f_bnd = collect_f_all_PODmodes.block(0, curr_trafo_iter, CLNS.curr_f_bnd.size(), 1);
		Eigen::VectorXd f_int = collect_f_all_PODmodes.block(CLNS.curr_f_bnd.size()+CLNS.curr_f_p.size(), curr_trafo_iter, CLNS.curr_f_int.size(), 1);
//		cout << "f_int " << f_int << endl;
		Array<OneD, MultiRegions::ExpListSharedPtr> fields = CLNS.UpdateFields(); 
	        Array<OneD, unsigned int> bmap, imap; 
		Array<OneD, double> field_0(CLNS.GetNcoeffs());
		Array<OneD, double> field_1(CLNS.GetNcoeffs());
		Array<OneD, double> curr_PhysBaseVec_x(CLNS.GetNpoints(), 0.0);
		Array<OneD, double> curr_PhysBaseVec_y(CLNS.GetNpoints(), 0.0);
	        int cnt = 0;
		int cnt1 = 0;
		int nvel = 2;
	        int nz_loc = 1;
	        int  nplanecoeffs = fields[0]->GetNcoeffs();
	        for(int i = 0; i < nel; ++i) // loop over elements
	        {
	            int eid  = i;
	            fields[0]->GetExp(eid)->GetBoundaryMap(bmap);
	            fields[0]->GetExp(eid)->GetInteriorMap(imap);
	            int nbnd   = bmap.num_elements();
	            int nint   = imap.num_elements();
	            int offset = fields[0]->GetCoeff_Offset(eid);
	            
	            for(int j = 0; j < nvel; ++j) // loop over velocity fields 
	            {
	                for(int n = 0; n < nz_loc; ++n)
	                {
	                    for(int k = 0; k < nbnd; ++k)
	                    {
	                        fields[j]->SetCoeff(n*nplanecoeffs + offset+bmap[k], f_bnd(cnt+k)); // could also just go through some "custom" Eigen Matrix
				if (j==0)
				{
					field_0[n*nplanecoeffs + offset+bmap[k]] = f_bnd(cnt+k);
				}
				if (j==1)
				{
					field_1[n*nplanecoeffs + offset+bmap[k]] = f_bnd(cnt+k);
				}
	                    }
	                    
	                    for(int k = 0; k < nint; ++k)
	                    {
	                        fields[j]->SetCoeff(n*nplanecoeffs + offset+imap[k], f_int(cnt1+k));
				if (j==0)
				{
					field_0[n*nplanecoeffs + offset+imap[k]] = f_int(cnt1+k);
				}
				if (j==1)
				{
					field_1[n*nplanecoeffs + offset+imap[k]] = f_int(cnt1+k);
				}
	                    }
	                    cnt  += nbnd;
	                    cnt1 += nint;
	                }
	            }
	        }
		Array<OneD, double> test_nn = fields[0]->GetCoeffs();
		fields[0]->BwdTrans_IterPerExp(fields[0]->GetCoeffs(), curr_PhysBaseVec_x);
		fields[1]->BwdTrans_IterPerExp(fields[1]->GetCoeffs(), curr_PhysBaseVec_y);
		PhysBaseVec_x[curr_trafo_iter] = curr_PhysBaseVec_x;
		PhysBaseVec_y[curr_trafo_iter] = curr_PhysBaseVec_y;
		
	}

	Eigen::MatrixXd eigen_phys_basis_x(CLNS.GetNpoints(), RBsize);
	Eigen::MatrixXd eigen_phys_basis_y(CLNS.GetNpoints(), RBsize);
	for (int index_phys_base=0; index_phys_base<CLNS.GetNpoints(); index_phys_base++)
	{
		for (int index_RBsize=0; index_RBsize<RBsize; index_RBsize++)
		{
			eigen_phys_basis_x(index_phys_base,index_RBsize) = PhysBaseVec_x[index_RBsize][index_phys_base];
			eigen_phys_basis_y(index_phys_base,index_RBsize) = PhysBaseVec_y[index_RBsize][index_phys_base];
		}
	}
	Eigen::VectorXd curr_col = eigen_phys_basis_x.col(0);
	double norm_curr_col = curr_col.norm();
	eigen_phys_basis_x.col(0) = curr_col / norm_curr_col;
	curr_col = eigen_phys_basis_y.col(0);
	norm_curr_col = curr_col.norm();
	eigen_phys_basis_y.col(0) = curr_col / norm_curr_col;
	Eigen::VectorXd orthogonal_complement;
	
	for (int orth_iter=1; orth_iter<RBsize; orth_iter++)
	{
		curr_col = eigen_phys_basis_x.col(orth_iter);
		Eigen::MatrixXd leftmostCols = eigen_phys_basis_x.leftCols(orth_iter);
		orthogonal_complement = curr_col - leftmostCols * leftmostCols.transpose() * curr_col;
		norm_curr_col = orthogonal_complement.norm();
		eigen_phys_basis_x.col(orth_iter) = orthogonal_complement / norm_curr_col;
		curr_col = eigen_phys_basis_y.col(orth_iter);
		leftmostCols = eigen_phys_basis_y.leftCols(orth_iter);
		orthogonal_complement = curr_col - leftmostCols * leftmostCols.transpose() * curr_col;
		norm_curr_col = orthogonal_complement.norm();
		eigen_phys_basis_y.col(orth_iter) = orthogonal_complement / norm_curr_col;
	}


	Array<OneD, Array<OneD, double> > orth_PhysBaseVec_x(RBsize); // also have it as an Eigen-matrix: eigen_phys_basis_ x/y
	Array<OneD, Array<OneD, double> > orth_PhysBaseVec_y(RBsize);
	for (int index_RBsize=0; index_RBsize<RBsize; index_RBsize++)
	{
		Array<OneD, double> curr_iter_x(CLNS.GetNpoints());
		Array<OneD, double> curr_iter_y(CLNS.GetNpoints());

		for (int index_phys_base=0; index_phys_base<CLNS.GetNpoints(); index_phys_base++)	
		{
			curr_iter_x[index_phys_base] = eigen_phys_basis_x(index_phys_base,index_RBsize);
			curr_iter_y[index_phys_base] = eigen_phys_basis_y(index_phys_base,index_RBsize);			
		}
		orth_PhysBaseVec_x[index_RBsize] = curr_iter_x;
		orth_PhysBaseVec_y[index_RBsize] = curr_iter_y;			

	}



	Array<OneD, Eigen::MatrixXd> adv_mats_proj_x(RBsize);
	Array<OneD, Eigen::MatrixXd> adv_mats_proj_y(RBsize);
	Array<OneD, Eigen::VectorXd> adv_vec_proj_x(RBsize);
	Array<OneD, Eigen::VectorXd> adv_vec_proj_y(RBsize);


	for(int trafo_iter = 0; trafo_iter < RBsize; trafo_iter++)
	{
		Array<OneD, double> curr_PhysBaseVec_x = orth_PhysBaseVec_x[trafo_iter];
		Array<OneD, double> curr_PhysBaseVec_y = orth_PhysBaseVec_y[trafo_iter];
		
		CLNS.InitObject();
//		CLNS.DoInitialise();
		CLNS.DoInitialiseAdv(curr_PhysBaseVec_x, PhysBase_zero); // call with parameter in phys state
		Eigen::MatrixXd adv_matrix = Eigen::MatrixXd::Zero(CLNS.RB_A.rows() + CLNS.RB_Dbnd.rows() + CLNS.RB_C.cols(), CLNS.RB_A.cols() + CLNS.RB_Dbnd.rows() + CLNS.RB_B.cols() );
		adv_matrix = CLNS.Get_advection_matrix();
		
		// compute add_to rhs
		Eigen::VectorXd add_to_rhs_adv(collect_f_all_PODmodes.rows()); // probably need this for adv and non-adv
		add_to_rhs_adv = adv_matrix * f_bnd_dbc_full_size;   

		// also implement the removal of dbc dof
		// the set of to-be-removed cols and rows is elem_loc_dbc
		Eigen::MatrixXd adv_matrix_simplified = Eigen::MatrixXd::Zero(adv_matrix.rows() - no_dbc_in_loc, adv_matrix.cols() - no_dbc_in_loc);
		Eigen::VectorXd adv_rhs_add = Eigen::VectorXd::Zero(adv_matrix.rows() - no_dbc_in_loc);
		// a naive algorithm (should test the timing here on a "real-world" example)
		int counter_row_simplified = 0;
		int counter_col_simplified = 0;
		for (int row_index=0; row_index < f_bnd_dbc_full_size.rows(); ++row_index)
		{
			for (int col_index=0; col_index < f_bnd_dbc_full_size.rows(); ++col_index)
			{
				if ((!elem_loc_dbc.count(row_index)) && (!elem_loc_dbc.count(col_index)))
				{
					adv_matrix_simplified(counter_row_simplified, counter_col_simplified) = adv_matrix(row_index, col_index);
					counter_col_simplified++;
				}
				
			}
			counter_col_simplified = 0;
			if (!elem_loc_dbc.count(row_index))
			{
				adv_rhs_add(counter_row_simplified) = add_to_rhs_adv(row_index);
				counter_row_simplified++;
			}
		
		}

		Eigen::MatrixXd adv_mat_proj = c_f_all_PODmodes_wo_dbc.transpose() * adv_matrix_simplified * c_f_all_PODmodes_wo_dbc;
		Eigen::VectorXd adv_rhs_proj = c_f_all_PODmodes_wo_dbc.transpose() * adv_rhs_add;

		adv_mats_proj_x[trafo_iter] = adv_mat_proj;
		adv_vec_proj_x[trafo_iter] = adv_rhs_proj;

		CLNS.DoInitialiseAdv(PhysBase_zero , curr_PhysBaseVec_y ); // call with parameter in phys state

		adv_matrix = CLNS.Get_advection_matrix();
		add_to_rhs_adv = adv_matrix * f_bnd_dbc_full_size;   
		counter_row_simplified = 0;
		counter_col_simplified = 0;
		for (int row_index=0; row_index < f_bnd_dbc_full_size.rows(); ++row_index)
		{
			for (int col_index=0; col_index < f_bnd_dbc_full_size.rows(); ++col_index)
			{
				if ((!elem_loc_dbc.count(row_index)) && (!elem_loc_dbc.count(col_index)))
				{
					adv_matrix_simplified(counter_row_simplified, counter_col_simplified) = adv_matrix(row_index, col_index);
					counter_col_simplified++;
				}
			}
			counter_col_simplified = 0;
			if (!elem_loc_dbc.count(row_index))
			{
				adv_rhs_add(counter_row_simplified) = add_to_rhs_adv(row_index);
				counter_row_simplified++;
			}
		
		}
		adv_mat_proj = c_f_all_PODmodes_wo_dbc.transpose() * adv_matrix_simplified * c_f_all_PODmodes_wo_dbc;
		adv_rhs_proj = c_f_all_PODmodes_wo_dbc.transpose() * adv_rhs_add;
		adv_mats_proj_y[trafo_iter] = adv_mat_proj;
		adv_vec_proj_y[trafo_iter] = adv_rhs_proj;

	}

	Eigen::MatrixXd mat_compare = Eigen::MatrixXd::Zero(f_bnd_dbc_full_size.rows(), 3);
//	get it all together:
	int current_index = 2; // i.e., \nu = 1, serves as reference parameter
	double current_nu = param_vector[current_index];
	CLNS.Set_m_kinvis( current_nu );
	CLNS.DoInitialiseAdv(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);
	Eigen::MatrixXd test_complete_matrix = CLNS.Get_complete_matrix();
	Eigen::MatrixXd the_const_one = CLNS.Get_no_advection_matrix_pressure();
	Eigen::MatrixXd the_ABCD_one = CLNS.Get_no_advection_matrix_ABCD();
	Eigen::MatrixXd the_const_one_simplified = CLNS.remove_cols_and_rows(the_const_one, elem_loc_dbc);
	Eigen::MatrixXd the_ABCD_one_simplified = CLNS.remove_cols_and_rows(the_ABCD_one, elem_loc_dbc);
	Eigen::MatrixXd the_const_one_proj = c_f_all_PODmodes_wo_dbc.transpose() * the_const_one_simplified * c_f_all_PODmodes_wo_dbc;
	Eigen::MatrixXd the_ABCD_one_proj = c_f_all_PODmodes_wo_dbc.transpose() * the_ABCD_one_simplified * c_f_all_PODmodes_wo_dbc;
	Eigen::VectorXd the_ABCD_one_rhs = the_ABCD_one * f_bnd_dbc_full_size;
	Eigen::VectorXd the_const_one_rhs = the_const_one * f_bnd_dbc_full_size;
	Eigen::VectorXd the_ABCD_one_rhs_simplified = CLNS.remove_rows(the_ABCD_one_rhs, elem_loc_dbc);
	Eigen::VectorXd the_const_one_rhs_simplified = CLNS.remove_rows(the_const_one_rhs, elem_loc_dbc);
	Eigen::VectorXd the_ABCD_one_rhs_proj = c_f_all_PODmodes_wo_dbc.transpose() * the_ABCD_one_rhs_simplified;
	Eigen::VectorXd the_const_one_rhs_proj = c_f_all_PODmodes_wo_dbc.transpose() * the_const_one_rhs_simplified;


	// start sweeping 
	for (int iter_index = 0; iter_index < Nmax; ++iter_index)
	{
		current_index = iter_index;
		current_nu = param_vector[current_index];

		Eigen::VectorXd c_snapshot_x(CLNS.GetNpoints());
		Eigen::VectorXd c_snapshot_y(CLNS.GetNpoints());
		for (int i = 0; i < CLNS.GetNpoints(); ++i)
		{
			c_snapshot_x(i) = snapshot_x_collection[current_index][i];
			c_snapshot_y(i) = snapshot_y_collection[current_index][i];
		}
		Eigen::VectorXd curr_x_proj = eigen_phys_basis_x.transpose() * c_snapshot_x;
		Eigen::VectorXd curr_y_proj = eigen_phys_basis_y.transpose() * c_snapshot_y;
		Eigen::MatrixXd recovered_affine_adv_mat_proj_xy = Eigen::MatrixXd::Zero(RBsize, RBsize);
		for (int i = 0; i < RBsize; ++i)
		{
			recovered_affine_adv_mat_proj_xy += adv_mats_proj_x[i] * curr_x_proj[i] + adv_mats_proj_y[i] * curr_y_proj[i];
		}
//		Eigen::MatrixXd affine_mat_proj = recovered_affine_adv_mat_proj_xy;
		Eigen::MatrixXd affine_mat_proj = the_const_one_proj + current_nu * the_ABCD_one_proj + recovered_affine_adv_mat_proj_xy;
//		cout << "sweep affine mat_proj " << affine_mat_proj << endl;

		Eigen::VectorXd recovered_affine_adv_rhs_proj_xy = Eigen::VectorXd::Zero(RBsize); 
		for (int i = 0; i < RBsize; ++i)
		{
			recovered_affine_adv_rhs_proj_xy += adv_vec_proj_x[i] * curr_x_proj[i] + adv_vec_proj_y[i] * curr_y_proj[i];
		}		
//		cout << "sweep affine rhs_proj " << the_const_one_rhs_proj + current_nu * the_ABCD_one_rhs_proj + recovered_affine_adv_rhs_proj_xy << endl;
		Eigen::VectorXd solve_affine = affine_mat_proj.colPivHouseholderQr().solve(the_const_one_rhs_proj + current_nu * the_ABCD_one_rhs_proj + recovered_affine_adv_rhs_proj_xy);
		cout << "solve_affine " << solve_affine << endl;
		Eigen::VectorXd repro_solve_affine = c_f_all_PODmodes_wo_dbc * solve_affine;
		Eigen::VectorXd reconstruct_solution = Eigen::VectorXd::Zero(f_bnd_dbc_full_size.rows());
		int counter_wo_dbc = 0;
		for (int row_index=0; row_index < f_bnd_dbc_full_size.rows(); ++row_index)
		{
			if (!elem_loc_dbc.count(row_index))
			{
				reconstruct_solution(row_index) = repro_solve_affine(counter_wo_dbc);
				counter_wo_dbc++;
			}
			else
			{
				reconstruct_solution(row_index) = -f_bnd_dbc_full_size(row_index);
			}
		}
		mat_compare.col(0) = collect_f_all.col(current_index);
		mat_compare.col(1) = reconstruct_solution; // sembra abbastanza bene
		mat_compare.col(2) = mat_compare.col(1) + mat_compare.col(0);
//		cout << mat_compare << endl;
		cout << "relative error norm: " << mat_compare.col(2).norm() / mat_compare.col(0).norm() << endl;

	}

        session->Finalise();
    }
    catch (const std::runtime_error&)
    {
        return 1;
    }
    catch (const std::string& eStr)
    {
        cout << "Error: " << eStr << endl;
    }
    
    return 0;
}
