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

	CLNS.collect_f_all = collect_f_all;
	CLNS.PODmodes = collect_f_all_PODmodes;
	CLNS.set_MtM();
	Eigen::MatrixXd MtM = CLNS.MtM;
	Eigen::VectorXd f_bnd_dbc_full_size = CLNS.f_bnd_dbc_full_size;
	// c_f_all_PODmodes_wo_dbc becomes CLNS.RB
	Eigen::MatrixXd c_f_all_PODmodes_wo_dbc = CLNS.RB;

	Array<OneD, double> PhysBase_zero(CLNS.GetNpoints(), 0.0);

	CLNS.gen_phys_base_vecs();
	Array<OneD, Array<OneD, double> > PhysBaseVec_x = CLNS.PhysBaseVec_x;
	Array<OneD, Array<OneD, double> > PhysBaseVec_y = CLNS.PhysBaseVec_y;
	Eigen::MatrixXd eigen_phys_basis_x = CLNS.eigen_phys_basis_x;
	Eigen::MatrixXd eigen_phys_basis_y = CLNS.eigen_phys_basis_y;
	Array<OneD, Array<OneD, double> > orth_PhysBaseVec_x = CLNS.orth_PhysBaseVec_x;
	Array<OneD, Array<OneD, double> > orth_PhysBaseVec_y = CLNS.orth_PhysBaseVec_y; 

	CLNS.gen_proj_adv_terms();
	Array<OneD, Eigen::MatrixXd> adv_mats_proj_x = CLNS.adv_mats_proj_x;
	Array<OneD, Eigen::MatrixXd> adv_mats_proj_y = CLNS.adv_mats_proj_y;
	Array<OneD, Eigen::VectorXd> adv_vec_proj_x = CLNS.adv_vec_proj_x;
	Array<OneD, Eigen::VectorXd> adv_vec_proj_y = CLNS.adv_vec_proj_y;

	Eigen::MatrixXd mat_compare = Eigen::MatrixXd::Zero(f_bnd_dbc_full_size.rows(), 3);
	CLNS.gen_reference_matrices();

	// start sweeping 
	for (int iter_index = 0; iter_index < Nmax; ++iter_index)
	{
		int current_index = iter_index;
		double current_nu = param_vector[current_index];
		Eigen::MatrixXd curr_xy_proj = CLNS.project_onto_basis(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);

		Eigen::MatrixXd affine_mat_proj = CLNS.gen_affine_mat_proj(current_nu);
		Eigen::VectorXd affine_vec_proj = CLNS.gen_affine_vec_proj(current_nu);

		Eigen::VectorXd solve_affine = affine_mat_proj.colPivHouseholderQr().solve(affine_vec_proj);
		cout << "solve_affine " << solve_affine << endl;
		Eigen::VectorXd repro_solve_affine = c_f_all_PODmodes_wo_dbc * solve_affine;
		Eigen::VectorXd reconstruct_solution = CLNS.reconstruct_solution_w_dbc(repro_solve_affine);
		mat_compare.col(0) = collect_f_all.col(current_index);
		mat_compare.col(1) = reconstruct_solution; // sembra abbastanza bene
		mat_compare.col(2) = mat_compare.col(1) - mat_compare.col(0);
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
