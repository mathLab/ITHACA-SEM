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
        drv = GetDriverFactory().CreateInstance(vDriverModule, session);  // why is that necessary?

	int snapshots_to_be_collected_aka_Nmax = 4;  
	int Nmax = snapshots_to_be_collected_aka_Nmax;
	Array<OneD, NekDouble> param_vector(snapshots_to_be_collected_aka_Nmax);
	// = [0.1 0.5 1 10];
	
	param_vector[0] = 0.1;
	param_vector[1] = 0.5;
	param_vector[2] = 1;
	param_vector[3] = 10;
	
	for (int i = 0; i<Nmax; i++)
	{
//		cout << "param_vector[i] " << param_vector[i] << endl;   
	}

	CoupledLinearNS_TT babyCLNS(session);
	babyCLNS.InitObject();


	// some IO checks
//	cout << "session->DefinesFunction(InitialConditions) " << session->DefinesFunction("InitialConditions") << endl;
//	cout << "session->DefinesFunction(TestSnap1) " << session->DefinesFunction("TestSnap1") << endl;
//	cout << "session->DefinesFunction(TestSnap2) " << session->DefinesFunction("TestSnap2") << endl;
//	cout << "session->DefinesFunction(TestSnap3) " << session->DefinesFunction("TestSnap3") << endl;
//	cout << "session->DefinesFunction(TestSnap4) " << session->DefinesFunction("TestSnap4") << endl;
//	cout << "session->DefinesFunction(TestSnap5) " << session->DefinesFunction("TestSnap5") << endl;
//	cout << "session->DefinesFunction(\"AdvectionVelocity\") " << session->DefinesFunction("AdvectionVelocity") << endl;

	int nvelo = 2;
        Array<OneD, Array<OneD, NekDouble> > test_load_snapshot1(nvelo); // for a 2D problem
        Array<OneD, Array<OneD, NekDouble> > test_load_snapshot2(nvelo); // for a 2D problem
        Array<OneD, Array<OneD, NekDouble> > test_load_snapshot3(nvelo); // for a 2D problem
        Array<OneD, Array<OneD, NekDouble> > test_load_snapshot4(nvelo); // for a 2D problem
	Array<OneD, Array<OneD, NekDouble> > snapshot_collection(snapshots_to_be_collected_aka_Nmax);
	Array<OneD, Array<OneD, NekDouble> > snapshot_x_collection(snapshots_to_be_collected_aka_Nmax);
	Array<OneD, Array<OneD, NekDouble> > snapshot_y_collection(snapshots_to_be_collected_aka_Nmax);
	// probably also need this as EigenMatrix

        for(int i = 0; i < nvelo; ++i)
        {
            test_load_snapshot1[i] = Array<OneD, NekDouble> (babyCLNS.GetNpoints(), 0.0);  // number of phys points
            test_load_snapshot2[i] = Array<OneD, NekDouble> (babyCLNS.GetNpoints(), 0.0);  // number of phys points
            test_load_snapshot3[i] = Array<OneD, NekDouble> (babyCLNS.GetNpoints(), 0.0);  // number of phys points
            test_load_snapshot4[i] = Array<OneD, NekDouble> (babyCLNS.GetNpoints(), 0.0);  // number of phys points
        }
               
        std::vector<std::string> fieldStr;
        for(int i = 0; i < nvelo; ++i)
        {
           fieldStr.push_back(session->GetVariable(i));
//           cout << "session->GetVariable(i) " << session->GetVariable(i) << endl;
        }
        babyCLNS.EvaluateFunction(fieldStr, test_load_snapshot1, "TestSnap1");
	snapshot_x_collection[0] = test_load_snapshot1[0];
	snapshot_y_collection[0] = test_load_snapshot1[1];
        babyCLNS.EvaluateFunction(fieldStr, test_load_snapshot2, "TestSnap2");
	snapshot_x_collection[1] = test_load_snapshot2[0];
	snapshot_y_collection[1] = test_load_snapshot2[1];
        babyCLNS.EvaluateFunction(fieldStr, test_load_snapshot3, "TestSnap3");
	snapshot_x_collection[2] = test_load_snapshot3[0];
	snapshot_y_collection[2] = test_load_snapshot3[1];
        babyCLNS.EvaluateFunction(fieldStr, test_load_snapshot4, "TestSnap4");
	snapshot_x_collection[3] = test_load_snapshot4[0];
	snapshot_y_collection[3] = test_load_snapshot4[1];


	// trafo the snapshots and see if they have really converged
	// so I need to capture the bnd / p / int as well as the local / global / phys
	// Domanda: how to change the k_invis for an already declared class -- Risposta: use Getter and Setter functions (to be implemented - done)



	session->SetSolverInfo("SolverType", "CoupledLinearisedNS_trafoP");
	CoupledLinearNS_trafoP babyCLNS_trafo(session);
	babyCLNS_trafo.InitObject();


	// Q: how to know the bnd / p / int sizes beforehand
	babyCLNS_trafo.DoInitialise(); 
	babyCLNS_trafo.DoSolve();
	Eigen::MatrixXd collect_f_bnd( babyCLNS_trafo.curr_f_bnd.size() , Nmax );
	Eigen::MatrixXd collect_f_p( babyCLNS_trafo.curr_f_p.size() , Nmax );
	Eigen::MatrixXd collect_f_int( babyCLNS_trafo.curr_f_int.size() , Nmax );

	for (int i=0; i<Nmax; i++)
	{
		babyCLNS_trafo.Set_m_kinvis( param_vector[i] );	
//		cout << "babyCLNS_trafo.Get_m_kinvis " << babyCLNS_trafo.Get_m_kinvis() << endl;
	//	babyCLNS_trafo.DoInitialise();
		babyCLNS_trafo.DoInitialiseAdv(snapshot_x_collection[i], snapshot_y_collection[i]); // replaces .DoInitialise();
		babyCLNS_trafo.DoSolve();

		// compare the accuracy
		Array<OneD, MultiRegions::ExpListSharedPtr> m_fields_t = babyCLNS_trafo.UpdateFields();
		m_fields_t[0]->BwdTrans(m_fields_t[0]->GetCoeffs(), m_fields_t[0]->UpdatePhys());
		m_fields_t[1]->BwdTrans(m_fields_t[1]->GetCoeffs(), m_fields_t[1]->UpdatePhys());
		Array<OneD, NekDouble> out_field_trafo_x(babyCLNS_trafo.GetNpoints(), 0.0);
		Array<OneD, NekDouble> out_field_trafo_y(babyCLNS_trafo.GetNpoints(), 0.0);

		Eigen::VectorXd csx0_trafo(babyCLNS_trafo.GetNpoints());
		Eigen::VectorXd csy0_trafo(babyCLNS_trafo.GetNpoints());
		Eigen::VectorXd csx0(babyCLNS_trafo.GetNpoints());
		Eigen::VectorXd csy0(babyCLNS_trafo.GetNpoints());


		babyCLNS_trafo.CopyFromPhysField(0, out_field_trafo_x); 
		babyCLNS_trafo.CopyFromPhysField(1, out_field_trafo_y);
		for( int index_conv = 0; index_conv < babyCLNS_trafo.GetNpoints(); ++index_conv)
		{
			csx0_trafo(index_conv) = out_field_trafo_x[index_conv];
			csy0_trafo(index_conv) = out_field_trafo_y[index_conv];
			csx0(index_conv) = snapshot_x_collection[i][index_conv];
			csy0(index_conv) = snapshot_y_collection[i][index_conv];
		}

//		cout << "csx0.norm() " << csx0.norm() << endl;
//		cout << "csx0_trafo.norm() " << csx0_trafo.norm() << endl;
//		cout << "csy0.norm() " << csy0.norm() << endl;
//		cout << "csy0_trafo.norm() " << csy0_trafo.norm() << endl;
		

		Eigen::VectorXd trafo_f_bnd = babyCLNS_trafo.curr_f_bnd;
		Eigen::VectorXd trafo_f_p = babyCLNS_trafo.curr_f_p;
		Eigen::VectorXd trafo_f_int = babyCLNS_trafo.curr_f_int;

		collect_f_bnd.col(i) = trafo_f_bnd;
		collect_f_p.col(i) = trafo_f_p;
		collect_f_int.col(i) = trafo_f_int;

	}

	Eigen::MatrixXd collect_f_all( babyCLNS_trafo.curr_f_bnd.size()+babyCLNS_trafo.curr_f_p.size()+babyCLNS_trafo.curr_f_int.size() , Nmax );
	collect_f_all.block(0,0,collect_f_bnd.rows(),collect_f_bnd.cols()) = collect_f_bnd;
	collect_f_all.block(collect_f_bnd.rows(),0,collect_f_p.rows(),collect_f_p.cols()) = collect_f_p;
	collect_f_all.block(collect_f_bnd.rows()+collect_f_p.rows(),0,collect_f_int.rows(),collect_f_int.cols()) = collect_f_int;

	Eigen::BDCSVD<Eigen::MatrixXd> svd_collect_f_all(collect_f_all, Eigen::ComputeThinU);
	cout << "svd_collect_f_all.singularValues() " << svd_collect_f_all.singularValues() << endl << endl;
	Eigen::MatrixXd collect_f_all_PODmodes = svd_collect_f_all.matrixU();

	// here probably limit to something like 99.99 percent of PODenergy, this will set RBsize
	int RBsize = Nmax; // the case of using all




/*	if (session->DefinesParameter("Kinvis"))
	{
		double testp = session->GetParameter("Kinvis");
		cout << "testp " << testp << endl;
	}
*/

	// now can get the phys as well as bnd / p / int solution :)
/*	Array<OneD, MultiRegions::ExpListSharedPtr> m_fields_t = myCLNS_trafo.UpdateFields();

	// map the loc field to the phys field
	m_fields_t[0]->BwdTrans(m_fields_t[0]->GetCoeffs(), m_fields_t[0]->UpdatePhys()); // needed if going with the equation system 


	Array<OneD, NekDouble> out_field_trafo_x(m_equ[0]->GetNpoints(), 0.0);
	Array<OneD, NekDouble> out_field_trafo_x2(m_equ[0]->GetNpoints(), 0.0);
	Array<OneD, NekDouble> out_field_trafo_y(m_equ[0]->GetNpoints(), 0.0);
	myCLNS_trafo.CopyFromPhysField(0, out_field_trafo_x); 
	myCLNS_trafo.CopyFromPhysField(1, out_field_trafo_y);
	out_field_trafo_x2 = m_fields_t[0]->UpdatePhys();
	// compare input and output
	Eigen::VectorXd csx0(collected_snapshots_x[0].num_elements());
	Eigen::VectorXd csx0_trafo(collected_snapshots_x[0].num_elements());
	Eigen::VectorXd csx0_trafo2(collected_snapshots_x[0].num_elements());
	Eigen::VectorXd csy0(collected_snapshots_x[0].num_elements());
	Eigen::VectorXd csy0_trafo(collected_snapshots_x[0].num_elements());
	for( int i = 0; i < collected_snapshots_x[0].num_elements(); ++i)
	{
		csx0(i) = collected_snapshots_x[0][i];
		csx0_trafo(i) = out_field_trafo_x[i];
		csx0_trafo2(i) = out_field_trafo_x2[i];
		csy0(i) = collected_snapshots_y[0][i];
		csy0_trafo(i) = out_field_trafo_y[i];
	}
	cout << "collected_snapshots_x[0] " << csx0.norm() << endl;
	cout << "trafo collected_snapshots_x[0] " << csx0_trafo.norm() << endl;
	cout << "trafo2 collected_snapshots_x[0] " << csx0_trafo2.norm() << endl;

//	Eigen::VectorXd trafo_f_bnd = myCLNS_trafo.curr_f_bnd;
//	Eigen::VectorXd trafo_f_p = myCLNS_trafo.curr_f_p;
//	Eigen::VectorXd trafo_f_int = myCLNS_trafo.curr_f_int;

	// do a full order looping
	int no_of_loops = -1;
	for ( int iter = 0; iter < no_of_loops; ++iter)
	{
//		myCLNS_trafo.InitObject();
		myCLNS_trafo.DoInitialiseAdv(out_field_trafo_x, out_field_trafo_y); // replaces .DoInitialise();
//		myCLNS_trafo.DoInitialise();
		myCLNS_trafo.DoSolve();
		myCLNS_trafo.Output();

		m_fields_t = myCLNS_trafo.UpdateFields();
		m_fields_t[0]->BwdTrans(m_fields_t[0]->GetCoeffs(), m_fields_t[0]->UpdatePhys()); 
		myCLNS_trafo.CopyFromPhysField(0, out_field_trafo_x); 
		myCLNS_trafo.CopyFromPhysField(1, out_field_trafo_y);
		for( int i = 0; i < collected_snapshots_x[0].num_elements(); ++i)
		{
			csx0_trafo(i) = out_field_trafo_x[i];
			csy0_trafo(i) = out_field_trafo_y[i];
		}

		cout << "trafo collected_snapshots_x[0] " << csx0_trafo.norm() << endl;

	}

	// do the bnd / p / int transform
	int no_snaps = snapshots_to_be_collected_aka_Nmax;
	Eigen::MatrixXd c_f_bnd( myCLNS_trafo.curr_f_bnd.size() , no_snaps );
	Eigen::MatrixXd c_f_p( myCLNS_trafo.curr_f_p.size() , no_snaps );
	Eigen::MatrixXd c_f_int( myCLNS_trafo.curr_f_int.size() , no_snaps );
	for(int trafo_iter = 0; trafo_iter < no_snaps; trafo_iter++)
	{
		myCLNS_trafo.InitObject();
		myCLNS_trafo.DoInitialiseAdv(collected_snapshots_x[trafo_iter], collected_snapshots_y[trafo_iter]); // replaces .DoInitialise();
		myCLNS_trafo.DoSolve();
		myCLNS_trafo.Output();
	
		c_f_bnd.col(trafo_iter) = myCLNS_trafo.curr_f_bnd;
		c_f_p.col(trafo_iter) = myCLNS_trafo.curr_f_p;
		c_f_int.col(trafo_iter) = myCLNS_trafo.curr_f_int;
		

	}
*/


	// identifiy the D-BC dofs --- can do this by checking the c_f_bnd dofs
	int no_dbc_in_loc = 0;
	int no_not_dbc_in_loc = 0;
	std::set<int> elem_loc_dbc;
	set<int> elem_not_loc_dbc;
	for ( int index_c_f_bnd = 0; index_c_f_bnd < collect_f_bnd.rows(); index_c_f_bnd++ )
	{
		if (collect_f_bnd(index_c_f_bnd,0) == collect_f_bnd(index_c_f_bnd,1))
		{
			no_dbc_in_loc++;
			elem_loc_dbc.insert(index_c_f_bnd);
		}
		else
		{
			no_not_dbc_in_loc++;
			elem_not_loc_dbc.insert(index_c_f_bnd);
		}
	}

//	cout << "no_dbc_in_loc " << no_dbc_in_loc << endl;
//	cout << "elem_loc_dbc.size() " << elem_loc_dbc.size() << endl;

//	cout << "no_not_dbc_in_loc " << no_not_dbc_in_loc << endl;
//	cout << "elem_not_loc_dbc.size() " << elem_not_loc_dbc.size() << endl;



	// deal with the Dirichlet BC
	// can this be simplified when not doing the bnd / p / int separation ?
	// can be seen a 2-step procedure, non_adv and adv separately



	// do a all = bnd/p/int POD
/*

	Eigen::MatrixXd c_f_all( myCLNS_trafo.curr_f_bnd.size()+myCLNS_trafo.curr_f_p.size()+myCLNS_trafo.curr_f_int.size() , no_snaps );
	c_f_all.block(0,0,c_f_bnd.rows(),c_f_bnd.cols()) = c_f_bnd;
	c_f_all.block(c_f_bnd.rows(),0,c_f_p.rows(),c_f_p.cols()) = c_f_p;
	c_f_all.block(c_f_bnd.rows()+c_f_p.rows(),0,c_f_int.rows(),c_f_int.cols()) = c_f_int;

	Eigen::BDCSVD<Eigen::MatrixXd> svd_c_f_all(c_f_all, Eigen::ComputeThinU);
	cout << svd_c_f_all.singularValues() << endl << endl;
	Eigen::MatrixXd c_f_all_PODmodes = svd_c_f_all.matrixU();
	// limit thyself to the 99.99%

*/

	// do the MM multiplication
	// need now a proper _TT object
	Array<OneD, CoupledLocalToGlobalC0ContMapSharedPtr> locToGloMap = babyCLNS_trafo.m_locToGloMap;
        const Array<OneD,const int>& loctoglobndmap = locToGloMap[0]->GetLocalToGlobalBndMap();
        const Array<OneD,const NekDouble>& loctoglobndsign = locToGloMap[0]->GetLocalToGlobalBndSign();


	session->SetSolverInfo("SolverType", "CoupledLinearisedNS_TT");
//	CoupledLinearNS_TT babyCLNS(session);  has happened already
//	babyCLNS.InitObject();   has happened already
	babyCLNS.DoInitialise();
	babyCLNS.DoSolve(); // a solve to get the dims right


	int nBndDofs = locToGloMap[0]->GetNumGlobalBndCoeffs();  // number of global bnd dofs
	Eigen::MatrixXd Mtrafo(babyCLNS.RB_A.rows(), nBndDofs);
	Array<OneD, MultiRegions::ExpListSharedPtr> m_fields = babyCLNS.UpdateFields();
        int  nel  = m_fields[0]->GetNumElmts(); // number of spectral elements
//	cout << "nel " << nel << endl;
//	cout << "Mtrafo.rows() " << Mtrafo.rows() << endl;
//	cout << "Mtrafo.cols() " << Mtrafo.cols() << endl;
	int nsize_bndry_p1 = loctoglobndmap.num_elements() / nel;
	int nsize_bndry = nsize_bndry_p1-1;
//	cout << "loctoglobndmap.num_elements() / nel -1 " << loctoglobndmap.num_elements() / nel -1 << endl;

	for (int curr_elem = 0; curr_elem < nel; curr_elem++)
	{
		int cnt = curr_elem*nsize_bndry_p1;
		int cnt_no_pp = curr_elem*nsize_bndry;
		for ( int index_ele = 0; index_ele < nsize_bndry_p1; index_ele++ )
		{
			int gid1 = loctoglobndmap[cnt+index_ele];
			int sign1 = loctoglobndsign[cnt+index_ele];
			if ((gid1 >= 0) && (index_ele < nsize_bndry))
			{
//				cout << "cnt_no_pp + index_ele " << cnt_no_pp + index_ele << endl;
//				cout << "gid1 " << gid1 << endl;
				Mtrafo(cnt_no_pp + index_ele, gid1) = sign1;
			}
		}
	}

	// do the MtM multiplication or try with an all-M approach
	Eigen::MatrixXd MtM(babyCLNS.RB_A.rows(), babyCLNS.RB_A.rows());
	MtM = Mtrafo * Mtrafo.transpose();
	babyCLNS.MtM = MtM;
	// assume have all the local sub-matrices in myCLNS
	// to break the bnd/p/int, could just form two matrices, "no_adv" and "adv", whereby "adv" has the trilinear projection
	// better initialise to zero
	Eigen::MatrixXd no_adv_matrix = Eigen::MatrixXd::Zero(babyCLNS.RB_A.rows() + babyCLNS.RB_Dbnd.rows() + babyCLNS.RB_C.cols(), babyCLNS.RB_A.cols() + babyCLNS.RB_Dbnd.rows() + babyCLNS.RB_B.cols() );
//	cout << "no_adv_matrix.rows() " << no_adv_matrix.rows() << endl;
//	cout << "no_adv_matrix.cols() " << no_adv_matrix.cols() << endl;
	// write Eigen::MatrixXd no_adv_matrix blockwise 
	no_adv_matrix.block(0, 0, babyCLNS.RB_A.rows(), babyCLNS.RB_A.cols()) = MtM * babyCLNS.RB_A_no_adv;
	no_adv_matrix.block(0, babyCLNS.RB_A.cols(), babyCLNS.RB_Dbnd.cols(), babyCLNS.RB_Dbnd.rows()) = -MtM * babyCLNS.RB_Dbnd.transpose();
	no_adv_matrix.block(0, babyCLNS.RB_A.cols() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_B.rows(), babyCLNS.RB_B.cols()) = MtM * babyCLNS.RB_B_no_adv;
	no_adv_matrix.block(babyCLNS.RB_A.rows(), 0, babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_Dbnd.cols()) = -babyCLNS.RB_Dbnd;
	no_adv_matrix.block(babyCLNS.RB_A.rows(), babyCLNS.RB_A.cols() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_Dint.rows(), babyCLNS.RB_Dint.cols()) = -babyCLNS.RB_Dint;	
	no_adv_matrix.block(babyCLNS.RB_A.rows() + babyCLNS.RB_Dbnd.rows(), 0, babyCLNS.RB_C.cols(), babyCLNS.RB_C.rows()) = babyCLNS.RB_C_no_adv.transpose();
	no_adv_matrix.block(babyCLNS.RB_A.rows() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_A.cols(), babyCLNS.RB_Dint.cols(), babyCLNS.RB_Dint.rows()) = -babyCLNS.RB_Dint.transpose();
	no_adv_matrix.block(babyCLNS.RB_A.rows() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_A.cols() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_D.rows(), babyCLNS.RB_D.cols()) = babyCLNS.RB_D_no_adv;


	// get the D-BC to the other side
	// remove from the c_f_all block the D_BC and put seperately
/*	cout << "the elements of elem_loc_dbc: ";
	for (std::set<int>::iterator it=elem_loc_dbc.begin(); it!=elem_loc_dbc.end(); ++it)
	{
		cout << " " << *it << endl;
		cout << "c_f_all @ iter " << c_f_all.row(*it) << endl;
		cout << "c_f_all_PODmodes @ iter " << c_f_all_PODmodes.row(*it) << endl;

	}
*/
//	cout << "all of c_f_all" << endl;
//	cout << c_f_all;
	
	Eigen::MatrixXd c_f_all_PODmodes_wo_dbc(collect_f_all_PODmodes.rows() - no_dbc_in_loc, collect_f_all_PODmodes.cols());
	Eigen::VectorXd f_bnd_dbc(no_dbc_in_loc);
	Eigen::VectorXd f_bnd_dbc_full_size(collect_f_all_PODmodes.rows());
	int counter_all = 0;
	int counter_dbc = 0;
	for (int index=0; index < collect_f_all_PODmodes.rows(); ++index)
	{
//		cout << index << endl;

		if (!elem_loc_dbc.count(index))
		{
			f_bnd_dbc_full_size(index) = 0;
			c_f_all_PODmodes_wo_dbc.row(counter_all) = collect_f_all_PODmodes.row(index);
			counter_all++;
//			cout << " counter all " << counter_all << endl;
		}
		else
		{
			f_bnd_dbc_full_size(index) = collect_f_all(index,0);
			f_bnd_dbc(counter_dbc) = collect_f_all(index,0);
			counter_dbc++;
//			cout << " counter dbc " << counter_dbc << endl;
//			cout << " no_dbc_in_loc " << no_dbc_in_loc << endl;
		}

	}

	// compute add_to rhs
	Eigen::VectorXd add_to_rhs(collect_f_all_PODmodes.rows()); // probably need this for adv and non-adv
	add_to_rhs = no_adv_matrix * f_bnd_dbc_full_size;   
	// also implement the removal of dbc dof
	// the set of to-be-removed cols and rows is elem_loc_dbc
	Eigen::MatrixXd no_adv_matrix_simplified = Eigen::MatrixXd::Zero(no_adv_matrix.rows() - no_dbc_in_loc, no_adv_matrix.cols() - no_dbc_in_loc);
	// a naive algorithm (should test the timing here on a "real-world" example)
	int counter_row_simplified = 0;
	int counter_col_simplified = 0;
	for (int row_index=0; row_index < no_adv_matrix.rows(); ++row_index)
	{
		for (int col_index=0; col_index < no_adv_matrix.cols(); ++col_index)
		{
			if ((!elem_loc_dbc.count(row_index)) && (!elem_loc_dbc.count(col_index)))
			{
				no_adv_matrix_simplified(counter_row_simplified, counter_col_simplified) = no_adv_matrix(row_index, col_index);
				counter_col_simplified++;
			}
			
		}
		counter_col_simplified = 0;
		if (!elem_loc_dbc.count(row_index))
		{
			counter_row_simplified++;
		}
		
	}

//	singular:       cout << "no_adv_matrix.eigenvalues() " << no_adv_matrix.eigenvalues() << endl;
//	singular:	cout << "no_adv_matrix_simplified.eigenvalues() " << no_adv_matrix_simplified.eigenvalues() << endl;
	

	Array<OneD, double> PhysBase_zero(babyCLNS.GetNpoints(), 0.0);

	Array<OneD, Array<OneD, double> > PhysBaseVec_x(RBsize); // or do I want it as an Eigen-matrix
	Array<OneD, Array<OneD, double> > PhysBaseVec_y(RBsize);

	// deal with the advection projection
	// first build the full advection matrix and then multiply with D-BC, add to right-hand-side and then remove the dbc dof
	// also need to have nicely orthonormalized (iteratively) basis vectors...
	// need to have every RB basis vector in phys for that, i.e., after the svd misch-masch
	// probably can use a separate object for that
	
	// first a mapping bnd / int --> loc  // this should just be a super-simple re-mapping of dofs... find a matrix to do that or at least put it in a function
        // Unpack solution from Bnd and F_int to v_coeffs 


	
	for (int curr_trafo_iter=0; curr_trafo_iter < RBsize; curr_trafo_iter++)
	{
		Eigen::VectorXd f_bnd = collect_f_all_PODmodes.block(0, curr_trafo_iter, babyCLNS_trafo.curr_f_bnd.size(), 1);
		Eigen::VectorXd f_int = collect_f_all_PODmodes.block(babyCLNS_trafo.curr_f_bnd.size()+babyCLNS_trafo.curr_f_p.size(), curr_trafo_iter, babyCLNS_trafo.curr_f_int.size(), 1);
//		cout << "f_int " << f_int << endl;
		Array<OneD, MultiRegions::ExpListSharedPtr> fields = babyCLNS.UpdateFields(); 
	        Array<OneD, unsigned int> bmap, imap; 
		Array<OneD, double> field_0(babyCLNS.GetNcoeffs());
		Array<OneD, double> field_1(babyCLNS.GetNcoeffs());
		Array<OneD, double> curr_PhysBaseVec_x(babyCLNS.GetNpoints(), 0.0);
		Array<OneD, double> curr_PhysBaseVec_y(babyCLNS.GetNpoints(), 0.0);
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
//		fields[0]->BwdTrans_IterPerExp(field_0, curr_PhysBaseVec_x);
//		fields[0]->BwdTrans_IterPerExp(field_1, curr_PhysBaseVec_y);
		Array<OneD, double> test_nn = fields[0]->GetCoeffs();
//		cout << test_nn[323] << endl;
		fields[0]->BwdTrans_IterPerExp(fields[0]->GetCoeffs(), curr_PhysBaseVec_x);
//		cout << curr_PhysBaseVec_x[babyCLNS.GetNpoints()-1] << endl;
		fields[1]->BwdTrans_IterPerExp(fields[1]->GetCoeffs(), curr_PhysBaseVec_y);
		PhysBaseVec_x[curr_trafo_iter] = curr_PhysBaseVec_x;
		PhysBaseVec_y[curr_trafo_iter] = curr_PhysBaseVec_y;
		
	}

	// orthonormalize the phys basis
	// iterative modified Gram-Schmidt would be nice

//	babyCLNS.Output();
//	babyCLNS_trafo.Output();

//	babyCLNS.InitObject();
//	babyCLNS.DoInitialise();

	// do this better with eigen-like data-structures

	Eigen::MatrixXd eigen_phys_basis_x(babyCLNS.GetNpoints(), RBsize);
	Eigen::MatrixXd eigen_phys_basis_y(babyCLNS.GetNpoints(), RBsize);
	for (int index_phys_base=0; index_phys_base<babyCLNS.GetNpoints(); index_phys_base++)
	{
		for (int index_RBsize=0; index_RBsize<RBsize; index_RBsize++)
		{
			eigen_phys_basis_x(index_phys_base,index_RBsize) = PhysBaseVec_x[index_RBsize][index_phys_base];
			eigen_phys_basis_y(index_phys_base,index_RBsize) = PhysBaseVec_y[index_RBsize][index_phys_base];
		}
	}

	// the actual orthonormalization
//	Eigen::Vector3d test1;
//	test1 << 1, 2, 3;
//	cout << "test1.norm() " << test1.norm() << endl; // seems to be the 2-norm

	// do the orthonormalization in-place on eigen_phys_basis
	// beware of aliasing!

//	cout << eigen_phys_basis_x << endl;
//	cout << "eigen_phys_basis_x.rows() " << eigen_phys_basis_x.rows() << endl;
//	cout << "eigen_phys_basis_x.cols() " << eigen_phys_basis_x.cols() << endl;
//	cout << "babyCLNS.GetNpoints() " << babyCLNS.GetNpoints() << endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////// causes the segfault at the end /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
//		cout << "leftmostCols.rows() " << leftmostCols.rows() << endl;
//		cout << "leftmostCols.cols() " << leftmostCols.cols() << endl;
//		cout << "curr_col.norm() " << curr_col.norm() << endl;
//		cout << "leftmostCols.norm() " << leftmostCols.norm() << endl;
		orthogonal_complement = curr_col - leftmostCols * leftmostCols.transpose() * curr_col;
		// now the orthogonal complement should be orthogonal to the leftmostCols
//		cout << "orthogonal_complement.transpose() * leftmostCols " << (orthogonal_complement.transpose() * leftmostCols) << endl;
//		cout << "orthogonal_complement.norm() " << (orthogonal_complement.norm()) << endl;
		norm_curr_col = orthogonal_complement.norm();
		eigen_phys_basis_x.col(orth_iter) = orthogonal_complement / norm_curr_col;
//		eigen_phys_basis_x.col(orth_iter) = orthogonal_complement;
//		cout << "eigen_phys_basis_x.col(orth_iter).transpose() * leftmostCols " << (eigen_phys_basis_x.col(orth_iter).transpose() * leftmostCols) << endl;
		
		
//		cout << "eigen_phys_basis_x.transpose() * eigen_phys_basis_x " << eigen_phys_basis_x.transpose() * eigen_phys_basis_x << endl;


		curr_col = eigen_phys_basis_y.col(orth_iter);
		leftmostCols = eigen_phys_basis_y.leftCols(orth_iter);
		orthogonal_complement = curr_col - leftmostCols * leftmostCols.transpose() * curr_col;
		norm_curr_col = orthogonal_complement.norm();
		eigen_phys_basis_y.col(orth_iter) = orthogonal_complement / norm_curr_col;
	}


	// move it back to Nektar array

	Array<OneD, Array<OneD, double> > orth_PhysBaseVec_x(RBsize); // also have it as an Eigen-matrix: eigen_phys_basis_ x/y
	Array<OneD, Array<OneD, double> > orth_PhysBaseVec_y(RBsize);
	for (int index_RBsize=0; index_RBsize<RBsize; index_RBsize++)
	{
		Array<OneD, double> curr_iter_x(babyCLNS.GetNpoints());
		Array<OneD, double> curr_iter_y(babyCLNS.GetNpoints());

		for (int index_phys_base=0; index_phys_base<babyCLNS.GetNpoints(); index_phys_base++)	
		{
			curr_iter_x[index_phys_base] = eigen_phys_basis_x(index_phys_base,index_RBsize);
			curr_iter_y[index_phys_base] = eigen_phys_basis_y(index_phys_base,index_RBsize);			
		}
		orth_PhysBaseVec_x[index_RBsize] = curr_iter_x;
		orth_PhysBaseVec_y[index_RBsize] = curr_iter_y;			

	}

//	cout << "orth_PhysBaseVec_x[0][123] " << orth_PhysBaseVec_x[0][523] << endl;
//	cout << "orth_PhysBaseVec_x[1][123] " << orth_PhysBaseVec_x[1][523] << endl;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// 	babyCLNS.DoInitialise();

	// python code
/*	phys_basis_x[:,0] = phys_basis_x[:,0] / np.linalg.norm(phys_basis_x[:,0])
	phys_basis_y[:,0] = phys_basis_y[:,0] / np.linalg.norm(phys_basis_y[:,0])

	for orth_iter in range(1,RBsize):
		curr_phys_basis_x = phys_basis_x[:,orth_iter]
		curr_phys_basis_y = phys_basis_y[:,orth_iter]
		curr_phys_basis_x = curr_phys_basis_x - np.dot(phys_basis_x[:,0:orth_iter], np.dot( np.transpose(phys_basis_x[:,0:orth_iter]), curr_phys_basis_x))
		curr_phys_basis_y = curr_phys_basis_y - np.dot(phys_basis_y[:,0:orth_iter], np.dot( np.transpose(phys_basis_y[:,0:orth_iter]), curr_phys_basis_y))
		phys_basis_x[:,orth_iter] = curr_phys_basis_x / np.linalg.norm(curr_phys_basis_x)
		phys_basis_y[:,orth_iter] = curr_phys_basis_y / np.linalg.norm(curr_phys_basis_y)
*/





/*	int curr_trafo_iter = 0;
	Eigen::VectorXd f_bnd = collect_f_all_PODmodes.block(0, curr_trafo_iter, babyCLNS_trafo.curr_f_bnd.size(), 1);
	Eigen::VectorXd f_int = collect_f_all_PODmodes.block(babyCLNS_trafo.curr_f_bnd.size()+babyCLNS_trafo.curr_f_p.size(), curr_trafo_iter, babyCLNS_trafo.curr_f_int.size(), 1);
	// I think this messes up the content of myCLNS local field, b.c. of pointer
	// if it gets initialized again before re-used it should be fine
	Array<OneD, MultiRegions::ExpListSharedPtr> fields = babyCLNS.UpdateFields(); 
        Array<OneD, unsigned int> bmap, imap; 
	Array<OneD, double> field_0(babyCLNS.GetNcoeffs());
	Array<OneD, double> field_1(babyCLNS.GetNcoeffs());
        int cnt = 0;
	int cnt1 = 0;
	int nvel = 2;
//	cout  << "  " << myCLNS_trafo.m_singleMode << endl; would need to implement a "getter" function
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
*/

	// then map the loc field to the phys field m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(), m_fields[0]->UpdatePhys()); // needed if going with the equation system 
/*	Array<OneD, double> PhysBaseVec_x(myCLNS.GetNpoints(),0.0);
	Array<OneD, double> PhysBaseVec_y(myCLNS.GetNpoints(),0.0);
	Array<OneD, double> PhysBase_zero(myCLNS.GetNpoints(),0.0);
//	fields[0]->BwdTrans(fields[0]->GetCoeffs(), PhysBaseVec_x);
//	fields[1]->BwdTrans(fields[1]->GetCoeffs(), PhysBaseVec_y);		

	fields[0]->BwdTrans(field_0, PhysBaseVec_x);
	fields[0]->BwdTrans(field_1, PhysBaseVec_y);

	cout << "PhysBaseVec_x.num_elements() " << PhysBaseVec_x.num_elements() << endl;


*/


	// do create the projected adv matrices
//	int RBsize = 1; // initially have just one for test purposes
//	Eigen::MatrixXd c_f_bnd( myCLNS_trafo.curr_f_bnd.size() , no_snaps );
//	Eigen::MatrixXd c_f_p( myCLNS_trafo.curr_f_p.size() , no_snaps );
//	Eigen::MatrixXd c_f_int( myCLNS_trafo.curr_f_int.size() , no_snaps );


//        session->Finalise();
//	CoupledLinearNS_TT babyCLNS2(session2);
//	babyCLNS2.InitObject();
//	babyCLNS2.DoInitialise();
//	babyCLNS2.DoSolve();


	// need to collect all the projected adv matrices
	// have in total for lhs and rhs: RBsize*2 (because of x and y seperately)
	// also need to do the MtM mults, already there :)
	Array<OneD, Eigen::MatrixXd> adv_mats_proj_x(RBsize);
	Array<OneD, Eigen::MatrixXd> adv_mats_proj_y(RBsize);
	Array<OneD, Eigen::VectorXd> adv_vec_proj_x(RBsize);
	Array<OneD, Eigen::VectorXd> adv_vec_proj_y(RBsize);


	for(int trafo_iter = 0; trafo_iter < RBsize; trafo_iter++)
	{
		Array<OneD, double> curr_PhysBaseVec_x = orth_PhysBaseVec_x[trafo_iter];
		Array<OneD, double> curr_PhysBaseVec_y = orth_PhysBaseVec_y[trafo_iter];
		
		babyCLNS.InitObject();
//		babyCLNS.DoInitialise();
		babyCLNS.DoInitialiseAdv(curr_PhysBaseVec_x, PhysBase_zero); // call with parameter in phys state
//		babyCLNS2.DoSolve(); // do I need this here ?
//		babyCLNS2.Output();  // do I need this here ?

	//	c_f_bnd.col(trafo_iter) = myCLNS_trafo.curr_f_bnd;
	//	c_f_p.col(trafo_iter) = myCLNS_trafo.curr_f_p;
	//	c_f_int.col(trafo_iter) = myCLNS_trafo.curr_f_int;

		// collect the advection matrices
		// separate into rhs and lhs part
		// save in projected form for each RB basis function

		Eigen::MatrixXd adv_matrix = Eigen::MatrixXd::Zero(babyCLNS.RB_A.rows() + babyCLNS.RB_Dbnd.rows() + babyCLNS.RB_C.cols(), babyCLNS.RB_A.cols() + babyCLNS.RB_Dbnd.rows() + babyCLNS.RB_B.cols() );
//		cout << "adv_matrix.rows() " << adv_matrix.rows() << endl;
//		cout << "adv_matrix.cols() " << adv_matrix.cols() << endl;
		// write Eigen::MatrixXd no_adv_matrix blockwise 
/*		adv_matrix.block(0, 0, babyCLNS.RB_A.rows(), babyCLNS.RB_A.cols()) = MtM * babyCLNS.RB_A_adv;
		adv_matrix.block(0, babyCLNS.RB_A.cols(), babyCLNS.RB_Dbnd.cols(), babyCLNS.RB_Dbnd.rows()) = -MtM * babyCLNS.RB_Dbnd.transpose();
		adv_matrix.block(0, babyCLNS.RB_A.cols() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_B.rows(), babyCLNS.RB_B.cols()) = MtM * babyCLNS.RB_B_adv;
		adv_matrix.block(babyCLNS.RB_A.rows(), 0, babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_Dbnd.cols()) = -babyCLNS.RB_Dbnd;
		adv_matrix.block(babyCLNS.RB_A.rows(), babyCLNS.RB_A.cols() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_Dint.rows(), babyCLNS.RB_Dint.cols()) = -babyCLNS.RB_Dint;	
		adv_matrix.block(babyCLNS.RB_A.rows() + babyCLNS.RB_Dbnd.rows(), 0, babyCLNS.RB_C.cols(), babyCLNS.RB_C.rows()) = babyCLNS.RB_C_adv.transpose();
		adv_matrix.block(babyCLNS.RB_A.rows() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_A.cols(), babyCLNS.RB_Dint.cols(), babyCLNS.RB_Dint.rows()) = -babyCLNS.RB_Dint.transpose();
		adv_matrix.block(babyCLNS.RB_A.rows() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_A.cols() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_D.rows(), babyCLNS.RB_D.cols()) = babyCLNS.RB_D_adv;
*/
		adv_matrix = babyCLNS.Get_advection_matrix();
		
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
		for (int row_index=0; row_index < no_adv_matrix.rows(); ++row_index)
		{
			for (int col_index=0; col_index < no_adv_matrix.cols(); ++col_index)
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

		// generate projected form 
		// using Ritz-Galerkin without any supremizers
		// the projector is 'c_f_all_PODmodes_wo_dbc'
		Eigen::MatrixXd adv_mat_proj = c_f_all_PODmodes_wo_dbc.transpose() * adv_matrix_simplified * c_f_all_PODmodes_wo_dbc;
//		cout << "adv_mat_proj.rows() " << adv_mat_proj.rows() << endl;
//		cout << "adv_mat_proj.cols() " << adv_mat_proj.cols() << endl;
		Eigen::VectorXd adv_rhs_proj = c_f_all_PODmodes_wo_dbc.transpose() * adv_rhs_add;


//		cout << "checking Mat to be saved: " << adv_mat_proj << endl;
		adv_mats_proj_x[trafo_iter] = adv_mat_proj;
		adv_vec_proj_x[trafo_iter] = adv_rhs_proj;

		babyCLNS.DoInitialiseAdv(PhysBase_zero , curr_PhysBaseVec_y ); // call with parameter in phys state
/*		adv_matrix.block(0, 0, babyCLNS.RB_A.rows(), babyCLNS.RB_A.cols()) = MtM * babyCLNS.RB_A_adv;
		adv_matrix.block(0, babyCLNS.RB_A.cols(), babyCLNS.RB_Dbnd.cols(), babyCLNS.RB_Dbnd.rows()) = -MtM * babyCLNS.RB_Dbnd.transpose();
		adv_matrix.block(0, babyCLNS.RB_A.cols() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_B.rows(), babyCLNS.RB_B.cols()) = MtM * babyCLNS.RB_B_adv;
		adv_matrix.block(babyCLNS.RB_A.rows(), 0, babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_Dbnd.cols()) = -babyCLNS.RB_Dbnd;
		adv_matrix.block(babyCLNS.RB_A.rows(), babyCLNS.RB_A.cols() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_Dint.rows(), babyCLNS.RB_Dint.cols()) = -babyCLNS.RB_Dint;	
		adv_matrix.block(babyCLNS.RB_A.rows() + babyCLNS.RB_Dbnd.rows(), 0, babyCLNS.RB_C.cols(), babyCLNS.RB_C.rows()) = babyCLNS.RB_C_adv.transpose();
		adv_matrix.block(babyCLNS.RB_A.rows() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_A.cols(), babyCLNS.RB_Dint.cols(), babyCLNS.RB_Dint.rows()) = -babyCLNS.RB_Dint.transpose();
		adv_matrix.block(babyCLNS.RB_A.rows() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_A.cols() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_D.rows(), babyCLNS.RB_D.cols()) = babyCLNS.RB_D_adv;
*/
		adv_matrix = babyCLNS.Get_advection_matrix();
		add_to_rhs_adv = adv_matrix * f_bnd_dbc_full_size;   
		counter_row_simplified = 0;
		counter_col_simplified = 0;
		for (int row_index=0; row_index < no_adv_matrix.rows(); ++row_index)
		{
			for (int col_index=0; col_index < no_adv_matrix.cols(); ++col_index)
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

	// do a bunch of check-ups to get an idea if all is good so far
	// primo: verify with the full order model all the snapshot data -- is this even possible?
	// should it give zero without the f_bnd_dbc ?

	// have the no_adv_matrix at the reference \nu = 1
//	cout << "no_adv_matrix.rows() " << no_adv_matrix.rows() << endl;
//	cout << "no_adv_matrix.cols() " << no_adv_matrix.cols() << endl;

	// have the snapshot in bnd / p / int format
	//	Eigen::MatrixXd collect_f_all( babyCLNS_trafo.curr_f_bnd.size()+babyCLNS_trafo.curr_f_p.size()+babyCLNS_trafo.curr_f_int.size() , Nmax );
//	cout << "collect_f_all.rows() " << collect_f_all.rows() << endl;
//	cout << "collect_f_all.cols() " << collect_f_all.cols() << endl;

	// have the adv_matrix corresponding to the snapshot 
	//	Array<OneD, Array<OneD, NekDouble> > snapshot_x_collection(snapshots_to_be_collected_aka_Nmax);
	//	Array<OneD, Array<OneD, NekDouble> > snapshot_y_collection(snapshots_to_be_collected_aka_Nmax);
//	cout << "snapshot_x_collection.num_elements() " << snapshot_x_collection.num_elements() << endl;
//	cout << "snapshot_x_collection[0].num_elements() " << snapshot_x_collection[0].num_elements() << endl;

	// have the full_dbc_vector for prior subtraction
	//	Eigen::VectorXd f_bnd_dbc_full_size(collect_f_all_PODmodes.rows());
//	cout << "f_bnd_dbc_full_size.rows() " << f_bnd_dbc_full_size.rows() << endl;
//	cout << "f_bnd_dbc_full_size.cols() " << f_bnd_dbc_full_size.cols() << endl;

	// the \nu = 1 case is the third snapshot
	//int check_index = 2;
	int check_index = 3;
	double current_nu = param_vector[check_index];
	babyCLNS.Set_m_kinvis( current_nu );	
//	cout << "babyCLNS.Get_m_kinvis " << babyCLNS.Get_m_kinvis() << endl;

	// gen the adv mat
//	babyCLNS.InitObject();
	babyCLNS.DoInitialiseAdv(snapshot_x_collection[check_index], snapshot_y_collection[check_index]); 
	Eigen::MatrixXd adv_matrix = Eigen::MatrixXd::Zero(babyCLNS.RB_A.rows() + babyCLNS.RB_Dbnd.rows() + babyCLNS.RB_C.cols(), babyCLNS.RB_A.cols() + babyCLNS.RB_Dbnd.rows() + babyCLNS.RB_B.cols() );
	Eigen::MatrixXd mat_totale = Eigen::MatrixXd::Zero(babyCLNS.RB_A.rows() + babyCLNS.RB_Dbnd.rows() + babyCLNS.RB_C.cols(), babyCLNS.RB_A.cols() + babyCLNS.RB_Dbnd.rows() + babyCLNS.RB_B.cols() );
	adv_matrix.block(0, 0, babyCLNS.RB_A.rows(), babyCLNS.RB_A.cols()) = MtM * babyCLNS.RB_A_adv;
	adv_matrix.block(0, babyCLNS.RB_A.cols(), babyCLNS.RB_Dbnd.cols(), babyCLNS.RB_Dbnd.rows()) = -MtM * babyCLNS.RB_Dbnd.transpose();
	adv_matrix.block(0, babyCLNS.RB_A.cols() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_B.rows(), babyCLNS.RB_B.cols()) = MtM * babyCLNS.RB_B_adv;
	adv_matrix.block(babyCLNS.RB_A.rows(), 0, babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_Dbnd.cols()) = -babyCLNS.RB_Dbnd;
	adv_matrix.block(babyCLNS.RB_A.rows(), babyCLNS.RB_A.cols() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_Dint.rows(), babyCLNS.RB_Dint.cols()) = -babyCLNS.RB_Dint;	
	adv_matrix.block(babyCLNS.RB_A.rows() + babyCLNS.RB_Dbnd.rows(), 0, babyCLNS.RB_C.cols(), babyCLNS.RB_C.rows()) = babyCLNS.RB_C_adv.transpose();
	adv_matrix.block(babyCLNS.RB_A.rows() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_A.cols(), babyCLNS.RB_Dint.cols(), babyCLNS.RB_Dint.rows()) = -babyCLNS.RB_Dint.transpose();
	adv_matrix.block(babyCLNS.RB_A.rows() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_A.cols() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_D.rows(), babyCLNS.RB_D.cols()) = babyCLNS.RB_D_adv;

	mat_totale.block(0, 0, babyCLNS.RB_A.rows(), babyCLNS.RB_A.cols()) = MtM * babyCLNS.RB_A;
	mat_totale.block(0, babyCLNS.RB_A.cols(), babyCLNS.RB_Dbnd.cols(), babyCLNS.RB_Dbnd.rows()) = -MtM * babyCLNS.RB_Dbnd.transpose();
	mat_totale.block(0, babyCLNS.RB_A.cols() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_B.rows(), babyCLNS.RB_B.cols()) = MtM * babyCLNS.RB_B;
	mat_totale.block(babyCLNS.RB_A.rows(), 0, babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_Dbnd.cols()) = -babyCLNS.RB_Dbnd;
	mat_totale.block(babyCLNS.RB_A.rows(), babyCLNS.RB_A.cols() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_Dint.rows(), babyCLNS.RB_Dint.cols()) = -babyCLNS.RB_Dint;	
	mat_totale.block(babyCLNS.RB_A.rows() + babyCLNS.RB_Dbnd.rows(), 0, babyCLNS.RB_C.cols(), babyCLNS.RB_C.rows()) = babyCLNS.RB_C.transpose();
	mat_totale.block(babyCLNS.RB_A.rows() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_A.cols(), babyCLNS.RB_Dint.cols(), babyCLNS.RB_Dint.rows()) = -babyCLNS.RB_Dint.transpose();
	mat_totale.block(babyCLNS.RB_A.rows() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_A.cols() + babyCLNS.RB_Dbnd.rows(), babyCLNS.RB_D.rows(), babyCLNS.RB_D.cols()) = babyCLNS.RB_D;

//	cout << "adv_matrix.rows() " << adv_matrix.rows() << endl;
//	cout << "adv_matrix.cols() " << adv_matrix.cols() << endl;

	Eigen::VectorXd snap_minus_dbc = collect_f_all.col(check_index) - f_bnd_dbc_full_size;
	Eigen::MatrixXd mat_totale2 = no_adv_matrix + adv_matrix;
	// might need to zero the dbc rows and cols
	for (int index=0; index < collect_f_all_PODmodes.rows(); ++index)
	{
		if (elem_loc_dbc.count(index))
		{
//			mat_totale.col(index) = Eigen::VectorXd::Zero(f_bnd_dbc_full_size.rows());
//			mat_totale.row(index) = Eigen::VectorXd::Zero(f_bnd_dbc_full_size.rows());
		}
	}
	Eigen::VectorXd zero_check = mat_totale * collect_f_all.col(check_index);
	Eigen::VectorXd zero_check2 = mat_totale * snap_minus_dbc;
	Eigen::VectorXd zero_check_dbc = mat_totale * f_bnd_dbc_full_size;
	Eigen::MatrixXd mat_compare = Eigen::MatrixXd::Zero(f_bnd_dbc_full_size.rows(),3);
	mat_compare.col(0) = zero_check;
	mat_compare.col(1) = zero_check2;
//	mat_compare.col(2) = zero_check_dbc;
//	cout << mat_compare << endl;
	// BUT: need to look at it under the "action of M", i.e., Mtrafo  -- maybe cannot use the MtM approach when writing the matrix ?????
//	cout << "Mtrafo.rows() " << Mtrafo.rows() << endl; // it is an eigen matrix
//	cout << "Mtrafo.cols() " << Mtrafo.cols() << endl;	
	Eigen::VectorXd check_more = zero_check.head(Mtrafo.rows());
	Eigen::VectorXd check_more_result = Mtrafo.transpose() * check_more;
//	cout << check_more_result << endl;
	
	// why not work -- try in the projected world
	// using directly the snapshots from collect_f_all as projection
	// use zero_check_dbc as r-h-s
	Eigen::MatrixXd mat_totale_proj = collect_f_all.transpose() * mat_totale * collect_f_all;
	Eigen::VectorXd r_h_s_proj = collect_f_all.transpose() *  zero_check_dbc;
	Eigen::VectorXd solve_test = mat_totale_proj.colPivHouseholderQr().solve(r_h_s_proj);
//	cout << "solve_test " << solve_test << endl;
	Eigen::VectorXd reproj_solve_test = collect_f_all * solve_test;
//	cout << "(reproj_solve_test-collect_f_all.col(check_index)).norm() " << (reproj_solve_test-collect_f_all.col(check_index)).norm() << endl;
	mat_compare.col(0) = collect_f_all.col(check_index);
	mat_compare.col(1) = reproj_solve_test;
//	cout << mat_compare << endl;

	// better check with the divided mats
	// be sure to have the proper mat_totale of a fixed snapshot point!
	// is the zero_check_dbc   :: add_to_rhs = mat_totale * f_bnd_dbc_full_size;   
	counter_row_simplified = 0;
	counter_col_simplified = 0;
	Eigen::MatrixXd mat_totale_simplified = Eigen::MatrixXd::Zero(adv_matrix.rows() - no_dbc_in_loc, adv_matrix.cols() - no_dbc_in_loc);
	Eigen::VectorXd adv_rhs_add = Eigen::VectorXd::Zero(adv_matrix.rows() - no_dbc_in_loc);
	for (int row_index=0; row_index < no_adv_matrix.rows(); ++row_index)
	{
		for (int col_index=0; col_index < no_adv_matrix.cols(); ++col_index)
		{
			if ((!elem_loc_dbc.count(row_index)) && (!elem_loc_dbc.count(col_index)))
			{
				mat_totale_simplified(counter_row_simplified, counter_col_simplified) = mat_totale(row_index, col_index);
				counter_col_simplified++;
			}
		}
		counter_col_simplified = 0;
		if (!elem_loc_dbc.count(row_index))
		{
			adv_rhs_add(counter_row_simplified) = zero_check_dbc(row_index);
			counter_row_simplified++;
		}
	}
	Eigen::MatrixXd mat_proj = c_f_all_PODmodes_wo_dbc.transpose() * mat_totale_simplified * c_f_all_PODmodes_wo_dbc;
	Eigen::VectorXd rhs_proj = c_f_all_PODmodes_wo_dbc.transpose() * adv_rhs_add;
	solve_test = mat_proj.colPivHouseholderQr().solve(rhs_proj);
//	cout << "solve_test " << solve_test << endl;
	Eigen::VectorXd agg_reproj_solve_test = c_f_all_PODmodes_wo_dbc * solve_test;
	//cout << agg_reproj_solve_test << endl;
	mat_compare.col(0) = collect_f_all.col(check_index);
	mat_compare.block(0,1,agg_reproj_solve_test.rows() ,1) = agg_reproj_solve_test; // sembra abbastanza bene
	// cout << mat_compare << endl;
	// reconstruct the full one :)
	Eigen::VectorXd reconstruct_solution = Eigen::VectorXd::Zero(adv_matrix.rows());
	int counter_wo_dbc = 0;
	for (int row_index=0; row_index < no_adv_matrix.rows(); ++row_index)
	{
		if (!elem_loc_dbc.count(row_index))
		{
			reconstruct_solution(row_index) = agg_reproj_solve_test(counter_wo_dbc);
			counter_wo_dbc++;
		}
		else
		{
			reconstruct_solution(row_index) = -f_bnd_dbc_full_size(row_index);
		}
	}
	mat_compare.col(0) = collect_f_all.col(check_index);
	mat_compare.col(1) = reconstruct_solution; // sembra abbastanza bene
	mat_compare.col(2) = mat_compare.col(1) + mat_compare.col(0);
	//cout << mat_compare << endl;
//	cout << "relative error norm: " << mat_compare.col(2).norm() / mat_compare.col(0).norm() << endl;
	


	// then should be good for the online solve - in regard to the Q_a term, keep the all structure
	// so what do I need? 
	// since at this point, I have all the projected parameter-independent parts,
	// collect the parameter-dependent matrices 


	// steps to do in particolare
	// given a current iterate, look at the projection coefficients
	// form the sum of the advection term with the projection coefficients
	// \nu dependent are the non-advection, non-pressure terms, so need to differentiate those

	// have an online solve
	// use the affine decomposition
	// for the matrix side, use the \nu-independent part re-use wise 
	// do the looping for all snapshot locations (need some more data to verify intermediate locations)
	int number_of_param_locations = 4;
	int current_index = 0;
	current_nu = param_vector[current_index];
	babyCLNS.Set_m_kinvis( current_nu );	
//	cout << "babyCLNS.Get_m_kinvis " << babyCLNS.Get_m_kinvis() << endl;
	babyCLNS.DoInitialiseAdv(snapshot_x_collection[check_index], snapshot_y_collection[check_index]);
//	have these functions:
//	Eigen::MatrixXd Get_no_advection_matrix(void);
//	Eigen::MatrixXd Get_advection_matrix(void);
//	Eigen::MatrixXd Get_complete_matrix(void);
	Eigen::MatrixXd no_adv_matrix1 = babyCLNS.Get_no_advection_matrix();
	Eigen::MatrixXd no_adv_matrix1_ABCD = babyCLNS.Get_no_advection_matrix_ABCD();
	Eigen::MatrixXd adv_matrix1 = babyCLNS.Get_advection_matrix();
	Eigen::MatrixXd matrix1 = babyCLNS.Get_complete_matrix();

	current_index = 2;
	current_nu = param_vector[current_index];
	babyCLNS.Set_m_kinvis( current_nu );	
//	cout << "babyCLNS.Get_m_kinvis " << babyCLNS.Get_m_kinvis() << endl;
//	babyCLNS.DoInitialiseAdv(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);
	// only checking the x-recovery a.t.m.
	babyCLNS.DoInitialiseAdv(snapshot_x_collection[current_index], PhysBase_zero);
	Eigen::MatrixXd no_adv_matrix2 = babyCLNS.Get_no_advection_matrix();
	Eigen::MatrixXd no_adv_matrix2_ABCD = babyCLNS.Get_no_advection_matrix_ABCD();
	Eigen::MatrixXd adv_matrix2 = babyCLNS.Get_advection_matrix();
	Eigen::MatrixXd matrix2 = babyCLNS.Get_complete_matrix();

	// the no_adv matrices should be the same
	
	// need to treat the Dint and Dbnd independently
//	cout << "(100*no_adv_matrix1_ABCD).norm() " << (100*no_adv_matrix1_ABCD).norm() << endl;
//	cout << "no_adv_matrix2_ABCD.norm() " << no_adv_matrix2_ABCD.norm() << endl;
//	cout << "(no_adv_matrix1_ABCD*100 - no_adv_matrix2_ABCD).norm() " << (no_adv_matrix1_ABCD*100 - no_adv_matrix2_ABCD).norm() << endl;

	// assume no body force for now
	// choose/sceglie a reference parameter configuration
	// the two no_adv terms are clear
	// how about the adv ones?  -- compute all the mats corresponding to all phys x/y bases and do the dbc kicking, then project
	// these things have been created as:
	// Array<OneD, Eigen::MatrixXd> adv_mats_proj_x(RBsize);
	// Array<OneD, Eigen::MatrixXd> adv_mats_proj_y(RBsize);
	// Array<OneD, Eigen::VectorXd> adv_vec_proj_x(RBsize);
	// Array<OneD, Eigen::VectorXd> adv_vec_proj_y(RBsize);
	// thus: project the current iterate onto the phys basis space:
	// in this case, init with a snapshot:
	// need to project snapshot_x_collection[check_index], snapshot_y_collection[check_index]
	// all involved are in phys so, should be no big deal, can also check the projection accuracy
	// probably need the init snapshot as an Eigen::VectorXd
	Eigen::VectorXd current_snapshot_x(babyCLNS.GetNpoints());
	Eigen::VectorXd current_snapshot_y(babyCLNS.GetNpoints());
	// or try the map functionality
	// maybe works only for matrices	Eigen::VectorXd mapped_current_snapshot_x = Eigen::Map<Eigen::VectorXd>(snapshot_x_collection[check_index],babyCLNS.GetNpoints());

	for (int i = 0; i < babyCLNS.GetNpoints(); ++i)
	{
		current_snapshot_x(i) = snapshot_x_collection[current_index][i];
		current_snapshot_y(i) = snapshot_y_collection[current_index][i];
	}
	// could do this in a function:	eigen_phys_basis_ x/y	
/*	cout << "current_snapshot_x.rows() " << current_snapshot_x.rows() << endl;
	cout << "current_snapshot_x.cols() " << current_snapshot_x.cols() << endl;
	cout << "current_snapshot_y.rows() " << current_snapshot_y.rows() << endl;
	cout << "current_snapshot_y.cols() " << current_snapshot_y.cols() << endl;
	cout << "eigen_phys_basis_x.rows() " << eigen_phys_basis_x.rows() << endl;
	cout << "eigen_phys_basis_x.cols() " << eigen_phys_basis_x.cols() << endl;
	cout << "eigen_phys_basis_y.rows() " << eigen_phys_basis_y.rows() << endl;
	cout << "eigen_phys_basis_y.cols() " << eigen_phys_basis_y.cols() << endl;
*/
	Eigen::VectorXd curr_x_proj = eigen_phys_basis_x.transpose() * current_snapshot_x;
	Eigen::VectorXd curr_y_proj = eigen_phys_basis_y.transpose() * current_snapshot_y;
//	cout << "curr_x_proj " << curr_x_proj << endl;
//	cout << "curr_y_proj " << curr_y_proj << endl;

//	cout << "eigen_phys_basis_x.transpose() * eigen_phys_basis_x " << eigen_phys_basis_x.transpose() * eigen_phys_basis_x << endl;
//	cout << "eigen_phys_basis_y.transpose() * eigen_phys_basis_y " << eigen_phys_basis_y.transpose() * eigen_phys_basis_y << endl;

	// check the accuracy of the recovered (affine-wise) projected advection matrix
	// how about the adv ones?  -- compute all the mats corresponding to all phys x/y bases and do the dbc kicking, then project
	// these things have been created as:
	// Array<OneD, Eigen::MatrixXd> adv_mats_proj_x(RBsize);
	// Array<OneD, Eigen::MatrixXd> adv_mats_proj_y(RBsize);
	// Array<OneD, Eigen::VectorXd> adv_vec_proj_x(RBsize);
	// Array<OneD, Eigen::VectorXd> adv_vec_proj_y(RBsize);
	Eigen::MatrixXd recovered_affine_adv_mat_proj_x = Eigen::MatrixXd::Zero(RBsize, RBsize);
	for (int i = 0; i < RBsize; ++i)
	{
		recovered_affine_adv_mat_proj_x += adv_mats_proj_x[i] * curr_x_proj[i];
//		cout << "adv_mats_proj_x[i] " << adv_mats_proj_x[i] << endl;
	}
//	cout << "recovered_affine_adv_mat_proj_x " << recovered_affine_adv_mat_proj_x << endl;
	// compare this to a projected adv_matrix2 
	// have to generate the zero_check_dbc beforehand 	Eigen::VectorXd zero_check_dbc = mat_totale * f_bnd_dbc_full_size;
	// do functionalise this!!
	counter_row_simplified = 0;
	counter_col_simplified = 0;
	mat_totale_simplified = Eigen::MatrixXd::Zero(adv_matrix.rows() - no_dbc_in_loc, adv_matrix.cols() - no_dbc_in_loc);
	adv_rhs_add = Eigen::VectorXd::Zero(adv_matrix.rows() - no_dbc_in_loc);
	for (int row_index=0; row_index < no_adv_matrix.rows(); ++row_index)
	{
		for (int col_index=0; col_index < no_adv_matrix.cols(); ++col_index)
		{
			if ((!elem_loc_dbc.count(row_index)) && (!elem_loc_dbc.count(col_index)))
			{
				mat_totale_simplified(counter_row_simplified, counter_col_simplified) = adv_matrix2(row_index, col_index);
				counter_col_simplified++;
			}
		}
		counter_col_simplified = 0;
		if (!elem_loc_dbc.count(row_index))
		{
//			adv_rhs_add(counter_row_simplified) = zero_check_dbc(row_index); have to generate the zero_check_dbc beforehand
			counter_row_simplified++;
		}
	}
	mat_proj = c_f_all_PODmodes_wo_dbc.transpose() * mat_totale_simplified * c_f_all_PODmodes_wo_dbc;
//	Eigen::VectorXd rhs_proj = c_f_all_PODmodes_wo_dbc.transpose() * adv_rhs_add;
//	cout << "mat_proj " << mat_proj << endl;

// 	now check the y-part
	Eigen::MatrixXd recovered_affine_adv_mat_proj_y = Eigen::MatrixXd::Zero(RBsize, RBsize);
	for (int i = 0; i < RBsize; ++i)
	{
		recovered_affine_adv_mat_proj_y += adv_mats_proj_y[i] * curr_y_proj[i];
//		cout << "adv_mats_proj_x[i] " << adv_mats_proj_x[i] << endl;
	}
//	cout << "recovered_affine_adv_mat_proj_y " << recovered_affine_adv_mat_proj_y << endl;

	babyCLNS.DoInitialiseAdv(PhysBase_zero, snapshot_y_collection[current_index]);
	Eigen::MatrixXd adv_matrix3 = babyCLNS.Get_advection_matrix();
	counter_row_simplified = 0;
	counter_col_simplified = 0;
	mat_totale_simplified = Eigen::MatrixXd::Zero(adv_matrix.rows() - no_dbc_in_loc, adv_matrix.cols() - no_dbc_in_loc);
	adv_rhs_add = Eigen::VectorXd::Zero(adv_matrix.rows() - no_dbc_in_loc);
	for (int row_index=0; row_index < no_adv_matrix.rows(); ++row_index)
	{
		for (int col_index=0; col_index < no_adv_matrix.cols(); ++col_index)
		{
			if ((!elem_loc_dbc.count(row_index)) && (!elem_loc_dbc.count(col_index)))
			{
				mat_totale_simplified(counter_row_simplified, counter_col_simplified) = adv_matrix3(row_index, col_index);
				counter_col_simplified++;
			}
		}
		counter_col_simplified = 0;
		if (!elem_loc_dbc.count(row_index))
		{
			counter_row_simplified++;
		}
	}
	mat_proj = c_f_all_PODmodes_wo_dbc.transpose() * mat_totale_simplified * c_f_all_PODmodes_wo_dbc;
//	cout << "mat_proj " << mat_proj << endl;

//	part totale
	babyCLNS.DoInitialiseAdv(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);
	Eigen::MatrixXd adv_matrix4 = babyCLNS.Get_advection_matrix();
/*	counter_row_simplified = 0;
	counter_col_simplified = 0;
	mat_totale_simplified = Eigen::MatrixXd::Zero(adv_matrix.rows() - no_dbc_in_loc, adv_matrix.cols() - no_dbc_in_loc);
	adv_rhs_add = Eigen::VectorXd::Zero(adv_matrix.rows() - no_dbc_in_loc);
	for (int row_index=0; row_index < no_adv_matrix.rows(); ++row_index)
	{
		for (int col_index=0; col_index < no_adv_matrix.cols(); ++col_index)
		{
			if ((!elem_loc_dbc.count(row_index)) && (!elem_loc_dbc.count(col_index)))
			{
				mat_totale_simplified(counter_row_simplified, counter_col_simplified) = adv_matrix4(row_index, col_index);
				counter_col_simplified++;
			}
		}
		counter_col_simplified = 0;
		if (!elem_loc_dbc.count(row_index))
		{
			counter_row_simplified++;
		}
	}
*/
	mat_totale_simplified = babyCLNS.remove_cols_and_rows(adv_matrix4, elem_loc_dbc);

	mat_proj = c_f_all_PODmodes_wo_dbc.transpose() * mat_totale_simplified * c_f_all_PODmodes_wo_dbc;
//	cout << "mat_proj " << mat_proj << endl;

	Eigen::MatrixXd recovered_affine_adv_mat_proj_xy = Eigen::MatrixXd::Zero(RBsize, RBsize);
	for (int i = 0; i < RBsize; ++i)
	{
		recovered_affine_adv_mat_proj_xy += adv_mats_proj_x[i] * curr_x_proj[i] + adv_mats_proj_y[i] * curr_y_proj[i];
//		cout << "adv_mats_proj_x[i] " << adv_mats_proj_x[i] << endl;
	}
//	cout << "recovered_affine_adv_mat_proj_xy " << recovered_affine_adv_mat_proj_xy << endl;


//	get it all together:
	current_index = 2; // i.e., \nu = 1
	current_nu = param_vector[current_index];
	babyCLNS.Set_m_kinvis( current_nu );
	babyCLNS.DoInitialiseAdv(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);
	Eigen::MatrixXd test_complete_matrix = babyCLNS.Get_complete_matrix();
	mat_totale_simplified = babyCLNS.remove_cols_and_rows(test_complete_matrix, elem_loc_dbc);
	mat_proj = c_f_all_PODmodes_wo_dbc.transpose() * mat_totale_simplified * c_f_all_PODmodes_wo_dbc;
//	cout << "non-affine mat_proj " << mat_proj << endl;
	// can re-use recovered_affine_adv_mat_proj_xy
	// two more terms const one depending pressure and \nu one
	Eigen::MatrixXd the_const_one = babyCLNS.Get_no_advection_matrix_pressure();
	Eigen::MatrixXd the_ABCD_one = babyCLNS.Get_no_advection_matrix_ABCD();
	Eigen::MatrixXd the_const_one_simplified = babyCLNS.remove_cols_and_rows(the_const_one, elem_loc_dbc);
	Eigen::MatrixXd the_ABCD_one_simplified = babyCLNS.remove_cols_and_rows(the_ABCD_one, elem_loc_dbc);
	Eigen::MatrixXd the_const_one_proj = c_f_all_PODmodes_wo_dbc.transpose() * the_const_one_simplified * c_f_all_PODmodes_wo_dbc;
	Eigen::MatrixXd the_ABCD_one_proj = c_f_all_PODmodes_wo_dbc.transpose() * the_ABCD_one_simplified * c_f_all_PODmodes_wo_dbc;
//	cout << "affine mat_proj " << recovered_affine_adv_mat_proj_xy + the_const_one_proj + the_ABCD_one_proj << endl;	
	Eigen::VectorXd the_ABCD_one_rhs = the_ABCD_one * f_bnd_dbc_full_size;
	Eigen::VectorXd the_const_one_rhs = the_const_one * f_bnd_dbc_full_size;
	Eigen::VectorXd the_ABCD_one_rhs_simplified = babyCLNS.remove_rows(the_ABCD_one_rhs, elem_loc_dbc);
	Eigen::VectorXd the_const_one_rhs_simplified = babyCLNS.remove_rows(the_const_one_rhs, elem_loc_dbc);
	Eigen::VectorXd the_ABCD_one_rhs_proj = c_f_all_PODmodes_wo_dbc.transpose() * the_ABCD_one_rhs_simplified;
	Eigen::VectorXd the_const_one_rhs_proj = c_f_all_PODmodes_wo_dbc.transpose() * the_const_one_rhs_simplified;


	// start sweeping 
	for (int iter_index = 0; iter_index < Nmax; ++iter_index)
	{
		current_index = iter_index;
		current_nu = param_vector[current_index];
		babyCLNS.Set_m_kinvis( current_nu );
		babyCLNS.DoInitialiseAdv(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);
		Eigen::MatrixXd test_complete_matrix = babyCLNS.Get_complete_matrix();
//		Eigen::MatrixXd test_complete_matrix = babyCLNS.Get_advection_matrix();
		mat_totale_simplified = babyCLNS.remove_cols_and_rows(test_complete_matrix, elem_loc_dbc);
		mat_proj = c_f_all_PODmodes_wo_dbc.transpose() * mat_totale_simplified * c_f_all_PODmodes_wo_dbc;
//		cout << "sweep non-affine mat_proj " << mat_proj << endl;
		Eigen::VectorXd test_complete_dbc = test_complete_matrix * f_bnd_dbc_full_size;
		Eigen::VectorXd test_complete_dbc_simplified = babyCLNS.remove_rows(test_complete_dbc, elem_loc_dbc);
		Eigen::VectorXd rhs_proj = c_f_all_PODmodes_wo_dbc.transpose() * test_complete_dbc_simplified;
//		cout << "sweep non-affine rhs_proj " << rhs_proj << endl;
		Eigen::VectorXd solve_non_affine = mat_proj.colPivHouseholderQr().solve(rhs_proj);
//		cout << "solve_non_affine " << solve_non_affine << endl;

		// affine part
		Eigen::VectorXd c_snapshot_x(babyCLNS.GetNpoints());
		Eigen::VectorXd c_snapshot_y(babyCLNS.GetNpoints());
		for (int i = 0; i < babyCLNS.GetNpoints(); ++i)
		{
			c_snapshot_x(i) = snapshot_x_collection[current_index][i];
			c_snapshot_y(i) = snapshot_y_collection[current_index][i];
		}
		curr_x_proj = eigen_phys_basis_x.transpose() * c_snapshot_x;
		curr_y_proj = eigen_phys_basis_y.transpose() * c_snapshot_y;
		recovered_affine_adv_mat_proj_xy = Eigen::MatrixXd::Zero(RBsize, RBsize);
		for (int i = 0; i < RBsize; ++i)
		{
			recovered_affine_adv_mat_proj_xy += adv_mats_proj_x[i] * curr_x_proj[i] + adv_mats_proj_y[i] * curr_y_proj[i];
		}
//		Eigen::MatrixXd affine_mat_proj = recovered_affine_adv_mat_proj_xy;
		Eigen::MatrixXd affine_mat_proj = the_const_one_proj + current_nu * the_ABCD_one_proj + recovered_affine_adv_mat_proj_xy;
//		cout << "sweep affine mat_proj " << affine_mat_proj << endl;

		// do deal with the right-hand-side
		// central correlation: 	Eigen::VectorXd zero_check_dbc = mat_totale * f_bnd_dbc_full_size;
		// already have computed:	Array<OneD, Eigen::VectorXd> adv_vec_proj_x(RBsize); 	Array<OneD, Eigen::VectorXd> adv_vec_proj_y(RBsize);
		Eigen::VectorXd recovered_affine_adv_rhs_proj_xy = Eigen::VectorXd::Zero(RBsize); 
		for (int i = 0; i < RBsize; ++i)
		{
			recovered_affine_adv_rhs_proj_xy += adv_vec_proj_x[i] * curr_x_proj[i] + adv_vec_proj_y[i] * curr_y_proj[i];
		}		
//		cout << "sweep affine rhs_proj " << the_const_one_rhs_proj + current_nu * the_ABCD_one_rhs_proj + recovered_affine_adv_rhs_proj_xy << endl;
		Eigen::VectorXd solve_affine = affine_mat_proj.colPivHouseholderQr().solve(the_const_one_rhs_proj + current_nu * the_ABCD_one_rhs_proj + recovered_affine_adv_rhs_proj_xy);
		cout << "solve_affine " << solve_affine << endl;
		Eigen::VectorXd repro_solve_affine = c_f_all_PODmodes_wo_dbc * solve_affine;
		Eigen::VectorXd reconstruct_solution = Eigen::VectorXd::Zero(adv_matrix.rows());
		int counter_wo_dbc = 0;
		for (int row_index=0; row_index < no_adv_matrix.rows(); ++row_index)
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
		//cout << mat_compare << endl;
		cout << "relative error norm: " << mat_compare.col(2).norm() / mat_compare.col(0).norm() << endl;

	}





//	babyCLNS.DoSolve(); // do I need this here ?
//	babyCLNS.Output();  // do I need this here ?


        session->Finalise();
//        session2->Finalise();


	// la fine:   /////////////////////////////
	return 0;


/*


	cout <<	"session->DefinesParameter(Kinvis_max) " << session->DefinesParameter("Kinvis_max") << endl;
	if (session->DefinesParameter("Kinvis_max"))
	{
		int testp = session->GetParameter("Kinvis_max");
		cout << "testp " << testp << endl;
	}

	Array<OneD, EquationSystemSharedPtr> m_equ;
	m_equ = drv->GetEqu();

	cout << "m_equ.num_elements() " << m_equ.num_elements() << endl;
	cout << "session->GetSolverInfo(SolverType) " << session->GetSolverInfo("SolverType") << endl;
        m_equ[0]->PrintSummary(cout);
//	cout << "boost::lexical_cast<std::string>(m_kinvis) " << boost::lexical_cast<std::string>(m_kinvis) << endl; 

        // Execute driver
// 	drv->Execute();

	m_equ[0]->InitObject();
	m_equ[0]->DoInitialise();
	m_equ[0]->DoSolve();
	m_equ[0]->Output();
	cout << "m_equ[0]->GetNvariables() " << m_equ[0]->GetNvariables() << endl;
	Array<OneD, NekDouble> exactsoln(m_equ[0]->GetTotPoints(), 0.0);
	m_equ[0]->EvaluateExactSolution(0, exactsoln, m_equ[0]->GetFinalTime());
	m_equ[0]->EvaluateExactSolution(1, exactsoln, m_equ[0]->GetFinalTime());

        // Finalise communications
        session->Finalise();

	Array<OneD, MultiRegions::ExpListSharedPtr> m_fields = m_equ[0]->UpdateFields();
	cout << "m_fields[0]->GetNcoeffs() " << m_fields[0]->GetNcoeffs() << endl; // this is nlc
	cout << "m_fields[0]->GetNpoints() " << m_fields[0]->GetNpoints() << endl;  // is this number of phys ??
	cout << "m_fields[0]->GetTotPoints() " << m_fields[0]->GetTotPoints() << endl;
	cout << "m_fields[0]->GetPhysState() " << m_fields[0]->GetPhysState() << endl;

	// map the loc field to the phys field
	m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(), m_fields[0]->UpdatePhys()); // needed if going with the equation system 

	Array<OneD, NekDouble> out_field_0_x(m_equ[0]->GetNpoints(), 0.0);
	Array<OneD, NekDouble> out_field_0_y(m_equ[0]->GetNpoints(), 0.0);
	Array<OneD, NekDouble> out_field_o(m_equ[0]->GetNpoints(), 0.0);
	m_equ[0]->CopyFromPhysField(0, out_field_0_x);  // how would I know the phys is updated ??
	m_equ[0]->CopyFromPhysField(1, out_field_0_y);
        Array<OneD, Array<OneD, NekDouble> > collected_snapshots_x(snapshots_to_be_collected_aka_Nmax);
        Array<OneD, Array<OneD, NekDouble> > collected_snapshots_y(snapshots_to_be_collected_aka_Nmax);
	collected_snapshots_x[0] = out_field_0_x;
	collected_snapshots_y[0] = out_field_0_y;

	Array<OneD, MultiRegions::ExpListSharedPtr> m_fields_o = m_equ[0]->UpdateFields();
	out_field_o = m_fields_o[0]->UpdatePhys();
	// compare input and output
	Eigen::VectorXd csx0_o(collected_snapshots_x[0].num_elements());
	Eigen::VectorXd csx0_oo(collected_snapshots_x[0].num_elements());
	for( int i = 0; i < collected_snapshots_x[0].num_elements(); ++i)
	{
		csx0_o(i) = out_field_o[i];
		csx0_oo(i) = collected_snapshots_x[0][i];
	}
	cout << "o collected_snapshots_x[0] " << csx0_o.norm() << endl;
	cout << "oo collected_snapshots_x[0] " << csx0_oo.norm() << endl;

	Array<OneD, NekDouble> loc_out_field_0_x(m_fields[0]->GetNcoeffs(), 0.0);
	loc_out_field_0_x = m_fields[0]->GetCoeffs();
	Eigen::VectorXd lof0x(loc_out_field_0_x.num_elements());
	for( int i = 0; i < loc_out_field_0_x.num_elements(); ++i)
	{
		lof0x(i) = loc_out_field_0_x[i];
	}
	cout << "norm of the loc field " << lof0x.norm() << endl;











//	for (int i_phys_dof = 0; i_phys_dof < m_equ[0]->GetNpoints(); i_phys_dof++)
//	{
//		cout << out_field_0[i_phys_dof] << "\t";
//	}


	
	cout << endl;
//	cout << "collected_snapshots.num_elements() " << collected_snapshots.num_elements() << endl;
//	cout << "collected_snapshots[0].num_elements() " << collected_snapshots[0].num_elements() << endl;
//	cout << "collected_snapshots[1].num_elements() " << collected_snapshots[1].num_elements() << endl;
//	cout << "collected_snapshots[2].num_elements() " << collected_snapshots[2].num_elements() << endl;

	


	// can I loop here somehow ????????
	const string name = "Kinvis";
	double value = 0.001;
	session->SetParameter(name, value);
        drv = GetDriverFactory().CreateInstance(vDriverModule, session);
	drv->Execute();
        session->Finalise();

	m_equ = drv->GetEqu();
	Array<OneD, NekDouble> out_field_1_x(m_equ[0]->GetNpoints(), 0.0);
	Array<OneD, NekDouble> out_field_1_y(m_equ[0]->GetNpoints(), 0.0);

	m_equ[0]->CopyFromPhysField(0, out_field_1_x);  // how would I know the phys is updated ??
	m_equ[0]->CopyFromPhysField(1, out_field_1_y);

//	for (int i_phys_dof = 0; i_phys_dof < m_equ[0]->GetNpoints(); i_phys_dof++)
//	{
//		cout << out_field_1[i_phys_dof] << "\t";
//	}

	collected_snapshots_x[1] = out_field_1_x;
	collected_snapshots_y[1] = out_field_1_y;

	value = 10000;
	session->SetParameter(name, value);
        drv = GetDriverFactory().CreateInstance(vDriverModule, session);
	drv->Execute();
        session->Finalise();

	m_equ = drv->GetEqu();
	Array<OneD, NekDouble> out_field_2_x(m_equ[0]->GetNpoints(), 0.0);
	Array<OneD, NekDouble> out_field_2_y(m_equ[0]->GetNpoints(), 0.0);
	m_equ[0]->CopyFromPhysField(0, out_field_2_x);  // how would I know the phys is updated ??
	m_equ[0]->CopyFromPhysField(1, out_field_2_y);

//	for (int i_phys_dof = 0; i_phys_dof < m_equ[0]->GetNpoints(); i_phys_dof++)
//	{
//		cout << out_field_2[i_phys_dof] << "\t";
//	}

	collected_snapshots_x[2] = out_field_2_x;
	collected_snapshots_y[2] = out_field_2_y;

	cout << "m_equ[0]->GetNcoeffs() " << m_equ[0]->GetNcoeffs() << endl; // this is nlc
	cout << "m_equ[0]->GetNpoints() " << m_equ[0]->GetNpoints() << endl;  // is this number of phys ??
	cout << "m_equ[0]->GetTotPoints() " << m_equ[0]->GetTotPoints() << endl;
	cout << "m_equ[0]->GetNvariables() " << m_equ[0]->GetNvariables() << endl;

//	Array<OneD, NekDouble> out_field_0(m_equ[0]->GetNpoints(), 0.0);
//	m_equ[0]->CopyFromPhysField(0, out_field_0);  // how would I know the phys is updated ??

//	for (int i_phys_dof = 0; i_phys_dof < m_equ[0]->GetNpoints(); i_phys_dof++)
//	{
//		cout << out_field_0[i_phys_dof] << "\t";
//	}


	// can also get the fields with:
//	Array<OneD, MultiRegions::ExpListSharedPtr> m_fields = m_equ[0]->UpdateFields();
	cout << "m_fields[0]->GetNcoeffs() " << m_fields[0]->GetNcoeffs() << endl; // this is nlc
	cout << "m_fields[0]->GetNpoints() " << m_fields[0]->GetNpoints() << endl;  // is this number of phys ??
	cout << "m_fields[0]->GetTotPoints() " << m_fields[0]->GetTotPoints() << endl;

	cout << "m_fields[0]->GetPhysState() " << m_fields[0]->GetPhysState() << endl;
//	cout << "collected_snapshots.num_elements() " << collected_snapshots.num_elements() << endl;
//	cout << "collected_snapshots[0].num_elements() " << collected_snapshots[0].num_elements() << endl;
//	cout << "collected_snapshots[1].num_elements() " << collected_snapshots[1].num_elements() << endl;
//	cout << "collected_snapshots[2].num_elements() " << collected_snapshots[2].num_elements() << endl;


	// now do an eigen-powered POD
	Eigen::MatrixXd SnapMatrix(snapshots_to_be_collected_aka_Nmax, m_equ[0]->GetNpoints());
	Eigen::VectorXd v(m_equ[0]->GetNpoints());
	for (int i_phys_dof = 0; i_phys_dof < m_equ[0]->GetNpoints(); i_phys_dof++)
	{
		SnapMatrix(0,i_phys_dof) = collected_snapshots_x[0][i_phys_dof]; // actually did it on v_bnd, f_p, v_int
		SnapMatrix(1,i_phys_dof) = collected_snapshots_x[1][i_phys_dof];
		SnapMatrix(2,i_phys_dof) = collected_snapshots_x[2][i_phys_dof];
	}
	// v << out_field_1;
	// cout << v << endl;

	// cout << SnapMatrix << endl << endl;
	// cout << SnapMatrix.bdcSvd(Eigen::ComputeThinU) << endl;
	Eigen::BDCSVD<Eigen::MatrixXd> svdTT(SnapMatrix, Eigen::ComputeThinV);
	cout << svdTT.singularValues() << endl << endl;
	Eigen::MatrixXd PODmodes = svdTT.matrixV();
	// cout << svdTT.matrixV() << endl << endl;
	// cout << PODmodes << endl << endl;

	// now to the projection business:
	// m_equ[0]->DefineRBspace;
	
	// two options: project in class derived from CoupledNS or here?
	// or in some new class?

	session->SetSolverInfo("SolverType", "CoupledLinearisedNS_TT");
	value = 1;
	session->SetParameter(name, value);
        drv = GetDriverFactory().CreateInstance(vDriverModule, session);

	Array<OneD, EquationSystemSharedPtr> my_equ;
	my_equ = drv->GetEqu();

//	load to hdd and then reload in the init??
// 	it is not the right bnd / p / int


//	CoupledLinearNS_TT *myCLNS;
//	EquationSystem *myCLNS;
//	boost::shared_ptr<EquationSystem> myCLNS;
//	boost::shared_ptr<CoupledLinearNS_TT> myCLNS;
//	myCLNS = my_equ[0];

	//drv->Execute();
        //session->Finalise();

	CoupledLinearNS_TT myCLNS(session);

	myCLNS.DefineRBspace(PODmodes);
	myCLNS.PrintSummary(cout);
	myCLNS.InitObject();
	myCLNS.DoInitialise();
	myCLNS.DoSolve();

	// now should have the A,B,C,D,Dbnd,Dint available
	// create the bnd / p / int snaphots

	Eigen::VectorXd trafo_snapshot;
	Eigen::VectorXd snapshot;
	cout << "PODmodes.cols() " << PODmodes.cols() << endl;
	cout << "PODmodes.rows() " << PODmodes.rows() << endl;
	snapshot = PODmodes.col(1);
	cout << "snapshot.cols() " << snapshot.cols() << endl;
	cout << "snapshot.rows() " << snapshot.rows() << endl;
	double parameter;
	parameter = 0.15;
	// can get the info by doing a solve for each: trafo_snapshot = myCLNS.trafoSnapshot(snapshot, parameter);
	// importante: need to set an appropriate advection field, i.e., the snapshot 'itself' quasi so-to-speak

	// getting the correct AdvField requires to go from local -> phys // can already capture the phys
	CoupledLinearNS_trafoP myCLNS_trafo(session);
	if (session->DefinesParameter("Kinvis"))
	{
		double testp = session->GetParameter("Kinvis");
		cout << "testp " << testp << endl;
	}


	myCLNS_trafo.InitObject();
	myCLNS_trafo.DoInitialiseAdv(collected_snapshots_x[2], collected_snapshots_y[2]); // replaces .DoInitialise();
//	myCLNS_trafo.DoInitialise();
	myCLNS_trafo.DoSolve();
	myCLNS_trafo.Output();


	// now can get the phys as well as bnd / p / int solution :)
	Array<OneD, MultiRegions::ExpListSharedPtr> m_fields_t = myCLNS_trafo.UpdateFields();

	// map the loc field to the phys field
	m_fields_t[0]->BwdTrans(m_fields_t[0]->GetCoeffs(), m_fields_t[0]->UpdatePhys()); // needed if going with the equation system 


	Array<OneD, NekDouble> out_field_trafo_x(m_equ[0]->GetNpoints(), 0.0);
	Array<OneD, NekDouble> out_field_trafo_x2(m_equ[0]->GetNpoints(), 0.0);
	Array<OneD, NekDouble> out_field_trafo_y(m_equ[0]->GetNpoints(), 0.0);
	myCLNS_trafo.CopyFromPhysField(0, out_field_trafo_x); 
	myCLNS_trafo.CopyFromPhysField(1, out_field_trafo_y);
	out_field_trafo_x2 = m_fields_t[0]->UpdatePhys();
	// compare input and output
	Eigen::VectorXd csx0(collected_snapshots_x[0].num_elements());
	Eigen::VectorXd csx0_trafo(collected_snapshots_x[0].num_elements());
	Eigen::VectorXd csx0_trafo2(collected_snapshots_x[0].num_elements());
	Eigen::VectorXd csy0(collected_snapshots_x[0].num_elements());
	Eigen::VectorXd csy0_trafo(collected_snapshots_x[0].num_elements());
	for( int i = 0; i < collected_snapshots_x[0].num_elements(); ++i)
	{
		csx0(i) = collected_snapshots_x[0][i];
		csx0_trafo(i) = out_field_trafo_x[i];
		csx0_trafo2(i) = out_field_trafo_x2[i];
		csy0(i) = collected_snapshots_y[0][i];
		csy0_trafo(i) = out_field_trafo_y[i];
	}
	cout << "collected_snapshots_x[0] " << csx0.norm() << endl;
	cout << "trafo collected_snapshots_x[0] " << csx0_trafo.norm() << endl;
	cout << "trafo2 collected_snapshots_x[0] " << csx0_trafo2.norm() << endl;

//	Eigen::VectorXd trafo_f_bnd = myCLNS_trafo.curr_f_bnd;
//	Eigen::VectorXd trafo_f_p = myCLNS_trafo.curr_f_p;
//	Eigen::VectorXd trafo_f_int = myCLNS_trafo.curr_f_int;

	// do a full order looping
	int no_of_loops = -1;
	for ( int iter = 0; iter < no_of_loops; ++iter)
	{
//		myCLNS_trafo.InitObject();
		myCLNS_trafo.DoInitialiseAdv(out_field_trafo_x, out_field_trafo_y); // replaces .DoInitialise();
//		myCLNS_trafo.DoInitialise();
		myCLNS_trafo.DoSolve();
		myCLNS_trafo.Output();

		m_fields_t = myCLNS_trafo.UpdateFields();
		m_fields_t[0]->BwdTrans(m_fields_t[0]->GetCoeffs(), m_fields_t[0]->UpdatePhys()); 
		myCLNS_trafo.CopyFromPhysField(0, out_field_trafo_x); 
		myCLNS_trafo.CopyFromPhysField(1, out_field_trafo_y);
		for( int i = 0; i < collected_snapshots_x[0].num_elements(); ++i)
		{
			csx0_trafo(i) = out_field_trafo_x[i];
			csy0_trafo(i) = out_field_trafo_y[i];
		}

		cout << "trafo collected_snapshots_x[0] " << csx0_trafo.norm() << endl;

	}

	// do the bnd / p / int transform
	int no_snaps = snapshots_to_be_collected_aka_Nmax;
	Eigen::MatrixXd c_f_bnd( myCLNS_trafo.curr_f_bnd.size() , no_snaps );
	Eigen::MatrixXd c_f_p( myCLNS_trafo.curr_f_p.size() , no_snaps );
	Eigen::MatrixXd c_f_int( myCLNS_trafo.curr_f_int.size() , no_snaps );
	for(int trafo_iter = 0; trafo_iter < no_snaps; trafo_iter++)
	{
		myCLNS_trafo.InitObject();
		myCLNS_trafo.DoInitialiseAdv(collected_snapshots_x[trafo_iter], collected_snapshots_y[trafo_iter]); // replaces .DoInitialise();
		myCLNS_trafo.DoSolve();
		myCLNS_trafo.Output();
	
		c_f_bnd.col(trafo_iter) = myCLNS_trafo.curr_f_bnd;
		c_f_p.col(trafo_iter) = myCLNS_trafo.curr_f_p;
		c_f_int.col(trafo_iter) = myCLNS_trafo.curr_f_int;
		

	}



*/



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
