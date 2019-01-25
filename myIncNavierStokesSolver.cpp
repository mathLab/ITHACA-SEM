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
	cout << "vDriverModule " << vDriverModule << endl;
        // Create driver
        drv = GetDriverFactory().CreateInstance(vDriverModule, session);  // why is that necessary?

	int snapshots_to_be_collected_aka_Nmax = 4;        

	CoupledLinearNS_TT babyCLNS(session);
	babyCLNS.InitObject();
//	babyCLNS.DoInitialise();

	// some IO checks
	cout << "session->DefinesFunction(InitialConditions) " << session->DefinesFunction("InitialConditions") << endl;
	cout << "session->DefinesFunction(TestSnap1) " << session->DefinesFunction("TestSnap1") << endl;
	cout << "session->DefinesFunction(TestSnap2) " << session->DefinesFunction("TestSnap2") << endl;
	cout << "session->DefinesFunction(TestSnap3) " << session->DefinesFunction("TestSnap3") << endl;
	cout << "session->DefinesFunction(TestSnap4) " << session->DefinesFunction("TestSnap4") << endl;
	cout << "session->DefinesFunction(TestSnap5) " << session->DefinesFunction("TestSnap5") << endl;
	cout << "session->DefinesFunction(\"AdvectionVelocity\") " << session->DefinesFunction("AdvectionVelocity") << endl;

	int nvelo = 2;
        Array<OneD, Array<OneD, NekDouble> > test_load_snapshot(nvelo); // for a 2D problem
	Array<OneD, Array<OneD, NekDouble> > snapshot_collection(snapshots_to_be_collected_aka_Nmax);
	Array<OneD, Array<OneD, NekDouble> > snapshot_x_collection(snapshots_to_be_collected_aka_Nmax);
	Array<OneD, Array<OneD, NekDouble> > snapshot_y_collection(snapshots_to_be_collected_aka_Nmax);
	// probably also need this as EigenMatrix

        for(int i = 0; i < nvelo; ++i)
        {
            test_load_snapshot[i] = Array<OneD, NekDouble> (babyCLNS.GetNpoints(), 0.0);  // number of phys points
        }
               
        std::vector<std::string> fieldStr;
        for(int i = 0; i < nvelo; ++i)
        {
           fieldStr.push_back(session->GetVariable(i));
           cout << "session->GetVariable(i) " << session->GetVariable(i) << endl;
        }
        babyCLNS.EvaluateFunction(fieldStr, test_load_snapshot, "TestSnap1");
	snapshot_x_collection[0] = test_load_snapshot[0];
	snapshot_y_collection[0] = test_load_snapshot[1];
        babyCLNS.EvaluateFunction(fieldStr, test_load_snapshot, "TestSnap2");
	snapshot_x_collection[1] = test_load_snapshot[0];
	snapshot_y_collection[1] = test_load_snapshot[1];
        babyCLNS.EvaluateFunction(fieldStr, test_load_snapshot, "TestSnap3");
	snapshot_x_collection[2] = test_load_snapshot[0];
	snapshot_y_collection[2] = test_load_snapshot[1];
        babyCLNS.EvaluateFunction(fieldStr, test_load_snapshot, "TestSnap4");
	snapshot_x_collection[3] = test_load_snapshot[0];
	snapshot_y_collection[3] = test_load_snapshot[1];


	session->SetSolverInfo("SolverType", "CoupledLinearisedNS_trafoP");


	return 0;


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



	// deal with the Dirichlet BC
	// can this be simplified when not doing the bnd / p / int separation ?
	// can be seen a 2-step procedure, non_adv and adv separately

	// identifiy the D-BC dofs --- can do this by checking the c_f_bnd dofs
	int no_dbc_in_loc = 0;
	int no_not_dbc_in_loc = 0;
	std::set<int> elem_loc_dbc;
	set<int> elem_not_loc_dbc;
	for ( int index_c_f_bnd = 0; index_c_f_bnd < c_f_bnd.rows(); index_c_f_bnd++ )
	{
		if (c_f_bnd(index_c_f_bnd,0) == c_f_bnd(index_c_f_bnd,1))
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

	cout << "no_dbc_in_loc " << no_dbc_in_loc << endl;
	cout << "elem_loc_dbc.size() " << elem_loc_dbc.size() << endl;

	cout << "no_not_dbc_in_loc " << no_not_dbc_in_loc << endl;
	cout << "elem_not_loc_dbc.size() " << elem_not_loc_dbc.size() << endl;

	// do a all = bnd/p/int POD


	Eigen::MatrixXd c_f_all( myCLNS_trafo.curr_f_bnd.size()+myCLNS_trafo.curr_f_p.size()+myCLNS_trafo.curr_f_int.size() , no_snaps );
	c_f_all.block(0,0,c_f_bnd.rows(),c_f_bnd.cols()) = c_f_bnd;
	c_f_all.block(c_f_bnd.rows(),0,c_f_p.rows(),c_f_p.cols()) = c_f_p;
	c_f_all.block(c_f_bnd.rows()+c_f_p.rows(),0,c_f_int.rows(),c_f_int.cols()) = c_f_int;

	Eigen::BDCSVD<Eigen::MatrixXd> svd_c_f_all(c_f_all, Eigen::ComputeThinU);
	cout << svd_c_f_all.singularValues() << endl << endl;
	Eigen::MatrixXd c_f_all_PODmodes = svd_c_f_all.matrixU();
	// limit thyself to the 99.99%



	// do the MM multiplication
	Array<OneD, CoupledLocalToGlobalC0ContMapSharedPtr> locToGloMap = myCLNS_trafo.m_locToGloMap;
        const Array<OneD,const int>& loctoglobndmap = locToGloMap[0]->GetLocalToGlobalBndMap();
        const Array<OneD,const NekDouble>& loctoglobndsign = locToGloMap[0]->GetLocalToGlobalBndSign();

	// original python code
/*	M_trafo_no_pp_incl_dbc = np.zeros([ num_elem*nsize_bndry , nBndDofs ])
	for curr_elem in range(0, num_elem):
		cnt = curr_elem*nsize_bndry_p1
		cnt_no_pp = curr_elem*nsize_bndry
		for i in range(0, nsize_bndry_p1):
			gid1 = int(LocGloBndMap[cnt + i])
			sign1 = int(LocGloBndSign[cnt + i])
			if (gid1 >= 0):
				if (i < nsize_bndry):
					M_trafo_no_pp_incl_dbc[ cnt_no_pp + i, gid1 ] = sign1
*/

	int nBndDofs = locToGloMap[0]->GetNumGlobalBndCoeffs();  // number of global bnd dofs
	Eigen::MatrixXd Mtrafo(myCLNS.RB_A.rows(), nBndDofs);
        int  nel  = m_fields[0]->GetNumElmts(); // number of spectral elements
	cout << "nel " << nel << endl;
	cout << "Mtrafo.rows() " << Mtrafo.rows() << endl;
	cout << "Mtrafo.cols() " << Mtrafo.cols() << endl;
	cout << "loctoglobndmap.num_elements() " << loctoglobndmap.num_elements() << endl;
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
	Eigen::MatrixXd MtM(myCLNS.RB_A.rows(), myCLNS.RB_A.rows());
	MtM = Mtrafo * Mtrafo.transpose();
	// assume have all the local sub-matrices in myCLNS
	// to break the bnd/p/int, could just form two matrices, "no_adv" and "adv", whereby "adv" has the trilinear projection
	// better initialise to zero
	Eigen::MatrixXd no_adv_matrix = Eigen::MatrixXd::Zero(myCLNS.RB_A.rows() + myCLNS.RB_Dbnd.rows() + myCLNS.RB_C.cols(), myCLNS.RB_A.cols() + myCLNS.RB_Dbnd.rows() + myCLNS.RB_B.cols() );
	cout << "no_adv_matrix.rows() " << no_adv_matrix.rows() << endl;
	cout << "no_adv_matrix.cols() " << no_adv_matrix.cols() << endl;
	// write Eigen::MatrixXd no_adv_matrix blockwise 
	no_adv_matrix.block(0, 0, myCLNS.RB_A.rows(), myCLNS.RB_A.cols()) = MtM * myCLNS.RB_A_no_adv;
	no_adv_matrix.block(0, myCLNS.RB_A.cols(), myCLNS.RB_Dbnd.cols(), myCLNS.RB_Dbnd.rows()) = -MtM * myCLNS.RB_Dbnd.transpose();
	no_adv_matrix.block(0, myCLNS.RB_A.cols() + myCLNS.RB_Dbnd.rows(), myCLNS.RB_B.rows(), myCLNS.RB_B.cols()) = MtM * myCLNS.RB_B_no_adv;
	no_adv_matrix.block(myCLNS.RB_A.rows(), 0, myCLNS.RB_Dbnd.rows(), myCLNS.RB_Dbnd.cols()) = -myCLNS.RB_Dbnd;
	no_adv_matrix.block(myCLNS.RB_A.rows(), myCLNS.RB_A.cols() + myCLNS.RB_Dbnd.rows(), myCLNS.RB_Dint.rows(), myCLNS.RB_Dint.cols()) = -myCLNS.RB_Dint;	
	no_adv_matrix.block(myCLNS.RB_A.rows() + myCLNS.RB_Dbnd.rows(), 0, myCLNS.RB_C.cols(), myCLNS.RB_C.rows()) = myCLNS.RB_C_no_adv.transpose();
	no_adv_matrix.block(myCLNS.RB_A.rows() + myCLNS.RB_Dbnd.rows(), myCLNS.RB_A.cols(), myCLNS.RB_Dint.cols(), myCLNS.RB_Dint.rows()) = -myCLNS.RB_Dint.transpose();
	no_adv_matrix.block(myCLNS.RB_A.rows() + myCLNS.RB_Dbnd.rows(), myCLNS.RB_A.cols() + myCLNS.RB_Dbnd.rows(), myCLNS.RB_D.rows(), myCLNS.RB_D.cols()) = myCLNS.RB_D_no_adv;


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

	Eigen::MatrixXd c_f_all_PODmodes_wo_dbc(c_f_all_PODmodes.rows() - no_dbc_in_loc, c_f_all_PODmodes.cols());
	Eigen::VectorXd f_bnd_dbc(no_dbc_in_loc);
	Eigen::VectorXd f_bnd_dbc_full_size(c_f_all_PODmodes.rows());
	int counter_all = 0;
	int counter_dbc = 0;
	for (int index=0; index < c_f_all_PODmodes.rows(); ++index)
	{
//		cout << index << endl;

		if (!elem_loc_dbc.count(index))
		{
			f_bnd_dbc_full_size(index) = 0;
			c_f_all_PODmodes_wo_dbc.row(counter_all) = c_f_all_PODmodes.row(index);
			counter_all++;
//			cout << " counter all " << counter_all << endl;
		}
		else
		{
			f_bnd_dbc_full_size(index) = c_f_all_PODmodes(index,0);
			f_bnd_dbc(counter_dbc) = c_f_all_PODmodes(index,0);
			counter_dbc++;
//			cout << " counter dbc " << counter_dbc << endl;
//			cout << " no_dbc_in_loc " << no_dbc_in_loc << endl;
		}

	}

	// compute add_to rhs
	Eigen::VectorXd add_to_rhs(c_f_all_PODmodes.rows()); // probably need this for adv and non-adv
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
	
	


	// deal with the advection projection
	// first build the full advection matrix and then multiply with D-BC, add to right-hand-side and then remove the dbc dof
	// also need to have nicely orthonormalized (iteratively) basis vectors...
	// need to have every RB basis vector in phys for that, i.e., after the svd misch-masch
	// probably can use a separate object for that
	
	// first a mapping bnd / int --> loc
        // Unpack solution from Bnd and F_int to v_coeffs 
	
	int curr_trafo_iter = 0;
	Eigen::VectorXd f_bnd = c_f_all_PODmodes.block(0, curr_trafo_iter, myCLNS_trafo.curr_f_bnd.size(), 1);
	Eigen::VectorXd f_int = c_f_all_PODmodes.block(myCLNS_trafo.curr_f_bnd.size()+myCLNS_trafo.curr_f_p.size(), curr_trafo_iter, myCLNS_trafo.curr_f_int.size(), 1);
	Array<OneD, MultiRegions::ExpListSharedPtr> fields = myCLNS.UpdateFields(); // I think this messes up the content of myCLNS local field, b.c. of pointer
        Array<OneD, unsigned int> bmap, imap; 
	Array<OneD, double> field_0(m_equ[0]->GetNcoeffs());
	Array<OneD, double> field_1(m_equ[0]->GetNcoeffs());
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


	// then map the loc field to the phys field m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(), m_fields[0]->UpdatePhys()); // needed if going with the equation system 
	Array<OneD, double> PhysBaseVec_x(myCLNS.GetNpoints(),0.0);
	Array<OneD, double> PhysBaseVec_y(myCLNS.GetNpoints(),0.0);
	Array<OneD, double> PhysBase_zero(myCLNS.GetNpoints(),0.0);
//	fields[0]->BwdTrans(fields[0]->GetCoeffs(), PhysBaseVec_x);
//	fields[1]->BwdTrans(fields[1]->GetCoeffs(), PhysBaseVec_y);		

	fields[0]->BwdTrans(field_0, PhysBaseVec_x);
	fields[0]->BwdTrans(field_1, PhysBaseVec_y);

	cout << "PhysBaseVec_x.num_elements() " << PhysBaseVec_x.num_elements() << endl;

	// orthonormalize the phys basis



	// do create the projected adv matrices
	int RBsize = 1; // initially have just one for test purposes
//	Eigen::MatrixXd c_f_bnd( myCLNS_trafo.curr_f_bnd.size() , no_snaps );
//	Eigen::MatrixXd c_f_p( myCLNS_trafo.curr_f_p.size() , no_snaps );
//	Eigen::MatrixXd c_f_int( myCLNS_trafo.curr_f_int.size() , no_snaps );
	for(int trafo_iter = 0; trafo_iter < RBsize; trafo_iter++)
	{
		myCLNS.InitObject();
		myCLNS.DoInitialiseAdv(PhysBaseVec_x, PhysBase_zero); // call with parameter in phys state
		myCLNS.DoSolve(); // do I need this here ?
		myCLNS.Output();  // do I need this here ?

		

	
	//	c_f_bnd.col(trafo_iter) = myCLNS_trafo.curr_f_bnd;
	//	c_f_p.col(trafo_iter) = myCLNS_trafo.curr_f_p;
	//	c_f_int.col(trafo_iter) = myCLNS_trafo.curr_f_int;

		// collect the advection matrices
		// separate into rhs and lhs part
		// save in projected form for each RB basis function

		Eigen::MatrixXd adv_matrix = Eigen::MatrixXd::Zero(myCLNS.RB_A.rows() + myCLNS.RB_Dbnd.rows() + myCLNS.RB_C.cols(), myCLNS.RB_A.cols() + myCLNS.RB_Dbnd.rows() + myCLNS.RB_B.cols() );
		cout << "adv_matrix.rows() " << adv_matrix.rows() << endl;
		cout << "adv_matrix.cols() " << adv_matrix.cols() << endl;
		// write Eigen::MatrixXd no_adv_matrix blockwise 
		adv_matrix.block(0, 0, myCLNS.RB_A.rows(), myCLNS.RB_A.cols()) = MtM * myCLNS.RB_A_adv;
		adv_matrix.block(0, myCLNS.RB_A.cols(), myCLNS.RB_Dbnd.cols(), myCLNS.RB_Dbnd.rows()) = -MtM * myCLNS.RB_Dbnd.transpose();
		adv_matrix.block(0, myCLNS.RB_A.cols() + myCLNS.RB_Dbnd.rows(), myCLNS.RB_B.rows(), myCLNS.RB_B.cols()) = MtM * myCLNS.RB_B_adv;
		adv_matrix.block(myCLNS.RB_A.rows(), 0, myCLNS.RB_Dbnd.rows(), myCLNS.RB_Dbnd.cols()) = -myCLNS.RB_Dbnd;
		adv_matrix.block(myCLNS.RB_A.rows(), myCLNS.RB_A.cols() + myCLNS.RB_Dbnd.rows(), myCLNS.RB_Dint.rows(), myCLNS.RB_Dint.cols()) = -myCLNS.RB_Dint;	
		adv_matrix.block(myCLNS.RB_A.rows() + myCLNS.RB_Dbnd.rows(), 0, myCLNS.RB_C.cols(), myCLNS.RB_C.rows()) = myCLNS.RB_C_adv.transpose();
		adv_matrix.block(myCLNS.RB_A.rows() + myCLNS.RB_Dbnd.rows(), myCLNS.RB_A.cols(), myCLNS.RB_Dint.cols(), myCLNS.RB_Dint.rows()) = -myCLNS.RB_Dint.transpose();
		adv_matrix.block(myCLNS.RB_A.rows() + myCLNS.RB_Dbnd.rows(), myCLNS.RB_A.cols() + myCLNS.RB_Dbnd.rows(), myCLNS.RB_D.rows(), myCLNS.RB_D.cols()) = myCLNS.RB_D_adv;
		
		// compute add_to rhs
		Eigen::VectorXd add_to_rhs_adv(c_f_all_PODmodes.rows()); // probably need this for adv and non-adv
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
		cout << "adv_mat_proj.rows() " << adv_mat_proj.rows() << endl;
		cout << "adv_mat_proj.cols() " << adv_mat_proj.cols() << endl;
		Eigen::VectorXd adv_rhs_proj = c_f_all_PODmodes_wo_dbc.transpose() * adv_rhs_add;



	}


	// then should be good for the online solve - in regard to the Q_a term, keep the all structure
	// so what do I need? 
	// since at this point, I have all the projected parameter-independent parts,
	// collect the parameter-dependent matrices 

	// steps to do in particolare
	// given a current iterate, look at the projection coefficients
	// form the sum of the advection term with the projection coefficients
	// \nu dependent are the non-advection, non-pressure terms, so need to differentiate those

	// have an online solve





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
