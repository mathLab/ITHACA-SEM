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

	int snapshots_to_be_collected_aka_Nmax = 3;        

        // Create driver
        session->LoadSolverInfo("Driver", vDriverModule, "Standard");
	cout << "vDriverModule " << vDriverModule << endl;
	session->SetSolverInfo("SolverType", "CoupledLinearisedNS_trafoP");
        drv = GetDriverFactory().CreateInstance(vDriverModule, session);

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
        drv->Execute();
        // Finalise communications
        session->Finalise();

	Array<OneD, MultiRegions::ExpListSharedPtr> m_fields = m_equ[0]->UpdateFields();
	cout << "m_fields[0]->GetNcoeffs() " << m_fields[0]->GetNcoeffs() << endl; // this is nlc
	cout << "m_fields[0]->GetNpoints() " << m_fields[0]->GetNpoints() << endl;  // is this number of phys ??
	cout << "m_fields[0]->GetTotPoints() " << m_fields[0]->GetTotPoints() << endl;
	cout << "m_fields[0]->GetPhysState() " << m_fields[0]->GetPhysState() << endl;

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
	CoupledLinearNS_TT myCLNS_trafo(session);
	if (session->DefinesParameter("Kinvis"))
	{
		double testp = session->GetParameter("Kinvis");
		cout << "testp " << testp << endl;
	}
	myCLNS_trafo.InitObject();
	myCLNS_trafo.DoInitialiseAdv(collected_snapshots_x[0], collected_snapshots_y[0]); // replaces .DoInitialise();

	myCLNS_trafo.DoSolve();
	myCLNS_trafo.Output();

	// now can get the phys as well as bnd / p / int solution :)
	Array<OneD, MultiRegions::ExpListSharedPtr> m_fields_t = myCLNS_trafo.UpdateFields();

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

	Eigen::VectorXd trafo_f_bnd = myCLNS_trafo.curr_f_bnd;
	Eigen::VectorXd trafo_f_p = myCLNS_trafo.curr_f_p;
	Eigen::VectorXd trafo_f_int = myCLNS_trafo.curr_f_int;

	// deal with the Dirichlet BC
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
