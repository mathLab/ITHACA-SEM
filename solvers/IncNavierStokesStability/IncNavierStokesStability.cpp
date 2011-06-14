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
// Description: Incompressible Navier Stokes Stability Solver
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <cmath> 

#include <ADRSolver/EquationSystem.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/LinearAlgebra/Arpack.hpp>

using namespace Nektar;

int main(int argc, char *argv[])
{
    
    if(argc != 2)
    {
        cerr << "\n \t Usage: IncNavierStokes  input.xml \n" << endl;
        exit(1);
    }

    string filename(argv[1]);
    time_t starttime, endtime;
    NekDouble CPUtime;
    string vCommModule("Serial");

    //----------------------------------------------------------------
    // Read the mesh and construct container class

    LibUtilities::SessionReaderSharedPtr session;
    LibUtilities::CommSharedPtr vComm;
    EquationSystemSharedPtr equ;

    // Record start time.
    time(&starttime);
    
    // Create session reader.
    session = MemoryManager<LibUtilities::SessionReader>::AllocateSharedPtr(filename);

    // Create communicator
    if (session->DefinesSolverInfo("Communication"))
    {
        vCommModule = session->GetSolverInfo("Communication");
    }
    else if (LibUtilities::GetCommFactory().ModuleExists("ParallelMPI"))
    {
        vCommModule = "ParallelMPI";
    }
    vComm = LibUtilities::GetCommFactory().CreateInstance(vCommModule, argc, argv);
    
    // Create instance of module to solve the equation specified in the session.
    try
    {
        equ = GetEquationSystemFactory().CreateInstance(session->GetSolverInfo("SOLVERTYPE"), vComm, session);
    }
    catch (int e)
    {
        ASSERTL0(e == -1, "No such solver class defined.");
    }


	//Evaluation of the time period
	NekDouble ts=session->GetParameter("TimeStep");
	NekDouble numstep=session->GetParameter("NumSteps");
    NekDouble period= ts*numstep;
	
	//Get the initial Arnoldi Vector
	std::string invecType = session->GetSolverInfo("InitialVector");
	
	//NekDouble info_ar= session->GetParameter("InitialVector");
    
	
    // Print a summary of solver and problem parameters and initialise
    // the solver.
    equ->PrintSummary(cout);
    equ->DoInitialise();
    //initialise force if it necessary
    equ->SetInitialForce(0.0);
	
	    
		int maxn=1000000; //Maximum size of the problem
        int maxnev=12;  //maximum number of eigenvalues requested
        int maxncv=30;  //Largest number of basis vector used in Implicitly Restarted Arnoldi

        int       nfields= equ->UpdateFields().num_elements();
        int       nq     = equ->UpdateFields()[0]->GetNpoints(); // Number of points in the mesh
        int       n      = nfields*nq; // Number of points in eigenvalue calculation
        NekDouble tol    = 1e-6; // determines the stopping criterion.
        int       ido    = 0;  //REVERSE COMMUNICATION parameter. At the first call must be initialised at 0
		int       info;        // do not set initial vector (info=0 random initial vector, info=1 read initial vector from session file)
        int       nev    = 2;  // Number of eigenvalues to be evaluated
        int       ncv    = 16; // Length of the Arnoldi factorisation
        int       lworkl = 3*ncv*(ncv+2); // Size of work array


			
        // Error alerts
        ASSERTL0(n   <= maxn,  "N is greater than MAXN");
        ASSERTL0(nev <= maxnev,"NEV is greater than MAXNEV");
        ASSERTL0(ncv <= maxncv,"NEV is greater than MAXNEV");
        ASSERTL0(2   <= ncv-nev,"NCV-NEV is less than 2");

        Array<OneD, NekDouble> resid (n);
        Array<OneD, NekDouble> v     (n*ncv);
        Array<OneD, NekDouble> workl (lworkl);
        Array<OneD, NekDouble> workd (3*n, 0.0);
        Array<OneD, MultiRegions::ExpListSharedPtr>& fields = equ->UpdateFields();

	

		if (session->DefinesFunction("InitialConditions"))
		{
			cout << " info=1: initial Arnoldi vector from session file" << endl;
			info=1;
			//Initialise resid to values specified in the initial conditions of the session file
			for (int k = 0; k < nfields; ++k)
			{
				Vmath::Vcopy(nq, &fields[k]->GetPhys()[0], 1, &resid[k*nq], 1);
			}
			 
		}
		else
		{
			cout << "info=0: random initial Arnoldi vector " << endl;
			info=0;
		}


        // Parameters
        int iparam[11];
        iparam[0] = 1;      // strategy for shift-invert
        iparam[1] = 0;      // (deprecated)
        iparam[2] = 400;    // maximum number of iterations allowed/taken
        iparam[3] = 1;      // blocksize to be used for recurrence
        iparam[4] = 0;      // number of converged ritz eigenvalues
        iparam[5] = 0;      // (deprecated)
        iparam[6] = 1;      // computation mode 1=> matrix-vector prod
        iparam[7] = 0;      // (for shift-invert)
        iparam[8] = 0;      // number of MV operations
        iparam[9] = 0;      // number of BV operations
        iparam[10]= 0;      // number of reorthogonalisation steps

        int ipntr[14];

        int cycle = 0;
        while(ido != 99)//ido==-1 || ido==1 || ido==0)
        {
            //Routine for eigenvalue evaluation for non-symmetric operators
            Arpack::Dnaupd( ido,        // reverse comm flag
                            "I",       // B='I' for std eval problem
                            n,          // problem size
                            "LM",      // type of problem (Largest Mag)
                            nev,        // number of eigenvalues to compute
                            tol,        // stopping tolerance
                            resid.get(),// array to store residuals
                            ncv,        // number of Lanczos vectors to use
                            v.get(),    // storage for lanczos vectors
                            n,
                            iparam,     // mode
                            ipntr,
                            workd.get(),
                            workl.get(),
                            lworkl,
                            info);

            cout << "Iteration " << cycle << ", output: " << info << ", ido=" << ido << endl;
            for(int k=0; k<=nev-1; ++k)
            {
             
				//Plotting of real and imaginary part of the eigenvalues from workl
                double r = workl[ipntr[5]-1+k];
                double i = workl[ipntr[6]-1+k];
                double res;
				

				cout << k << ": Mag " << sqrt(r*r+i*i) << ", angle " << atan2(i,r) << " growth " << log(sqrt(r*r+i*i))/period << 
				" Frequency " << atan2(i,r)/period << endl;
            }

            cycle++;

            if (ido == 99) break;

            ASSERTL0(ido == 1, "Unexpected reverse communication request.");

			
            //fields are copied in workd[ipntr[0]-1] and following
            for (int k = 0; k < nfields; ++k)
            {
                Vmath::Vcopy(nq, &workd[ipntr[0]-1+k*nq], 1, &fields[k]->UpdatePhys()[0] , 1);
            }

            equ->DoSolve();

            //Evoluted fields are copied into workd[inptr[1]-1] and following
            for (int k = 0; k < nfields; ++k)
            {
                Vmath::Vcopy(nq, &fields[k]->GetPhys()[0], 1, &workd[ipntr[1]-1+k*nq], 1);
            }
			
       }

        cout<< "Converged in " << iparam[8] << " iterations" << endl;

        ASSERTL0(info >= 0," Error with Dnaupd");

        Array<OneD, int> ritzSelect(ncv,0);
        Array<OneD, NekDouble> dr(nev+1,0.0);
        Array<OneD, NekDouble> di(nev+1,0.0);
        Array<OneD, NekDouble> workev(3*ncv);
        Array<OneD, NekDouble> z(n*(nev+1));
        double sigmar = 0.0, sigmai = 0.0;

        //Setting 'A', Ritz vectors are computed. 'S' for Shur vectors
        Arpack::Dneupd(1, "A", ritzSelect.get(), dr.get(), di.get(), z.get(), n,
                       sigmar, sigmai, workev.get(), "I", n, "LM", nev, tol,
                       resid.get(), ncv, v.get(), n, iparam, ipntr, workd.get(), workl.get(),
                       lworkl,info);

        ASSERTL0(info == 0, " Error with Dneupd");

        //number of "converged" Ritz values. This represents the number of Ritz values that satisfy
        //the convergence criterion.
        int nconv=iparam[4];

        cout << "===========================================" << endl;
        cout << "Stability Analysis " << endl;
        cout << "1)Size of the problem: " << n << " points" << endl;
        cout << "2)The number of Ritz values requested: " << nev << endl;
        cout << "3)The number of Arnoldi vectors generated: " << ncv << endl;
        cout << "4)What portion of the spectrum: " << "LM" << endl;
        cout << "5)The number of the converged Ritz value: " << nconv << endl;
        cout << "6)The number of OP*x: "<< iparam[8] <<endl;
        cout << "7)The number of Implicit Arnoldi update iterations taken: "<< iparam[2] <<endl;
        cout << "8)The convergence criterion: "<< tol << endl;

        //Printing of eigenvalues
        cout << "~~~~~~~~Eigenvalues of the problem~~~~~~~~" <<endl;

	for (int i = nev; i >= 0; --i)
	{
		cout << "Eigenvalue n. " << i+1 << " Re= "
		<< dr[i]<< "   Im=" << di[i]
		<< " Growth= " << log(sqrt(dr[i]*dr[i]+di[i]*di[i]))/period 
		<< " Frequency=" <<atan2(di[i],dr[i])/period <<endl;
		
		for (int k = 0; k < nfields; ++k)
		{
			Vmath::Vcopy(nq, &z[k*nq+i*n], 1, &fields[k]->UpdatePhys()[0] , 1);
			fields[k]->SetPhysState(true);
		}
		std::string file = session->GetFilename().substr(0,session->GetFilename().find_last_of('.'))
		+ ".eig." + boost::lexical_cast<std::string>(nev-i);
		equ->WriteFld(file);
	}

 
    // Record end time.
    time(&endtime);
    CPUtime = (1.0/60.0/60.0)*difftime(endtime,starttime);
 
    // Write output to .fld file
    equ->Output();
    
    // Evaluate and output computation time and solution accuracy.
    // The specific format of the error output is essential for the
    // regression tests to work.
    cout << "-------------------------------------------" << endl;
    cout << "Total Computation Time = " << CPUtime << " hr." << endl;
    for(int i = 0; i < equ->GetNvariables(); ++i)
    {
        // Get Exact solution
        Array<OneD, NekDouble> exactsoln(equ->GetTotPoints(),0.0);
        equ->EvaluateExactSolution(i,exactsoln,equ->GetFinalTime());
        
        cout << "L 2 error (variable " << equ->GetVariable(i)  << "): " << equ->L2Error(i,exactsoln) << endl;
        cout << "L inf error (variable " << equ->GetVariable(i)  << "): " << equ->LinfError(i,exactsoln) << endl;
    }
}

/**
 * $Log $
**/
