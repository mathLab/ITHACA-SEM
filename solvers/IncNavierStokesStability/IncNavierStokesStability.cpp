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
#include <iomanip>

#include <ADRSolver/EquationSystem.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/LinearAlgebra/Arpack.hpp>

using namespace Nektar;

/// Generates a new vector in the sequence by applying the linear operator.
void EV_update(
        EquationSystemSharedPtr &equ,
        Array<OneD, NekDouble> &src,
        Array<OneD, NekDouble> &tgt);

/// Generates the upper Hessenberg matrix H and computes its eigenvalues.
void EV_small(
        Array<OneD, Array<OneD, NekDouble> > &Kseq,
        const int ntot,
        const Array<OneD, NekDouble> &alpha,
        const int kdim,
        Array<OneD, NekDouble> &zvec,
        Array<OneD, NekDouble> &wr,
        Array<OneD, NekDouble> &wi,
        NekDouble &resnorm);

/// Tests for convergence of eigenvalues of H.
int EV_test(
        const int itrn,
        const int kdim,
        Array<OneD, NekDouble> &zvec,
        Array<OneD, NekDouble> &wr,
        Array<OneD, NekDouble> &wi,
        const NekDouble resnorm,
        const NekDouble evtol,
        const int nvec,
        ofstream &evlout);

/// Sorts a sequence of eigenvectors/eigenvalues by magnitude.
void EV_sort(
        Array<OneD, NekDouble> &evec,
        Array<OneD, NekDouble> &wr,
        Array<OneD, NekDouble> &wi,
        Array<OneD, NekDouble> &test,
        const int dim);


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
	bool useModifiedArnoldi;
	
	//Get the initial Arnoldi Vector
	std::string invecType = session->GetSolverInfo("InitialVector");
	
	//NekDouble info_ar= session->GetParameter("InitialVector");
    session->MatchSolverInfo("EigenvalueAlgorithm","ModifiedArnoldi",useModifiedArnoldi,true);
	
    // Print a summary of solver and problem parameters and initialise
    // the solver.
    equ->PrintSummary(cout);
    equ->DoInitialise();
    //initialise force if it necessary
    equ->SetInitialForce(0.0);
	
    // --- ORIGINAL ARNOLDI ALGORITHM (USING ARPACK) ---
    if (!useModifiedArnoldi)
    {
		int maxn=1000000; //Maximum size of the problem
        int maxnev=12;  //maximum number of eigenvalues requested
        int maxncv=30;  //Largest number of basis vector used in Implicitly Restarted Arnoldi

        int       nfields= equ->UpdateFields().num_elements()-1;
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
    }

    // --- MODIFIED ARNOLDI ALGORITHM ---
    else
    {
        // Declare variables
        Array<OneD, MultiRegions::ExpListSharedPtr>& fields
                                = equ->UpdateFields();
        int nfields             = fields.num_elements();
        int nq                  = fields[0]->GetNpoints();
        int ntot                = nfields*nq;
        int converged           = 0;
        NekDouble resnorm;
        std::string evlFile     = session->GetFilename() + ".evl";
        int kdim                = 0;
        int nvec                = 0;
        int nits                = 0;
        NekDouble evtol         = 1e-06;

        // Initialise progress output and load session parameters
        ofstream evlout(evlFile.c_str());
        session->LoadParameter("kdim",  kdim,  8);
        session->LoadParameter("nvec",  nvec,  1);
        session->LoadParameter("nits",  nits,  500);
        session->LoadParameter("evtol", evtol, 1e-06);

        // Allocate memory
        Array<OneD, NekDouble> alpha                (kdim + 1,      0.0);
        Array<OneD, NekDouble> wr                   (kdim,          0.0);
        Array<OneD, NekDouble> wi                   (kdim,          0.0);
        Array<OneD, NekDouble> zvec                 (kdim * kdim,   0.0);
        Array<OneD, Array<OneD, NekDouble> > Kseq   (kdim + 1);
        Array<OneD, Array<OneD, NekDouble> > Tseq   (kdim + 1);
        for (int i = 0; i < kdim + 1; ++i)
        {
            Kseq[i] = Array<OneD, NekDouble>(ntot, 0.0);
            Tseq[i] = Array<OneD, NekDouble>(ntot, 0.0);
        }

        // Print session parameters
        cout << "Krylov-space dimension: " << kdim << endl;
        cout << "Number of vectors:      " << nvec << endl;
        cout << "Max iterations:         " << nits << endl;
        cout << "Eigenvalue tolerance:   " << evtol << endl;

        // Copy starting vector into second sequence element (temporary).
        for (int k = 0; k < nfields; ++k)
        {
            Vmath::Vcopy(nq, &fields[k]->GetPhys()[0], 1, &Kseq[1][0] + k*nq, 1);
        }

        // Perform one iteration to enforce boundary conditions.
        // Set this as the initial value in the sequence.
        EV_update(equ, Kseq[1], Kseq[0]);

        // Normalise first vector in sequence
        alpha[0] = std::sqrt(Vmath::Dot(ntot, &Kseq[0][0], 1, &Kseq[0][0], 1));
        alpha[0] = std::sqrt(alpha[0]);
        Vmath::Smul(ntot, 1.0/alpha[0], Kseq[0], 1, Kseq[0], 1);

        // Fill initial krylov sequence
        for (int i = 1; !converged && i <= kdim; ++i)
        {
            // Compute next vector
            EV_update(equ, Kseq[i-1], Kseq[i]);

            // Normalise
            alpha[i] = std::sqrt(Vmath::Dot(ntot, &Kseq[i][0], 1, &Kseq[i][0], 1));
            alpha[i] = std::sqrt(alpha[i]);
            Vmath::Smul(ntot, 1.0/alpha[i], Kseq[i], 1, Kseq[i], 1);

            // Copy Krylov sequence into temporary storage
            for (int k = 0; k < i + 1; ++k)
            {
                Vmath::Vcopy(ntot, Kseq[k], 1, Tseq[k], 1);
            }

            // Generate Hessenberg matrix and compute eigenvalues of it.
            EV_small(Tseq, ntot, alpha, i, zvec, wr, wi, resnorm);

            // Test for convergence.
            converged = EV_test(i,i,zvec,wr,wi,resnorm,evtol,std::min(i,nvec),evlout);
        }

        // Continue with full sequence
        for (int i = kdim + 1; !converged && i <= nits; ++i)
        {
            // Shift all the vectors in the sequence.
            // First vector is removed.
            for (int j = 1; j <= kdim; ++j)
            {
                alpha[j-1] = alpha[j];
                Vmath::Vcopy(ntot, Kseq[j], 1, Kseq[j-1], 1);
            }

            // Compute next vector
            EV_update(equ, Kseq[kdim - 1], Kseq[kdim]);

            // Compute new scale factor
            alpha[kdim] = std::sqrt(Vmath::Dot(ntot, &Kseq[kdim][0], 1, &Kseq[kdim][0], 1));
            alpha[kdim] = std::sqrt(alpha[kdim]);
            Vmath::Smul(ntot, 1.0/alpha[kdim], Kseq[kdim], 1, Kseq[kdim], 1);

            // Copy Krylov sequence into temporary storage
            for (int k = 0; k < kdim + 1; ++k)
            {
                Vmath::Vcopy(ntot, Kseq[k], 1, Tseq[k], 1);
            }

            // Generate Hessenberg matrix and compute eigenvalues of it
            EV_small(Tseq, ntot, alpha, kdim, zvec, wr, wi, resnorm);

            // Test for convergence.
            converged = EV_test(i,kdim,zvec,wr,wi,resnorm,evtol,nvec,evlout);
        }

        // Close the runtime info file.
        evlout.close();
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
 *
 */
void EV_update(
        EquationSystemSharedPtr &equ,
        Array<OneD, NekDouble> &src,
        Array<OneD, NekDouble> &tgt)
{
    Array<OneD, MultiRegions::ExpListSharedPtr>& fields = equ->UpdateFields();
    int nfields = fields.num_elements();
    int nq = fields[0]->GetNpoints();

    // Copy starting vector into first sequence element.
    for (int k = 0; k < nfields; ++k)
    {
        Vmath::Vcopy(nq, &src[0] + k*nq, 1, &fields[k]->UpdatePhys()[0], 1);
    }
    equ->DoSolve();
    // Copy starting vector into first sequence element.
    for (int k = 0; k < nfields; ++k)
    {
        Vmath::Vcopy(nq, &fields[k]->GetPhys()[0], 1, &tgt[0] + k*nq, 1);
    }

}


/**
 *
 */
void EV_small(
        Array<OneD, Array<OneD, NekDouble> > &Kseq,
        const int ntot,
        const Array<OneD, NekDouble> &alpha,
        const int kdim,
        Array<OneD, NekDouble> &zvec,
        Array<OneD, NekDouble> &wr,
        Array<OneD, NekDouble> &wi,
        NekDouble &resnorm)
{
    int kdimp = kdim + 1;
    int lwork = 10*kdim;
    int ier;
    Array<OneD, NekDouble> R(kdimp * kdimp, 0.0);
    Array<OneD, NekDouble> H(kdimp * kdim, 0.0);
    Array<OneD, NekDouble> rwork(lwork, 0.0);

    // Modified G-S orthonormalisation
    for (int i = 0; i < kdimp; ++i)
    {
        NekDouble gsc = std::sqrt(Vmath::Dot(ntot, &Kseq[i][0], 1, &Kseq[i][0], 1));
        ASSERTL0(gsc != 0.0, "Vectors are linearly independent.");
        R[i*kdimp+i] = gsc;
        Vmath::Smul(ntot, 1.0/gsc, Kseq[i], 1, Kseq[i], 1);
        for (int j = i + 1; j < kdimp; ++j)
        {
            gsc = Vmath::Dot(ntot, Kseq[i], 1, Kseq[j], 1);
            Vmath::Svtvp(ntot, -gsc, Kseq[i], 1, Kseq[j], 1, Kseq[j], 1);
            R[j*kdimp+i] = gsc;
        }
    }

    // Compute H matrix
    for (int i = 0; i < kdim; ++i)
    {
        for (int j = 0; j < kdim; ++j)
        {
            H[j*kdim+i] = alpha[j+1] * R[(j+1)*kdimp+i]
                          - Vmath::Dot(j, &H[0] + i, kdim, &R[0] + j*kdimp, 1);
            H[j*kdim+i] /= R[j*kdimp+j];
        }
    }
    H[(kdim-1)*kdim+kdim] = alpha[kdim]
                * std::fabs(R[kdim*kdimp+kdim] / R[(kdim-1)*kdimp + kdim-1]);

    Lapack::dgeev_('N','V',kdim,&H[0],kdim,&wr[0],&wi[0],0,1,&zvec[0],kdim,&rwork[0],lwork,ier);

    ASSERTL0(!ier, "Error with dgeev");

    resnorm = H[(kdim-1)*kdim + kdim];
}


/**
 *
 */
int EV_test(
        const int itrn,
        const int kdim,
        Array<OneD, NekDouble> &zvec,
        Array<OneD, NekDouble> &wr,
        Array<OneD, NekDouble> &wi,
        const NekDouble resnorm,
        const NekDouble evtol,
        const int nvec,
        ofstream &evlout)
{
    int idone = 0;
    NekDouble re_ev, im_ev, abs_ev, ang_ev, re_Aev, im_Aev;
    NekDouble period = 0.1;
    Array<OneD, NekDouble> resid(kdim);
    for (int i = 0; i < kdim; ++i)
    {
        resid[i] = resnorm * std::fabs(zvec[kdim - 1 + i*kdim]) /
                std::sqrt(Vmath::Dot(kdim, &zvec[0] + i*kdim, 1, &zvec[0] + i*kdim, 1));
        if (wi[i] < 0.0) resid[i-1] = resid[i] = hypot(resid[i-1], resid[i]);
    }
    EV_sort(zvec, wr, wi, resid, kdim);

    if (resid[nvec-1] < evtol) idone = nvec;

    evlout << "-- Iteration = " << itrn << ", H(k+1, k) = " << resnorm << endl;
    evlout.precision(4);
    evlout.setf(ios::scientific, ios::floatfield);
    evlout << "EV  Magnitude   Angle       Growth      Frequency   Residual"
            << endl;
    for (int i = 0; i < kdim; i++) {
      re_ev  = wr[i];
      im_ev  = wi[i];
      abs_ev = hypot (re_ev, im_ev);
      ang_ev = atan2 (im_ev, re_ev);
      re_Aev = log (abs_ev) / period;
      im_Aev = ang_ev       / period;
      evlout << setw(2)  << i
           << setw(12) << abs_ev
           << setw(12) << ang_ev
       << setw(12) << re_Aev
       << setw(12) << im_Aev
       << setw(12) << resid[i]
       << endl;
    }

    return idone;
}


/**
 *
 */
void EV_sort(
        Array<OneD, NekDouble> &evec,
        Array<OneD, NekDouble> &wr,
        Array<OneD, NekDouble> &wi,
        Array<OneD, NekDouble> &test,
        const int dim)
{
    Array<OneD, NekDouble> z_tmp(dim,0.0);
    NekDouble wr_tmp, wi_tmp, te_tmp;
    for (int j = 1; j < dim; ++j)
    {
        wr_tmp = wr[j];
        wi_tmp = wi[j];
        te_tmp = test[j];
        Vmath::Vcopy(dim, &evec[0] + j*dim, 1, &z_tmp[0], 1);
        int i = j - 1;
        while (i >= 0 && test[i] > te_tmp)
        {
            wr[i+1] = wr[i];
            wi[i+1] = wi[i];
            test[i+1] = test[i];
            Vmath::Vcopy(dim, &evec[0] + i*dim, 1, &evec[0] + (i+1)*dim, 1);
            i--;
        }
        wr[i+1] = wr_tmp;
        wi[i+1] = wi_tmp;
        test[i+1] = te_tmp;
        Vmath::Vcopy(dim, &z_tmp[0], 1, &evec[0] + (i+1)*dim, 1);
    }
}

/**
 * $Log $
**/
