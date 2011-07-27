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

#include <LibUtilities/LinearAlgebra/Arpack.hpp>

#include <Auxiliary/DriverArnoldi.h>

namespace Nektar

{
    string DriverArnoldi::className = GetDriverFactory().RegisterCreatorFunction("Arnoldi", DriverArnoldi::create);
    
    /**
     *
     */
    DriverArnoldi::DriverArnoldi(LibUtilities::CommSharedPtr          pComm,
                                 LibUtilities::SessionReaderSharedPtr pSession)
        : Driver(pComm,pSession)
    {
    }
    
    //Destructor
    DriverArnoldi::~DriverArnoldi()
    {
    }
    
    void DriverArnoldi::v_InitObject()
    {
        try
        {
            ASSERTL0(m_session->DefinesSolverInfo("EqType"),
                     "EqType SolverInfo tag must be defined.");
            std::string vEquation = m_session->GetSolverInfo("EqType");
            if (m_session->DefinesSolverInfo("SolverType"))
            {
                vEquation = m_session->GetSolverInfo("SolverType");
            }
            ASSERTL0(GetEquationSystemFactory().ModuleExists(vEquation),
                     "EquationSystem '" + vEquation + "' is not defined.\n"
                     "Ensure equation name is correct and module is compiled.\n");
            
            m_equ = Array<OneD, EquationSystemSharedPtr>(1);
            m_equ[0] = GetEquationSystemFactory().CreateInstance(vEquation, m_comm, m_session);
        }
        catch (int e)
            
        {
            ASSERTL0(e == -1, "No such class class defined.");
        }
	
        //Evaluation of the time period
        NekDouble ts=m_session->GetParameter("TimeStep");
        NekDouble numstep=m_session->GetParameter("NumSteps");
        period= ts*numstep;
	
	
        //Initialisation of Arnoldi parameters
        maxn=1000000; //Maximum size of the problem
        maxnev=12;  //maximum number of eigenvalues requested
        maxncv=200;  //Largest number of basis vector used in Implicitly Restarted Arnoldi		
	
        nfields= m_equ[0]->UpdateFields().num_elements()-1;
        nq     = m_equ[0]->UpdateFields()[0]->GetNpoints(); // Number of points in the mesh
        n      = nfields*nq; // Number of points in eigenvalue calculation
        tol    = 1e-6; // determines the stopping criterion.
        ido    = 0;  //REVERSE COMMUNICATION parameter. At the first call must be initialised at 0
        nev    = 20;  // Number of eigenvalues to be evaluated
        ncv=     100;          // Length of the Arnoldi factorisation
        lworkl = 3*ncv*(ncv+2); // Size of work array
		
        //Load values from session file if defined 
        m_session->LoadParameter("kdim",  ncv, 16);
        m_session->LoadParameter("nvec",  nev, 2);
        m_session->LoadParameter("evtol", tol, 1e-6);
		
        // Error alerts
        ASSERTL0(n   <= maxn,  "N is greater than MAXN");
        ASSERTL0(nev <= maxnev,"NEV is greater than MAXNEV");
        ASSERTL0(ncv <= maxncv,"NEV is greater than MAXNEV");
        ASSERTL0(2   <= ncv-nev,"NCV-NEV is less than 2");
	
        m_equ[0]->PrintSummary(cout);
        m_equ[0]->DoInitialise();
	
        resid= Array<OneD, NekDouble> (n);
        v    = Array<OneD, NekDouble> (n*ncv);
        workl =Array<OneD, NekDouble> (lworkl);
        workd =Array<OneD, NekDouble> (3*n, 0.0);
        //CHECK AGAIN
        fields = m_equ[0]->UpdateFields();
	
        
    }
    
    void DriverArnoldi::v_Execute()
        
    {
        
        if (m_session->DefinesFunction("InitialConditions"))
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

        int cycle = 0;
		
        FILE *pFile;
        std::string name = m_session->GetFilename().substr(0,m_session->GetFilename().find_last_of('.'))+".evl";
        pFile= fopen (name.c_str(), "w");
	
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
            
            fprintf (pFile, "Iteration: %i\t ", cycle);
            fprintf (pFile, "\n");
            
            for(int k=0; k<=nev-1; ++k)
            {
                
                //Plotting of real and imaginary part of the eigenvalues from workl
                double r = workl[ipntr[5]-1+k];
                double i = workl[ipntr[6]-1+k];
                double res;
		
				
                cout << k << ": Mag " << sqrt(r*r+i*i) << ", angle " << atan2(i,r) << " growth " << log(sqrt(r*r+i*i))/period << 
                    " Frequency " << atan2(i,r)/period << endl;
		
                fprintf (pFile, "EV: %i\t , Mag: %f\t, angle:  %f\t, growth:  %f\t, Frequency:  %f\t ",k, sqrt(r*r+i*i), atan2(i,r),log(sqrt(r*r+i*i))/period, atan2(i,r)/period );
                fprintf (pFile, "\n");
		
		
            }
            
            cycle++;
            
            if (ido == 99) break;
            
            
            
            ASSERTL0(ido == 1, "Unexpected reverse communication request.");
            
			
            //fields are copied in workd[ipntr[0]-1] and following
            for (int k = 0; k < nfields; ++k)
            {
                Vmath::Vcopy(nq, &workd[ipntr[0]-1+k*nq], 1, &fields[k]->UpdatePhys()[0] , 1);
            }
			
            m_equ[0]->DoSolve();
			
            //Evoluted fields are copied into workd[inptr[1]-1] and following
            for (int k = 0; k < nfields; ++k)
            {
                Vmath::Vcopy(nq, &fields[k]->GetPhys()[0], 1, &workd[ipntr[1]-1+k*nq], 1);
            }
			
        }
		
		fclose (pFile);
        
		cout<< "Converged in " << iparam[8] << " iterations" << endl;
		
        ASSERTL0(info >= 0," Error with Dnaupd");
	    
        ritzSelect = Array<OneD, int> (ncv,0);
        dr         = Array<OneD, NekDouble>(nev+1,0.0);
        di         = Array<OneD, NekDouble>(nev+1,0.0);
        workev     = Array<OneD, NekDouble> (3*ncv);
        z          = Array<OneD, NekDouble> (n*(nev+1));
        sigmar = 0.0; 
        sigmai = 0.0;
	
        //Setting 'A', Ritz vectors are computed. 'S' for Shur vectors
        Arpack::Dneupd(1, "A", ritzSelect.get(), dr.get(), di.get(), z.get(), n, sigmar, sigmai, workev.get(), "I", n, "LM", nev, tol, resid.get(), ncv, v.get(), n, iparam, ipntr, workd.get(), workl.get(),lworkl,info);
		
        ASSERTL0(info == 0, " Error with Dneupd");
	
	
        bool VelCorrectionScheme; 
        m_session->MatchSolverInfo("SolverType","VelocityCorrectionScheme",VelCorrectionScheme, false);
	
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
        
        for(int i= 0; i< nconv; ++i)
        {
            if(VelCorrectionScheme==false)
            {
                cout << "Eigenvalue n. " << i+1 << " Re= "
                     << dr[i]<< "   Im=" << di[i]
                     << " Growth= " << log(sqrt(dr[i]*dr[i]+di[i]*di[i]))/period
                     << " Frequency=" <<atan2(di[i],dr[i])/period <<endl;
            }
            else
            {
                NekDouble invmag = 1.0/(dr[i]*dr[i]+di[i]*di[i]);
                cout << "Eigenvalue n. " << i+1 << " Re= "
                     << dr[i]<< "   Im=" << di[i]
                     << " Inverse real= " << -dr[i]*invmag
                     << " Inverse imag=" <<  di[i]*invmag << endl;
            }
            
            for (int k = 0; k < nfields; ++k)
            {
                Vmath::Vcopy(nq, &z[k*nq+i*n], 1, &fields[k]->UpdatePhys()[0] , 1);
                fields[k]->SetPhysState(true);
			}
            
            std::string file = m_session->GetFilename().substr(0,m_session->GetFilename().find_last_of('.')) + "_eig_" + boost::lexical_cast<std::string>(i);
            m_equ[0]->WriteFld(file);
        }
        
	
    };
}
	
	



/**
 * $Log $
**/
