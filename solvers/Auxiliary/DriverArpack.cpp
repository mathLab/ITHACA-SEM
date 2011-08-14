///////////////////////////////////////////////////////////////////////////////
//
// File DriverArpack.cpp
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
// Description:  Arnoldi solver using Arpack
//
///////////////////////////////////////////////////////////////////////////////

#include <Auxiliary/DriverArpack.h>

namespace Nektar

{
    string DriverArpack::className = GetDriverFactory().RegisterCreatorFunction("Arpack", DriverArpack::create);
    
    
    /**
     *
     */
    DriverArpack::DriverArpack(LibUtilities::SessionReaderSharedPtr pSession)
        : DriverArnoldi(pSession)
    {
    }
    
    //Destructor
    DriverArpack::~DriverArpack()
    {
    }
    
    // Arpack problem type character string mappings
    int ArpackProbLen = 6;
    std::string ArpackProbTypes[] = 
        {
            "LargestReal",  "SmallestReal",
            "LargestImag",  "SmallestImag",
            "LargestMag",   "SmallestMag"
        };
    std::string ArpackProbTypesTrans[] = 
        { "LR", "SR",
          "LI", "SI",
          "LM", "SM"
        };
    
    void DriverArpack::v_InitObject()
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
            
		//	std:: string vAdvectionForm = m_session->GetSolverInfo("AdvectionForm");
             m_session->SetTag("AdvectiveType","Linearised");
                
            m_equ = Array<OneD, EquationSystemSharedPtr>(1);
            m_equ[0] = GetEquationSystemFactory().CreateInstance(vEquation, m_comm, m_session);
        }
        catch (int e)
            
        {
            ASSERTL0(e == -1, "No such class class defined.");
        }
	
	
        /// @todo This should be an independent Arnoldi call 
        m_session->MatchSolverInfo("SolverType","VelocityCorrectionScheme",m_TimeSteppingAlgorithm, false);
        
        //Initialisation of Arnoldi parameters
        m_maxn   = 1000000; // Maximum size of the problem
        m_maxnev = 12;      // maximum number of eigenvalues requested
        m_maxncv = 200;     // Largest number of basis vector used in Implicitly Restarted Arnoldi		
	
        if(m_TimeSteppingAlgorithm)
        {
            //Evaluation of the time period
            NekDouble ts      = m_session->GetParameter("TimeStep");
            NekDouble numstep = m_session->GetParameter("NumSteps");
            m_period          = ts*numstep;
            
            m_nfields = m_equ[0]->UpdateFields().num_elements()-1;
        }
        else
        {
            ASSERTL0(m_session->DefinesFunction("BodyForce"),"A BodyForce section needs to be defined for this solver type");
            m_nfields = m_equ[0]->UpdateFields().num_elements();
        }

        //Load values from session file if defined 
        
        // Length of the Arnoldi factorisation
        m_session->LoadParameter("kdim",  m_kdim, 16);
        // Number of eigenvalues to be evaluated 
        m_session->LoadParameter("nvec",  m_nvec,  2);
        // maximum number of iterations. 
        m_session->LoadParameter("nits",  m_nits,  500);
        // determines the stopping criterion.
        m_session->LoadParameter("evtol", m_evtol, 1e-6); 
        
        m_session->LoadParameter("realShift", m_realShift, 0.0);
        
        m_equ[0]->SetLambda(m_realShift);
                
        bool IsProbType;
        int i;
        for(i = 0; i < ArpackProbLen; ++i)
        {        
            m_session->MatchSolverInfo("ArpackProblemType",ArpackProbTypes[i].c_str(), IsProbType,false);
            if(IsProbType)
            {
                m_arpackProblemType = ArpackProbTypesTrans[i];                
                break;
            }
        }

        ASSERTL0(i  < ArpackProbLen,"Cannot determine the Arpack Problem Type defiend in ArpackProblemType")


        // Error alerts
        ASSERTL0(m_nvec <= m_maxnev,"NEV is greater than MAXNEV");
        ASSERTL0(m_kdim <= m_maxncv,"NEV is greater than MAXNEV");
        ASSERTL0(2      <= m_kdim-m_nvec,"NCV-NEV is less than 2");
	
        m_equ[0]->PrintSummary(cout);
        
        ArpackSummary(cout);
        
        m_equ[0]->DoInitialise();
    }
    
    void DriverArpack::ArpackSummary(std::ostream &out)
    {
        // Print session parameters
        out << "\tArnoldi solver type    : Arpack" << endl;

        out << "\tArpack problem type    : ";
        for(int i = 0; i < ArpackProbLen; ++i)
        {
            if(m_arpackProblemType == ArpackProbTypesTrans[i])
            {
                out << ArpackProbTypes[i] << endl;
            }
        }

        if(m_session->DefinesSolverInfo("SingleMode"))
        {
            out << "\tSingle Fourier mode    : true " << endl;
            ASSERTL0(m_session->DefinesSolverInfo("Homogeneous"),"Expected a homogeneous expansion to be defined with single mode");
        }
        else
        {
            out << "\tSingle Fourier mode    : false " << endl;
        }
        if(m_session->DefinesSolverInfo("BetaZero"))           
        {
            out << "\tBeta set to Zero       : true (overrides LHom)" << endl;
        }
        else
        {
            out << "\tBeta set to Zero       : false " << endl;
        }
        out << "\tReal Shift             : " << m_realShift << endl;
        out << "\tKrylov-space dimension : " << m_kdim << endl;
        out << "\tNumber of vectors      : " << m_nvec << endl;
        out << "\tMax iterations         : " << m_nits << endl;
        out << "\tEigenvalue tolerance   : " << m_evtol << endl;
	out << "=======================================================================" << endl;
    }

    void DriverArpack::v_Execute()
        
    {
        Array<OneD, NekDouble> tmpworkd;
        bool random;

        int  nq     = m_equ[0]->UpdateFields()[0]->GetNpoints(); // Number of points in the mesh
        int  n      = m_nfields*nq;    // Number of points in eigenvalue calculation
        int lworkl  = 3*m_kdim*(m_kdim+2); // Size of work array
        int       ido ;		//REVERSE COMMUNICATION parameter. At the first call must be initialised at 0
        int       info;     // do not set initial vector (info=0 random initial vector, info=1 read initial vector from session file)
		
        int iparam[11];
        int ipntr[14];
        
        Array<OneD, int> ritzSelect;
        Array<OneD, NekDouble> dr;
        Array<OneD, NekDouble> di;
        Array<OneD, NekDouble> workev;
        Array<OneD, NekDouble> z;
        NekDouble sigmar, sigmai;

        Array<OneD, NekDouble> resid(n);
        Array<OneD, NekDouble> v(n*m_kdim);
        Array<OneD, NekDouble> workl(lworkl);
        Array<OneD, NekDouble> workd(3*n, 0.0);

        ASSERTL0(n <= m_maxn,  "N is greater than   MAXN");

        m_session->MatchSolverInfo("InitialVector","Random",random,false);

        if(random)
        {
            cout << "\tInital vector       : random  " << endl;
            info = 0;
        }
        else
        {
            cout << "\tInital vector       : input file  " << endl;
            info = 1;
            CopyFieldToArnoldiArray(resid);
        }

	
        iparam[0] = 1;      // strategy for shift-invert
        iparam[1] = 0;      // (deprecated)
        iparam[2] = m_nits; // maximum number of iterations allowed/taken
        iparam[3] = 1;      // blocksize to be used for recurrence
        iparam[4] = 0;      // number of converged ritz eigenvalues
        iparam[5] = 0;      // (deprecated)
        if(fabs(m_realShift) > NekConstants::kNekZeroTol) // use shift if m_realShift > 1e-12
        {
            iparam[6] = 3;
        }
        else
        {
            iparam[6] = 1;      // computation mode 1=> matrix-vector prod
        }
        iparam[7] = 0;      // (for shift-invert)
        iparam[8] = 0;      // number of MV operations
        iparam[9] = 0;      // number of BV operations
        iparam[10]= 0;      // number of reorthogonalisation steps

        int cycle = 0;
		
        FILE *pFile;
        std::string name = m_session->GetFilename().substr(0,m_session->GetFilename().find_last_of('.'))+".evl";
        pFile= fopen (name.c_str(), "w");
        
        ido     = 0;    //At the first call must be initialisedat 0

        while(ido != 99)//ido==-1 || ido==1 || ido==0)
        {
            //Routine for eigenvalue evaluation for non-symmetric operators
            Arpack::Dnaupd( ido, "I",       // B='I' for std eval problem
                            n, m_arpackProblemType.c_str(),  m_nvec,
                            m_evtol, &resid[0], m_kdim, 
                            &v[0], n, iparam, ipntr, &workd[0],
                            &workl[0], lworkl, info);
            
            //Plotting of real and imaginary part of the eigenvalues from workl
            cout << "Iteration " << cycle << ", output: " << info << ", ido=" << ido << endl;
            if(cycle >= m_kdim)
            {
                fprintf (pFile, "Krylov spectrum at iteration: %i\t \n", cycle);
                for(int k=0; k<=m_kdim-1; ++k)
                {                    
                    // write m_nvec eigs to screen
                    if(k < m_nvec)
                    {
                        WriteEvs(stdout,k, workl[ipntr[5]-1+k],workl[ipntr[6]-1+k]);
                    }
                    // write m_kdim eigs to screen
                    WriteEvs(pFile,k, workl[ipntr[5]-1+k],workl[ipntr[6]-1+k]);
                }
            }
            
            cycle++;
            
            if (ido == 99) break;
                        
            ASSERTL0(ido == 1, "Unexpected reverse communication request.");

            //workd[inptr[0]-1] copied into operator fields
            CopyArnoldiArrayToField(tmpworkd = workd + (ipntr[0]-1));

            m_equ[0]->DoSolve();

            // operated fields are copied into workd[inptr[1]-1] 
            CopyFieldToArnoldiArray(tmpworkd = workd + (ipntr[1]-1));
            
        }
		
        cout<< "Converged in " << iparam[8] << " iterations" << endl;
	
        ASSERTL0(info >= 0," Error with Dnaupd");
	
        ritzSelect = Array<OneD, int>       (m_kdim,0);
        dr         = Array<OneD, NekDouble> (m_nvec+1,0.0);
        di         = Array<OneD, NekDouble> (m_nvec+1,0.0);
        workev     = Array<OneD, NekDouble> (3*m_kdim);
        z          = Array<OneD, NekDouble> (n*(m_nvec+1));
        
        sigmar     = m_realShift; 
        sigmai     = 0.0;
	
        //Setting 'A', Ritz vectors are computed. 'S' for Shur vectors
        Arpack::Dneupd(1, "A", ritzSelect.get(), dr.get(), di.get(), z.get(), n, sigmar, sigmai, workev.get(), "I", n, m_arpackProblemType.c_str(), m_nvec, m_evtol, resid.get(), m_kdim, v.get(), n, iparam, ipntr, workd.get(), workl.get(),lworkl,info);
		
        ASSERTL0(info == 0, " Error with Dneupd");
	       	
        int nconv=iparam[4];	
        Array<OneD, MultiRegions::ExpListSharedPtr>  fields = m_equ[0]->UpdateFields();
        
        cout << "Converged Eigenvalues: " << nconv << endl;
        fprintf(pFile,"Converged Eigenvalues: %d\n:",nconv);
        for(int i= 0; i< nconv; ++i)
        {
            WriteEvs(stdout,i,dr[i],di[i]);
            WriteEvs(pFile,i,dr[i],di[i]);
            
            for (int k = 0; k < m_nfields; ++k)
            {
                Vmath::Vcopy(nq, &z[k*nq+i*n], 1, &fields[k]->UpdatePhys()[0] , 1);
                fields[k]->SetPhysState(true);
            }
            
            std::string file = m_session->GetFilename().substr(0,m_session->GetFilename().find_last_of('.')) + "_eig_" + boost::lexical_cast<std::string>(i);

            m_equ[0]->WriteFld(file);
        }

        fclose (pFile);
    };
    
    void DriverArpack::WriteEvs(FILE *fp, const int k,  const NekDouble real, const NekDouble imag)
    {
        if(m_TimeSteppingAlgorithm)
        {
            fprintf (fp, "EV: %i\t , Mag: %10.6lf\t, angle:  %10.6lf\t, growth:  %10.6le\t, Frequency:  %10.6le \n",k, sqrt(real*real+imag*imag), atan2(imag,real),log(sqrt(real*real+imag*imag))/m_period, atan2(imag,real)/m_period );
        }
        else
        {
            NekDouble invmag = 1.0/(real*real + imag*imag);
            fprintf (fp, "EV: %i\t , Re: %10.6lf\t Imag:  %10.6lf\t inverse real:  %10.6le\t, inverse imag:  %10.6le\n",k, real, imag,-real*invmag, imag*invmag);
        }
    }

}
	
       

/**
 * $Log $
**/
