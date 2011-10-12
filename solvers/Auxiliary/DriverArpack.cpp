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
    std::string DriverArpack::arpackProblemTypeLookupIds[6] = {
            LibUtilities::SessionReader::RegisterEnumValue("ArpackProblemType","LargestReal"    ,0),
            LibUtilities::SessionReader::RegisterEnumValue("ArpackProblemType","SmallestReal"   ,1),
            LibUtilities::SessionReader::RegisterEnumValue("ArpackProblemType","LargestImag"    ,2),
            LibUtilities::SessionReader::RegisterEnumValue("ArpackProblemType","SmallestImag"   ,3),
            LibUtilities::SessionReader::RegisterEnumValue("ArpackProblemType","LargestMag"     ,4),
            LibUtilities::SessionReader::RegisterEnumValue("ArpackProblemType","SmallestMag"    ,5),
    };
    std::string DriverArpack::arpackProblemTypeDefault = LibUtilities::SessionReader::RegisterDefaultSolverInfo("ArpackProblemType","LargestMag");
    std::string DriverArpack::driverLookupId = LibUtilities::SessionReader::RegisterEnumValue("Driver","Arpack",0);

    std::string DriverArpack::className = GetDriverFactory().RegisterCreatorFunction("Arpack", DriverArpack::create);

    std::string DriverArpack::ArpackProblemTypeTrans[6] =
            { "LR", "SR", "LI", "SI", "LM", "SM" };


    /**
     *
     */
    DriverArpack::DriverArpack(const LibUtilities::SessionReaderSharedPtr pSession)
        : DriverArnoldi(pSession)
    {
    }
    

    /**
     *
     */
    DriverArpack::~DriverArpack()
    {
    }
    
    /**
     *
     */
    void DriverArpack::v_InitObject()
    {
        DriverArnoldi::v_InitObject();
        
        //Initialisation of Arnoldi parameters
        m_maxn   = 1000000; // Maximum size of the problem
        m_maxnev = 200;      // maximum number of eigenvalues requested
        m_maxncv = 500;     // Largest number of basis vector used in Implicitly Restarted Arnoldi		
	
        m_session->LoadParameter("realShift", m_realShift, 0.0);
        
        m_equ[0]->SetLambda(m_realShift);
                
        // Error alerts
        ASSERTL0(m_nvec <= m_maxnev,"NEV is greater than MAXNEV");
        ASSERTL0(m_kdim <= m_maxncv,"NEV is greater than MAXNEV");
        ASSERTL0(2      <= m_kdim-m_nvec,"NCV-NEV is less than 2");
	
        m_equ[0]->PrintSummary(cout);
        
        ArpackSummary(cout);
        
        m_equ[m_nequ - 1]->DoInitialise();
		
		//FwdTrans Initial conditions to be in Coefficient Space
		m_equ[m_nequ-1] ->TransPhysToCoeff();

    }
    
    void DriverArpack::ArpackSummary(std::ostream &out)
    {
        // Print session parameters
        out << "\tArnoldi solver type    : Arpack" << endl;

        out << "\tArpack problem type    : ";
        out << ArpackProblemTypeTrans[m_session->GetSolverInfoAsEnum<int>("ArpackProblemType")] << endl;

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
        out << "\tEvolution operator     : " << m_session->GetSolverInfo("EvolutionOperator") << endl;
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

        int  nq     = m_equ[0]->UpdateFields()[0]->GetNcoeffs(); // Number of points in the mesh
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


			if(m_session->DefinesFunction("InitialConditions"))
			{
				
				cout << "\tInital vector       : input file  " << endl;
				info = 1;
				CopyFieldToArnoldiArray(resid);
				
			}
			else
			{
				cout << "\tInital vector       : random  " << endl;
				info = 0;
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
		const char* problem = ArpackProblemTypeTrans[m_session->GetSolverInfoAsEnum<int>("ArpackProblemType")].c_str();

        FILE *pFile;
        std::string name = m_session->GetFilename().substr(0,m_session->GetFilename().find_last_of('.'))+".evl";
        pFile= fopen (name.c_str(), "w");
        
        ido     = 0;    //At the first call must be initialisedat 0

        while(ido != 99)//ido==-1 || ido==1 || ido==0)
        {
            //Routine for eigenvalue evaluation for non-symmetric operators
            Arpack::Dnaupd( ido, "I",       // B='I' for std eval problem
                            n, problem,  m_nvec,
                            m_evtol, &resid[0], m_kdim, 
                            &v[0], n, iparam, ipntr, &workd[0],
                            &workl[0], lworkl, info);
            
            //Plotting of real and imaginary part of the eigenvalues from workl
            cout << "\rIteration " << cycle << ", output: " << info << ", ido=" << ido << " " << std::flush <<endl;

            if(!((cycle-1)%m_kdim)&&(cycle> m_kdim))
            {
                cout << endl;
                fprintf (pFile, "Krylov spectrum at iteration: %i\t \n", cycle);
                for(int k=0; k<=m_kdim-1; ++k)
                {                    
                    // write m_nvec eigs to screen
                    if(m_kdim-1-k < m_nvec)
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

	    m_equ[0]->TransCoeffToPhys();

            m_equ[0]->DoSolve();

            if(!(cycle%m_infosteps))
            {
                cout << endl;
                m_equ[0]->Output();
            }

            if(m_EvolutionOperator == eTransientGrowth)
            {
				//start Adjoint with latest fields of direct 
				CopyFwdToAdj();
				
				m_equ[1]->TransCoeffToPhys();
				
				m_equ[1]->DoSolve();
            }
        
            // operated fields are copied into workd[inptr[1]-1] 
            CopyFieldToArnoldiArray(tmpworkd = workd + (ipntr[1]-1));
            
        }

        cout<< endl << "Converged in " << iparam[8] << " iterations" << endl;
	
        ASSERTL0(info >= 0," Error with Dnaupd");
	
        ritzSelect = Array<OneD, int>       (m_kdim,0);
        dr         = Array<OneD, NekDouble> (m_nvec+1,0.0);
        di         = Array<OneD, NekDouble> (m_nvec+1,0.0);
        workev     = Array<OneD, NekDouble> (3*m_kdim);
        z          = Array<OneD, NekDouble> (n*(m_nvec+1));
        
        sigmar     = m_realShift; 
        sigmai     = 0.0;
	
        //Setting 'A', Ritz vectors are computed. 'S' for Shur vectors
        Arpack::Dneupd(1, "A", ritzSelect.get(), dr.get(), di.get(), z.get(), n, sigmar, sigmai, workev.get(), "I", n, 
        	problem, m_nvec, m_evtol, resid.get(), m_kdim, v.get(), n, iparam, ipntr, workd.get(), workl.get(),lworkl,info);
		
        ASSERTL0(info == 0, " Error with Dneupd");
		int nconv=iparam[4];	
        Array<OneD, MultiRegions::ExpListSharedPtr>  fields = m_equ[0]->UpdateFields();
        
        cout << "Converged Eigenvalues: " << nconv << endl;
        fprintf(pFile,"Converged Eigenvalues: %d\n:",nconv);
		
		
        for(int i= 0; i< nconv; ++i)
        {
            WriteEvs(stdout,i,dr[i],di[i]);
            WriteEvs(pFile,i,dr[i],di[i]);
			
            std::string file = m_session->GetFilename().substr(0,m_session->GetFilename().find_last_of('.')) + "_eig_" + boost::lexical_cast<std::string>(i);
            
            WriteFld(file,z + k*nq);
        }

        fclose (pFile);
        
        if(m_EvolutionOperator != eTransientGrowth)
        {
            cout<<"Dump eigenvector: "<<nconv-1<<endl;
            CopyArnoldiArrayToField(z);               
            m_equ[0]->DoSolve(); 
            
            for (int k = 0; k < m_nfields; ++k)
            {
                Vmath::Vcopy(nq, &z[k*nq+(nconv-1)*n], 1, &fields[k]->UpdateCoeffs()[0] , 1);				
            }			
            
            for (int k = 0; k < m_nfields; ++k)
            {
                //Backward transformation in the physical space for plotting eigenmodes
                fields[k]->BwdTrans_IterPerExp(fields[k]->GetCoeffs(),fields[k]->UpdatePhys());
                fields[k]->SetPhysState(true);
            } 
            
            m_equ[0]->Output();
        }
        
	for(int j = 0; j < m_equ[0]->GetNvariables(); ++j)
        {
            NekDouble vL2Error = m_equ[0]->L2Error(j,false);
            NekDouble vLinfError = m_equ[0]->LinfError(j);
            if (m_comm->GetRank() == 0)
            {
                cout << "L 2 error (variable " << m_equ[0]->GetVariable(j) << ") : " << vL2Error << endl;
                cout << "L inf error (variable " << m_equ[0]->GetVariable(j) << ") : " << vLinfError << endl;
            }
			
        }
    }

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
