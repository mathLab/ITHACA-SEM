/////////////////////////////////////////////////////////////////////////////
//
// File MMFEMFHN.cpp
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
//ㅡ_
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: MMFEMFHN solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <boost/algorithm/string.hpp>

#include <MMFSolver/EquationSystems/MMFEMFHN.h>

namespace Nektar
{
    Gs::string MMFEMFHN::className = SolverUtils::GetEquationSystemFactory().
        RegisterCreatorFunction("MMFEMFHN",
                                MMFEMFHN::create,
                                "MMFEMFHN equation.");

    MMFEMFHN::MMFEMFHN(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession),
          MMFSystem(pSession)
    {
    }

    void MMFEMFHN::v_InitObject()
    {
	// Variable assign
	m_E1  = 0;
	m_E2  = 1;
	m_Hz  = 2;
	m_phi = 3;
	m_psi = 4;
	m_chi = 5;
	m_Fs  = 6;
	m_Gr  = 7;
	m_dFdt = 8;
      
      UnsteadySystem::v_InitObject();

      int nq = m_fields[0]->GetNpoints();
      // int shapedim = m_fields[0]->GetShapeDimension();

      m_session->LoadParameter("ElemtGroup1", m_ElemtGroup1, 1);
      m_session->LoadParameter("StimulusPeriod", m_StimulusPeriod, -1.0);
      m_session->LoadParameter("radiusofinit", m_radiusofinit, 10.0);

      m_session->LoadParameter("InitPtx", m_InitPtx, 0.0);
      m_session->LoadParameter("InitPty", m_InitPty, 0.0);
      m_session->LoadParameter("InitPtz", m_InitPtz, 0.0);

      m_varepsilon = Array<OneD, NekDouble>(m_spacedim);
      m_session->LoadParameter("varepsilon1", m_varepsilon[0], 1.0);
      m_session->LoadParameter("varepsilon2", m_varepsilon[1], 1.0);
      m_session->LoadParameter("varepsilon3", m_varepsilon[2], 1.0);

      m_mu = Array<OneD, NekDouble>(m_spacedim);
      m_session->LoadParameter("mu1", m_mu[0], 1.0);
      m_session->LoadParameter("mu2", m_mu[1], 1.0);
      m_session->LoadParameter("mu3", m_mu[2], 1.0);

      m_sigma = Array<OneD, NekDouble>(m_spacedim);
      m_session->LoadParameter("sigma1", m_sigma[0], 1.0);
      m_session->LoadParameter("sigma2", m_sigma[1], 1.0);
      m_session->LoadParameter("sigma3", m_sigma[2], 1.0);

      m_session->LoadParameter("PoissonTau", m_PoissonTau, 1.0);
      m_session->LoadParameter("PoissonNeumann", m_PoissonNeumann, 1);
      m_session->LoadParameter("kp", m_kp, 10.0);
      /*
      for(int n = 0; n <  m_fields[m_phi]->GetBndConditions().num_elements(); ++n)
	{				
	  if ((m_fields[m_phi]->GetBndConditions())[n]->GetBoundaryConditionType() != SpatialDomains::eNeumann)
	    {
	      m_PoissonNeumann = 0;     
	    }
	}
      */
      
      Array<OneD, Array<OneD, NekDouble> > Anisotropy(m_spacedim);
      for(int j=0; j<m_spacedim; ++j)
	{
	  Anisotropy[j] = Array<OneD, NekDouble>(nq,1.0);
	}
      
      MMFSystem::MMFInitObject(Anisotropy);

      // Compute m_Vn
      int nTraceNumPoints = GetTraceNpoints();
	
      m_Vn = Array<OneD, Array<OneD, NekDouble> >(m_shapedim);
      
      Array<OneD, NekDouble> Vc(nq);
      Array<OneD, NekDouble> VcFwd(nTraceNumPoints);
      for (int j=0; j < m_shapedim; ++j)
	{
	  m_Vn[j] = Array<OneD, NekDouble> (nTraceNumPoints,0.0);
	  for(int i = 0; i < m_spacedim; ++i)
	    {		
	      Vmath::Vcopy(nq, &m_movingframes[j][i*nq], 1, &Vc[0], 1);
	      
	      m_fields[0]->ExtractTracePhys(Vc, VcFwd);
	      Vmath::Vvtvp(nTraceNumPoints,&m_traceNormals[i][0],1,&VcFwd[0],1,&m_Vn[j][0],1,&m_Vn[j][0],1);
	    }
	}
      
      // Sigma Block
      m_SigmaBlock = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
      for (int j=0; j<m_spacedim; j++)
	{
	  m_SigmaBlock[j] = Array<OneD, NekDouble>(nq);
	  m_fields[0]->GenerateElementVector(m_ElemtGroup1, 1.0, m_sigma[j], m_SigmaBlock[j]);
	}

      Gs::cout << "*** SigmaBlock = ( " << RootMeanSquare(m_SigmaBlock[0]) << " , " << RootMeanSquare(m_SigmaBlock[1]) << " ) " << Gs::endl;

      // Define TestType
      if(m_session->DefinesSolverInfo("TESTTYPE"))
	{
	  std::string TestTypeStr = m_session->GetSolverInfo("TESTTYPE");
	  for(int i = 0; i < (int) SIZE_TestType; ++i)
            {
	      if(NoCaseStringCompare(TestTypeMap[i],TestTypeStr) == 0)
                {
                    m_TestType = (TestType)i;
                    break;
                }
            }
	  }
        else
	  {
            m_TestType = (TestType)0;
	  }
	
	// Initwave type
        if(m_session->DefinesSolverInfo("INITWAVETYPE"))
	  {
            std::string InitWaveTypeStr = m_session->GetSolverInfo("INITWAVETYPE");
            for(int i = 0; i < (int) SIZE_TestType; ++i)
            {
                if(NoCaseStringCompare(InitWaveTypeMap[i],InitWaveTypeStr) == 0)
                {
                    m_InitWaveType = (InitWaveType)i;
                    break;
                }
            }
	  }
        else
	  {
            m_InitWaveType = (InitWaveType)0;
	  }

	// Define FentonKarmaType
        if(m_session->DefinesSolverInfo("FENTONKARMATYPE"))
	  {
            std::string FentonKarmaTypeStr = m_session->GetSolverInfo("FENTONKARMATYPE");
            for(int i = 0; i < (int) SIZE_FentonKarmaType; ++i)
            {
                if(NoCaseStringCompare(FentonKarmaTypeMap[i],FentonKarmaTypeStr) == 0)
                {
                    m_FentonKarmaType = (FentonKarmaType)i;
                    break;
                }
            }
	  }
	
        else
	  {
            m_FentonKarmaType = (FentonKarmaType)0;
	  }

	// Define TestMaxwellType
        if(m_session->DefinesSolverInfo("TESTMAXWELLTYPE"))
	  {
            std::string TestMaxwellTypeStr = m_session->GetSolverInfo("TESTMAXWELLTYPE");
            for(int i = 0; i < (int) SolverUtils::SIZE_TestMaxwellType; ++i)
            {
                if(NoCaseStringCompare(SolverUtils::TestMaxwellTypeMap[i],TestMaxwellTypeStr) == 0)
                {
                    m_TestMaxwellType = (SolverUtils::TestMaxwellType)i;
                    break;
                }
            }
	  }
	
        else
	  {
            m_TestMaxwellType = (SolverUtils::TestMaxwellType)0;
        }
	
	// Define Polarization
        if(m_session->DefinesSolverInfo("POLTYPE"))
	  {
            std::string PolTypeStr = m_session->GetSolverInfo("POLTYPE");
            for(int i = 0; i < (int) SolverUtils::SIZE_PolType; ++i)
            {
	      if(NoCaseStringCompare(SolverUtils::PolTypeMap[i],PolTypeStr) == 0)
                {
                    m_PolType = (SolverUtils::PolType)i;
                    break;
                }
            }
        }
        else
        {
	  m_PolType = (SolverUtils::PolType)0;
        }
	
	Array<OneD, NekDouble> One(nq, 1.0);	
	m_Jac = m_fields[0]->PhysIntegral(One);
	
	// Compute varcoeff for Diffusion
	ComputeVarcoeffDiffusion();
	
	// Compute the cross producted MF
	DeriveCrossProductMF(m_CrossProductMF);
	
	// Compute n_timesMFFwd and m_times_timesMFFwd
	ComputeNtimesMF();

	// Compute vaepsilon and mu vector (m_epsveci, m_muvec0);
	m_epsvec = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
	m_muvec = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
	for (int k=0; k<m_spacedim; ++k)
	  {
	    m_epsvec[k] = Array<OneD, NekDouble>(nq,1.0);
	    m_muvec[k] = Array<OneD, NekDouble>(nq,1.0);
	  }
	
	// MMEMFHN Special: Consider m_epsvec = m_muvec = 1
	ComputeZimYim(m_epsvec, m_muvec);

        // If explicit it computes RHS and PROJECTION for the time integration
	if (!m_explicitDiffusion)
        {
	  m_ode.DefineImplicitSolve (&MMFEMFHN::DoImplicitSolve, this);
        }

	m_ode.DefineOdeRhs(&MMFEMFHN::DoOdeRhs, this);
	m_ode.DefineProjection (&MMFEMFHN::DoOdeProjection, this);
    }

    MMFEMFHN::~MMFEMFHN()
    {
    }

  
  void MMFEMFHN::v_DoSolve()
  {
    ASSERTL0(m_intScheme != 0, "No time integration scheme.");
    
    int i, nchk = 1;
    int nq         = GetTotPoints();
    int nvariables;
    int nfields = m_fields.num_elements();
    
    if (m_intVariables.empty())
      {
	for (i = 0; i < nfields; ++i)
	  {
	    m_intVariables.push_back(i);
	  }
	nvariables = nfields;
      }
    else
      {
	nvariables = m_intVariables.size();
      }

    // Set up wrapper to fields data storage.
    Array<OneD, Array<OneD, NekDouble> > fields(nvariables);
    Array<OneD, Array<OneD, NekDouble> > tmp   (nvariables);

    // Order storage to list time-integrated fields first.
    for(i = 0; i < nvariables; ++i)
      {
	fields[i] = m_fields[m_intVariables[i]]->GetPhys();
	m_fields[m_intVariables[i]]->SetPhysState(false);
      }

    // Initialise time integration scheme
    m_intSoln = m_intScheme->InitializeScheme(m_timestep, fields, m_time, m_ode);
    
    // Check uniqueness of checkpoint output
    ASSERTL0((m_checktime == 0.0 && m_checksteps == 0) ||
	     (m_checktime >  0.0 && m_checksteps == 0) || 
	     (m_checktime == 0.0 && m_checksteps >  0),
	     "Only one of IO_CheckTime and IO_CheckSteps "
	     "should be set!");
    
    Timer     timer;
    bool      doCheckTime   = false;
    int       step          = 0;
    NekDouble intTime       = 0.0;
    NekDouble cpuTime       = 0.0;
    NekDouble elapsed       = 0.0;

    int      NoStim = 1 ;
    
    while (step   < m_steps || m_time < m_fintime - NekConstants::kNekZeroTol)
      {
	timer.Start();
	fields = m_intScheme->TimeIntegrate(step, m_timestep, m_intSoln, m_ode);
	timer.Stop();
	
	m_time  += m_timestep;
	elapsed  = timer.TimePerTest(1);
	intTime += elapsed;
	cpuTime += elapsed;

	if(m_StimulusPeriod>0)
	  {
	    if(!((step+1) % m_StimulusPeriod))
	      {
		Array<OneD, NekDouble> Initpulse(nq,0.0);
		PlaneWaveForPhi(Initpulse);

		Gs::cout << "*********************************" << Gs::endl;
		Gs::cout << " Stimulus is initiated at time = " << m_time << Gs::endl;
		Gs::cout << "*********************************" << Gs::endl;
		
		Vmath::Vadd(nq, &Initpulse[0], 1, &fields[m_phi][0], 1, &fields[m_phi][0], 1);

		/*
		int cnt=0;
		if(m_TestType==eFentonKarma)
		  {
		    for (int i=0; i<nq; i++)
		      {
     			if(fabs(Initpulse[i])>0.001)
			  {
			    // Initpsi[i] = 1.0;
			    // fields[m_psi][i] = 1.0;
			    // cnt++;
			  }
		      }
		    
		    // Vmath::Vadd(nq, &Initpsi[0], 1, &fields[m_psi][0], 1, &fields[m_psi][0], 1);
		  }
		*/
		NoStim++;
	      }
	  }

	// Input for m_Fs
	Vmath::Vcopy(nq, m_fields[m_Fs]->GetPhys(), 1, fields[m_Fs], 1);
	Vmath::Vcopy(nq, m_fields[m_Gr]->GetPhys(), 1, fields[m_Gr], 1);
	Vmath::Vcopy(nq, m_fields[m_dFdt]->GetPhys(), 1, fields[m_dFdt], 1);
	
	// Write out status information
	if (m_session->GetComm()->GetRank() == 0 && !((step+1) % m_infosteps))
	  {
	    Gs::cout << "Steps: " << std::setw(8)  << Gs::left << step+1 << " "
		 << "Time: "  << std::setw(12) << Gs::left << m_time;
	    
	    Gs::stringstream ss;
	    ss << cpuTime/60.0 << " min.";
	    Gs::cout << " CPU Time: " << std::setw(8) << Gs::left
		 << ss.str() << Gs::endl;

	    Array<OneD, Array<OneD, NekDouble> >Energy;
	    // ComputeEnergy(fields, Energy);
	    
	    Array<OneD, NekDouble> APDlength(10,0.0);
	    // ComputeAPD(fields[m_phi], APDlength);
	    
	    Gs::cout << "APD0 = " << APDlength[0] << ", APD1 = " << APDlength[1] << ", APD2 = " << APDlength[2] << Gs::endl;
	    Gs::cout << "phi = " << AbsIntegral(fields[m_phi]) << ", Fs = " << AbsIntegral(fields[m_Fs])
		 << ", dFdt = " << AbsIntegral(fields[m_dFdt]) << ", Gravity = " << AbsIntegral(fields[m_Gr]) << Gs::endl;

	    cpuTime = 0.0;
	  }
       
	// Transform data into coefficient space
	for (i = 0; i < nvariables; ++i)
	  {
	    m_fields[i]->SetPhys(fields[i]);
	    m_fields[i]->FwdTrans_IterPerExp(fields[i],m_fields[i]->UpdateCoeffs());
	    m_fields[i]->SetPhysState(false);
	  }

	// Write out checkpoint files
	if ((m_checksteps && step && !((step + 1) % m_checksteps)) || doCheckTime)
	  {
	    Checkpoint_Output(nchk);
	    // Checkpoint_XYZOutput(nchk,fields);
	    // Checkpoint_EnergyOutput(nchk,fields);
	    nchk++;
	    
	    /*
	    if(nvariables>7)
	    {
		Checkpoint_XYZOutput(nchk,fields);
		Checkpoint_EnergyOutput(nchk++,fields);
	    }

	    else
	      {
		nchk++;
	      }
	    */
	    doCheckTime = false;
	  }
	
	// Step advance
	++step;
      }
    
    // Print out summary statistics
    if (m_session->GetComm()->GetRank() == 0)
      {
	Gs::cout << "Time-integration  : " << intTime/60.0  << " min"   << Gs::endl;
      }
    
    for(i = 0; i < nvariables; ++i)
      {
	m_fields[m_intVariables[i]]->SetPhys(fields[i]);
	m_fields[m_intVariables[i]]->SetPhysState(true);
      }

    for(i = 0; i < nvariables; ++i)
      {
	m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
      }
  } 

  
  void MMFEMFHN::DoOdeRhs(const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
			      Array<OneD,        Array<OneD, NekDouble> >&outarray,
			      const NekDouble time)
    {
    int i, j;
    int nvar    = inarray.num_elements();
    int ncoeffs    = GetNcoeffs();
    int nq         = GetTotPoints();
    
    Array<OneD, Array<OneD, NekDouble> > physarray(nvar);
    Array<OneD, Array<OneD, NekDouble> > modarray(nvar);
    for (i = 0; i < nvar; ++i)
      {
	physarray[i] = Array<OneD, NekDouble>(nq);
	modarray[i]  = Array<OneD, NekDouble>(ncoeffs,0.0);
	
	Vmath::Vcopy(nq, &inarray[i][0], 1, &physarray[i][0], 1);
      }
    
    for (i=0; i<nvar; i++)
      {
	outarray[i] = Array<OneD, NekDouble>(nq,0.0);
	m_fields[i]->SetPhysState(true);
      }

    // modarray[m_E1var] = d H^z / d \xi^1 
    // modarray[m_E2var] = - d H^z / d \xi^2 
    // modarray[m_H3var] = d E^1 / d \xi^2 - d E^2 / d \xi^1
    WeakDGMaxwellDirDeriv(physarray, modarray, time);
    m_fields[m_E1]->SetPhysState(false);
    m_fields[m_E2]->SetPhysState(false);
    m_fields[m_Hz]->SetPhysState(false);

    /*
    if(nvar>7)
      {
	Vmath::Vcopy(ncoeffs, modarray[m_E1], 1, modarray[m_D1], 1);
	Vmath::Vcopy(ncoeffs, modarray[m_E2], 1, modarray[m_D2], 1);
	m_fields[m_D1]->SetPhysState(false);
	m_fields[m_D2]->SetPhysState(false);
      }
    */
	
    // To obtain \nabla phi
    // Compute the reactions for \phi and \psi
    // modarray[m_phi] = F(\phi,\psi)
    // modarray[m_psi] = G(\phi, \psi)

    Array<OneD, NekDouble> phireaction(nq,0.0);
    Array<OneD, NekDouble> psireaction(nq,0.0);
    Array<OneD, NekDouble> chireaction(nq,0.0);
    
    // Implement Boundary Conditinos
    SetBoundaryConditions(m_phi,time);
    SetBoundaryConditions(m_psi,time);
    SetBoundaryConditions(m_chi,time);

    if(m_TestType==eTestEMFHN)
      {
	phireaction = FlowTestEMFHN(time, 100);
      }

    else
      {
	ComputeDRReaction(physarray[m_phi], physarray[m_psi], physarray[m_chi], phireaction, psireaction, chireaction);
      }
	
    m_fields[m_phi]->IProductWRTBase(phireaction,modarray[m_phi]);
    m_fields[m_psi]->IProductWRTBase(psireaction,modarray[m_psi]);
    m_fields[m_chi]->IProductWRTBase(chireaction,modarray[m_chi]);
    m_fields[m_phi]->SetPhysState(false);
    m_fields[m_psi]->SetPhysState(false);
    m_fields[m_chi]->SetPhysState(false);

    
    // outarray[m_phi] = nabla \phi
    // dphidt = nabla \phi + F(phi,psi)
    Array<OneD, NekDouble> dphidt(nq);
    Vmath::Vadd(nq, outarray[m_phi], 1, phireaction, 1, dphidt, 1);

    // Solve Esource such that \nabla \cdot \sigma \nabla Esource = dFdt such that \int dFdt = 0
    Array<OneD, NekDouble> dFdt(nq);
    // Compute \dot{F} = \partial F/ \partial t
    if(m_TestType==eTestEMFHN)
      {
	dFdt = FlowTestEMFHN(time,40);
      }
    
    else
      {
	// ComputedFdt(physarray[m_phi], physarray[m_psi], dphidt, psireaction, dFdt);
	ComputedFdt(physarray[m_phi], physarray[m_psi], physarray[m_chi], dphidt, psireaction, chireaction, dFdt);
      }
    
    // Solve \nabla \phi_s \nabla = dFdt
    Vmath::Vcopy(nq, dFdt, 1, m_fields[m_dFdt]->UpdatePhys(), 1);
    SolvePoisson(dFdt);
    
    // dD/dt = \nabla \times H - \hat{\sigma} \nabla \phi
    // dE/dt = \nabla \times H - \nabla phi_s
    Array<OneD, NekDouble> tmpc(ncoeffs);
    Array<OneD, NekDouble> tmp(nq);

    Array<OneD, NekDouble> Grmag(nq,0.0);
    for (j=0; j<m_shapedim; j++)
      {
	/*
	if(nvar>7)
	  {
	    WeakDGAdv(m_D1+j, j, m_fields[m_phi]->GetPhys(), tmpc);
	    
	    // Multiplying conductivity tensor, \hat{sigma}
	    m_fields[m_D1+j]->BwdTrans(tmpc,tmp);
	    Vmath::Vmul(nq, m_SigmaBlock[j], 1, tmp, 1, tmp, 1);
	    m_fields[m_D1+j]->FwdTrans_IterPerExp(tmp, tmpc);
	    
	     Vmath::Vadd(ncoeffs, tmpc, 1, modarray[m_D1+j], 1, modarray[m_D1+j], 1);
	  }
	*/
	
	WeakDGAdv(j, j, m_fields[m_Fs]->GetPhys(), tmpc);
	Vmath::Neg(ncoeffs, tmpc, 1);

	m_fields[0]->BwdTrans(tmpc,tmp);

	Vmath::Vvtvp(nq, tmp, 1, tmp, 1, Grmag, 1, Grmag, 1);

	// Vmath::Vcopy(ncoeffs, tmpc, 1, m_fields[m_Gr1+j]->UpdateCoeffs(), 1);
	// m_fields[m_Gr1+j]->BwdTrans(m_fields[m_Gr1+j]->GetCoeffs(),m_fields[m_Gr1+j]->UpdatePhys());
				  
	Vmath::Vadd(ncoeffs, tmpc, 1, modarray[m_E1+j], 1, modarray[m_E1+j], 1);
      }

    Vmath::Vsqrt(nq, Grmag, 1, m_fields[m_Gr]->UpdatePhys(), 1);
    
    // Compensation for Covariation differentiation 
    AddGreenDerivCompensate(physarray,modarray);
     
    for(i = 0; i < nvar; ++i)
      {
	m_fields[i]->MultiplyByElmtInvMass(modarray[i],modarray[i]);
	m_fields[i]->BwdTrans(modarray[i],outarray[i]);
      }
    }


  void MMFEMFHN::WeakDGMaxwellDirDeriv(const Array<OneD, const Array<OneD, NekDouble> >& InField,
				       Array<OneD, Array<OneD, NekDouble> >& OutField,
				       const NekDouble time)
    {
        int i;
        int nq              = GetNpoints();
        int ncoeffs         = GetNcoeffs();
        int nTracePointsTot = GetTraceNpoints();
	int nvar = 3;

        Array<OneD, Array<OneD, NekDouble> > fluxvector(m_shapedim);
        for(i = 0; i < m_shapedim; ++i)
        {
            fluxvector[i]    = Array<OneD, NekDouble>(nq);
        }

        Array<OneD, Array<OneD, NekDouble> > physfield (nvar);
	for(i = 0; i < nvar; ++i)
	  {
	    physfield[i] = InField[i];
	  }
 
	Array<OneD, NekDouble> tmpc(ncoeffs);
	for(i = 0; i < nvar; ++i)
	  {	    
	    GetMaxwellFluxVector(i, physfield, fluxvector);
	    
	    OutField[i] = Array<OneD, NekDouble>(ncoeffs,0.0);
	    for(int j=0; j<m_shapedim;++j)
	      {		    
		// Directional derivation with respect to the j'th moving frame
		// tmp_j = ( \nabla \phi, fluxvector[j] \mathbf{e}^j )
		m_fields[i]->IProductWRTDirectionalDerivBase(m_CrossProductMF[j], fluxvector[j], tmpc);
		Vmath::Vadd(ncoeffs, &tmpc[0], 1, &OutField[i][0], 1, &OutField[i][0], 1);
	      }
	  }

	
        // V the numerical flux and add to the modal coeffs
        // if the NumericalFlux function does not include the
        // normal in the output
	Array<OneD, Array<OneD, NekDouble> > numfluxFwd   (nvar);
	Array<OneD, Array<OneD, NekDouble> > numfluxBwd   (nvar);
	
	for(i = 0; i < nvar; ++i)
	  {
	    numfluxFwd[i]   = Array<OneD, NekDouble>(nTracePointsTot,0.0);
	    numfluxBwd[i]   = Array<OneD, NekDouble>(nTracePointsTot,0.0);
	  }
	
	// Evaluate numerical flux in physical space which may in
	// general couple all component of vectors
	NumericalMaxwellFlux(physfield, numfluxFwd, numfluxBwd, time);
	
	// Evaulate  <\phi, \hat{F}\cdot n> - OutField[i]
	for(i = 0; i < nvar; ++i)
	  {
	    Vmath::Neg(ncoeffs,OutField[i],1);
	    m_fields[i]->AddFwdBwdTraceIntegral(numfluxFwd[i], numfluxBwd[i], OutField[i]);
	    m_fields[i]->SetPhysState(false);
	  }
    }


  void MMFEMFHN::DoImplicitSolve(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
				 Array<OneD, Array<OneD, NekDouble> >&outarray,
				 const NekDouble time,
				 const NekDouble lambda)
    {
        int nvariables = inarray.num_elements();
        int nq = m_fields[0]->GetNpoints();

	for(int n=0; n<nvariables; ++n)
	  {
	    Vmath::Vcopy(nq, &inarray[n][0], 1, &outarray[n][0], 1);
	  }
	
        StdRegions::ConstFactorMap factors;
        factors[StdRegions::eFactorTau] = 1.0;
	
	Array<OneD, Array< OneD, NekDouble> > F(nvariables);
	F[0] = Array<OneD, NekDouble> (nq*nvariables);
        for(int n = 1; n < nvariables; ++n)
        {
            F[n] = F[n-1] + nq;
        }

        // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
        // inarray = input: \hat{rhs} -> output: \hat{Y}
        // outarray = output: nabla^2 \hat{Y} where \hat = modal coeffs
	
	//Setting boundary conditions
	SetBoundaryConditions(m_phi, time);

	// factors[StdRegions::eFactorLambda] = 1.0/lambda/m_epsu[i];
	factors[StdRegions::eFactorLambda] = 1.0/lambda;
	
	// Multiply 1.0/timestep/lambda
	Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], inarray[m_phi], 1, F[m_phi], 1);
	
	// d phi / dt = \naba \sigma \nabla phi + F
	m_fields[m_phi]->HelmSolve(F[m_phi], m_fields[m_phi]->UpdateCoeffs(), NullFlagList, factors, m_varcoeff);
	m_fields[m_phi]->BwdTrans(m_fields[m_phi]->GetCoeffs(),outarray[m_phi]);
    }

  void MMFEMFHN::ComputeDRReaction(const Array<OneD, const NekDouble> &uarray,
				   const Array<OneD, const NekDouble> &varray,
				   const Array<OneD, const NekDouble> &warray,
				   Array<OneD, NekDouble> &uReaction,
				   Array<OneD, NekDouble> &vReaction,
				   Array<OneD, NekDouble> &wReaction)
  {
    int i;
    int nq = m_fields[0]->GetNpoints();

    const NekDouble *u = &uarray[0];
    const NekDouble *v = &varray[0];
    const NekDouble *w = &warray[0];

    NekDouble *u_new = &uReaction[0];
    NekDouble *v_new = &vReaction[0];
    NekDouble *w_new = &wReaction[0];

    switch(m_TestType)
      {
      case eTestEMFHN:
	{
	  vReaction = Array<OneD, NekDouble>(nq, 0.0);
	  wReaction = Array<OneD, NekDouble>(nq, 0.0);
	}
	break;

	// Reaction for \phi = c1 \phi ( \phi - a)*(1 - \phi) - c2 v
	// Reaction for \psi = b (\phi - d \psi )
      case eFHNStandard:
	{
	  NekDouble a = 0.12;
	  NekDouble b = 0.011;
	  NekDouble c1 = 0.175;
	  NekDouble c2 = 0.03;
	  NekDouble d = 0.55;
	  
	  for (i=0; i<nq; ++i)
	    {
	      *u_new = c1*(*u)*(*u - a)*(1 -(*u)) - c2*(*v);
	      *v_new = b*(*u - d*(*v));
	      *w_new = 0.0;

	      ++u, ++v, ++w, ++u_new, ++v_new, ++w_new;
	    }
	}
	break;

      case eRogers:
	{
	  NekDouble a = 0.13;
	  NekDouble b = 0.013;
	  NekDouble c1 = 0.26;
	  NekDouble c2 = 0.1;
	  NekDouble d = 1.0;

	  for (i=0; i<nq; ++i)
	    {
	      *u_new = c1*(*u)*(*u - a)*(1 -(*u)) - c2*(*u)*(*v);
	      *v_new = b*(*u - d*(*v));
	      *w_new = 0.0;
	      
	      ++u, ++v, ++w, ++u_new, ++v_new, ++w_new;
	    }
	}
	break;	
	
      case eAlievPanf:
	{
	  NekDouble a = 0.15;
	  NekDouble c1 = 8.0;
	  NekDouble c2 = 1.0;
	  NekDouble c0 = 0.002;
	  NekDouble mu1 = 0.2;
	  NekDouble mu2 = 0.3;

	  for (i=0; i<nq; ++i)
	    {
	      *u_new = c1*(*u)*(*u - a)*(1 - (*u)) - c2*(*u)*(*v);
	      *v_new = (c0 + (mu1*(*v)/(mu2+*u)))*(-1.0*(*v) - c1*(*u)*(*u - a -1.0) );
	      *w_new = 0.0;

	      ++u, ++v, ++w, ++u_new, ++v_new, ++w_new;
	    }
	}
	break;

      case eFentonKarma:
	{
	  NekDouble J_fi, J_so, J_si, hc, hv, hr, tau_v_minus;
	  //  NekDouble alpha, beta;
	  
	  for (i=0; i<nq; ++i)
	    {
	      // Heavyside functions
	      hc = (*u < m_u_c) ? 0.0 : 1.0;
	      hv = (*u < m_u_v) ? 0.0 : 1.0;
	      hr = (*u < m_u_r) ? 0.0 : 1.0;
	      
	      // u-gate
	      J_fi = -(*v)*hc*(1 - *u)*(*u - m_u_fi)/m_tau_d;
	
	      // added extra (1-k2*v) term from Cherry&Fenton 2004
	      // J_so = (*u)*(1-hr)*(1-k2*(*v))/tau_0 + (isCF3 ? h3*(*u)*(*y)/tau_r : h3/tau_r);
	      J_so = (*u)*(1-hr)*(1-m_k2*(*v))/m_tau_0 + hr/m_tau_r;
	      
	      J_si = -(*w)*(1 + tanh(m_k1*(*u - m_u_csi)))/(2.0*m_tau_si);
	      
	      *u_new = -J_fi - J_so - J_si;
	      
	      // v-gate
	      // tau_v_minus = hv*m_tau_v1_minus + (1-hv)*m_tau_v2_minus;
	      // alpha = (1-hc)/tau_v_minus;
	      // beta = hc/m_tau_v_plus;
	      // *v_new = alpha / (alpha + beta);
	      
	      tau_v_minus = (1-hv)*m_tau_v1_minus + hv*m_tau_v2_minus;
	      *v_new = (1-hc)*(1-*v)/tau_v_minus - hc*(*v)/m_tau_v_plus;
	      
	      // w-gate
	      // alpha = (1.0-hc)*(1.0-*w)/m_tau_w_minus;
	      // beta = hc/m_tau_w_plus;
	      // *w_new = alpha / (alpha + beta);
	      *w_new = (1-hc)*(1-*w)/m_tau_w_minus - hc*(*w)/m_tau_w_plus;
	      
	      ++u, ++v, ++w, ++u_new, ++v_new, ++w_new;
	    }
	}
	break;

      default:
	break;
      }
  }

  void MMFEMFHN::ComputedFdt(const Array<OneD, const NekDouble> &uarray,
			     const Array<OneD, const NekDouble> &varray,
			     const Array<OneD, const NekDouble> &warray,
			     const Array<OneD, const NekDouble> &dudtarray,
			     const Array<OneD, const NekDouble> &dvdtarray,
			     const Array<OneD, const NekDouble> &dwdtarray,
			     Array<OneD, NekDouble> &outarray)
  {
    int i;
    int nq = m_fields[0]->GetNpoints();

    const NekDouble *u = &uarray[0];
    const NekDouble *v = &varray[0];
    const NekDouble *w = &varray[0];

    const NekDouble *dudt = &dudtarray[0];
    const NekDouble *dvdt = &dvdtarray[0];
    const NekDouble *dwdt = &dwdtarray[0];

    NekDouble *dFdt = &outarray[0];

    switch(m_TestType)
      {
	// Compute (-3 c1 \phi^2 + 2 c_1 (a+1) \phi - a*c_1 ) \dot{\phi} - c_2 \dot{\psi}
      case eFHNStandard:
	{
	  NekDouble a = 0.12;
	  NekDouble c1 = 0.175;
	  NekDouble c2 = 0.03;

	  for (i=0; i<nq; ++i)
	    {
	      *dFdt = (-3.0*c1*(*u)*(*u) + 2.0*c1*(a+1)*(*u) - a*c1)*(*dudt) - c2*(*dvdt);
	      ++u, ++v, ++dudt, ++dvdt, ++dFdt;
	    }
	}
	break;
	
	// Compute (-3 c1 \phi^2 + 2 c_1 (a+1) \phi - a*c_1 ) \dot{\phi} - c2( \phi \dot{\psi} + \psi \dot{\phi})
      case eRogers:
	{
	  NekDouble a = 0.13;
	  NekDouble c1 = 0.26;
	  NekDouble c2 = 0.1;

	  for (i=0; i<nq; ++i)
	    {
	      *dFdt = c1*(-3.0*(*u)*(*u) + 2.0*(a+1)*(*u) - a)*(*dudt) - c2*((*u)*(*dvdt) + (*v)*(*dudt));
	      ++u, ++v, ++dudt, ++dvdt, ++dFdt;
	    }
	}
	break;

	// Compute -K * ( 3 \phi^2 - 2 (a+1) \phi + a  ) \dot{\phi} - ( \phi \dot{\psi} + \psi \dot{\phi})
      case eAlievPanf:
	{
	  NekDouble a = 0.15;
	  NekDouble c1 = 8.0;
	  NekDouble c2 = 1.0;

	  for (i=0; i<nq; ++i)
	    {
	      *dFdt = c1*(-3.0*(*u)*(*u) + 2.0*(a+1)*(*u) - a)*(*dudt) - c2*((*u)*(*dvdt) + (*v)*(*dudt));
	      ++u, ++v, ++dudt, ++dvdt, ++dFdt;
	    }
	}
	break;

      case eFentonKarma:
	{
	  NekDouble Jdot_fi1, Jdot_fi2, Jdot_fi;
	  NekDouble Jdot_so1, Jdot_so2, Jdot_so, Jdot_si;
	  NekDouble hc, hv, hr, csh, cshuc, cshur, tanhc, tanhr;

	  for (i=0; i<nq; ++i)
	    {
	      // Heavyside functions
	      hc = (*u < m_u_c) ? 0.0 : 1.0;
	      hv = (*u < m_u_v) ? 0.0 : 1.0;
	      hr = (*u < m_u_r) ? 0.0 : 1.0;
	      
	      // u-gate
	      // Jdot_fi = -(*v)*hc*(1 - *u)*(*u - m_u_fi)/m_tau_d;

	      if(m_kp>0)
		{
		  tanhc = 0.5*(1.0+tanh(m_kp*((*u)-m_u_c)));
		  tanhr =  0.5*(1.0+tanh(m_kp*((*u)-m_u_r)));
		}

	      else
		{
		  tanhc = hc;
		  tanhr = hr;
		}

	      cshuc = cosh(m_kp*((*u)-m_u_c));
	      cshur = cosh(m_kp*((*u)-m_u_r));
	      
	      
	      Jdot_fi1 = tanhc*( (*v)*(m_u_fi+1-2*(*u))*(*dudt) + ((*u)-m_u_fi)*(1-(*u))*(*dvdt) )/m_tau_d;
	      Jdot_fi2 = 0.5*m_kp/cshuc/cshuc*(*v)*((*u)-m_u_fi)*(1-(*u))*(*dudt)/m_tau_d;
	      Jdot_fi = Jdot_fi1+Jdot_fi2;
	
	      // added extra (1-k2*v) term from Cherry&Fenton 2004
	      // J_so = (*u)*(1-hr)*(1-k2*(*v))/tau_0 + (isCF3 ? h3*(*u)*(*y)/tau_r : h3/tau_r);
	      // Jdot_so = (*u)*(1-hr)*(1-m_k2*(*v))/m_tau_0 + hr/m_tau_r;

	      Jdot_so1 = (1-tanhr)*( -(1 - m_k2*(*v))*(*dudt) + m_k2*(*u)*(*dvdt) )/m_tau_0;
	      Jdot_so2 = 0.5*m_kp/cshur/cshur*((*u)*(1-(*v)*m_k2)/m_tau_0 - 1/m_tau_r)*(*dudt);
	      Jdot_so = Jdot_so1 + Jdot_so2;

	      // Jdot_si = -(*w)*(1 + tanh(m_k1*(*u - m_u_csi)))/(2.0*m_tau_si);
	      csh = cosh(m_k1*(*u - m_u_csi));	
	      Jdot_si = ( m_k1*(*w)/csh/csh*(*dudt) + (1 + tanh(m_k1*(*u - m_u_csi)))*(*dwdt) )/(2.0*m_tau_si);

	      *dFdt = Jdot_fi + Jdot_so + Jdot_si; 
	      
	      ++u, ++v, ++w, ++dudt, ++dvdt, ++dwdt, ++dFdt; 
	    }
	}
	break;

      default:
	break;
	
      }
  }
  
  void MMFEMFHN::AddGreenDerivCompensate(const Array<OneD, const Array<OneD, NekDouble> > &physarray,
					   Array<OneD, Array<OneD, NekDouble> > &outarray)
  {    
    // routine works for both primitive and conservative formulations
    int ncoeffs = outarray[0].num_elements();
    int nq      = physarray[0].num_elements();
    
    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> tmpc(ncoeffs);

    Array<OneD, Array<OneD, NekDouble> > fluxvector (m_shapedim);
    for(int j=0; j<m_shapedim; ++j)
      {
	fluxvector[j] = Array<OneD, NekDouble>(nq);
      }

    // m_CurlMF[0][0] = e^3 \cdot (\nabla \times e^1) [ NEW m_CurlMF[0][2] ]
    // m_CurlMF[0][1] = 0.0                           
    // m_CurlMF[1][0] = 0.0,
    // m_CurlMF[1][1] = e^3 \cdot (\nabla \times e^2) [ NEW m_CurlMF[1][2]  ]
    // m_CurlMF[2][0] = e^1 \cdot (\nabla \times e^3) [ NEW m_CurlMF[2][0]  ]
    // m_CurlMF[2][1] = e^2 \cdot (\nabla \times e^3) [ NEW m_CurlMF[2][1]  ]
    
    int var;
    
    var = 0;
    GetMaxwellFluxVector(var, physarray, fluxvector);
    Vmath::Vmul(nq, &fluxvector[0][0], 1, &m_CurlMF[0][2][0], 1, &tmp[0], 1);
    m_fields[var]->IProductWRTBase(tmp,tmpc);
    Vmath::Vadd(ncoeffs, tmpc, 1, outarray[var], 1, outarray[var], 1);

    var = 1;
    GetMaxwellFluxVector(var, physarray, fluxvector);
    Vmath::Vmul(nq, &fluxvector[1][0], 1, &m_CurlMF[1][2][0], 1, &tmp[0], 1);
    Vmath::Neg(nq, tmp, 1);
    m_fields[var]->IProductWRTBase(tmp,tmpc);
    Vmath::Vadd(ncoeffs, tmpc, 1, outarray[var], 1, outarray[var], 1);

    var = 2;
    GetMaxwellFluxVector(var, physarray, fluxvector);
    Vmath::Vmul(nq, &fluxvector[0][0], 1, &m_CurlMF[2][0][0], 1, &tmp[0], 1);
    Vmath::Vvtvm(nq, &fluxvector[1][0], 1, &m_CurlMF[2][1][0], 1, &tmp[0], 1, &tmp[0], 1);
    m_fields[var]->IProductWRTBase(tmp,tmpc);
    Vmath::Vadd(ncoeffs, tmpc, 1, outarray[var], 1, outarray[var], 1);
  }

  
    /**
     * @brief Compute the projection for the linear advection equation.
     *
     * @param inarray    Given fields.
     * @param outarray   Calculated solution.
     * @param time       Time.
     */
    void MMFEMFHN::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> >&inarray,
              Array<OneD,       Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
    {
      //  int nvar = inarray.num_elements();
      
      // for (int i=m_D1; i<nvar; ++i)
      //	{
      //   SetBoundaryConditions(i,time);
      //	}
    }

  void MMFEMFHN::v_EvaluateExactSolution(unsigned int field,
					     Array<OneD, NekDouble> &outfield,
					     const NekDouble time)
  
    {
      int nq = m_fields[0]->GetNpoints();

       switch(m_TestType)
	 {
	 case eTestEMFHN:
	   {
	     outfield = FlowTestEMFHN(time,field);
	   }
	   break;

	 default:
	   {
	     outfield = Array<OneD, NekDouble>(nq, 0.0);
	   }
	   break;
	 }
    }



  void MMFEMFHN::v_SetInitialConditions(const NekDouble initialtime,
					    bool dumpInitialConditions,
					    const int domain)
  {
    int nq  = GetTotPoints();
    int nvar = m_fields.num_elements();
    int ncoeffs    = GetNcoeffs();
    
    // Test SolvePoissson
    Array<OneD, NekDouble> tmp(nq,0.0);
    Array<OneD, NekDouble> tmpc(ncoeffs,0.0);	    
    
    NekDouble testtime = m_pi/6;
    
    Array<OneD, NekDouble> phisctmp(nq, 0.0);
    Array<OneD, NekDouble> phiscexact(nq, 0.0);

    Array<OneD, Array<OneD, NekDouble> > Diffphisc(m_shapedim);
    Array<OneD, Array<OneD, NekDouble> > Directgradphisc(m_shapedim);
    Array<OneD, Array<OneD, NekDouble> > gradphisc(m_shapedim);
    Array<OneD, Array<OneD, NekDouble> > gradphiscexact(m_shapedim);
    
    tmp = FlowTestEMFHN(testtime,40);
    SolvePoisson(tmp, testtime);

    Vmath::Vcopy(nq, m_fields[m_Fs]->GetPhys(), 1, phisctmp, 1);

    for (int j=0; j<m_shapedim; ++j)
      {
	Diffphisc[j] = Array<OneD, NekDouble>(nq);
	Directgradphisc[j] = Array<OneD, NekDouble>(nq);
	gradphisc[j] = Array<OneD, NekDouble>(nq);
	gradphiscexact[j] = Array<OneD, NekDouble>(nq);
	
	WeakDGAdv(j, j, phisctmp, tmpc);
	
	m_fields[j]->MultiplyByElmtInvMass(tmpc,tmpc);
	m_fields[j]->BwdTrans(tmpc, gradphisc[j]);

	m_fields[j]->PhysDirectionalDeriv(m_movingframes[j], phisctmp, Directgradphisc[j]);
      }

    phiscexact        = FlowTestEMFHN(testtime, 30);
    gradphiscexact[0] = FlowTestEMFHN(testtime, 31);
    gradphiscexact[1] = FlowTestEMFHN(testtime, 32);
    
    Vmath::Vsub(nq, phisctmp, 1, phiscexact, 1, phisctmp, 1);
    Vmath::Sadd(nq, -1.0*phisctmp[0], phisctmp, 1, phisctmp, 1);

    Vmath::Vsub(nq, Directgradphisc[0], 1, gradphisc[0], 1, Diffphisc[0], 1);
    Vmath::Vsub(nq, Directgradphisc[1], 1, gradphisc[1], 1, Diffphisc[1], 1);
    
    Vmath::Vsub(nq, gradphiscexact[0], 1, gradphisc[0], 1, gradphisc[0], 1);
    Vmath::Vsub(nq, gradphiscexact[1], 1, gradphisc[1], 1, gradphisc[1], 1);

    Vmath::Vsub(nq, gradphiscexact[0], 1, Directgradphisc[0], 1, Directgradphisc[0], 1);
    Vmath::Vsub(nq, gradphiscexact[1], 1, Directgradphisc[1], 1, Directgradphisc[1], 1);

    if(m_TestType==eTestEMFHN)
      {
	Gs::cout << "Error of SolvePoisson = " << AbsIntegral(phisctmp) << Gs::endl;
	Gs::cout << "CoeffsDeriv = ( " << AbsIntegral(gradphisc[0]) << ", " << AbsIntegral(gradphisc[1]) << " ) " << Gs::endl;
	Gs::cout << "DirectDeriv = ( " << AbsIntegral(Directgradphisc[0]) << ", " << AbsIntegral(Directgradphisc[1]) << " ) " << Gs::endl;
	Gs::cout << "Difference = ( " << AbsIntegral(Diffphisc[0]) << ", " << AbsIntegral(Diffphisc[1]) << " ) " << Gs::endl << Gs::endl;
      }

    switch(m_TestType)
      {
      case eFHNStandard:
      case eRogers:
      case eAlievPanf:
	{
	  Array<OneD, NekDouble> Zero(nq,0.0);
	  Array<OneD, NekDouble> Initpulse(nq,0.0);

	  /*
	  // Assign  E = - \nabla \phi
	  // Assign  D = - \nabla \phi
	  Array<OneD, NekDouble> NablaInitpulse(nq,0.0);
	  for (int j=0; j<m_shapedim; j++)
	    {
	      WeakDGAdv(m_phi+j, j, Initpulse, tmpc);
	      
	      m_fields[m_phi]->MultiplyByElmtInvMass(tmpc, tmpc);
	      m_fields[m_phi]->BwdTrans(tmpc,NablaInitpulse);
	    }
	  */

	    for(int i = 0; i < m_fields.num_elements(); ++i)
	      {
		m_fields[i]->SetPhys(Zero);
	      }

	    // Assign \phi
	    PlaneWaveForPhi(Initpulse);
	    m_fields[m_phi]->SetPhys(Initpulse);
	  }
	  break;

      case eFentonKarma:
	{
	  Array<OneD, NekDouble> Zero(nq,0.0);
	  Array<OneD, NekDouble> Ones(nq,1.0);
	  Array<OneD, NekDouble> Initpulse(nq,0.0);

	  SetUpParametersFentonKarma();

	  for(int i = 0; i < m_fields.num_elements(); ++i)
	    {
	      m_fields[i]->SetPhys(Zero);
	    }
	  m_fields[m_psi]->SetPhys(Ones);
	  m_fields[m_chi]->SetPhys(Ones);
	  
	  // Assign \phi
	  PlaneWaveForPhi(Initpulse);
	  m_fields[m_phi]->SetPhys(Initpulse);
	}
	break;
	
      case eTestEMFHN:
	  {
	    for(int i = 0; i < m_fields.num_elements(); ++i)
	      {
		m_fields[i]->SetPhys(FlowTestEMFHN(initialtime, i));
	      } 
	  }
	  break;

	default:
	  break;
        }

    // forward transform to fill the modal coeffs
    for(int i = 0; i < m_fields.num_elements(); ++i)
      {
	m_fields[i]->SetPhysState(true);
	m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
      }       
    
    Array<OneD, Array<OneD, NekDouble> > fields(nvar);
    for (int i=0; i<nvar; i++)
      {
	fields[i] = Array<OneD, NekDouble>(nq);
	Vmath::Vcopy(nq, &(m_fields[i]->GetPhys())[0], 1, &fields[i][0], 1);
      }
    
  // 0: |E|
  // 1: |H|
  // 2: |U|
  // 3: Sx
  // 4: Sy
  // 5: Jx
  // 6: Jy
  // 7: DTheta
  // 8: \nabla \times J
    Array<OneD, Array<OneD, NekDouble> >Energy;
    ComputeEnergy(fields, Energy);
    
    Gs::cout << "Initial: |E| = " << AbsIntegral(Energy[0]) << ", |H| = " << AbsIntegral(Energy[1]) << ", U = " << AbsIntegral(Energy[2])
	 << ", phi = " << AbsIntegral(fields[m_phi]) << ", Fs = " << AbsIntegral(m_fields[m_Fs]->GetPhys()) << Gs::endl;

    // dump initial conditions to file
    std::string outname = m_sessionName + "_initial.chk";
    WriteFld(outname);
  }

  
  // WeakDGAdv(indx, dir, inarray, outarray)
  // outarray = \nabla_{dir} inarray
  void MMFEMFHN::WeakDGAdv(const int indx,
			   const int dir,
			   const Array<OneD, const NekDouble> &inarray,
			   Array<OneD, NekDouble> &outarray)

  {
    int ncoeffs         = GetNcoeffs();
    int nTracePointsTot = GetTraceNpoints();

    outarray = Array<OneD, NekDouble>(ncoeffs, 0.0);
    
    m_fields[indx]->IProductWRTDirectionalDerivBase(m_movingframes[dir], inarray, outarray);

    // Assign Boundary Conditions
    Array<OneD, NekDouble> Fwd(nTracePointsTot);
    Array<OneD, NekDouble> Bwd(nTracePointsTot);
    m_fields[indx]->GetFwdBwdTracePhys(inarray,Fwd,Bwd);

    Array<OneD, NekDouble> numflux(nTracePointsTot);
    m_fields[indx]->GetTrace()->Upwind(m_Vn[dir],Fwd,Bwd,numflux);

    // calculate numflux = (n \cdot MF)*flux
    Vmath::Vmul(nTracePointsTot, &numflux[0], 1, &m_ncdotMFFwd[dir][0], 1, &Fwd[0], 1);
    Vmath::Vmul(nTracePointsTot, &numflux[0], 1, &m_ncdotMFBwd[dir][0], 1, &Bwd[0], 1);
    
    Vmath::Neg(ncoeffs,outarray,1);	    
    m_fields[indx]->AddFwdBwdTraceIntegral(Fwd, Bwd, outarray);
  }

  
  Array<OneD, NekDouble> MMFEMFHN::FlowTestEMFHN(const NekDouble time, const int var)

  {
        int nq  = GetTotPoints();

        Array<OneD, NekDouble> x(nq), y(nq), z(nq);
        m_fields[0]->GetCoords(x,y,z);

	Array<OneD, NekDouble> outfield(nq);

	NekDouble freqm, freqn, omega;
	NekDouble mpi, npi, xp, yp;
	
	freqm = 1.0;
	freqn = 1.0;
	mpi = freqm*m_pi;
	npi = freqn*m_pi;
	
	omega = m_pi*sqrt(freqm*freqm+freqn*freqn);

	Array<OneD, NekDouble> dphidt(nq), phi(nq), phisc(nq);
	Array<OneD, NekDouble> F(nq), Fdot(nq), phisc_dx(nq), phisc_dy(nq);
	
	Array<OneD, NekDouble> A1(nq),  A2(nq), E1(nq), E2(nq), Hz(nq), dHzdt(nq), D1(nq), D2(nq);
	Array<OneD, NekDouble> E1Force(nq), E2Force(nq), HzForce(nq);
	
	NekDouble Ex, Ey, ExForce, EyForce, Ax, Ay, Dx, Dy, phix, phiy, phiscx, phiscy;

	switch(m_TestType)
	  {
	  case eTestEMFHN:
	    {
	      for (int i=0; i<nq; ++i)
		{
		  mpi = freqm*m_pi;
		  npi = freqn*m_pi;
		  
		  xp = x[i];
		  yp = y[i];

		  // modarray[m_E1var] = d H^z / d \xi^1 -  d phi_sc / dx
		  // modarray[m_E2var] = - d H^z / d \xi^2 - d phi_sc / dy
		  ExForce = -(m_pi/omega)*sin(mpi*xp)*cos(npi*yp)*( -omega*cos(omega*time) + sin(omega*time) );
		  EyForce = -(m_pi/omega)*cos(mpi*xp)*sin(npi*yp)*( -omega*cos(omega*time) + sin(omega*time) );
		  		  
		  // Ex = -(m_pi/omega)*cos(m_pi*xp)*sin(m_pi*yp)*sin(omega*time);
		  // Ey = (m_pi/omega)*sin(m_pi*xp)*cos(m_pi*yp)*sin(omega*time);
		  // Hz[i] = cos(m_pi*xp)*cos(m_pi*yp)*cos(omega*time);
		  Ex =  (m_pi/omega)*sin(mpi*xp)*cos(npi*yp)*( (1-omega)*cos(omega*time) + sin(omega*time) );
		  Ey =  (m_pi/omega)*cos(mpi*xp)*sin(npi*yp)*( -1.0*(1+omega)*cos(omega*time) + sin(omega*time) );
		  Hz[i] = -sin(mpi*xp)*sin(npi*yp)*sin(omega*time);		  
		  dHzdt[i] = omega*sin(mpi*xp)*sin(npi*yp)*cos(omega*time);

		  // Satisfy d \phi dt = \nabla \cdot \sigma \nabla \phi + F
		  phi[i] = cos(mpi*xp)*cos(npi*yp)*cos(omega*time);
		  phix = -mpi/omega*sin(mpi*xp)*cos(npi*yp)*sin(omega*time);
		  phiy = -npi/omega*cos(mpi*xp)*sin(npi*yp)*sin(omega*time);
		  
		  dphidt[i] = -1.0*omega*cos(mpi*xp)*cos(npi*yp)*sin(omega*time);
		  F[i] = dphidt[i] + m_pi*m_pi*(m_SigmaBlock[0][i]+m_SigmaBlock[1][i])*cos(mpi*xp)*cos(npi*yp)*cos(omega*time);

		  // phi[i] = exp(-mpi*mpi*time)*cos(mpi*xp)*cos(npi*yp);
		  // dphidt[i] = mpi*mpi*phi[i];
		  // F[i] = (m_SigmaBlock[0][i]+m_SigmaBlock[1][i] - 1)*mpi*mpi*phi[i];
		  
		  // \nabla \cdot \sigma \nabla phisc = Fdot
		  Fdot[i] = -m_pi*m_pi*(m_SigmaBlock[0][i]+m_SigmaBlock[1][i])*(omega*sin(omega*time)+cos(omega*time))*cos(mpi*xp)*cos(npi*yp);
		  phisc[i] = (omega*sin(omega*time)+cos(omega*time))*cos(mpi*xp)*cos(npi*yp);

		  phiscx = -1.0*mpi*(omega*sin(omega*time)+cos(omega*time))*sin(mpi*xp)*cos(npi*yp);
		  phiscy = -1.0*npi*(omega*sin(omega*time)+cos(omega*time))*cos(mpi*xp)*sin(npi*yp);

		  // modarray[m_D1var] = d H_z/dy + \sigma_x( \nabla_x phi  )
		  // modarray[m_D2var] = - d H_z/dx + \sigma_y ( \nabla_y phi )
		  Dx = phix*m_SigmaBlock[0][i] + (m_pi/omega)*sin(m_pi*xp)*cos(m_pi*yp)*cos(omega*time);
		  Dy = phiy*m_SigmaBlock[1][i] - (m_pi/omega)*cos(m_pi*xp)*sin(m_pi*yp)*cos(omega*time);

		  // Satisfy d A_x / dt = - d \phi / d x - E_x
		  Ax = m_SigmaBlock[0][i]*(-phix - (m_pi/omega/omega)*cos(m_pi*xp)*sin(m_pi*yp)*cos(omega*time));
		  Ay = m_SigmaBlock[0][i]*(-phiy + (m_pi/omega/omega)*sin(m_pi*xp)*cos(m_pi*yp)*cos(omega*time));

		  phisc_dx[i] = phiscx*m_movingframes[0][i] + phiscy*m_movingframes[0][i+nq];
		  phisc_dy[i] = phiscx*m_movingframes[1][i] + phiscy*m_movingframes[1][i+nq];
		  
		  E1[i] = Ex*m_movingframes[0][i] + Ey*m_movingframes[0][i+nq];
		  E2[i] = Ex*m_movingframes[1][i] + Ey*m_movingframes[1][i+nq];

		  E1Force[i] = ExForce*m_movingframes[0][i] + EyForce*m_movingframes[0][i+nq];
		  E2Force[i] = ExForce*m_movingframes[1][i] + EyForce*m_movingframes[1][i+nq];
		  
		  D1[i] = Dx*m_movingframes[0][i] + Dy*m_movingframes[0][i+nq];
		  D2[i] = Dx*m_movingframes[1][i] + Dy*m_movingframes[1][i+nq];
		  
		  A1[i] = Ax*m_movingframes[0][i] + Ay*m_movingframes[0][i+nq];
		  A2[i] = Ax*m_movingframes[1][i] + Ay*m_movingframes[1][i+nq];
		}
	    }
	    break;

	  default:
	    {
	      NekDouble xsize = Vmath::Vmax(nq, x, 1);
	      NekDouble ysize = Vmath::Vmax(nq, y, 1);
	      
	      for (int i=0; i<nq; ++i)
		{
		  mpi = freqm*m_pi;
		  npi = freqn*m_pi;
		  
		  xp = x[i];
		  yp = y[i];
		  
		  Fdot[i] = -0.5*omega*omega*(1/xsize/xsize+1/ysize/ysize)*cos(mpi*xp/xsize)*cos(npi*yp/ysize)*(omega*sin(omega*time)+cos(omega*time));
		  
		  phisc[i] = (omega*sin(omega*time)+cos(omega*time))*cos(mpi*xp/xsize)*cos(npi*yp/ysize);
		  phisc_dx[i] =  -1.0*mpi/xsize*(omega*sin(omega*time)+cos(omega*time))*sin(mpi*xp/xsize)*cos(npi*yp/ysize);
		  phisc_dy[i] =  -1.0*npi/ysize*(omega*sin(omega*time)+cos(omega*time))*cos(mpi*xp/xsize)*sin(npi*yp/ysize);
		}
	    }
	    break;
	  }
	
	switch(var)
	  {	
	  case(0):
	    {
	      outfield = E1;
	    }
	    break;
	    
	  case(1):
	    {
	      outfield = E2;
	    }
	    break;

	  case(7):
	    {
	      outfield = D1;
	    }
	    break;
	    
	  case(8):
	    {
	      outfield = D2;
	    }
	    break;

	  case(10):
	    {
	      outfield = E1Force;
	    }
	    break;
	    
	  case(11):
	    {
	      outfield = E2Force;
	    }
	    break;

	  case(2):
	    {
	      outfield = Hz;
	    }
	    break;

	  case(102):
	    {
	      outfield = dHzdt;
	    }
	    break;
	    
	  case(3):
	    {
	      outfield = phi;
	    }
	    break;

	  case(100):
	    {
	      outfield = F;
	    }
	    break;
	    
	  case(20):
	    {
	      outfield = dphidt;
	    }
	    break;

	  case(30):
	    {
	      outfield = phisc;
	    }
	    break;

	  case(31):
	    {
	      outfield = phisc_dx;
	    }
	    break;
	    
	  case(32):
	    {
	      outfield = phisc_dy;
	    }
	    break;

	  case(40):
	    {
	      outfield = Fdot;
	    }
	    break;
	    
	  default:
	    {
	      outfield = Array<OneD, NekDouble>(nq,0.0);
	    }
	    break;						     
	  }

	return outfield;
  }

  void MMFEMFHN::SolvePoisson(const Array<OneD, const NekDouble> &inarray,
			      const NekDouble time)
  {
    int nq = m_fields[0]->GetTotPoints();
    int index = m_Fs;
    
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorLambda] = 0.0;
    factors[StdRegions::eFactorTau] = m_PoissonTau;

    m_fields[index]->SetPhys(inarray);
    SetBoundaryConditions(index, time);
    
    if(m_PoissonNeumann)
      {
	// Adjust force such that \int f = 0
	NekDouble avg = -1.0*(m_fields[0]->PhysIntegral(m_fields[index]->GetPhys()))/m_Jac;
	Vmath::Sadd(nq, avg, m_fields[index]->GetPhys(), 1, m_fields[index]->UpdatePhys(), 1);
      }

    m_fields[index]->HelmSolve(m_fields[index]->GetPhys(), m_fields[index]->UpdateCoeffs(), NullFlagList, factors, m_varcoeff);
    m_fields[index]->BwdTrans(m_fields[index]->GetCoeffs(), m_fields[index]->UpdatePhys());
  }

  
  void MMFEMFHN::ComputeVarcoeffDiffusion()
  {
    int MFdim = 3;
    int nq = m_fields[0]->GetTotPoints();
    
    StdRegions::VarCoeffType MMFCoeffs[15] = {StdRegions::eVarCoeffMF1x,
					       StdRegions::eVarCoeffMF1y,
					       StdRegions::eVarCoeffMF1z,
					       StdRegions::eVarCoeffMF1Div,
					       StdRegions::eVarCoeffMF1Mag,
					       StdRegions::eVarCoeffMF2x,
					       StdRegions::eVarCoeffMF2y,
					       StdRegions::eVarCoeffMF2z,
					       StdRegions::eVarCoeffMF2Div,
					       StdRegions::eVarCoeffMF2Mag,
					       StdRegions::eVarCoeffMF3x,
					       StdRegions::eVarCoeffMF3y,
					       StdRegions::eVarCoeffMF3z,
					       StdRegions::eVarCoeffMF3Div,
					       StdRegions::eVarCoeffMF3Mag};

      int indx;
      Array<OneD, NekDouble> tmp(nq);
      Array<OneD, NekDouble> Dtmp(nq);
	for (int k=0; k<MFdim; ++k)
	  {
	    // For Moving Frames
	    indx = 5*k;

	    for (int i=0; i<m_spacedim; ++i)
	      {
		m_varcoeff[MMFCoeffs[indx+i]] = Array<OneD, NekDouble>(nq, 0.0);
		Vmath::Vcopy(nq, &m_movingframes[k][i*nq], 1, &m_varcoeff[MMFCoeffs[indx+i]][0], 1);

		Vmath::Vsqrt(nq, m_SigmaBlock[i], 1, tmp, 1);
		Vmath::Vmul(nq, &tmp[0], 1, &m_movingframes[k][i*nq], 1, &m_varcoeff[MMFCoeffs[indx+i]][0], 1);
	      }
	    
	    
	    // m_DivMF
	    m_varcoeff[MMFCoeffs[indx+3]] = Array<OneD, NekDouble>(nq, 0.0);
	    for (int i=0; i<m_spacedim; ++i)
	      {	    	      
		Vmath::Vcopy(nq, &m_varcoeff[MMFCoeffs[indx+i]][0], 1, &tmp[0], 1);
		
		m_fields[0]->PhysDeriv(i, tmp, Dtmp);
		Vmath::Vadd(nq, &Dtmp[0], 1, &m_varcoeff[MMFCoeffs[indx+3]][0], 1, &m_varcoeff[MMFCoeffs[indx+3]][0], 1);
	      }

	    
	    // \| e^k \|
	    m_varcoeff[MMFCoeffs[indx+4]] = Array<OneD, NekDouble>(nq,0.0);
	    tmp = Array<OneD, NekDouble>(nq,0.0);
	    for (int i=0; i<m_spacedim; ++i)
	      {
		Vmath::Vvtvp(nq, &m_varcoeff[MMFCoeffs[indx+i]][0], 1, &m_varcoeff[MMFCoeffs[indx+i]][0], 1, &tmp[0], 1, &tmp[0], 1);
	      }
	    
	    Vmath::Vcopy(nq, &tmp[0], 1, &m_varcoeff[MMFCoeffs[indx+4]][0], 1);
	  }
  }
  

  void MMFEMFHN::PlaneWaveForPhi(Array<OneD, NekDouble> &outarray)
  {
    int nq  = GetTotPoints();
    
    outarray = Array<OneD, NekDouble>(nq,0.0);
    
    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Array<OneD, NekDouble> z(nq);
    
    m_fields[0]->GetCoords(x,y,z);
    
    NekDouble xmin, ymin, xmax, ymax;

    xmin = Vmath::Vmin(nq, x, 1);
    xmax = Vmath::Vmax(nq, x, 1);
    ymin = Vmath::Vmin(nq, y, 1);
    ymax = Vmath::Vmax(nq, y, 1);

    NekDouble radiusofinit, frontstiff;
    switch(m_TestType)
      {
      case eFHNStandard:
	{
	  radiusofinit = 6.0;
	  frontstiff = 0.01;
	}
      break;
      
      case eAlievPanf:
	{
	  radiusofinit = 20.0;
	  frontstiff = 0.01;
	}
      break;
	      
      case eFentonKarma:
	{
	  //radiusofinit = m_radiusofinit;
	  // frontstiff = 0.1;
	  
	  radiusofinit = 20.0;
	  frontstiff = 0.5;
	  
	}
      break;
	      
      default:
	{
	  radiusofinit = 4.0;
	  frontstiff = 0.1;
	}

      }
    
    NekDouble xp, yp, xp2;
    for (int i=0; i<nq; i++)
      {
	switch(m_InitWaveType)
	  {
	  case eLeft:
	    {
	      xp = x[i] - xmin;
	      outarray[i] = 2.0/( 1.0 + exp( ( xp - radiusofinit)/frontstiff ) );

	    }
	    break;
	
	  case eBothEnds:
	    {
	      NekDouble radiusofinit = 3.0;
	      NekDouble frontstiff = 0.1;
	      
	      xp = x[i] - xmin;
	      xp2 = x[i] - xmax;
	      
	      // outarray[i] = 2.0/(1.0 + exp( 0.5*( sqrt(xp*xp) - 0.1) ) ) + 2.0/(1.0 + exp( 0.5*( sqrt(xp2*xp2) - 0.1) ) );
	      outarray[i] = 1.0/( 1.0 + exp( ( sqrt(xp*xp) - radiusofinit)/frontstiff ) ) + 1.0/( 1.0 + exp( ( sqrt(xp2*xp2) - radiusofinit)/frontstiff ) );
	    }
	    break;

	  case eCenter:
	    {
	      // NekDouble radiusofinit = 6.0;
	      NekDouble frontstiff = 0.1;

	      // NekDouble xc;
	      // xc = 0.5*(Vmath::Vmax(nq, x, 1) + Vmath::Vmin(nq, x, 1));

	      xp = x[i] - xmin;
	      outarray[i] =1.0/( 1.0 + exp( ( xp - m_radiusofinit)/frontstiff ) );
	    }
	    break;

	  case eLeftTopCorner:
	    {
	      NekDouble radiusofinit = 6.0;
	      NekDouble frontstiff = 0.1;
	      NekDouble bs = 2.0;
	      
	      xp = x[i] - xmin;
	      yp = y[i] - ymax;
	      outarray[i] = 1.0/( 1.0 + exp( ( sqrt(xp*xp+yp*yp)/bs - radiusofinit)/frontstiff ) );
	    }
	    break;

	  case ePoint:
	    {
	      NekDouble xloc, yloc, zloc, rad;
	      // NekDouble radiusofinit = 10.0;

	      xloc = x[i]-m_InitPtx;
	      yloc = y[i]-m_InitPty;
	      zloc = z[i]-m_InitPtz;
	      
	      rad = sqrt(xloc*xloc + yloc*yloc + zloc*zloc);
	      
	      xloc = xloc/m_radiusofinit;
	      yloc = yloc/m_radiusofinit;
	      zloc = zloc/m_radiusofinit;
	      
	      if(rad<m_radiusofinit)
		{
		  outarray[i] = exp( -(1.0/2.0)*( xloc*xloc + yloc*yloc + zloc*zloc) ) ;
		}
	      
	      else
		{
		  outarray[i] = 0.0;
		}
	    }
	    break;

	  case eSpiralDock:
	    {
	      NekDouble radiusofinit = 3.0;
	      NekDouble frontstiff = 0.1;
	      xp = x[i] - 4.0;
	      yp = y[i];
	      outarray[i] = (1.0/(1.0+exp(2.0*yp)))*(1.0/(1.0+exp(-2.0*xp)))*( 1.0/( 1.0 + exp( ( xp - radiusofinit)/frontstiff ) ) );	      
	    }
	    break;

	  default:
	    break;
	  }

      }
  }

  
  // 0: |E|
  // 1: |H|
  // 2: |U|
  // 3: Sx
  // 4: Sy
  // 5: Jx
  // 6: Jy
  // 7: DTheta
  // 8: \nabla \times J

  void MMFEMFHN::ComputeEnergy(const Array<OneD, const Array<OneD, NekDouble> > &fields, 
			       Array<OneD, Array<OneD, NekDouble> > &outarray)
  {
    int nq  = GetTotPoints();
    int nvar = m_fields.num_elements();
    
    outarray = Array<OneD, Array<OneD, NekDouble> >(nvar);
    for (int i=0; i<nvar; ++i)
      {
	outarray[i] = Array<OneD, NekDouble>(nq,0.0);
      }

    // Ex
    Array<OneD, NekDouble> Ex(nq);
    Vmath::Vmul(nq, &fields[m_E1][0], 1, &m_movingframes[0][0], 1, &Ex[0], 1);
    Vmath::Vvtvp(nq, &fields[m_E2][0], 1, &m_movingframes[1][0], 1, &Ex[0], 1, &Ex[0], 1);

    // Ey
    Array<OneD, NekDouble> Ey(nq);
    Vmath::Vmul(nq, &fields[m_E1][0], 1, &m_movingframes[0][nq], 1, &Ey[0], 1);
    Vmath::Vvtvp(nq, &fields[m_E2][0], 1, &m_movingframes[1][nq], 1, &Ey[0], 1, &Ey[0], 1);

    /*
    // Dx
    Array<OneD, NekDouble> Dx(nq);
    Vmath::Vmul(nq, &fields[m_D1][0], 1, &m_movingframes[0][0], 1, &Dx[0], 1);
    Vmath::Vvtvp(nq, &fields[m_D2][0], 1, &m_movingframes[1][0], 1, &Dx[0], 1, &Dx[0], 1);

    // Dy
    Array<OneD, NekDouble> Dy(nq);
    Vmath::Vmul(nq, &fields[m_D1][0], 1, &m_movingframes[0][nq], 1, &Dy[0], 1);
    Vmath::Vvtvp(nq, &fields[m_D2][0], 1, &m_movingframes[1][nq], 1, &Dy[0], 1, &Dy[0], 1);
    */

    
    // 0) |E| = sqrt(E1*E1 + E2*E2)
    Vmath::Vmul(nq, Ex, 1, Ex, 1, outarray[0], 1);
    Vmath::Vvtvp(nq, Ey, 1, Ey, 1, outarray[0], 1, outarray[0], 1);
    Vmath::Vsqrt(nq, outarray[0], 1, outarray[0], 1);

    // 1) |H| = |Hz|
    Vmath::Vabs(nq, fields[m_Hz], 1, outarray[1], 1);

    // Classical 2) Ue = E1*D1 + E2*D2 + H3^2
    // 2) Ue = E1*E1 + E2*E2 + H3^2
    Vmath::Vmul(nq, Ex, 1, Ex, 1, outarray[2], 1);
    Vmath::Vvtvp(nq, Ey, 1, Ey, 1, outarray[2], 1, outarray[2], 1);
    Vmath::Vvtvp(nq, fields[m_Hz], 1, fields[m_Hz], 1, outarray[2], 1, outarray[2], 1);

    
    // 3) Sex = Ey*Hz
    Vmath::Vmul(nq, Ey, 1, fields[m_Hz], 1, outarray[3], 1);
		
    // 4) Sey = -Ex*Hz
    Vmath::Vmul(nq, Ex, 1, fields[m_Hz], 1, outarray[4], 1);
    Vmath::Neg(nq, outarray[4], 1);
    
    // 5) Jx = \sigma_x Ex
    Vmath::Vmul(nq, m_SigmaBlock[0], 1, Ex, 1, outarray[5], 1);

    // 6) Jy = \sigma_y Ey
    Vmath::Vmul(nq, m_SigmaBlock[1], 1, Ey, 1, outarray[6], 1);

    /*
    // 7) Angle between J and D (Theta) = acos( (J \cdot D)/|J| /|D| )
    Vmath::Vmul(nq, Dx, 1, outarray[5], 1, outarray[7], 1);
    Vmath::Vvtvp(nq, Dy, 1, outarray[6], 1, outarray[7], 1, outarray[7], 1);
    Vmath::Vabs(nq, outarray[7], 1, outarray[7], 1);
    
    NekDouble Jmag, Dmag, theta;
    NekDouble Tol = 0.01;
    for (int i=0; i<nq; i++)
      {
	Jmag = sqrt(outarray[5][i]*outarray[5][i] + outarray[6][i]*outarray[6][i]);
	Dmag = sqrt(Dx[i]*Dx[i] + Dy[i]*Dy[i]);

	if( (fabs(Dmag)>Tol) && (fabs(Jmag)>Tol) )
	  {
	    theta = acos(outarray[7][i]/Jmag/Dmag);
	    if(theta>0.5*m_pi)
	      {
		outarray[7][i] = (180.0/m_pi)*(m_pi - theta);
	      }

	    else
	      {
		outarray[7][i] = (180.0/m_pi)*theta;
	      }
	  }
      }
    */

    /*
    // 8) \nabla \times J = \nabla_1 J2 - \nabla_2 J1
    Array<OneD, Array<OneD, NekDouble> > Jvec(m_spacedim);
    Array<OneD, Array<OneD, NekDouble> > Jcurl(m_spacedim);
    
    for(int k=0; k<m_spacedim; ++k)
      {
	Jvec[k] = Array<OneD, NekDouble>(nq,0.0);
	Jcurl[k] = Array<OneD, NekDouble>(nq,0.0);
	Vmath::Vmul(nq, &outarray[5][0], 1, &m_movingframes[0][k*nq], 1, &Jvec[k][0], 1);
	Vmath::Vvtvp(nq, &outarray[6][0], 1, &m_movingframes[1][k*nq], 1, &Jvec[k][0], 1, &Jvec[k][0], 1);
      }

    ComputeCurl(Jvec, Jcurl);

    Vmath::Vcopy(nq, Jcurl[2], 1, outarray[8], 1);
    Vmath::Vabs(nq, outarray[8], 1, outarray[8], 1);
    */

    /*
    // 9) Graviational field magnitude
    // Grx
    Array<OneD, NekDouble> Grx(nq);
    Vmath::Vmul(nq, &fields[m_Gr1][0], 1, &m_movingframes[0][0], 1, &Grx[0], 1);
    Vmath::Vvtvp(nq, &fields[m_Gr2][0], 1, &m_movingframes[1][0], 1, &Grx[0], 1, &Grx[0], 1);

    // Gry
    Array<OneD, NekDouble> Gry(nq);
    Vmath::Vmul(nq, &fields[m_Gr1][0], 1, &m_movingframes[0][nq], 1, &Gry[0], 1);
    Vmath::Vvtvp(nq, &fields[m_Gr2][0], 1, &m_movingframes[1][nq], 1, &Gry[0], 1, &Gry[0], 1);
    
    Vmath::Vmul(nq, Grx, 1, Grx, 1, outarray[9], 1);
    Vmath::Vvtvp(nq, Gry, 1, Gry, 1, outarray[9], 1, outarray[9], 1);
    Vmath::Vsqrt(nq, outarray[9], 1, outarray[9], 1);
    */
  }


  void MMFEMFHN::SetBoundaryConditions(const int var, NekDouble time)
  { 
      std::string varName;

      // loop over Boundary Regions
      for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
      {	
	varName = m_session->GetVariable(var);
	m_fields[var]->EvaluateBoundaryConditions(time, varName);
      }
  }
  
  /*
  void MMFEMFHN::SetBoundaryConditions(const int var, NekDouble time)
  {
    int cnt = 0;

    // loop over Boundary Regions
    for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
      {	
	m_fields[var]->EvaluateBoundaryConditions(time);
	cnt +=m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
      }
  }
*/

  
  void MMFEMFHN::ComputeAPD(const Array<OneD, const NekDouble> &inarray,
			    Array<OneD, NekDouble> &outarray)
  {
    int i,j;
    int nq         = GetTotPoints();
    
    NekDouble APDpick = 0.2;
    NekDouble Tol = 0.1;
    
    Array<OneD, NekDouble> x(nq), y(nq), z(nq);
    m_fields[0]->GetCoords(x,y,z);
    
    Array<OneD, NekDouble> xnew(nq);
    Array<OneD, NekDouble> phinew(nq);

    int nqnew=0;
    for (i=0; i<nq; ++i)
      {
	if ((y[i]-5.0)<Tol)
	  {
	    xnew[nqnew] = x[i];
	    phinew[nqnew] = inarray[i];
	    nqnew++;
	  }
      }
    
    BubbleSort(xnew,phinew);
    
    int cnt=0;
    i = nqnew-1;
    while (i>0)
      {	
	// Depolarization points
	if(phinew[i] >= APDpick)
	  {	    
	    for (j=i-1; j>0; j--)
	      {		
		if( (phinew[j] < APDpick ) && (fabs(xnew[i]-xnew[j])>10.0) )
		  {
		    outarray[cnt] = fabs(xnew[j]-xnew[i]);
		    break;
		  }
	      }
	    cnt++;
	    i=j-3;
	  }

	else
	  {
	    i--;
	  }
      }
  }

  void MMFEMFHN::Checkpoint_XYZOutput(const int n,
				      const Array<OneD, const Array<OneD, NekDouble> > &fieldphys)
  {
    int nvar = m_fields.num_elements();
    int nq = m_fields[0]->GetTotPoints();
    int ncoeffs = m_fields[0]->GetNcoeffs();

    std::string outname =  m_sessionName +  "XYZ_" + boost::lexical_cast<std::string>(n) + ".chk";
    
    std::vector<std::string> variables(nvar);
    variables[0] = "Ex";
    variables[1] = "Ey";
    variables[2] = "Hz";
    variables[3] = "phi";
    variables[4] = "psi";
    variables[5] = "chi";
    variables[6] = "Fs";
    // variables[7] = "Dx";
    // variables[8] = "Dy";
    // variables[9] = "dFdt";
    
    std::vector<Array<OneD, NekDouble> > fieldcoeffs(nvar);
    std::vector<Array<OneD, NekDouble> > fieldtmp(nvar);
    for(int i=0; i<nvar; ++i)
      {
	fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
	fieldtmp[i] = Array<OneD, NekDouble>(nq);
	Vmath::Vcopy(nq, fieldphys[i], 1, fieldtmp[i], 1);
      }
    
    Array<OneD, NekDouble> tmp1(nq);
    Array<OneD, NekDouble> tmp2(nq);

    // Ex
    Vmath::Vmul(nq, &fieldphys[m_E1][0], 1, &m_movingframes[0][0], 1, &tmp1[0], 1);
    Vmath::Vvtvp(nq, &fieldphys[m_E2][0], 1, &m_movingframes[1][0], 1, &tmp1[0], 1, &tmp1[0], 1);

    // Ey
    Vmath::Vmul(nq, &fieldphys[m_E1][0], 1, &m_movingframes[0][nq], 1, &tmp2[0], 1);
    Vmath::Vvtvp(nq, &fieldphys[m_E2][0], 1, &m_movingframes[1][nq], 1, &tmp2[0], 1, &tmp2[0], 1);

    // Put them into fieldphys
    Vmath::Vcopy(nq, tmp1, 1, fieldtmp[m_E1], 1);
    Vmath::Vcopy(nq, tmp2, 1, fieldtmp[m_E2], 1);


    /*
    // Dx
    Vmath::Vmul(nq, &fieldphys[m_D1][0], 1, &m_movingframes[0][0], 1, &tmp1[0], 1);
    Vmath::Vvtvp(nq, &fieldphys[m_D2][0], 1, &m_movingframes[1][0], 1, &tmp1[0], 1, &tmp1[0], 1);
    
    // Dy
    Vmath::Vmul(nq, &fieldphys[m_D1][0], 1, &m_movingframes[0][nq], 1, &tmp2[0], 1);
    Vmath::Vvtvp(nq, &fieldphys[m_D2][0], 1, &m_movingframes[1][nq], 1, &tmp2[0], 1, &tmp2[0], 1);
    */

    // Put them into fieldphys
    //   Vmath::Vcopy(nq, tmp1, 1, fieldtmp[m_D1], 1);
    //   Vmath::Vcopy(nq, tmp2, 1, fieldtmp[m_D2], 1);
    
    for(int i=0; i<nvar; ++i)
      {
	m_fields[i]->FwdTrans(fieldtmp[i], fieldcoeffs[i]);
      }
    
    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
  }


    
  // 0: |E|
  // 1: |H|
  // 2: |U|
  // 3: Sx
  // 4: Sy
  // 5: Jx
  // 6: Jy
  // 7: DTheta
  // 8: CurlJ
    void MMFEMFHN::Checkpoint_EnergyOutput(const int n,
					   const Array<OneD, const Array<OneD, NekDouble> > &fieldphys)
  {
    int nvar = m_fields.num_elements();
    int ncoeffs = m_fields[0]->GetNcoeffs();
    
    Array<OneD, Array<OneD, NekDouble> >Energy;
    ComputeEnergy(fieldphys,Energy);

    std::string outname =  m_sessionName +  "Energy_" + boost::lexical_cast<std::string>(n) + ".chk";
    
    std::vector<std::string> variables(nvar);
    variables[0] = "E";
    variables[1] = "H";
    variables[2] = "U";
    variables[3] = "Sx";
    variables[4] = "Sy";
    variables[5] = "Jx";
    variables[6] = "Jy";
    // variables[7] = "DTheta";
    // variables[8] = "CurlJ";
    // variables[9] = "Gr";
    
    std::vector<Array<OneD, NekDouble> > fieldcoeffs(nvar);
    for(int i=0; i<nvar; ++i)
      {
	fieldcoeffs[i] = Array<OneD, NekDouble>(ncoeffs);
      }
    
    for(int i=0; i<nvar; ++i)
      {
	m_fields[i]->FwdTrans(Energy[i], fieldcoeffs[i]);
      }
    
    WriteFld(outname, m_fields[0], fieldcoeffs, variables);
  }

  

  void MMFEMFHN::SetUpParametersFentonKarma()
  {
    m_Cm  =  1; // picoF
    m_V0  = -85;
    
    switch (m_FentonKarmaType)
      {
      case eBR:
	m_g_fi_max     = 4;
	m_tau_r        = 33.33;
	m_tau_si       = 29;
	m_tau_0        = 12.5;
	m_tau_v_plus   = 3.33;
	m_tau_v1_minus = 1250;
	m_tau_v2_minus = 19.6;
	m_tau_w_plus   = 870;
	m_tau_w_minus  = 41;
	m_u_c          = 0.13;
	m_u_v          = 0.04;
	m_u_r          = 0.13;
	m_u_fi         = 0.13;
	m_u_csi        = 0.85;
	m_k1           = 10;
	m_k2           = 0.0;
	break;
	
      case eMBR:
	m_g_fi_max     = 4;
	m_tau_r        = 50;
	// m_tau_si       = 44.84;
	m_tau_si       = 45;
	m_tau_0        = 8.3;
	m_tau_v_plus   = 3.33;
	m_tau_v1_minus = 1000;
	m_tau_v2_minus = 19.2;
	m_tau_w_plus   = 667;
	m_tau_w_minus  = 11;
	m_u_c          = 0.13;
	// m_u_v          = 0.04;
	m_u_v          = 0.055;
	m_u_r          = 0.13;
	m_u_fi         = 0.13;
	m_u_csi        = 0.85;
	m_k1           = 10;
	m_k2           = 0.0;
	break;
	
      case eMLR1:
	m_g_fi_max     = 5.8;
	m_tau_r        = 130;
	m_tau_si       = 127;
	m_tau_0        = 12.5;
	m_tau_v_plus   = 10;
	m_tau_v1_minus = 18.2;
	m_tau_v2_minus = 18.2;
	m_tau_w_plus   = 1020;
	m_tau_w_minus  = 80;
	m_u_c          = 0.13;
	m_u_v          = 0.0;
	m_u_r          = 0.13;
	m_u_fi         = 0.13;
	m_u_csi        = 0.85;
	m_k1           = 10;
	m_k2           = 0.0;
	break;
	
      case eGP:
	m_g_fi_max     = 8.7;
	m_tau_r        = 25;
	m_tau_si       = 22.22;
	m_tau_0        = 12.5;
	m_tau_v_plus   = 10;
	m_tau_v1_minus = 333;
	m_tau_v2_minus = 40;
	m_tau_w_plus   = 1000;
	m_tau_w_minus  = 65;
	m_u_c          = 0.13;
	m_u_v          = 0.025;
	m_u_r          = 0.13;
	m_u_fi         = 0.13;
	m_u_csi        = 0.85;
	m_k1           = 10;
	m_k2           = 0.0;
	break;
	
      case eCF1:
	m_g_fi_max     = 6.6666;
	m_tau_r        = 12.5;
	m_tau_si       = 10;
	m_tau_0        = 1.5;
	m_tau_v_plus   = 10;
	m_tau_v1_minus = 350;
	m_tau_v2_minus = 80;
	m_tau_w_plus   = 562;
	m_tau_w_minus  = 48.5;
	m_u_c          = 0.25;
	m_u_v          = 0.001;
	m_u_r          = 0.25;
	m_u_fi         = 0.15;
	m_u_csi        = 0.2;
	m_k1           = 15;
	m_k2           = 0;
	break;
	
      case eCF2a:
	m_g_fi_max     = 6.6666;
	m_tau_r        = 31;
	m_tau_si       = 26.5;
	m_tau_0        = 1.5;
	m_tau_v_plus   = 10;
	m_tau_v1_minus = 20;
	m_tau_v2_minus = 20;
	m_tau_w_plus   = 800;
	m_tau_w_minus  = 45;
	m_u_c          = 0.25;
	m_u_v          = 0.05;
	m_u_r          = 0.6;
	m_u_fi         = 0.11;
	m_u_csi        = 0.7;
	m_k1           = 10;
	m_k2           = 1;
	break;
	
      case eCF2b:
	m_g_fi_max     = 6.6666;
	m_tau_r        = 31;
	m_tau_si       = 26.5;
	m_tau_0        = 1.5;
	m_tau_v_plus   = 10;
	m_tau_v1_minus = 100;
	m_tau_v2_minus = 20;
	m_tau_w_plus   = 800;
	m_tau_w_minus  = 45;
	m_u_c          = 0.25;
	m_u_v          = 0.05;
	m_u_r          = 0.6;
	m_u_fi         = 0.11;
	m_u_csi        = 0.7;
	m_k1           = 10;
	m_k2           = 1;
	break;
	
      case eCF2c:
	m_g_fi_max     = 6.6666;
	m_tau_r        = 31;
	m_tau_si       = 26.5;
	m_tau_0        = 1.5;
	m_tau_v_plus   = 10;
	m_tau_v1_minus = 150;
	m_tau_v2_minus = 20;
	m_tau_w_plus   = 800;
	m_tau_w_minus  = 45;
	m_u_c          = 0.25;
	m_u_v          = 0.05;
	m_u_r          = 0.6;
	m_u_fi         = 0.11;
	m_u_csi        = 0.7;
	m_k1           = 10;
	m_k2           = 1;
	break;
	
      case eCF3a:
	m_g_fi_max     = 13.3333;
	m_tau_r        = 38;
	m_tau_si       = 127;
	m_tau_0        = 8.3;
	m_tau_v_plus   = 3.33;
	m_tau_v1_minus = 45;
	m_tau_v2_minus = 300;
	m_tau_w_plus   = 600;
	m_tau_w_minus  = 40;
	m_tau_y_plus   = 1000;
	m_tau_y_minus  = 230;
	m_u_c          = 0.25;
	m_u_v          = 0.5;
	m_u_r          = 0.25;
	m_u_fi         = 0.25;
	m_u_csi        = 0.7;
	m_k1           = 60;
	m_k2           = 0;
	break;
	
      case eCF3b:
	m_g_fi_max     = 13.3333;
	m_tau_r        = 38;
	m_tau_si       = 127;
	m_tau_0        = 8.3;
	m_tau_v_plus   = 3.33;
	m_tau_v1_minus = 20;
	m_tau_v2_minus = 300;
	m_tau_w_plus   = 600;
	m_tau_w_minus  = 40;
	m_tau_y_plus   = 1000;
	m_tau_y_minus  = 230;
	m_u_c          = 0.25;
	m_u_v          = 0.5;
	m_u_r          = 0.25;
	m_u_fi         = 0.25;
	m_u_csi        = 0.7;
	m_k1           = 60;
	m_k2           = 0;
	break;

      default:
	break;
      }

     m_tau_d  = m_Cm/m_g_fi_max;

     /*
     Gs::cout << "m_g_fi_max = " << m_g_fi_max << ", tau_r = " << m_tau_r << ", tau_si = " << m_tau_si
	  << ", tau_0 = " << m_tau_0 << ", tau_v_plus = " << m_tau_v_plus << ", tau_v1_minus = " << m_tau_v1_minus
	  << ", tau_v2_minus = " << m_tau_v2_minus << ", tau_w_plus = " << m_tau_w_plus << ", tau_w_minus = " << m_tau_w_minus
	  << ", u_c = " << m_u_c << ", u_v = " << m_u_v << ", u_fi = " << m_u_fi << ", u_csi = " << m_u_csi
	  << ", k1 = " << m_k1 << ", k2 = " << m_k2 << Gs::endl << Gs::endl;
     */
  }
  
    
  void MMFEMFHN::v_GenerateSummary(SolverUtils::SummaryList& s)
  {
    MMFSystem::v_GenerateSummary(s);
    SolverUtils::AddSummaryItem(s, "TestMaxwellType", SolverUtils::TestMaxwellTypeMap[m_TestMaxwellType]);
    SolverUtils::AddSummaryItem(s, "PolType", SolverUtils::PolTypeMap[m_PolType]);
    SolverUtils::AddSummaryItem(s, "InitWaveType", InitWaveTypeMap[m_InitWaveType]);
    SolverUtils::AddSummaryItem(s, "Radius of Init", m_radiusofinit);
    SolverUtils::AddSummaryItem(s, "PoissonTau", m_PoissonTau);
    SolverUtils::AddSummaryItem(s, "kp", m_kp);

    if(fabs(m_sigma[0]-1.0)>0.000001)
      {
	SolverUtils::AddSummaryItem(s, "Sigma1", m_sigma[0]);
      }

    if(fabs(m_sigma[1]-1.0)>0.000001)
      {
	SolverUtils::AddSummaryItem(s, "Sigma2", m_sigma[1]);
      }
      
    if(m_TestType==eFentonKarma)
      {
	SolverUtils::AddSummaryItem(s, "FentonKarmaType", FentonKarmaTypeMap[m_FentonKarmaType]);
      }

    if(m_StimulusPeriod>0)
      {
	SolverUtils::AddSummaryItem(s, "StimulusTime", m_StimulusPeriod*m_timestep);
      }
  }
  
}
