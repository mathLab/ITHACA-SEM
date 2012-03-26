///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadyAdvection.cpp
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
// Description: CFL tester solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <complex>

#include <ADRSolver/EquationSystems/CFLtester.h>

namespace Nektar
{
    string CFLtester::className = GetEquationSystemFactory().RegisterCreatorFunction("CFLtester", CFLtester::create, "Testing CFL restriction");

    CFLtester::CFLtester(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession)
    {
    }

    void CFLtester::v_InitObject()
    {
        UnsteadySystem::v_InitObject();

        m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        int nq = m_fields[0]->GetNpoints();
        std::vector<std::string> vel;
        vel.push_back("Vx");
        vel.push_back("Vy");
        vel.push_back("Vz");
        vel.resize(m_spacedim);

        EvaluateFunction(vel, m_velocity, "AdvectionVelocity");

        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs        (&CFLtester::DoOdeRhs,        this);
            m_ode.DefineProjection (&CFLtester::DoOdeProjection, this);
        }
        else
        {
            ASSERTL0(false, "Implicit unsteady Advection not set up.");
        }
		
		//m_max_time_step = CalculateMaximumTimeStep();
    }

    CFLtester::~CFLtester()
    {
    }

    void CFLtester::DoOdeRhs(const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
							Array<OneD,        Array<OneD, NekDouble> >&outarray,
							const NekDouble time)
    {
        int j;
        int nvariables = inarray.num_elements();
        int npoints = GetNpoints();

        switch (m_projectionType)
        {
        case MultiRegions::eDiscontinuousGalerkin:
            {
                int ncoeffs    = inarray[0].num_elements();
                Array<OneD, Array<OneD, NekDouble> > WeakAdv(nvariables);

                WeakAdv[0] = Array<OneD, NekDouble>(ncoeffs*nvariables);
                for(j = 1; j < nvariables; ++j)
                {
                    WeakAdv[j] = WeakAdv[j-1] + ncoeffs;
                }

                WeakDGAdvection(inarray, WeakAdv,true,true);

                for(j = 0; j < nvariables; ++j)
                {
                    m_fields[j]->MultiplyByElmtInvMass(WeakAdv[j],
                                                       WeakAdv[j]);
                    m_fields[j]->BwdTrans(WeakAdv[j],outarray[j]);
                    Vmath::Neg(npoints,outarray[j],1);
                }

                break;
            }
            case MultiRegions::eGalerkin:
            {
                // Calculate -V\cdot Grad(u);
                for(j = 0; j < nvariables; ++j)
                {
                    AdvectionNonConservativeForm(m_velocity,
                                                 inarray[j],
                                                 outarray[j]);
                    Vmath::Neg(npoints,outarray[j],1);
                }
                break;
            }
        }
    }

    void CFLtester::DoOdeProjection(const Array<OneD,
									const Array<OneD, NekDouble> >&inarray,
									Array<OneD,       Array<OneD, NekDouble> >&outarray,
									const NekDouble time)
    {
        int j;
        int nvariables = inarray.num_elements();
        SetBoundaryConditions(time);

        switch(m_projectionType)
        {
        case MultiRegions::eDiscontinuousGalerkin:
            {
                // Just copy over array
                int npoints = GetNpoints();

                for(j = 0; j < nvariables; ++j)
                {
                    Vmath::Vcopy(npoints,inarray[j],1,outarray[j],1);
                }
            }
            break;
        case MultiRegions::eGalerkin:
            {
                Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());

                for(j = 0; j < nvariables; ++j)
                {
                    m_fields[j]->FwdTrans(inarray[j],coeffs,false);
                    m_fields[j]->BwdTrans_IterPerExp(coeffs,outarray[j]);
                }
                break;
            }
        default:
            ASSERTL0(false,"Unknown projection scheme");
            break;
        }
    }


    void CFLtester::v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        ASSERTL1(flux.num_elements() == m_velocity.num_elements(),"Dimension of flux array and velocity array do not match");

        for(int j = 0; j < flux.num_elements(); ++j)
        {
            Vmath::Vmul(GetNpoints(),physfield[i],1,
                m_velocity[j],1,flux[j],1);
        }
    }

    void CFLtester::v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        int i;

        int nTraceNumPoints = GetTraceNpoints();
        int nvel = m_spacedim;

        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
        Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);

        // Get Edge Velocity - Could be stored if time independent
        for(i = 0; i < nvel; ++i)
        {
            m_fields[0]->ExtractTracePhys(m_velocity[i], Fwd);
            Vmath::Vvtvp(nTraceNumPoints,m_traceNormals[i],1,Fwd,1,Vn,1,Vn,1);
        }

        for(i = 0; i < numflux.num_elements(); ++i)
        {
            m_fields[i]->GetFwdBwdTracePhys(physfield[i],Fwd,Bwd);
            //evaulate upwinded m_fields[i]
            m_fields[i]->GetTrace()->Upwind(Vn,Fwd,Bwd,numflux[i]);
            // calculate m_fields[i]*Vn
            Vmath::Vmul(nTraceNumPoints,numflux[i],1,Vn,1,numflux[i],1);
        }
    }

    void CFLtester::v_PrintSummary(std::ostream &out)
    {
        UnsteadySystem::v_PrintSummary(out);
    }
	
	NekDouble CFLtester::v_GetTimeStep(const Array<OneD,int> ExpOrder, const Array<OneD,NekDouble> CFL, NekDouble timeCFL)
	{ 
		
		int nvariables = m_fields.num_elements();               // Number of variables in the mesh
		int nTotQuadPoints  = GetTotPoints();
		int n_element  = m_fields[0]->GetExpSize(); 
		Array<OneD, NekDouble> tstep(n_element,0.0);
		const NekDouble minLengthStdTri  = 1.414213;
		const NekDouble minLengthStdQuad = 2.0;
		const NekDouble cLambda = 0.2; // Spencer book pag. 317
		Array<OneD, NekDouble> stdVelocity(n_element,0.0);
		stdVelocity = GetStdVelocity(m_velocity);
		
		for(int el = 0; el < n_element; ++el)
		{
			int npoints = m_fields[0]->GetExp(el)->GetTotPoints();
			Array<OneD, NekDouble> one2D(npoints, 1.0);
			NekDouble Area = m_fields[0]->GetExp(el)->Integral(one2D);
			if(boost::dynamic_pointer_cast<LocalRegions::TriExp>(m_fields[0]->GetExp(el)))
			{
				//tstep[el] =  timeCFL/(stdVelocity[el]*cLambda*(ExpOrder[el]-1)*(ExpOrder[el]-1));
				//tstep[el] =  timeCFL*minLengthStdTri/(stdVelocity[el]*cLambda*(ExpOrder[el]-1)*(ExpOrder[el]-1));
				tstep[el] = CFL[el]/(stdVelocity[el]);
			}
			else if(boost::dynamic_pointer_cast<LocalRegions::QuadExp>(m_fields[0]->GetExp(el)))
			{ 
				//tstep[el] =  timeCFL/(stdVelocity[el]*cLambda*(ExpOrder[el]-1)*(ExpOrder[el]-1));
				//tstep[el] =  timeCFL*minLengthStdQuad/(stdVelocity[el]*cLambda*(ExpOrder[el]-1)*(ExpOrder[el]-1));
				tstep[el] = CFL[el]/(stdVelocity[el]);
			}
		}
		
		NekDouble TimeStep = Vmath::Vmin(n_element,tstep,1);
		
		return TimeStep;
	}
	
	NekDouble CFLtester::v_GetTimeStep(NekDouble CFL)
	{
		return CFL*m_max_time_step;
	}
	
	NekDouble CFLtester::CalculateMaximumTimeStep()
	{
		NekDouble max_time_step;
		int j;
		Array<OneD, Array<OneD, NekDouble> > inarray(1);
		Array<OneD, Array<OneD, NekDouble> > tmp(1);
		Array<OneD, Array<OneD, NekDouble> > outarray(1);
		Array<OneD, Array<OneD, NekDouble> > WeakAdv(1);
		int npoints = GetNpoints();
		int ncoeffs = GetNcoeffs();
		inarray[0]  = Array<OneD, NekDouble>(npoints,0.0);
		outarray[0] = Array<OneD, NekDouble>(npoints,0.0);
		tmp[0] = Array<OneD, NekDouble>(npoints,0.0);		
		WeakAdv[0]  = Array<OneD, NekDouble>(ncoeffs,0.0);
		Array<OneD, NekDouble> MATRIX(npoints*npoints,0.0);
		
		for (j = 0; j < npoints; j++)
		{
			inarray[0][j] = 1.0;
			switch (m_projectionType)
			{
				case MultiRegions::eDiscontinuousGalerkin:
				{
					WeakDGAdvection(inarray, WeakAdv,true,true,1);
					m_fields[0]->MultiplyByElmtInvMass(WeakAdv[0],WeakAdv[0]);
					m_fields[0]->BwdTrans(WeakAdv[0],outarray[0]);
					Vmath::Neg(npoints,outarray[0],1);
					break;
				}
				case MultiRegions::eGalerkin:
				{
					m_fields[0]->FwdTrans(inarray[0],WeakAdv[0]);
					m_fields[0]->BwdTrans_IterPerExp(WeakAdv[0],tmp[0]);
					AdvectionNonConservativeForm(m_velocity,tmp[0],outarray[0]);
					Vmath::Neg(npoints,outarray[0],1);
					m_fields[0]->FwdTrans(outarray[0],WeakAdv[0]);
					m_fields[0]->BwdTrans_IterPerExp(WeakAdv[0],outarray[0]);
				    break;
				}
			}

			Vmath::Vcopy(npoints,&(outarray[0][0]),1,&(MATRIX[j]),npoints);
			inarray[0][j] = 0.0;
		}
		
		char jobvl = 'N';
		char jobvr = 'N';
		int info = 0, lwork = 3*npoints;
		NekDouble dum;
		Array<OneD, NekDouble> EIG_R(npoints);
		Array<OneD, NekDouble> EIG_I(npoints);
		Array<OneD, NekDouble> work(lwork);
		Lapack::Dgeev(jobvl,jobvr,npoints,MATRIX.get(),npoints,EIG_R.get(),EIG_I.get(),&dum,1,&dum,1,&work[0],lwork,info);
				
		LibUtilities::TimeIntegrationMethod  timeIntMethod;
		
		// Determine TimeIntegrationMethod to use.
        ASSERTL0(m_session->DefinesSolverInfo("TIMEINTEGRATIONMETHOD"),
				 "No TIMEINTEGRATIONMETHOD defined in session.");
        for (j = 0; j < (int)LibUtilities::SIZE_TimeIntegrationMethod; ++j)
        {
            bool match;
            m_session->MatchSolverInfo("TIMEINTEGRATIONMETHOD",
									   LibUtilities::TimeIntegrationMethodMap[j], match, false);
            if (match)
            {
                timeIntMethod = (LibUtilities::TimeIntegrationMethod) j;
                break;
            }
        }
        ASSERTL0(j != (int) LibUtilities::SIZE_TimeIntegrationMethod,
				 "Invalid time integration type.");
		
		NekDouble phase,module;
		NekDouble tol = 0.00001;
		int nteta = 10000;
		NekDouble dteta = 2*M_PI/(nteta-1);
		Array<OneD, NekDouble> time_steps(npoints,1.0);
		Array<OneD, NekDouble> teta(nteta,0.0);
		
		for(j = 1; j < nteta; j++)
		{
			teta[j] = teta[j-1] + dteta;
		}
		
		Array<OneD, complex<NekDouble> > z (nteta);
		Array<OneD, complex<NekDouble> > r (nteta);
		Array<OneD, complex<NekDouble> > s (nteta);
		Array<OneD, complex<NekDouble> > w (nteta);
		
		for(j = 0; j < nteta; j++)
		{
			complex<NekDouble> iteta(0.0,teta[j]);
			z[j] = exp(iteta);
		}
		
		switch(m_timeIntMethod)
        {
			case LibUtilities::eAdamsBashforthOrder2:
			{
				if(m_projectionType == MultiRegions::eDiscontinuousGalerkin)
				{
					for(j = 0; j < nteta; j++)
					{
						r[j] = z[j]-1.0;
						s[j] = 1.5-0.5/z[j];
						w[j] = r[j]/s[j];
					}
					
					for(j = 0; j < npoints; j++)
					{
						phase  = atan2(EIG_I[j],EIG_R[j]);
						module = sqrt(EIG_R[j]*EIG_R[j] + EIG_I[j]*EIG_I[j]);
						
						for(int k = 0; k < nteta; k++)
						{
							if(abs((arg(w[k])-phase))<=tol)
							{
								time_steps[j] = abs(w[k])/module;
								k = nteta+1;
							}
							if(k == nteta-1)
							{
								tol = 10*tol;
								k = 0;
							}
						}
					}
				}
				else 
				{
					ASSERTL0(false,"AB2 is theoretically unstable using a CG projection");
				}
				break;
			}
			case LibUtilities::eRungeKutta2_ImprovedEuler:
			{
				if(m_projectionType == MultiRegions::eDiscontinuousGalerkin)
				{
					for(j = 0; j < nteta; j++)
					{
						w[j] = z[j]-1.0;
						for(int k = 0; k < 3; k++)
						{
							w[j] = w[j] - (1.0 + w[j] + (w[j]*w[j])/2.0 - z[j]*z[j])/(1.0+w[j]);
						}
					}
					for(j = 0; j < npoints; j++)
					{
						phase  = atan2(EIG_I[j],EIG_R[j]);
						module = sqrt(EIG_R[j]*EIG_R[j] + EIG_I[j]*EIG_I[j]);
						
						for(int k = 0; k < nteta; k++)
						{
							if(abs((arg(w[k])-phase))<=tol)
							{
								time_steps[j] = abs(w[k])/module;
								k = nteta+1;
							}
							if(k == nteta-1)
							{
								tol = 10*tol;
								k = 0;
							}
						}
						
					}
				}
				else 
				{
					ASSERTL0(false,"RK2 is theoretically unstable using a CG projection");
				}
				break;
			}
			case LibUtilities::eClassicalRungeKutta4:
			{
				if(m_projectionType == MultiRegions::eDiscontinuousGalerkin)
				{
					for(j = 0; j < nteta; j++)
					{
						w[j] = z[j]-1.0;
						for(int k = 0; k < 3; k++)
						{
							w[j] = w[j] - (1.0 + w[j] + (w[j]*w[j])/2.0 - z[j]*z[j])/(1.0+w[j]);
						}
						for(int k = 0; k < 4; k++)
						{
							w[j] = w[j] - (1.0 + w[j] + (w[j]*w[j])/2.0 + (w[j]*w[j]*w[j])/6.0 - z[j]*z[j]*z[j])/(1.0 + w[j] + (w[j]*w[j])/2.0);
						}
						
						for(int k = 0; k < 4; k++)
						{
							w[j] = w[j] - (1.0 + w[j] + (w[j]*w[j])/2.0 + (w[j]*w[j]*w[j])/6.0 + (w[j]*w[j]*w[j]*w[j])/24.0 - z[j]*z[j]*z[j]*z[j])/(1.0 + w[j] + (w[j]*w[j])/2.0 + (w[j]*w[j]*w[j])/6.0);
						}
					}
					for(j = 0; j < npoints; j++)
					{
						phase  = atan2(EIG_I[j],EIG_R[j]);
						module = sqrt(EIG_R[j]*EIG_R[j] + EIG_I[j]*EIG_I[j]);
						
						for(int k = 0; k < nteta; k++)
						{
							if(abs((arg(w[k])-phase))<=tol)
							{
								time_steps[j] = abs(w[k])/module;
								k = nteta+1;
							}
							if(k == nteta-1)
							{
								tol = 10*tol;
								k = 0;
							}
						}
					}
				}
				else 
				{
					NekDouble max_imag = Vmath::Vmax(npoints,EIG_I,1);
					
					time_steps[0] = 2.828/max_imag;
				}
				break;
			}
			default:
            {
                ASSERTL0(false,"Exact CFL limit not implemented for this time-integration scheme");
            }
		}
		
		max_time_step = Vmath::Vmin(npoints,time_steps,1);
		
		cout << "max_time_step = " <<  max_time_step << endl;
		
		return max_time_step;
	}
	
	
}
