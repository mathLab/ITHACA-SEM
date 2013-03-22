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

#include <ADRSolver/EquationSystems/CFLtester.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>

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
        std::vector<std::string> vel;
        vel.push_back("Vx");
        vel.push_back("Vy");
        vel.push_back("Vz");
        vel.resize(m_spacedim);

        EvaluateFunction(vel, m_velocity, "AdvectionVelocity");

        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs     (&CFLtester::DoOdeRhs,        this);
            m_ode.DefineProjection (&CFLtester::DoOdeProjection, this);
        }
        else
        {
            ASSERTL0(false, "Implicit unsteady Advection not set up.");
        }
    }

    CFLtester::~CFLtester()
    {
    }

    void CFLtester::DoOdeRhs(
        const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
              Array<OneD,        Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
    {
        int j;
        int nvariables = inarray.num_elements();
        int npoints = GetNpoints();

        switch (m_projectionType)
        {
            case MultiRegions::eDiscontinuous:
            {
                int ncoeffs    = inarray[0].num_elements();
                Array<OneD, Array<OneD, NekDouble> > WeakAdv(nvariables);

                WeakAdv[0] = Array<OneD, NekDouble>(ncoeffs*nvariables);
                for(j = 1; j < nvariables; ++j)
                {
                    WeakAdv[j] = WeakAdv[j-1] + ncoeffs;
                }

                WeakDGAdvection(inarray, WeakAdv, true, true);

                for(j = 0; j < nvariables; ++j)
                {
                    m_fields[j]->MultiplyByElmtInvMass(WeakAdv[j], WeakAdv[j]);
                    m_fields[j]->BwdTrans(WeakAdv[j],outarray[j]);
                    Vmath::Neg(npoints,outarray[j],1);
                }
                break;
            }
            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                // Calculate -V\cdot Grad(u);
                for(j = 0; j < nvariables; ++j)
                {
                    AdvectionNonConservativeForm(m_velocity, 
                                                 inarray[j],
                                                 outarray[j]);
                    
                    Vmath::Neg(npoints, outarray[j], 1);
                }
                break;
            }
        }
    }

    void CFLtester::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> >&inarray,
              Array<OneD,       Array<OneD, NekDouble> >&outarray,
    const NekDouble time)
    {
        int j;
        int nvariables = inarray.num_elements();
        SetBoundaryConditions(time);

        switch(m_projectionType)
        {
            case MultiRegions::eDiscontinuous:
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
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());

                for(j = 0; j < nvariables; ++j)
                {
                    m_fields[j]->FwdTrans(inarray[j],coeffs);
                    m_fields[j]->BwdTrans_IterPerExp(coeffs,outarray[j]);
                }
                break;
            }
            default:
                ASSERTL0(false,"Unknown projection scheme");
                break;
        }
    }


    void CFLtester::v_GetFluxVector(
        const int i, 
        Array<OneD, Array<OneD, NekDouble> > &physfield,
        Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        ASSERTL1(flux.num_elements() == m_velocity.num_elements(),
                 "Dimension of flux array and velocity array do not match");

        for (int j = 0; j < flux.num_elements(); ++j)
        {
            Vmath::Vmul(GetNpoints(),
                        physfield[i], 1,
                        m_velocity[j], 1,
                        flux[j], 1);
        }
    }

    void CFLtester::v_NumericalFlux(
        Array<OneD, Array<OneD, NekDouble> > &physfield, 
        Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        int i;
        int nTraceNumPoints = GetTraceNpoints();
        int nvel = m_spacedim;

        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
        Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);

        // Get Edge Velocity - Could be stored if time independent
        for (i = 0; i < nvel; ++i)
        {
            m_fields[0]->ExtractTracePhys(m_velocity[i], Fwd);
            Vmath::Vvtvp(nTraceNumPoints,
                         m_traceNormals[i], 1,
                         Fwd, 1, Vn, 1, Vn, 1);
        }

        for (i = 0; i < numflux.num_elements(); ++i)
        {
            m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd, Bwd);
            //evaulate upwinded m_fields[i]
            m_fields[i]->GetTrace()->Upwind(Vn, Fwd, Bwd, numflux[i]);
            // calculate m_fields[i]*Vn
            Vmath::Vmul(nTraceNumPoints, numflux[i], 1, Vn, 1, numflux[i], 1);
        }
    }

    
    
    void CFLtester::v_PrintSummary(std::ostream &out)
    {
        UnsteadySystem::v_PrintSummary(out);
    }
	
    
    
	NekDouble CFLtester::v_GetTimeStep(
        const Array<OneD,int> ExpOrder, 
        const Array<OneD,NekDouble> CFL, 
        NekDouble timeCFL)
	{ 
		
		int n_element       = m_fields[0]->GetExpSize();
        
		//const NekDouble minLengthStdTri  = 1.414213;
		//const NekDouble minLengthStdQuad = 2.0;
		//const NekDouble cLambda          = 0.2; // Spencer book pag. 317

        Array<OneD, NekDouble> tstep      (n_element, 0.0);
		Array<OneD, NekDouble> stdVelocity(n_element, 0.0);
		stdVelocity = GetStdVelocity(m_velocity);
		
		for (int el = 0; el < n_element; ++el)
		{
			int npoints = m_fields[0]->GetExp(el)->GetTotPoints();
			Array<OneD, NekDouble> one2D(npoints, 1.0);
			//NekDouble Area = m_fields[0]->GetExp(el)->Integral(one2D);
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
		
		NekDouble TimeStep = Vmath::Vmin(n_element, tstep, 1);
		
		return TimeStep;
	}
	
    
    
	NekDouble CFLtester::v_GetTimeStep(
        int ExpOrder, 
        NekDouble CFL, 
        NekDouble TimeStability)
	{
		//================================================================
		// This function has been created just to test specific problems, hence is not general
		// and it has been implemented in a rude fashion, as the full CFLtester class.
		// For real CFL calculations refer to the general implementation above. (A.Bolis)
		//================================================================
		
		NekDouble TimeStep;
		NekDouble SpatialStability;
		int n_elements = m_fields[0]->GetExpSize();
		
		//solve ambiguity in windows
		NekDouble n_elem = n_elements;
		NekDouble DH     = sqrt(n_elem);
		
		int H = (int)DH;
		int P = ExpOrder-1;
		
		//================================================================
		// Regular meshes
		
		//SpatialStability = EigenvaluesRegMeshes[H-1][P-1];
		
		//================================================================
		// Anisotropic meshes
		
		if (TimeStability == 1.0) 
		{
			SpatialStability = EigenvaluesAnaMeshesAB2[H/2][P-1];
		}
		else if (TimeStability == 2.0) 
		{
			SpatialStability = EigenvaluesAnaMeshesRK2[H/2][P-1];
		}
		else if (TimeStability == 2.784) 
		{
			SpatialStability = EigenvaluesAnaMeshesRK4[H/2][P-1];
		}
		else 
		{
			ASSERTL0(false,"Dominant eigenvalues database not present for this time-stepping scheme")
		}
		
		//================================================================
		
		TimeStep = (TimeStability/SpatialStability)*CFL;
		
		//================================================================
		
		return TimeStep;
	}
    
    
    
    Array<OneD, NekDouble> CFLtester::GetStdVelocity(
        const Array<OneD, Array<OneD,NekDouble> > inarray)
    {
        // Checking if the problem is 2D
        ASSERTL0(m_expdim >= 2, "Method not implemented for 1D");
        
        int nTotQuadPoints  = GetTotPoints();
        int n_element       = m_fields[0]->GetExpSize();       
        int nvel            = inarray.num_elements();
        
        NekDouble pntVelocity;
        
        // Getting the standard velocity vector on the 2D normal space
        Array<OneD, Array<OneD, NekDouble> > stdVelocity(nvel);
        Array<OneD, NekDouble> stdV(n_element, 0.0);
        
        for (int i = 0; i < nvel; ++i)
        {
            stdVelocity[i] = Array<OneD, NekDouble>(nTotQuadPoints);
        }
		
        if (nvel == 2)
        {
            for (int el = 0; el < n_element; ++el)
            { 
                
                int n_points = m_fields[0]->GetExp(el)->GetTotPoints();
                
                Array<OneD, const NekDouble> jac  = 
                m_fields[0]->GetExp(el)->GetGeom2D()->GetJac();
                Array<TwoD, const NekDouble> gmat = 
                m_fields[0]->GetExp(el)->GetGeom2D()->GetGmat();
                
                if (m_fields[0]->GetExp(el)->GetGeom2D()->GetGtype() 
                    == SpatialDomains::eDeformed)
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][i]*inarray[0][i] 
                        + gmat[2][i]*inarray[1][i];
                        
                        stdVelocity[1][i] = gmat[1][i]*inarray[0][i] 
                        + gmat[3][i]*inarray[1][i];
                    }
                }
                else
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][0]*inarray[0][i] 
                        + gmat[2][0]*inarray[1][i];
                        
                        stdVelocity[1][i] = gmat[1][0]*inarray[0][i] 
                        + gmat[3][0]*inarray[1][i];
                    }
                }
                
                
                for (int i = 0; i < n_points; i++)
                {
                    pntVelocity = sqrt(stdVelocity[0][i]*stdVelocity[0][i] 
                                       + stdVelocity[1][i]*stdVelocity[1][i]);
                    
                    if (pntVelocity>stdV[el])
                    {
                        stdV[el] = pntVelocity;
                    }
                }
            }
        }
        else
        {
            for (int el = 0; el < n_element; ++el)
            { 
                
                int n_points = m_fields[0]->GetExp(el)->GetTotPoints();
                
                Array<OneD, const NekDouble> jac =
                m_fields[0]->GetExp(el)->GetGeom3D()->GetJac();
                Array<TwoD, const NekDouble> gmat =
                m_fields[0]->GetExp(el)->GetGeom3D()->GetGmat();
                
                if (m_fields[0]->GetExp(el)->GetGeom3D()->GetGtype() 
                    == SpatialDomains::eDeformed)
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][i]*inarray[0][i] 
                        + gmat[3][i]*inarray[1][i] 
                        + gmat[6][i]*inarray[2][i];
                        
                        stdVelocity[1][i] = gmat[1][i]*inarray[0][i] 
                        + gmat[4][i]*inarray[1][i] 
                        + gmat[7][i]*inarray[2][i];
                        
                        stdVelocity[2][i] = gmat[2][i]*inarray[0][i] 
                        + gmat[5][i]*inarray[1][i] 
                        + gmat[8][i]*inarray[2][i];
                    }
                }
                else
                {
                    Array<OneD, const NekDouble> jac =
                    m_fields[0]->GetExp(el)->GetGeom3D()->GetJac();
                    Array<TwoD, const NekDouble> gmat = 
                    m_fields[0]->GetExp(el)->GetGeom3D()->GetGmat();
                    
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][0]*inarray[0][i] 
                        + gmat[3][0]*inarray[1][i] 
                        + gmat[6][0]*inarray[2][i];
                        
                        stdVelocity[1][i] = gmat[1][0]*inarray[0][i] 
                        + gmat[4][0]*inarray[1][i] 
                        + gmat[7][0]*inarray[2][i];
                        
                        stdVelocity[2][i] = gmat[2][0]*inarray[0][i] 
                        + gmat[5][0]*inarray[1][i] 
                        + gmat[8][0]*inarray[2][i];
                    }
                }
                
                for (int i = 0; i < n_points; i++)
                {
                    pntVelocity = sqrt(stdVelocity[0][i]*stdVelocity[0][i] 
                                       + stdVelocity[1][i]*stdVelocity[1][i] 
                                       + stdVelocity[2][i]*stdVelocity[2][i]);
                    
                    if (pntVelocity > stdV[el])
                    {
                        stdV[el] = pntVelocity;
                    }
                }
            }
        }
		
        return stdV;
	}
}
