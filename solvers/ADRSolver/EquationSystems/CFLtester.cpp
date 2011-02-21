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
// Description: Unsteady Advection  solve
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <ADRSolver/EquationSystems/CFLtester.h>

namespace Nektar
{
    string CFLtester::className = EquationSystemFactory::RegisterCreatorFunction("CFLtester", CFLtester::create, "Testing CFL restriction");

    CFLtester::CFLtester(SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession)
    {
        // Define Velocity fields
        m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        int nq = m_fields[0]->GetNpoints();
        std::string velStr[3] = {"Vx","Vy","Vz"};

        for(int i = 0; i < m_spacedim; ++i)
        {
            m_velocity[i] = Array<OneD, NekDouble> (nq,0.0);

            SpatialDomains::ConstUserDefinedEqnShPtr ifunc
                = m_boundaryConditions->GetUserDefinedEqn(velStr[i]);

            EvaluateFunction(m_velocity[i],ifunc);
        }

        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs        (&CFLtester::DoOdeRhs,        this);
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
        int i;
        int nvariables = inarray.num_elements();
        int npoints = GetNpoints();

        switch (m_projectionType)
        {
        case eDiscontinuousGalerkin:
            {
                int ncoeffs    = inarray[0].num_elements();
                Array<OneD, Array<OneD, NekDouble> > WeakAdv(nvariables);

                WeakAdv[0] = Array<OneD, NekDouble>(ncoeffs*nvariables);
                for(i = 1; i < nvariables; ++i)
                {
                    WeakAdv[i] = WeakAdv[i-1] + ncoeffs;
                }

                WeakDGAdvection(inarray, WeakAdv,true,true);

                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->MultiplyByElmtInvMass(WeakAdv[i],
                                                       WeakAdv[i]);
                    m_fields[i]->BwdTrans(WeakAdv[i],outarray[i]);
                    Vmath::Neg(npoints,outarray[i],1);
                }

                break;
            }
            case eGalerkin:
            {
                // Calculate -V\cdot Grad(u);
                for(i = 0; i < nvariables; ++i)
                {
                    AdvectionNonConservativeForm(m_velocity,
                                                 inarray[i],
                                                 outarray[i]);
                    Vmath::Neg(npoints,outarray[i],1);
                }
                break;
            }
        }
    }



    /**
     *
     */
    void CFLtester::DoOdeProjection(const Array<OneD,
                                            const Array<OneD, NekDouble> >&inarray,
                                            Array<OneD,       Array<OneD, NekDouble> >&outarray,
                                            const NekDouble time)
    {
        int i;
        int nvariables = inarray.num_elements();
        SetBoundaryConditions(time);

        switch(m_projectionType)
        {
        case eDiscontinuousGalerkin:
            {
                // Just copy over array
                int npoints = GetNpoints();

                for(i = 0; i < nvariables; ++i)
                {
                    Vmath::Vcopy(npoints,inarray[i],1,outarray[i],1);
                }
            }
            break;
        case eGalerkin:
            {
                Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());

                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->FwdTrans(inarray[i],coeffs,false);
                    m_fields[i]->BwdTrans_IterPerExp(coeffs,outarray[i]);
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
        int nvel = m_spacedim; //m_velocity.num_elements();

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
		const NekDouble minLengthStdTri  = 0.7072*0.5;
		const NekDouble minLengthStdQuad = 0.5;
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
				//tstep[el] =  timeCFL*minLengthStdTri/(stdVelocity[el]*cLambda*(ExpOrder[el]-1)*(ExpOrder[el]-1));
				tstep[el] = CFL[el]*minLengthStdTri/(stdVelocity[el]);
			}
			else if(boost::dynamic_pointer_cast<LocalRegions::QuadExp>(m_fields[0]->GetExp(el)))
			{ 
				//tstep[el] =  timeCFL*minLengthStdQuad/(stdVelocity[el]*cLambda*(ExpOrder[el]-1)*(ExpOrder[el]-1));
				tstep[el] = CFL[el]*minLengthStdQuad/(stdVelocity[el]);
			}
		}
		
		NekDouble TimeStep = Vmath::Vmin(n_element,tstep,1);
		
		return TimeStep;
	}


}
