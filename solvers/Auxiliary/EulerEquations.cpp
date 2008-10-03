///////////////////////////////////////////////////////////////////////////////
//
// File EulerEquations.cpp
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
// Description: Two-dimensional Euler equations 
//
///////////////////////////////////////////////////////////////////////////////

#include <../solvers/Auxiliary/EulerEquations.h>

namespace Nektar
{
    EulerEquations::EulerEquations(void):
        AdvectionDiffusionReaction()
    {
    }
    
    EulerEquations::EulerEquations(SpatialDomains::MeshGraph2D &graph2D,
				   SpatialDomains::BoundaryConditions &bcs):
        AdvectionDiffusionReaction(graph2D,bcs,4)
    {
    }


    //  void EulerEquations::ConservativeToPrimitive(void)
    //     {
      

    //     }

    //     void EulerEquations::PrimitiveToConservative(void)
    //     {


    //     }


    void EulerEquations::GetFluxVector(Array<OneD, Array<OneD, NekDouble> >&FX,
				       Array<OneD, Array<OneD, NekDouble> >&FY)
    { 
      
        NekDouble gamma = GetGamma();

        // compute the pressure
        Array<OneD, NekDouble> p( m_fields[0]->GetPointsTot() );
      
        for (int i = 0; i < m_fields[0]->GetPointsTot(); ++i)
	{
            p[i] = (gamma - 1.0)*(m_fields[3]->GetPhys()[i] - 
                                  0.5*(m_fields[1]->GetPhys()[i]*m_fields[1]->GetPhys()[i]/ m_fields[0]->GetPhys()[i] +
                                       m_fields[2]->GetPhys()[i]*m_fields[2]->GetPhys()[i]/ m_fields[0]->GetPhys()[i]));
	}
      
        // flux function for the \rho equation
        FX[0]  =  m_fields[1]->GetPhys();
        FY[0]  =  m_fields[2]->GetPhys();
      
        // Fill the flux vector using the old time step
        for (int i = 0; i < m_fields[0]->GetPointsTot(); ++i)
	{
            // flux function for the \rho u equation
            FX[1][i] = m_fields[1]->GetPhys()[i]*m_fields[1]->GetPhys()[i]/m_fields[0]->GetPhys()[i] + p[i];
            FY[1][i] = m_fields[1]->GetPhys()[i]*m_fields[2]->GetPhys()[i]/m_fields[0]->GetPhys()[i];
	  
            // flux function for the \rho v equation
            FX[2][i] = m_fields[1]->GetPhys()[i]*m_fields[2]->GetPhys()[i]/m_fields[0]->GetPhys()[i];
            FY[2][i] = m_fields[2]->GetPhys()[i]*m_fields[2]->GetPhys()[i]/m_fields[0]->GetPhys()[i] + p[i];
	  
            //flux function for the E equation
            FX[3][i] = (m_fields[1]->GetPhys()[i]/m_fields[0]->GetPhys()[i])*(m_fields[3]->GetPhys()[i] + p[i]);
            FY[3][i] = (m_fields[2]->GetPhys()[i]/m_fields[0]->GetPhys()[i])*(m_fields[3]->GetPhys()[i] + p[i]);
	}
    }


    void EulerEquations::NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &upX, 
				       Array<OneD, Array<OneD, NekDouble> > &upY)
    {
      
        // get temporary arrays
        Array<OneD, Array<OneD, NekDouble> > Fwd(GetNvariables());
        Array<OneD, Array<OneD, NekDouble> > Bwd(GetNvariables());

        for (int i = 0; i < GetNvariables(); ++i)
	{
            Fwd[i] = Array<OneD, NekDouble>(GetNpoints());
            Bwd[i] = Array<OneD, NekDouble>(GetNpoints());
	}
      
        // get the physical values at the trace
        for (int i = 0; i < GetNvariables(); ++i)
	{
            GetFwdBwdTracePhys(Fwd[i],Bwd[i],i);
	}

        // get the normal vector of the trace
        Array<OneD, Array<OneD, NekDouble> > normals(2);
        normals[0] = Array<OneD, NekDouble>(GetNpoints());
        normals[1] = Array<OneD, NekDouble>(GetNpoints());

        GetTraceNormals(normals);

        // rotate the values to the normal direction
        NekDouble tmpX, tmpY;
        for (int i = 0; i < GetNpoints(); ++i)
	{
            tmpX =  Fwd[1][i]*normals[0][i]+Fwd[2][i]*normals[1][i];
            tmpY = -Fwd[1][i]*normals[1][i]+Fwd[2][i]*normals[0][i];
            Fwd[1][i] = tmpX;
            Fwd[2][i] = tmpY;

            tmpX =  Bwd[1][i]*normals[0][i]+Bwd[2][i]*normals[1][i];
            tmpY = -Bwd[1][i]*normals[1][i]+Bwd[2][i]*normals[0][i];
            Bwd[1][i] = tmpX;
            Bwd[2][i] = tmpY;
	}

        // Solve the Riemann problem
        NekDouble rhoflux, rhouflux, rhovflux, Eflux;
      
        for (int i = 0; i < GetNpoints(); ++i)
	{
            RiemannSolver(Fwd[0][i],Fwd[1][i],Fwd[2][i],Fwd[3][i],
                          Bwd[0][i],Bwd[1][i],Bwd[2][i],Bwd[3][i],
                          rhoflux, rhouflux, rhovflux, Eflux );
	  
            // rotate back to Cartesian
            upX[0][i] =  rhoflux*normals[0][i];
            upY[0][i] =  rhoflux*normals[1][i];
            upX[1][i] = (rhouflux*normals[0][i] - rhovflux*normals[1][i]) * normals[0][i];
            upY[1][i] = (rhouflux*normals[0][i] - rhovflux*normals[1][i]) * normals[1][i];
            upX[2][i] = (rhouflux*normals[1][i] + rhovflux*normals[0][i]) * normals[0][i];
            upY[2][i] = (rhouflux*normals[1][i] + rhovflux*normals[0][i]) * normals[1][i];
            upX[3][i] =  Eflux*normals[0][i];
            upY[3][i] =  Eflux*normals[1][i];
	}
    }

    
    void EulerEquations::SetBoundaryConditions(void)
    {
        // loop over Boundary Regions
        for(int n = 0; n < m_fields[0]->GetBndCondExpansions().num_elements(); ++n)
	{	  
            // check if UserSpecified Boundary 
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() != "NoUserDefined")
	    {
                // Wall Boundary Condition
                if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() == "Wall")
		{
                    WallBoundary(n);
		}
                // Transparent Boundary Condition
                else if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() == "Transparent")
		{
                    // TODO:: initial values stored in Exp.
		}
	      
                // Timedependent Boundary Condition
                else if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() == "Timedependent")
		{
                    IsenTropicVortexBoundary(n);

                    //    m_fields[0]->EvaluateBoundaryConditions(time,m_fields[0]->GetBndCondExpansions()[n],
                    //				      m_fields[0]->GetBndConditions()[n]);
		  
                    // TODO: evaluate the bc equations
		}
	      
                // No matching User Defined Type
                else
		{
                    ASSERTL0(0, "No matching User Defined Boundary Condition.");
		}
	    }
	}
      
    }
    //----------------------------------------------------
    
    
    void EulerEquations::WallBoundary(int bcRegion)
    { 
        // get physical values of h, hu, hv for the forward trace
        Array<OneD, NekDouble> rho(GetNpoints());
        Array<OneD, NekDouble> rhou(GetNpoints());
        Array<OneD, NekDouble> rhov(GetNpoints());
        Array<OneD, NekDouble> E(GetNpoints());

        ExtractTracePhys(rho,0);
        ExtractTracePhys(rhou,1);
        ExtractTracePhys(rhov,2);
        ExtractTracePhys(E,3);
	 
        // get trace normals
        Array<OneD, Array<OneD, NekDouble> > normals(2);
        normals[0] = Array<OneD, NekDouble>(GetNpoints());
        normals[1] = Array<OneD, NekDouble>(GetNpoints());
        GetTraceNormals(normals);
               
        // Adjust the physical values of the trace to take 
        // user defined boundaries into account
        int e, id1, id2, npts, cnt = 0; 
      
        for(e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
	{
            npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetNumPoints(0);
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e) ;
            id2  = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));//GetBndExpToTraceExpMap(cnt+e));
	  
            Array<OneD, NekDouble> tmp_n(npts);
            Array<OneD, NekDouble> tmp_t(npts);
	  
            // rotate to compute the normal and tangential flux components
            Vmath::Vmul(npts,&rhou[id2],1,&normals[0][id2],1,&tmp_n[0],1);
            Vmath::Vvtvp(npts,&rhov[id2],1,&normals[1][id2],1,&tmp_n[0],1,&tmp_n[0],1);
	  
            Vmath::Vmul(npts,&rhou[id2],1,&normals[1][id2],1,&tmp_t[0],1);
            Vmath::Vvtvm(npts,&rhov[id2],1,&normals[0][id2],1,&tmp_t[0],1,&tmp_t[0],1);
	  
            // negate the normal flux
            Vmath::Neg(npts,tmp_n,1);		      
	  
            // rotate back to Cartesian
            Vmath::Vmul(npts,&tmp_t[0],1,&normals[1][id2],1,&rhou[id2],1);
            Vmath::Vvtvm(npts,&tmp_n[0],1,&normals[0][id2],1,&rhou[id2],1,&rhou[id2],1);
	  
            Vmath::Vmul(npts,&tmp_t[0],1,&normals[0][id2],1,&rhov[id2],1);
            Vmath::Vvtvp(npts,&tmp_n[0],1,&normals[1][id2],1,&rhov[id2],1,&rhov[id2],1);
	  
            // copy boundary adjusted values into the boundary expansion
            Vmath::Vcopy(npts,&rho[id2], 1,&(m_fields[0]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
            Vmath::Vcopy(npts,&rhou[id2],1,&(m_fields[1]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
            Vmath::Vcopy(npts,&rhov[id2],1,&(m_fields[2]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
            Vmath::Vcopy(npts,&E[id2],1,&(m_fields[3]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
	}
        cnt +=e;
    }

    void EulerEquations::IsenTropicVortexBoundary(int bcRegion)
    { 
        // get physical values of h, hu, hv for the forward trace
        Array<OneD, NekDouble> rho(GetNpoints());
        Array<OneD, NekDouble> rhou(GetNpoints());
        Array<OneD, NekDouble> rhov(GetNpoints());
        Array<OneD, NekDouble> E(GetNpoints());
        Array<OneD, NekDouble> x(GetNpoints());
        Array<OneD, NekDouble> y(GetNpoints());
        Array<OneD, NekDouble> z(GetNpoints());
  
        m_fields[0]->GetTrace()->GetCoords(x,y,z);
      
        for (int i = 0; i < GetNpoints(); ++i)
	{
            GetIsenTropicVortex(x[i], y[i], GetTime(), rho[i], rhou[i], rhov[i], E[i]);
	}
      
        // Adjust the physical values of the trace to take 
        // user defined boundaries into account
        int e, id1, id2, npts, cnt = 0; 
      
        for(e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
	{
            npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetNumPoints(0);
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e) ;
            id2  = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));//GetBndExpToTraceExpMap(cnt+e));
	  
            // copy the boundary  values into the boundary expansion
            Vmath::Vcopy(npts,&rho[id2], 1,&(m_fields[0]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
            Vmath::Vcopy(npts,&rhou[id2],1,&(m_fields[1]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
            Vmath::Vcopy(npts,&rhov[id2],1,&(m_fields[2]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
            Vmath::Vcopy(npts,&E[id2],1,&(m_fields[3]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
	}
        cnt +=e;
    }
    
    void EulerEquations::GetIsenTropicVortex(NekDouble x, NekDouble y, NekDouble time, NekDouble &rho, 
                                             NekDouble &rhou, NekDouble &rhov, NekDouble &E)
    {
      
        //---------------------------------
        // flow parameters
        NekDouble x0   = 5.0;
        NekDouble y0   = 0.0;
        NekDouble beta  = 5.0;
        NekDouble u0    = 1.0;
        NekDouble v0    = 0.0;
        NekDouble gamma = m_gamma;
        NekDouble r;
      
        r    = sqrt( pow(x-u0*time-x0, 2.0) + pow(y-v0*time-y0, 2.0));
        rho  = pow( (1.0-((gamma-1.0)/(16.0*gamma*M_PI*M_PI))*beta*beta*exp(2.0*(1.0-r*r))), (1.0/(gamma-1.0)) );
        rhou = rho * (1.0 - beta*exp(1.0-r*r)*((y-y0)/(2.0*M_PI)));
        rhov = rho * (beta*exp(1.0-r*r)*((x-x0)/(2.0*M_PI)));
        E    = (pow(rho,gamma)/(gamma-1.0)) + 0.5*rho*(pow(rhou/rho,2.0)+pow(rhov/rho,2.0));
    }
    


    void EulerEquations::RiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
				       NekDouble rhoR, NekDouble rhouR, NekDouble rhovR, NekDouble ER,
				       NekDouble &rhoflux, NekDouble &rhouflux, NekDouble &rhovflux, NekDouble &Eflux)
    {
  
        NekDouble gamma = GetGamma();
      
        NekDouble uL = rhouL/rhoL;
        NekDouble vL = rhovL/rhoL;
        NekDouble uR = rhouR/rhoR;
        NekDouble vR = rhovR/rhoR;
        NekDouble pL = (gamma - 1.0) * (EL - 0.5 * (rhouL*uL + rhovL*vL));
        NekDouble pR = (gamma - 1.0) * (ER - 0.5 * (rhouR*uR + rhovR*vR));
        NekDouble cL = sqrt(gamma * pL / rhoL);
        NekDouble cR = sqrt(gamma * pR / rhoR);
        NekDouble hL = (EL+pL)/rhoL;
        NekDouble hR = (ER+pR)/rhoR;
 
        // compute Roe averages
        NekDouble rhoRoe = sqrt(rhoL) * sqrt(rhoR);
        NekDouble uRoe   = (sqrt(rhoL)*uL + sqrt(rhoR)*uR)/(sqrt(rhoL)+sqrt(rhoR));
        NekDouble vRoe   = (sqrt(rhoL)*vL + sqrt(rhoR)*vR)/(sqrt(rhoL)+sqrt(rhoR));
        NekDouble hRoe   = (sqrt(rhoL)*hL + sqrt(rhoR)*hR)/(sqrt(rhoL)+sqrt(rhoR));
        NekDouble cRoe   = sqrt( (gamma - 1.0)*(hRoe - 0.5*(uRoe*uRoe + vRoe*vRoe)) );

        // compute the wave speeds
        NekDouble SL = min(uL-cL, uRoe-cRoe);
        NekDouble SR = max(uR+uR, uRoe+cRoe);

        // compute the HLL flux
        if (SL >= 0)
	{
            rhoflux  = rhoL*uL;
            rhouflux = rhoL*uL*uL + pL;
            rhovflux = rhoL*uL*vL;
            Eflux    = uL*(EL + pL);
	}
        else if (SR <= 0)
	{
            rhoflux  = rhoR*uR;
            rhouflux = rhoR*uR*uR + pR;
            rhovflux = rhoR*uR*vR;
            Eflux    = uR*(ER + pR);
	}
        else
	{
            rhoflux  = (( SR*(rhoL*uL) - SL*(rhoR*uR)+SR*SL*(rhoR-rhoL))/(SR-SL) );
            rhouflux = (( SR*(rhoL*uL*uL + pL) - SL*(rhoR*uR*uR + pR) +
                          SR*SL*(rhouR-rhouL)) / (SR-SL) );
            rhovflux = (( SR*(rhoL*uL*vL) - SL*(rhoR*uR*vR) +
                          SR*SL*(rhovR-rhovL)) / (SR-SL) );
            Eflux    = (( SR*(uL*EL+uL*pL) - SL*(uR*ER+uR*pR) +
                          SR*SL*(ER-EL)) / (SR-SL) );
	}

    }

} //end of namespace

/**
* $Log: EulerEquations.cpp,v $
* Revision 1.1  2008/08/22 09:48:23  pvos
* Added Claes' AdvectionDiffusionReaction, ShallowWater and Euler solver
*
**/
