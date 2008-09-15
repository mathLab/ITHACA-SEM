///////////////////////////////////////////////////////////////////////////////
//
// File ShallowWaterEquations.cpp
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
// Description: Shallow water equations routines
//
///////////////////////////////////////////////////////////////////////////////

#include <../solvers/Auxiliary/ShallowWaterEquations.h>

namespace Nektar
{
    ShallowWaterEquations::ShallowWaterEquations(void):
        AdvectionDiffusionReaction()
    {
    }
     
 //  ShallowWaterEquations::ShallowWaterEquations(ShallowWaterEquations &In):
//         AdvectionDiffusionReaction(In)
//     {
//     }
  
    ShallowWaterEquations::ShallowWaterEquations(SpatialDomains::MeshGraph2D &graph2D,
                                                 SpatialDomains::BoundaryConditions &bcs,
						 int &nVariables):
        AdvectionDiffusionReaction(graph2D,bcs,nVariables)
    {
        // here we must check that all read userdefined boundaries
        // are INDEED implemented...
      
    }


  void ShallowWaterEquations::ConservativeToPrimitive(void)
  {
    
    // here we should then add a check that the 
    // water depth is positive
    for (int i = 0; i < m_fields[0]->GetPointsTot(); ++i)
      {
	
	// Hu->u
	(m_fields[1]->UpdatePhys())[i] = ((m_fields[1]->GetPhys())[i]/
					  (m_fields[0]->GetPhys())[i]);
	
	// Hv->v
	(m_fields[2]->UpdatePhys())[i] = ((m_fields[2]->GetPhys())[i]/
					  (m_fields[0]->GetPhys())[i]);
	
	// H -> eta
	(m_fields[0]->UpdatePhys())[i] = ((m_fields[0]->GetPhys())[i] - 
					  (m_fields[4]->GetPhys())[i]);

      }

  }

  void ShallowWaterEquations::PrimitiveToConservative(void)
  {   
    // here we should then add a check that the 
    // water depth is positive
    for (int i = 0; i < m_fields[0]->GetPointsTot(); ++i)
      {
	// eta->H
	(m_fields[0]->UpdatePhys())[i] = ((m_fields[0]->GetPhys())[i] + 
					  (m_fields[4]->GetPhys())[i]);
	
	
	// u->Hu
	(m_fields[1]->UpdatePhys())[i] = ((m_fields[1]->GetPhys())[i] *
					  (m_fields[0]->GetPhys())[i]);
	
	// v->Hv
	(m_fields[2]->UpdatePhys())[i] = ((m_fields[2]->GetPhys())[i] *
					  (m_fields[0]->GetPhys())[i]);
      }
  }

  void ShallowWaterEquations::SetCoriolis(SpatialDomains::BoundaryConditions &bcs)
  {
    
    // Hardcoded and horrible...
    // Should be changed to be read as a function...
    
        if (bcs.GetParameter("Coriolis"))
	{
            int nTotQuadPoints  = GetPointsTot();
	  
            Array<OneD,NekDouble> x1(nTotQuadPoints,0.0);
            Array<OneD,NekDouble> x2(nTotQuadPoints,0.0);
	  
            GetCoords(x1,x2);
	  
            m_coriolis = Array<OneD, NekDouble>(nTotQuadPoints);
	  
            for(int i = 0; i < nTotQuadPoints; ++i)
	    {
                m_coriolis[i]   = 0.0+1.0*x2[i];
	    }      
	}
    }


  void ShallowWaterEquations::GetFluxVector(Array<OneD, Array<OneD, NekDouble> >&FX,
					    Array<OneD, Array<OneD, NekDouble> >&FY)
  { 
    
    NekDouble g = m_g;
        
    // Fill the flux vector using the old time step
    for (int i = 0; i < m_fields[0]->GetPointsTot(); ++i)
      {
	// flux function for the h equation
	FX[0][i]  =  m_fields[1]->GetPhys()[i];
	FY[0][i]  =  m_fields[2]->GetPhys()[i];
	
	// flux function for the hu equation
	FX[1][i] = m_fields[1]->GetPhys()[i]*m_fields[1]->GetPhys()[i]/m_fields[0]->GetPhys()[i] +
	  0.5*g*m_fields[0]->GetPhys()[i]*m_fields[0]->GetPhys()[i];
	FY[1][i] = m_fields[1]->GetPhys()[i]*m_fields[2]->GetPhys()[i]/m_fields[0]->GetPhys()[i];
	
	// flux function for the hv equation
	FX[2][i] = m_fields[1]->GetPhys()[i]*m_fields[2]->GetPhys()[i]/m_fields[0]->GetPhys()[i];
	FY[2][i] = m_fields[2]->GetPhys()[i]*m_fields[2]->GetPhys()[i]/m_fields[0]->GetPhys()[i]+
	  0.5*g*m_fields[0]->GetPhys()[i]*m_fields[0]->GetPhys()[i];
      }
  }
  
  void ShallowWaterEquations::GetFluxVectorPrimitive(Array<OneD, Array<OneD, NekDouble> >&FX,
						     Array<OneD, Array<OneD, NekDouble> >&FY)
  { 
    
    NekDouble g = m_g;
    
    NekDouble d = m_d;
    
    // Fill the flux vector using the old time step
    for (int i = 0; i < m_fields[0]->GetPointsTot(); ++i)
      {
	// flux function for the h equation
	FX[0][i]  =  (d + m_fields[0]->GetPhys()[i]) * m_fields[1]->GetPhys()[i];
	FY[0][i]  =  (d + m_fields[0]->GetPhys()[i]) * m_fields[2]->GetPhys()[i];
	
	// flux function for the hu equation
	FX[1][i] = g * m_fields[0]->GetPhys()[i] + m_fields[1]->GetPhys()[i]*m_fields[1]->GetPhys()[i]
	  + m_fields[2]->GetPhys()[i]*m_fields[2]->GetPhys()[i];
	FY[1][i] = 0.0;
	
	// flux function for the hv equation
	FX[2][i] = 0.0;
	FY[2][i] = g * m_fields[0]->GetPhys()[i] + m_fields[1]->GetPhys()[i]*m_fields[1]->GetPhys()[i]
	  + m_fields[2]->GetPhys()[i]*m_fields[2]->GetPhys()[i];
      }
  }

  void ShallowWaterEquations::GetFluxVectorPrimitiveLinear(Array<OneD, Array<OneD, NekDouble> >&FX,
							   Array<OneD, Array<OneD, NekDouble> >&FY)
  { 
    
    NekDouble g = m_g;
    
    NekDouble d = m_d;

    // Fill the flux vector using the old time step
    for (int i = 0; i < m_fields[0]->GetPointsTot(); ++i)
      {
	// flux function for the h equation
	FX[0][i]  =  (d) * m_fields[1]->GetPhys()[i];
	FY[0][i]  =  (d) * m_fields[2]->GetPhys()[i];

	// flux function for the hu equation
	FX[1][i] = g * m_fields[0]->GetPhys()[i];
	FY[1][i] = 0.0;
	
	// flux function for the hv equation
	FX[2][i] = 0.0;
	FY[2][i] = g * m_fields[0]->GetPhys()[i];
      }
  }

  

    void ShallowWaterEquations::NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &upX, 
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
        NekDouble hflux, huflux, hvflux;
      
        for (int i = 0; i < GetNpoints(); ++i)
	{
            RiemannSolver(Fwd[0][i],Fwd[1][i],Fwd[2][i],
                          Bwd[0][i],Bwd[1][i],Bwd[2][i],
                          hflux, huflux, hvflux );
	  
            // rotate back to Cartesian
            upX[0][i]  = hflux*normals[0][i];
            upY[0][i]  = hflux*normals[1][i];
            upX[1][i] = (huflux*normals[0][i] - hvflux*normals[1][i]) * normals[0][i];
            upY[1][i] = (huflux*normals[0][i] - hvflux*normals[1][i]) * normals[1][i];
            upX[2][i] = (huflux*normals[1][i] + hvflux*normals[0][i]) * normals[0][i];
            upY[2][i] = (huflux*normals[1][i] + hvflux*normals[0][i]) * normals[1][i];
	}
      
      
    }
  
  void ShallowWaterEquations::NumericalFluxPrimitiveLinear(Array<OneD, Array<OneD, NekDouble> > &upX, 
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
            
      // Solve the Riemann problem
      NekDouble eta, u, v;
      NekDouble g = m_g;
      NekDouble d = m_d;
      
      // averaging
      for (int i = 0; i < GetNpoints(); ++i)
	{
	  eta = 0.5*(Fwd[0][i] + Bwd[0][i]);
	  u   = 0.5*(Fwd[1][i] + Bwd[1][i]);
	  v   = 0.5*(Fwd[2][i] + Bwd[2][i]);
	  
	  upX[0][i]  = d * u;
	  upY[0][i]  = d * v;
	  upX[1][i]  = g * eta;
	  upY[1][i]  = 0.0;
	  upX[2][i]  = 0.0;
	  upY[2][i]  = g * eta;
	}
    }
  
  
  //----------------------------------------------------
    void ShallowWaterEquations::SetBoundaryConditions(void)
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
    
    
    void ShallowWaterEquations::WallBoundary(int bcRegion)
    { 
        // get physical values of h, hu, hv for the forward trace
        Array<OneD, NekDouble> h(GetNpoints());
        Array<OneD, NekDouble> hu(GetNpoints());
        Array<OneD, NekDouble> hv(GetNpoints());
        ExtractTracePhys(h,0);
        ExtractTracePhys(hu,1);
        ExtractTracePhys(hv,2);

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
            id2  = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndExpToTraceExpMap(cnt+e));
	  
            Array<OneD, NekDouble> tmp_n(npts);
            Array<OneD, NekDouble> tmp_t(npts);
	  
            // rotate to compute the normal and tangential flux components
            Vmath::Vmul(npts,&hu[id2],1,&normals[0][id2],1,&tmp_n[0],1);
            Vmath::Vvtvp(npts,&hv[id2],1,&normals[1][id2],1,&tmp_n[0],1,&tmp_n[0],1);
	  
            Vmath::Vmul(npts,&hu[id2],1,&normals[1][id2],1,&tmp_t[0],1);
            Vmath::Vvtvm(npts,&hv[id2],1,&normals[0][id2],1,&tmp_t[0],1,&tmp_t[0],1);
	  
            // negate the normal flux
            Vmath::Neg(npts,tmp_n,1);		      
	  
            // rotate back to Cartesian
            Vmath::Vmul(npts,&tmp_t[0],1,&normals[1][id2],1,&hu[id2],1);
            Vmath::Vvtvm(npts,&tmp_n[0],1,&normals[0][id2],1,&hu[id2],1,&hu[id2],1);
	  
            Vmath::Vmul(npts,&tmp_t[0],1,&normals[0][id2],1,&hv[id2],1);
            Vmath::Vvtvp(npts,&tmp_n[0],1,&normals[1][id2],1,&hv[id2],1,&hv[id2],1);
	  
            // copy boundary adjusted values into the boundary expansion
            Vmath::Vcopy(npts,&h[id2], 1,&(m_fields[0]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
            Vmath::Vcopy(npts,&hu[id2],1,&(m_fields[1]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
            Vmath::Vcopy(npts,&hv[id2],1,&(m_fields[2]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
	}
        cnt +=e;
    }
    
    
    void ShallowWaterEquations::RiemannSolver(NekDouble hL,NekDouble huL,NekDouble hvL,NekDouble hR,NekDouble huR, 
                                              NekDouble hvR, NekDouble &hflux, NekDouble &huflux,NekDouble &hvflux )
    {
  
        NekDouble hC,huC,hvC,SL,SR,hstar,Sstar;
      
        NekDouble g = GetGravity();
      
        NekDouble uL = huL/hL;
        NekDouble vL = hvL/hL;
        NekDouble uR = huR/hR;
        NekDouble vR = hvR/hR;
        NekDouble cL = sqrt(g * hL);
        NekDouble cR = sqrt(g * hR);
      
        // the two-rarefaction wave assumption
        hstar = 0.5*(cL + cR) + 0.25*(uL - uR);
        hstar *= hstar;
        hstar *= (1.0/g);
      
        // Compute SL
        if (hstar > hL)
            SL = uL - cL * sqrt(0.5*((hstar*hstar + hstar*hL)/(hL*hL)));
        else
            SL = uL - cL;
      
        // Compute SR
        if (hstar > hR)
            SR = uR + cR * sqrt(0.5*((hstar*hstar + hstar*hR)/(hR*hR)));
        else
            SR = uR + cR;
      
      
        if (fabs(hR*(uR-SR)-hL*(uL-SL)) <= 1.0e-15)
            Sstar = 0.0;
        else
            Sstar = (SL*hR*(uR-SR)-SR*hL*(uL-SL))/(hR*(uR-SR)-hL*(uL-SL));
      
        if (SL >= 0)
	{
            hflux  = hL * uL;
            huflux  = uL * uL * hL + 0.5 * g * hL * hL;
            hvflux  = hL * uL * vL;
	}
        else if (SR <= 0)
	{
            hflux  = hR * uR;
            huflux  = uR * uR * hR + 0.5 * g * hR * hR;
            hvflux  = hR * uR *vR;
	}
        else
	{
            if ((SL < 0) && (Sstar >= 0))
	    {
                hC  = hL * ((SL - uL) / (SL - Sstar));
                huC = hC * Sstar;
                hvC = hC * vL;
	      
                hflux = hL*uL + SL * (hC - hL);
                huflux = (uL*uL*hL+0.5*g*hL*hL) + SL * (huC - hL*uL);
                hvflux = (uL*vL*hL) + SL * (hvC - hL*vL);
	    }
            else
	    {
                hC  = hR * ((SR - uR) / (SR - Sstar));
                huC = hC * Sstar;
                hvC = hC * vR;
	      
                hflux = hR*uR + SR * (hC - hR);
                huflux = (uR*uR*hR+0.5*g*hR*hR) + SR * (huC - hR*uR);
                hvflux = (uR*vR*hR) + SR * (hvC - hR*vR);
	    }
	}
      
    }
} //end of namespace

/**
* $Log: ShallowWaterEquations.cpp,v $
* Revision 1.1  2008/08/22 09:48:23  pvos
* Added Claes' AdvectionDiffusionReaction, ShallowWater and Euler solver
*
**/
