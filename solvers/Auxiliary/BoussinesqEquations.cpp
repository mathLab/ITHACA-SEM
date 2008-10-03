///////////////////////////////////////////////////////////////////////////////
//
// File BoussinesqEquations.cpp
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
// Description: Boussinesq equations routines
//
///////////////////////////////////////////////////////////////////////////////

#include <../solvers/Auxiliary/BoussinesqEquations.h>

namespace Nektar
{
  BoussinesqEquations::BoussinesqEquations(void):
        ShallowWaterEquations()
  {
  }
  
//   BoussinesqEquations::BoussinesqEquations(const BoussinesqEquations &In):
//     ShallowWaterEquations(In)
//   {
//   }
  
  BoussinesqEquations::BoussinesqEquations(SpatialDomains::MeshGraph2D &graph2D,
					   SpatialDomains::BoundaryConditions &bcs,
					   int &nVariables):
    ShallowWaterEquations(graph2D,bcs,nVariables)
  {
    // here we must check that all read userdefined boundaries
    // are INDEED implemented...
    
  }
  
  void BoussinesqEquations::SetBoundaryConditionsWaveCont(void)
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
		  WallBoundaryWaveCont(n);
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
  
  
  void BoussinesqEquations::WallBoundaryWaveCont(int bcRegion)
  { 
    
    //std::cout << " WallBoundaryWaveCont" << std::endl;
    
    // get physical values of h, hu, hv for the forward trace
    Array<OneD, NekDouble> f1(GetNpoints());
    Array<OneD, NekDouble> f2(GetNpoints());
    ExtractTracePhys(f1,1);
    ExtractTracePhys(f2,2);
    
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
	Vmath::Vmul(npts,&f1[id2],1,&normals[0][id2],1,&tmp_n[0],1);
	Vmath::Vvtvp(npts,&f2[id2],1,&normals[1][id2],1,&tmp_n[0],1,&tmp_n[0],1);
	
	Vmath::Vmul(npts,&f1[id2],1,&normals[1][id2],1,&tmp_t[0],1);
	Vmath::Vvtvm(npts,&f2[id2],1,&normals[0][id2],1,&tmp_t[0],1,&tmp_t[0],1);
	
	// negate the normal flux
	Vmath::Neg(npts,tmp_n,1);		      
	
	// rotate back to Cartesian
	Vmath::Vmul(npts,&tmp_t[0],1,&normals[1][id2],1,&f1[id2],1);
	Vmath::Vvtvm(npts,&tmp_n[0],1,&normals[0][id2],1,&f1[id2],1,&f1[id2],1);
	
	Vmath::Vmul(npts,&tmp_t[0],1,&normals[0][id2],1,&f2[id2],1);
	Vmath::Vvtvp(npts,&tmp_n[0],1,&normals[1][id2],1,&f2[id2],1,&f2[id2],1);
	
	// copy boundary adjusted values into the boundary expansion
	Vmath::Vcopy(npts,&f1[id2],1,&(m_fields[1]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
	Vmath::Vcopy(npts,&f2[id2],1,&(m_fields[2]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
      }
    cnt +=e;
  }
  

  void BoussinesqEquations::NumericalFluxWaveCont(Array<OneD, NekDouble> &outX, Array<OneD, NekDouble> &outY)
  {
    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble> > Fwd(2);
    Array<OneD, Array<OneD, NekDouble> > Bwd(2);
    
    for (int i = 0; i < 2; ++i)
      {
	Fwd[i] = Array<OneD, NekDouble>(GetNpoints());
	Bwd[i] = Array<OneD, NekDouble>(GetNpoints());
      }
    
    // get the physical values at the trace
    GetFwdBwdTracePhys(Fwd[0],Bwd[0],1);
    GetFwdBwdTracePhys(Fwd[1],Bwd[1],2);
    
    for (int i = 0; i < GetNpoints(); ++i)
      {
	outX[i]  = 0.5*(Fwd[0][i] + Bwd[0][i]);
	outY[i]  = 0.5*(Fwd[1][i] + Bwd[1][i]);
      }
    
  }
  
  void BoussinesqEquations::SetBoundaryConditionsSolve(void)
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
		WallBoundarySolve(n);
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
  
  void BoussinesqEquations::WallBoundarySolve(int bcRegion)
  { 
    
    //cout << " WallBoundarySolve" << endl;
    
    // get physical values of h, hu, hv for the forward trace
    Array<OneD, NekDouble> z(GetNpoints(),0.0);
    ExtractTracePhys(z,3);
    
    // Adjust the physical values of the trace to take 
    // user defined boundaries into account
    int e, id1, id2, npts, cnt = 0; 
    
    for(e = 0; e < m_fields[3]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
      {
	npts = m_fields[3]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetNumPoints(0);
	id1  = m_fields[3]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e) ;
	id2  = m_fields[3]->GetTrace()->GetPhys_Offset(m_fields[3]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));//GetBndExpToTraceExpMap(cnt+e));
	
	// copy boundary adjusted values into the boundary expansion
	Vmath::Vcopy(npts,&z[id2],1,&(m_fields[3]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
      }
    cnt +=e;
  }
  
  void BoussinesqEquations::WaveContSolve(Array<OneD, NekDouble> &fce, NekDouble lambda)
  {
    
    MultiRegions::DisContField2DSharedPtr Fce;
    Fce = MemoryManager<MultiRegions::DisContField2D>::AllocateSharedPtr(*m_fields[3]);
    
    Fce->SetPhys(fce);
    
    m_fields[3]->HelmSolve(*Fce, lambda);
    
    m_fields[3]->BwdTrans(*m_fields[3]);
    
  }


  void BoussinesqEquations::SetBoundaryConditionsContVariables(void)
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
		WallBoundaryContVariables(n);
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
  
  void BoussinesqEquations::WallBoundaryContVariables(int bcRegion)
  { 
    
    //cout << " WallBoundaryContVariables" << endl;
    
    // get physical values of h, hu, hv for the forward trace
    Array<OneD, NekDouble> z(GetNpoints(),0.0);
    ExtractTracePhys(z,3);
    // Vmath::Neg(GetNpoints(),z,1);
    
    // Adjust the physical values of the trace to take 
    // user defined boundaries into account
    int e, id1, id2, npts, cnt = 0; 
    
    for(e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
      {
	npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetNumPoints(0);
	id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e) ;
	id2  = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));//GetBndExpToTraceExpMap(cnt+e));
	
	// copy boundary adjusted values into the boundary expansion
	Vmath::Vcopy(npts,&z[id2],1,&(m_fields[3]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
	
      }
    cnt +=e;
  }
  
  void BoussinesqEquations::NumericalFluxConsVariables(Array<OneD, NekDouble> &outX, Array<OneD, NekDouble> &outY)
  {
    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble> > Fwd(1);
    Array<OneD, Array<OneD, NekDouble> > Bwd(1);
    
    Fwd[0] = Array<OneD, NekDouble> (GetNpoints(),0.0); 
    Bwd[0] = Array<OneD, NekDouble> (GetNpoints(),0.0);
    
    // get the physical values at the trace
    GetFwdBwdTracePhys(Fwd[0],Bwd[0],3);
    
    for (int i = 0; i < GetNpoints(); ++i)
      {
 	outX[i]  = 0.5*(Fwd[0][i] + Bwd[0][i]);
	outY[i]  = 0.5*(Fwd[0][i] + Bwd[0][i]);
      }
    
  }




  void BoussinesqEquations::Madsen92SpatialTerms(Array<OneD, NekDouble> &outX, Array<OneD, NekDouble> &outY)
  {

    //--------------------------------------
    // local parameters
    int nTotQuadPoints = GetPointsTot();
    int nTotCoeffs = GetNcoeffs();
    
    NekDouble d = GetConstantDepth();
    NekDouble g = GetGravity();
    NekDouble B = GetAlpha1();
    //--------------------------------------

    
    //--------------------------------------
    // local arrays
    Array<OneD, Array<OneD, NekDouble> > b(2);
    b[0] = Array<OneD, NekDouble>(nTotQuadPoints);
    b[1] = Array<OneD, NekDouble>(nTotQuadPoints);
    Array<OneD, NekDouble > a(nTotQuadPoints);

    Array<OneD, NekDouble > store_u(nTotQuadPoints);
    Array<OneD, NekDouble > store_v(nTotQuadPoints);

    Array<OneD, NekDouble> tmpX(GetNcoeffs());
    Array<OneD, NekDouble> tmpY(GetNcoeffs());
    
    Array<OneD, NekDouble> upwindX(GetNpoints());
    Array<OneD, NekDouble> upwindY(GetNpoints());
    //--------------------------------------

    
    //--------------------------------------
    // Store the physical values of u and v
    
    Vmath::Vcopy(nTotQuadPoints,GetPhys(1),1,store_u,1);
    Vmath::Vcopy(nTotQuadPoints,GetPhys(2),1,store_v,1);
    //--------------------------------------

    //--------------------------------------
    // solve {\bf b} = d \nabla \eta
    // (we store b in u and v)
    
    IProductWRTDerivBase(0,GetPhys(0),tmpX,0);
    IProductWRTDerivBase(1,GetPhys(0),tmpY,0);
    
    Vmath::Neg(nTotCoeffs,tmpX,1);
    Vmath::Neg(nTotCoeffs,tmpY,1);
    
//     // Evaluate  upwind numerical flux (physical space)
//     NumericalFluxGradient(upwindX,upwindY,0);
    
//     Array<OneD, NekDouble> uptemp(GetNpoints(),0.0);
    
//     AddTraceIntegral(upwindX,uptemp,tmpX,0);
//     AddTraceIntegral(uptemp,upwindY,tmpY,0);

    Vmath::Smul(nTotCoeffs,d,tmpX,1,tmpX,1);
    Vmath::Smul(nTotCoeffs,d,tmpY,1,tmpY,1);

    MultiplyByElmtInvMass(tmpX,UpdateCoeffs(1),0);
    MultiplyByElmtInvMass(tmpY,UpdateCoeffs(2),0);
    
    BwdTrans(1);
    BwdTrans(2);
    //--------------------------------------
    

    //--------------------------------------
    // solve a = \nabla \cdot {\bf b}
    // (we store a in u)
    
    IProductWRTDerivBase(0,GetPhys(1),tmpX,0);
    IProductWRTDerivBase(1,GetPhys(2),tmpY,0);
    
    Vmath::Vadd(nTotCoeffs,tmpX,1,tmpY,1,tmpX,1);

    Vmath::Neg(nTotCoeffs,tmpX,1);
    
   //  // Evaluate  upwind numerical flux (physical space)
//     NumericalFluxDivergence(upwindX,upwindY,1,2);
    
//     AddTraceIntegral(upwindX,upwindY,tmpX,0);
    
    MultiplyByElmtInvMass(tmpX,UpdateCoeffs(1),0);
    
    BwdTrans(1);
    //--------------------------------------
    

    //--------------------------------------
    // solve {\bf D^s} = - g B d \nabla a
    // (we store D^s in u and v)
    
    IProductWRTDerivBase(0,GetPhys(1),tmpX,0);
    IProductWRTDerivBase(1,GetPhys(1),tmpY,0);
    
    Vmath::Neg(nTotCoeffs,tmpX,1);
    Vmath::Neg(nTotCoeffs,tmpY,1);
    
    // Evaluate  upwind numerical flux (physical space)
  //   NumericalFluxGradient(upwindX,upwindY,1);
    
//     AddTraceIntegral(upwindX,uptemp,tmpX,0);
//     AddTraceIntegral(uptemp,upwindY,tmpY,0);

    Vmath::Smul(nTotCoeffs,g*B*d,tmpX,1,tmpX,1);
    Vmath::Smul(nTotCoeffs,g*B*d,tmpY,1,tmpY,1);

    //MultiplyByElmtInvMass(tmpX,UpdateCoeffs(1),0);
    //MultiplyByElmtInvMass(tmpY,UpdateCoeffs(2),0);
    
    //BwdTrans(1);
    //BwdTrans(2);
    //--------------------------------------
    
    
    //--------------------------------------
    // Add D^s to the RHS terms 
    
    Vmath::Vadd(nTotCoeffs,outX,1,tmpX,1,outX,1);
    Vmath::Vadd(nTotCoeffs,outY,1,tmpY,1,outY,1);
    //--------------------------------------

    
    //--------------------------------------
    // fill u and v with the ingoing values
    
    Vmath::Vcopy(nTotQuadPoints,store_u,1,UpdatePhys(1),1);
    Vmath::Vcopy(nTotQuadPoints,store_v,1,UpdatePhys(2),1);
    //--------------------------------------
    

  }

  /**
   * \brief This function evaluates the spatial linear dispersive terms 
   *
   * \f[ \Lambda_{20} = -\nabla \left( \nabla \cdot \left(H_t {\bf u}\right) \right) \f]
   *
   **/ 
  
  void BoussinesqEquations::Lambda20Spatial(BoussinesqEquations &Uf, Array<OneD, NekDouble> &outX, 
					    Array<OneD, NekDouble> &outY, Array<OneD, NekDouble> &Ht)
  {
    
    
    //--------------------------------------
    // local parameters
    int nTotQuadPoints = GetPointsTot();
    int nTotCoeffs = GetNcoeffs();
    
    NekDouble d = GetConstantDepth();
    NekDouble g = GetGravity();
    NekDouble alpha1 = GetAlpha1();
    NekDouble alpha2 = GetAlpha2();
    //--------------------------------------
    
    
    //--------------------------------------
    // local arrays
    Array<OneD, NekDouble> tmpX(GetNcoeffs());
    Array<OneD, NekDouble> tmpY(GetNcoeffs());
    
    Array<OneD, NekDouble> upwindX(GetNpoints());
    Array<OneD, NekDouble> upwindY(GetNpoints());
    Array<OneD, NekDouble> upwindZero(GetNpoints(),0.0);

    //--------------------------------------------------
    
    
    //---------------------------------------------------
    // Add the 
    // g H d^2 (\alpha_2 - \alpha_1)\nabla(\nabla^2 \eta)
    // term to the rhs 
    

    // solve {\bf b} =  \nabla \eta
    // [we store b in Uf(1) and Uf(2)]
    
    // eta in Uf(0)
    Vmath::Vsub(nTotQuadPoints,GetPhys(0),1,GetPhys(4),1,Uf.UpdatePhys(0),1);
    
    IProductWRTDerivBase(0,Uf.GetPhys(0),tmpX,0);
    IProductWRTDerivBase(1,Uf.GetPhys(0),tmpY,0);
    
    Vmath::Neg(nTotCoeffs,tmpX,1);
    Vmath::Neg(nTotCoeffs,tmpY,1);
    
    Uf.NumericalFluxGradient(upwindX,upwindY,0);
    
    AddTraceIntegral(upwindX,upwindZero,tmpX,0);
    AddTraceIntegral(upwindZero,upwindY,tmpY,0);
    
    MultiplyByElmtInvMass(tmpX,Uf.UpdateCoeffs(1),0);
    MultiplyByElmtInvMass(tmpY,Uf.UpdateCoeffs(2),0);
    
    Uf.BwdTrans(1);
    Uf.BwdTrans(2);
    

    // solve a = \nabla \cdot {\bf b}
    // [we store a in Uf(1)]
     
     IProductWRTDerivBase(0,Uf.GetPhys(1),tmpX,0);
     IProductWRTDerivBase(1,Uf.GetPhys(2),tmpY,0);
     
     Vmath::Vadd(nTotCoeffs,tmpX,1,tmpY,1,tmpX,1); 
     Vmath::Neg(nTotCoeffs,tmpX,1);
     
     Uf.NumericalFluxDivergence(upwindX,upwindY,1,2);
     
     AddTraceIntegral(upwindX,upwindY,tmpX,0);
     
     MultiplyByElmtInvMass(tmpX,Uf.UpdateCoeffs(1),0);
     
     Uf.BwdTrans(1);
     
     
     // solve {\bf c} =  \nabla a
     // [we store c in Uf(1) and Uf(2)]
     
     IProductWRTDerivBase(0,Uf.GetPhys(1),tmpX,0);
     IProductWRTDerivBase(1,Uf.GetPhys(1),tmpY,0);
     
     Vmath::Neg(nTotCoeffs,tmpX,1);
     Vmath::Neg(nTotCoeffs,tmpY,1);
    
     Uf.NumericalFluxGradient(upwindX,upwindY,1);
     
     AddTraceIntegral(upwindX,upwindZero,tmpX,0);
     AddTraceIntegral(upwindZero,upwindY,tmpY,0);
     
     MultiplyByElmtInvMass(tmpX,Uf.UpdateCoeffs(1),0);
     MultiplyByElmtInvMass(tmpY,Uf.UpdateCoeffs(2),0);
     
     Uf.BwdTrans(1);
     Uf.BwdTrans(2);
     
    
     // compute {\bf Lambda_20^1  = g d^2 H (alpha2-alpha1){\bf c}
    
     Vmath::Vmul(nTotQuadPoints,GetPhys(0),1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,GetPhys(0),1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);
          
     Vmath::Smul(nTotQuadPoints,g*(alpha2-alpha1),Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Smul(nTotQuadPoints,g*(alpha2-alpha1),Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     // due to move to rhs
     Vmath::Smul(nTotQuadPoints,-1.0,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Smul(nTotQuadPoints,-1.0,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     IProductWRTBase(Uf.GetPhys(1),tmpX,0);
     IProductWRTBase(Uf.GetPhys(2),tmpY,0);
     
     // Add to the RHS terms 
     
     Vmath::Vadd(nTotCoeffs,outX,1,tmpX,1,outX,1);
     Vmath::Vadd(nTotCoeffs,outY,1,tmpY,1,outY,1);
     //--------------------------------------
 
     
     //---------------------------------------------------
     // Add the 
     // -g H d \alpha_2 \nabla(\nabla \cdot (d \nabla \eta))
     // term to the rhs 
    

    // solve {\bf b} =  \nabla \eta
    // [we store b in Uf(1) and Uf(2)]
    
    // eta in Uf(0)
     Vmath::Vsub(nTotQuadPoints,GetPhys(0),1,GetPhys(4),1,Uf.UpdatePhys(0),1);

    IProductWRTDerivBase(0,Uf.GetPhys(0),tmpX,0);
    IProductWRTDerivBase(1,Uf.GetPhys(0),tmpY,0);
    
    Vmath::Neg(nTotCoeffs,tmpX,1);
    Vmath::Neg(nTotCoeffs,tmpY,1);
    
    Uf.NumericalFluxGradient(upwindX,upwindY,0);
    
    AddTraceIntegral(upwindX,upwindZero,tmpX,0);
    AddTraceIntegral(upwindZero,upwindY,tmpY,0);
    
    MultiplyByElmtInvMass(tmpX,Uf.UpdateCoeffs(1),0);
    MultiplyByElmtInvMass(tmpY,Uf.UpdateCoeffs(2),0);
    
    Uf.BwdTrans(1);
    Uf.BwdTrans(2);
    

    // solve a = \nabla \cdot d{\bf b}
    // [we store a in Uf(1)]

    Vmath::Vmul(nTotQuadPoints,Uf.GetPhys(1),1,GetPhys(4),1,Uf.UpdatePhys(1),1);
    Vmath::Vmul(nTotQuadPoints,Uf.GetPhys(2),1,GetPhys(4),1,Uf.UpdatePhys(2),1);
    
    IProductWRTDerivBase(0,Uf.GetPhys(1),tmpX,0);
    IProductWRTDerivBase(1,Uf.GetPhys(2),tmpY,0);
     
    Vmath::Vadd(nTotCoeffs,tmpX,1,tmpY,1,tmpX,1); 
    Vmath::Neg(nTotCoeffs,tmpX,1);
    
    Uf.NumericalFluxDivergence(upwindX,upwindY,1,2);
    
    AddTraceIntegral(upwindX,upwindY,tmpX,0);
    
    MultiplyByElmtInvMass(tmpX,Uf.UpdateCoeffs(1),0);
    
    Uf.BwdTrans(1);
    
     
    // solve {\bf c} =  \nabla a
    // [we store c in Uf(1) and Uf(2)]
     
    IProductWRTDerivBase(0,Uf.GetPhys(1),tmpX,0);
    IProductWRTDerivBase(1,Uf.GetPhys(1),tmpY,0);
    
    Vmath::Neg(nTotCoeffs,tmpX,1);
    Vmath::Neg(nTotCoeffs,tmpY,1);
    
    Uf.NumericalFluxGradient(upwindX,upwindY,1);
    
    AddTraceIntegral(upwindX,upwindZero,tmpX,0);
    AddTraceIntegral(upwindZero,upwindY,tmpY,0);
    
    MultiplyByElmtInvMass(tmpX,Uf.UpdateCoeffs(1),0);
    MultiplyByElmtInvMass(tmpY,Uf.UpdateCoeffs(2),0);
    
    Uf.BwdTrans(1);
    Uf.BwdTrans(2);
    
    
    // compute {\bf Lambda_20^2  = - g d H alpha2 {\bf c}
    
    Vmath::Vmul(nTotQuadPoints,GetPhys(0),1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
    Vmath::Vmul(nTotQuadPoints,GetPhys(0),1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);
    
     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     Vmath::Smul(nTotQuadPoints,-g*alpha2,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Smul(nTotQuadPoints,-g*alpha2,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     // due to move to rhs
     Vmath::Smul(nTotQuadPoints,-1.0,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Smul(nTotQuadPoints,-1.0,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     IProductWRTBase(Uf.GetPhys(1),tmpX,0);
     IProductWRTBase(Uf.GetPhys(2),tmpY,0);
     
     // Add to the RHS terms 
     
     Vmath::Vadd(nTotCoeffs,outX,1,tmpX,1,outX,1);
     Vmath::Vadd(nTotCoeffs,outY,1,tmpY,1,outY,1);
     //--------------------------------------

     

     //---------------------------------------------------
     // Add the 
     // -d^2 (1/6 + \alpha_2 - \alpha_1)\nabla(\nabla \cdot (H_t {\bf u}))
     // term to the rhs 

     
     // solve a =  \nabla \cdot (H_t {\bf u})
     // [we store a in Uf(1)]
     
     // get u and v 
     // can't get Vmath::Vdiv to work ...
     for (int j = 0; j < nTotQuadPoints; ++j)
       {
	 Uf.UpdatePhys(1)[j] = GetPhys(1)[j]/GetPhys(0)[j];
	 Uf.UpdatePhys(2)[j] = GetPhys(2)[j]/GetPhys(0)[j];
       }
     
     // get uHt and vHt
     Vmath::Vmul(nTotQuadPoints,Ht,1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,Ht,1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);
     
     IProductWRTDerivBase(0,Uf.GetPhys(1),tmpX,0);
     IProductWRTDerivBase(1,Uf.GetPhys(2),tmpY,0);
     
     Vmath::Neg(nTotCoeffs,tmpX,1);
     Vmath::Neg(nTotCoeffs,tmpY,1);
     
     Uf.NumericalFluxDivergence(upwindX,upwindY,1,2);
     
     AddTraceIntegral(upwindX,upwindY,tmpX,0);
     
     MultiplyByElmtInvMass(tmpX,Uf.UpdateCoeffs(1),0);
     
     Uf.BwdTrans(1);
     
     
     // solve {\bf c} =  \nabla a
     // [we store c in Uf(1) and Uf(2)]
     
     IProductWRTDerivBase(0,Uf.GetPhys(1),tmpX,0);
     IProductWRTDerivBase(1,Uf.GetPhys(1),tmpY,0);
     
     Vmath::Neg(nTotCoeffs,tmpX,1);
     Vmath::Neg(nTotCoeffs,tmpY,1);
    
     Uf.NumericalFluxGradient(upwindX,upwindY,1);
     
     AddTraceIntegral(upwindX,upwindZero,tmpX,0);
     AddTraceIntegral(upwindZero,upwindY,tmpY,0);
     
     MultiplyByElmtInvMass(tmpX,Uf.UpdateCoeffs(1),0);
     MultiplyByElmtInvMass(tmpY,Uf.UpdateCoeffs(2),0);
     
     Uf.BwdTrans(1);
     Uf.BwdTrans(2);
    

     // compute \Lambda_20^3 =  -d^2 (1/6 + \alpha_2 - \alpha_1) {\bf c}
    
     Vmath::Smul(nTotQuadPoints,-((1.0/6.0)+alpha2-alpha1),Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Smul(nTotQuadPoints,-((1.0/6.0)+alpha2-alpha1),Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     // due to move to rhs
     Vmath::Smul(nTotQuadPoints,-1.0,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Smul(nTotQuadPoints,-1.0,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);
     
     IProductWRTBase(Uf.GetPhys(1),tmpX,0);
     IProductWRTBase(Uf.GetPhys(2),tmpY,0);
     
     
     // Add to the RHS terms 
     
     Vmath::Vadd(nTotCoeffs,outX,1,tmpX,1,outX,1);
     Vmath::Vadd(nTotCoeffs,outY,1,tmpY,1,outY,1);
     //--------------------------------------



     //---------------------------------------------------
     // Add the 
     // d (1/2 + \alpha_2)\nabla(\nabla \cdot (H_t {\bf u}))
     // term to the rhs 

     
     // solve a =  \nabla \cdot (H_t d {\bf u})
     // [we store a in Uf(1)]
     
     // get u and v 
     // can't get Vmath::Vdiv to work ...
     for (int j = 0; j < nTotQuadPoints; ++j)
       {
	 Uf.UpdatePhys(1)[j] = GetPhys(1)[j]/GetPhys(0)[j];
	 Uf.UpdatePhys(2)[j] = GetPhys(2)[j]/GetPhys(0)[j];
       }
     
     // get u d Ht and v d Ht
     Vmath::Vmul(nTotQuadPoints,Ht,1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,Ht,1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);
     
     IProductWRTDerivBase(0,Uf.GetPhys(1),tmpX,0);
     IProductWRTDerivBase(1,Uf.GetPhys(2),tmpY,0);
     
     Vmath::Neg(nTotCoeffs,tmpX,1);
     Vmath::Neg(nTotCoeffs,tmpY,1);
     
     Uf.NumericalFluxDivergence(upwindX,upwindY,1,2);
     
     AddTraceIntegral(upwindX,upwindY,tmpX,0);
     
     MultiplyByElmtInvMass(tmpX,Uf.UpdateCoeffs(1),0);
     
     Uf.BwdTrans(1);
     
     
     // solve {\bf c} =  \nabla a
     // [we store c in Uf(1) and Uf(2)]
     
     IProductWRTDerivBase(0,Uf.GetPhys(1),tmpX,0);
     IProductWRTDerivBase(1,Uf.GetPhys(1),tmpY,0);
     
     Vmath::Neg(nTotCoeffs,tmpX,1);
     Vmath::Neg(nTotCoeffs,tmpY,1);
    
     Uf.NumericalFluxGradient(upwindX,upwindY,1);
     
     AddTraceIntegral(upwindX,upwindZero,tmpX,0);
     AddTraceIntegral(upwindZero,upwindY,tmpY,0);
     
     MultiplyByElmtInvMass(tmpX,Uf.UpdateCoeffs(1),0);
     MultiplyByElmtInvMass(tmpY,Uf.UpdateCoeffs(2),0);
     
     Uf.BwdTrans(1);
     Uf.BwdTrans(2);
    

     // compute \Lambda_20^4 =  d (1/2 + \alpha_2) {\bf c}
    
     Vmath::Smul(nTotQuadPoints,((1.0/2.0)+alpha2),Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Smul(nTotQuadPoints,((1.0/2.0)+alpha2),Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     // due to move to rhs
     Vmath::Smul(nTotQuadPoints,-1.0,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Smul(nTotQuadPoints,-1.0,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);
     
     IProductWRTBase(Uf.GetPhys(1),tmpX,0);
     IProductWRTBase(Uf.GetPhys(2),tmpY,0);
     
     
     // Add to the RHS terms 
     
     Vmath::Vadd(nTotCoeffs,outX,1,tmpX,1,outX,1);
     Vmath::Vadd(nTotCoeffs,outY,1,tmpY,1,outY,1);
     //--------------------------------------


     
     //---------------------------------------------------
     // Add the 
     // -d^2 (1/6 + \alpha_2 - \alpha_1)\nabla(\nabla H \cdot {\bf u}_t)
     // term to the rhs 
   

     // {\bf b} = \nabla H
     // [we store b in Uf(1) and Uf(2)]
     IProductWRTDerivBase(0,GetPhys(0),tmpX,0);
     IProductWRTDerivBase(1,GetPhys(0),tmpY,0);
     
     Vmath::Neg(nTotCoeffs,tmpX,1);
     Vmath::Neg(nTotCoeffs,tmpY,1);
    
     NumericalFluxGradient(upwindX,upwindY,0);
     
     AddTraceIntegral(upwindX,upwindZero,tmpX,0);
     AddTraceIntegral(upwindZero,upwindY,tmpY,0);
     
     MultiplyByElmtInvMass(tmpX,Uf.UpdateCoeffs(1),0);
     MultiplyByElmtInvMass(tmpY,Uf.UpdateCoeffs(2),0);
     
     Uf.BwdTrans(1);
     Uf.BwdTrans(2);

     
     // a = {\bf b} \cdot {\bf u}_t
     // [we store a in Uf(1)]
     
     Vmath::Vmul(nTotQuadPoints,GetPhys(5),1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,GetPhys(6),1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     Vmath::Vadd(nTotQuadPoints,Uf.GetPhys(1),1,Uf.GetPhys(2),1,Uf.UpdatePhys(1),1);

     // {\bf c} = nabla a
     // [we store c in Uf(1) and Uf(2)]
     IProductWRTDerivBase(0,Uf.GetPhys(1),tmpX,0);
     IProductWRTDerivBase(1,Uf.GetPhys(1),tmpY,0);
     
     Vmath::Neg(nTotCoeffs,tmpX,1);
     Vmath::Neg(nTotCoeffs,tmpY,1);
    
     Uf.NumericalFluxGradient(upwindX,upwindY,1);
     
     AddTraceIntegral(upwindX,upwindZero,tmpX,0);
     AddTraceIntegral(upwindZero,upwindY,tmpY,0);
     
     MultiplyByElmtInvMass(tmpX,Uf.UpdateCoeffs(1),0);
     MultiplyByElmtInvMass(tmpY,Uf.UpdateCoeffs(2),0);
     
     Uf.BwdTrans(1);
     Uf.BwdTrans(2);
     
     // compute \Lambda_{20}^5 = -d^2 (1/6 + \alpha_2 - \alpha_1){\bf c}
    
     Vmath::Smul(nTotQuadPoints,-((1.0/6.0)+alpha2-alpha1),Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Smul(nTotQuadPoints,-((1.0/6.0)+alpha2-alpha1),Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);
 
     // due to move to rhs
     Vmath::Smul(nTotQuadPoints,-1.0,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Smul(nTotQuadPoints,-1.0,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);
     
     IProductWRTBase(Uf.GetPhys(1),tmpX,0);
     IProductWRTBase(Uf.GetPhys(2),tmpY,0);
     
     
     // Add to the RHS terms 
     
     Vmath::Vadd(nTotCoeffs,outX,1,tmpX,1,outX,1);
     Vmath::Vadd(nTotCoeffs,outY,1,tmpY,1,outY,1);
     //--------------------------------------


     //---------------------------------------------------
     // Add the 
     // d (1/2 + \alpha_2)\nabla(\nabla H \cdot (d{\bf u}_t))
     // term to the rhs 
   

     // {\bf b} = \nabla H
     // [we store b in Uf(1) and Uf(2)]
     IProductWRTDerivBase(0,GetPhys(0),tmpX,0);
     IProductWRTDerivBase(1,GetPhys(0),tmpY,0);
     
     Vmath::Neg(nTotCoeffs,tmpX,1);
     Vmath::Neg(nTotCoeffs,tmpY,1);
    
     NumericalFluxGradient(upwindX,upwindY,0);
     
     AddTraceIntegral(upwindX,upwindZero,tmpX,0);
     AddTraceIntegral(upwindZero,upwindY,tmpY,0);
     
     MultiplyByElmtInvMass(tmpX,Uf.UpdateCoeffs(1),0);
     MultiplyByElmtInvMass(tmpY,Uf.UpdateCoeffs(2),0);
     
     Uf.BwdTrans(1);
     Uf.BwdTrans(2);

     
     // a = {\bf b} \cdot (d{\bf u}_t)
     // [we store a in Uf(1)]
     
     Vmath::Vmul(nTotQuadPoints,GetPhys(5),1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,GetPhys(6),1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     Vmath::Vadd(nTotQuadPoints,Uf.GetPhys(1),1,Uf.GetPhys(2),1,Uf.UpdatePhys(1),1);

     // {\bf c} = nabla a
     // [we store c in Uf(1) and Uf(2)]
     IProductWRTDerivBase(0,Uf.GetPhys(1),tmpX,0);
     IProductWRTDerivBase(1,Uf.GetPhys(1),tmpY,0);
     
     Vmath::Neg(nTotCoeffs,tmpX,1);
     Vmath::Neg(nTotCoeffs,tmpY,1);
    
     Uf.NumericalFluxGradient(upwindX,upwindY,1);
     
     AddTraceIntegral(upwindX,upwindZero,tmpX,0);
     AddTraceIntegral(upwindZero,upwindY,tmpY,0);
     
     MultiplyByElmtInvMass(tmpX,Uf.UpdateCoeffs(1),0);
     MultiplyByElmtInvMass(tmpY,Uf.UpdateCoeffs(2),0);
     
     Uf.BwdTrans(1);
     Uf.BwdTrans(2);
     
     // compute \Lambda_{20}^6 = d (1/2 + \alpha_2){\bf c}
    
     Vmath::Smul(nTotQuadPoints,((1.0/2.0)+alpha2),Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Smul(nTotQuadPoints,((1.0/2.0)+alpha2),Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);
 
     // due to move to rhs
     Vmath::Smul(nTotQuadPoints,-1.0,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Smul(nTotQuadPoints,-1.0,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);
     
     IProductWRTBase(Uf.GetPhys(1),tmpX,0);
     IProductWRTBase(Uf.GetPhys(2),tmpY,0);
     
     
     // Add to the RHS terms 
     
     Vmath::Vadd(nTotCoeffs,outX,1,tmpX,1,outX,1);
     Vmath::Vadd(nTotCoeffs,outY,1,tmpY,1,outY,1);
     //--------------------------------------


         
     //---------------------------------------------------
     // Add the 
     // -d^2 (1/6 + \alpha_2 - \alpha_1)\nabla H (\nabla \cdot {\bf u}_t)
     // term to the rhs 
   

     // {\bf b} = \nabla H
     // [we store b in Uf(1) and Uf(2)]
     IProductWRTDerivBase(0,GetPhys(0),tmpX,0);
     IProductWRTDerivBase(1,GetPhys(0),tmpY,0);
     
     Vmath::Neg(nTotCoeffs,tmpX,1);
     Vmath::Neg(nTotCoeffs,tmpY,1);
    
     NumericalFluxGradient(upwindX,upwindY,0);
     
     AddTraceIntegral(upwindX,upwindZero,tmpX,0);
     AddTraceIntegral(upwindZero,upwindY,tmpY,0);
     
     MultiplyByElmtInvMass(tmpX,Uf.UpdateCoeffs(1),0);
     MultiplyByElmtInvMass(tmpY,Uf.UpdateCoeffs(2),0);
     
     Uf.BwdTrans(1);
     Uf.BwdTrans(2);

     // c = \nabla \cdot {\bf u}_t
     // [we store c in Uf(3)]
     
     IProductWRTDerivBase(0,GetPhys(5),tmpX,0);
     IProductWRTDerivBase(1,GetPhys(6),tmpY,0);

     Vmath::Vadd(nTotCoeffs,tmpX,1,tmpY,1,tmpX,1);
     Vmath::Neg(nTotCoeffs,tmpX,1);
     
     NumericalFluxDivergence(upwindX,upwindY,5,6);
     
     AddTraceIntegral(upwindX,upwindY,tmpX,0);
     
     MultiplyByElmtInvMass(tmpX,Uf.UpdateCoeffs(3),0);
     
     Uf.BwdTrans(3);
     
     
     // compute Lambda_{20}^7 = - d^2(1/6+alpha2-alpha1){\bf b} c

     Vmath::Vmul(nTotQuadPoints,Uf.GetPhys(1),1,Uf.GetPhys(3),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,Uf.GetPhys(2),1,Uf.GetPhys(3),1,Uf.UpdatePhys(2),1);
     
     Vmath::Smul(nTotQuadPoints,-((1.0/6.0)+alpha2-alpha1),Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Smul(nTotQuadPoints,-((1.0/6.0)+alpha2-alpha1),Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);
     
     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     // due to move to rhs
     Vmath::Smul(nTotQuadPoints,-1.0,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Smul(nTotQuadPoints,-1.0,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);
     
     IProductWRTBase(Uf.GetPhys(1),tmpX,0);
     IProductWRTBase(Uf.GetPhys(2),tmpY,0);
     
     
     // Add to the RHS terms 
     
     Vmath::Vadd(nTotCoeffs,outX,1,tmpX,1,outX,1);
     Vmath::Vadd(nTotCoeffs,outY,1,tmpY,1,outY,1);
     //--------------------------------------


     //---------------------------------------------------
     // Add the 
     // d (1/2 + \alpha_2)\nabla H (\nabla \cdot (d{\bf u}_t))
     // term to the rhs 
   

     // {\bf b} = \nabla H
     // [we store b in Uf(1) and Uf(2)]
     IProductWRTDerivBase(0,GetPhys(0),tmpX,0);
     IProductWRTDerivBase(1,GetPhys(0),tmpY,0);
     
     Vmath::Neg(nTotCoeffs,tmpX,1);
     Vmath::Neg(nTotCoeffs,tmpY,1);
    
     NumericalFluxGradient(upwindX,upwindY,0);
     
     AddTraceIntegral(upwindX,upwindZero,tmpX,0);
     AddTraceIntegral(upwindZero,upwindY,tmpY,0);
     
     MultiplyByElmtInvMass(tmpX,Uf.UpdateCoeffs(1),0);
     MultiplyByElmtInvMass(tmpY,Uf.UpdateCoeffs(2),0);
     
     Uf.BwdTrans(1);
     Uf.BwdTrans(2);

     // c = \nabla \cdot (d{\bf u}_t)
     // [we store c in Uf(3)]

     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,GetPhys(5),1,Uf.UpdatePhys(3),1);
     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,GetPhys(6),1,Uf.UpdatePhys(4),1);
     
     IProductWRTDerivBase(0,Uf.GetPhys(3),tmpX,0);
     IProductWRTDerivBase(1,Uf.GetPhys(4),tmpY,0);

     Vmath::Vadd(nTotCoeffs,tmpX,1,tmpY,1,tmpX,1);
     Vmath::Neg(nTotCoeffs,tmpX,1);
     
     Uf.NumericalFluxDivergence(upwindX,upwindY,3,4);
     
     AddTraceIntegral(upwindX,upwindY,tmpX,0);
     
     MultiplyByElmtInvMass(tmpX,Uf.UpdateCoeffs(3),0);
     
     Uf.BwdTrans(3);
     
     
     // compute Lambda_{20}^8 = d(1/2+alpha2){\bf b} c

     Vmath::Vmul(nTotQuadPoints,Uf.GetPhys(1),1,Uf.GetPhys(3),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,Uf.GetPhys(2),1,Uf.GetPhys(3),1,Uf.UpdatePhys(2),1);
     
     Vmath::Smul(nTotQuadPoints,((1.0/2.0)+alpha2),Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Smul(nTotQuadPoints,((1.0/2.0)+alpha2),Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);
     
     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);

     // due to move to rhs
     Vmath::Smul(nTotQuadPoints,-1.0,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
     Vmath::Smul(nTotQuadPoints,-1.0,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);
     
     IProductWRTBase(Uf.GetPhys(1),tmpX,0);
     IProductWRTBase(Uf.GetPhys(2),tmpY,0);
     
     
     // Add to the RHS terms 
     
     Vmath::Vadd(nTotCoeffs,outX,1,tmpX,1,outX,1);
     Vmath::Vadd(nTotCoeffs,outY,1,tmpY,1,outY,1);
     //--------------------------------------

     
  }

  void BoussinesqEquations::FullyNonlinearSpatial(BoussinesqEquations &Uf, Array<OneD, NekDouble> &outX, 
						  Array<OneD, NekDouble> &outY, Array<OneD, NekDouble> &Ht)
  {
    //--------------------------------------
    // local parameters
    int nTotQuadPoints = GetPointsTot();
    int nTotCoeffs     = GetNcoeffs();
    int nTracePoints   = GetNpoints();
    
    NekDouble g        = GetGravity();
    NekDouble alpha1   = GetAlpha1();
    NekDouble alpha2   = GetAlpha2();
    //--------------------------------------
  
    Array<OneD, NekDouble> tmpX(GetNcoeffs());
    Array<OneD, NekDouble> tmpY(GetNcoeffs());
    
    Array<OneD, NekDouble> upwindX(GetNpoints());
    Array<OneD, NekDouble> upwindY(GetNpoints());
    Array<OneD, NekDouble> upwindZero(GetNpoints(),0.0);
    
    Array<OneD, NekDouble> eta(nTotQuadPoints);
    Array<OneD, NekDouble> u(nTotQuadPoints);
    Array<OneD, NekDouble> v(nTotQuadPoints);
    Array<OneD, NekDouble> physX(nTotQuadPoints);
    Array<OneD, NekDouble> physY(nTotQuadPoints);
    
    //--------------------------------------
    // Get primitive variables
    for (int j = 0; j < nTotQuadPoints; ++j)
      {
	eta[j] = GetPhys(0)[j]-GetPhys(4)[j];
	u[j]   = GetPhys(1)[j]/GetPhys(0)[j];
	v[j]   = GetPhys(2)[j]/GetPhys(0)[j];
      }
    //--------------------------------------
    
    
    //----------------------------------------
    // compute a_1 = \nabla \cdot {\bf u}
    
    // hey hey ... those boundary conditions...
    Array<OneD, NekDouble > a1(nTotQuadPoints);
    Vmath::Vcopy(nTotQuadPoints,u,1,Uf.UpdatePhys(1),1);
    Vmath::Vcopy(nTotQuadPoints,v,1,Uf.UpdatePhys(2),1);
    Uf.DivergenceFluxTerms(Uf.GetPhys(1),Uf.GetPhys(2),a1,1,2);
    //----------------------------------------
    

    //----------------------------------------
    // compute a_2 = \nabla \cdot {\bf u}_t
     
    // hey hey ... those boundary conditions...
    Array<OneD, NekDouble > a2(nTotQuadPoints);
    Vmath::Vcopy(nTotQuadPoints,GetPhys(5),1,Uf.UpdatePhys(1),1);
    Vmath::Vcopy(nTotQuadPoints,GetPhys(6),1,Uf.UpdatePhys(2),1);
    Uf.DivergenceFluxTerms(Uf.GetPhys(1),Uf.GetPhys(2),a2,1,2);
    //----------------------------------------


    //----------------------------------------
    // compute a_3 = \nabla \cdot (h{\bf u})
     
    // hey hey ... those boundary conditions...
    Array<OneD, NekDouble > a3(nTotQuadPoints);
    Vmath::Vmul(nTotQuadPoints,u,1,GetPhys(4),1,Uf.UpdatePhys(1),1);
    Vmath::Vmul(nTotQuadPoints,v,1,GetPhys(4),1,Uf.UpdatePhys(2),1);
    Uf.DivergenceFluxTerms(Uf.GetPhys(1),Uf.GetPhys(2),a3,1,2);
    //----------------------------------------

    
    //----------------------------------------
    // compute a_4 = \nabla \cdot (h {\bf u}_t)
     
    // hey hey ... those boundary conditions...
    Array<OneD, NekDouble > a4(nTotQuadPoints);
    Vmath::Vmul(nTotQuadPoints,GetPhys(5),1,GetPhys(4),1,Uf.UpdatePhys(1),1);
    Vmath::Vmul(nTotQuadPoints,GetPhys(6),1,GetPhys(4),1,Uf.UpdatePhys(2),1);
    Uf.DivergenceFluxTerms(Uf.GetPhys(1),Uf.GetPhys(2),a4,1,2);
    //----------------------------------------


    //----------------------------------------
    // compute a_5 = \nabla \cdot (H_t {\bf u})
     
    // hey hey ... those boundary conditions...
    Array<OneD, NekDouble > a5(nTotQuadPoints);
    Vmath::Vmul(nTotQuadPoints,Ht,1,u,1,Uf.UpdatePhys(1),1);
    Vmath::Vmul(nTotQuadPoints,Ht,1,v,1,Uf.UpdatePhys(2),1);
    Uf.DivergenceFluxTerms(Uf.GetPhys(1),Uf.GetPhys(2),a5,1,2);
    //----------------------------------------


    //----------------------------------------
    // compute a_6 = \nabla \cdot (h H_t {\bf u})
     
    // hey hey ... those boundary conditions...
    Array<OneD, NekDouble > a6(nTotQuadPoints);

    Vmath::Vmul(nTotQuadPoints,Ht,1,u,1,Uf.UpdatePhys(1),1);
    Vmath::Vmul(nTotQuadPoints,Ht,1,v,1,Uf.UpdatePhys(2),1);
    Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
    Vmath::Vmul(nTotQuadPoints,GetPhys(4),1,Uf.GetPhys(2),1,Uf.UpdatePhys(2),1);
    Uf.DivergenceFluxTerms(Uf.GetPhys(1),Uf.GetPhys(2),a6,1,2);
    //----------------------------------------


    //----------------------------------------
    // compute a_7 = \nabla \cdot (\eta {\bf u})
     
    // hey hey ... those boundary conditions...
    Array<OneD, NekDouble > a7(nTotQuadPoints);

    Vmath::Vmul(nTotQuadPoints,eta,1,u,1,Uf.UpdatePhys(1),1);
    Vmath::Vmul(nTotQuadPoints,eta,1,v,1,Uf.UpdatePhys(2),1);
    Uf.DivergenceFluxTerms(Uf.GetPhys(1),Uf.GetPhys(2),a7,1,2);
    //----------------------------------------

    
    //----------------------------------------
    // compute {\bf b}_1 = \nabla a_1
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > b1(2);
    b1[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    b1[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    Vmath::Vcopy(nTotQuadPoints,a1,1,Uf.UpdatePhys(1),1);
    Uf.GradientFluxTerms(Uf.GetPhys(1),b1,1);
    //----------------------------------------
    
    
    //----------------------------------------
    // compute {\bf b}_2 = \nabla a_2
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > b2(2);
    b2[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    b2[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    Vmath::Vcopy(nTotQuadPoints,a2,1,Uf.UpdatePhys(1),1);
    Uf.GradientFluxTerms(Uf.GetPhys(1),b2,1);
    //----------------------------------------

    
    //----------------------------------------
    // compute {\bf b}_3 = \nabla a_3
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > b3(2);
    b3[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    b3[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    Vmath::Vcopy(nTotQuadPoints,a3,1,Uf.UpdatePhys(1),1);
    Uf.GradientFluxTerms(Uf.GetPhys(1),b3,1);
    //----------------------------------------


    //----------------------------------------
    // compute {\bf b}_4 = \nabla a_4
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > b4(2);
    b4[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    b4[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    Vmath::Vcopy(nTotQuadPoints,a4,1,Uf.UpdatePhys(1),1);
    Uf.GradientFluxTerms(Uf.GetPhys(1),b4,1);
    //----------------------------------------


    //----------------------------------------
    // compute \Gamma   = (1/6)h{\bf b}_1-(1/2){\bf b}_3
    // compute \Gamma_t = (1/6)h{\bf b}_2-(1/2){\bf b}_4
    
    Array<OneD, Array<OneD, NekDouble> > Gamma(2);
    Gamma[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    Gamma[1] = Array<OneD, NekDouble> (nTotQuadPoints);

    Array<OneD, Array<OneD, NekDouble> > Gammat(2);
    Gammat[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    Gammat[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    
    for (int i = 0; i < nTotQuadPoints; ++i)
      {
	Gamma[0][i] = (1.0/6.0)*(GetPhys(4))[i]*b1[0][i] -
	  0.5 *b3[0][i];
	Gamma[1][i] = (1.0/6.0)*(GetPhys(4))[i]*b1[1][i] -
	  0.5 *b3[1][i];
	Gammat[0][i] = (1.0/6.0)*(GetPhys(4))[i]*b2[0][i] -
	  0.5 *b4[0][i];
	Gammat[1][i] = (1.0/6.0)*(GetPhys(4))[i]*b2[1][i] -
	  0.5 *b4[1][i];
      }
    //----------------------------------------
	

    //----------------------------------------
    // compute {\bf b}_5 = \nabla a_5
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > b5(2);
    b5[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    b5[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    Vmath::Vcopy(nTotQuadPoints,a5,1,Uf.UpdatePhys(1),1);
    Uf.GradientFluxTerms(Uf.GetPhys(1),b5,1);
    //----------------------------------------


    //----------------------------------------
    // compute {\bf b}_6 = \nabla a_6
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > b6(2);
    b6[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    b6[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    Vmath::Vcopy(nTotQuadPoints,a6,1,Uf.UpdatePhys(1),1);
    Uf.GradientFluxTerms(Uf.GetPhys(1),b6,1);
    //----------------------------------------

    
    //----------------------------------------
    // compute {\bf b}_7 = \nabla a_7
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > b7(2);
    b7[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    b7[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    Vmath::Vcopy(nTotQuadPoints,a7,1,Uf.UpdatePhys(1),1);
    Uf.GradientFluxTerms(Uf.GetPhys(1),b7,1);
    //----------------------------------------


    //----------------------------------------
    // compute {\bf b}_8 = \nabla (\eta a_4)
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > b8(2);
    b8[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    b8[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    Vmath::Vmul(nTotQuadPoints,eta,1,a4,1,Uf.UpdatePhys(1),1);
    Uf.GradientFluxTerms(Uf.GetPhys(1),b8,1);
    //----------------------------------------


    //----------------------------------------
    // compute {\bf b}_9 = \nabla (a3 a_3)
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > b9(2);
    b9[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    b9[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    Vmath::Vmul(nTotQuadPoints,a3,1,a3,1,Uf.UpdatePhys(1),1);
    Uf.GradientFluxTerms(Uf.GetPhys(1),b9,1);
    //----------------------------------------

    
    //----------------------------------------
    // compute {\bf b}_10 = \nabla (\eta \eta a_2)
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > b10(2);
    b10[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    b10[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    Vmath::Vmul(nTotQuadPoints,eta,1,eta,1,Uf.UpdatePhys(1),1);
    Vmath::Vmul(nTotQuadPoints,a2,1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
    Uf.GradientFluxTerms(Uf.GetPhys(1),b10,1);
    //----------------------------------------


    //----------------------------------------
    // compute {\bf b}_11 = \nabla (\eta a_1 a_3)
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > b11(2);
    b11[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    b11[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    Vmath::Vmul(nTotQuadPoints,eta,1,a1,1,Uf.UpdatePhys(1),1);
    Vmath::Vmul(nTotQuadPoints,a3,1,Uf.GetPhys(1),1,Uf.UpdatePhys(1),1);
    Uf.GradientFluxTerms(Uf.GetPhys(1),b11,1);
    //----------------------------------------

    
    //----------------------------------------
    // compute {\bf b}_12 = \nabla ({\bf u} \cdot (h \Gamma))
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > b12(2);
    b12[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    b12[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    for (int i = 0; i < nTotQuadPoints; ++i)
      {
	Uf.UpdatePhys(1)[i] = GetPhys(4)[i]*(u[i]*Gamma[0][i] + 
					     v[i]*Gamma[1][i]);
      }
    Uf.GradientFluxTerms(Uf.GetPhys(1),b12,1);
    //----------------------------------------

    
    //----------------------------------------
    // compute {\bf b}_13 = \nabla (\eta {\bf u} \cdot \Gamma)
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > b13(2);
    b13[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    b13[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    for (int i = 0; i < nTotQuadPoints; ++i)
      {
	Uf.UpdatePhys(1)[i] = eta[i]*(u[i]*Gamma[0][i] + 
				      v[i]*Gamma[1][i]);
      }
    Uf.GradientFluxTerms(Uf.GetPhys(1),b13,1);
    //----------------------------------------
    
    
    //----------------------------------------
    // compute {\bf b}_14 = \nabla (\eta {\bf u} \cdot {\bf b}_3)
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > b14(2);
    b14[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    b14[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    for (int i = 0; i < nTotQuadPoints; ++i)
      {
	Uf.UpdatePhys(1)[i] = eta[i]*(u[i]*b3[0][i] + 
				      v[i]*b3[1][i]);
      }
    Uf.GradientFluxTerms(Uf.GetPhys(1),b14,1);
    //----------------------------------------


    //----------------------------------------
    // compute {\bf b}_15 = \nabla (\eta \eta {\bf u} \cdot {\bf b}_1)
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > b15(2);
    b15[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    b15[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    for (int i = 0; i < nTotQuadPoints; ++i)
      {
	Uf.UpdatePhys(1)[i] = eta[i]* eta[i] * (u[i]*b1[0][i] + 
						v[i]*b1[1][i]);
      }
    Uf.GradientFluxTerms(Uf.GetPhys(1),b15,1);
    //----------------------------------------

    
    //----------------------------------------
    // compute {\bf b}_16 = \nabla (\eta \eta a_1 a_1)
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > b16(2);
    b16[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    b16[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    for (int i = 0; i < nTotQuadPoints; ++i)
      {
	Uf.UpdatePhys(1)[i] = eta[i]* eta[i] * a1[i] * a1[i];	
      }
    Uf.GradientFluxTerms(Uf.GetPhys(1),b16,1);
    //----------------------------------------


    //----------------------------------------
    // compute {\bf c}_1 = \nabla H
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > c1(2);
    c1[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    c1[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    Vmath::Vcopy(nTotQuadPoints,GetPhys(0),1,Uf.UpdatePhys(1),1);
    Uf.GradientFluxTerms(Uf.GetPhys(1),c1,1);
    //----------------------------------------

    
    //----------------------------------------
    // compute {\bf l}_1 = \nabla ({\bf c}_1 \cdot {\bf u}_t)
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > l1(2);
    l1[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    l1[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    for (int i = 0; i < nTotQuadPoints; ++i)
      {
	Uf.UpdatePhys(1)[i] = c1[0][i] * GetPhys(5)[i] +
	  c1[1][i] * GetPhys(6)[i];
      }
    Uf.GradientFluxTerms(Uf.GetPhys(1),l1,1);
    //----------------------------------------
    
    
    //----------------------------------------
    // compute {\bf l}_2 = \nabla ({\bf c}_1 \cdot (h{\bf u}_t))
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > l2(2);
    l2[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    l2[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    for (int i = 0; i < nTotQuadPoints; ++i)
      {
	Uf.UpdatePhys(1)[i] = GetPhys(4)[i] * (c1[0][i] * GetPhys(5)[i] +
					       c1[1][i] * GetPhys(6)[i]);
      }
    Uf.GradientFluxTerms(Uf.GetPhys(1),l2,1);
    //----------------------------------------
    

    //----------------------------------------
    // compute {\bf e}_1 = \nabla \eta
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > e1(2);
    e1[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    e1[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    Vmath::Vcopy(nTotQuadPoints,eta,1,Uf.UpdatePhys(1),1);
    Uf.GradientFluxTerms(Uf.GetPhys(1),e1,1);
    //----------------------------------------


    //----------------------------------------
    // compute g_1 = \nabla \cdot {\bf e}_1
    
    // hey hey ... those boundary conditions...
    Array<OneD, NekDouble>  g1(nTotQuadPoints); 
    Vmath::Vcopy(nTotQuadPoints,e1[0],1,Uf.UpdatePhys(1),1);
    Vmath::Vcopy(nTotQuadPoints,e1[1],1,Uf.UpdatePhys(2),1);
    Uf.DivergenceFluxTerms(Uf.GetPhys(1),Uf.GetPhys(2),g1,1,2);
    //----------------------------------------
    
    
    //----------------------------------------
    // compute g_2 = \nabla \cdot (h{\bf e}_1)
    
    // hey hey ... those boundary conditions...
    Array<OneD, NekDouble>  g2(nTotQuadPoints); 
    Vmath::Vmul(nTotQuadPoints,e1[0],1,GetPhys(4),1,Uf.UpdatePhys(1),1);
    Vmath::Vmul(nTotQuadPoints,e1[1],1,GetPhys(4),1,Uf.UpdatePhys(2),1);
    Uf.DivergenceFluxTerms(Uf.GetPhys(1),Uf.GetPhys(2),g2,1,2);
    //----------------------------------------


    //----------------------------------------
    // compute {\bf k}_1 = \nabla g1
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > k1(2);
    k1[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    k1[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    Vmath::Vcopy(nTotQuadPoints,g1,1,Uf.UpdatePhys(1),1);
    Uf.GradientFluxTerms(Uf.GetPhys(1),k1,1);
    //----------------------------------------
    
    
    //----------------------------------------
    // compute {\bf k}_2 = \nabla g2
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > k2(2);
    k2[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    k2[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    Vmath::Vcopy(nTotQuadPoints,g2,1,Uf.UpdatePhys(1),1);
    Uf.GradientFluxTerms(Uf.GetPhys(1),k2,1);
    //----------------------------------------


    //----------------------------------------
    // compute {\bf m}_1 = \nabla u
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > m1(2);
    m1[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    m1[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    Vmath::Vcopy(nTotQuadPoints,u,1,Uf.UpdatePhys(1),1);
    Uf.GradientFluxTerms(Uf.GetPhys(1),m1,1);
    //----------------------------------------
    

    //----------------------------------------
    // compute {\bf m}_2 = \nabla v
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > m2(2);
    m2[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    m2[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    Vmath::Vcopy(nTotQuadPoints,v,1,Uf.UpdatePhys(1),1);
    Uf.GradientFluxTerms(Uf.GetPhys(1),m2,1);
    //----------------------------------------


    //----------------------------------------
    // compute n_1 = \partial_x ({\bf u} \cdot {\bf m}_1)
    
    Array<OneD, NekDouble>  n1(nTotQuadPoints); 
    for (int i = 0; i < nTotQuadPoints; ++i)
      {
	Uf.UpdatePhys(1)[i] = u[i]*m1[0][i] + v[i]*m1[1][i];
      }
    
    IProductWRTDerivBase(0,Uf.GetPhys(1),tmpX,0);
    Vmath::Neg(nTotCoeffs,tmpX,1);
    Uf.NumericalFluxGradient(upwindX,upwindY,1);
    AddTraceIntegral(upwindX,upwindZero,tmpX,0);
    MultiplyByElmtInvMass(tmpX,tmpX,0);
    BwdTrans(tmpX,n1,1);
    //----------------------------------------

    
    //----------------------------------------
    // compute n_3 = \partial_x (h {\bf u} \cdot {\bf m}_1)
    
    Array<OneD, NekDouble>  n3(nTotQuadPoints); 
    for (int i = 0; i < nTotQuadPoints; ++i)
      {
	Uf.UpdatePhys(1)[i] = GetPhys(4)[i]*(u[i]*m1[0][i] + v[i]*m1[1][i]);
      }
    
    IProductWRTDerivBase(0,Uf.GetPhys(1),tmpX,0);
    Vmath::Neg(nTotCoeffs,tmpX,1);
    Uf.NumericalFluxGradient(upwindX,upwindY,1);
    AddTraceIntegral(upwindX,upwindZero,tmpX,0);
    MultiplyByElmtInvMass(tmpX,tmpX,0);
    BwdTrans(tmpX,n3,1);
    //----------------------------------------
   
    
    //----------------------------------------
    // compute n_2 = \partial_y ({\bf u} \cdot {\bf m}_2)
    
    Array<OneD, NekDouble>  n2(nTotQuadPoints); 
    for (int i = 0; i < nTotQuadPoints; ++i)
      {
	Uf.UpdatePhys(1)[i] = u[i]*m2[0][i] + v[i]*m2[1][i];
      }
    
    IProductWRTDerivBase(1,Uf.GetPhys(1),tmpY,0);
    Vmath::Neg(nTotCoeffs,tmpY,1);
    Uf.NumericalFluxGradient(upwindX,upwindY,1);
    AddTraceIntegral(upwindZero,upwindY,tmpY,0);
    MultiplyByElmtInvMass(tmpY,tmpY,0);
    BwdTrans(tmpY,n2,1);
    //----------------------------------------

    
    //----------------------------------------
    // compute n_4 = \partial_y (h{\bf u} \cdot {\bf m}_2)
    
    Array<OneD, NekDouble>  n4(nTotQuadPoints); 
    for (int i = 0; i < nTotQuadPoints; ++i)
      {
	Uf.UpdatePhys(1)[i] = GetPhys(4)[i] *(u[i]*m2[0][i] + v[i]*m2[1][i]);
      }
    
    IProductWRTDerivBase(1,Uf.GetPhys(1),tmpY,0);
    Vmath::Neg(nTotCoeffs,tmpY,1);
    Uf.NumericalFluxGradient(upwindX,upwindY,1);
    AddTraceIntegral(upwindZero,upwindY,tmpY,0);
    MultiplyByElmtInvMass(tmpY,tmpY,0);
    BwdTrans(tmpY,n4,1);
    //----------------------------------------


    //----------------------------------------
    // compute {\bf p}_1 = \nabla (n_1+n_2)
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > p1(2);
    p1[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    p1[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    Vmath::Vadd(nTotQuadPoints,n1,1,n2,1,Uf.UpdatePhys(1),1);
    Uf.GradientFluxTerms(Uf.GetPhys(1),p1,1);
    //----------------------------------------

    
    //----------------------------------------
    // compute {\bf p}_2 = \nabla (n_3+n_4)
    
    // hey hey ... those boundary conditions...
    Array<OneD, Array<OneD, NekDouble> > p2(2);
    p2[0] = Array<OneD, NekDouble> (nTotQuadPoints); 
    p2[1] = Array<OneD, NekDouble> (nTotQuadPoints);
    Vmath::Vadd(nTotQuadPoints,n3,1,n4,1,Uf.UpdatePhys(1),1);
    Uf.GradientFluxTerms(Uf.GetPhys(1),p2,1);
    //----------------------------------------


    //----------------------------------------
    // compute q_1 = \partial_y (h \Gamma^(x))
    
    Array<OneD, NekDouble>  q1(nTotQuadPoints); 
    for (int i = 0; i < nTotQuadPoints; ++i)
      {
	Uf.UpdatePhys(1)[i] = GetPhys(4)[i]*Gamma[0][i];
      }
    
    IProductWRTDerivBase(1,Uf.GetPhys(1),tmpY,0);
    Vmath::Neg(nTotCoeffs,tmpY,1);
    Uf.NumericalFluxGradient(upwindX,upwindY,1);
    AddTraceIntegral(upwindZero,upwindY,tmpY,0);
    MultiplyByElmtInvMass(tmpY,tmpY,0);
    BwdTrans(tmpY,q1,1);
    //----------------------------------------

    
    //----------------------------------------
    // compute q_2 = \partial_x (h \Gamma^(y))
    
    Array<OneD, NekDouble>  q2(nTotQuadPoints); 
    for (int i = 0; i < nTotQuadPoints; ++i)
      {
	Uf.UpdatePhys(1)[i] = GetPhys(4)[i]*Gamma[1][i];
      }
    
    IProductWRTDerivBase(0,Uf.GetPhys(1),tmpX,0);
    Vmath::Neg(nTotCoeffs,tmpX,1);
    Uf.NumericalFluxGradient(upwindX,upwindY,1);
    AddTraceIntegral(upwindX,upwindZero,tmpX,0);
    MultiplyByElmtInvMass(tmpX,tmpX,0);
    BwdTrans(tmpX,q2,1);
    //----------------------------------------


    //----------------------------------------
    // compute q_3 = \partial_y (\eta \Gamma^(x))
    
    Array<OneD, NekDouble>  q3(nTotQuadPoints); 
    for (int i = 0; i < nTotQuadPoints; ++i)
      {
	Uf.UpdatePhys(1)[i] = eta[i]*Gamma[0][i];
      }
    
    IProductWRTDerivBase(1,Uf.GetPhys(1),tmpY,0);
    Vmath::Neg(nTotCoeffs,tmpY,1);
    Uf.NumericalFluxGradient(upwindX,upwindY,1);
    AddTraceIntegral(upwindZero,upwindY,tmpY,0);
    MultiplyByElmtInvMass(tmpY,tmpY,0);
    BwdTrans(tmpY,q3,1);
    //----------------------------------------

    
    //----------------------------------------
    // compute q_4 = \partial_x (\eta \Gamma^(y))
    
    Array<OneD, NekDouble>  q4(nTotQuadPoints); 
    for (int i = 0; i < nTotQuadPoints; ++i)
      {
	Uf.UpdatePhys(1)[i] = eta[i]*Gamma[1][i];
      }
    
    IProductWRTDerivBase(0,Uf.GetPhys(1),tmpX,0);
    Vmath::Neg(nTotCoeffs,tmpX,1);
    Uf.NumericalFluxGradient(upwindX,upwindY,1);
    AddTraceIntegral(upwindX,upwindZero,tmpX,0);
    MultiplyByElmtInvMass(tmpX,tmpX,0);
    BwdTrans(tmpX,q4,1);
    //----------------------------------------


   //----------------------------------------
    // compute q_5 = \partial_y (\eta \eta {\bf b}_1^(x))
    
    Array<OneD, NekDouble>  q5(nTotQuadPoints); 
    for (int i = 0; i < nTotQuadPoints; ++i)
      {
	Uf.UpdatePhys(1)[i] = eta[i]*eta[i]*b1[0][i];
      }
    
    IProductWRTDerivBase(1,Uf.GetPhys(1),tmpY,0);
    Vmath::Neg(nTotCoeffs,tmpY,1);
    Uf.NumericalFluxGradient(upwindX,upwindY,1);
    AddTraceIntegral(upwindZero,upwindY,tmpY,0);
    MultiplyByElmtInvMass(tmpY,tmpY,0);
    BwdTrans(tmpY,q5,1);
    //----------------------------------------

    
    //----------------------------------------
    // compute q_6 = \partial_x (\eta \eta {\bf b}_1^(y))
    
    Array<OneD, NekDouble>  q6(nTotQuadPoints); 
    for (int i = 0; i < nTotQuadPoints; ++i)
      {
	Uf.UpdatePhys(1)[i] = eta[i]*eta[i]*b1[1][i];
      }
    
    IProductWRTDerivBase(0,Uf.GetPhys(1),tmpX,0);
    Vmath::Neg(nTotCoeffs,tmpX,1);
    Uf.NumericalFluxGradient(upwindX,upwindY,1);
    AddTraceIntegral(upwindX,upwindZero,tmpX,0);
    MultiplyByElmtInvMass(tmpX,tmpX,0);
    BwdTrans(tmpX,q6,1);
    //----------------------------------------

    

    for (int i = 0; i < nTotQuadPoints; ++i)
      {
	// Lambda 20 terms (spatial terms that is)
	physX[i] = g*GetPhys(0)[i]*GetPhys(4)[i]*(GetPhys(4)[i]*(alpha2-alpha1)*k1[0][i] - alpha2*k2[0][i])
	  + GetPhys(4)[i]*GetPhys(4)[i]*((1.0/6.0)+alpha2-alpha1)*(-c1[0][i]*a2[i]-l1[0][i]-b5[0][i])
	  - GetPhys(4)[i]*(0.5+alpha2)*(-c1[0][i]*a4[i] - l2[0][i] - b6[0][i]);
	physY[i] = g*GetPhys(0)[i]*GetPhys(4)[i]*(GetPhys(4)[i]*(alpha2-alpha1)*k1[1][i] - alpha2*k2[1][i])
	  + GetPhys(4)[i]*GetPhys(4)[i]*((1.0/6.0)+alpha2-alpha1)*(-c1[1][i]*a2[i]-l1[1][i]-b5[1][i])
	  - GetPhys(4)[i]*(0.5+alpha2)*(-c1[1][i]*a4[i] - l2[1][i] - b6[1][i]);
	
 	// Lambda 21 terms 
 	physX[i] += GetPhys(0)[i]*(GetPhys(4)[i]*GetPhys(4)[i]*(alpha2-alpha1)*p1[0][i]-GetPhys(4)[i]*alpha2*p2[0][i]
				  -eta[i]*Gammat[0][i]+a3[i]*Gamma[0][i] + b12[0][i] - b8[0][i] + 0.5*b9[0][i]);
	physY[i] += GetPhys(0)[i]*(GetPhys(4)[i]*GetPhys(4)[i]*(alpha2-alpha1)*p1[1][i]-GetPhys(4)[i]*alpha2*p2[1][i]
				  -eta[i]*Gammat[1][i]+a3[i]*Gamma[1][i] + b12[1][i] - b8[1][i] + 0.5*b9[1][i]);

	// Lambda 21 curl terms (missing)
	physX[i] += -(-v[i]*q1[i]+v[i]*q2[i]);
	physX[i] += -( u[i]*q1[i]-u[i]*q2[i]);

	// Lambda 22 terms 
	physX[i] += GetPhys(0)[i]*((1.0/6.0)*eta[i]*eta[i]*b2[0][i]-(1.0/3.0)*eta[i]*a3[i]*b1[0][i]+a7[i]*Gamma[0][i]
				  -b13[0][i]-0.5*b10[0][i]+b11[0][i]-b14[0][i]);
	physY[i] += GetPhys(0)[i]*((1.0/6.0)*eta[i]*eta[i]*b2[1][i]-(1.0/3.0)*eta[i]*a3[i]*b1[1][i]+a7[i]*Gamma[1][i]
				  -b13[1][i]-0.5*b10[1][i]+b11[1][i]-b14[1][i]);

	// Lambda 22 curl terms (missing)
	physX[i] += (-v[i]*q3[i]+v[i]*q4[i]);
	physX[i] += ( u[i]*q3[i]-u[i]*q4[i]);

	// Lambda 23 terms 
	physX[i] += GetPhys(0)[i]*(-(1.0/3.0)*eta[i]*a7[i]*b1[0][i]-(1.0/3.0)*b15[0][i]+0.5*b16[0][i]);
	physY[i] += GetPhys(0)[i]*(-(1.0/3.0)*eta[i]*a7[i]*b1[1][i]-(1.0/3.0)*b15[1][i]+0.5*b16[1][i]);
	
	// Lambda 23 curl terms (missing)
	physX[i] += -(1.0/6.0)*(-v[i]*q5[i]+v[i]*q6[i]);
	physX[i] += -(1.0/6.0)*( u[i]*q5[i]-u[i]*q6[i]);
	
      }

    
    // Get modal values
    IProductWRTBase(physX,tmpX,0);
    IProductWRTBase(physY,tmpY,0);
    MultiplyByElmtInvMass(tmpX,tmpX,0);
    MultiplyByElmtInvMass(tmpY,tmpY,0);


    // add to rhs
    Vmath::Vsub(nTotCoeffs,outX,1,tmpX,1,outX,1);
    Vmath::Vsub(nTotCoeffs,outY,1,tmpY,1,outY,1);
    
    

  } 


  // Term should be in Uf.m_fields[1]

  void BoussinesqEquations::GradientFluxTerms(const Array<OneD, const NekDouble> &in, 
					      Array<OneD, Array<OneD, NekDouble> > &out, int field_no)
  {
     
    Array<OneD, NekDouble> tmpX(GetNcoeffs());
    Array<OneD, NekDouble> tmpY(GetNcoeffs());
    
    Array<OneD, NekDouble> upwindX(GetNpoints());
    Array<OneD, NekDouble> upwindY(GetNpoints());
    Array<OneD, NekDouble> upwindZero(GetNpoints(),0.0);

    IProductWRTDerivBase(0,GetPhys(field_no),tmpX,0);
    IProductWRTDerivBase(1,GetPhys(field_no),tmpY,0);
    
    Vmath::Neg(GetNcoeffs(),tmpX,1);
    Vmath::Neg(GetNcoeffs(),tmpY,1);

    NumericalFluxGradient(upwindX,upwindY,field_no);
    
    AddTraceIntegral(upwindX,upwindZero,tmpX,0);
    AddTraceIntegral(upwindZero,upwindY,tmpY,0);
    
    MultiplyByElmtInvMass(tmpX,tmpX,0);
    MultiplyByElmtInvMass(tmpY,tmpY,0);
    
    BwdTrans(tmpX,out[0],1);
    BwdTrans(tmpY,out[1],2);
  }
  
  void BoussinesqEquations::DivergenceFluxTerms(const Array<OneD, const NekDouble> &inX, 
						const Array<OneD, const NekDouble> &inY, 
						Array<OneD, NekDouble> &out,  int field_no_1, int field_no_2)
  {
       
    Array<OneD, NekDouble> tmpX(GetNcoeffs());
    Array<OneD, NekDouble> tmpY(GetNcoeffs());
    
    Array<OneD, NekDouble> upwindX(GetNpoints());
    Array<OneD, NekDouble> upwindY(GetNpoints());
 
    IProductWRTDerivBase(0,GetPhys(field_no_1),tmpX,0);
    IProductWRTDerivBase(1,GetPhys(field_no_2),tmpY,0);
    
    Vmath::Vadd(GetNcoeffs(),tmpX,1,tmpY,1,tmpX,1);
    Vmath::Neg(GetNcoeffs(),tmpX,1);
    
    NumericalFluxDivergence(upwindX,upwindY,field_no_1,field_no_2);
    
    AddTraceIntegral(upwindX,upwindY,tmpX,0);
    
    MultiplyByElmtInvMass(tmpX,tmpX,0);

    BwdTrans(tmpX,out,1);
    
  }
  
  
  /**
   * Computes the \hat{f} term in the  \int_{\partial \Omega^e} \phi \hat{f} {\bf n} dS 
   * integral. Using averaged fluxes.
   **/
  void BoussinesqEquations::NumericalFluxGradient(Array<OneD, NekDouble> &outX, Array<OneD, NekDouble> &outY, int field_no)
  {
    
    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble> > Fwd(1);
    Array<OneD, Array<OneD, NekDouble> > Bwd(1);
    
    Fwd[0] = Array<OneD, NekDouble> (GetNpoints(),0.0); 
    Bwd[0] = Array<OneD, NekDouble> (GetNpoints(),0.0);
    
    // get the physical values at the trace
    GetFwdBwdTracePhys(Fwd[0],Bwd[0],field_no);
    
    for (int i = 0; i < GetNpoints(); ++i)
      {
 	outX[i]  = 0.5*(Fwd[0][i] + Bwd[0][i]);
	outY[i]  = 0.5*(Fwd[0][i] + Bwd[0][i]);
      }

  }

  
  /**
   * Computes the \hat{\bf f} term in the  \int_{\partial \Omega^e} \phi \hat{\bf f} \cdot {\bf n} dS 
   * integral. Using averaged fluxes.
   **/
  void BoussinesqEquations::NumericalFluxDivergence(Array<OneD, NekDouble> &outX, Array<OneD, NekDouble> &outY,
						    int field_no_1, int field_no_2)
  {

    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble> > Fwd(2);
    Array<OneD, Array<OneD, NekDouble> > Bwd(2);
    
    for (int i = 0; i < 2; ++i)
      {
	Fwd[i] = Array<OneD, NekDouble>(GetNpoints());
	Bwd[i] = Array<OneD, NekDouble>(GetNpoints());
      }
    
    // get the physical values at the trace
    GetFwdBwdTracePhys(Fwd[0],Bwd[0],field_no_1);
    GetFwdBwdTracePhys(Fwd[1],Bwd[1],field_no_2);
    
    for (int i = 0; i < GetNpoints(); ++i)
      {
	outX[i]  = 0.5*(Fwd[0][i] + Bwd[0][i]);
	outY[i]  = 0.5*(Fwd[1][i] + Bwd[1][i]);
      }
  }

} //end of namespace

