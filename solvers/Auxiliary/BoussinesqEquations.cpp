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
	id2  = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndExpToTraceExpMap(cnt+e));
	
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
	id2  = m_fields[3]->GetTrace()->GetPhys_Offset(m_fields[3]->GetTraceMap()->GetBndExpToTraceExpMap(cnt+e));
	
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
	id2  = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndExpToTraceExpMap(cnt+e));
	
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

    
    // compute DivU = \nabla \cdot {\bf u}
    {
    Array<OneD, NekDouble> DivU(nTotQuadPoints);

    
    
    

    }

    // compute GradDivU = \nabla DivU
    Array<OneD, Array<OneD, NekDouble> > GradDivU(2);
    GradDivU[0] = Array<OneD, NekDouble>(nTotQuadPoints);
    GradDivU[1] = Array<OneD, NekDouble>(nTotQuadPoints);


    


    // compute DivdU = \nabla \cdot (d {\bf u})
    Array<OneD, NekDouble> DivdU(nTotQuadPoints);

    

    // compute GradDivdU = \nabla DivdU
    Array<OneD, Array<OneD, NekDouble> > GradDivdU(2);
    GradDivdU[0] = Array<OneD, NekDouble>(nTotQuadPoints);
    GradDivdU[1] = Array<OneD, NekDouble>(nTotQuadPoints);


    // compute DivUt = \nabla \cdot {\bf ut}



    // compute GradDivUt = \nabla DivUt

    


    // compute DivdUt = \nabla \cdot (d {\bf ut})

    

    // compute GradDivdUt = \nabla DivdUt


    // compute Gamma = (d/6)GradDivU - (1/2)GradDivdU


    
    
    // compute Gammat = (d/6)GradDivUt - (1/2)GradDivdUt

    

    // Add Lambda_20 terms



    // Add Lambda_21 terms



    // Add Lambda_22 terms



    // Add Lambda_23 terms



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

