//////////////////////////////////////////////////////////////////////////////
//
// File BoussinesqSolver.cpp
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
// Description: Boussinesq equation solver 
// (comes in three flavours: Peregrine, enhanced and fully nonlinear)
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <cmath> 

#include <../solvers/Auxiliary/BoussinesqEquations.h>


using namespace Nektar;


void rhsFunction(BoussinesqEquationsSharedPtr U,
		 BoussinesqEquationsSharedPtr Uf,
		 Array<OneD, Array<OneD, NekDouble> > &rhs,
		 Array<OneD, NekDouble> f,
		 const NekDouble time);

int main(int argc, char *argv[])
{

  if(argc != 2)
    {
      fprintf(stderr,"Usage: BoussinesqSolver  meshfile \n");
      exit(1);
    }

   cout << "====================================================" <<endl;
   cout << "=============== BOUSSINESQ EQUATIONS ===============" <<endl;
   cout << "====================================================" <<endl;

   cout << "--------------------------------" <<endl;
   cout << "--------- DG version -----------" <<endl;
   cout << "--------------------------------" <<endl;

   {
     //-----------------------------------------
     // Read the mesh
     string fileNameString(argv[1]);
     SpatialDomains::MeshGraph2D mesh; 
     // Both the geometry and th expansion information should be read
     mesh.ReadGeometry(fileNameString);
     mesh.ReadExpansions(fileNameString);
     // Also read the boundary conditions
     SpatialDomains::BoundaryConditions boundaryConds(&mesh); 
     boundaryConds.Read(fileNameString);
     //-----------------------------------------


     //-----------------------------------------
     // Construct objects from the class BoussDisCont2D.
     int nVariables = 7;

     // Construct U: for the primary variables
     // U[0]: H  (total water depth)
     // U[1]: Hu (depth-averaged flux in the x-direction)
     // U[2]: Hv (depth-averaged flux in the y-direction)
     // U[3]: z  (=\nabla \cdot (H{\bf u)_t)
     // U[4]: d  (still water depth)
     // U[5]: ut (time derivative of u)
     // U[6]: vt (time derivative of v)
     
     BoussinesqEquationsSharedPtr U =
       MemoryManager<BoussinesqEquations>::AllocateSharedPtr(mesh,boundaryConds,nVariables);
     
     // Construct Uf: forcing storage for the variables in U 
     // (should be changed to copy constructor)
     BoussinesqEquationsSharedPtr Uf =
       MemoryManager<BoussinesqEquations>::AllocateSharedPtr(mesh,boundaryConds,nVariables);
     //-----------------------------------------
     
 
     //-----------------------------------------------
     // Store input parameters
     U->SetTime(boundaryConds.GetParameter("T0"));
     U->SetTimeStep(boundaryConds.GetParameter("Dt"));
     U->SetSteps(boundaryConds.GetParameter("Steps"));
     U->SetCheckSteps(boundaryConds.GetParameter("Check"));
     U->SetGravity(boundaryConds.GetParameter("G")); 
     U->SetConstantDepth(boundaryConds.GetParameter("D"));
     U->SetAlpha1(boundaryConds.GetParameter("Alpha1"));
     U->SetAlpha2(boundaryConds.GetParameter("Alpha2"));
     //-----------------------------------------------
     
    
     //-----------------------------------------
     // Read and evaluate the initial conditions
     // (stored in U)
     U->SetInitialConditions(boundaryConds, U->GetTime());

     // Write initial conditions to file
     stringstream outfileName_h;
     outfileName_h << "h_initial.dat";
     ofstream outfile((outfileName_h.str()).data());      
     U->WriteToFile(outfile,eTecplot,0);
     stringstream outfileName_hu;
     outfileName_hu << "hu_initial.dat";
     ofstream outfile1((outfileName_hu.str()).data());      
     U->WriteToFile(outfile1,eTecplot,1);
     stringstream outfileName_hv;
     outfileName_hv << "hv_initial.dat";
     ofstream outfile2((outfileName_hv.str()).data());      
     U->WriteToFile(outfile2,eTecplot,2);
     //-----------------------------------------

     
     //-----------------------------------------
     // Set Coriolis parameter
     U->SetCoriolis(boundaryConds);

     
     int nTotQuadPoints  = U->GetPointsTot();
     
     Array<OneD,NekDouble> x1(nTotQuadPoints,0.0);
     Array<OneD,NekDouble> x2(nTotQuadPoints,0.0);
     Array<OneD,NekDouble> x3(nTotQuadPoints,0.0);

     Uf->GetCoords(x1,x2,x3);
     
     Array<OneD, NekDouble> f(nTotQuadPoints);
     
     for(int i = 0; i < nTotQuadPoints; ++i)
	{
	  f[i]   = 0.0+1.0*x2[i];
	}      
     //-----------------------------------------------
     

     //-----------------------------------------------
     // MultiStep Arrays (should be "automatically" defined 
     // by the timeintegratorkey

     Array<OneD, Array<OneD, NekDouble> > U0(4);
     U0[0] = Array<OneD, NekDouble> (nTotQuadPoints);
     U0[1] = Array<OneD, NekDouble> (nTotQuadPoints);
     U0[2] = Array<OneD, NekDouble> (nTotQuadPoints);
     U0[3] = Array<OneD, NekDouble> (nTotQuadPoints);

     Array<OneD, Array<OneD, NekDouble> > U1(4);
     U1[0] = Array<OneD, NekDouble> (nTotQuadPoints);
     U1[1] = Array<OneD, NekDouble> (nTotQuadPoints);
     U1[2] = Array<OneD, NekDouble> (nTotQuadPoints);
     U1[3] = Array<OneD, NekDouble> (nTotQuadPoints);

     Array<OneD, Array<OneD, NekDouble> > U2(4);
     U2[0] = Array<OneD, NekDouble> (nTotQuadPoints);
     U2[1] = Array<OneD, NekDouble> (nTotQuadPoints);
     U2[2] = Array<OneD, NekDouble> (nTotQuadPoints);
     U2[3] = Array<OneD, NekDouble> (nTotQuadPoints);


     Array<OneD, Array<OneD,NekDouble> > f1(4);
     f1[0] = Array<OneD, NekDouble> (nTotQuadPoints,0.0);
     f1[1] = Array<OneD, NekDouble> (nTotQuadPoints,0.0);
     f1[2] = Array<OneD, NekDouble> (nTotQuadPoints,0.0);
     f1[3] = Array<OneD, NekDouble> (nTotQuadPoints,0.0);

     Array<OneD, Array<OneD,NekDouble> > f2(4);
     f2[0] = Array<OneD, NekDouble> (nTotQuadPoints,0.0);
     f2[1] = Array<OneD, NekDouble> (nTotQuadPoints,0.0);
     f2[2] = Array<OneD, NekDouble> (nTotQuadPoints,0.0);
     f2[3] = Array<OneD, NekDouble> (nTotQuadPoints,0.0);

     Array<OneD, Array<OneD,NekDouble> > f3(4);
     f3[0] = Array<OneD, NekDouble> (nTotQuadPoints,0.0);
     f3[1] = Array<OneD, NekDouble> (nTotQuadPoints,0.0);
     f3[2] = Array<OneD, NekDouble> (nTotQuadPoints,0.0);
     f3[3] = Array<OneD, NekDouble> (nTotQuadPoints,0.0);

     //-----------------------------------------------


     //-----------------------------------------------
     // start the time loop
     int chk = 0;
     for (int i = 1; i < U->GetSteps() + 1; ++i)
	{
	  

	  //-----------------------------------------------
	  // Third-order TVD Runge-Kutta scheme
	  
	  // Store the starting value in U0
	  //U->GetPhys(U0);	  
	  U0[0] = U->GetPhys(0);
	  U0[1] = U->GetPhys(1);
	  U0[2] = U->GetPhys(2);
	  U0[3] = U->GetPhys(3);

	  //RK step one
	  rhsFunction(U,Uf,f1,f,U->GetTime());
	  Vmath::Svtvp(nTotQuadPoints,U->GetTimeStep(),f1[0],1,U0[0],1,U1[0],1);
	  Vmath::Svtvp(nTotQuadPoints,U->GetTimeStep(),f1[1],1,U0[1],1,U1[1],1);
	  Vmath::Svtvp(nTotQuadPoints,U->GetTimeStep(),f1[2],1,U0[2],1,U1[2],1);
	  U->SetPhys(U1,0);
	  U->SetPhys(U1,1);
	  U->SetPhys(U1,2);


	  //RK step two
	  rhsFunction(U,Uf,f2,f,U->GetTime()+U->GetTimeStep());
	 //  for(int j = 0; j < nTotQuadPoints; ++j)
// 	    {
// 	      U1[0][j] = U0[0][j] + U->GetTimeStep() * (0.5*f1[0][j]+0.5*f2[0][j]);
// 	      U1[1][j] = U0[1][j] + U->GetTimeStep() * (0.5*f1[1][j]+0.5*f2[1][j]);
// 	      U1[2][j] = U0[2][j] + U->GetTimeStep() * (0.5*f1[2][j]+0.5*f2[2][j]);

// 	    }

	  Vmath::Svtvp(nTotQuadPoints,U->GetTimeStep(),f2[0],1,U1[0],1,U2[0],1);
	  Vmath::Svtvp(nTotQuadPoints,3.0,U0[0],1,U2[0],1,U2[0],1);
	  Vmath::Smul(nTotQuadPoints,0.25,U2[0],1,U2[0],1);

	  Vmath::Svtvp(nTotQuadPoints,U->GetTimeStep(),f2[1],1,U1[1],1,U2[1],1);
	  Vmath::Svtvp(nTotQuadPoints,3.0,U0[1],1,U2[1],1,U2[1],1);
	  Vmath::Smul(nTotQuadPoints,0.25,U2[1],1,U2[1],1);
	  
	  Vmath::Svtvp(nTotQuadPoints,U->GetTimeStep(),f2[2],1,U1[2],1,U2[2],1);
	  Vmath::Svtvp(nTotQuadPoints,3.0,U0[2],1,U2[2],1,U2[2],1);
	  Vmath::Smul(nTotQuadPoints,0.25,U2[2],1,U2[2],1);
 
	  U->SetPhys(U2,0);
	  U->SetPhys(U2,1);
	  U->SetPhys(U2,2);

	  //RK step three
	  rhsFunction(U,Uf,f3,f,U->GetTime()+U->GetTimeStep());

	  Vmath::Svtvp(nTotQuadPoints,U->GetTimeStep(),f3[0],1,U2[0],1,U1[0],1);
	  Vmath::Svtvp(nTotQuadPoints,2.0,U1[0],1,U0[0],1,U1[0],1);
	  Vmath::Smul(nTotQuadPoints,(1.0/3.0),U1[0],1,U1[0],1);

	  Vmath::Svtvp(nTotQuadPoints,U->GetTimeStep(),f3[1],1,U2[1],1,U1[1],1);
	  Vmath::Svtvp(nTotQuadPoints,2.0,U1[1],1,U0[1],1,U1[1],1);
	  Vmath::Smul(nTotQuadPoints,(1.0/3.0),U1[1],1,U1[1],1);
	  
	  Vmath::Svtvp(nTotQuadPoints,U->GetTimeStep(),f3[2],1,U2[2],1,U1[2],1);
	  Vmath::Svtvp(nTotQuadPoints,2.0,U1[2],1,U0[2],1,U1[2],1);
	  Vmath::Smul(nTotQuadPoints,(1.0/3.0),U1[2],1,U1[2],1);

	  U->SetPhys(U1,0);
	  U->SetPhys(U1,1);
	  U->SetPhys(U1,2);
	  
	  U->SetTime(U->GetTime()+U->GetTimeStep());
	  //-----------------------------------------------

	  
	  //-----------------------------------------------
	  // Write chk files
	  if (i%U->GetCheckSteps() == 0)
	    {  
	      cout << "Time: " << U->GetTime() << endl;
	      stringstream outfileName_h;
	      outfileName_h << "h_Quad_"<< chk <<".dat";
	      ofstream outfile((outfileName_h.str()).data());      
	      U->WriteToFile(outfile,eTecplot,0);
	     
	      stringstream outfileName_hu;
	      outfileName_hu << "hu_Quad_"<< chk <<".dat";
	      ofstream outfile1((outfileName_hu.str()).data());      
	      U->WriteToFile(outfile1,eTecplot,1);
	      
	      stringstream outfileName_hv;
	      outfileName_hv << "hv_Quad_"<< chk <<".dat";
	      ofstream outfile2((outfileName_hv.str()).data());      
	      U->WriteToFile(outfile2,eTecplot,2);

	      stringstream outfileName_z;
	      outfileName_z << "z_Quad_"<< chk <<".dat";
	      ofstream outfile3((outfileName_z.str()).data());      
	      U->WriteToFile(outfile3,eTecplot,5);
	      ++chk;
	    }
	  //-----------------------------------------------
	  
	  
	  //-----------------------------------------------
	  // Compute the L2 error
	  
	  if (i%U->GetSteps() == 0)
	    {
                Array<OneD, NekDouble> exactSolution(nTotQuadPoints);
                SpatialDomains::ConstForcingFunctionShPtr exactSolutionEquation 
                    = boundaryConds.GetExactSolution(boundaryConds.GetVariable(0));
                for(int j = 0; j < nTotQuadPoints; ++j)
		{
		  exactSolution[j] = exactSolutionEquation->Evaluate(x1[j],x2[j],x3[j],0.0);//U->GetTime());
		}      
	      
                MultiRegions::ExpList2DSharedPtr exactSolutionExp =
                    MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(mesh);
                exactSolutionExp->SetPhys(exactSolution);
	      
                NekDouble error = U->L2(*exactSolutionExp,0);
	      
                // Display the output              
                cout << "eta: Ndof = " << U->GetNcoeffs() << " => L2 error = " << error << endl;
	    }
            //-----------------------------------------------
	}

   }
   
}


void rhsFunction(BoussinesqEquationsSharedPtr U,
		 BoussinesqEquationsSharedPtr Uf,
		 Array<OneD, Array<OneD, NekDouble> > &rhs,
		 Array<OneD, NekDouble> f,
		 const NekDouble time)
{

    int nTotQuadPoints  = U->GetPointsTot();
    int nTotCoeffs      = U->GetNcoeffs();
    int nTotTracePoints = U->GetNpoints();
    int nVariables      = U->GetNvariables();
    
    //---------------------------------------
    // Define local arrays  
    
    Array<OneD, Array<OneD, NekDouble> > upwindX(nVariables);
    Array<OneD, Array<OneD, NekDouble> > upwindY(nVariables);
    Array<OneD, Array<OneD, NekDouble> > iprod_rhs(nVariables);
    Array<OneD, Array<OneD, NekDouble> > fluxVectorX(nVariables);
    Array<OneD, Array<OneD, NekDouble> > fluxVectorY(nVariables);
    
    for (int i = 0; i < 4; ++i)
      {
	upwindX[i] = Array<OneD, NekDouble>(nTotTracePoints);
	upwindY[i] = Array<OneD, NekDouble>(nTotTracePoints);
	iprod_rhs[i] = Array<OneD, NekDouble>(nTotCoeffs);
	fluxVectorX[i] = Array<OneD, NekDouble>(nTotQuadPoints);
	fluxVectorY[i] = Array<OneD, NekDouble>(nTotQuadPoints);
      }
    //---------------------------------------
    
    
    //---------------------------------------
    // Now we start by computing the advection part
    // corresponding to the Shallow Water Equations
    
    // Update user defined boundary conditions					
    //U->SetBoundaryConditions();
    
    // Compute the flux vector {\bf F} (physical space)
    U->GetFluxVector(fluxVectorX,fluxVectorY);
    
    // Compute the innerproduct (Grad \phi \cdot {\bf F}) (modal space)
    Array<OneD, NekDouble> tmp(nTotCoeffs,0.0);
    
    for (int i = 0; i < 3; ++i)
      {
	U->IProductWRTDerivBase(0,fluxVectorX[i],iprod_rhs[i],0);
	U->IProductWRTDerivBase(1,fluxVectorY[i],tmp,0);
	Vmath::Vadd(nTotCoeffs,tmp,1,iprod_rhs[i],1,iprod_rhs[i],1);
      }
    
    // Evaluate  upwind numerical flux (physical space)
    U->NumericalFlux(upwindX,upwindY);
    
    // Compute the boundary integral (modal space)
    for (int i = 0; i < 3; ++i)
      {
	Vmath::Neg(nTotTracePoints,upwindX[i],1);
	Vmath::Neg(nTotTracePoints,upwindY[i],1);
	
	U->AddTraceIntegral(upwindX[i],upwindY[i],iprod_rhs[i],i);
      }
    
    // ok: SWE advection... done!
    //---------------------------------------
    

    //---------------------------------------
    // As no more terms is required for the
    // continuity equations we evaluate the 
    // physical values for H_t here.
    // we use it later in the spatial
    // dispersive terms
    
    // Solve the block-diagonal system
    Uf->MultiplyByElmtInvMass(iprod_rhs[0],Uf->UpdateCoeffs(0),0);
    
    // Go from modal to physical space
    Uf->BwdTrans(0);
    
    // copy into the rhs arrays
    Vmath::Vcopy(nTotQuadPoints,Uf->GetPhys(0),1,rhs[0],1);
    
    // ok: H_t in rhs[0]... done!
    //---------------------------------------

        
  //   //---------------------------------------
//     // Compute source terms (in this case the Coriolis force)
//     Array<OneD, NekDouble> tmp1(nTotQuadPoints,0.0);

//     Vmath::Vmul(nTotQuadPoints,U->GetCoriolis(),1,U->GetPhys(2),1,tmp1,1);
//     U->IProductWRTBase(tmp1,tmp,0);
//     Vmath::Vadd(nTotCoeffs,tmp,1,iprod_rhs[1],1,iprod_rhs[1],1);

//     Vmath::Vmul(nTotQuadPoints,U->GetCoriolis(),1,U->GetPhys(1),1,tmp1,1);
//     Vmath::Neg(nTotQuadPoints,tmp1,1);
//     U->IProductWRTBase(tmp1,tmp,0);
//     Vmath::Vadd(nTotCoeffs,tmp,1,iprod_rhs[2],1,iprod_rhs[2],1);

    // Add source terms
    //U->SourceTerms(iprod_rhs);

//     //---------------------------------------


    //---------------------------------------
    // Add the spatial dispersive terms
    
    //U->Madsen92SpatialTerms(iprod_rhs[1],iprod_rhs[2]);
    U->Lambda20Spatial(*Uf,iprod_rhs[1],iprod_rhs[2],rhs[0]);
    //--------------------------

    
    //--------------------------
    // store f1 and f2 for later use
    
    Array<OneD, NekDouble> f1(nTotCoeffs);
    Array<OneD, NekDouble> f2(nTotCoeffs);
    Vmath::Vcopy(nTotCoeffs,iprod_rhs[1],1,f1,1); // f1
    Vmath::Vcopy(nTotCoeffs,iprod_rhs[2],1,f2,1); // f2
    //--------------------------



    //---------------------------------------
    // Solve the block-diagonal system
    
    U->MultiplyByElmtInvMass(iprod_rhs[1],Uf->UpdateCoeffs(1),1);
    U->MultiplyByElmtInvMass(iprod_rhs[2],Uf->UpdateCoeffs(2),2);

    //---------------------------------------

    
    //---------------------------------------
    // Go from modal to physical space
    
    Uf->BwdTrans(1);
    Uf->BwdTrans(2);
    //---------------------------------------
    

    //---------------------------------------
    // copy into the rhs arrays
    
    Vmath::Vcopy(nTotQuadPoints,Uf->GetPhys(1),1,rhs[1],1);
    Vmath::Vcopy(nTotQuadPoints,Uf->GetPhys(2),1,rhs[2],1);
    //---------------------------------------

 

    //---------------------------------------
    // Start for solve of mixed dispersive terms
    // using the 'scalar method' 
    // (Eskilsson & Sherwin, JCP 2006)
    
    // Set boundary condidtions for z
    //U->SetBoundaryConditionsWaveCont(); // wall like bc
    
    // Compute the forcing function for the
    // wave continuity equation
    
    U->IProductWRTDerivBase(0,rhs[1],iprod_rhs[3],0);
    U->IProductWRTDerivBase(1,rhs[2],tmp,0);
    Vmath::Vadd(nTotCoeffs,tmp,1,iprod_rhs[3],1,iprod_rhs[3],1);
    
    // Evaluate  upwind numerical flux (physical space)
    Uf->NumericalFluxWaveCont(upwindX[3],upwindY[3]);
    
    Vmath::Neg(nTotTracePoints,upwindX[3],1);
    Vmath::Neg(nTotTracePoints,upwindY[3],1);
	
    U->AddTraceIntegral(upwindX[3],upwindY[3],iprod_rhs[3],0);
    
    U->MultiplyByElmtInvMass(iprod_rhs[3],U->UpdateCoeffs(3),0);
    U->BwdTrans(3); 

    Vmath::Vcopy(nTotQuadPoints,U->GetPhys(3),1,rhs[3],1);
    
    // ok: forcing function for HelmSolve... done!
    //--------------------------------------
    
    
    //--------------------------------------
    // Solve the Helmhotz-type equation
    // for the wave continuity equation
    
    NekDouble gamma = (U->GetConstantDepth() * U->GetConstantDepth())*(U->GetAlpha1()+1.0/3.0);
    NekDouble invgamma = 1.0/gamma;
    
//      // U->SetBoundaryConditionsSolve(); // equal zero

    Vmath::Smul(nTotQuadPoints,invgamma,rhs[3],1,rhs[3],1);
    U->WaveContSolve(rhs[3],invgamma);
    
    Vmath::Vcopy(nTotQuadPoints,U->GetPhys(3),1,rhs[3],1);
    
    // ok: Wave Continuity Equation... done! 
    //------------------------------------
    

    //------------------------------------
    // Return to the primary variables
    
    // Set boundary conditions 
    //U->SetBoundaryConditionsContVariables(); // -interior = exteroir

    U->IProductWRTDerivBase(0,rhs[3],iprod_rhs[1],0);
    U->IProductWRTDerivBase(1,rhs[3],iprod_rhs[2],0);

    Vmath::Neg(nTotCoeffs,iprod_rhs[1],1);
    Vmath::Neg(nTotCoeffs,iprod_rhs[2],1);
    
    // Evaluate  upwind numerical flux (physical space)
    U->NumericalFluxConsVariables(upwindX[3],upwindY[3]);
    
    {
      Array<OneD, NekDouble> uptemp(nTotTracePoints,0.0);
      
      U->AddTraceIntegral(upwindX[3],uptemp,iprod_rhs[1],0);
      U->AddTraceIntegral(uptemp,upwindY[3],iprod_rhs[2],0);
    }
    
    Vmath::Smul(nTotCoeffs,gamma,iprod_rhs[1],1,iprod_rhs[1],1);
    Vmath::Smul(nTotCoeffs,gamma,iprod_rhs[2],1,iprod_rhs[2],1);

    Vmath::Vadd(nTotCoeffs,f1,1,iprod_rhs[1],1,iprod_rhs[1],1);
    Vmath::Vadd(nTotCoeffs,f2,1,iprod_rhs[2],1,iprod_rhs[2],1);
    
    Uf->MultiplyByElmtInvMass(iprod_rhs[1],Uf->UpdateCoeffs(1),0);
    Uf->MultiplyByElmtInvMass(iprod_rhs[2],Uf->UpdateCoeffs(2),0);
    
    Uf->BwdTrans(1);
    Uf->BwdTrans(2);

    Vmath::Vcopy(nTotQuadPoints,Uf->GetPhys(1),1,rhs[1],1);
    Vmath::Vcopy(nTotQuadPoints,Uf->GetPhys(2),1,rhs[2],1);

    // ok: returned to conservative variables... done!
    //---------------------------------
  

    //---------------------------------
    // Store u_t and v_t
    
    // Note: these values should be time averaged but	
    // as for now we just keep them as is
    {
      Array<OneD, NekDouble> u(nTotQuadPoints);
      Array<OneD, NekDouble> v(nTotQuadPoints);
      
      
      for (int j = 0; j < nTotQuadPoints; ++j)
	{
	  u[j] = U->GetPhys(1)[j]/U->GetPhys(0)[j];
	  v[j] = U->GetPhys(2)[j]/U->GetPhys(0)[j];
	}
      
      // u_t =  ( (Hu)_t - H_t u ) / H
      // v_t =  ( (Hv)_t - H_t v ) / H
      // (Hu)_t is stored in rhs[1]
      // (Hu)_t is stored in rhs[2]
      //  H_t is stored in rhs[0]
      
      for (int j = 0; j < nTotQuadPoints; ++j)
 	{
 	  U->UpdatePhys(5)[j] = (rhs[1][j] - rhs[0][j]*u[j]) / (U->GetPhys(0)[j]);
 	  U->UpdatePhys(6)[j] = (rhs[2][j] - rhs[0][j]*v[j]) / (U->GetPhys(0)[j]);
 	}
      

    }
    
    // ok: u_t and v_t... done!
    //---------------------------------
}




