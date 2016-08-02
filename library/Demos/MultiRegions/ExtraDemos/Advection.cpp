#include <cstdio>
#include <cstdlib>
#include <cmath> 

#include <MultiRegions/DisContField2D.h>

using namespace std;
using namespace Nektar;


void rhsFunction(MultiRegions::DisContField2DSharedPtr hpExp, 
                 Array<OneD, const NekDouble> v1, 
		 Array<OneD, const NekDouble> v2,
                 Array<OneD, const NekDouble> v1Trace, 
		 Array<OneD, const NekDouble> v2Trace, 
                 Array<OneD, const NekDouble> u0,
		 Array<OneD, NekDouble> &rhs, 
                 const NekDouble time);

int main(int argc, char *argv[])
{

 if(argc != 2)
   {
     fprintf(stderr,"Usage: AdvectionDis  meshfile \n");
     exit(1);
   }

  cout << "=====================================================" <<endl;
  cout << "============= LINEAR ADVECTION EQUATION==============" <<endl;
  cout << "=====================================================" <<endl;

  cout << "-------------------------------" <<endl;
  cout << "--- Rotating Gaussian hump ----" <<endl;
  cout << "--------- DG version ----------" <<endl;
  cout << "-------------------------------" <<endl;

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
    // Construct an object from the class DisContField2D.
    MultiRegions::DisContField2DSharedPtr u =
        MemoryManager<MultiRegions::DisContField2D>::AllocateSharedPtr(mesh,boundaryConds,"u");
    //-----------------------------------------

    //-----------------------------------------
    // Read and evaluate the initial conditions
    int nTotQuadPoints  = u->GetTotPoints();

    Array<OneD,NekDouble> x1(nTotQuadPoints,0.0);
    Array<OneD,NekDouble> x2(nTotQuadPoints,0.0);

    u->GetCoords(x1,x2);

    Array<OneD, NekDouble> u0(nTotQuadPoints);
    Array<OneD, NekDouble> u1(nTotQuadPoints);
    Array<OneD, NekDouble> v1(nTotQuadPoints,0.0);
    Array<OneD, NekDouble> v2(nTotQuadPoints,0.0);

    SpatialDomains::ConstForcingFunctionShPtr u0SolutionEquation 
        = boundaryConds.GetExactSolution(boundaryConds.GetVariable(0));

    u0SolutionEquation->Evaluate(x1,x2, u0);
    for(int i = 0; i < nTotQuadPoints; ++i)
    {
      v1[i] =  1.0;//2.0*M_PI*x2[i];
      v2[i] =  0.0;//-2.0*M_PI*x1[i];
    }

    u->SetPhys(u0);
    //-----------------------------------------------

    //-----------------------------------------------
    // Initiate the velocities on the trace
    int nTotTracePoints = u->GetTrace()->GetNpoints();

    Array<OneD, NekDouble> v1Trace(nTotTracePoints,0.0);
    Array<OneD, NekDouble> v2Trace(nTotTracePoints,0.0);

    u->SetPhys(v1);
    u->ExtractTracePhys(v1Trace);

    u->SetPhys(v2);
    u->ExtractTracePhys(v2Trace);

    u->SetPhys(u0);
    //-----------------------------------------------

    //-----------------------------------------------
    // Write initial conditions to file
    stringstream outfileName;
    outfileName << "Advection_Quad_P9_initial.dat";
    ofstream outfile((outfileName.str()).data());      
    u->WriteToFile(outfile,eTecplot);
    //-----------------------------------------------

    //-----------------------------------------------
    // RK vectors
    Array<OneD, NekDouble> f1(nTotQuadPoints,0.0);
    Array<OneD, NekDouble> f2(nTotQuadPoints,0.0);
    //-----------------------------------------------

    //-----------------------------------------------
    // Time stepping parameters
    NekDouble dt = boundaryConds.GetParameter("Dt");
    int timeStep = boundaryConds.GetParameter("Steps");
    int chkStep  = boundaryConds.GetParameter("Check");
    //-----------------------------------------------

    //-----------------------------------------------
    // start the time loop
    int chk = 0;
    for (int i = 1; i < timeStep + 1; ++i)
	{
	  
	  // Store the starting value
	  u0 = u->GetPhys();
	  
	  //RK step one
	  rhsFunction(u,v1,v2,v1Trace,v2Trace,u0,f1,i*dt);
	  Vmath::Svtvp(nTotQuadPoints,dt,f1,1,u0,1,u1,1);
	  u->SetPhys(u1);
	
	  
	  //RK step two
	  rhsFunction(u,v1,v2,v1Trace,v2Trace,u1,f2,i*dt+dt);
	  for(int j = 0; j < nTotQuadPoints; ++j)
	    {
	      u1[j] = u0[j] + dt *(0.5*f1[j]+0.5*f2[j]);
	    }
	  u->SetPhys(u1);


	  //-----------------------------------------------
	  // Write chk files
	  if (i%chkStep == 0)
	    {  
	      cout << "Time: " << i*dt << endl;
	      stringstream outfileName;
	      outfileName << "Advection_Quad_P9_"<< chk <<".dat";
	      ofstream outfile((outfileName.str()).data());      
	      u->WriteToFile(outfile,eTecplot);
	      ++chk;
	    }
	  //-----------------------------------------------
	  	  
	  //-----------------------------------------------
	  // Compute the L2 error
	  if (i == timeStep)
          {
	      Array<OneD, NekDouble> exactSolution(nTotQuadPoints);
	      SpatialDomains::ConstForcingFunctionShPtr exactSolutionEquation 
            = boundaryConds.GetExactSolution(boundaryConds.GetVariable(0));

        exactSolutionEquation->Evaluate(x1,x2, exactSolution);
        MultiRegions::DisContField2DSharedPtr exactSolutionExp =
            MemoryManager<MultiRegions::DisContField2D>::AllocateSharedPtr(*u);
        exactSolutionExp->SetPhys(exactSolution);
	      
	      NekDouble error = u->L2(exactSolutionExp->GetPhys());
	      
	      // Display the output              
	      cout << "Ndof = " << u->GetNcoeffs() << " => Error = " << error << endl;
	    }
	  //-----------------------------------------------

	}
  }
}


void rhsFunction(MultiRegions::DisContField2DSharedPtr hpExp,
                Array<OneD, const NekDouble> v1, 
                Array<OneD, const NekDouble> v2, 
                Array<OneD, const NekDouble> v1Trace, 
		 Array<OneD, const NekDouble> v2Trace, 
                Array<OneD, const NekDouble> u0, 
		 Array<OneD, NekDouble> &rhs, const NekDouble time)
{

   int nTotQuadPoints  = hpExp->GetTotPoints();
   int nTotCoeffs      = hpExp->GetNcoeffs();
   int nTotTracePoints = hpExp->GetTrace()->GetNpoints();

   //---------------------------------------
   // get the trace values

   Array<OneD,NekDouble> Fwd(nTotTracePoints,0.0);
   Array<OneD,NekDouble> Bwd(nTotTracePoints,0.0);
   Array<OneD,NekDouble> upwind1(nTotTracePoints,0.0);
   Array<OneD,NekDouble> upwind2(nTotTracePoints,0.0);

   Array<OneD, Array<OneD, NekDouble> > normals;

   hpExp->GetFwdBwdTracePhys(Fwd,Bwd); 
   //---------------------------------------


   //---------------------------------------
   // Compute the flux vector

   Array<OneD, NekDouble> fluxVector0(nTotQuadPoints,0.0);
   Array<OneD, NekDouble> fluxVector1(nTotQuadPoints,0.0);

   // Fill the flux vector using the old time step
   for (int i = 0; i < nTotQuadPoints; ++i)
   {
       fluxVector0[i] = v1[i]*u0[i];
       fluxVector1[i] = v2[i]*u0[i];
   }
   //---------------------------------------


   Array<OneD, NekDouble> iprod_rhs(nTotCoeffs,0.0);
   Array<OneD, NekDouble> tmp(nTotCoeffs,0.0);
   Array<OneD, NekDouble> iprod_rhs_test(nTotCoeffs,0.0);
   //--------------------------------------


   //---------------------------------------
   // Compute the (Grad \phi \cdot F)

   hpExp->IProductWRTDerivBase(0,fluxVector0,iprod_rhs);
   hpExp->IProductWRTDerivBase(1,fluxVector1,tmp);
   Vmath::Vadd(nTotCoeffs,tmp,1,iprod_rhs,1,iprod_rhs,1);
   //---------------------------------------


   //---------------------------------------
   // Evaluate  upwind numerical flux
   // note: upwind1 contains  a_1 * upwind
   //       upwind2 contains  a_2 * upwind

#if 1
   Array<OneD, Array<OneD, NekDouble> > Vel(2);
   Vel[0] = v1Trace;
   Vel[1] = v2Trace;

   hpExp->GetTrace()->Upwind(Vel,Fwd,Bwd,upwind1);

   Vmath::Vmul(nTotTracePoints,upwind1,1,v2Trace,1,upwind2,1);
   Vmath::Vmul(nTotTracePoints,upwind1,1,v1Trace,1,upwind1,1);
#else
   for (int i = 0; i < nTotTracePoints; ++i)
     {
       if (v1Trace[i]*normals[0][i]+v2Trace[i]*normals[1][i] > 0.0)
	  {
	    upwind1[i] = v1Trace[i]*Fwd[i];
	    upwind2[i] = v2Trace[i]*Fwd[i];
	  }
       else
	  {
	    upwind1[i] = v1Trace[i]*Bwd[i];
	    upwind2[i] = v2Trace[i]*Bwd[i];
	  }
     }
#endif

   Vmath::Neg(nTotTracePoints,upwind1,1);
   Vmath::Neg(nTotTracePoints,upwind2,1);

   // this should then be changed 
   hpExp->AddTraceIntegral(upwind1, upwind2, iprod_rhs);

   hpExp->MultiplyByElmtInvMass(iprod_rhs,hpExp->UpdateCoeffs());


   //---------------------------------------
   // Go to physical space
   hpExp->BwdTrans(hpExp->GetCoeffs(),rhs);
   //---------------------------------------
}
