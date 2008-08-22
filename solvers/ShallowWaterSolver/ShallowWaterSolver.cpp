///////////////////////////////////////////////////////////////////////////////
//
// File ShallowWaterSolver.cpp
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
// Description: Shallow water equations solver
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <cmath> 

#include <../solvers/Auxiliary/ShallowWaterEquations.h>


using namespace Nektar;


void rhsFunction(ShallowWaterEquationsSharedPtr U,
		 Array<OneD, Array<OneD, NekDouble> > &rhs,
		 Array<OneD, NekDouble> f,
		 const NekDouble time);

int main(int argc, char *argv[])
{

    if(argc != 2)
    {
        fprintf(stderr,"Usage: ShallowWaterSolver meshfile \n");
        exit(1);
    }

    cout << "====================================================" <<endl;
    cout << "============= SHALLOW WATER EQUATIONS ==============" <<endl;
    cout << "====================================================" <<endl;

    cout << "--------------------------------" <<endl;
    cout << "--- Equatorial Rossby Modon ----" <<endl;
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
        // Construct an object from the class ShallowWaterEquations.
        ShallowWaterEquationsSharedPtr U =
            MemoryManager<ShallowWaterEquations>::AllocateSharedPtr(mesh,boundaryConds);
        //-----------------------------------------
     
     
        //-----------------------------------------------
        // Store input parameters
        U->SetTime(boundaryConds.GetParameter("T0"));
        U->SetTimeStep(boundaryConds.GetParameter("Dt"));
        U->SetSteps(boundaryConds.GetParameter("Steps"));
        U->SetCheckSteps(boundaryConds.GetParameter("Check"));
        U->SetGravity(boundaryConds.GetParameter("G"));
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
     
        U->GetCoords(x1,x2);
     
        Array<OneD, NekDouble> f(nTotQuadPoints);
     
        for(int i = 0; i < nTotQuadPoints; ++i)
	{
            f[i]   = 0.0+1.0*x2[i];
	}      
        //-----------------------------------------------
     

        //-----------------------------------------------
        // MultiStep Arrays (should be "automatically" defined 
        // by the timeintegratorkey

        Array<OneD, Array<OneD, NekDouble> > U0(3);
        U0[0] = Array<OneD, NekDouble> (nTotQuadPoints);
        U0[1] = Array<OneD, NekDouble> (nTotQuadPoints);
        U0[2] = Array<OneD, NekDouble> (nTotQuadPoints);

        Array<OneD, Array<OneD, NekDouble> > U1(3);
        U1[0] = Array<OneD, NekDouble> (nTotQuadPoints);
        U1[1] = Array<OneD, NekDouble> (nTotQuadPoints);
        U1[2] = Array<OneD, NekDouble> (nTotQuadPoints);
        
        Array<OneD, Array<OneD,NekDouble> > f1(3);
        f1[0] = Array<OneD, NekDouble> (nTotQuadPoints,0.0);
        f1[1] = Array<OneD, NekDouble> (nTotQuadPoints,0.0);
        f1[2] = Array<OneD, NekDouble> (nTotQuadPoints,0.0);
     
        Array<OneD, Array<OneD,NekDouble> > f2(3);
        f2[0] = Array<OneD, NekDouble> (nTotQuadPoints,0.0);
        f2[1] = Array<OneD, NekDouble> (nTotQuadPoints,0.0);
        f2[2] = Array<OneD, NekDouble> (nTotQuadPoints,0.0);
        //-----------------------------------------------


        //-----------------------------------------------
        // Time stepping parameters
        // NekDouble dt = boundaryConds.GetParameter("Dt");
        // int timeStep = boundaryConds.GetParameter("Steps");
        // int chkStep  = boundaryConds.GetParameter("Check");
        //-----------------------------------------------

        //-----------------------------------------------
        // start the time loop
        int chk = 0;
        for (int i = 1; i < U->GetSteps() + 1; ++i)
	{
	  

            //-----------------------------------------------
            // Runge-Kutta scheme
	  
            // Store the starting value in U0
            U->GetPhys(U0);
	  
            //RK step one
            rhsFunction(U,f1,f,U->GetTime());
            Vmath::Svtvp(nTotQuadPoints,U->GetTimeStep(),f1[0],1,U0[0],1,U1[0],1);
            Vmath::Svtvp(nTotQuadPoints,U->GetTimeStep(),f1[1],1,U0[1],1,U1[1],1);
            Vmath::Svtvp(nTotQuadPoints,U->GetTimeStep(),f1[2],1,U0[2],1,U1[2],1);
            U->SetPhys(U1);
	  
            //RK step two
            rhsFunction(U,f2,f,U->GetTime()+U->GetTimeStep());
            for(int j = 0; j < nTotQuadPoints; ++j)
	    {
                U1[0][j] = U0[0][j] + U->GetTimeStep() * (0.5*f1[0][j]+0.5*f2[0][j]);
                U1[1][j] = U0[1][j] + U->GetTimeStep() * (0.5*f1[1][j]+0.5*f2[1][j]);
                U1[2][j] = U0[2][j] + U->GetTimeStep() * (0.5*f1[2][j]+0.5*f2[2][j]);

	    }
            U->SetPhys(U1);
	  
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
                ++chk;
	    }
            //-----------------------------------------------
	  
	  
            // 	  //-----------------------------------------------
            // 	  // Compute the L2 error
            // 	  if (i == timeStep)
            // 	    {
            // 	      Array<OneD, NekDouble> exactSolution(nTotQuadPoints);
            // 	      SpatialDomains::ConstForcingFunctionShPtr exactSolutionEquation 
            // 		= boundaryConds.GetExactSolution(boundaryConds.GetVariable(0));
            // 	      for(int j = 0; j < nTotQuadPoints; ++j)
            // 		{
            // 		  exactSolution[j] = exactSolutionEquation->Evaluate(x1[j],x2[j]);
            // 		}      
	      
            // 	      MultiRegions::DisContField2DSharedPtr exactSolutionExp =
            // 		MemoryManager<MultiRegions::DisContField2D>::AllocateSharedPtr(*u);
            // 	      exactSolutionExp->SetPhys(exactSolution);
	      
            // 	      NekDouble error = u->L2(*exactSolutionExp);
	      
            // 	      // Display the output              
            // 	      cout << "Ndof = " << u->GetNcoeffs() << " => Error = " << error << endl;
            // 	    }
            // 	  //-----------------------------------------------
	}

    }
   
}


void rhsFunction(ShallowWaterEquationsSharedPtr U,
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
    
    for (int i = 0; i < nVariables; ++i)
    {
	upwindX[i] = Array<OneD, NekDouble>(nTotTracePoints);
	upwindY[i] = Array<OneD, NekDouble>(nTotTracePoints);
	iprod_rhs[i] = Array<OneD, NekDouble>(nTotCoeffs);
	fluxVectorX[i] = Array<OneD, NekDouble>(nTotQuadPoints);
	fluxVectorY[i] = Array<OneD, NekDouble>(nTotQuadPoints);
    }
     
    //---------------------------------------
    // Update user defined boundary conditions					
    U->SetBoundaryConditions();
    //---------------------------------------
    
    
    //---------------------------------------
    // Compute the advection part
      
    // Compute the flux vector {\bf F} (physical space)
    U->GetFluxVector(fluxVectorX,fluxVectorY);
    
    // Compute the innerproduct (Grad \phi \cdot {\bf F}) (modal space)
    Array<OneD, NekDouble> tmp(nTotCoeffs,0.0);
    
    for (int i = 0; i < nVariables; ++i)
    {
	U->IProductWRTDerivBase(0,fluxVectorX[i],iprod_rhs[i],0);
	U->IProductWRTDerivBase(1,fluxVectorY[i],tmp,0);
	Vmath::Vadd(nTotCoeffs,tmp,1,iprod_rhs[i],1,iprod_rhs[i],1);
    }
    
    // Evaluate  upwind numerical flux (physical space)
    U->NumericalFlux(upwindX,upwindY);
    
    // Compute the boundary integral (modal space)
    for (int i = 0; i < nVariables; ++i)
    {
	Vmath::Neg(nTotTracePoints,upwindX[i],1);
	Vmath::Neg(nTotTracePoints,upwindY[i],1);
	
	U->AddTraceIntegral(upwindX[i],upwindY[i],iprod_rhs[i],i);
    }
    //---------------------------------------
    

    //---------------------------------------
    // Compute source terms (in this case the Coriolis force)
    Array<OneD, NekDouble> tmp1(nTotQuadPoints,0.0);

    Vmath::Vmul(nTotQuadPoints,U->GetCoriolis(),1,U->GetPhys(2),1,tmp1,1);
    U->IProductWRTBase(tmp1,tmp,0);
    Vmath::Vadd(nTotCoeffs,tmp,1,iprod_rhs[1],1,iprod_rhs[1],1);

    Vmath::Vmul(nTotQuadPoints,U->GetCoriolis(),1,U->GetPhys(1),1,tmp1,1);
    Vmath::Neg(nTotQuadPoints,tmp1,1);
    U->IProductWRTBase(tmp1,tmp,0);
    Vmath::Vadd(nTotCoeffs,tmp,1,iprod_rhs[2],1,iprod_rhs[2],1);
    //---------------------------------------

    
    //---------------------------------------
    // Solve the block-diagonal system
    for (int i = 0; i < nVariables; ++i)
    {
	U->MultiplyByElmtInvMass(iprod_rhs[i],U->UpdateCoeffs(i),i);
    }
    
    //---------------------------------------
    // Go from modal to physical space
    U->BwdTrans(rhs);
    
    //---------------------------------------
}

/**
* $Log: $
**/
