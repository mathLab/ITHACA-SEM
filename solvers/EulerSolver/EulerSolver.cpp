///////////////////////////////////////////////////////////////////////////////
//
// File EulerSolver.cpp
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
// Description: Euler equations solver
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <cmath> 

#include <../solvers/Auxiliary/EulerEquations.h>

using namespace Nektar;

void rhsFunction(EulerEquationsSharedPtr U,
		 Array<OneD, Array<OneD, NekDouble> > &rhs,
		 const NekDouble time);

void SetIsenTropicVortex(EulerEquationsSharedPtr U);
void GetExactIsenTropicVortex(EulerEquationsSharedPtr U, Array<OneD, NekDouble> &outarray, int field_no);
int main(int argc, char *argv[])
{

    if(argc != 2)
    {
        cerr << "Usage: EulerSolver meshfile" << endl;
        exit(1);
    }

    cout << "===============================================" <<endl;
    cout << "============= 2D EULER EQUATIONS ==============" <<endl;
    cout << "===============================================" <<endl;

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
        // Construct an object from the class SWEDisCont2D.
        EulerEquationsSharedPtr U =
            MemoryManager<EulerEquations>::AllocateSharedPtr(mesh,boundaryConds);
        //-----------------------------------------
     
     
        //-----------------------------------------------
        // Store input parameters
        U->SetTime(boundaryConds.GetParameter("T0"));
        U->SetTimeStep(boundaryConds.GetParameter("Dt"));
        U->SetSteps((int)boundaryConds.GetParameter("Steps"));
        U->SetCheckSteps((int)boundaryConds.GetParameter("Check"));
        U->SetGamma(boundaryConds.GetParameter("Gamma"));
        //-----------------------------------------------
     
    
        //-----------------------------------------
        // Read and evaluate the initial conditions
        // (stored in U)
        U->SetInitialConditions(boundaryConds, U->GetTime());

        SetIsenTropicVortex(U);


        // Write initial conditions to file
        stringstream outfileName_rho;
        outfileName_rho << "rho_initial.dat";
        ofstream outfile((outfileName_rho.str()).data());      
        U->WriteToFile(outfile,eTecplot,0);
        stringstream outfileName_rhou;
        outfileName_rhou << "rhou_initial.dat";
        ofstream outfile1((outfileName_rhou.str()).data());      
        U->WriteToFile(outfile1,eTecplot,1);
        stringstream outfileName_rhov;
        outfileName_rhov << "rhov_initial.dat";
        ofstream outfile2((outfileName_rhov.str()).data());      
        U->WriteToFile(outfile2,eTecplot,2);
        stringstream outfileName_E;
        outfileName_E << "E_initial.dat";
        ofstream outfile3((outfileName_E.str()).data());      
        U->WriteToFile(outfile3,eTecplot,3);
        //-----------------------------------------

     
        //-----------------------------------------------
        // MultiStep Arrays (should be "automatically" defined 
        // by the timeintegratorkey
     
        int nTotQuadPoints  = U->GetPointsTot();
     
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
            rhsFunction(U,f1,U->GetTime());
            Vmath::Svtvp(nTotQuadPoints,U->GetTimeStep(),f1[0],1,U0[0],1,U1[0],1);
            Vmath::Svtvp(nTotQuadPoints,U->GetTimeStep(),f1[1],1,U0[1],1,U1[1],1);
            Vmath::Svtvp(nTotQuadPoints,U->GetTimeStep(),f1[2],1,U0[2],1,U1[2],1);
            Vmath::Svtvp(nTotQuadPoints,U->GetTimeStep(),f1[3],1,U0[3],1,U1[3],1);
            U->SetPhys(U1);
	  
            //RK step two
            rhsFunction(U,f2,U->GetTime()+U->GetTimeStep());
            for(int j = 0; j < nTotQuadPoints; ++j)
	    {
                U1[0][j] = U0[0][j] + U->GetTimeStep() * (0.5*f1[0][j]+0.5*f2[0][j]);
                U1[1][j] = U0[1][j] + U->GetTimeStep() * (0.5*f1[1][j]+0.5*f2[1][j]);
                U1[2][j] = U0[2][j] + U->GetTimeStep() * (0.5*f1[2][j]+0.5*f2[2][j]);
                U1[3][j] = U0[3][j] + U->GetTimeStep() * (0.5*f1[3][j]+0.5*f2[3][j]);

	    }
            U->SetPhys(U1);
	  
            U->SetTime(U->GetTime()+U->GetTimeStep());
            //-----------------------------------------------

	  
            //-----------------------------------------------
            // Write chk files
            if (i%U->GetCheckSteps() == 0)
	    {  
                cout << "Time: " << U->GetTime() << endl;
                stringstream outfileName_rho;
                outfileName_rho << "rho_Quad_"<< chk <<".dat";
                ofstream outfile((outfileName_rho.str()).data());      
                U->WriteToFile(outfile,eTecplot,0);
	     
                stringstream outfileName_rhou;
                outfileName_rhou << "rhou_Quad_"<< chk <<".dat";
                ofstream outfile1((outfileName_rhou.str()).data());      
                U->WriteToFile(outfile1,eTecplot,1);
                stringstream outfileName_rhov;
                outfileName_rhov << "rhov_Quad_"<< chk <<".dat";
                ofstream outfile2((outfileName_rhov.str()).data());      
                U->WriteToFile(outfile2,eTecplot,2);
                stringstream outfileName_E;
                outfileName_E << "E_Quad_"<< chk <<".dat";
                ofstream outfile3((outfileName_E.str()).data());      
                U->WriteToFile(outfile3,eTecplot,3);
                ++chk;
	    }
            //-----------------------------------------------
	  
	  
            //-----------------------------------------------
            // Compute the L2 error
            if (i == U->GetSteps())
	    {
                Array<OneD, NekDouble> exactSolution(nTotQuadPoints);
                GetExactIsenTropicVortex(U,exactSolution,0);
	      
                MultiRegions::DisContField2DSharedPtr exactSolutionExp =
                    MemoryManager<MultiRegions::DisContField2D>::AllocateSharedPtr(mesh,boundaryConds);
	      
                exactSolutionExp->SetPhys(exactSolution);
	      
                NekDouble error = U->L2(*exactSolutionExp,0);
	      
                // Display the output              
                cout << "Ndof = " << U->GetNcoeffs() << " => Error = " << error << endl;
	    }
            //-----------------------------------------------
  
	}
     
    }
}

void rhsFunction(EulerEquationsSharedPtr U,
		 Array<OneD, Array<OneD, NekDouble> > &rhs,
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
    // U->SetBoundaryConditions(time);
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
    // Solve the block-diagonal system
    U->MultiplyByElmtInvMass(iprod_rhs[0],U->UpdateCoeffs(0),0);
    U->MultiplyByElmtInvMass(iprod_rhs[1],U->UpdateCoeffs(1),1);
    U->MultiplyByElmtInvMass(iprod_rhs[2],U->UpdateCoeffs(2),2);
    U->MultiplyByElmtInvMass(iprod_rhs[3],U->UpdateCoeffs(3),3);
	
    //---------------------------------------
    // Go from modal to physical space
    U->BwdTrans(rhs);
    
    //---------------------------------------
}

void SetIsenTropicVortex(EulerEquationsSharedPtr U)
{
    int nTotQuadPoints  = U->GetPointsTot();
  
    Array<OneD, NekDouble> rho(nTotQuadPoints,100.0);
    Array<OneD, NekDouble> rhou(nTotQuadPoints);
    Array<OneD, NekDouble> rhov(nTotQuadPoints);
    Array<OneD, NekDouble> E(nTotQuadPoints);
    Array<OneD, NekDouble> x(nTotQuadPoints);
    Array<OneD, NekDouble> y(nTotQuadPoints);
    Array<OneD, NekDouble> z(nTotQuadPoints);
  
    U->GetCoords(x,y,z);

    //---------------------------------
    // flow parameters
    NekDouble x0   = 5.0;
    NekDouble y0   = 0.0;
    NekDouble beta  = 5.0;
    NekDouble u0    = 1.0;
    NekDouble v0    = 0.0;
    NekDouble gamma = U->GetGamma();
    NekDouble time  = U->GetTime();
    NekDouble r;

    for (int i = 0; i < nTotQuadPoints; ++i)
    {
        r       = sqrt( pow(x[i]-u0*time-x0, 2.0) + pow(y[i]-v0*time-y0, 2.0));
        rho[i]  = pow( (1.0-((gamma-1.0)/(16.0*gamma*M_PI*M_PI))*beta*beta*exp(2.0*(1.0-r*r))), (1.0/(gamma-1.0)) );
        rhou[i] = rho[i] * (1.0 - beta*exp(1.0-r*r)*((y[i]-y0)/(2.0*M_PI)));
        rhov[i] = rho[i] * (beta*exp(1.0-r*r)*((x[i]-x0)/(2.0*M_PI)));
        E[i]    = (pow(rho[i],gamma)/(gamma-1.0)) + 0.5*rho[i]*(pow(rhou[i]/rho[i],2.0)+pow(rhov[i]/rho[i],2.0));
    }

    U->SetPhys(rho,0);
    U->SetPhys(rhou,1);
    U->SetPhys(rhov,2);
    U->SetPhys(E,3);
  
}


void GetExactIsenTropicVortex(EulerEquationsSharedPtr U, Array<OneD, NekDouble> &outarray, int field_no)
{
    int nTotQuadPoints  = U->GetPointsTot();
  
    Array<OneD, NekDouble> rho(nTotQuadPoints,100.0);
    Array<OneD, NekDouble> rhou(nTotQuadPoints);
    Array<OneD, NekDouble> rhov(nTotQuadPoints);
    Array<OneD, NekDouble> E(nTotQuadPoints);
    Array<OneD, NekDouble> x(nTotQuadPoints);
    Array<OneD, NekDouble> y(nTotQuadPoints);
    Array<OneD, NekDouble> z(nTotQuadPoints);
  
    U->GetCoords(x,y,z);
  
    //---------------------------------
    // flow parameters
    NekDouble x0   = 5.0;
    NekDouble y0   = 0.0;
    NekDouble beta  = 5.0;
    NekDouble u0    = 1.0;
    NekDouble v0    = 0.0;
    NekDouble gamma = U->GetGamma();
    NekDouble time  = U->GetTime();
    NekDouble r;

    for (int i = 0; i < nTotQuadPoints; ++i)
    {
        r       = sqrt( pow(x[i]-u0*time-x0, 2.0) + pow(y[i]-v0*time-y0, 2.0));
        rho[i]  = pow( (1.0-((gamma-1.0)/(16.0*gamma*M_PI*M_PI))*beta*beta*exp(2.0*(1.0-r*r))), (1.0/(gamma-1.0)) );
        rhou[i] = rho[i] * (1.0 - beta*exp(1.0-r*r)*((y[i]-y0)/(2.0*M_PI)));
        rhov[i] = rho[i] * (beta*exp(1.0-r*r)*((x[i]-x0)/(2.0*M_PI)));
        E[i]    = (pow(rho[i],gamma)/(gamma-1.0)) + 0.5*rho[i]*(pow(rhou[i]/rho[i],2.0)+pow(rhov[i]/rho[i],2.0));
    }

    switch (field_no){
    case 0:
        outarray = rho;
        break;
    case 1:
        outarray = rhou;
        break;
    case 2:
        outarray = rhov;
        break;
    case 3:
        outarray = E;
        break;
    }
  
}

/**
* $Log: EulerSolver.cpp,v $
* Revision 1.1  2008/08/22 09:48:23  pvos
* Added Claes' AdvectionDiffusionReaction, ShallowWater and Euler solver
*
**/
