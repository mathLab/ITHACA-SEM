///////////////////////////////////////////////////////////////////////////////
//
// File AdvectionDiffusionReactionSolver.cpp
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
// Description: Advection Diffusion Reaction solver
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <cmath> 

#include <../solvers/Auxiliary/AdvectionDiffusionReaction.h>

using namespace Nektar;

//----------------------------------------------
// local functions
void rhsFunction(AdvectionDiffusionReactionSharedPtr hpExp, 
		 Array<OneD, Array<OneD, NekDouble> > v,
		 Array<OneD, Array<OneD, NekDouble> > vTrace,
		 Array<OneD, const NekDouble> u0,
		 Array<OneD, NekDouble> &rhs, NekDouble time);
//----------------------------------------------

int main(int argc, char *argv[])
{

    if(argc != 2)
    {
        fprintf(stderr,"Usage: AdvectionDiffusionReactionSolver  meshfile \n");
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
        AdvectionDiffusionReactionSharedPtr u =
            MemoryManager<AdvectionDiffusionReaction>::AllocateSharedPtr(mesh,boundaryConds);
        //-----------------------------------------
    
        //-----------------------------------------
        // Create arrays of arrays
    
        int nTotQuadPoints  = u->GetPointsTot();
        int nTotTracePoints = u->GetNpoints();
    
        Array<OneD, Array<OneD, NekDouble> > v(2);
        Array<OneD, Array<OneD, NekDouble> > vTrace(2);
    
        for (int i = 0; i < 2; ++i)
        {
            v[i] = Array<OneD, NekDouble>(nTotQuadPoints);
            vTrace[i] = Array<OneD, NekDouble>(nTotTracePoints);
        }
      
        //-----------------------------------------
        // Read and evaluate the initial conditions
        Array<OneD,NekDouble> x1(nTotQuadPoints,0.0);
        Array<OneD,NekDouble> x2(nTotQuadPoints,0.0);

        u->GetCoords(x1,x2);
    
        Array<OneD, NekDouble> u0(nTotQuadPoints);
        Array<OneD, NekDouble> u1(nTotQuadPoints);
    
        SpatialDomains::ConstForcingFunctionShPtr u0SolutionEquation 
            = boundaryConds.GetExactSolution(boundaryConds.GetVariable(0));

        for(int i = 0; i < nTotQuadPoints; ++i)
	{
            u0[i] = u0SolutionEquation->Evaluate(x1[i],x2[i]);
            v[0][i] =  2.0*M_PI*x2[i];
            v[1][i] = -2.0*M_PI*x1[i];
	}      

        u->SetPhys(u0,0);
        //-----------------------------------------------

    
        //-----------------------------------------------
        // Initiate the velocities on the trace
        u->SetPhys(v[0],0);
        u->ExtractTracePhys(vTrace[0],0);

        u->SetPhys(v[1],0);
        u->ExtractTracePhys(vTrace[1],0);
    
        u->SetPhys(u0,0);
        //-----------------------------------------------


        //-----------------------------------------------
        // Write initial conditions to file
        stringstream outfileName;
        outfileName << "Advection_Quad_P9_initial.dat";
        ofstream outfile((outfileName.str()).data());      
        u->WriteToFile(outfile,eTecplot,0);
        //-----------------------------------------------
   

        //-----------------------------------------------
        // RK vectors
        Array<OneD, NekDouble> f1(nTotQuadPoints,0.0);
        Array<OneD, NekDouble> f2(nTotQuadPoints,0.0);
        //-----------------------------------------------

   
        //-----------------------------------------------
        // Time stepping parameters
        NekDouble dt = boundaryConds.GetParameter("Dt");
        int timeStep = (int)boundaryConds.GetParameter("Steps");
        int chkStep  = (int)boundaryConds.GetParameter("Check");
        //-----------------------------------------------


        //-----------------------------------------------
        // start the time loop
        int chk = 0;
        for (int i = 1; i < timeStep + 1; ++i)
	{
	  
            // Store the starting value
            u0 = u->GetPhys(0);
	  
	  
            //RK step one
            rhsFunction(u,v,vTrace,u0,f1,i*dt);
            Vmath::Svtvp(nTotQuadPoints,dt,f1,1,u0,1,u1,1);
            u->SetPhys(u1,0);
	
	 
            //RK step two
            rhsFunction(u,v,vTrace,u1,f2,i*dt+dt);
            for(int j = 0; j < nTotQuadPoints; ++j)
	    {
                u1[j] = u0[j] + dt *(0.5*f1[j]+0.5*f2[j]);
	    }
            u->SetPhys(u1,0);


            //-----------------------------------------------
            // Write chk files
            if (i%chkStep == 0)
	    {  
                cout << "Time: " << i*dt << endl;
                stringstream outfileName;
                outfileName << "Advection_Quad_P9_"<< chk <<".dat";
                ofstream outfile((outfileName.str()).data());      
                u->WriteToFile(outfile,eTecplot,0);
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
                for(int j = 0; j < nTotQuadPoints; ++j)
		{
                    exactSolution[j] = exactSolutionEquation->Evaluate(x1[j],x2[j]);
		}      
	      
                MultiRegions::ExpList2DSharedPtr exactSolutionExp =
                    MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(mesh);
                exactSolutionExp->SetPhys(exactSolution);
	      
                NekDouble error = u->L2(*exactSolutionExp,0);
	      
                // Display the output              
                cout << "Ndof = " << u->GetNcoeffs() << " => L2 error = " << error << endl;
	    }
            //-----------------------------------------------

	}

    }

}


//-----------------------------------------
void rhsFunction(const AdvectionDiffusionReactionSharedPtr hpExp,
            	 Array<OneD, Array<OneD, NekDouble> > v,
		 Array<OneD, Array<OneD, NekDouble> > vTrace,
		 Array<OneD, const NekDouble> u0, 
		 Array<OneD, NekDouble> &rhs, NekDouble time)
{
  
    int nTotQuadPoints  = hpExp->GetPointsTot();
    int nTotCoeffs      = hpExp->GetNcoeffs();
    int nTotTracePoints = hpExp->GetNpoints();
 
    //---------------------------------------
    // get the trace values

    Array<OneD,NekDouble> Fwd(nTotTracePoints,0.0);
    Array<OneD,NekDouble> Bwd(nTotTracePoints,0.0);
    Array<OneD,NekDouble> upwind1(nTotTracePoints,0.0);
    Array<OneD,NekDouble> upwind2(nTotTracePoints,0.0);

    Array<OneD, Array<OneD, NekDouble> > normals(2);
    normals[0] = Array<OneD, NekDouble>(nTotTracePoints);
    normals[1] = Array<OneD, NekDouble>(nTotTracePoints);
   

    hpExp->GetFwdBwdTracePhys(Fwd,Bwd,0);
    hpExp->GetTraceNormals(normals); 

    //---------------------------------------


    //---------------------------------------
    // Compute the flux vector

    Array<OneD, NekDouble> fluxVector0(nTotQuadPoints,0.0);
    Array<OneD, NekDouble> fluxVector1(nTotQuadPoints,0.0);

    // Fill the flux vector using the old time step
    for (int i = 0; i < nTotQuadPoints; ++i)
    {
        fluxVector0[i] = v[0][i]*u0[i];
        fluxVector1[i] = v[1][i]*u0[i];
    }
    //---------------------------------------


    Array<OneD, NekDouble> iprod_rhs(nTotCoeffs,0.0);
    Array<OneD, NekDouble> tmp(nTotCoeffs,0.0);
    Array<OneD, NekDouble> iprod_rhs_test(nTotCoeffs,0.0);
    //--------------------------------------


    //---------------------------------------
    // Compute the (Grad \phi \cdot F)

    hpExp->IProductWRTDerivBase(0,fluxVector0,iprod_rhs,0);
    hpExp->IProductWRTDerivBase(1,fluxVector1,tmp,0);
    Vmath::Vadd(nTotCoeffs,tmp,1,iprod_rhs,1,iprod_rhs,1);
    //---------------------------------------


    //---------------------------------------
    // Evaluate  upwind numerical flux
    // note: upwind1 contains  a_1 * upwind
    //       upwind2 contains  a_2 * upwind

#if 1
    Array<OneD, Array<OneD, const NekDouble> > Vel(2);
    Vel[0] = vTrace[0];
    Vel[1] = vTrace[1];

    hpExp->UpwindTrace(Vel,Fwd,Bwd,upwind1);

    Vmath::Vmul(nTotTracePoints,upwind1,1,vTrace[1],1,upwind2,1);
    Vmath::Vmul(nTotTracePoints,upwind1,1,vTrace[0],1,upwind1,1);
#else
    for (int i = 0; i < nTotTracePoints; ++i)
    {
        if (vTrace[0][i]*normals[0][i]+vTrace[1][i]*normals[1][i] > 0.0)
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

    hpExp->AddTraceIntegral(upwind1, upwind2, iprod_rhs, 0);

    hpExp->MultiplyByElmtInvMass(iprod_rhs,hpExp->UpdateCoeffs(0),0);


    //---------------------------------------
    // Go to physical space
    hpExp->BwdTrans(hpExp->GetCoeffs(0),rhs,0);
    //---------------------------------------
}

/**
* $Log: AdvectionDiffusionReactionSolver.cpp,v $
* Revision 1.1  2008/08/22 09:48:23  pvos
* Added Claes' AdvectionDiffusionReaction, ShallowWater and Euler solver
*
**/
