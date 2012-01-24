///////////////////////////////////////////////////////////////////////////////
//
// File Advection1DFR.cpp
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
// Description: Demo to test the Flux Reconstruction 1D
// Unsteady advection (linear or spatially non linear advection) 
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <FluxReconstruction/Advection1DFR.h>

namespace Nektar
{
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// CONSTRUCTOR
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Advection1DFR::Advection1DFR(LibUtilities::SessionReaderSharedPtr &vSession)
	{
		// Read in mesh from input file and create an object of class MeshGraph1D
		// to encaplusate the mesh
		graph1D = MemoryManager<SpatialDomains::MeshGraph1D>::AllocateSharedPtr(vSession);
		
		/*const SpatialDomains::ExpansionMap &expansions = graph1D->GetExpansions(); 
		LibUtilities::BasisKey bkey0 = expansions.begin()->second->m_basisKeyVector[0];
		int Npoints = bkey0.GetNumModes();
		
		const LibUtilities::PointsKey FRpoints(Npoints,LibUtilities::eGaussGaussLegendre);
		const LibUtilities::BasisKey  FRBase(LibUtilities::eLagrange,Npoints,FRpoints);*/
		
		// Feed our spatial discretisation object with the information coming from the session file
		// i.e. initialise all the memory
		Domain = MemoryManager<MultiRegions::DisContField1D>::AllocateSharedPtr(vSession,graph1D,vSession->GetVariable(0));
		
		//number of total Gauss point
		nq = Domain->GetTotPoints();
		
		x = Array<OneD,NekDouble>(nq);
		y = Array<OneD,NekDouble>(nq);
		z = Array<OneD,NekDouble>(nq);

		Domain->GetCoords(x,y,z);
		
		// interface points
		ni = graph1D->GetNvertices();
		ne = Domain->GetExpSize();
		
		K = Array<OneD,NekDouble>(ni);
		
		xi = Array<OneD,NekDouble>(ni);
		yi = Array<OneD,NekDouble>(ni);
		zi = Array<OneD,NekDouble>(ni);
		
		for (int i = 0; i < ni; i++)
		{
			graph1D->GetVertex(i)->GetCoords(xi[i],yi[i],zi[i]);
		}
		
		//filling the correction function array
		GFtype = vSession->GetSolverInfo("Gfunctions");
		dgL = Array<OneD,NekDouble>(nq/ne);
		dgR = Array<OneD,NekDouble>(nq/ne);
		GFunctionsGrad(dgL,dgR);
		
		// Loading from the session file the functions describing the advection term,
		// the initial condition and the exact solution
		AdveFunc = vSession->GetFunction("Advection",0);
		InitCond = vSession->GetFunction("InitialConditions",0);
		ExSol = vSession->GetFunction("ExactSolution",0);
		
		// Loading time-stepping parameters
		vSession->LoadParameter("InitialTime",Time,0.0);
		vSession->LoadParameter("TimeStep",TimeStep,0.01);
		vSession->LoadParameter("NumSteps",NumSteps,0);
		
		RiemSol = vSession->GetSolverInfo("ReimannSolver");
		
		// Filling two vectors with the value of the prescibed functions at the quadrature points
		Adv  = Array<OneD, Array<OneD, NekDouble> >(1);
		Ini  = Array<OneD, Array<OneD, NekDouble> >(1);
		Exac = Array<OneD, Array<OneD, NekDouble> >(1);
		
		Adv[0]  = Array<OneD,NekDouble>(nq);
		Ini[0]  = Array<OneD,NekDouble>(nq);
		Exac[0] = Array<OneD,NekDouble>(nq);
		
		for(int i = 0; i < nq; ++i)
		{
			Adv[0][i] = AdveFunc->Evaluate(x[i],y[i],z[i],Time);
			Ini[0][i] = InitCond->Evaluate(x[i],y[i],z[i],Time);
		}
		
		// Setting the physical value of the initial condition inside the physical space
		Domain->UpdatePhys() = Ini[0];
		// Setting also the coefficients consequently
		//Domain->FwdTrans(Ini[0],Domain->UpdateCoeffs());
	}
	
	//disctructor
	Advection1DFR::~Advection1DFR()
	{
	}
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// FUNCTIONs IMPLEMENTATION
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	void Advection1DFR::EvaluateAdvectionTerm(const Array<OneD, Array<OneD, NekDouble> > & inarray,
											  Array<OneD, Array<OneD, NekDouble> > & outarray,
											  const NekDouble time)
	{
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// inarray is intended to be the solution U at time n, so this RHS calculation we will be used to work out U at time n+1
		// which will then used to feed this fuction at the next step
		// outarray is our advection term (A*dU/dx) corrected with the flux recostruction technique
		// Note: The value of the Advection coefficients are stored in the vector Adv[i].
		// This is because we can have non-constant advection (some spatial non lineraity).
		// The advection coefficients could also be function of time, then the evaluation of Adv[i]
		// must be redone here as we did in the constructor with the function "Evaluate".
		// For a full non-linear advection where A=U, we need just to substitute in the Vmath function "Adv" with "inarray"
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// switching to Peter nomenclature where inarray = Ud
		
		Array<OneD,NekDouble> fd(nq,0.0);
		Array<OneD,NekDouble> gradfd(nq,0.0);
		
		Array<OneD,NekDouble> udi(2*ne,0.0);
		Array<OneD,NekDouble> fdi(2*ne,0.0);
		Array<OneD,NekDouble> fi(ni,0.0);
		
		Array<OneD,NekDouble> tmp1,tmp2;
		
		//calculating the discontinous flux fd = Adv*Ud
		Vmath::Vmul(nq,Adv[0],1,inarray[0],1,fd,1);
		
		LibUtilities::BasisSharedPtr Basis;
		
		Basis = Domain->GetExp(0)->GetBasis(0);
		
		StdRegions::StdSegExp StdSeg(Basis->GetBasisKey());
		
		for(int i = 0; i < ne ; i++)
		{
			StdSeg.PhysDeriv(tmp1 = fd + i*nq/ne,tmp2 = gradfd + i*nq/ne);
		}
		
		//calculate the value of the fd at the interface points ==> fdi
		InterpToInterface(fd,fdi);
		
		//calculate the value of ud at the interface points ==> udi
		InterpToInterface(inarray[0],udi);
		
		//calculate the interface fluxes fi
		ReimannSolver(udi,fi);
		
		//calculate correction final flux gradient df/dx and putting it into the output vector ==> outarray
		FluxesReconstruction(gradfd,fdi,fi,outarray[0]);
		
		// negate advection term
		Vmath::Neg(nq,outarray[0],1);
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Advection1DFR::Projection(const Array<OneD, Array<OneD, NekDouble> > & inarray,
								   Array<OneD, Array<OneD, NekDouble> > & outarray,
								   const NekDouble time) const
	{
		// For DG it is just a copy
		Vmath::Vcopy(nq,inarray[0],1,outarray[0],1);
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Advection1DFR::ReimannSolver(const Array<OneD, NekDouble> & inarray, 
									  Array<OneD, NekDouble> & F)
	{
		Array<OneD,NekDouble> a,sumInt,difInt;
		a = Array<OneD,NekDouble>(ni);
		sumInt = Array<OneD,NekDouble>(ni);
		difInt = Array<OneD,NekDouble>(ni);
		
		sumInt[0] = inarray[0]+inarray[2*ne-1];
		sumInt[ni-1] = sumInt[0];
		difInt[0] = inarray[0]-inarray[2*ne-1];
		difInt[ni-1] = difInt[0];
		
		for(int i = 1; i < ni-1; i++)
		{
			sumInt[i] = inarray[2*i] + inarray[2*i-1];
			difInt[i] = inarray[2*i] - inarray[2*i-1];
		}
	
		for(int i = 0; i < ni; ++i)
		{
			a[i] = AdveFunc->Evaluate(xi[i],yi[i],zi[i],Time);
			
			if(RiemSol == "Up-Wind"){K[i] = 0.0;}
			else if(RiemSol == "Euler-Centered"){K[i] = 1.0;}
			else if(RiemSol == "Lax-Friedrichs"){K[i] = 1.0 - 1.0/(TimeStep*abs(a[i]));}
			else if(RiemSol == "Lax-Wendroff"){K[i] = 1.0 - (TimeStep*abs(a[i]));}
			else{ASSERTL0(false,"Reimann solver not implemented");}
			
			F[i] = 0.5*a[i]*(sumInt[i]) - 0.5*abs(a[i])*(1.0-K[i])*(difInt[i]);
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Advection1DFR::FluxesReconstruction(const Array<OneD, NekDouble> & fd,
											 const Array<OneD, NekDouble> & fdi,
											 const Array<OneD, NekDouble> & fi,
											 Array<OneD, NekDouble> & outarray)
	{
		Array<OneD,NekDouble> jL(ni-1,0.0);  // jumps on the left interfaces
		Array<OneD,NekDouble> jR(ni-1,0.0);  // jumps on the righ interfaces
		
		Array<OneD,NekDouble> tmpDGL(nq/ne,0.0); // temp arrays containing j*dg/dx
		Array<OneD,NekDouble> tmpDGR(nq/ne,0.0);
		Array<OneD,NekDouble> tmp,tmparray;
		
		// calculating the jumps on the left and right side
		for(int i = 0; i< ne; i++)
		{
			jL[i] = fi[i] - fdi[2*i];
			jR[i] = fi[i+1] - fdi[2*i+1];
		}

		for(int i = 0; i < ne; i++)
		{
			Vmath::Smul(nq/ne,jL[i],tmp = dgL,1,tmpDGL,1);
			Vmath::Smul(nq/ne,jR[i],tmp = dgR,1,tmpDGR,1);
			Vmath::Vadd(nq/ne,tmp = fd+i*nq/ne,1,tmpDGL,1,tmparray = outarray +i*nq/ne,1);
			Vmath::Vadd(nq/ne,tmparray = outarray +i*nq/ne,1,tmpDGR,1,tmparray = outarray +i*nq/ne,1);
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Advection1DFR::GFunctionsGrad(Array<OneD, NekDouble> & dGL,
									   Array<OneD, NekDouble> & dGR)
	{	
		LibUtilities::BasisSharedPtr Basis;
		
		Basis = Domain->GetExp(0)->GetBasis(0);
		
		StdRegions::StdSegExp StdSeg(Basis->GetBasisKey());
		
		int k  = Basis->GetNumModes();
		int np = Basis->GetNumPoints();
		
		Array<OneD,NekDouble> zeros(np,0.0);
		
		zeros = Basis->GetZ();
		
		NekDouble sign = pow(-1.0,double(k));
		NekDouble etak;
		
		if(GFtype == "DG"){ etak = 0.0; }
		else if(GFtype == "SD"){ etak = k/(1.0+k); }
		else if(GFtype == "HU"){ etak = (1.0+k)/k; }
		else {ASSERTL0(false,"options for the g functions are DG, SD and HU");}
		
		NekDouble overeta = 1.0/(1.0+etak);
						
		Array<OneD,NekDouble> Lk(np,0.0);
		Array<OneD,NekDouble> dLk(np,0.0);
		Array<OneD,NekDouble> Lkp1(np,0.0);
		Array<OneD,NekDouble> dLkp1(np,0.0);
		Array<OneD,NekDouble> Lkm1(np,0.0);
		Array<OneD,NekDouble> dLkm1(np,0.0);
		Array<OneD,NekDouble> GL(np,0.0);
		Array<OneD,NekDouble> GR(np,0.0);
		
		Polylib::jacobfd(np,&(zeros[0]),&(Lk[0]),&(dLk[0]),k,0.0,0.0);
		Polylib::jacobfd(np,&(zeros[0]),&(Lkm1[0]),&(dLkm1[0]),k-1,0.0,0.0);
		Polylib::jacobfd(np,&(zeros[0]),&(Lkp1[0]),&(dLkp1[0]),k+1,0.0,0.0);
	
		//dgL
		Vmath::Smul(np,etak,Lkm1,1,GL,1);
		Vmath::Vadd(np,GL,1,Lkp1,1,GL,1);
		Vmath::Smul(np,overeta,GL,1,GL,1);
		Vmath::Vsub(np,Lk,1,GL,1,GL,1);
		Vmath::Smul(np,0.5*sign,GL,1,GL,1);
		//dgR
		Vmath::Smul(np,etak,Lkm1,1,GR,1);
		Vmath::Vadd(np,GR,1,Lkp1,1,GR,1);
		Vmath::Smul(np,overeta,GR,1,GR,1);
		Vmath::Vadd(np,Lk,1,GR,1,GR,1);
		Vmath::Smul(np,0.5,GR,1,GR,1);
		
		StdSeg.PhysDeriv(GL,dGL);
		StdSeg.PhysDeriv(GR,dGR);
		

		np = 1000;
		NekDouble dx = 2.0/(np-1);
		Array<OneD,NekDouble> points(np);
		points[0] = -1.0;
		for (int i=1; i < np ; i++)
		{
			points[i] = points[i-1] + dx;
		}
			
		Array<OneD,NekDouble> plot_Lk(np,0.0);
		Array<OneD,NekDouble> plot_dLk(np,0.0);
		Array<OneD,NekDouble> plot_Lkp1(np,0.0);
		Array<OneD,NekDouble> plot_dLkp1(np,0.0);
		Array<OneD,NekDouble> plot_Lkm1(np,0.0);
		Array<OneD,NekDouble> plot_dLkm1(np,0.0);
		Array<OneD,NekDouble> plot_GL(np,0.0);
		Array<OneD,NekDouble> plot_GR(np,0.0);
			
		Polylib::jacobfd(np,&(points[0]),&(plot_Lk[0]),&(plot_dLk[0]),k,0.0,0.0);
		Polylib::jacobfd(np,&(points[0]),&(plot_Lkm1[0]),&(plot_dLkm1[0]),k-1,0.0,0.0);
		Polylib::jacobfd(np,&(points[0]),&(plot_Lkp1[0]),&(plot_dLkp1[0]),k+1,0.0,0.0);
			
			//dgL
		Vmath::Smul(np,etak,plot_Lkm1,1,plot_GL,1);
		Vmath::Vadd(np,plot_GL,1,plot_Lkp1,1,plot_GL,1);
		Vmath::Smul(np,overeta,plot_GL,1,plot_GL,1);
		Vmath::Vsub(np,plot_Lk,1,plot_GL,1,plot_GL,1);
		Vmath::Smul(np,0.5*sign,plot_GL,1,plot_GL,1);
		
		ofstream outfile;
		
		outfile.open("GL.dat");
		for(int i = 0; i < np; i++)
		{
			outfile << scientific 
			<< setw (17) 
			<< setprecision(10) 
			<< points[i]
			<< "  " 
			<< plot_GL[i]  
			<< endl;
		}
		outfile << endl << endl;
		outfile.close();
		
	}
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Advection1DFR::InterpToInterface(const Array<OneD, NekDouble> & total,
										  Array<OneD, NekDouble> & interfaceValue)
	{
		LibUtilities::BasisSharedPtr Basis;
		Basis = Domain->GetExp(0)->GetBasis(0);
		
		StdRegions::StdSegExp StdSeg(Basis->GetBasisKey());
		
		Array<OneD,NekDouble> coordL(3,0.0);
		Array<OneD,NekDouble> coordR(3,0.0);
		Array<OneD,NekDouble> tmp(nq/ne,0.0);
		coordL[0] = -1.0;
		coordR[0] = 1.0;
		int cnt = 0;
		for(int i=0; i < ne; i++)
		{
			tmp = total + i*nq/ne;
			interfaceValue[cnt] = StdSeg.PhysEvaluate(coordL,tmp);cnt++;
			interfaceValue[cnt] = StdSeg.PhysEvaluate(coordR,tmp);cnt++;
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Advection1DFR::AppendOutput(const Array<OneD, NekDouble> & approx,
									 const Array<OneD, NekDouble> & exact) const
	{
		ofstream outfile;
		outfile.open("Solution.dat");
		for(int i = 0; i < nq; i++)
		{
			outfile << scientific 
			<< setw (17) 
			<< setprecision(10) 
			<< x[i]
			<< "  " 
			<< approx[i] 
			<< "  " 
			<< exact[i] 
			<< endl;
		}
		outfile << endl << endl;
		outfile.close();
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Advection1DFR::GenerateGnuplotScript() const
	{
		ofstream outfile;
		outfile.open("FRSolution.p");
			
		outfile << "# Gnuplot script file" << endl;
		outfile << "set   autoscale" << endl;                       
		outfile << "unset log" << endl;                           
		outfile << "unset label" << endl;                          
		outfile << "set xtic auto" << endl;                    
		outfile << "set ytic auto" << endl;                        
		outfile << "set xlabel \"x\"" << endl;
		outfile << "set ylabel \"u\"" << endl;

		outfile << "plot    \"Solution.dat\" ";
		outfile << "using 1:2 index 0 ";
		outfile << " title 'FR solution (t = " << Time << ")' with linespoints lt 3, ";
		outfile << "\"Solution.dat\" ";
		outfile << "using 1:3 index 0 ";
		outfile << " title 'Exact Solution (t = " << Time << ")' with linespoints lt 1" << endl;

		outfile.close();
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Advection1DFR::SolutionPrint()
	{
		/////////////////////////////////////////////////////////////
		// Evaluating the error
		NekDouble L2, Linf;
		for(int i = 0; i < nq; ++i)
		{
			Exac[0][i] = ExSol->Evaluate(x[i],y[i],z[i],Time);
		}
		
		L2   = Domain->L2(Exac[0]);
		Linf = Domain->Linf(Exac[0]);
		
		/////////////////////////////////////////////////////////////
		// Write solution to file
		AppendOutput(Domain->GetPhys(),Exac[0]);
		GenerateGnuplotScript();
		/////////////////////////////////////////////////////////////
		// Print summary of solution details
		const SpatialDomains::ExpansionMap &expansions = graph1D->GetExpansions(); 
		LibUtilities::BasisKey bkey0 = expansions.begin()->second->m_basisKeyVector[0];
		
		cout << endl;
		cout << "Solving Unsteady 1D Advection with Flux Recostruction"  << endl;
		
		cout << endl;
		cout << "Expansion : " << LibUtilities::BasisTypeMap[bkey0.GetBasisType()] << endl;
		cout << "No. modes : " << bkey0.GetNumModes() << endl;
		cout << "Advection : " << AdveFunc->GetEquation() << endl;
		
		cout << endl;
		
		cout << "              Time Step : " << TimeStep << endl;
		cout << "        Number of Steps : " << NumSteps << endl;
		cout << "      Initial Condition : " << InitCond->GetEquation() << endl;
		
		cout << endl;
		cout << "Exact Solution : " << ExSol->GetEquation() << endl;
		cout << "    Linf error : " << L2 << endl;
		cout << "    L2   error : " << Linf << endl;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	MultiRegions::DisContField1DSharedPtr Advection1DFR::GetDomain() const
	{
		return Domain;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	NekDouble Advection1DFR::GetTime() const
	{
		return Time;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Advection1DFR::UpdateTime()
	{
		Time = Time + TimeStep;
		
		return;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	NekDouble Advection1DFR::GetTimeStep() const
	{
		return TimeStep;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int Advection1DFR::GetNumSteps() const
	{
		return NumSteps;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int Advection1DFR::GetNumPoints() const
	{
		return nq;
	}
	//////////////////////////////////////////////////////////////////////////////	
}