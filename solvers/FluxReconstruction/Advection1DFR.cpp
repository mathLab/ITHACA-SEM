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
// Description: Demo to test the Flux Reconstruction scheme 1D
// Unsteady advection problem (linear or spatially non linear advection) 
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <FluxReconstruction/Advection1DFR.h>

namespace Nektar
{
    /// Constructor
	Advection1DFR::Advection1DFR(LibUtilities::SessionReaderSharedPtr &vSession)
	{
		/// Reading mesh from the session file
		graph1D = MemoryManager<SpatialDomains::MeshGraph1D>::AllocateSharedPtr(vSession);
		
		/// Feeding spatial discretisation using the information from the session file
		Domain = MemoryManager<MultiRegions::DisContField1D>::AllocateSharedPtr(vSession,graph1D, vSession->GetVariable(0));

        /// Number of elements
		ne = Domain->GetExpSize();
        
        /// Total number of quadrature points
        /// For nodal expansion the number of quadrature points is equal to the number of coefficients
		nq = Domain->GetTotPoints();

        /// Coordinates of the quadrature points in the real physical space
		x = Array<OneD,NekDouble>(nq);
		y = Array<OneD,NekDouble>(nq);
		z = Array<OneD,NekDouble>(nq);
		Domain->GetCoords(x, y, z);
        
		/// Number of interface points
		ni = graph1D->GetNvertices();
        
        /// Coefficient for the Riemann solver
		K = Array<OneD,NekDouble>(ni);
		
        /// Coordinates of the interface points in the real physical space
		xi = Array<OneD,NekDouble>(ni);
		yi = Array<OneD,NekDouble>(ni);
		zi = Array<OneD,NekDouble>(ni);
		for (int i = 0; i < ni; i++)
		{
			graph1D->GetVertex(i)->GetCoords(xi[i], yi[i], zi[i]);
		}

		/// Filling the derivatives of the correction functions arrays
		GFtype = vSession->GetSolverInfo("Gfunctions");
		dgL = Array<OneD,NekDouble>(nq/ne);
		dgR = Array<OneD,NekDouble>(nq/ne);
		GFunctionsGrad(dgL, dgR);
		
		/// Loading from the session file the advection term
		AdveFunc    = vSession->GetFunction("Advection",        0);
        
		/// Loading from the session file the initial conditions
		InitCond    = vSession->GetFunction("InitialConditions",0);
        
		/// Loading from the session file the exact solution
		ExSol       = vSession->GetFunction("ExactSolution",    0);
		
		/// Loading from the session file the time-stepping parameters
		vSession->LoadParameter("InitialTime",  Time,       0.0);
		vSession->LoadParameter("TimeStep",     TimeStep,   0.01);
		vSession->LoadParameter("NumSteps",     NumSteps,   0);
		
        /// Loading from the session file the Riemann solver
		RiemSol = vSession->GetSolverInfo("RiemannSolver");
		
		/// Filling the Adv vector with the value of the advection function at the quadrature points
		Adv  = Array<OneD, Array<OneD, NekDouble> >(1);
        Adv[0]  = Array<OneD,NekDouble>(nq);
        AdveFunc->Evaluate(x, y, z, Time, Adv[0]);

		/// Filling the Ini vector with the value of the initial condition at the quadrature points
		Ini  = Array<OneD, Array<OneD, NekDouble> >(1);
        Ini[0]  = Array<OneD,NekDouble>(nq);
        InitCond->Evaluate(x, y, z, Time, Ini[0]);

        /// Filling the Exac vector with the value of the exact solution at the quadrature points
		Exac = Array<OneD, Array<OneD, NekDouble> >(1);
		Exac[0] = Array<OneD,NekDouble>(nq);
		
		/// Setting the physical value of the initial condition inside the physical space
		Domain->UpdatePhys() = Ini[0];
        
		/// Setting the coefficients of the initial condition consequently
		Domain->FwdTrans(Ini[0], Domain->UpdateCoeffs());
	}
	
	/// Destructor
	Advection1DFR::~Advection1DFR()
	{
	}
	

    
    //! 1. METHOD 1: EvaluateAdvectionTerm
    void Advection1DFR::EvaluateAdvectionTerm(const Array<OneD, Array<OneD, NekDouble> > & inarray,
											  Array<OneD, Array<OneD, NekDouble> > & outarray,
											  const NekDouble time)
	{
		/// *****************************************************************************************
		/// inarray is intended to be the solution U at time n, so this RHS calculation we will be *
		/// used to work out U at time n+1 which will then used to feed this fuction at the next   *
		/// step outarray is our advection term (A*dU/dx) corrected with the flux recostruction    *
		/// technique. Note: The value of the Advection coefficients are stored in the vector      *
		/// Adv[i]. This is because we can have non-constant advection (some spatial nonlinearity).*
		/// The advection coefficients could also be function of time, then the evaluation of      *
		/// Adv[i] must be redone here as we did in the constructor with the function "Evaluate".  *
		/// For a full non-linear advection where A=U, we need just to substitute in the Vmath     *
        /// function "Adv" with "inarray".                                                         *
        /// ****************************************************************************************
		
		
        
        /// 1.1) Switching to Peter nomenclature where inarray = Ud 
        /// fd: discountinuous flux (A * Ud) ---> see P.Vincent papers for nomenclature 		
		Array<OneD,NekDouble> fd(nq,0.0);
		
        /// gradfd: \grad(fd), gradient of the discountinuous flux
        Array<OneD,NekDouble> gradfd(nq,0.0);
		
        /// udi: discountinuous solution at the interfaces
		Array<OneD,NekDouble> udi(2*ne,0.0);
		
        /// fdi: discontinuous flux fd at the interfaces
        Array<OneD,NekDouble> fdi(2*ne,0.0);
		
        /// fi = interface flux computed via Riemann solver
        Array<OneD,NekDouble> fi(ni,0.0);
		
        /// tmp1, tmp2: auxiliary variables
		Array<OneD,NekDouble> tmp1,tmp2;
		
        
		/// 1.3) Calculating the discontinous flux fd = (Adv * Ud) = (A * Ud)
		Vmath::Vmul(nq, Adv[0], 1, inarray[0], 1, fd, 1);
		
        
		/// 1.4) Taking the gradient of gradfd = d(fd)/dx on the reference element
		LibUtilities::BasisSharedPtr Basis;
		Basis = Domain->GetExp(0)->GetBasis(0);
		StdRegions::StdSegExp StdSeg(Basis->GetBasisKey());
		for(int i=0; i < ne; i++)
		{
			StdSeg.PhysDeriv(tmp1 = fd + i*nq/ne, tmp2 = gradfd + i*nq/ne);
		}
		

		/// 1.5) Calculate the value of the fd at the interface points ==> fdLi and fdRi 
		InterpToInterface(fd,fdi);
        
        
		/// 1.6) Calculate the value of ud at the interface points ==> udi 
		InterpToInterface(inarray[0], udi);
        
        
		/// 1.7) Calculate the interface fluxes fi
		RiemannSolver(udi, fi);
        		
        //SpatialDomains::GeomFactors1D::GeomFactors1D(eRegular, 1, xi, Basis); 
        
		/// 1.8) Calculate the final correction flux gradient ==> (df/dx) and putting 
        /// it into the output vector ==> outarray
		FluxesReconstruction(gradfd, fdi, fi, outarray[0]);
		
        Vmath::Smul(nq, 10.0, outarray[0], 1.0, outarray[0], 1.0);
		
        /// 1.9) Negate advection term
		Vmath::Neg(nq, outarray[0], 1);
	}


    
    //! 2. METHOD 2: Projection
	void Advection1DFR::Projection(const Array<OneD, Array<OneD, NekDouble> > & inarray,
								   Array<OneD, Array<OneD, NekDouble> > & outarray,
								   const NekDouble time) const
	{
		/// 2.1) For DG it is just a copy
		Vmath::Vcopy(nq,inarray[0],1,outarray[0],1);
	}
    	
    
    
    //! 3. METHOD 3: RiemannSolver
    void Advection1DFR::RiemannSolver(const Array<OneD, NekDouble> & inarray, 
									  Array<OneD, NekDouble> & F)
	{
		Array<OneD,NekDouble> a,sumInt,difInt;
        
		a       = Array<OneD,NekDouble>(ni);
		sumInt  = Array<OneD,NekDouble>(ni);
		difInt  = Array<OneD,NekDouble>(ni);
		
		sumInt[0]       = inarray[0]+inarray[2*ne-1];
		sumInt[ni-1]    = sumInt[0];
        
		difInt[0]       = inarray[0]-inarray[2*ne-1];
		difInt[ni-1]    = difInt[0];
		
		for(int i = 1; i < ni-1; i++)
		{
			sumInt[i] = inarray[2*i] + inarray[2*i-1];
			difInt[i] = inarray[2*i] - inarray[2*i-1];
		}
	
        AdveFunc->Evaluate(xi,yi,zi,Time,a);
		for(int i = 0; i < ni; ++i)
		{
			
			if      (RiemSol == "Up-Wind")       {K[i] = 0.0;}
			else if (RiemSol == "Euler-Centered"){K[i] = 1.0;}
			else if (RiemSol == "Lax-Friedrichs"){K[i] = 1.0 - 0.07/(TimeStep*abs(a[i]));}
			else if (RiemSol == "Lax-Wendroff")  {K[i] = 1.0 - (TimeStep*abs(a[i]))/(0.07);}
			else                                 {ASSERTL0(false,"Riemann solver not implemented");}
			
			F[i] = 0.5*a[i]*(sumInt[i]) - 0.5*abs(a[i])*(1.0-K[i])*(difInt[i]);
		}
	}

    
    
    //! 4. METHOD 4: FluxesReconstruction
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
		
		/// Calculating the jumps on the left and right side
		for(int i = 0; i< ne; i++)
		{
			jL[i] = fi[i]   - fdi[2*i];
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

    
    
    //! 5. METHOD 5: GFunctionsGrad
    void Advection1DFR::GFunctionsGrad(Array<OneD, NekDouble> & dGL,
									   Array<OneD, NekDouble> & dGR)
	{	
		LibUtilities::BasisSharedPtr Basis;
		LibUtilities::BasisSharedPtr BasisFR_Left;
        LibUtilities::BasisSharedPtr BasisFR_Right;
		Basis = Domain->GetExp(0)->GetBasis(0);
		
        /// Number of modes
		int k  = Basis->GetNumModes();
        
        /// Total number of quadrature points
		int np = Basis->GetNumPoints();
        
        /// Type of points
        const LibUtilities::PointsKey FRpoints = Basis->GetPointsKey();
        
        /// Computing the derivatives of the Radau correction functions to recover the DG scheme
        if      (GFtype == "DG")
        {
            std::cout << "\n======= Scheme recovered: DG ========" << std::endl;
            const LibUtilities::BasisKey  FRBase_Left (LibUtilities::eDG_DG_Left,  np, FRpoints);
            const LibUtilities::BasisKey  FRBase_Right(LibUtilities::eDG_DG_Right, np, FRpoints);
            
            BasisFR_Left  = LibUtilities::BasisManager()[FRBase_Left];
            BasisFR_Right = LibUtilities::BasisManager()[FRBase_Right];
            
            dGL = BasisFR_Left ->GetBdata();
            dGR = BasisFR_Right->GetBdata();
        }
        
        /// Computing the derivatives of the correction functions to recover the SD scheme
		else if (GFtype == "SD")
        {   
            std::cout << "\n======= Scheme recovered: SD ========" << std::endl;
            const LibUtilities::BasisKey  FRBase_Left (LibUtilities::eDG_SD_Left,  np, FRpoints);
            const LibUtilities::BasisKey  FRBase_Right(LibUtilities::eDG_SD_Right, np, FRpoints);
            
            BasisFR_Left  = LibUtilities::BasisManager()[FRBase_Left];
            BasisFR_Right = LibUtilities::BasisManager()[FRBase_Right];
            
            dGL = BasisFR_Left ->GetBdata();
            dGR = BasisFR_Right->GetBdata();
        }
        
        /// Computing the derivatives of the correction functions to recover the HU scheme
		else if (GFtype == "HU")
        { 
            std::cout << "\n======= Scheme recovered: HU ========" << std::endl;
            const LibUtilities::BasisKey  FRBase_Left (LibUtilities::eDG_HU_Left,  np, FRpoints);
            const LibUtilities::BasisKey  FRBase_Right(LibUtilities::eDG_HU_Right, np, FRpoints);
            
            BasisFR_Left  = LibUtilities::BasisManager()[FRBase_Left];
            BasisFR_Right = LibUtilities::BasisManager()[FRBase_Right];
            
            dGL = BasisFR_Left ->GetBdata();
            dGR = BasisFR_Right->GetBdata();
        }
        
		else                    {ASSERTL0(false,"options for the g functions are DG, SD and HU");}

        /*
        /// Definition of a file for derivatives inspection
        ofstream outfile;
        
        /// Printing a file with the data of the left derivative
        outfile.open("dGL.dat");
		for(int i = 0; i < np; i++)
		{
			outfile << scientific 
			<< setw (17) 
			<< setprecision(10) 
			<< points[i]
			<< "  " 
			<< dGL[i]  
			<< endl;
		}
		outfile << endl << endl;
		outfile.close();
        
        /// Printing a file with the data of the right derivative
        outfile.open("dGR.dat");
		for(int i = 0; i < np; i++)
		{
			outfile << scientific 
			<< setw (17) 
			<< setprecision(10) 
			<< points[i]
			<< "  " 
			<< dGR[i]  
			<< endl;
		}
		outfile << endl << endl;
		outfile.close();
         */
    }

    
    
    //! 6. METHOD 6: InterpToInterface
	void Advection1DFR::InterpToInterface(const Array<OneD, NekDouble> & total,
										  Array<OneD, NekDouble> & interfaceValue)
	{
        /// Getting the basis
		LibUtilities::BasisSharedPtr Basis;
		Basis = Domain->GetExp(0)->GetBasis(0);
		StdRegions::StdSegExp StdSeg(Basis->GetBasisKey());
		
        /// Definition of the coordinates of the standard segment (-1, 1)
		Array<OneD,NekDouble> coordL(3,0.0);
		Array<OneD,NekDouble> coordR(3,0.0);
		Array<OneD,NekDouble> tmp(nq/ne,0.0);
		coordL[0] = -1.0;
		coordR[0] = 1.0;
        
        /// Interpolation of the input variable total at the left and right interfaces of each element
		int cnt = 0;
		for(int i=0; i < ne; i++)
		{
			tmp = total + i*nq/ne;
			interfaceValue[cnt] = StdSeg.PhysEvaluate(coordL,tmp);cnt++;
			interfaceValue[cnt] = StdSeg.PhysEvaluate(coordR,tmp);cnt++;
		}
	}

    
    
    //! 7. OUTPUT METHODS
    /// a) AppendOutput 
	void Advection1DFR::AppendOutput(const Array<OneD, NekDouble> & approx,
									 const Array<OneD, NekDouble> & exact) const
	{
		ofstream outfile;
		outfile.open("Solution.dat");
		for(int i = 0; i < nq; i++)
		{
			outfile << scientific 
			<< setw (17) 
			<< setprecision(16) 
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
    
    /// b) GenerateGnuplotScript
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
    
    /// b) GenerateMatlabScript
	void Advection1DFR::GenerateMatlabScript() const
	{
		ofstream outfile;
		outfile.open("Matlab_FR_Solution.m");
        
		outfile << "% ========= Matlab Script File =========" << endl;
		outfile << "clc; clear all; close all;"               << endl;                       
		outfile << "solution    = load('Solution.dat');"      << endl;                           
		outfile << "x           = solution(1:end,1);"         << endl;                    
		outfile << "approx      = solution(1:end,2);"         << endl;                        
		outfile << "exact       = solution(1:end,3);"         << endl;
		outfile << "figure()"                                 << endl;
		outfile << "plot(x,approx,'--or')"                      << endl;
		outfile << "hold on"                                  << endl;
		outfile << "plot(x,exact,'-o')"                      << endl;
		outfile << "legend('Approximated','Exact')"           << endl;
		outfile << "hold off"                                 << endl;

		outfile.close();
	}
    
    /// c) SolutionPrint
    void Advection1DFR::SolutionPrint()
	{
		/// Evaluating the error
		NekDouble L2, Linf;
        ExSol->Evaluate(x,y,z,Time,Exac[0]);
		
		L2   = Domain->L2(Exac[0]);
		Linf = Domain->Linf(Exac[0]);
		
		/// Writing solution to file
		AppendOutput(Domain->GetPhys(),Exac[0]);
		GenerateGnuplotScript();
        GenerateMatlabScript();
        
		/// Printing summary of solution details
		const SpatialDomains::ExpansionMap &expansions = graph1D->GetExpansions(); 
		LibUtilities::BasisKey bkey0 = expansions.begin()->second->m_basisKeyVector[0];
		
		cout << endl;
		cout << "Solving Unsteady 1D Advection problem with Flux Recostruction scheme"  << endl;
		
		cout << endl;
        cout << "Number of elements          : " << Domain->GetExpSize()                            << endl;
        cout << "Number of quadrature points : " << Domain->GetTotPoints()                          << endl;
        cout << "Number of interface points  : " << graph1D->GetNvertices()                         << endl;

        cout << endl;
		cout << "Expansion                   : " << LibUtilities::BasisTypeMap[bkey0.GetBasisType()]<< endl;
		cout << "No. modes                   : " << bkey0.GetNumModes()                             << endl;
		cout << "Advection                   : " << AdveFunc->GetExpression()                       << endl;
	
        cout << endl;
		cout << "Time Step                   : " << TimeStep                    << endl;
		cout << "Number of Steps             : " << NumSteps                    << endl;
		cout << "Initial Condition           : " << InitCond->GetExpression()   << endl;

		cout << endl;
		cout << "Exact Solution              : " << ExSol->GetExpression()      << endl;
		cout << "Linf error                  : " << setprecision(16) << L2      << endl;
		cout << "L2   error                  : " << setprecision(16) << Linf    << endl;
	}

    
    
    //! 8. GETTER and AUXILIARY METHODS    
    /// a) GetDomain 
    MultiRegions::DisContField1DSharedPtr Advection1DFR::GetDomain() const
	{
		return Domain;
	}
    
    /// b) GetTime 
	NekDouble Advection1DFR::GetTime() const
	{
		return Time;
	}
  
    /// c) UpdateTime
	void Advection1DFR::UpdateTime()
	{
		Time = Time + TimeStep;
	}
    
    /// d) GetTimeStep
	NekDouble Advection1DFR::GetTimeStep() const
	{
		return TimeStep;
	}
    
    /// e) GetNumSteps
	int Advection1DFR::GetNumSteps() const
	{
		return NumSteps;
	}
    
    /// f) GetNumPoints
	int Advection1DFR::GetNumPoints() const
	{
		return nq;
	}
}