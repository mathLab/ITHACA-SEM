///////////////////////////////////////////////////////////////////////////////
//
// File FHN.cpp
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
// Description: FHN class definition built on
// ADRBase class
//
///////////////////////////////////////////////////////////////////////////////

#include <FitzHugh-Nagumo/FHN.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>

namespace Nektar

{
    /**
     * Basic construnctor
     */
    FHN::FHN(void):
        ADRBase(),
        m_infosteps(100)
    {     
    }
    
    int nocase_cmp(const string & s1, const string& s2);

    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */
    FHN::FHN(string &fileNameString):
        ADRBase(fileNameString,true),
        m_infosteps(10)
    {

        int i;

        // Set up equation type enum using kEquationTypeStr
        std::string typeStr = m_boundaryConditions->GetSolverInfo("EQTYPE");

        for(i = 0; i < (int) eEquationTypeSize; ++i)
        {
            if(nocase_cmp(kEquationTypeStr[i],typeStr) == 0 )
            {
                m_equationType = (EquationType)i; 
                break;
            }
        }        

        // Equation Setups 
  
            m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        
            for(int i = 0; i < m_spacedim; ++i)
            {
                m_velocity[i] = Array<OneD, NekDouble> (GetNpoints());
            }
            
            EvaluateAdvectionVelocity();
            
            if(m_boundaryConditions->CheckForParameter("IO_InfoSteps") == true)
            {
                m_infosteps =  m_boundaryConditions->GetParameter("IO_InfoSteps");
            }
            
            // check that any user defined boundary condition is indeed implemented
            for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
            {	
                // Time Dependent Boundary Condition (if no use defined then this is empty)
                if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() != "")
                {
                    if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() != "TimeDependent")
                    {
                        ASSERTL0(false,"Unknown USERDEFINEDTYPE boundary condition");
                    }
                }
            }
    }

    void FHN::Evaluateepsilon()
    {
        m_epsilon = m_boundaryConditions->GetParameter("epsilon");
    }

    void FHN::Evaluatebeta()
    {
        m_beta = m_boundaryConditions->GetParameter("beta");
    }

    void FHN::ReadTimemarchingwithmass()
    {
        m_Timemarchingwithmass = m_boundaryConditions->GetParameter("TimemarchingWithmass");
    }


    void FHN::SetUSERInitialConditions(NekDouble initialtime)
    {
        int nq = m_fields[0]->GetNpoints();
      
        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);
      
        NekDouble unew, uinit, vinit, f, fd, Tol, fval, gval;

        // get the coordinates (assuming all fields have the same
        // discretisation)
        m_fields[0]->GetCoords(x0,x1,x2);

        // Get the initial constant u and v

        uinit=0.0; unew=100.0; Tol=0.00000001;
        for(int j=0; j<1000; j++)
        {
            f = fvalue(uinit);
            fd = fderiv(uinit);

            unew = uinit - f/fd;

            if(abs(unew-uinit) < Tol)
            {
                break;
            }

            uinit = unew;
        }
        vinit = uinit - (1.0/3.0)*uinit*uinit*uinit;

        fval = uinit - (1.0/3.0)*uinit*uinit*uinit - vinit;
        gval = uinit + m_beta - 0.5*vinit;

        cout <<"uinit = " << uinit << ", vinit = " << vinit << endl;
        cout << "fval = " << fval << ", gval = " << gval << endl;
        cout << " " << endl;

        for(int j = 0; j < nq; j++)
             {
                  if( x0[j]<=3.50 )
                   {
                      (m_fields[0]->UpdatePhys())[j] = 2.0;
                      (m_fields[1]->UpdatePhys())[j] = vinit;
                   }
                  else
                  {
                     (m_fields[0]->UpdatePhys())[j] = uinit;
                     (m_fields[1]->UpdatePhys())[j] = vinit;
                   }
             }

        for(int i = 0 ; i < m_fields.num_elements(); i++)
	{
            m_fields[i]->SetPhysState(true);
            m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
	}


	// dump initial conditions to file
	for(int i = 0; i < m_fields.num_elements(); ++i)
	  {
	    std::string outname = m_sessionName +"_" + m_boundaryConditions->GetVariable(i) + "_initial.chk";
            ofstream outfile(outname.c_str());
	    m_fields[i]->WriteToFile(outfile,eTecplot);
	  }       
    }


        NekDouble FHN::fvalue(const NekDouble u)
            {
                NekDouble value;

                value = u*u*u + 3.0*u + 6.0*m_beta;
                return value;
            }

        NekDouble FHN::fderiv(const NekDouble u)
            {
                NekDouble value;

                value = 3.0*u*u + 3.0;
                return value;
            }


    void FHN::EvaluateAdvectionVelocity()
    {
        int nq = m_fields[0]->GetNpoints();
        
        std::string velStr[3] = {"Vx","Vy","Vz"};

        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);
      
        // get the coordinates (assuming all fields have the same
        // discretisation)
        m_fields[0]->GetCoords(x0,x1,x2);

        for(int i = 0 ; i < m_velocity.num_elements(); i++)
	{
            SpatialDomains::ConstUserDefinedEqnShPtr ifunc = m_boundaryConditions->GetUserDefinedEqn(velStr[i]);
            
            for(int j = 0; j < nq; j++)
	    {
                m_velocity[i][j] = ifunc->Evaluate(x0[j],x1[j],x2[j]);
	    }
	}
    }
    
    void FHN::ODETest_rhs_u(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
			   Array<OneD, Array<OneD, NekDouble> >&outarray, 
			   const NekDouble time)
	
	{
                NekDouble PI = 3.14159265;

                int nvariables = inarray.num_elements();
                int ncoeffs    = inarray[0].num_elements();        
                int npoints = m_fields[0]->GetNpoints();

		Array<OneD,NekDouble> x0(npoints,0.0);
                Array<OneD,NekDouble> x1(npoints,0.0);
                Array<OneD,NekDouble> x2(npoints,0.0);  

                Array<OneD,NekDouble> physfield(npoints);
                Array<OneD,NekDouble> ft(npoints,0.0);

                MassMultiply(inarray[0], outarray[0], -1, m_Timemarchingwithmass);

                //Compute f(u) in physical field
                m_fields[0]->BwdTrans(outarray[0],physfield);
                m_fields[0]->GetCoords(x0,x1,x2);

                for (int i =0; i< npoints; ++i)
                {
                    ft[i] = (2*PI*PI)*exp(-1.0*time)*cos(PI*x0[i])*cos(PI*x1[i]);
                }
                    
                Vmath::Vsub(npoints, ft, 1, physfield, 1, physfield, 1);

                // Back to modal field
                m_fields[0]->FwdTrans(physfield,outarray[0]);

                MassMultiply(outarray[0], outarray[0], 1, m_Timemarchingwithmass);
	}


    void FHN::ODETest_rhs_u2(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
			   Array<OneD, Array<OneD, NekDouble> >&outarray, 
			   const NekDouble time)
	
	{
                NekDouble PI = 3.14159265;

                int nvariables = inarray.num_elements();
                int ncoeffs    = inarray[0].num_elements();        
                int npoints = m_fields[0]->GetNpoints();

		Array<OneD,NekDouble> x0(npoints,0.0);
                Array<OneD,NekDouble> x1(npoints,0.0);
                Array<OneD,NekDouble> x2(npoints,0.0);  

                Array<OneD,NekDouble> physfield(npoints);
                Array<OneD,NekDouble> ft(npoints,0.0);

                NekDouble csx,csy;

                MassMultiply(inarray[0], outarray[0], -1, m_Timemarchingwithmass);

                //Compute f(u) in physical field
                m_fields[0]->BwdTrans(outarray[0],physfield);

                Vmath::Vmul(npoints,physfield, 1, physfield, 1, physfield, 1);

                m_fields[0]->GetCoords(x0,x1,x2);

                for (int i =0; i< npoints; ++i)
                {
                    csx = cos(PI*x0[i]);
                    csy = cos(PI*x1[i]);
                    ft[i] = exp(-2.0*time)*csx*csx*csy*csy + (2*PI*PI-1.0)*exp(-1.0*time)*csx*csy; 
                }
                    
                Vmath::Vsub(npoints, ft, 1, physfield, 1, physfield, 1);

                // Back to modal field
                m_fields[0]->FwdTrans(physfield,outarray[0]);

                MassMultiply(outarray[0], outarray[0], 1, m_Timemarchingwithmass);
	}


    void FHN::ODEFHNtype_v1(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
			   Array<OneD, Array<OneD, NekDouble> >&outarray, 
			   const NekDouble time)
	
	{
                int nvariables = inarray.num_elements();
                int ncoeffs    = inarray[0].num_elements();
                int npoints = m_fields[0]->GetNpoints();

		const NekDouble coeff = 2.0/m_epsilon;

                Array<OneD, NekDouble> physfield(npoints);
                Array<OneD, NekDouble> temp2(npoints,0.0);
                Array<OneD, NekDouble> temp3(npoints,0.0);
                Array<OneD, NekDouble> temp(npoints,0.0);
					
		for (int i = 0; i < nvariables; ++i)
		{  

                    MassMultiply(inarray[i], outarray[i], -1, m_Timemarchingwithmass);
                    m_fields[i]->BwdTrans(outarray[i],physfield);


                    // temp3 = (2/epsilon)*(u*u - u*u*u)
                    Vmath::Vmul(npoints, physfield, 1, physfield, 1, temp, 1);
                    Vmath::Vcopy(npoints, temp, 1, temp2, 1);
                    Vmath::Vmul(npoints, physfield, 1, temp, 1, temp3, 1);
                    Vmath::Vsub(npoints, temp2, 1, temp3, 1, physfield, 1);
		    Vmath::Smul(npoints, coeff, physfield, 1, physfield, 1);

                    m_fields[i]->FwdTrans(physfield,outarray[i]);
                    MassMultiply(outarray[i], outarray[i], 1, m_Timemarchingwithmass);
		}
	}


    void FHN::ODEFHNtype_v2(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
			   Array<OneD, Array<OneD, NekDouble> >&outarray, 
			   const NekDouble time)
	
	{
                int nvariables = inarray.num_elements();
                int ncoeffs    = inarray[0].num_elements();
                int npoints = m_fields[0]->GetNpoints();

		Array<OneD,NekDouble> x0(npoints,0.0);
                Array<OneD,NekDouble> x1(npoints,0.0);
                Array<OneD,NekDouble> x2(npoints,0.0);  

		const NekDouble coeff = 2.0/m_epsilon;
                const NekDouble Tol = 0.0000001;

                Array<OneD, NekDouble> physfield(npoints);
                Array<OneD, NekDouble> dist(npoints,0.0);
                Array<OneD, NekDouble> temp(npoints,0.0);
					
		for (int i = 0; i < nvariables; ++i)
		{  

                    MassMultiply(inarray[i], outarray[i], -1, m_Timemarchingwithmass);
                    m_fields[i]->BwdTrans(outarray[i],physfield);

                    // temp = u - u*u
                    Vmath::Vmul(npoints, physfield, 1, physfield, 1, temp, 1);
                    Vmath::Vsub(npoints, physfield, 1, temp, 1, temp, 1);

                    m_fields[0]->GetCoords(x0,x1,x2);

                    for (int j =0; j< npoints; ++j)
                    {  
                        dist[j] = 1.0/sqrt(x0[j]*x0[j]+x1[j]*x1[j]+Tol);
                    }

                    // v = (2/epsilon)*u + 1/|x|
                    Vmath::Svtvp(npoints, coeff, physfield, 1, dist, 1, physfield, 1);

                    // f(u) = temp*v = u(1-u)*( (2/epsilon)*u + 1/|x| )
                    Vmath::Vmul(npoints, temp, 1, physfield, 1, physfield, 1);

                    m_fields[i]->FwdTrans(physfield,outarray[i]);
                    MassMultiply(outarray[i], outarray[i], 1, m_Timemarchingwithmass);
		}
	}

    void FHN::ODEFHN_Reaction(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
			   Array<OneD, Array<OneD, NekDouble> >&outarray, 
			   const NekDouble time)
	
	{
                NekDouble Rlambda = -100.0;
                NekDouble theta = 0.25;
                NekDouble alpha = 0.16875;
                NekDouble beta = 1.0;

                int nvariables = inarray.num_elements();
                int ncoeffs    = inarray[0].num_elements();
                int npoints = m_fields[0]->GetNpoints();

                Array<OneD, NekDouble> physfield(npoints);
                Array<OneD, NekDouble> temp1(npoints,0.0);
                Array<OneD, NekDouble> temp2(npoints,0.0);
                Array<OneD, NekDouble> temp3(npoints,0.0);

                MassMultiply(inarray[0], outarray[0], -1, m_Timemarchingwithmass);
                m_fields[0]->BwdTrans(outarray[0],physfield);

                // For v: lambda*(q + v*(v-1)*(v-theta) )
                Vmath::Sadd(npoints, -1.0, physfield, 1, temp1, 1);
                Vmath::Sadd(npoints, -1.0*theta, physfield, 1, temp2, 1);
                Vmath::Vmul(npoints, temp1, 1, temp2, 1, temp3, 1);
                Vmath::Vmul(npoints, temp3, 1, physfield, 1, physfield, 1);

                m_fields[0]->FwdTrans(physfield,outarray[0]);
                MassMultiply(outarray[0], outarray[0], 1, m_Timemarchingwithmass);

                Vmath::Vadd(ncoeffs, inarray[1], 1, outarray[0], 1, outarray[0], 1);
		Vmath::Smul(ncoeffs, Rlambda, outarray[0], 1, outarray[0], 1);

                // For q: alpha*v - q
                Vmath::Svtvp(ncoeffs, -1.0*alpha, inarray[0], 1, inarray[1], 1, outarray[1], 1);
                Vmath::Neg(ncoeffs, outarray[1], 1);
	}


    void FHN::ODEFHN_Spiral_Reaction(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
			   Array<OneD, Array<OneD, NekDouble> >&outarray, 
			   const NekDouble time)
	
	{

                NekDouble m_gamma = 0.5;

                int nvariables = inarray.num_elements();
                int ncoeffs    = inarray[0].num_elements();
                int npoints    = m_fields[0]->GetNpoints();

                Array<OneD, NekDouble> physfield(npoints);
                Array<OneD, NekDouble> temp2(npoints,0.0);
                Array<OneD, NekDouble> temp3(npoints,0.0);

                Array<OneD, NekDouble> temp(ncoeffs, 0.0);

                MassMultiply(inarray[0], outarray[0], -1, m_Timemarchingwithmass);
                m_fields[0]->BwdTrans(outarray[0],physfield);

                // For v: (1/m_epsilon)*( u*(1-u*u/3) - q )
                // physfield = u - (1.0/3.0)*u*u*u
                Vmath::Vmul(npoints, physfield, 1, physfield, 1, temp2, 1);
                Vmath::Vmul(npoints, physfield, 1, temp2, 1, temp3, 1);
                Vmath::Svtvp(npoints, (-1.0/3.0), temp3, 1, physfield, 1, physfield, 1);

                m_fields[0]->FwdTrans(physfield,outarray[0]);
                MassMultiply(outarray[0], outarray[0], 1, m_Timemarchingwithmass);

                Vmath::Vsub(ncoeffs, inarray[1], 1, outarray[0], 1, outarray[0], 1);
		Vmath::Smul(ncoeffs, -1.0/m_epsilon, outarray[0], 1, outarray[0], 1);

                // For q: m_epsilon*( v + m_beta - m_gamma*q )
                Vmath::Smul(ncoeffs, -1.0*m_gamma, inarray[1], 1, temp, 1);
                Vmath::Svtvp(ncoeffs, 1.0, inarray[0], 1, temp, 1, outarray[1], 1);
                Vmath::Sadd(ncoeffs, m_beta, outarray[1], 1, outarray[1], 1);
                Vmath::Smul(ncoeffs, m_epsilon, outarray[1], 1, outarray[1], 1);
	}

  
    void FHN::ODEhelmSolve(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
			   Array<OneD, Array<OneD, NekDouble> >&outarray,
			   const NekDouble time, 
                           const NekDouble lambda)
       {
                int nvariables = inarray.num_elements();
                int ncoeffs    = inarray[0].num_elements();

                NekDouble kappa = 1.0/(lambda*m_epsilon);
									
		for (int i = 0; i < nvariables; ++i)
		{

                         MassMultiply(inarray[i], outarray[i], -1, m_Timemarchingwithmass);

			// Multiply rhs[i] with -1.0/gamma/timestep
                         Vmath::Smul(ncoeffs, -1.0*kappa, outarray[i], 1, outarray[i], 1);
			
			// Update coeffs to m_fields
			 m_fields[i]->UpdateCoeffs() = outarray[i];
			
			// Backward Transformation to nodal coefficients
			  m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), m_fields[i]->UpdatePhys());

	    	        // Solve a system of equations with Helmholtz solver
                          m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs(),kappa);
			
			// The solution is Y[i]
			  outarray[i] = m_fields[i]->GetCoeffs();

                         MassMultiply(outarray[i], outarray[i], 1, m_Timemarchingwithmass);

		}
	}
	
    void FHN::ODEFHN_helmSolve(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
			   Array<OneD, Array<OneD, NekDouble> >&outarray,
			   const NekDouble time, 
                           const NekDouble lambda)
       {
                int nvariables = inarray.num_elements();
                int ncoeffs    = inarray[0].num_elements();

                NekDouble kappa = 1.0/(lambda*m_epsilon);

                // For v: ==============================

                   MassMultiply(inarray[0], outarray[0], -1, m_Timemarchingwithmass);
				
		// Multiply rhs[i] with -1.0/gamma/timestep
                   Vmath::Smul(ncoeffs, -1.0*kappa, outarray[0], 1, outarray[0], 1);
			
		// Update coeffs to m_fields
		   m_fields[0]->UpdateCoeffs() = outarray[0];
			
		// Backward Transformation to nodal coefficients
		   m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(), m_fields[0]->UpdatePhys());

	    	// Solve a system of equations with Helmholtz solver
                   m_fields[0]->HelmSolve(m_fields[0]->GetPhys(),m_fields[0]->UpdateCoeffs(),kappa);
			
		// The solution is Y[i]
		   outarray[0] = m_fields[0]->GetCoeffs();	  
             						  
                   MassMultiply(outarray[0], outarray[0], 1, m_Timemarchingwithmass);

                // For q: No helmholtz solver is needed=============================
                   Vmath::Vcopy(ncoeffs, inarray[1], 1, outarray[1], 1);
	}
	

    void FHN::ODEFHN_Spiral_helmSolve(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
			   Array<OneD, Array<OneD, NekDouble> >&outarray,
			   const NekDouble time, 
                           const NekDouble lambda)
       {
                int nvariables = inarray.num_elements();
                int ncoeffs    = inarray[0].num_elements();

                NekDouble kappa = 1.0/lambda;

                // For v: ==============================

                   MassMultiply(inarray[0], outarray[0], -1, m_Timemarchingwithmass);
				
		// Multiply rhs[i] with -1.0/gamma/timestep
                   Vmath::Smul(ncoeffs, -1.0*kappa, outarray[0], 1, outarray[0], 1);
			
		// Update coeffs to m_fields
		   m_fields[0]->UpdateCoeffs() = outarray[0];
			
		// Backward Transformation to nodal coefficients
		   m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(), m_fields[0]->UpdatePhys());

	    	// Solve a system of equations with Helmholtz solver
                   m_fields[0]->HelmSolve(m_fields[0]->GetPhys(), m_fields[0]->UpdateCoeffs(), kappa);
			
		// The solution is Y[i]
		   outarray[0] = m_fields[0]->GetCoeffs();	  
             						  
                   MassMultiply(outarray[0], outarray[0], 1, m_Timemarchingwithmass);

                // For q: No helmholtz solver is needed=============================
                   Vmath::Vcopy(ncoeffs, inarray[1], 1, outarray[1], 1);
	}


    void FHN::MassMultiply(const Array<OneD, NekDouble> &inarray, 
                                 Array<OneD, NekDouble> &outarray, 
                           const int direction, const int turnon )
    {
        int ncoeffs = inarray.num_elements();
    
          if(turnon)
             {
                 if(direction == -1)
                 {
                   m_fields[0]->MultiplyByInvMassMatrix(inarray,outarray,false);
                 }
                 else if(direction == 1)
                 {
		   MultiRegions::GlobalLinSysKey key(StdRegions::eMass);
                   m_fields[0]->MultiRegions::ExpList::GeneralMatrixOp(key,inarray,outarray);
                 }
             }

          else
          {
              Vmath::Vcopy(ncoeffs, inarray, 1, outarray, 1);
          }
    }



    void FHN::GeneralTimeIntegration(int nsteps, 
	                              LibUtilities::TimeIntegrationMethod IntMethod,
				      LibUtilities::TimeIntegrationSchemeOperators ode)
    {
        int i,n,nchk = 0;
        int ncoeffs = m_fields[0]->GetNcoeffs();
        int nvariables = m_fields.num_elements();

        // Set up wrapper to fields data storage. 
        Array<OneD, Array<OneD, NekDouble> >   fields(nvariables);
        Array<OneD, Array<OneD, NekDouble> >   tmp(nvariables);
        
        for(i = 0; i < nvariables; ++i)
        {
            fields[i]  = m_fields[i]->UpdateCoeffs();
        }

        if( (m_projectionType==eGalerkin) && (m_Timemarchingwithmass==1) )
        {
            // calculate the variable u* = Mu
            // we are going to TimeIntegrate this new variable u*
            MultiRegions::GlobalLinSysKey key(StdRegions::eMass);
            for(int i = 0; i < nvariables; ++i)
            {
                tmp[i] = Array<OneD, NekDouble>(ncoeffs);
                m_fields[i]->MultiRegions::ExpList::GeneralMatrixOp(key,fields[i],fields[i]);
            }
        }

        // Declare an array of TimeIntegrationSchemes
        // For multi-stage methods, this array will have just one entry containing
        // the actual multi-stage method...
        // For multi-steps method, this can have multiple entries
        //  - the first scheme will used for the first timestep (this is an initialization scheme)
        //  - the second scheme will used for the first timestep (this is an initialization scheme)
        //  - ...
        //  - the last scheme will be used for all other time-steps (this will be the actual scheme)
        Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> IntScheme;
        LibUtilities::TimeIntegrationSolutionSharedPtr u;
        int numMultiSteps;

        switch(IntMethod)
        {
	case LibUtilities::eIMEXdirk_3_4_3:
	case LibUtilities::eDIRKOrder3:
        case LibUtilities::eBackwardEuler:      
        case LibUtilities::eForwardEuler:      
        case LibUtilities::eClassicalRungeKutta4:
            {
                numMultiSteps = 1;

                IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

                LibUtilities::TimeIntegrationSchemeKey IntKey(IntMethod);
                IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey];

                u = IntScheme[0]->InitializeScheme(m_timestep,fields,m_time,ode);
            }
            break;
        default:
            {
                ASSERTL0(false,"populate switch statement for integration scheme");
            }
        }
					          
        for(n = 0; n < nsteps; ++n)
        {
            //----------------------------------------------
            // Perform time step integration
            //----------------------------------------------
            if( n < numMultiSteps-1)
            {
                // Use initialisation schemes
                fields = IntScheme[n]->TimeIntegrate(m_timestep,u,ode);
            }
            else
            {
                fields = IntScheme[numMultiSteps-1]->TimeIntegrate(m_timestep,u,ode);
            }

            m_time += m_timestep;

         if( (m_projectionType==eGalerkin) && (m_Timemarchingwithmass==1) )
            {
                // Project the solution u* onto the boundary conditions to
                // obtain the actual solution
                SetBoundaryConditions(m_time);
                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->MultiplyByInvMassMatrix(fields[i],tmp[i],false);
                    fields[i] = tmp[i];	   		    
                }
            }

            //----------------------------------------------
            // Dump analyser information
            //----------------------------------------------
            if(!((n+1)%m_infosteps))
            {
                     for(i = 0; i < nvariables; ++i)
                	{
		          (m_fields[i]->UpdateCoeffs()) = fields[i];
                          m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),m_fields[i]->UpdatePhys());
                        }

                cout << "Steps: " << n+1 << "\t Time: " << m_time 
                     << "\t L2err: " << this->L2Error(0) << endl;
            }
            
            if(n&&(!((n+1)%m_checksteps)))
            {
                     for(i = 0; i < nvariables; ++i)
                	{
		             (m_fields[i]->UpdateCoeffs()) = fields[i];
                	}
                  Checkpoint_Output(nchk++);
            }
        }
        
        for(i = 0; i < nvariables; ++i)
        {
	  (m_fields[i]->UpdateCoeffs()) = fields[i];
        }
    }
	
    
  //----------------------------------------------------
  void FHN::SetBoundaryConditions(NekDouble time)
  {
    int nvariables = m_fields.num_elements();
    for (int i = 0; i < nvariables; ++i)
    {
        m_fields[i]->EvaluateBoundaryConditions(time);
    }
    
//     // loop over Boundary Regions
//     for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
//       {	
	
// 	// Time Dependent Boundary Condition
// 	if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() == "TimeDependent")
// 	  {
// 	    for (int i = 0; i < nvariables; ++i)
// 	      {
// 		m_fields[i]->EvaluateBoundaryConditions(time);
// 	      }
// 	  }
//       }
  }
  
 

    void FHN::Summary(std::ostream &out)
    {   
        cout << "=======================================================================" << endl;
        cout << "\tEquation Type   : "<< kEquationTypeStr[m_equationType] << endl;
        out << "\tSession Name    : " << m_sessionName << endl;
	out << "\tExp. Dimension  : " << m_expdim << endl;
        out << "\tMax Exp. Order  : " << m_fields[0]->EvalBasisNumModesMax() << endl;
        if(m_projectionType == eGalerkin)
        {
            out << "\tProjection Type : Galerkin" <<endl;
        }
        else
        {
            out << "\tProjection Type : Discontinuous Galerkin" <<endl;
        }
        out << "\tTime Step       : " << m_timestep << endl;
        out << "\tNo. of Steps    : " << m_steps << endl;
        out << "\tCheckpoints     : " << m_checksteps <<" steps" <<endl;
        out << "\tepsilon         : " << m_epsilon << endl;
        out << "\tbeta            : " << m_beta << endl;
      cout << "=======================================================================" << endl;

    }

    
} //end of namespace

/**
* $Log: FHN.cpp,v $
* Revision 1.3  2009/03/16 14:41:22  sehunchun
* FHN model update
*
* Revision 1.2  2009/03/07 21:18:00  sehunchun
* FHN updated
*
* Revision 1.1  2009/03/06 16:02:55  sehunchun
* FitzHugh-Nagumo modeling
*
*
* Revision 1.1  2008/08/22 09:48:23  pvos
* Added Sehun' FHN solver
*
**/
