///////////////////////////////////////////////////////////////////////////////
//
// File EigenValuesAdvection.cpp
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
// Description: 
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <ADRSolver/EquationSystems/EigenValuesAdvection.h>

namespace Nektar
{
    string EigenValuesAdvection::className = GetEquationSystemFactory().RegisterCreatorFunction("EigenValuesAdvection", EigenValuesAdvection::create, "Eigenvalues of the weak advection operator.");

    EigenValuesAdvection::EigenValuesAdvection(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : EquationSystem(pSession)
    {
    }

    void EigenValuesAdvection::v_InitObject()
    {
        EquationSystem::v_InitObject();

        // Define Velocity fields
        m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        std::vector<std::string> vel;
        vel.push_back("Vx");
        vel.push_back("Vy");
        vel.push_back("Vz");
        vel.resize(m_spacedim);

        EvaluateFunction(vel, m_velocity, "AdvectionVelocity");
    }
	
	void EigenValuesAdvection::v_DoInitialise()
    {
       
    }

    EigenValuesAdvection::~EigenValuesAdvection()
    {

    }

    void EigenValuesAdvection::v_DoSolve()
    {
        int nvariables = 1;
		int i,dofs;
		//bool UseContCoeffs = false;
		
		Array<OneD, Array<OneD, NekDouble> > inarray(nvariables);
		Array<OneD, Array<OneD, NekDouble> > tmp(nvariables);
		Array<OneD, Array<OneD, NekDouble> > outarray(nvariables);
		Array<OneD, Array<OneD, NekDouble> > WeakAdv(nvariables);
		
		int npoints = GetNpoints();
		int ncoeffs = GetNcoeffs();
		
		switch (m_projectionType)
		{
                case MultiRegions::eDiscontinuous:
                    {
                        dofs = ncoeffs;
                        break;
                    }
                case MultiRegions::eGalerkin:
                case MultiRegions::eMixed_CG_Discontinuous:
                    {
                        //dofs = GetContNcoeffs();
                        //UseContCoeffs = true;
                        break;
                    }
		}
		
		cout << endl;
		cout << "Num Phys Points = " << npoints << endl; // phisical points
		cout << "Num Coeffs      = " << ncoeffs << endl; //
		cout << "Num Cont Coeffs = " << dofs << endl;
		
		inarray[0]  = Array<OneD, NekDouble>(npoints,0.0);
		outarray[0] = Array<OneD, NekDouble>(npoints,0.0);
		tmp[0] = Array<OneD, NekDouble>(npoints,0.0);
		
		WeakAdv[0]  = Array<OneD, NekDouble>(ncoeffs,0.0);
		Array<OneD, NekDouble> MATRIX(npoints*npoints,0.0);
		
		for (int j = 0; j < npoints; j++)
		{
		
		inarray[0][j] = 1.0;
       
	    /// Feeding the weak Advection oprator with  a vector (inarray)
        /// Looping on inarray and changing the position of the only non-zero entry
		/// we simulate the multiplication by the identity matrix.
		/// The results stored in outarray is one of the columns of the weak advection oprators
		/// which are then stored in MATRIX for the futher eigenvalues calculation.

        switch (m_projectionType)
        {
        case MultiRegions::eDiscontinuous:
            {
                WeakDGAdvection(inarray, WeakAdv,true,true,1);
                
                m_fields[0]->MultiplyByElmtInvMass(WeakAdv[0],WeakAdv[0]);
		
                m_fields[0]->BwdTrans(WeakAdv[0],outarray[0]);
                
                Vmath::Neg(npoints,outarray[0],1);
                break;
            }
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
            {
                // Calculate -V\cdot Grad(u);
                for(i = 0; i < nvariables; ++i)
                {
                    //Projection
                    m_fields[i]->FwdTrans(inarray[i],WeakAdv[i]);
                    
                    m_fields[i]->BwdTrans_IterPerExp(WeakAdv[i],tmp[i]);
                    
                    //Advection operator
                    AdvectionNonConservativeForm(m_velocity,tmp[i],outarray[i]);
                    
                    Vmath::Neg(npoints,outarray[i],1);
                    
                    //m_fields[i]->MultiplyByInvMassMatrix(WeakAdv[i],WeakAdv[i]);
                    //Projection
                    m_fields[i]->FwdTrans(outarray[i],WeakAdv[i]);
                    
                    m_fields[i]->BwdTrans_IterPerExp(WeakAdv[i],outarray[i]);
                }
                break;
            }
        }
	
        /// The result is stored in outarray (is the j-th columns of the weak advection operator).
        /// We now store it in MATRIX(j)
        Vmath::Vcopy(npoints,&(outarray[0][0]),1,&(MATRIX[j]),npoints);
	
        /// Set the j-th entry of inarray back to zero
        inarray[0][j] = 0.0;
		}
                
		////////////////////////////////////////////////////////////////////////////////
		/// Calulating the eigenvalues of the weak advection operator stored in (MATRIX)
		/// using Lapack routines
		
		char jobvl = 'N';
		char jobvr = 'N';
		int info = 0, lwork = 3*npoints;
		NekDouble dum;
		
		Array<OneD, NekDouble> EIG_R(npoints);
		Array<OneD, NekDouble> EIG_I(npoints);
		
		Array<OneD, NekDouble> work(lwork);
		
		Lapack::Dgeev(jobvl,jobvr,npoints,MATRIX.get(),npoints,EIG_R.get(),EIG_I.get(),&dum,1,&dum,1,&work[0],lwork,info);
		
		////////////////////////////////////////////////////////
		//Print Matrix
		FILE *mFile;
		
		mFile = fopen ("WeakAdvMatrix.txt","w");
		for(int j = 0; j<npoints; j++)
		{
			for(int k = 0; k<npoints; k++)
			{
				fprintf(mFile,"%e ",MATRIX[j*npoints+k]);
			}
			fprintf(mFile,"\n");
		}
		fclose (mFile);
		
		////////////////////////////////////////////////////////
		//Output of the EigenValues
		FILE *pFile;
		
		pFile = fopen ("Eigenvalues.txt","w");
		for(int j = 0; j<npoints; j++)
		{
			fprintf(pFile,"%e %e\n",EIG_R[j],EIG_I[j]);
		}
		fclose (pFile);
		
		cout << "\nEigenvalues : " << endl;
		for(int j = 0; j<npoints; j++)
		{
			cout << EIG_R[j] << "\t" << EIG_I[j] << endl;
		}
		cout << endl;
    }

    void EigenValuesAdvection::v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield,
											Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        ASSERTL1(flux.num_elements() == m_velocity.num_elements(),"Dimension of flux array and velocity array do not match");
		
        for(int j = 0; j < flux.num_elements(); ++j)
        {
            Vmath::Vmul(GetNpoints(),physfield[i],1,
						m_velocity[j],1,flux[j],1);
        }
    }
	
    void EigenValuesAdvection::v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        int i;
		
        int nTraceNumPoints = GetTraceNpoints();
        int nvel = m_spacedim; //m_velocity.num_elements();
		
        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
        Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);		
        
        //Get Edge Velocity - Could be stored if time independent
        for(i = 0; i < nvel; ++i)
        {
            m_fields[0]->ExtractTracePhys(m_velocity[i], Fwd);
            Vmath::Vvtvp(nTraceNumPoints,m_traceNormals[i],1,Fwd,1,Vn,1,Vn,1);
        }
        
        for(i = 0; i < numflux.num_elements(); ++i)
        {
            m_fields[i]->GetFwdBwdTracePhys(physfield[i],Fwd,Bwd);
            m_fields[i]->GetTrace()->Upwind(Vn,Fwd,Bwd,numflux[i]);
            // calculate m_fields[i]*Vn
            Vmath::Vmul(nTraceNumPoints,numflux[i],1,Vn,1,numflux[i],1);
        }
    }
}
