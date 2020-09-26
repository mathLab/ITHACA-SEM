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

using namespace std;

namespace Nektar
{
    string EigenValuesAdvection::className = GetEquationSystemFactory().RegisterCreatorFunction("EigenValuesAdvection", EigenValuesAdvection::create, "Eigenvalues of the weak advection operator.");

    EigenValuesAdvection::EigenValuesAdvection(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : EquationSystem(pSession, pGraph)
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

        GetFunction( "AdvectionVelocity")->Evaluate(vel,  m_velocity);

        // Type of advection class to be used
        switch(m_projectionType)
        {
            // Discontinuous field
            case MultiRegions::eDiscontinuous:
            {
                // Define the normal velocity fields
                if (m_fields[0]->GetTrace())
                {
                    m_traceVn = Array<OneD, NekDouble>(GetTraceNpoints());
                }

                string advName;
                string riemName;
                m_session->LoadSolverInfo(
                    "AdvectionType", advName, "WeakDG");
                m_advObject = SolverUtils::
                    GetAdvectionFactory().CreateInstance(advName, advName);
                m_advObject->SetFluxVector(
                    &EigenValuesAdvection::GetFluxVector, this);
                m_session->LoadSolverInfo(
                    "UpwindType", riemName, "Upwind");
                m_riemannSolver = SolverUtils::
                    GetRiemannSolverFactory().CreateInstance(
                        riemName, m_session);
                m_riemannSolver->SetScalar(
                    "Vn", &EigenValuesAdvection::GetNormalVelocity, this);

                m_advObject->SetRiemannSolver(m_riemannSolver);
                m_advObject->InitObject(m_session, m_fields);
                break;
            }
            // Continuous field
            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                string advName;
                m_session->LoadSolverInfo(
                    "AdvectionType", advName, "NonConservative");
                m_advObject = SolverUtils::
                    GetAdvectionFactory().CreateInstance(advName, advName);
                m_advObject->SetFluxVector(
                    &EigenValuesAdvection::GetFluxVector, this);
                break;
            }
            default:
            {
                ASSERTL0(false, "Unsupported projection type.");
                break;
            }
        }
    }

    /**
     * @brief Get the normal velocity
     */
    Array<OneD, NekDouble> &EigenValuesAdvection::GetNormalVelocity()
    {
        // Number of trace (interface) points
        int nTracePts = GetTraceNpoints();

        // Auxiliary variable to compute the normal velocity
        Array<OneD, NekDouble> tmp(nTracePts);

        // Reset the normal velocity
        Vmath::Zero(nTracePts, m_traceVn, 1);

        for (int i = 0; i < m_velocity.size(); ++i)
        {
            m_fields[0]->ExtractTracePhys(m_velocity[i], tmp);

            Vmath::Vvtvp(nTracePts,
                         m_traceNormals[i], 1,
                         tmp,               1,
                         m_traceVn,         1,
                         m_traceVn,         1);
        }

        return m_traceVn;
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
        int dofs = GetNcoeffs();
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
        m_advObject->Advect(nvariables, m_fields, m_velocity, inarray,
                            outarray, 0.0);
        Vmath::Neg(npoints,outarray[0],1);
        switch (m_projectionType)
        {
        case MultiRegions::eDiscontinuous:
            {
                break;
            }
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
            {
                for(int i = 0; i < nvariables; ++i)
                {
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

    void EigenValuesAdvection::GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        ASSERTL1(flux[0].size() == m_velocity.size(),
                 "Dimension of flux array and velocity array do not match");

        int nq = physfield[0].size();

        for (int i = 0; i < flux.size(); ++i)
        {
            for (int j = 0; j < flux[0].size(); ++j)
            {
                Vmath::Vmul(nq, physfield[i], 1, m_velocity[j], 1,
                            flux[i][j], 1);
            }
        }
    }
}
