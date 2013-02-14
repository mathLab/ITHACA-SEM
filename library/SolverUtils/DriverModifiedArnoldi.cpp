///////////////////////////////////////////////////////////////////////////////
//
// File DriverModifiedArnoldi.cpp
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
// Description: Driver for eigenvalue analysis using the modified Arnoldi
//              method.
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <SolverUtils/DriverModifiedArnoldi.h>

namespace Nektar
{
    namespace SolverUtils
    {
        string DriverModifiedArnoldi::className = GetDriverFactory().RegisterCreatorFunction("ModifiedArnoldi", DriverModifiedArnoldi::create);
        string DriverModifiedArnoldi::driverLookupId = LibUtilities::SessionReader::RegisterEnumValue("Driver","ModifiedArnoldi",0);
    
        /**
         *
         */
        DriverModifiedArnoldi::DriverModifiedArnoldi(
            const LibUtilities::SessionReaderSharedPtr        pSession)
            : DriverArnoldi(pSession)
        {
        }
    
    
        /**
         *
         */
        DriverModifiedArnoldi::~DriverModifiedArnoldi()
        {
        }
    
    
        /**
         *
         */
        void DriverModifiedArnoldi::v_InitObject(ostream &out)
        {
            DriverArnoldi::v_InitObject(out);
        
            m_equ[0]->PrintSummary(out);
            
            // Print session parameters
            out << "\tArnoldi solver type   : Modified Arnoldi" << endl;

            DriverArnoldi::ArnoldiSummary(out);

            m_equ[m_nequ - 1]->DoInitialise();
	
            //FwdTrans Initial conditions to be in Coefficient Space
            m_equ[m_nequ-1] ->TransPhysToCoeff();

        }
    
	
        /**
         *
         */
        void DriverModifiedArnoldi::v_Execute(ostream &out)
        {
            int i                   = 0;
            int j                   = 0;
            int nq                  = m_equ[0]->UpdateFields()[0]->GetNcoeffs();
            int ntot                = m_nfields*nq;
            int converged           = 0;
            NekDouble resnorm       = 0.0;
            std::string evlFile     = m_session->GetFilename().substr(0,m_session->GetFilename().find_last_of('.')) + ".evl";
            ofstream evlout(evlFile.c_str());
        
            // Allocate memory
            Array<OneD, NekDouble> alpha     = Array<OneD, NekDouble> (m_kdim+1,      0.0);
            Array<OneD, NekDouble> wr        = Array<OneD, NekDouble> (m_kdim,        0.0);
            Array<OneD, NekDouble> wi        = Array<OneD, NekDouble> (m_kdim,        0.0);
            Array<OneD, NekDouble> zvec      = Array<OneD, NekDouble> (m_kdim*m_kdim, 0.0);
        
            Array<OneD, Array<OneD, NekDouble> > Kseq
                = Array<OneD, Array<OneD, NekDouble> > (m_kdim + 1);
            Array<OneD, Array<OneD, NekDouble> > Tseq
                = Array<OneD, Array<OneD, NekDouble> > (m_kdim + 1);
            for (i = 0; i < m_kdim + 1; ++i)
            {
                Kseq[i] = Array<OneD, NekDouble>(ntot, 0.0);
                Tseq[i] = Array<OneD, NekDouble>(ntot, 0.0);
            }

		
            // Copy starting vector into second sequence element (temporary).
            if(m_session->DefinesFunction("InitialConditions"))
            {
                out << "\tInital vector       : specified in input file " << endl;
                m_equ[0]->SetInitialConditions(0.0,false);

                CopyFieldToArnoldiArray(Kseq[1]);
            }
            else 
            {
                out << "\tInital vector       : random  " << endl;
             
                NekDouble eps=1;
             
                Vmath::FillWhiteNoise(ntot, eps , &Kseq[1][0], 1);
             
            }
         
            // Perform one iteration to enforce boundary conditions.
            // Set this as the initial value in the sequence.
            EV_update(Kseq[1], Kseq[0]);
            out << "Iteration: " << 0 <<  endl;
         
            // Normalise first vector in sequence
            alpha[0] = std::sqrt(Vmath::Dot(ntot, &Kseq[0][0], 1, &Kseq[0][0], 1));
            //alpha[0] = std::sqrt(alpha[0]);
            Vmath::Smul(ntot, 1.0/alpha[0], Kseq[0], 1, Kseq[0], 1);

            // Fill initial krylov sequence
            NekDouble resid0;
            for (i = 1; !converged && i <= m_kdim; ++i)
            {
                // Compute next vector
                EV_update(Kseq[i-1], Kseq[i]);
             
                // Normalise
                alpha[i] = std::sqrt(Vmath::Dot(ntot, &Kseq[i][0], 1, &Kseq[i][0], 1));
                //alpha[i] = std::sqrt(alpha[i]);
                Vmath::Smul(ntot, 1.0/alpha[i], Kseq[i], 1, Kseq[i], 1);
             
                // Copy Krylov sequence into temporary storage
                for (int k = 0; k < i + 1; ++k)
                {
                    Vmath::Vcopy(ntot, Kseq[k], 1, Tseq[k], 1);
                }
             
                // Generate Hessenberg matrix and compute eigenvalues of it.
                EV_small(Tseq, ntot, alpha, i, zvec, wr, wi, resnorm);
             
                // Test for convergence.
                converged = EV_test(i,i,zvec,wr,wi,resnorm,std::min(i,m_nvec),evlout,resid0); 
                converged = max (converged, 0);
                out << "Iteration: " <<  i << " (residual : " << resid0 << ")" <<endl;
            }
         
            // Continue with full sequence
            if (!converged)
            {
                for (i = m_kdim + 1; !converged && i <= m_nits; ++i)
                {
                    // Shift all the vectors in the sequence.
                    // First vector is removed.
                    //NekDouble invnorm = 1.0/sqrt(Blas::Ddot(ntot,Kseq[1],1,Kseq[1],1));
                    for (int j = 1; j <= m_kdim; ++j)
                    {
                        alpha[j-1] = alpha[j];
                        //Vmath::Smul(ntot,invnorm,Kseq[j],1,Kseq[j],1);
                        Vmath::Vcopy(ntot, Kseq[j], 1, Kseq[j-1], 1);
                    }
                 
                    // Compute next vector
                    EV_update(Kseq[m_kdim - 1], Kseq[m_kdim]);
                 
                    // Compute new scale factor
                    alpha[m_kdim] = std::sqrt(Vmath::Dot(ntot, &Kseq[m_kdim][0], 1, &Kseq[m_kdim][0], 1));
                    //alpha[m_kdim] = std::sqrt(alpha[m_kdim]);
                    Vmath::Smul(ntot, 1.0/alpha[m_kdim], Kseq[m_kdim], 1, Kseq[m_kdim], 1);
                 
                    // Copy Krylov sequence into temporary storage
                    for (int k = 0; k < m_kdim + 1; ++k)
                    {
                        Vmath::Vcopy(ntot, Kseq[k], 1, Tseq[k], 1);
                    }
                 
                    // Generate Hessenberg matrix and compute eigenvalues of it
                    EV_small(Tseq, ntot, alpha, m_kdim, zvec, wr, wi, resnorm);
                 
                    // Test for convergence.
                    converged = EV_test(i,m_kdim,zvec,wr,wi,resnorm,m_nvec,evlout,resid0);
                    out << "Iteration: " <<  i << " (residual : " << resid0 << ")" <<endl;
                }
            }
         
            m_equ[0]->Output();
         
            // Evaluate and output computation time and solution accuracy.
            // The specific format of the error output is essential for the
            // regression tests to work.
            // Evaluate L2 Error
            for(j = 0; j < m_equ[0]->GetNvariables(); ++j)
            {
                NekDouble vL2Error = m_equ[0]->L2Error(j,false);
                NekDouble vLinfError = m_equ[0]->LinfError(j);
                if (m_comm->GetRank() == 0)
                {
                    out << "L 2 error (variable " << m_equ[0]->GetVariable(j) << ") : " << vL2Error << endl;
                    out << "L inf error (variable " << m_equ[0]->GetVariable(j) << ") : " << vLinfError << endl;
                }
            }
         
            // Process eigenvectors and write out.
            EV_post(Tseq, Kseq, ntot, min(--i, m_kdim), m_nvec, zvec, wr, wi, converged);

            // store eigenvalues so they can be access from driver class
            m_real_evl = wr;
            m_imag_evl = wi;
         
            // Close the runtime info file.
            evlout.close();
        }
    
    
        /**
         *
         */
        void DriverModifiedArnoldi::EV_update(
            Array<OneD, NekDouble> &src,
            Array<OneD, NekDouble> &tgt)
        {
            // Copy starting vector into first sequence element.
            CopyArnoldiArrayToField(src);
		
            m_equ[0]->TransCoeffToPhys();

		
            m_equ[0]->DoSolve();

            if(m_EvolutionOperator == eTransientGrowth)
            {
                Array<OneD, MultiRegions::ExpListSharedPtr> fields;
                fields = m_equ[0]->UpdateFields();
		
                //start Adjoint with latest fields of direct 
                CopyFwdToAdj();

                m_equ[1]->TransCoeffToPhys();
			
                m_equ[1]->DoSolve();
            }

            // Copy starting vector into first sequence element.
            CopyFieldToArnoldiArray(tgt);
        }
	
	
        /**
         *
         */
        void DriverModifiedArnoldi::EV_small(
            Array<OneD, Array<OneD, NekDouble> > &Kseq,
            const int ntot,
            const Array<OneD, NekDouble> &alpha,
            const int kdim,
            Array<OneD, NekDouble> &zvec,
            Array<OneD, NekDouble> &wr,
            Array<OneD, NekDouble> &wi,
            NekDouble &resnorm)
        {
            int kdimp = kdim + 1;
            int lwork = 10*kdim;
            int ier;
            Array<OneD, NekDouble> R(kdimp * kdimp, 0.0);
            Array<OneD, NekDouble> H(kdimp * kdim, 0.0);
            Array<OneD, NekDouble> rwork(lwork, 0.0);
	
            // Modified G-S orthonormalisation
            for (int i = 0; i < kdimp; ++i)
            {
                NekDouble gsc = std::sqrt(Vmath::Dot(ntot, &Kseq[i][0], 1, &Kseq[i][0], 1));
                ASSERTL0(gsc != 0.0, "Vectors are linearly independent.");
                R[i*kdimp+i] = gsc;
                Vmath::Smul(ntot, 1.0/gsc, Kseq[i], 1, Kseq[i], 1);
                for (int j = i + 1; j < kdimp; ++j)
                {
                    gsc = Vmath::Dot(ntot, Kseq[i], 1, Kseq[j], 1);
                    Vmath::Svtvp(ntot, -gsc, Kseq[i], 1, Kseq[j], 1, Kseq[j], 1);
                    R[j*kdimp+i] = gsc;
                }
            }
	
            // Compute H matrix
            for (int i = 0; i < kdim; ++i)
            {
                for (int j = 0; j < kdim; ++j)
                {
                    H[j*kdim+i] = alpha[j+1] * R[(j+1)*kdimp+i]
                        - Vmath::Dot(j, &H[0] + i, kdim, &R[0] + j*kdimp, 1);
                    H[j*kdim+i] /= R[j*kdimp+j];
                }
            }
	
            H[(kdim-1)*kdim+kdim] = alpha[kdim]
                * std::fabs(R[kdim*kdimp+kdim] / R[(kdim-1)*kdimp + kdim-1]);
	
            Lapack::dgeev_('N','V',kdim,&H[0],kdim,&wr[0],&wi[0],0,1,&zvec[0],kdim,&rwork[0],lwork,ier);
	
            ASSERTL0(!ier, "Error with dgeev");
	
            resnorm = H[(kdim-1)*kdim + kdim];
        }
    
	
        /**
         *
         */
        int DriverModifiedArnoldi::EV_test(
            const int itrn,
            const int kdim,
            Array<OneD, NekDouble> &zvec,
            Array<OneD, NekDouble> &wr,
            Array<OneD, NekDouble> &wi,
            const NekDouble resnorm,
            const int nvec,
            ofstream &evlout,
            NekDouble &resid0)
        {
            int idone = 0;
            // NekDouble period = 0.1;
	
            Array<OneD, NekDouble> resid(kdim);
            for (int i = 0; i < kdim; ++i)
            {
                resid[i] = resnorm * std::fabs(zvec[kdim - 1 + i*kdim]) /
                    std::sqrt(Vmath::Dot(kdim, &zvec[0] + i*kdim, 1, &zvec[0] + i*kdim, 1));
                if (wi[i] < 0.0) resid[i-1] = resid[i] = hypot(resid[i-1], resid[i]);
            }
            EV_sort(zvec, wr, wi, resid, kdim);
	
            if (resid[nvec-1] < m_evtol) idone = nvec;
	
            evlout << "-- Iteration = " << itrn << ", H(k+1, k) = " << resnorm << endl;
		
            evlout.precision(4);
            evlout.setf(ios::scientific, ios::floatfield);
            if(m_timeSteppingAlgorithm)
            {
                evlout << "EV  Magnitude   Angle       Growth      Frequency   Residual"
                       << endl;
            }
            else
            {
                evlout << "EV  Real        Imaginary   inverse real  inverse imag  Residual"
                       << endl;
            }
        
            for (int i = 0; i < kdim; i++) 
            {
                WriteEvs(evlout,i,wr[i],wi[i],resid[i]);
            }
            
	
            resid0 = resid[nvec-1];
            return idone;
        }
    
    
        /**
         *
         */
        void DriverModifiedArnoldi::EV_sort(
            Array<OneD, NekDouble> &evec,
            Array<OneD, NekDouble> &wr,
            Array<OneD, NekDouble> &wi,
            Array<OneD, NekDouble> &test,
            const int dim)
        {
            Array<OneD, NekDouble> z_tmp(dim,0.0);
            NekDouble wr_tmp, wi_tmp, te_tmp;
            for (int j = 1; j < dim; ++j)
            {
                wr_tmp = wr[j];
                wi_tmp = wi[j];
                te_tmp = test[j];
                Vmath::Vcopy(dim, &evec[0] + j*dim, 1, &z_tmp[0], 1);
                int i = j - 1;
                while (i >= 0 && test[i] > te_tmp)
                {
                    wr[i+1] = wr[i];
                    wi[i+1] = wi[i];
                    test[i+1] = test[i];
                    Vmath::Vcopy(dim, &evec[0] + i*dim, 1, &evec[0] + (i+1)*dim, 1);
                    i--;
                }
                wr[i+1] = wr_tmp;
                wi[i+1] = wi_tmp;
                test[i+1] = te_tmp;
                Vmath::Vcopy(dim, &z_tmp[0], 1, &evec[0] + (i+1)*dim, 1);
            }	
        };
    

        /**
         *
         */
        void DriverModifiedArnoldi::EV_post(
            Array<OneD, Array<OneD, NekDouble> > &Tseq,
            Array<OneD, Array<OneD, NekDouble> > &Kseq,
            const int ntot,
            const int kdim,
            const int nvec,
            Array<OneD, NekDouble> &zvec,
            Array<OneD, NekDouble> &wr,
            Array<OneD, NekDouble> &wi,
            const int icon)
        {
            if (icon == 0)
            {
                // Not converged, write final Krylov vector
                ASSERTL0(false, "Convergence was not achieved within the prescribed number of iterations.");
            }
            else if (icon < 0)
            {
                // Minimum residual reached
                ASSERTL0(false, "Minimum residual reached.");
            }
            else if (icon == nvec)
            {
                // Converged, write out eigenvectors
                EV_big(Tseq, Kseq, ntot, kdim, icon, zvec, wr, wi);
                Array<OneD, MultiRegions::ExpListSharedPtr> fields
                    = m_equ[0]->UpdateFields();

                for (int j = 0; j < icon; ++j)
                {
                    std::string file = m_session->GetFilename().substr(0,m_session->GetFilename().find_last_of('.')) + "_eig_" + boost::lexical_cast<std::string>(j);
                
                    WriteFld(file,Kseq[j]);
                }
            }
            else
            {
                // Not recognised
                ASSERTL0(false, "Unrecognised value.");
            }
        }


        /**
         *
         */
        void DriverModifiedArnoldi::EV_big(
            Array<OneD, Array<OneD, NekDouble> > &bvecs,
            Array<OneD, Array<OneD, NekDouble> > &evecs,
            const int ntot,
            const int kdim,
            const int nvec,
            Array<OneD, NekDouble> &zvec,
            Array<OneD, NekDouble> &wr,
            Array<OneD, NekDouble> &wi)
        {
            NekDouble wgt, norm;

            // Generate big eigenvectors
            for (int j = 0; j < nvec; ++j)
            {
                Vmath::Zero(ntot, evecs[j], 1);
                for (int i = 0; i < kdim; ++i)
                {
                    wgt = zvec[i + j*kdim];
                    Vmath::Svtvp(ntot, wgt, bvecs[i], 1, evecs[j], 1, evecs[j], 1);
                }
            }

            // Normalise the big eigenvectors
            for (int i = 0; i < nvec; ++i)
            {
                if (wi[i] == 0.0)   // Real mode
                {
                    norm = std::sqrt(Vmath::Dot(ntot, evecs[i], evecs[i]));
                    Vmath::Smul(ntot, 1.0/norm, evecs[i], 1, evecs[i], 1);
                }
                else
                {
                    norm = Vmath::Dot(ntot, evecs[i], 1, evecs[i], 1);
                    norm += Vmath::Dot(ntot, evecs[i+1], 1, evecs[i+1], 1);
                    norm = std::sqrt(norm);
                    Vmath::Smul(ntot, 1.0/norm, evecs[i], 1, evecs[i], 1);
                    Vmath::Smul(ntot, 1.0/norm, evecs[i+1], 1, evecs[i+1], 1);
                    i++;
                }
            }
        }
    }
}
