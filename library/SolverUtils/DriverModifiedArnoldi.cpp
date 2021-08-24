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

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/DriverModifiedArnoldi.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{

string DriverModifiedArnoldi::className =
        GetDriverFactory().RegisterCreatorFunction("ModifiedArnoldi",
                                    DriverModifiedArnoldi::create);
string DriverModifiedArnoldi::driverLookupId =
        LibUtilities::SessionReader::RegisterEnumValue("Driver",
                                    "ModifiedArnoldi",0);

/**
 *
 */
DriverModifiedArnoldi::DriverModifiedArnoldi(
    const LibUtilities::SessionReaderSharedPtr pSession,
    const SpatialDomains::MeshGraphSharedPtr pGraph)
    : DriverArnoldi(pSession, pGraph)
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
    if (m_comm->GetRank() == 0)
    {
        out << "\tArnoldi solver type    : Modified Arnoldi" << endl;
    }

    DriverArnoldi::ArnoldiSummary(out);

    for( int i = 0; i < m_nequ; ++i)
    {
        m_equ[i]->DoInitialise();
    }

    //FwdTrans Initial conditions to be in Coefficient Space
    m_equ[m_nequ-1] ->TransPhysToCoeff();

}


/**
 *
 */
void DriverModifiedArnoldi::v_Execute(ostream &out)
{
    int i               = 0;
    int j               = 0;
    int nq              = m_equ[0]->UpdateFields()[0]->GetNcoeffs();
    int ntot            = m_nfields*nq;
    int converged       = 0;
    NekDouble resnorm   = 0.0;
    ofstream evlout;
    std::string evlFile = m_session->GetSessionName() + ".evl";

    if (m_comm->GetRank() == 0)
    {
        evlout.open(evlFile.c_str());
    }

    // Allocate memory
    Array<OneD, NekDouble> alpha = Array<OneD, NekDouble> (m_kdim+1,      0.0);
    Array<OneD, NekDouble> wr    = Array<OneD, NekDouble> (m_kdim,        0.0);
    Array<OneD, NekDouble> wi    = Array<OneD, NekDouble> (m_kdim,        0.0);
    Array<OneD, NekDouble> zvec  = Array<OneD, NekDouble> (m_kdim*m_kdim, 0.0);

    Array<OneD, Array<OneD, NekDouble> > Kseq
            = Array<OneD, Array<OneD, NekDouble> > (m_kdim + 1);
    Array<OneD, Array<OneD, NekDouble> > Kseqcopy
            = Array<OneD, Array<OneD, NekDouble> > (m_kdim + 1);
    Array<OneD, Array<OneD, NekDouble> > Tseq
            = Array<OneD, Array<OneD, NekDouble> > (m_kdim + 1);
    for (i = 0; i < m_kdim + 1; ++i)
    {
        Kseq[i] = Array<OneD, NekDouble>(ntot, 0.0);
        if (m_useMask)
        {
            Kseqcopy[i] = Array<OneD, NekDouble>(ntot, 0.0);
        }
        else
        {
            Kseqcopy[i] = Kseq[i];
        }
        Tseq[i] = Array<OneD, NekDouble>(ntot, 0.0);
    }

    // Copy starting vector into second sequence element (temporary).
    if(m_session->DefinesFunction("InitialConditions"))
    {
        if (m_comm->GetRank() == 0)
        {
            out << "\tInital vector       : specified in input file " << endl;
        }
        m_equ[0]->SetInitialConditions(0.0,false);

        CopyFieldToArnoldiArray(Kseq[1]);
    }
    else
    {
        if (m_comm->GetRank() == 0)
        {
            out << "\tInital vector       : random  " << endl;
        }

        NekDouble eps=0.0001;
        Vmath::FillWhiteNoise(ntot, eps , &Kseq[1][0], 1);
        if (m_useMask)
        {
            Vmath::Vmul(ntot, Kseq[1], 1, GetMaskCoeff(), 1, Kseq[1], 1);
        }
    }

    // Perform one iteration to enforce boundary conditions.
    // Set this as the initial value in the sequence.
    EV_update(Kseq[1], Kseq[0]);
    if (m_comm->GetRank() == 0)
    {
        out << "Iteration: " << 0 <<  endl;
    }

    // Normalise first vector in sequence
    if (m_useMask)
    {
        Vmath::Vmul(ntot, Kseq[0], 1, GetMaskCoeff(), 1, Kseqcopy[0], 1);
    }
    alpha[0] = Blas::Ddot(ntot, &Kseqcopy[0][0], 1, &Kseqcopy[0][0], 1);
    m_comm->AllReduce(alpha[0], Nektar::LibUtilities::ReduceSum);
    alpha[0] = std::sqrt(alpha[0]);
    Vmath::Smul(ntot, 1.0/alpha[0], Kseq[0], 1, Kseq[0], 1);

    // Fill initial krylov sequence
    NekDouble resid0;
    for (i = 1; !converged && i <= m_kdim; ++i)
    {
        // Compute next vector
        EV_update(Kseq[i-1], Kseq[i]);

        // Normalise
        if (m_useMask)
        {
            Vmath::Vmul(ntot, Kseq[i], 1, GetMaskCoeff(), 1, Kseqcopy[i], 1);
        }
        alpha[i] = Blas::Ddot(ntot, &Kseqcopy[i][0], 1, &Kseqcopy[i][0], 1);
        m_comm->AllReduce(alpha[i], Nektar::LibUtilities::ReduceSum);
        alpha[i] = std::sqrt(alpha[i]);

        //alpha[i] = std::sqrt(alpha[i]);
        Vmath::Smul(ntot, 1.0/alpha[i], Kseq[i], 1, Kseq[i], 1);

        // Copy Krylov sequence into temporary storage
        for (int k = 0; k < i + 1; ++k)
        {
            if (m_useMask)
            {
                Vmath::Vmul(ntot, Kseq[k], 1, GetMaskCoeff(), 1, Tseq[k], 1);
                Vmath::Vcopy(ntot, Kseq[k], 1, Kseqcopy[k], 1);
            }
            else
            {
                Vmath::Vcopy(ntot, Kseq[k], 1, Tseq[k], 1);;
            }
        }

        // Generate Hessenberg matrix and compute eigenvalues of it.
        EV_small(Tseq, Kseqcopy, ntot, alpha, i, zvec, wr, wi, resnorm);

        // Test for convergence.
        converged = EV_test(i, i, zvec, wr, wi, resnorm,
            std::min(i, m_nvec), evlout, resid0);
        if ( i >= m_nvec)
        {
            converged = max (converged, 0);
        }
        else
        {
            converged = 0;
        }

        if (m_comm->GetRank() == 0)
        {
            out << "Iteration: " <<  i << " (residual : " << resid0
            << ")" <<endl;
        }
    }

    // Continue with full sequence
    if (!converged)
    {
        for (i = m_kdim + 1; !converged && i <= m_nits; ++i)
        {
            // Shift all the vectors in the sequence.
            // First vector is removed.
            for (int j = 1; j <= m_kdim; ++j)
            {
                alpha[j-1] = alpha[j];
                Vmath::Vcopy(ntot, Kseq[j], 1, Kseq[j-1], 1);
            }

            // Compute next vector
            EV_update(Kseq[m_kdim - 1], Kseq[m_kdim]);

            // Compute new scale factor
            if (m_useMask)
            {
                Vmath::Vmul(ntot, Kseq[m_kdim], 1, GetMaskCoeff(), 1,
                                  Kseqcopy[m_kdim], 1);
            }
            alpha[m_kdim] = Blas::Ddot(ntot, &Kseqcopy[m_kdim][0], 1,
                                             &Kseqcopy[m_kdim][0], 1);
            m_comm->AllReduce(alpha[m_kdim], Nektar::LibUtilities::ReduceSum);
            alpha[m_kdim] = std::sqrt(alpha[m_kdim]);
            Vmath::Smul(ntot, 1.0/alpha[m_kdim], Kseq[m_kdim], 1,
                                                 Kseq[m_kdim], 1);

            // Copy Krylov sequence into temporary storage
            for (int k = 0; k < m_kdim + 1; ++k)
            {
                if (m_useMask)
                {
                    Vmath::Vmul(ntot, Kseq[k], 1, GetMaskCoeff(), 1, Tseq[k], 1);
                    Vmath::Vcopy(ntot, Kseq[k], 1, Kseqcopy[k], 1);
                }
                else
                {
                    Vmath::Vcopy(ntot, Kseq[k], 1, Tseq[k], 1);;
                }
            }

            // Generate Hessenberg matrix and compute eigenvalues of it
            EV_small(Tseq, Kseqcopy, ntot, alpha, m_kdim, zvec, wr, wi, resnorm);

            // Test for convergence.
            converged = EV_test(i, m_kdim, zvec, wr, wi, resnorm,
                                m_nvec, evlout, resid0);

            if (m_comm->GetRank() == 0)
            {
                out << "Iteration: " <<  i << " (residual : "
                    << resid0 << ")" <<endl;
            }
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
            out << "L 2 error (variable " << m_equ[0]->GetVariable(j)
            << ") : " << vL2Error << endl;
            out << "L inf error (variable " << m_equ[0]->GetVariable(j)
            << ") : " << vLinfError << endl;
        }
    }

    // Process eigenvectors and write out.
    m_nvec = converged;
    EV_post(Tseq, Kseqcopy, ntot, min(--i, m_kdim), m_nvec, zvec, wr, wi,
            converged);

    WARNINGL0(m_imagShift == 0,"Complex Shift applied. "
              "Need to implement Ritz re-evaluation of"
              "eigenvalue. Only one half of complex "
              "value will be correct");

    // store eigenvalues so they can be accessed from driver class
    m_real_evl = wr;
    m_imag_evl = wi;

    // Close the runtime info file.
    if (m_comm->GetRank() == 0)
    {
        evlout.close();
    }
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

    m_equ[0]->SetTime(0.);
    m_equ[0]->DoSolve();

    if(m_EvolutionOperator == eTransientGrowth)
    {
        Array<OneD, MultiRegions::ExpListSharedPtr> fields;
        fields = m_equ[0]->UpdateFields();

        //start Adjoint with latest fields of direct
        CopyFwdToAdj();
        m_equ[1]->TransCoeffToPhys();

        m_equ[1]->SetTime(0.);
        m_equ[1]->DoSolve();
    }

    // Copy starting vector into first sequence element.
    CopyFieldToArnoldiArray(tgt);
}


/**
 * Computes the Ritz eigenvalues and eigenvectors of the Hessenberg matrix
 * constructed from the Krylov sequence.
 */
void DriverModifiedArnoldi::EV_small(
    Array<OneD, Array<OneD, NekDouble> > &Kseq,
    Array<OneD, Array<OneD, NekDouble> > &Kseqcopy,
    const int                             ntot,
    const Array<OneD, NekDouble>         &alpha,
    const int                             kdim,
    Array<OneD, NekDouble>               &zvec,
    Array<OneD, NekDouble>               &wr,
    Array<OneD, NekDouble>               &wi,
    NekDouble                            &resnorm)
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
        NekDouble gsc = Blas::Ddot(ntot, &Kseq[i][0], 1, &Kseq[i][0], 1);
        m_comm->AllReduce(gsc, Nektar::LibUtilities::ReduceSum);
        gsc = std::sqrt(gsc);
        ASSERTL0(gsc != 0.0, "Vectors are linearly independent.");

        R[i*kdimp+i] = gsc;
        Vmath::Smul(ntot, 1.0/gsc, Kseq[i], 1, Kseq[i], 1);
        if (m_useMask)
        {
            Vmath::Smul(ntot, 1.0/gsc, Kseqcopy[i], 1, Kseqcopy[i], 1);
        }

        for (int j = i + 1; j < kdimp; ++j)
        {
            gsc = Blas::Ddot(ntot, &Kseq[i][0], 1, &Kseq[j][0], 1);
            m_comm->AllReduce(gsc, Nektar::LibUtilities::ReduceSum);
            Vmath::Svtvp(ntot, -gsc, Kseq[i], 1, Kseq[j], 1, Kseq[j], 1);
            if (m_useMask)
            {
                Vmath::Svtvp(ntot, -gsc, Kseqcopy[i], 1, Kseqcopy[j], 1, Kseqcopy[j], 1);
            }
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

    Lapack::dgeev_('N', 'V', kdim, &H[0], kdim, &wr[0], &wi[0], 0, 1,
                   &zvec[0], kdim, &rwork[0], lwork, ier);

    ASSERTL0(!ier, "Error with dgeev");

    resnorm = H[(kdim-1)*kdim + kdim];
}


/**
 * Check for convergence of the residuals of the eigenvalues of H.
 */
int DriverModifiedArnoldi::EV_test(
    const int               itrn,
    const int               kdim,
    Array<OneD, NekDouble> &zvec,
    Array<OneD, NekDouble> &wr,
    Array<OneD, NekDouble> &wi,
    const NekDouble         resnorm,
    int                     nvec,
    ofstream               &evlout,
    NekDouble              &resid0)
{
    int idone = 0;
    // NekDouble period = 0.1;

    Array<OneD, NekDouble> resid(kdim);
    for (int i = 0; i < kdim; ++i)
    {
        NekDouble tmp = std::sqrt(Vmath::Dot(kdim, &zvec[0] + i*kdim, 1,
                                             &zvec[0] + i*kdim, 1));
        resid[i] = resnorm * std::fabs(zvec[kdim - 1 + i*kdim]) / tmp;
        if (wi[i] < 0.0)
        {
            resid[i-1] = resid[i] = hypot(resid[i-1], resid[i]);
        }
    }

    EV_sort(zvec, wr, wi, resid, kdim);

    while (nvec <= kdim && resid[nvec-1] < m_evtol)
    {
        idone = nvec;
        ++nvec;
    }
    nvec -= (idone > 0);

    if (m_comm->GetRank() == 0)
    {
        evlout << "-- Iteration = " << itrn << ", H(k+1, k) = "
               << resnorm << endl;
        evlout.precision(4);
        evlout.setf(ios::scientific, ios::floatfield);
        if(m_timeSteppingAlgorithm)
        {
            evlout << "        Magnitude   Angle       Growth      "
                   << "Frequency   Residual" << endl;
        }
        else
        {
            evlout << "        Real        Imaginary   inverse real  "
                   << "inverse imag  Residual" << endl;
        }

        for (int i = 0; i < kdim; i++)
        {
            WriteEvs(evlout,i,wr[i],wi[i],resid[i]);
        }
    }

    resid0 = resid[nvec-1];
    return idone;
}


/**
 * Sorts the computed eigenvalues by smallest residual first.
 */
void DriverModifiedArnoldi::EV_sort(
    Array<OneD, NekDouble> &evec,
    Array<OneD, NekDouble> &wr,
    Array<OneD, NekDouble> &wi,
    Array<OneD, NekDouble> &test,
    const int               dim)
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
}


/**
 * Post-process the Ritz eigenvalues/eigenvectors of the matrix H, to compute
 * estimations of the leading eigenvalues and eigenvectors of the original
 * matrix.
 */
void DriverModifiedArnoldi::EV_post(
    Array<OneD, Array<OneD, NekDouble> > &Tseq,
    Array<OneD, Array<OneD, NekDouble> > &Kseq,
    const int                             ntot,
    const int                             kdim,
    const int                             nvec,
    Array<OneD, NekDouble>               &zvec,
    Array<OneD, NekDouble>               &wr,
    Array<OneD, NekDouble>               &wi,
    const int                             icon)
{
    if (icon == 0)
    {
        // Not converged, write final Krylov vector
        ASSERTL0(false, "Convergence was not achieved within the "
                        "prescribed number of iterations.");
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
            std::string file = m_session->GetSessionName() + "_eig_"
                + boost::lexical_cast<std::string>(j)
                + ".fld";

            if (m_comm->GetRank() == 0)
            {
                WriteEvs(cout, j, wr[j], wi[j]);
            }
            WriteFld(file,Kseq[j]);
            if (m_useMask)
            {
                std::string fileunmask = m_session->GetSessionName() + "_eig_masked_"
                    + boost::lexical_cast<std::string>(j)
                    + ".fld";
                WriteFld(fileunmask,Tseq[j]);
            }
        }
    }
    else
    {
        // Not recognised
        ASSERTL0(false, "Unrecognised value.");
    }
}


/**
 * Compute estimates of the eigenvalues/eigenvectors of the original system.
 */
void DriverModifiedArnoldi::EV_big(
    Array<OneD, Array<OneD, NekDouble> > &bvecs,
    Array<OneD, Array<OneD, NekDouble> > &evecs,
    const int                             ntot,
    const int                             kdim,
    const int                             nvec,
    Array<OneD, NekDouble>               &zvec,
    Array<OneD, NekDouble>               &wr,
    Array<OneD, NekDouble>               &wi)
{
    boost::ignore_unused(wr);

    NekDouble wgt, norm;
    Array<OneD, Array<OneD, NekDouble> > btmp(nvec);
    Array<OneD, Array<OneD, NekDouble> > etmp(nvec);
    for(int i=0; i<nvec; ++i)
    {
        if (m_useMask)
        {
            btmp[i] = Array<OneD, NekDouble>(ntot, 0.);
            etmp[i] = Array<OneD, NekDouble>(ntot, 0.);
        }
        else
        {
            btmp[i] = evecs[i];
        }
    }

    // Generate big eigenvectors
    for (int j = 0; j < nvec; ++j)
    {
        Vmath::Zero(ntot, btmp[j], 1);
        if (m_useMask)
        {
            Vmath::Zero(ntot, etmp[j], 1);
        }
        for (int i = 0; i < kdim; ++i)
        {
            wgt = zvec[i + j*kdim];
            Vmath::Svtvp(ntot, wgt, bvecs[i], 1, btmp[j], 1, btmp[j], 1);
            if (m_useMask)
            {
                Vmath::Svtvp(ntot, wgt, evecs[i], 1, etmp[j], 1, etmp[j], 1);
            }
        }
    }

    // Normalise the big eigenvectors
    for (int i = 0; i < nvec; ++i)
    {
        if (wi[i] == 0.0)   // Real mode
        {
            if (m_session->GetComm()->GetRank() == 0)
            {
                cout << "eigenvalue " << i << ": real mode" << endl;
            }
            norm = Blas::Ddot(ntot, &btmp[i][0], 1, &btmp[i][0], 1);
            m_comm->AllReduce(norm, Nektar::LibUtilities::ReduceSum);
            norm = std::sqrt(norm);
            if (m_useMask)
            {
                Vmath::Smul(ntot, 1.0/norm, btmp[i], 1, bvecs[i], 1);
                Vmath::Smul(ntot, 1.0/norm, etmp[i], 1, evecs[i], 1);
            }
            else
            {
                Vmath::Smul(ntot, 1.0/norm, btmp[i], 1, evecs[i], 1);
            }
        }
        else
        {
            if (m_session->GetComm()->GetRank() == 0)
            {
                cout << "eigenvalues " << i << ", " << i+1
                     <<  ": complex modes" << endl;
            }
            norm  = Blas::Ddot(ntot, &btmp[i][0],   1, &btmp[i][0],   1);
            norm += Blas::Ddot(ntot, &btmp[i+1][0], 1, &btmp[i+1][0], 1);
            m_comm->AllReduce(norm, Nektar::LibUtilities::ReduceSum);
            norm = std::sqrt(norm);

            if (m_useMask)
            {
                Vmath::Smul(ntot, 1.0/norm, btmp[i],   1, bvecs[i],   1);
                Vmath::Smul(ntot, 1.0/norm, btmp[i+1], 1, bvecs[i+1], 1);
                Vmath::Smul(ntot, 1.0/norm, etmp[i],   1, evecs[i],   1);
                Vmath::Smul(ntot, 1.0/norm, etmp[i+1], 1, evecs[i+1], 1);
            }
            else
            {
                Vmath::Smul(ntot, 1.0/norm, btmp[i],   1, evecs[i],   1);
                Vmath::Smul(ntot, 1.0/norm, btmp[i+1], 1, evecs[i+1], 1);
            }

            i++;
        }
    }
}

}
}
