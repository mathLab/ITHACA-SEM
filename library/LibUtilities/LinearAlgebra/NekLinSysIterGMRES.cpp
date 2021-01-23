///////////////////////////////////////////////////////////////////////////////
//
// File:  NekLinSysIterGMRES.cpp
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
// Description:  NekLinSysIterGMRES definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/LinearAlgebra/NekLinSysIterGMRES.h>

using namespace std;

namespace Nektar
{
namespace LibUtilities
{
/**
 * @class  NekLinSysIterGMRES
 *
 * Solves a linear system using iterative methods.
 */
string NekLinSysIterGMRES::className =
    LibUtilities::GetNekLinSysIterFactory().RegisterCreatorFunction(
        "GMRES", NekLinSysIterGMRES::create, "NekLinSysIterGMRES solver.");

NekLinSysIterGMRES::NekLinSysIterGMRES(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const LibUtilities::CommSharedPtr &vComm, const int nDimen,
    const NekSysKey &pKey)
    : NekLinSysIter(pSession, vComm, nDimen, pKey)
{
    std::vector<std::string> variables(1);
    variables[0]    = pSession->GetVariable(0);
    string variable = variables[0];

    pSession->MatchSolverInfo("GMRESLeftPrecond", "True", m_NekLinSysLeftPrecond,
                              pKey.m_NekLinSysLeftPrecond);
    pSession->MatchSolverInfo("GMRESRightPrecond", "False", m_NekLinSysRightPrecond,
                              pKey.m_NekLinSysRightPrecond);

    if (pSession->DefinesGlobalSysSolnInfo(variable, "GMRESMaxHessMatBand"))
    {
        m_KrylovMaxHessMatBand = boost::lexical_cast<int>(
            pSession->GetGlobalSysSolnInfo(variable, "GMRESMaxHessMatBand")
                .c_str());
    }
    else
    {
        pSession->LoadParameter("GMRESMaxHessMatBand", m_KrylovMaxHessMatBand,
                                m_LinSysMaxStorage + 1);
    }

    m_maxrestart = ceil(NekDouble(m_maxiter) / NekDouble(m_LinSysMaxStorage));
    m_LinSysMaxStorage = min(m_maxiter, m_LinSysMaxStorage);

    // cout 
    //     << " m_maxiter = " << m_maxiter
    //     << " m_maxrestart = " << m_maxrestart
    //     << " m_LinSysMaxStorage = " << m_LinSysMaxStorage
    //     << endl;

    int GMRESCentralDifference = 0;
    pSession->LoadParameter("GMRESCentralDifference", GMRESCentralDifference,
                            0);

    switch (GMRESCentralDifference)
    {
        case 1:
            m_DifferenceFlag0 = true;
            m_DifferenceFlag1 = false;
            break;
        case 2:
            m_DifferenceFlag0 = true;
            m_DifferenceFlag1 = true;
            break;

        default:
            m_DifferenceFlag0 = false;
            m_DifferenceFlag1 = false;
            break;
    }
}

void NekLinSysIterGMRES::v_InitObject()
{
    NekLinSysIter::v_InitObject();
}

NekLinSysIterGMRES::~NekLinSysIterGMRES()
{
}

/**
 *
 */
int NekLinSysIterGMRES::v_SolveSystem(
    const int nGlobal, const TensorOfArray1D<NekDouble> &pInput,
    Array<OneD, NekDouble> &pOutput, const int nDir, const NekDouble tol,
    const NekDouble factor)
{
    boost::ignore_unused(tol);

    m_tolerance     = max(tol, 1.0E-16);
    m_prec_factor   = factor;
    int niterations = DoGMRES(nGlobal, pInput, pOutput, nDir);

    return niterations;
}

/**  
 * Solve a global linear system using the Gmres 
 * We solve only for the non-Dirichlet modes. The operator is evaluated  
 * using an auxiliary function v_DoMatrixMultiply defined by the  
 * specific solver. Distributed math routines are used to support  
 * parallel execution of the solver.  
 *  
 * @param       pInput      Input residual  of all DOFs.  
 * @param       pOutput     Solution vector of all DOFs.  
 */

int NekLinSysIterGMRES::DoGMRES(const int nGlobal,
                                const TensorOfArray1D<NekDouble> &pInput,
                                Array<OneD, NekDouble> &pOutput, const int nDir)
{
    m_prec_factor = NekConstants::kNekUnsetDouble;

    if (m_rhs_magnitude == NekConstants::kNekUnsetDouble)
    {
        NekVector<NekDouble> inGlob(nGlobal, pInput, eWrapper);
        Set_Rhs_Magnitude(inGlob);
    }
    
    m_rhs_magnitude = 1.0;

    // Get vector sizes
    NekDouble eps = 0.0;
    int nNonDir   = nGlobal - nDir;

    Array<OneD, NekDouble> tmp;

    // zero homogeneous out array ready for solution updates
    // Should not be earlier in case input vector is same as
    // output and above copy has been peformed
    Vmath::Zero(nNonDir, tmp = pOutput + nDir, 1);

    m_totalIterations = 0;
    m_converged       = false;

    bool restarted = false;
    bool truncted  = false;

    if (m_KrylovMaxHessMatBand > 0)
    {
        truncted = true;
    }

    int nwidthcolm = 13;

    for (int nrestart = 0; nrestart < m_maxrestart; ++nrestart)
    {
        eps =
            DoGmresRestart(restarted, truncted, nGlobal, pInput, pOutput, nDir);

        if (m_converged)
        {
            break;
        }
        restarted = true;
    }

    if (m_verbose)
    {
        Array<OneD, NekDouble> r0(nGlobal, 0.0);
        m_operator.DoNekSysLhsEval(pOutput, r0, m_DifferenceFlag0);
        Vmath::Svtvp(nNonDir, -1.0, &r0[0] + nDir, 1, &pInput[0] + nDir, 1,
                     &r0[0] + nDir, 1);
        NekDouble vExchange = Vmath::Dot2(nNonDir, &r0[0] + nDir, &r0[0] + nDir,
                                          &m_map[0] + nDir);
        m_Comm->AllReduce(vExchange, LibUtilities::ReduceSum);
        NekDouble eps1 = vExchange;

        if (m_root)
        {
            cout << std::scientific << std::setw(nwidthcolm)
                 << std::setprecision(nwidthcolm - 8)
                 << "       GMRES iterations made = " << m_totalIterations
                 << " using tolerance of " << m_tolerance
                 << " (error = " << sqrt(eps * m_prec_factor / m_rhs_magnitude)
                 << ")";

            cout << " WITH (GMRES eps = " << eps << " REAL eps= " << eps1
                 << ")";

            if (m_converged)
            {
                cout << " CONVERGED" << endl;
            }
            else
            {
                cout << " WARNING: Exceeded maxIt" << endl;
            }
        }
    }

    WARNINGL1(m_converged, "GMRES did not converged !!");
    return m_totalIterations;
}

NekDouble NekLinSysIterGMRES::DoGmresRestart(
    const bool restarted, const bool truncted, const int nGlobal,
    const TensorOfArray1D<NekDouble> &pInput, Array<OneD, NekDouble> &pOutput,
    const int nDir)
{
    int nNonDir = nGlobal - nDir;

    // Allocate array storage of coefficients
    // Hessenburg matrix
    Array<OneD, Array<OneD, NekDouble>> hes(m_LinSysMaxStorage);
    for (int i = 0; i < m_LinSysMaxStorage; ++i)
    {
        hes[i] = Array<OneD, NekDouble>(m_LinSysMaxStorage + 1, 0.0);
    }
    // Hesseburg matrix after rotation
    Array<OneD, Array<OneD, NekDouble>> Upper(m_LinSysMaxStorage);
    for (int i = 0; i < m_LinSysMaxStorage; ++i)
    {
        Upper[i] = Array<OneD, NekDouble>(m_LinSysMaxStorage + 1, 0.0);
    }
    // Total search directions
    Array<OneD, Array<OneD, NekDouble>> V_total(m_LinSysMaxStorage + 1);
    for (int i = 0; i < m_LinSysMaxStorage + 1; ++i)
    {
        V_total[i] = Array<OneD, NekDouble>(nGlobal, 0.0);
    }
    // Residual
    Array<OneD, NekDouble> eta(m_LinSysMaxStorage + 1, 0.0);
    // Givens rotation c
    Array<OneD, NekDouble> cs(m_LinSysMaxStorage, 0.0);
    // Givens rotation s
    Array<OneD, NekDouble> sn(m_LinSysMaxStorage, 0.0);
    // Total coefficients, just for check
    Array<OneD, NekDouble> y_total(m_LinSysMaxStorage, 0.0);
    // Residual
    NekDouble eps;
    // Search direction order
    Array<OneD, int> id(m_LinSysMaxStorage, 0);
    Array<OneD, int> id_start(m_LinSysMaxStorage, 0);
    Array<OneD, int> id_end(m_LinSysMaxStorage, 0);
    // Temporary variables
    int idtem;
    int starttem;
    int endtem;

    NekDouble beta, alpha;
    NekDouble vExchange = 0;
    // Temporary Array
    Array<OneD, NekDouble> r0(nGlobal, 0.0);
    Array<OneD, NekDouble> tmp1;
    Array<OneD, NekDouble> tmp2;
    Array<OneD, NekDouble> Vsingle1;
    Array<OneD, NekDouble> Vsingle2;
    Array<OneD, NekDouble> hsingle1;
    Array<OneD, NekDouble> hsingle2;

    if (restarted)
    {
        // This is tmp2=Ax
        m_operator.DoNekSysLhsEval(pOutput, r0, m_DifferenceFlag0);

        // The first search direction
        beta = -1.0;
        // This is r0=b-AX
        Vmath::Svtvp(nNonDir, beta, &r0[0] + nDir, 1, &pInput[0] + nDir, 1,
                     &r0[0] + nDir, 1);
    }
    else
    {
        // If not restarted, x0 should be zero
        Vmath::Vcopy(nNonDir, &pInput[0] + nDir, 1, &r0[0] + nDir, 1);
    }

    if (m_NekLinSysLeftPrecond)
    {
        tmp1 = r0 + nDir;
        tmp2 = r0 + nDir;
        m_operator.DoNekSysPrecond(tmp1, tmp2);
    }

    // Norm of (r0)
    // The m_map tells how to connect
    vExchange =
        Vmath::Dot2(nNonDir, &r0[0] + nDir, &r0[0] + nDir, &m_map[0] + nDir);
    m_Comm->AllReduce(vExchange, LibUtilities::ReduceSum);
    eps = vExchange;

    if (!restarted)
    {
        if (m_prec_factor == NekConstants::kNekUnsetDouble)
        {
            if (m_NekLinSysLeftPrecond)
            {
                // Evaluate initial residual error for exit check
                vExchange = Vmath::Dot2(nNonDir, &pInput[0] + nDir,
                                        &pInput[0] + nDir, &m_map[0] + nDir);
                m_Comm->AllReduce(vExchange, LibUtilities::ReduceSum);
                m_prec_factor = vExchange / eps;
            }
            else
            {
                m_prec_factor = 1.0;
            }
        }
    }

    tmp2 = r0 + nDir;
    Vmath::Smul(nNonDir, sqrt(m_prec_factor), tmp2, 1, tmp2, 1);
    eps    = eps * m_prec_factor;
    eta[0] = sqrt(eps);

    // Give an order for the entries in Hessenburg matrix
    for (int nd = 0; nd < m_LinSysMaxStorage; ++nd)
    {
        id[nd]     = nd;
        id_end[nd] = nd + 1;
        starttem   = id_end[nd] - m_KrylovMaxHessMatBand;
        if (truncted && (starttem) > 0)
        {
            id_start[nd] = starttem;
        }
        else
        {
            id_start[nd] = 0;
        }
    }

    // Normlized by r0 norm V(:,1)=r0/norm(r0)
    alpha = 1.0 / eta[0];
    // Scalar multiplication
    Vmath::Smul(nNonDir, alpha, &r0[0] + nDir, 1, &V_total[0][0] + nDir, 1);

    // restarted Gmres(m) process
    int nswp = 0;
    if (m_NekLinSysRightPrecond)
    {
        Vsingle1 = Array<OneD, NekDouble>(nGlobal, 0.0);
    }

    // cout 
    //     << " m_tolerance = " << m_tolerance 
    //     << " m_rhs_magnitude = " << m_rhs_magnitude 
    //     << endl;

    for (int nd = 0; nd < m_LinSysMaxStorage; ++nd)
    {
        Vsingle2 = V_total[nd + 1];
        hsingle1 = hes[nd];

        if (m_NekLinSysRightPrecond)
        {
            tmp1 = V_total[nd] + nDir;
            tmp2 = Vsingle1 + nDir;
            m_operator.DoNekSysPrecond(tmp1, tmp2);
        }
        else
        {
            Vsingle1 = V_total[nd];
        }
        // w here is no need to add nDir due to temporary Array
        idtem    = id[nd];
        starttem = id_start[idtem];
        endtem   = id_end[idtem];

        DoArnoldi(starttem, endtem, nGlobal, nDir, V_total, Vsingle1, Vsingle2,
                  hsingle1);

        if (starttem > 0)
        {
            starttem = starttem - 1;
        }

        hsingle2 = Upper[nd];
        Vmath::Vcopy(m_LinSysMaxStorage + 1, &hsingle1[0], 1, &hsingle2[0], 1);
        DoGivensRotation(starttem, endtem, nGlobal, nDir, cs, sn, hsingle2,
                         eta);

        eps = eta[nd + 1] * eta[nd + 1];
        // This Gmres merge truncted Gmres to accelerate.
        // If truncted, cannot jump out because
        // the last term of eta is not residual
        if ((!truncted) || (nd < m_KrylovMaxHessMatBand))
        {
            // If (eps * m_prec_factor < m_tolerance *
            // m_tolerance * m_rhs_magnitude )
            if ((eps < m_tolerance * m_tolerance * m_rhs_magnitude) && nd > 0)
            {
                m_converged = true;
            }
            NekDouble tolmin = 1.0E-15;
            if (eps < tolmin * tolmin * m_rhs_magnitude)
            {
                m_converged = true;
            }
        }
        nswp++;
        m_totalIterations++;

        if (m_converged)
        {
            break;
        }
    }

    DoBackward(nswp, Upper, eta, y_total);
    // calculate output y_total*V_total
    Array<OneD, NekDouble> solution(nGlobal, 0.0);
    for (int i = 0; i < nswp; ++i)
    {
        beta = y_total[i];
        Vmath::Svtvp(nNonDir, beta, &V_total[i][0] + nDir, 1,
                     &solution[0] + nDir, 1, &solution[0] + nDir, 1);
    }

    if (m_NekLinSysRightPrecond)
    {
        tmp1 = solution + nDir;
        tmp2 = solution + nDir;
        m_operator.DoNekSysPrecond(tmp1, tmp2);
    }
    Vmath::Vadd(nNonDir, &solution[0] + nDir, 1, &pOutput[0] + nDir, 1,
                &pOutput[0] + nDir, 1);

    return eps;
}

// Arnoldi Subroutine
void NekLinSysIterGMRES::DoArnoldi(const int starttem, const int endtem,
                                   const int nGlobal, const int nDir,
                                   // V_total(:,1:nd)
                                   Array<OneD, Array<OneD, NekDouble>> &V_local,
                                   // V[nd]
                                   Array<OneD, NekDouble> &Vsingle1,
                                   // V[nd+1]
                                   Array<OneD, NekDouble> &Vsingle2,
                                   // h
                                   Array<OneD, NekDouble> &hsingle)
{
    // To notice, V_local's order not certainly equal to starttem:endtem
    // starttem:endtem is the entry position in Hessenburg matrix
    NekDouble alpha, beta;
    Array<OneD, NekDouble> tmp1, tmp2;
    int numbertem;
    int nNonDir = nGlobal - nDir;
    // Later for parallel
    NekDouble vExchange = 0.0;
    // w=AV(:,nd)
    Array<OneD, NekDouble> w(nGlobal, 0.0);

    m_operator.DoNekSysLhsEval(Vsingle1, w, m_DifferenceFlag1);

    tmp1 = w + nDir;
    tmp2 = w + nDir;
    if (m_NekLinSysLeftPrecond)
    {
        m_operator.DoNekSysPrecond(tmp1, tmp2);
    }

    Vmath::Smul(nNonDir, sqrt(m_prec_factor), tmp2, 1, tmp2, 1);

    // Modified Gram-Schmidt
    // The pointer not certainly equal to starttem.
    // Like initially, Gmres-deep need to use numbertem=0
    numbertem = starttem;
    for (int i = starttem; i < endtem; ++i)
    {
        vExchange =
            Vmath::Dot2(nNonDir, &w[0] + nDir, &V_local[numbertem][0] + nDir,
                        &m_map[0] + nDir);
        m_Comm->AllReduce(vExchange, LibUtilities::ReduceSum);

        hsingle[i] = vExchange;

        beta = -1.0 * vExchange;
        Vmath::Svtvp(nNonDir, beta, &V_local[numbertem][0] + nDir, 1,
                     &w[0] + nDir, 1, &w[0] + nDir, 1);
        numbertem = numbertem + 1;
    }
    // end of Modified Gram-Schmidt

    // calculate the L2 norm and normalize
    vExchange =
        Vmath::Dot2(nNonDir, &w[0] + nDir, &w[0] + nDir, &m_map[0] + nDir);

    m_Comm->AllReduce(vExchange, LibUtilities::ReduceSum);

    hsingle[endtem] = sqrt(vExchange);

    alpha = 1.0 / hsingle[endtem];
    Vmath::Smul(nNonDir, alpha, &w[0] + nDir, 1, &Vsingle2[0] + nDir, 1);
}

// QR factorization through Givens rotation
void NekLinSysIterGMRES::DoGivensRotation(const int starttem, const int endtem,
                                          const int nGlobal, const int nDir,
                                          Array<OneD, NekDouble> &c,
                                          Array<OneD, NekDouble> &s,
                                          Array<OneD, NekDouble> &hsingle,
                                          Array<OneD, NekDouble> &eta)
{
    boost::ignore_unused(nGlobal, nDir);
    NekDouble temp_dbl;
    NekDouble dd;
    NekDouble hh;
    int idtem = endtem - 1;
    // The starttem and endtem are beginning and ending order of Givens rotation
    // They usually equal to the beginning position and ending position of
    // Hessenburg matrix But sometimes starttem will change, like if it is
    // initial 0 and becomes nonzero because previous Givens rotation See Yu
    // Pan's User Guide
    for (int i = starttem; i < idtem; ++i)
    {
        temp_dbl       = c[i] * hsingle[i] - s[i] * hsingle[i + 1];
        hsingle[i + 1] = s[i] * hsingle[i] + c[i] * hsingle[i + 1];
        hsingle[i]     = temp_dbl;
    }
    dd = hsingle[idtem];
    hh = hsingle[endtem];
    if (hh == 0.0)
    {
        c[idtem] = 1.0;
        s[idtem] = 0.0;
    }
    else if (abs(hh) > abs(dd))
    {
        temp_dbl = -dd / hh;
        s[idtem] = 1.0 / sqrt(1.0 + temp_dbl * temp_dbl);
        c[idtem] = temp_dbl * s[idtem];
    }
    else
    {
        temp_dbl = -hh / dd;
        c[idtem] = 1.0 / sqrt(1.0 + temp_dbl * temp_dbl);
        s[idtem] = temp_dbl * c[idtem];
    }

    hsingle[idtem]  = c[idtem] * hsingle[idtem] - s[idtem] * hsingle[endtem];
    hsingle[endtem] = 0.0;

    temp_dbl    = c[idtem] * eta[idtem] - s[idtem] * eta[endtem];
    eta[endtem] = s[idtem] * eta[idtem] + c[idtem] * eta[endtem];
    eta[idtem]  = temp_dbl;
}

// Backward calculation
// To notice, Hesssenburg matrix's column
// and row changes due to use Array<OneD,Array<OneD,NekDouble>>
void NekLinSysIterGMRES::DoBackward(const int number,
                                    Array<OneD, Array<OneD, NekDouble>> &A,
                                    const TensorOfArray1D<NekDouble> &b,
                                    Array<OneD, NekDouble> &y)
{
    // Number is the entry number
    // but C++'s order need to be one smaller
    int maxid = number - 1;
    NekDouble sum;
    y[maxid] = b[maxid] / A[maxid][maxid];

    for (int i = maxid - 1; i > -1; --i)
    {
        sum = b[i];

        for (int j = i + 1; j < number; ++j)
        {
            // i and j changes due to use Array<OneD,Array<OneD,NekDouble>>
            sum = sum - y[j] * A[j][i];
        }
        y[i] = sum / A[i][i];
    }
}
} // namespace LibUtilities
} // namespace Nektar
