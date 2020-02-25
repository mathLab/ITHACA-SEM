///////////////////////////////////////////////////////////////////////////////
//
// File:  NekLinSysIteratGMRES.cpp
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
// Description:  NekLinSysIteratGMRES definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/NekLinSysIteratGMRES.h>
#include <LibUtilities/BasicUtils/Timer.h>

using namespace std;

namespace Nektar
{      
    namespace LibUtilities
    {  
        /**
         * @class  NekLinSysIteratGMRES
         *
         * Solves a linear system using iterative methods.
         */
        string NekLinSysIteratGMRES::className =
        LibUtilities::GetNekLinSysIteratFactory().RegisterCreatorFunction(
            "GMRES", NekLinSysIteratGMRES::create,
            "NekLinSysIteratGMRES solver.");

        /// Constructor for full direct matrix solve.
        NekLinSysIteratGMRES::NekLinSysIteratGMRES(
            const LibUtilities::SessionReaderSharedPtr  &pSession,
            const LibUtilities::CommSharedPtr           &vComm,
            const int                                   nDimen)
            : NekLinSysIterat(pSession, vComm, nDimen)
        {
            m_maxstorage        =   30;
            m_maxhesband        =   30;
            // m_maxiter           =   60;
            std::vector<std::string>  variables(1);
            variables[0] =  pSession->GetVariable(0);
            string variable = variables[0];

            pSession->MatchSolverInfo(
                "flag_LeftPrecond", "True", m_flag_LeftPrecond, false);
            pSession->MatchSolverInfo(
                "flag_RightPrecond", "False", m_flag_RightPrecond, true);
            
            if(pSession->DefinesGlobalSysSolnInfo(variable,
                                                "MaxStorage"))
            {
                m_maxstorage = boost::lexical_cast<int>(
                        pSession->GetGlobalSysSolnInfo(variable,
                                "MaxStorage").c_str());
            }
            else
            {
                pSession->LoadParameter("MaxStorage",
                                        m_maxstorage,
                                        30);
            }
            if(pSession->DefinesGlobalSysSolnInfo(variable , 
                                                "MaxHesband"))
            {
                m_maxhesband = boost::lexical_cast<int>(
                        pSession->GetGlobalSysSolnInfo(variable,
                                "MaxHesband").c_str());
            }
            else
            {
                pSession->LoadParameter("MaxHesband",
                                        m_maxhesband,
                                        m_maxstorage+1);
            }

            m_maxrestart = ceil(NekDouble(m_maxiter)/NekDouble(m_maxstorage));

            int flaguseCentralDifference = 0;
            pSession->LoadParameter("flaguseCentralDifference",
                                        flaguseCentralDifference,
                                        0);
            
            switch (flaguseCentralDifference)
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

        void NekLinSysIteratGMRES::v_InitObject()
        {
            NekLinSysIterat::v_InitObject();
        }


        NekLinSysIteratGMRES::~NekLinSysIteratGMRES()
        {
        }

        /**
         *
         */
        int NekLinSysIteratGMRES::v_SolveSystem(
            const int                           nGlobal,
            const Array<OneD, const NekDouble>  &pInput,
            Array<OneD,      NekDouble>         &pOutput,
            const int                           nDir,
            const NekDouble                     tol,
            const NekDouble                     factor)
        {
            boost::ignore_unused(tol);

            m_tolerance = max(tol,1.0E-16);
            m_prec_factor = factor;
            int niterations = DoGMRES(nGlobal, pInput, pOutput, nDir);

            return niterations;
        }

        bool NekLinSysIteratGMRES::v_ConvergenceCheck(
                const int                           nIteration,
                const Array<OneD, const NekDouble>  &Residual,
                const NekDouble                     tol         )
        {
            bool converged = false;
            boost::ignore_unused(nIteration,Residual,tol);
            return converged;
        }

        /**  
        * Solve a global linear system using the Gmres 
        * We solve only for the non-Dirichlet modes. The operator is evaluated  
        * using an auxiliary function v_DoMatrixMultiply defined by the  
        * specific solver. Distributed math routines are used to support  
        * parallel execution of the solver.  
        *  
        * The implemented algorithm uses a reduced-communication reordering of  
        * the standard PCG method (Demmel, Heath and Vorst, 1993)  
        *  
        * @param       pInput      Input residual  of all DOFs.  
        * @param       pOutput     Solution vector of all DOFs.  
        */

        int  NekLinSysIteratGMRES::DoGMRES(
            const int                          nGlobal,
            const Array<OneD, const NekDouble> &pInput,
            Array<OneD,      NekDouble> &pOutput,
            const int                          nDir)
        {
            m_prec_factor = NekConstants::kNekUnsetDouble;
            // m_rhs_magnitude = NekConstants::kNekUnsetDouble;

            if(m_rhs_magnitude == NekConstants::kNekUnsetDouble)
            {
                NekVector<NekDouble> inGlob (nGlobal, pInput, eWrapper);
                Set_Rhs_Magnitude(inGlob);
            }

            // Get vector sizes
            NekDouble eps = 0.0;
            int nNonDir = nGlobal - nDir;

            // zero homogeneous out array ready for solution updates
            // Should not be earlier in case input vector is same as
            // output and above copy has been peformed
            Vmath::Zero(nNonDir, &pOutput[nDir], 1);
            //Vmath::Zero(nGlobal, &qk_a[iqk0],1);


            m_totalIterations = 0;
            m_converged       = false;

            bool restarted = false;
            bool truncted = false;

            // m_maxstorage = 5;
            // m_maxrestart = m_maxiter / m_maxstorage + 1;
            // default truncted Gmres(m) closed
            // m_maxhesband = 0;
            if(m_maxhesband > 0)
            {
                truncted = true;
            }

            int nwidthcolm = 13;

            for(int nrestart = 0; nrestart < m_maxrestart; ++nrestart)
            {
                eps = DoGmresRestart(restarted,
                                    truncted,
                                    nGlobal,
                                    pInput,
                                    pOutput,
                                    nDir);

                if(m_converged)
                {
                    break;
                }
                restarted = true;
            }


            if(m_verbose)
            {
                Array<OneD, NekDouble> r0(nGlobal, 0.0);
                m_operator.DoNonlinLinSysLhsEval(pOutput, r0,m_DifferenceFlag0);
                Vmath::Svtvp(nNonDir, -1.0, &r0[0] + nDir, 1, &pInput[0] + nDir, 1, &r0[0] + nDir, 1);
                NekDouble vExchange    = Vmath::Dot2(nNonDir,
                                            &r0[0] + nDir,
                                            &r0[0] + nDir,
                                            &m_map[0] + nDir);
                m_Comm->AllReduce(vExchange, LibUtilities::ReduceSum);
                NekDouble eps1 = vExchange;

                if (m_root)
                {
                    cout <<std::scientific<<std::setw(nwidthcolm)<<std::setprecision(nwidthcolm-8) 
                            << "       GMRES iterations made = " << m_totalIterations
                            << " using tolerance of "  << m_tolerance
                            << " (error = " << sqrt(eps * m_prec_factor / m_rhs_magnitude) << ")";
                    
                    cout << " WITH (GMRES eps = " << eps << " REAL eps= "<<eps1<<")";

                    if(m_converged)
                    {
                        cout <<" CONVERGED"<<endl;
                    }
                    else
                    {
                        cout <<" WARNING: Exceeded maxIt"<<endl;
                    }
                }
            }

            WARNINGL0(m_converged, "GMRES did not converged !!");
            return m_totalIterations;
        }

        NekDouble  NekLinSysIteratGMRES::DoGmresRestart(
            const bool                          restarted,
            const bool                          truncted,
            const int                           nGlobal,
            const Array<OneD, const NekDouble>  &pInput,
            Array<OneD,      NekDouble>         &pOutput,
            const int                           nDir)
        {
            int nNonDir = nGlobal - nDir;

            // Allocate array storage of coefficients
            // Hessenburg matrix
            Array<OneD, Array<OneD, NekDouble> > hes    (m_maxstorage);
            for (int i = 0; i < m_maxstorage; ++i)
            {
                hes[i] = Array<OneD, NekDouble>(m_maxstorage + 1, 0.0);
            }
            // Hesseburg matrix after rotation
            Array<OneD, Array<OneD, NekDouble> >  Upper  (m_maxstorage);
            for (int i = 0; i < m_maxstorage; ++i)
            {
                Upper[i] = Array<OneD, NekDouble>(m_maxstorage + 1, 0.0);
            }
            // Total search directions
            Array<OneD, Array<OneD, NekDouble> >  V_total(m_maxstorage + 1);
            for (int i = 0; i < m_maxstorage + 1; ++i)
            {
                V_total[i] = Array<OneD, NekDouble>(nGlobal, 0.0);
            }
            //Residual
            Array<OneD, NekDouble> eta    (m_maxstorage + 1, 0.0);
            //Givens rotation c
            Array<OneD, NekDouble> cs     (m_maxstorage, 0.0);
            //Givens rotation s
            Array<OneD, NekDouble> sn     (m_maxstorage, 0.0);
            //Total coefficients, just for check
            Array<OneD, NekDouble> y_total     (m_maxstorage, 0.0);
            // Residual
            NekDouble eps;
            //Search direction order
            Array<OneD, int> id (m_maxstorage, 0);
            Array<OneD, int> id_start (m_maxstorage, 0);
            Array<OneD, int> id_end (m_maxstorage, 0);
            // temporary variables
            int idtem, starttem, endtem;
            NekDouble beta, alpha;
            NekDouble vExchange = 0;
            // temporary Array
            Array<OneD, NekDouble> r0(nGlobal, 0.0);
            Array<OneD, NekDouble> tmp1;
            Array<OneD, NekDouble> tmp2;
            Array<OneD, NekDouble> Vsingle1;
            Array<OneD, NekDouble> Vsingle2;
            Array<OneD, NekDouble> hsingle1;
            Array<OneD, NekDouble> hsingle2;
            ///////////////////////////////////////////////////////////////////////////////
            // // tmp2 for preconditioner multiplication, later consider it

            if(restarted)
            {
                // tmp2=Ax
                m_operator.DoNonlinLinSysLhsEval(pOutput, r0, m_DifferenceFlag0);

                //The first search direction
                beta = -1.0;
                //PYT: r0=b-AX
                Vmath::Svtvp(nNonDir, beta, &r0[0] + nDir, 1, &pInput[0] + nDir, 1, &r0[0] + nDir, 1);

            }
            else
            {
                // If not restarted, x0 should be zero
                Vmath::Vcopy(nNonDir, &pInput[0] + nDir, 1, &r0[0] + nDir, 1);
            }
            
            if(m_flag_LeftPrecond)
            {
                tmp1 = r0 + nDir;
                tmp2 = r0 + nDir;
                m_operator.DoNonlinLinPrecond(tmp1, tmp2);
            }
            
            // norm of (r0)
            // m_map tells how to connect
            vExchange    = Vmath::Dot2(nNonDir,
                                    &r0[0] + nDir,
                                    &r0[0] + nDir,
                                    &m_map[0] + nDir);
            m_Comm->AllReduce(vExchange, LibUtilities::ReduceSum);
            eps = vExchange;

            if(!restarted)
            {
                if(m_prec_factor == NekConstants::kNekUnsetDouble)
                {
                    if(m_flag_LeftPrecond)
                    {
                        // evaluate initial residual error for exit check
                        vExchange    = Vmath::Dot2(nNonDir,
                                                &pInput[0] + nDir,
                                                &pInput[0] + nDir,
                                                &m_map[0] + nDir);
                        m_Comm->AllReduce(vExchange, LibUtilities::ReduceSum);
                        m_prec_factor = vExchange / eps;
                    }
                    else
                    {
                        // cout << "Right precondtioning"<<endl;
                        m_prec_factor = 1.0;
                    }
                }
            }

            tmp2 = r0 + nDir;
            Vmath::Smul(nNonDir,sqrt(m_prec_factor),tmp2,1,tmp2,1);
            eps     =   eps*m_prec_factor;
            eta[0] = sqrt(eps);

            // // If input residual is less than tolerance skip solve.
            // // if (eps * m_prec_factor < m_tolerance * m_tolerance * m_rhs_magnitude)
            // if (eps < m_tolerance * m_tolerance * m_rhs_magnitude)
            // {
            //     m_converged = true;
            //     return eps;
            // }

            // Give an order for the entries in Hessenburg matrix
            for(int nd = 0; nd < m_maxstorage; ++nd)
            {
                id[nd] = nd;
                id_end[nd] = nd + 1;
                starttem = id_end[nd] - m_maxhesband;
                if(truncted && (starttem) > 0)
                {
                    id_start[nd] = starttem;
                }
                else
                {
                    id_start[nd] = 0;
                }
            }

            //Normlized by r0 norm V(:,1)=r0/norm(r0)
            alpha = 1.0 / eta[0];
            //Scalar multiplication
            Vmath::Smul(nNonDir, alpha, &r0[0] + nDir, 1, &V_total[0][0] + nDir, 1);

            // restarted Gmres(m) process
            int nswp = 0;
            if(m_flag_RightPrecond)
            {
                Vsingle1    =   Array<OneD, NekDouble>(nGlobal, 0.0);
            }

            for (int nd = 0; nd < m_maxstorage; ++nd)
            {
                Vsingle2 = V_total[nd + 1];
                hsingle1 = hes[nd];

                if(m_flag_RightPrecond)
                {
                    tmp1 = V_total[nd] + nDir;
                    tmp2 = Vsingle1 + nDir;
                    m_operator.DoNonlinLinPrecond(tmp1, tmp2);
                }
                else
                {
                    Vsingle1 = V_total[nd];
                }
                // w here is no need to add nDir due to temporary Array
                idtem = id[nd];
                starttem = id_start[idtem];
                endtem = id_end[idtem];

                DoArnoldi(starttem, endtem, nGlobal, nDir, V_total, Vsingle1, Vsingle2, hsingle1);

                if(starttem > 0)
                {
                    starttem = starttem - 1;
                }

                hsingle2 = Upper[nd];
                Vmath::Vcopy(m_maxstorage + 1, &hsingle1[0], 1, &hsingle2[0], 1);
                DoGivensRotation(starttem, endtem, nGlobal, nDir, cs, sn, hsingle2, eta);

                eps = eta[nd + 1] * eta[nd + 1];
                // This Gmres merge truncted Gmres to accelerate.
                // If truncted, cannot jump out because the last term of eta is not residual
                if((!truncted) || (nd < m_maxhesband))
                {
                    // if (eps * m_prec_factor < m_tolerance * m_tolerance * m_rhs_magnitude )
                    if ((eps < m_tolerance * m_tolerance * m_rhs_magnitude)&&nd>1 )
                    {
                        m_converged = true;
                    }
                }
                nswp++;
                m_totalIterations++;

                if(m_converged)
                {
                    break;
                }
            }

            DoBackward(nswp, Upper, eta, y_total);
            // calculate output y_total*V_total
            Array<OneD, NekDouble> solution(nGlobal, 0.0);
            for(int i = 0; i < nswp; ++i)
            {
                beta = y_total[i];
                // Vmath::Svtvp(nNonDir, beta, &V_total[i][0] + nDir, 1, &pOutput[0] + nDir, 1, &pOutput[0] + nDir, 1);
                Vmath::Svtvp(nNonDir, beta, &V_total[i][0] + nDir, 1, &solution[0] + nDir, 1, &solution[0] + nDir, 1);
            }

            if(m_flag_RightPrecond)
            {
                tmp1 = solution + nDir;
                tmp2 = solution + nDir;
                m_operator.DoNonlinLinPrecond(tmp1, tmp2);
            }
            Vmath::Vadd(nNonDir, &solution[0] + nDir, 1, &pOutput[0] + nDir, 1, &pOutput[0] + nDir, 1);

            return eps;
        }

        // Arnoldi Subroutine
        void  NekLinSysIteratGMRES::DoArnoldi(
            const int starttem,
            const int endtem,
            const int nGlobal,
            const int nDir,
            // V_total(:,1:nd)
            Array<OneD, Array<OneD,  NekDouble> > &V_local,
            // V[nd]
            Array<OneD, NekDouble> &Vsingle1,
            // V[nd+1]
            Array<OneD, NekDouble> &Vsingle2,
            //h
            Array<OneD, NekDouble> &hsingle)
        {
            // // Get the communicator for performing data exchanges
            // LibUtilities::CommSharedPtr m_Comm
            //     = m_expList.lock()->GetComm()->GetRowComm();

            // To notice, V_local's order not certainly equal to starttem:endtem
            // starttem:endtem is the entry position in Hessenburg matrix
            NekDouble alpha, beta;
            Array<OneD, NekDouble> tmp1, tmp2;
            int numbertem;
            int nNonDir = nGlobal - nDir;
            //later for parallel
            NekDouble vExchange = 0.0;
            // w=AV(:,nd)
            Array<OneD, NekDouble> w(nGlobal, 0.0);

            m_operator.DoNonlinLinSysLhsEval(Vsingle1, w,m_DifferenceFlag1);

            tmp1 = w + nDir;
            tmp2 = w + nDir;
            if(m_flag_LeftPrecond)
            {
                m_operator.DoNonlinLinPrecond(tmp1, tmp2);
            }
            
            Vmath::Smul(nNonDir,sqrt(m_prec_factor),tmp2,1,tmp2,1);

            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Modified Gram-Schmidt
            // The pointer not certainly equal to starttem.
            // Like initially, Gmres-deep need to use numbertem=0
            numbertem = starttem;
            for(int i = starttem; i < endtem; ++i)
            {
                vExchange = Vmath::Dot2(nNonDir,
                                        &w[0] + nDir,
                                        &V_local[numbertem][0] + nDir,
                                        &m_map[0] + nDir);
                m_Comm->AllReduce(vExchange, LibUtilities::ReduceSum);

                hsingle[i] = vExchange;

                beta = -1.0 * vExchange;
                Vmath::Svtvp(nNonDir, beta, &V_local[numbertem][0] + nDir, 1, &w[0] + nDir, 1, &w[0] + nDir, 1);
                numbertem = numbertem + 1;
            }
            // end of Modified Gram-Schmidt
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // // Classical Gram-Schmidt
            // int number = endtem - starttem;
            // Array<OneD,NekDouble> vExchangeArray(number, 0.0);
            // numbertem=0;
            // for (int i = starttem; i < endtem; ++i)
            // {
            //     vExchangeArray[numbertem] =
            //         Vmath::Dot2(nNonDir, &w[0] + nDir, &V_local[i][0] + nDir,
            //                     &m_map[0] + nDir);
            //     numbertem = numbertem + 1;
            // }
            // m_Comm->AllReduce(vExchangeArray,  LibUtilities::ReduceSum);
            // numbertem = 0;
            // for (int i = starttem; i < endtem; ++i)
            // {
            //     hsingle[i] = vExchangeArray[numbertem];
            //     beta       = -1.0 * vExchangeArray[numbertem];
            //     Vmath::Svtvp(nNonDir, beta, &V_local[i][0] + nDir, 1,
            //                  &w[0] + nDir, 1, &w[0] + nDir, 1);
            //     numbertem = numbertem + 1;
            // }
            // // end of Classical Gram-Schmidt
            // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // calculate the L2 norm and normalize
            vExchange = Vmath::Dot2(nNonDir,
                                    &w[0] + nDir,
                                    &w[0] + nDir,
                                    &m_map[0] + nDir);

            m_Comm->AllReduce(vExchange, LibUtilities::ReduceSum);

            hsingle[endtem] = sqrt(vExchange);

            alpha = 1.0 / hsingle[endtem];
            Vmath::Smul(nNonDir, alpha, &w[0] + nDir, 1, &Vsingle2[0] + nDir, 1);
        }

        // QR factorization through Givens rotation
        void  NekLinSysIteratGMRES::DoGivensRotation(
            const int starttem,
            const int endtem,
            const int nGlobal,
            const int nDir,
            Array<OneD, NekDouble> &c,
            Array<OneD, NekDouble> &s,
            Array<OneD, NekDouble> &hsingle,
            Array<OneD, NekDouble> &eta)
        {
            boost::ignore_unused(nGlobal, nDir);
            NekDouble temp_dbl, dd, hh;
            int idtem = endtem - 1;
            // The starttem and endtem are beginning and ending order of Givens rotation
            // They usually equal to the beginning position and ending position of Hessenburg matrix
            // But sometimes starttem will change, like if it is initial 0 and becomes nonzero because previous Givens rotation
            // See Yu Pan's User Guide
            for(int i = starttem; i < idtem; ++i)
            {
                temp_dbl = c[i] * hsingle[i] - s[i] * hsingle[i + 1];
                hsingle[i + 1] = s[i] * hsingle[i] + c[i] * hsingle[i + 1];
                hsingle[i] = temp_dbl;
            }
            dd = hsingle[idtem];
            hh = hsingle[endtem];
            if(hh == 0.0)
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

            hsingle[idtem] = c[idtem] * hsingle[idtem] - s[idtem] * hsingle[endtem];
            hsingle[endtem] = 0.0;

            temp_dbl = c[idtem] * eta[idtem] - s[idtem] * eta[endtem];
            eta[endtem]       = s[idtem] * eta[idtem] + c[idtem] * eta[endtem];
            eta[idtem] = temp_dbl;
        }

        // Backward calculation
        // to notice, Hesssenburg matrix's column and row changes due to use Array<OneD,Array<OneD,NekDouble>>
        void  NekLinSysIteratGMRES::DoBackward(
            const int  number,
            Array<OneD, Array<OneD, NekDouble> > &A,
            const Array<OneD, const NekDouble> &b,
            Array <OneD, NekDouble> &y
        )
        {
            // number is the entry number, but C++'s order need to be one smaller
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
    }
}

