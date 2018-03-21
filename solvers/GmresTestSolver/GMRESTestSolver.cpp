///////////////////////////////////////////////////////////////////////////////
//
// File ADRSolver.cpp
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
// Description: Advection Diffusion Reaction framework solver
//
///////////////////////////////////////////////////////////////////////////////
#include <LibUtilities/BasicUtils/VmathArray.hpp> 
#include <MultiRegions/AssemblyMap/AssemblyMap.h>


#define GMRES_DEBUG

using namespace std;
using namespace Nektar;
//using namespace MultiRegions;

#ifdef GMRES_DEBUG        
            
            NekDouble       m_rhs_magnitude;
            NekDouble       m_tolerance; 
            bool            m_converged;
            bool            m_verbose=true;
            bool            m_root=true;       
            unsigned int    m_maxrestart;
            unsigned int    m_totalIterations;
            unsigned int    m_maxdirction;

            // stores the dimension of the linear sysytem
            unsigned int    m_nlinsys;
            
            // stores the A of the linear system Ax = f
            Array<OneD, NekDouble> m_mat;

            // stores the f of the linear system Ax = f
            Array<OneD,       NekDouble> m_rhs;

#endif

        void v_DoMatrixMultiply(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {
            //ASSERTL1(m_nlinsys > pInput.num_elements()+pInput.GetOffset(),"Array out of bounds");
            //ASSERTL1(m_nlinsys > pInput.num_elements()+pInput.GetOffset(),"Array out of bounds");
            for(int i=0;i<m_nlinsys;++i)
            {
                pOutput[i]  =   Vmath::Dot(m_nlinsys,&m_mat[i*m_nlinsys],&pInput[0]);
                //pOutput[i]  =   sqrt(pOutput[i]);
            }
            return;
        }


        // void initializeLinSys(const LibUtilities::SessionReaderSharedPtr &session)
        void initializeLinSys()
        {
            
            // m_nlinsys                   = session->GetParameter("LinSysDimens");
            // m_maxdirction               = session->GetParameter("MaxDirction" );
            // m_maxrestart                = session->GetParameter("MaxRestart" );
            // m_tolerance                 = session->GetParameter("Tolerance" );

            m_nlinsys                   = 5;
            m_maxdirction               = m_nlinsys-4;
            m_maxrestart                = 1;
            m_tolerance                 = 1.0E-8;
            int nSubMatrix              = 0;
            
            m_rhs       =  Array<OneD, NekDouble>(m_nlinsys);
            m_mat       =  Array<OneD, NekDouble>(m_nlinsys*m_nlinsys);
            Vmath::Fill(m_nlinsys,1.0,&m_rhs[0],1);
            Vmath::Zero(m_nlinsys*m_nlinsys,&m_mat[0],1);

            // // initial the matrix(A) in Ax = f
            // for(int i=0; i<m_nlinsys*m_nlinsys; ++i);
            // {
            //     m_mat[i]    =   1.0+i;
            // }
            // for(int i=0; i<m_nlinsys; ++i)
            // {
            //     m_mat[i*m_nlinsys+i]    =   m_mat[i*m_nlinsys+i]*10.0;
            // }
            
            int       nnn=0;
            for(int i=0; i<m_nlinsys-nSubMatrix; ++i)
            {
                for(int j=0; j<m_nlinsys-nSubMatrix; ++j)
                {
                    m_mat[i*m_nlinsys+j] = 1.0+nnn;
                    nnn++;
                    if (i==j)
                    {
                        m_mat[i*m_nlinsys+j]    =   m_mat[i*m_nlinsys+j]*10.0;
                    }
                }
            }

            for(int i=m_nlinsys-nSubMatrix; i<m_nlinsys; ++i)
            {
                m_mat[i*m_nlinsys+i]    =   1.0;
            }
            
            NekDouble tmp;
            cout <<"The Matrix A is:"<<endl;
            for(int i=0; i<m_nlinsys; ++i)
            {
                cout <<"i="<<i<<"   :";
                for(int j=0; j<m_nlinsys; ++j)
                {
                    tmp =   m_mat[i*m_nlinsys+j];
                    cout<< tmp<<"    ";
                }
                cout <<endl;
            }
            return;
        }

        


        /**  
         * Solve a global linear system(Ax=f, r =f-Ax) using the GMRES method.  
         * We solve only for the non-Dirichlet modes. The operator is evaluated  
         * using an auxiliary function v_DoMatrixMultiply defined by the  
         * specific solver. Distributed math routines are used to support  
         * parallel execution of the solver.  
         *  
         * The implemented algorithm uses a reduced-communication reordering of  
         * the standard PCG method (Demmel, Heath and Vorst, 1993)  
         *  
         * @param       pInput      Input residual(f)  of all DOFs.  
         * @param       pOutput     Solution vector(x) of all DOFs.  
         */
        /**  
         * Solve a global linear system(Ax=f, r =f-Ax) using the conjugate gradient method.  
         * We solve only for the non-Dirichlet modes. The operator is evaluated  
         * using an auxiliary function v_DoMatrixMultiply defined by the  
         * specific solver. Distributed math routines are used to support  
         * parallel execution of the solver.  
         *  
         * The implemented algorithm uses a reduced-communication reordering of  
         * the standard PCG method (Demmel, Heath and Vorst, 1993)  
         *  
         * @param       pInput      Input residual(f)  of all DOFs.  
         * @param       pOutput     Solution vector(x) of all DOFs.  
         */
        NekDouble DoGmresRestart(
            const bool                         rested,
            const int                          nGlobal,
            const Array<OneD,const NekDouble> &pInput,
                  Array<OneD,      NekDouble> &pOutput,
            const int                          nDir)
        {

#ifndef GMRES_DEBUG        
            // Get the communicator for performing data exchanges
            LibUtilities::CommSharedPtr vComm
                = m_expList.lock()->GetComm()->GetRowComm();
#endif
            int nNonDir = nGlobal - nDir;

            Array<TwoD, NekDouble> han    (m_maxdirction+1,m_maxdirction, 0.0);
            Array<OneD, NekDouble> eta    (m_maxdirction+1, 0.0);
            Array<OneD, NekDouble> cs     (m_maxdirction+1, 0.0);
            Array<OneD, NekDouble> sn     (m_maxdirction+1, 0.0);
            Array<OneD, NekDouble> yk     (m_maxdirction+1, 0.0);
            Array<OneD, NekDouble> tmp0;
            Array<OneD, NekDouble> tmp1;

            // Allocate array storage
            int nqk =   (m_maxdirction+1)*nGlobal;
            int iqk0 =0;
            int iqk1 =0;
            Array<OneD, NekDouble> qk_a   (nqk, 0.0);
            Array<OneD, NekDouble> tm_a   (nGlobal, 0.0);
            Array<OneD, NekDouble> tm_b   (nGlobal, 0.0);


            NekDouble alpha, beta, eps, dd, hh,temp_dbl;
            //Array<OneD, NekDouble> vExchange(3,0.0);
            NekDouble   vExchange=0.0;
            Array<OneD, NekDouble> vExchange_a(m_maxdirction,0.0);
            for(int nd=0;nd<m_maxdirction+1;++nd)
            {
                iqk0 = nd*nGlobal;
                Vmath::Zero(nGlobal,&qk_a[iqk0],1);
                //Vmath::Vcopy(nDir,pOutput,1,qk_a[nd],1);
            }            
             
            if(rested)
            {
                // qk_a[0] = A*x0

                iqk0 = 0;
                v_DoMatrixMultiply(pOutput, tm_b);
                tmp0 = qk_a+iqk0;
                for(int k =0;k<nGlobal;++k)
                {
                    tmp0[k] = tm_b[k]; 
                }
            }

            beta = -1.0;
            // q_k[0] = f-A*x0
            iqk0 = 0;
            Vmath::Svtvp(nNonDir, beta, &qk_a[iqk0]+nDir, 1, &pInput[0]+nDir, 1, &qk_a[iqk0]+nDir, 1);

#ifndef GMRES_DEBUG        
            // evaluate initial residual error for exit check
            vExchange    = Vmath::Dot2(nNonDir,
                                       &qk_a[iqk0]+nDir,
                                       &qk_a[iqk0]+nDir,
                                       &m_map[0] + nDir);
            vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);
#else
            // evaluate initial residual error for exit check
            vExchange    = Vmath::Dot(nNonDir,
                                       &qk_a[iqk0]+nDir,
                                       &qk_a[iqk0]+nDir);
#endif            
            eps          = vExchange;
         
            // If input residual is less than tolerance skip solve.
            if (eps < m_tolerance * m_tolerance * m_rhs_magnitude)
            {
                m_converged = true;
                return eps;
            }

            eta[0]       = sqrt(eps);
            alpha        = 1.0/eta[0] ;
            iqk0 = 0;
            Vmath::Smul(nNonDir,alpha,&qk_a[iqk0]+nDir,1,&qk_a[iqk0]+nDir,1);


            int nswp = 0;
            // cout <<"(m_maxdirction-1)="<< (m_maxdirction-1)<<endl;
            // cout <<endl;
            for(int nd=0; nd<(m_maxdirction); ++nd)
            {
                iqk0 = nd*nGlobal;
                iqk1 = (nd+1)*nGlobal;
                
                //tmp1 = qk_a[nd+1][0];
                tmp0 = qk_a+iqk0 + nDir;
                //tmp0 = pInput + nDir;
                tmp1 = tm_a + nDir;
#ifndef GMRES_DEBUG        
                m_precon->DoPreconditioner(tmp0, tmp1);
#else
                for(int k =0;k<nNonDir;++k)
                {
                    tmp1[k] = tmp0[k]; 
                }
#endif

                tmp0 = qk_a+iqk0;
                tmp1 = qk_a+iqk1;
                v_DoMatrixMultiply(tm_a, tm_b);
                for(int k =0;k<nGlobal;++k)
                {
                    tmp1[k] = tm_b[k]; 
                }

#ifdef GMRES_DEBUG        
                cout<< nd<<"th search direction"<<endl;
                // cout<<"debug output"<<endl;
                // for(int jj=0;jj<m_maxdirction; ++jj)
                // {
                //     cout << tmp1[jj]<<endl;
                // }
                
#endif            
                for(int i=0;i<nd+1;++i)
                {
                    // evaluate initial residual error for exit check
                    // tmp0 = qk_a[nDir];
                    int iqki = i*nGlobal;
#ifndef GMRES_DEBUG        
                    vExchange    = Vmath::Dot2(nNonDir,
                                               &qk_a[iqki]+nDir,
                                               &qk_a[iqk1]+nDir,
                                               &m_map[0] + nDir);
#else
                    // evaluate initial residual error for exit check
                    vExchange    = Vmath::Dot(nNonDir,
                                               &qk_a[iqki]+nDir,
                                               &qk_a[iqk1]+nDir);
#endif 
                    han[i][nd] = vExchange;
                }

#ifndef GMRES_DEBUG        
                
                for(int i=0;i<nd+1;++i)
                {
                    vExchange_a[i] = han[i][nd];
                } 
                vComm->AllReduce(vExchange_a, Nektar::LibUtilities::ReduceSum);
                for(int i=0;i<nd+1;++i)
                {
                    han[i][nd] = vExchange_a[i];
                } 
#endif
                
                for(int i=0;i<nd+1;++i)
                {
                    // evaluate initial residual error for exit check
                    // tmp0 = qk_a[nDir];
                    int iqki = i*nGlobal;
                    vExchange   =   han[i][nd];
                    beta = -1.0*vExchange;
                    Vmath::Svtvp(nNonDir, beta, &qk_a[iqki]+nDir, 1, &qk_a[iqk1]+nDir, 1, &qk_a[iqk1]+nDir, 1);
                }


#ifndef GMRES_DEBUG        
                vExchange    = Vmath::Dot2(nNonDir,
                                           &qk_a[iqk1]+nDir,
                                           &qk_a[iqk1]+nDir,
                                           &m_map[0] + nDir);
#else
                // cout<<"debug output"<<endl;
                // for(int jj=0;jj<m_maxdirction; ++jj)
                // {
                //     cout << qk_a[iqk1+jj]<<endl;
                // }
                
                // evaluate initial residual error for exit check
                vExchange    = Vmath::Dot(nNonDir,
                                           &qk_a[iqk1]+nDir,
                                           &qk_a[iqk1]+nDir);
#endif 
#ifndef GMRES_DEBUG        
                vComm->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);
#endif
                han[nd+1][nd] = sqrt(vExchange);
                // q_k[0] = f-A*x0
                alpha        = 1.0/han[nd+1][nd] ;
                Vmath::Smul(nNonDir,alpha,&qk_a[iqk1]+nDir,1,&qk_a[iqk1]+nDir,1);

                for(int i=0;i<nd;++i)
                {
                    temp_dbl        = cs[i]*han[i][nd] - sn[i]*han[i+1][nd];
                    han[i+1][nd]    = sn[i]*han[i][nd] + cs[i]*han[i+1][nd];
                    han[i][nd]      = temp_dbl;
                }
                dd = han[nd][nd];
                hh = han[nd+1][nd];
                if(hh==0.0)
                {
                    cs[nd] = 1.0;
                    sn[nd] = 0.0;
                }
                else if (abs(hh) > abs(dd))
                {
                    temp_dbl = -dd/hh;
                    sn[nd] = 1.0 / sqrt(1.0 + temp_dbl*temp_dbl);
                    cs[nd] = temp_dbl * sn[nd];
                }
                else
                {
                    temp_dbl = -hh/dd;
                    cs[nd] = 1.0 / sqrt(1.0 + temp_dbl*temp_dbl);
                    sn[nd] = temp_dbl * cs[nd];
                }

                temp_dbl = sn[nd]*han[nd][nd]+cs[nd]*han[nd+1][nd];

                han[nd][nd] = cs[nd]*han[nd][nd]-sn[nd]*han[nd+1][nd];
                han[nd+1][nd] = 0.0;

                //dd = cs[nd]*eta[nd] ;
                //eta[nd+1] = -sn[nd]*eta[nd] ;

                temp_dbl        = cs[nd]*eta[nd] - sn[nd]*eta[nd+1];
                eta[nd+1]       = sn[nd]*eta[nd] + cs[nd]*eta[nd+1];
                eta[nd]         = temp_dbl;



                eps          = eta[nd+1]*eta[nd+1];

                nswp++;
                m_totalIterations++;


                // If input residual is less than tolerance skip solve.
                if (eps < m_tolerance * m_tolerance * m_rhs_magnitude)
                {
                    m_converged = true;
                    break;
                }
            }
#ifdef GMRES_DEBUG        
            for(int i=0;i<m_maxdirction+1;++i)
            {
                cout << "eta[" <<i<<"]="<< eta[i]<<endl;
            }
            for(int i=0;i<m_maxdirction+1;++i)
            {
                for(int j=0;j<m_maxdirction;++j)
                {
                    cout << "han[" <<i<<"][" <<j<<"]="<< han[i][j]<<"   ";
                }

                cout << endl;
            }
#endif
            for(int i=0;i<nswp;++i)
            {
                yk[i] = eta[i];
            }

            for(int i=nswp-1;i>-1;--i)
            {
                //cout <<"yk"<<i<<"="<<yk[i]<<endl;
                //cout <<"han"<<i<<"="<<han[i][i]<<endl;
                yk[i] = yk[i]/han[i][i];
                //cout <<"yk"<<i<<"="<<yk[i]<<endl;
                for(int j=0;j<i;++j)
                {
                    yk[j] = yk[j]-han[j][i]*yk[i];
                }
                //cout <<"yk"<<i<<"="<<yk[i]<<endl;
            }



            //tmp1 = qk_a[nDir][0];
            for(int i=0;i<nswp;++i)
            {
                // q_k[0] = f-A*x0
                beta = yk[i];
                //cout << beta <<endl;
                int iqki = i*nGlobal;
                //tmp0 = qk_a[nDir][0];
                Vmath::Svtvp(nNonDir, beta, &qk_a[iqki]+nDir, 1, &pOutput[0]+nDir, 1, &pOutput[0]+nDir, 1);
            }
#ifdef GMRES_DEBUG        
            for(int i=0;i<m_maxdirction+1;++i)
            {
                cout << "yk[" <<i<<"]="<< yk[i]<<endl;
            }
#endif
            return eps;
            
        }


        void DoGMRES(
            const int                          nGlobal,
            const Array<OneD,const NekDouble> &pInput,
                  Array<OneD,      NekDouble> &pOutput,
            const int                          nDir)
        {
#ifndef GMRES_DEBUG        
            if (!m_precon)
            {
                v_UniqueMap();
                m_precon = CreatePrecon(plocToGloMap);
                m_precon->BuildPreconditioner();
            }
#endif

            /*
            // Get the communicator for performing data exchanges
            LibUtilities::CommSharedPtr vComm
                = m_expList.lock()->GetComm()->GetRowComm();
            */
            
            // Get vector sizes
            NekDouble eps = 0.0;
            int nNonDir = nGlobal - nDir;

            // zero homogeneous out array ready for solution updates
            // Should not be earlier in case input vector is same as
            // output and above copy has been peformed
            Vmath::Zero(nNonDir, &pOutput[nDir],1);
            //Vmath::Zero(nGlobal, &qk_a[iqk0],1);


#ifndef GMRES_DEBUG        
            if(m_rhs_magnitude == NekConstants::kNekUnsetDouble)
            {
                NekVector<NekDouble> inGlob (nGlobal, pInput, eWrapper);
                Set_Rhs_Magnitude(inGlob);
            }
#else
            m_rhs_magnitude = 1.0;
#endif

            m_totalIterations = 0;
            m_converged       = false;
            
            bool restarted = false;
            for(int nrestart=0;nrestart<m_maxrestart;++nrestart)
            {
                eps = DoGmresRestart(restarted, nGlobal,pInput,pOutput,nDir);
   
                if(m_converged)
                {   
                    cout << "Solution converged!!!"<<endl;
                    if (m_verbose && m_root)
                    {
                        cout << "GMRES iterations made = " << m_totalIterations 
                             << " using tolerance of "  << m_tolerance 
                             << " (error = " << sqrt(eps/m_rhs_magnitude) 
                             << ", rhs_mag = " << sqrt(m_rhs_magnitude) <<  ")" 
                             << endl;
                    }
                    return;
                }
                restarted = true;
            }
            
            if(m_root)
            {
                cout << "GMRES iterations made = " << m_totalIterations 
                     << " using tolerance of "  << m_tolerance 
                     << " (error = " << sqrt(eps/m_rhs_magnitude)
                     << ", rhs_mag = " << sqrt(m_rhs_magnitude) <<  ")"
                     << endl;
            }
#ifndef GMRES_DEBUG        
            ROOTONLY_NEKERROR(ErrorUtil::efatal,
                              "Exceeded maximum number of iterations");
#endif

            return;
        }


int main(int argc, char *argv[])
{
    // LibUtilities::SessionReaderSharedPtr pppsession;
    
        // Create session reader.
        // pppsession = LibUtilities::SessionReader::CreateInstance(argc, argv);

        // // Get some information about the session
        // string       sessionName    = pppsession->GetSessionName();
        // string       outFile        = sessionName + ".fld";
        


        initializeLinSys();

        int ndim = m_nlinsys;

        Array<OneD, NekDouble> soltn_cg(ndim);

        //DoConjugateGradient(ndim, , soltn_cg, NoUsePtr, 0);

        Array<OneD, NekDouble> soltn_gmres(ndim);

        DoGMRES(ndim, m_rhs, soltn_gmres, 0);

#ifdef GMRES_DEBUG        
        cout<<"the gmres solution is :"<<endl;
        for(int jj=0;jj<ndim; ++jj)
        {
            cout << soltn_gmres[jj]<<endl;
        }
#endif   

        NekDouble l2error=0.0;
        NekDouble tmp;
        for(int i; i< ndim; ++i)
        {   
            tmp = soltn_cg[i]-soltn_gmres[i];
            l2error += tmp*tmp;
        }

        cout << "The L2 difference between CG&GMRES is"<<l2error<<endl;

        // Finalise session
        // pppsession->Finalise();


    return 0;
}
